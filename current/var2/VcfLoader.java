package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ListNum;

/**
 * Loads VAR records or VcfToVar-compatible VCF into a newly allocated VarMap.
 * A producer handles headers and queues batches of up to200 data lines; parser
 * workers merge equivalent variants under the result map's monitor. Parsing/merging
 * order need not follow input order. VCF requires a supplied scaffold registry and the custom
 * INFO blocks requested by the flags; this is not a general VCF-normalization layer.
 *
 * Static load methods join all workers and throw with the first recorded failure.
 * Queue waits periodically check shared failure state so a failed parser cannot
 * strand the producer or other parsers (VAR2-029/030). Reader closure runs on all
 * producer exits and its error flag is propagated (VAR2-031). Blank lines are ignored.
 * The producer also changes global ByteFile mode selection without restoring it.
 * @author Brian Bushnell
 */
public class VcfLoader{

	/**
	 * Allocates result and queue state without starting the loader workers.
	 * Captures max(1, Shared.threads()) parser workers and the same queue capacity;
	 * DEFAULT_THREADS is not consulted. The input descriptor defaults to TXT while
	 * vcfMode selects record parsing. The scaffold registry is borrowed, not populated
	 * from contig headers by this loader.
	 * @param fname_ Input path
	 * @param scafMap_ Registry required by VCF parsing; VAR records carry numeric scaffold IDs
	 * @param vcfMode_ True for VcfToVar records, false for tab-delimited VAR
	 */
	public VcfLoader(String fname_, ScafMap scafMap_, boolean vcfMode_){
		fname=fname_;
		scafMap=scafMap_;
		varMap=new VarMap(scafMap);
		threads=Tools.max(1, Shared.threads());
		inq=new ArrayBlockingQueue<ListNum<byte[]>>(threads);
		vcfMode=vcfMode_;
		ffin=FileFormat.testInput(fname, FileFormat.TXT, null, true, false);
	}

	/**
	 * Dispatches VAR descriptors to loadVarFile and all other descriptors to loadVcfFile.
	 * @param ff Non-null input descriptor; only its name and var flag are used here
	 * @param scafMap Scaffold registry for VCF parsing
	 * @param loadCoverage Load the VCF custom read-count/coverage block; ignored for VAR
	 * @param extendedInfo Load the VCF extended-statistics block; ignored for VAR
	 * @return Newly allocated map, after successful worker completion
	 */
	public static VarMap loadFile(FileFormat ff, ScafMap scafMap, boolean loadCoverage, boolean extendedInfo){
		final VarMap varMap;
		if(ff.var()){
			varMap=VcfLoader.loadVarFile(ff.name(), scafMap);
		}else{
			varMap=VcfLoader.loadVcfFile(ff.name(), scafMap, loadCoverage, extendedInfo);
		}
		return varMap;
	}

	/**
	 * Starts VAR loading and waits for all workers before returning their map.
	 * Headers may populate dataset statistics; equivalent variants merge read evidence.
	 * @param fname Input path
	 * @param scafMap Borrowed map associated with the result; may be null for numeric VAR records
	 * @return Result map after all workers terminate
	 * @throws RuntimeException if parsing, reading, closure or worker completion records an error
	 */
	public static VarMap loadVarFile(String fname, ScafMap scafMap){
		VcfLoader loader=new VcfLoader(fname, scafMap, false);
		ArrayList<ProcessThread> alpt=loader.spawnThreads(false, false);
		loader.waitForFinish(alpt);
		if(loader.errorState){throw new RuntimeException("VcfLoader: failed to load "+fname, loader.failure);}
		return loader.varMap;
	}

	/**
	 * Starts VCF loading and waits for all workers. Uses VcfToVar's single-ALT and
	 * historical coordinate conventions. Scalar dataset headers may populate the map;
	 * contig headers are retained internally but do not construct the scaffold registry.
	 * @param fname Input path
	 * @param scafMap Non-null registry containing the records' scaffold names
	 * @param loadCoverage Require/parse the custom coverage block
	 * @param extendedInfo Require/parse the extended-statistics block
	 * @return Result map after successful worker completion
	 * @throws RuntimeException for recorded worker failures after all workers terminate
	 */
	public static VarMap loadVcfFile(String fname, ScafMap scafMap, boolean loadCoverage, boolean extendedInfo){
		VcfLoader loader=new VcfLoader(fname, scafMap, true);
		ArrayList<ProcessThread> alpt=loader.spawnThreads(loadCoverage, extendedInfo);
		loader.waitForFinish(alpt);
		if(loader.errorState){throw new RuntimeException("VcfLoader: failed to load "+fname, loader.failure);}
		return loader.varMap;
	}

	/**
	 * Constructs the complete thread list, then starts one producer and the captured
	 * number of parser workers. The producer is tid0. This method does not join them.
	 * @param loadCoverage Request the VCF coverage block; ignored for VAR records
	 * @param extendedInfo Request extended VCF fields; ignored for VAR records
	 * @return List of started workers, producer first
	 */
	private ArrayList<ProcessThread> spawnThreads(boolean loadCoverage, boolean extendedInfo){

		//Do anything necessary prior to processing

		//Determine how many threads may be used
		final int threads=this.threads+1;//1 extra for thread 0, which is different than the others.
		assert(threads>=2);

		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, alpt, loadCoverage, extendedInfo));
		}
		if(verbose){outstream.println("Spawned threads.");}

		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}

		//Do anything necessary after processing
		return alpt;
	}

	/**
	 * Joins every worker. Interruption records failure to release queue waits, then
	 * joining continues until resources are released; interrupt status is restored.
	 * Sets errorState if any worker failed or a shared failure was recorded.
	 * @param alpt Worker list from spawnThreads
	 */
	private void waitForFinish(ArrayList<ProcessThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true, interrupted=false;
		for(ProcessThread pt : alpt){
			if(verbose){outstream.println("Waiting for thread "+pt.tid);}
			while(pt.getState()!=Thread.State.TERMINATED){
				try{
					//Attempt a join operation
					pt.join();
				}catch(InterruptedException e){
					interrupted=true;
					recordFailure(e);
				}
			}

			allSuccess&=pt.success;
		}

		//Track whether any threads failed
		if(!allSuccess || failure!=null){errorState=true;}
		if(interrupted){Thread.currentThread().interrupt();}
	}

	/** Publishes the first failure; queue waiters observe it without taking this monitor. */
	private synchronized void recordFailure(Throwable cause){
		assert(cause!=null) : "A nonnull cause releases failed loader queue waits (VAR2-030)";
		if(failure==null){failure=cause;}
	}

	/**
	 * Parses numeric-ID VAR data through Var's tab-delimited constructor.
	 * For headers, splits tab fields with the default charset and asserts at least two
	 * fields. Recognizes ploidy, pairingRate, totalQualityAvg, mapqAvg and readLengthAvg
	 * case-insensitively, updating the result under its monitor. Retains each header's
	 * original byte array; unrecognized headers do not alter statistics.
	 * @param line Raw line, possibly null or empty
	 * @return New variant, or null for a header/null/empty line
	 */
	private Var loadVarLine(byte[] line){
		if(line==null || line.length<1){return null;}
		if(line[0]!='#'){
			return new Var(line, (byte)'\t');
		}else{
			String[] split=new String(line).split("\t");
			String a=split[0], b=(split.length>1 ? split[1] : null);
			assert(split.length>1) : new String(line);
			if(a.equalsIgnoreCase("#ploidy")){
				synchronized(varMap){varMap.ploidy=Integer.parseInt(b);}
			}else if(a.equalsIgnoreCase("#pairingRate")){
				synchronized(varMap){varMap.properPairRate=Double.parseDouble(b);}
			}else if(a.equalsIgnoreCase("#totalQualityAvg")){
				synchronized(varMap){varMap.totalQualityAvg=Double.parseDouble(b);}
			}else if(a.equalsIgnoreCase("#mapqAvg")){
				synchronized(varMap){varMap.totalMapqAvg=Double.parseDouble(b);}
			}else if(a.equalsIgnoreCase("#readLengthAvg")){
				synchronized(varMap){varMap.readLengthAvg=Double.parseDouble(b);}
			}
			synchronized(header){header.add(line);}
			return null;
		}
	}

	/**
	 * Parses data through VcfToVar; malformed data throws with the original cause.
	 * Worker failure handling publishes Exceptions and Errors to the load caller.
	 * Headers split on '=' using the default charset; exactly two parts are needed to
	 * recognize ploidy, properPairRate, totalQualityAvg, mapqAvg or readLengthAvg.
	 * Retains header bytes internally; it does not interpret contig definitions.
	 * @param line Raw line, possibly null or empty
	 * @param loadCoverage Require/parse custom coverage fields on data records
	 * @param loadExtended Require/parse custom extended fields on data records
	 * @return New variant, or null for header/empty input
	 */
	private Var loadVcfLine(byte[] line, boolean loadCoverage, boolean loadExtended){
		if(line==null || line.length<1){return null;}
		if(line[0]!='#'){
			try{
				return VcfToVar.fromVCF(line, scafMap, loadCoverage, loadExtended);
			}catch(Exception e){
				throw new RuntimeException("Unable to parse VCF line: '"+new String(line)+"'", e);
			}
		}else{
			String[] split=new String(line).split("=");
			if(split.length==2){
				String a=split[0], b=split[1];
				if(a.equalsIgnoreCase("##ploidy")){
					synchronized(varMap){varMap.ploidy=Integer.parseInt(b);}
				}else if(a.equalsIgnoreCase("##properPairRate")){
					synchronized(varMap){varMap.properPairRate=Double.parseDouble(b);}
				}else if(a.equalsIgnoreCase("##totalQualityAvg")){
					synchronized(varMap){varMap.totalQualityAvg=Double.parseDouble(b);}
				}else if(a.equalsIgnoreCase("##mapqAvg")){
					synchronized(varMap){varMap.totalMapqAvg=Double.parseDouble(b);}
				}else if(a.equalsIgnoreCase("##readLengthAvg")){
					synchronized(varMap){varMap.readLengthAvg=Double.parseDouble(b);}
				}
			}
			synchronized(header){header.add(line);}
			return null;
		}
	}

	/**
	 * Producer for tid0 or parser for a positive tid, sharing the bounded input queue.
	 * Normal completion forwards the shared terminal. A failed role publishes its
	 * cause instead; producer and parser queue operations then stop waiting.
	 */
	private class ProcessThread extends Thread{

		/**
		 * Stores the role and parsing flags. Only the producer retains alpt, which is
		 * not otherwise used by the current worker implementation.
		 * @param tid_ Zero for producer, positive for parser
		 * @param alpt_ Complete worker-list object assembled before starting threads
		 * @param loadCoverage_ Request VCF coverage parsing
		 * @param extendedInfo_ Request extended VCF parsing
		 */
		ProcessThread(int tid_, ArrayList<ProcessThread> alpt_, boolean loadCoverage_, boolean extendedInfo_){
			tid=tid_;
			alpt=(tid==0 ? alpt_ : null);
			loadCoverage=loadCoverage_;
			extendedInfo=extendedInfo_;
		}

		/** Runs one role, forwards normal completion, and publishes any failure. */
		@Override
		public void run(){
			try{
				if(tid==0){processBytes();}
				else{makeVars();}
				putBytes(POISON_BYTES);
			}catch(Throwable t){
				recordFailure(t);
			}
			success=(failure==null);
		}

		/** Runs the producer; success is determined after terminal forwarding. */
		void processBytes(){processBytes0();}

		/**
		 * Opens a ByteFile, handles headers on this producer and enqueues data batches.
		 * Sets global FORCE_MODE_BF2 without restoring it; higher-priority force flags can
		 * still select another reader. Batch IDs increase from zero but do not order merges.
		 * Ignores blank lines before indexing (VAR2-032). Closes the reader on every exit,
		 * preserving a primary failure if closure also fails (VAR2-031).
		 */
		public final void processBytes0(){
			if(verbose){outstream.println("tid "+tid+" started processBytes.");}

			ByteFile.FORCE_MODE_BF2=true;
			ByteFile bf=ByteFile.makeByteFile(ffin);

			try{
				long number=0;
				ArrayList<byte[]> list=new ArrayList<byte[]>(LIST_SIZE);
				for(byte[] line=bf.nextLine(); line!=null && failure==null; line=bf.nextLine()){
					if(line.length==0){continue;}
					if(line[0]=='#'){
						if(vcfMode){loadVcfLine(line, false, false);}
						else{loadVarLine(line);}
					}else{
						list.add(line);
						if(list.size()>=LIST_SIZE){
							if(!putBytes(new ListNum<byte[]>(list, number++))){return;}
							list=new ArrayList<byte[]>(LIST_SIZE);
						}
					}
				}
				if(list.size()>0 && failure==null){putBytes(new ListNum<byte[]>(list, number));}
			}catch(Throwable t){
				recordFailure(t);
			}finally{
				if(bf.close()){recordFailure(new RuntimeException("VcfLoader: read/close failed for "+fname));}
			}
			if(verbose){outstream.println("tid "+tid+" closed stream.");}
		}

		/**
		 * Transfers a batch without copying it, or returns false after shared failure.
		 * Callers must stop mutating transferred data. Interruptions fail the load.
		 * Periodic failure checks also release a full queue after its consumers fail.
		 * @param list Nonnull data batch or the shared terminal identity
		 * @return True if the batch was enqueued
		 */
		final boolean putBytes(ListNum<byte[]> list){
			assert(list!=null) : "Only real batches or POISON_BYTES are valid queue handoffs";
			while(failure==null){
				try{
					if(inq.offer(list, 100, TimeUnit.MILLISECONDS)){return true;}
				}catch(InterruptedException e){
					recordFailure(e);
					Thread.currentThread().interrupt();
					return false;
				}
			}
			return false;
		}

		/**
		 * Waits for data or normal completion, checking shared failure between waits.
		 * @return Dequeued batch, or the terminal identity when loading has failed
		 */
		final ListNum<byte[]> takeBytes(){
			while(failure==null){
				try{
					final ListNum<byte[]> list=inq.poll(100, TimeUnit.MILLISECONDS);
					if(list!=null){return list;}
				}catch(InterruptedException e){
					recordFailure(e);
					Thread.currentThread().interrupt();
					return POISON_BYTES;
				}
			}
			return POISON_BYTES;
		}

		/**
		 * Parses batches until the exact shared poison object arrives. Accumulates each
		 * batch's parsed variants, then locks varMap while merging non-null results.
		 * An uncaught parse failure can discard earlier parsed results from the same batch.
		 * Equivalent variants merge evidence; completion order need not match file order.
		 * Normal producer routing excludes headers. Shared failure stops further batches.
		 */
		void makeVars(){
			if(verbose){outstream.println("tid "+tid+" started makeVars.");}

			ListNum<byte[]> list=takeBytes();
			ArrayList<Var> vars=new ArrayList<Var>(LIST_SIZE);
			while(list!=POISON_BYTES && failure==null){
				vars.clear();
				if(vcfMode){
					for(byte[] line : list){
						assert(line[0]!='#') : new String(line);
						Var v=loadVcfLine(line, loadCoverage, extendedInfo);
						vars.add(v);
					}
				}else{
					for(byte[] line : list){
						assert(line[0]!='#') : new String(line);
						Var v=loadVarLine(line);
						vars.add(v);
					}
				}
				synchronized(varMap){
					for(Var v : vars){
						if(v!=null){varMap.addUnsynchronized(v);}
					}
				}
				list=takeBytes();
			}
			if(verbose){outstream.println("tid "+tid+" done making vars.");}

			//Normal terminal forwarding is handled by run(); failures cancel queue waits.
			if(verbose){outstream.println("tid "+tid+" done poisoning bytes.");}
		}

		/*--------------------------------------------------------------*/
		/*----------------           Fields             ----------------*/
		/*--------------------------------------------------------------*/

		/** Worker-list reference retained only by producer; currently unused after construction. */
		final ArrayList<ProcessThread> alpt;
		/** Role selector: zero reads bytes, positive values parse. */
		final int tid;
		/** Whether this parser requests the VCF custom coverage block. */
		final boolean loadCoverage;
		/** Whether this parser requests VCF extended statistics. */
		final boolean extendedInfo;
		/** Role and terminal handoff finished without observing failure; checked after join. */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Input name retained from construction. */
	final String fname;
	/** Input descriptor; record parsing is controlled separately by vcfMode. */
	final FileFormat ffin;

	/** Bounded shared queue; capacity equals the captured parser-worker count. */
	final ArrayBlockingQueue<ListNum<byte[]>> inq;

	/** Captured number of parser workers; producer is additional. */
	final int threads;

	/** Original header arrays retained under the list monitor; not exposed by static helpers. */
	ArrayList<byte[]> header=new ArrayList<byte[]>();

	/** Borrowed status stream; this loader does not close it. */
	protected PrintStream outstream=System.err;

	/** Borrowed VCF scaffold registry; not built from headers here. */
	final ScafMap scafMap;
	/** New result map; worker merges and header updates lock this object. */
	final VarMap varMap;
	/** Select VcfToVar parsing instead of numeric-ID VAR parsing. */
	final boolean vcfMode;

	/** Final status set after all worker joins; failed loads never return their partial map. */
	boolean errorState=false;

	/** First failure, published to producer/parser queue waits and retained as the caller cause. */
	private volatile Throwable failure;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared terminal identity, forwarded through the same bounded queue. */
	static final ListNum<byte[]> POISON_BYTES=new ListNum<byte[]>(null, Long.MAX_VALUE, true, false);
	/** Maximum number of data lines per producer batch. */
	static final int LIST_SIZE=200;
	/** Legacy setting unused by this loader; construction uses Shared.threads(). */
	public static int DEFAULT_THREADS=3;
	/** Diagnostic logging switch; not a cancellation or progress mechanism. */
	static boolean verbose=false;
}
