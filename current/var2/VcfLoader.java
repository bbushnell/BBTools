package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ListNum;

/**
 * Multithreaded loader for VCF and VAR format files.
 * Uses producer-consumer pattern with one thread reading lines and
 * multiple worker threads parsing variants to maximize I/O and CPU efficiency.
 *
 * Supports both custom VAR format and standard VCF format with optional
 * coverage and extended information parsing.
 *
 * @author Brian Bushnell
 */
public class VcfLoader {
	
	public VcfLoader(String fname_, ScafMap scafMap_, boolean vcfMode_){
		fname=fname_;
		scafMap=scafMap_;
		varMap=new VarMap(scafMap);
		threads=Tools.max(1, Shared.threads());
		inq=new ArrayBlockingQueue<ListNum<byte[]>>(threads);
		vcfMode=vcfMode_;
		ffin=FileFormat.testInput(fname, FileFormat.TXT, null, true, false);
	}
	
	public static VarMap loadFile(FileFormat ff, ScafMap scafMap, boolean loadCoverage, boolean extendedInfo){
		final VarMap varMap;
		if(ff.var()){
			varMap=VcfLoader.loadVarFile(ff.name(), scafMap);
		}else{
			varMap=VcfLoader.loadVcfFile(ff.name(), scafMap, loadCoverage, extendedInfo);
		}
		return varMap;
	}
	
	public static VarMap loadVarFile(String fname, ScafMap scafMap){
		VcfLoader loader=new VcfLoader(fname, scafMap, false);
		ArrayList<ProcessThread> alpt=loader.spawnThreads(false, false);
		loader.waitForFinish(alpt);
		return loader.varMap;
	}
	
	public static VarMap loadVcfFile(String fname, ScafMap scafMap, boolean loadCoverage, boolean extendedInfo){
		VcfLoader loader=new VcfLoader(fname, scafMap, true);
		ArrayList<ProcessThread> alpt=loader.spawnThreads(loadCoverage, extendedInfo);
		loader.waitForFinish(alpt);
		return loader.varMap;
	}
	
	/**
	 * Spawns producer and consumer threads for parallel file processing.
	 * Creates one reader thread (tid=0) and multiple parser threads.
	 *
	 * @param loadCoverage Whether to parse coverage information
	 * @param extendedInfo Whether to parse extended statistical fields
	 * @return List of spawned ProcessThread objects
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
	
	private void waitForFinish(ArrayList<ProcessThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true;
		for(ProcessThread pt : alpt){
			if(verbose){outstream.println("Waiting for thread "+pt.tid);}
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}

			allSuccess&=pt.success;
		}
		
		//Track whether any threads failed
		if(!allSuccess){errorState=true;}
	}
	
	/**
	 * Parses a single line from VAR format file.
	 * Handles both data lines (variants) and header lines (metadata).
	 * @param line Raw line bytes
	 * @return Parsed Var object or null for header lines
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
	 * Parses a single line from VCF format file.
	 * Handles both data lines (variants) and header lines (metadata).
	 *
	 * @param line Raw line bytes
	 * @param loadCoverage Whether to parse coverage information
	 * @param loadExtended Whether to parse extended statistical fields
	 * @return Parsed Var object or null for header lines
	 */
	private Var loadVcfLine(byte[] line, boolean loadCoverage, boolean loadExtended){
		if(line==null || line.length<1){return null;}
		if(line[0]!='#'){
			try {
				return VcfToVar.fromVCF(line, scafMap, loadCoverage, loadExtended);
			} catch (Exception e) {
				System.err.println("Unable to parse VCF line: '"+new String(line)+"'");
				return null;
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
	
	/** Worker thread that either reads file bytes (tid=0) or parses variants (tid>0).
	 * Implements producer-consumer pattern using ArrayBlockingQueue for coordination. */
	private class ProcessThread extends Thread{
		
		ProcessThread(int tid_, ArrayList<ProcessThread> alpt_, boolean loadCoverage_, boolean extendedInfo_){
			tid=tid_;
			alpt=(tid==0 ? alpt_ : null);
			loadCoverage=loadCoverage_;
			extendedInfo=extendedInfo_;
		}
		
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			if(tid==0){
				processBytes();
			}else{
				makeVars();
			}
			
			//Indicate successful exit status
			success=true;
		}
		
		void processBytes(){
			processBytes0();
			
			success=true;
		}
		
		/** Reads file line-by-line and distributes batches to parser threads.
		 * Handles header lines directly and queues data lines for parsing. */
		public final void processBytes0(){
			if(verbose){outstream.println("tid "+tid+" started processBytes.");}

			ByteFile.FORCE_MODE_BF2=true;
			ByteFile bf=ByteFile.makeByteFile(ffin);
			
			long number=0;
			
			ArrayList<byte[]> list=new ArrayList<byte[]>(LIST_SIZE);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				assert(line!=null);
				if(line[0]=='#'){
					if(vcfMode){loadVcfLine(line, false, false);}
					else{loadVarLine(line);}
				}else{
					list.add(line);
					if(list.size()>=LIST_SIZE){
						putBytes(new ListNum<byte[]>(list, number));
						number++;
						list=new ArrayList<byte[]>(LIST_SIZE);
					}
				}
			}
			if(verbose){outstream.println("tid "+tid+" ran out of input.");}
			if(list.size()>0){
				putBytes(new ListNum<byte[]>(list, number));
				number++;
				list=null;
			}
			if(verbose){outstream.println("tid "+tid+" done reading bytes.");}
			putBytes(POISON_BYTES);
			if(verbose){outstream.println("tid "+tid+" done poisoning.");}
			bf.close();
			if(verbose){outstream.println("tid "+tid+" closed stream.");}
		}
		
		/**
		 * Adds a batch of lines to the processing queue.
		 * Blocks if queue is full.
		 * @param list Batch of lines to queue
		 */
		final void putBytes(ListNum<byte[]> list){
			if(verbose){outstream.println("tid "+tid+" putting blist size "+list.size());}
			while(list!=null){
				try {
					inq.put(list);
					list=null;
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("tid "+tid+" done putting blist");}
		}
		
		/**
		 * Takes a batch of lines from the processing queue.
		 * Blocks if queue is empty.
		 * @return Batch of lines to process
		 */
		final ListNum<byte[]> takeBytes(){
			if(verbose){outstream.println("tid "+tid+" taking blist");}
			ListNum<byte[]> list=null;
			while(list==null){
				try {
					list=inq.take();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("tid "+tid+" took blist size "+list.size());}
			return list;
		}
		
		/** Parser thread main loop - processes batches of variant lines.
		 * Continues until poison pill is received, then forwards it to next parser. */
		void makeVars(){
			if(verbose){outstream.println("tid "+tid+" started makeVars.");}
			
			ListNum<byte[]> list=takeBytes();
			ArrayList<Var> vars=new ArrayList<Var>(LIST_SIZE);
			while(list!=POISON_BYTES){
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
						varMap.addUnsynchronized(v);
					}
				}
				list=takeBytes();
			}
			if(verbose){outstream.println("tid "+tid+" done making vars.");}

			putBytes(POISON_BYTES);
			if(verbose){outstream.println("tid "+tid+" done poisoning bytes.");}
		}
		
		/*--------------------------------------------------------------*/
		/*----------------           Fields             ----------------*/
		/*--------------------------------------------------------------*/

		final ArrayList<ProcessThread> alpt;
		final int tid;
		final boolean loadCoverage;
		final boolean extendedInfo;
		boolean success=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final String fname;
	final FileFormat ffin;

	final ArrayBlockingQueue<ListNum<byte[]>> inq;
	
	final int threads;
	
	ArrayList<byte[]> header=new ArrayList<byte[]>();
	
	protected PrintStream outstream=System.err;
	
	final ScafMap scafMap;
	final VarMap varMap;
	final boolean vcfMode;
	
	boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	static final ListNum<byte[]> POISON_BYTES=new ListNum<byte[]>(null, Long.MAX_VALUE, true, false);
	static final int LIST_SIZE=200;
	public static int DEFAULT_THREADS=3;
	static boolean verbose=false;
}