package tax;

import java.io.IOException;
import java.io.PrintStream;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import cardinality.DynamicDemiLog;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Builds disk-backed binary tables for the taxonomy server.
 * Reads NCBI accession-to-taxid files and GI tables, writes
 * memory-mappable binary files.
 *
 * Usage:
 *   BuildDiskTables accession=auto pattern=auto table=auto tree=null \
 *     accout=accession_disk.bin giout=gi_disk.bin
 *
 * @author Brian Bushnell, Chloe
 */
public class BuildDiskTables{

	/** Runs disk-table construction from command-line arguments.
	 * @param args BBTools-style arguments selecting input and output tables */
	public static void main(String[] args){
		Timer t=new Timer();
		BuildDiskTables x=new BuildDiskTables(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/** Parses input-table, output-table and buffering options.
	 * @param args BBTools-style command-line arguments */
	public BuildDiskTables(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_UNPIGZ=true;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			final String rawValue=split.length>1 ? split[1] : null;
			final String b=rawValue==null || rawValue.isEmpty() ? null : rawValue;

			if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=parseTreeFile(rawValue);
			}else if(a.equals("pattern")){
				patternFile=b;
			}else if(a.equals("accout") || a.equals("accessionout")){
				accOut=b;
			}else if(a.equals("giout")){
				giOut=b;
			}else if(a.equals("bufferentries")){
				bufferEntries=Parse.parseIntKMG(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		if("auto".equalsIgnoreCase(patternFile)){patternFile=TaxTree.defaultPatternFile();}
	}

	/** Builds requested disk tables after loading optional pattern and tree inputs.
	 * @param t Timer for whole-program reporting */
	void process(Timer t){
		if(taxTreeFile!=null){
			outstream.println("Loading tax tree.");
			tree=TaxTree.loadTaxTree(taxTreeFile, outstream, false, false);
		}

		if(patternFile!=null){
			outstream.println("Loading pattern table.");
			AnalyzeAccession.loadCodeMap(patternFile);
		}

		if(giTableFile!=null && giOut!=null){
			buildGiTable();
		}

		if(accessionFile!=null && accOut!=null){
			if(AnalyzeAccession.codeMap==null){
				throw new RuntimeException("Building a disk accession table requires pattern=...");
			}
			buildAccessionTable();
		}

		t.stop();
		outstream.println("Total time: "+t);
	}

	/*--------------------------------------------------------------*/

	/** Loads the configured GI source table and writes the memory-mappable image. */
	private void buildGiTable(){
		Timer t=new Timer();
		outstream.println("\n--- Building GI disk table ---");

		GiToTaxid.initialize(giTableFile);

		try{
			DiskGiTable.build(giOut, GiToTaxid.array, GiToTaxid.maxGiLoaded);
		}catch(IOException e){
			throw new RuntimeException(e);
		}

		t.stop();
		outstream.println("GI table built in "+t);
		GiToTaxid.unload();
	}

	/*--------------------------------------------------------------*/

	/** Counts digitizable accession keys, then builds the disk accession table. */
	private void buildAccessionTable(){
		Timer t=new Timer();
		outstream.println("\n--- Building accession disk table ---");

		final String[] files=accessionFile.split(",");

		outstream.println("Pass 1: counting digitizable entries.");
		final CountStats stats=CountStats.sum(spawnCountThreads(files));
		outstream.println("Total lines: "+stats.total);
		outstream.println("Digitizable: "+stats.digitizable);
		outstream.println("Non-digitizable: "+stats.nonDigitizable);
		outstream.println("Estimated Unique: "+stats.estimatedUnique);

		outstream.println("Pass 2: building disk table.");
		long collected=0, added=0;
		try(DiskAccessionTable.Builder builder=DiskAccessionTable.build(accOut, stats.estimatedUnique)){
			final DiskAccessionTable.BufferSet buffers=builder.bufferSet(bufferEntries());
			final ProcessThread[] threads=spawnThreads(files, buffers);
			// All started workers must finish before the shared builder closes.
			collected=waitForThreads(threads);
			buffers.flush();
			added=buffers.added();
		}catch(IOException e){
			throw new RuntimeException(e);
		}

		outstream.println("Collected "+collected+" entries.");
		outstream.println("Added "+added+" entries.");
		t.stop();
		outstream.println("Accession table built in "+t);
	}

	/** Reports whether an accession input row is a header or comment.
	 * @param line Raw input row
	 * @return True if builders should skip the row */
	private static boolean skipAccessionLine(byte[] line){
		return line.length<1 || Tools.startsWith(line, "accession") || Tools.startsWith(line, "#");
	}

	/** Closes one accession input and preserves an earlier processing failure.
	 * @param bf Input stream to close
	 * @param fname Input file name
	 * @param failure Earlier failure to keep primary, or null
	 * @throws IOException If close reports an error and there was no earlier failure */
	private static void closeByteFile(final ByteFile bf, final String fname,
			final Throwable failure) throws IOException{
		try{
			if(bf.close()){throw new IOException("Error reading "+fname);}
		}catch(IOException|RuntimeException|Error e){
			if(failure==null){throw e;}
			if(e!=failure){failure.addSuppressed(e);}
		}
	}

	/** Parses optional tree arguments.
	 * @param s Raw value: null for a bare flag, empty for an explicit empty/null option
	 * @return auto for enabled default trees, null for disabled trees, otherwise the path */
	private static String parseTreeFile(final String s){
		if(s==null){return "auto";}
		if(s.length()<1){return null;}
		if(Parse.isBoolean(s) || s.equals("1") || s.equals("0") ||
				s.equalsIgnoreCase("null") || s.equalsIgnoreCase("none")){
			return Parse.parseBoolean(s) ? "auto" : null;
		}
		return s;
	}

	/** Starts first-pass counting workers.
	 * If construction or startup fails after any worker starts, the caller joins
	 * every started worker before the original startup failure is rethrown.
	 * @param files Input file names
	 * @return Started workers sharing the file-name cursor */
	private CountThread[] spawnCountThreads(String[] files){
		final int threads=Tools.mid(1, files.length, Shared.threads());
		final CountThread[] array=new CountThread[threads];
		final AtomicInteger next=new AtomicInteger(0);
		final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
		try{
			for(int i=0; i<array.length; i++){
				array[i]=new CountThread(files, next, error);
				startWorker(array[i]);
			}
		}catch(RuntimeException|Error failure){
			error.compareAndSet(null, failure);
			joinStartedAfterStartupFailure(array, failure);
		}
		return array;
	}

	/** Starts second-pass table-building workers.
	 * If construction or startup fails after any worker starts, the caller joins
	 * every started worker before the original startup failure is rethrown.
	 * @param files Input file names
	 * @param buffers Shared partition buffers for the output table
	 * @return Started workers sharing the file-name cursor */
	private ProcessThread[] spawnThreads(String[] files, DiskAccessionTable.BufferSet buffers){
		final int threads=Tools.mid(1, files.length, Shared.threads());
		final ProcessThread[] array=new ProcessThread[threads];
		final AtomicInteger next=new AtomicInteger(0);
		final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
		try{
			for(int i=0; i<array.length; i++){
				array[i]=new ProcessThread(files, next, error, buffers);
				startWorker(array[i]);
			}
		}catch(RuntimeException|Error failure){
			error.compareAndSet(null, failure);
			joinStartedAfterStartupFailure(array, failure);
		}
		return array;
	}

	/** Starts a worker thread.
	 * @param worker Worker to start */
	protected void startWorker(final Thread worker){
		worker.start();
	}

	/** Joins workers after construction or startup fails, then rethrows the
	 * original startup failure.
	 * @param threads Partly populated worker array
	 * @param failure Construction or startup failure to rethrow */
	private static void joinStartedAfterStartupFailure(final Thread[] threads,
			final Throwable failure){
		boolean interrupted=false;
		for(int i=0; i<threads.length; i++){
			final Thread worker=threads[i];
			if(worker!=null){
				while(worker.isAlive()){
					try{
						worker.join();
					}catch(InterruptedException e){
						failure.addSuppressed(e);
						interrupted=true;
					}
				}
			}
		}
		if(interrupted){Thread.currentThread().interrupt();}
		rethrowStartupFailure(failure);
	}

	/** Rethrows a startup failure preserving its exact unchecked type.
	 * @param failure Failure to rethrow */
	private static void rethrowStartupFailure(final Throwable failure){
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		throw (Error)failure;
	}

	/** Joins every second-pass worker before reporting the first observed failure.
	 * A caller interruption is recorded and restored after all joins, so failure
	 * cannot close the shared builder while a worker is still using its buffers.
	 * @param threads Started workers to join, in error-reporting order
	 * @return Number of mappings sent to the output table buffers
	 * @throws RuntimeException After all workers finish, if a worker or join failed */
	private static long waitForThreads(ProcessThread[] threads){
		long collected=0;
		Throwable failure=null;
		boolean interrupted=false;
		for(int i=0; i<threads.length; i++){
			while(threads[i].isAlive()){
				try{
					threads[i].join();
				}catch(InterruptedException e){
					if(failure==null){failure=e;}
					interrupted=true;
				}
			}
			if(failure==null && threads[i].error!=null){failure=threads[i].error;}
			collected+=threads[i].collected;
		}
		if(interrupted){Thread.currentThread().interrupt();}
		if(failure!=null){throw new RuntimeException(failure);}
		return collected;
	}

	/** Chooses per-partition buffering from heap size when not set explicitly.
	 * @return Buffer entries per disk-table partition */
	private int bufferEntries(){
		if(bufferEntries>0){return bufferEntries;}
		final long memory=Runtime.getRuntime().maxMemory();
		final long entries=memory/(4L*BUFFER_ENTRY_BYTES*DiskAccessionTable.BUILD_PARTITIONS);
		return (int)Tools.mid(64, MAX_BUFFER_ENTRIES, entries);
	}

	/** First-pass worker that counts digitizable accession keys and estimates
	 * their distinct cardinality before the disk table is allocated. */
	private static class CountThread extends Thread{

		/** Creates one first-pass worker.
		 * @param files_ File list shared by all workers
		 * @param next_ Atomic next-file cursor
		 * @param error_ Shared worker-failure slot */
		CountThread(String[] files_, AtomicInteger next_, AtomicReference<Throwable> error_){
			files=files_;
			next=next_;
			errorRef=error_;
		}

		@Override
		public void run(){
			try{
				for(int i=next.getAndIncrement(); i<files.length && errorRef.get()==null; i=next.getAndIncrement()){
					processFile(files[i]);
				}
			}catch(Throwable e){
				error=e;
				errorRef.compareAndSet(null, e);
			}
		}

		/** Counts digitizable keys in one accession file, closing the input on
		 * normal completion and exceptional row-processing failures.
		 * @param fname Input file name
		 * @throws IOException If the input stream reports an error */
		private void processFile(String fname) throws IOException{
			final ByteFile bf=ByteFile.makeByteFile(fname, true);
			Throwable failure=null;
			try{
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					if(skipAccessionLine(line)){continue;}
					total++;
					final long number=AnalyzeAccession.digitize(line);
					if(number>0){
						unique.add(number);
						digitizable++;
					}else{
						nonDigitizable++;
					}
				}
			}catch(RuntimeException|Error e){
				failure=e;
				throw e;
			}finally{
				closeByteFile(bf, fname, failure);
			}
		}

		/** Shared accession input path list. */
		private final String[] files;
		/** Atomic next-file cursor shared by first-pass workers. */
		private final AtomicInteger next;
		/** Shared slot recording a worker failure. */
		private final AtomicReference<Throwable> errorRef;
		/** Approximate distinct digitized-key counter for this worker. */
		private final DynamicDemiLog unique=DynamicDemiLog.create(UNIQUE_BUCKETS, 31, -1, 0, false);
		/** Source rows processed by this worker. */
		private long total=0;
		/** Rows that digitized to a storable accession key. */
		private long digitizable=0;
		/** Rows skipped because the accession did not digitize. */
		private long nonDigitizable=0;
		/** First failure thrown by this worker. */
		private Throwable error=null;

	}

	/** Aggregate first-pass accession counts and distinct-key estimate. */
	private static class CountStats{

		/** Creates an empty accumulator for first-pass worker totals. */
		private CountStats(){}

		/** Joins all first-pass workers before reporting the first failure.
		 * Accumulator allocation and merging happen only after every worker has
		 * terminated, so failures cannot leave later workers unjoined. Caller
		 * interruption is restored before the recorded failure is reported.
		 * @param threads Workers to join and merge in error-reporting order
		 * @return Summed counters and unique-cardinality estimate
		 * @throws RuntimeException After all workers finish, if a worker or join failed */
		private static CountStats sum(CountThread[] threads){
			Throwable failure=null;
			boolean interrupted=false;
			for(int i=0; i<threads.length; i++){
				while(threads[i].getState()!=Thread.State.TERMINATED){
					try{
						threads[i].join();
					}catch(InterruptedException e){
						if(failure==null){failure=e;}
						interrupted=true;
					}
				}
				if(failure==null && threads[i].error!=null){failure=threads[i].error;}
			}
			if(interrupted){Thread.currentThread().interrupt();}
			if(failure!=null){throw new RuntimeException(failure);}

			CountStats stats=new CountStats();
			DynamicDemiLog unique=DynamicDemiLog.create(UNIQUE_BUCKETS, 31, -1, 0, false);
			for(int i=0; i<threads.length; i++){
				stats.total+=threads[i].total;
				stats.digitizable+=threads[i].digitizable;
				stats.nonDigitizable+=threads[i].nonDigitizable;
				unique.add(threads[i].unique);
			}
			stats.estimatedUnique=Tools.max(1, unique.cardinality());
			return stats;
		}

		/** Total source rows across first-pass workers. */
		private long total=0;
		/** Total rows that digitized to a storable accession key. */
		private long digitizable=0;
		/** Total rows skipped because the accession did not digitize. */
		private long nonDigitizable=0;
		/** Approximate number of distinct digitized keys. */
		private long estimatedUnique=0;

	}

	/** Second-pass worker that buffers digitized accession to TaxID mappings for
	 * the partitioned disk table builder. */
	private static class ProcessThread extends Thread{

		/** Creates one second-pass worker.
		 * @param files_ File list shared by all workers
		 * @param next_ Atomic next-file cursor
		 * @param error_ Shared worker-failure slot
		 * @param buffers_ Shared partition buffers for the output table */
		ProcessThread(String[] files_, AtomicInteger next_, AtomicReference<Throwable> error_,
				DiskAccessionTable.BufferSet buffers_){
			files=files_;
			next=next_;
			errorRef=error_;
			buffers=buffers_;
		}

		@Override
		public void run(){
			try{
				for(int i=next.getAndIncrement(); i<files.length && errorRef.get()==null; i=next.getAndIncrement()){
					processFile(files[i]);
				}
			}catch(Throwable e){
				error=e;
				errorRef.compareAndSet(null, e);
			}
		}

		/** Parses one accession file and buffers valid digitized mappings,
		 * closing the input on normal completion and exceptional row-processing
		 * or buffer-insertion failures.
		 * @param fname Input file name
		 * @throws IOException If input reading or buffer insertion fails */
		private void processFile(String fname) throws IOException{
			final ByteFile bf=ByteFile.makeByteFile(fname, true);
			Throwable failure=null;
			try{
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					if(skipAccessionLine(line)){continue;}
					final long number=AnalyzeAccession.digitize(line);
					if(number>0){
						final int taxid=AccessionToTaxid.parseLineToTaxid(line, (byte)'\t');
						if(taxid>0){
							buffers.add(number, taxid);
							collected++;
						}
					}
				}
			}catch(IOException|RuntimeException|Error e){
				failure=e;
				throw e;
			}finally{
				closeByteFile(bf, fname, failure);
			}
		}

		/** Shared accession input path list. */
		private final String[] files;
		/** Atomic next-file cursor shared by second-pass workers. */
		private final AtomicInteger next;
		/** Shared slot recording a worker failure. */
		private final AtomicReference<Throwable> errorRef;
		/** Shared disk-table partition buffers. */
		private final DiskAccessionTable.BufferSet buffers;
		/** Mappings sent to the disk-table buffers by this worker. */
		private long collected=0;
		/** First failure thrown by this worker. */
		private Throwable error=null;

	}

	/*--------------------------------------------------------------*/

	/** Comma-delimited accession2taxid input paths. */
	private String accessionFile=null;
	/** Serialized or source GI table input path. */
	private String giTableFile=null;
	/** Optional taxonomy tree input path. */
	private String taxTreeFile=null;
	/** Accession-pattern table required for disk accession builds. */
	private String patternFile=null;
	/** Disk-backed accession table output path. */
	private String accOut=null;
	/** Disk-backed GI table output path. */
	private String giOut=null;
	/** Buffered entries per accession-table partition, or zero for heap-based sizing. */
	private int bufferEntries=0;
	/** Optional loaded tree retained for future tree-dependent build steps. */
	private TaxTree tree=null;
	/** Destination for status messages. */
	private PrintStream outstream=System.err;

	/** Approximate bytes per buffered accession-table entry. */
	private static final long BUFFER_ENTRY_BYTES=20;
	/** Upper bound on per-partition buffered entries. */
	private static final int MAX_BUFFER_ENTRIES=1<<24;
	/** DynamicDemiLog buckets used to estimate distinct digitized accessions. */
	private static final int UNIQUE_BUCKETS=1<<16;

}
