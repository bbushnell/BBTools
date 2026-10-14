package template;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;
import tracker.ReadStats;
import var2.Realigner;
import var2.SamFilter;
import var2.ScafMap;

/**
 * Incomplete SAM/BAM processing template with worker-local read processing.
 * Copy it and resolve the parameter, reference, output-format and processing
 * TODOs before use; assertions deliberately stop the untouched skeleton.
 * Workers share a Streamer and optional legacy output stream. The supplied
 * processing assumes unlinked SAM-record roots, filters short/absent sequences
 * and SAM flags, and leaves the actual transformation to the custom hook.
 * Reference scaffolds and optional Realigner objects are prepared, but the
 * default hook does not perform realignment.
 * <p>Created September 6, 2019.
 * @author Brian Bushnell, Shinobu
 */
public class A_SampleSamStreamer implements Accumulator<A_SampleSamStreamer.ProcessThread>{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point that constructs an A_SampleSamStreamer and processes reads.
	 * @param args Command line arguments */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		A_SampleSamStreamer x=new A_SampleSamStreamer(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses arguments and initializes the processing environment.
	 * Configures descriptors and global I/O/filter settings. Intentional TODO assertions
	 * require decisions about validation, reference requirements and output format.
	 * @param args Command line arguments for configuration
	 */
	public A_SampleSamStreamer(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
//		samFilter.includeUnmapped=false;
//		samFilter.includeSupplementary=false;
//		samFilter.includeDuplicate=false;
//		samFilter.includeNonPrimary=false;
//		samFilter.includeQfail=false;
//		samFilter.minMapq=4;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in=parser.in1;
			extin=parser.extin;

			out=parser.out1;
			extout=parser.extout;
		}
		
		{
//			if("auto".equalsIgnoreCase(atomic)){Scaffold.setCA3A(Shared.threads()>8);}
//			else{Scaffold.setCA3A(Parse.parseBoolean(atomic));}

			if(ploidy<1){System.err.println("WARNING: ploidy not set; assuming ploidy=1."); ploidy=1;}
			samFilter.setSamtoolsFilter();
			
			streamerThreads=Tools.max(1, Tools.min(streamerThreads, Shared.threads()));
			assert(streamerThreads>0) : streamerThreads;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout=FileFormat.testOutput(out, FileFormat.SAM, extout, true, overwrite, append, ordered);
		assert(false) : "TODO: Default output format might be fasta.";

		//Create input FileFormat objects
		ffin=FileFormat.testInput(in, FileFormat.SAM, extin, true, true);
		ffref=FileFormat.testInput(ref, FileFormat.FASTA, null, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses command-line arguments into configuration parameters.
	 * Handles verbose mode, reference file, ordering, realignment, ploidy, SAM filtering, and other options.
	 * @param args Command line arguments to parse
	 * @return Configured Parser object with standard parameters set
	 */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			//Preserve the complete value after the first equals sign.
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("realign")){
				realign=Parse.parseBoolean(b);
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("clearfilters")){
				if(Parse.parseBoolean(b)){
					samFilter.clear();
				}
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(samFilter.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/**
	 * Adds or removes .gz or .bz2 extensions as needed for input and reference files.
	 */
	private void fixExtensions(){
		in=Tools.fixExtension(in);
		ref=Tools.fixExtension(ref);
	}
	
	/** Checks primary paths after the reference-policy TODO is resolved.
	 * The skeleton currently requires both input and reference, tests access and
	 * rejects duplicate input/reference/output path strings. */
	private void checkFileExistence(){

		//Ensure there is an input file
		if(in==null){throw new RuntimeException("Error - an input file is required.");}

		//Ensure there is an input file
		assert(false) : "TODO: Check.";
		if(ref==null){throw new RuntimeException("Error - a reference file is required.");}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in, ref, out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Enables BF2 when BF1/BF2 are not forced and more than two threads are available.
	 * Also asserts the current FASTA configuration. */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/** Intentional parameter-validation stub; replace the TODO before use.
	 * @return True after validation is implemented; currently only reachable with assertions off */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Starts input, loads reference scaffolds, creates output and waits for all workers.
	 * Finalizes input/output, collects their error states and restores constructor
	 * validation on the normal path before reporting any accumulated error.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final Streamer ss=makeStreamer(ffin);
		
		//Load reference, if desired (and if present);
		loadScafMapFromReference();
//		loadReferenceCustom();
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(ss, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Include input completion errors as well as output completion errors.
		errorState|=ReadWrite.closeStream(ss);
		errorState|=ReadWrite.closeStream(ros);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/** Loads reference scaffolds once and supplies the global Realigner map when requested. */
	private void loadScafMapFromReference(){
		if(loadedRef){return;}
		assert(ref!=null);
		scafMap=ScafMap.loadReference(ref, scafMap, samFilter, true);
		if(realign){Realigner.map=scafMap;}
		loadedRef=true;
	}

	/** Dormant alternative reference hook; its call is commented out.
	 * Before enabling it, implement the TODO, return fetched lists and close the stream. */
	private void loadReferenceCustom(){
		ConcurrentReadInputStream cris=makeRefCris();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			//Do something
			assert(false) : "TODO";
		}
	}
	
	/** Starts the unpaired reference stream used only by the dormant custom hook.
	 * The main input limit is also forwarded here. Caller owns list return and closure.
	 * @return Started reference input stream */
	private ConcurrentReadInputStream makeRefCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffref, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		assert(!paired) : "References should not be paired.";
		return cris;
	}
	
	/** Starts a Read-producing reader with shared SAM header retention.
	 * @param ff Input descriptor; the intended input is SAM/BAM
	 * @return Started reader, or null for a null descriptor */
	private Streamer makeStreamer(FileFormat ff){
		if(ff==null){return null;}
		Streamer ss=StreamerFactory.makeSamOrBamStreamer(ff, streamerThreads, true, false, maxReads, true);
		ss.start(); //Start the stream
		if(verbose){outstream.println("Started Streamer");}
		return ss;
	}
	
	/** Starts optional legacy output with the requested batch ordering.
	 * Shared-header forwarding is disabled here; SAM output uses a generated header.
	 * A customized tool must choose its header/reference policy before processing alignments.
	 * @return Started output stream, or null when output is disabled */
	private ConcurrentReadOutputStream makeCros(){
		if(ffout==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Starts Shared.threads workers and collects each through ThreadWaiter.
	 * @param ss Shared started input reader
	 * @param ros Shared started output stream, or null */
	private void spawnThreads(final Streamer ss, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(ss, ros, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
		errorState|=!success;
		
		//Do anything necessary after processing
		
	}
	
	/** Collects one completed worker under the caller's accumulator lock.
	 * @param pt Worker to collect once */
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		readsOut+=pt.readsOutT;
		basesOut+=pt.basesOutT;
		errorState|=(!pt.success);
	}
	
	/** Reports the accumulated worker/input/output/statistics status.
	 * @return True when no accumulated error has been recorded */
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Processes each fetched batch before handing it to the shared output. */
	class ProcessThread extends Thread{
		
		/** Binds shared I/O and optionally creates a realigner for custom processing.
		 * @param ss_ Started reader
		 * @param ros_ Optional started output
		 * @param tid_ Worker index */
		ProcessThread(final Streamer ss_, final ConcurrentReadOutputStream ros_, final int tid_){
			ss=ss_;
			ros=ros_;
			tid=tid_;
			realigner=(realign ? new Realigner() : null);
		}
		
		//Called by start()
		/** Main execution method for the worker thread.
		 * Processes reads and sets the success flag upon completion. */
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/**
		 * Consumes through terminal null, passing empty batches through with their IDs.
		 */
		void processInner(){
			
			//Grab and process all lists
			for(ListNum<Read> ln=ss.nextList(); ln!=null; ln=ss.nextList()){
//				if(verbose){outstream.println("Got list of size "+list.size());} //Disabled due to non-static access
				
				processList(ln);
			}
			
		}
		
		/** Validates roots, counts input before customization and nulls rejected entries.
		 * Counters assume unlinked SAM-record roots: bases use root length and retained
		 * reads count roots, while the input read count adds pairCount. Adapt these
		 * together if the template is changed to accept linked pairs.
		 * @param ln Owned batch until output handoff; do not mutate it after submission */
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r=reads.get(idx);
				
				//Validate reads in worker threads
				if(!r.validated()){r.validate(true);}

				//Track the initial length for statistics
				final int initialLength=r.length();

				//Increment counters
				readsProcessedT+=r.pairCount();
				basesProcessedT+=initialLength;
				
				{
					//Reads are processed in this block.
					boolean keep=processRead(r);
					
					if(!keep){reads.set(idx, null);}
					else{
						readsOutT++;
						basesOutT+=r.length();
					}
				}
			}

			//Output reads to the output stream
			if(ros!=null){ros.add(reads, ln.id);}
		}
		
		/**
		 * Rejects missing/short bases and records failing the shared SAM filter.
		 * Passing records reach the intentional TODO assertion; implement the actual
		 * transformation there before use. The optional realigner is not invoked here.
		 * @param r Read to process
		 * @return true if the read should be retained, false if it should be discarded
		 */
		boolean processRead(final Read r){
			if(r.bases==null || r.length()<=1){return false;}
			final SamLine sl=r.samline;
			if(samFilter!=null && !samFilter.passesFilter(sl)){return false;}
			
//			System.err.println("A: "+sl);
			
//			final SamLine oldSL=new SamLine(sl);
//			final Read oldRead=r.clone();
			
			assert(false) : "TODO";
			return true;
		}

		/** Sum of pairCount for input roots before filtering. */
		protected long readsProcessedT=0;
		/** Input root bases before filtering; does not add mate bases. */
		protected long basesProcessedT=0;
		
		/** Retained roots, one per accepted SAM record. */
		protected long readsOutT=0;
		/** Retained root bases after customization. */
		protected long basesOutT=0;
		
		/** Flag indicating successful completion of this thread. */
		boolean success=false;
		
		/** Shared input stream for reading SAM data. */
		private final Streamer ss;
		/** Shared output stream for writing processed reads. */
		private final ConcurrentReadOutputStream ros;
		/** Worker index, currently unused by the supplied hook. */
		final int tid;
		/**
		 * Optional realigner provided for customization; the supplied hook never uses it.
		 */
		final Realigner realigner;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input path. */
	private String in=null;
	/** Reference file path for scaffold mapping. */
	private String ref=null;

	/** Optional primary output path. */
	private String out=null;
	
	/** Input extension override. */
	private String extin=null;
	/** Output extension override. */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Total number of reads processed across all threads. */
	protected long readsProcessed=0;
	/** Total number of bases processed across all threads. */
	protected long basesProcessed=0;

	/** Total number of reads retained after processing. */
	protected long readsOut=0;
	/** Total number of bases retained after processing. */
	protected long basesOut=0;

	/** Requested input limit forwarded to the reader; negative for unlimited. */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	/** Reader-selection hint; the default -1 is clamped to one by the constructor. */
	private int streamerThreads=-1;
	
	/** Whether the reference-loading step completed. */
	private boolean loadedRef=false;
	
	/** Whether to provision a realigner for each worker; the default hook does not use it. */
	private boolean realign=false;
	
	/** Reserved customization setting, normalized to at least one. */
	private int ploidy=1;
	
	public ScafMap scafMap;
	/** Shared filter configured before worker startup. */
	public final SamFilter samFilter=new SamFilter();
	
	/** Supplies the lock used by ThreadWaiter for accumulation, not all object access.
	 * @return Accumulator coordination lock */
	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	/** Lock used for ThreadWaiter accumulation. */
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file format specification. */
	private final FileFormat ffin;
	/** Reference input file format specification. */
	private final FileFormat ffref;
	
	/** Primary output file format specification. */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages. */
	private PrintStream outstream=System.err;
	/** Whether to print selected lifecycle diagnostics. */
	public static boolean verbose=false;
	/** Accumulated worker/input/output/statistics error flag. */
	public boolean errorState=false;
	/** Whether output may replace an existing file. */
	private boolean overwrite=true;
	/** Whether output requests append mode. */
	private boolean append=false;
	/** Whether legacy output reorders batches by their IDs. */
	private boolean ordered=false;
	
}
