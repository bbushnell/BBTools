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
 * Template class for processing SAM/BAM files in multiple threads.
 * Designed to be easily modified by filling in the processRead method.
 * Provides a framework for multi-threaded read processing with SAM filtering, realignment, and reference loading.
 * @author Brian Bushnell
 * @date September 6, 2019
 */
public class A_SampleSamStreamer implements Accumulator<A_SampleSamStreamer.ProcessThread> {
	
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
	 * Sets up input/output files, SAM filtering, threading parameters, and validates configuration.
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
			String[] split=arg.split("=");
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
	
	/** Validates that all required files exist and can be accessed.
	 * Ensures input and reference files are readable, output files are writable, and no file paths are duplicated. */
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
	
	/** Adjusts static configuration for optimal file I/O performance.
	 * Sets ByteFile modes based on thread count and validates FASTA stream settings. */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/** Validates that parameter ranges are within bounds and required parameters are set.
	 * @return true if all parameters are valid */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing method that orchestrates the read processing pipeline.
	 * Creates input/output streams, loads reference data, spawns worker threads, and reports timing statistics.
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
		//Close the read streams
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

	private void loadScafMapFromReference(){
		if(loadedRef){return;}
		assert(ref!=null);
		scafMap=ScafMap.loadReference(ref, scafMap, samFilter, true);
		if(realign){Realigner.map=scafMap;}
		loadedRef=true;
	}

	private void loadReferenceCustom(){
		ConcurrentReadInputStream cris=makeRefCris();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()) {
			//Do something
			assert(false) : "TODO";
		}
	}
	
	private ConcurrentReadInputStream makeRefCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffref, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		assert(!paired) : "References should not be paired.";
		return cris;
	}
	
	private Streamer makeStreamer(FileFormat ff){
		if(ff==null){return null;}
		Streamer ss=StreamerFactory.makeSamOrBamStreamer(ff, streamerThreads, true, false, maxReads, true);
		ss.start(); //Start the stream
		if(verbose){outstream.println("Started Streamer");}
		return ss;
	}
	
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
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		readsOut+=pt.readsOutT;
		basesOut+=pt.basesOutT;
		errorState|=(!pt.success);
	}
	
	/** Reports overall processing success status.
	 * @return true if no errors occurred during processing */
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	class ProcessThread extends Thread {
		
		//Constructor
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
		 * Core processing loop that retrieves and processes read lists until input is exhausted.
		 */
		void processInner(){
			
			//Grab and process all lists
			for(ListNum<Read> ln=ss.nextList(); ln!=null; ln=ss.nextList()){
//				if(verbose){outstream.println("Got list of size "+list.size());} //Disabled due to non-static access
				
				processList(ln);
			}
			
		}
		
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
		 * Template method for processing individual reads.
		 * Currently contains basic validation and SAM filtering logic and is intended for customization.
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

		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		
		protected long readsOutT=0;
		protected long basesOutT=0;
		
		/** Flag indicating successful completion of this thread. */
		boolean success=false;
		
		/** Shared input stream for reading SAM data. */
		private final Streamer ss;
		/** Shared output stream for writing processed reads. */
		private final ConcurrentReadOutputStream ros;
		final int tid;
		/**
		 * Realigner instance for adjusting read alignments (null if realignment disabled).
		 */
		final Realigner realigner;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in=null;
	/** Reference file path for scaffold mapping. */
	private String ref=null;

	private String out=null;
	
	private String extin=null;
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

	/** Maximum number of reads to process; -1 for no limit. */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	/** Number of threads dedicated to reading the SAM file. */
	private int streamerThreads=-1;
	
	private boolean loadedRef=false;
	
	private boolean realign=false;
	
	private int ploidy=1;
	
	public ScafMap scafMap;
	public final SamFilter samFilter=new SamFilter();
	
	/** Provides read-write lock for thread-safe access.
	 * @return ReadWriteLock instance for synchronization */
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
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
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean ordered=false;
	
}
