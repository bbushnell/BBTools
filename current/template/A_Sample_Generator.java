package template;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

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
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Template for multithreaded synthetic read generators that preload template reads.
 * Concrete generators supply the read-construction logic and paired-output policy.
 *
 * <p>Created June 8, 2019.
 * @author Brian Bushnell
 */
public class A_Sample_Generator{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		A_Sample_Generator x=new A_Sample_Generator(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses arguments, validates configuration, and creates file-format descriptors.
	 * @param args command-line arguments
	 */
	public A_Sample_Generator(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			extin=parser.extin;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			extout=parser.extout;
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses tool-specific flags and delegates standard flags to {@link Parser}.
	 * @param args command-line arguments to parse
	 * @return parser containing shared BBTools settings
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
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
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
	 * Applies paired-output # replacement and validates output pairing.
	 * Concrete generators can require template input here if needed.
	 */
	private void doPoundReplacement(){

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		//Ensure there is an input file
		assert(false) : "May or may not be required...";
//		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}
	
	/** Adds or removes .gz and .bz2 extensions from input names as needed. */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
	}
	
	/**
	 * Ensures inputs are readable, outputs are writable, and file names are unique.
	 * @throws RuntimeException if any file path is unusable
	 */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjusts static file I/O settings for this tool. */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/**
	 * Validates generator-specific parameter ranges and required settings.
	 * @return true if all parameters are valid
	 */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads template reads, runs generator workers, and reports final statistics.
	 * @param t timer for tracking execution time
	 * @throws RuntimeException if a worker, statistics writer, or stream close fails
	 */
	void process(Timer t){
		
		//Force validation while loading template reads.
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=true;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Optionally create a read output stream
		assert(false) : "TODO: Determine whether output should be paired.";
		final ConcurrentReadOutputStream ros=makeCros(false);
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Load all template reads before starting generator workers.
		data=loadData(cris);
		
		//Process the reads in separate threads
		spawnThreads(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros);
		
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
	
	/**
	 * Creates and starts a concurrent read input stream.
	 * @return stream for loading template reads
	 */
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/**
	 * Creates and starts a concurrent read output stream.
	 * @param pairedInput true when generated reads should be written as paired
	 * @return output stream, or null if no output file was specified
	 */
	private ConcurrentReadOutputStream makeCros(boolean pairedInput){
		if(ffout1==null){return null;}

		//Set output buffer size
		final int buff=4;

		//Notify user of output mode
		if(pairedInput && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
			outstream.println("Writing interleaved.");
		}

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(
				ffout1, ffout2, qfout1, qfout2, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/**
	 * Starts generator workers, waits for completion, and aggregates counters.
	 * @param cris input stream used to load template reads
	 * @param ros output stream for generated reads
	 */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(ros, i, maxReads, nextReadID));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		waitForThreads(alpt);
		
		//Do anything necessary after processing
		
	}
	
	/**
	 * Waits for workers to complete and aggregates their output counters.
	 * @param alpt worker threads to wait for
	 */
	private void waitForThreads(ArrayList<ProcessThread> alpt){
		
		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try{
					//Attempt a join operation
					pt.join();
				}catch(InterruptedException e){
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads all input reads into memory for generator-specific sampling.
	 * @param cris stream of template reads
	 * @return loaded template reads
	 */
	private ArrayList<Read> loadData(ConcurrentReadInputStream cris){
		
		ArrayList<Read> input=new ArrayList<Read>();
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();

		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){
//			if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
			
			for(Read r : ln){
				
				//Optional filter criteria
				if(true){
					input.add(r);
				}
				
				//Increment counters
				readsProcessed+=r.pairCount();
				basesProcessed+=r.pairLength();
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
//			if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
			
			//Fetch a new list
			ln=cris.nextList();
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		
		return input;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Worker that generates batches of synthetic reads from a shared ID counter. */
	private static class ProcessThread extends Thread{
		
		//Constructor
		ProcessThread(final ConcurrentReadOutputStream ros_, final int tid_, 
				final long maxReads_, final AtomicLong nextReadID_){
			ros=ros_;
			tid=tid_;
			maxReads=maxReads_;
			atomicReadID=nextReadID_;
		}
		
		//Called by start()
		/** Generates reads and records successful completion. */
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Claims ID ranges from the shared counter and generates each claimed batch. */
		void processInner(){

			//As long as there is a nonempty read list...
			for(long generated=atomicReadID.getAndAdd(readsPerList); generated<maxReads; 
					generated=atomicReadID.getAndAdd(readsPerList)){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

				long toGenerate=Tools.min(readsPerList, maxReads-generated);
				
				ArrayList<Read> reads=generateList((int)toGenerate, generated);

				//Output reads to the output stream
				if(ros!=null){ros.add(reads, 0);}
			}
		}
		
		/**
		 * Generates one output batch and updates per-thread counters.
		 * @param toGenerate number of reads to generate in this batch
		 * @param nextID first read ID to assign in the batch
		 * @return generated reads
		 */
		private ArrayList<Read> generateList(int toGenerate, long nextID){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=new ArrayList<Read>(toGenerate);
			
			//Loop through each read in the list
			for(int i=0; i<toGenerate; i++, nextID++){
				Read r=generateRead(nextID);
				readsOutT+=r.pairCount();
				basesOutT+=r.pairLength();
				reads.add(r);
			}

			return reads;
		}
		
		/**
		 * Generates one read with the specified ID.
		 * @param nextID read ID to assign
		 * @return generated read
		 */
		private Read generateRead(final long nextID){
//			Read r=new Read(null, null, nextID);
			throw new RuntimeException("TODO: Implement this method."); //TODO
		}
		
		/** Number of reads generated by this thread */
		protected long readsOutT=0;
		/** Number of bases generated by this thread */
		protected long basesOutT=0;
		
		/** True only if this thread completed successfully without errors */
		boolean success=false;
		
		private final AtomicLong atomicReadID;
		private final long maxReads;
		private final int readsPerList=Shared.bufferLen();
		
		/** Output stream for writing generated reads */
		private final ConcurrentReadOutputStream ros;
		/** Thread identifier for debugging and logging */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file path for template sequences */
	private String in1=null;

	/** Primary output file path for generated reads */
	private String out1=null;
	/** Secondary output file path for paired-end reads */
	private String out2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	/** Override extension for input file format detection */
	private String extin=null;
	/** Override extension for output file format */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Total number of input reads processed from template files */
	protected long readsProcessed=0;
	/** Total number of input bases processed from template files */
	protected long basesProcessed=0;

	/** Total number of reads generated and written to output */
	protected long readsOut=0;
	/** Total number of bases generated and written to output */
	protected long basesOut=0;

	/** Maximum number of reads to generate; -1 for unlimited */
	private long maxReads=-1;
	
	/** Template reads loaded into memory for concrete generators */
	private ArrayList<Read> data=new ArrayList<Read>();
	
	private AtomicLong nextReadID=new AtomicLong(0);
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** File format specification for primary input file */
	private final FileFormat ffin1;
	
	/** File format specification for primary output file */
	private final FileFormat ffout1;
	/** File format specification for secondary output file */
	private final FileFormat ffout2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging and progress tracking */
	public static boolean verbose=false;
	/** True if an error occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	
}
