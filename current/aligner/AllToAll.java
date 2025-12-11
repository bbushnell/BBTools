package aligner;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.SketchObject;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Aligns all sequences to all sequences and produces an identity matrix.
 * Performs pairwise sequence alignment between all input sequences using SketchObject.align()
 * and writes identity scores to a similarity matrix with optional multithreading.
 * @author Brian Bushnell
 * @date January 27, 2020
 */
public class AllToAll implements Accumulator<AllToAll.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Program entry point.
	 * Constructs an AllToAll instance, runs processing with a Timer, and closes output streams.
	 * @param args Command line arguments for input/output files and parameters
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		AllToAll x=new AllToAll(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses command line arguments and initializes file formats.
	 * Sets up input/output handlers, validates parameters, and prepares for sequence loading and alignment.
	 * @param args Command line arguments containing file paths and options
	 */
	public AllToAll(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			qfin1=parser.qfin1;
			extin=parser.extin;

			out1=parser.out1;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses command line arguments using the standard BBTools Parser framework.
	 * Processes verbose/ordered flags and standard parser fields.
	 * @param args Command line arguments to parse
	 * @return Configured Parser object with parsed parameters
	 */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		
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
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
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
	 * Adds or removes compression extensions (e.g., .gz or .bz2) as needed for input files.
	 */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		qfin1=Tools.fixExtension(qfin1);
	}
	
	/** Validates that input files exist, output files can be written, and no duplicate paths are specified.
	 * Throws RuntimeException if file access requirements are not met. */
	private void checkFileExistence(){
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/**
	 * Adjusts static settings for optimal file I/O performance, including ByteFile mode selection.
	 */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/**
	 * Validates that required parameters are set and within acceptable ranges.
	 * Ensures at least one input file is specified and that extension overrides are consistent.
	 * @return true if validation passes
	 * @throws RuntimeException if required parameters are missing
	 */
	private boolean validateParams(){
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing method that performs all-to-all sequence alignment.
	 * Loads all sequences, spawns worker threads for pairwise alignments, mirrors the matrix, and writes results.
	 * @param t Timer for tracking execution time and reporting performance
	 */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Reset counters
		readsProcessed=alignments=0;
		basesProcessed=0;
		
		//Fetch data
		reads=ConcurrentReadInputStream.getReads(maxReads, true, ffin1, null, qfin1, null); //TODO:  Note that this does not return the error state
		results=new float[reads.size()][];
		
		outstream.println("Loaded "+reads.size()+" sequences.");
		
		//Process the reads in separate threads
		spawnThreads();
		mirrorMatrix(results);
		if(verbose){outstream.println("Finished alignment.");}
		
		printResults();
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.number("Alignments:", alignments, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Mirrors an alignment matrix across the diagonal to fill the upper triangle.
	 * Sets diagonal elements to 1.0 and copies lower-triangle values to corresponding upper-triangle positions.
	 * @param matrix Square similarity matrix with computed lower triangle values
	 */
	private static void mirrorMatrix(float[][] matrix){
		for(int i=0; i<matrix.length; i++) {
			for(int j=i; j<matrix.length; j++) {
				assert(matrix[i][j]==0) : matrix[i][j];
				matrix[i][j]=(i==j ? 1 : matrix[j][i]);
			}
		}
	}
	
	/** Outputs the similarity matrix to the specified output file.
	 * Writes a header line of sequence names followed by tab-separated identity percentages with two decimal places. */
	private void printResults(){
		if(ffout1==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		final int max=reads.size();
		bsw.print("Name");
		for(int rnum=0; rnum<max; rnum++){
			bsw.tab().print(reads.get(rnum).id);
		}
		bsw.println();
		for(int qnum=0; qnum<max; qnum++){
			bsw.print(reads.get(qnum).id);
			final float[] scores=results[qnum];
			for(int rnum=0; rnum<max; rnum++){
				bsw.tab().print(100*scores[rnum], 2);
			}
			bsw.println();
		}
		bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates and manages worker threads for parallel alignment processing.
	 * Uses an AtomicInteger for work distribution and ThreadWaiter to start and join all ProcessThreads. */
	private void spawnThreads(){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		AtomicInteger atom=new AtomicInteger(0);
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(reads, results, atom, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	/**
	 * Accumulates statistics from a completed worker thread.
	 * Aggregates read counts, base counts, alignment counts, and error status.
	 * @param pt Completed ProcessThread containing execution statistics
	 */
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		alignments+=pt.alignmentsT;
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
	
	/** Worker thread for processing query sequences against all reference sequences.
	 * Each thread claims query indices atomically and aligns each query against all earlier sequences to avoid redundant work. */
	static class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ArrayList<Read> reads_, float[][] results_, final AtomicInteger atom_, final int tid_){
			reads=reads_;
			results=results_;
			atom=atom_;
			tid=tid_;
		}
		
		//Called by start()
		/** Thread execution entry point.
		 * Calls processInner() to perform alignment work and marks success status on completion. */
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Claims and processes query sequences using atomic work distribution.
		 * Continues until all sequences have been processed by some ProcessThread. */
		void processInner(){
			
			for(int next=atom.getAndIncrement(); next<reads.size(); next=atom.getAndIncrement()){
				processQuery(next);
			}
			
		}
		
		void processQuery(final int qnum){
			final Read query=reads.get(qnum);
			final float[] scores=new float[reads.size()];
			readsProcessedT++;
			basesProcessedT+=query.length();
			for(int rnum=0; rnum<qnum; rnum++){
				final Read ref=reads.get(rnum);
				float identity=SketchObject.align(query.bases, ref.bases);
				scores[rnum]=identity;
				alignmentsT++;
			}
			synchronized(results){
				results[qnum]=scores;
			}
		}
		
		boolean processReadPair(final Read r1, final Read r2){
			throw new RuntimeException("TODO: Implement this method."); //TODO
//			return true;
		}

		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		
		/** Number of alignments computed by this thread */
		protected long alignmentsT=0;
		protected long basesOutT=0;
		
		boolean success=false;
		
		final int tid;
		
		final ArrayList<Read> reads;
		final float[][] results;
		final AtomicInteger atom;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	
	private String qfin1=null;

	private String out1=null;
	
	private String extin=null;
	
	/*--------------------------------------------------------------*/

	ArrayList<Read> reads;
	float[][] results;
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;

	/** Number of pairwise alignments computed */
	protected long alignments=0;

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file format handler */
	private final FileFormat ffin1;
	
	/** Primary output file format handler */
	private final FileFormat ffout1;
	
	/** Returns the read-write lock for thread synchronization.
	 * @return ReadWriteLock instance for coordinating access */
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean ordered=false;
	
}
