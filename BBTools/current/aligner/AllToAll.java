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
 * 
 * @author Brian Bushnell
 * @date January 27, 2020
 *
 */
public class AllToAll implements Accumulator<AllToAll.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
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
	 * Constructor.
	 * @param args Command line arguments
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
	
	/** Parse arguments from the command line */
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
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		qfin1=Tools.fixExtension(qfin1);
	}
	
	/** Ensure files can be read and written */
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
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all reads */
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
	
	private static void mirrorMatrix(float[][] matrix){
		for(int i=0; i<matrix.length; i++) {
			for(int j=i; j<matrix.length; j++) {
				assert(matrix[i][j]==0) : matrix[i][j];
				matrix[i][j]=(i==j ? 1 : matrix[j][i]);
			}
		}
	}
	
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
	
	/** Spawn process threads */
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
	
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		alignments+=pt.alignmentsT;
		errorState|=(!pt.success);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	static class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ArrayList<Read> reads_, float[][] results_, final AtomicInteger atom_, final int tid_){
			reads=reads_;
			results=results_;
			atom=atom_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
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
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		boolean processReadPair(final Read r1, final Read r2){
			throw new RuntimeException("TODO: Implement this method."); //TODO
//			return true;
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long alignmentsT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Thread ID */
		final int tid;
		
		final ArrayList<Read> reads;
		final float[][] results;
		final AtomicInteger atom;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	
	private String qfin1=null;

	/** Primary output file path */
	private String out1=null;
	
	/** Override input file extension */
	private String extin=null;
	
	/*--------------------------------------------------------------*/

	ArrayList<Read> reads;
	float[][] results;
	
	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long alignments=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Primary output file */
	private final FileFormat ffout1;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	/** Reads are output in input order */
	private boolean ordered=false;
	
}
