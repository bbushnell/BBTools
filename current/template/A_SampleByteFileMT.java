package template;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Template for tools that read one small text file before processing a second
 * text file with worker threads. The single-threaded prepass can write valid
 * prefixes; dormant invalid-line output scaffolding remains disabled because
 * outInvalid is private and never assigned. The threaded pass parses valid
 * lines without writing output.
 *
 * <p>Created February 6, 2023.
 * @author Brian Bushnell
 */
public class A_SampleByteFileMT implements Accumulator<A_SampleByteFileMT.ProcessThread>{

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
		A_SampleByteFileMT x=new A_SampleByteFileMT(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses arguments, validates paths and initializes file descriptors.
	 * @param args command-line arguments
	 */
	public A_SampleByteFileMT(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.setZipThreads(Shared.threads());

		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;

			in1=parser.in1;
			in2=parser.in2;
			out1=parser.out1;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.TXT, null, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses command-line arguments and delegates standard flags to Parser.
	 * @param args command-line arguments
	 * @return Parser with recognized standard fields populated
	 */
	private Parser parse(String[] args){

		//Create a parser object
		Parser parser=new Parser();
		//		parser.out1="stdout";

		//Set any necessary Parser defaults here
		//parser.foo=bar;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Split on only the first '=' so argument values may contain '='.
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
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		return parser;
	}

	/**
	 * Normalizes input compression suffixes and requires both input paths.
	 */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		if(in1==null || in2==null){throw new RuntimeException("Error - two input files are required.");}
	}

	/**
	 * Verifies read/write path accessibility and rejects duplicate path strings
	 * across the two inputs and optional outputs.
	 */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, outInvalid)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, outInvalid)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Enables ByteFile2 when neither BF1 nor BF2 is forced and more than two
	 * threads are available. */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}

	/**
	 * Placeholder parameter validator; replace before using this template.
	 * @return true if parameters are valid
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
	 * Runs the single-threaded first-input prepass, then the threaded second-input
	 * pass, and reports aggregate nonempty-line statistics.
	 * @param t timer started at program launch
	 */
	void process(Timer t){

		//Reset counters
		linesProcessed=linesOut=0;
		bytesProcessed=bytesOut=0;

		processSinglethreaded(ffin1);

		//Process the reads in separate threads
		spawnThreads(ffin2);

		if(verbose){outstream.println("Finished; closing streams.");}

		//Report timing and results
		t.stop();

		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));

		outstream.println();
		outstream.println("Valid Lines:       \t"+linesOut);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesOut));

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Starts worker threads that share one ByteFile and accumulates their counters.
	 * @param ffin input descriptor to process with workers
	 */
	private void spawnThreads(FileFormat ffin){

		//Do anything necessary prior to processing
		ByteFile bf=ByteFile.makeByteFile(ffin);

		//Determine how many threads may be used
		final int threads=Shared.threads();

		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(bf, i));
		}

		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);

		//Do anything necessary after processing
		boolean bfErr=bf.close();
		errorState|=!success || bfErr;

	}

	/**
	 * Adds one completed worker's nonempty-line counters to the shared totals.
	 * @param pt completed worker containing counters to accumulate
	 */
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt){
			linesProcessed+=pt.linesProcessedT;
			bytesProcessed+=pt.bytesProcessedT;
			linesOut+=pt.linesOutT;
			bytesOut+=pt.bytesOutT;
			errorState|=(!pt.success);
			errorState|=(pt.errorStateT);
		}
	}

	/**
	 * Reports aggregate prepass, worker, reader and writer success to ThreadWaiter.
	 * @return true if no error has been recorded
	 */
	@Override
	public final boolean success(){return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Processes the first input on the main thread, counting nonempty lines and
	 * routing valid prefixes to optional output. Invalid-line output is dormant:
	 * outInvalid is private, always null and never parsed or assigned.
	 * @param ff input descriptor for the prepass
	 */
	private void processSinglethreaded(FileFormat ff){
		ByteFile bf=ByteFile.makeByteFile(ff);
		ByteStreamWriter bsw=makeBSW(ffout1);
		ByteStreamWriter bswInvalid=makeBSW(ffoutInvalid);

		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();

		while(line!=null){
			if(line.length>0){
				linesProcessed++;
				bytesProcessed+=(line.length+1);

				final boolean valid=(line[0]!='#');

				if(valid){
					linesOut++;
					bytesOut+=(line.length+1);
					for(int i=0; i<line.length && line[i]!='\t'; i++){bb.append(line[i]);}
					if(bsw!=null){bsw.print(bb.nl().toBytes());}
					bb.clear();
				}else{
					if(bswInvalid!=null){bswInvalid.println(line);}
				}
			}
			line=bf.nextLine();
		}

		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Worker that pulls line batches from a shared ByteFile and parses valid
	 * non-comment lines.
	 */
	static class ProcessThread extends Thread{

		//Constructor
		ProcessThread(final ByteFile bf_, final int tid_){
			bf=bf_;
			tid=tid_;
		}

		/** Pulls and processes lines until the shared ByteFile is exhausted. */
		@Override
		public synchronized void run(){
			//Do anything necessary prior to processing

			//Process the reads
			processInner();

			//Do anything necessary after processing

			//Indicate successful exit status
			success=true;
		}

		/**
		 * Fetches byte-line batches, counts nonempty and valid lines, and
		 * delegates valid lines to processLine().
		 */
		void processInner(){

			//bf.nextList() is synchronized
			for(ListNum<byte[]> ln=bf.nextList(); ln!=null; ln=bf.nextList()){
				for(byte[] line : ln){
					if(line!=null && line.length>0){
						linesProcessedT++;
						bytesProcessedT+=(line.length+1);

						final boolean valid=(line[0]!='#');
						if(valid){
							linesOutT++;
							bytesOutT+=(line.length+1);
							processLine(line);
						}
					}
				}
			}
		}

		/**
		 * Placeholder valid-line processing hook; replace before using this
		 * template.
		 * @param line tab-delimited input line
		 * @return first parsed integer, as a placeholder result
		 */
		int processLine(final byte[] line){
			lp.set(line); //Zero-allocation parsing
			int value=lp.parseInt(0);
			assert(false) : KillSwitch.assertDie("TODO: Implement this method.");//TODO
			return value;
		}

		/** Number of nonempty input lines processed by this thread. */
		protected long linesProcessedT=0;
		/** Input bytes processed by this thread, with one newline byte counted per nonempty line. */
		protected long bytesProcessedT=0;

		/** Number of processed input lines that did not start with '#'. */
		protected long linesOutT=0;
		/** Full valid input-line bytes, with one newline byte counted per valid line. */
		protected long bytesOutT=0;

		/** True if this worker hit a local error. */
		protected boolean errorStateT=false;

		/** True after normal completion. */
		boolean success=false;

		/** Reused tab-delimited parser. */
		private final LineParser1 lp=new LineParser1(delimiter);
		/** Shared input reader owned by spawnThreads. */
		private final ByteFile bf;
		/** Thread identifier for diagnostics. */
		final int tid;
	}

	/**
	 * Creates and starts a line writer.
	 * @param ff output descriptor
	 * @return started writer, or null when output is disabled
	 */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** First input text file, processed on the main thread. */
	private String in1=null;

	/** Second input text file, processed by worker threads. */
	private String in2=null;

	/** Valid-line output file for the first input prepass. */
	private String out1=null;

	/** Dormant invalid-line output path; currently never assigned. */
	private String outInvalid=null;

	/*--------------------------------------------------------------*/

	/** Total nonempty input lines processed across both inputs. */
	private long linesProcessed=0;
	/** Total processed input lines that did not start with '#'. */
	private long linesOut=0;
	/** Total input bytes processed, with one newline byte counted per nonempty line. */
	private long bytesProcessed=0;
	/** Total full valid input-line bytes, with one newline byte counted per valid line. */
	private long bytesOut=0;

	/** Delimiter for worker line parsing. */
	private static byte delimiter='\t';

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Validated first input descriptor. */
	private final FileFormat ffin1;

	/** Validated second input descriptor. */
	private final FileFormat ffin2;
	/** Validated valid-line output descriptor, or null when output is disabled. */
	private final FileFormat ffout1;
	/** Validated invalid-line output descriptor, or null when output is disabled. */
	private final FileFormat ffoutInvalid;

	/**
	 * Supplies ThreadWaiter with this accumulator's lock.
	 * @return read-write lock for worker accumulation
	 */
	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Diagnostic and status output stream, possibly redirected by PreParser. */
	private PrintStream outstream=System.err;
	/** Print verbose diagnostic messages. */
	public static boolean verbose=false;
	/** True if any worker, reader or writer reported an error. */
	public boolean errorState=false;
	/** Overwrite existing output files. */
	private boolean overwrite=true;
	/** Append to existing output files. */
	private boolean append=false;

}
