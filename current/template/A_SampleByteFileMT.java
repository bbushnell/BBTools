package template;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Multithreaded text file processing template with file input/output handling.
 * Loads text files, spawns worker threads to process input files, and manages
 * concurrent file reading. Supports multiple input/output file formats with configurable
 * thread management for file processing.
 *
 * @author Brian Bushnell
 * @date February 6, 2023
 */
public class A_SampleByteFileMT implements Accumulator<A_SampleByteFileMT.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point that creates an instance and processes input files.
	 * @param args Command line arguments */
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
	 * Constructor that initializes the processor from command line arguments.
	 * Parses arguments, validates parameters, and sets up input/output file formats.
	 * @param args Command line arguments containing input/output paths and options
	 */
	public A_SampleByteFileMT(String[] args){
		
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
	 * Parses command line arguments into configuration parameters.
	 * Supports verbose mode and line processing limits.
	 *
	 * @param args Command line arguments array
	 * @return Configured Parser object with parsed settings
	 * @throws RuntimeException If unknown parameters are encountered
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
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
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
	 * Adjusts file extensions by adding or removing compression suffixes as needed.
	 * Validates that at least one input file is specified.
	 * @throws RuntimeException If no input file is provided
	 */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/**
	 * Validates that input files can be read and output files can be written.
	 * Checks for duplicate file specifications and proper file permissions.
	 * @throws RuntimeException If files cannot be accessed or duplicates exist
	 */
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
	
	/** Adjusts static file-related settings for optimal performance.
	 * Forces ByteFile mode BF2 when using more than 2 threads. */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/**
	 * Validates parameter ranges and ensures required parameters are set.
	 * Currently contains placeholder assertion that needs implementation.
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
	 * Main processing method that orchestrates file reading and thread spawning.
	 * Processes primary input file sequentially, then spawns threads for secondary file.
	 * Reports timing and processing statistics upon completion.
	 *
	 * @param t Timer for tracking execution time
	 * @throws RuntimeException If processing encounters errors
	 */
	void process(Timer t){
		
		//Reset counters
		linesProcessed=linesOut=0;
		bytesProcessed=bytesOut=0;
		
		processFF1(ffin1);
		
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
	 * Creates and manages worker threads for concurrent file processing.
	 * Spawns threads equal to system thread count and waits for completion.
	 * @param ffin Input file format to be processed by threads
	 */
	private void spawnThreads(FileFormat ffin){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(ffin, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	/**
	 * Accumulates processing statistics from completed worker threads.
	 * Thread-safe method that combines counters from individual threads.
	 * @param pt ProcessThread containing processing statistics to accumulate
	 */
	@Override
	public final void accumulate(ProcessThread pt){
		linesProcessed+=pt.linesProcessedT;
		bytesProcessed+=pt.bytesProcessedT;
		linesOut+=pt.linesOutT;
		bytesOut+=pt.bytesOutT;
		errorState|=(!pt.success);
		errorState|=(pt.errorStateT);
	}
	
	/** Reports whether processing completed successfully.
	 * @return true if no errors were encountered during processing */
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void processFF1(FileFormat ff){
		ByteFile bf=ByteFile.makeByteFile(ff);
		ByteStreamWriter bsw=makeBSW(ffout1);
		ByteStreamWriter bswInvalid=makeBSW(ffoutInvalid);
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				final boolean valid=(line[0]!='#');
				
				if(valid){
					linesOut++;
					bytesOut+=(line.length+1);
					for(int i=0; i<line.length && line[i]!='\t'; i++){
						bb.append(line[i]);
					}
					bb.nl();
					bsw.print(bb.toBytes());
					bb.clear();
				}else{
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
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
	 * Worker thread for processing input files concurrently.
	 * Static class to prevent accidental access to shared variables.
	 * Each thread processes lines from the input file independently.
	 */
	static class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final FileFormat ffin_, final int tid_){
			ffin=ffin_;
			tid=tid_;
		}
		
		//Called by start()
		/** Thread execution method that processes input and sets success status.
		 * Calls processInner() to handle the actual file processing work. */
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
		 * Iterates through input file lines and processes valid entries.
		 * Reads lines from ByteFile, applies validation, and delegates to processLine().
		 * Respects maxLines limit if configured.
		 */
		void processInner(){
			ByteFile bf=ByteFile.makeByteFile(ffin);
			byte[] line=bf.nextLine();
			ByteBuilder bb=new ByteBuilder();
			
			while(line!=null){
				if(line.length>0){
					if(maxLines>0 && linesProcessedT>=maxLines){break;}
					linesProcessedT++;
					bytesProcessedT+=(line.length+1);
					
					final boolean valid=(line[0]!='#');
					
					if(valid){
						processLine(line);
					}else{
						//Do something else
					}
				}
				line=bf.nextLine();
			}
			errorStateT|=bf.close();
		}
		
		/**
		 * Processes a single line from the input file.
		 * Currently contains placeholder assertion that needs implementation.
		 * @param line The input line to process as byte array
		 * @return true if the line should be retained, false if discarded
		 */
		boolean processLine(final byte[] line){
//			throw new RuntimeException("TODO: Implement this method."); //TODO
			assert(false);
			return true;
		}

		/** Number of lines processed by this thread */
		protected long linesProcessedT=0;
		/** Number of bytes processed by this thread */
		protected long bytesProcessedT=0;
		
		/** Number of lines output by this thread */
		protected long linesOutT=0;
		/** Number of bytes output by this thread */
		protected long bytesOutT=0;
		
		protected boolean errorStateT=false;
		
		/** Flag indicating successful completion of this thread */
		boolean success=false;
		
		/** Input file format for this thread */
		private final FileFormat ffin;
		/** Thread identifier for this worker thread */
		final int tid;
	}
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;

	private String in2=null;

	private String out1=null;

	/** Output file path for invalid/filtered lines */
	private String outInvalid=null;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	private static long maxLines=Long.MAX_VALUE;//TODO: Static to compile; make non-static if threads are non-static
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** FileFormat wrapper for primary input file */
	private final FileFormat ffin1;
	
	/** FileFormat wrapper for secondary input file */
	private final FileFormat ffin2;
	/** FileFormat wrapper for primary output file */
	private final FileFormat ffout1;
	/** FileFormat wrapper for invalid line output file */
	private final FileFormat ffoutInvalid;
	
	/** Returns the read-write lock for thread synchronization.
	 * @return ReadWriteLock instance for coordinating thread access */
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages */
	private PrintStream outstream=System.err;
	/** Flag to enable verbose output messages */
	public static boolean verbose=false;
	/** Flag indicating whether an error was encountered during processing */
	public boolean errorState=false;
	/** Flag to allow overwriting existing output files */
	private boolean overwrite=true;
	/** Flag to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
