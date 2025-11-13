package stream;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile4;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Tests and benchmarks the Streamer and Writer interfaces.
 * 
 * Reads data from input files using StreamerFactory,
 * optionally processes reads/SamLines,
 * and writes output using WriterFactory.
 * 
 * Supports:
 * - FASTQ, FASTA, SAM, and BAM input/output
 * - Paired and interleaved files
 * - Subsampling with samplerate
 * - Multithreaded streaming and writing
 * - Ordered or unordered output
 * 
 * Usage examples:
 *   StreamerWrapper in=reads.fq.gz out=filtered.fq.gz samplerate=0.1
 *   StreamerWrapper in=mapped.bam out=reads.fq.gz
 *   StreamerWrapper in1=r1.fq in2=r2.fq out=interleaved.fq
 * 
 * @author Brian Bushnell, Isla
 * @date November 4, 2025
 */
public class StreamerWrapper{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args) {
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		StreamerWrapper x=new StreamerWrapper(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public StreamerWrapper(String[] args){
		
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
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			out1=parser.out1;
			out2=parser.out2;
		}

		doPoundReplacement();
		fixExtensions();
		adjustInterleaving();
		checkFileExistence();
		checkStatics();
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, true);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, true);
		
		final boolean samIn=(ffin1!=null && ffin1.samOrBam());
		final boolean samOut=(ffout1!=null && ffout1.samOrBam());
		SamLine.SET_FROM_OK=samIn;
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=samIn && samOut;
		//Determine if we need to parse SAM fields or can skip for performance
		if(!forceParse && samIn && !samOut){
			SamLine.PARSE_2=false;
			SamLine.PARSE_5=false;
			SamLine.PARSE_6=false;
			SamLine.PARSE_7=false;
			SamLine.PARSE_8=false;
			SamLine.PARSE_OPTIONAL=false;
		}
	}

	/** 
	 * Parse arguments from the command line.
	 * @param args Command line arguments
	 * @return Parser object with standard flags processed
	 */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
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
			}else if(a.equals("threadsin") || a.equals("tin")){
				threadsIn=Integer.parseInt(b);
			}else if(a.equals("threadsout") || a.equals("tout")){
				threadsOut=Integer.parseInt(b);
			}else if(a.equals("samplerate") || a.equals("sample")){
				samplerate=Float.parseFloat(b);
			}else if(a.equals("sampleseed") || a.equals("seed")){
				sampleseed=Long.parseLong(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("forceparse")){
				forceParse=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(i==0 && !arg.contains("=") && parser.in1==null &&
				FileFormat.isSequence(arg) && new File(arg).isFile()){
				parser.in1=arg;
			}else if(i==1 && !arg.contains("=") && parser.out1==null && parser.in1!=null &&
				FileFormat.isSequence(arg)){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	

	/** 
	 * Replace # with 1 and 2 in file names.
	 * Example: reads_#.fq becomes reads_1.fq and reads_2.fq
	 */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}
	
	/** Add or remove .gz or .bz2 extensions as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}
	
	/** Ensure input files can be read and output files can be written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** 
	 * Adjust interleaved mode based on number of input and output files.
	 * Two input files forces non-interleaved mode.
	 * Two output files with one input forces interleaved mode.
	 */
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Empty
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Primary Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Create Streamer and Writer, then process all data.
	 * @param t Timer for tracking elapsed time
	 */
	private void process(Timer t) {
		final boolean inputReads=(ffin1!=null && !ffin1.samOrBam());
		final boolean inputSam=(ffin1!=null && ffin1.samOrBam());
		final boolean outputReads=(ffout1!=null && !ffout1.samOrBam());
		final boolean outputSam=(ffout1!=null && ffout1.samOrBam());
		final boolean saveHeader=inputSam && outputSam;
		
		Streamer st=StreamerFactory.makeStreamer(ffin1, ffin2, ordered, maxReads, 
			saveHeader, outputReads, threadsIn);
		st.setSampleRate(samplerate, sampleseed);
		Writer fw=WriterFactory.makeWriter(ffout1, ffout2, threadsOut, null, saveHeader);
		
		process(st, fw, t, inputReads || outputReads);
	}
	
	/**
	 * Main processing loop - reads from Streamer, processes, writes to Writer.
	 * @param st Input Streamer
	 * @param fw Output Writer (may be null)
	 * @param t Timer for tracking elapsed time
	 * @param readMode True for Read objects, false for SamLine objects
	 */
	private void process(Streamer st, Writer fw, Timer t, boolean readMode) {
		if(threadsIn==0 || threadsIn==1) {Read.VALIDATE_IN_CONSTRUCTOR=false;}
//		System.err.println("VIC="+Read.VALIDATE_IN_CONSTRUCTOR);
		st.start();
		if(fw!=null) {fw.start();}
		if(readMode) {
			for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()) {
				for(Read r : ln) {
					processReadPair(r, r.mate);
				}
				if(fw!=null) {fw.addReads(ln);}
			}
		}else {
			for(ListNum<SamLine> ln=st.nextLines(); ln!=null; ln=st.nextLines()) {
				final ArrayList<SamLine> list=ln.list;
				for(int i=0, len=list.size(); i<len; i++) {
					SamLine sl=list.get(i);
					boolean keep=processSamLine(sl);
					if(!keep) {list.set(i, null);}
				}
				if(fw!=null) {fw.addLines(ln);}
			}
		}
		if(fw!=null) {
			fw.poisonAndWait();
			assert(!readMode || readsIn==fw.readsWritten());
			assert(!readMode || basesIn==fw.basesWritten());
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}
		t.stop();
		System.err.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		st.close();//Prevents a BF4 hang with limited reads
//		try {
//		    Thread.sleep(100);
//		} catch (InterruptedException e) {
//		    // Ignore
//		}
//		Shared.listThreads2();
	}
	
	/**
	 * Process a read pair - override this method to add custom processing.
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 */
	private void processReadPair(Read r1, Read r2) {
		readsIn+=r1.pairCount();
		basesIn+=r1.pairLength();
		if(!r1.validated()) {r1.validate(true);}
		if(r2!=null && !r2.validated()) {r2.validate(true);}
	}
	
	/**
	 * Process a SamLine - override this method to add custom processing.
	 * @param sl SamLine to process
	 * @return keep;
	 */
	private boolean processSamLine(SamLine sl) {
		final int len=sl.lengthOrZero();
		readsIn++;
		basesIn+=len;

		//Apply filter if present
		boolean keep=(len>0 || ffout1==null || ffout1.samOrBam());
		return keep;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;
	
	/** Primary input file format */
	private FileFormat ffin1;
	/** Secondary input file format */
	private FileFormat ffin2;
	/** Primary output file format */
	private FileFormat ffout1;
	/** Secondary output file format */
	private FileFormat ffout2;
	
	/** Number of threads for input streaming (-1 = auto) */
	private int threadsIn=-1;
	/** Number of threads for output writing (-1 = auto) */
	private int threadsOut=-1;
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	/** Fraction of reads to keep (1.0 = all reads) */
	private float samplerate=1f;
	/** Random seed for subsampling */
	private long sampleseed=17;
	/** Force parsing of all SAM fields even if not needed */
	private boolean forceParse=false;

	/** Number of reads processed */
	private long readsIn=0;
	/** Number of bases processed */
	private long basesIn=0;
	/** Number of reads written */
	private long readsOut=0;
	/** Number of bases written */
	private long basesOut=0;
	
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	/** Maintain input order in output */
	private boolean ordered=true;
	/** Whether interleaved was explicitly set */
	private boolean setInterleaved=false;
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print verbose messages */
	public static boolean verbose=false;
	
}