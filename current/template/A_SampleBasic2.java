package template;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Reads a text file.
 * Does something.
 * @author Brian Bushnell
 * @date April 9, 2020
 *
 */
public class A_SampleBasic2 {
	
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
		A_SampleBasic2 x=new A_SampleBasic2(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public A_SampleBasic2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		
		{//Parse the arguments
			final Parser parser=parse(args);
			parser.out1="stdout.txt";
			overwrite=parser.overwrite;
			append=parser.append;

			out=parser.out1;
		}
		
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, append, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that coordinates the workflow.
	 * Creates output writer, processes data, handles errors, and reports timing.
	 * @param t Timer for execution time tracking
	 */
	void process(Timer t){
		
		ByteStreamWriter bsw=makeBSW(ffout);
		
//		assert(false) : "Header goes here.";
		if(bsw!=null){
//			assert(false) : "Header goes here.";
		}
		
		processInner(bsw);
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesOut, bytesOut, 8));
		outstream.println();
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Inner processing loop that performs the main work.
	 * Executes cycles of processing using the provided output writer.
	 * @param bsw Output stream writer for results
	 */
	private void processInner(ByteStreamWriter bsw){
		ByteBuilder bb=new ByteBuilder();
		
		for(long cycle=0; cycle<maxCycles; cycle++){
			doSomething(bsw, bb, cycle);
		}
	}
	
	/**
	 * Creates and starts a ByteStreamWriter for the given file format.
	 * @param ff File format specification
	 * @return Started ByteStreamWriter, or null if format is null
	 */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/**
	 * Placeholder method for specific processing logic.
	 * Currently contains commented-out template code for line processing.
	 *
	 * @param bsw Output stream writer
	 * @param bb Byte buffer for building output
	 * @param cycle Current processing cycle number
	 * @return true indicating successful processing
	 */
	private boolean doSomething(ByteStreamWriter bsw, ByteBuilder bb, long cycle){

//		if(line.length>0){
//			linesProcessed++;
//			bytesProcessed+=(line.length+1);
//
//			if(true){
//				linesOut++;
//				bytesOut+=(line.length+1);
//				for(int i=0; i<line.length && line[i]!='\t'; i++){
//					bb.append(line[i]);
//				}
//				bb.nl();
//				bsw.print(bb.toBytes());
//				bb.clear();
//			}
//		}
//		line=bf.nextLine();
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output file path */
	private String out=null;
	/** Maximum number of processing cycles to execute */
	private long maxCycles=100;
	
	/*--------------------------------------------------------------*/
	
	/** Count of output lines written */
	private long linesOut=0;
	/** Count of output bytes written */
	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** File format specification for output */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and errors */
	private PrintStream outstream=System.err;
	/** Enable verbose output and debugging information */
	public static boolean verbose=false;
	/** Tracks whether an error has occurred during processing */
	public boolean errorState=false;
	/** Allow overwriting existing output files */
	private boolean overwrite=true;
	/** Append to existing output files instead of overwriting */
	private boolean append=false;
	
}
