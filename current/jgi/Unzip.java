package jgi;

import java.io.IOException;
import java.io.InputStream;
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

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class Unzip {
	
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
		Unzip x=new Unzip(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Unzip(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
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

			out1=parser.out1;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		parser.overwrite=true;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
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
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
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
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that executes the decompression operation.
	 * Opens input stream, creates output writer, and processes data.
	 * Handles cleanup, error reporting, and performance statistics.
	 * @param t Timer for tracking execution performance
	 */
	void process(Timer t){
		
		InputStream is=ReadWrite.getInputStream(ffin1.name(), false, true);
		ByteStreamWriter bsw=makeBSW(ffout1);
		
//		assert(false) : "Header goes here.";
		if(bsw!=null){
//			assert(false) : "Header goes here.";
		}
		
		processInner(is, bsw);
		
		try {
			is.close();
		} catch (IOException e) {
			errorState=true;
			e.printStackTrace();
		}
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		
		if(showSpeed){
			outstream.println(Tools.timeLinesBytesProcessed(t, 1, bytesProcessed, 8));

			outstream.println();
			outstream.println("Bytes Processed:   \t"+bytesProcessed);
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Core decompression loop that streams data from input to output.
	 * Uses buffered reading for efficient processing of large files.
	 * Tracks bytes processed and handles I/O exceptions.
	 *
	 * @param is Input stream containing compressed data
	 * @param bsw Output writer for decompressed data
	 */
	private void processInner(InputStream is, ByteStreamWriter bsw){
		final byte[] buffer=new byte[65536<<2];
		
		int len=0;
		try {
			len = is.read(buffer);
		} catch (IOException e) {
			errorState=true;
			e.printStackTrace();
		}
		while(len>0){
			bytesProcessed+=len;
			bsw.print(buffer, len);
			len=0;
			try {
				len=is.read(buffer);
			} catch (IOException e) {
				errorState=true;
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Creates and initializes a ByteStreamWriter for output.
	 * Returns null if no file format is specified.
	 * @param ff File format specification for output
	 * @return Initialized ByteStreamWriter or null
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
	
	/** Primary input file path for decompression */
	private String in1=null;
	/** Primary output file path for decompressed data */
	private String out1=null;
	/** Output file path for invalid or corrupted data */
	private String outInvalid=null;
	
	/*--------------------------------------------------------------*/
	
	/** Total number of bytes processed during decompression */
	private long bytesProcessed=0;
	
	/** Maximum number of lines to process before stopping */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** File format configuration for primary input */
	private final FileFormat ffin1;
	/** File format configuration for primary output */
	private final FileFormat ffout1;
	/** File format configuration for invalid data output */
	private final FileFormat ffoutInvalid;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and error reporting */
	private PrintStream outstream=System.err;
	/** Enables verbose logging and detailed progress reporting */
	public static boolean verbose=false;
	/** Tracks whether any errors occurred during processing */
	public boolean errorState=false;
	/** Controls display of processing speed and performance statistics */
	public boolean showSpeed=false;
	/** Allows overwriting existing output files */
	private boolean overwrite=true;
	/** Appends output to existing files rather than overwriting */
	private boolean append=false;
	
}
