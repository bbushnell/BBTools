package fileIO;

import java.io.File;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import structures.ListNum;

/**
 * Accepts multiple input files.
 * Reads them each sequentially, and outputs everything to a single output file.
 * Generically, it can be used to concatenate files while recompressing them
 * and avoiding the use of stdio.
 * @author Brian Bushnell
 * @date January 21, 2025
 *
 */
public class Concatenate2 {

	/**
	 * Program entry point for file concatenation.
	 * Creates a timer, instantiates Concatenate, processes files, and closes streams.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		Concatenate2 x=new Concatenate2(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs a Concatenate instance and parses command-line arguments.
	 * Handles input file specification, output file configuration, and parser setup.
	 * Supports multiple input files via 'in' parameter or direct file arguments.
	 * @param args Command-line arguments including input/output file specifications
	 */
	public Concatenate2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		parser.out1=out1;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("in")){
				in.clear();
				String[] b2=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
				for(String b3 : b2){in.add(b3);}
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
	}
	
	/**
	 * Main processing method that concatenates all input files to output.
	 * Creates output stream writer, validates file names don't conflict,
	 * and processes each input file sequentially.
	 * @param t Timer for tracking execution time and performance metrics
	 */
	void process(Timer t){

		final FileFormat ffout=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, true, false, true);
		if(out1!=null){
			for(String s : in) {
				assert(!out1.equalsIgnoreCase(s)) : "Input file and output file have same name.";
			}
		}
		final OutputStream os=(ffout!=null ? ReadWrite.getOutputStream(out1, false, false, true) : null);
		
		for(String s : in) {
			processInner(s, os);
		}
		
		if(out1!=null) {ReadWrite.finishWriting(null, os, out1, ffout.subprocess);}
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		if(verbose) {
			outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 12));
			outstream.println(bytesProcessed);
		}
	}
	
	/**
	 * Processes a single input file by reading all lines and writing to output stream.
	 *
	 * @param fname Input filename to process
	 * @param os OutputStream for output (may be null for stdout)
	 */
	void processInner(String fname, OutputStream os) {
		final FileFormat ffin=FileFormat.testInput(fname, FileFormat.TXT, null, true, true);
		final byte[] buffer=new byte[1024*256];
		final InputStream is=ReadWrite.getInputStream(fname, false, true);
		
		try {
			for(int r=is.read(buffer); r>0; r=is.read(buffer)) {
				os.write(buffer, 0, r);
				bytesProcessed+=r;
				linesProcessed=Vector.countSymbols(buffer, 0, r, (byte)'\n');
			}
		}catch(Exception e) {
			KillSwitch.exceptionKill(e);
		}
		ReadWrite.finishReading(is, fname, ffin.subprocess);
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** List of input filenames to concatenate */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output filename for concatenated result (default: stdout.txt) */
	private String out1="stdout.txt";
	
	/** FileFormat object for output file format detection and handling */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Counter for total lines processed across all input files */
	private long linesProcessed=0, bytesProcessed=0;
	
	/*--------------------------------------------------------------*/
	
	/** Print stream for status messages and error output (default: stderr) */
	private java.io.PrintStream outstream=System.err;
	/** Controls verbose output messages during processing */
	public static boolean verbose=false;
	
}
