package fileIO;

import java.io.File;
import java.util.ArrayList;

import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
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
public class Concatenate {

	/**
	 * Program entry point for file concatenation.
	 * Creates a timer, instantiates Concatenate, processes files, and closes streams.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		Concatenate x=new Concatenate(args);
		
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
	public Concatenate(String[] args){
		
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

		final ByteStreamWriter bsw;
		if(out1!=null){
			final int buff=4;

			for(String s : in) {
				assert(!out1.equalsIgnoreCase(s)) : "Input file and output file have same name.";
			}
			
			bsw=ByteStreamWriter.makeBSW(ffout1);
		}else{bsw=null;}
		
		for(String s : in) {
			processInner(s, bsw);
		}
		
		bsw.poisonAndWait();
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		if(verbose) {
			outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 12));
			outstream.println(bytesProcessed);
		}
	}
	
	/**
	 * Processes a single input file by reading all lines and writing to output stream.
	 * Creates ByteFile from input filename, reads line batches, and tracks statistics.
	 * Updates line and byte counters for performance reporting.
	 *
	 * @param fname Input filename to process
	 * @param bsw ByteStreamWriter for output (may be null for stdout)
	 */
	void processInner(String fname, ByteStreamWriter bsw) {
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.TXT, null, true, true);
		
		final ByteFile bf=ByteFile.makeByteFile(ffin);
		
		ListNum<byte[]> ln;
		
		while((ln=bf.nextList())!=null) {
			for(byte[] line : ln.list) {
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				if(bsw!=null) {bsw.println(line);}
			}
		}
		bf.close();
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
