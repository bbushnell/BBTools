package jgi;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Generates some stats from multiple files.
 * Uses the new Assembly class.
 * @author Brian Bushnell
 * @date January 21, 2025
 *
 */
public class AssemblyStats3 {

	/**
	 * Program entry point. Creates AssemblyStats3 instance and executes processing.
	 * Manages timer initialization and output stream cleanup.
	 * @param args Command-line arguments for assembly statistics generation
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		AssemblyStats3 x=new AssemblyStats3(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs AssemblyStats3 instance with command-line argument parsing.
	 * Processes input files, output configuration, and parser settings.
	 * Supports single files or comma-separated file lists as input.
	 * @param args Command-line arguments including input files and parameters
	 */
	public AssemblyStats3(String[] args){
		
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
	 * Main processing method that generates assembly statistics for all input files.
	 * Creates output writer, processes each input file, and writes tabular results.
	 * Tracks processing statistics and execution timing.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){

		final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout1);
		if(bsw!=null) {bsw.println(makeHeader());}
		
		for(String s : in) {
			processInner(s, bsw);
		}
		
		if(bsw!=null) {bsw.poisonAndWait();}
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 12);
	}
	
	/**
	 * Creates tab-delimited header row for assembly statistics output.
	 * Includes columns for filename, size, contigs, GC content, max contig,
	 * and length thresholds (5k+, 10k+, 25k+, 50k+).
	 * @return Formatted header string with tab separators
	 */
	public static String makeHeader() {
		ByteBuilder bb=new ByteBuilder();
		bb.append("fname");
		bb.tab().append("size");
		bb.tab().append("contigs");
		bb.tab().append("gc");
		bb.tab().append("maxContig");
		bb.tab().append("5kplus");
		bb.tab().append("10kplus");
		bb.tab().append("25kplus");
		bb.tab().append("50kplus");
		return bb.toString();
	}
	
	/**
	 * Processes a single assembly file and writes statistics to output stream.
	 * Creates Assembly object, extracts metrics, and formats results as
	 * tab-delimited row with filename, size, contig count, GC content,
	 * maximum contig length, and length distribution statistics.
	 *
	 * @param fname Assembly filename to process
	 * @param bsw Output stream writer for results
	 */
	void processInner(String fname, ByteStreamWriter bsw) {
		Assembly a=new Assembly(fname);
		if(bsw==null) {return;}
		int max=a.contigs.size>0 ? a.contigs.get(0) : 0;
		bsw.print(fname).tab().print(a.length).tab().print(a.contigs.size);
		bsw.tab().print(a.gc(), 3).tab().print(max);
		bsw.tab().print(a.lengthAtLeast(5000));
		bsw.tab().print(a.lengthAtLeast(10000));
		bsw.tab().print(a.lengthAtLeast(25000));
		bsw.tab().print(a.lengthAtLeast(50000));
		bsw.nl();
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** List of input assembly filenames to process */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output filename for assembly statistics results */
	private String out1="stdout.txt";
	
	/** File format handler for output file configuration */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Counter for total bytes processed across all input files */
	/** Counter for total lines processed across all input files */
	private long linesProcessed=0, bytesProcessed=0;
	
	/*--------------------------------------------------------------*/
	
	/** Output print stream for status messages and logging */
	private java.io.PrintStream outstream=System.err;
	/** Flag to enable verbose output and detailed logging */
	public static boolean verbose=false;
	
}
