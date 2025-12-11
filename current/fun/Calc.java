package fun;

import java.io.PrintStream;

import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tracker.ReadStats;

/**
 * Command-line tool to compute cumulative distributions of bit-encoded integer combinations.
 * Iterates over 2^(numStats*5) combinations, counts bucket sums, and reports cumulative percentages.
 * @author Brian Bushnell
 */
public class Calc {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Entry point: parses args, runs processing, closes redirected streams.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		Calc x=new Calc(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and configures processing
	 * parameters. Supports verbose mode, numstats configuration, and standard
	 * parser flags for output file handling.
	 * @param args Command-line arguments to parse
	 */
	public Calc(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(a.equals("num") || a.equals("numstats")){
				numStats=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out1=parser.out1;
		}
	}
	

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing method that executes statistical calculation pipeline.
	 * Calls processInner for computation, reports timing information, and
	 * checks for error states.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		//Process the read stream
		processInner(numStats);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Report timing and results
		{
			t.stop();
			outstream.println("Time: \t"+t);
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Core statistical calculation engine that generates cumulative distribution
	 * statistics. Iterates through 2^(numStats*5) combinations, calculates
	 * frequency distribution using sum() method, builds cumulative counts,
	 * and outputs percentage statistics to stdout.
	 *
	 * @param numStats Number of statistics to process (affects bit combinations)
	 */
	void processInner(int numStats){
		int bits=numStats*5;
		final int iters=1<<bits;
		final int buckets=1+31*numStats;
		int[] counts=new int[buckets];
		for(int i=0; i<iters; i++){
			counts[sum(i)]++;
		}
		int[] cumulative=new int[buckets];
		cumulative[0]=counts[0];
		for(int i=1; i<buckets; i++){
			cumulative[i]=cumulative[i-1]+counts[i];
		}
		//StringBuilder sb=new StringBuilder();
		final double mult=100.0/iters;
		for(int i=0; i<buckets; i++){
			String s=(Tools.format("%d\t%.4f%%\n", i, cumulative[i]*mult));
			System.out.print(s);
		}
	}
	
	/**
	 * Calculates the sum of 5-bit chunks (mask 0x1F) in an integer for bucket assignment.
	 * @param stats Bit-encoded integer
	 * @return Sum of 5-bit segments
	 */
	int sum(int stats){
		int sum=0;
		while(stats>0){
			sum+=(stats&0x1F);
			stats>>>=5;
		}
		return sum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output path (default stdout.txt). */
	private String out1="stdout.txt";
	
	/** Number of statistics to process (controls 2^(numStats*5) iterations). */
	private int numStats=6;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print stream for status messages and verbose output */
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	/** Controls whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Controls whether to append to existing output files */
	private boolean append=false;
	
}
