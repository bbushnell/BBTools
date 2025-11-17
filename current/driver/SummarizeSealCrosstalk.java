package driver;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import tracker.SealStats;
import tracker.SealStats.SealStatsLine;

/**
 * @author Brian Bushnell
 * @date August 15, 2023
 *
 */
public class SummarizeSealCrosstalk {

	/**
	 * Program entry point for SEAL crosstalk analysis.
	 * Creates instance, processes files, and manages output stream cleanup.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		SummarizeSealCrosstalk x=new SummarizeSealCrosstalk(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and initializes processing parameters.
	 * Handles file inputs, output settings, and primary sequence usage configuration.
	 * @param args Command-line arguments to parse
	 */
	public SummarizeSealCrosstalk(String[] args){
		
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
			}else if(a.equals("primary") || a.equals("useprimary")){
				usePrimary=Parse.parseBoolean(b);
			}else if(arg.indexOf('=')<0 && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();

			overwrite=parser.overwrite;
			append=parser.append;
			if(append) {overwrite=false;}
			if(parser.in1!=null){
				Tools.addFiles(parser.in1, in);
			}
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
	}
	
	/**
	 * Main processing method that analyzes SEAL statistics files for crosstalk.
	 * Reads SealStats from input files, calculates contamination rates, and outputs
	 * summary results in PPM format sorted by contamination level.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){

		ArrayList<byte[]> header=new ArrayList<byte[]>();
//		ArrayList<byte[]> results=new ArrayList<byte[]>();
		ArrayList<Result> results=new ArrayList<Result>();
		ByteBuilder bb=new ByteBuilder();

		bb.append("#All values are in read PPM.\n");
		bb.append("#Files");
		for(String s : in) {bb.tab().append(s);}
		//header.add(bb.nl().toBytes());
		bb.clear();
		
		Result totals=new Result("#Totals");
		
		for(String fname : in) {
			SealStats ss=new SealStats(fname);
			String core=ss.fnamePrefix();
			core=core.replace("stats_", "").replace("_stats", "");
			SealStatsLine primary=(usePrimary ? ss.primary() : ss.map.get(core));
			String name=(usePrimary && primary!=null ? primary.name : core);
			SealStatsLine contam=ss.countNonmatching(name);

			long primaryReads=(primary==null ? 0 : primary.reads);
			long contamReads=(contam==null ? 0 : contam.reads);
			Result r=new Result(core, ss.totalReads, ss.matchedReads, primaryReads, contamReads);
			totals.add(r);
			results.add(r);
		}
		
		if(verbose){outstream.println("Finished reading data.");}
		
		{
			totals.appendTo(bb);
			header.add(bb.nl().toBytes());
			bb.clear();
			
			bb.append("#Name\tCorrect\tContam\tAmbig");
			header.add(bb.nl().toBytes());
			bb.clear();
		}
		
		outputResults(header, results);
		
		t.stop();
		outstream.println("Time:                         \t"+t);
//		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
//		assert(!errorState) : "An error was encountered.";
	}
	
	/**
	 * Writes header and sorted results to the output stream.
	 * Results are sorted by contamination rate and formatted for tab-delimited output.
	 * @param header Header lines to write before results
	 * @param results Result objects to sort and output
	 */
	private void outputResults(ArrayList<byte[]> header, ArrayList<Result> results){
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();

		ByteBuilder bb=new ByteBuilder();
		for(byte[] line : header) {bsw.print(line);}
		
		Collections.sort(results);
//		Collections.reverse(results);
		for(Result r : results) {
			r.appendTo(bb).nl();
			bsw.print(bb);
			bb.clear();
		}

		errorState=bsw.poisonAndWait() | errorState;
	}
	
	/**
	 * Represents contamination analysis results for a single sample.
	 * Stores read counts for total, matched, primary, contamination, and ambiguous reads.
	 * Implements comparison based on contamination rates for sorting output.
	 */
	private static class Result implements Comparable<Result> {
		
		/** Creates Result with name only, initializing all counts to zero.
		 * @param name_ Sample name identifier */
		Result(String name_){
			this(name_, 0, 0, 0, 0);
		}
		
		/**
		 * Creates Result with specified read counts and calculates ambiguous reads.
		 * Ambiguous reads are calculated as total minus matched reads.
		 *
		 * @param name_ Sample name identifier
		 * @param total_ Total number of reads processed
		 * @param matched_ Number of reads that matched any reference
		 * @param primary_ Number of reads matching the primary reference
		 * @param contam_ Number of contaminating reads (non-primary matches)
		 */
		Result(String name_, long total_, long matched_, long primary_, long contam_){
			name=name_;
			total=total_;
			matched=matched_;
			primary=primary_;
			contam=contam_;
			ambig=total-matched;
		}
		
		/**
		 * Adds another Result's counts to this Result for totaling across samples.
		 * All count fields are summed together.
		 * @param r Result to add to this one
		 * @return this Result with updated counts
		 */
		Result add(Result r) {
			total+=r.total;
			matched+=r.matched;
			primary+=r.primary;
			contam+=r.contam;
			ambig+=r.ambig;
			return this;
		}
		
		/**
		 * Appends tab-delimited result data to ByteBuilder in PPM format.
		 * Converts read counts to parts-per-million for standardized comparison.
		 * Output format: name, primary_ppm, contamination_ppm, ambiguous_ppm.
		 *
		 * @param bb ByteBuilder to append formatted data to
		 * @return the same ByteBuilder for method chaining
		 */
		ByteBuilder appendTo(ByteBuilder bb) {
			double inv=1.0/(Tools.max(1, total));
			double ppmMult=inv*1000000;
			bb.append(name).tab();//Could use name here instead
			bb.append(primary*ppmMult, 2).tab();
			bb.append(contam*ppmMult, 2).tab();
			bb.append(ambig*ppmMult, 2);
			return bb;
		}

		@Override
		public int compareTo(Result o) {
			double inv1=1.0/(Tools.max(1, total));
			double inv2=1.0/(Tools.max(1, o.total));
			if(contam*inv1!=o.contam*inv2) {return contam*inv1>o.contam*inv2 ? 1 : -1;}
			if(primary*inv1!=o.primary*inv2) {return primary*inv1>o.primary*inv2 ? 1 : -1;}
			if(ambig*inv1!=o.ambig*inv2) {return ambig*inv1>o.ambig*inv2 ? 1 : -1;}
			return name.compareTo(o.name);
		}
		
		/** Sample name identifier */
		String name;
		/** Total number of reads processed for this sample */
		long total;
		/** Number of reads that matched any reference sequence */
		long matched;
		/** Number of reads matching the primary/expected reference */
		long primary;
		/** Number of contaminating reads (matched non-primary references) */
		long contam;
		/** Number of ambiguous reads (total - matched) */
		long ambig;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** List of input SEAL statistics files to process */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output file path for summary results */
	private String out1="stdout.txt";
	
	/** Output file format handler */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	
	/** Tracks whether any errors occurred during processing */
	private boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
	/**
	 * Whether to use primary reference from SEAL stats instead of filename-based matching
	 */
	private boolean usePrimary=false;
	
	/*--------------------------------------------------------------*/
	
	/** Print stream for status messages and logging */
	private java.io.PrintStream outstream=System.err;
	/** Controls verbose output for debugging and progress tracking */
	public static boolean verbose=false;
	
}
