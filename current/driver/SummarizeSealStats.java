package driver;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date May 8, 2015
 *
 */
public class SummarizeSealStats {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Create a new SummarizeSealStats instance
		SummarizeSealStats sss=new SummarizeSealStats(args);
		
		///And run it
		sss.summarize();
	}
	
	/**
	 * Constructs a SummarizeSealStats instance with command-line arguments.
	 * Parses configuration options and input file specifications.
	 * Processes file patterns and expands directories into individual files.
	 * @param args Command-line arguments including input files and options
	 */
	public SummarizeSealStats(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		ArrayList<String> names=new ArrayList<String>();
		Parser parser=new Parser();
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("printtotal") || a.equals("pt")){
				printTotal=Parse.parseBoolean(b);
			}else if(a.equals("ignoresametaxa")){
				ignoreSameTaxa=Parse.parseBoolean(b);
			}else if(a.equals("ignoresamebarcode") || a.equals("ignoresameindex")){
				ignoreSameBarcode=Parse.parseBoolean(b);
			}else if(a.equals("ignoresamelocation") || a.equals("ignoresameloc")){
				ignoreSameLocation=Parse.parseBoolean(b);
			}else if(a.equals("usetotal") || a.equals("totaldenominator") || a.equals("totald") || a.equals("td")){
				totalDenominator=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(!arg.contains("=")){
				String[] x=(new File(arg).exists() ? new String[] {arg} : arg.split(","));
				for(String x2 : x){names.add(x2);}
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		
		{//Process parser fields
			out=(parser.out1==null ? "stdout" : parser.out1);
			if(parser.in1!=null){
				String[] x=(new File(parser.in1).exists() ? new String[] {parser.in1} : parser.in1.split(","));
				for(String x2 : x){names.add(x2);}
			}
		}

		in=new ArrayList<String>();
		for(String s : names){
			Tools.getFileOrFiles(s, in, false, false, false, true);
		}
	}
	
	/**
	 * Executes the main summarization process for all input files.
	 * Creates SealSummary objects for each input file, aggregates totals,
	 * and writes formatted output with primary taxa and contamination statistics.
	 */
	public void summarize(){
		ArrayList<SealSummary> list=new ArrayList<SealSummary>();
		
		SealSummary total=new SealSummary(null);
		total.pname="TOTAL";
		for(String fname : in){
			SealSummary ss=new SealSummary(fname);
			list.add(ss);
			total.add(ss);
		}
		
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		tsw.start();
		tsw.print("#File\tPrimary_Name\tPrimary_Count\tOther_Count\tPrimary_Bases\tOther_Bases\tOther_ppm\n");
		if(printTotal){
			tsw.println(total.toString());
		}
		for(SealSummary ss : list){
			tsw.println(ss.toString());
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Inner class representing statistics summary for a single Seal output file.
	 * Tracks primary and other taxa counts, bases, and calculates contamination metrics.
	 * Supports filtering modes to ignore same taxa, barcode, or location matches.
	 */
	private class SealSummary {
		
		/**
		 * Constructs a SealSummary for the specified file.
		 * Automatically processes the file using appropriate summarization method
		 * based on configured filtering options.
		 * @param fname_ Path to the Seal statistics file to process
		 */
		SealSummary(String fname_){
			fname=fname_;
			if(fname!=null){
				if(ignoreSameTaxa || ignoreSameBarcode || ignoreSameLocation){
					cleanAndSummarize();
				}else{
					summarize();
				}
			}
		}
		
		/**
		 * Aggregates statistics from another SealSummary into this one.
		 * Adds counts and bases, then recalculates contamination ppm.
		 * @param ss The SealSummary to add to this summary
		 */
		public void add(SealSummary ss){
			pcount+=ss.pcount;
			ocount+=ss.ocount;
			tcount+=ss.tcount;
			pbases+=ss.pbases;
			obases+=ss.obases;
			tbases+=ss.tbases;
			
			if(totalDenominator && tbases>0){
				ppm=obases*1000000.0/tbases;
			}else{
				ppm=(obases==0 ? 0 : obases*1000000.0/(obases+pbases));
			}
		}
		
		@Override
		public String toString(){
			return Tools.format("%s\t%s\t%d\t%d\t%d\t%d\t%.2f", fname, pname, pcount, ocount, pbases, obases, ppm);
		}
		
		/**
		 * Basic summarization method that processes Seal statistics file.
		 * Identifies primary taxa (highest base count) and calculates other taxa statistics.
		 * Parses total counts from header lines and individual taxa from data lines.
		 */
		private void summarize(){
			TextFile tf=new TextFile(fname);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(line.startsWith("#")){
					if(line.startsWith("#Total")){
						String[] split=line.split("\t");
						tcount=Long.parseLong(split[1]);
						tbases=Long.parseLong(split[2]);
					}
				}else{
					String[] split=line.split("\t");
					long count=Long.parseLong(split[1]);
					long bases=Long.parseLong(split[3]);
					if(pcount==0 || bases>pbases || (bases==pbases && count>pcount)){
						pname=split[0];
						ocount+=pcount;
						obases+=pbases;
						pcount=count;
						pbases=bases;
					}else{
						ocount+=count;
						obases+=bases;
					}
				}
			}
			tf.close();
			if(totalDenominator && tbases>0){
				ppm=obases*1000000.0/tbases;
			}else{
				ppm=(obases==0 ? 0 : obases*1000000.0/(obases+pbases));
			}
		}
		
		/**
		 * Advanced summarization with filtering for same taxa, barcode, or location.
		 * Parses taxa names and barcodes to apply filtering logic before counting.
		 * Ignores contamination from specified same-category matches based on configuration.
		 */
		public void cleanAndSummarize(){
			TextFile tf=new TextFile(fname);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(line.startsWith("#")){
					if(line.startsWith("#Total")){
						String[] split=line.split("\t");
						tcount=Long.parseLong(split[1]);
						tbases=Long.parseLong(split[2]);
					}
				}else{
					String[] split=line.split("\t");
					String[] name=split[0].toLowerCase().split(",");
					String[] barcode=name[0].split("-");
					
					long count=Long.parseLong(split[1]);
					long bases=Long.parseLong(split[3]);
					if(pcount==0 || bases>pbases || (bases==pbases && count>pcount)){
						name0=name;
						barcode0=barcode;
						pname=split[0];
						ocount+=pcount;
						obases+=pbases;
						pcount=count;
						pbases=bases;
					}else{
						boolean process=true;
						if(ignoreSameTaxa){
							if(name[2].contains(name0[2]) || name0[2].contains(name[2])){
								process=false;
							}
						}
						if(ignoreSameBarcode){
							if(barcode[0].equals(barcode0[0]) || barcode[1].equals(barcode0[1])){
								process=false;
							}
						}
						if(ignoreSameLocation){
							assert(name.length==4) : "Too many delimiters: "+name.length+"\n"+line+"\n";
							if(name[3].equals(name0[3])){
								process=false;
							}
						}
						if(process){
							ocount+=count;
							obases+=bases;
						}
					}
				}
			}
			tf.close();
			if(totalDenominator && tbases>0){
				ppm=obases*1000000.0/tbases;
			}else{
				ppm=(obases==0 ? 0 : obases*1000000.0/(obases+pbases));
			}
		}
		
		/** Path to the input Seal statistics file */
		final String fname;
		/** Name of the primary taxa with highest base count */
		String pname=null;
		/** Total read count from file header */
		/** Combined read count for all other (non-primary) taxa */
		/** Read count for the primary taxa */
		long pcount=0, ocount=0, tcount=0;
		/** Total base count from file header */
		/** Combined base count for all other (non-primary) taxa */
		/** Base count for the primary taxa */
		long pbases=0, obases=0, tbases=0;
		/**
		 * Parts per million contamination level (other bases relative to total or primary+other)
		 */
		double ppm;
		/** Parsed barcode components of the primary taxa for filtering comparisons */
		/** Parsed components of the primary taxa name for filtering comparisons */
		String[] name0=null, barcode0=null;
		
	}
	
	/** List of input file paths to process */
	final ArrayList<String> in;
	/** Output file path or "stdout" for console output */
	final String out;
	/** Whether to ignore contamination from taxa with similar names */
	boolean ignoreSameTaxa=false;
	/** Whether to ignore contamination from sequences with matching barcodes */
	boolean ignoreSameBarcode=false;
	/** Whether to ignore contamination from samples from the same location */
	boolean ignoreSameLocation=false;
	/**
	 * Whether to use total bases as denominator for ppm calculation instead of primary+other
	 */
	boolean totalDenominator=false;
	/** Whether to include aggregated totals in the output */
	boolean printTotal=true;
	
}
