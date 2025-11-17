package ml;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date June 1, 2016
 *
 */
public class ProcessBBMergeHeaders {

	/** Program entry point for processing BBMerge headers.
	 * @param args Command-line arguments specifying input/output files and parameters */
	public static void main(String[] args){
		Timer t=new Timer();
		ProcessBBMergeHeaders x=new ProcessBBMergeHeaders(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs ProcessBBMergeHeaders instance and parses command-line arguments.
	 * Sets up input/output file formats and initializes processing parameters.
	 * @param args Command-line arguments containing file paths and options
	 */
	public ProcessBBMergeHeaders(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.FORCE_INTERLEAVED=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
	}
	
	/**
	 * Main processing method that reads BBMerge output and extracts header statistics.
	 * Processes reads containing merge statistics, parses header information, and outputs
	 * tabular data for machine learning analysis. Adds column headers and processes
	 * each read to extract overlap metrics, insert size predictions, and quality scores.
	 *
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ffin1, null);
			cris.start();
		}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		long readsProcessed=0;
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			ArrayList<Read> keep=new ArrayList<Read>(reads.size());
			keep.add(new Read(null, null, headerString(), 0));

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					Header h=makeHeader(r1.id);
					if(h!=null){
						r1.id=h.toString();
						keep.add(r1);
					}
					
					readsProcessed++;
				}
				
				if(ros!=null){ros.add(keep, ln.id);}
				keep=new ArrayList<Read>();

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		ReadWrite.closeStreams(cris, ros);
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses a BBMerge header line into a structured Header object.
	 * Validates that the line contains required BBMerge statistics markers
	 * and attempts to parse the encoded merge statistics.
	 *
	 * @param line The read identifier line containing merge statistics
	 * @return Header object containing parsed statistics, or null if invalid
	 */
	private Header makeHeader(String line){
		if(!line.startsWith("insert=")){return null;}
		if(!line.contains(" mo=")){return null;}
		Header h=new Header(line);
		return h.valid ? h : null;
	}
	
	/**
	 * Returns the column header string for tabular output.
	 * Defines the columns that will be output for each processed header,
	 * including correctness, overlap metrics, error rates, and probabilities.
	 * @return Tab-separated header string for output table
	 */
	public String headerString(){
		return "#Correct\tminOverlap\tbestOverlap\tbestBadInt\tsecondBestOverlap\tsecondBestBadInt\t"
				+ "expectedErrors\tbestExpectedErrors\tbestRatio\tbestBad\tsecondBestRatio\tsecondBestBad\tprobability";
	}
	
	/**
	 * Represents parsed BBMerge header statistics for a single read pair.
	 * Contains overlap predictions, insert size estimates, error rates, and
	 * correctness validation for machine learning analysis of merge performance.
	 */
	private class Header {
	
		//mo=14_r1ee=5.2728_r2ee=3.4856_bi=202_bo=98_bb=5.3063_br=0.0598_bbi=6_sbi=270_sbo=30_sbb=12.4775_sbr=0.4343_sbbi=14_be=6.5990_pr=0.0007
		
		/**
		 * Parses a BBMerge header line to extract all merge statistics.
		 * Decodes embedded statistics including insert sizes, overlap lengths,
		 * error rates, and probability scores from the encoded header format.
		 * @param line_ The header line containing encoded merge statistics
		 */
		Header(String line_){
			line=line_;
			String[] split=line.split(" ");
			trueInsert=Integer.parseInt(split[0].split("=")[1]);
			
			split=split[2].split("_");
			for(String s : split){
				String[] split2=s.split("=");
				String a=split2[0], b=split2[1];
				if(a.equals("mo")){
					minOverlap=Integer.parseInt(b);
				}else if(a.equals("bi")){
					bestInsert=Integer.parseInt(b);
				}else if(a.equals("bo")){
					bestOverlap=Integer.parseInt(b);
				}else if(a.equals("bbi")){
					bestBadInt=Integer.parseInt(b);
				}else if(a.equals("sbi")){
					secondBestInsert=Integer.parseInt(b);
				}else if(a.equals("sbo")){
					secondBestOverlap=Integer.parseInt(b);
				}else if(a.equals("sbbi")){
					secondBestBadInt=Integer.parseInt(b);
				}else if(a.equals("r1ee")){
					expectedErrors1=Float.parseFloat(b);
				}else if(a.equals("r2ee")){
					expectedErrors2=Float.parseFloat(b);
				}else if(a.equals("be")){
					bestExpected=Float.parseFloat(b);
				}else if(a.equals("pr")){
					probability=Float.parseFloat(b);
				}else if(a.equals("br")){
					bestRatio=Float.parseFloat(b);
				}else if(a.equals("bb")){
					bestBad=Float.parseFloat(b);
				}else if(a.equals("sbr")){
					secondBestRatio=Float.parseFloat(b);
				}else if(a.equals("sbb")){
					secondBestBad=Float.parseFloat(b);
				}else{
					throw new RuntimeException(s);
				}
			}
			
			correct=(bestInsert==trueInsert);
			valid=(split.length==15 && bestInsert>0 && secondBestInsert>0);
			assert(!valid || bestOverlap>0) : bestOverlap+", "+bestInsert;
		}
		
		@Override
		public String toString(){
			return Tools.format("%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.8f",
					correct ? 1 : 0, minOverlap, bestOverlap, bestBadInt, secondBestOverlap, secondBestBadInt,
							expectedErrors1+expectedErrors2, bestExpected, bestRatio, bestBad, secondBestRatio, secondBestBad, probability);
		}
		
		/** True insert size from the read identifier */
		int trueInsert;
		
		/** Minimum overlap length required for merging */
		int minOverlap;
		/** Expected errors in read 1 based on quality scores */
		float expectedErrors1;
		/** Expected errors in read 2 based on quality scores */
		float expectedErrors2;
		
		/** Expected errors for the best overlap prediction */
		float bestExpected;
		/** Probability score for the merge prediction */
		float probability;
		
		/** Best predicted insert size from merge analysis */
		int bestInsert;
		/** Best overlap length found during merge analysis */
		int bestOverlap;
		/** Quality ratio for the best merge prediction */
		float bestRatio;
		/** Bad score (mismatch/error metric) for best prediction */
		float bestBad;
		/** Integer version of bad score for best prediction */
		int bestBadInt;
		
		/** Second best predicted insert size */
		int secondBestInsert;
		/** Second best overlap length found */
		int secondBestOverlap;
		/** Quality ratio for the second best prediction */
		float secondBestRatio;
		/** Bad score for second best prediction */
		float secondBestBad;
		/** Integer version of bad score for second best prediction */
		int secondBestBadInt;

		/** Whether the best prediction matches the true insert size */
		boolean correct;
		/** Whether the header contains valid, complete merge statistics */
		boolean valid=false;
		
		/** Original header line containing the encoded statistics */
		String line;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Input file path */
	private String in1=null;
	/** Output file path */
	private String out1=null;
	
	/** Input file format specification */
	private final FileFormat ffin1;
	/** Output file format specification */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process, -1 for unlimited */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and results */
	private java.io.PrintStream outstream=System.err;
	/** Whether to print verbose status messages during processing */
	public static boolean verbose=false;
	
}
