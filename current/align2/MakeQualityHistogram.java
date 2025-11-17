package align2;

import java.util.ArrayList;

import fileIO.ReadWrite;
import shared.Tools;
import stream.ConcurrentLegacyReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.SiteScore;
import structures.ListNum;

/**
 * Generates quality score histograms for mapped vs unmapped and paired vs single reads.
 * Analyzes quality score distributions to assess mapping and pairing success rates across
 * different quality ranges, providing diagnostic information for sequencing data evaluation.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class MakeQualityHistogram {
	
	/**
	 * Program entry point for quality histogram generation.
	 * Expects one or two file arguments for single-end or paired-end reads.
	 * @param args Command-line arguments: input file(s)
	 */
	public static void main(String[] args){
		
		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
		
		long maxReads=0;
		RTextInputStream rtis=new RTextInputStream(fname1, fname2, maxReads);
		ConcurrentLegacyReadInputStream cris=new ConcurrentLegacyReadInputStream(rtis, maxReads);
		
		int[][][] counts=process(cris);
		printMappedHistogram(counts[0]);
		System.out.println();
		printPairedHistogram(counts[1]);
//		System.out.println("*** main() finished ***");
	}
	
	/**
	 * Prints histogram showing mapping success rates by quality score.
	 * Displays quality score bins with counts of mapped vs unmapped reads
	 * and calculates mapping percentage for each quality level.
	 * @param mapped Array containing [mapped_counts, unmapped_counts] by quality
	 */
	public static void printMappedHistogram(int[][] mapped){
		System.out.println("#Error Quality Histogram");
		System.out.println("Quality\tMapped\tUnmapped\tPercent Mapped");
		for(int i=0; i<mapped[0].length; i++){
			int e=mapped[0][i];
			int m=mapped[1][i];
			float percent=e*100f/(e+m);
			System.out.println(i+"\t"+e+"\t"+m+"\t"+Tools.format("%.3f", percent));
		}
	}
	
	/**
	 * Prints histogram showing pairing success rates by quality score.
	 * Displays quality score bins with counts of paired vs single reads
	 * and calculates pairing percentage for each quality level.
	 * @param paired Array containing [paired_counts, single_counts] by quality
	 */
	public static void printPairedHistogram(int[][] paired){
		System.out.println("#Error Quality Histogram");
		System.out.println("Quality\tPaired\tSingle\tPercent Paired");
		for(int i=0; i<paired[0].length; i++){
			int e=paired[0][i];
			int m=paired[1][i];
			float percent=e*100f/(e+m);
			System.out.println(i+"\t"+e+"\t"+m+"\t"+Tools.format("%.3f", percent));
		}
	}
	
	/**
	 * Main processing method that generates quality histograms from read stream.
	 * Iterates through all reads, categorizing by mapping status and pairing status,
	 * then bins by average quality score to create distribution histograms.
	 *
	 * @param cris Concurrent read input stream
	 * @return 3D array containing [mapped_histogram, paired_histogram]
	 */
	public static int[][][] process(ConcurrentLegacyReadInputStream cris){
		
		cris.start();

		int[][] mapped=new int[2][50];
		int[][] paired=new int[2][50];
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> readlist=ln.list;
		while(!readlist.isEmpty()){
			
			processList(readlist, mapped, paired);
			
			cris.returnList(ln.id, readlist.isEmpty());
			//System.err.println("Waiting on a list...");
			ln=cris.nextList();
			readlist=ln.list;
		}
		
		//System.err.println("Returning a list... (final)");
		assert(readlist.isEmpty());
		cris.returnList(ln.id, readlist.isEmpty());
		ReadWrite.closeStream(cris);
		
		return new int[][][] {mapped, paired};
	}

	/**
	 * Processes a batch of reads and updates histogram counters.
	 * Iterates through the read list and categorizes each read by mapping
	 * and pairing status for histogram generation.
	 *
	 * @param list List of reads to process
	 * @param mapped Histogram counters for mapped vs unmapped reads
	 * @param paired Histogram counters for paired vs single reads
	 */
	private static void processList(ArrayList<Read> list, int[][] mapped, int[][] paired) {
		for(Read r : list){
			processRead(r, mapped, paired);
//			if(r.mate!=null){
//				processRead(r.mate, mapped, paired);
//			}
		}
	}

	/**
	 * Categorizes a single read by mapping and pairing status for histograms.
	 * Updates mapping position from top site if needed, calculates average quality,
	 * then increments appropriate histogram bins based on mapping and pairing status.
	 *
	 * @param r Read to analyze and categorize
	 * @param mapped Histogram counters for mapped vs unmapped reads
	 * @param paired Histogram counters for paired vs single reads
	 */
	private static void processRead(Read r, int[][] mapped, int[][] paired) {
		
		if(r.chrom<1 && r.numSites()>0){
			SiteScore ss=r.topSite(); //Should not be necessary
			r.start=ss.start;
			r.stop=ss.stop;
			r.chrom=ss.chrom;
			r.setStrand(ss.strand);
		}
		
		int avgQ=r.avgQualityInt(true, 0);
		if(r.chrom>0){
			mapped[0][avgQ]++;
		}else{
			mapped[1][avgQ]++;
		}
		if(r.paired()){
			paired[0][avgQ]++;
		}else{
			paired[1][avgQ]++;
		}
		
	}
	
}
