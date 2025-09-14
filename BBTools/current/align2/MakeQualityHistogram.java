package align2;

import java.util.ArrayList;

import fileIO.ReadWrite;
import shared.Tools;
import stream.ConcurrentLegacyReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.SiteScore;
import structures.ListNum;

public class MakeQualityHistogram {
	
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

	private static void processList(ArrayList<Read> list, int[][] mapped, int[][] paired) {
		for(Read r : list){
			processRead(r, mapped, paired);
//			if(r.mate!=null){
//				processRead(r.mate, mapped, paired);
//			}
		}
	}

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
