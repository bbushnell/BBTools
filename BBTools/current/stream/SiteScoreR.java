package stream;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import shared.Shared;



/**
 * @author Brian Bushnell
 * @date Jul 16, 2012
 *
 */
public final class SiteScoreR implements Comparable<SiteScoreR>{
	
	public SiteScoreR(SiteScore ss, int readlen_, long numericID_, byte pairnum_){
		this(ss.chrom, ss.strand, ss.start, ss.stop, readlen_, numericID_, pairnum_, ss.score, ss.pairedScore, ss.perfect, ss.semiperfect);
	}
	
	public SiteScoreR(int chrom_, byte strand_, int start_, int stop_, int readlen_, long numericID_, byte pairnum_, int score_, int pscore_, boolean perfect_, boolean semiperfect_){
		chrom=chrom_;
		strand=strand_;
		start=start_;
		stop=stop_;
		readlen=readlen_;
		numericID=numericID_;
		pairnum=pairnum_;
		score=score_;
		pairedScore=pscore_;
		perfect=perfect_;
		semiperfect=semiperfect_|perfect_;
		assert(start_<=stop_) : this.toText();
	}
	
	@Override
	public int compareTo(SiteScoreR other) {
		int x=other.score-score;
		if(x!=0){return x;}
		
		x=other.pairedScore-pairedScore;
		if(x!=0){return x;}
		
		x=chrom-other.chrom;
		if(x!=0){return x;}
		
		x=strand-other.strand;
		if(x!=0){return x;}
		
		x=start-other.start;
		return x;
	}
	
	@Override
	public boolean equals(Object other){
		return compareTo((SiteScoreR)other)==0;
	}
	
	@Override
	public int hashCode() {
		assert(false) : "This class should not be hashed.";
		return super.hashCode();
	}
	
	public boolean equals(SiteScore other){
		if(other.start!=start){return false;}
		if(other.stop!=stop){return false;}
		if(other.chrom!=chrom){return false;}
		if(other.strand!=strand){return false;}
		return true;
	}
	
	public boolean equals(SiteScoreR other){
		return compareTo(other)==0;
	}
	
	@Override
	public String toString(){
//		StringBuilder sb=new StringBuilder();
//		sb.append('\t');
//		sb.append(start);
//		int spaces=10-sb.length();
//		for(int i=0; i<spaces; i++){
//			sb.append(" ");
//		}
//		sb.append('\t');
//		sb.append(quickScore);
//		sb.append('\t');
//		sb.append(score);
//
//		return "chr"+chrom+"\t"+Gene.strandCodes[strand]+sb;
		return toText().toString();
	}
	
//	9+2+1+9+9+1+1+4+4+4+4+gaps
	public StringBuilder toText(){
		StringBuilder sb=new StringBuilder(50);
		if(correct){sb.append('*');}
		sb.append(chrom);
		sb.append(',');
		sb.append(strand);
		sb.append(',');
		sb.append(start);
		sb.append(',');
		sb.append(stop);
		sb.append(',');
		sb.append(readlen);
		sb.append(',');
		sb.append(numericID);
		sb.append(',');
		sb.append(pairnum);
		sb.append(',');
		sb.append((semiperfect ? 1 : 0));
		sb.append((perfect ? 1 : 0));
		sb.append(',');
		sb.append(pairedScore);
		sb.append(',');
		sb.append(score);
//		sb.append(',');
//		sb.append((long)normalizedScore);
		return sb;
//		chrom+","+strand+","+start+","+stop+","+(rescued ? 1 : 0)+","+
//		(perfect ? 1 : 0)+","+quickScore+","+slowScore+","+pairedScore+","+score;
	}
	
	public final boolean overlaps(SiteScoreR ss){
		return chrom==ss.chrom && strand==ss.strand && overlap(start, stop, ss.start, ss.stop);
	}
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	public static String header() {
		return "chrom,strand,start,stop,readlen,numericID,pairnum,semiperfect+perfect,quickScore,slowScore,pairedScore,score";
	}
	
	public static SiteScoreR fromText(String s){
//		System.err.println("Trying to make a SS from "+s);
		String line[]=s.split(",");
		
		SiteScoreR ss;

		assert(line.length==10 || line.length==11) : "\n"+line.length+"\n"+s+"\n"+Arrays.toString(line)+"\n";
		boolean correct=false;
		if(line[0].charAt(0)=='*'){
			correct=true;
			line[0]=line[0].substring(1);
		}
		int chrom=Byte.parseByte(line[0]);
		byte strand=Byte.parseByte(line[1]);
		int start=Integer.parseInt(line[2]);
		int stop=Integer.parseInt(line[3]);
		int readlen=Integer.parseInt(line[4]);
		long numericID=Long.parseLong(line[5]);
		byte pairnum=Byte.parseByte(line[6]);
		int p=Integer.parseInt(line[7], 2);
		boolean perfect=(p&1)==1;
		boolean semiperfect=(p&2)==2;
		int pairedScore=Integer.parseInt(line[8]);
		int score=Integer.parseInt(line[9]);
		ss=new SiteScoreR(chrom, strand, start, stop, readlen, numericID, pairnum, score, pairedScore, perfect, semiperfect);
		ss.correct=correct;
		
		return ss;
	}
	
	public static SiteScoreR[] fromTextArray(String s){
		String[] split=s.split("\t");
		SiteScoreR[] out=new SiteScoreR[split.length];
		for(int i=0; i<split.length; i++){out[i]=fromText(split[i]);}
		return out;
	}
	
	public boolean positionalMatch(SiteScoreR b){
//		return chrom==b.chrom && strand==b.strand && start==b.start && stop==b.stop;
		if(chrom!=b.chrom || strand!=b.strand || start!=b.start || stop!=b.stop){
			return false;
		}
		return true;
	}
	
	public static class PositionComparator implements Comparator<SiteScoreR>{
		
		private PositionComparator(){}
		
		@Override
		public int compare(SiteScoreR a, SiteScoreR b) {
			if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
			if(a.start!=b.start){return a.start-b.start;}
			if(a.stop!=b.stop){return a.stop-b.stop;}
			if(a.strand!=b.strand){return a.strand-b.strand;}
			if(a.score!=b.score){return b.score-a.score;}
			if(a.perfect!=b.perfect){return a.perfect ? -1 : 1;}
			return 0;
		}
		
		public void sort(List<SiteScoreR> list){
			if(list==null || list.size()<2){return;}
			Collections.sort(list, this);
		}
		
		public void sort(SiteScoreR[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}
	
	public static class NormalizedComparator implements Comparator<SiteScoreR>{
		
		private NormalizedComparator(){}
		
		@Override
		public int compare(SiteScoreR a, SiteScoreR b) {
			if((int)a.normalizedScore!=(int)b.normalizedScore){return (int)b.normalizedScore-(int)a.normalizedScore;}
			if(a.score!=b.score){return b.score-a.score;}
			if(a.pairedScore!=b.pairedScore){return b.pairedScore-a.pairedScore;}
			if(a.retainVotes!=b.retainVotes){return b.retainVotes-a.retainVotes;}
			if(a.perfect!=b.perfect){return a.perfect ? -1 : 1;}
			if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
			if(a.start!=b.start){return a.start-b.start;}
			if(a.stop!=b.stop){return a.stop-b.stop;}
			if(a.strand!=b.strand){return a.strand-b.strand;}
			return 0;
		}
		
		public void sort(List<SiteScoreR> list){
			if(list==null || list.size()<2){return;}
			Collections.sort(list, this);
		}
		
		public void sort(SiteScoreR[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}
	
	public static class IDComparator implements Comparator<SiteScoreR>{
		
		private IDComparator(){}
		
		@Override
		public int compare(SiteScoreR a, SiteScoreR b) {
			if(a.numericID!=b.numericID){return a.numericID>b.numericID ? 1 : -1;}
			if(a.pairnum!=b.pairnum){return a.pairnum-b.pairnum;}
			
			if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
			if(a.start!=b.start){return a.start-b.start;}
			if(a.stop!=b.stop){return a.stop-b.stop;}
			if(a.strand!=b.strand){return a.strand-b.strand;}
			if(a.score!=b.score){return b.score-a.score;}
			if(a.perfect!=b.perfect){return a.perfect ? -1 : 1;}
			return 0;
		}
		
		public void sort(ArrayList<SiteScoreR> list){
			if(list==null || list.size()<2){return;}
			Shared.sort(list, this);
		}
		
		public void sort(SiteScoreR[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}

	public static final PositionComparator PCOMP=new PositionComparator();
	public static final NormalizedComparator NCOMP=new NormalizedComparator();
	public static final IDComparator IDCOMP=new IDComparator();
	
	public int reflen(){return stop-start+1;}
	
	public int start;
	public int stop;
	public int readlen;
	public int score;
	public int pairedScore;
	public final int chrom;
	public final byte strand;
	public boolean perfect;
	public boolean semiperfect;
	public final long numericID;
	public final byte pairnum;
	public float normalizedScore;
//	public int weight=0; //Temp variable, for calculating normalized score
	public boolean correct=false;
	public int retainVotes=0;
	
}
