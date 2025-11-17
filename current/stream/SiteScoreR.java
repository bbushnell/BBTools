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
	
	/**
	 * Constructs a SiteScoreR from an existing SiteScore with additional read metadata.
	 * Copies all scoring information from the SiteScore and adds read-specific data.
	 *
	 * @param ss The SiteScore to copy alignment data from
	 * @param readlen_ Length of the aligned read in bases
	 * @param numericID_ Unique numeric identifier for the read
	 * @param pairnum_ Pair number (1 or 2 for paired reads, 0 for single)
	 */
	public SiteScoreR(SiteScore ss, int readlen_, long numericID_, byte pairnum_){
		this(ss.chrom, ss.strand, ss.start, ss.stop, readlen_, numericID_, pairnum_, ss.score, ss.pairedScore, ss.perfect, ss.semiperfect);
	}
	
	/**
	 * Constructs a SiteScoreR with all alignment and scoring parameters.
	 * Creates a complete alignment site record with position, scoring, and read metadata.
	 *
	 * @param chrom_ Chromosome or reference sequence identifier
	 * @param strand_ Strand orientation (0=forward, 1=reverse)
	 * @param start_ Start position on reference (0-based, inclusive)
	 * @param stop_ Stop position on reference (0-based, inclusive)
	 * @param readlen_ Length of the aligned read in bases
	 * @param numericID_ Unique numeric identifier for the read
	 * @param pairnum_ Pair number (1 or 2 for paired reads, 0 for single)
	 * @param score_ Alignment score for this site
	 * @param pscore_ Paired alignment score when applicable
	 * @param perfect_ True if alignment is perfect (no mismatches)
	 * @param semiperfect_ True if alignment is semiperfect (few mismatches)
	 */
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
	
	/**
	 * Tests positional equality with a SiteScore object.
	 * Compares chromosome, strand, start, and stop positions only.
	 * @param other SiteScore to compare positions with
	 * @return true if positions match exactly
	 */
	public boolean equals(SiteScore other){
		if(other.start!=start){return false;}
		if(other.stop!=stop){return false;}
		if(other.chrom!=chrom){return false;}
		if(other.strand!=strand){return false;}
		return true;
	}
	
	/**
	 * Tests equality with another SiteScoreR using compareTo semantics.
	 * @param other SiteScoreR to compare with
	 * @return true if objects are equal according to compareTo
	 */
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
	/**
	 * Creates comma-separated text representation of alignment site data.
	 * Includes position, scoring, and metadata in a parseable format.
	 * Prefixes with '*' if marked as correct for validation purposes.
	 * @return StringBuilder containing formatted alignment data
	 */
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
	
	/**
	 * Tests if this alignment site overlaps with another on the same reference.
	 * Requires same chromosome and strand, with overlapping position ranges.
	 * @param ss Other SiteScoreR to test overlap with
	 * @return true if sites overlap on same chromosome and strand
	 */
	public final boolean overlaps(SiteScoreR ss){
		return chrom==ss.chrom && strand==ss.strand && overlap(start, stop, ss.start, ss.stop);
	}
	/**
	 * Tests if two inclusive ranges overlap.
	 * Range [a1,b1] overlaps [a2,b2] if they share any positions.
	 *
	 * @param a1 Start of first range (inclusive)
	 * @param b1 End of first range (inclusive)
	 * @param a2 Start of second range (inclusive)
	 * @param b2 End of second range (inclusive)
	 * @return true if ranges overlap
	 */
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	/**
	 * Returns column headers for text output format.
	 * Matches the field order produced by toText() method.
	 * @return Comma-separated column header string
	 */
	public static String header() {
		return "chrom,strand,start,stop,readlen,numericID,pairnum,semiperfect+perfect,quickScore,slowScore,pairedScore,score";
	}
	
	/**
	 * Parses a SiteScoreR from comma-separated text format.
	 * Handles both 10 and 11 column formats, with optional correctness flag.
	 * Reconstructs all alignment data from the serialized string.
	 *
	 * @param s Comma-separated string in toText() format
	 * @return Reconstructed SiteScoreR object
	 */
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
	
	/**
	 * Parses multiple SiteScoreR objects from tab-separated text.
	 * Each tab-separated field should contain a comma-separated SiteScoreR.
	 * @param s Tab-separated string of SiteScoreR text representations
	 * @return Array of parsed SiteScoreR objects
	 */
	public static SiteScoreR[] fromTextArray(String s){
		String[] split=s.split("\t");
		SiteScoreR[] out=new SiteScoreR[split.length];
		for(int i=0; i<split.length; i++){out[i]=fromText(split[i]);}
		return out;
	}
	
	/**
	 * Tests if this alignment has identical position to another.
	 * Compares chromosome, strand, start, and stop coordinates only.
	 * @param b Other SiteScoreR to compare positions with
	 * @return true if positions match exactly
	 */
	public boolean positionalMatch(SiteScoreR b){
//		return chrom==b.chrom && strand==b.strand && start==b.start && stop==b.stop;
		if(chrom!=b.chrom || strand!=b.strand || start!=b.start || stop!=b.stop){
			return false;
		}
		return true;
	}
	
	/**
	 * Comparator for sorting SiteScoreR objects by genomic position.
	 * Orders by chromosome, start position, stop position, strand, then score.
	 * Used for position-based analysis and coordinate-sorted output.
	 */
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
		
		/**
		 * Sorts a list of SiteScoreR objects by genomic position.
		 * No-op for null or single-element lists.
		 * @param list List to sort in-place by position
		 */
		public void sort(List<SiteScoreR> list){
			if(list==null || list.size()<2){return;}
			Collections.sort(list, this);
		}
		
		/**
		 * Sorts an array of SiteScoreR objects by genomic position.
		 * No-op for null or single-element arrays.
		 * @param list Array to sort in-place by position
		 */
		public void sort(SiteScoreR[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}
	
	/**
	 * Comparator for sorting SiteScoreR objects by normalized alignment scores.
	 * Primary sort by normalized score, then regular score, paired score, and retain votes.
	 * Used for score-based ranking after normalization adjustments.
	 */
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
		
		/**
		 * Sorts a list of SiteScoreR objects by normalized scores.
		 * No-op for null or single-element lists.
		 * @param list List to sort in-place by normalized score
		 */
		public void sort(List<SiteScoreR> list){
			if(list==null || list.size()<2){return;}
			Collections.sort(list, this);
		}
		
		/**
		 * Sorts an array of SiteScoreR objects by normalized scores.
		 * No-op for null or single-element arrays.
		 * @param list Array to sort in-place by normalized score
		 */
		public void sort(SiteScoreR[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}
	
	/**
	 * Comparator for sorting SiteScoreR objects by read identifier.
	 * Primary sort by numeric ID and pair number, then by position and score.
	 * Used for grouping alignments by read for paired-end processing.
	 */
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
		
		/**
		 * Sorts an ArrayList of SiteScoreR objects by read ID.
		 * Uses Shared.sort() for ArrayList-specific optimization.
		 * @param list ArrayList to sort in-place by read ID
		 */
		public void sort(ArrayList<SiteScoreR> list){
			if(list==null || list.size()<2){return;}
			Shared.sort(list, this);
		}
		
		/**
		 * Sorts an array of SiteScoreR objects by read ID.
		 * No-op for null or single-element arrays.
		 * @param list Array to sort in-place by read ID
		 */
		public void sort(SiteScoreR[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}

	/** Singleton instance of PositionComparator for position-based sorting */
	public static final PositionComparator PCOMP=new PositionComparator();
	/** Singleton instance of NormalizedComparator for score-based sorting */
	public static final NormalizedComparator NCOMP=new NormalizedComparator();
	/** Singleton instance of IDComparator for read ID-based sorting */
	public static final IDComparator IDCOMP=new IDComparator();
	
	/** Calculates the reference length spanned by this alignment.
	 * @return Length in bases from start to stop (inclusive) */
	public int reflen(){return stop-start+1;}
	
	/** Start position of alignment on reference sequence (0-based, inclusive) */
	public int start;
	/** Stop position of alignment on reference sequence (0-based, inclusive) */
	public int stop;
	/** Length of the aligned read in bases */
	public int readlen;
	/** Alignment score for this site */
	public int score;
	/** Combined score when considering paired-end alignment */
	public int pairedScore;
	/** Chromosome or reference sequence identifier */
	public final int chrom;
	/** Strand orientation (0=forward, 1=reverse) */
	public final byte strand;
	/** True if alignment is perfect with no mismatches */
	public boolean perfect;
	/** True if alignment is semiperfect with few mismatches */
	public boolean semiperfect;
	/** Unique numeric identifier for the read */
	public final long numericID;
	/** Pair number (1 or 2 for paired reads, 0 for single reads) */
	public final byte pairnum;
	/** Score normalized for comparison across different alignment contexts */
	public float normalizedScore;
//	public int weight=0; //Temp variable, for calculating normalized score
	/** Flag indicating if this alignment is marked as correct for validation */
	public boolean correct=false;
	/** Number of votes to retain this alignment during filtering */
	public int retainVotes=0;
	
}
