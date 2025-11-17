package tracker;

import java.util.Arrays;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;
import structures.LongList;

/**
 * Tracks the number of homopolymers observed of given lengths.
 * Only the longest homopolymer for a given base is counted per sequence.
 * @author Brian Bushnell
 * @date August 27, 2018
 *
 */
public class PolymerTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a new PolymerTracker and initializes data structures */
	public PolymerTracker(){
		reset();
	}
	
	/** Resets all counters and reinitializes data structures to empty state */
	public void reset(){
		Arrays.fill(maxACGTN, 0);
		for(int i=0; i<5; i++){
			countsACGTN[i]=new LongList();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Public Add Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds a paired read and its mate to the homopolymer tracking.
	 * Processes both reads in the pair if they exist.
	 * @param r The read pair to analyze
	 */
	public void addPair(Read r){
		if(r==null){return;}
		add(r.bases);
		add(r.mate);
	}
	
	/** Adds a single read to the homopolymer tracking.
	 * @param r The read to analyze for homopolymers */
	public void add(Read r){
		if(r==null){return;}
		add(r.bases);
	}
	
	/**
	 * Merges counts from another PolymerTracker into this one.
	 * Combines homopolymer counts for all bases (A, C, G, T, N).
	 * @param pt The PolymerTracker to merge into this one
	 */
	public void add(PolymerTracker pt){
		for(int i=0; i<5; i++){
			LongList list=countsACGTN[i];
			LongList ptList=pt.countsACGTN[i];
			for(int len=0; len<ptList.size; len++){
				long count=ptList.get(len);
				list.increment(len, count);
			}
		}
	}
	
	/**
	 * Analyzes a sequence for homopolymers and updates counts.
	 * Uses either per-sequence or per-polymer counting based on PER_SEQUENCE flag.
	 * @param bases The DNA sequence bases to analyze
	 */
	public void add(byte[] bases){
		if(bases==null || bases.length<1){return;}
		if(PER_SEQUENCE){
			addPerSequence(bases);
		}else{
			addPerPolymer(bases);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Computes reverse-cumulative sums for all base counts.
	 * Returns cached result if already computed.
	 * @return Array of LongLists containing cumulative counts for A,C,G,T,N
	 */
	public LongList[] accumulate(){
		if(cumulativeACGTN!=null){return cumulativeACGTN;}
		LongList[] sums=new LongList[5];
		for(int i=0; i<5; i++){//Make reverse-cumulative version
			LongList list=countsACGTN[i];
			LongList sumList=new LongList(list.size);
			for(int len=list.size-1; len>=0; len--){
				sumList.set(len, sumList.get(len+1)+list.get(len));
			}
			sums[i]=sumList;
		}
		cumulativeACGTN=sums;
		return sums;
	}
	
	//Non-cumulative
	/**
	 * Generates a tab-delimited histogram of homopolymer counts.
	 * Shows raw counts (non-cumulative) for each length and base type.
	 * @return String containing histogram with header and data rows
	 */
	public String toHistogram(){
		StringBuilder sb=new StringBuilder();
		sb.append("#Length\tA\tC\tG\tT\tN\n");
		
		final int maxIndex=longest();
		for(int len=0; len<maxIndex; len++){
			sb.append(len);
			for(int i=0; i<5; i++){
				long count=countsACGTN[i].get(len);
				sb.append('\t').append(count);
			}
			sb.append('\n');
		}
		return sb.toString();
	}
	
	//Cumulative
	/**
	 * Generates a tab-delimited histogram of cumulative homopolymer counts.
	 * Shows cumulative counts (length X or longer) for each base type.
	 * @return String containing cumulative histogram with header and data rows
	 */
	public String toHistogramCumulative(){
		LongList[] sums=accumulate();
		
		StringBuilder sb=new StringBuilder();
		sb.append("#Length\tA\tC\tG\tT\tN\n");

		final int maxIndex=longest();
		for(int len=0; len<maxIndex; len++){
			sb.append(len);
			for(int i=0; i<5; i++){
				long count=sums[i].get(len);
				sb.append('\t').append(count);
			}
			sb.append('\n');
		}
		return sb.toString();
	}
	
	/**
	 * Calculates the ratio of homopolymer counts between two bases at a given length.
	 *
	 * @param base1 First base for ratio calculation
	 * @param base2 Second base for ratio calculation (denominator)
	 * @param length Homopolymer length to compare
	 * @return Ratio of base1 count to base2 count
	 */
	public double calcRatio(byte base1, byte base2, int length){
		long count1=getCount(base1, length);
		long count2=getCount(base2, length);
		return count1/Tools.max(1.0, count2);
	}
	
	/**
	 * Gets the raw count of homopolymers for a specific base and length.
	 * @param base The base to query (A, C, G, T, or N)
	 * @param length The homopolymer length to query
	 * @return Number of homopolymers of the specified base and length
	 */
	public long getCount(byte base, int length) {
		int x=AminoAcid.baseToNumberACGTN[base];
		return countsACGTN[x].get(length);
	}
	
	/**
	 * Calculates the ratio of cumulative homopolymer counts between two bases.
	 * Uses cumulative counts (length X or longer).
	 *
	 * @param base1 First base for ratio calculation
	 * @param base2 Second base for ratio calculation (denominator)
	 * @param length Minimum homopolymer length for cumulative count
	 * @return Ratio of base1 cumulative count to base2 cumulative count
	 */
	public double calcRatioCumulative(byte base1, byte base2, int length){
		long count1=getCountCumulative(base1, length);
		long count2=getCountCumulative(base2, length);
		return count1/Tools.max(1.0, count2);
	}
	
	/**
	 * Gets the cumulative count of homopolymers for a base at or above a length.
	 * @param base The base to query (A, C, G, T, or N)
	 * @param length Minimum homopolymer length for cumulative count
	 * @return Number of homopolymers of the specified base at length or longer
	 */
	public long getCountCumulative(byte base, int length) {
		int x=AminoAcid.baseToNumberACGTN[base];
		return accumulate()[x].get(length);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Increment once per sequence, 
	 * for the longest homopolymer of each base */
	private void addPerSequence(byte[] bases){
		Arrays.fill(maxACGTN, 0);
		byte prev=-1;
		int current=0;
		for(byte b : bases){
			if(b==prev){
				current++;
			}else{
				recordMax(prev, current);
				prev=b;
				current=1;
			}
		}
		recordMax(prev, current);
		
		for(int i=0; i<maxACGTN.length; i++){
			countsACGTN[i].increment(maxACGTN[i], 1);
		}
	}
	
	/** Increment once per homopolymer */
	private void addPerPolymer(byte[] bases){
		byte prev=-1;
		int current=0;
		for(byte b : bases){
			if(b==prev){
				current++;
			}else{
				recordCounts(prev, current);
				prev=b;
				current=1;
			}
		}
		recordCounts(prev, current);
	}
	
	/**
	 * Records the maximum homopolymer length seen for a base in current sequence.
	 * Updates the running maximum if the new length is longer.
	 * @param base The base type
	 * @param len The homopolymer length to potentially record
	 */
	private void recordMax(byte base, int len){
		if(base<0){return;}
		int x=AminoAcid.baseToNumberACGTN[base];
		if(x<0){return;}
		maxACGTN[x]=Tools.max(maxACGTN[x], len);
	}
	
	/**
	 * Increments the count for a specific base and homopolymer length.
	 * @param base The base type
	 * @param len The homopolymer length to increment
	 */
	private void recordCounts(byte base, int len){
		if(base<0){return;}
		int x=AminoAcid.baseToNumberACGTN[base];
		if(x<0){return;}
		countsACGTN[x].increment(len, 1);
	}
	
	/** Finds the longest homopolymer length recorded across all base types.
	 * @return Maximum homopolymer length in any count list */
	private int longest(){
		int max=0;
		for(LongList list : countsACGTN){
			max=Tools.max(list.size(), max);
		}
		return max;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Array tracking maximum homopolymer length per base (A,C,G,T,N) in current sequence
	 */
	private final int[] maxACGTN=KillSwitch.allocInt1D(5);
	/**
	 * Array of count lists for each base type (A,C,G,T,N) indexed by homopolymer length
	 */
	final LongList[] countsACGTN=new LongList[5];
	/** Cached cumulative count arrays, computed on demand by accumulate() */
	private LongList[] cumulativeACGTN;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * If true, count only longest homopolymer per base per sequence; if false, count all
	 */
	public static boolean PER_SEQUENCE=true;
	/** Controls whether cumulative statistics are computed and used in output */
	public static boolean CUMULATIVE=true;
	
}
