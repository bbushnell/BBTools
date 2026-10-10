package tracker;

import java.util.Arrays;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;
import structures.LongList;

/**
 * Tracks the number of homopolymers observed of given lengths.
 * Per-sequence mode counts the longest run of each base, including a zero bin
 * for absent bases. Per-polymer mode counts every observed run. Runs compare
 * literal bytes: a case change breaks a run even though both cases are recognized.
 * Unknown symbols break runs and are omitted. Instances are unsynchronized;
 * keep PER_SEQUENCE fixed during collection/merging. Add/reset invalidates cached totals.
 *
 * @author Brian Bushnell
 * @date August 27, 2018
 */
public class PolymerTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a PolymerTracker and initializes internal counters. */
	public PolymerTracker(){
		reset();
	}

	/**
	 * Resets all homopolymer counters and longest-length trackers to empty state.
	 */
	public void reset(){
		cumulativeACGTN=null;
		Arrays.fill(maxACGTN, 0);
		for(int i=0; i<5; i++){
			countsACGTN[i]=new LongList();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Public Add Methods      ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds a paired read (and mate) to the homopolymer statistics if present.
	 * @param r Read (with optional mate) to analyze */
	public void addPair(Read r){
		if(r==null){return;}
		add(r.bases);
		add(r.mate);
	}

	/** Adds a single read to the homopolymer statistics.
	 * @param r Read to analyze */
	public void add(Read r){
		if(r==null){return;}
		add(r.bases);
	}

	/** Merges counts from another PolymerTracker into this one.
	 * @param pt Tracker to merge */
	public void add(PolymerTracker pt){
		cumulativeACGTN=null;
		for(int i=0; i<5; i++){
			LongList list=countsACGTN[i];
			LongList ptList=pt.countsACGTN[i];
			for(int len=0; len<ptList.size; len++){
				long count=ptList.get(len);
				list.increment(len, count);
			}
		}
	}

	/** Records a sequence under PER_SEQUENCE; ignores null and empty input. */
	public void add(byte[] bases){
		if(bases==null || bases.length<1){return;}
		cumulativeACGTN=null;
		if(PER_SEQUENCE){
			addPerSequence(bases);
		}else{
			addPerPolymer(bases);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns borrowed histograms whose bin n counts lengths at least n.
	 * Do not mutate the lists. Subsequent tracker mutation invalidates the cache;
	 * an earlier returned array remains a snapshot and does not update itself. */
	public LongList[] accumulate(){
		if(cumulativeACGTN!=null){return cumulativeACGTN;}
		LongList[] sums=new LongList[5];
		for(int i=0; i<5; i++){//Make reverse-cumulative version
			LongList list=countsACGTN[i];
			LongList sumList=new LongList(list.size);
			//n comprehension: reverse-cumulative sum[len]=sum[len+1]+count[len]. Correct because (a) descending len means
			//n sum[len+1] was set on the PRIOR iteration, and (b) the base case sum[list.size] reads past sumList's live size
			//n so LongList.get returns 0 (default) — no seed needed. sums[len] thus = total count of homopolymers of length>=len.
			for(int len=list.size-1; len>=0; len--){
				sumList.set(len, sumList.get(len+1)+list.get(len));
			}
			sums[i]=sumList;
		}
		cumulativeACGTN=sums;
		return sums;
	}

	/** Formats exact-length counts through the longest observed length, including zero bins. */
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

	/** Formats counts of lengths at least each threshold, refreshing the cache if needed. */
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

	public double calcRatio(byte base1, byte base2, int length){
		long count1=getCount(base1, length);
		long count2=getCount(base2, length);
		return count1/Tools.max(1.0, count2);
	}

	public long getCount(byte base, int length){
		//The mapping accepts both cases of A/C/G/T/N and U as T; other query keys are invalid.
		int x=AminoAcid.baseToNumberACGTN[base];
		return countsACGTN[x].get(length);
	}

	public double calcRatioCumulative(byte base1, byte base2, int length){
		long count1=getCountCumulative(base1, length);
		long count2=getCountCumulative(base2, length);
		return count1/Tools.max(1.0, count2);
	}

	public long getCountCumulative(byte base, int length){
		int x=AminoAcid.baseToNumberACGTN[base];
		return accumulate()[x].get(length);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

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

	private void recordMax(byte base, int len){
		if(base<0){return;}
		int x=AminoAcid.baseToNumberACGTN[base];
		if(x<0){return;}
		maxACGTN[x]=Tools.max(maxACGTN[x], len);
	}

	private void recordCounts(byte base, int len){
		if(base<0){return;}
		int x=AminoAcid.baseToNumberACGTN[base];
		if(x<0){return;}
		countsACGTN[x].increment(len, 1);
	}

	/** Number of bins, one more than the greatest stored length. */
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

	private final int[] maxACGTN=KillSwitch.allocInt1D(5);
	final LongList[] countsACGTN=new LongList[5];
	/** Borrowed query cache, invalidated by every count mutation (JT001). */
	private LongList[] cumulativeACGTN;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean PER_SEQUENCE=true;
	public static boolean CUMULATIVE=true;

}
