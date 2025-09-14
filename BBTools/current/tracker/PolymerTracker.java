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
	
	public PolymerTracker(){
		reset();
	}
	
	public void reset(){
		Arrays.fill(maxACGTN, 0);
		for(int i=0; i<5; i++){
			countsACGTN[i]=new LongList();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Public Add Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	public void addPair(Read r){
		if(r==null){return;}
		add(r.bases);
		add(r.mate);
	}
	
	public void add(Read r){
		if(r==null){return;}
		add(r.bases);
	}
	
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
	
	public long getCount(byte base, int length) {
		int x=AminoAcid.baseToNumberACGTN[base];
		return countsACGTN[x].get(length);
	}
	
	public double calcRatioCumulative(byte base1, byte base2, int length){
		long count1=getCountCumulative(base1, length);
		long count2=getCountCumulative(base2, length);
		return count1/Tools.max(1.0, count2);
	}
	
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
	private LongList[] cumulativeACGTN;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean PER_SEQUENCE=true;
	public static boolean CUMULATIVE=true;
	
}
