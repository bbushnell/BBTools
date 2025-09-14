package tracker;

import repeat.Palindrome;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Tracks palindrome stats to determine which kind occur in a given feature.
 * 
 * @author Brian Bushnell
 * @date Sept 3, 2023
 *
 */
public class PalindromeTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void add(final Palindrome p, final int a0, final int b0) {
		int tail1=p.a-a0, tail2=b0-p.a;
		if(tail1>tail2) {
			int x=tail1;
			tail1=tail2;
			tail2=x;
		}
		int tailDif=tail2-tail1;
		int rlen=b0-a0+1;
		plenList.increment(p.plen());
		loopList.increment(p.loop());
		tailList.increment(tail1);
		tailList.increment(tail2);
		tailDifList.increment(tailDif);
		matchList.increment(p.matches);
		mismatchList.increment(p.mismatches);
		rlenList.increment(rlen);
		found++;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public PalindromeTracker add(PalindromeTracker p) {
		for(int i=0; i<lists.length; i++) {
			lists[i].incrementBy(p.lists[i]);
		}
		found+=p.found;
		return this;
	}
	
	public ByteBuilder appendTo(ByteBuilder bb) {
		return append(bb, "#Value\tplen\tloop\ttail\ttaildif\tmatch\tmismtch\trlen", lists, histmax);
	}
	
	@Override
	public String toString() {
		return appendTo(new ByteBuilder()).toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Can be used to make generic histograms */
	public static ByteBuilder append(ByteBuilder bb, String header, LongList[] lists, int histmax) {
		int maxSize=1;
		for(LongList ll : lists) {
			ll.capHist(histmax);
			maxSize=Tools.max(maxSize, ll.size);
		}
		
		bb.append(header).nl();
		
		for(int i=0; i<maxSize; i++) {
			bb.append(i);
			for(LongList ll : lists) {
				bb.tab().append(ll.get(i));
			}
			bb.nl();
		}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public long found=0;
	
	public LongList plenList=new LongList();
	public LongList loopList=new LongList();
	public LongList tailList=new LongList();
	public LongList tailDifList=new LongList();
	public LongList matchList=new LongList();
	public LongList mismatchList=new LongList();
	public LongList rlenList=new LongList();//Region of interest length
	
	public final LongList[] lists={plenList, loopList, tailList, 
			tailDifList, matchList, mismatchList, rlenList};
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int histmax=50;
	
}
