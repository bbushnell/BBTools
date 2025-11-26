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
	
	/**
	 * Adds a palindrome to the tracking statistics.
	 * Records palindrome length, loop size, tail lengths, matches, mismatches,
	 * and region length for histogram generation.
	 *
	 * @param p The palindrome to add to statistics
	 * @param a0 Start position of the region of interest
	 * @param b0 End position of the region of interest
	 */
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
	
	/**
	 * Merges statistics from another PalindromeTracker into this one.
	 * Combines all histogram counts and total palindrome counts.
	 * @param p The PalindromeTracker to merge into this one
	 * @return This PalindromeTracker for method chaining
	 */
	public PalindromeTracker add(PalindromeTracker p) {
		for(int i=0; i<lists.length; i++) {
			lists[i].incrementBy(p.lists[i]);
		}
		found+=p.found;
		return this;
	}
	
	/**
	 * Appends formatted histogram data to a ByteBuilder.
	 * Creates a tab-separated table with headers for all tracked statistics.
	 * @param bb The ByteBuilder to append the formatted data to
	 * @return The same ByteBuilder for method chaining
	 */
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

	/** Total number of palindromes found and tracked */
	public long found=0;
	
	/** Histogram of palindrome lengths */
	public LongList plenList=new LongList();
	/** Histogram of loop sizes in palindromes */
	public LongList loopList=new LongList();
	/**
	 * Histogram of tail lengths (distances from palindrome to region boundaries)
	 */
	public LongList tailList=new LongList();
	/** Histogram of differences between left and right tail lengths */
	public LongList tailDifList=new LongList();
	/** Histogram of number of matching bases in palindromes */
	public LongList matchList=new LongList();
	/** Histogram of number of mismatching bases in palindromes */
	public LongList mismatchList=new LongList();
	/** Histogram of region of interest lengths */
	public LongList rlenList=new LongList();//Region of interest length
	
	/** Array containing all histogram lists for bulk operations */
	public final LongList[] lists={plenList, loopList, tailList, 
			tailDifList, matchList, mismatchList, rlenList};
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Maximum value for histogram display capping */
	public static int histmax=50;
	
}
