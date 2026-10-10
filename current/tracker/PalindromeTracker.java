package tracker;

import repeat.Palindrome;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Tracks palindrome statistics to determine which types occur in a given feature.
 * Collects counts and histograms for palindrome lengths, loop sizes, tail lengths,
 * matches, mismatches, and region lengths for analysis and reporting.
 * Coordinates are inclusive. Instances are unsynchronized. Reporting caps the
 * original histograms in place; collect and merge before final reporting.
 *
 * @author Brian Bushnell
 * @date Sept 3, 2023
 */
public class PalindromeTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Records a palindrome contained in the inclusive region [a0,b0].
	 * Tail lengths exclude the palindrome, matching Palindrome.appendTo.
	 * @param p Palindrome with nonnegative arm, loop and mismatch counts
	 * @param a0 Region start, no greater than p.a
	 * @param b0 Region end, no less than p.b */
	public void add(final Palindrome p, final int a0, final int b0){
		int tail1=p.a-a0, tail2=b0-p.b;
		if(tail1>tail2){
			int x=tail1;
			tail1=tail2;
			tail2=x;
		}
		int tailDif=tail2-tail1;
		int rlen=b0-a0+1;
		//[tracker/PalindromeTracker#001] Histogram indices require a0<=p.a<=p.b<=b0.
		//JT002: right tail uses p.b, as in Palindrome.toString(a0,b0), not p.a.
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

	/** Adds p's histograms and count, leaving p unchanged unless p is this. */
	public PalindromeTracker add(PalindromeTracker p){
		for(int i=0; i<lists.length; i++){
			lists[i].incrementBy(p.lists[i]);
		}
		found+=p.found;
		return this;
	}

	public ByteBuilder appendTo(ByteBuilder bb){
		return append(bb, "#Value\tplen\tloop\ttail\ttaildif\tmatch\tmismtch\trlen", lists, histmax);
	}

	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Generic method to create histogram tables from LongList arrays.
	 * Caps the original histograms in place and formats a tab-separated table.
	 * The final bin sums all counts at or above histmax; those original indices
	 * are lost. Avoid resuming collection into lists after reporting them.
	 *
	 * @param bb The ByteBuilder to append formatted data to
	 * @param header The header row for the table
	 * @param lists Array of LongList histograms to format
	 * @param histmax Maximum value for histogram capping
	 * @return The same ByteBuilder for method chaining
	 */
	public static ByteBuilder append(ByteBuilder bb, String header, LongList[] lists, int histmax){
		int maxSize=1;
		for(LongList ll : lists){
			ll.capHist(histmax);
			maxSize=Tools.max(maxSize, ll.size);
		}

		bb.append(header).nl();

		for(int i=0; i<maxSize; i++){
			bb.append(i);
			for(LongList ll : lists){
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
