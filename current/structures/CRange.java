package structures;

import java.util.ArrayList;
import java.util.Collections;

import shared.Tools;

/** A numeric range, assuming 0-based, base-centered numbering,
 * including a contig number. */
public class CRange implements Comparable<CRange>{
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a CRange with specified contig and coordinates.
	 * @param c_ Contig number
	 * @param a_ Left coordinate (inclusive)
	 * @param b_ Right coordinate (inclusive)
	 */
	public CRange(long c_, int a_, int b_){
		this(c_, a_, b_, null);
	}
	
	/**
	 * Constructs a CRange with specified contig, coordinates, and attached object.
	 *
	 * @param c_ Contig number
	 * @param a_ Left coordinate (inclusive)
	 * @param b_ Right coordinate (inclusive)
	 * @param o_ Optional object to attach to this range
	 */
	public CRange(long c_, int a_, int b_, Object o_){
		a=a_;
		b=b_;
		c=c_;
		obj=o_;
		assert(a<=b) : a+">"+b;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Tests if a point is within this range.
	 * @param p Point to test
	 * @return true if point is within range [a, b]
	 */
	public boolean includes(int p){
		return p>=a && p<=b;
	}
	
	/**
	 * Tests if this range intersects with the specified coordinate range.
	 * @param p1 Start coordinate of test range
	 * @param p2 End coordinate of test range
	 * @return true if ranges overlap
	 */
	public boolean intersects(int p1, int p2){
		return overlap(a, b, p1, p2);
	}
	
	/**
	 * Tests if this range is adjacent to the specified coordinate range.
	 * Adjacent means ranges touch but don't overlap.
	 *
	 * @param p1 Start coordinate of test range
	 * @param p2 End coordinate of test range
	 * @return true if ranges are adjacent
	 */
	public boolean adjacent(int p1, int p2) {
		return adjacent(a, b, p1, p2);
	}
	
	/**
	 * Tests if this range touches the specified coordinate range.
	 * Touching includes both intersection and adjacency.
	 *
	 * @param p1 Start coordinate of test range
	 * @param p2 End coordinate of test range
	 * @return true if ranges touch or overlap
	 */
	public boolean touches(int p1, int p2) {
		return touch(a, b, p1, p2);
	}
	
	/**
	 * Tests if this range completely contains the specified coordinate range.
	 * @param p1 Start coordinate of test range
	 * @param p2 End coordinate of test range
	 * @return true if this range completely contains the test range
	 */
	public boolean includes(int p1, int p2){
		return include(a, b, p1, p2);
	}
	/*--------------------------------------------------------------*/
	
	/**
	 * Tests if this range intersects with another CRange.
	 * Only ranges on the same contig can intersect.
	 * @param r Range to test for intersection
	 * @return true if ranges are on same contig and overlap
	 */
	public boolean intersects(CRange r){
		return c==r.c && intersects(r.a, r.b);
	}
	
	/**
	 * Tests if this range touches another CRange.
	 * Only ranges on the same contig can touch.
	 * @param r Range to test for touching
	 * @return true if ranges are on same contig and touch or overlap
	 */
	public boolean touches(CRange r){
		return c==r.c && (intersects(r.a, r.b) || adjacent(r.a, r.b));
	}
	
	/**
	 * Tests if this range completely contains another CRange.
	 * Only ranges on the same contig can have inclusion relationships.
	 * @param r Range to test for inclusion
	 * @return true if this range is on same contig and completely contains r
	 */
	public boolean includes(CRange r){
		return c==r.c && includes(r.a, r.b);
	}
	
	/** Returns the length of this range in bases */
	public int length() {
		return b-a+1;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a new CRange that encompasses both this range and another.
	 * The ranges must touch (intersect or be adjacent).
	 * @param r Range to merge with
	 * @return New CRange spanning both input ranges
	 */
	public CRange merge(CRange r){
		assert(touches(r));
		CRange r2=new CRange(c, min(a, r.a), max(b, r.b), obj);
		
		assert(r2.includes(this));
		assert(r2.includes(r));
		assert(r2.length()<=length()+r.length());
		return r2;
	}
	
	/**
	 * Expands this range to encompass another range.
	 * The ranges must touch (intersect or be adjacent).
	 * Modifies this range in place.
	 * @param r Range to absorb into this range
	 */
	public void absorb(CRange r){
		assert(touches(r));
		a=min(a, r.a);
		b=max(b, r.b);
	}
	
	@Override
	public int hashCode(){
		return Integer.rotateLeft(~a, 16)^b^Long.hashCode(Long.rotateRight(c, 8));
	}
	
	@Override
	public boolean equals(Object r){
		return equals((CRange)r);
	}
	
	/**
	 * Tests equality with another CRange.
	 * Ranges are equal if they have same contig and coordinates.
	 * @param r Range to compare
	 * @return true if ranges have identical contig and coordinates
	 */
	public boolean equals(CRange r){
		return c==r.c &&a==r.a && b==r.b;
	}
	
	@Override
	public int compareTo(CRange r) {
		if(c!=r.c) {return c>r.c ? 1 : -1;}
		if(a!=r.a) {return a-r.a;}
		return b-r.b;
	}
	
	@Override
	public String toString(){
		return "(c"+c+":"+a+(a==b ? "" : (" - "+b))+")";
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Merges overlapping or adjacent ranges in a list.
	 * Modifies the input list by absorbing touching ranges and removing nulls.
	 *
	 * @param ranges List of ranges to merge
	 * @param sort Whether to sort the list before merging
	 * @return Number of ranges that were merged and removed
	 */
	public static int mergeList(ArrayList<CRange> ranges, boolean sort) {
		if(ranges.size()<2){return 0;}
		if(sort){Collections.sort(ranges);}
		
		CRange current=ranges.get(0);
		int removed=0;
		for(int i=1; i<ranges.size(); i++) {
			CRange r=ranges.get(i);
			if(current.touches(r)) {
				current.absorb(r);
				ranges.set(i, null);
				removed++;
			}else{
				current=r;
			}
		}
		if(removed>0){
			Tools.condenseStrict(ranges);
		}
		return removed;
	}
	
	/**
	 * Tests if one coordinate range completely contains another.
	 *
	 * @param a1 Start of first range
	 * @param b1 End of first range
	 * @param a2 Start of second range
	 * @param b2 End of second range
	 * @return true if first range completely contains second range
	 */
	public static boolean include(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2>=a1 && b2<=b1;
	}
	
	/**
	 * Tests if two coordinate ranges overlap.
	 *
	 * @param a1 Start of first range
	 * @param b1 End of first range
	 * @param a2 Start of second range
	 * @param b2 End of second range
	 * @return true if ranges overlap
	 */
	public static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	/**
	 * Tests if two coordinate ranges touch (overlap or are adjacent).
	 *
	 * @param a1 Start of first range
	 * @param b1 End of first range
	 * @param a2 Start of second range
	 * @param b2 End of second range
	 * @return true if ranges touch or overlap
	 */
	public static boolean touch(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1+1 && b2>=a1-1;
	}
	
	/**
	 * Tests if two coordinate ranges are adjacent but do not overlap.
	 *
	 * @param a1 Start of first range
	 * @param b1 End of first range
	 * @param a2 Start of second range
	 * @param b2 End of second range
	 * @return true if ranges are adjacent but don't overlap
	 */
	public static boolean adjacent(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2==b1+1 && b2==a1-1;
	}
	
	/** Returns the smaller of two integers */
	private static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger of two integers */
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/*--------------------------------------------------------------*/
	
	/** Left point, inclusive */
	public int a;
	/** Right point, inclusive */
	public int b;
	/** Contig or sequence number */
	public final long c;
	
	/** Optional object attachment for storing associated data */
	public Object obj; //For attaching things
}
