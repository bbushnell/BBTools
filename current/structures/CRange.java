package structures;

import java.util.ArrayList;
import java.util.Collections;

import shared.Tools;

/**
 * Mutable inclusive range on one contig, using 0-based, base-centered coordinates.
 * The contig ID is fixed; endpoints and the optional attached object are public.
 * Equality and hashing use contig and endpoints, ignoring the attached object.
 * Callers must serialize mutations and avoid changing endpoints while a range
 * is used as a hash key or as an element of a sorted collection.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class CRange implements Comparable<CRange>{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an inclusive range without an attached object. */
	public CRange(long c_, int a_, int b_){
		this(c_, a_, b_, null);
	}

	/** Creates a range with ordered endpoints and a borrowed attachment reference. */
	public CRange(long c_, int a_, int b_, Object o_){
		a=a_;
		b=b_;
		c=c_;
		obj=o_;
		assert(a<=b) : a+">"+b;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Tests whether p lies within both inclusive endpoints. */
	public boolean includes(int p){
		return p>=a && p<=b;
	}

	/** Tests coordinate overlap with the inclusive interval p1 through p2. */
	public boolean intersects(int p1, int p2){
		return overlap(a, b, p1, p2);
	}

	/** Tests immediate coordinate adjacency using the static int-arithmetic helper. */
	public boolean adjacent(int p1, int p2){
		return adjacent(a, b, p1, p2);
	}

	/** Tests coordinate overlap or immediate adjacency, with widened endpoint arithmetic. */
	public boolean touches(int p1, int p2){
		return touch(a, b, p1, p2);
	}

	/** Tests whether this range contains both endpoints of an ordered interval. */
	public boolean includes(int p1, int p2){
		return include(a, b, p1, p2);
	}
	/*--------------------------------------------------------------*/

	/** Tests overlap on the same contig; r must be nonnull. */
	public boolean intersects(CRange r){
		return c==r.c && intersects(r.a, r.b);
	}

	/** Tests same-contig overlap or int-helper adjacency; r must be nonnull. */
	public boolean touches(CRange r){
		return c==r.c && (intersects(r.a, r.b) || adjacent(r.a, r.b));
	}

	/** Tests containment on the same contig; r must be nonnull. */
	public boolean includes(CRange r){
		return c==r.c && includes(r.a, r.b);
	}

	/** Returns the inclusive length as an int; spans above Integer.MAX_VALUE overflow. */
	public int length(){
		//TODO: Confirm supported span bounds; [0, Integer.MAX_VALUE] produces a negative length.
		return b-a+1;
	}

	/*--------------------------------------------------------------*/

	/**
	 * Returns a new covering range for touching same-contig inputs.
	 * Neither input is mutated; this range's attachment reference is retained.
	 * Assertions check touching, containment and an int-arithmetic length bound.
	 */
	public CRange merge(CRange r){
		assert(touches(r));
		CRange r2=new CRange(c, min(a, r.a), max(b, r.b), obj);

		assert(r2.includes(this));
		assert(r2.includes(r));
		assert(r2.length()<=length()+r.length());
		return r2;
	}

	/** Expands this range to cover a touching same-contig range, retaining this attachment. */
	public void absorb(CRange r){
		assert(touches(r));
		a=min(a, r.a);
		b=max(b, r.b);
	}

	/** Hashes the contig and current endpoints, excluding the attachment. */
	@Override
	public int hashCode(){
		return Integer.rotateLeft(~a, 16)^b^Long.hashCode(Long.rotateRight(c, 8));
	}

	/** Tests contig/endpoints equality; null and other types are unequal. */
	@Override
	public boolean equals(Object r){
		return r instanceof CRange && equals((CRange)r);
	}

	/** Tests contig/endpoints equality; null is unequal and attachments are ignored. */
	public boolean equals(CRange r){
		return r!=null && c==r.c && a==r.a && b==r.b;
	}

	/** Compares contig, start, then end; endpoint differences must fit in an int. */
	@Override
	public int compareTo(CRange r){
		//TODO: Endpoint subtraction can overflow when signed coordinates span the int range.
		if(c!=r.c){return c>r.c ? 1 : -1;}
		if(a!=r.a){return a-r.a;}
		return b-r.b;
	}

	/** Formats the contig and inclusive endpoints, omitting a repeated endpoint. */
	@Override
	public String toString(){
		return "(c"+c+":"+a+(a==b ? "" : (" - "+b))+")";
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Merges touching ranges in place and returns the number of removed list entries.
	 * Inputs must be nonnull and ordered by compareTo, or sort must be true to sort
	 * them here. Each merged group mutates its first range and retains that range's
	 * attachment; later group entries are removed from the list.
	 */
	public static int mergeList(ArrayList<CRange> ranges, boolean sort){
		if(ranges.size()<2){return 0;}
		if(sort){Collections.sort(ranges);}

		CRange current=ranges.get(0);
		int removed=0;
		for(int i=1; i<ranges.size(); i++){
			CRange r=ranges.get(i);
			if(current.touches(r)){
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

	/** Tests whether ordered inclusive interval 1 contains ordered interval 2. */
	public static boolean include(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2>=a1 && b2<=b1;
	}

	/** Tests overlap between two ordered inclusive intervals. */
	public static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}

	/** Tests overlap or immediate adjacency, widening endpoint arithmetic to avoid int wrap. */
	public static boolean touch(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=(long)b1+1 && b2>=(long)a1-1;
	}

	/** Tests end-plus-one adjacency with int arithmetic; extreme signed endpoints can wrap. */
	public static boolean adjacent(int a1, int b1, int a2, int b2){
		//TODO: Widen endpoint arithmetic; a maximum end plus one wraps to a negative start.
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2==b1+1 || a1==b2+1;
	}

	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Inclusive start coordinate. */
	public int a;
	/** Inclusive end coordinate, at least a for a valid range. */
	public int b;
	/** Contig identifier used for equality, ordering and range-to-range operations. */
	public final long c;

	/** Optional borrowed attachment, excluded from equality and hashing. */
	public Object obj;
}
