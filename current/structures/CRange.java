package structures;

import java.util.ArrayList;
import java.util.Collections;

import shared.Tools;

/** A numeric range, assuming 0-based, base-centered numbering,
 * including a contig number. */
public class CRange implements Comparable<CRange>{
	/*--------------------------------------------------------------*/
	
	public CRange(long c_, int a_, int b_){
		this(c_, a_, b_, null);
	}
	
	public CRange(long c_, int a_, int b_, Object o_){
		a=a_;
		b=b_;
		c=c_;
		obj=o_;
		assert(a<=b) : a+">"+b;
	}
	
	/*--------------------------------------------------------------*/
	
	public boolean includes(int p){
		return p>=a && p<=b;
	}
	
	public boolean intersects(int p1, int p2){
		return overlap(a, b, p1, p2);
	}
	
	public boolean adjacent(int p1, int p2) {
		return adjacent(a, b, p1, p2);
	}
	
	public boolean touches(int p1, int p2) {
		return touch(a, b, p1, p2);
	}
	
	public boolean includes(int p1, int p2){
		return include(a, b, p1, p2);
	}
	/*--------------------------------------------------------------*/
	
	public boolean intersects(CRange r){
		return c==r.c && intersects(r.a, r.b);
	}
	
	public boolean touches(CRange r){
		return c==r.c && (intersects(r.a, r.b) || adjacent(r.a, r.b));
	}
	
	public boolean includes(CRange r){
		return c==r.c && includes(r.a, r.b);
	}
	
	public int length() {
		return b-a+1;
	}
	
	/*--------------------------------------------------------------*/
	
	public CRange merge(CRange r){
		assert(touches(r));
		CRange r2=new CRange(c, min(a, r.a), max(b, r.b), obj);
		
		assert(r2.includes(this));
		assert(r2.includes(r));
		assert(r2.length()<=length()+r.length());
		return r2;
	}
	
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
	
	public static boolean include(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2>=a1 && b2<=b1;
	}
	
	public static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	public static boolean touch(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1+1 && b2>=a1-1;
	}
	
	public static boolean adjacent(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2==b1+1 && b2==a1-1;
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/*--------------------------------------------------------------*/
	
	/** Left point, inclusive */
	public int a;
	/** Right point, inclusive */
	public int b;
	/** Contig or sequence number */
	public final long c;
	
	public Object obj; //For attaching things
}
