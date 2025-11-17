package structures;
import java.util.ArrayList;
import java.util.Collections;

import shared.Tools;


/** A numeric range, assuming 0-based, base-centered numbering. */
public class Range implements Comparable<Range>, Cloneable{
	
	/**
	 * Constructs a range from start to end coordinates.
	 * @param aa Start coordinate (inclusive)
	 * @param bb End coordinate (inclusive)
	 */
	public Range(int aa, int bb){
		assert(aa<=bb) : aa+">"+bb;
		a=aa;
		b=bb;
	}
	
	/**
	 * Constructs a range with an associated object.
	 * @param aa Start coordinate (inclusive)
	 * @param bb End coordinate (inclusive)
	 * @param oo Object to associate with this range
	 */
	public Range(int aa, int bb, Object oo){
		assert(aa<=bb) : aa+">"+bb;
		a=aa;
		b=bb;
		obj1=oo;
	}
	
	/** Creates a shallow copy of this range.
	 * @return Cloned Range instance */
	public Range clone() {
		try {
			return (Range) super.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Tests if this range includes a specific position.
	 * @param p Position to test
	 * @return true if position is within range bounds (inclusive)
	 */
	public boolean includes(int p){
		return p>=a && p<=b;
	}
	
	/**
	 * Tests if this range intersects with another range defined by coordinates.
	 * @param p1 Start coordinate of test range
	 * @param p2 End coordinate of test range
	 * @return true if ranges overlap
	 */
	public boolean intersects(int p1, int p2){
		return overlaps(a, b, p1, p2);
	}
	
	/**
	 * Tests if this range overlaps with a Feature.
	 * @param f Feature to test against
	 * @return true if this range overlaps the feature's coordinates
	 */
	public <K extends Feature> boolean overlaps(K f) {
		return intersects(f.start(), f.stop());
	}
	
	/**
	 * Tests if this range overlaps with another range.
	 * @param f Range to test against
	 * @return true if ranges overlap
	 */
	public boolean overlaps(Range f) {
		return intersects(f.a, f.b);
	}
	
	/**
	 * Calculates the overlap length between this range and another.
	 * @param f Range to calculate overlap with
	 * @return Number of overlapping bases
	 */
	public int overlap(Range f) {
		return overlap(a, b, f.a, f.b);
	}
	
	/**
	 * Tests if this range is adjacent to another range.
	 * @param p1 Start coordinate of test range
	 * @param p2 End coordinate of test range
	 * @return true if ranges touch but don't overlap
	 */
	public boolean adjacent(int p1, int p2) {
		return adjacent(a, b, p1, p2);
	}
	
	/**
	 * Tests if this range is adjacent to a specific position.
	 * @param p1 Position to test
	 * @return true if position is adjacent to range boundary
	 */
	public boolean adjacent(int p1) {
		return adjacent(a, b, p1);
	}
	
	/**
	 * Tests if this range touches another range (overlaps or is adjacent).
	 * @param p1 Start coordinate of test range
	 * @param p2 End coordinate of test range
	 * @return true if ranges touch or overlap
	 */
	public boolean touches(int p1, int p2) {
		return touch(a, b, p1, p2);
	}
	
	/**
	 * Tests if this range touches a specific position.
	 * @param p1 Position to test
	 * @return true if position touches or overlaps the range
	 */
	public boolean touches(int p1) {
		return touch(a, b, p1);
	}
	
	/**
	 * Tests if this range completely includes another range.
	 * @param p1 Start coordinate of test range
	 * @param p2 End coordinate of test range
	 * @return true if test range is completely contained within this range
	 */
	public boolean includes(int p1, int p2){
		return include(a, b, p1, p2);
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Tests if this range intersects with another range.
	 * @param r Range to test against
	 * @return true if ranges overlap
	 */
	public boolean intersects(Range r){
		return intersects(r.a, r.b);
	}
	
	/**
	 * Tests if this range touches another range (overlaps or is adjacent).
	 * @param r Range to test against
	 * @return true if ranges touch or overlap
	 */
	public boolean touches(Range r){
		return (intersects(r.a, r.b) || adjacent(r.a, r.b));
	}
	
	/**
	 * Tests if this range completely includes another range.
	 * @param r Range to test
	 * @return true if test range is completely contained within this range
	 */
	public boolean includes(Range r){
		return includes(r.a, r.b);
	}
	
	/** Returns the length of this range in bases (inclusive of both endpoints) */
	public int length() {
		return b-a+1;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a new range that encompasses both this range and another.
	 * Ranges must touch (overlap or be adjacent) for merging to be valid.
	 * @param r Range to merge with
	 * @return New range spanning both input ranges
	 */
	public Range merge(Range r){
		assert(touches(r));
		Range r2=new Range(min(a, r.a), max(b, r.b), obj1);
		
		assert(r2.includes(this));
		assert(r2.includes(r));
		assert(r2.length()<=length()+r.length());
		return r2;
	}
	
	/**
	 * Extends this range to encompass another range.
	 * Modifies this range in-place rather than creating a new one.
	 * @param r Range to absorb
	 */
	public void absorb(Range r){
		assert(touches(r));
		a=min(a, r.a);
		b=max(b, r.b);
	}
	
	/** Extends this range to include a specific position.
	 * @param p Position to absorb into the range */
	public void absorb(int p){
		assert(touches(p));
		a=min(a, p);
		b=max(b, p);
	}
	
	/**
	 * Constrains range boundaries to valid sequence coordinates.
	 * Ensures start is non-negative and end doesn't exceed sequence length.
	 * @param length Maximum sequence length (exclusive upper bound)
	 */
	public void fixBounds(int length) {
		a=Tools.max(0, a);
		b=Tools.min(length-1, b);
	}
	
	@Override
	public int hashCode(){
		return Integer.rotateLeft(~a, 16)^b;
	}
	
	@Override
	public boolean equals(Object r){
		return equals((Range)r);
	}
	
	/**
	 * Tests equality with another range based on coordinates.
	 * @param r Range to compare with
	 * @return true if ranges have identical start and end coordinates
	 */
	public boolean equals(Range r){
		return a==r.a && b==r.b;
	}
	
	@Override
	public int compareTo(Range r) {
		if(a!=r.a) {return a-r.a;}
		return b-r.b;
	}
	
	@Override
	public String toString(){
		return a+(a==b ? "" : ("-"+b));
	}
	
	/**
	 * Appends range representation to a ByteBuilder.
	 * @param bb ByteBuilder to append to
	 * @return The modified ByteBuilder
	 */
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append(a);
		return (b==a ? bb : bb.dash().append(b));
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Merges overlapping and adjacent ranges in a list.
	 * Modifies the input list by consolidating touching ranges and removing nulls.
	 *
	 * @param ranges List of ranges to merge
	 * @param sort Whether to sort the list first
	 * @return Number of ranges that were merged and removed
	 */
	public static int mergeList(ArrayList<Range> ranges, boolean sort) {
		if(ranges.size()<2){return 0;}
		if(sort){Collections.sort(ranges);}
		
		Range current=ranges.get(0);
		int removed=0;
		for(int i=1; i<ranges.size(); i++) {
			Range r=ranges.get(i);
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
	 * Tests if a coordinate range includes a specific position.
	 *
	 * @param a1 Range start coordinate
	 * @param b1 Range end coordinate
	 * @param p Position to test
	 * @return true if position is within range bounds (inclusive)
	 */
	public static boolean include(int a1, int b1, int p){
		assert(a1<=b1) : a1+", "+b1+", "+p;
		return p>=a1 && p<=b1;
	}
	
	/**
	 * Tests if one coordinate range completely includes another.
	 *
	 * @param a1 First range start
	 * @param b1 First range end
	 * @param a2 Second range start
	 * @param b2 Second range end
	 * @return true if second range is completely contained within first
	 */
	public static boolean include(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2>=a1 && b2<=b1;
	}
	
	/**
	 * Tests if two coordinate ranges overlap.
	 *
	 * @param a1 First range start
	 * @param b1 First range end
	 * @param a2 Second range start
	 * @param b2 Second range end
	 * @return true if ranges overlap
	 */
	public static boolean overlaps(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	/**
	 * Calculates overlap length between two coordinate ranges.
	 * Returns negative value if ranges don't overlap.
	 *
	 * @param a1 First range start
	 * @param b1 First range end
	 * @param a2 Second range start
	 * @param b2 Second range end
	 * @return Number of overlapping bases, or negative if no overlap
	 */
	public static int overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		if(a1>a2) {return overlap(a2, b2, a1, b1);}
		assert(a1<=a2);
//		System.err.println("Overlap("+a1+", "+b1+", "+a2+", "+b2+")");
		if(b1<a2) {return b1-a2+1;}//a1a2b1b2 no overlap
		else if(b2<=b1) {return b2-a2+1;}//a1a2b2b1 contained
		return b1-a2+1;//a1a2b1b2 normal overlap
		
	}
	
	/**
	 * Calculates asymmetry between two ranges when aligned.
	 * Measures how much one range extends beyond the other on each side.
	 *
	 * @param a1 First range start
	 * @param b1 First range end
	 * @param a2 Second range start
	 * @param b2 Second range end
	 * @return Measure of range asymmetry
	 */
	public static int lopsidedness(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		if(a1>a2) {return lopsidedness(a2, b2, a1, b1);}
		assert(a1<=a2);
		int left=a2-a1;
		int right=b2-b1;
		return Tools.mid(0, left, right);
	}
	
	/**
	 * Tests if a coordinate range touches a position (overlaps or is adjacent).
	 *
	 * @param a1 Range start coordinate
	 * @param b1 Range end coordinate
	 * @param p Position to test
	 * @return true if position touches or overlaps the range
	 */
	public static boolean touch(int a1, int b1, int p){
		assert(a1<=b1) : a1+", "+b1+", "+p;
		return p<=b1+1 && p>=a1-1;
	}
	
	//Adjacent or overlapping
	/**
	 * Tests if two coordinate ranges touch (overlap or are adjacent).
	 *
	 * @param a1 First range start
	 * @param b1 First range end
	 * @param a2 Second range start
	 * @param b2 Second range end
	 * @return true if ranges touch or overlap
	 */
	public static boolean touch(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1+1 && b2>=a1-1;
	}
	
	/**
	 * Tests if a position is adjacent to a coordinate range.
	 *
	 * @param a1 Range start coordinate
	 * @param b1 Range end coordinate
	 * @param p Position to test
	 * @return true if position is exactly adjacent to range boundary
	 */
	public static boolean adjacent(int a1, int b1, int p){
		assert(a1<=b1) : a1+", "+b1+", "+p;
		return p==b1+1 || p==a1-1;
	}
	
	//Touch but don't overlap
	/**
	 * Tests if two coordinate ranges are adjacent but don't overlap.
	 *
	 * @param a1 First range start
	 * @param b1 First range end
	 * @param a2 Second range start
	 * @param b2 Second range end
	 * @return true if ranges touch but don't overlap
	 */
	public static boolean adjacent(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2==b1+1 || b2==a1-1;
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/*--------------------------------------------------------------*/
	
	/** This should produce the minimal number of non-overlapping ranges,
	 * such that each range has a list of all features that fully contain the range. */
	public static <K extends Feature> Range[] toRanges(ArrayList<K> features){
		if(features==null || features.isEmpty()) {
			assert(false) : "No features found.";//TODO: OK for some sequences; disable after testing.
			return null;
		}
		IntList starts=new IntList(features.size());
		IntList stops=new IntList(features.size());
		for(K f : features) {
			if(f.stop()>=f.start()) {//If it has nonzero length
				starts.add(f.start());
				stops.add(f.stop());
			}else {
				assert(false);
			}
		}
		starts.sort();
		stops.sort();
		starts.shrinkToUnique();
		stops.shrinkToUnique();
		ArrayList<Range> ranges=new ArrayList<Range>(starts.size+stops.size);
		int i=0, j=0;
		int openFeatures=0;//This does not work correctly due to shrinkToUnique; simpler to ignore it.
		int lastPoint=-1;
		boolean lastPointWasStart=false;
		int a=-1, b=-1;
		while(i<starts.size && j<stops.size) {
			a=starts.get(i);
			b=stops.get(j);
			if(a<=b) {
				//New start; make a new range as needed
				if(lastPoint>=0) {
					Range r=null;
					if(lastPointWasStart) {
						if(lastPoint<=a-1) {
//							System.err.println("A: Adding "+(lastPoint)+", "+(a-1)+"; a="+a+", b="+b+", last="+lastPoint+", start="+lastPointWasStart);
							r=new Range(lastPoint, a-1);//right exclusive
						}
					}else {
						if(lastPoint+1<=a-1) {
//							System.err.println("B: Adding "+(lastPoint+1)+", "+(a-1)+"; a="+a+", b="+b+", last="+lastPoint+", start="+lastPointWasStart+"; \tnon-feature");
							r=new Range(lastPoint+1, a-1);//fully exclusive
						}
					}
					if(r!=null) {ranges.add(r);}
				}
				lastPoint=a;
				lastPointWasStart=true;
				openFeatures++;
				i++;
			}else {
				//New stop; always make a new range
				assert(lastPoint>=0);
				Range r=null;
				if(lastPointWasStart) {
					if(lastPoint<=b) {
//						System.err.println("C: Adding "+(lastPoint)+", "+(b)+"; a="+a+", b="+b+", last="+lastPoint+", start="+lastPointWasStart+";  \tfeature");
						r=new Range(lastPoint, b);//fully inclusive
					}
				}else{
					if(lastPoint+1<=b) {
//						System.err.println("D: Adding "+(lastPoint+1)+", "+(b)+"; a="+a+", b="+b+", last="+lastPoint+", start="+lastPointWasStart);
						r=new Range(lastPoint+1, b);//left exclusive
					}
				}
				if(r!=null) {ranges.add(r);}
				lastPoint=b;
				lastPointWasStart=false;
				openFeatures--;
				j++;
			}
		}
		assert(i>=starts.size); //Starts should run out first.
		while(j<stops.size) {
			b=stops.get(j);
			//New stop; always make a new range
			assert(lastPoint>=0);
			Range r;
			if(lastPointWasStart) {
				r=new Range(lastPoint, b);//fully inclusive
			}else{
				try {
					r=new Range(lastPoint+1, b);//left exclusive
				} catch (Throwable e) {
					System.err.println("Error! lastPointWasStart="+lastPointWasStart+", lastPoint="+lastPoint+"\nstarts="+starts+"\nstops="+stops);
					throw new RuntimeException(e);
				}
			}
			ranges.add(r);
			lastPoint=b;
			lastPointWasStart=false;
			openFeatures--;
			j++;
		}
		return ranges.isEmpty() ? null : ranges.toArray(new Range[ranges.size()]);//These are still empty and need to be populated.
	}
	
	/**
	 * Associates features with ranges that overlap them.
	 * Populates each range's obj1 field with a list of overlapping features.
	 * @param ranges Array of ranges to populate
	 * @param features List of features to associate with ranges
	 */
	public static <K extends Feature> void populateRanges(Range[] ranges, ArrayList<K> features) {
		for(K f : features) {
			//TODO: I could maintain a start that occasionally gets incremented rather than
			//binary search each time, since the list is sorted.
			int index=findIndex(f.start(), ranges);
			if(index>-1) {
				for(int i=index; i<ranges.length; i++) {
					Range r=ranges[i];
					if(r.overlaps(f)) {
						ArrayList<K> list=(ArrayList<K>) r.obj1;
						if(list==null) {r.obj1=list=new ArrayList<K>(2);}//Must be pre-allocated (?) ...seems not...
						list.add(f);
					}else if(r.a>f.stop()) {
						break;
					}
				}
			}
		}
	}
	
	/**
	 * Finds the index of the range containing a specific position.
	 * Uses binary search for efficiency on sorted range arrays.
	 *
	 * @param p Position to search for
	 * @param ranges Sorted array of ranges
	 * @return Index of containing range, or negative insertion point if not found
	 */
	public static int findIndex(int p, Range[] ranges) {
		final int x=findIndexBinary(p, ranges);
//		final int y=findIndexLinear(p, ranges);
//		assert(x==y) : p+", "+x+", "+y;
//		assert(x==findIndexLinear(p, ranges)) : p+", "+x+", "+findIndexLinear(p, ranges)+
//			"\n"+Arrays.toString(ranges);
//		assert(x!=y) : p+", "+x+", "+y;//Reminder to remove y.
		return x;
	}
	
	/**
	 * Linear search implementation for finding range index.
	 * Used for debugging and verification against binary search.
	 *
	 * @param p Position to search for
	 * @param ranges Array of ranges
	 * @return Index of containing range, or negative insertion point if not found
	 */
	private static int findIndexLinear(int p, Range[] ranges) {
		if(verbose) {System.err.println("findIndexLinear("+p+", ranges)");}
//		if(ranges==null) {return -1;}
		for(int i=0; i<ranges.length; i++){
			Range r=ranges[i];
			if(verbose) {System.err.println("p="+p+", r="+r);}
			
			if(p<r.a) {
				if(verbose) {System.err.println("returning "+(-i-1));}
				return -i-1;
			}else if(p<=r.b) {
				if(verbose) {System.err.println("returning "+i);}
				return i;
			}
		}
		if(verbose) {System.err.println("returning "+(-ranges.length-1));}
		return -ranges.length-1;
	}
	
	/**
	 * Binary search implementation for finding range index.
	 * More efficient than linear search for large range arrays.
	 *
	 * @param p Position to search for
	 * @param ranges Sorted array of ranges
	 * @return Index of containing range, or negative insertion point if not found
	 */
	private static int findIndexBinary(int p, Range[] ranges) {
//		final int idxl=findIndexLinear(p, ranges);
//		if(ranges==null) {return -1;}
		if(verbose) {System.err.println("findIndexBinary("+p+", ranges)");}
		int a=0, b=ranges.length-1;
		while(a<=b) {
			final int mid=(a+b)/2;
			final Range r=ranges[mid];
			if(verbose) {System.err.println("a="+a+", b="+b+", p="+p+", r="+r);}
			if(p<r.a) {b=mid-1;}
			else if(p>r.b) {a=mid+1;}
			else {
				if(verbose) {System.err.println("Returning "+mid);}
//				assert(mid==idxl) : mid+", "+idxl;
				return mid;
			}
		}
		assert(a>b);//Not found
		final int ret=-max(a, b)-1;
		if(verbose) {System.err.println("Returning "+ret);}
//		assert(ret==idxl) : ret+", "+idxl;
		return ret;//Return negative insert index minus 1
	}
	
	/** List is sorted by start point; may be overlapping.
	 * Returns the leftmost feature starting AFTER this point. */
	public static <K extends Feature> int findIndexAfterBinary(int p, ArrayList<K> list) {
//		if(list==null) {return -1;}
		if(list.isEmpty()) {return -1;}
		int a=0, b=list.size()-1;
		while(a<b) {
			int mid=(a+b)/2;
			K f=list.get(mid);
			if(p<f.start()) {b=mid-1;}
			else{a=mid+1;}
		}
		assert(a==b);
		//Now a should be the index of rightmost element starting with or before p,
		//or the leftmost element after p
		K fm1=(a==0 ? null : list.get(a-1));
		K f=list.get(a);
		K fp1=(a==list.size()-1 ? null : list.get(a+1));
		assert((f.start()>p && (fm1==null || fm1.start()<=p)) || 
				(f.start()<=p && (fp1==null || fp1.start()>p)));//This is the goal at this point
		while(a<list.size() && list.get(a).start()<=p) {a++;}
		return a;
	}
		
//		/** List is sorted by start point; may be overlapping.
//		 * Returns the leftmost feature starting at or after this point. */
//		private static <K extends Feature> int findIndexBinary(int p, ArrayList<K> list) {
////			if(list==null) {return -1;}
//			if(list.isEmpty()) {return -1;}
//			int a=0, b=list.size()-1;
//			while(a<b) {
//				int mid=(a+b)/2;
//				K f=list.get(mid);
//				if(p<=f.start()) {b=mid-1;}
//				else if(p>f.start()) {a=mid+1;}
//				else {return a=b=mid;}
//			}
//			assert(a==b);
//			//Now a should be the index before the first element with index at least p, or 0
//			K f=list.get(a);
//			if(p<f.start()) {
//				a++;
//				if(a>list.size()) {return 
//			}
//			assert(a>b);//Not found
//			return -max(a, b)-1;//Return negative insert index minus 1
//		}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses a string representation into a Range object.
	 * Accepts various formats including single positions and ranges with separators.
	 * @param s String to parse (e.g., "100", "100-200", "[100-200]")
	 * @return Range object representing the parsed coordinates
	 */
	public static Range toRange(String s){
		String[] s2=s.replace("[","").replace("]","").replace("(","").replace(")","").replace(",","").split("-");
		
		int a, b;
		if(s2.length==1){
			a=b=Integer.parseInt(s2[0]);
		}else{
			a=Integer.parseInt(s2[0]);
			b=Integer.parseInt(s2[1]);
		}
		return new Range(a, b);
	}
	
//	@Override
//	public int compareTo(Range other) {
//		if(a<other.a){return -1;}
//		if(a>other.a){return 1;}
//		
//		if(b<other.b){return -1;}
//		if(b>other.b){return 1;}
//		
//		return 0;
//	}
//	
//	public boolean includes(int p){
//		return p>=a && p<=b;
//	}
//	
//	public boolean intersects(int p1, int p2){
//		return overlap(a, b, p1, p2);
//	}
//	
//	public boolean includes(int p1, int p2){
//		assert(p1<=p2);
//		return p1>=a && p2<=b;
//	}
//	
//	public boolean intersects(Range other){
//		return intersects(other.a, other.b);
//	}
//	
//	public boolean touches(Range other){
//		if(intersects(other.a, other.b)){return true;}
//		return b==other.a-1 || a==other.b+1;
//	}
//	
//	public boolean includes(Range other){
//		return includes(other.a, other.b);
//	}
//	
//	@Override
//	public boolean equals(Object other){
//		return equals((Range)other);
//	}
//	
//	public Range merge(Range other){
//		assert(touches(other));
//		Range r=new Range(min(a, other.a), max(b, other.b));
//		
//		assert(r.includes(this));
//		assert(r.includes(other));
//		assert(r.length()<=length()+other.length());
//		return r;
//	}
//	
//	public boolean equals(Range other){
//		return a==other.a && b==other.b;
//	}
//	
//	@Override
//	public int hashCode(){
//		return Long.valueOf(Long.rotateLeft(a, 16)^b).hashCode();
//	}
//	
//	@Override
//	public String toString(){
//		return a+(a==b ? "" : ("-"+b));
//	}
//	
//	public static boolean overlap(int a1, int b1, int a2, int b2){
//		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
//		return a2<=b1 && b2>=a1;
//	}
//	
//	public static Range[] toRanges(int[] ...arrays){
//		int len=0;
//		int[] combined=null;
//
//		if(arrays.length==1){
//			combined=arrays[0];
//			len=combined.length;
//		}else{
//			for(int i=0; i<arrays.length; i++){
//				len+=arrays[i].length;
//			}
//			combined=new int[len];
//			for(int i=0, index=0; i<arrays.length; i++){
//				for(int j=0; j<arrays[i].length; j++){
//					combined[index]=arrays[i][j];
//					index++;
//				}
//			}
//			Arrays.sort(combined);
//		}
//		
//		ArrayList<Range> list=new ArrayList<Range>(16);
//		int start=combined[0], last=combined[0];
//		
////		System.out.println(Arrays.toString(combined));
//		
//		for(int i=0; i<len; i++){
//			int x=combined[i];
//			if(x>last+1){
//				list.add(new Range(start, last));
//				start=last=x;
//			}else{
//				last=x;
//			}
//		}
//		list.add(new Range(start, last));
//		return list.toArray(new Range[list.size()]);
//	}
//	
//	public int length() {
//		return b-a+1;
//	}
//	
//	private static final int min(int x, int y){return x<y ? x : y;}
//	private static final int max(int x, int y){return x>y ? x : y;}
	
	/*--------------------------------------------------------------*/
	
	/** Start coordinate of the range (inclusive) */
	public int a;
	/** End coordinate of the range (inclusive) */
	public int b;

	/** Optional object associated with this range */
	public Object obj1=null;
	/** Debug flag for verbose output during range operations */
	public static final boolean verbose=false;
	
}
