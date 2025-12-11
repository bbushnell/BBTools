package align2;

import java.util.ArrayList;
import java.util.Arrays;

import shared.Shared;
import shared.Tools;
import stream.SiteScore;

/**
 * Utility class for handling gaps in sequence alignments within the BBTools framework.
 * Provides methods for gap array processing, validation, merging, and length calculations.
 * Gaps are represented as integer arrays with alternating start/stop coordinates.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class GapTools {
	
	/**
	 * Fixes gap array inconsistencies within a SiteScore object.
	 * Updates the SiteScore's gaps array to ensure proper coordinate ordering and boundaries.
	 * @param ss SiteScore containing gap array to fix
	 * @return Fixed gap array with corrected coordinates
	 */
	public static int[] fixGaps(SiteScore ss){
		int[] r=fixGaps(ss.start(), ss.stop(), ss.gaps, Shared.MINGAP);
		ss.gaps=r;
		return r;
	}
	
	/**
	 * Converts a gap array to string representation with tilde separators.
	 * Returns null if input array is null.
	 * @param gaps Gap array to convert
	 * @return String representation with coordinates separated by tildes, or null
	 */
	public static String toString(int[] gaps){
		if(gaps==null){return null;}
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<gaps.length; i++){
			if(i>0){sb.append('~');}
			sb.append(gaps[i]);
		}
		return sb.toString();
	}
	
	/**
	 * Fixes gap array inconsistencies within specified coordinate boundaries.
	 * Ensures gap coordinates are properly ordered, within bounds, and removes invalid gaps.
	 * Constrains all coordinates to the range [a, b] and maintains monotonic ordering.
	 *
	 * @param a Start coordinate boundary
	 * @param b End coordinate boundary
	 * @param gaps Gap array with alternating start/stop coordinates
	 * @param minGap Minimum gap size threshold
	 * @return Fixed gap array or null if no valid gaps remain
	 */
	public static int[] fixGaps(int a, int b, int[] gaps, int minGap){
//		System.err.println("fixGaps Input: "+a+", "+b+", "+Arrays.toString(gaps)+", "+minGap);
//		assert(false) : "fixGaps called!";
		if(verbose){System.err.println("fixGaps a: "+Arrays.toString(gaps));}
		assert(b>a);
		if(gaps==null){return null;}
		assert(gaps.length>=4);
		if(verbose){System.err.println("fixGaps b: "+Arrays.toString(gaps));}
		
		int g0=gaps[0];
		int gN=gaps[gaps.length-1];
		if(!Tools.overlap(a, b, g0, gN)){return null;}

		int changed=0;
		if(gaps[0]!=a){gaps[0]=a; changed++;}
		if(gaps[gaps.length-1]!=b){gaps[gaps.length-1]=b; changed++;}
		for(int i=0; i<gaps.length; i++){
			if(gaps[i]<a){gaps[i]=a; changed++;}
			else if(gaps[i]>b){gaps[i]=b; changed++;}
		}
		
		if(verbose){System.err.println("fixGaps c0: "+Arrays.toString(gaps));}
		
		for(int i=1; i<gaps.length; i++){
			if(gaps[i-1]>gaps[i]){gaps[i]=gaps[i-1]; changed++;}
		}
		
		if(changed==0){return gaps;}
		
		if(verbose){System.err.println("fixGaps c1: "+Arrays.toString(gaps));}
		
		gaps[0]=a;
		gaps[gaps.length-1]=b;
		if(verbose){System.err.println("fixGaps d: "+Arrays.toString(gaps));}
		
		int remove=0;
		for(int i=0; i<gaps.length; i+=2){
			gaps[i]=Tools.constrict(gaps[i], a, b);
			gaps[i+1]=Tools.constrict(gaps[i+1], a, b);
			if(gaps[i]==gaps[i+1]){remove++;}
		}
		if(verbose){System.err.println("fixGaps e: "+Arrays.toString(gaps));}
		if(remove==0){return gaps;}
		if(verbose){System.err.println("fixGaps f: "+Arrays.toString(gaps));}
		
		return fixGaps2(a, b, gaps, minGap);
	}
	
	/**
	 * Calculates genome reference length accounting for gaps in a SiteScore.
	 * Delegates to coordinate-based calculation method.
	 * @param ss SiteScore containing alignment coordinates and gaps
	 * @return Reference length adjusted for gap compression
	 */
	public static final int calcGrefLen(SiteScore ss){
		return calcGrefLen(ss.start(), ss.stop(), ss.gaps);
	}
	
	/**
	 * Calculates genome reference length accounting for gap compression.
	 * Computes total length minus gap symbol savings based on GAPLEN compression ratio.
	 * May have off-by-one errors as noted in source comments.
	 *
	 * @param a Start coordinate
	 * @param b End coordinate
	 * @param gaps Gap array with alternating coordinates
	 * @return Reference length adjusted for gap compression
	 */
	public static final int calcGrefLen(int a, int b, int[] gaps){
		int total=b-a+1;
		if(gaps==null){return total;}
		for(int i=2; i<gaps.length; i+=2){
			int b1=gaps[i-1];
			int b2=gaps[i];
			int syms=calcNumGapSymbols(b1, b2);
			total=total-syms*(Shared.GAPLEN-1);
		}
		assert(total>0) : "total="+total+", a="+a+", b="+b+", gaps="+Arrays.toString(gaps);
		return total;
	}
	
	/**
	 * Calculates buffer space needed for alignment with gaps.
	 * Accounts for gap compression savings and required buffer padding.
	 * Includes GAPBUFFER2 padding for each gap region.
	 *
	 * @param a Start coordinate
	 * @param b End coordinate
	 * @param gaps Gap array with alternating coordinates
	 * @return Buffer size needed including gap compression and padding
	 */
	public static final int calcBufferNeeded(int a, int b, int[] gaps){
		int total=b-a+1;
		if(gaps==null){return total;}
		for(int i=2; i<gaps.length; i+=2){
			int b1=gaps[i-1];
			int b2=gaps[i];
			int syms=calcNumGapSymbols(b1, b2);
			total=total-syms*(Shared.GAPLEN-1)+Shared.GAPBUFFER2;
		}
		assert(total>0) : a+", "+b+", "+Arrays.toString(gaps);
		return total;
	}
	
	/**
	 * Calculates compressed gap length between two coordinates.
	 * Uses GAPLEN compression for gaps exceeding MINGAP threshold.
	 * Includes GAPBUFFER2 padding plus compressed representation.
	 *
	 * @param a Start coordinate
	 * @param b End coordinate (must be greater than a)
	 * @return Compressed gap length including buffer and symbol compression
	 */
	public static int calcGapLen(int a, int b){
		assert(b>a);
		int gap=b-a;
		if(gap<Shared.MINGAP){return gap;}
		int len=Shared.GAPBUFFER2;
		gap-=Shared.GAPBUFFER2;
		int div=gap/Shared.GAPLEN;
		int rem=gap%Shared.GAPLEN;
		len+=(div+rem);
		return len;
	}
	
	/**
	 * Calculates number of gap symbols needed for coordinate span.
	 * Subtracts GAPBUFFER2 padding and divides by GAPLEN compression ratio.
	 *
	 * @param a Start coordinate
	 * @param b End coordinate (must be greater than a)
	 * @return Number of gap symbols required, minimum 0
	 */
	public static int calcNumGapSymbols(int a, int b){
		assert(b>a);
		int gap=b-a-Shared.GAPBUFFER2;
		return Tools.max(0, gap/Shared.GAPLEN);
	}
	
	/**
	 * Advanced gap fixing algorithm that merges overlapping or closely spaced gaps.
	 * Converts gaps to Range objects, merges ranges separated by less than minGap,
	 * then converts back to gap array format. Handles null range cleanup.
	 *
	 * @param a Start coordinate boundary
	 * @param b End coordinate boundary
	 * @param gaps Gap array to process
	 * @param minGap Minimum spacing between gaps before merging
	 * @return Merged and fixed gap array, or null if insufficient gaps remain
	 */
	public static final int[] fixGaps2(int a, int b, int[] gaps, int minGap){
		if(verbose){System.err.println("Input: "+a+", "+b+", "+Arrays.toString(gaps)+", "+minGap);}
		ArrayList<Range> list=toList(gaps);
		if(verbose){System.err.println("Before fixing: "+list);}
		assert(list.size()>1);
		for(int i=1; i<list.size(); i++){
			Range r1=list.get(i-1);
			Range r2=list.get(i);
			
			if(verbose){
				System.err.println("\nRound "+i);
				System.err.println("r1="+r1);
				System.err.println("r2="+r2);
			}
			
			if(r1!=null){
				if(r2.a-r1.b<=minGap){
					r2.a=Tools.min(r1.a, r2.a);
					r2.b=Tools.max(r1.b, r2.b);
					list.set(i-1, null);
				}
			}
			
			if(verbose){
				System.err.println("->");
				System.err.println(list.get(i-1));
				System.err.println(list.get(i));
			}
			
		}
		if(verbose){System.err.println("After fixing: "+list);}
		Tools.condenseStrict(list);
		if(verbose){System.err.println("After condensing: "+list);}
		
		if(list.size()<2){return null;}
		
		int[] gaps2;
		if(gaps.length==list.size()*2){
			gaps2=gaps;
		}else{
			gaps2=new int[list.size()*2];
		}
		for(int i=0, j=0; i<list.size(); i++, j+=2){
			Range r=list.get(i);
			gaps2[j]=r.a;
			gaps2[j+1]=r.b;
		}
		if(verbose){System.err.println("Final gaps: "+Arrays.toString(gaps2));}
		return gaps2;
	}
	
	/**
	 * Converts gap array to list of Range objects.
	 * Creates Range objects from alternating start/stop coordinates in gaps array.
	 * @param gaps Gap array with alternating start/stop coordinates
	 * @return ArrayList of Range objects representing gaps
	 */
	public static final ArrayList<Range> toList(int[] gaps){
		ArrayList<Range> list=new ArrayList<Range>(gaps.length/2);
		for(int i=0; i<gaps.length; i+=2){list.add(new Range(gaps[i], gaps[i+1]));}
		return list;
	}
	
	/**
	 * Represents a coordinate range with start and end positions.
	 * Used internally for gap processing and merging operations.
	 * Implements Comparable for sorting by start then end coordinates.
	 */
	public static class Range implements Comparable<Range>{
		
		/**
		 * Constructs a Range with specified start and end coordinates.
		 * Asserts that end coordinate is not less than start coordinate.
		 * @param a_ Start coordinate
		 * @param b_ End coordinate (must be >= a_)
		 */
		public Range(int a_, int b_){
			assert(b_>=a_);
			a=a_;
			b=b_;
		}
		
		/**
		 * Compares this Range to another Range for ordering.
		 * Primary sort by start coordinate, secondary sort by end coordinate.
		 * @param r Range to compare against
		 * @return Negative if this < r, positive if this > r, zero if equal
		 */
		@Override
		public int compareTo(Range r){
			int x;
			x=a-r.a;
			if(x!=0){return x;}
			return b-r.b;
		}
		
		/** Returns string representation of Range as (start,end).
		 * @return String in format "(a,b)" where a and b are coordinates */
		@Override
		public String toString(){
			return "("+a+","+b+")";
		}

		/**
		 * Tests equality with another object by delegating to Range-specific equals.
		 * Casts other object to Range before comparison.
		 * @param other Object to compare against
		 * @return true if ranges have identical coordinates, false otherwise
		 */
		@Override
		public boolean equals(Object other){return equals((Range)other);}
		public boolean equals(Range other){return compareTo(other)==0;}
		
		/**
		 * Hash code method that should not be used.
		 * Throws assertion error to prevent Range objects from being hashed.
		 * @return Never returns; always throws AssertionError
		 */
		@Override
		public int hashCode() {
			assert(false) : "This class should not be hashed.";
			return super.hashCode();
		}

		public int a;
		public int b;
	}
	
	public static boolean verbose=false;
	
}
