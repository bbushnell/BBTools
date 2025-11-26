package repeat;

import shared.Tools;
import structures.ByteBuilder;

/**
 * Represents a palindromic sequence within DNA/RNA with position and match information.
 * Tracks start/stop positions and counts of matching and mismatching bases.
 * Supports comparison and cloning operations for analysis and sorting.
 * @author Brian Bushnell
 */
public class Palindrome implements Comparable<Palindrome>, Cloneable {
	
	/** Creates an empty palindrome with all values initialized to zero */
	public Palindrome(){}
	
	/**
	 * Creates a palindrome with specified position and match information.
	 *
	 * @param a_ Start position of the palindrome
	 * @param b_ Stop position of the palindrome
	 * @param matches_ Number of matching bases in the palindromic region
	 * @param mismatches_ Number of mismatching bases in the palindromic region
	 */
	public Palindrome(int a_, int b_, int matches_, int mismatches_){
		set(a_, b_, matches_, mismatches_);
	}
	
	/**
	 * Sets the position and match information for this palindrome.
	 *
	 * @param a_ Start position of the palindrome
	 * @param b_ Stop position of the palindrome
	 * @param matches_ Number of matching bases in the palindromic region
	 * @param mismatches_ Number of mismatching bases in the palindromic region
	 */
	public void set(int a_, int b_, int matches_, int mismatches_){
		a=a_;
		b=b_;
		matches=matches_;
		mismatches=mismatches_;
	}
	
	/** Copies all values from another palindrome to this one.
	 * @param p Source palindrome to copy from (must not be this instance) */
	public void setFrom(Palindrome p){
		assert(p!=this);
		set(p.a, p.b, p.matches, p.mismatches);
	}
	
	/**
	 * Returns a string representation showing position range and match statistics.
	 * Format: (start-stop,matches=X,mismatches=Y)
	 * @return String representation of this palindrome
	 */
	public String toString() {
		return "("+a+"-"+b+",matches="+matches+",mismatches="+mismatches+")";
	}
	
	/**
	 * Returns detailed string representation including tail length calculations.
	 * Shows position range, matches, mismatches, loop size, and tail statistics.
	 *
	 * @param a0 Reference start position for tail calculation
	 * @param b0 Reference stop position for tail calculation
	 * @return Detailed string representation with tail analysis
	 */
	public String toString(final int a0, final int b0) {
		int tail1=a-a0, tail2=b0-b;
		return "("+a+"-"+b+",matches="+matches+",mismatches="+mismatches+
				",loop="+loop()+",tail1="+Tools.min(tail1,tail2)+
				",tail2="+Tools.max(tail1,tail2)+",taildif="+Tools.absdif(tail1,tail2)+")";
	}
	
//	public ByteBuilder appendTo(ByteBuilder bb, final int a0, final int b0) {
//		int tail1=a-a0, tail2=b0-b;
//		int tmin=Tools.min(tail1, tail2), tmax=Tools.max(tail1, tail2);
//		int tdif=tail2-tail1;
//		bb.append('(').append(a).dash().append(b).comma();
//		bb.append('m','=').append(matches).comma();
//		bb.append('m','m','=').append(mismatches).comma();
//		bb.append('l','=').append(loop()).comma();
//		bb.append('t','1','=').append(tmin).comma();
//		bb.append('t','2','=').append(tmax).comma();
//		bb.append('t','d','=').append(tdif);
//		return bb.append(')');
//	}
	
	/**
	 * Appends palindrome information to a ByteBuilder in compact format.
	 * Includes position, palindrome length, matches, loop size, and tail info.
	 *
	 * @param bb ByteBuilder to append to
	 * @param a0 Reference start position for tail calculation
	 * @param b0 Reference stop position for tail calculation
	 * @return The ByteBuilder with appended palindrome data
	 */
	public ByteBuilder appendTo(ByteBuilder bb, final int a0, final int b0) {
		int tail1=a-a0, tail2=b0-b;
		int tmin=Tools.min(tail1, tail2), tmax=Tools.max(tail1, tail2);
//		int tdif=tail2-tail1;
		bb.append('(').append(a).dash().append(b).comma();
		bb.append('P','=').append(matches+mismatches).comma();
		bb.append('M','=').append(matches).comma();
		bb.append('L','=').append(loop()).comma();
		bb.append('T','=').append(tmin).plus().append(tmax);
		return bb.append(')');
	}
	
	/** Returns the palindromic length (matches plus mismatches).
	 * @return Total number of bases in the palindromic region */
	public int plen() {return matches+mismatches;}
	
	/** Returns the total length of the palindrome region including loop.
	 * @return Distance from start to stop position plus one */
	public int length() {return b-a+1;}
	
	/**
	 * Calculates the loop size between palindromic arms.
	 * Loop size is total length minus twice the matching region length.
	 * @return Number of bases in the central non-palindromic loop
	 */
	public int loop() {return length()-2*matches;}
	
	/** Resets all fields to zero and returns this instance for chaining.
	 * @return This palindrome instance after clearing */
	public Palindrome clear() {
		a=b=matches=mismatches=0;
		return this;
	}
	
	/** Creates a shallow copy of this palindrome.
	 * @return New palindrome instance with identical field values */
	public Palindrome clone() {
		Palindrome p=null;
		try {p=(Palindrome) super.clone();} 
		catch (CloneNotSupportedException e) {e.printStackTrace();}
		return p;
	}
	
	/** 
	 * The greater of the two will have more matches, or fewer mismatches, 
	 * or be longer, or more to the left.
	 */
	public int compareTo(Palindrome p) {
		if(p==null) {return 1;}
		if(matches!=p.matches) {return matches-p.matches;}
		if(mismatches!=p.mismatches) {return p.mismatches-mismatches;}
		int lenDif=length()-p.length();
		if(lenDif!=0) {return lenDif;}
		return p.a-a;
	}
	
	/** Start location */
	public int a=0;
	/** Stop location */
	public int b=0;
	/** Length of the palindromic sequence, excluding the loop */
	public int matches=0;
	/** Number of mismatches in the palindrome */
	public int mismatches=0;
	
}
