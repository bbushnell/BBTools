package repeat;

import shared.Tools;
import structures.ByteBuilder;

/** Mutable inclusive palindrome bounds and paired-arm counts.
 * One arm contains matches+mismatches positions; the central loop is excluded.
 * Callers own coordinate/count consistency and synchronize shared mutation.
 * Quality comparison is positive for a better candidate; identity equality is retained.
 */
public class Palindrome implements Comparable<Palindrome>, Cloneable{

	/** Creates a zero-score scratch value; its zero bounds do not denote an empty interval. */
	public Palindrome(){}

	/** Stores inclusive bounds and matching/mismatching position counts for one arm. */
	public Palindrome(int a_, int b_, int matches_, int mismatches_){
		set(a_, b_, matches_, mismatches_);
	}

	/** Replaces all four scalar fields without validating coordinates or counts. */
	public void set(int a_, int b_, int matches_, int mismatches_){
		a=a_;
		b=b_;
		matches=matches_;
		mismatches=mismatches_;
	}

	/** Copies the four scalar fields from a nonnull, distinct source. */
	public void setFrom(Palindrome p){
		assert(p!=this);
		set(p.a, p.b, p.matches, p.mismatches);
	}

	/** Formats inclusive bounds and paired-arm match/mismatch counts. */
	public String toString(){
		return "("+a+"-"+b+",matches="+matches+",mismatches="+mismatches+")";
	}

	/** Adds loop length and sorted left/right tail lengths relative to inclusive a0..b0. */
	public String toString(final int a0, final int b0){
		int tail1=a-a0, tail2=b0-b;
		return "("+a+"-"+b+",matches="+matches+",mismatches="+mismatches+
				",loop="+loop()+",tail1="+Tools.min(tail1, tail2)+
				",tail2="+Tools.max(tail1, tail2)+",taildif="+Tools.absdif(tail1, tail2)+")";
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

	/** Appends bounds, arm length P, matches M, loop L and sorted tails T; returns bb.
	 * Tail coordinates are relative to the caller's inclusive containing region. */
	public ByteBuilder appendTo(ByteBuilder bb, final int a0, final int b0){
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

	/** Returns the length of one arm, including mismatches but excluding the loop/tails. */
	public int plen(){return matches+mismatches;}

	/** Returns the inclusive span; the caller maintains ordered bounds for real candidates. */
	public int length(){return b-a+1;}

	//FIXED [repeat/Palindrome#001] (Brian-greenlit): was length()-2*matches, which overcounted the loop by 2*mismatches on the
	//imperfect/CRISPR path (arm=plen=matches+mismatches per bbcrisprfinder help "plen excludes the loop or tail"). Now length()-2*plen():
	//identical on the perfect path (mismatches=0), and corrects the CRISPR minLoop/maxLoop accept-filter for interior-mismatch palindromes.
	/** Returns the span remaining after both arms; valid candidate counts make it nonnegative. */
	public int loop(){return length()-2*plen();}

	/** Restores the all-zero scratch sentinel and returns this instance. */
	public Palindrome clear(){
		a=b=matches=mismatches=0;
		return this;
	}

	/** Returns an independent copy of these scalar fields using Object.clone. */
	public Palindrome clone(){
		Palindrome p=null;
		try{p=(Palindrome) super.clone();}
		catch(CloneNotSupportedException e){e.printStackTrace();}
		return p;
	}

	/**
	 * Compares quality; natural ascending order places lower-quality candidates first.
	 * Priority order: more matches, fewer mismatches, longer length, leftward position.
	 * @param p Palindrome to compare against
	 * @return Positive if this palindrome is better, negative if worse, zero if equal
	 */
	public int compareTo(Palindrome p){
		if(p==null){return 1;}
		if(matches!=p.matches){return matches-p.matches;}
		if(mismatches!=p.mismatches){return p.mismatches-mismatches;}
		int lenDif=length()-p.length();
		if(lenDif!=0){return lenDif;}
		return p.a-a;
	}

	/** Start position of the palindrome in the sequence */
	public int a=0;
	/** Stop position of the palindrome in the sequence */
	public int b=0;
	/** Number of matching bases in the palindromic arms */
	public int matches=0;
	/** Number of mismatching bases in the palindromic arms */
	public int mismatches=0;

}
