package repeat;

import shared.Tools;
import structures.ByteBuilder;

public class Palindrome implements Comparable<Palindrome>, Cloneable {
	
	public Palindrome(){}
	
	public Palindrome(int a_, int b_, int matches_, int mismatches_){
		set(a_, b_, matches_, mismatches_);
	}
	
	public void set(int a_, int b_, int matches_, int mismatches_){
		a=a_;
		b=b_;
		matches=matches_;
		mismatches=mismatches_;
	}
	
	public void setFrom(Palindrome p){
		assert(p!=this);
		set(p.a, p.b, p.matches, p.mismatches);
	}
	
	public String toString() {
		return "("+a+"-"+b+",matches="+matches+",mismatches="+mismatches+")";
	}
	
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
	
	public int plen() {return matches+mismatches;}
	
	public int length() {return b-a+1;}
	
	public int loop() {return length()-2*matches;}
	
	public Palindrome clear() {
		a=b=matches=mismatches=0;
		return this;
	}
	
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
