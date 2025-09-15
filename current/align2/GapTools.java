package align2;

import java.util.ArrayList;
import java.util.Arrays;

import shared.Shared;
import shared.Tools;
import stream.SiteScore;

public class GapTools {
	
	public static int[] fixGaps(SiteScore ss){
		int[] r=fixGaps(ss.start(), ss.stop(), ss.gaps, Shared.MINGAP);
		ss.gaps=r;
		return r;
	}
	
	public static String toString(int[] gaps){
		if(gaps==null){return null;}
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<gaps.length; i++){
			if(i>0){sb.append('~');}
			sb.append(gaps[i]);
		}
		return sb.toString();
	}
	
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
	
	/** This may have some off-by-one errors... */
	public static final int calcGrefLen(SiteScore ss){
		return calcGrefLen(ss.start(), ss.stop(), ss.gaps);
	}
	
	/** This may have some off-by-one errors... */
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
	
	/** TODO: Verify. */
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
	
	/** TODO: Verify. */
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
	
	public static int calcNumGapSymbols(int a, int b){
		assert(b>a);
		int gap=b-a-Shared.GAPBUFFER2;
		return Tools.max(0, gap/Shared.GAPLEN);
	}
	
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
	
	public static final ArrayList<Range> toList(int[] gaps){
		ArrayList<Range> list=new ArrayList<Range>(gaps.length/2);
		for(int i=0; i<gaps.length; i+=2){list.add(new Range(gaps[i], gaps[i+1]));}
		return list;
	}
	
	public static class Range implements Comparable<Range>{
		
		public Range(int a_, int b_){
			assert(b_>=a_);
			a=a_;
			b=b_;
		}
		
		@Override
		public int compareTo(Range r){
			int x;
			x=a-r.a;
			if(x!=0){return x;}
			return b-r.b;
		}
		
		@Override
		public String toString(){
			return "("+a+","+b+")";
		}

		@Override
		public boolean equals(Object other){return equals((Range)other);}
		public boolean equals(Range other){return compareTo(other)==0;}
		
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
