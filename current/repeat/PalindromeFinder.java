package repeat;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.IntList;
import tracker.PalindromeTracker;

/**
 * Finds the longest palindrome in a specified region.
 * Designed for hairpins, so it allows both a bounded loop size,
 * number of mismatches, and tail length.
 * 
 * @author Brian Bushnell
 * @date August 30, 2023
 *
 */
public class PalindromeFinder {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		
		ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		
		int minPLen=1;
		int maxMismatches=0;
		int minMatches=0;
		int minLoop=0;//must be>=0
		int maxLoop=40;
		int minTail=0;
		int maxTail=Integer.MAX_VALUE;
		int maxTailDif=Integer.MAX_VALUE;
		
		String fname=null;
		
		for(String s : args) {
			if(Character.isDigit(s.charAt(0))) {
				maxMismatches=Integer.parseInt(s);
			}else if(s.startsWith("mismatch") || s.startsWith("maxmismatch")) {
				maxMismatches=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("minmatch")) {
				minMatches=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("loop") || s.startsWith("maxloop")) {
				maxLoop=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("minloop")) {
				minLoop=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("plen") || s.startsWith("minplen")) {
				minPLen=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("maxtaildif")) {
				maxTailDif=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("mintail")) {
				minTail=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("maxtail")) {
				maxTail=Integer.parseInt(s.split("=")[1]);
			}else if(new File(s).exists()) {
				fname=s;
			}else {
				sequences.add(s.getBytes());
			}
		}
		
		if(fname!=null) {
			ArrayList<Read> reads=ReadInputStream.toReads(FileFormat.testInput(fname, null, true), -1);
			for(Read r : reads) {sequences.add(r.bases);}
			
		}

		final int histlen=30;
		long found=0;
//		long[] plenHist=new long[histlen+1];
//		long[] loopHist=new long[histlen+1];
//		long[] tailHist=new long[histlen+1];
//		long[] tailDifHist=new long[histlen+1];
//		long[] mismatchHist=new long[histlen+1];
		long symmetric=0;
		Palindrome best=null;
		PalindromeFinder pf=new PalindromeFinder(minPLen, minLoop, maxLoop, minMatches, maxMismatches,
				minTail, maxTail, maxTailDif);
//		assert(false):pf;
//		PalindromeFinder pf2=new PalindromeFinder(minPLen, minLoop-1, maxLoop, maxMismatches);
		for(byte[] s : sequences) {
			Palindrome p=pf.longestPalindrome(s);
			if(p!=null) {
//				Palindrome p2=pf2.longestPalindrome(s);
//				assert(p2!=null) : "\n"+new String(s)+"\n"+p+"\n"+p2+"\n"+pf+"\n"+pf2;
				found++;
				int tail1=p.a, tail2=s.length-p.b-1;
				if(tail1==tail2) {symmetric++;}
//				int tailDif=Tools.absdif(tail1, tail2);
//				assert(tailDif<=maxTailDif) : tail1+", "+tail2+", "+tailDif+", "+maxTailDif+"\n"
//						+ p+", len="+s.length;
//				int loop=p.loop();
//				int plen=p.plen();
//				plenHist[Tools.min(plen, histlen)]++;
//				loopHist[Tools.min(loop, histlen)]++;
//				mismatchHist[Tools.min(p.mismatches, histlen)]++;
//				tailHist[Tools.min(Tools.max(tail1, tail2), histlen)]++;
//				tailDifHist[Tools.min(tailDif, histlen)]++;
				pf.tracker.add(p, 0, s.length-1);
				if(p.compareTo(best)>0) {best=p;}
			}
			
		}
		System.out.println("Longest palindrome is: "+best);
		System.out.println("Palindromes found:     "+found+"/"+sequences.size()+
				" = "+String.format("%.2f%%", found*100.0/sequences.size()));
		System.out.println("Symmetric:             "+symmetric+"/"+found);
		System.out.println("Histogram:");
//		System.out.println("Length\tplen\tloop\ttail\ttaildif\tmismatches");
//		ByteBuilder bb=new ByteBuilder();
//		for(int i=0; i<plenHist.length; i++) {
//			bb.append(i).tab().append(plenHist[i]).tab().append(loopHist[i]).tab();
//			bb.append(tailHist[i]).tab().append(tailDifHist[i]).tab().append(mismatchHist[i]);
//			System.out.println(bb);
//			bb.clear();
//		}
		System.out.println(pf.tracker.toString());
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public PalindromeFinder(int minPLen_, int minLoop_, int maxLoop_, 
			int minMatches_, int maxMismatches_, int minTail_, int maxTail_, int maxTailDif_) {
		minLoop=minLoop_;
		maxLoop=maxLoop_;
		minMatches=Tools.max(1, minMatches_);
		minPLen=Tools.max(minMatches, minPLen_);
		maxMismatches=maxMismatches_;
		minTail=minTail_;
		maxTail=maxTail_;
		maxTailDif=maxTailDif_;
		assert(minLoop>=0);
		assert(maxLoop>=minLoop);
		halfMinLoopEven=minLoop/2;
		halfMinLoopOdd=((minLoop|1)+1)/2;
		int minLength=minLoop+2*minPLen;
		halfMinLength=Tools.max(0, minLength/2);
//		System.err.println("minPLen="+minPLen+", minLoop="+minLoop+", maxLoop="+maxLoop+", maxMismatches="+
//				maxMismatches+", halfMinLoopEven="+halfMinLoopEven+", halfMinLoopOdd="+halfMinLoopOdd);
//		assert(maxMismatches==0) : "maxMismatches!=0 is not yet supported.";
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public Palindrome longestPalindrome(byte[] s){
		return longestPalindrome(s, 0, s.length-1);
	}
	
	public Palindrome longestPalindrome(byte[] s, final int minPos, final int maxPos){
		Palindrome best=tempB.clear(), p;
//		System.err.println("longestPalindrome()");
		final int minStart=minPos+halfMinLength;
		final int maxStop=maxPos-halfMinLength;
		for(int i=minStart; i<=maxStop; i++){
//			System.err.println();
//			System.err.println("longestPalindrome cycle "+i+" odd");
			p=longestPalindromeOdd(s, i, minPos, maxPos);
			if(p!=null && p.compareTo(best)>0){
				best.setFrom(p);
//				System.err.println("New best: "+best);
			}
//			System.err.println("longestPalindrome cycle "+i+" even");
			p=longestPalindromeEven(s, i, minPos, maxPos);
			if(p!=null && p.compareTo(best)>0){
				best.setFrom(p);
//				System.err.println("New best: "+best);
			}
		}
		if(best==null) {return null;}
		int loop=best.loop(), plen=best.plen();
		Palindrome ret=(plen<minPLen || loop<minLoop || loop>maxLoop || best.mismatches>maxMismatches) ? null : best.clone();
//		System.err.println("plen="+plen+", minPLen="+minPLen+", loop="+loop+
//				", minLoop="+minLoop+", maxLoop="+maxLoop+", mismatches="+best.mismatches+", maxMismatches="+maxMismatches);
		return ret;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	//Ignores palindromes in the loop which may violate minLoop
	private Palindrome longestPalindromeOddIgnoringLoop(byte[] s, int middle, final int minPos, final int maxPos){
		return (maxMismatches<1 ? longestPerfectPalindrome(s, middle-halfMinLoopOdd, middle+halfMinLoopOdd, minPos, maxPos)
				: longestImperfectPalindrome(s, middle-halfMinLoopOdd, middle+halfMinLoopOdd, minPos, maxPos));
	}

	//Ignores palindromes in the loop which may violate minLoop
	private Palindrome longestPalindromeEvenIgnoringLoop(byte[] s, int middle, final int minPos, final int maxPos){
		return (maxMismatches<1 ? longestPerfectPalindrome(s, middle-halfMinLoopEven, middle+halfMinLoopEven+1, minPos, maxPos)
				: longestImperfectPalindrome(s, middle-halfMinLoopEven, middle+halfMinLoopEven+1, minPos, maxPos));
	}
	
	private Palindrome longestPalindromeOdd(byte[] s, int middle, final int minPos, final int maxPos){
		return (maxMismatches<1 ? longestPerfectPalindrome(s, middle-1, middle+1, minPos, maxPos)
				: longestImperfectPalindrome(s, middle-1, middle+1, minPos, maxPos));
	}

	private Palindrome longestPalindromeEven(byte[] s, int middle, final int minPos, final int maxPos){
		return (maxMismatches<1 ? longestPerfectPalindrome(s, middle-1, middle+2, minPos, maxPos)
				: longestImperfectPalindrome(s, middle-1, middle+2, minPos, maxPos));
	}
	
	private Palindrome longestPerfectPalindrome(final byte[] s, final int a0, final int b0, final int minPos, final int maxPos){
//		System.err.println("longestPerfectPalindrome("+a0+", "+b0+")");
		assert(b0>a0) : a0+", "+b0;
		if(a0<minPos || b0>maxPos) {
//			System.err.println("Out of bounds: minPos="+minPos+", maxPos="+maxPos);
			return null;
		}
		final Palindrome p=tempC.clear(), best=tempD.clear();
		int a=a0, b=b0;
		int matches=0, mismatches=0;
//		int lastMismatch=-1;
		for(; a>=minPos && b<=maxPos; a--, b++){
//			System.err.println("a="+a+", b="+b+", matches="+matches);
			if(matches(s[a], s[b])){
				matches++;
			}else{
				if(matches>=minMatches && matches>=best.matches) {

					p.set(a+1, b-1, matches, mismatches);//This is so I can use the compareTo method
					int tail1=p.a-minPos, tail2=maxPos-p.b;
					int tailDif=Tools.absdif(tail1, tail2);
					int loop=p.loop();
					int plen=p.plen();
					if(tailDif<=maxTailDif && tail1<=maxTail && tail2<=maxTail
							&& tail1>=minTail && tail2>=minTail && loop>=minLoop && loop<=maxLoop
							&& plen>=minPLen && p.mismatches<=maxMismatches) {

						//					System.err.println("Considering "+p);
						if(p.compareTo(best)>0) {
							best.setFrom(p);
							assert(tailDif<=maxTailDif);
							//						System.err.println("Set to "+best);
						}
					}else {
//						System.err.println(tail1+", "+tail2+", "+
//								tailDif+", "+loop+", "+plen+", "+p.mismatches);
					}
				}
				matches=0;
//				mismatches|=1;
//				lastMismatch=b;
				if(b-a+1>maxLoop) {break;}
			}
		}
//		System.err.println("Exit loop: a="+a+", b="+b+", matches="+matches);
		while(a<minPos || b>maxPos) {a++; b--;}
		if(matches>=minMatches && matches>=best.matches) {
			p.set(a, b, matches, mismatches);
			//			System.err.println("Considering "+p);
			int tail1=p.a-minPos, tail2=maxPos-p.b;
			int tailDif=Tools.absdif(tail1, tail2);
			int loop=p.loop();
			int plen=p.plen();
			if(tailDif<=maxTailDif && tail1<=maxTail && tail2<=maxTail
					&& tail1>=minTail && tail2>=minTail && loop>=minLoop && loop<=maxLoop
					&& plen>=minPLen && p.mismatches<=maxMismatches) {
				if(p.compareTo(best)>0) {
					assert(tailDif<=maxTailDif);
					best.setFrom(p);
					//				System.err.println("Set to "+best);
				}
			}
		}

//		System.err.println("Returning "+best);
		return best;
	}
	
	/** TODO: Make a new dynamic programming version with +10 for match and -9 for mismatch */
	private Palindrome longestImperfectPalindrome(final byte[] s, final int a0, final int b0, final int minPos, final int maxPos){
//		System.err.println("longestPerfectPalindrome("+a0+", "+b0+")");
		assert(b0>a0) : a0+", "+b0;
		if(a0<minPos || b0>maxPos) {
//			System.err.println("Out of bounds: minPos="+minPos+", maxPos="+maxPos);
			return null;
		}
		final Palindrome p=tempC.clear(), best=tempD.clear();
		int a=a0, b=b0;
		int a2=a0, b2=b0;
		int matches=0, mismatches=0;

		for(; a>=minPos && b<=maxPos; a--, b++){
//			System.err.println("a="+a+", b="+b+", matches="+matches);
			if(matches(s[a], s[b])){
				matches++;
			}else{
//				System.err.println(matches+", "+mismatches+", "+a+"-"+a2+", "+b+"-"+b2);
			
				//Shrink palindrome inner bounds until mismatches is under thresh
				for(; mismatches>maxMismatches; a2--, b2++) {
					if(matches(s[a2], s[b2])) {
						matches--;
					}else {
						mismatches--;
					}
				}
				//Shrink palindrome inner bounds until the outermost base is a match
				for(; mismatches>0 && !matches(s[a2], s[b2]); a2--, b2++) {
					mismatches--;
				}
				assert(mismatches<=maxMismatches);
				assert(matches+mismatches==a2-a) : matches+", "+mismatches+", "+a+"-"+a2+", "+b+"-"+b2;
				assert(matches==0 || matches(s[a2],s[b2]));
				if(matches>=minMatches && matches>=best.matches) {

					p.set(a+1, b-1, matches, mismatches);//This is so I can use the compareTo method
					while(p.mismatches>0 && !matches(s[p.a], s[p.b])) {
						p.mismatches--; p.a++; p.b--;
					}
					assert(matches(s[p.a],s[p.b])) : matches+", "+minMatches;
					int tail1=p.a-minPos, tail2=maxPos-p.b;
					int tailDif=Tools.absdif(tail1, tail2);
					int loop=p.loop();
					int plen=p.plen();
					if(tailDif<=maxTailDif && tail1<=maxTail && tail2<=maxTail
							&& tail1>=minTail && tail2>=minTail && loop>=minLoop && loop<=maxLoop
							&& plen>=minPLen && p.mismatches<=maxMismatches) {

						//					System.err.println("Considering "+p);
						if(p.compareTo(best)>0) {
							best.setFrom(p);
							//						System.err.println("Set to "+best);
						}
					}else {
//						System.err.println(tail1+", "+tail2+", "+
//								tailDif+", "+loop+", "+plen+", "+p.mismatches);
					}
				}
//				mismatches|=1;
//				lastMismatch=b;
				mismatches++;
				assert(matches==0 || matches+mismatches==a2-a+1) : matches+", "+mismatches+", "+a+"-"+a2+", "+b2+"-"+b;
				
				if(b-a+1>maxLoop) {break;}
			}
		}
//		System.err.println("Exit loop: a="+a+", b="+b+", matches="+matches);
		while(a<minPos || b>maxPos) {a++; b--;}
		assert(matches==0 || matches+mismatches==a2-a+1) : matches+", "+mismatches+", "+a+"-"+a2+", "+b2+"-"+b;
		
		
		//Shrink palindrome inner bounds until mismatches is under thresh
		for(; mismatches>maxMismatches; a2--, b2++) {
			if(matches(s[a2], s[b2])) {
				matches--;
			}else {
				mismatches--;
			}
		}
		//Shrink palindrome inner bounds until the outermost base is a match
		for(; mismatches>0 && !matches(s[a2], s[b2]); a2--, b2++) {
			mismatches--;
		}
		while(mismatches>0 && !matches(s[a], s[b])) {
			a++; b--; mismatches--;
		}
		assert(mismatches<=maxMismatches);
		assert(matches==0 || matches+mismatches==a2-a+1) : matches+", "+mismatches+", "+a+"-"+a2+", "+b2+"-"+b;
		assert(matches==0 || matches(s[a],s[b]));
		assert(matches==0 || matches(s[a2],s[b2]));
		
		if(matches>=minMatches && matches>=best.matches) {
			p.set(a, b, matches, mismatches);
			//			System.err.println("Considering "+p);
			int tail1=p.a-minPos, tail2=maxPos-p.b;
			int tailDif=Tools.absdif(tail1, tail2);
			int loop=p.loop();
			int plen=p.plen();
			if(tailDif<=maxTailDif && tail1<=maxTail && tail2<=maxTail
					&& tail1>=minTail && tail2>=minTail && loop>=minLoop && loop<=maxLoop
					&& plen>=minPLen && p.mismatches<=maxMismatches) {
				//			System.err.println("Considering "+p);
				if(p.compareTo(best)>0) {
					best.setFrom(p);
					//				System.err.println("Set to "+best);
				}
			}
		}

//		System.err.println("Returning "+best);
		return best;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString() {
		return "minMatches="+minMatches+", maxMismatches="+maxMismatches+", minLoop="+minLoop+
				", maxLoop="+maxLoop+", minPLen="+minPLen+
				", halfMinLoopOdd="+halfMinLoopOdd+", halfMinLoopEven="+halfMinLoopEven+
				", halfMinLength="+halfMinLength+
				", minTail="+minTail+", maxTail="+maxTail+", maxTailDif="+maxTailDif;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	static boolean matches(byte a, byte b) {
		return a==(rcomp ? baseToComp[b] : b);
	}
	
	static byte[] makeBaseToComp() {
		byte[] array=new byte[128];
		Arrays.fill(array, (byte)'~');//Nothing should match this since it is invalid
		array['A']=array['a']='T';
		array['C']=array['c']='G';
		array['G']=array['g']='C';
		array['T']=array['t']=array['U']=array['u']='A';
		return array;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final Palindrome tempB=new Palindrome();
	private final Palindrome tempC=new Palindrome();
	private final Palindrome tempD=new Palindrome();
	private final IntList mismatchList=new IntList();

	public PalindromeTracker tracker=new PalindromeTracker();
	public PalindromeTracker trackerFull=new PalindromeTracker();

	public final int maxMismatches;
	public final int minMatches;//Should be >=1
	public final int minLoop;//must be>=0
	public final int maxLoop;
	public final int minPLen;
	
	public final int minTail;
	public final int maxTail;
	public final int maxTailDif;
	
	public final int halfMinLoopOdd;
	public final int halfMinLoopEven;
	public final int halfMinLength;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	static final boolean rcomp=true;
	static final byte[] baseToComp=makeBaseToComp();
	
}
