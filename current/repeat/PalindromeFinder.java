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
 * Searches inclusive sequence regions for reverse-complement stems with bounded
 * loops, mismatches and tails. Candidate quality follows Palindrome.compareTo:
 * matches first, then mismatches, span and position; this is not simply maximum span.
 * An instance reuses mutable scratch and is confined to one caller at a time.
 * Searches borrow but do not modify input bases and return independent result copies.
 * Complement lookup returns uppercase bases, so matching is not generally
 * case-insensitive for lowercase input.
 *
 * @author Brian Bushnell
 * @date August 30, 2023
 *
 */
public class PalindromeFinder{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Diagnostic entry point accepting literal sequences, numeric/flag constraints and files.
	 * Loads the last named existing file fully, then reports candidates and histograms.
	 * This is an in-memory diagnostic, not the streaming RepeatFinder launcher. */
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

		for(String s : args){
			if(Character.isDigit(s.charAt(0))){
				maxMismatches=Integer.parseInt(s);
			}else if(s.startsWith("mismatch") || s.startsWith("maxmismatch")){
				maxMismatches=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("minmatch")){
				minMatches=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("loop") || s.startsWith("maxloop")){
				maxLoop=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("minloop")){
				minLoop=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("plen") || s.startsWith("minplen")){
				minPLen=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("maxtaildif")){
				maxTailDif=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("mintail")){
				minTail=Integer.parseInt(s.split("=")[1]);
			}else if(s.startsWith("maxtail")){
				maxTail=Integer.parseInt(s.split("=")[1]);
			}else if(new File(s).exists()){
				fname=s;
			}else{
				sequences.add(s.getBytes());
			}
		}

		if(fname!=null){
			ArrayList<Read> reads=ReadInputStream.toReads(FileFormat.testInput(fname, null, true), -1);
			for(Read r : reads){sequences.add(r.bases);}

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
		for(byte[] s : sequences){
			Palindrome p=pf.longestPalindrome(s);
			if(p!=null){
//				Palindrome p2=pf2.longestPalindrome(s);
//				assert(p2!=null) : "\n"+new String(s)+"\n"+p+"\n"+p2+"\n"+pf+"\n"+pf2;
				found++;
				int tail1=p.a, tail2=s.length-p.b-1;
				if(tail1==tail2){symmetric++;}
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
				if(p.compareTo(best)>0){best=p;}
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

	/**
	 * Constructs a palindrome finder with specified constraints.
	 * Clamps minimum matches to at least one and arm length to at least that count.
	 * Other constraints are retained; assertions check nonnegative ordered loop bounds.
	 *
	 * @param minPLen_ Minimum length of one arm, including mismatches and excluding loop/tails
	 * @param minLoop_ Minimum loop size allowed
	 * @param maxLoop_ Maximum loop size allowed
	 * @param minMatches_ Minimum number of matching base pairs required
	 * @param maxMismatches_ Maximum number of mismatches allowed
	 * @param minTail_ Minimum tail length on each side
	 * @param maxTail_ Maximum tail length on each side
	 * @param maxTailDif_ Maximum difference between left and right tail lengths
	 */
	public PalindromeFinder(int minPLen_, int minLoop_, int maxLoop_,
			int minMatches_, int maxMismatches_, int minTail_, int maxTail_, int maxTailDif_){
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
		halfMinLoopEven=minLoop/2+(minLoop&1);
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

	/**
	 * Finds the highest-ranked qualifying palindrome in the entire sequence.
	 * @param s Nonnull borrowed sequence; no bases are modified
	 * @return Independent result copy, or null if no candidate meets the constraints
	 */
	public Palindrome longestPalindrome(byte[] s){
		return longestPalindrome(s, 0, s.length-1);
	}

	/**
	 * Finds the highest-ranked palindrome within an inclusive sequence region.
	 * Searches for both odd-length and even-length palindromes at each position,
	 * applying all configured constraints including loop size, tail requirements,
	 * and mismatch tolerance.
	 *
	 * @param s The sequence to search for palindromes
	 * @param minPos Minimum position to consider (inclusive)
	 * @param maxPos Maximum position to consider (inclusive)
	 * @return Independent result copy, or null if no scanned candidate meets criteria;
	 * tracker and trackerFull are not updated by this method
	 */
	public Palindrome longestPalindrome(byte[] s, final int minPos, final int maxPos){
		Palindrome best=tempB.clear(), p;
//		System.err.println("longestPalindrome()");
		//REPEAT-004: an even center lies between bases and can start one position before the odd-center bound.
		final int minStart=minPos+halfMinLength-1;
		final int maxStop=maxPos-halfMinLength;
		for(int i=minStart; i<=maxStop; i++){
//			System.err.println();
//			System.err.println("longestPalindrome cycle "+i+" odd");
			//The shorter side bounds all possible pairs outside the smallest allowed loop.
			//Keep equal-match candidates: compareTo also ranks mismatches, span and position.
			final int oddCapacity=Tools.min(i-minPos, maxPos-i)-halfMinLoopOdd+1;
			p=(oddCapacity<Tools.max(minPLen, best.matches) ? null : longestPalindromeOdd(s, i, minPos, maxPos));
			if(p!=null && p.compareTo(best)>0){
				best.setFrom(p);
//				System.err.println("New best: "+best);
			}
//			System.err.println("longestPalindrome cycle "+i+" even");
			final int evenCapacity=Tools.min(i-minPos+1, maxPos-i)-halfMinLoopEven;
			p=(evenCapacity<Tools.max(minPLen, best.matches) ? null : longestPalindromeEven(s, i, minPos, maxPos));
			if(p!=null && p.compareTo(best)>0){
				best.setFrom(p);
//				System.err.println("New best: "+best);
			}
		}
		if(best==null){return null;}
		int loop=best.loop(), plen=best.plen();
		Palindrome ret=(plen<minPLen || loop<minLoop || loop>maxLoop || best.mismatches>maxMismatches) ? null : best.clone();
//		System.err.println("plen="+plen+", minPLen="+minPLen+", loop="+loop+
//				", minLoop="+minLoop+", maxLoop="+maxLoop+", mismatches="+best.mismatches+", maxMismatches="+maxMismatches);
		return ret;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	//Start outside the smallest allowed odd loop, even if its interior happens to match.
	/**
	 * Finds longest odd-length palindrome centered at specified position.
	 * Ignores palindromes in the loop region that may violate minimum loop constraints.
	 *
	 * @param s The sequence to search
	 * @param middle Center position for odd-length palindrome
	 * @param minPos Minimum position boundary
	 * @param maxPos Maximum position boundary
	 * @return Longest odd palindrome found, or null if none
	 */
	private Palindrome longestPalindromeOddIgnoringLoop(byte[] s, int middle, final int minPos, final int maxPos){
		return (maxMismatches<1 ? longestPerfectPalindrome(s, middle-halfMinLoopOdd, middle+halfMinLoopOdd, minPos, maxPos)
				: longestImperfectPalindrome(s, middle-halfMinLoopOdd, middle+halfMinLoopOdd, minPos, maxPos));
	}

	//Start outside the smallest allowed even loop; the offset rounds an odd minimum upward.
	/**
	 * Finds longest even-length palindrome centered at specified position.
	 * Ignores palindromes in the loop region that may violate minimum loop constraints.
	 *
	 * @param s The sequence to search
	 * @param middle Center position for even-length palindrome
	 * @param minPos Minimum position boundary
	 * @param maxPos Maximum position boundary
	 * @return Longest even palindrome found, or null if none
	 */
	private Palindrome longestPalindromeEvenIgnoringLoop(byte[] s, int middle, final int minPos, final int maxPos){
		return (maxMismatches<1 ? longestPerfectPalindrome(s, middle-halfMinLoopEven, middle+halfMinLoopEven+1, minPos, maxPos)
				: longestImperfectPalindrome(s, middle-halfMinLoopEven, middle+halfMinLoopEven+1, minPos, maxPos));
	}

	/**
	 * Finds longest odd-length palindrome centered at specified position.
	 * Uses perfect matching for zero mismatches, imperfect matching otherwise.
	 *
	 * @param s The sequence to search
	 * @param middle Center position for odd-length palindrome
	 * @param minPos Minimum position boundary
	 * @param maxPos Maximum position boundary
	 * @return Longest odd palindrome found, or null if none meets criteria
	 */
	private Palindrome longestPalindromeOdd(byte[] s, int middle, final int minPos, final int maxPos){
		return longestPalindromeOddIgnoringLoop(s, middle, minPos, maxPos);
	}

	/**
	 * Finds longest even-length palindrome centered at specified position.
	 * Uses perfect matching for zero mismatches, imperfect matching otherwise.
	 *
	 * @param s The sequence to search
	 * @param middle Center position for even-length palindrome
	 * @param minPos Minimum position boundary
	 * @param maxPos Maximum position boundary
	 * @return Longest even palindrome found, or null if none meets criteria
	 */
	private Palindrome longestPalindromeEven(byte[] s, int middle, final int minPos, final int maxPos){
		return longestPalindromeEvenIgnoringLoop(s, middle, minPos, maxPos);
	}

	/**
	 * Finds longest perfect palindrome (no mismatches) extending from initial positions.
	 * Expands outward from starting positions, validating all constraints including
	 * tail lengths, loop size, and minimum match requirements.
	 *
	 * @param s The sequence to search
	 * @param a0 Initial left position
	 * @param b0 Initial right position
	 * @param minPos Minimum allowed position
	 * @param maxPos Maximum allowed position
	 * @return Borrowed best scratch candidate (possibly zero-score), or null for an
	 * out-of-bounds initial pair; the next inner search overwrites the scratch
	 */
	private Palindrome longestPerfectPalindrome(final byte[] s, final int a0, final int b0, final int minPos, final int maxPos){
//		System.err.println("longestPerfectPalindrome("+a0+", "+b0+")");
		assert(b0>a0) : a0+", "+b0;
		if(a0<minPos || b0>maxPos){
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
				if(matches>=minMatches && matches>=best.matches){

					p.set(a+1, b-1, matches, mismatches);//This is so I can use the compareTo method
					int tail1=p.a-minPos, tail2=maxPos-p.b;
					int tailDif=Tools.absdif(tail1, tail2);
					int loop=p.loop();
					int plen=p.plen();
					if(tailDif<=maxTailDif && tail1<=maxTail && tail2<=maxTail
							&& tail1>=minTail && tail2>=minTail && loop>=minLoop && loop<=maxLoop
							&& plen>=minPLen && p.mismatches<=maxMismatches){

						//					System.err.println("Considering "+p);
						if(p.compareTo(best)>0){
							best.setFrom(p);
							assert(tailDif<=maxTailDif);
							//						System.err.println("Set to "+best);
						}
					}else{
//						System.err.println(tail1+", "+tail2+", "+
//								tailDif+", "+loop+", "+plen+", "+p.mismatches);
					}
				}
				matches=0;
//				mismatches|=1;
//				lastMismatch=b;
				if(b-a+1>maxLoop){break;}
			}
		}
//		System.err.println("Exit loop: a="+a+", b="+b+", matches="+matches);
		while(a<minPos || b>maxPos){a++; b--;}
		if(matches>=minMatches && matches>=best.matches){
			p.set(a, b, matches, mismatches);
			//			System.err.println("Considering "+p);
			int tail1=p.a-minPos, tail2=maxPos-p.b;
			int tailDif=Tools.absdif(tail1, tail2);
			int loop=p.loop();
			int plen=p.plen();
			if(tailDif<=maxTailDif && tail1<=maxTail && tail2<=maxTail
					&& tail1>=minTail && tail2>=minTail && loop>=minLoop && loop<=maxLoop
					&& plen>=minPLen && p.mismatches<=maxMismatches){
				if(p.compareTo(best)>0){
					assert(tailDif<=maxTailDif);
					best.setFrom(p);
					//				System.err.println("Set to "+best);
				}
			}
		}

//		System.err.println("Returning "+best);
		return best;
	}

	/** Extends an imperfect stem while moving its inner boundary to enforce mismatch limits.
	 * Uses caller-owned bounds and the same borrowed scratch/sentinel convention as
	 * the perfect scanner. Mismatches are interior; accepted terminal pairs match.
	 * TODO: Consider a separate dynamic-programming version with +10/-9 scoring. */
	//Comprehension: each stem arm spans matches+mismatches aligned positions (invariant asserted below: matches+mismatches==a2-a),
	//with mismatches interior and outermost/innermost forced to matches. Thus true loop=length-2*plen, NOT length-2*matches -
	//[repeat/Palindrome#001] corrected the old length-2*matches loop calculation; this remains the default CRISPR path (maxMismatches=2).
	private Palindrome longestImperfectPalindrome(final byte[] s, final int a0, final int b0, final int minPos, final int maxPos){
//		System.err.println("longestPerfectPalindrome("+a0+", "+b0+")");
		assert(b0>a0) : a0+", "+b0;
		if(a0<minPos || b0>maxPos){
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
				for(; mismatches>maxMismatches; a2--, b2++){
					if(matches(s[a2], s[b2])){
						matches--;
					}else{
						mismatches--;
					}
				}
				//Shrink palindrome inner bounds until the outermost base is a match
				for(; mismatches>0 && !matches(s[a2], s[b2]); a2--, b2++){
					mismatches--;
				}
				assert(mismatches<=maxMismatches);
				assert(matches+mismatches==a2-a) : matches+", "+mismatches+", "+a+"-"+a2+", "+b+"-"+b2;
				assert(matches==0 || matches(s[a2], s[b2]));
				if(matches>=minMatches && matches>=best.matches){

					p.set(a+1, b-1, matches, mismatches);//This is so I can use the compareTo method
					while(p.mismatches>0 && !matches(s[p.a], s[p.b])){
						p.mismatches--; p.a++; p.b--;
					}
					assert(matches(s[p.a], s[p.b])) : matches+", "+minMatches;
					int tail1=p.a-minPos, tail2=maxPos-p.b;
					int tailDif=Tools.absdif(tail1, tail2);
					int loop=p.loop();
					int plen=p.plen();
					if(tailDif<=maxTailDif && tail1<=maxTail && tail2<=maxTail
							&& tail1>=minTail && tail2>=minTail && loop>=minLoop && loop<=maxLoop
							&& plen>=minPLen && p.mismatches<=maxMismatches){

						//					System.err.println("Considering "+p);
						if(p.compareTo(best)>0){
							best.setFrom(p);
							//						System.err.println("Set to "+best);
						}
					}else{
//						System.err.println(tail1+", "+tail2+", "+
//								tailDif+", "+loop+", "+plen+", "+p.mismatches);
					}
				}
//				mismatches|=1;
//				lastMismatch=b;
				mismatches++;
				assert(matches==0 || matches+mismatches==a2-a+1) : matches+", "+mismatches+", "+a+"-"+a2+", "+b2+"-"+b;

				//REPEAT-010: an allowed interior mismatch need not enlarge the loop.
				//Only inner bounds move the loop outward; once that gap is too large it cannot shrink again.
				if(b2-a2-1>maxLoop){break;}
			}
		}
//		System.err.println("Exit loop: a="+a+", b="+b+", matches="+matches);
		while(a<minPos || b>maxPos){a++; b--;}
		assert(matches==0 || matches+mismatches==a2-a+1) : matches+", "+mismatches+", "+a+"-"+a2+", "+b2+"-"+b;


		//Shrink palindrome inner bounds until mismatches is under thresh
		for(; mismatches>maxMismatches; a2--, b2++){
			if(matches(s[a2], s[b2])){
				matches--;
			}else{
				mismatches--;
			}
		}
		//Shrink palindrome inner bounds until the outermost base is a match
		for(; mismatches>0 && !matches(s[a2], s[b2]); a2--, b2++){
			mismatches--;
		}
		while(mismatches>0 && !matches(s[a], s[b])){
			a++; b--; mismatches--;
		}
		assert(mismatches<=maxMismatches);
		assert(matches==0 || matches+mismatches==a2-a+1) : matches+", "+mismatches+", "+a+"-"+a2+", "+b2+"-"+b;
		assert(matches==0 || matches(s[a], s[b]));
		assert(matches==0 || matches(s[a2], s[b2]));

		if(matches>=minMatches && matches>=best.matches){
			p.set(a, b, matches, mismatches);
			//			System.err.println("Considering "+p);
			int tail1=p.a-minPos, tail2=maxPos-p.b;
			int tailDif=Tools.absdif(tail1, tail2);
			int loop=p.loop();
			int plen=p.plen();
			if(tailDif<=maxTailDif && tail1<=maxTail && tail2<=maxTail
					&& tail1>=minTail && tail2>=minTail && loop>=minLoop && loop<=maxLoop
					&& plen>=minPLen && p.mismatches<=maxMismatches){
				//			System.err.println("Considering "+p);
				if(p.compareTo(best)>0){
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

	/** Formats configured constraints and derived center-search bounds. */
	@Override
	public String toString(){
		return "minMatches="+minMatches+", maxMismatches="+maxMismatches+", minLoop="+minLoop+
				", maxLoop="+maxLoop+", minPLen="+minPLen+
				", halfMinLoopOdd="+halfMinLoopOdd+", halfMinLoopEven="+halfMinLoopEven+
				", halfMinLength="+halfMinLength+
				", minTail="+minTail+", maxTail="+maxTail+", maxTailDif="+maxTailDif;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Tests if two bases form a matching pair in a palindrome.
	 * Uses the fixed reverse-complement mode; the alternate direct comparison is retained.
	 * Right-side lookup accepts DNA/RNA letters in either case but returns uppercase
	 * complements; left-side comparison is literal. Inputs must be valid lookup indices.
	 *
	 * @param a First base
	 * @param b Second base
	 * @return true if bases match according to current matching mode
	 */
	static boolean matches(byte a, byte b){
		return a==(rcomp ? baseToComp[b] : b);
	}

	/**
	 * Creates lookup table for base complement conversion.
	 * Maps upper/lowercase A/T and C/G to uppercase complements, with U treated as T.
	 * Other ASCII entries use the invalid sentinel '~'; this is not an IUPAC matcher.
	 * @return Array mapping base characters to complement characters
	 */
	static byte[] makeBaseToComp(){
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

	/** Temporary palindrome object for best palindrome tracking */
	private final Palindrome tempB=new Palindrome();
	/** Temporary palindrome object for current palindrome evaluation */
	private final Palindrome tempC=new Palindrome();
	/** Temporary palindrome object for intermediate best palindrome storage */
	private final Palindrome tempD=new Palindrome();
	/** Reserved mismatch-list scratch; current scanners maintain counters instead. */
	private final IntList mismatchList=new IntList();

	/** Caller-populated statistics for selected palindromes; searches do not add automatically. */
	public PalindromeTracker tracker=new PalindromeTracker();
	/** Caller-populated statistics for the caller's full/internal-region subset. */
	public PalindromeTracker trackerFull=new PalindromeTracker();

	/** Maximum number of mismatches allowed in palindromic sequence */
	public final int maxMismatches;
	/** Minimum number of matching base pairs required */
	public final int minMatches;//Should be >=1
	/** Minimum loop size allowed between palindromic arms */
	public final int minLoop;//must be>=0
	/** Maximum loop size allowed between palindromic arms */
	public final int maxLoop;
	/** Minimum arm length including mismatches, excluding the central loop and tails. */
	public final int minPLen;

	/** Minimum tail length required on each side of palindrome */
	public final int minTail;
	/** Maximum tail length allowed on each side of palindrome */
	public final int maxTail;
	/** Maximum difference allowed between left and right tail lengths */
	public final int maxTailDif;

	/** Offset from an odd center to the first pair outside the smallest permitted odd loop. */
	public final int halfMinLoopOdd;
	/** Left offset from an even center to the first pair outside the smallest permitted even loop. */
	public final int halfMinLoopEven;
	/** Half of minimum total palindrome length for search bounds */
	public final int halfMinLength;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Flag indicating whether to use reverse complement matching */
	static final boolean rcomp=true;
	/** Lookup table for converting bases to their complements */
	static final byte[] baseToComp=makeBaseToComp();

}
