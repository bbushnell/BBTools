package align2;

import java.util.ArrayList;
import dna.Data;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

/** Observes paired candidates to decide whether a wider search is useful and
 * whether its resulting deletion geometry is compatible with the mate.
 * Native matches are expanded operations, not CIGAR strings. Candidate positions
 * use chromosome-array coordinates; temporary intervals are half-open.
 * Observations clone the selected sites but reuse the owner's aligner and cache.
 * This class neither selects an attempt nor changes the index search limits.
 * @author Collei */
public final class HybridPairPolicy {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private HybridPairPolicy(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Observes cloned top sites using the legacy/diagnostic terminal policy. */
	public static Observation observe(BBMapThread owner, Read first, byte[] minus1, byte[] minus2,
			int imperfect1, int max1, int imperfect2, int max2, boolean geometryOnly){
		return observe(owner, first, minus1, minus2, imperfect1, max1, imperfect2, max2, geometryOnly, false);
	}

	/** Generates each mate's native match without modifying its selected SiteScore.
	 * May replace the owner's two match-cache slots and alter aligner scratch state.
	 * @param owner Worker whose aligner generates the disposable matches
	 * @param first Read with a non-null mate
	 * @param minus1 Reverse-complement bases for the first read
	 * @param minus2 Reverse-complement bases for its mate
	 * @param imperfect1 Maximum imperfect alignment score for the first read
	 * @param max1 Maximum alignment score for the first read
	 * @param imperfect2 Maximum imperfect alignment score for the mate
	 * @param max2 Maximum alignment score for the mate
	 * @param geometryOnly True evaluates geometry and suppresses retry reasons
	 * @param terminalHalfWindow Enables terminal windows up to half the read length
	 * @return Per-mate retry flags, or the pair geometry when geometryOnly is true */
	public static Observation observe(BBMapThread owner, Read first, byte[] minus1, byte[] minus2,
			int imperfect1, int max1, int imperfect2, int max2, boolean geometryOnly, boolean terminalHalfWindow){
		assert(first!=null && first.mate!=null) : "Paired retry observation requires both reads";
		final Read second=first.mate;
		final SiteScore a=generate(owner, first, minus1, imperfect1, max1, 0);
		final SiteScore b=generate(owner, second, minus2, imperfect2, max2, 1);
		final int flags1=geometryOnly ? 0 : flags(a, first.length(), terminalHalfWindow);
		final int flags2=geometryOnly ? 0 : flags(b, second.length(), terminalHalfWindow);
		return new Observation(flags1, flags2,
				geometryOnly ? geometry(a, first.length(), b, second.length()) : "NOT_RUN");
	}

	/** Generates a disposable match without populating an observation cache slot. */
	static SiteScore generate(BBMapThread owner, Read read, byte[] minus, int imperfect, int max){
		return generate(owner, read, minus, imperfect, max, -1);
	}
	/** Deep-copies mutable candidate arrays because SiteScore.clone is shallow.
	 * A nonnegative slot admits eligible unchanged results to HybridMatchCache. */
	private static SiteScore generate(BBMapThread owner, Read read, byte[] minus, int imperfect, int max, int slot){
		if(slot>=0 && owner.hybridMatchCache!=null){owner.hybridMatchCache.clear(slot);}
		final SiteScore original=read.topSite();
		if(original==null){return null;}
		final SiteScore copy=original.clone();
		copy.gaps=original.gaps==null ? null : original.gaps.clone();
		copy.match=original.match==null ? null : original.match.clone();
		owner.genMatchStringForSite(read.numericID, copy, read.bases, minus, imperfect, max, read.mate, false);
		if(slot>=0 && owner.hybridMatchCache!=null){
			owner.hybridMatchCache.remember(slot, read.numericID, original, copy, read.bases, minus, imperfect, max);
		}
		return copy;
	}

	/** Applies the legacy/diagnostic policy to an expanded native match. */
	static int flags(SiteScore site, int length){
		return flags(site, length, false);
	}

	/** Returns UNKNOWN before classifying clusters if the match contains N or B. */
	static int flags(SiteScore site, int length, boolean terminalHalfWindow){
		if(site==null){return NO_SITE;}
		if(site.match==null){return UNKNOWN;}
		final int[] counts=new int[5];
		count(site.match, length, !site.plus(), counts);
		if(counts[4]>0){return UNKNOWN;}
		return (signal(counts[2], counts[3], length) ? HALF_ERRORS : 0) |
				((terminalHalfWindow || TERMINAL_RETRY) && terminalSignal(site.match, length, counts[2]+counts[3],
						(terminalHalfWindow || TERMINAL_HALF_WINDOW) ? length/2 : Math.min(32, length/2)) ? TERMINAL_ERRORS : 0);
	}

	/** Counts discrepancies before applying the diagnostic terminal-window policy. */
	static boolean terminalSignal(byte[] match, int length){
		assert(match!=null && length>=2) : "HybridPairPolicy.flags requires a generated native match and two query halves";
		int q=0, total=0;
		for(byte op:match){
			if(op=='D'){continue;}
			if("SVIXYC".indexOf((char)op)>=0){total++;} q++;
		}
		assert(q==length) : "Native operations must consume the query exactly before terminal windows are located: "+q+" versus "+length;
		return terminalSignal(match, length, total, TERMINAL_HALF_WINDOW ? length/2 : Math.min(32, length/2));
	}

	/** Finds a terminal window with at least 4 discrepancies, at least 40% inside,
	 * and at most 10% outside. D consumes reference only and is skipped.
	 * Reuses flags' discrepancy count; no additional full pass or scratch array.
	 * @param length Query length in bases, at least 2
	 * @param total Query-consuming discrepancies across the entire match
	 * @param limit Maximum terminal window length, at most half the query */
	static boolean terminalSignal(byte[] match, int length, int total, int limit){
		assert(match!=null && length>=2 && total>=0 && total<=length && limit>=0 && limit<=length/2) :
				"Native count and terminal limit must describe disjoint query-end windows: length="+length+", errors="+total+", limit="+limit;
		if(total<4){return false;}
		int first=0, last=0, left=0, right=match.length-1;
		for(int n=1; n<=limit; n++){
			while(left<match.length && match[left]=='D'){left++;}
			while(right>=0 && match[right]=='D'){right--;}
			assert(left<=right) : "Reference-only deletions must not consume a query position in either terminal window";
			first+=("SVIXYC".indexOf((char)match[left++])>=0 ? 1 : 0);
			last+=("SVIXYC".indexOf((char)match[right--])>=0 ? 1 : 0);
			if((first>=4 && 5*first>=2*n && 10L*(total-first)<=length-n) ||
					(last>=4 && 5*last>=2*n && 10L*(total-last)<=length-n)){return true;}
		}
		return false;
	}

	/** Adds native-match counts into out; callers must zero it before a new read.
	 * Slots 0/1 count substitutions, 2/3 all discrepancies, and 4 unknown query bases.
	 * Halves are in original query orientation; odd lengths put the extra base in
	 * the second half. Expanded D consumes no query base; all other supported
	 * operations consume one. This method does not clear caller-owned counters. */
	static void count(byte[] match, int length, boolean reverse, int[] out){
		assert(length>=2 && out.length==5) : "NativeHalfSignalProbe contract: two nonempty halves and five counters";
		int q=0;
		for(byte op:match){
			assert("msSVIDXYCNB".indexOf((char)op)>=0) : "Unsupported native operation in paired retry observation: "+op;
			if(op=='D'){continue;}
			assert(q<length) : "Native match cannot consume more than the query length";
			final int position=reverse ? length-1-q : q, half=position<length/2 ? 0 : 1;
			if(op=='S' || op=='V'){out[half]++;}
			if("SVIXYC".indexOf((char)op)>=0){out[2+half]++;}
			if(op=='N' || op=='B'){out[4]++;}
			q++;
		}
		assert(q==length) : "Native match must consume exactly the query for half-read classification";
	}

	/** Tests asymmetric half-read discrepancy rates with at least 12 total errors. */
	static boolean signal(int a, int b, int length){
		assert(length>=2) : "Half-error ratios require two nonempty query halves";
		return a+b>=12 && Math.max(a/(double)(length/2), b/(double)(length-length/2))>=.40
				&& Math.min(a/(double)(length/2), b/(double)(length-length/2))<=.10;
	}

	/** Conservatively rejects long deletions overlapping aligned mate bases.
	 * Unevaluable matches return their reason before any interval comparison. */
	static String geometry(SiteScore a, int length1, SiteScore b, int length2){
		final Shape one=shape(a, length1), two=shape(b, length2);
		if(!one.state.equals("OK")){return one.state;}
		if(!two.state.equals("OK")){return two.state;}
		if(one.chrom!=two.chrom || one.scaffold!=two.scaffold){return "UNEVALUABLE_REFERENCE";}
		return overlaps(one.deletions, two.aligned)+overlaps(two.deletions, one.aligned)>0 ? "CONFLICT" : "NO_CONFLICT";
	}

	/** Converts a generated match to half-open reference intervals. Only deletions
	 * of 101 bases or more enter the conflict policy; shorter deletions still move
	 * the reference cursor. SiteScore stop is inclusive, so consumption ends at stop+1. */
	private static Shape shape(SiteScore site, int length){
		if(site==null){return new Shape("UNEVALUABLE_UNMAPPED");}
		if(site.match==null){return new Shape("UNEVALUABLE_NATIVE_MATCH");}
		if(Data.scaffoldLocs==null || Data.scaffoldLengths==null){return new Shape("UNEVALUABLE_NATIVE_METADATA");}
		final int idx=Data.scaffoldIndex(site.chrom, site.start);
		final int left=Data.scaffoldLocs[site.chrom][idx], size=Data.scaffoldLengths[site.chrom][idx];
		if(site.start<left || site.stop<site.start || (long)site.stop>=(long)left+size){return new Shape("UNEVALUABLE_NATIVE_BOUNDS");}
		final Shape out=new Shape("OK"); out.chrom=site.chrom; out.scaffold=idx;
		int r=site.start, q=0;
		for(byte op:site.match){if("CXY".indexOf((char)op)>=0){return new Shape("UNEVALUABLE_NATIVE_CLIP");}}
		for(int i=0; i<site.match.length; ){
			final byte op=site.match[i];
			assert("msSVNBID".indexOf((char)op)>=0) : "Unsupported native geometry operation: "+op;
			int n=1; while(i+n<site.match.length && site.match[i+n]==op){n++;} i+=n;
			if(op=='D'){
				if(n>SamLine.INTRON_LIMIT){return new Shape("UNEVALUABLE_NATIVE_INTRON");}
				if(n>=101){out.deletions.add(new int[]{r, Math.addExact(r, n)});}
				r=Math.addExact(r, n);
			}else if(op=='I'){q=Math.addExact(q, n);}
			else{out.aligned.add(new int[]{r, Math.addExact(r, n)}); r=Math.addExact(r, n); q=Math.addExact(q, n);}
		}
		assert(q==length && r==(long)site.stop+1) : "Native geometry consumption must match query length and candidate limits";
		return out;
	}

	/** Counts each deletion once if it intersects at least one mate interval. */
	private static int overlaps(ArrayList<int[]> gaps, ArrayList<int[]> blocks){
		int n=0;
		for(int[] gap:gaps){for(int[] block:blocks){if(gap[0]<block[1] && block[0]<gap[1]){n++;break;}}}
		return n;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable observation; flags combines both mates without deciding acceptance. */
	public static final class Observation {
		Observation(int f1, int f2, String g){flags1=f1; flags2=f2; flags=f1|f2; geometry=g;}
		public final int flags1;
		public final int flags2;
		public final int flags;
		public final String geometry;
		/** Legacy retry predicate; general hybrid mode adds its pairing policy in BBMapThread. */
		public boolean route(){return (flags&(NO_SITE|HALF_ERRORS|TERMINAL_ERRORS))!=0;}
	}
	/** Disposable geometry in one chromosome array and one scaffold. */
	private static final class Shape {
		Shape(String s){state=s;}
		final String state;
		int chrom, scaffold;
		final ArrayList<int[]> aligned=new ArrayList<int[]>(), deletions=new ArrayList<int[]>();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Independent retry reasons; UNKNOWN is a guard, not evidence of an error cluster. */
	public static final int NO_SITE=1, HALF_ERRORS=2, UNKNOWN=4, PSEUDO_HALF=8, TERMINAL_ERRORS=16;
	/** Diagnostic compatibility switches; public hybrid mode supplies an explicit policy. */
	private static final boolean TERMINAL_RETRY=Boolean.getBoolean("bbmap3.hybridTerminalRetry");
	private static final boolean TERMINAL_HALF_WINDOW=Boolean.getBoolean("bbmap3.hybridTerminalHalfWindow");
}
