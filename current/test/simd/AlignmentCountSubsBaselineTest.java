package test.simd;

import java.nio.charset.StandardCharsets;

import align2.QuantumRanker;
import idaligner.IDAlignerStatics;
import shared.Shared;
import simd.SIMDAlignByte;
import simd.Vector;
import structures.IntList;

/** Baselines SIMD clipped-offset enumeration and prealignment offset selection.
 * Oracles are scalar and independent of the production helpers under test.
 * This finite fixture is not a throughput benchmark.
 * @author Jean */
public final class AlignmentCountSubsBaselineTest{

	/** Runs bounded alignment/countSubs/QuantumRanker checks with SIMD enabled.
	 * @param args Unused
	 */
	public static void main(final String[] args){
		if(args.length>0){throw new IllegalArgumentException("No arguments expected");}
		if(!Vector.simd256){
			throw new IllegalStateException("Alignment SIMD baseline requires 256-bit Vector support");
		}
		final boolean oldSimd=Shared.SIMD;
		int failures=0, cases=0;
		try{
			Shared.SIMD=true;
			failures+=checkAlignDiagonal("paid_left_overhang",
					bytes("AAAAAAAAAA"), bytes("AAAA"), 5, 1, "-6,-5,-4,-3,-2,-1,0"); cases++;
			failures+=checkAlignDiagonal("lane_boundary_31",
					repeat('A', 4), repeat('A', 31), 0, 0,
					"0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27"); cases++;
			failures+=checkAlignDiagonal("lane_boundary_32",
					repeat('A', 4), repeat('A', 32), 0, 0,
					"0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28"); cases++;
			failures+=checkAlignDiagonal("lane_boundary_33",
					repeat('A', 4), repeat('A', 33), 0, 0,
					"0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29"); cases++;
			failures+=checkAlignDiagonal("undefined_query_base",
					bytes("A$A"), bytes("A$A"), 0, 0, ""); cases++;
			failures+=checkAlignDiagonal("ref_n_is_mismatch",
					bytes("AAA"), bytes("ANACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"), 0, 1, ""); cases++;
			failures+=checkAlignDiagonal("full_budget_128",
					repeat('A', 128), repeat('C', 160), 128, 0, rangeCsv(-127, 159)); cases++;
			failures+=checkAlignDiagonal("normal_150",
					repeat('A', 150), repeat('C', 182), 150, 0, rangeCsv(-149, 181)); cases++;
			failures+=checkAlignDiagonal("full_budget_255",
					repeat('A', 255), repeat('C', 287), 255, 0, rangeCsv(-254, 286)); cases++;
			failures+=checkAlignDiagonal("wrap_127",
					repeat('A', 256), repeat('C', 288), 127, 0, ""); cases++;
			failures+=checkAlignDiagonal("wrap_128",
					repeat('A', 256), repeat('C', 288), 128, 0, ""); cases++;
			failures+=checkAlignDiagonal("wrap_255",
					repeat('A', 256), repeat('C', 288), 255, 0, ""); cases++;

			failures+=checkCountSubs("short_reference_left_edge",
					bytes("AAAA"), bytes("AAAAAA"), -1, 10, -1, 1); cases++;
			failures+=checkCountSubs("short_reference_right_edge",
					bytes("AAAA"), bytes("CCAAAA"), 2, 10, 2, 0); cases++;
			failures+=checkCountSubs("vector_right_edge", bytes("AAAA"),
					bytes("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAA"), 36, 10, 36, 0); cases++;
			failures+=checkCountSubs("min_value_fallback", bytes("ACGT"), bytes("ACGT"),
					Integer.MIN_VALUE, 3, 0, 0); cases++;
			failures+=checkCountSubs("suffix_changes_best_lane",
					bytes("CCCCCACCAAAAACAAAACACACAACCAAACACCCCCCAACCCCACACAC"),
					bytes("AAACAAACCACAACCCACAAACAACCACCAAAACAACACA"),
					8, 40, 1, 29); cases++;
			failures+=checkCountSubs("full_window_argmin_128",
					repeat('A', 128), offsetRef(128, 17, 31), 8, 128, 17, 0); cases++;
			failures+=checkCountSubs("normal_150_boundary",
					repeat('A', 150), offsetRef(150, 20, 31), 28, 150, 20, 0); cases++;
			failures+=checkCountSubs("full_window_argmin_255",
					repeat('A', 255), offsetRef(255, 17, 31), 8, 255, 17, 0); cases++;
			failures+=checkCountSubs("fully_outside_left_window",
					repeat('A', 150), repeat('A', 64), -512, 150, -520, 150); cases++;
			failures+=checkCountSubs("fully_outside_right_window",
					repeat('A', 150), repeat('A', 64), 512, 150, 504, 150); cases++;

			failures+=checkQuantumBandwidth(); cases++;
			System.out.println("ALIGNMENT_COUNTSUBS_RESULT cases="+cases+" failures="+failures);
			if(failures!=0){
				throw new AssertionError("Alignment/countSubs baselines failed: "+failures);
			}
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Compares clipped candidate enumeration with an independent scalar oracle.
	 * @param name Diagnostic case name
	 * @param query Query bases
	 * @param ref Reference bases
	 * @param maxSubs Maximum substitutions or paid clips
	 * @param maxClips Free clipping allowance
	 * @param expectedCsv Literal expected accepted offsets
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkAlignDiagonal(final String name, final byte[] query,
			final byte[] ref, final int maxSubs, final int maxClips,
			final String expectedCsv){
		final IntList actualList=SIMDAlignByte.alignDiagonal(query, ref, maxSubs, maxClips);
		final String actual=toCsv(actualList);
		final String expected=toCsv(expectedAlignments(query, ref, maxSubs, maxClips));
		if(!expected.equals(expectedCsv)){
			throw new AssertionError("Internal literal drift for "+name+": expected literal "+
					expectedCsv+" but oracle produced "+expected);
		}
		if(actual.equals(expected)){return 0;}
		System.out.println("FAIL alignDiagonal "+name+" maxSubs="+maxSubs+
				" maxClips="+maxClips+" expected="+expected+" actual="+actual);
		return 1;
	}

	/**
	 * Builds all scalar accepted offsets over the full paid-overhang domain.
	 * @param query query bases to align
	 * @param ref reference bases to scan
	 * @param maxSubs maximum accepted substitutions or paid clips
	 * @param maxClips free clipping allowance
	 * @return accepted offsets, or null when no offsets are accepted
	 */
	private static IntList expectedAlignments(final byte[] query, final byte[] ref,
			final int maxSubs, final int maxClips){
		IntList list=null;
		if(query.length==0 || ref.length==0){return null;}
		final long budget=(long)maxSubs+maxClips;
		final int first=(int)Math.max(1L-query.length, -budget);
		final int last=(int)Math.min((long)ref.length-1,
				(long)ref.length-query.length+budget);
		for(int pos=first; pos<=last; pos++){
			if(scoreClipped(query, ref, pos, maxClips)<=maxSubs){
				if(list==null){list=new IntList();}
				list.add(pos);
			}
		}
		return list;
	}

	/**
	 * Scores one clipped offset without calling IndelFreeAligner helpers.
	 * @param query query bases to align
	 * @param ref reference bases to scan
	 * @param rStart reference start offset for query base zero
	 * @param maxClips free clipping allowance
	 * @return total substitutions plus paid overhang bases
	 */
	private static int scoreClipped(final byte[] query, final byte[] ref,
			final int rStart, final int maxClips){
		final int rStop=rStart+query.length;
		final int leftClip=Math.max(0, -rStart);
		final int rightClip=Math.max(0, rStop-ref.length);
		final int clips=leftClip+rightClip;
		if(clips>=query.length){return query.length;}
		int subs=Math.max(0, clips-maxClips);
		for(int q=leftClip, r=rStart+leftClip; q<query.length && r<ref.length; q++, r++){
			final byte qb=query[q];
			if(qb!=ref[r] || !isDefined(qb)){subs++;}
		}
		return subs;
	}

	/**
	 * Compares prealignment against the 32-lane expected-offset window.
	 * @param name diagnostic case name
	 * @param query query bases to compare
	 * @param ref reference bases to compare
	 * @param start expected start offset
	 * @param maxSubs maximum substitutions to count before an offset saturates
	 * @param expectedOffset expected selected offset
	 * @param expectedSubs expected substitution count
	 * @return one on mismatch, otherwise zero
	 */
	private static int checkCountSubs(final String name, final byte[] query,
			final byte[] ref, final int start, final int maxSubs,
			final int expectedOffset, final int expectedSubs){
		final int[] pos={start, 0};
		final int[] expected=expectedCountSubs(query, ref, start, maxSubs);
		if(expected[0]!=expectedOffset || expected[1]!=expectedSubs){
			throw new AssertionError("Internal literal drift for "+name+": expected literal "+
					expectedOffset+"/"+expectedSubs+" but oracle produced "+expected[0]+"/"+expected[1]);
		}
		final int actualSubs=IDAlignerStatics.countSubs(query, ref, pos, maxSubs);
		if(pos[0]==expectedOffset && actualSubs==expectedSubs){return 0;}
		System.out.println("FAIL countSubs "+name+" start="+start+" maxSubs="+maxSubs+
				" expected="+expectedOffset+"/"+expectedSubs+
				" actual="+pos[0]+"/"+actualSubs);
		return 1;
	}

	/**
	 * Returns fixed short-reference offsets or the best searched 32-lane offset.
	 * @param query query bases to compare
	 * @param ref reference bases to compare
	 * @param start expected start offset
	 * @param maxSubs maximum substitutions to count before an offset saturates
	 * @return selected offset and substitution count
	 */
	private static int[] expectedCountSubs(final byte[] query, final byte[] ref,
			final int start, final int maxSubs){
		if(ref.length<32 || query.length==0 || start==Integer.MIN_VALUE){
			final int offset=(start==Integer.MIN_VALUE ? 0 : start);
			return new int[] {offset, scoreOffset(query, ref, offset, maxSubs)};
		}
		int bestOffset=start-8, bestSubs=Integer.MAX_VALUE;
		for(int i=0; i<32; i++){
			final int offset=start-8+i;
			final int subs=scoreOffset(query, ref, offset, maxSubs);
			if(subs<bestSubs){
				bestSubs=subs;
				bestOffset=offset;
			}
		}
		return new int[] {bestOffset, bestSubs};
	}

	/**
	 * Scores one offset, counting out-of-window query bases as mismatches.
	 * @param query query bases to compare
	 * @param ref reference bases to compare
	 * @param offset reference offset for query base zero
	 * @param maxSubs maximum substitutions to count before early exit
	 * @return substitutions at the selected offset
	 */
	private static int scoreOffset(final byte[] query, final byte[] ref,
			final int offset, final int maxSubs){
		int subs=0;
		for(int q=0; q<query.length && subs<maxSubs; q++){
			final int r=q+offset;
			if(r<0 || r>=ref.length || query[q]!=ref[r]){subs++;}
		}
		return subs;
	}

	/**
	 * Checks QuantumRanker prealignment metadata; this is not a score-parity oracle.
	 * @return one on mismatch, otherwise zero
	 */
	private static int checkQuantumBandwidth(){
		final byte[] query=bytes("ACGTTGCAAGTCCTAGGACTTGCA");
		final byte[] ref=bytes("TTTTTTTTTTTTACGTTGCCAGTCCTAGGACTTGCATTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		final QuantumRanker.Result scalar=new QuantumRanker().align(query, ref, 0,
				ref.length-1, 8, false);
		final QuantumRanker.Result simd=new QuantumRanker().align(query, ref, 0,
				ref.length-1, 8, true);
		if(scalar.bandwidth==5 && !scalar.vectorized &&
				simd.bandwidth==2 && simd.vectorized){return 0;}
		System.out.println("FAIL QuantumRanker expected scalarBandwidth=5 simdBandwidth=2 "+
				"scalarVectorized=false simdVectorized=true actual="+scalar.bandwidth+"/"+
				simd.bandwidth+"/"+scalar.vectorized+"/"+simd.vectorized);
		return 1;
	}

	/**
	 * Renders accepted offsets for diagnostics and literal checks.
	 * @param list accepted offsets
	 * @return comma-separated offset list
	 */
	private static String toCsv(final IntList list){
		if(list==null || list.size<1){return "";}
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<list.size; i++){
			if(i>0){sb.append(',');}
			sb.append(list.get(i));
		}
		return sb.toString();
	}

	/**
	 * Renders an inclusive offset range for literal checks.
	 * @param first first offset, inclusive
	 * @param last last offset, inclusive
	 * @return comma-separated inclusive offset range
	 */
	private static String rangeCsv(final int first, final int last){
		final StringBuilder sb=new StringBuilder();
		for(int i=first; i<=last; i++){
			if(i>first){sb.append(',');}
			sb.append(i);
		}
		return sb.toString();
	}

	/**
	 * Encodes an ASCII literal as bytes.
	 * @param s ASCII literal
	 * @return US-ASCII bytes for the literal
	 */
	private static byte[] bytes(final String s){return s.getBytes(StandardCharsets.US_ASCII);}

	/**
	 * Builds an ASCII repeated-base array.
	 * @param c base character to repeat
	 * @param length array length
	 * @return byte array filled with the selected base
	 */
	private static byte[] repeat(final char c, final int length){
		final byte[] bases=new byte[length];
		for(int i=0; i<length; i++){bases[i]=(byte)c;}
		return bases;
	}

	/**
	 * Builds a reference with one exact query-sized A block at a searched offset.
	 * @param queryLength query length to match
	 * @param goodOffset offset containing the exact A block
	 * @param lastOffset final offset in the 32-lane search window
	 * @return synthetic C reference with one A match block
	 */
	private static byte[] offsetRef(final int queryLength, final int goodOffset,
			final int lastOffset){
		final byte[] ref=repeat('C', queryLength+lastOffset+1);
		for(int i=0; i<queryLength; i++){
			ref[goodOffset+i]=(byte)'A';
		}
		return ref;
	}

	/**
	 * Returns whether b is one of the ordinary uppercase DNA bases.
	 * @param b base byte
	 * @return true for A, C, G or T
	 */
	private static boolean isDefined(final byte b){
		return b=='A' || b=='C' || b=='G' || b=='T';
	}

	/** Prevents construction of this standalone test. */
	private AlignmentCountSubsBaselineTest(){}
}
