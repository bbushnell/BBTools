package test.simd;

import java.util.Arrays;

import shared.Shared;
import simd.Vector;

/** Checks fused-reference alignment padding budgets in scalar and SIMD facade modes.
 * @author Sucrose */
public final class VectorAlignFusedTest{

	/** Runs independent fused-reference alignment checks.
	 * @param args Optional simd=t or simd=f selector */
	public static void main(final String[] args){
		final boolean simd=parseSimd(args);
		if(simd && !Vector.simd128){
			throw new IllegalStateException("Fused alignment SIMD test requires 128/256-bit Vector support");
		}
		final boolean oldSimd=Shared.SIMD;
		int cases=0, failures=0;
		try{
			Shared.SIMD=simd;
			for(final int length : LENGTHS){
				final int[] anchors=anchors(length);
				failures+=checkCase("free-n-"+length, length, 0, FREE_NS,
						take(anchors, FREE_NS), EMPTY); cases++;
				failures+=checkCase("paid-n-"+length, length, MAX_SUBS, FREE_NS,
						take(anchors, FREE_NS+MAX_SUBS), EMPTY); cases++;
				failures+=checkCase("n-overflow-"+length, length, MAX_SUBS, FREE_NS,
						anchors, EMPTY); cases++;
				failures+=checkCase("sub-overflow-"+length, length, MAX_SUBS, FREE_NS,
						EMPTY, take(anchors, MAX_SUBS+1)); cases++;
				failures+=checkConservativeFalsePositive("clip-reject-"+length,
						length, anchors); cases++;
			}
			failures+=checkCase("budget-overflow", 32, MAX_SUBS, Integer.MAX_VALUE,
					EMPTY, EMPTY); cases++;
			System.out.println("ALIGN_FUSED_RESULT simd="+simd+" cases="+cases+
					" failures="+failures);
			if(failures>0){throw new AssertionError("alignFused failures="+failures);}
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Parses the optional SIMD selector.
	 * @param args Command-line arguments
	 * @return True for simd=t */
	private static boolean parseSimd(final String[] args){
		boolean simd=false;
		for(final String arg : args){
			if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		return simd;
	}

	/** Checks one Vector.alignFused case against an independent oracle.
	 * @param label Diagnostic case label
	 * @param length Query length
	 * @param maxSubs Substitution budget
	 * @param maxNs Free reference-N budget
	 * @param nPositions Query positions whose reference bases become N
	 * @param subPositions Query positions whose reference bases become mismatches
	 * @return One on mismatch, otherwise zero */
	private static int checkCase(final String label, final int length, final int maxSubs,
			final int maxNs, final int[] nPositions, final int[] subPositions){
		final byte[] query=filled(length, (byte)'A');
		final byte[] ref=reference(length, nPositions, subPositions);
		final int expected=expectedAlignFused(query, ref, maxSubs, maxNs, RSTART);
		final int actual=Vector.alignFused(query, ref, maxSubs, maxNs, RSTART);
		if(actual==expected){return 0;}
		System.out.println("FAIL alignFused "+label+" maxSubs="+maxSubs+
				" maxNs="+maxNs+" nPositions="+Arrays.toString(nPositions)+
				" subPositions="+Arrays.toString(subPositions)+
				" expected="+expected+" actual="+actual);
		return 1;
	}

	/** Checks a candidate that must pass the prefilter and fail downstream clipping.
	 * @param label Diagnostic case label
	 * @param length Query length
	 * @param anchors Representative positions across vector chunks and the tail
	 * @return One on mismatch, otherwise zero */
	private static int checkConservativeFalsePositive(final String label,
			final int length, final int[] anchors){
		final int[] nPositions=take(anchors, FREE_NS+MAX_SUBS);
		final int[] subPositions=take(reverse(anchors), MAX_SUBS);
		final int failures=checkCase(label, length, MAX_SUBS, FREE_NS,
				nPositions, subPositions);
		final byte[] query=filled(length, (byte)'A');
		final byte[] ref=reference(length, nPositions, subPositions);
		final int subs=expectedAlignFused(query, ref, MAX_SUBS, FREE_NS, RSTART);
		final int clipPenalty=Math.max(0, nPositions.length-FREE_NS);
		final boolean downstreamAccepts=(subs+clipPenalty)<=MAX_SUBS;
		if(!downstreamAccepts){return failures;}
		System.out.println("FAIL clip-reject "+label+" subs="+subs+
				" clipPenalty="+clipPenalty+" maxSubs="+MAX_SUBS);
		return failures+1;
	}

	/** Computes the fused-reference prefilter contract independently.
	 * @param query Query bases
	 * @param ref Reference bases
	 * @param maxSubs Maximum allowed true substitutions
	 * @param maxNs Free reference-N padding count
	 * @param rStart Reference offset for query base zero
	 * @return True substitutions, or maxSubs+1 on prefilter rejection */
	private static int expectedAlignFused(final byte[] query, final byte[] ref,
			final int maxSubs, final int maxNs, final int rStart){
		final long maxPadding=(long)maxNs+maxSubs;
		int subs=0, ns=0;
		for(int q=0, r=rStart; q<query.length; q++, r++){
			final byte rb=ref[r];
			if(rb=='N'){ns++;}
			else if(query[q]!=rb){subs++;}
			if(subs>maxSubs || ns>maxPadding){return maxSubs+1;}
		}
		return subs;
	}

	/** Creates an all-one-base byte array.
	 * @param length Array length
	 * @param base Fill value
	 * @return New filled array */
	private static byte[] filled(final int length, final byte base){
		final byte[] array=new byte[length];
		Arrays.fill(array, base);
		return array;
	}

	/** Creates a padded reference with requested Ns and substitutions.
	 * @param length Query length
	 * @param nPositions Positions to overwrite with N
	 * @param subPositions Positions to overwrite with a non-N mismatch
	 * @return Reference with guard bases around the aligned interval */
	private static byte[] reference(final int length, final int[] nPositions,
			final int[] subPositions){
		final byte[] ref=filled(length+2*RSTART, (byte)'A');
		for(final int pos : subPositions){ref[RSTART+pos]='C';}
		for(final int pos : nPositions){ref[RSTART+pos]='N';}
		return ref;
	}

	/** Returns representative positions for the requested length.
	 * @param length Query length
	 * @return Six positions ordered across the vector and tail boundaries */
	private static int[] anchors(final int length){
		if(length==7){return new int[] {0, 1, 2, 3, 5, 6};}
		if(length==16){return new int[] {0, 1, 7, 8, 14, 15};}
		if(length==17){return new int[] {0, 7, 8, 15, 16};}
		if(length==24){return new int[] {0, 7, 8, 15, 16, 23};}
		if(length==25){return new int[] {0, 8, 15, 16, 23, 24};}
		if(length==32){return new int[] {0, 7, 16, 23, 30, 31};}
		if(length==41){return new int[] {0, 8, 31, 32, 39, 40};}
		throw new IllegalArgumentException("Unexpected length: "+length);
	}

	/** Copies the first n elements from an array.
	 * @param array Source array
	 * @param n Number of elements to copy
	 * @return New prefix copy */
	private static int[] take(final int[] array, final int n){
		return Arrays.copyOf(array, n);
	}

	/** Reverses an int array.
	 * @param array Source array
	 * @return New reversed array */
	private static int[] reverse(final int[] array){
		final int[] copy=array.clone();
		for(int left=0, right=copy.length-1; left<right; left++, right--){
			final int temp=copy[left];
			copy[left]=copy[right];
			copy[right]=temp;
		}
		return copy;
	}

	/** Prevents construction of this standalone test. */
	private VectorAlignFusedTest(){}

	/** Short-scalar, 128-bit, 256-bit and vector-plus-tail query lengths. */
	private static final int[] LENGTHS={7, 16, 17, 24, 25, 32, 41};
	/** Reference start offset used to catch accidental zero-offset assumptions. */
	private static final int RSTART=3;
	/** True-substitution budget. */
	private static final int MAX_SUBS=2;
	/** Free N-padding budget. */
	private static final int FREE_NS=2;
	/** Empty position vector. */
	private static final int[] EMPTY={};
}
