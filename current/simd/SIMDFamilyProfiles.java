package simd;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/** Counts equal nonzero registers in corresponding family-profile buckets.
 * Compares complete stored bit patterns; zero is the empty sentinel, so an
 * empty/empty pair contributes nothing. The caller must supply matching sketch
 * geometry: this class checks array widths, not the meaning of the registers.
 * Arrays are borrowed, never retained or modified. Concurrent calls may share
 * arrays only while their contents remain stable; accumulators are call-local.
 * No converted arrays or result objects are created.
 * Call only behind Shared.SIMD with the Vector API available. This class does
 * not select a scalar fallback or check the runtime gate itself.
 * @author Sayu
 */
public final class SIMDFamilyProfiles{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Static comparison kernels only. */
	private SIMDFamilyProfiles(){}

	/** Sixteen 16-bit registers per vector. */
	private static final VectorSpecies<Short> S=ShortVector.SPECIES_256;
	/** Eight 32-bit registers or widened weights per vector. */
	private static final VectorSpecies<Integer> I=IntVector.SPECIES_256;
	/** Empty-register value and the masked-out weight value. */
	private static final ShortVector ZERO=ShortVector.zero(S);

	/*--------------------------------------------------------------*/
	/*----------------       Static Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts matching nonzero uint16 registers directly from char arrays.
	 * Full vectors and the scalar remainder cover each bucket exactly once.
	 * @param query Nonnull query registers; zero denotes an empty bucket.
	 * @param reference Nonnull reference registers with the same length and geometry.
	 * @return Match count from zero through query.length; empty arrays return zero.
	 */
	public static int matches(final char[] query, final char[] reference){
		assert(query!=null && reference!=null && query.length==reference.length) :
			"Live-tail SIMD compares identical buckets of equal-width immutable sketches";
		final int limit=S.loopBound(query.length);
		int score=0;
		for(int b=0; b<limit; b+=S.length()){
			final ShortVector q=ShortVector.fromCharArray(S, query, b), r=ShortVector.fromCharArray(S, reference, b);
			score+=q.compare(VectorOperators.EQ, r).and(q.compare(VectorOperators.NE, ZERO)).trueCount();
		}
		for(int b=limit; b<query.length; b++){if(query[b]!=0 && query[b]==reference[b]){score++;}}
		assert(score>=0 && score<=query.length) : "Each nonzero uint16 register match contributes once, including scalar tail lanes";
		return score;
	}

	/** Counts full uint32 matches, including equal negative Java int bit patterns.
	 * Prefix agreement alone does not count. Zero remains the empty sentinel.
	 * @param query Nonnull query registers.
	 * @param reference Nonnull reference registers with the same length and geometry.
	 * @return Match count from zero through query.length; empty arrays return zero.
	 */
	public static int matches(int[] query, int[] reference){
		assert(query!=null && reference!=null && query.length==reference.length) : "Wide SIMD compares the same bucket in equal-width profiles";
		final int limit=I.loopBound(query.length);
		int score=0;
		for(int b=0; b<limit; b+=I.length()){
			final IntVector q=IntVector.fromArray(I, query, b), r=IntVector.fromArray(I, reference, b);
			score+=q.compare(VectorOperators.EQ, r).and(q.compare(VectorOperators.NE, 0)).trueCount();
		}
		for(int b=limit; b<query.length; b++){if(query[b]!=0 && query[b]==reference[b]){score++;}}
		assert(score>=0 && score<=query.length) : "Each nonzero full32 match contributes once, including high-bit values and scalar tails";
		return score;
	}

	/** Counts matching nonzero uint16 bit patterns stored in short arrays.
	 * Equality needs no signedness bias; negative shorts retain their high bits.
	 * @param query Nonnull query registers.
	 * @param reference Nonnull reference registers with the same length and geometry.
	 * @return Match count from zero through query.length; empty arrays return zero.
	 */
	public static int matches(short[] query, short[] reference){
		assert(query.length==reference.length) : "Each lane must compare the same bucket in both profiles";
		final int limit=S.loopBound(query.length);
		int score=0;
		for(int b=0; b<limit; b+=S.length()){
			final ShortVector q=ShortVector.fromArray(S, query, b), r=ShortVector.fromArray(S, reference, b);
			score+=q.compare(VectorOperators.EQ, r).and(q.compare(VectorOperators.NE, ZERO)).trueCount();
		}
		for(int b=limit; b<query.length; b++){if(query[b]!=0 && query[b]==reference[b]){score++;}}
		assert(score>=0 && score<=query.length) : "Each nonempty matching bucket contributes exactly one";
		return score;
	}

	/** Sums unsigned16 reference weights at matching nonzero registers.
	 * Every weight is interpreted as weights[b]&amp;65535; a zero weight adds nothing.
	 * Nonmatches and empty/empty pairs add nothing regardless of their weight.
	 * Accumulates in int lanes and reduces once; short-lane sums would overflow.
	 * The length bound keeps the worst-case sum, 32768*65535=2147450880, within int.
	 * @param query Nonnull query registers, at most 32768 buckets.
	 * @param reference Nonnull reference registers with the same length and geometry.
	 * @param weights Nonnull per-bucket reference weights with the same length.
	 * @return Unsigned-weight sum from zero through query.length*65535; empty arrays return zero.
	 */
	public static int weighted(short[] query, short[] reference, short[] weights){
		assert(query.length==reference.length && weights.length==query.length && query.length<=32768) :
			"Weights address the same buckets;32768*65535 is the largest supported signed-int score";
		final int limit=S.loopBound(query.length);
		IntVector low=IntVector.zero(I), high=IntVector.zero(I);
		for(int b=0; b<limit; b+=S.length()){
			final ShortVector q=ShortVector.fromArray(S, query, b), r=ShortVector.fromArray(S, reference, b);
			final VectorMask<Short> match=q.compare(VectorOperators.EQ, r).and(q.compare(VectorOperators.NE, ZERO));
			final ShortVector w=ZERO.blend(ShortVector.fromArray(S, weights, b), match);
			// Widen both halves, then remove sign extension before adding unsigned weights.
			low=low.add(((IntVector)w.convertShape(VectorOperators.S2I, I, 0)).and(65535));
			high=high.add(((IntVector)w.convertShape(VectorOperators.S2I, I, 1)).and(65535));
		}
		int score=low.add(high).reduceLanes(VectorOperators.ADD);
		for(int b=limit; b<query.length; b++){
			if(query[b]!=0 && query[b]==reference[b]){score+=weights[b]&65535;}
		}
		assert(score>=0 && (long)score<=(long)query.length*65535) : "Unsigned weights must not wrap during SIMD reduction";
		return score;
	}
}
