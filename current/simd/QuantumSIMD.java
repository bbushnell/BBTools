package simd;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/** Vectorizes QuantumRanker diagonal/up candidates for independent sparse columns.
 * The caller keeps the dependent quantum-jump/frontier phase scalar.
 * {@code active[from..to)} must hold one-based reference columns, and each
 * referenced column {@code j} must have readable {@code ref[j-1]},
 * {@code prevScore[j-1]}, {@code prevScore[j]}, {@code prevMeta[j-1]}, and
 * {@code prevMeta[j]}. Output arrays are indexed by active-list index and must
 * be writable through {@code to-1}. Output arrays must be mutually distinct and
 * must not alias any input array they could overwrite before a later gathered
 * lane reads it, including {@code ref} or {@code active}. Tied diagonal/up
 * scores keep the diagonal predecessor. An {@code N} in either sequence scores
 * {@value #N_SCORE} and is not a match.
 * @contributor Sucrose
 */
public class QuantumSIMD{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Allows construction for historical public no-arg compatibility. */
	public QuantumSIMD(){}

	/*--------------------------------------------------------------*/
	/*----------------       Static Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Calculates independent diagonal/up candidates for active sparse columns.
	 * Results are contiguous by active-list index for the caller's scalar
	 * quantum-jump/frontier pass.
	 * @param q Encoded query base for this row
	 * @param ref Encoded reference bases with ref[column-1] addressing
	 * @param active Dense list of active reference columns
	 * @param from First active-list index to compute; must skip column zero
	 * @param to Exclusive active-list index limit
	 * @param prevScore Previous dynamic-programming score row
	 * @param prevMeta Previous metadata row
	 * @param scoreOut Output scores by active-list index
	 * @param metaOut Output metadata by active-list index
	 */
	public static void diagonalUp(final int q, final int[] ref, final int[] active,
			final int from, final int to, final int[] prevScore, final int[] prevMeta,
			final int[] scoreOut, final int[] metaOut){
		assert(from>=1 && from<=to && to<=active.length) :
			"Active range must exclude column zero and remain in scratch: "+from+"-"+to;
		final int alignedEnd=from+((to-from)/WIDTH)*WIDTH;
		final IntVector qVector=IntVector.broadcast(SPECIES, q);
		final IntVector matchVector=IntVector.broadcast(SPECIES, MATCH);
		final IntVector subVector=IntVector.broadcast(SPECIES, SUB);
		final IntVector nVector=IntVector.broadcast(SPECIES, N_SCORE);
		final IntVector insVector=IntVector.broadcast(SPECIES, INS);
		final boolean qIsN=(q=='N');

		int idx=from;
		for(; idx<alignedEnd; idx+=WIDTH){
			final IntVector refVector=IntVector.fromArray(SPECIES, ref, -1, active, idx);
			final VectorMask<Integer> matchMask=(qIsN ? SPECIES.maskAll(false) :
					qVector.compare(VectorOperators.EQ, refVector));
			final VectorMask<Integer> nMask=(qIsN ? SPECIES.maskAll(true) :
					refVector.compare(VectorOperators.EQ, (int)'N'));
			final IntVector addVector=subVector.blend(matchVector, matchMask).blend(nVector, nMask);
			final IntVector diagonal=IntVector.fromArray(SPECIES, prevScore, -1,
					active, idx).add(addVector);
			final IntVector up=IntVector.fromArray(SPECIES, prevScore, 0,
					active, idx).add(insVector);
			final VectorMask<Integer> useUp=up.compare(VectorOperators.GT, diagonal);
			diagonal.blend(up, useUp).intoArray(scoreOut, idx);

			final IntVector diagonalMeta=IntVector.fromArray(SPECIES, prevMeta, -1,
					active, idx);
			final IntVector upMeta=IntVector.fromArray(SPECIES, prevMeta, 0, active, idx);
			diagonalMeta.blend(upMeta, useUp).intoArray(metaOut, idx);
		}

		for(; idx<to; idx++){
			final int j=active[idx], r=ref[j-1];
			final boolean match=(q==r && q!='N');
			final boolean hasN=(q=='N' || r=='N');
			final int add=match ? MATCH : (hasN ? N_SCORE : SUB);
			final int diagonal=prevScore[j-1]+add;
			final int up=prevScore[j]+INS;
			if(up>diagonal){scoreOut[idx]=up; metaOut[idx]=prevMeta[j];}
			else{scoreOut[idx]=diagonal; metaOut[idx]=prevMeta[j-1];}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** 256-bit integer vector species used by the gather loads. */
	private static final VectorSpecies<Integer> SPECIES=IntVector.SPECIES_256;
	/** Number of integer lanes per vector. */
	private static final int WIDTH=SPECIES.length();
	/** Match score. */
	private static final int MATCH=1;
	/** Substitution score. */
	private static final int SUB=-1;
	/** Insertion score. */
	private static final int INS=-1;
	/** Score for an N in either sequence. */
	private static final int N_SCORE=0;
}
