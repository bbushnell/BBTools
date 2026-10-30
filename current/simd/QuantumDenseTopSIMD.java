package simd;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/** Vectorizes QuantumRanker diagonal/up candidates for contiguous dense top rows.
 * The caller keeps the dependent deletion prefix scalar. {@code ref} must hold
 * {@code columns} encoded bases. Previous-row arrays are read from zero through
 * {@code columns}; current-row arrays are written from one through
 * {@code columns}. Column zero is caller-owned boundary state. Current-row
 * output arrays must be mutually distinct and must not alias any input array they
 * could overwrite before a later vector or scalar lane reads it, including
 * {@code ref}. Tied diagonal/up scores keep the diagonal predecessor. An
 * {@code N} in either sequence scores {@value #N_SCORE} and contributes one edit.
 * @contributor Sucrose
 */
public class QuantumDenseTopSIMD{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Allows construction for historical public no-arg compatibility. */
	public QuantumDenseTopSIMD(){}

	/*--------------------------------------------------------------*/
	/*----------------       Static Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Calculates independent diagonal/up candidates for a dense contiguous row.
	 * Outputs start at column one because column zero is owned by the caller's
	 * scalar boundary handling.
	 * @param q Encoded query base for this row
	 * @param ref Encoded reference bases with ref[column-1] addressing
	 * @param columns Number of reference columns to compute
	 * @param prevScore Previous dynamic-programming score row
	 * @param prevMeta Previous metadata row
	 * @param prevEdits Previous edit-count row
	 * @param currScore Current score row
	 * @param currMeta Current metadata row
	 * @param currEdits Current edit-count row
	 */
	public static void diagonalUp(final int q, final int[] ref, final int columns,
			final int[] prevScore, final int[] prevMeta, final int[] prevEdits,
			final int[] currScore, final int[] currMeta, final int[] currEdits){
		assert(columns>0 && columns<=ref.length) :
			"Dense Quantum columns must fit encoded reference: "+columns+" / "+ref.length;
		final IntVector qVector=IntVector.broadcast(SPECIES, q);
		final boolean qIsN=(q=='N');
		final int vectorEnd=(columns/WIDTH)*WIDTH;
		int offset=0;
		for(; offset<vectorEnd; offset+=WIDTH){
			final IntVector refVector=IntVector.fromArray(SPECIES, ref, offset);
			final VectorMask<Integer> matchMask=(qIsN ? SPECIES.maskAll(false) :
					qVector.compare(VectorOperators.EQ, refVector));
			final VectorMask<Integer> nMask=(qIsN ? SPECIES.maskAll(true) :
					refVector.compare(VectorOperators.EQ, (int)'N'));
			final IntVector add=SUB_VECTOR.blend(MATCH_VECTOR, matchMask).blend(N_VECTOR, nMask);
			final IntVector diagonal=IntVector.fromArray(SPECIES, prevScore, offset).add(add);
			final IntVector up=IntVector.fromArray(SPECIES, prevScore, offset+1).add(INS_VECTOR);
			final VectorMask<Integer> useUp=up.compare(VectorOperators.GT, diagonal);
			diagonal.blend(up, useUp).intoArray(currScore, offset+1);

			final IntVector diagonalMeta=IntVector.fromArray(SPECIES, prevMeta, offset);
			final IntVector upMeta=IntVector.fromArray(SPECIES, prevMeta, offset+1);
			diagonalMeta.blend(upMeta, useUp).intoArray(currMeta, offset+1);

			final IntVector diagonalEdits=IntVector.fromArray(SPECIES, prevEdits, offset)
					.add(ONE_VECTOR.blend(ZERO_VECTOR, matchMask));
			final IntVector upEdits=IntVector.fromArray(SPECIES, prevEdits, offset+1)
					.add(ONE_VECTOR);
			diagonalEdits.blend(upEdits, useUp).intoArray(currEdits, offset+1);
		}
		for(int j=offset+1; j<=columns; j++){
			final int r=ref[j-1];
			final boolean match=(q==r && q!='N');
			final boolean hasN=(q=='N' || r=='N');
			final int diagonal=prevScore[j-1]+(match ? MATCH : (hasN ? N_SCORE : SUB));
			final int up=prevScore[j]+INS;
			final boolean useUp=(up>diagonal);
			currScore[j]=(useUp ? up : diagonal);
			currMeta[j]=(useUp ? prevMeta[j] : prevMeta[j-1]);
			currEdits[j]=(useUp ? prevEdits[j]+1 : prevEdits[j-1]+(match ? 0 : 1));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** 256-bit integer vector species used by contiguous row loads. */
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
	/** Vector of match scores. */
	private static final IntVector MATCH_VECTOR=IntVector.broadcast(SPECIES, MATCH);
	/** Vector of substitution scores. */
	private static final IntVector SUB_VECTOR=IntVector.broadcast(SPECIES, SUB);
	/** Vector of insertion scores. */
	private static final IntVector INS_VECTOR=IntVector.broadcast(SPECIES, INS);
	/** Vector of N scores. */
	private static final IntVector N_VECTOR=IntVector.broadcast(SPECIES, N_SCORE);
	/** Vector of zero edit increments. */
	private static final IntVector ZERO_VECTOR=IntVector.zero(SPECIES);
	/** Vector of one edit increment. */
	private static final IntVector ONE_VECTOR=IntVector.broadcast(SPECIES, 1);
}
