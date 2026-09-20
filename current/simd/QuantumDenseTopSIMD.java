package simd;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/** Contiguous int diagonal/up pass for QuantumRanker's dense top rows.
 * The dependent deletion prefix remains scalar in the caller.
 */
public class QuantumDenseTopSIMD {

	public static void diagonalUp(final int q, final int[] ref, final int columns,
			final int[] prevScore, final int[] prevMeta, final int[] prevEdits,
			final int[] currScore, final int[] currMeta, final int[] currEdits){
		assert(columns>0 && columns<=ref.length) :
			"Dense Quantum columns must fit encoded reference: "+columns+" / "+ref.length;
		final IntVector qVector=IntVector.broadcast(SPECIES,q);
		final boolean qIsN=(q=='N');
		final int vectorEnd=(columns/WIDTH)*WIDTH;
		int offset=0;
		for(; offset<vectorEnd; offset+=WIDTH){
			final IntVector refVector=IntVector.fromArray(SPECIES,ref,offset);
			final VectorMask<Integer> matchMask=(qIsN ? SPECIES.maskAll(false) :
				qVector.compare(VectorOperators.EQ,refVector));
			final VectorMask<Integer> nMask=(qIsN ? SPECIES.maskAll(true) :
				refVector.compare(VectorOperators.EQ,(int)'N'));
			final IntVector add=SUB_VECTOR.blend(MATCH_VECTOR,matchMask).blend(N_VECTOR,nMask);
			final IntVector diagonal=IntVector.fromArray(SPECIES,prevScore,offset).add(add);
			final IntVector up=IntVector.fromArray(SPECIES,prevScore,offset+1).add(INS_VECTOR);
			final VectorMask<Integer> useUp=up.compare(VectorOperators.GT,diagonal);
			diagonal.blend(up,useUp).intoArray(currScore,offset+1);
			final IntVector diagonalMeta=IntVector.fromArray(SPECIES,prevMeta,offset);
			final IntVector upMeta=IntVector.fromArray(SPECIES,prevMeta,offset+1);
			diagonalMeta.blend(upMeta,useUp).intoArray(currMeta,offset+1);
			final IntVector diagonalEdits=IntVector.fromArray(SPECIES,prevEdits,offset)
				.add(ONE_VECTOR.blend(ZERO_VECTOR,matchMask));
			final IntVector upEdits=IntVector.fromArray(SPECIES,prevEdits,offset+1).add(ONE_VECTOR);
			diagonalEdits.blend(upEdits,useUp).intoArray(currEdits,offset+1);
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

	private static final VectorSpecies<Integer> SPECIES=IntVector.SPECIES_256;
	private static final int WIDTH=SPECIES.length();
	private static final int MATCH=1, SUB=-1, INS=-1, N_SCORE=0;
	private static final IntVector MATCH_VECTOR=IntVector.broadcast(SPECIES,MATCH);
	private static final IntVector SUB_VECTOR=IntVector.broadcast(SPECIES,SUB);
	private static final IntVector INS_VECTOR=IntVector.broadcast(SPECIES,INS);
	private static final IntVector N_VECTOR=IntVector.broadcast(SPECIES,N_SCORE);
	private static final IntVector ZERO_VECTOR=IntVector.zero(SPECIES);
	private static final IntVector ONE_VECTOR=IntVector.broadcast(SPECIES,1);
}
