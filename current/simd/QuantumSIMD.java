package simd;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/** Gather-vectorized diagonal/up phase for align2.QuantumRanker. */
public class QuantumSIMD {

	/**
	 * Calculates independent diagonal/up candidates for active sparse columns.
	 * Results are contiguous by active-list index for the caller's scalar
	 * quantum-jump/frontier pass.
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
			final IntVector diagonal=IntVector.fromArray(SPECIES, prevScore, -1, active, idx).add(addVector);
			final IntVector up=IntVector.fromArray(SPECIES, prevScore, 0, active, idx).add(insVector);
			final VectorMask<Integer> useUp=up.compare(VectorOperators.GT, diagonal);
			diagonal.blend(up, useUp).intoArray(scoreOut, idx);

			final IntVector diagonalMeta=IntVector.fromArray(SPECIES, prevMeta, -1, active, idx);
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

	private static final VectorSpecies<Integer> SPECIES=IntVector.SPECIES_256;
	private static final int WIDTH=SPECIES.length();
	private static final int MATCH=1;
	private static final int SUB=-1;
	private static final int INS=-1;
	private static final int N_SCORE=0;
}
