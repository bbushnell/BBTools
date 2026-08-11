package prot;

import shared.Tools;

/**
 * A single column node in an {@link AAGraph} protein consensus graph.
 *
 * <p>Holds a per-residue weight and count histogram over the 20 standard amino
 * acids plus X, so the plurality residue (and its allele fraction) can be read
 * off at consensus time. Three node types mirror the nucleotide
 * {@code consensus.BaseNode}: REF (a pivot column), INS (an inserted column
 * hanging off a REF/INS node via {@link #insEdge}), and DEL (a deletion column,
 * which stores only summary counts, no residues).</p>
 *
 * @author Eru
 */
public final class AAGraphNode {

	/** Reference (pivot) column. */
	public static final int REF=0;
	/** Inserted column (relative to the pivot). */
	public static final int INS=1;
	/** Deletion column (no residue, summary counts only). */
	public static final int DEL=2;

	/**
	 * Residue-histogram length: 20 standard amino acids at indices 0-19 and X at
	 * {@link Blosum62#X_CODE} (21). Index 20 (stop) is left unused; proteins never
	 * carry it because {@link Blosum62#encode} rejects '*'.
	 */
	public static final int NAA=Blosum62.X_CODE+1;

	/**
	 * Constructs a graph node.
	 * @param refResidue_ Encoded pivot residue at this column (unused for DEL).
	 * @param type_ One of REF, INS, DEL.
	 * @param rpos_ Reference position this node is anchored to.
	 */
	public AAGraphNode(byte refResidue_, int type_, int rpos_){
		type=type_;
		rpos=rpos_;
		refResidue=refResidue_;
		weight=(type==DEL ? null : new int[NAA]);
		count=(type==DEL ? null : new int[NAA]);
	}

	/**
	 * Adds one aligned member observation to this column.
	 * @param enc Encoded residue (0-19 or X_CODE); ignored for DEL nodes.
	 * @param w Weight of the observation (identity-inverted weight from AAGraph.add).
	 */
	public void add(byte enc, int w){
		if(type!=DEL && enc>=0 && enc<NAA){
			weight[enc]+=w;
			count[enc]++;
		}
		countSum++;
		weightSum+=w;
	}

	/**
	 * Consensus letter for this column: the plurality over REAL residues (0-19). X is never a
	 * valid consensus letter -- an emitted X would mean "members aligned here but their votes were
	 * discarded," which is nonsensical (if nothing aligned here the column would not exist). REF
	 * columns keep the real pivot residue when no real variant clears {@code maf}/{@code minDepth}.
	 * The only way X is returned is the vanishing case where every aligned member had X here and the
	 * pivot residue is not real (a padding column) -- which the depth trim normally removes anyway.
	 *
	 * @param maf Minimum allele fraction (by count) for a real variant to override the pivot.
	 * @param minDepth Minimum count for the winning residue to override the pivot.
	 * @return Encoded consensus residue (a real residue 0-19, or X_CODE only in the all-X case).
	 */
	public byte consensus(float maf, int minDepth){
		assert(type!=DEL) : "consensus() called on a DEL node (no residue histogram).";
		//Plurality over REAL residues (0-19) only; the X seed/observations never win the column.
		int maxPos=-1, maxWeight=0, maxDepth=0;
		for(int i=0; i<=19; i++){
			final int x=weight[i], y=count[i];
			if(x>maxWeight || (x==maxWeight && y>maxDepth)){maxWeight=x; maxDepth=y; maxPos=i;}
		}
		if(maxPos<0){//No real residue observed at all here -- vanishingly rare (all members had X).
			return (refResidue>=0 && refResidue<=19) ? refResidue : Blosum62.X_CODE;
		}
		if(type==REF && refResidue>=0 && refResidue<=19){
			//Keep the real pivot residue unless a real variant is confident enough to replace it.
			final float af=maxDepth/(float)Tools.max(1, countSum);
			if(af<maf || maxDepth<minDepth){return refResidue;}
		}
		//Logically certain: a real residue was found, so the consensus letter is real (0-19), never X.
		assert(maxPos>=0 && maxPos<=19) : "consensus letter must be a real residue; got "+maxPos;
		return (byte)maxPos;
	}

	/**
	 * Column agreement: (max - 0.5*second)/depth over the residue counts.
	 * @return A rough 0-1 measure of how dominant the plurality residue is.
	 */
	public float alleleDif(){
		int max=0, second=0, sum=0;
		for(int x : count){
			if(x>max){second=max; max=x;}
			else if(x>second){second=x;}
			sum+=x;
		}
		return (max-0.5f*second)/(float)Tools.max(1, sum);
	}

	/** Node type: REF, INS, or DEL. */
	public final int type;
	/** Reference position this node is anchored to. */
	public final int rpos;
	/** Encoded pivot residue at this column (unused for DEL). */
	public final byte refResidue;
	/** Per-residue weight histogram (null for DEL). */
	public final int[] weight;
	/** Per-residue count histogram (null for DEL). */
	public final int[] count;
	/** Total observations at this column. */
	public int countSum;
	/** Total weight at this column. */
	public int weightSum;
	/** Next inserted column after this one (the only edge the consensus port needs). */
	public AAGraphNode insEdge;
}
