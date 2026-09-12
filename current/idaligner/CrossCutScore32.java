package idaligner;

import java.util.Arrays;

/**
 * Prototype score-only CrossCut aligner using three reusable int diagonals.
 * The query is global and the reference ends are free; no traceback or
 * alignment coordinates are retained. One instance is intended per worker
 * thread and is not thread-safe.
 */
public final class CrossCutScore32 {

	/** Returns the best raw score over the final query row. */
	public int score(final byte[] query, final byte[] ref) {
		return score(query, ref, 0, ref.length, true);
	}

	/** Returns the best raw score in an inclusive reference window. */
	public int score(final byte[] query, final byte[] ref, final int refStart,
			final int refEnd) {
		validateRange(ref, refStart, refEnd);
		return score(query, ref, refStart, refEnd-refStart+1, true);
	}

	/** Returns the lower-right matrix cell for equivalence testing. */
	int scoreLastCell(final byte[] query, final byte[] ref) {
		return score(query, ref, 0, ref.length, false);
	}

	private int score(final byte[] query, final byte[] ref, final int refStart,
			final int rLen, final boolean freeRefSuffix) {
		assert(query!=null) : "A query is required.";
		assert(ref!=null) : "A reference is required.";
		final int qLen=query.length;
		if(qLen==0) {return 0;}
		if(rLen==0) {return -qLen;}

		ensureCapacity(rLen);
		int[] diagKm2=diagonalA;
		int[] diagKm1=diagonalB;
		int[] diagK=diagonalC;
		int best=Integer.MIN_VALUE;

		for(int k=2; k<=qLen+rLen; k++) {
			final int minCol=Math.max(1, k-qLen);
			final int maxCol=Math.min(rLen, k-1);
			for(int col=minCol; col<=maxCol; col++) {
				final int row=k-col;
				final int diagValue=(row==1 ? 0 : (col==1 ? -(row-1) : diagKm2[col-2]));
				final int upValue=(row==1 ? 0 : diagKm1[col-1]);
				final int leftValue=(col==1 ? -row : diagKm1[col-2]);
				final int diagonal=diagValue+baseScore(query[row-1], ref[refStart+col-1]);
				final int insertion=upValue-1;
				final int deletion=leftValue-1;
				final int value=Math.max(diagonal, Math.max(insertion, deletion));
				diagK[col-1]=value;
				if(row==qLen && value>best) {best=value;}
			}

			final int[] old=diagKm2;
			diagKm2=diagKm1;
			diagKm1=diagK;
			diagK=old;
		}
		return freeRefSuffix ? best : diagKm1[rLen-1];
	}

	private static void validateRange(final byte[] ref, final int refStart, final int refEnd) {
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		if(refStart<0 || refEnd<refStart || refEnd>=ref.length) {
			throw new IndexOutOfBoundsException("Invalid inclusive reference range "+
					refStart+".."+refEnd+" for length "+ref.length);
		}
	}

	private void ensureCapacity(final int length) {
		if(diagonalA.length>=length) {return;}
		diagonalA=new int[length];
		diagonalB=new int[length];
		diagonalC=new int[length];
	}

	/** CrossCut scoring: canonical match +1, canonical mismatch -1, ambiguity 0. */
	static int baseScore(final byte queryBase, final byte refBase) {
		final int q=baseCode(queryBase);
		final int r=baseCode(refBase);
		return q==r && q<=8 ? 1 : (q>8 || r>8 ? 0 : -1);
	}

	static int baseCode(final byte base) {
		assert(base>=0 && base<CODES.length) : "Base is outside ASCII: "+base;
		return CODES[base];
	}

	private static byte[] makeCodes() {
		final byte[] codes=new byte[128];
		Arrays.fill(codes, (byte)31);
		codes['A']=codes['a']=1;
		codes['C']=codes['c']=2;
		codes['G']=codes['g']=4;
		codes['T']=codes['t']=8;
		codes['U']=codes['u']=8;
		return codes;
	}

	private int[] diagonalA=new int[0];
	private int[] diagonalB=new int[0];
	private int[] diagonalC=new int[0];
	private static final byte[] CODES=makeCodes();
}
