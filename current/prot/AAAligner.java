package prot;

/**
 * Local (Smith-Waterman) protein aligner with affine gaps and BLOSUM62 scoring.
 *
 * <p>Implements Gotoh's three-state affine-gap recurrence with a local
 * zero-floor and full traceback. Given two BBTools-encoded amino-acid arrays
 * (see {@link Blosum62#encode}), it returns the single best-scoring local
 * alignment (HSP) as an {@link AAAlignment}, or null if no positive-scoring
 * alignment exists.</p>
 *
 * <p>This is a correctness-first implementation: matrices are allocated per
 * call (O(m*n) memory). It is not tuned for speed or huge sequences; the MVP
 * contract prioritizes correctness over speed.</p>
 *
 * @author Eru
 */
public final class AAAligner {

	/** Utility class; no instances. */
	private AAAligner(){}

	/** Large negative sentinel for forbidden matrix states (no overflow risk). */
	private static final int NEG=Integer.MIN_VALUE/4;

	//Traceback states.
	private static final int STATE_M=0;//diagonal match/mismatch
	private static final int STATE_IX=1;//gap in target (query residue vs gap)
	private static final int STATE_IY=2;//gap in query (target residue vs gap)

	/**
	 * Aligns an encoded query against an encoded target and returns the best
	 * local alignment, or null if the best score is not positive.
	 *
	 * @param q Encoded query residues (0-19 or X_CODE).
	 * @param t Encoded target residues (0-19 or X_CODE).
	 * @return Best HSP as an AAAlignment, or null if none scores above zero.
	 */
	public static final AAAlignment align(final byte[] q, final byte[] t){
		final int m=q.length, n=t.length;
		if(m==0 || n==0){return null;}
		final int GO=Blosum62.GAP_OPEN, GE=Blosum62.GAP_EXTEND;
		final int openExtend=GO+GE;//cost of the first position of a gap

		//M=diagonal state; E=gap-in-target (consumes query); F=gap-in-query (consumes target).
		final int[][] M=new int[m+1][n+1];
		final int[][] E=new int[m+1][n+1];
		final int[][] F=new int[m+1][n+1];

		//Borders: local alignment starts fresh (M=0); gap states cannot start an alignment.
		for(int j=0; j<=n; j++){M[0][j]=0; E[0][j]=NEG; F[0][j]=NEG;}
		for(int i=0; i<=m; i++){M[i][0]=0; E[i][0]=NEG; F[i][0]=NEG;}

		int best=0, bestI=0, bestJ=0;
		for(int i=1; i<=m; i++){
			final byte qi=q[i-1];
			final int[] Mi=M[i], Mp=M[i-1], Ei=E[i], Ep=E[i-1], Fi=F[i], Fp=F[i-1];
			for(int j=1; j<=n; j++){
				final int s=Blosum62.score(qi, t[j-1]);
				int diag=Mp[j-1];
				if(Ep[j-1]>diag){diag=Ep[j-1];}
				if(Fp[j-1]>diag){diag=Fp[j-1];}
				if(diag<0){diag=0;}//local zero-floor: allow a fresh start
				Mi[j]=diag+s;

				//Gap in target: open from M[i-1][j] or extend E[i-1][j].
				int e=Mp[j]-openExtend;
				final int e2=Ep[j]-GE;
				if(e2>e){e=e2;}
				Ei[j]=e;

				//Gap in query: open from M[i][j-1] or extend F[i][j-1].
				int f=Mi[j-1]-openExtend;
				final int f2=Fi[j-1]-GE;
				if(f2>f){f=f2;}
				Fi[j]=f;

				if(Mi[j]>best){best=Mi[j]; bestI=i; bestJ=j;}
			}
		}

		if(best<=0){return null;}
		return traceback(q, t, M, E, F, best, bestI, bestJ, openExtend, GE);
	}

	/**
	 * Reconstructs the best alignment by walking the DP matrices from the
	 * maximum-scoring cell back to a local start, tallying identities,
	 * mismatches, gap openings, and end coordinates.
	 */
	private static final AAAlignment traceback(final byte[] q, final byte[] t,
			final int[][] M, final int[][] E, final int[][] F,
			final int best, final int bestI, final int bestJ,
			final int openExtend, final int GE){

		int i=bestI, j=bestJ, state=STATE_M;
		final int qStop=bestI-1, tStop=bestJ-1;//0-based inclusive ends
		int qStart=bestI-1, tStart=bestJ-1;

		int identities=0, mismatches=0, gapOpens=0, length=0;
		boolean prevQGap=false, prevTGap=false;//track gap-run boundaries (right-to-left)

		while(true){
			if(state==STATE_M){
				final byte qa=q[i-1], ta=t[j-1];
				length++;
				if(qa==ta && Blosum62.isStandard(qa)){identities++;}
				else{mismatches++;}
				prevQGap=false; prevTGap=false;
				qStart=i-1; tStart=j-1;//leftmost M column processed = alignment start

				final int s=Blosum62.score(qa, ta);
				final int pred=M[i][j]-s;
				if(pred<=0){break;}//local start reached
				i--; j--;
				if(pred==M[i][j]){state=STATE_M;}
				else if(pred==E[i][j]){state=STATE_IX;}
				else{assert(pred==F[i][j]) : "Traceback inconsistency in M"; state=STATE_IY;}
			}else if(state==STATE_IX){
				//query residue aligned to a gap in target (gap columns are neither
				//identities nor mismatches, so only length and gapOpens change)
				length++;
				if(!prevQGap){gapOpens++;}//new gap run in target (query has a residue here)
				prevQGap=true; prevTGap=false;
				final boolean open=(E[i][j]==M[i-1][j]-openExtend);
				i--;
				state=open ? STATE_M : STATE_IX;
			}else{//STATE_IY: target residue aligned to a gap in query
				length++;
				if(!prevTGap){gapOpens++;}
				prevTGap=true; prevQGap=false;
				final boolean open=(F[i][j]==M[i][j-1]-openExtend);
				j--;
				state=open ? STATE_M : STATE_IY;
			}
		}

		//Gap runs are counted by their left boundary during right-to-left walk; the
		//prevXGap flags reset at each M column, so each maximal run is counted once.
		return new AAAlignment(best, qStart, qStop, tStart, tStop,
			identities, mismatches, gapOpens, length);
	}
}
