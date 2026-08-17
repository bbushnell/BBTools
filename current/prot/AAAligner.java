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
 * <p>M/E/F matrices are stored flat (row-major, one array per state) in
 * thread-local scratch buffers reused across calls and grown (never shrunk)
 * to the largest sequence pair seen so far on that thread -- zero heap
 * allocation per call once the buffers reach steady-state size. This is a
 * storage-strategy change only; the recurrence and traceback logic are
 * unmodified, so output is byte-identical to the original per-call-allocated
 * implementation.</p>
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
	 * Per-thread reused DP scratch: flat M/E/F buffers (row-major, row length
	 * varies per call and is tracked externally via rowStride) grown to the
	 * largest (m+1)*(n+1) needed so far, plus the best-cell result fields that
	 * used to be a fresh {@code new int[]{best,bestI,bestJ}} per call.
	 */
	private static final class Scratch{
		int[] M=new int[0];
		int[] E=new int[0];
		int[] F=new int[0];
		int capacity=0;
		int best, bestI, bestJ;

		/** Grows M/E/F to at least {@code needed} ints, never shrinking. */
		void ensureCapacity(final int needed){
			if(capacity<needed){
				M=new int[needed];
				E=new int[needed];
				F=new int[needed];
				capacity=needed;
			}
		}
	}
	private static final ThreadLocal<Scratch> SCRATCH=ThreadLocal.withInitial(Scratch::new);

	/**
	 * Aligns an encoded query against an encoded target and returns the best
	 * local alignment, or null if the best score is not positive.
	 *
	 * @param q Encoded query residues (0-19 or X_CODE).
	 * @param t Encoded target residues (0-19 or X_CODE).
	 * @return Best HSP as an AAAlignment, or null if none scores above zero.
	 */
	public static final AAAlignment align(final byte[] q, final byte[] t){
		return align(q, t, false);
	}

	/**
	 * As {@link #align(byte[],byte[])}, but when {@code recordPath} is true the returned
	 * alignment's {@link AAAlignment#match} is populated with the per-column m/D/I op string,
	 * which an alignment/consensus graph needs to add the query as a member.
	 *
	 * @param q Encoded query residues (0-19 or X_CODE).
	 * @param t Encoded target residues (0-19 or X_CODE).
	 * @param recordPath Whether to record the per-column op path into the result's match field.
	 * @return Best HSP as an AAAlignment, or null if none scores above zero.
	 */
	public static final AAAlignment align(final byte[] q, final byte[] t, final boolean recordPath){
		final int m=q.length, n=t.length;
		if(m==0 || n==0){return null;}

		//M=diagonal state; E=gap-in-target (consumes query); F=gap-in-query (consumes target).
		final int rowStride=n+1;
		final Scratch s=SCRATCH.get();
		s.ensureCapacity((m+1)*rowStride);
		final int[] M=s.M, E=s.E, F=s.F;
		fillMatrices(q, t, M, E, F, rowStride, true, s);
		final int best=s.best, bestI=s.bestI, bestJ=s.bestJ;

		if(best<=0){return null;}
		final int openExtend=Blosum62.GAP_OPEN+Blosum62.GAP_EXTEND;
		return traceback(q, t, M, E, F, rowStride, best, bestI, bestJ, openExtend, Blosum62.GAP_EXTEND, recordPath);
	}

	/**
	 * Semi-global ("glocal") alignment: the whole query is placed (global in the query) while the
	 * reference ends are free (the query aligns to any substring of the reference). This is what a
	 * consensus graph needs to add a member: with a padded reference, a member overhanging the seed
	 * maps into the padding and can grow the consensus in both directions. Unlike {@link #align}, it
	 * does not require a positive score (cluster members may be divergent) and returns an alignment
	 * covering the entire query. When {@code recordPath} is true the result's {@link AAAlignment#match}
	 * holds the m/D/I ops.
	 *
	 * @param q Encoded query (member) residues.
	 * @param t Encoded reference residues (typically an X-padded pivot).
	 * @param recordPath Whether to record the per-column op path.
	 * @return Alignment covering the full query, or null only for empty input.
	 */
	public static final AAAlignment alignGlocal(final byte[] q, final byte[] t, final boolean recordPath){
		final int m=q.length, n=t.length;
		if(m==0 || n==0){return null;}
		final int rowStride=n+1;
		final Scratch s=SCRATCH.get();
		s.ensureCapacity((m+1)*rowStride);
		final int[] M=s.M, E=s.E, F=s.F;
		fillMatrices(q, t, M, E, F, rowStride, false, s);
		final int best=s.best, bestI=s.bestI, bestJ=s.bestJ;
		final int openExtend=Blosum62.GAP_OPEN+Blosum62.GAP_EXTEND;
		return glocalTraceback(q, t, M, E, F, rowStride, best, bestI, bestJ, openExtend, recordPath);
	}

	/**
	 * Runs the Gotoh affine-gap local DP, filling the caller-provided flat M/E/F buffers
	 * (row-major, each row {@code rowStride} ints wide, valid region [0..m]x[0..n]).
	 * Factored out so the recurrence lives in exactly one place, shared by
	 * {@link #align(byte[],byte[])} and its path-recording overload. Writes the best cell
	 * and its coordinates into {@code s.best}/{@code s.bestI}/{@code s.bestJ} instead of
	 * allocating a result tuple.
	 *
	 * @param q Encoded query residues.
	 * @param t Encoded target residues.
	 * @param M Diagonal-state matrix (filled in place).
	 * @param E Gap-in-target matrix (filled in place).
	 * @param F Gap-in-query matrix (filled in place).
	 * @param rowStride Row width (n+1) for flat (i,j)->i*rowStride+j indexing.
	 * @param s Scratch to receive best/bestI/bestJ.
	 */
	private static final void fillMatrices(final byte[] q, final byte[] t,
			final int[] M, final int[] E, final int[] F, final int rowStride, final boolean local, final Scratch s){
		final int m=q.length, n=t.length;
		final int GO=Blosum62.GAP_OPEN, GE=Blosum62.GAP_EXTEND;
		final int openExtend=GO+GE;//cost of the first position of a gap

		//Borders. Row 0 (empty query) = 0 across all reference columns, so an alignment may begin at
		//any reference position (free reference start, both modes). Column 0 (empty reference): local
		//can start fresh at any query position (0); glocal is GLOBAL in the query, so leading query
		//with no reference is forbidden (NEG) -- the query must be placed against the (padded) reference.
		for(int j=0; j<=n; j++){M[j]=0; E[j]=NEG; F[j]=NEG;}
		for(int i=0; i<=m; i++){final int idx=i*rowStride; M[idx]=(local || i==0) ? 0 : NEG; E[idx]=NEG; F[idx]=NEG;}

		for(int i=1; i<=m; i++){
			final byte qi=q[i-1];
			final int cur=i*rowStride, prev=cur-rowStride;
			for(int j=1; j<=n; j++){
				final int sc=Blosum62.score(qi, t[j-1]);
				int diag=M[prev+j-1];
				if(E[prev+j-1]>diag){diag=E[prev+j-1];}
				if(F[prev+j-1]>diag){diag=F[prev+j-1];}
				if(local && diag<0){diag=0;}//local zero-floor: allow a fresh start (glocal has none)
				M[cur+j]=diag+sc;

				//Gap in target: open from M[i-1][j] or extend E[i-1][j].
				int e=M[prev+j]-openExtend;
				final int e2=E[prev+j]-GE;
				if(e2>e){e=e2;}
				E[cur+j]=e;

				//Gap in query: open from M[i][j-1] or extend F[i][j-1].
				int f=M[cur+j-1]-openExtend;
				final int f2=F[cur+j-1]-GE;
				if(f2>f){f=f2;}
				F[cur+j]=f;
			}
		}

		//Best cell: local = highest M anywhere. glocal (global query, free reference ends) = highest M
		//in the LAST query row (all query consumed), so the reference may end anywhere.
		int best, bestI=0, bestJ=0;
		if(local){
			best=0;
			for(int i=1; i<=m; i++){
				final int cur=i*rowStride;
				for(int j=1; j<=n; j++){if(M[cur+j]>best){best=M[cur+j]; bestI=i; bestJ=j;}}
			}
		}else{
			best=NEG; bestI=m;
			final int cur=m*rowStride;
			for(int j=1; j<=n; j++){if(M[cur+j]>best){best=M[cur+j]; bestJ=j;}}
		}
		s.best=best; s.bestI=bestI; s.bestJ=bestJ;
	}

	/**
	 * Reconstructs the best alignment by walking the DP matrices from the
	 * maximum-scoring cell back to a local start, tallying identities,
	 * mismatches, gap openings, and end coordinates.
	 */
	private static final AAAlignment traceback(final byte[] q, final byte[] t,
			final int[] M, final int[] E, final int[] F, final int rowStride,
			final int best, final int bestI, final int bestJ,
			final int openExtend, final int GE, final boolean recordPath){

		int i=bestI, j=bestJ, state=STATE_M;
		final int qStop=bestI-1, tStop=bestJ-1;//0-based inclusive ends
		int qStart=bestI-1, tStart=bestJ-1;

		int identities=0, mismatches=0, gapOpens=0, length=0;
		boolean prevQGap=false, prevTGap=false;//track gap-run boundaries (right-to-left)

		//Op path is written right-to-left into the tail of ops[], so ops[op..end] reads
		//left-to-right in alignment order; total columns cannot exceed q.length+t.length.
		final byte[] ops=(recordPath ? new byte[q.length+t.length] : null);
		int op=(ops==null ? 0 : ops.length);

		while(true){
			if(state==STATE_M){
				final byte qa=q[i-1], ta=t[j-1];
				length++;
				if(qa==ta && Blosum62.isStandard(qa)){identities++;}
				else{mismatches++;}
				prevQGap=false; prevTGap=false;
				qStart=i-1; tStart=j-1;//leftmost M column processed = alignment start
				if(ops!=null){ops[--op]='m';}

				final int sc=Blosum62.score(qa, ta);
				final int pred=M[i*rowStride+j]-sc;
				if(pred<=0){break;}//local start reached
				i--; j--;
				final int idx=i*rowStride+j;
				if(pred==M[idx]){state=STATE_M;}
				else if(pred==E[idx]){state=STATE_IX;}
				else{assert(pred==F[idx]) : "Traceback inconsistency in M"; state=STATE_IY;}
			}else if(state==STATE_IX){
				//query residue aligned to a gap in target (gap columns are neither
				//identities nor mismatches, so only length and gapOpens change)
				length++;
				if(!prevQGap){gapOpens++;}//new gap run in target (query has a residue here)
				prevQGap=true; prevTGap=false;
				if(ops!=null){ops[--op]='I';}//insertion: query residue vs a reference gap
				final boolean open=(E[i*rowStride+j]==M[(i-1)*rowStride+j]-openExtend);
				i--;
				state=open ? STATE_M : STATE_IX;
			}else{//STATE_IY: target residue aligned to a gap in query
				length++;
				if(!prevTGap){gapOpens++;}
				prevTGap=true; prevQGap=false;
				if(ops!=null){ops[--op]='D';}//deletion: reference residue vs a query gap
				final boolean open=(F[i*rowStride+j]==M[i*rowStride+j-1]-openExtend);
				j--;
				state=open ? STATE_M : STATE_IY;
			}
		}

		//Gap runs are counted by their left boundary during right-to-left walk; the
		//prevXGap flags reset at each M column, so each maximal run is counted once.
		final byte[] match=(ops==null ? null : java.util.Arrays.copyOfRange(ops, op, ops.length));
		return new AAAlignment(best, qStart, qStop, tStart, tStop,
			identities, mismatches, gapOpens, length, match);
	}

	/**
	 * Traceback for the glocal (global-query, free-reference-ends) alignment: walks from the best
	 * last-row cell until the entire query is consumed (i reaches 0). Reference outside the query's
	 * placement is free and not recorded, so every query residue appears as an 'm' or 'I' op (the
	 * op counts m+I always equal the query length). Records the path when requested.
	 */
	private static final AAAlignment glocalTraceback(final byte[] q, final byte[] t,
			final int[] M, final int[] E, final int[] F, final int rowStride,
			final int best, final int bestI, final int bestJ,
			final int openExtend, final boolean recordPath){

		int i=bestI, j=bestJ, state=STATE_M;
		final int qStop=bestI-1, tStop=bestJ-1;
		int qStart=bestI-1, tStart=bestJ-1;

		int identities=0, mismatches=0, gapOpens=0, length=0;
		boolean prevQGap=false, prevTGap=false;

		final byte[] ops=(recordPath ? new byte[q.length+t.length] : null);
		int op=(ops==null ? 0 : ops.length);

		while(i>0){
			if(state==STATE_M){
				final byte qa=q[i-1], ta=t[j-1];
				length++;
				if(qa==ta && Blosum62.isStandard(qa)){identities++;}
				else{mismatches++;}
				prevQGap=false; prevTGap=false;
				qStart=i-1; tStart=j-1;
				if(ops!=null){ops[--op]='m';}
				final int sc=Blosum62.score(qa, ta);
				final int pred=M[i*rowStride+j]-sc;
				i--; j--;
				if(i==0){break;}//entire query consumed (reached the free reference start)
				final int idx=i*rowStride+j;
				if(pred==M[idx]){state=STATE_M;}
				else if(pred==E[idx]){state=STATE_IX;}
				else{assert(pred==F[idx]) : "Glocal traceback inconsistency in M"; state=STATE_IY;}
			}else if(state==STATE_IX){//query residue vs a reference gap = insertion (consumes query)
				length++;
				if(!prevQGap){gapOpens++;}
				prevQGap=true; prevTGap=false;
				if(ops!=null){ops[--op]='I';}
				final boolean open=(E[i*rowStride+j]==M[(i-1)*rowStride+j]-openExtend);
				i--;
				state=open ? STATE_M : STATE_IX;
			}else{//STATE_IY: reference residue vs a query gap = deletion (consumes reference only)
				length++;
				if(!prevTGap){gapOpens++;}
				prevTGap=true; prevQGap=false;
				if(ops!=null){ops[--op]='D';}
				final boolean open=(F[i*rowStride+j]==M[i*rowStride+j-1]-openExtend);
				j--;
				state=open ? STATE_M : STATE_IY;
			}
		}

		//Global query: the entire query is placed, so it always begins at residue 0. (The local
		//qStart is only updated on 'm' columns, so a path ending in a leading insertion would leave
		//it stale; report 0 explicitly, which the graph relies on to index the member from residue 0.)
		final byte[] match=(ops==null ? null : java.util.Arrays.copyOfRange(ops, op, ops.length));
		return new AAAlignment(best, 0, qStop, tStart, tStop,
			identities, mismatches, gapOpens, length, match);
	}
}
