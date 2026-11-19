package prot;

import java.util.Arrays;

/** Explicit opt-in stop-anchored protein scoring. Existing callers are untouched.
 * GLOCAL delegates to the historical query-global/reference-local aligner.
 * QUERY_SUFFIX consumes the entire family and a nonempty suffix of the query;
 * only the upstream query prefix is free. Internal and terminal gaps cost GAP.
 * Both encoded input arrays are borrowed read-only. Per-thread scratch is never
 * exposed; returned alignments and optional paths belong to the caller.
 * @author Keqing
 */
public final class StopAnchoredAminoLinear{

	private StopAnchoredAminoLinear(){}

	/** The no-mode API keeps the old scoring convention; suffix scoring is opt-in. */
	public static Score score(final ProteinSequence query, final ProteinSequence family){
		return score(query, family, Mode.GLOCAL, false);
	}

	/** Normalizes by the self score of the consumed query region, excluding a free
	 * prefix in QUERY_SUFFIX mode. Raw negative scores remain negative; no clipping
	 * or current family thresholds are applied. Nonpositive denominators fail. */
	public static Score score(final ProteinSequence query, final ProteinSequence family,
			final Mode mode, final boolean recordPath){
		require(query!=null && family!=null, "Scoring requires both query and family sequences");
		final AAAlignment a=align(query.enc, family.enc, mode, recordPath);
		int self=0;
		for(int i=a.qStart; i<=a.qStop; i++){self=Math.addExact(self, Blosum62.score(query.enc[i], query.enc[i]));}
		require(self>0, "Consumed query self score must be positive for a finite ratio: query="+query.id+" self="+self);
		return new Score(mode, a, self);
	}

	/** Defaults to exactly the existing glocal implementation, without a path. */
	public static AAAlignment align(final byte[] query, final byte[] family){
		return align(query, family, Mode.GLOCAL, false);
	}

	/** Encoded residues must be standard amino acids or the existing X code.
	 * Coordinates are zero-based inclusive in the original arrays. In suffix mode,
	 * family coverage is global and both right ends are consumed. "Anchored" does
	 * not require the last residues to be identical or forbid charged terminal gaps. */
	public static AAAlignment align(final byte[] query, final byte[] family,
			final Mode mode, final boolean recordPath){
		validate(query, "query"); validate(family, "family");
		require(mode!=null, "Alignment mode must be explicit or use the default overload");
		if(mode==Mode.GLOCAL){return GlocalAminoLinear.align(query, family, recordPath);}
		final int m=query.length, n=family.length, stride=Math.addExact(n, 1);
		final int cells=Math.multiplyExact(Math.addExact(m, 1), stride);
		// BLOSUM62 has maximum11; a whole path has at most m+n columns. Leave
		// headroom for recurrence arithmetic instead of silently overflowing scores.
		require(11L*(m+(long)n)<Integer.MAX_VALUE/2, "Alignment dimensions exceed safe integer-score bounds");
		final Scratch s=SCRATCH.get(); s.ensure(cells);
		fill(query, family, s.matrix, stride);
		return trace(query, family, s.matrix, stride, recordPath);
	}

	/** Same linear recurrence as GlocalAminoLinear, with reversed free-end policy:
	 * H[0,j]=-j*GAP consumes the family prefix; H[i,0]=0 skips query prefix.
	 * The endpoint is H[m,n], never the best internal column or a zero reset. */
	private static void fill(final byte[] q, final byte[] t, final int[] h, final int stride){
		assert(stride==t.length+1 && h.length>=(long)(q.length+1)*stride) : "Scratch must fit every traceback cell";
		for(int j=0; j<=t.length; j++){h[j]=-j*GAP;}
		for(int i=1; i<=q.length; i++){
			final int current=i*stride, previous=current-stride;
			final int[] scores=GlocalAminoScoreOnly.ROWS[q[i-1]];
			int diagonal=h[previous], left=0; h[current]=0;
			for(int j=1; j<=t.length; j++){
				final int above=h[previous+j];
				int best=diagonal+scores[t[j-1]];
				final int up=above-GAP, across=left-GAP;
				if(up>best){best=up;} if(across>best){best=across;}
				h[current+j]=best; diagonal=above; left=best;
			}
		}
	}

	/** Deterministic ties are diagonal, query insertion, then family insertion.
	 * Stop at family column0: remaining query prefix is outside the alignment. */
	private static AAAlignment trace(final byte[] q, final byte[] t, final int[] h,
			final int stride, final boolean recordPath){
		int i=q.length, j=t.length, identical=0, mismatch=0, gaps=0, columns=0, traceScore=0;
		final int best=h[i*stride+j];
		boolean queryGap=false, familyGap=false;
		final byte[] ops=recordPath ? new byte[Math.addExact(q.length, t.length)] : null;
		int op=ops==null ? 0 : ops.length;
		while(j>0){
			if(i>0 && h[i*stride+j]==h[(i-1)*stride+j-1]+Blosum62.score(q[i-1], t[j-1])){
				if(q[i-1]==t[j-1] && Blosum62.isStandard(q[i-1])){identical++;}else{mismatch++;}
				traceScore+=Blosum62.score(q[i-1], t[j-1]); columns++;
				queryGap=familyGap=false; if(ops!=null){ops[--op]='m';} i--; j--;
			}else if(i>0 && h[i*stride+j]==h[(i-1)*stride+j]-GAP){
				if(!familyGap){gaps++;} familyGap=true; queryGap=false;
				traceScore-=GAP; columns++; if(ops!=null){ops[--op]='I';} i--;
			}else{
				assert(h[i*stride+j]==h[i*stride+j-1]-GAP) : "Suffix traceback must reproduce a recorded predecessor at "+i+","+j;
				if(!queryGap){gaps++;} queryGap=true; familyGap=false;
				traceScore-=GAP; columns++; if(ops!=null){ops[--op]='D';} j--;
			}
		}
		// With BLOSUM62 minimum -4 and GAP=4, a terminal diagonal is never worse
		// than deleting the entire family; diagonal-first ties consume query data.
		require(i<q.length, "A scored suffix must consume query residues; check BLOSUM/gap policy if this fails");
		assert(traceScore==best && columns>=t.length && identical+mismatch<=columns) :
			"Free prefix contributes zero; consumed columns must reconstruct the endpoint score and whole family";
		return new AAAlignment(best, i, q.length-1, 0, t.length-1, identical, mismatch, gaps, columns,
			ops==null ? null : Arrays.copyOfRange(ops, op, ops.length));
	}

	private static void validate(final byte[] sequence, final String role){
		require(sequence!=null && sequence.length>0, "Nonempty encoded "+role+" required");
		for(int i=0; i<sequence.length; i++){
			final byte residue=sequence[i];
			require(Blosum62.isStandard(residue) || residue==Blosum62.X_CODE,
				"Invalid encoded "+role+" residue at "+i+": "+residue);
		}
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	public enum Mode{GLOCAL, QUERY_SUFFIX}

	/** Values retain the mode that produced them; do not mix with old calibrated
	 * cutoffs. QUERY_SUFFIX requires its own measured acceptance policy. */
	public static final class Score{
		private Score(final Mode mode_, final AAAlignment alignment_, final int self_){
			mode=mode_; alignment=alignment_; querySelfScore=self_;
			identityPct=alignment.pident(); scoreRatio=alignment.rawScore/(double)self_;
			require(Double.isFinite(identityPct) && identityPct>=0 && identityPct<=100 && Double.isFinite(scoreRatio),
				"Alignment API must return finite metrics in the stated coordinates");
		}
		public final Mode mode;
		public final AAAlignment alignment;
		public final int querySelfScore;
		public final double identityPct, scoreRatio;
	}
	private static final class Scratch{
		void ensure(final int cells){assert(cells>0) : "A nonempty DP matrix must contain cells"; if(matrix.length<cells){matrix=new int[cells];}}
		int[] matrix=new int[0];
	}
	public static final int GAP=GlocalAminoScoreOnly.GAP;
	private static final ThreadLocal<Scratch> SCRATCH=ThreadLocal.withInitial(Scratch::new);
}
