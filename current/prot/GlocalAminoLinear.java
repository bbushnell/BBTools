package prot;

/**
 * Detailed-pass glocal (query-global, reference-local) BLOSUM62 aligner with a LINEAR gap and full
 * traceback — the D55/D56 canonical detailed alignment. Reconstructs the SAME optimum the score-only
 * fast pass ({@link GlocalAminoScoreOnly}/{@link GlocalAminoSimd}) finds (identical {@code rawScore}
 * and reference end) and additionally reports the statistics D50 accepts on: percent identity
 * (identical residue columns / total alignment columns INCLUDING gap columns, matching
 * {@link AAAlignment#pident}) and the aligned spans for D55's span-based overlap
 * (min(query span / qLen, reference span / tLen)).
 *
 * <p><b>Recurrence</b> (single matrix, LINEAR gap, no gap-open — D55; GAP is the shared constant
 * {@link GlocalAminoScoreOnly#GAP}): {@code H[i][j]=max(H[i-1][j-1]+S(q_i,r_j), H[i-1][j]-GAP,
 * H[i][j-1]-GAP)}, {@code H[0][j]=0} (free reference start), {@code H[i][0]=-i*GAP} (query-global:
 * leading query is gap-charged, not free), score = first (leftmost) maximum over the LAST query row
 * (reference-local end). The BOUNDARY (query-global/reference-local) is shared with both
 * {@link GlocalAminoScoreOnly} and {@link AAAligner#alignGlocal}; the RECURRENCE is LINEAR (single
 * gap, no open) like {@link GlocalAminoScoreOnly}, NOT affine like {@link AAAligner#alignGlocal}.
 * Score and reference-end PARITY is therefore with the linear score-only path
 * ({@link GlocalAminoScoreOnly}/{@link GlocalAminoSimd}), verified exactly in
 * {@code GlocalAminoLinearTest} — NOT with the affine {@code AAAligner.alignGlocal}, whose scores
 * differ by gap model (Sayu's review, 2026-09-09).</p>
 *
 * <p><b>Traceback</b>: a single-matrix walk-back — a linear gap has no open/extend distinction, so
 * there are no M/E/F states and this is structurally simpler than {@link AAAligner}'s affine
 * {@code glocalTraceback} (Eru's observation, 2026-09-09). From the best last-row cell the walk
 * re-derives each step by which predecessor reproduces {@code H[i][j]}, with a FIXED deterministic
 * tie order DIAG &gt; up (query-residue vs reference-gap, 'I') &gt; left (reference-residue vs
 * query-gap, 'D') — matching AAAligner's M &gt; IX &gt; IY priority — so the reconstructed path, and
 * therefore identity and the spans, are deterministic across runs and implementations even when
 * several optimal paths share the score. The walk ends when the whole query is consumed (i==0);
 * query-global ⇒ {@code qStart=0}, {@code qStop=qLen-1}, and query coverage is always 1, so D55's
 * span-based overlap reduces to reference coverage {@code (tStop-tStart+1)/tLen}.</p>
 *
 * <p>Zero heap allocation at steady state on the {@link #align(byte[],byte[])} path (one reused
 * thread-local matrix, grown never shrunk); {@code recordPath} and {@link #alignReference} allocate
 * and exist for path consumers and the gate.</p>
 *
 * @author UMP45
 */
public final class GlocalAminoLinear {

	private GlocalAminoLinear(){}

	/** Shared D55 linear gap cost (per gap symbol, no gap-open); one source of truth with the fast pass. */
	public static final int GAP=GlocalAminoScoreOnly.GAP;

	/** Per-thread reused FULL DP matrix (traceback needs the whole matrix, not two rolling rows); grown never shrunk. */
	private static final class Scratch{
		int[] H=new int[0];
		int capacity=0;
		void ensureCapacity(final int needed){if(capacity<needed){H=new int[needed]; capacity=needed;}}
	}
	private static final ThreadLocal<Scratch> SCRATCH=ThreadLocal.withInitial(Scratch::new);

	/** Aligns PRE-ENCODED query (query-global) vs reference (reference-local); never null for nonempty input. */
	public static AAAlignment align(final byte[] q, final byte[] t){return align(q, t, false);}

	/**
	 * As {@link #align(byte[],byte[])}; when {@code recordPath} is true the result's
	 * {@link AAAlignment#match} holds the per-column m/D/I ops.
	 */
	public static AAAlignment align(final byte[] q, final byte[] t, final boolean recordPath){
		final int m=q.length, n=t.length;
		assert(m>0 && n>0) : "empty sequence: q="+m+" t="+n;
		final int rowStride=n+1;
		final Scratch s=SCRATCH.get();
		s.ensureCapacity((m+1)*rowStride);
		final int[] H=s.H;
		fill(q, t, H, rowStride);
		final int lastRow=m*rowStride;
		int best=H[lastRow+1], bestJ=1;//leftmost maximum over the last query row (matches GlocalAminoScoreOnly)
		for(int j=2; j<=n; j++){if(H[lastRow+j]>best){best=H[lastRow+j]; bestJ=j;}}
		return walk(q, t, H, rowStride, best, m, bestJ, recordPath);
	}

	/**
	 * Fills the flat single matrix H (row-major, rowStride=n+1) with the linear glocal recurrence.
	 * Hot-loop shape borrowed from {@link GlocalAminoScoreOnly}: the BLOSUM row for the current query
	 * residue is hoisted ({@code row=ROWS[q_i]}, one 1-D load per cell instead of a 2-D
	 * {@code Blosum62.score} call), and the diagonal predecessor and left neighbour are carried in
	 * registers rather than re-read from H — the only extra cost over the score-only scorer is that
	 * the whole matrix is written (traceback needs it), not two rolling rows.
	 */
	private static void fill(final byte[] q, final byte[] t, final int[] H, final int rowStride){
		final int m=q.length, n=t.length;
		for(int j=0; j<=n; j++){H[j]=0;}//row 0: free reference start
		for(int i=1; i<=m; i++){
			final int[] row=GlocalAminoScoreOnly.ROWS[q[i-1]];//hoisted BLOSUM row
			final int cur=i*rowStride, prev=cur-rowStride;
			int diag=H[prev];//H[i-1][0], the diagonal predecessor for j=1
			int left=H[cur]=-i*GAP;//column 0: query-global gap ramp; also the left neighbour
			for(int j=1; j<=n; j++){
				final int up=H[prev+j];
				int v=diag+row[t[j-1]];
				final int u=up-GAP; if(u>v){v=u;}
				final int l=left-GAP; if(l>v){v=l;}
				H[cur+j]=v;
				left=v; diag=up;
			}
		}
	}

	/**
	 * Walks the filled matrix from the best last-row cell to i==0, tallying identity/mismatch/gap
	 * stats and the aligned spans under the fixed tie order DIAG &gt; up('I') &gt; left('D').
	 */
	private static AAAlignment walk(final byte[] q, final byte[] t, final int[] H, final int rowStride,
			final int best, final int m, final int bestJ, final boolean recordPath){
		int i=m, j=bestJ;
		final int qStop=m-1, tStop=bestJ-1;
		int tStart=bestJ-1;//smallest reference column consumed by an 'm' or 'D'; shrinks as we walk left
		int identities=0, mismatches=0, gapOpens=0, length=0;
		boolean prevQGap=false, prevTGap=false;//gap-run boundaries, tracked right-to-left
		final byte[] ops=(recordPath ? new byte[q.length+t.length] : null);
		int op=(ops==null ? 0 : ops.length);
		while(i>0){
			if(j>0 && H[i*rowStride+j]==H[(i-1)*rowStride+(j-1)]+Blosum62.score(q[i-1], t[j-1])){
				//DIAG (highest tie priority): match/mismatch column, consumes one query and one reference residue.
				final byte qa=q[i-1], ta=t[j-1];
				length++;
				if(qa==ta && Blosum62.isStandard(qa)){identities++;}else{mismatches++;}
				prevQGap=false; prevTGap=false;
				tStart=j-1;
				if(ops!=null){ops[--op]='m';}
				i--; j--;
			}else if(H[i*rowStride+j]==H[(i-1)*rowStride+j]-GAP){
				//up: query residue vs a reference gap = insertion (consumes query only). Also the sole
				//predecessor along column 0 (H[i][0]=-i*GAP came from H[i-1][0]-GAP), so j==0 lands here.
				length++;
				if(!prevQGap){gapOpens++;}
				prevQGap=true; prevTGap=false;
				if(ops!=null){ops[--op]='I';}
				i--;
			}else{
				//left: reference residue vs a query gap = deletion (consumes reference only).
				assert(j>0 && H[i*rowStride+j]==H[i*rowStride+(j-1)]-GAP)
					: "linear glocal traceback inconsistency at i="+i+" j="+j+" H="+H[i*rowStride+j];
				length++;
				if(!prevTGap){gapOpens++;}
				prevTGap=true; prevQGap=false;
				tStart=j-1;
				if(ops!=null){ops[--op]='D';}
				j--;
			}
		}
		final byte[] match=(ops==null ? null : java.util.Arrays.copyOfRange(ops, op, ops.length));
		//Query-global: the whole query is placed, so qStart is 0 (leading 'I' columns do not move it).
		return new AAAlignment(best, 0, qStop, tStart, tStop, identities, mismatches, gapOpens, length, match);
	}

	/**
	 * Span-based overlap per D55: {@code min(query aligned span / qLen, reference aligned span / tLen)}.
	 * Under query-global geometry the query span always equals qLen (coverage 1), so this reduces to
	 * reference coverage; computed generally so the formula stays correct if the geometry ever changes.
	 */
	public static double overlap(final AAAlignment a, final int qLen, final int tLen){
		assert(qLen>0 && tLen>0) : "qLen="+qLen+" tLen="+tLen;
		final double covQ=(a.qStop-a.qStart+1)/(double)qLen;
		final double covT=(a.tStop-a.tStart+1)/(double)tLen;
		return Math.min(covQ, covT);
	}

	/**
	 * A SEPARATE allocating full-matrix reference IMPLEMENTATION (2-D {@code int[][]}, separate best-scan
	 * and walk) with the SAME recurrence and SAME tie order, for the gate. Different memory layout and
	 * index arithmetic from {@link #align}, so it catches index/scratch bugs — but it is NOT a wholly
	 * independent algorithmic oracle: it shares the recurrence and tie order by design and must agree
	 * exactly. The cross-implementation checks are the score/endpoint parity against the separate class
	 * {@link GlocalAminoScoreOnly} and the op-path structural invariants (Sayu's review, 2026-09-09).
	 */
	static AAAlignment alignReference(final byte[] q, final byte[] t, final boolean recordPath){
		final int m=q.length, n=t.length;
		if(m==0 || n==0){return null;}
		final int[][] H=new int[m+1][n+1];
		for(int j=0; j<=n; j++){H[0][j]=0;}
		for(int i=1; i<=m; i++){H[i][0]=-i*GAP;}
		for(int i=1; i<=m; i++){
			for(int j=1; j<=n; j++){
				final int d=H[i-1][j-1]+Blosum62.score(q[i-1], t[j-1]);
				final int u=H[i-1][j]-GAP;
				final int l=H[i][j-1]-GAP;
				H[i][j]=Math.max(d, Math.max(u, l));
			}
		}
		int best=H[m][1], bestJ=1;
		for(int j=2; j<=n; j++){if(H[m][j]>best){best=H[m][j]; bestJ=j;}}
		int i=m, j=bestJ;
		final int qStop=m-1, tStop=bestJ-1;
		int tStart=bestJ-1;
		int identities=0, mismatches=0, gapOpens=0, length=0;
		boolean prevQGap=false, prevTGap=false;
		final byte[] ops=(recordPath ? new byte[m+n] : null);
		int op=(ops==null ? 0 : ops.length);
		while(i>0){
			if(j>0 && H[i][j]==H[i-1][j-1]+Blosum62.score(q[i-1], t[j-1])){
				final byte qa=q[i-1], ta=t[j-1];
				length++;
				if(qa==ta && Blosum62.isStandard(qa)){identities++;}else{mismatches++;}
				prevQGap=false; prevTGap=false; tStart=j-1;
				if(ops!=null){ops[--op]='m';}
				i--; j--;
			}else if(H[i][j]==H[i-1][j]-GAP){
				length++;
				if(!prevQGap){gapOpens++;}
				prevQGap=true; prevTGap=false;
				if(ops!=null){ops[--op]='I';}
				i--;
			}else{
				length++;
				if(!prevTGap){gapOpens++;}
				prevTGap=true; prevQGap=false; tStart=j-1;
				if(ops!=null){ops[--op]='D';}
				j--;
			}
		}
		final byte[] match=(ops==null ? null : java.util.Arrays.copyOfRange(ops, op, ops.length));
		return new AAAlignment(best, 0, qStop, tStart, tStop, identities, mismatches, gapOpens, length, match);
	}
}
