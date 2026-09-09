package prot;

/**
 * Score-only glocal (query-global, reference-local) protein aligner: BLOSUM62 substitution, linear
 * gap, two rolling {@code int[]} rows, no traceback, no packed position bits. Built for the
 * shortlist scoring pass of the family assigner (Brian, 2026-09-03: "1 matrix instead of 3 and 2
 * arrays instead of a full matrix"): score every shortlisted family cheaply, then trace only the
 * winner with the full aligner.
 *
 * <p>Why this is faster than {@link GlocalAminoBlosumFast#scoreOnly}: the BLOSUM row for the
 * current query residue is hoisted out of the inner loop ({@code row=ROWS[q]}, one 1-D lookup per
 * cell instead of a 2-D {@code MATRIX[a][b]} through a method call), the left neighbour and the
 * diagonal predecessor are carried in registers instead of re-read from {@code curr}/{@code prev},
 * the max is branch-free ({@code Math.max}), and the last-row scan happens once after the loop
 * instead of a {@code lastRow &&} test per cell.
 *
 * <p>Boundary semantics (glocal): row 0 is free (reference may start anywhere: {@code prev[j]=0});
 * column 0 charges a gap per consumed query residue ({@code -i*GAP}: the query may begin before
 * the reference starts, at gap cost, exactly as idaligner.GlocalAligner's {@code curr[0]=i*INS});
 * the score is the maximum over the LAST row (the reference may end anywhere after the query
 * ends). {@code X} scores 0 against everything via the BLOSUM table (Blosum62.MATRIX row for X).
 *
 * @author UMP45
 */
public final class GlocalAminoScoreOnly {

	private GlocalAminoScoreOnly(){}

	/** Linear per-residue gap cost; same value as GlocalAminoBlosum/Fast so scores are comparable. */
	public static final int GAP=Blosum62.GAP_OPEN;

	/** Alphabet size of Blosum62.encode's output (0-19 standard + X_CODE). */
	private static final int ALPHA=Blosum62.X_CODE+1;

	/** ROWS[a][b] == Blosum62.score(a,b); a flat copy so the inner loop does one 1-D load. */
	private static final int[][] ROWS=buildRows();

	/** Codes Blosum62.encode can emit: the 20 standard residues and X_CODE (21). Code 20 is never emitted and
	 *  Blosum62.score asserts on unmapped pairs, so only emitted codes are copied; the rest stay at a
	 *  poison value that would make any accidental use obvious (score far below the all-gap floor). */
	private static int[][] buildRows(){
		final int[][] r=new int[ALPHA][ALPHA];
		final int poison=Integer.MIN_VALUE/8;
		for(int a=0; a<ALPHA; a++){java.util.Arrays.fill(r[a], poison);}
		final int[] codes=new int[21];
		for(int i=0; i<20; i++){codes[i]=i;}
		codes[20]=Blosum62.X_CODE;
		for(int a : codes){for(int b : codes){r[a][b]=Blosum62.score((byte)a, (byte)b);}}
		return r;
	}

	/** Two rolling rows per thread, grown never shrunk (BBTools zero-allocation idiom). */
	private static final class Scratch{
		int[] prev=new int[0], curr=new int[0];
		void ensureCapacity(final int needed){
			if(prev.length<needed){prev=new int[needed]; curr=new int[needed];}
		}
	}
	private static final ThreadLocal<Scratch> SCRATCH=ThreadLocal.withInitial(Scratch::new);

	/**
	 * Scores a glocal alignment of PRE-ENCODED query against PRE-ENCODED reference.
	 * Zero allocation at steady state.
	 * @param qEnc Encoded query (Blosum62.encode), length >= 1.
	 * @param rEnc Encoded reference, length >= 1.
	 * @param out Receives {score, refEndPos (0-based inclusive, of the best last-row cell)}; length >= 2.
	 */
	public static final void scoreOnly(final byte[] qEnc, final byte[] rEnc, final int[] out){
		final int qLen=qEnc.length, rLen=rEnc.length;
		assert(qLen>0 && rLen>0) : "empty sequence: q="+qLen+" r="+rLen;
		assert(out.length>=2) : "out needs 2 slots";
		final Scratch s=SCRATCH.get();
		s.ensureCapacity(rLen+1);
		int[] prev=s.prev, curr=s.curr;
		for(int j=0; j<=rLen; j++){prev[j]=0;}//free reference start

		for(int i=1; i<=qLen; i++){
			final int[] row=ROWS[qEnc[i-1]];//hoisted: one 1-D lookup per cell below
			int left=-i*GAP;//column 0: i query residues consumed before the reference starts
			curr[0]=left;
			int pd=prev[0];//diagonal predecessor prev[j-1], carried in a register
			for(int j=1; j<=rLen; j++){
				final int pu=prev[j];
				final int d=pd+row[rEnc[j-1]];
				final int v=Math.max(Math.max(d, pu-GAP), left-GAP);
				curr[j]=v;
				left=v;
				pd=pu;
			}
			final int[] t=prev; prev=curr; curr=t;
		}
		//best over the last row (reference may end anywhere after the query ends)
		int best=prev[1], bestJ=0;
		for(int j=2; j<=rLen; j++){final int v=prev[j]; if(v>best){best=v; bestJ=j-1;}}
		assert(best>=-(qLen)*GAP-(rLen)*GAP) : "score below the all-gap floor: "+best;
		out[0]=best;
		out[1]=bestJ;
	}

	/**
	 * Reference implementation with a full matrix, same recurrence, for tests only (allocates).
	 * @return {score, refEndPos}.
	 */
	static int[] scoreOnlyReference(final byte[] q, final byte[] r){
		final int[][] m=new int[q.length+1][r.length+1];
		for(int j=0; j<=r.length; j++){m[0][j]=0;}
		for(int i=1; i<=q.length; i++){
			m[i][0]=-i*GAP;
			for(int j=1; j<=r.length; j++){
				final int d=m[i-1][j-1]+Blosum62.score(q[i-1], r[j-1]);
				final int u=m[i-1][j]-GAP;
				final int l=m[i][j-1]-GAP;
				m[i][j]=Math.max(d, Math.max(u, l));
			}
		}
		int best=m[q.length][1], bestJ=0;
		for(int j=2; j<=r.length; j++){if(m[q.length][j]>best){best=m[q.length][j]; bestJ=j-1;}}
		return new int[]{best, bestJ};
	}
}
