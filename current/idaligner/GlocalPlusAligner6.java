package idaligner;

import java.util.Arrays;
import shared.Tools;
import java.util.concurrent.atomic.AtomicLong;

import shared.Shared;
import structures.LongList;

/**
 * GlocalPlusAligner5 + traceback. Same full-matrix (unbanded) fill and bit-packed
 * scoring as GlocalPlusAligner5 -- alignStatic() below is that class's algorithm,
 * unchanged, kept here so this class is a self-contained IDAligner implementation.
 * Adds alignAndTraceStatic(), which records the same fill into a Tracer-compatible
 * LongList trace and returns the operation string via AlignmentStats.
 *
 * Chosen over adding traceback to QuantumAligner (Neptune/Noire/Brian, Aug 19 2026):
 * Quantum's value is its sparse teleportation, which a traced/dense fill would defeat
 * entirely -- there'd be no reason to use Quantum over a full-matrix aligner at that
 * point. GlocalPlusAligner5 ALREADY fills the whole matrix every row (bandStart=1,
 * bandEnd=rLen, no banding at all -- see its alignStatic below), so it is already
 * robust to a padded/shifted window start the way Quantum's teleportation was meant to
 * provide, AND its trace is trivially one contiguous block per row (Tracer's format
 * assumption) with no teleportation-gap contiguity problem to solve at all. Traceback
 * uses Tracer's sequence-aware overload: the original byte[] query0/ref0 params are
 * still in scope even though the fill itself works on Factory.encodeLong()-packed
 * long[] sequences (an earlier version of this comment claimed no byte[] was
 * available -- that was wrong, corrected Aug 19 2026 when the blind score-derivative
 * variant was found to be non-injective and disabled in Tracer.java).
 *
 * @author Neptune
 * @date August 19, 2026
 */
public class GlocalPlusAligner6 implements IDAligner{

	public static <C extends IDAligner> void main(String[] args) throws Exception {
	    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public GlocalPlusAligner6() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "Glocal+6";}
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}

	/** Overrides the interface default (which asserts !doTrace) to actually support tracing. */
	@Override
	public final float align(byte[] a, byte[] b, AlignmentStats stats){
		return alignAndTraceStatic(a, b, stats);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Core static alignment method using dynamic programming. Unchanged copy of
	 * GlocalPlusAligner5.alignStatic -- kept here so this class doesn't depend on
	 * GlocalPlusAligner5 staying binary-compatible, and so untraced callers of Glocal+6
	 * get the exact same identity/positions as Glocal+5 would.
	 */
	public static final float alignStatic(byte[] query0, byte[] ref0, int[] posVector) {
		if(posVector==null && query0.length>ref0.length) {
			byte[] temp=query0;
			query0=ref0;
			ref0=temp;
		}

		final long[] query=Factory.encodeLong(query0, (byte)15);
		final long[] ref=Factory.encodeLong(ref0, (byte)31);

		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;
		long mloops=0;
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));

		final int bandStart=1, bandEnd=rLen;

		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

		{
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}

		for(int i=1; i<=qLen; i++){
			curr[0]=i*INS;
			final long q=query[i-1];

			if(Shared.SIMD) {
				simd.SIMDAlign.alignBandVectorDel(q, ref, bandStart, bandEnd, prev, curr);
			}else {
				for(int j=bandStart; j<=bandEnd; j++){
					final long r=ref[j-1];
					final boolean isMatch=(q==r);
					final boolean hasN=((q|r)>=15);
					final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);
					final long pj1=prev[j-1], pj=prev[j];
					final long diagScore=pj1+scoreAdd;
					final long upScore=pj+INS;
					final long maxDiagUp=Math.max(diagScore, upScore);
					curr[j]=maxDiagUp;
				}
			}

			if(!Shared.SIMD) {
				long leftCell=curr[bandStart-1];
				for(int j=bandStart; j<=bandEnd; j++){
					final long maxDiagUp=curr[j];
					final long leftScore=leftCell+DEL_INCREMENT;
					leftCell=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
					curr[j]=leftCell;
				}
			}

			if(viz!=null) {viz.print(curr, bandStart, bandEnd, rLen);}
			mloops+=(bandEnd-bandStart+1);

			long[] temp=prev;
			prev=curr;
			curr=temp;
		}
		if(viz!=null) {viz.shutdown();}
		loops.addAndGet(mloops);
		return postprocess(prev, qLen, bandStart, bandEnd, posVector);
	}

	/**
	 * Aligns sequences and returns identity, recording a Tracer-compatible trace when
	 * stats.doTrace is set. Same full-matrix fill as alignStatic above (SIMD-or-scalar
	 * diag/up pass, then a scalar deletion tail-pass), plus a LongList trace recorded
	 * after each row's FINAL values (post deletion-pass) are known -- recording the
	 * preliminary diag/up-only values would store the wrong cell contents.
	 *
	 * @param query0 Query sequence
	 * @param ref0 Reference sequence
	 * @param stats Optional container for scores, counts, and backtrace
	 * @return Identity value between 0.0 and 1.0
	 */
	public static final float alignAndTraceStatic(byte[] query0, byte[] ref0, AlignmentStats stats){
		final boolean swapped, doTrace=(stats!=null && stats.doTrace);
		if(stats==null && query0.length>ref0.length) {
			byte[] temp=query0;
			query0=ref0;
			ref0=temp;
			swapped=true;
		}else{swapped=false;}

		final long[] query=Factory.encodeLong(query0, (byte)15);
		final long[] ref=Factory.encodeLong(ref0, (byte)31);

		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;

		final int bandStart=1, bandEnd=rLen;

		// Trace storage: same header format as ScrabbleAligner.alignAndTraceStatic
		// (sign | row(21) | col(21) | distToPrevHeader(21)). bandStart is always 1 here
		// (full matrix, no banding), so contiguity holds trivially every row.
		final LongList trace=(doTrace ? new LongList((qLen+1)*(rLen+2)) : null);
		int lastHeaderIdx=0;

		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

		{
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			if(doTrace){
				trace.add(0x8000000000000000L); // Row 0, Col 0, Dist 0
				lastHeaderIdx=0;
			}
			for(int j=0; j<=rLen; j++){
				prev[j]=j*mult;
				if(doTrace){trace.add(prev[j]);}
			}
		}

		for(int i=1; i<=qLen; i++){
			curr[0]=i*INS;
			final long q=query[i-1];

			if(doTrace){
				final int dist=trace.size-lastHeaderIdx;
				assert(dist<=(int)POSITION_MASK);
				long header=0x8000000000000000L | ((long)i<<42) | ((long)bandStart<<21) | dist;
				lastHeaderIdx=trace.size;
				trace.add(header);
			}

			if(Shared.SIMD) {
				simd.SIMDAlign.alignBandVectorDel(q, ref, bandStart, bandEnd, prev, curr);
			}else {
				for(int j=bandStart; j<=bandEnd; j++){
					final long r=ref[j-1];
					final boolean isMatch=(q==r);
					final boolean hasN=((q|r)>=15);
					final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);
					final long pj1=prev[j-1], pj=prev[j];
					final long diagScore=pj1+scoreAdd;
					final long upScore=pj+INS;
					final long maxDiagUp=Math.max(diagScore, upScore);
					curr[j]=maxDiagUp;
				}
			}

			{
				long leftCell=curr[bandStart-1];
				for(int j=bandStart; j<=bandEnd; j++){
					final long maxDiagUp=curr[j];
					final long leftScore=leftCell+DEL_INCREMENT;
					leftCell=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
					curr[j]=leftCell;
				}
			}

			if(doTrace) {trace.add(curr, bandStart, bandEnd+1);}

			long[] temp=prev;
			prev=curr;
			curr=temp;
		}

		long maxScore=Long.MIN_VALUE;
		int maxPos=bandEnd;
		if(GLOBAL){
			maxScore=prev[bandEnd];
		}else{
			for(int j=bandStart; j<=bandEnd; j++){
				if(prev[j]>maxScore){maxScore=prev[j]; maxPos=j;}
			}
		}
		// Tracer uses the identical 21/21-bit packed-score layout as this class (same
		// POSITION_BITS/DEL_BITS/MATCH/SUB/INS/DEL_INCREMENT constants), so its shared
		// postprocess does the matches/subs/ins/dels system-of-equations solve and fills
		// stats directly -- no need to duplicate that math here.
		float identity=Tracer.postprocess(maxScore, maxPos, qLen, rLen, null, stats);

		if(stats!=null && stats.doTrace){
			// Sequence-aware overload -- see Tracer.java's disabled blind traceback for why.
			// query0/ref0 (the original byte[] params, pre-encodeLong) are still valid and
			// unswapped here: the swap branch above requires stats==null, but doTrace
			// requires stats!=null, so on this path the swap never fires.
			final byte[] matchString=Tracer.traceback(trace, query0, ref0, qLen, maxPos, null);
			if(swapped){Tracer.invertMatchString(matchString);}//Should never happen (stats!=null forces no swap)
			stats.setFromMatchString(matchString);
		}
		return identity;
	}

	/**
	 * Post-processes alignment matrix to extract identity score and positions.
	 * Unchanged copy of GlocalPlusAligner5's private postprocess, kept for
	 * alignStatic()'s independence from GlocalPlusAligner5.
	 */
	private static final float postprocess(long[] prev, int qLen, int bandStart, int bandEnd, int[] posVector) {
		long maxScore=Long.MIN_VALUE;
		int maxPos=bandEnd;
		if(GLOBAL){
			maxScore=prev[bandEnd];
		}else{
			for(int j=bandStart; j<=bandEnd; j++){
				long score=prev[j];
				if(score>maxScore){
					maxScore=score;
					maxPos=j;
				}
			}
		}

		final int originPos=(int)(maxScore&POSITION_MASK);
		final int endPos=maxPos;
		if(posVector!=null){
			posVector[0]=originPos;
			posVector[1]=endPos-1;
		}

		final int deletions=(int)((maxScore & DEL_MASK) >> POSITION_BITS);
		final int refAlnLength=(endPos-originPos);
		final int rawScore=(int)(maxScore >> SCORE_SHIFT);

		final int insertions=Math.max(0, qLen+deletions-refAlnLength);
		final float matches=((rawScore+qLen+deletions)/2f);
		final float substitutions=Math.max(0, qLen-matches-insertions);
	    final float identity=matches/(matches+substitutions+insertions+deletions);

		if(PRINT_OPS) {
			System.err.println("originPos="+originPos);
			System.err.println("endPos="+endPos);
			System.err.println("qLen="+qLen);
			System.err.println("matches="+matches);
			System.err.println("refAlnLength="+refAlnLength);
			System.err.println("rawScore="+rawScore);
			System.err.println("deletions="+deletions);
			System.err.println("substitutions="+substitutions);
			System.err.println("insertions="+insertions);
			System.err.println("identity="+identity);
		}

		return identity;
	}

	/**
	 * Aligns query to a specified window of the reference sequence. Unchanged copy of
	 * GlocalPlusAligner5's windowed alignStatic.
	 */
	public static final float alignStatic(final byte[] query, final byte[] ref,
			final int[] posVector, int refStart, int refEnd) {
		refStart=Math.max(refStart, 0);
		refEnd=Math.min(refEnd, ref.length-1);
		final int rlen=refEnd-refStart+1;
		final byte[] region=(rlen==ref.length ? ref : Arrays.copyOfRange(ref, refStart, refEnd+1));
		if(region!=ref){Tools.toUpperCase(region);}
		final float id=alignStatic(query, region, posVector);
		if(posVector!=null) {
			posVector[0]+=refStart;
			posVector[1]+=refStart;
		}
		return id;
	}

	/*--------------------------------------------------------------*/

	private static AtomicLong loops=new AtomicLong(0);
	public long loops() {return loops.get();}
	public void setLoops(long x) {loops.set(x);}
	public static String output=null;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final int POSITION_BITS=21;
	private static final int DEL_BITS=21;
	private static final int SCORE_SHIFT=POSITION_BITS+DEL_BITS;

	private static final long POSITION_MASK=(1L << POSITION_BITS)-1;
	private static final long DEL_MASK=((1L << DEL_BITS)-1) << POSITION_BITS;
	private static final long SCORE_MASK=~(POSITION_MASK | DEL_MASK);

	private static final long MATCH=1L << SCORE_SHIFT;
	private static final long SUB=(-1L) << SCORE_SHIFT;
	private static final long INS=(-1L) << SCORE_SHIFT;
	private static final long DEL=(-1L) << SCORE_SHIFT;
	private static final long N_SCORE=0L;
	private static final long BAD=Long.MIN_VALUE/2;
	private static final long DEL_INCREMENT=(1L<<POSITION_BITS)+DEL;

	private static final boolean PRINT_OPS=false;
	public static boolean GLOBAL=false;

}
