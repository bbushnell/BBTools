package idaligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import structures.LongList;

/**
 * Performs sequence alignment to calculate Average Nucleotide Identity (ANI).
 * Uses a space-efficient algorithm with only 2 arrays and no traceback.
 * Calculates exact alignment scores and reference coordinates without full
 * traceback matrix reconstruction. Limited to sequences up to 2Mbp with
 * 21 position bits. Supports optional traceback via sparse LongList.
 *
 * @author Brian Bushnell
 * @contributor Isla, Neptune
 * @date April 23, 2025
 */
public class GlocalAligner implements IDAligner{

	public static <C extends IDAligner> void main(String[] args) throws Exception {
	    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public GlocalAligner() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "GlocalFull";}
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}

	@Override
	public final float align(byte[] a, byte[] b, AlignmentStats stats){
		return alignAndTraceStatic(a, b, stats);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public static final float alignStatic(byte[] query, byte[] ref, int[] posVector) {
		if(posVector==null) {return alignAndTraceStatic(query, ref, null);}
		AlignmentStats as=new AlignmentStats();
		as.doTrace=false;
		float id=alignAndTraceStatic(query, ref, as);
		posVector[0]=as.rStart;
		posVector[1]=as.rStop;
		if(posVector.length>2) {posVector[2]=as.score;}
		if(posVector.length>3) {posVector[3]=as.dels;}
		return id;
	}

	public static final float alignAndTraceStatic(byte[] query, byte[] ref, AlignmentStats stats){
		final boolean swapped, doTrace=(stats!=null && stats.doTrace);
		if(stats==null && query.length>ref.length){
			byte[] temp=query;
			query=ref;
			ref=temp;
			swapped=true;
		}else{swapped=false;}

		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;
		long mloops=0;
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));

		final LongList trace=(doTrace ? new LongList(qLen*rLen) : null);
		int lastHeaderIdx=0;

		long[] prev=new long[rLen+1], curr=new long[rLen+1];

		{
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			if(doTrace){
				trace.add(0x8000000000000000L);
				lastHeaderIdx=0;
			}
			for(int j=0; j<=rLen; j++){
				prev[j]=j*mult;
				if(doTrace){trace.add(prev[j]);}
			}
		}

		int maxPos=0;
		long maxScore=2*SUB;

		for(int i=1; i<=qLen; i++){
			curr[0]=i*INS;

			if(doTrace){
				final int dist=trace.size-lastHeaderIdx;
				assert(dist<=(int)POSITION_MASK);
				long header=0x8000000000000000L | ((long)i<<42) | ((long)1<<21) | dist;
				lastHeaderIdx=trace.size;
				trace.add(header);
			}

			final byte q=query[i-1];

			maxScore=BAD;
			maxPos=0;

			for(int j=1; j<=rLen; j++){
				final byte r=ref[j-1];

				final boolean isMatch=(q==r && q!='N');
				final boolean hasN=(q=='N' || r=='N');
				final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;
				final long upScore=pj+INS;
				final long leftScore=cj1+DEL_INCREMENT;

				final long maxDiagUp=Math.max(diagScore, upScore);
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;

				curr[j]=maxValue;

				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;
				maxPos=better ? j : maxPos;
			}
			if(viz!=null) {viz.print(curr, 1, rLen, rLen);}
			if(doTrace) {trace.add(curr, 1, rLen+1);}
			mloops+=rLen;

			long[] temp=prev;
			prev=curr;
			curr=temp;
		}
		if(viz!=null) {viz.shutdown();}
		if(GLOBAL) {maxPos=rLen;maxScore=prev[rLen-1]+DEL_INCREMENT;}
		loops.addAndGet(mloops);

		float identity=Tracer.postprocess(maxScore, maxPos, qLen, rLen, null, stats);
		if(stats!=null && stats.doTrace){
			final byte[] matchString=Tracer.traceback(trace, qLen, maxPos, null);
			if(swapped){Tracer.invertMatchString(matchString);}
			stats.setFromMatchString(matchString);
		}
		return identity;
	}

	public static final float alignStatic(final byte[] query, final byte[] ref,
			final int[] posVector, int refStart, int refEnd) {
		refStart=Math.max(refStart, 0);
		refEnd=Math.min(refEnd, ref.length-1);
		final int rlen=refEnd-refStart+1;
		final byte[] region=(rlen==ref.length ? ref : Arrays.copyOfRange(ref, refStart, refEnd));
		final float id=alignStatic(query, region, posVector);
		if(posVector!=null) {
			posVector[0]+=refStart;
			posVector[1]+=refStart;
		}
		return id;
	}

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
