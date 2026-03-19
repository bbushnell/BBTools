package idaligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import parse.PreParser;
import shared.Tools;
import structures.LongList;

/**
 * Banded sequence aligner for computing average nucleotide identity (ANI).
 * Uses dynamic programming with a restricted band around the diagonal to align sequences.
 * Optimized for high-identity alignments with minimal indels, avoiding full traceback.
 * Limited to sequences up to 2Mbp using 21-bit position encoding.
 * Supports SIMD prealignment for optimal band centering and optional traceback.
 *
 * @author Brian Bushnell
 * @contributor Isla, Neptune
 * @date April 24, 2025
 */
public class BandedAligner implements IDAligner{

	public static <C extends IDAligner> void main(String[] args) throws Exception {
		args=new PreParser(args, System.err, null, false, true, false).args;
	    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public BandedAligner() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "Banded";}
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

	private static int decideBandwidth(byte[] query, byte[] ref, int[] pos) {
		int qLen=query.length, rLen=ref.length;
		int bandwidth=Math.min(60+(int)Math.sqrt(rLen), 4+Math.max(qLen, rLen)/8);
		return IDAlignerStatics.decideBandwidth(query, ref, pos, bandwidth);
	}

	public static final float alignStatic(byte[] query, byte[] ref, int[] posVector) {
		// Swap to ensure query is not longer than ref when no posVector needed
		if(posVector==null && query.length>ref.length) {
			byte[] temp=query;
			query=ref;
			ref=temp;
		}
		AlignmentStats as=(posVector==null ? null : new AlignmentStats());
		if(as!=null){as.doTrace=false;}
		float id=alignAndTraceStatic(query, ref, as);
		if(posVector!=null){
			posVector[0]=as.rStart;
			posVector[1]=as.rStop;
			if(posVector.length>2) {posVector[2]=as.score;}
			if(posVector.length>3) {posVector[3]=as.dels;}
		}
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

		final LongList trace=(doTrace ? new LongList(qLen*20) : null);
		int lastHeaderIdx=0;

		// Banding parameters — SIMD prealignment finds optimal offset
		int[] bwPos=new int[2];
		if(stats!=null){bwPos[0]=stats.rStart;}
		final int bandWidth=decideBandwidth(query, ref, bwPos);
		final int offset=bwPos[0];

		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

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

		int bandStart=1, bandEnd=rLen-1;

		for(int i=1; i<=qLen; i++){
			final int center=i+offset;
			bandStart=Math.max(1, Math.min(rLen, center-bandWidth));
			bandEnd=Math.min(rLen, center+bandWidth);

			if(doTrace){
				final int dist=trace.size-lastHeaderIdx;
				assert(dist<=(int)POSITION_MASK);
				long header=0x8000000000000000L | ((long)i<<42) | ((long)bandStart<<21) | dist;
				lastHeaderIdx=trace.size;
				trace.add(header);
			}

			if(bandStart-1>=0) {curr[bandStart-1]=BAD;}
			curr[0]=i*INS;

			final byte q=query[i-1];

			maxScore=BAD;
			maxPos=-1;

			for(int j=bandStart; j<=bandEnd; j++){
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
			if(viz!=null) {viz.print(curr, bandStart, bandEnd, rLen);}
			if(doTrace) {trace.add(curr, bandStart, bandEnd+1);}
			mloops+=(bandEnd-bandStart+1);

			long[] temp=prev;
			prev=curr;
			curr=temp;
		}
		if(viz!=null) {viz.shutdown();}
		if(GLOBAL){maxPos=rLen;maxScore=prev[rLen-1]+DEL_INCREMENT;}
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
		assert(posVector[1]>0) : id+", "+Arrays.toString(posVector)+", "+refStart;
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
	private static final long DEL_INCREMENT=DEL+(1L<<POSITION_BITS);

	private static final boolean PRINT_OPS=false;
	public static final boolean GLOBAL=false;

}
