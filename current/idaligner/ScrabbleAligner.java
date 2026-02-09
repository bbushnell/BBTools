package idaligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Shared;
import shared.Tools;
import structures.LongList;

/**
 * Aligns two sequences to return ANI using dynamic banded alignment.
 * Uses only 2 arrays and avoids traceback for memory efficiency.
 * Calculates rstart and rstop without traceback.
 * Limited to length 2Mbp with 21 position bits.
 * Center of band drifts toward highest score.
 * Band starts wide and narrows to allow glocal alignments.
 * Band dynamically widens and narrows in response to sequence identity.
 * Like Wobble but trades the ring buffer for scalars.
 *
 * @author Brian Bushnell
 * @contributor Opus
 * @date May 31, 2025
 */
public class ScrabbleAligner implements IDAligner{

	/**
	 * Program entry point that delegates to Test class to avoid redundant code.
	 * @param <C> Type parameter extending IDAligner
	 * @param args Command-line arguments
	 * @throws Exception If testing fails
	 */
	public static <C extends IDAligner> void main(String[] args) throws Exception{
		StackTraceElement[] stackTrace=Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public ScrabbleAligner(){}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name(){return "Scrabble";}

	@Override
	public final float align(byte[] a, byte[] b){return alignStatic(a, b, null);}

	@Override
	public final float align(byte[] a, byte[] b, int[] pos){return alignStatic(a, b, pos);}

	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore){return alignStatic(a, b, pos);}

	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop){return alignStatic(a, b, pos, rStart, rStop);}

	@Override
	public final float align(byte[] a, byte[] b, AlignmentStats stats){
		return alignAndTraceStatic(a, b, stats);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Determines optimal bandwidth for banded alignment by testing for high-identity
	 * indel-free alignments. Calculates bandwidth based on sequence lengths and
	 * early substitution counts.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @return Optimal bandwidth for alignment
	 */
	private static int decideBandwidth(byte[] query, byte[] ref){
		int subs=0, qLen=query.length, rLen=ref.length;
		int bandwidth=Tools.mid(7, 1+Math.max(qLen, rLen)/32, 20+(int)Math.sqrt(rLen)/8);
		for(int i=0, minlen=Math.min(qLen, rLen); i<minlen && subs<bandwidth; i++){
			subs+=(query[i]!=ref[i] ? 1 : 0);}
		return Math.min(subs+1, bandwidth);
	}

	private static int decideBandwidth(byte[] query, byte[] ref, int[] pos) {
		int qLen=query.length, rLen=ref.length;
		int bandwidth=Tools.mid(7, 1+Math.max(qLen, rLen)/32, 20+(int)Math.sqrt(rLen)/8);
		return IDAlignerStatics.decideBandwidth(query, ref, pos, bandwidth);
	}

	/**
	 * Static alignment method that performs banded dynamic programming alignment.
	 * May swap sequences to ensure query is not longer than reference if posVector is null.
	 * Uses adaptive banding that widens and narrows based on alignment quality.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning{rStart, rStop}of optimal alignment.
	 * If null, sequences may be swapped for efficiency.
	 * @return Identity(0.0-1.0)
	 */
	public static final float alignStatic(byte[] query, byte[] ref, int[] posVector){
		// Swap to ensure query is not longer than ref
		if(posVector==null && query.length>ref.length){
			byte[] temp=query;
			query=ref;
			ref=temp;
		}
		//		assert(false); //123
		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;
		long mloops=0;

		//Create a visualizer if an output file is defined
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));

		// Banding parameters
		final int bandWidth0=decideBandwidth(query, ref);
		final int maxDrift=2, maxDynamic=(bandWidth0*12)/4;

		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}

		// Initialize band limits for use outside main loop
		int bandStart=1, bandEnd=rLen-1;
		int center=0;

		// Best scoring position
		int maxPos=0;
		long maxScore=2*SUB;
		int dynamicBW=0;
		int deltaBW=0;//Bandwidth Change

		// Fill alignment matrix
		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];// Cache the query

			// Calculate bonus bandwidth due to low local alignment quality
			final boolean nextMatch=(q==ref[Math.min(rLen-1, maxPos)]);
			if(nextMatch){// Reduce bandwidth cautiously
				deltaBW=(deltaBW<0 ? Math.max(-maxDynamic, deltaBW*2) : -2);
			}else{// Increase bandwidth rapidly
				deltaBW=Tools.mid(1,(maxDynamic-dynamicBW)/2, 8);
			}
			dynamicBW=Tools.mid(0, dynamicBW+deltaBW, maxDynamic);

			// Add dynamic bandwidth from score and near the top row
			final int bandWidth=bandWidth0+Math.max(16+bandWidth0*12-maxDrift*i, dynamicBW);
			final int quarterBand=bandWidth/4;
			// Center drift for this round
			final int drift=Tools.mid(-1, maxPos-center, maxDrift);
			// New band center
			center=center+1+drift;
			bandStart=Math.max(bandStart, center-bandWidth+quarterBand);
			bandEnd=Math.min(rLen, center+bandWidth+quarterBand);

			//Clear stale data to the left of the band
			curr[bandStart-1]=BAD;

			// Clear first column score
			curr[0]=i*INS;

			//Swap row best scores
			maxScore=BAD;
			maxPos=-1;

			// Process only cells within the band
			for(int j=bandStart; j<=bandEnd; j++){
				final byte r=ref[j-1];

				// Branchless score calculation
				final boolean isMatch=(q==r && q!='N');
				final boolean hasN=(q=='N' || r=='N');
				final long scoreAdd=isMatch ? MATCH :(hasN ? N_SCORE : SUB);

				// Read adjacent scores
				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;// Match/Sub
				final long upScore=pj+INS;
				final long leftScore=cj1+DEL_INCREMENT;

				// Find max using conditional expressions
				final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
				// Changing this conditional to max or removing the mask causes a slowdown.
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;

				// Write score to current cell
				curr[j]=maxValue;

				// Track best score in row
				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;
				maxPos=better ? j : maxPos;
			}
			if(viz!=null){viz.print(curr, bandStart, bandEnd, rLen);}
			mloops+=(bandEnd-bandStart+1);

			// Swap rows
			long[] temp=prev;
			prev=curr;
			curr=temp;
		}
		if(viz!=null){viz.shutdown();}// Terminate visualizer
		if(GLOBAL){maxPos=rLen;maxScore=prev[rLen-1]+DEL_INCREMENT;}//The last cell may be empty 
		loops.addAndGet(mloops);
		return Tracer.postprocess(maxScore, maxPos, qLen, rLen, posVector, null);
	}

	/**
	 * Lightweight wrapper for aligning to a window of the reference sequence.
	 * Extracts reference region and adjusts coordinates in result.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning{rStart, rStop}of optimal alignment
	 * @param refStart Alignment window start position
	 * @param refEnd Alignment window end position
	 * @return Identity(0.0-1.0)
	 */
	public static final float alignStatic(final byte[] query, final byte[] ref, 
		final int[] posVector, int refStart, int refEnd){
		refStart=Math.max(refStart, 0);
		refEnd=Math.min(refEnd, ref.length-1);
		final int rlen=refEnd-refStart+1;
		final byte[] region=(rlen==ref.length ? ref : Arrays.copyOfRange(ref, refStart, refEnd));
		final float id=alignStatic(query, region, posVector);
		assert(posVector[1]>0) : id+", "+Arrays.toString(posVector)+", "+refStart;
		if(posVector!=null){
			posVector[0]+=refStart;
			posVector[1]+=refStart;
		}
		return id;
	}

	private static AtomicLong loops=new AtomicLong(0);
	public long loops(){return loops.get();}
	public void setLoops(long x){loops.set(x);}
	public static String output=null;

	/*--------------------------------------------------------------*/
	/*----------------          Traceback           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Aligns sequences and returns the match string(CIGAR-like but with m/S/I/D/N).
	 * Records the alignment trace using a sparse LongList to minimize memory overhead.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param stats Optional container for scores, counts, and backtrace
	 * @return Byte array containing the match string(e.g. "mmmmmSmmDmmm")
	 */
	public static final float alignAndTraceStatic(byte[] query, byte[] ref, AlignmentStats stats){
		// Swap to ensure query is not longer than ref
		final boolean swapped, doTrace=(stats!=null && stats.doTrace);
		if(stats==null && query.length>ref.length){
			byte[] temp=query;
			query=ref;
			ref=temp;
			swapped=true;
		}else{swapped=false;}
		//		assert(false); //123

		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;

		// Initialize Trace Storage
		// Header format: 1(sign) | row(21) | col(21) | distToPrevHeader(21)
		final LongList trace=(doTrace ? new LongList(qLen*20) : null); 
		int lastHeaderIdx=0; // Tracks index of the previous row's header

		final int bandWidth0=decideBandwidth(query, ref);
		final int maxDrift=2, maxDynamic=(bandWidth0*12)/4;

		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

		{// Initialize first row(Row 0)
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			if(doTrace){
				// Row 0, Col 0, Dist 0
				trace.add(0x8000000000000000L); 
				lastHeaderIdx=0;
			}
			for(int j=0; j<=rLen; j++){
				prev[j]=j*mult;
				if(doTrace){trace.add(prev[j]);}
			}
		}

		int bandStart=1, bandEnd=rLen-1;
		int center=0;
		int maxPos=0;
		long maxScore=2*SUB;
		int dynamicBW=0;
		int deltaBW=0;

		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];

			final boolean nextMatch=(q==ref[Math.min(rLen-1, maxPos)]);
			if(nextMatch){deltaBW=(deltaBW<0 ? Math.max(-maxDynamic, deltaBW*2) : -2);}
			else{deltaBW=Tools.mid(1, (maxDynamic-dynamicBW)/2, 8);}
			dynamicBW=Tools.mid(0, dynamicBW+deltaBW, maxDynamic);

			final int bandWidth=bandWidth0+Math.max(16+bandWidth0*12-maxDrift*i, dynamicBW);
			final int quarterBand=bandWidth/4;
			final int drift=Tools.mid(-1, maxPos-center, maxDrift);
			center=center+1+drift;
			bandStart=Math.max(bandStart, center-bandWidth+quarterBand);
			bandEnd=Math.min(rLen, center+bandWidth+quarterBand);

			if(doTrace){
				final int dist=trace.size-lastHeaderIdx;
				assert(dist<=(int)POSITION_MASK);
				// Header: Sign | Row | Col | Dist
				long header=0x8000000000000000L | ((long)i<<42) | ((long)bandStart<<21) | dist;
				lastHeaderIdx=trace.size;
				trace.add(header);
			}

			curr[bandStart-1]=BAD;
			curr[0]=i*INS;
			maxScore=BAD;
			maxPos=-1;

			for(int j=bandStart; j<=bandEnd; j++){
				final byte r=ref[j-1];

				final boolean isMatch=(q==r && q!='N');
				final boolean hasN=(q=='N' || r=='N');
				final long scoreAdd=isMatch ? MATCH :(hasN ? N_SCORE : SUB);

				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long maxDiagUp=Math.max(pj1+scoreAdd, pj+INS);
				final long maxValue=(maxDiagUp&SCORE_MASK)>=(cj1+DEL_INCREMENT) ? maxDiagUp : (cj1+DEL_INCREMENT);

				curr[j]=maxValue;

				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;
				maxPos=better ? j : maxPos;
			}
			if(doTrace) {trace.add(curr, bandStart, bandEnd+1);}
			long[] temp=prev; prev=curr; curr=temp;
		}

		if(GLOBAL){maxPos=rLen;}
		float identity=Tracer.postprocess(maxScore, maxPos, qLen, rLen, null, stats);
		if(stats!=null && stats.doTrace){
			final byte[] matchString=Tracer.traceback(trace, qLen, maxPos, null);
			if(swapped){Tracer.invertMatchString(matchString);}//Should never happen
			stats.setFromMatchString(matchString);
		}
		return identity;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	// Bit field definitions
	private static final int POSITION_BITS=21;
	private static final int DEL_BITS=21;
	private static final int SCORE_SHIFT=POSITION_BITS+DEL_BITS;

	// Masks
	private static final long POSITION_MASK=(1L << POSITION_BITS)-1;
	private static final long DEL_MASK=((1L << DEL_BITS)-1) << POSITION_BITS;
	private static final long SCORE_MASK=~(POSITION_MASK | DEL_MASK);

	// Scoring constants
	private static final long MATCH=1L << SCORE_SHIFT;
	private static final long SUB=(-1L) << SCORE_SHIFT;
	private static final long INS=(-1L) << SCORE_SHIFT;
	private static final long DEL=(-1L) << SCORE_SHIFT;
	private static final long N_SCORE=0L;
	private static final long BAD=Long.MIN_VALUE/2;
	private static final long DEL_INCREMENT=DEL+(1L<<POSITION_BITS);

	// Run modes
	private static final boolean PRINT_OPS=false;
	public static final boolean GLOBAL=false;

}
