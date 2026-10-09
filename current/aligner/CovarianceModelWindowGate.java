package aligner;

import java.util.Arrays;

/** Exact CYK over all nonempty target intervals in a supplied window.
 * Global model by default; explicit local mode still uses CYK, not Inside. GA remains a
 * provisional per-model threshold. The caller's coordinates are not mutated.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelWindowGate {
	public CovarianceModelWindowGate(CovarianceModel model_, float threshold_, long maxCells){
		this(model_, threshold_, maxCells, false);
	}
	public CovarianceModelWindowGate(CovarianceModel model_, float threshold_, long maxCells, boolean local_){
		require(model_!=null && Float.isFinite(threshold_) && maxCells>0,
			"Window verification requires a parsed model, finite bit cutoff and positive matrix-cell bound");
		model=model_;threshold=threshold_;local=local_;scorer=new CovarianceModelCyk(model, maxCells, local);
	}
	/** Inclusive coordinates on an already oriented sequence. Search is bounded
	 * by the entire supplied window, never silently capped at the model's W. */
	public Decision evaluate(byte[] bases, int start, int stop){
		return evaluate(bases, start, stop, -1, -1, false);
	}
	/** Candidate association uses the same original oriented-array coordinates. */
	public Decision evaluate(byte[] bases, int start, int stop, int candidateStart, int candidateStop){
		require(bases!=null && candidateStart>=0 && candidateStart<=candidateStop && candidateStop<bases.length,
			"Candidate bounds must fit the source array before window-coordinate conversion");
		return evaluate(bases, start, stop, candidateStart, candidateStop, true);
	}
	private Decision evaluate(byte[] bases, int start, int stop, int candidateStart, int candidateStop, boolean associate){
		return evaluate(bases,start,stop,candidateStart,candidateStop,associate,false);
	}
	/** Explicit experiment: retain real extracted hits instead of selecting an
	 * arbitrary root interval merely because its bounding span touches a call. */
	public Decision evaluateGreedy(byte[] bases,int start,int stop,int candidateStart,int candidateStop){
		require(bases!=null&&candidateStart>=0&&candidateStop>=candidateStart&&candidateStop<bases.length,"Greedy candidate bounds must fit the original oriented sequence");
		return evaluate(bases,start,stop,candidateStart,candidateStop,true,true);
	}
	private Decision evaluate(byte[] bases, int start, int stop, int candidateStart, int candidateStop, boolean associate, boolean greedy){
		require(bases!=null && start>=0 && start<=stop && stop<bases.length,
			"CM window must be a nonempty existing interval: "+start+".."+stop);
		final long begin=System.nanoTime();
		final byte[] sequence=Arrays.copyOfRange(bases, start, stop+1);
		final CovarianceModelCyk.Result result=greedy ? scorer.searchGreedyOverlapping(sequence,sequence.length,candidateStart-start+1,candidateStop-start+1,threshold) : associate
			? scorer.searchOverlapping(sequence, sequence.length, candidateStart-start+1, candidateStop-start+1)
			: scorer.search(sequence, sequence.length);
		final double seconds=(System.nanoTime()-begin)*1e-9;
		require(Float.isFinite(result.score) || result.score==Float.NEGATIVE_INFINITY,
			"Invalid arithmetic cannot be interpreted as a biological CM rejection");
		final boolean finite=Float.isFinite(result.score);
		final int from=finite ? start+result.targetFrom-1 : -1, to=finite ? start+result.targetTo-1 : -1;
		require(!finite || from>=start && from<=to && to<=stop,
			"Best-interval coordinates must remain within the supplied caller window");
		return new Decision(result.score, finite && result.score>=threshold, sequence.length, from, to,
			result.rootStarts.size, seconds, result.cells, result.unrestrictedScore,
			Float.isFinite(result.unrestrictedScore) ? start+result.unrestrictedFrom-1 : -1,
			Float.isFinite(result.unrestrictedScore) ? start+result.unrestrictedTo-1 : -1, local,result.greedyHits,start-1);
	}
	public static final class Decision {
		Decision(float bits_, boolean accepted_, int length_, int from_, int to_, int ties_, double seconds_, long cells_,
				float unrestrictedBits_, int unrestrictedFrom_, int unrestrictedTo_, boolean local_){
			this(bits_,accepted_,length_,from_,to_,ties_,seconds_,cells_,unrestrictedBits_,unrestrictedFrom_,unrestrictedTo_,local_,null,0);
		}
		Decision(float bits_, boolean accepted_, int length_, int from_, int to_, int ties_, double seconds_, long cells_,
				float unrestrictedBits_, int unrestrictedFrom_, int unrestrictedTo_, boolean local_,CovarianceModelCyk.GreedyHits hits_,int hitOffset_){
			assert(length_>0 && seconds_>=0 && cells_>0):"A window decision describes a completed nonempty M3 search";
			bits=bits_;accepted=accepted_;windowLength=length_;from=from_;to=to_;ties=ties_;seconds=seconds_;cells=cells_;
			unrestrictedBits=unrestrictedBits_;unrestrictedFrom=unrestrictedFrom_;unrestrictedTo=unrestrictedTo_;
			// CovarianceModelCyk.score allocates one extra int per root cell to
			// retain local-entry choices; it is separate from per-state cells.
			localBeginCells=local_ ? Math.multiplyExact((long)length_+1, (long)length_+1) : 0;
			greedyHits=hits_;greedyOffset=hitOffset_;
		}
		public final float bits;
		public final boolean accepted;
		/** Zero-based inclusive on the supplied oriented array; -1/-1 if impossible. */
		public final int from, to, windowLength, ties;
		/** Best anywhere in the window; kept separately from candidate-associated result. */
		public final float unrestrictedBits;
		public final int unrestrictedFrom, unrestrictedTo;
		/** Greedy hit positions plus greedyOffset give original zero-based bounds. */
		public final CovarianceModelCyk.GreedyHits greedyHits;
		public final int greedyOffset;
		public final double seconds;
		/** Exact CYK allocates one float score and one int traceback choice per cell;
		 * excludes array headers and result metadata, and is not JVM heap/RSS. */
		public final long cells;
		public final long localBeginCells;
		public long matrixBytes(){return Math.addExact(Math.multiplyExact(cells, 8L), Math.multiplyExact(localBeginCells, 4L));}
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	public final CovarianceModel model;
	public final float threshold;
	public final boolean local;
	private final CovarianceModelCyk scorer;
}
