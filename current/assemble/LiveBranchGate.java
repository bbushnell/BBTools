package assemble;

import java.util.ArrayList;

import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;
import shared.Shared;

/**
 * Optional interior-only rescue of an immutable node's preferred route.
 * Each worker owns inference state; native accepts and hard stops never score.
 * The prototype, borrowed count table and supplied population moments must stay
 * stable while workers use this gate. Keep Shared.SIMD_FMA disabled throughout
 * scoring, and join workers before reading aggregate counters or scan timing.
 * @author Fischl
 */
public final class LiveBranchGate{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads inference without changing legacy global network dispatch. */
	LiveBranchGate(final Tadpole tad_, final String path, final float cutoff_){
		this(tad_, CellNetParser.loadInferenceFromLines(ByteFile.toLines(path)), tad_.k(), tad_.minCountExtend, cutoff_, null);
	}

	/** Supplies fixed moments without scanning an assembler's count table.
	 * Borrows the model and moments unchanged for this gate's lifetime; each scoring
	 * thread lazily copies the model into private dense inference state. The caller
	 * keeps Shared.SIMD_FMA disabled until all scoring threads have finished.
	 * @param model Nonnull 12-input, one-output model supporting dense inference copies
	 * @param k_ Word length in bases, from 24 through 127
	 * @param minExtend_ Positive minimum extension count
	 * @param cutoff_ Finite inclusive acceptance threshold in [0,1]
	 * @param moments_ Nonnull completed population moments, never mutated during use
	 */
	public LiveBranchGate(final CellNet model, final int k_, final int minExtend_, final float cutoff_,
			final SameKJoinFeatures.Moments moments_){
		this(null, model, k_, minExtend_, cutoff_, moments_);
	}

	/** Validates the fixed schema and scalar experimental scoring contract once. */
	private LiveBranchGate(final Tadpole tad_, final CellNet model, final int k_, final int minExtend_,
			final float cutoff_, final SameKJoinFeatures.Moments moments_){
		if(model==null || model.numInputs()!=12 || model.numOutputs()!=1 || k_<24 || k_>127 || minExtend_<1 ||
				!Float.isFinite(cutoff_) || cutoff_<0 || cutoff_>1 || (tad_==null && moments_==null) || Shared.SIMD_FMA){
			throw new IllegalArgumentException("Live branch gate requires12 inputs, one output, K24-127, mce>=1, cutoff in[0,1], moments and simdfma=f.");
		}
		tad=tad_; prototype=model; k=k_; minExtend=minExtend_; cutoff=cutoff_; moments=moments_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * The caller excludes hidden branches before this method. Requiring the actual
	 * predecessor to be the unique preferred one also protects callers with their
	 * native hidden-branch flag disabled. This does not override ownership/cycles.
	 * Borrows both count arrays for this call without modifying them; callers keep
	 * them stable during the call. Scoring state and counters belong to each thread.
	 * @param count Nonnegative depth of the current node
	 * @param left Four nonnegative predecessor counts in A/C/G/T order
	 * @param right Four nonnegative successor counts in A/C/G/T order
	 * @param leftBranch Native branch flag for the predecessor side
	 * @param rightBranch Native branch flag for the successor side
	 * @param arrivalBase Actual predecessor base index, from 0 through 3
	 * @return True only for an eligible route whose score reaches the inclusive cutoff;
	 * native accepts, hard stops and other predecessors return false without scoring
	 */
	public boolean rescue(final int count, final int[] left, final int[] right,
			final boolean leftBranch, final boolean rightBranch, final int arrivalBase){
		if(arrivalBase<0 || arrivalBase>3){throw new IllegalArgumentException("Interior rescue requires its actual preceding base.");}
		if(LiveBranchPolicy.level(count, left, right, leftBranch, rightBranch, minExtend, 0.5)!=LiveBranchPolicy.GREY_ZONE ||
				LiveBranchPolicy.preferredBase(left)!=arrivalBase){return false;}
		final State state=local.get();
		final SameKJoinFeatures.Moments population=moments();
		LiveBranchPolicy.fill(state.vector, k, count, left, right, population.uniqueMean(), population.observationMean());
		final float score=score(state, state.vector);
		state.evaluations++;
		final boolean accepted=score>=cutoff;
		if(accepted){state.rescues++;}
		return accepted;
	}

	/** Scores with this thread's private inference state, without eligibility checks
	 * or changes to rescue/evaluation counters. Keep Shared.SIMD_FMA disabled.
	 * @param vector Twelve finite nonnegative features in LiveBranchPolicy.fill order,
	 * held stable by the caller during this call
	 * @return A finite score in [0,1]; an invalid model result throws
	 */
	public float score(final float[] vector){return score(local.get(), vector);}

	/** Immutable count moments are computed once, without a second table or read pass. */
	private SameKJoinFeatures.Moments moments(){
		SameKJoinFeatures.Moments result=moments;
		if(result==null){
			synchronized(this){
				result=moments;
				if(result==null){
					assert(tad!=null) : "Lazy moments borrow the already loaded, immutable contig count table.";
					final long start=System.nanoTime();
					result=SameKJoinFeatures.moments(tad);
					momentNanos=System.nanoTime()-start;
					moments=result;
				}
			}
		}
		return result;
	}

	/** Reports only the one-time scan, after workers join; it adds no per-node clock calls. */
	long momentNanos(){return momentNanos;}

	/** Reuses the immutable population after workers join; the postpass must not reload the table. */
	SameKJoinFeatures.Moments populationMoments(){return moments();}

	/** Merges worker counters after the caller has joined every scoring worker.
	 * The registry lock does not make concurrent counter writes safe to observe.
	 * @return A new two-element array: completed eligible evaluations, then rescues;
	 * direct score calls are excluded, and counters are not reset
	 */
	public long[] totals(){
		long evaluations=0, rescues=0;
		synchronized(states){
			for(State state : states){evaluations+=state.evaluations; rescues+=state.rescues;}
		}
		assert(rescues<=evaluations) : "Only evaluated grey routes can be rescued; callers must join workers before reporting.";
		return new long[]{evaluations, rescues};
	}

	/** Uses the same scalar dense arithmetic as the frozen validation curves. */
	private static float score(final State state, final float[] vector){
		assert(!Shared.SIMD_FMA) : "The initial cutoff was validated with simdfma=f; do not change it during worker inference.";
		if(vector==null || vector.length!=12){throw new IllegalArgumentException("Expected12 live branch features.");}
		for(float value : vector){if(!Float.isFinite(value) || value<0){throw new IllegalArgumentException("Invalid live branch feature.");}}
		final float score=state.net.applyInput(vector).feedForward();
		if(!Float.isFinite(score) || score<0 || score>1){throw new IllegalStateException("Live branch score outside[0,1]: "+score);}
		return score;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Worker-private model arrays, reusable feature buffer and unsynchronized counters. */
	private static final class State{
		/** Creates this worker's dense inference copy with private arrays and normalization. */
		State(final CellNet prototype){
			net=prototype.copyDenseForInference();
			net.simdFF=false;
		}
		/** Private inference copy; per-network batched SIMD is disabled. */
		final CellNet net;
		/** Reused live-branch feature vector, distinct from the network's input scratch. */
		final float[] vector=new float[12];
		/** This worker's scored and accepted eligible routes; merged only after joining. */
		long evaluations, rescues;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed assembler for lazy population scanning; null with supplied test moments. */
	private final Tadpole tad;
	/** Stable source for per-worker inference copies; never scored directly here. */
	private final CellNet prototype;
	/** Word length and minimum extension count used by the fixed eligibility policy. */
	private final int k, minExtend;
	/** Inclusive model-score acceptance boundary. */
	private final float cutoff;
	/** Published once after lazy accumulation, or borrowed unchanged from the constructor. */
	private volatile SameKJoinFeatures.Moments moments;
	/** Lazy scan duration in nanoseconds; zero when moments were supplied. */
	private long momentNanos;
	/** Synchronized worker registry; counter reads additionally require completed workers. */
	private final ArrayList<State> states=new ArrayList<State>();
	/** Allocates one private model/vector/counter set on each worker's first scoring use. */
	private final ThreadLocal<State> local=new ThreadLocal<State>(){
		@Override protected State initialValue(){
			final State state=new State(prototype);
			synchronized(states){states.add(state);}
			return state;
		}
	};
}
