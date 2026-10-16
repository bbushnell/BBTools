package assemble;

import java.util.ArrayList;

import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;
import shared.Shared;

/**
 * Optional interior-only rescue of an immutable node's preferred route.
 * Each worker owns inference state; native accepts and hard stops never score.
 * @author Fischl
 */
final class LiveBranchGate {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads inference without changing legacy global network dispatch. */
	LiveBranchGate(final Tadpole tad_, final String path, final float cutoff_){
		this(tad_, CellNetParser.loadInferenceFromLines(ByteFile.toLines(path)), tad_.k(), tad_.minCountExtend, cutoff_, null);
	}

	/** Supplies fixed moments for deterministic policy/inference tests. */
	LiveBranchGate(final CellNet model, final int k_, final int minExtend_, final float cutoff_,
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
	 */
	boolean rescue(final int count, final int[] left, final int[] right,
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

	/** Tests exact inference independently of native walk eligibility. */
	float score(final float[] vector){return score(local.get(), vector);}

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

	/** Merges primitive worker counters only after the assembly workers have joined. */
	long[] totals(){
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
	private static final class State {
		State(final CellNet prototype){
			net=prototype.copyDenseForInference();
			net.simdFF=false;
		}
		final CellNet net;
		final float[] vector=new float[12];
		long evaluations, rescues;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Tadpole tad;
	private final CellNet prototype;
	private final int k, minExtend;
	private final float cutoff;
	private volatile SameKJoinFeatures.Moments moments;
	private long momentNanos;
	private final ArrayList<State> states=new ArrayList<State>();
	private final ThreadLocal<State> local=new ThreadLocal<State>(){
		@Override protected State initialValue(){
			final State state=new State(prototype);
			synchronized(states){states.add(state);}
			return state;
		}
	};
}
