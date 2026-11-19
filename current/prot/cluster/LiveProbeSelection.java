package prot.cluster;

import java.util.Arrays;

/** Deterministic greedy selection of independent local probes with overlapping rescues.
 * This helper selects records only: it never mutates models or the frozen pool.
 * A consumed parent invalidates every remaining operation that uses that parent.
 * Each newly claimed residual tail reduces every remaining claimant's gain by one.
 * @author Brian Bushnell, Fischl
 */
public final class LiveProbeSelection{

	/*--------------------------------------------------------------*/
	/*----------------           Selection          ----------------*/
	/*--------------------------------------------------------------*/

	/** Selects the greatest current gain strictly above the external pass margin.
	 * Equal gains prefer the smaller request ID. Inputs and their arrays remain
	 * unchanged. Gains are original measured gains minus previously claimed rescues;
	 * they do not represent a new alignment or reevaluation of the changed models.
	 * Cost is O(P log P + R log R + (N+R) log N), storage O(N+P+R), where P and R
	 * count parent and rescue postings. IDs need not be dense.
	 */
	public static Result select(final Operation[] operations, final double margin){
		return select(operations, margin, null, null);
	}

	/** Adds a transactional callback after a measured operation clears the shared margin. */
	public static Result select(final Operation[] operations,final double margin,final AcceptanceGate gate){
		require(gate!=null,"Gated measured selection requires an explicit gate");
		return select(operations,margin,null,null,gate);
	}

	/** Ranks and accepts measured marginal gain plus a fixed, signed operation offset.
	 * Rescue overlap discounts only the measured gain. Returned coverage gains remain
	 * exact integers; the offsets are objective terms, not extra rescued tails.
	 */
	public static Result selectAdjusted(final Operation[] operations,final double[] offsets,final double margin){
		return selectAdjusted(operations,offsets,margin,null);
	}

	/** Ranks by measured gain plus objective offset, then transactionally vets exact post-overlap winners. */
	public static Result selectAdjusted(final Operation[] operations,final double[] offsets,final double margin,
			final AcceptanceGate gate){
		require(operations!=null && offsets!=null && offsets.length==operations.length && Double.isFinite(margin),
			"Adjusted selection requires one finite offset per operation and a finite signed margin");
		final double[] scores=new double[operations.length],limits=new double[operations.length];
		boolean nonzero=false;
		for(int i=0;i<operations.length;i++){
			require(operations[i]!=null && Double.isFinite(offsets[i]),"Objective offsets require present operations and finite values");
			scores[i]=operations[i].measuredNetGain+offsets[i];limits[i]=margin-offsets[i];
			require(Double.isFinite(scores[i]) && Double.isFinite(limits[i]),"Adjusted scores and equivalent raw-gain thresholds must remain finite");
			nonzero|=offsets[i]!=0;
		}
		if(!nonzero){return select(operations,margin,null,null,gate);}
		final Result result=select(operations,0,scores,limits,gate);
		return new Result(result.ids,result.gains,null);
	}

	/** Ranks exact measured marginal gains with per-operation thresholds and a transactional veto.
	 * Acceptance requires measured gain strictly above the operation's signed threshold.
	 * A failed threshold or veto consumes no parents or rescues; other thresholds
	 * may still pass, so rejection never stops the mixed-phase selection early.
	 */
	public static Result selectMeasured(final Operation[] operations,final double[] thresholds,final AcceptanceGate gate){
		require(thresholds!=null && gate!=null,"Gated measured selection requires explicit thresholds and a gate");
		return select(operations,0,null,thresholds,gate);
	}

	/** Ranks predicted marginal gains, but returns exact measured marginal gains.
	 * Predictions align with operation array order, not request ID. Acceptance needs
	 * measured gain strictly above the signed threshold; predictions affect order only.
	 * Both gain types lose one per uniquely claimed rescue; prediction is always
	 * recomputed from its initial value and exact claim count, avoiding roundoff drift.
	 */
	public static Result selectPredicted(final Operation[] operations, final double[] initialPredictions,
			final double threshold){
		require(initialPredictions!=null, "Predicted selection requires explicit initial predictions");
		return select(operations, threshold, initialPredictions, null);
	}

	/** Per-operation final thresholds allow mixed phase policy without coupling it to this helper.
	 * Threshold failures discard only that operation: a lower-ranked operation may
	 * have a lower threshold. No rejection claims tails or invalidates shared parents.
	 */
	public static Result selectPredicted(final Operation[] operations, final double[] initialPredictions,
			final double[] thresholds){
		require(initialPredictions!=null && thresholds!=null, "Predicted selection requires explicit predictions and thresholds");
		return select(operations, 0, initialPredictions, thresholds);
	}

	/** Adds a transactional callback after the measured threshold passes.
	 * Rejected candidates are skipped for this pass, even if a later commit could
	 * make them admissible. Neither their tails nor parents are consumed.
	 */
	public static Result selectPredicted(final Operation[] operations, final double[] initialPredictions,
			final double[] thresholds, final AcceptanceGate gate){
		require(initialPredictions!=null && thresholds!=null && gate!=null, "Gated prediction requires explicit scores, thresholds and a gate");
		return select(operations, 0, initialPredictions, thresholds, gate);
	}

	/** Owns a private acceptance ledger. A true result commits only that ledger;
	 * false must leave it unchanged. The index refers to the original operation
	 * array. Implementations must not mutate operations or publish model state.
	 */
	@FunctionalInterface
	public interface AcceptanceGate{
		boolean tryAccept(int operationIndex, long measuredMarginalGain);
	}

	/** Shares overlap bookkeeping while keeping predicted ordering separate from exact coverage. */
	private static Result select(final Operation[] operations, final double margin,
			final double[] initialPredictions, final double[] thresholds){
		return select(operations, margin, initialPredictions, thresholds, null);
	}

	private static Result select(final Operation[] operations, final double margin,
			final double[] initialPredictions, final double[] thresholds, final AcceptanceGate gate){
		require(operations!=null && Double.isFinite(margin), "Selection needs operations and a finite external margin");
		final Operation[] input=operations.clone();
		final double[] predictions=initialPredictions==null ? null : checkedScores(initialPredictions, input.length);
		final double[] limits=thresholds==null ? null : checkedScores(thresholds, input.length);
		final long[] requestOrder=new long[input.length];
		long parents=0, rescues=0;
		for(int i=0; i<input.length; i++){
			require(input[i]!=null, "A missing operation cannot be silently omitted from selection");
			requestOrder[i]=pair(input[i].requestId, i);
			parents+=input[i].parents.length; rescues+=input[i].rescues.length;
		}
		Arrays.sort(requestOrder);
		for(int i=1; i<requestOrder.length; i++){
			require((requestOrder[i]>>>32)!=(requestOrder[i-1]>>>32), "Request IDs must be unique for stable ties and output identity");
		}
		require(parents<=Integer.MAX_VALUE-8 && rescues<=Integer.MAX_VALUE-8,
			"Parent and rescue reverse postings exceed primitive array capacity");
		final long[] parentPosts=postings(input, true, (int)parents);
		final long[] rescuePosts=postings(input, false, (int)rescues);
		// Only a group's first posting is marked, avoiding storage indexed by tail ID.
		final boolean[] claimed=new boolean[rescuePosts.length];
		final Heap heap=new Heap(input, predictions);
		final int[] accepted=new int[input.length];
		final long[] gains=new long[input.length];
		final double[] predictedGains=predictions==null ? null : new double[input.length];
		int count=0;
		while(heap.size>0){
			final int best=heap.order[0];
			final long gain=heap.gains[best];
			if(predictions==null){
				if(limits==null){
					if(!aboveMargin(gain,margin)){break;}
				}else if(!aboveMargin(gain,limits[best])){
					heap.remove(best); continue;
				}
			}else{
				final double predicted=heap.predictedGain(best), threshold=limits==null ? margin : limits[best];
				if(!aboveMargin(gain,threshold)){
					// Ranking by prediction does not imply that other measured gains or
					// per-operation thresholds fail. Remove without claiming anything.
					heap.remove(best); continue;
				}
				predictedGains[count]=predicted;
			}
			if(gate!=null && !gate.tryAccept(best,gain)){heap.remove(best); continue;}
			accepted[count]=input[best].requestId; gains[count++]=gain;
			heap.remove(best);
			for(int parent:input[best].parents){
				for(int at=lowerBound(parentPosts, parent); sameId(parentPosts, at, parent); at++){
					heap.remove((int)parentPosts[at]);
				}
			}
			for(int tail:input[best].rescues){
				final int first=lowerBound(rescuePosts, tail);
				assert(sameId(rescuePosts, first, tail)) : "Selected rescue IDs were included in the immutable reverse postings";
				if(claimed[first]){continue;}
				claimed[first]=true;
				for(int at=first; sameId(rescuePosts, at, tail); at++){
					heap.decrement((int)rescuePosts[at]);
				}
			}
		}
		return new Result(Arrays.copyOf(accepted, count), Arrays.copyOf(gains, count),
			predictedGains==null ? null : Arrays.copyOf(predictedGains, count));
	}

	/** Predictions and thresholds must be finite and bound to the exact input roster. */
	private static double[] checkedScores(final double[] scores, final int size){
		require(scores.length==size, "Prediction and threshold counts must match the operation roster");
		final double[] copy=scores.clone();
		for(double score:copy){require(Double.isFinite(score), "Predictions and thresholds must be finite; missing scores are not NaN");}
		return copy;
	}

	/** Integer comparison avoids rounding a large long gain when the margin is a double. */
	private static boolean aboveMargin(final long gain, final double margin){
		if(margin>=0x1p63){return false;}
		if(margin< -0x1p63){return true;}
		return gain>(long)Math.floor(margin);
	}

	/** Packs nonnegative IDs and operation offsets so ordinary signed sorting groups postings. */
	private static long pair(final int id, final int operation){return ((long)id<<32)|(operation&0xffffffffL);}

	/** Builds one primitive reverse index without per-posting objects or boxed keys. */
	private static long[] postings(final Operation[] input, final boolean parents, final int count){
		final long[] result=new long[count];
		int next=0;
		for(int i=0; i<input.length; i++){
			for(int id:parents ? input[i].parents : input[i].rescues){result[next++]=pair(id, i);}
		}
		assert(next==count) : "Posting allocation must match every immutable operation's exact ID count";
		Arrays.sort(result);
		return result;
	}

	/** Finds the first posting for an ID; unused integer key space costs no memory. */
	private static int lowerBound(final long[] posts, final int id){
		final long key=pair(id, 0);
		int low=0, high=posts.length;
		while(low<high){
			final int mid=low+((high-low)>>>1);
			if(posts[mid]<key){low=mid+1;}else{high=mid;}
		}
		return low;
	}
	private static boolean sameId(final long[] posts, final int at, final int id){
		return at<posts.length && (posts[at]>>>32)==id;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Immutable records      ----------------*/
	/*--------------------------------------------------------------*/

	/** A measured operation on a single frozen baseline; arrays are defensively copied. */
	public static final class Operation{
		public Operation(final int requestId_, final int[] sortedUniqueParentIds,
				final long measuredNetGain_, final int[] sortedUniqueRescuedTailOrdinals){
			require(requestId_>=0, "Request IDs must be nonnegative stable identities");
			parents=checkedCopy(sortedUniqueParentIds, true);
			rescues=checkedCopy(sortedUniqueRescuedTailOrdinals, false);
			requestId=requestId_; measuredNetGain=measuredNetGain_;
		}
		public int[] parentIds(){return parents.clone();}
		public int[] rescuedTailOrdinals(){return rescues.clone();}
		public final int requestId;
		public final long measuredNetGain;
		private final int[] parents, rescues;
	}

	/** Accepted request IDs and their actual overlap-adjusted gains, in commit order. */
	public static final class Result{
		private Result(final int[] ids_, final long[] gains_, final double[] predictions_){
			assert(ids_.length==gains_.length && (predictions_==null || predictions_.length==ids_.length)) :
				"Each accepted request needs one measured gain and, only in predicted mode, one predicted gain";
			ids=ids_; gains=gains_; predictions=predictions_;
		}
		public int size(){return ids.length;}
		public int requestId(final int at){return ids[at];}
		public long marginalGain(final int at){return gains[at];}
		public boolean hasPredictions(){return predictions!=null;}
		/** Exposes the accepted predicted score without changing the exact measured-gain accessor. */
		public double predictedMarginalGain(final int at){
			if(predictions==null){throw new IllegalStateException("Measured selection has no predicted marginal gains");}
			return predictions[at];
		}
		private final int[] ids;
		private final long[] gains;
		private final double[] predictions;
	}

	/** Validates once at construction, keeping all later selection inputs immutable. */
	private static int[] checkedCopy(final int[] ids, final boolean parents){
		require(ids!=null && (!parents || ids.length>0), "An operation needs nonnull ID lists and at least one input parent");
		final int[] copy=ids.clone();
		for(int i=0; i<copy.length; i++){
			require(copy[i]>=0 && (i==0 || copy[i]>copy[i-1]), "Parent and rescue IDs must be sorted, unique and nonnegative");
		}
		return copy;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Indexed heap         ----------------*/
	/*--------------------------------------------------------------*/

	/** Mutable selection scratch only. Each active operation appears exactly once. */
	private static final class Heap{
		Heap(final Operation[] input_, final double[] predictions_){
			input=input_; size=input.length;
			predictions=predictions_; claims=predictions==null ? null : new int[size];
			order=new int[size]; position=new int[size]; gains=new long[size];
			for(int i=0; i<size; i++){
				order[i]=position[i]=i; gains[i]=input[i].measuredNetGain;
			}
			for(int i=(size>>>1)-1; i>=0; i--){down(i);}
		}
		/** Invalidates a consumed parent or removes the selected winner in logarithmic time. */
		void remove(final int operation){
			final int at=position[operation];
			if(at<0){return;}
			position[operation]=-1;
			final int last=order[--size];
			if(at==size){return;}
			order[at]=last; position[last]=at;
			if(at>0 && better(last, order[(at-1)>>>1])){up(at);}else{down(at);}
		}
		/** A tail can discount an active operation only once; overflow is a loud invalid gain. */
		void decrement(final int operation){
			final int at=position[operation];
			if(at<0){return;}
			gains[operation]=Math.subtractExact(gains[operation], 1L);
			if(claims!=null){
				claims[operation]++;
				assert(claims[operation]<=input[operation].rescues.length) : "A unique rescued-tail posting can be discounted at most once";
			}
			down(at);
		}
		/** Accumulates exact small discounts before rounding, rather than losing each decrement separately. */
		double predictedGain(final int operation){
			assert(predictions!=null) : "Only predicted selection has a predicted ranking score";
			return predictions[operation]-claims[operation];
		}
		private boolean better(final int a, final int b){
			if(predictions!=null){
				final double left=predictedGain(a), right=predictedGain(b);
				return left>right || left==right && input[a].requestId<input[b].requestId;
			}
			return gains[a]>gains[b] || gains[a]==gains[b] && input[a].requestId<input[b].requestId;
		}
		/** Repairs a decreased key or a replacement that belongs below its new position. */
		private void down(int at){
			while(at<size/2){
				int child=at*2+1;
				if(child+1<size && better(order[child+1], order[child])){child++;}
				if(!better(order[child], order[at])){break;}
				swap(at, child); at=child;
			}
		}
		/** Repairs arbitrary deletion when the replacement outranks its new parent. */
		private void up(int at){
			while(at>0){
				final int parent=(at-1)>>>1;
				if(!better(order[at], order[parent])){break;}
				swap(at, parent); at=parent;
			}
		}
		private void swap(final int a, final int b){
			final int x=order[a], y=order[b];
			order[a]=y; order[b]=x; position[x]=b; position[y]=a;
		}
		final Operation[] input;
		final int[] order, position;
		final long[] gains;
		final double[] predictions;
		final int[] claims;
		int size;
	}
	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}
	private LiveProbeSelection(){}
}
