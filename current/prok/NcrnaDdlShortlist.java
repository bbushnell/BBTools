package prok;

/** Worker-owned RF00174 ordering scratch over one immutable profile resource.
 * Counts, scores, proposal order and emitted IDs are reused for every window.
 * The borrowed emitted array has a live prefix given by size(), never its capacity.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlShortlist{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlShortlist(final NcrnaDdlProfileResource profiles, final int cap_){
		validateCap(profiles, cap_);
		unlimited=cap_==Integer.MAX_VALUE;
		selected=profiles.selected();
		cap=Math.min(cap_, profiles.models());
		dual=profiles.dual();
		counts=dual ? null : new char[profiles.models()];
		scores=new float[profiles.models()];
		order=new NcrnaDdlInterleave(profiles.models());
		scorer=profiles.newScorer();
		models=new int[cap];
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Legacy resources keep their reviewed caps; only the selected format admits all models. */
	public static void validateCap(final NcrnaDdlProfileResource profiles, final int cap){
		if(profiles==null || profiles.models()<1 || cap!=50 && cap!=384
			&& !(cap==Integer.MAX_VALUE && profiles.selected())){
			throw new IllegalArgumentException("RF00174 ordering requires cap50 or384; cap=all additionally requires selected-format profiles");
		}
	}

	/** Unlimited is explicit, rather than inferred from a small library fitting a finite cap. */
	public boolean unlimited(){return unlimited;}
	public boolean selected(){return selected;}
	public boolean gcAdjusted(){return scorer.eligibility()!=null;}

	/** Dual resources use the complete sketch-only order and never touch an index.
	 * Hybrid resources retain their native-first filtered interleave unchanged. */
	public int[] build(final byte[] bases, final TrnaKmerIndex nativeIndex){
		if(bases==null || !dual && nativeIndex==null){throw new IllegalArgumentException("Ranking needs window bytes and the hybrid alone needs a native index");}
		if(!dual){nativeIndex.countShared(bases, counts);}
		scorer.score(bases, scores);
		final boolean[] eligible=scorer.eligibility();
		final int size=eligible!=null ? order.buildFilteredSketchCaller(scores, eligible, cap)
			: dual ? order.buildSketchCaller(scores, cap) : order.buildCaller(counts, scores, 2, cap);
		for(int position=0; position<size; position++){models[position]=order.modelAt(position);}
		assert(size<=models.length)
			: shared.KillSwitch.assertDie("Only the emitted prefix may reach alignment; unused model-array capacity is stale scratch");
		return models;
	}

	/** Actual number of distinct proposals in the most recent query. */
	public int size(){return order.size();}

	/** Only this identity denotes the borrowed live-prefix array. */
	public boolean owns(final int[] proposed){return proposed==models;}

	/** Native support is -1 when deliberately unmeasured by sketch-only ordering. */
	public int nativeCount(final int model){
		if(model<0 || model>=scores.length || order.rank(model)==0){
			throw new IllegalArgumentException("Proposal provenance requires an emitted original model");
		}
		return dual ? -1 : counts[model];
	}
	public boolean sketchFirst(final int model){return order.sketchFirst(model);}
	public int rank(final int model){return order.rank(model);}
	public float score(final int model){
		if(order.rank(model)==0){throw new IllegalArgumentException("A proposal score must belong to the current emitted order");}
		return scores[model];
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int cap;
	private final boolean unlimited, selected;
	private final boolean dual;
	private final char[] counts;
	private final float[] scores;
	private final int[] models;
	private final NcrnaDdlInterleave order;
	private final NcrnaDdlProfileResource.Scorer scorer;
}
