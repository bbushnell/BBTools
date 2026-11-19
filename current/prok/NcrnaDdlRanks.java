package prok;

import cardinality.DynamicDemiLog;
import ddl.DDLComparison;
import simd.Vector;

/** Worker-owned scores for one member-profile geometry and one captured window.
 * Profiles and their keys are borrowed read-only in current consensus-model order.
 * Primitive score arrays and one native comparison record survive replacement;
 * no sorted model objects are retained. Rank and tie spans use two linear walks.
 * This is experiment-only ranking, not model verification or caller policy.
 * The driver binds the profiles' k, orientation, middle mask, seed12345 and
 * exponent4 to this geometry; both member and query modes must agree.
 * One native query sketch and four comparison counters are owned by this worker
 * and reused after their preceding scores have been consumed.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlRanks{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlRanks(final NcrnaDdlProfile[] profiles_, final int k_){
		this(profiles_, k_, true, false);
	}

	public NcrnaDdlRanks(final NcrnaDdlProfile[] profiles_, final int k_, final boolean canonical_, final boolean middle_){
		if(profiles_==null || profiles_.length<1 || k_<1 || k_>31){
			throw new IllegalArgumentException("Ranking requires current-model-order profiles and k1..31");
		}
		if(middle_ && (canonical_ || k_<3)){
			throw new IllegalArgumentException("Middle-masked query scoring requires forward words with k at least3");
		}
		canonical=canonical_;
		middle=middle_;
		profiles=profiles_.clone();
		if(profiles[0]==null){throw new IllegalArgumentException("Every current model needs a member profile");}
		buckets=profiles[0].buckets();
		for(final NcrnaDdlProfile profile:profiles){
			if(profile==null || profile.buckets()!=buckets){
				throw new IllegalArgumentException("Every model must use this geometry's independently voted profile width");
			}
		}
		k=k_;
		query=DynamicDemiLog.create(buckets, k, 12345, 0, false, false);
		raw=new NcrnaDdlQueryKeys(buckets, k);
		scores=new float[METHODS][profiles.length];
	}

	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces all three score vectors from exactly the same whole-window bytes.
	 * The input is borrowed read-only. A failed replacement invalidates old scores.
	 */
	public void score(final byte[] window){
		ready=false;
		if(window==null || DynamicDemiLog.exponentBits()!=4 || DynamicDemiLog.blacklistExists()){
			throw new IllegalArgumentException("Ranking requires captured bytes, fixed exponent4 and no DDL blacklist");
		}
		query.clear();
		if(canonical){
			query.hash(window, null);
			raw.build(window);
		}else{raw.build(window, query, false, middle);}
		for(int model=0; model<profiles.length; model++){
			final char[] reference=profiles[model].keys();
			Vector.compareDDL(query.maxArray(), reference, counts);
			comparison.set(counts, null, null, k);
			scores[SSU][model]=comparison.cachedScore;
			scores[EQUAL][model]=counts[1];
			scores[ALL_KEYS][model]=raw.equalBuckets(reference);
		}
		ready=true;
	}

	/** Returns one score without exposing mutable worker scratch. */
	public float score(final int method, final int model){
		checkReady(method);
		if(model<0 || model>=profiles.length){throw new IllegalArgumentException("Unknown current model: "+model);}
		return scores[method][model];
	}

	/** Writes {one-based rank, tie-first rank, tie-last rank, model index}.
	 * Acceptability is the independently verified, locus-specific oracle result.
	 * No acceptable model writes four -1 values; zero scores retain their full tie.
	 * Neither oracle bits nor scores are modified; output must hold four integers.
	 */
	public void best(final int method, final boolean[] acceptable, final int[] output){
		checkReady(method);
		best(scores[method], acceptable, output);
	}

	/** Ranks by descending score then ascending original model index, without sort.
	 * Public for literal independent rank/tie fixtures; all inputs are borrowed.
	 */
	public static void best(final float[] scores, final boolean[] acceptable, final int[] output){
		if(scores==null || acceptable==null || scores.length!=acceptable.length || scores.length<1 || output==null || output.length<4){
			throw new IllegalArgumentException("Ranks require aligned nonempty score/oracle arrays and four output fields");
		}
		int selected=-1;
		for(int model=0; model<scores.length; model++){
			final float score=scores[model];
			if(!Float.isFinite(score) || score<0){throw new IllegalArgumentException("Ranking score must be finite and nonnegative: model="+model);}
			if(acceptable[model] && (selected<0 || score>scores[selected])){selected=model;}
		}
		if(selected<0){
			output[0]=output[1]=output[2]=output[3]=-1;
			return;
		}
		final float target=scores[selected];
		int higher=0, equal=0, earlier=0;
		for(int model=0; model<scores.length; model++){
			if(scores[model]>target){higher++;}
			else if(scores[model]==target){
				equal++;
				if(model<selected){earlier++;}
			}
		}
		output[0]=higher+earlier+1;
		output[1]=higher+1;
		output[2]=higher+equal;
		output[3]=selected;
		assert(output[1]<=output[0] && output[0]<=output[2] && output[2]<=scores.length)
			: "The selected acceptable model's deterministic rank must lie inside its complete score tie";
	}

	private void checkReady(final int method){
		if(!ready || method<0 || method>=METHODS){throw new IllegalStateException("Score a complete query before requesting one of the three ranking methods");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final NcrnaDdlProfile[] profiles;
	private final int buckets, k;
	private final boolean canonical, middle;
	private final NcrnaDdlQueryKeys raw;
	private final DynamicDemiLog query;
	private final int[] counts=new int[4];
	private final DDLComparison comparison=new DDLComparison();
	private final float[][] scores;
	private boolean ready;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	public static final int SSU=0, EQUAL=1, ALL_KEYS=2, METHODS=3;
}
