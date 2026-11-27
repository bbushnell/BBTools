package prok;

import java.util.Arrays;

import shared.KillSwitch;

/** Worker-owned native eligibility OR a minimum raw maximum-DDL hit count.
 * Only membership changes: ranking, caps and stopping remain the caller's job.
 * Copies membership into reused primitive scratch before overwriting the borrowed
 * model-ID buffer, so admitting an earlier ID cannot clobber a later native ID.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlAdmission{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlAdmission(final int models, final int minimumHits_, final int maximumHits_){
		if(models<1 || minimumHits_<1 || maximumHits_<minimumHits_){
			throw new IllegalArgumentException("Admission needs a positive model axis and raw-hit cutoff within its profile bins: models="
				+models+", minimum="+minimumHits_+", maximum="+maximumHits_);
		}
		nativeModels=new boolean[models];
		minimumHits=minimumHits_;
		maximumHits=maximumHits_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Forms the caller's existing maximum and ascending-ID union in one pass.
	 * Both complete score vectors already exist; no query is scored again.
	 * memberScores becomes the same maximum-score vector used by final ordering.
	 * Zero-hit native members, including an existing fallback, remain members.
	 */
	public int select(final float[] memberScores, final float[] consensusScores, final int[] eligible, final int count){
		assert(memberScores!=null && memberScores.length==nativeModels.length
			&& consensusScores!=null && consensusScores.length==nativeModels.length
			&& eligible!=null && eligible.length>=nativeModels.length && count>=0 && count<=nativeModels.length)
			: KillSwitch.assertDie("OR admission requires complete current score/model scratch; a capped buffer would silently lose candidates");
		Arrays.fill(nativeModels, false);
		for(int position=0; position<count; position++){
			final int model=eligible[position];
			assert(model>=0 && model<nativeModels.length && !nativeModels[model])
				: KillSwitch.assertDie("The native member gate must supply unique original model IDs before OR admission: model="+model);
			nativeModels[model]=true;
		}
		int written=0;
		for(int model=0; model<nativeModels.length; model++){
			final float a=memberScores[model], b=consensusScores[model];
			assert(Float.isFinite(a) && Float.isFinite(b) && a>=0 && b>=0
				&& a<=maximumHits && b<=maximumHits && a==(int)a && b==(int)b)
				: KillSwitch.assertDie("OR admission uses both complete raw integer hit vectors, never normalized or imputed scores: model="+model+", member="+a+", consensus="+b);
			final float score=memberScores[model]=Math.max(a, b);
			if(nativeModels[model] || score>=minimumHits){eligible[written++]=model;}
		}
		assert(written>=count && written<=nativeModels.length)
			: KillSwitch.assertDie("An eligibility union must retain every native member and fit the original model axis");
		return written;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final boolean[] nativeModels;
	private final int minimumHits, maximumHits;
}
