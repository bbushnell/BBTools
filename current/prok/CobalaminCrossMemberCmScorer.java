package prok;

import java.util.Arrays;

/** Immutable clan resources with mutable CM/query state owned by one worker.
 * The caller and saved-data tools share the same member-only anchor selection.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminCrossMemberCmScorer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public CobalaminCrossMemberCmScorer(final NcrnaFamily[] families,
			final NcrnaBoundaryRefiner.Factory[] factories_, final NcrnaDdlProfileResource profiles_){
		require(families!=null && families.length==3 && factories_!=null && factories_.length==3
			&& profiles_!=null && profiles_.orderOnly() && profiles_.buckets()==1024 && !profiles_.consensus(),
			"Cross-member anchors require all three ordered libraries and their existing 1024-bin member-only index");
		profiles=profiles_;
		this.families=families.clone();
		for(int member=0; member<3; member++){
			require(families[member]!=null && FAMILIES[member].equals(families[member].name)
				&& families[member].library!=null && families[member].library.length>0
				&& factories_[member] instanceof CobalaminCmEndpointRefiner.Resources,
				"Each target must bind its own immutable CM/map factory in main, second, third order");
			factories[member]=factories_[member];
			thresholds[member]=((CobalaminCmEndpointRefiner.Resources)factories[member]).gatheringThreshold();
			require(Float.isFinite(thresholds[member]), "Cross-member comparisons require each CM's finite gathering threshold");
			offsets[member+1]=Math.addExact(offsets[member], families[member].library.length);
		}
		require(offsets[3]==profiles.models(), "Merged profile ordinals must conserve the three complete target libraries");
	}

	/** Explicit workers suit standalone tools; lazy thread ownership suits caller
	 * inputs whose setup may run on a different thread from candidate processing. */
	public Worker newWorker(){return new Worker(this);}
	public Worker currentWorker(){return workers.get();}
	public float ga(final int member){
		require(member>=0 && member<3, "Gathering thresholds retain the explicit main, second, third target axis");
		return thresholds[member];
	}

	/*--------------------------------------------------------------*/
	/*----------------        Worker State          ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Worker{
		private Worker(final CobalaminCrossMemberCmScorer owner_){
			owner=owner_;
			ranker=owner.profiles.newScorer();
			ranks=new float[owner.offsets[3]];
			for(int member=0; member<3; member++){
				scorers[member]=(CobalaminCmEndpointRefiner)owner.factories[member].create(owner.families[member]);
			}
		}

		/** Build one sense-oriented query and rank all targets without repeating
		 * the immutable index or allocating a query object for each candidate. */
		public void prepare(final Orf candidate, final byte[] bases){
			bindThread();
			require(candidate!=null && bases!=null && candidate.start>=0 && candidate.stop>=candidate.start
				&& candidate.stop<bases.length, "Cross-member ranking consumes the unchanged sense-oriented candidate span");
			final long begin=NcrnaCpuProfile.cpuNow();
			final int length=candidate.stop-candidate.start+1;
			if(length>=windows.length){windows=Arrays.copyOf(windows, Math.max(length+1, windows.length*2));}
			if(windows[length]==null){windows[length]=new byte[length];}
			final byte[] window=windows[length];
			for(int i=0; i<length; i++){
				final byte base=bases[candidate.start+i];
				window[i]=base>='a' && base<='z' ? (byte)(base-32) : base;
			}
			ranker.score(window, ranks);
			for(int member=0; member<3; member++){
				models[member]=bestModel(ranks, owner.offsets[member], owner.offsets[member+1])-owner.offsets[member];
			}
			prepared=candidate;
			preparedStart=candidate.start;
			preparedStop=candidate.stop;
			rankingCpu=elapsed(begin);
			rankingPending=true;
		}

		/** The first additional target bears the single shared ranking cost. */
		public float scoreChosen(final int member, final Orf candidate, final byte[] bases){
			bindThread();
			require(member>=0 && member<3 && candidate!=null && prepared==candidate && candidate.start==preparedStart
				&& candidate.stop==preparedStop, "Every additional CM must use the same prepared candidate, never stale window ranks");
			final long ranking=rankingPending ? rankingCpu : 0;
			rankingPending=false;
			return measure(member, models[member], candidate, bases, ranking);
		}

		/** Own-member measurements use the emitted model, not a newly ranked one. */
		public float scoreExact(final int member, final int model, final Orf candidate, final byte[] bases){
			bindThread();
			return measure(member, model, candidate, bases, 0);
		}

		/** A real ordinary alignment chooses both the target model and its sense
		 * location. No second DDL ranking is performed, and the emitted candidate
		 * remains the CM overlap constraint rather than a neighboring interval. */
		public float scoreAnchored(final int member, final int model, final Orf candidate, final byte[] bases,
				final int anchorStart, final int anchorStop){
			bindThread();
			require(bases!=null && anchorStart>=0 && anchorStop>=anchorStart && anchorStop<bases.length,
				"A qualifying target must supply its actual sense-oriented genomic alignment span");
			return measure(member, model, candidate, bases, 0, anchorStart, anchorStop);
		}

		private float measure(final int member, final int model, final Orf candidate, final byte[] bases, final long ranking){
			return measure(member, model, candidate, bases, ranking, -1, -1);
		}
		private float measure(final int member, final int model, final Orf candidate, final byte[] bases,
				final long ranking, final int anchorStart, final int anchorStop){
			require(member>=0 && member<3 && model>=0 && model<owner.offsets[member+1]-owner.offsets[member],
				"A local target model ordinal cannot be interpreted in another member's CM/map axis");
			require(candidate!=null && bases!=null && candidate.start>=0 && candidate.stop>=candidate.start
				&& candidate.stop<bases.length, "Target scoring requires one valid unchanged oriented contig span");
			final int start=candidate.start, stop=candidate.stop;
			final CobalaminCmEndpointRefiner scorer=scorers[member];
			lastModel=model;
			lastGa=owner.thresholds[member];
			lastRankingCpu=ranking;
			lastAnchorCpu=lastCmCpu=lastFills=0;
			lastScore=Float.NaN;
			lastStatus=UNMAPPED_TARGET_MODEL;
			lastWindowStart=lastWindowStop=-1;
			lastProposalStart=lastProposalStop=-1;
			final long begin=NcrnaCpuProfile.cpuNow();
			if(scorer.applies(model)){
				final long anchor=scorer.count(CobalaminCmEndpointRefiner.ANCHOR_CPU);
				final long cm=scorer.searchCpu(), fills=scorer.searchCalls();
				lastScore=anchorStart>=0 ? scorer.scoreOnlyAnchored(candidate, bases, model, anchorStart, anchorStop)
					: scorer.scoreOnly(candidate, bases, model);
				lastWindowStart=scorer.windowStartFor(candidate, model);
				lastWindowStop=scorer.windowStopFor(candidate, model);
				lastProposalStart=scorer.proposalStartFor(candidate, model);
				lastProposalStop=scorer.proposalStopFor(candidate, model);
				lastAnchorCpu=scorer.count(CobalaminCmEndpointRefiner.ANCHOR_CPU)-anchor;
				lastCmCpu=scorer.searchCpu()-cm;
				lastFills=scorer.searchCalls()-fills;
				lastStatus=Float.isNaN(lastScore) ? UNAVAILABLE_TARGET_ANCHOR
					: anchorStart>=0 ? QUALIFIED_ORDINARY_TARGET_ANCHOR : INDEPENDENT_TARGET_ANCHOR;
			}
			lastTotalCpu=ranking+elapsed(begin);
			require(candidate.start==start && candidate.stop==stop
				&& (Float.isNaN(lastScore) || Float.isFinite(lastScore) || lastScore==Float.NEGATIVE_INFINITY)
				&& Float.isFinite(lastGa) && lastAnchorCpu>=0 && lastCmCpu>=0 && lastFills>=0,
				"Score-only work preserves endpoints, missing/impossible states and monotonic CPU/fill counters");
			return lastScore;
		}

		/** Do not retain a previous input's candidate after its decision finishes. */
		public void clearPrepared(){
			prepared=null;
			rankingPending=false;
		}
		public int model(){return lastModel;}
		public float ga(){return lastGa;}
		public int windowStart(){return lastWindowStart;}
		public int windowStop(){return lastWindowStop;}
		/** Raw proposal from the last measured target; application guards remain
		 * the consumer's explicit responsibility, and no candidate is changed. */
		public int proposalStart(){return lastProposalStart;}
		public int proposalStop(){return lastProposalStop;}
		public long rankingCpu(){return lastRankingCpu;}
		public long anchorCpu(){return lastAnchorCpu;}
		public long cmCpu(){return lastCmCpu;}
		public long fills(){return lastFills;}
		public long totalCpu(){return lastTotalCpu;}
		public String status(){return lastStatus;}

		private void bindThread(){
			final long current=Thread.currentThread().getId();
			if(thread==0){thread=current;}
			require(thread==current, "Mutable CM/query state belongs to the candidate-processing thread, not the input setup thread");
		}
		private final CobalaminCrossMemberCmScorer owner;
		private final CobalaminCmEndpointRefiner[] scorers=new CobalaminCmEndpointRefiner[3];
		private final NcrnaDdlProfileResource.Scorer ranker;
		private final float[] ranks;
		private final int[] models=new int[3];
		private byte[][] windows=new byte[16][];
		private Orf prepared;
		private int preparedStart, preparedStop, lastModel, lastWindowStart, lastWindowStop;
		private int lastProposalStart=-1, lastProposalStop=-1;
		private float lastScore, lastGa;
		private long thread, rankingCpu, lastRankingCpu, lastAnchorCpu, lastCmCpu, lastFills, lastTotalCpu;
		private boolean rankingPending;
		private String lastStatus;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	public static int bestModel(final float[] scores, final int from, final int to){
		require(scores!=null && from>=0 && to>from && to<=scores.length,
			"Target ranking must cover one complete nonempty member partition");
		int best=from;
		for(int model=from; model<to; model++){
			require(Float.isFinite(scores[model]) && scores[model]>=0, "Member-only bucket counts must be finite and nonnegative");
			if(scores[model]>scores[best]){best=model;}
		}
		return best;
	}
	private static long elapsed(final long begin){
		if(begin<0){return 0;}
		final long value=NcrnaCpuProfile.cpuNow()-begin;
		require(value>=0, "A supported current-thread CPU clock must remain monotonic");
		return value;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final NcrnaDdlProfileResource profiles;
	private final NcrnaFamily[] families;
	private final int[] offsets=new int[4];
	private final float[] thresholds=new float[3];
	private final NcrnaBoundaryRefiner.Factory[] factories=new NcrnaBoundaryRefiner.Factory[3];
	private final ThreadLocal<Worker> workers=new ThreadLocal<Worker>(){
		@Override
		protected Worker initialValue(){return newWorker();}
	};
	private static final String[] FAMILIES={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
	public static final String UNMAPPED_TARGET_MODEL="UNMAPPED_TARGET_MODEL";
	public static final String UNAVAILABLE_TARGET_ANCHOR="UNAVAILABLE_TARGET_ANCHOR";
	public static final String INDEPENDENT_TARGET_ANCHOR="INDEPENDENT_TARGET_ANCHOR";
	public static final String QUALIFIED_ORDINARY_TARGET_ANCHOR="QUALIFIED_ORDINARY_TARGET_ANCHOR";
}
