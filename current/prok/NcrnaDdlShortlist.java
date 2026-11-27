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
		this.profiles=profiles;
		unlimited=cap_==Integer.MAX_VALUE;
		selected=profiles.selected();
		orderOnly=profiles.orderOnly();
		cap=Math.min(cap_, profiles.models());
		dual=profiles.dual();
		counts=dual && !orderOnly ? null : new char[profiles.models()];
		nativeEligible=orderOnly ? new int[profiles.models()] : null;
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
		if(profiles!=null && profiles.orderOnly()){
			if(cap<1){throw new IllegalArgumentException("Ordering-only profiles preserve the caller's positive cap");}
			return;
		}
		if(profiles==null || profiles.models()<1 || cap!=50 && cap!=384
			&& !(cap==Integer.MAX_VALUE && profiles.selected())){
			throw new IllegalArgumentException("RF00174 ordering requires cap50 or384; cap=all additionally requires selected-format profiles");
		}
	}

	/** Unlimited is explicit, rather than inferred from a small library fitting a finite cap. */
	public boolean unlimited(){return unlimited;}
	public boolean selected(){return selected;}
	public boolean orderOnly(){return orderOnly;}
	/** Immutable resource binding for cold compatibility checks. */
	public NcrnaDdlProfileResource profiles(){return profiles;}
	/** The owning caller installs its existing observer; null is the default. */
	void setModelAttemptSink(final NcrnaModelAttemptInstrumentSink sink){modelAttemptSink=sink;}
	/** Both resources were loaded against the same caller model names/sequences.
	 * This check keeps the primary member-only and the rejected-window order bonus-only. */
	public static void validateBonusRescue(final NcrnaDdlProfileResource primary,
			final NcrnaDdlProfileResource rescue){
		if(primary==null || rescue==null || primary.models()!=rescue.models()
			|| !primary.orderOnly() || primary.orderingBonus()!=0 || primary.buckets()!=1024
			|| !rescue.orderOnly() || rescue.orderingBonus()!=100
			|| rescue.orderingBonusRule()!=NcrnaDdlProfileResource.BONUS_MAXIMA
			|| !primary.sameOrderingQuery(rescue)){
			throw new IllegalArgumentException("Bonus-only rescue needs member-only primary profiles and maxima100/1024 rescue profiles on the same model/query axis and blacklist");
		}
	}
	public boolean gcAdjusted(){return scorer.eligibility()!=null;}
	/** One validated query axis, either ranked independently or by maximum hits.
	 * Configuration is cold; query buffers and both model orders are reused. */
	public void configureConsensusPair(final NcrnaDdlProfileResource.ConsensusPair pair, final boolean merged){
		if(pair==null || pair.memberOnly!=profiles || !orderOnly || ties!=null
			|| memberPrefix!=null || consensusPair!=null){
			throw new IllegalArgumentException("Paired ranking requires its original member resource without another ordering experiment");
		}
		consensusPair=pair;
		replacementWalk=merged ? null : new NcrnaDdlShortlist(pair.replacement, unlimited ? Integer.MAX_VALUE : cap);
		replacementScores=merged ? new float[scores.length] : replacementWalk.scores;
	}
	/** Borrowed second order for the last build; null denotes a single merged walk. */
	public NcrnaDdlShortlist replacementWalk(){return replacementWalk;}
	/** Isolate member-local admission within the selected single merged walk. */
	public void configureMemberEligibility(final NcrnaFamily[] families){
		if(!orderOnly || consensusPair==null || replacementWalk!=null || ties!=null
			|| memberPrefix!=null || memberEligibility!=null || ddlAdmission!=null){
			throw new IllegalStateException("Member-native admission requires one configured merged consensus walk without another eligibility experiment");
		}
		final NcrnaMemberNativeEligibility eligibility=new NcrnaMemberNativeEligibility(families);
		if(eligibility.models()!=scores.length){throw new IllegalArgumentException("Member-native admission must bind the complete merged model axis");}
		memberEligibility=eligibility;
	}
	/** Explicit merged-member experiment; ordinary paths never allocate admission scratch. */
	public void configureDdlAdmission(final int minimumHits){
		if(!orderOnly || consensusPair==null || replacementWalk!=null || memberEligibility==null
			|| ties!=null || memberPrefix!=null || ddlAdmission!=null || profiles.buckets()!=1024){
			throw new IllegalStateException("DDL admission requires the complete1024-bin merged maximum-score walk with member-native eligibility, configured once before scanning");
		}
		ddlAdmission=new NcrnaDdlAdmission(scores.length, minimumHits, profiles.buckets());
	}
	/** Explicit replacement of native admission, not the historical native OR DDL rule. */
	public void configureDdlOnly(final int minimumHits){
		if(!orderOnly || consensusPair==null || replacementWalk!=null || memberEligibility!=null
			|| ties!=null || memberPrefix!=null || ddlAdmission!=null || modelAttemptSink!=null || profiles.buckets()!=1024){
			throw new IllegalStateException("DDL-only admission requires one complete merged 1024-bin pair without a native gate or another ordering experiment");
		}
		ddlAdmission=new NcrnaDdlAdmission(scores.length, minimumHits, profiles.buckets());
		ddlOnly=true;
	}
	public boolean ddlOnly(){return ddlOnly;}
	/** Frozen experimental reordering after raw-hit admission, before cap/patience.
	 * Only the order changes; score(model) retains its raw-count provenance. */
	public void configureWindowDimers(final NcrnaRankingFeatures.Profiles windows){
		if(!ddlOnly || consensusPair==null || replacementWalk!=null || ddlAdmission==null
			|| ties!=null || memberPrefix!=null || windowRanking!=null || windows==null
			|| windows.models()!=scores.length || profiles.buckets()!=1024){
			throw new IllegalStateException("Window dimers require one complete index-free merged 1024-bin model axis, configured once before scanning");
		}
		windowRanking=new NcrnaRankingFeatures.WindowRanking(windows);
	}
	public void configureTies(final NcrnaWindowGcProfile profile, final byte[][] consensus, final double gc, final double posting){
		if(!orderOnly || ties!=null || memberPrefix!=null || consensusPair!=null){throw new IllegalArgumentException("Configure ties once on the native-eligible ordering path without another ranking experiment");}
		scorer.enableSpecificity();
		ties=new NcrnaDdlTieBreak(profile, consensus, gc, posting);
	}
	public NcrnaDdlTieBreak ties(){return ties;}
	/** Explicit shared-clan experiment; default ordering remains untouched. */
	public void configureMemberPrefix(final int[] sizes){
		if(!orderOnly || ties!=null || memberPrefix!=null || consensusPair!=null || sizes==null || sizes.length!=3
			|| sizes[0]+sizes[1]+sizes[2]!=scores.length){
			throw new IllegalArgumentException("Member prefix requires one complete shared-clan ordering-only axis, without tie experiments");
		}
		memberPrefix=new NcrnaMemberNativePrefix(sizes);
		prefixOrder=new int[models.length];
	}
	public int extendNativeEligibility(final int[] destination, final int count, final TrnaKmerIndex index){
		if(memberEligibility!=null){return memberEligibility.select(index, destination);}
		return memberPrefix==null ? count : memberPrefix.extendEligibility(destination, count);
	}
	public int copyMemberPrefix(final int[] destination){
		if(memberPrefix==null || destination==null || destination.length<memberPrefix.size()){
			throw new IllegalStateException("Prefix telemetry requires the active native-prefix experiment and sufficient scratch");
		}
		for(int i=0; i<memberPrefix.size(); i++){destination[i]=memberPrefix.modelAt(i);}
		return memberPrefix.size();
	}
	public void appendPostingRuns(final structures.ByteBuilder out){
		if(ties==null){throw new IllegalStateException("Posting runs require the explicit tie experiment");}
		final char[] keys=scorer.queryKeys();
		final int[] lengths=scorer.postingLengths();
		boolean comma=false;
		for(int bin=0; bin<keys.length; bin++){
			if(lengths[bin]==0){continue;}
			if(comma){out.append(',');}
			out.append(bin).append(':').append((int)keys[bin]).append(':').append(lengths[bin]);
			comma=true;
		}
	}

	/** Dual resources use the complete sketch-only order and never touch an index.
	 * Hybrid resources retain their native-first filtered interleave unchanged. */
	public int[] build(final byte[] bases, final TrnaKmerIndex nativeIndex){
		if(bases==null || (ddlOnly ? nativeIndex!=null : (!dual || orderOnly) && nativeIndex==null)){
			throw new IllegalArgumentException("Ranking needs window bytes and exactly its declared index mode; DDL-only forbids a native index");
		}
		if(orderOnly){
			int count=ddlOnly ? 0 : nativeIndex.shortlistInto(bases, Integer.MAX_VALUE, -1, nativeEligible);
			if(memberEligibility!=null){count=memberEligibility.select(nativeIndex, nativeEligible);}
			if(ddlAdmission==null){
				for(int position=0; position<count; position++){
					final int model=nativeEligible[position];
					counts[model]=(char)nativeIndex.lastSharedCount(model);
				}
			}
			if(consensusPair==null){scorer.score(bases, scores);}
			else{
				scorer.scorePair(bases, consensusPair, scores, replacementScores);
				if(modelAttemptSink!=null && !ddlOnly){
					try{modelAttemptSink.rankingInputs(nativeEligible, count, scores, replacementScores, nativeIndex);}
					catch(RuntimeException | AssertionError failure){
						shared.KillSwitch.assertDie("Paired ranking-input observation failed on caller worker: "+failure);
						throw failure;
					}
				}
				if(replacementWalk==null){
					if(ddlAdmission==null){
						for(int model=0; model<scores.length; model++){scores[model]=Math.max(scores[model], replacementScores[model]);}
					}else{count=ddlAdmission.select(scores, replacementScores, nativeEligible, count);}
				}else{
					final int secondSize=replacementWalk.order.buildEligibleCaller(replacementScores, nativeEligible, count, cap);
					for(int position=0; position<secondSize; position++){
						final int model=replacementWalk.order.modelAt(position);
						replacementWalk.models[position]=model;
						replacementWalk.counts[model]=counts[model];
					}
				}
			}
			if(ddlAdmission!=null && !ddlOnly){
				// New entrants can have zero native hits; never expose a count left
				// from another query on their first visit to the emitted order.
				for(int position=0; position<count; position++){
					final int model=nativeEligible[position];
					counts[model]=(char)nativeIndex.lastSharedCount(model);
				}
			}
			final int[] ranked=windowRanking!=null ? windowRanking.build(bases, scores, nativeEligible, count)
				: ties==null ? nativeEligible : ties.build(bases, scores, scorer.specificitySums(), nativeEligible, count);
			final int size=order.buildEligibleCaller(scores, ranked, count, cap, ties!=null || windowRanking!=null);
			for(int position=0; position<size; position++){models[position]=order.modelAt(position);}
			if(memberPrefix!=null){
				memberPrefix.select(nativeIndex);
				final int emitted=memberPrefix.merge(models, size, prefixOrder, cap);
				order.buildEligibleCaller(scores, prefixOrder, emitted, cap, true);
				for(int position=0; position<emitted; position++){
					final int model=order.modelAt(position);
					models[position]=model;
					counts[model]=(char)nativeIndex.lastSharedCount(model);
				}
			}
			return models;
		}
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
	/** Original model ID in the current emitted prefix, including a paired order. */
	public int modelAt(final int position){
		if(position<0 || position>=size()){
			throw new IllegalArgumentException("Read only the emitted model prefix, not stale shortlist capacity");
		}
		return models[position];
	}

	/** Only this identity denotes the borrowed live-prefix array. */
	public boolean owns(final int[] proposed){return proposed==models;}

	/** Native support is -1 when deliberately unmeasured by sketch-only ordering. */
	public int nativeCount(final int model){
		if(model<0 || model>=scores.length || order.rank(model)==0){
			throw new IllegalArgumentException("Proposal provenance requires an emitted original model");
		}
		return ddlOnly || dual && !orderOnly ? -1 : counts[model];
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
	private NcrnaModelAttemptInstrumentSink modelAttemptSink;
	private final NcrnaDdlProfileResource profiles;
	private final boolean unlimited, selected, orderOnly;
	private final boolean dual;
	private final char[] counts;
	private final float[] scores;
	private final int[] models;
	private final int[] nativeEligible;
	private final NcrnaDdlInterleave order;
	private final NcrnaDdlProfileResource.Scorer scorer;
	private NcrnaDdlTieBreak ties;
	private NcrnaMemberNativePrefix memberPrefix;
	private int[] prefixOrder;
	private NcrnaDdlProfileResource.ConsensusPair consensusPair;
	private NcrnaDdlShortlist replacementWalk;
	private float[] replacementScores;
	private NcrnaMemberNativeEligibility memberEligibility;
	private NcrnaDdlAdmission ddlAdmission;
	private boolean ddlOnly;
	private NcrnaRankingFeatures.WindowRanking windowRanking;
}
