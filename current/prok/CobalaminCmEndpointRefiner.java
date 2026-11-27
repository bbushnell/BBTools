package prok;

import java.util.Arrays;
import java.util.ArrayList;

import aligner.CovarianceModel;
import aligner.CovarianceModelBandedCyk;
import aligner.CovarianceModelConsensusMap;
import aligner.CovarianceModelConsensusMaps;
import aligner.CovarianceModelGate;
import aligner.CovarianceModelLengthBands;
import aligner.CovarianceModelParser;
import aligner.CovarianceModelPlacementBands;
import aligner.CovarianceModelStateCorridors;
import idaligner.AlignmentStats;
import idaligner.EndClippedAligner;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;

/** Optional post-trim endpoint search. Detection identity, model and path score
 * remain the caller's. Models/maps are immutable; all search state is worker-owned.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminCmEndpointRefiner implements NcrnaBoundaryRefiner{

	public CobalaminCmEndpointRefiner(final CovarianceModel cm, final CovarianceModelConsensusMap[] maps_,
			final byte[][] library_, final int policy_){
		this(cm, maps_, library_, policy_, false);
	}
	public CobalaminCmEndpointRefiner(final CovarianceModel cm, final CovarianceModelConsensusMap[] maps_,
			final byte[][] library_, final int policy_, final boolean gaHold_){
		this(cm, maps_, library_, policy_, gaHold_, null);
	}
	public CobalaminCmEndpointRefiner(final CovarianceModel cm, final CovarianceModelConsensusMap[] maps_,
			final byte[][] library_, final int policy_, final boolean gaHold_, final CovarianceModelStateCorridors[] corridors_){
		this(cm, maps_, library_, policy_, gaHold_, corridors_, -1);
	}
	/** Standalone learned-region comparison may add one narrow glocal component.
	 * Ordinary factories keep -1 and retain their existing serving policy. */
	public CobalaminCmEndpointRefiner(final CovarianceModel cm, final CovarianceModelConsensusMap[] maps_,
			final byte[][] library_, final int policy_, final boolean gaHold_, final CovarianceModelStateCorridors[] corridors_,
			final int unionRadius_){
		this(cm, maps_, library_, policy_, gaHold_, corridors_, unionRadius_,
			policy_==POOLED ? CovarianceModelStateCorridors.pool(cm, corridors_) : null, null, PADDING, false);
	}
	/** Explicit separate anchor context and per-model search flanks. The decoded
	 * resources were bound to actual consensus bases before worker construction. */
	public CobalaminCmEndpointRefiner(final CovarianceModel cm, final CovarianceModelConsensusMap[] maps_,
			final byte[][] library_, final int policy_, final boolean gaHold_,
			final NcrnaCmCorridorResource.Data context, final NcrnaCmInputPadding.Data padding){
		this(cm, maps_, library_, policy_, gaHold_, context==null ? null : context.profiles, -1, null,
			padding, context==null ? PADDING : context.anchorPadding, context!=null && context.unclippedAnchors);
		require(context!=null && padding!=null, "Separate-context serving requires both decoded corridor and input-padding resources");
	}
	/** Factories share the immutable member pool instead of rebuilding it per input. */
	private CobalaminCmEndpointRefiner(final CovarianceModel cm, final CovarianceModelConsensusMap[] maps_,
			final byte[][] library_, final int policy_, final boolean gaHold_, final CovarianceModelStateCorridors[] corridors_,
			final int unionRadius_, final CovarianceModelStateCorridors pooled_, final NcrnaCmInputPadding.Data inputPadding_,
			final int anchorPadding_, final boolean unclippedAnchors_){
		require(cm!=null && maps_!=null && library_!=null && maps_.length==library_.length
			&& maps_.length>0 && policy_>=CONTACT && policy_<=POOLED,
			"CM endpoints require one complete member model/map axis and an explicit operating point");
		require(usesLearned(policy_) ? corridors_!=null && corridors_.length==maps_.length : corridors_==null,
			"Learned serving requires its complete model axis; old operating points cannot silently ignore supplied corridors");
		require(unionRadius_==-1 || policy_==LEARNED && unionRadius_>=0 && unionRadius_<=32,
			"Only an explicit learned replay may add a glocal union component within the unseen-state radius32 fallback");
		require(policy_==POOLED ? pooled_!=null && pooled_.belongsTo(cm) : pooled_==null,
			"Sparse pooling needs observed bounds from this exact member; other policies must not acquire a pool");
		require(anchorPadding_>=0 && (inputPadding_!=null)==unclippedAnchors_
			&& (inputPadding_==null ? anchorPadding_==PADDING : policy_==LEARNED && unionRadius_==-1
				&& inputPadding_.size()>=maps_.length
				&& (inputPadding_.fixed || inputPadding_.minCalls>=5 && inputPadding_.minGroups>=2)),
			"Input flanks require learned separate-context corridors and declared fixed or sufficiently supported padding on the complete mapped prefix; legacy serving keeps its original context");
		for(int m=0; m<maps_.length; m++){
			require(maps_[m]!=null && library_[m]!=null && maps_[m].consensusLength==library_[m].length,
				"Every selected consensus must have its own complete CM map");
			require(corridors_==null || corridors_[m]==null || corridors_[m].belongsTo(cm),
				"A fitted profile must bind this exact member CM; null profiles retain sparse support and select the declared fallback");
		}
		maps=maps_;
		library=library_;
		policy=policy_;
		learned=usesLearned(policy_);
		pooled=pooled_;
		unionRadius=unionRadius_;
		inputPadding=inputPadding_;
		anchorPadding=anchorPadding_;
		unclippedAnchors=unclippedAnchors_;
		corridors=corridors_==null ? null : corridors_.clone();
		cmLength=cm.clen;
		threshold=cm.header("GA")==null ? Float.NaN : CovarianceModelGate.gatheringThreshold(cm);
		gaHold=gaHold_;
		require(!gaHold || Float.isFinite(threshold), "An enabled CM score hold needs the selected model's finite GA");
		scorer=new CovarianceModelBandedCyk(cm, new CovarianceModelLengthBands(cm, 2), 500000000L, 1);
	}
	@Override
	public boolean applies(final int model){return model>=0 && model<maps.length;}
	@Override
	public boolean afterAlignmentTrim(){return true;}
	@Override
	public boolean needsAlignmentTrace(){return !learned;}
	@Override
	public boolean usesContigGc(){return false;}
	@Override
	public void setProposalTrace(final NcrnaBoundaryTrace.Family trace){
		require(trace==null || trace.invocationsOnly(), "CM proposals require invocation diagnostics, never NN score replay");
		proposalTrace=trace;
	}
	@Override
	public void alignment(final AlignmentStats stats, final int origin, final boolean accepting){
		require(origin>=0, "Borrowed trace origins are positions on the complete oriented contig");
		if(accepting){
			acceptedTrace=stats;
			acceptedOrigin=origin;
			trimTrace=null;
		}else{
			trimTrace=stats;
			trimOrigin=origin;
		}
	}
	@Override
	public boolean refine(final Orf orf, final byte[] bases, final int model, final float gc, final float identity,
			final boolean quantumIdentity, final int fixedMask, final int minLen, final int maxLen){
		return search(orf, bases, model, fixedMask, minLen, maxLen, true);
	}
	/** Runs the same anchored search without applying its proposed endpoints.
	 * Used by acceptance-only consumers; no tip/coordinate policy is replaced. */
	public float scoreOnly(final Orf orf, final byte[] bases, final int model){
		search(orf, bases, model, 0, 1, Integer.MAX_VALUE, false);
		return scoreFor(orf, model);
	}
	/** Rebuilds the selected target consensus trace at its actual qualifying
	 * location. Search overlap and padding still belong to the emitted candidate,
	 * so a neighboring RNA cannot supply an unrelated passing CM interval.
	 * Reconstruction is measured work; ordinary identity/endpoints are not a CIGAR.
	 */
	public float scoreOnlyAnchored(final Orf orf, final byte[] bases, final int model,
			final int anchorStart, final int anchorStop){
		search(orf, bases, model, 0, 1, Integer.MAX_VALUE, false, anchorStart, anchorStop, true);
		return scoreFor(orf, model);
	}
	/** A completed unavailable ordinary search is also reusable. Explicit target
	 * anchors are decision-local and cannot stand in for the ordinary own result. */
	public boolean hasResult(final Orf orf, final int model){
		return scoredCandidate==orf && scoredModel==model && !candidateExplicitAnchor;
	}
	/** Metadata for the last request, including an explicit anchored request.
	 * True includes measured impossibility; false means no measurement. Only
	 * hasResult() decides whether ordinary own scoring may reuse that request. */
	public boolean scored(final Orf orf, final int model){
		return scoredCandidate==orf && scoredModel==model && candidateScored;
	}
	public float scoreFor(final Orf orf, final int model){
		require(scoredCandidate==orf && scoredModel==model, "A cached covariance score belongs only to its exact candidate/model, never a preceding call");
		return candidateScored ? candidateScore : Float.NaN;
	}
	public float gatheringThreshold(){return threshold;}
	/** Sense-oriented inclusive target window; -1 means no window was prepared. */
	public int windowStartFor(final Orf orf, final int model){
		scoreFor(orf, model);
		return candidateWindowStart;
	}
	public int windowStopFor(final Orf orf, final int model){
		scoreFor(orf, model);
		return candidateWindowStop;
	}
	/** Raw sense-oriented inclusive proposal, before length/fixed-end/GA guards.
	 * -1 means no finite proposal. Reading it never applies or repeats a search. */
	public int proposalStartFor(final Orf orf, final int model){
		scoreFor(orf, model);
		return candidateProposalStart;
	}
	public int proposalStopFor(final Orf orf, final int model){
		scoreFor(orf, model);
		return candidateProposalStop;
	}
	/** Length, physical-end and GA failures for one finite raw proposal.
	 * Ordinary application and saved prediction use the same guards; this
	 * primitive query never scores, changes coordinates or updates counters. */
	public static int proposalFailures(final int start, final int stop,
			final int proposedStart, final int proposedStop, final float score,
			final float ga, final boolean gaHold, final int fixedMask, final int minLen, final int maxLen){
		require(start>=0 && stop>=start && proposedStart>=0 && proposedStop>=proposedStart
			&& (fixedMask&~3)==0 && minLen>0 && maxLen>=minLen
			&& Float.isFinite(score) && (!gaHold || Float.isFinite(ga)),
			"Proposal guards require valid spans, actual caller constraints and a finite measured score; an enabled GA hold needs its finite threshold");
		final long length=(long)proposedStop-proposedStart+1;
		int failures=length<minLen || length>maxLen ? PROPOSAL_LENGTH : 0;
		if((fixedMask&1)!=0 && proposedStart!=start || (fixedMask&2)!=0 && proposedStop!=stop){
			failures|=PROPOSAL_FIXED_END;
		}
		if(gaHold && score<ga){failures|=PROPOSAL_SCORE;}
		return failures;
	}
	/** Actual DP fills and their charged stages; retries count separately. */
	public long searchCalls(){return counters[RADIUS16_CALLS]+counters[RADIUS24_CALLS]+counters[RADIUS32_CALLS]
		+counters[RADIUS48_CALLS]+counters[RADIUS64_CALLS]+counters[EXACT_CALLS]+counters[LEARNED_CALLS];}
	public long searchCpu(){return counters[RADIUS16_CPU]+counters[RADIUS24_CPU]+counters[RADIUS32_CPU]
		+counters[RADIUS48_CPU]+counters[RADIUS64_CPU]+counters[EXACT_CPU]+counters[LEARNED_CPU];}

	private boolean search(final Orf orf, final byte[] bases, final int model,
			final int fixedMask, final int minLen, final int maxLen, final boolean apply){
		return search(orf, bases, model, fixedMask, minLen, maxLen, apply,
			orf==null ? -1 : orf.start, orf==null ? -1 : orf.stop, false);
	}
	private boolean search(final Orf orf, final byte[] bases, final int model,
			final int fixedMask, final int minLen, final int maxLen, final boolean apply,
			final int targetAnchorStart, final int targetAnchorStop, final boolean explicitAnchor){
		require(orf!=null && bases!=null && applies(model) && orf.start>=0 && orf.stop>=orf.start
			&& orf.stop<bases.length && (fixedMask&~3)==0 && minLen>0 && maxLen>=minLen,
			"CM refinement must bind an accepted model/span and its unchanged caller length/physical-end constraints");
		require(targetAnchorStart>=0 && targetAnchorStop>=targetAnchorStart && targetAnchorStop<bases.length
			&& (!explicitAnchor || !apply),
			"A separately qualified target anchor must be a real sense interval used only for score-only candidate acceptance");
		counters[CALLS]++;
		if(searchProfiles!=null){
			for(final CovarianceModelBandedCyk.SearchProfile profile:searchProfiles){profile.lastResult=null;}
		}
		final int start=orf.start, stop=orf.stop;
		int proposedStart=-1, proposedStop=-1, width=0, contacts16=0, contacts32=0, anchorKind=0, anchorOrigin=-1;
		int searchStart=-1, searchStop=-1;
		boolean unavailable=false, scored=false;
		float score=Float.NaN;
		String decision=null;
		try{
			if(fixedMask==3){
				decision="FIXED_BOTH";
				counters[FIXED_BOTH]++;
				return false;
			}
			final int pad5=inputPadding==null ? PADDING : inputPadding.pad5(model);
			final int pad3=inputPadding==null ? PADDING : inputPadding.pad3(model);
			final int from=(int)Math.max(0L, (long)start-pad5);
			final int to=(int)Math.min(bases.length-1L, (long)stop+pad3);
			searchStart=from;
			searchStop=to;
			final int anchorFrom=(int)Math.max(0L, (long)targetAnchorStart-anchorPadding);
			final int anchorTo=(int)Math.min(bases.length-1L, (long)targetAnchorStop+anchorPadding);
			final byte[] anchorWindow=fill(bases, anchorFrom, anchorTo);
			final long begin=NcrnaCpuProfile.cpuNow();
			final AlignmentStats anchor;
			final int origin;
			if(!explicitAnchor && !learned && usable(acceptedTrace, library[model].length)){
				anchor=acceptedTrace;
				origin=acceptedOrigin;
				counters[ACCEPTING_TRACE]++;
				anchorKind=1;
			}else if(!explicitAnchor && !learned && usable(trimTrace, library[model].length)){
				anchor=trimTrace;
				origin=trimOrigin;
				counters[TRIM_TRACE]++;
				anchorKind=2;
			}else if(library[model].length<=anchorWindow.length){
				reconstructed.clear();
				final float rebuiltIdentity=NcrnaAlignmentWork.quantumTrace(library[model], anchorWindow, reconstructed, library[model].length, 1);
				require(Float.isFinite(rebuiltIdentity), "Reconstruction arithmetic must be finite before composing an anchor");
				anchor=usable(reconstructed, library[model].length) ? reconstructed : null;
				origin=anchorFrom;
				counters[RECONSTRUCTED_TRACE]++;
				anchorKind=3;
			}else{
				anchor=null;
				origin=anchorFrom;
			}
			charge(ANCHOR_CPU, begin);
			// Equal-sized contexts can have different origins. Only reuse their
			// bytes when both bounds match; otherwise fill after tracing the anchor.
			final byte[] window=from==anchorFrom && to==anchorTo ? anchorWindow : fill(bases, from, to);
			anchorOrigin=origin;
			CovarianceModelBandedCyk.EndpointResult result=null;
			if(anchor==null){
				counters[UNAVAILABLE_ANCHOR]++;
				unavailable=true;
			}
			else{
				if(learned){
					result=learned(model, anchor, origin, from, to, window, start-from+1, stop-from+1);
					width=lastLearnedWidth;
					contacts32=initialUniformContacts;
				}else{
					final int first=policy==CONTACT ? 16 : policy==RETRY24_64 || policy==RETRY24_48 ? 24 : 32;
					result=band(model, anchor, origin, from, to, first, window, start-from+1, stop-from+1);
					width=first;
					if(first==16){contacts16=result.contacts;}else if(first==32){contacts32=result.contacts;}
					if(policy==CONTACT && (result==null || !Float.isFinite(result.score) || result.contacts!=0)){
						result=band(model, anchor, origin, from, to, 32, window, start-from+1, stop-from+1);
						width=32;
						contacts32=result.contacts;
					}
				}
			}
			final boolean wouldExact=result==null || !Float.isFinite(result.score) || result.contacts!=0;
			if(policy!=CONTACT && wouldExact){counters[WOULD_EXACT]++;}
			if(result!=null && (policy==RETRY_64 && result.contacts!=0
				|| (policy==RETRY24_64 || policy==RETRY24_48) && wouldExact)){
				final int retry=policy==RETRY24_48 ? 48 : 64;
				result=band(model, anchor, origin, from, to, retry, window, start-from+1, stop-from+1);
				width=retry;
			}
			if(policy==FIXED_32 && wouldExact){
				counters[EXACT_CALLS]++;
				final long exactStart=NcrnaCpuProfile.cpuNow();
				result=scorer.searchOverlapping(window, window.length, start-from+1, stop-from+1, null);
				width=-1;
				charge(EXACT_CPU, exactStart);
			}
			if(result==null){
				decision="UNAVAILABLE";
				counters[UNREFINED_UNAVAILABLE]++;
				return false;
			}
			require(Float.isFinite(result.score) || result.score==Float.NEGATIVE_INFINITY,
				"Invalid CM arithmetic cannot become an ordinary retained-original outcome");
			score=result.score;
			scored=true;
			if(!Float.isFinite(result.score)){
				decision="IMPOSSIBLE";
				counters[UNREFINED_IMPOSSIBLE]++;
				return false;
			}
			final int left=from+result.targetFrom-1, right=from+result.targetTo-1;
			proposedStart=left;
			proposedStop=right;
			require(left>=from && right>=left && right<=to, "Finite CM endpoints must remain inside the actual padded window");
			if(!apply){
				decision="SCORE_ONLY";
				return false;
			}
			final int failures=proposalFailures(start, stop, left, right, score, threshold, gaHold, fixedMask, minLen, maxLen);
			final boolean badLength=(failures&PROPOSAL_LENGTH)!=0;
			final boolean fixedMoved=(failures&PROPOSAL_FIXED_END)!=0;
			if(badLength || fixedMoved){
				decision=badLength ? fixedMoved ? "LENGTH_AND_FIXED_END_REJECTED" : "LENGTH_REJECTED" : "FIXED_END_REJECTED";
				counters[REJECTED_PROPOSAL]++;
				if(badLength){counters[REJECTED_LENGTH]++;}
				if(fixedMoved){counters[REJECTED_FIXED_END]++;}
				return false;
			}
			if((failures&PROPOSAL_SCORE)!=0){
				decision="SCORE_REJECTED";
				counters[REJECTED_PROPOSAL]++;
				counters[REJECTED_SCORE]++;
				return false;
			}
			orf.start=left;
			orf.stop=right;
			final boolean changed=left!=start || right!=stop;
			decision=changed ? "APPLIED" : "UNCHANGED";
			if(changed){counters[CHANGED_CALLS]++;}
			return changed;
		}finally{
			scoredCandidate=orf;
			scoredModel=model;
			candidateExplicitAnchor=explicitAnchor;
			candidateScored=scored;
			candidateScore=score;
			candidateWindowStart=searchStart;
			candidateWindowStop=searchStop;
			candidateProposalStart=proposedStart;
			candidateProposalStop=proposedStop;
			acceptedTrace=null;
			trimTrace=null;
			if(apply && proposalTrace!=null && decision!=null){
				proposalTrace.cmProposal(orf, start, stop, proposedStart, proposedStop, library[model].length, cmLength,
					scored, score, threshold, width, contacts16, contacts32, anchorKind, anchorOrigin, unavailable, decision);
			}
		}
	}
	/** Reconstructs the teacher's independent anchor before this method. A supported
	 * fit gets one learned fill; contact OR impossibility escalates to32 then64.
	 * Sparse models begin at32, except in explicit POOLED mode, where the immutable
	 * member pool is tried first. Sparse starts remain counted as sparse even then;
	 * LEARNED_CALLS includes that pooled search. No route here calls exact search.
	 */
	private CovarianceModelBandedCyk.EndpointResult learned(final int model, final AlignmentStats anchor, final int origin,
			final int from, final int to, final byte[] window, final int candidateFrom, final int candidateTo){
		assert(learned && corridors!=null && anchor!=null) :
			"Only the explicit learned policy supplies fitted bounds and an independently rebuilt usable anchor";
		initialUniformContacts=0;
		CovarianceModelBandedCyk.EndpointResult result=null;
		CovarianceModelStateCorridors profile=corridors[model];
		if(profile==null){
			counters[SPARSE_MODEL_CALLS]++;
			profile=pooled;
		}
		if(profile!=null){
			final long begin=NcrnaCpuProfile.cpuNow();
			final CovarianceModelPlacementBands placement=placement(model, anchor, origin, from, to, 32);
			require(placement!=null, "A rebuilt full-query anchor must compose the same placement used during training");
			if(searchProfiles!=null){compositionCpu[0]+=NcrnaCpuProfile.cpuNow()-begin;}
			result=scorer.searchOverlapping(window, window.length, candidateFrom, candidateTo, placement, profile, unionRadius,
				searchProfiles==null ? null : searchProfiles[0]);
			require(Float.isFinite(result.score) || result.score==Float.NEGATIVE_INFINITY,
				"Invalid learned arithmetic cannot be hidden by a uniform retry");
			counters[LEARNED_CALLS]++;
			if(result.contacts!=0){counters[LEARNED_CONTACTS]++;}
			if(!Float.isFinite(result.score)){counters[LEARNED_IMPOSSIBLE]++;}
			charge(LEARNED_CPU, begin);
			lastLearnedWidth=LEARNED_WIDTH;
		}
		if(retry(result)){
			result=band(model, anchor, origin, from, to, 32, window, candidateFrom, candidateTo);
			lastLearnedWidth=32;
			initialUniformContacts=result.contacts;
			if(retry(result)){
				result=band(model, anchor, origin, from, to, 64, window, candidateFrom, candidateTo);
				lastLearnedWidth=64;
			}
		}
		return result;
	}
	private static boolean retry(final CovarianceModelBandedCyk.EndpointResult result){
		return result==null || !Float.isFinite(result.score) || result.contacts!=0;
	}
	private CovarianceModelBandedCyk.EndpointResult band(final int model, final AlignmentStats trace, final int origin,
			final int from, final int to, final int radius, final byte[] window, final int candidateFrom, final int candidateTo){
		final long begin=NcrnaCpuProfile.cpuNow();
		final CovarianceModelPlacementBands placement=placement(model, trace, origin, from, to, radius);
		require(placement!=null, "A validated full-query anchor must produce a placement layout");
		final int route=radius==32 ? 1 : 2;
		if(searchProfiles!=null){
			assert(radius==32 || radius==64) : "The explicitly profiled learned/uniform64 policies have only32and64uniform searches";
			compositionCpu[route]+=NcrnaCpuProfile.cpuNow()-begin;
		}
		final CovarianceModelBandedCyk.EndpointResult result=scorer.searchOverlapping(window, window.length,
			candidateFrom, candidateTo, placement, null, -1, searchProfiles==null ? null : searchProfiles[route]);
		require(Float.isFinite(result.score) || result.score==Float.NEGATIVE_INFINITY,
			"Invalid CM arithmetic must fail before widening or fallback could hide it");
		final int calls, contacts, impossible, cpu;
		switch(radius){
			case 16:
				calls=RADIUS16_CALLS;
				contacts=RADIUS16_CONTACTS;
				impossible=RADIUS16_IMPOSSIBLE;
				cpu=RADIUS16_CPU;
				break;
			case 24:
				calls=RADIUS24_CALLS;
				contacts=RADIUS24_CONTACTS;
				impossible=RADIUS24_IMPOSSIBLE;
				cpu=RADIUS24_CPU;
				break;
			case 32:
				calls=RADIUS32_CALLS;
				contacts=RADIUS32_CONTACTS;
				impossible=RADIUS32_IMPOSSIBLE;
				cpu=RADIUS32_CPU;
				break;
			case 48:
				calls=RADIUS48_CALLS;
				contacts=RADIUS48_CONTACTS;
				impossible=RADIUS48_IMPOSSIBLE;
				cpu=RADIUS48_CPU;
				break;
			case 64:
				calls=RADIUS64_CALLS;
				contacts=RADIUS64_CONTACTS;
				impossible=RADIUS64_IMPOSSIBLE;
				cpu=RADIUS64_CPU;
				break;
			default: throw new IllegalArgumentException("Band counters must describe the actual requested search radius: "+radius);
		}
		counters[calls]++;
		if(result.contacts!=0){counters[contacts]++;}
		if(!Float.isFinite(result.score)){counters[impossible]++;}
		charge(cpu, begin);
		return result;
	}
	/** The resource declares whether residuals precede clipping to the search. */
	private CovarianceModelPlacementBands placement(final int model, final AlignmentStats trace, final int origin,
			final int from, final int to, final int radius){
		assert(model>=0 && model<maps.length && trace!=null) : "The selected mapped consensus supplies the anchor projection";
		return unclippedAnchors ? maps[model].composeWindowUnclipped(trace, origin, from, to, radius)
			: maps[model].composeWindow(trace, origin, from, to, radius);
	}
	private byte[] fill(final byte[] bases, final int from, final int to){
		assert(bases!=null && from>=0 && to>=from && to<bases.length) : "Only real sense-oriented contig bases enter an anchor or search window";
		final byte[] result=window(to-from+1);
		for(int i=0; i<result.length; i++){
			final byte b=bases[from+i];
			result[i]=(byte)(b>='a' && b<='z' ? b-32 : b);
		}
		return result;
	}
	/** Pool exact-length windows so the existing full-array Quantum API needs no
	 * padded-capacity reinterpretation and ordinary repeated lengths allocate once. */
	private byte[] window(final int length){
		require(length>0, "An accepted padded interval must contain real sequence");
		if(length>=windows.length){windows=Arrays.copyOf(windows, Math.max(length+1, windows.length*2));}
		if(windows[length]==null){windows[length]=new byte[length];}
		return windows[length];
	}
	private static boolean usable(final AlignmentStats stats, final int length){
		if(stats==null || stats.matchString==null || stats.qLen!=length){return false;}
		if(stats instanceof EndClippedAligner.Result){
			final EndClippedAligner.Result clipped=(EndClippedAligner.Result)stats;
			if(clipped.qStart!=0 || clipped.qStop!=length-1){return false;}
		}
		if(stats.rStop==stats.rStart-1){
			require(stats.rStart>=0 && stats.rStart<=stats.rLen && stats.ins==stats.qLen
				&& stats.matches==0 && stats.subs==0 && stats.dels==0 && stats.ns==0,
				"An empty reference span is usable only as an explicit all-insertion missing anchor");
			return false;
		}
		require(stats.rStart>=0 && stats.rStop>=stats.rStart && stats.rStop<stats.rLen,
			"A retained full-query trace must have valid reference endpoints before map composition");
		return true;
	}
	private void charge(final int field, final long begin){
		if(begin<0){return;}
		final long elapsed=NcrnaCpuProfile.cpuNow()-begin;
		require(elapsed>=0, "Profiled CM stages use one monotonic worker CPU clock");
		counters[field]+=elapsed;
	}
	/** Reads each member CM/map once. No payload digest or unrelated tip resource is opened. */
	public static Resources load(final NcrnaFamily family, final FamilyResourceResolver resolver,
			final FamilyResourceResolver.Overrides overrides, final int policy){
		return load(family, resolver, overrides, policy, false);
	}
	public static Resources load(final NcrnaFamily family, final FamilyResourceResolver resolver,
			final FamilyResourceResolver.Overrides overrides, final int policy, final boolean gaHold){
		require(family!=null && NcrnaFamilyRegistry.isCobalamin(family.name), "Only the explicitly selected cobalamin group uses this refiner");
		final CovarianceModel cm=CovarianceModelParser.read(resolver.resolve(family.name, "cm", overrides).path);
		require(!gaHold || cm.header("GA")!=null && Float.isFinite(CovarianceModelGate.gatheringThreshold(cm)),
			"An enabled endpoint score hold requires a finite GA in every selected member CM");
		require(cm.accession.equals(NcrnaFamilyRegistry.get(family.name).member), "The selected member must bind its own Rfam CM");
		final FamilyResourceResolver.Handle unmapped=resolver.optional(family.name, "cm_unmapped_models", overrides);
		final int mapped=unmapped==null ? family.library.length : mappedPrefix(unmapped.path, family.modelNames);
		final String[] mappedNames=mapped==family.modelNames.length ? family.modelNames : Arrays.copyOf(family.modelNames, mapped);
		final byte[][] mappedLibrary=mapped==family.library.length ? family.library : Arrays.copyOf(family.library, mapped);
		final CovarianceModelConsensusMap[] maps=CovarianceModelConsensusMaps.read(cm,
			resolver.resolve(family.name, "cm_map", overrides).path, mappedNames, mappedLibrary);
		for(int m=0; m<maps.length; m++){
			require(maps[m]!=null, "CM endpoint serving requires every selected consensus map; missing "+family.modelNames[m]);
		}
		CovarianceModelStateCorridors[] profiles=null;
		NcrnaCmCorridorResource.Data context=null;
		final FamilyResourceResolver.Handle padding=resolver.optional(family.name, "cm_input_padding", overrides);
		require(padding==null || policy==LEARNED, "Per-model input padding is an explicit learned-context experiment, not an ignored option on another CM policy");
		final NcrnaCmInputPadding.Data inputPadding=padding==null ? null
			: NcrnaCmInputPadding.read(padding.path, family.modelNames, family.library);
		if(usesLearned(policy)){
			final String[] names=new String[mapped], members=new String[names.length];
			final CovarianceModel[] models=new CovarianceModel[names.length];
			for(int m=0; m<names.length; m++){
				names[m]=NcrnaSupport.recordName(family.modelNames[m]);
				members[m]=cm.accession;
				models[m]=cm;
			}
			context=NcrnaCmCorridorResource.read(
				resolver.resolve(family.name, "cm_corridors", overrides).path, names, members, models, mappedLibrary);
			profiles=context.profiles;
			require(context.unclippedAnchors==(inputPadding!=null),
				"A separate-context corridor and per-model input-padding resource must be selected together; legacy bounds cannot be transplanted");
			for(int m=0; m<profiles.length; m++){
				require(profiles[m]==null || context.retainedCalls[m]>=100 && context.resolvedGroups[m]>=10,
					"Learned caller bounds require at least100retained usable winners from10resolved assembly groups");
			}
		}
		return new Resources(family, mappedLibrary, cm, maps, policy, gaHold, profiles, inputPadding,
			context==null ? PADDING : context.anchorPadding, context!=null && context.unclippedAnchors);
	}
	/** An explicit suffix alone may lack maps. The existing prefix still requires
	 * complete maps/corridors, and new models use ordinary alignment endpoints.
	 * One-column TSV: model header, then exact record names in library suffix order.
	 * Loaded once, never consulted or allocated per candidate.
	 */
	private static int mappedPrefix(final String path, final String[] names){
		assert(names!=null && names.length>0) : "The loaded family supplies the complete unique model axis";
		final ArrayList<String> added=new ArrayList<String>();
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, "model", path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==1 && line.length>0 && added.size()<names.length-1,
					"Unmapped models must be a nonempty suffix retaining at least one mapped consensus: "+path);
				added.add(p.parseString(0));
			}
		}
		require(!added.isEmpty(), "A declared unmapped-model resource must name the added suffix: "+path);
		final int mapped=names.length-added.size();
		for(int i=0; i<added.size(); i++){
			require(added.get(i).equals(NcrnaSupport.recordName(names[mapped+i])),
				"Unmapped models must match the exact appended library suffix in order; unexpected "+added.get(i));
		}
		return mapped;
	}
	public static final class Resources implements NcrnaBoundaryRefiner.Factory{
		private Resources(final NcrnaFamily family, final byte[][] mappedLibrary_, final CovarianceModel model_, final CovarianceModelConsensusMap[] maps_,
				final int policy_, final boolean gaHold_, final CovarianceModelStateCorridors[] corridors_,
				final NcrnaCmInputPadding.Data inputPadding_, final int anchorPadding_, final boolean unclippedAnchors_){
			name=family.name;
			library=family.library;
			mappedLibrary=mappedLibrary_;
			model=model_;
			threshold=model.header("GA")==null ? Float.NaN : CovarianceModelGate.gatheringThreshold(model);
			maps=maps_;
			policy=policy_;
			gaHold=gaHold_;
			corridors=corridors_;
			inputPadding=inputPadding_;
			anchorPadding=anchorPadding_;
			unclippedAnchors=unclippedAnchors_;
			pooled=policy_==POOLED ? CovarianceModelStateCorridors.pool(model, corridors) : null;
			require(policy_!=POOLED || pooled!=null, "Explicit sparse pooling needs at least one observed supported-model envelope before workers start");
		}
		@Override
		public NcrnaBoundaryRefiner create(final NcrnaFamily family){
			require(family!=null && family.name.equals(name) && family.library==library,
				"A worker factory must retain the exact loaded member/library used to validate its maps");
			return new CobalaminCmEndpointRefiner(model, maps, mappedLibrary, policy, gaHold, corridors, -1, pooled,
				inputPadding, anchorPadding, unclippedAnchors);
		}
		/** Read once from the immutable model, without constructing a scoring worker. */
		public float gatheringThreshold(){return threshold;}
		private final String name;
		private final byte[][] library, mappedLibrary;
		private final CovarianceModel model;
		private final float threshold;
		private final CovarianceModelConsensusMap[] maps;
		private final int policy;
		private final boolean gaHold;
		private final CovarianceModelStateCorridors[] corridors;
		private final CovarianceModelStateCorridors pooled;
		private final NcrnaCmInputPadding.Data inputPadding;
		private final int anchorPadding;
		private final boolean unclippedAnchors;
	}
	/** Scalar inspection for fixtures and joined-worker reporting. */
	public long count(final int field){
		require(field>=0 && field<COUNTER_COUNT, "Endpoint counter must belong to the published scalar schema");
		return counters[field];
	}
	/** Standalone measurement only; ordinary factories never enable these clocks. */
	public void enableSearchProfiling(){
		enableSearchProfiling(false);
	}
	/** Detailed state accounting is confined to the saved-call measurement path. */
	public void enableSearchProfiling(final boolean detailed){
		require(searchProfiles==null && (policy==LEARNED || policy==RETRY_64) && unionRadius==-1
			&& NcrnaCpuProfile.enabled(), "Profile unchanged learned or uniform32->64 replay before any algorithm experiment");
		searchProfiles=new CovarianceModelBandedCyk.SearchProfile[3];
		compositionCpu=new long[3];
		for(int route=0; route<searchProfiles.length; route++){
			searchProfiles[route]=new CovarianceModelBandedCyk.SearchProfile(detailed);
		}
	}
	/** Routes are learned, uniform32, uniform64; all counts belong to this worker. */
	public CovarianceModelBandedCyk.SearchProfile searchProfile(final int route){
		require(searchProfiles!=null && route>=0 && route<3, "Search-work inspection requires an enabled existing measurement route");
		return searchProfiles[route];
	}
	public long compositionCpu(final int route){
		searchProfile(route);
		return compositionCpu[route];
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private static boolean usesLearned(final int policy){return policy==LEARNED || policy==POOLED;}
	private final CovarianceModelConsensusMap[] maps;
	private final byte[][] library;
	private final int policy;
	private final int unionRadius;
	private final CovarianceModelStateCorridors[] corridors;
	private final CovarianceModelStateCorridors pooled;
	private final boolean learned;
	private final boolean gaHold;
	private final NcrnaCmInputPadding.Data inputPadding;
	private final int anchorPadding;
	private final boolean unclippedAnchors;
	private final int cmLength;
	private final float threshold;
	private NcrnaBoundaryTrace.Family proposalTrace;
	private Orf scoredCandidate;
	private int scoredModel=-1;
	private boolean candidateScored;
	private boolean candidateExplicitAnchor;
	private float candidateScore=Float.NaN;
	private int candidateWindowStart=-1, candidateWindowStop=-1;
	private int candidateProposalStart=-1, candidateProposalStop=-1;
	private final CovarianceModelBandedCyk scorer;
	private CovarianceModelBandedCyk.SearchProfile[] searchProfiles;
	private long[] compositionCpu;
	private final AlignmentStats reconstructed=new AlignmentStats(true);
	private AlignmentStats acceptedTrace, trimTrace;
	private int acceptedOrigin, trimOrigin;
	private int lastLearnedWidth, initialUniformContacts;
	private byte[][] windows=new byte[16][];
	final long[] counters=new long[COUNTER_COUNT];
	public static final int COUNTER_COUNT=42;
	public static final int OFF=0, CONTACT=1, FIXED_32=2, ONLY_32=3, RETRY_64=4, RETRY24_64=5, RETRY24_48=6, LEARNED=7, POOLED=8, PADDING=20;
	public static final int LEARNED_WIDTH=-2;
	public static final int PROPOSAL_LENGTH=1, PROPOSAL_FIXED_END=2, PROPOSAL_SCORE=4;
	public static final int CALLS=0, ACCEPTING_TRACE=1, TRIM_TRACE=2, RECONSTRUCTED_TRACE=3, UNAVAILABLE_ANCHOR=4,
		RADIUS16_CALLS=5, RADIUS32_CALLS=6, RADIUS16_CONTACTS=7, RADIUS32_CONTACTS=8,
		RADIUS16_IMPOSSIBLE=9, RADIUS32_IMPOSSIBLE=10, EXACT_CALLS=11, UNREFINED_UNAVAILABLE=12,
		UNREFINED_IMPOSSIBLE=13, REJECTED_PROPOSAL=14, FIXED_BOTH=15, ANCHOR_CPU=16,
		RADIUS16_CPU=17, RADIUS32_CPU=18, EXACT_CPU=19, CHANGED_CALLS=20, REJECTED_LENGTH=21, REJECTED_FIXED_END=22,
		REJECTED_SCORE=23, RADIUS64_CALLS=24, RADIUS64_CONTACTS=25, RADIUS64_IMPOSSIBLE=26, RADIUS64_CPU=27, WOULD_EXACT=28,
		RADIUS24_CALLS=29, RADIUS24_CONTACTS=30, RADIUS24_IMPOSSIBLE=31, RADIUS24_CPU=32,
		RADIUS48_CALLS=33, RADIUS48_CONTACTS=34, RADIUS48_IMPOSSIBLE=35, RADIUS48_CPU=36,
		LEARNED_CALLS=37, LEARNED_CONTACTS=38, LEARNED_IMPOSSIBLE=39, LEARNED_CPU=40, SPARSE_MODEL_CALLS=41;
	// Under POOLED only, sparse starts are also counted among learned searches.
	// The resource still has no per-model fit for them; pooled is a shared prior.
	static final String COUNTER_HEADER="calls\taccepting_trace\ttrim_trace\treconstructed_trace\tunavailable_anchor"
		+"\tradius16_calls\tradius32_calls\tradius16_contacts\tradius32_contacts\tradius16_impossible\tradius32_impossible"
		+"\texact_calls\tunrefined_unavailable\tunrefined_impossible\trejected_proposal\tfixed_both"
		+"\tanchor_cpu_ns\tradius16_cpu_ns\tradius32_cpu_ns\texact_cpu_ns\tchanged_calls\trejected_length\trejected_fixed_end\trejected_score"
		+"\tradius64_calls\tradius64_contacts\tradius64_impossible\tradius64_cpu_ns\twould_exact"
		+"\tradius24_calls\tradius24_contacts\tradius24_impossible\tradius24_cpu_ns"
		+"\tradius48_calls\tradius48_contacts\tradius48_impossible\tradius48_cpu_ns"
		+"\tlearned_calls\tlearned_contacts\tlearned_impossible\tlearned_cpu_ns\tsparse_model_calls";
}
