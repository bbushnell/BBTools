package prok;

import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import idaligner.AlignmentStats;
import idaligner.EndClippedAligner;
import map.LongHashSet;

/** Combined cobalamin caller: one ordered model axis, frozen pass schedules,
 * shared claims and patience. Resource-only member workers retain their own
 * endpoint inference and local model ordinals, but never rank or scavenge.
 * Mutable ranking, refinement and accounting scratch belongs to one caller
 * thread. Family resources and consensus/HBM arrays are borrowed read-only.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaClanScavenger extends NcrnaScavenger{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaClanScavenger(final NcrnaFamily[] members, final NcrnaScavenger[] refiners){
		this(members, refiners, true);
	}

	/** false constructs no clan nine-mer index; it must later bind DDL-only admission. */
	public NcrnaClanScavenger(final NcrnaFamily[] members, final NcrnaScavenger[] refiners,
			final boolean nativeRanking){
		this(new Models(members), refiners, nativeRanking);
	}
	private NcrnaClanScavenger(final Models axis_, final NcrnaScavenger[] refiners_, final boolean nativeRanking){
		super(axis_.library, axis_.models, axis_.names, axis_.families[0].kmerSet,
			17, axis_.minimum, axis_.families[0].windowPad, 9, 384, true, 2, 0, 0, 0,
			axis_.families[0].scoreA, axis_.families[0].scoreB, axis_.families[0].idPass, axis_.families[0].idBorderline,
			null, null, null, null, -1, -1, -1, -1, 0,
			NcrnaFamily.LEGACY_START_OFFSETS, NcrnaFamily.LEGACY_STOP_OFFSETS, 0, 0, nativeRanking);
		require(refiners_!=null && refiners_.length==3, "A joint caller needs three existing member refinement contexts");
		axis=axis_;
		refiners=refiners_.clone();
		boolean traced=false;
		for(int member=0; member<3; member++){
			final NcrnaScavenger refiner=refiners[member];
			require(refiner!=null && !refiner.hasNativeRankingIndex() && refiner.modelCount()==axis.families[member].library.length
				&& axis.families[member].name.equals(refiner.family),
				"Member refinement contexts must bind local resources without building unused ranking indexes");
			traced|=refiner.hasBoundaryTrace();
		}
		memberBoundaryTrace=traced;
		final NcrnaFamily main=axis.families[0];
		family="cobalamin_clan";
		hbmPass=main.hbmPass;
		collapseFrac=main.collapseFrac;
		minKmerHits=main.seedMinHits;
		seedDistinct=main.seedDistinct;
		quantumThresh=main.quantumThresh;
		scavengePass2=main.scavengePass2;
		trimAlignmentExtent=main.trimAlignmentExtent;
		physicalEdgeMinLen=main.physicalEdgeMinLen==0 ? 0 : axis.minimum;
		refreshClaimedWindows=false;
		setModelPatience(8);
		setModelThresholds(NcrnaModelThresholds.forModels(axis.names, axis.pass, axis.borderline, axis.minima, axis.maxima));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Shared Scheduling     ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	protected boolean requiresWindowSeedKeys(){return overlapBuilders!=null;}

	/** Shared hit stream and claims. The ordinary route combines exact coordinates;
	 * the explicit overlap experiment pools union bounds across members instead.
	 * Historical separately accepted claims are never borrowed as shared state. */
	@Override
	public ArrayList<int[]> prepareCandidateWindows(final int[] positions, final long[] keys, final int seqLen){
		require(!refreshClaimedWindows, "Joint passes must be frozen before alignment, without a snapshot conflict resolver");
		proposalsStarted=true;
		final ArrayList<int[]> proposals=new ArrayList<int[]>();
		for(int member=0; member<3; member++){
			final ArrayList<int[]> windows;
			if(overlapBuilders!=null){
				final NcrnaFamily resource=axis.families[member];
				final NcrnaOverlapWindowBuilder builder=overlapBuilders[member];
				windows=builder.build(positions, keys, seqLen, proposalPads[member], withinMemberUnion,
					resource.minLen, resource.physicalEdgeMinLen, keyCaps==null ? null : keyCaps[member],
					proposalCeilings==null ? 0 : proposalCeilings[member]);
				if(memberGeometry!=null){memberGeometry[member].add(builder.geometry);}
			}else{
				windows=proposalPads==null ? refiners[member].prepareCandidateWindows(positions, keys, seqLen)
					: refiners[member].preparePaddedCandidateWindows(positions, seqLen, proposalPads[member]);
			}
			for(final int[] window:windows){
				assert(window.length==(proposalCeilings==null ? 3 : 4))
					: "Member windows must retain the optional chain ceiling before shared pooling";
				proposals.add(proposalCeilings==null ? new int[]{window[0], window[1], window[2], 1<<member}
					: new int[]{window[0], window[1], window[2], 1<<member, window[3]});
			}
		}
		return overlapBuilders==null ? deduplicate(proposals) : NcrnaOverlapWindowBuilder.pool(proposals);
	}

	/** Claims may split proposals. Preserve source membership and each member's
	 * scheduling floor, then apply the selected shared pooling rule again. */
	@Override
	protected ArrayList<int[]> subtractClaimed(final ArrayList<int[]> windows,
			final ArrayList<int[]> claimed, final int seqLen){
		final ArrayList<int[]> proposals=new ArrayList<int[]>();
		final ArrayList<int[]> one=new ArrayList<int[]>(1);
		for(final int[] window:windows){
			one.clear();
			one.add(window);
			for(final int[] remainder:super.subtractClaimed(one, claimed, seqLen)){
				int mask=0;
				for(int member=0; member<3; member++){
					if((window[3]&(1<<member))!=0 && retains(remainder, seqLen, axis.families[member])){mask|=1<<member;}
				}
				if(mask!=0){
					if(proposalCeilings==null){proposals.add(new int[]{remainder[0], remainder[1], remainder[2], mask});}
					else{
						int limit=remainder[1]-remainder[0]+1;
						for(int member=0; member<3; member++){
							if((mask&(1<<member))!=0){limit=Math.max(limit, proposalCeilings[member]);}
						}
						proposals.add(new int[]{remainder[0], remainder[1], remainder[2], mask, limit});
					}
				}
			}
		}
		return overlapBuilders==null ? deduplicate(proposals) : NcrnaOverlapWindowBuilder.pool(proposals);
	}

	private static boolean retains(final int[] window, final int seqLen, final NcrnaFamily member){
		assert(window[0]>=0 && window[1]<seqLen && window[1]>=window[0])
			: "Clan claim remainders retain original inclusive contig coordinates";
		return window[1]-window[0]>=member.minLen || member.physicalEdgeMinLen>0
			&& (window[0]==0 || window[1]==seqLen-1) && window[1]-window[0]+1>=member.physicalEdgeMinLen;
	}

	/** Input arrays are private proposal state; only exact coordinates combine. */
	private static ArrayList<int[]> deduplicate(final ArrayList<int[]> proposals){
		proposals.sort((a, b)->{
			final int start=Integer.compare(a[0], b[0]);
			return start==0 ? Integer.compare(a[1], b[1]) : start;
		});
		final ArrayList<int[]> result=new ArrayList<int[]>();
		for(final int[] window:proposals){
			assert(window.length==4 && window[3]>0 && window[3]<8)
				: "A clan proposal must retain a nonempty three-member provenance mask";
			final int[] last=result.isEmpty() ? null : result.get(result.size()-1);
			if(last!=null && last[0]==window[0] && last[1]==window[1]){last[3]|=window[3];}
			else{result.add(window);}
		}
		return result;
	}

	@Override
	protected void scheduledWindow(final int[] window){
		assert(window.length==(proposalCeilings==null ? 4 : 5) && window[3]>0 && window[3]<8
			&& (proposalCeilings==null || window[1]-(long)window[0]+1<=window[4]))
			: "Every verified clan window must retain its exact member proposal provenance";
		windowMembers=window[3];
	}

	public int windowMemberMask(){return windowMembers;}
	public int member(final int model){return axis.owner[model];}
	public int localModel(final int model){return axis.local[model];}
	public String modelName(final int model){return axis.names[model];}
	/** Missing covariance evidence may retain an ordinary pass, never invent a
	 * newly rescued lower-identity call. This does not change the walk's gate. */
	public float nativeIdentityPass(final int model){
		assert(model>=0 && model<axis.pass.length) : "Fallback acceptance must bind the chosen global model's original gate";
		return candidateNativePass(model);
	}
	/** Only explicit drivers allocate geometry/model counters. */
	public void setWindowExperiment(final NcrnaWindowCounters counters){
		require(queryGeometry==null && !proposalsStarted && counters!=null && counters.requests.length==axis.names.length*8,
			"Install window observations once before the shared caller prepares any proposals");
		memberGeometry=new WindowGeometry[3];
		for(int member=0; member<3; member++){
			memberGeometry[member]=overlapBuilders==null ? refiners[member].enableWindowGeometry() : new WindowGeometry();
		}
		queryGeometry=counters;
	}
	public WindowGeometry memberWindowGeometry(final int member){return memberGeometry[member];}
	public NcrnaWindowCounters windowExperiment(){return queryGeometry;}
	/** Member proposal padding. Equal resource geometry remains checked;
	 * these explicit proposal bounds do not change shared model admission. */
	public void setMemberProposalPads(final int[] padding){
		require(proposalPads==null && !proposalsStarted && padding!=null && padding.length==3,
			"Install three explicit member proposal pads exactly once, before preparing any windows");
		for(int member=0; member<3; member++){
			require(padding[member]>=0 && padding[member]<=axis.families[member].windowPad,
				"A member proposal experiment may only narrow its existing padding: "+axis.families[member].name);
		}
		proposalPads=padding.clone();
	}
	/** Explicit driver experiment. The member stage chooses nonempty intersection
	 * or union; the shared stage ALWAYS unions overlapping proposals and masks.
	 * Caps are immutable and loaded once outside workers; geometry scratch is local. */
	public void setOverlapWindows(final boolean union, final NcrnaOverlapWindowBuilder.Caps[] caps){
		require(overlapBuilders==null && !proposalsStarted && memberGeometry==null && proposalPads!=null
			&& (caps==null || caps.length==3 && union),
			"Configure overlap clustering once after member pads and before observations; per-key bounds require union");
		withinMemberUnion=union;
		keyCaps=caps==null ? null : caps.clone();
		overlapBuilders=new NcrnaOverlapWindowBuilder[3];
		for(int member=0; member<3; member++){
			require(keyCaps==null || keyCaps[member]!=null, "All three members need their own immutable per-key bounds");
			overlapBuilders[member]=new NcrnaOverlapWindowBuilder();
		}
	}
	/** Bound only chain growth in the explicit per-key union experiment. A raw
	 * interval wider than its member's limit survives without being clipped. */
	public void setWindowCeilings(final int[] ceilings){
		require(proposalCeilings==null && !proposalsStarted && memberGeometry==null
			&& overlapBuilders!=null && withinMemberUnion && keyCaps!=null
			&& ceilings!=null && ceilings.length==3,
			"Install chain ceilings after per-key union configuration and before workers prepare windows");
		for(final int ceiling:ceilings){require(ceiling>0, "Every member needs a positive window ceiling in bases");}
		proposalCeilings=ceilings.clone();
	}
	/** Explicit driver-only experiment; member refiners retain their ordinary rule. */
	public void setMatchNumeratorWinner(final boolean enabled){enableMatchedBaseWinner(enabled);}
	/** Configure after the explicit score rule and finite patience value. */
	public void setScorePatience(final int shortlist){enableScorePatience(shortlist);}
	/** Configure after shared ordering, score winner and primary patience. */
	public void setNativeRescue(final int count, final double minimumScore){enableNativeRescue(count, minimumScore);}
	/** Separate immutable profiles reorder only gated, rejected-window rescue. */
	public void setBonusRescueProfiles(final NcrnaDdlProfileResource profiles){enableBonusRescue(profiles);}
	/** Existing rejected-only gates, with the runtime replacement index as its order. */
	public void setReplacementRescueProfiles(final NcrnaDdlProfileResource.ConsensusPair pair){enableReplacementRescue(pair);}
	/** Two independent patient walks, or one maximum-hit-count walk. */
	public void setConsensusPolicy(final NcrnaDdlProfileResource.ConsensusPair pair, final boolean twoWalks){enableConsensusPolicy(pair, twoWalks);}
	/** Explicit admission-only experiment; ranking and patience remain configured separately. */
	public void setMemberNativeEligibility(){enableMemberNativeEligibility(axis.families);}
	/** Add raw maximum-DDL-supported models without changing native provenance or ordering. */
	public void setDdlAdmission(final int minimumHits){enableDdlAdmission(minimumHits);}
	/** Bind the index-free clan's ordering resource before configuring its paired walk. */
	public void useDdlOnlyProfiles(final NcrnaDdlProfileResource profiles){useMemberDdlProfiles(profiles, true);}
	public void setDdlOnlyAdmission(final int minimumHits){enableDdlOnlyAdmission(minimumHits);}
	public void setWindowDimers(final NcrnaRankingFeatures.Profiles windows){enableWindowDimers(windows);}
	public void setWinnerRule(final int rule){
		setWinnerRule(rule, 0.54);
	}
	public void setWinnerRule(final int rule, final double offset){
		if(rule==NcrnaMatchedBaseWinner.IDENTITY){enableMatchedBaseWinner(false);}
		else if(rule==NcrnaMatchedBaseWinner.MATCHES){enableMatchedBaseWinner(true);}
		else if(rule==NcrnaMatchedBaseWinner.CIGAR_SCORE){enableAdjustedWinner(false);}
		else if(rule==NcrnaMatchedBaseWinner.MODEL_SCORE){enableAdjustedWinner(true);}
		else{throw new IllegalArgumentException("Unknown shared winner rule: "+rule);}
		if(rule==NcrnaMatchedBaseWinner.CIGAR_SCORE || rule==NcrnaMatchedBaseWinner.MODEL_SCORE){setAdjustedWinnerOffset(offset);}
	}

	/*--------------------------------------------------------------*/
	/*----------------     Member Finalization      ----------------*/
	/*--------------------------------------------------------------*/

	/** Absolute path coefficients; installed before calling.
	 * This changes only the already chosen model's path score, never its ranking,
	 * verifier, endpoints or member resources. The cell is borrowed immutable.
	 */
	public void setPathScoreCell(final NcrnaPathScoreConfig.Cell cell){
		require(pathScoreCell==null && cell!=null, "Install exactly one explicit runtime score cell before the clan is used");
		cell.validateRuntime(axis.families);
		pathScoreCell=cell;
	}

	/** Each accepted model keeps its own calibrated path coefficients. The model
	 * winner is already chosen; length/identity and float arithmetic stay unchanged. */
	@Override
	protected float pathScore(final int model, final int length, final float identity){
		assert(model>=0 && model<axis.owner.length && length>0 && Float.isFinite(identity))
			: "A shared path score must bind the actual global model and its aligned span/identity";
		final NcrnaFamily member=axis.families[axis.owner[model]];
		if(pathScoreCell!=null){return pathScoreCell.runtimeScore(axis.owner[model], axis.local[model], length, identity);}
		return member.scoreA+member.scoreB*length*identity*identity;
	}
	@Override
	protected String pathScoreFamily(final int model){return axis.families[axis.owner[model]].name;}
	@Override
	protected float pathScoreA(final int model){
		require(pathScoreCell==null, "The historical score-basis schema cannot capture an offset runtime score cell");
		return axis.families[axis.owner[model]].scoreA;
	}
	@Override
	protected float pathScoreB(final int model){
		require(pathScoreCell==null, "The historical score-basis schema cannot capture an offset runtime score cell");
		return axis.families[axis.owner[model]].scoreB;
	}

	@Override
	public boolean finishAcceptedBoundary(final Orf orf, final byte[] bases, final int model,
			final int start, final int stop, final int source, final float identity, final int alignedLength,
			final float locusAni, final boolean quantum, final int[] positions, final long[] keys,
			final int fixedMask, final EndClippedAligner.Result partial){
		return finishAcceptedBoundary(orf, bases, model, start, stop, source, identity, alignedLength,
			locusAni, quantum, positions, keys, fixedMask, partial, partial, start);
	}
	@Override
	public boolean finishAcceptedBoundary(final Orf orf, final byte[] bases, final int model,
			final int start, final int stop, final int source, final float identity, final int alignedLength,
			final float locusAni, final boolean quantum, final int[] positions, final long[] keys,
			final int fixedMask, final EndClippedAligner.Result partial, final AlignmentStats trace, final int origin){
		final int member=axis.owner[model];
		attribute(model);
		orf.ncrnaFamily=axis.families[member].name;
		try{
			return refiners[member].finishAcceptedBoundary(orf, bases, axis.local[model], start, stop, source,
				identity, alignedLength, locusAni, quantum, positions, keys, fixedMask, partial, trace, origin);
		}finally{flushWork();}
	}
	@Override
	protected boolean tracedBoundary(){return refiners[0].tracedBoundary();}

	/** The member that finalized this object also owns its candidate-bound record. */
	@Override
	protected void traceReturned(final Orf orf){
		if(!memberBoundaryTrace){return;}
		for(int member=0; member<3; member++){
			if(axis.families[member].name.equals(orf.ncrnaFamily)){
				refiners[member].traceReturned(orf);
				return;
			}
		}
		throw new IllegalStateException("A returned clan call must retain its actual member identity for boundary tracing");
	}
	/** Original member refinement contexts own the actual CM work in a shared walk. */
	long[] cmEndpointCounters(final int member){return refiners[member].cmEndpointCounters();}
	/** The native main search is reused by acceptance; minor consumers score only. */
	CobalaminCmEndpointRefiner cmEndpointRefiner(final int member){return refiners[member].cmEndpointRefiner();}
	@Override
	public NcrnaBoundaryTrace.Family candidateTrace(final int model){return refiners[axis.owner[model]].candidateTrace(axis.local[model]);}
	public void requireAcceptanceTrace(){
		for(final NcrnaScavenger refiner:refiners){
			require(refiner.candidateTrace(0)!=null, "Real acceptance collection requires the complete three-member candidate trace before any genome is called");
		}
	}

	@Override
	protected void observeAlignedModel(final int model, final int stage){
		attribute(model);
		super.observeAlignedModel(model, stage);
	}
	@Override
	protected void observeAcceptedModel(final int model){
		if(queryGeometry!=null){queryGeometry.accepted(model, windowMembers);}
	}

	@Override
	protected void completedWindow(){
		flushWork();
		activeMember=-1;
		if(queryGeometry!=null){
			final NcrnaAlignmentWork work=NcrnaAlignmentWork.current();
			assert(work!=null) : shared.KillSwitch.assertDie("Observed window exit must retain its enclosing clan profile");
			work.clearQueryObserver(queryGeometry);
		}
	}

	/** Components are subdivisions of the active clan scope, never nested spans. */
	private void attribute(final int model){
		flushWork();
		activeMember=axis.owner[model];
		if(queryGeometry!=null){
			final NcrnaAlignmentWork work=NcrnaAlignmentWork.current();
			assert(work!=null) : shared.KillSwitch.assertDie("Explicit query geometry requires the actual enclosing clan profile");
			queryGeometry.activate(model, windowMembers);
			work.setQueryObserver(queryGeometry);
		}
	}

	private void flushWork(){
		final NcrnaAlignmentWork current=NcrnaAlignmentWork.current();
		if(current==null){return;}
		if(activeMember>=0){
			final NcrnaAlignmentWork target=memberWork[activeMember];
			target.logical+=current.logical-previous.logical;
			target.physical+=current.physical-previous.physical;
			target.hits+=current.hits-previous.hits;
			target.modelBases+=current.modelBases-previous.modelBases;
			target.queryBases+=current.queryBases-previous.queryBases;
			target.physicalModelBases+=current.physicalModelBases-previous.physicalModelBases;
			target.cpuNanos+=current.cpuNanos-previous.cpuNanos;
			target.hbmCalls+=current.hbmCalls-previous.hbmCalls;
			target.hbmCpuNanos+=current.hbmCpuNanos-previous.hbmCpuNanos;
		}
		previous.logical=current.logical;
		previous.physical=current.physical;
		previous.hits=current.hits;
		previous.modelBases=current.modelBases;
		previous.queryBases=current.queryBases;
		previous.physicalModelBases=current.physicalModelBases;
		previous.cpuNanos=current.cpuNanos;
		previous.hbmCalls=current.hbmCalls;
		previous.hbmCpuNanos=current.hbmCpuNanos;
	}

	/** Read only after the owning caller joins. These counters omit outer CPU. */
	public NcrnaAlignmentWork memberWork(final int member){return memberWork[member];}
	@Override
	public long alignmentCount(){
		long result=super.alignmentCount();
		for(final NcrnaScavenger refiner:refiners){result+=refiner.alignmentCount();}
		return result;
	}
	@Override
	public long alignedBases(){
		long result=super.alignedBases();
		for(final NcrnaScavenger refiner:refiners){result+=refiner.alignedBases();}
		return result;
	}
	@Override
	public long alignmentCpuNanos(){
		long result=super.alignmentCpuNanos();
		for(final NcrnaScavenger refiner:refiners){result+=refiner.alignmentCpuNanos();}
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Immutable Model Axis    ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Models{
		Models(final NcrnaFamily[] families_){
			require(families_!=null && families_.length==3, "The combined cobalamin axis requires exactly three loaded members");
			families=families_.clone();
			final String[] expected={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
			int size=0, shortest=Integer.MAX_VALUE;
			for(int member=0; member<3; member++){
				final NcrnaFamily family=families[member];
				require(family!=null && expected[member].equals(family.name) && family.library!=null
					&& family.models!=null && family.modelNames!=null && family.library.length>0
					&& family.models.length==family.library.length && family.modelNames.length==family.library.length,
					"Global ordinals require ordered, complete member consensus/HBM/name bindings");
				validateCommon(families[0], family);
				size=Math.addExact(size, family.library.length);
				shortest=Math.min(shortest, family.minLen);
			}
			minimum=shortest;
			library=new byte[size][];
			models=new BaseGraph[size];
			names=new String[size];
			owner=new int[size];
			local=new int[size];
			minima=new int[size];
			maxima=new int[size];
			pass=new float[size];
			borderline=new float[size];
			int next=0;
			for(int member=0; member<3; member++){
				final NcrnaFamily family=families[member];
				final NcrnaModelThresholds thresholds=family.modelThresholds;
				if(thresholds!=null){thresholds.validate(family.modelNames);}
				for(int model=0; model<family.library.length; model++, next++){
					library[next]=family.library[model];
					models[next]=family.models[model];
					names[next]=family.modelNames[model];
					owner[next]=member;
					local[next]=model;
					pass[next]=thresholds==null ? family.idPass : thresholds.pass(model);
					borderline[next]=thresholds==null ? family.idBorderline : thresholds.borderline(model);
					minima[next]=Math.max(family.minLen, thresholds==null ? 0 : thresholds.minimum(model));
					maxima[next]=Math.min(family.maxLen, thresholds==null ? Integer.MAX_VALUE : thresholds.maximum(model));
					require(pass[next]==borderline[next], "This joint experiment retains equal pass/borderline cutoffs without an unmeasured member HBM rescue policy");
				}
			}
			assert(next==size) : "Combined consensus, owner and local model axes must conserve every member row";
		}
		final NcrnaFamily[] families;
		final byte[][] library;
		final BaseGraph[] models;
		final String[] names;
		final int[] owner, local, minima, maxima;
		final float[] pass, borderline;
		final int minimum;
	}

	/** Fail on unsupported geometry rather than borrowing one member's defaults. */
	private static void validateCommon(final NcrnaFamily main, final NcrnaFamily member){
		require(main.kLong==17 && member.kLong==17 && main.kmerSet!=null && member.kmerSet!=null
			&& main.windowPad==member.windowPad && main.collapseFrac==member.collapseFrac
			&& main.seedMinHits==member.seedMinHits && main.seedDistinct==member.seedDistinct
			&& main.scavengePass2==member.scavengePass2 && main.quantumThresh==member.quantumThresh
			&& main.trimAlignmentExtent==member.trimAlignmentExtent
			&& main.hbmPass==member.hbmPass && member.hbmPass>1 && member.maxLen==Integer.MAX_VALUE
			&& member.minLen>60 && member.outputType==ProkObject.RNA && member.indexK==9
			&& !member.voteWindows && !member.voteEnds && !member.voteBeforePadding
			&& !member.preserveContainingWindows && member.collapsedWindowPad==0
			&& !member.reuseConsensusAlignment && !member.pacBioConsensusAlignment && !member.modelEndClipping
			&& !member.modelClipRescue && member.modelClipRescueId==-1 && member.prok16sPartialMinLen==0
			&& member.joinedProposals==null && member.joinedPositions==null && !member.deferBoundaryReport
			&& member.cdsStopConstraintFactory==null && member.rrnaEndpointFeatures==null
			&& !member.rankedModelFallback && !member.strictIndexCutoff && member.indexScoreMargin==-1
			&& (main.physicalEdgeMinLen==0 ? member.physicalEdgeMinLen==0 : member.physicalEdgeMinLen==member.minLen),
			"Joint cobalamin requires common ordinary seed/window rules and supported physical-edge controls");
		final long[] expected=main.kmerSet.toArray(), actual=member.kmerSet.toArray();
		Arrays.sort(expected);
		Arrays.sort(actual);
		require(Arrays.equals(expected, actual), "The shared clan scheduler requires equal actual loaded seed words");
	}

	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private final Models axis;
	private final NcrnaScavenger[] refiners;
	private int[] proposalPads;
	private int[] proposalCeilings;
	private NcrnaOverlapWindowBuilder[] overlapBuilders;
	private NcrnaOverlapWindowBuilder.Caps[] keyCaps;
	private boolean withinMemberUnion;
	private boolean proposalsStarted;
	private WindowGeometry[] memberGeometry;
	private NcrnaWindowCounters queryGeometry;
	private NcrnaPathScoreConfig.Cell pathScoreCell;
	private final boolean memberBoundaryTrace;
	private final NcrnaAlignmentWork[] memberWork={new NcrnaAlignmentWork(), new NcrnaAlignmentWork(), new NcrnaAlignmentWork()};
	private final NcrnaAlignmentWork previous=new NcrnaAlignmentWork();
	private int activeMember=-1, windowMembers;
}
