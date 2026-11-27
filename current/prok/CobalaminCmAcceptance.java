package prok;

import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Worker-owned covariance acceptance after native ends, before claims/path.
 * Capture admits the declared lower-floor candidates without a score cutoff;
 * filter applies an explicit, externally calibrated three-member policy.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminCmAcceptance implements NcrnaCandidateAcceptance{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public CobalaminCmAcceptance(final Resources resources_, final NcrnaClanScavenger clan, final int ordinal_){
		require(resources_!=null && clan!=null && ordinal_>=0, "Acceptance must bind one input's real clan and immutable resources");
		resources=resources_;
		targets=resources.clanOr ? new long[9*TARGET_FIELDS] : null;
		ordinal=ordinal_;
		axis=clan;
		alignments=resources.aligned ? axis.enableAcceptanceAlignments() : null;
		scorers[0]=clan.cmEndpointRefiner(0);
		require(scorers[0]!=null, "The main acceptance route reuses the actual endpoint scorer, never a second search worker");
		for(int member=1; !resources.clanOr && member<3; member++){
			scorers[member]=(CobalaminCmEndpointRefiner)resources.factories[member].create(resources.families[member]);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public float floor(){return resources.floor;}

	@Override
	public boolean accept(final Orf candidate, final byte[] bases, final int model, final float identity){
		require(candidate!=null && bases!=null && Float.isFinite(identity) && identity>=0 && identity<=1
			&& identity>=resources.floor
			&& candidate.start>=0 && candidate.stop>=candidate.start && candidate.stop<bases.length,
			"Acceptance consumes one real sense-oriented candidate and its actual verifier identity");
		final int member=axis.member(model), local=axis.localModel(model);
		final CobalaminCmEndpointRefiner scorer=scorers[member];
		final boolean identityOnly=!resources.capture && identity>=resources.upper[member];
		final long decisionBegin=NcrnaCpuProfile.cpuNow();
		float score=Float.NaN;
		boolean reused=false;
		boolean continued=false;
		long cpu=0;
		long ownAnchor=0, ownCm=0, ownFills=0, ownCpu=0;
		int requests=0, searchStart=-1, searchStop=-1;
		if(!identityOnly){
			if(resources.clanOr && (member>0 || !scorer.applies(local))){
				// A declared unmapped suffix keeps native ends but has no own CM score.
				// Other members may still supply a finite OR pass.
				final CobalaminCrossMemberCmScorer.Worker worker=crossWorker();
				score=worker.scoreExact(member, local, candidate, bases);
				ownCm=worker.cmCpu();
				ownAnchor=worker.anchorCpu();
				ownFills=worker.fills();
				ownCpu=worker.totalCpu();
				extraCmCpu[member]+=ownCm;
				extraAnchorCpu[member]+=ownAnchor;
				extraFills[member]+=ownFills;
				searchStart=worker.windowStart();
				searchStop=worker.windowStop();
				requests++;
			}else{
				reused=scorer.hasResult(candidate, local);
				if(reused){score=scorer.scoreFor(candidate, local);}
				else{
					final long begin=NcrnaCpuProfile.cpuNow();
					final long cmBefore=scorer.searchCpu(), fillsBefore=scorer.searchCalls();
					final long anchorBefore=scorer.count(CobalaminCmEndpointRefiner.ANCHOR_CPU);
					score=scorer.scoreOnly(candidate, bases, local);
					cpu=elapsed(begin);
					ownCm=scorer.searchCpu()-cmBefore;
					ownAnchor=scorer.count(CobalaminCmEndpointRefiner.ANCHOR_CPU)-anchorBefore;
					ownFills=scorer.searchCalls()-fillsBefore;
					ownCpu=cpu;
					extraCmCpu[member]+=ownCm;
					extraAnchorCpu[member]+=ownAnchor;
					extraFills[member]+=ownFills;
					requests++;
				}
				searchStart=scorer.windowStartFor(candidate, local);
				searchStop=scorer.windowStopFor(candidate, local);
			}
		}
		final float ga=resources.clanOr ? resources.cross.ga(member) : scorer.gatheringThreshold();
		require(Float.isFinite(ga) && (Float.isNaN(score) || Float.isFinite(score) || score==Float.NEGATIVE_INFINITY),
			"A missing score stays absent; malformed CM arithmetic or GA cannot become an ordinary acceptance");
		if(resources.clanOr && !identityOnly){observe(member, member, score, reused, ownFills, 0, ownAnchor, ownCm, ownCpu);}
		boolean accepted=resources.clanOr ? identityOnly || passes(score, ga, resources.margin[member])
			: resources.capture || identityOnly || (Float.isNaN(score)
				? identity>=axis.nativeIdentityPass(model) : (double)score-ga>=resources.margin[member]);
		final int crossTargets=resources.crossTargets[member];
		if(resources.aligned && !accepted && identity>=axis.nativeIdentityPass(model) && alignments.active()){
			require(alignments!=null && alignments.active(),
				"Aligned fallback consumes the still-active ordinary PRIMARY query");
			continued=true;
			final int origin=alignments.senseOrigin();
			require((long)origin+alignments.query().length<=bases.length,
				"The original sense-oriented query must remain within this unchanged contig");
			final CobalaminCrossMemberCmScorer.Worker worker=crossWorker();
			int tried=1<<member;
			for(int phase=0; phase<2 && !accepted; phase++){
				final int limit=phase==0 ? alignments.rankedCount() : Math.min(resources.depth, alignments.rankedCount());
				for(int rank=0; rank<limit && !accepted; rank++){
					final int targetModel=alignments.rankedModel(rank);
					final int target=axis.member(targetModel);
					if((tried&(1<<target))!=0){continue;}
					int state=alignments.state(targetModel);
					if(phase==1 && state==NcrnaAcceptanceAlignments.UNVISITED){
						axis.measureAcceptanceModel(targetModel);
						state=alignments.state(targetModel);
					}
					if(state!=NcrnaAcceptanceAlignments.ALIGNED_VALID
						|| alignments.identity(targetModel)<axis.nativeIdentityPass(targetModel)){continue;}
					// Mark before scoring: even an unavailable or impossible result
					// consumes this target's one CM request for this candidate only.
					tried|=1<<target;
					final float other=worker.scoreAnchored(target, axis.localModel(targetModel), candidate, bases,
						origin+alignments.localStart(targetModel), origin+alignments.localStop(targetModel));
					extraCmCpu[member]+=worker.cmCpu();
					extraAnchorCpu[member]+=worker.anchorCpu();
					extraFills[member]+=worker.fills();
					observe(member, target, other, false, worker.fills(), worker.rankingCpu(), worker.anchorCpu(), worker.cmCpu(), worker.totalCpu());
					requests++;
					accepted=passes(other, worker.ga(), resources.margin[target]);
				}
			}
		}else if(resources.clanOr && !resources.aligned && !accepted && crossTargets!=0){
			final CobalaminCrossMemberCmScorer.Worker worker=crossWorker();
			worker.prepare(candidate, bases);
			try{
				for(int position=0; position<3 && !accepted; position++){
					final int target=resources.priority[position];
					if((crossTargets&(1<<target))==0){continue;}
					final float other=worker.scoreChosen(target, candidate, bases);
					extraCmCpu[member]+=worker.cmCpu();
					extraAnchorCpu[member]+=worker.anchorCpu();
					extraFills[member]+=worker.fills();
					observe(member, target, other, false, worker.fills(), worker.rankingCpu(), worker.anchorCpu(), worker.cmCpu(), worker.totalCpu());
					requests++;
					accepted=passes(other, worker.ga(), resources.margin[target]);
				}
			}finally{worker.clearPrepared();}
		}
		if(resources.clanOr && (requests>0 || continued)){
			cpu=elapsed(decisionBegin);
		}
		candidates[member]++;
		if(accepted){kept[member]++;}
		if(!Float.isNaN(score)){measured[member]++;}
		if(reused){reuses[member]++;}
		extra[member]+=requests;
		extraCpu[member]+=cpu;
		if(resources.output!=null){
			write(candidate, bases.length, member, local, identity, score, ga, accepted,
				identityOnly ? "IDENTITY_ONLY" : reused ? "REUSED_ENDPOINT" : "SCORE_ONLY",
				searchStart, searchStop, axis.candidateTrace(model));
		}
		return accepted;
	}
	private CobalaminCrossMemberCmScorer.Worker crossWorker(){
		assert(resources.clanOr && resources.cross!=null) : "Only explicit OR admission uses the shared cross-member worker";
		if(cross==null){cross=resources.cross.currentWorker();}
		return cross;
	}
	private static boolean passes(final float score, final float ga, final double cutoff){
		return Float.isFinite(score) && (double)score-ga>=cutoff;
	}
	private void observe(final int emitted, final int target, final float score, final boolean reused,
			final long fills, final long ranking, final long anchor, final long cm, final long total){
		assert(targets!=null && emitted>=0 && emitted<3 && target>=0 && target<3) :
			"OR target counters retain the actual emitted/scoring member axes for one input";
		final int offset=(emitted*3+target)*TARGET_FIELDS;
		targets[offset]++;
		if(reused){targets[offset+1]++;}
		if(Float.isNaN(score)){targets[offset+2]++;}
		targets[offset+3]+=fills;
		targets[offset+4]+=ranking;
		targets[offset+5]+=anchor;
		targets[offset+6]+=cm;
		targets[offset+7]+=total;
	}

	/** Record genomic coordinates before the native caller flips its path nodes.
	 * These are candidate events, not final selected loci or one-to-one grades. */
	private void write(final Orf candidate, final int length, final int member, final int local,
			final float identity, final float score, final float ga, final boolean accepted, final String source,
			final int searchStart, final int searchStop, final NcrnaBoundaryTrace.Family trace){
		assert(length>0 && member>=0 && member<3) : "Acceptance records bind the same three-member model axis as the caller";
		final boolean minus=candidate.strand==Shared.MINUS;
		final int start=minus ? length-1-candidate.stop : candidate.start;
		final int stop=minus ? length-1-candidate.start : candidate.stop;
		row.clear().append(ordinal).tab();
		for(int i=0; i<candidate.scafName.length(); i++){
			final char c=candidate.scafName.charAt(i);
			if(c==' ' || c=='\t'){break;}
			row.append(c);
		}
		row.tab().append(minus ? '-' : '+').tab().append(MEMBERS[member]).tab()
			.append(resources.names[member][local]).tab()
			.append(length).tab().append(start).tab().append(stop).tab().append(Float.floatToIntBits(identity)).tab();
		if(!Float.isNaN(score)){row.append(Float.floatToIntBits(score));}
		row.tab().append(Float.floatToIntBits(ga)).tab().append(source).tab().append(accepted).tab();
		if(searchStart>=0){row.append(minus ? length-1-searchStop : searchStart);}
		row.tab();
		if(searchStop>=0){row.append(minus ? length-1-searchStart : searchStop);}
		row.tab();
		if(trace!=null){
			require(trace.inputOrdinal()==ordinal, "Acceptance and native trace must share the actual input ordinal, not merely a matching contig name");
			row.append(trace.workerId());
		}
		row.tab();
		if(trace!=null){row.append(trace.candidateId(candidate));}
		row.nl();
		synchronized(resources.output){resources.output.writer.print(row);}
	}

	/** Append scalar costs after this input's caller has finished. */
	public void finish(final Output output){
		require(output!=null, "Acceptance cost rows belong to the driver's existing joined output owner");
		for(int member=0; member<3; member++){
			output.writer.print(row.clear().append(ordinal).tab().append(MEMBERS[member]).tab()
				.append(candidates[member]).tab().append(kept[member]).tab().append(measured[member]).tab()
				.append(reuses[member]).tab().append(extra[member]).tab().append(extraFills[member]).tab()
				.append(extraCpu[member]).tab().append(extraCmCpu[member]).tab().append(extraAnchorCpu[member]).nl());
		}
	}
	/** Separate target accounting preserves the existing eleven-column receipt.
	 * An unavailable attempted target counts as evaluated even with zero DP fills. */
	public void finishTargets(final Output output){
		require(resources.clanOr && targets!=null && output!=null,
			"Target accounting belongs to explicit OR admission after the input's candidate worker finishes");
		for(int emitted=0; emitted<3; emitted++){
			for(int target=0; target<3; target++){
				row.clear().append(ordinal).tab().append(MEMBERS[emitted]).tab().append(MEMBERS[target]);
				final int offset=(emitted*3+target)*TARGET_FIELDS;
				for(int field=0; field<TARGET_FIELDS; field++){row.tab().append(targets[offset+field]);}
				output.writer.print(row.nl());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Shared Resources       ----------------*/
	/*--------------------------------------------------------------*/

	/** CM/map loading is once per batch, on the existing family consensus axes.
	 * Only minor score factories are loaded; their native tip factories stay intact.
	 * Filter output is optional; capture requires an output for its measured events. */
	public static final class Resources{
		public Resources(final NcrnaFamily[] families_, final String root, final String mode,
				final float floor_, final String policy, final Output output_){
			this(families_, root, mode, floor_, policy, output_, null);
		}
		public Resources(final NcrnaFamily[] families_, final String root, final String mode,
				final float floor_, final String policy, final Output output_, final NcrnaDdlProfileResource memberProfiles){
			this(families_, root, mode, floor_, policy, output_, memberProfiles, null);
		}
		public Resources(final NcrnaFamily[] families_, final String root, final String mode,
				final float floor_, final String policy, final Output output_, final NcrnaDdlProfileResource memberProfiles,
				final String targetOrder){
			this(families_, root, mode, floor_, policy, output_, memberProfiles, targetOrder, null);
		}
		public Resources(final NcrnaFamily[] families_, final String root, final String mode,
				final float floor_, final String policy, final Output output_, final NcrnaDdlProfileResource memberProfiles,
				final String targetOrder, final String directions){
			this(families_, root, mode, floor_, policy, output_, memberProfiles, targetOrder, directions, -1);
		}
		public Resources(final NcrnaFamily[] families_, final String root, final String mode,
				final float floor_, final String policy, final Output output_, final NcrnaDdlProfileResource memberProfiles,
				final String targetOrder, final String directions, final int acceptanceDepth){
			this(families_, root, mode, floor_, policy, output_, memberProfiles, targetOrder, directions, acceptanceDepth, false);
		}
		/** Ordinary serving consumes the selected band directly, without a legacy-table override. */
		public static Resources serving(final NcrnaFamily[] families, final String root, final String policy,
				final NcrnaDdlProfileResource profiles, final int depth){
			return new Resources(families, root, "filter_aligned", .60f, policy, null, profiles, null, null, depth, true);
		}
		private Resources(final NcrnaFamily[] families_, final String root, final String mode,
				final float floor_, final String policy, final Output output_, final NcrnaDdlProfileResource memberProfiles,
				final String targetOrder, final String directions, final int acceptanceDepth, final boolean exactBand){
			require(families_!=null && families_.length==3 && root!=null
				&& ("capture".equals(mode) || "filter".equals(mode) || "filter_or".equals(mode) || "filter_aligned".equals(mode))
				&& (!"capture".equals(mode) || output_!=null)
				&& Float.isFinite(floor_) && floor_>=0 && floor_<=1,
				"Acceptance requires three loaded families, an explicit score-resource root/mode and a finite floor; zero retains native admission");
			families=families_.clone();
			floor=floor_;
			capture=mode.equals("capture");
			aligned=mode.equals("filter_aligned");
			clanOr=mode.equals("filter_or") || aligned;
			depth=acceptanceDepth;
			validateAlignedOptions(mode, floor, depth, targetOrder, directions);
			require(targetOrder==null || clanOr && !aligned, "Explicit cross-target priority belongs only to filter_or; aligned mode follows consensus ranks");
			priority=targetPriority(targetOrder);
			require(directions==null || clanOr, "Explicit cross directions belong only to filter_or; own scoring is never pruned");
			final int allowed=CobalaminCrossDirectionMasks.parse(directions);
			for(int member=0; member<3; member++){
				crossTargets[member]=CobalaminCrossDirectionMasks.targets(allowed, member);
			}
			require(clanOr==(memberProfiles!=null), "Explicit cross-member modes require the caller's loaded member-only index; legacy modes do not consume it");
			output=output_;
			require(capture ? policy==null : policy!=null, "Capture fits no cutoff; filter requires an explicit three-member calibration table");
			if(policy!=null){
				if(clanOr){loadOrPolicy(policy);}
				else{loadPolicy(policy);}
			}
			if(aligned){
				for(int member=0; member<3; member++){
					require(upper[member]==.75f || !exactBand && upper[member]==.778f,
						"Packaged acceptance requires an explicit .75 policy; only development callers may apply the historical .778-table bypass");
					upper[member]=.75f;
				}
			}
			final FamilyResourceResolver resolver=new FamilyResourceResolver();
			final FamilyResourceResolver.Overrides overrides=new FamilyResourceResolver.Overrides();
			overrides.root=root;
			for(int member=0; member<3; member++){
				require(families[member]!=null && families[member].name.equals(FAMILIES[member]),
					"Covariance acceptance binds ordered main/second/third families, never a guessed model prefix");
				names[member]=new String[families[member].modelNames.length];
				for(int model=0; model<names[member].length; model++){
					names[member][model]=NcrnaSupport.recordName(families[member].modelNames[model]);
				}
				if(member>0){factories[member]=CobalaminCmEndpointRefiner.load(families[member], resolver, overrides,
					CobalaminCmEndpointRefiner.RETRY_64, false);}
			}
			if(clanOr){
				factories[0]=families[0].boundaryRefinerFactory;
				cross=new CobalaminCrossMemberCmScorer(families, factories, memberProfiles);
			}else{cross=null;}
			if(output!=null){output.writer.println(HEADER);}
		}
		private void loadOrPolicy(final String path){
			final LineParser1 parser=new LineParser1('\t');
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, OR_POLICY_HEADER, path);
				for(int member=0; member<3; member++){
					final byte[] line=input.reader.nextLine();
					require(line!=null, "OR calibration needs one finite rule for every target CM");
					parser.set(line);
					require(parser.terms()==6 && parser.termEquals("CLAN_OR_ANY_FINITE_PASS", 0)
						&& parser.termEquals(MEMBERS[member], 1) && parser.termEquals("FAIL_THIS_CM", 5),
						"The explicit OR schema fixes ordered CM targets and unavailable-score rejection");
					final float lower=parser.parseFloat(2);
					upper[member]=parser.parseFloat(3);
					margin[member]=parser.parseDouble(4);
					require(Float.isFinite(lower) && lower==floor && Float.isFinite(upper[member])
						&& upper[member]>lower && upper[member]<=1 && Double.isFinite(margin[member])
						&& (member==0 || upper[member]==upper[0]),
						"OR cutoffs must match the explicit capture floor, one shared upper identity, and finite score-minus-own-GA margins");
				}
				require(input.reader.nextLine()==null, "Unexpected OR target rows cannot redefine the calibrated three-member rule");
			}
		}
		private void loadPolicy(final String path){
			final LineParser1 parser=new LineParser1('\t');
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, "member\tupper_identity\tcm_minus_ga", path);
				for(int member=0; member<3; member++){
					final byte[] line=input.reader.nextLine();
					require(line!=null, "Acceptance calibration must contain every member, including family fallback rows");
					parser.set(line);
					require(parser.terms()==3 && parser.termEquals(MEMBERS[member], 0), "Calibration rows must retain the complete ordered three-member axis");
					upper[member]=parser.parseFloat(1);
					margin[member]=parser.parseDouble(2);
					require(Float.isFinite(upper[member]) && upper[member]>=floor && upper[member]<=1
						&& Double.isFinite(margin[member]), "Calibrated identity bounds and our own score margins must be finite and ordered");
				}
				require(input.reader.nextLine()==null, "Unexpected calibration rows must not be silently ignored");
			}
		}
		final NcrnaFamily[] families;
		final float floor;
		final boolean capture;
		final boolean clanOr;
		final boolean aligned;
		final int depth;
		final int[] priority;
		final int[] crossTargets=new int[3];
		final CobalaminCrossMemberCmScorer cross;
		final String[][] names=new String[3][];
		final float[] upper=new float[3];
		final double[] margin=new double[3];
		final NcrnaBoundaryRefiner.Factory[] factories=new NcrnaBoundaryRefiner.Factory[3];
		final Output output;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Aligned mode is an explicit depth experiment with frozen .60 admission.
	 * All cross directions are available; target priority comes from real ranks.
	 * Other modes keep their prior options and reject an inapplicable depth.
	 */
	public static void validateAlignedOptions(final String mode, final float floor, final int depth,
			final String targetOrder, final String directions){
		if("filter_aligned".equals(mode)){
			require(floor==.60f && depth>=1 && depth<=100 && targetOrder==null
				&& CobalaminCrossDirectionMasks.parse(directions)==63,
				"filter_aligned requires floor0.60, explicit depth1..100, all cross directions and no acceptanceorder");
		}else{require(depth==-1, "acceptancedepth belongs only to explicit filter_aligned mode");}
	}

	/** Parses only the cross-target priority. Acceptance always tests the emitted
	 * own model first, then skips that member in this complete permutation. */
	public static int[] targetPriority(final String order){
		if(order==null){return new int[]{0, 1, 2};}
		final String[] names=order.split(",", -1);
		require(names.length==3, "acceptanceorder requires each of the three explicit RFAM member IDs exactly once");
		final int[] priority=new int[3];
		int seen=0;
		for(int position=0; position<3; position++){
			int member=-1;
			for(int m=0; m<3; m++){if(MEMBERS[m].equals(names[position])){member=m;}}
			require(member>=0 && (seen&(1<<member))==0, "acceptanceorder cannot omit, duplicate or guess a target member: "+order);
			priority[position]=member;
			seen|=1<<member;
		}
		assert(seen==7) : "A complete priority permutation prevents an OR target from disappearing from acceptance";
		return priority;
	}

	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	/** Timing is optional observation; absence of a clock cannot change a decision. */
	private static long elapsed(final long begin){
		if(begin<0){return 0;}
		final long cpu=NcrnaCpuProfile.cpuNow()-begin;
		require(cpu>=0, "An enabled current-thread CPU clock must stay monotonic across acceptance work");
		return cpu;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Resources resources;
	private final NcrnaClanScavenger axis;
	private final NcrnaAcceptanceAlignments alignments;
	private final int ordinal;
	private final CobalaminCmEndpointRefiner[] scorers=new CobalaminCmEndpointRefiner[3];
	private CobalaminCrossMemberCmScorer.Worker cross;
	private final long[] candidates=new long[3], kept=new long[3], measured=new long[3], reuses=new long[3], extra=new long[3], extraCpu=new long[3];
	private final long[] extraCmCpu=new long[3], extraAnchorCpu=new long[3], extraFills=new long[3];
	private final ByteBuilder row=new ByteBuilder();
	private final long[] targets;
	private static final int TARGET_FIELDS=8;
	public static final String HEADER="input_ordinal0\tcontig\tstrand\tmember\tmodel\tcontig_length\tstart0\tstop0\tidentity_bits\tcm_score_bits\tga_bits\tscore_source\taccepted\tcm_window_start0\tcm_window_stop0\ttrace_worker\ttrace_id";
	public static final String COST_HEADER="input_ordinal0\tmember\tcandidates\tkept\tmeasured\treused_endpoint_results\textra_score_requests\textra_dp_fills\textra_acceptance_cpu_ns\textra_cm_cpu_ns\textra_anchor_cpu_ns";
	public static final String OR_POLICY_HEADER="acceptance_policy\tmember\tlower_identity\tupper_identity\tcm_minus_ga\tunavailable_score";
	public static final String TARGET_COST_HEADER="input_ordinal0\temitted_member\tscoring_member\tevaluated_events\treused\tunavailable\textra_dp_fills\textra_ranking_cpu_ns\textra_anchor_cpu_ns\textra_cm_cpu_ns\textra_total_cpu_ns";
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] FAMILIES={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
}
