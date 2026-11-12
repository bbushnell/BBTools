package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceCalls.Call;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Iterates one generic RNA caller, or an explicit complete group, in the whole-input caller.
 * Adapted from template.A_SampleBasic2 for argument/output ownership; CallGenes
 * owns sequence I/O, resources, PGMs and workers. Each policy executes the real
 * caller, including claims, second passes, endpoints and final path selection.
 * Exact identity kernels and same-contig shared seed scans may be reused.
 * Optional candidate reuse also retains unchanged pre-path CDS/non-target RNA.
 *
 * References and query-family projection use the ordinary native codecs and
 * matcher. Results retain fixed reference denominators, zero-call inputs and
 * paired lost/gained reference IDs. No biological admission or cutoff selection
 * occurs here. The first policy is the configured caller baseline.
 *
 * Without candidate reuse, CDS and non-target RNA calling repeat for each
 * policy. Cached sweep time is not production caller time. Legacy target-scavenger counts
 * remain separate from optional actual-kernel CPU/model-base accounting.
 * @author Brian Bushnell, Raiden
 */
public final class SweepNcrnaShortlist implements CallGenesIteration{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final SweepNcrnaShortlist tool=new SweepNcrnaShortlist(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	public SweepNcrnaShortlist(final String[] args){
		final PreParser pp=new PreParser(args, SweepNcrnaShortlist.class, false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("caller")){callerConfig=b;}
			else if(a.equals("reference")){reference=b;}
			else if(a.equals("familymap")){familymap=b;}
			else if(a.equals("target")){targetName=b;}
			else if(a.equals("grading")){gradingView=b;}
			else if(a.equals("policies")){policies=b;}
			else if(a.equals("joint")){joint=Parse.parseBoolean(b);}
			else if(a.equals("cachewindows")){cacheWindows=Parse.parseIntKMG(b);}
			else if(a.equals("cachelength")){cacheLength=Parse.parseIntKMG(b);}
			else if(a.equals("savegff")){saveGff=Parse.parseBoolean(b);}
			else if(a.equals("savestages")){saveStages=Parse.parseBoolean(b);}
			else if(a.equals("reusecandidates")){reuseCandidates=Parse.parseBoolean(b);}
			else if(a.equals("coverage")){coverage=Parse.parseBoolean(b);}
			else if(a.equals("members")){membersPath=b;}
			else if(a.equals("loci")){lociPath=b;}
			else if(a.equals("overlap")){overlap=Float.parseFloat(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		prefix=parser.out1;
		require(callerConfig!=null && reference!=null && familymap!=null && policies!=null && prefix!=null,
			"Supply caller=config reference=manifest familymap=rules policies=table out=prefix");
		require(!parser.overwrite && !parser.append, "Sweep outputs must be fresh");
		require(cacheWindows>=0 && cacheLength>0 && overlap>0 && overlap<=1, "Invalid cache bounds or reciprocal overlap");
		require((membersPath==null)==(lociPath==null), "Member coverage requires both members= and loci= metadata");
		if(membersPath!=null){
			coverage=true;
			require(Tools.testInputFiles(false, true, membersPath, lociPath), "Member coverage metadata must be readable");
			require(Tools.testOutputFiles(false, false, false, prefix+".coverage.tsv"), "Member matched-count output must be fresh");
		}
		require(Tools.testInputFiles(false, true, callerConfig, reference, familymap, policies), "Sweep inputs must be readable");
		require(Tools.testOutputFiles(false, false, false, prefix+".summary.tsv", prefix+".inputs.tsv", prefix+".changes.tsv", prefix+".target.tsv"),
			"Sweep reports must be fresh and writable");
		targetSelection=new NcrnaShortlistTarget(targetName, gradingView, joint);
		if(joint){
			final NcrnaJointPolicy[] points=NcrnaJointPolicy.load(policies);
			policy=new Policy[points.length];
			for(int i=0; i<points.length; i++){policy[i]=new Policy(points[i]);}
		}else{policy=loadPolicies(policies);}
		rules=GradeNcrnaReference.loadRules(familymap);
		for(final Rule rule:rules){
			if(rule.include){
				require(scoringFamily==null || scoringFamily.equals(rule.family), "This sweep grades one explicit family, excluding other generic RNA families");
				scoringFamily=rule.family;
			}
		}
		require(scoringFamily!=null, "At least one query projection must be included");
	}

	/** Load one ordered policy list; row zero is the original caller configuration. */
	private Policy[] loadPolicies(final String path){
		final ArrayList<Policy> list=new ArrayList<Policy>();
		final ObjectIntMap<String> names=new ObjectIntMap<String>(String.class);
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(path)){
			final byte[] header=input.reader.nextLine();
			require(header!=null, "Policy table must have a header: "+path);
			parser.set(header);
			require(parser.terms()==3 && parser.termEquals("policy", 0) && parser.termEquals("topn", 1)
				&& (parser.termEquals("floor", 2) || parser.termEquals("minhits", 2)),
				"Policy header must be policy, topn, and floor (adaptive) or minhits (fixed)");
			adaptiveColumn=parser.termEquals("floor", 2);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==3, "Each policy needs name, maximum model count and its mode-specific absolute cutoff");
				final String name=parser.parseString(0);
				final int cap=parser.parseInt(1);
				final float floor=adaptiveColumn ? parser.parseFloat(2) : 0;
				parser.setBounds(2);
				final int minHits=adaptiveColumn ? -1 : NcrnaShortlistTarget.fixedHits(line, parser.a(), parser.b());
				require(!name.isEmpty() && names.get(name)<0, "Policy names must be nonempty and unique");
				require(cap>0 && Float.isFinite(floor) && floor>=0, "Policy cap must be positive and floor finite/nonnegative");
				names.put(name, list.size());
				list.add(new Policy(name, cap, floor, minHits));
			}
		}
		require(!list.isEmpty(), "A sweep needs at least its baseline policy");
		return list.toArray(new Policy[list.size()]);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		final String[] expanded=Parser.parseConfig(new String[]{"config="+NcrnaCallerProcess.configurationPath(new java.io.File(callerConfig))});
		final ArrayList<String> arguments=new ArrayList<String>();
		for(final String arg:expanded){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			if(saveStages && (a.equals("ncrnadiag") || a.equals("r58diag"))){
				require(split.length>1 && (split[1].isEmpty() || split[1].equalsIgnoreCase("null")),
					"savestages owns per-input/policy diagnostics; omit the caller's shared ncrnadiag/r58diag output");
			}
			require(!a.endsWith("root") && !a.equals("ncrnaconsensus")
				&& !a.equals("ncrnamodels") && !a.equals("ncrnakmers"), "Shortlist iteration requires packaged resource roots and libraries: "+a);
			arguments.add(arg);
		}
		// The sweep owns its reports. Native aggregate output, if left enabled,
		// would look like a baseline while its work counters include every policy.
		arguments.add("out=null");
		arguments.add("outstats=null");
		final CallGenes nativeCaller=new CallGenes(arguments.toArray(new String[arguments.size()]));
		nativeCaller.run(timer, () -> {
			bind(nativeCaller.inputFiles());
			targetIndices=targetSelection.bindTargets(GeneCaller.ncrnaFamilies, rules, scoringFamily);
			if(joint){
				require(policy[0].joint.size()==targetIndices.length, "Joint policies must include every selected group member");
				for(int i=0; i<targetIndices.length; i++){
					policy[0].joint.member(i).validateBaseline(GeneCaller.ncrnaFamilies.get(targetIndices[i]), gradingView);
				}
			}else{
				targetSelection.validate(GeneCaller.ncrnaFamilies.get(targetIndices[0]), adaptiveColumn, policy[0].topN, policy[0].floor, policy[0].minHits);
			}
			writeTarget();
			nativeCaller.setIteration(this);
		});
		writeReports();
		outstream.println("Completed "+results.length+" inputs and "+policy.length+" policies. Cached sweep elapsed: "+timer);
		outstream.println("Legacy target-scavenger counts are retained separately. Kernel CPU/model-base columns require enabled ncRNA CPU profiling; cached CPU is physical work only.");
	}

	/** Bind exactly the configured FASTAs to original IDs; off-panel references
	 * are omitted explicitly and their GFFs are never read. No bare-contig joins.
	 */
	private void bind(final String[] inputs){
		require(inputs.length>0, "The native caller must have at least one input");
		final ObjectMap<String, Binding> references=new ObjectMap<String, Binding>(String.class, Binding.class);
		for(final Binding binding:GradeNcrnaReference.loadManifest(reference)){references.put(binding.fasta, binding);}
		results=new InputResult[inputs.length];
		bindings=new Binding[inputs.length];
		for(int i=0; i<inputs.length; i++){
			final Binding binding=references.get(inputs[i]);
			require(binding!=null && inputIndex.get(inputs[i])<0, "Every exact caller input must bind one unique reference FASTA: "+inputs[i]);
			inputIndex.put(inputs[i], i);
			bindings[i]=binding;
			if(saveGff){for(int p=0; p<policy.length; p++){
				require(Tools.testOutputFiles(false, false, false, gffPath(i, p)), "Bounded GFF evidence must be fresh");
			}}
			if(saveStages){
				for(int p=0; p<policy.length; p++){
					require(Tools.testOutputFiles(false, false, false, stagePath(i, p)), "Bound per-input/policy stage evidence must be fresh");
				}
			}
		}
		if(membersPath!=null){memberSupport=new NcrnaMemberCoverage(membersPath, lociPath, scoringFamily, bindings);}
	}

	/** Legacy tables keep their mode; joint tables apply the complete explicit tuple. */
	private void applyPolicy(final NcrnaScavenger scavenger, final Policy point, final int member){
		if(joint){point.joint.member(member).apply(scavenger);}
		else if(adaptiveColumn){scavenger.setShortlistPolicy(point.topN, point.floor);}
		else{scavenger.setFixedShortlistPolicy(point.topN, point.minHits);}
	}

	@Override
	public Worker create(final String input, final GeneCaller caller, final GeneModel model){
		final int index=inputIndex.get(input);
		require(index>=0 && results[index]==null, "Iteration worker must own one declared, unfinished input: "+input);
		return new PolicyWorker(index, caller, model);
	}

	/** Worker-local references, caller state and bounded cache. Only the completed
	 * scalar result and changed-reference rows survive whole-input processing.
	 */
	private final class PolicyWorker implements Worker{
		PolicyWorker(final int index_, final GeneCaller caller_, final GeneModel model_){
			index=index_;
			caller=caller_;
			model=model_;
			require(targetIndices!=null && targetIndices.length>0, "Target preflight must finish before workers start");
			scavengers=new NcrnaScavenger[targetIndices.length];
			work=new NcrnaAlignmentWork[targetIndices.length];
			caches=new NcrnaAlignmentCache[targetIndices.length];
			for(int i=0; i<targetIndices.length; i++){
				final int target=targetIndices[i];
				final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(target);
				scavengers[i]=caller.ncrnaScavenger(target);
				work[i]=caller.ncrnaAlignmentWork(target);
				caches[i]=new NcrnaAlignmentCache(family.library, cacheWindows, cacheLength);
				scavengers[i].setAlignmentCache(caches[i]);
			}
			for(int p=0; p<policy.length; p++){
				queries[p]=reader.result(bindings[index].id);
				stats[p]=new Stats();
				if(saveGff){gffs[p]=new Output(gffPath(index, p));}
				if(saveStages){
					stages[p]=new Output(stagePath(index, p));
					stages[p].writer.println(NcrnaStageDiagLogSink.HEADER);
					stageSinks[p]=new NcrnaStageDiagLogSink(stages[p].writer);
				}
			}
		}

		@Override
		public ArrayList<Orf> call(final Read read){
			ArrayList<Orf> baseline=null;
			for(final NcrnaAlignmentCache cache:caches){cache.clear();}
			if(reuseCandidates){caller.beginNcrnaIteration(read, targetIndices);}
			else{caller.beginNcrnaIteration(read);}
			try{
				for(int p=0; p<policy.length; p++){
					if(saveStages){caller.setNcrnaDiagSink(stageSinks[p]);}
					for(int i=0; i<scavengers.length; i++){
						applyPolicy(scavengers[i], policy[p], i);
						stats[p].accumulateTarget(scavengers[i], caches[i], work[i], -1);
					}
					final long cpuStart=NcrnaCpuProfile.cpuNow();
					final ArrayList<Orf> path=caller.callGenes(read, model, true);
					if(cpuStart>=0){stats[p].callingCpu+=NcrnaCpuProfile.cpuNow()-cpuStart;}
					for(int i=0; i<scavengers.length; i++){
						stats[p].accumulateTarget(scavengers[i], caches[i], work[i], 1);
					}
					if(p==0){baseline=path;}
					for(final Orf orf:path){
						reader.add(queries[p], orf);
						if(gffs[p]!=null){
							buffer.clear();
							gffs[p].writer.print(orf.appendGff(buffer).nl());
						}
					}
				}
			}finally{
				if(saveStages){caller.setNcrnaDiagSink(null);}
				caller.endNcrnaIteration();
				for(int i=0; i<scavengers.length; i++){
					caches[i].clear();
					applyPolicy(scavengers[i], policy[0], i);
				}
			}
			// A later path reused some baseline candidates. Restore their reporting
			// bounds before handing the first policy to ordinary output consumers.
			if(reuseCandidates){NcrnaReportedOrf.applyTo(baseline);}
			return baseline;
		}

		@Override
		public void finish(){
			if(saveStages){
				try{
					for(final Output output:gffs){if(output!=null){output.close();}}
					for(final Output output:stages){output.close();}
				}
				catch(RuntimeException | Error failure){
					shared.KillSwitch.kill("Per-input/policy stage output failed on caller worker: "+failure);
					throw failure;
				}
			}else{for(final Output output:gffs){if(output!=null){output.close();}}}
			final ArrayList<Row> allRows=new NcrnaReferenceGff().read(bindings[index].gff, bindings[index].id);
			final ArrayList<Row> rows=new ArrayList<Row>();
			final ArrayList<Locus> truth=new ArrayList<Locus>();
			for(final Row row:allRows){if(scoringFamily.equals(row.locus.family)){
				rows.add(row);
				truth.add(row.locus);
			}}
			final boolean[] masks=new boolean[rows.size()];
			for(int i=0; i<rows.size(); i++){masks[i]=rows.get(i).mask;}
			final int[] truthInputs=ones(rows.size());
			final long[] support=memberSupport==null ? null : memberSupport.bind(index, rows);
			final long[] memberRefs=support==null ? null : new long[memberSupport.size()];
			if(support!=null){
				for(final long bits:support){
					for(int m=0; m<memberRefs.length; m++){if((bits&(1L<<m))!=0){memberRefs[m]++;}}
				}
			}
			NcrnaReferenceMatcher.Result baseline=null;
			final ByteBuilder changes=new ByteBuilder();
			for(int p=0; p<policy.length; p++){
				final ArrayList<Locus> calls=new ArrayList<Locus>();
				for(final Call call:queries[p].calls){if(call.included){calls.add(call.locus);}}
				final NcrnaReferenceMatcher.Result match=NcrnaReferenceMatcher.match(calls, ones(calls.size()), truth, truthInputs, masks, overlap);
				if(p==0){baseline=match;}
				stats[p].grade(match, calls, truth);
				if(coverage){
					stats[p].lastLossInputs=stats[0].tp>0 && stats[p].tp==0 ? 1 : 0;
					stats[p].firstGainInputs=stats[0].tp==0 && stats[p].tp>0 ? 1 : 0;
				}
				if(support!=null){
					stats[p].memberReferences=memberRefs;
					stats[p].memberMatched=new long[memberRefs.length];
					for(int i=0; i<support.length; i++){
						if(match.referenceToCall[i]<0 || masks[i]){continue;}
						for(int m=0; m<memberRefs.length; m++){
							if((support[i]&(1L<<m))!=0){stats[p].memberMatched[m]++;}
						}
					}
				}
				stats[p].features=queries[p].features;
				stats[p].cds=queries[p].cds;
				stats[p].trna=queries[p].trna;
				stats[p].rrna=queries[p].rrna;
				stats[p].excluded=queries[p].excluded;
				for(int i=0; i<rows.size(); i++){
					if(masks[i]){continue;}
					final boolean was=baseline.referenceToCall[i]>=0, now=match.referenceToCall[i]>=0;
					if(was==now){continue;}
					if(was){stats[p].lost++;}
					else{stats[p].gained++;}
					final Row row=rows.get(i);
					field(changes, bindings[index].id).tab();
					field(changes, policy[p].name).tab().append(was ? "LOST" : "GAINED").tab();
					field(changes, row.id).tab();
					field(changes, row.locus.seqid).tab().append(row.locus.start)
						.tab().append(row.locus.stop).tab().append(row.locus.strand).nl();
				}
				require(stats[p].tp==stats[0].tp-stats[p].lost+stats[p].gained, "Paired reference changes must explain the net TP change");
			}
			results[index]=new InputResult(stats, changes);
		}

		private final int index;
		private final GeneCaller caller;
		private final GeneModel model;
		private final NcrnaScavenger[] scavengers;
		private final NcrnaAlignmentCache[] caches;
		private final NcrnaAlignmentWork[] work;
		private final NcrnaReferenceCalls reader=new NcrnaReferenceCalls(rules);
		private final NcrnaReferenceCalls.Result[] queries=new NcrnaReferenceCalls.Result[policy.length];
		private final Stats[] stats=new Stats[policy.length];
		private final Output[] gffs=new Output[policy.length];
		/** Allocated only for explicit tracing; one worker owns every bound sink. */
		private final Output[] stages=saveStages ? new Output[policy.length] : null;
		private final NcrnaStageDiagSink[] stageSinks=saveStages ? new NcrnaStageDiagSink[policy.length] : null;
		private final ByteBuilder buffer=new ByteBuilder();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	/** Declares every swept member separately from the reference-grading population.
	 * These rows describe the immutable root; the policy input declares each tuple.
	 */
	private void writeTarget(){
		final ByteBuilder b=new ByteBuilder();
		try(Output output=new Output(prefix+".target.tsv")){
			output.writer.println("target\tgrading\treference_family\tcutoff_mode\tmap_k\ttop_fraction\tquery_fraction\tmin_length\twindow_pad\tmodel_length_bounds");
			for(final int index:targetIndices){
				final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(index);
				b.clear();
				field(b, family.name).tab();
				field(b, gradingView).tab();
				field(b, scoringFamily).tab();
				b.append(family.adaptive ? "adaptive" : "fixed").tab().append(family.indexK).tab()
					.append(family.adaptTopFrac, 8).tab().append(family.adaptQFrac, 8).tab().append(family.minLen)
					.tab().append(family.windowPad).tab().append(family.modelThresholds!=null && family.modelThresholds.hasLengths()).nl();
				output.writer.print(b);
			}
		}
	}

	private String cutoffColumn(){return adaptiveColumn ? "floor" : "minhits";}
	private String policyColumns(){return joint ? "policy" : "policy\ttopn\t"+cutoffColumn();}
	private String metricColumns(){
		final String base=joint ? METRICS+"\tall_features\tCDS\ttRNA\trRNA\texcluded_generic_RNA" : METRICS;
		return coverage ? base+"\tlast_loss_inputs\tfirst_gain_inputs" : base;
	}

	private void writeReports(){
		final Stats[] totals=new Stats[policy.length];
		for(int p=0; p<totals.length; p++){totals[p]=new Stats();}
		final ByteBuilder buffer=new ByteBuilder();
		try(Output inputs=new Output(prefix+".inputs.tsv"); Output changes=new Output(prefix+".changes.tsv")){
			inputs.writer.println("input_id\t"+policyColumns()+"\t"+metricColumns());
			changes.writer.println("input_id\tpolicy\tchange\treference_id\tseqid\tstart\tstop\tstrand");
			for(int i=0; i<results.length; i++){
				require(results[i]!=null, "Every declared input must finish, including zero-call inputs");
				for(int p=0; p<policy.length; p++){
					final Stats stats=results[i].stats[p];
					totals[p].add(stats);
					buffer.clear();
					field(buffer, bindings[i].id).tab();
					row(buffer, p, stats);
					inputs.writer.print(buffer);
				}
				changes.writer.print(results[i].changes);
			}
		}
		try(Output output=new Output(prefix+".summary.tsv")){
			output.writer.println("# Cached sweep; legacy target-scavenger counts are separate from actual kernel counts. CPU/model-base columns are NA without profiling. No production-runtime inference.");
			output.writer.println(policyColumns()+"\t"+metricColumns());
			for(int p=0; p<totals.length; p++){
				buffer.clear();
				row(buffer, p, totals[p]);
				output.writer.print(buffer);
			}
		}
		if(memberSupport!=null){writeMemberCoverage();}
	}

	/** Complete matched-parent counts include unchanged detections and explicit zeroes. */
	private void writeMemberCoverage(){
		final ByteBuilder buffer=new ByteBuilder();
		try(Output output=new Output(prefix+".coverage.tsv")){
			output.writer.println("# Nonadditive qualified truth support under the union matcher; marginal support is excluded.");
			output.writer.println("input_id\tpolicy\tmember\tmatched_scored_parents\tqualified_scored_parents");
			for(int i=0; i<results.length; i++){
				for(int p=0; p<policy.length; p++){
					final Stats s=results[i].stats[p];
					assert(s.memberMatched!=null && s.memberReferences!=null)
						: "Every member-enabled input must finish its live matcher counts before report publication";
					for(int m=0; m<memberSupport.size(); m++){
						field(buffer.clear(), bindings[i].id).tab();
						field(buffer, policy[p].name).tab();
						field(buffer, memberSupport.name(m)).tab().append(s.memberMatched[m]).tab().append(s.memberReferences[m]);
						output.writer.print(buffer.nl());
					}
				}
			}
		}
	}

	private void row(final ByteBuilder buffer, final int p, final Stats s){
		assert(s.tp+s.fn==s.scored) : "Recall denominator must conserve scored reference loci";
		field(buffer, policy[p].name);
		if(!joint){
			buffer.tab().append(policy[p].topN).tab();
			if(adaptiveColumn){buffer.append(policy[p].floor, 6);}
			else{buffer.append(policy[p].minHits);}
		}
		buffer.tab().append(s.tp).tab().append(s.fp).tab().append(s.fn).tab().append(s.scored).tab().append(s.masked)
			.tab().append(s.maskCalls).tab().append(s.duplicates).tab().append(s.lost).tab().append(s.gained)
			.tab().append(s.exact5).tab().append(s.exact3).tab().append(s.exactBoth)
			.tab().append(s.alignments).tab().append(s.fills).tab().append(s.hits).tab();
		ratio(buffer, s.tp, s.tp+s.fp).tab();
		ratio(buffer, s.tp, s.scored).tab();
		ratio(buffer, s.fp, s.tp+s.fp);
		if(NcrnaCpuProfile.enabled()){
			buffer.tab().append(s.kernelCount).tab().append(s.modelBases).tab().append(s.physicalCount)
				.tab().append(s.physicalModelBases).tab().append(s.kernelCpu).tab().append(s.hbmCalls)
				.tab().append(s.hbmCpu).tab().append(s.shortlistCpu).tab().append(s.callingCpu);
		}else{buffer.append("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");}
		if(joint){buffer.tab().append(s.features).tab().append(s.cds).tab().append(s.trna).tab().append(s.rrna).tab().append(s.excluded);}
		if(coverage){buffer.tab().append(s.lastLossInputs).tab().append(s.firstGainInputs);}
		buffer.nl();
	}

	private static ByteBuilder field(final ByteBuilder b, final String text){return BuildNcrnaReference.field(b, text);}
	private String gffPath(final int input, final int p){return prefix+".input"+input+".policy"+p+".gff";}
	private String stagePath(final int input, final int p){return prefix+".input"+input+".policy"+p+".stages.tsv";}
	private static ByteBuilder ratio(final ByteBuilder b, final long n, final long d){return d==0 ? b.append("NA") : b.append(n/(double)d, 8);}
	private static int[] ones(final int length){
		final int[] array=new int[length];
		Arrays.fill(array, 1);
		return array;
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Policy{
		Policy(final String name_, final int topN_, final float floor_, final int minHits_){
			name=name_;
			topN=topN_;
			floor=floor_;
			minHits=minHits_;
			joint=null;
		}
		Policy(final NcrnaJointPolicy joint_){
			joint=joint_;
			name=joint.name;
			topN=-1;
			floor=Float.NaN;
			minHits=-1;
		}
		final String name;
		final NcrnaJointPolicy joint;
		final int topN, minHits;
		final float floor;
	}
	private static final class InputResult{
		InputResult(final Stats[] stats_, final ByteBuilder changes_){
			stats=stats_;
			changes=changes_;
		}
		final Stats[] stats;
		final ByteBuilder changes;
	}
	private static final class Stats{
		/** Sum work over all selected members; cached physical cost stays separate. */
		void accumulateTarget(final NcrnaScavenger scavenger, final NcrnaAlignmentCache cache,
				final NcrnaAlignmentWork work, final int sign){
			assert(sign==-1 || sign==1) : "Before/after snapshots must cancel cumulative counters outside this policy";
			accumulateWork(work, sign);
			shortlistCpu+=sign*scavenger.shortlistCpuNanos();
			alignments+=sign*scavenger.alignmentCount();
			fills+=sign*cache.fills();
			hits+=sign*cache.hits();
		}
		void grade(final NcrnaReferenceMatcher.Result result, final ArrayList<Locus> calls, final ArrayList<Locus> truth){
			tp=result.matched;
			fp=result.falsePositives;
			fn=result.falseNegatives;
			scored=result.scoredReferences;
			masked=result.maskedReferences;
			maskCalls=result.maskedCalls;
			duplicates=result.duplicates;
			for(int i=0; i<calls.size(); i++){
				final int ref=result.callToReference[i];
				if(ref<0){continue;}
				final Locus call=calls.get(i), reference=truth.get(ref);
				final boolean start=call.start==reference.start, stop=call.stop==reference.stop;
				if(call.strand=='+' ? start : stop){exact5++;}
				if(call.strand=='+' ? stop : start){exact3++;}
				if(start && stop){exactBoth++;}
			}
		}
		/** Add/subtract snapshots around one policy without per-contig metric objects. */
		void accumulateWork(final NcrnaAlignmentWork w, final int sign){
			if(w==null){return;}
			assert(sign==1 || sign==-1) : "A policy uses exactly one before and one after snapshot";
			kernelCount+=sign*w.logical;
			modelBases+=sign*w.modelBases;
			physicalCount+=sign*w.physical;
			physicalModelBases+=sign*w.physicalModelBases;
			kernelCpu+=sign*w.cpuNanos;
			hbmCalls+=sign*w.hbmCalls;
			hbmCpu+=sign*w.hbmCpuNanos;
		}
		void add(final Stats s){
			tp+=s.tp;
			fp+=s.fp;
			fn+=s.fn;
			scored+=s.scored;
			masked+=s.masked;
			maskCalls+=s.maskCalls;
			duplicates+=s.duplicates;
			lost+=s.lost;
			gained+=s.gained;
			exact5+=s.exact5;
			exact3+=s.exact3;
			exactBoth+=s.exactBoth;
			alignments+=s.alignments;
			fills+=s.fills;
			hits+=s.hits;
			kernelCount+=s.kernelCount;
			modelBases+=s.modelBases;
			physicalCount+=s.physicalCount;
			physicalModelBases+=s.physicalModelBases;
			kernelCpu+=s.kernelCpu;
			hbmCalls+=s.hbmCalls;
			hbmCpu+=s.hbmCpu;
			shortlistCpu+=s.shortlistCpu;
			callingCpu+=s.callingCpu;
			features+=s.features;
			cds+=s.cds;
			trna+=s.trna;
			rrna+=s.rrna;
			excluded+=s.excluded;
				lastLossInputs+=s.lastLossInputs;
				firstGainInputs+=s.firstGainInputs;
		}
		long tp, fp, fn, scored, masked, maskCalls, duplicates, lost, gained;
		long exact5, exact3, exactBoth, alignments, fills, hits;
		long kernelCount, modelBases, physicalCount, physicalModelBases, kernelCpu, hbmCalls, hbmCpu, shortlistCpu, callingCpu;
		long features, cds, trna, rrna, excluded;
			long lastLossInputs, firstGainInputs;
			long[] memberMatched, memberReferences;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String callerConfig, reference, familymap, policies, prefix, scoringFamily;
	private String targetName="tpp", gradingView="target";
	private final NcrnaShortlistTarget targetSelection;
	private int[] targetIndices;
	private boolean joint;
	private boolean adaptiveColumn;
	private int cacheWindows=128, cacheLength=4096;
	private float overlap=.5f;
	private boolean saveGff, saveStages;
	private final Policy[] policy;
	private final Rule[] rules;
	private final PrintStream outstream;
	private final ObjectIntMap<String> inputIndex=new ObjectIntMap<String>(String.class);
	private Binding[] bindings;
	private boolean reuseCandidates=false;
	private boolean coverage;
	private String membersPath, lociPath;
	private NcrnaMemberCoverage memberSupport;
	private InputResult[] results;
	private static final String METRICS="TP\tFP\tFN\tscored_references\tmasked_references\tmasked_calls\tduplicates\tlost_baseline_TP\tgained_TP\texact_5prime\texact_3prime\texact_both\ttarget_scavenger_alignments\tidentity_cache_fills\tidentity_cache_hits\tprecision\trecall\tFDR\tkernel_logical_alignments\taligned_model_bases\tphysical_kernel_alignments\tphysical_model_bases\tkernel_cpu_ns\thbm_calls\thbm_cpu_ns\tshortlist_cpu_ns\tcalling_cpu_ns";
}
