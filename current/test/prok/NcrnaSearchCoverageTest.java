package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaJointPolicy;
import prok.NcrnaMemberCoverage;
import prok.NcrnaPolicySearch.Measurement;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaSearchCoverage;
import prok.NcrnaSearchCoverage.Counts;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.OptimizeNcrnaShortlist;
import shared.Shared;
import structures.ByteBuilder;

/** Literal opt-in coverage fixtures. Adapts A_SampleBasic2's CLI/output ownership;
 * all files are tiny artificial metadata. No real caller, aligner or genome runs.
 * The injected evaluator is explicitly synthetic, not biological verification.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSearchCoverageTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSearchCoverageTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown test parameter: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Supply a fresh out=fixture_directory");
			final NcrnaSearchCoverageTest test=new NcrnaSearchCoverageTest(new File(parser.out1).getAbsoluteFile());
			test.metadata();
			test.support();
			test.driver();
			pp.outstream.println("NCRNA_SEARCH_COVERAGE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaSearchCoverageTest(final File root_){
		root=root_;
		check(root.mkdir(), "Coverage fixtures require a fresh directory under an existing parent");
		write("a.fa", ">a\nACGT\n");
		write("b.fa", ">b\nACGT\n");
		caller=new String[]{"in="+path("a.fa"), "in="+path("b.fa")};
		write("caller.config", caller[0]+"\n"+caller[1]+"\n");
		write("reference.tsv", "input_id\tfasta\tgff\na\t"+path("a.fa")+"\t"+path("a.gff")
			+"\nb\t"+path("b.fa")+"\t"+path("b.gff")
			+"\noff_panel\t"+path("absent.fa")+"\t"+path("absent.gff")+"\n");
		write("members.tsv", "family\tmember\treference_model\nTPP\tmajor\tRF00001\nTPP\tminor%α\tRF00002\n");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Metadata Tests       ----------------*/
	/*--------------------------------------------------------------*/

	private void metadata(){
		final String prefix=path("literal");
		counts(prefix, false);
		final NcrnaSearchCoverage coverage=new NcrnaSearchCoverage(caller, path("reference.tsv"), new String[]{"major", "minor%α"});
		coverage.bindBaseline(path("baseline.tsv"), false);
		final Counts[] result=coverage.read(prefix, POLICIES, measurements());
		check(result[1].lastLoss()==1 && result[1].firstGain()==1,
			"One lost input and one first gain do not cancel even with unchanged total TP");
		check(result[1].memberLastLoss(0)==1 && result[1].memberLastLoss(1)==0
			&& result[1].memberFirstGain(1)==1, "Member counters remain separate and use unchanged qualified denominators");
		coverage.read(prefix, POLICIES, measurements());
		final NcrnaSearchCoverage resumed=new NcrnaSearchCoverage(caller, path("reference.tsv"), new String[]{"major", "minor%α"});
		resumed.bindBaseline(path("baseline.tsv"), true);
		resumed.read(prefix, POLICIES, measurements());
		check(new File(path("baseline.tsv")).isFile(), "The immutable shipped vector must survive driver reconstruction");
		// Same detection evidence, deliberately different cost units in the two columns.
		final String weighted=path("modelbases");
		write(weighted+".inputs.tsv", read(prefix+".inputs.tsv")
			.replace("kernel_logical_alignments", "kernel_logical_alignments\taligned_model_bases")
			.replace("\t10\n", "\t10\t1000\n").replace("\t5\n", "\t5\t800\n"));
		write(weighted+".summary.tsv", read(prefix+".summary.tsv")
			.replace("kernel_logical_alignments", "kernel_logical_alignments\taligned_model_bases")
			.replace("\t20\t0\t0\n", "\t20\t2000\t0\t0\n")
			.replace("\t10\t1\t1\n", "\t10\t1600\t1\t1\n"));
		final Measurement[] weightedMeasurements={new Measurement(1, 0, 2000, 0), new Measurement(1, 0, 1600, 1)};
		check(new NcrnaSearchCoverage(caller, path("reference.tsv"), new String[0], true)
			.read(weighted, POLICIES, weightedMeasurements)[1].lastLoss()==1,
			"Model-base objective validation must conserve its own units while retaining coverage counts");
		rejected(() -> new NcrnaSearchCoverage(caller, path("reference.tsv"), new String[0], true)
			.read(weighted, POLICIES, measurements()), "Alignment-count values cannot masquerade as model-base objective measurements");

		final String hidden=path("hidden");
		counts(hidden, true);
		final Counts[] memberOnly=new NcrnaSearchCoverage(caller, path("reference.tsv"), new String[]{"major", "minor%α"})
			.read(hidden, POLICIES, measurements());
		check(memberOnly[1].lastLoss()==0 && memberOnly[1].memberLastLoss(0)==1,
			"A different member can preserve union coverage while the original member disappears");
		multicopy();

		final String original=read(prefix+".inputs.tsv");
		write(prefix+".inputs.tsv", original+original.substring(original.lastIndexOf("b\tselected")));
		rejected(() -> coverage.read(prefix, POLICIES, measurements()), "An extra duplicate input/policy cannot be ignored");
		write(prefix+".inputs.tsv", original.replace("b\tselected\t1\t0\t1\t2\t0\t0\t1\t5\n", ""));
		rejected(() -> coverage.read(prefix, POLICIES, measurements()), "Missing zero/positive input-policy accounting must fail");
		write(prefix+".inputs.tsv", original.replace("a\tselected\t0", "unknown\tselected\t0"));
		rejected(() -> coverage.read(prefix, POLICIES, measurements()), "Unknown input cannot acquire a count");
		write(prefix+".inputs.tsv", original.replace("a\tselected\t0\t0\t2", "a\tselected\t-1\t0\t3"));
		rejected(() -> coverage.read(prefix, POLICIES, measurements()), "Negative matched counts must fail");
		write(prefix+".inputs.tsv", original.replace("a\tselected\t0", "a\tselected\t18446744073709551616"));
		rejected(() -> coverage.read(prefix, POLICIES, measurements()), "Decimal overflow cannot wrap into a plausible zero count");
		write(prefix+".inputs.tsv", original);
		final String summary=read(prefix+".summary.tsv");
		write(prefix+".summary.tsv", summary.replace("selected\t1\t0\t3", "selected\t2\t0\t2"));
		rejected(() -> coverage.read(prefix, POLICIES, measurements()), "Equal denominators do not excuse inconsistent summary totals");
		write(prefix+".summary.tsv", summary);
		final String member=read(prefix+".coverage.tsv");
		write(prefix+".coverage.tsv", member.replace("a\tselected\tmajor\t0\t1", "a\tselected\tmajor\t1\t1"));
		rejected(() -> coverage.read(prefix, POLICIES, measurements()), "A member cannot report more matches than its zero union count");
		write(prefix+".coverage.tsv", member);

		// Swap shipped detections between original IDs while retaining the aggregate.
		write(prefix+".inputs.tsv", original.replace("a\tshipped\t1\t0\t1", "a\tshipped\t0\t0\t2")
			.replace("b\tshipped\t0\t0\t2", "b\tshipped\t1\t0\t1"));
		rejected(() -> resumed.read(prefix, POLICIES, measurements()), "Persisted per-input baseline detects same-total coverage drift");
		write(prefix+".inputs.tsv", original);
	}

	private void multicopy(){
		final String prefix=path("multicopy");
		final String input="input_id\tpolicy\t"+METRICS+"\na\tshipped\t3\t0\t1\t4\t0\t0\t0\t10\n"
			+"a\tselected\t2\t0\t2\t4\t0\t1\t0\t5\nb\tshipped\t0\t0\t2\t2\t0\t0\t0\t10\n"
			+"b\tselected\t0\t0\t2\t2\t0\t0\t0\t5\n";
		write(prefix+".inputs.tsv", input);
		write(prefix+".summary.tsv", "policy\t"+METRICS+"\tlast_loss_inputs\tfirst_gain_inputs\n"
			+"shipped\t3\t0\t3\t6\t0\t0\t0\t20\t0\t0\nselected\t2\t0\t4\t6\t0\t1\t0\t10\t0\t0\n");
		check(new NcrnaSearchCoverage(caller, path("reference.tsv"), new String[0]).read(prefix, POLICIES,
			new Measurement[]{new Measurement(3, 0, 20, 0), new Measurement(2, 0, 10, 1)})[1].lastLoss()==0,
			"A three-to-two input loses a locus without losing coverage");
		write(prefix+".inputs.tsv", input.replace("a\tselected\t2\t0\t2\t4\t0\t1", "a\tselected\t0\t0\t4\t4\t0\t3"));
		write(prefix+".summary.tsv", "policy\t"+METRICS+"\tlast_loss_inputs\tfirst_gain_inputs\n"
			+"shipped\t3\t0\t3\t6\t0\t0\t0\t20\t0\t0\nselected\t0\t0\t6\t6\t0\t3\t0\t10\t1\t0\n");
		check(new NcrnaSearchCoverage(caller, path("reference.tsv"), new String[0]).read(prefix, POLICIES,
			new Measurement[]{new Measurement(3, 0, 20, 0), new Measurement(0, 0, 10, 3)})[1].lastLoss()==1,
			"A three-to-zero input spends one coverage exception, not three");
	}

	/** Every expected table cell is literal, not derived by the implementation. */
	private void counts(final String prefix, final boolean hidden){
		final String input=hidden ? "a\tshipped\t1\t0\t1\t2\t0\t0\t0\t10\n"
			+"a\tselected\t1\t0\t1\t2\t0\t1\t1\t5\n"
			+"b\tshipped\t0\t0\t2\t2\t0\t0\t0\t10\n"
			+"b\tselected\t0\t0\t2\t2\t0\t0\t0\t5\n"
			: "a\tshipped\t1\t0\t1\t2\t0\t0\t0\t10\n"
			+"a\tselected\t0\t0\t2\t2\t0\t1\t0\t5\n"
			+"b\tshipped\t0\t0\t2\t2\t0\t0\t0\t10\n"
			+"b\tselected\t1\t0\t1\t2\t0\t0\t1\t5\n";
		write(prefix+".inputs.tsv", "input_id\tpolicy\t"+METRICS+"\n"+input);
		write(prefix+".summary.tsv", "policy\t"+METRICS+"\tlast_loss_inputs\tfirst_gain_inputs\n"
			+"shipped\t1\t0\t3\t4\t0\t0\t0\t20\t0\t0\n"
			+"selected\t1\t0\t3\t4\t0\t1\t1\t10\t"+(hidden ? "0\t0\n" : "1\t1\n"));
		final String member="a\tshipped\tmajor\t1\t1\na\tshipped\tminor%25α\t0\t1\n"
			+"a\tselected\tmajor\t0\t1\na\tselected\tminor%25α\t"+(hidden ? 1 : 0)+"\t1\n"
			+"b\tshipped\tmajor\t0\t1\nb\tshipped\tminor%25α\t0\t1\n"
			+"b\tselected\tmajor\t0\t1\nb\tselected\tminor%25α\t"+(hidden ? 0 : 1)+"\t1\n";
		write(prefix+".coverage.tsv", "input_id\tpolicy\tmember\tmatched_scored_parents\tqualified_scored_parents\n"+member);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Qualified Support     ----------------*/
	/*--------------------------------------------------------------*/

	private void support(){
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ArrayList<Row> rows=new ArrayList<Row>();
		rows.add(codec.create("a", new Locus("seq%α", "TPP", 10, 30, '+'), false, null, null, "literal_fixture"));
		rows.add(codec.create("a", new Locus("seq%α", "TPP", 40, 60, '+'), true, "marginal", null, "literal_fixture"));
		final ByteBuilder inventory=new ByteBuilder("input_id\tparent_id\tfamily\tseqid\tstart\tstop\tstrand\tstatus\tscored_members\tmarginal_members\n");
		for(int i=0; i<rows.size(); i++){
			final Row row=rows.get(i);
			prok.BuildNcrnaReference.field(inventory, "a").tab();
			prok.BuildNcrnaReference.field(inventory, row.id).append("\tTPP\t");
			prok.BuildNcrnaReference.field(inventory, "seq%α").tab().append(row.locus.start).tab().append(row.locus.stop)
				.append(i==0 ? "\t+\tscored\tmajor," : "\t+\tdont_care\t.\t");
			// Use the inventory producer's UTF-8 field codec; append(String) narrows chars.
			prok.BuildNcrnaReference.field(inventory, "minor%α").append(i==0 ? "\t.\n" : "\n");
		}
		final String text=new String(inventory.toBytes(), StandardCharsets.UTF_8);
		write("support.tsv", text);
		final GradeNcrnaReference.Binding[] bindings=GradeNcrnaReference.loadManifest(path("reference.tsv"))
			.subList(0, 2).toArray(new GradeNcrnaReference.Binding[2]);
		final NcrnaMemberCoverage support=new NcrnaMemberCoverage(path("members.tsv"), path("support.tsv"), "TPP", bindings);
		final long[] bits=support.bind(0, rows);
		check(bits[0]==3 && bits[1]==0 && support.size()==2,
			"One scored parent can support both members; marginal masked evidence protects neither");
		support.bind(1, new ArrayList<Row>());
		write("support.tsv", text.replace("\t10\t30\t", "\t11\t30\t"));
		rejected(() -> new NcrnaMemberCoverage(path("members.tsv"), path("support.tsv"), "TPP", bindings).bind(0, rows),
			"Qualified support cannot join mismatched parent geometry");
		write("support.tsv", text+text.substring(text.indexOf('\n')+1));
		rejected(() -> new NcrnaMemberCoverage(path("members.tsv"), path("support.tsv"), "TPP", bindings),
			"Repeated parent evidence cannot inflate member coverage");
		write("support.tsv", text.replace("major,minor%25α", "unknown"));
		rejected(() -> new NcrnaMemberCoverage(path("members.tsv"), path("support.tsv"), "TPP", bindings),
			"Undeclared support must fail rather than becoming missing coverage");
		write("support.tsv", text);
		final String originalMap=read(path("members.tsv"));
		write("members.tsv", originalMap.replace("RF00002", "RF00001"));
		rejected(() -> new NcrnaMemberCoverage(path("members.tsv"), path("support.tsv"), "TPP", bindings),
			"Two reference members cannot share a reference model");
		write("members.tsv", originalMap.replace("RF00002", "."));
		rejected(() -> new NcrnaMemberCoverage(path("members.tsv"), path("support.tsv"), "TPP", bindings),
			"The missing-model marker cannot identify a qualified reference member");
		write("members.tsv", originalMap.replace("RF00002", ""));
		rejected(() -> new NcrnaMemberCoverage(path("members.tsv"), path("support.tsv"), "TPP", bindings),
			"An empty reference model is invalid");
		write("members.tsv", originalMap);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Driver Lifecycle      ----------------*/
	/*--------------------------------------------------------------*/

	private void driver(){
		write("identity.tsv", "source_commit\t0000000000000000000000000000000000000000\nbuild_seal\t00000000000000000000\n"
			+"resources_sha80\t00000000000000000000\ninputs_sha80\t00000000000000000000\nreferences_sha80\t00000000000000000000\n");
		write("familymap.tsv", "synthetic evaluator does not load caller resources\n");
		write("initial.tsv", NcrnaJointPolicy.HEADER+"\nroot\ttpp\tadaptive\t1\t1\t1\t0\t0\t-1\n"
			+"start\ttpp\tadaptive\t1\t1\t1\t0\t0\t0\n");
		write("grid.tsv", "member\tknob\tvalues\ntpp\ttopn\t1\ntpp\tfloor\t1\ntpp\ttopfrac\t0\ntpp\tqfrac\t0\ntpp\tpatience\t0\n");
		final Synthetic runner=new Synthetic();
		new OptimizeNcrnaShortlist(arguments("disabled", -1, -1, false), runner).process();
		final String legacy=read(path("disabled/result.tsv"));
		check(legacy.equals("backend\tstatus\tTP\tFP\tlogical_alignments\tlost_root\ttarget_TP\n"
			+"synthetic_test\tLOCAL_POINT_MET_TARGET\t1\t0\t10\t1\t1\n"), "Disabled budgets retain the exact original result schema and choice");
		check(!new File(path("disabled/coverage_baseline.tsv")).exists(), "Disabled guards add no baseline sidecar");
		// Inspect option keys at line starts, not matching text inside fixture paths.
		check(!read(path("disabled/contract.txt")).contains("\ncoverage_")
			&& !read(path("disabled/round0/sweep.config")).contains("\ncoverage="), "Disabled protection retains the original contract and child arguments");
		final String[] guarded=arguments("guarded", 0, -1, false);
		new OptimizeNcrnaShortlist(guarded, runner).process();
		check(read(path("guarded/result.tsv")).contains("LOCAL_STOP_CONSTRAINTS_UNMET")
			&& read(path("guarded/result.tsv")).endsWith("\tfalse\ttrue\n"),
			"An infeasible center is not a successful point even when total TP meets target");
		final int completed=runner.calls;
		new OptimizeNcrnaShortlist(arguments("guarded", 0, -1, true), runner).process();
		check(runner.calls==completed, "Coverage-enabled resume reuses saved tables without another evaluator call");
		rejected(() -> new OptimizeNcrnaShortlist(arguments("guarded", 1, -1, true), runner).process(),
			"Changing a coverage budget cannot silently resume old choices");
		check(runner.calls==completed, "Rejected resume must not execute another round");
		new OptimizeNcrnaShortlist(arguments("permitted", 1, -1, false), runner).process();
		check(read(path("permitted/result.tsv")).contains("LOCAL_POINT_MET_TARGET"), "An explicit one-input allowance can admit the same point");
		new OptimizeNcrnaShortlist(arguments("member_guarded", 1, 0, false), runner).process();
		check(read(path("member_guarded/result.tsv")).contains("LOCAL_STOP_CONSTRAINTS_UNMET"),
			"Member protection composes with the independently permissive union budget");
		new OptimizeNcrnaShortlist(arguments("member_only", -1, 0, false), runner).process();
		final int memberCompleted=runner.calls;
		new OptimizeNcrnaShortlist(arguments("member_only", -1, 0, true), runner).process();
		check(runner.calls==memberCompleted, "Member-only resume validates saved counts without invoking the evaluator");
		final String originalMap=read(path("members.tsv"));
		write("members.tsv", originalMap.replace("RF00001", "RF00003"));
		rejected(() -> new OptimizeNcrnaShortlist(arguments("member_only", -1, 0, true), runner).process(),
			"Changing literal member-map bytes must reject resume even with unchanged declared pins");
		write("members.tsv", originalMap);
		rejected(() -> new OptimizeNcrnaShortlist(with(arguments("member_only", -1, 0, true),
			"locisha80", "11111111111111111111"), runner).process(), "Changing the frozen support identity must reject resume");
		rejected(() -> new OptimizeNcrnaShortlist(arguments("member_only", -1, 1, true), runner).process(),
			"Changing the member budget cannot reuse old choices");
		check(runner.calls==memberCompleted, "Member-map, support-identity and budget rejection must all precede evaluator work");
		write("grid.tsv", "member\tknob\tvalues\ntpp\ttopn\t1,2\ntpp\tfloor\t1\ntpp\ttopfrac\t0\ntpp\tqfrac\t0\ntpp\tpatience\t0\n");
		final ArrayList<String> limited=new ArrayList<String>();
		for(final String arg:arguments("limited", 0, -1, false)){limited.add(arg);}
		limited.add("maxrounds=1");
		final Synthetic neighbor=new Synthetic();
		neighbor.addNeighbor=true;
		new OptimizeNcrnaShortlist(limited.toArray(new String[limited.size()]), neighbor).process();
		check(read(path("limited/result.tsv")).contains("ROUND_LIMIT\t1\t0\t8\t0\t1")
			&& read(path("limited/result.tsv")).endsWith("\ttrue\ttrue\n"),
			"Round-limit results retain the measured feasible neighbor and explicit constraint/target flags");
		write("grid.tsv", "member\tknob\tvalues\ntpp\ttopn\t1,2\ntpp\tfloor\t1\ntpp\ttopfrac\t0\ntpp\tqfrac\t0\ntpp\tpatience\t0,1\n");
		final Neighborhood landscape=new Neighborhood();
		new OptimizeNcrnaShortlist(arguments("neighborhood", 0, -1, false), landscape).process();
		check(landscape.calls==3 && read(path("neighborhood/result.tsv")).contains("LOCAL_POINT_MET_TARGET\t2\t0\t60\t0\t2"),
			"The literal landscape takes a free saving, then a paid recovery to the shipped target");
		check(read(path("neighborhood/round0/coverage.tsv")).contains("0_0_0_0_1_\tnew\tfalse\tconstraints_pass")
			&& read(path("neighborhood/round2/coverage.tsv")).contains("0_0_0_0_1_\tcached\tfalse\ttarget_TP;"),
			"A cached below-target neighbor needs today's rejection reason after the center attains target");
		new OptimizeNcrnaShortlist(arguments("neighborhood", 0, -1, true), landscape).process();
		check(landscape.calls==3, "Full-neighborhood reason replay must not reevaluate cached points or completed rounds");
	}

	private String[] arguments(final String name, final long maxLast, final long maxMember, final boolean resume){
		final ArrayList<String> args=new ArrayList<String>();
		for(final String option:new String[]{"caller=caller.config", "reference=reference.tsv", "familymap=familymap.tsv",
			"initial=initial.tsv", "grid=grid.tsv", "identity=identity.tsv"}){
			final int equal=option.indexOf('=');
			args.add(option.substring(0, equal+1)+path(option.substring(equal+1)));
		}
		args.add("target=tpp");
		args.add("grading=target");
		args.add("targetrecall=1");
		args.add("out="+path(name));
		args.add("maxlastloss="+maxLast);
		args.add("maxmemberlastloss="+maxMember);
		args.add("resume="+resume);
		if(maxMember>=0){
			args.add("members="+path("members.tsv"));
			args.add("loci="+path("support.tsv"));
			args.add("membersha80=00000000000000000000");
			args.add("locisha80=00000000000000000000");
		}
		return args.toArray(new String[args.size()]);
	}

	private final class Synthetic implements OptimizeNcrnaShortlist.RoundRunner{
		@Override
		public void run(final File prefix, final ArrayList<NcrnaJointPolicy> policies){
			calls++;
			check(policies.size()==(addNeighbor ? 3 : 2), "The fixture evaluates only its declared root, center and optional neighbor");
			counts(prefix.getPath(), false);
			for(final String suffix:new String[]{".inputs.tsv", ".summary.tsv", ".coverage.tsv"}){
				write(prefix+suffix, read(prefix+suffix).replace("selected", "center"));
			}
			if(addNeighbor){
				final String input=read(prefix+".inputs.tsv");
				write(prefix+".inputs.tsv", input.replace("a\tcenter\t0\t0\t2\t2\t0\t1\t0\t5\n",
					"a\tcenter\t0\t0\t2\t2\t0\t1\t0\t5\na\tneighbor1\t1\t0\t1\t2\t0\t0\t0\t4\n")
					.replace("b\tcenter\t1\t0\t1\t2\t0\t0\t1\t5\n",
					"b\tcenter\t1\t0\t1\t2\t0\t0\t1\t5\nb\tneighbor1\t0\t0\t2\t2\t0\t0\t0\t4\n"));
				write(prefix+".summary.tsv", read(prefix+".summary.tsv")+"neighbor1\t1\t0\t3\t4\t0\t0\t0\t8\t0\t0\n");
			}
		}
		int calls;
		boolean addNeighbor;
	}

	/** Only these four literal outcomes exist; current-center ranking is the SUT. */
	private final class Neighborhood implements OptimizeNcrnaShortlist.RoundRunner{
		@Override
		public void run(final File prefix, final ArrayList<NcrnaJointPolicy> policies){
			calls++;
			if(calls==3){check(policies.size()==2, "The third round must not reevaluate its already measured neighbor");}
			final ByteBuilder input=new ByteBuilder("input_id\tpolicy\t"+METRICS+"\n");
			final ByteBuilder summary=new ByteBuilder("policy\t"+METRICS+"\tlast_loss_inputs\tfirst_gain_inputs\n");
			for(int i=0; i<2; i++){
				for(final NcrnaJointPolicy policy:policies){
					final NcrnaJointPolicy.Member member=policy.member(0);
					final int tp, cost;
					if(member.patience<0){
						tp=2;
						cost=100;
					}else{
						check(member.topN>=1 && member.topN<=2 && member.patience>=0 && member.patience<=1,
							"Unknown literal neighborhood coordinate");
						tp=NEIGHBOR_TP[member.topN-1][member.patience];
						cost=NEIGHBOR_COST[member.topN-1][member.patience];
					}
					input.append(i==0 ? "a\t" : "b\t").append(policy.name).tab().append(i==0 ? tp : 0)
						.append("\t0\t").append(i==0 ? 2-tp : 2).append("\t2\t0\t").append(i==0 ? 2-tp : 0)
						.append("\t0\t").append(cost/2).nl();
					if(i==0){
						summary.append(policy.name).tab().append(tp).append("\t0\t").append(4-tp)
							.append("\t4\t0\t").append(2-tp).append("\t0\t").append(cost).append("\t0\t0\n");
					}
				}
			}
			write(prefix+".inputs.tsv", new String(input.toBytes(), StandardCharsets.UTF_8));
			write(prefix+".summary.tsv", new String(summary.toBytes(), StandardCharsets.UTF_8));
		}
		int calls;
	}

	private static Measurement[] measurements(){return new Measurement[]{new Measurement(1, 0, 20, 0), new Measurement(1, 0, 10, 1)};}
	private static String[] with(final String[] args, final String key, final String value){
		final String[] changed=args.clone();
		for(int i=0; i<changed.length; i++){
			if(changed[i].startsWith(key+"=")){
				changed[i]=key+"="+value;
				return changed;
			}
		}
		throw new AssertionError("Missing fixture option: "+key);
	}
	private String path(final String name){return new File(root, name).getPath();}
	private void write(final String name, final String text){
		final String file=new File(name).isAbsolute() ? name : path(name);
		try(Output output=new Output(FileFormat.testOutput(file, FileFormat.TXT, null, true, true, false, false))){
			output.writer.print(new ByteBuilder(text.getBytes(StandardCharsets.UTF_8)));
		}
	}
	private static String read(final String path){
		final ByteBuilder text=new ByteBuilder();
		try(Input input=new Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){text.append(line).nl();}
		}
		return new String(text.toBytes(), StandardCharsets.UTF_8);
	}
	private static void rejected(final Runnable action, final String reason){
		try{action.run();}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError(reason);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
	private final String[] caller;
	private static final String[] POLICIES={"shipped", "selected"};
	private static final String METRICS="TP\tFP\tFN\tscored_references\tmasked_references\tlost_baseline_TP\tgained_TP\tkernel_logical_alignments";
	private static final int[][] NEIGHBOR_TP={{1, 1}, {1, 2}}, NEIGHBOR_COST={{40, 50}, {30, 60}};
}
