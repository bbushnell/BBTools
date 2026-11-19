package test.prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

import ddl.FamilySketchBuilder;
import prot.ProteinSequence;
import prot.cluster.ClusteringInput;
import prot.cluster.ClusteringShardWriter;
import prot.cluster.ClusterSketchCodec;
import prot.cluster.DdlClustStageB;
import prot.cluster.Family;
import prot.cluster.LiveProbeMetrics;
import prot.cluster.LiveProbeState;
import prot.cluster.LocalProbeBaseline;
import prot.cluster.MemberRankInputs;
import prot.cluster.StructuralProbeBatch;

/** Fixtures for Stage-B consensus collection and original-member union rebuilds. */
public final class DdlClustStageBTest{

	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Require NEW fixture directory");}
		run(Paths.get(args[0]));
	}

	public static void run(Path root) throws Exception{
		Files.createDirectory(root);
		final String a=repeat("ACDEFGHIKLMNPQRSTVWY", 4);
		final String b=repeat("RRRRACDEFGHIKLMNPQ", 5);
		final ProteinSequence[] originals=new ProteinSequence[]{
			new ProteinSequence("p0", a), new ProteinSequence("p1", b),
			new ProteinSequence("p2", a), new ProteinSequence("p3", repeat("M", 80))};
		check(decode(new ProteinSequence("roundtrip", "ILWD").enc).equals("ILWD"),
			"Test consensus decoder preserves BBTools amino-acid codes");
		final FamilySketchBuilder builder=new FamilySketchBuilder(64, 2, 7, false, 0, 0, null);
		final Path stage0=root.resolve("stage0.pddl"), stage1=root.resolve("stage1.pddl");
		ClusterSketchCodec.write(stage0.toString(), false,
			new Family[]{family(10, "s0_f10", originals, builder, 0),
				family(11, "s0_f11", originals, builder, 1)}, 2, 64, 7, false, 0, 0);
		ClusterSketchCodec.write(stage1.toString(), false,
			new Family[]{family(12, "s1_f12", originals, builder, 2)}, 2, 64, 7, false, 0, 0);

		final Path collectFasta=root.resolve("stageb_input.faa"), roster=root.resolve("roster.tsv");
		new DdlClustStageB(new String[]{
			"mode=collect", "codecs="+stage0+","+stage1, "out="+collectFasta, "roster="+roster, "ow=f"}).process();
		final List<String> rosterLines=Files.readAllLines(roster, StandardCharsets.US_ASCII);
		check(rosterLines.size()==4 && rosterLines.get(1).equals("0\ts0_f10\t0\t10\t1") &&
			rosterLines.get(3).equals("2\ts1_f12\t1\t12\t1"), "Collect roster preserves shard/family ordinals");
		final String collected=new String(Files.readAllBytes(collectFasta), StandardCharsets.US_ASCII);
		check(collected.contains(">s0_f10 ") && collected.contains(">s1_f12 "), "Collect FASTA uses canonical Stage-A IDs");

		final Path stageB=root.resolve("stageb.pddl");
		final Family merged=new Family(20, "stageb_20", 0, originals[0].enc);
		merged.addMember(2); merged.refreshStats(originals, builder);
		ClusterSketchCodec.write(stageB.toString(), false,
			new Family[]{merged, family(21, "stageb_21", originals, builder, 1)}, 2, 64, 7, false, 0, 0);
		final Path originalFasta=root.resolve("originals.faa");
		fasta(originalFasta, originals);
		final Path map=root.resolve("parent.tsv"), stage0Members=root.resolve("stage0.members.tsv"),
			stage1Members=root.resolve("stage1.members.tsv"), stageBMembers=root.resolve("stageb.members.tsv"),
			out=root.resolve("union.pddl"), report=root.resolve("union.tsv");
		write(map, "parent_input_index\tshard\tshard_input_index\tinput_id\n0\t0\t0\tp0\n1\t0\t1\tp1\n2\t1\t0\tp2\n3\t1\t1\tp3\n");
		write(stage0Members, "input_index\tinput_id\tfamily_id\tstatus\n0\tp0\t10\tassigned\n1\tp1\t11\tassigned\n");
		write(stage1Members, "input_index\tinput_id\tfamily_id\tstatus\n0\tp2\t12\tassigned\n1\tp3\t-1\tunassigned\n");
		write(stageBMembers, "input_index\tinput_id\tfamily_id\tstatus\n0\ts0_f10\t20\tassigned\n1\ts0_f11\t21\tassigned\n2\ts1_f12\t20\tassigned\n");

		final Path out4=root.resolve("union.t4.pddl"), report4=root.resolve("union.t4.tsv");
		new DdlClustStageB(new String[]{"mode=union", "in="+originalFasta, "out="+out, "roster="+roster,
			"partitionmap="+map, "stageacodecs="+stage0+","+stage1,
			"stageamembers="+stage0Members+","+stage1Members, "stagebcodec="+stageB,
			"stagebmembers="+stageBMembers, "report="+report, "minspan=5", "t=1", "ow=f"}).process();
		new DdlClustStageB(new String[]{"mode=union", "in="+originalFasta, "out="+out4, "roster="+roster,
			"partitionmap="+map, "stageacodecs="+stage0+","+stage1,
			"stageamembers="+stage0Members+","+stage1Members, "stagebcodec="+stageB,
			"stagebmembers="+stageBMembers, "report="+report4, "minspan=5", "t=4", "ow=f"}).process();
		check(Arrays.equals(Files.readAllBytes(out), Files.readAllBytes(out4)),
			"Union codec output is byte-identical with one versus four workers");
		check(Arrays.equals(Files.readAllBytes(report), Files.readAllBytes(report4)),
			"Union report output is byte-identical with one versus four workers");
		final Path explicitUnion=root.resolve("union.explicit_one.pddl"), explicitReport=root.resolve("union.explicit_one.tsv");
		final String[] serialUnion=unionArgs(originalFasta, explicitUnion, roster, map, stage0, stage1,
			stage0Members, stage1Members, stageB, stageBMembers, explicitReport);
		final String[] serialUnionOption=Arrays.copyOf(serialUnion,serialUnion.length+2);
		serialUnionOption[serialUnion.length]="t=4"; serialUnionOption[serialUnion.length+1]="frozenindexthreads=1";
		new DdlClustStageB(serialUnionOption).process();
		check(Arrays.equals(Files.readAllBytes(out4),Files.readAllBytes(explicitUnion)) &&
			Arrays.equals(Files.readAllBytes(report4),Files.readAllBytes(explicitReport)),
			"Explicit serial frozen index leaves union mode bytes unchanged");
		final ClusterSketchCodec.Loaded union=ClusterSketchCodec.read(out.toString());
		check(union.families.length==2 && union.families[0].familyId==20 && union.families[0].members==2 &&
			union.families[1].familyId==21 && union.families[1].members==1, "Union codec retains original-member counts");
		final List<String> reportLines=Files.readAllLines(report, StandardCharsets.US_ASCII);
		check(reportLines.size()==3 && reportLines.get(1).contains("\tunion_"), "Union report records rebuild status");

		final byte[] stageBBytes=Files.readAllBytes(stageB);
		final Path stageBAlias=root.resolve("stageb.alias.pddl");
		Files.createLink(stageBAlias, stageB);
		expectFail("Hardlinked output aliases a Stage-B codec",
			"Stage-B input/output paths alias",
			unionArgs(originalFasta, stageBAlias, roster, map, stage0, stage1, stage0Members,
				stage1Members, stageB, stageBMembers, report));
		check(Arrays.equals(stageBBytes, Files.readAllBytes(stageB)), "Alias preflight preserves the linked input");

		final Path duplicate=root.resolve("stageb.members.duplicate.tsv");
		write(duplicate, "input_index\tinput_id\tfamily_id\tstatus\n0\ts0_f10\t20\tassigned\n0\ts0_f10\t20\tassigned\n2\ts1_f12\t20\tassigned\n");
		expectFail("Stage-B duplicate ordinal",
			"Duplicate or out-of-range Stage-A ordinal",
			unionArgs(originalFasta, out, roster, map, stage0, stage1, stage0Members, stage1Members, stageB, duplicate, report));
		final Path missing=root.resolve("stageb.members.missing.tsv");
		write(missing, "input_index\tinput_id\tfamily_id\tstatus\n0\ts0_f10\t20\tassigned\n2\ts1_f12\t20\tassigned\n");
		expectFail("Stage-B missing ordinal",
			"Stage-B omitted Stage-A ordinal",
			unionArgs(originalFasta, out, roster, map, stage0, stage1, stage0Members, stage1Members, stageB, missing, report));
		final Path countMismatch=root.resolve("stageb.members.count_mismatch.tsv");
		write(countMismatch, "input_index\tinput_id\tfamily_id\tstatus\n0\ts0_f10\t20\tassigned\n1\ts0_f11\t21\tassigned\n2\ts1_f12\t21\tassigned\n");
		expectFail("Stage-B codec/member count mismatch",
			"Stage-B codec member count disagrees",
			unionArgs(originalFasta, out, roster, map, stage0, stage1, stage0Members, stage1Members,
				stageB, countMismatch, report));
		final Path badRoster=root.resolve("roster.bad.tsv");
		write(badRoster, "stagea_ordinal\tstagea_id\tshard\tstagea_family_id\tsource_members\n"+
			"0\ts0_f10\t0\t10\t1\n1\ts0_f11\t0\t11\t2\n2\ts1_f12\t1\t12\t1\n");
		expectFail("Roster no longer matches Stage-A codec",
			"Roster no longer matches the ordered Stage-A codec list",
			unionArgs(originalFasta, out, badRoster, map, stage0, stage1, stage0Members,
				stage1Members, stageB, stageBMembers, report));
		final Path badMap=root.resolve("parent.bad.tsv");
		write(badMap, "parent_input_index\tshard\tshard_input_index\tinput_id\n0\t0\t0\tp0\n1\t0\t1\tp1\n2\t1\t0\tp2\n3\t1\t4\tp3\n");
		expectFail("Partition local ordinal is bounded before allocation",
			"Out-of-range shard or local ordinal",
			unionArgs(originalFasta, out, roster, badMap, stage0, stage1, stage0Members,
				stage1Members, stageB, stageBMembers, report));
		expectFail("Reads-limited parser options are rejected",
			"reads= is unsupported",
			new String[]{"mode=union", "in="+originalFasta, "out="+out, "roster="+roster,
				"partitionmap="+map, "stageacodecs="+stage0+","+stage1,
				"stageamembers="+stage0Members+","+stage1Members, "stagebcodec="+stageB,
				"stagebmembers="+stageBMembers, "report="+report, "reads=1", "ow=t"});
		expectFail("Output pseudo-stream aliases are rejected",
			"Stage-B requires on-disk files",
			unionArgs(originalFasta, Paths.get("stdout.pddl"), roster, map, stage0, stage1,
				stage0Members, stage1Members, stageB, stageBMembers, report));

		final CaseResult accepted=singleUnionCase(root.resolve("accepted"), repeat("I", 80),
			new String[]{repeat("L", 80), repeat("L", 80)}, 0.35);
		check(accepted.report.contains("\tunion_accepted\tall_members_pass\t"),
			"Accepted Stage-B union is reported");
		check(accepted.consensus.equals(repeat("L", 80)), "Accepted union publishes the rebuilt consensus");
		final CaseResult strictRejected=singleUnionCase(root.resolve("strict_reject"), repeat("I", 80),
			new String[]{repeat("I", 80), repeat("L", 80)}, 0.95);
		check(strictRejected.report.contains("\tunion_rejected\tstrict_new_validation_failed\t"),
			"Strict validation failure is reported");
		check(strictRejected.consensus.equals(repeat("I", 80)), "Strict rejection preserves the Stage-B seed");
		final CaseResult noPositive=singleUnionCase(root.resolve("no_positive"), repeat("W", 80),
			new String[]{repeat("D", 80)}, 0.35);
		check(noPositive.report.contains("\tunion_rejected\tsnapshot_no_positive_placement\t"),
			"No-positive inclusive placement is reported");
		check(noPositive.consensus.equals(repeat("W", 80)), "No-positive rejection preserves the Stage-B seed");
		mutualCase(root.resolve("mutual"));
		mutualParallelPairs(root.resolve("mutual_parallel"));
		frozenIndexInvalid(root);
		fissionProbe(root.resolve("fission_probe"));
		liveProbe(root.resolve("live_probe"));
		liveProbeLoss(root.resolve("live_loss"));
		liveCycle(root.resolve("live_cycle"));
		localPlanRanking(root.resolve("local_plan_ranking"));
		mutualRanking(root.resolve("mutual_ranking"));
		System.out.println("DDLCLUST_STAGEB_TEST_PASS checks="+checks);
	}

	/** Identical models isolate current-member ranking from consensus and stale source descriptors. */
	private static void localPlanRanking(final Path root)throws Exception{
		Files.createDirectory(root);
		final ProteinSequence[] input={new ProteinSequence("a",repeat("A",99)),new ProteinSequence("b",repeat("I",100)),
			new ProteinSequence("c",repeat("D",40)),new ProteinSequence("d",repeat("A",99)),
			new ProteinSequence("e",repeat("I",100)),new ProteinSequence("f",repeat("R",160))};
		final FamilySketchBuilder builder=new FamilySketchBuilder(64,2,7,false,0,0,null);
		final MemberRankInputs cache=new MemberRankInputs(input);
		// Source membership is deliberately different from the completed current assignment.
		// Reading serialized nonempty measurements without rehydration would choose10--50.
		final Family[] models={localPlanModel(90,builder,2),localPlanModel(10,builder,0,1),
			localPlanModel(50,builder,3,4),localPlanModel(70,builder,5)};
		final Path plain=root.resolve("plain.pddl"),ranked=root.resolve("ranked.pddl");
		ClusterSketchCodec.write(plain.toString(),false,models,2,64,7,false,0,0);
		for(Family model:models){model.enableMemberRanking(cache,3,null);}
		ClusterSketchCodec.write(ranked.toString(),false,models,2,64,7,false,0,0);
		fasta(root.resolve("input.faa"),input);
		write(root.resolve("selected.tsv"),"family_id\tsource_members\n90\t1\n10\t2\n50\t2\n70\t1\n");
		write(root.resolve("baseline.tsv"),"input_index\tinput_id\tfamily_id\tstatus\n"+
			"0\ta\t90\tassigned\n1\tb\t90\tassigned\n2\tc\t10\tassigned\n"+
			"3\td\t50\tassigned\n4\te\t50\tassigned\n5\tf\t-1\tunassigned\n");
		final String[] inputs={"plain.pddl","ranked.pddl","input.faa","selected.tsv","baseline.tsv"};
		final byte[][] before=new byte[inputs.length][];
		for(int i=0;i<inputs.length;i++){before[i]=Files.readAllBytes(root.resolve(inputs[i]));}
		new DdlClustStageB(localPlanArgs(root,"off_plain",plain,1,"memberranking=f")).process();
		new DdlClustStageB(localPlanArgs(root,"off_ranked",ranked,1,"memberranking=f")).process();
		new DdlClustStageB(localPlanArgs(root,"on",ranked,1,"memberranking=t")).process();
		new DdlClustStageB(localPlanArgs(root,"implicit",ranked,1)).process();
		new DdlClustStageB(localPlanArgs(root,"parallel",ranked,4)).process();
		for(String suffix:new String[]{".requests.tsv",".metadata.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("off_plain"+suffix)),Files.readAllBytes(root.resolve("off_ranked"+suffix))),
				"Localplan explicit OFF preserves legacy bytes on the same model payload with or without descriptors");
			for(String other:new String[]{"implicit","parallel"}){
				check(Arrays.equals(Files.readAllBytes(root.resolve("on"+suffix)),Files.readAllBytes(root.resolve(other+suffix))),
					"Localplan omitted ranking equals explicit ON and requested-worker parity: "+other+suffix);
			}
		}
		final List<String> legacy=Files.readAllLines(root.resolve("off_plain.requests.tsv"),StandardCharsets.US_ASCII);
		final List<String> current=Files.readAllLines(root.resolve("on.requests.tsv"),StandardCharsets.US_ASCII);
		check(legacy.size()==2 && legacy.get(0).equals(StructuralProbeBatch.REQUEST_HEADER) &&
			legacy.get(1).equals("0\tfusion\t10\t50\t17"),"Legacy identical-payload tie chooses literal pair10--50");
		check(current.size()==2 && current.get(1).equals("0\tfusion\t50\t90\t17"),
			"Actual assigned A99/I100 members choose literal pair50--90, not consensus or stale source measurements");
		final List<String> metadata=Files.readAllLines(root.resolve("on.metadata.tsv"),StandardCharsets.US_ASCII);
		check(metadata.contains("member_ranking\ttrue") && metadata.contains("member_ranking_layer\t4") &&
			metadata.contains("selected_models\t4") && metadata.contains("positive_models\t3"),
			"Ranked metadata records hydration provenance and preserves the empty selected template");
		check(!Files.readAllLines(root.resolve("off_ranked.metadata.tsv"),StandardCharsets.US_ASCII).contains("member_ranking\ttrue"),
			"Legacy metadata retains its old schema");
		expectFail("Default localplan requires descriptors","Ranked hydration requires serialized member descriptors",
			localPlanArgs(root,"missing",plain,1));
		check(!Files.exists(root.resolve("missing.requests.tsv")) && !Files.exists(root.resolve("missing.metadata.tsv")),
			"Missing descriptors fail before publishing requests or metadata");
		// One empty current template still needs its own serialized last-nonempty provenance.
		final Family[] mixed=models.clone();mixed[3]=localPlanModel(70,builder,5);
		final Path mixedPath=root.resolve("mixed.pddl");
		ClusterSketchCodec.write(mixedPath.toString(),false,mixed,2,64,7,false,0,0);
		expectFail("Empty localplan template cannot use consensus as descriptor fallback",
			"Ranked hydration requires serialized member descriptors",localPlanArgs(root,"missing_empty",mixedPath,1));
		for(int i=0;i<inputs.length;i++){
			check(Arrays.equals(before[i],Files.readAllBytes(root.resolve(inputs[i]))),"Localplan leaves its input unchanged: "+inputs[i]);
		}
		localProbeRankingChoice(root,ranked);
	}

	/** Driver selection uses current assigned members, despite identical source model payloads. */
	private static void localProbeRankingChoice(final Path root,final Path source)throws Exception{
		final Path requests=root.resolve("local_choices.tsv");
		write(requests,StructuralProbeBatch.REQUEST_HEADER+"\n0\tfusion\t10\t50\t17\n1\tfusion\t50\t90\t17\n");
		for(String mode:new String[]{"off","on","default","parallel"}){
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("mode=localprobe","in="+root.resolve("input.faa"),
				"codecs="+source,"selectedseeds="+root.resolve("selected.tsv"),"baselinemembers="+root.resolve("baseline.tsv"),
				"probes="+requests,"out="+root.resolve("choices_"+mode+".tsv"),"selectedcount=4","minlen=40","minspan=6",
				"maxalignlen=200","trimleft=f","shortlist=1","mergeminid=0","mergeminratio=0","passhits=0","maxhydrodelta=9",
				"t="+(mode.equals("parallel") ? 4 : 1),"ow=f"));
			if(mode.equals("off") || mode.equals("on")){args.add("memberranking="+mode.equals("on"));}
			new DdlClustStageB(args.toArray(new String[args.size()])).process();
			final List<String> rows=Files.readAllLines(root.resolve("choices_"+mode+".tsv"),StandardCharsets.US_ASCII);
			check(rows.size()==3,"Both local candidate requests retain their independent response rows");
			for(int row=1;row<=2;row++){
				final String reason=rows.get(row).split("\t",-1)[7];
				final boolean reciprocal=mode.equals("off") ? row==1 : row==2;
				check(reason.equals(reciprocal ? "input_quality_unavailable" : "fusion_not_reciprocal"),
					"Current-member ranking changes the reciprocal endpoint pair before independent construction eligibility");
			}
		}
		cycleTableParity(root.resolve("choices_on.tsv"),root.resolve("choices_default.tsv"),28,3);
		cycleTableParity(root.resolve("choices_on.tsv"),root.resolve("choices_parallel.tsv"),28,3);
		System.out.println("LOCALPROBE_RANK_PAIR legacy=10,50 ranked=50,90 proposal=input_quality_unavailable");
	}

	/** Supplies equal model bytes while preserving distinct original ordinal membership. */
	private static Family localPlanModel(final int id,final FamilySketchBuilder builder,final int... ordinals){
		assert(ordinals.length>0) : "Fixture source codec requires positive source membership before current reassignment";
		final byte[] common=new ProteinSequence("common",repeat("A",80)).enc;
		final Family model=new Family(id,"plan"+id,ordinals[0],common);
		for(int i=1;i<ordinals.length;i++){model.addMember(ordinals[i]);}
		final char[] sketch=new char[64];builder.sketchEncoded(common,sketch);model.setSketch(sketch);
		final float[] monomers=new float[20];Arrays.fill(monomers,.05f);model.setMonomers(monomers);
		return model;
	}

	/** One fusion drawn from a one-pair pool gives a literal oracle independent of random ordering. */
	private static String[] localPlanArgs(final Path root,final String prefix,final Path source,
			final int threads,final String... extra){
		final String[] base={"mode=localplan","in="+root.resolve("input.faa"),"codecs="+source,
			"selectedseeds="+root.resolve("selected.tsv"),"baselinemembers="+root.resolve("baseline.tsv"),
			"out="+root.resolve(prefix+".requests.tsv"),"report="+root.resolve(prefix+".metadata.tsv"),
			"selectedcount=4","minlen=40","minspan=6","maxalignlen=200","trimleft=f","shortlist=1",
			"mergeminid=0","mergeminratio=0","fissionprobes=0","fusionprobes=1","probeseed=17",
			"t="+threads,"frozenindexthreads="+threads,"ow=f"};
		final String[] args=Arrays.copyOf(base,base.length+extra.length);
		System.arraycopy(extra,0,args,base.length,extra.length);
		return args;
	}

	/** Multiple shards and repeated strict-veto unions retain original-member statistics. */
	private static void mutualRanking(final Path root)throws Exception{
		Files.createDirectory(root);
		final ProteinSequence[] input={new ProteinSequence("a",repeat("A",99)),new ProteinSequence("b",repeat("I",100)),
			new ProteinSequence("c",repeat("D",40)),new ProteinSequence("d",repeat("A",99)),
			new ProteinSequence("e",repeat("I",100))};
		final FamilySketchBuilder builder=new FamilySketchBuilder(64,2,7,false,0,0,null);
		final MemberRankInputs cache=new MemberRankInputs(input);
		final Family[] models={localPlanModel(90,builder,0,1),localPlanModel(10,builder,2),localPlanModel(50,builder,3,4)};
		for(String tag:new String[]{"plain","ranked"}){
			if(tag.equals("ranked")){
				for(int i=0;i<models.length;i++){models[i].enableMemberRanking(cache,5+i,null);}
			}
			ClusterSketchCodec.write(root.resolve(tag+"0.pddl").toString(),false,new Family[]{models[0],models[1]},2,64,7,false,0,0);
			ClusterSketchCodec.write(root.resolve(tag+"1.pddl").toString(),false,new Family[]{models[2]},2,64,7,false,0,0);
		}
		fasta(root.resolve("input.faa"),input);
		write(root.resolve("map.tsv"),ClusteringShardWriter.MAPPING_HEADER+"\n"+
			"0\t0\t2\ta\n1\t0\t0\tb\n2\t0\t1\tc\n3\t1\t1\td\n4\t1\t0\te\n");
		write(root.resolve("members0.tsv"),"input_index\tinput_id\tfamily_id\tstatus\n"+
			"0\tb\t90\tassigned\n1\tc\t10\tassigned\n2\ta\t90\tassigned\n");
		write(root.resolve("members1.tsv"),"input_index\tinput_id\tfamily_id\tstatus\n"+
			"0\te\t50\tassigned\n1\td\t50\tassigned\n");
		new DdlClustStageB(new String[]{"mode=collect","codecs="+root.resolve("plain0.pddl")+","+root.resolve("plain1.pddl"),
			"out="+root.resolve("collected.faa"),"roster="+root.resolve("roster.tsv"),"t=1","ow=f"}).process();
		final String[] names={"input.faa","map.tsv","members0.tsv","members1.tsv","plain0.pddl","plain1.pddl",
			"ranked0.pddl","ranked1.pddl","roster.tsv"};
		final byte[][] before=new byte[names.length][];
		for(int i=0;i<names.length;i++){before[i]=Files.readAllBytes(root.resolve(names[i]));}
		new DdlClustStageB(mutualRankingArgs(root,"off_plain","plain",1,"memberranking=f")).process();
		new DdlClustStageB(mutualRankingArgs(root,"off_ranked","ranked",1,"memberranking=f")).process();
		new DdlClustStageB(mutualRankingArgs(root,"on","ranked",1,"memberranking=t")).process();
		new DdlClustStageB(mutualRankingArgs(root,"implicit","ranked",1)).process();
		new DdlClustStageB(mutualRankingArgs(root,"parallel","ranked",4)).process();
		new DdlClustStageB(mutualRankingArgs(root,"legacy_source","plain",1)).process();
		for(String[] pair:new String[][]{{"off_plain","off_ranked"},{"on","implicit"},{"on","parallel"}}){
			for(String suffix:new String[]{".pddl",".pairs.tsv",".members.tsv"}){
				check(Arrays.equals(Files.readAllBytes(root.resolve(pair[0]+suffix)),Files.readAllBytes(root.resolve(pair[1]+suffix))),
					"Mutual default/explicit/worker/legacy payload parity: "+pair[0]+"/"+pair[1]+suffix);
			}
			cycleTableParity(root.resolve(pair[0]+".passes.tsv"),root.resolve(pair[1]+".passes.tsv"),10,4);
		}
		final List<String> off=Files.readAllLines(root.resolve("off_plain.pairs.tsv"),StandardCharsets.US_ASCII);
		final List<String> on=Files.readAllLines(root.resolve("on.pairs.tsv"),StandardCharsets.US_ASCII);
		check(off.size()==3 && on.size()==3 && off.get(1).startsWith("1\t0\t1\t2\t1\t") &&
			on.get(1).startsWith("1\t0\t2\t2\t2\t"),"Member statistics change the literal first reciprocal pair on identical models");
		for(int i=1;i<on.size();i++){
			final String[] row=on.get(i).split("\t",-1);
			check(Integer.parseInt(row[8])==0 && Integer.parseInt(row[9])>0,
				"Both mutual consensus proposals are rejected while all construction members are retained");
		}
		final ClusterSketchCodec.Entry result=ClusterSketchCodec.read(root.resolve("on.pddl").toString()).families[0];
		check(result.members==5 && result.memberRanking.members==5 && result.memberRanking.layer==10 &&
			result.memberRanking.medianLength==99 && Math.abs(result.memberRanking.meanKd-1.82)<1e-6,
			"Two unions retain all five original measurements with literal median99, meanKD1.82 and imported-layer advance8--10");
		final ClusterSketchCodec.Entry legacySource=ClusterSketchCodec.read(root.resolve("legacy_source.pddl").toString()).families[0];
		check(legacySource.memberRanking.layer==2 && legacySource.memberRanking.members==5 &&
			legacySource.memberRanking.medianLength==99 && Math.abs(legacySource.memberRanking.meanKd-1.82)<1e-6,
			"Legacy codecs get actual membership-derived measurements, starting at layer0 without invented descriptors");
		final List<String> owners=Files.readAllLines(root.resolve("on.members.tsv"),StandardCharsets.US_ASCII);
		check(owners.size()==6,"Every original mapped ordinal survives the two unions exactly once");
		for(int i=0;i<input.length;i++){check(owners.get(i+1).equals(i+"\t"+input[i].id+"\t0"),"Shard-local remapping preserves original identity");}
		new DdlClustStageB(mutualRankingArgs(root,"blocker","ranked",1,"mergemaxsize=3")).process();
		final ClusterSketchCodec.Loaded blocked=ClusterSketchCodec.read(root.resolve("blocker.pddl").toString());
		check(blocked.families.length==2 && blocked.families[0].members==4 && blocked.families[0].memberRanking.layer==9 &&
			blocked.families[1].members==1 && blocked.families[1].memberRanking.layer==8,
			"A ranked union above the size limit still blocks merging, while its unchanged neighbor retains its layer");
		final Family overflow=localPlanModel(90,builder,0,1);
		overflow.enableMemberRanking(cache,Integer.MAX_VALUE,null);
		ClusterSketchCodec.write(root.resolve("overflow0.pddl").toString(),false,new Family[]{overflow,models[1]},2,64,7,false,0,0);
		expectFail("Imported measurement layer cannot wrap","Mutual hydration cannot advance beyond",
			mutualRankingArgs(root,"overflow","ranked",1,"stageacodecs="+root.resolve("overflow0.pddl")+","+root.resolve("ranked1.pddl")));
		check(!Files.exists(root.resolve("overflow.pddl")) && !Files.exists(root.resolve("overflow.pairs.tsv")),
			"Invalid imported provenance fails before opening publication outputs");
		final Family lastLayer=localPlanModel(90,builder,0,1);
		lastLayer.enableMemberRanking(cache,Integer.MAX_VALUE-1,null);
		ClusterSketchCodec.write(root.resolve("last0.pddl").toString(),false,new Family[]{lastLayer,models[1]},2,64,7,false,0,0);
		expectFail("A new union cannot wrap after valid last-layer hydration","Mutual union cannot advance beyond",
			mutualRankingArgs(root,"last","ranked",1,"stageacodecs="+root.resolve("last0.pddl")+","+root.resolve("ranked1.pddl")));
		check(!Files.exists(root.resolve("last.pddl")),"Provenance overflow cannot publish a union codec");
		for(int i=0;i<names.length;i++){check(Arrays.equals(before[i],Files.readAllBytes(root.resolve(names[i]))),"Mutual leaves source input unchanged: "+names[i]);}
	}

	/** Fixed gates force consensus rejection, independently of the reciprocal merge gate. */
	private static String[] mutualRankingArgs(final Path root,final String prefix,final String source,
			final int threads,final String... extra){
		final String[] base={"mode=mutual","in="+root.resolve("input.faa"),"out="+root.resolve(prefix+".pddl"),
			"roster="+root.resolve("roster.tsv"),"partitionmap="+root.resolve("map.tsv"),
			"stageacodecs="+root.resolve(source+"0.pddl")+","+root.resolve(source+"1.pddl"),
			"stageamembers="+root.resolve("members0.tsv")+","+root.resolve("members1.tsv"),
			"report="+root.resolve(prefix+".pairs.tsv"),"passesreport="+root.resolve(prefix+".passes.tsv"),
			"members="+root.resolve(prefix+".members.tsv"),"minspan=6","maxalignlen=200","trimleft=f","consensustrimleft=f",
			"minid=1","minblosumratio=100","mergeminid=0","mergeminratio=0","shortlist=1","mergemaxsize=20",
			"t="+threads,"frozenindexthreads="+threads,"ow=f"};
		final String[] args=Arrays.copyOf(base,base.length+extra.length);
		System.arraycopy(extra,0,args,base.length,extra.length);
		return args;
	}

	/** Profitable fission must update the actual input to fusion within the same invocation. */
	private static void liveCycle(final Path root) throws Exception{
		Files.createDirectory(root);
		final ProteinSequence[] input=new ProteinSequence[10];
		final StringBuilder membership=new StringBuilder("input_index\tinput_id\tfamily_id\tstatus\n");
		for(int i=0;i<input.length;i++){
			final String sequence=repeat(i%4<2 && i!=9 ? "R" : "K",100);
			input[i]=new ProteinSequence("cycle"+i,sequence);
			membership.append(i).append("\tcycle").append(i).append(i<4 ? "\t100\tassigned\n" :
				i<8 ? "\t300\tassigned\n" : "\t-1\tunassigned\n");
		}
		final FamilySketchBuilder builder=new FamilySketchBuilder(64,2,7,false,0,0,null);
		final Family first=new Family(100,"first",0,input[0].enc),second=new Family(300,"second",4,input[4].enc);
		for(int i=1;i<4;i++){first.addMember(i);second.addMember(i+4);}
		first.refreshStats(input,builder);second.refreshStats(input,builder);
		ClusterSketchCodec.write(root.resolve("source.pddl").toString(),false,new Family[]{first,second},2,64,7,false,0,0);
		fasta(root.resolve("input.faa"),input);write(root.resolve("baseline.tsv"),membership.toString());
		write(root.resolve("selected.tsv"),"family_id\tsource_members\n100\t4\n300\t4\n");
		for(int threads:new int[]{1,4}){new DdlClustStageB(cycleArgs(root,"cycle"+threads,threads)).process();}
		for(String suffix:new String[]{".pddl",".members.tsv",".report.tsv",".fusion.tsv",".sizes.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("cycle1"+suffix)),Files.readAllBytes(root.resolve("cycle4"+suffix))),
				"Cycle serial/parallel deterministic export: "+suffix);
		}
		cycleTableParity(root.resolve("cycle1.passes.tsv"),root.resolve("cycle4.passes.tsv"),8,3);
		cycleTableParity(root.resolve("cycle1.metrics.tsv"),root.resolve("cycle4.metrics.tsv"),17,4);
		final List<String> phases=Files.readAllLines(root.resolve("cycle1.passes.tsv"),StandardCharsets.US_ASCII);
		check(phases.get(0).equals("phase\tcandidates\tprobed\tunavailable\taccepted\tslot_cost\tthreshold\tnet_gain\telapsed_ns\tcovered_before\tcovered_after\tresidual_before\tresidual_after"),
			"Cycle phase report uses the exact thirteen-column contract");
		final String[] split=phases.get(1).split("\t",-1),merge=phases.get(2).split("\t",-1);
		check(split.length==13 && merge.length==13 && split[0].equals("fission") && merge[0].equals("fusion"),"One cycle executes exactly fission then fusion");
		check(split[1].equals("2") && split[2].equals("2") && split[4].equals("1") && split[5].equals("0") && split[7].equals("2") &&
			split[9].equals("8") && split[10].equals("10") && split[11].equals("2") && split[12].equals("0"),
			"Fission resolves shared rescues before updating coverage and residuals");
		check(merge[9].equals(split[10]) && merge[11].equals(split[12]) && merge[5].equals("3"),
			"Fusion consumes committed fission state and recomputes B from its three active families");
		check(merge[1].equals("1") && merge[4].equals("1") && merge[7].equals("0") && merge[10].equals("10"),
			"Fusion evaluates current mutual parents and accepts a zero-gain merge with slot credit");
		final List<String> metrics=Files.readAllLines(root.resolve("cycle1.metrics.tsv"),StandardCharsets.US_ASCII);
		check(metrics.get(0).equals(LiveProbeMetrics.SUMMARY_HEADER),"Cycle metrics use the shared summary schema");
		final String[] stages={"before","after_fission","after_fusion"};
		final int[] covered={8,10,10},models={2,3,2},rescored={2,2,1},reused={0,1,1};
		for(int i=0;i<stages.length;i++){
			final String[] row=metrics.get(i+1).split("\t",-1);
			check(row.length==18 && row[0].equals(stages[i]) && Integer.parseInt(row[1])==models[i] &&
				Integer.parseInt(row[3])==covered[i] && Integer.parseInt(row[4])==covered[i],"Cycle before/after metrics match published populations");
			check(Integer.parseInt(row[14])==rescored[i] && Integer.parseInt(row[15])==reused[i],
				"Cycle quality measurements reuse unchanged Family objects and rescore changed ones");
			check(Long.parseLong(row[6])+Long.parseLong(row[7])==covered[i],"Available and missing quality denominators partition current coverage");
		}
		final List<String> histogram=Files.readAllLines(root.resolve("cycle1.sizes.tsv"),StandardCharsets.US_ASCII);
		check(histogram.get(0).equals(LiveProbeMetrics.HISTOGRAM_HEADER),"Cycle size histogram has the shared schema");
		for(int stage=0;stage<stages.length;stage++){
			long familyCount=0,tails=0;
			for(int row=1;row<histogram.size();row++){
				final String[] fields=histogram.get(row).split("\t",-1);
				check(fields.length==4,"Exact size histogram has four fields");
				if(fields[0].equals(stages[stage])){
					final long size=Long.parseLong(fields[1]),count=Long.parseLong(fields[2]),total=Long.parseLong(fields[3]);
					check(size*count==total,"Histogram tails equal exact size times count");familyCount+=count;tails+=total;
				}
			}
			check(familyCount==models[stage] && tails==covered[stage],"Histogram agrees with each metric snapshot");
		}
		final List<String> members=Files.readAllLines(root.resolve("cycle1.members.tsv"),StandardCharsets.US_ASCII);
		check(members.size()==11,"Final cycle publishes every original ordinal");
		for(int i=0;i<10;i++){check(members.get(i+1).startsWith(i+"\tcycle"+i+"\t") && members.get(i+1).endsWith("\tassigned"),"Final cycle retains both rescued originals");}
		new DdlClustStageB(cycleArgs(root,"cap2",1,"mergemaxsize=2")).process();
		final List<String> capped=Files.readAllLines(root.resolve("cap2.passes.tsv"),StandardCharsets.US_ASCII);
		check(capped.get(2).split("\t",-1)[1].equals("0") && ClusterSketchCodec.read(root.resolve("cap2.pddl").toString()).families.length==3,
			"Fusion cap filters current post-fission endpoint sizes without changing fission");
		final String[] cycleOutputs={"fusionreport","passesreport","metrics","sizehist"};
		for(int i=0;i<cycleOutputs.length;i++){
			final String prefix="missing"+i;
			expectFail("Missing cycle output","Live cycle requires",cycleArgs(root,prefix,1,cycleOutputs[i]+"=null"));
			expectFail("New cycle output aliases input","Stage-B input/output paths alias",cycleArgs(root,prefix,1,cycleOutputs[i]+"="+root.resolve("source.pddl"),"ow=t"));
			check(!Files.exists(root.resolve(prefix+".pddl")) && !Files.exists(root.resolve(prefix+".members.tsv")),"Cycle output preflight precedes all writers");
		}
		expectFail("Cycle output aliases another new output","Stage-B input/output paths alias",
			cycleArgs(root,"alias",1,"sizehist="+root.resolve("alias.metrics.tsv")));
	}

	private static String[] cycleArgs(final Path root,final String prefix,final int threads,final String...extra){
		final String[] fixed={"phase=cycle","selectedcount=2","minlen=100","minspan=70","minid=0","minblosumratio=0",
			"mergeminid=0","mergeminratio=0","fissionminsize=4","familybudget=3","mergemaxsize=19",
			"fusionreport="+root.resolve(prefix+".fusion.tsv"),"passesreport="+root.resolve(prefix+".passes.tsv"),
			"metrics="+root.resolve(prefix+".metrics.tsv"),"sizehist="+root.resolve(prefix+".sizes.tsv")};
		final String[] options=Arrays.copyOf(fixed,fixed.length+extra.length);
		System.arraycopy(extra,0,options,fixed.length,extra.length);
		return liveArgs(root,prefix,threads,options);
	}

	private static void cycleTableParity(final Path one,final Path four,final int timeColumn,final int expectedRows) throws Exception{
		final List<String> a=Files.readAllLines(one,StandardCharsets.US_ASCII),b=Files.readAllLines(four,StandardCharsets.US_ASCII);
		check(a.size()==expectedRows && b.size()==expectedRows && a.get(0).equals(b.get(0)),"Cycle table row/header parity");
		for(int i=1;i<a.size();i++){
			final String[] x=a.get(i).split("\t",-1),y=b.get(i).split("\t",-1);
			check(x.length==y.length && timeColumn<x.length,"Cycle metric column count");
			for(int col=0;col<x.length;col++){
				if(col==timeColumn){check(Long.parseLong(x[col])>=0 && Long.parseLong(y[col])>=0,"Cycle elapsed time is nonnegative");}
				else{check(x[col].equals(y[col]),"Cycle scientific fields are worker-invariant: "+col);}
			}
		}
	}

	/** Current-state exports preserve zero-member models, duplicate names and all original statuses. */
	private static void liveProbe(final Path root) throws Exception{
		Files.createDirectory(root);
		final String sequence=repeat("ACDEFGHIKLMNPQRSTVWY",4);
		final ProteinSequence[] input={new ProteinSequence("duplicate",sequence),new ProteinSequence("duplicate",sequence),
			new ProteinSequence("rescue",sequence),new ProteinSequence("residual",repeat("W",80)),new ProteinSequence("short","ACDE")};
		final FamilySketchBuilder builder=new FamilySketchBuilder(64,2,7,false,0,0,null);
		final Family[] families={family(10,"a",input,builder,0),family(20,"b",input,builder,1),
			seedFamily(30,"empty",repeat("C",80),1,builder)};
		final Path source=root.resolve("source.pddl"),selected=root.resolve("selected.tsv"),baseline=root.resolve("baseline.tsv");
		ClusterSketchCodec.write(source.toString(),false,families,2,64,7,false,0,0);
		final Family[] rankedFamilies={family(10,"a",input,builder,0),family(20,"b",input,builder,1),
			seedFamily(30,"empty",repeat("C",80),1,builder)};
		final MemberRankInputs rankingInputs=new MemberRankInputs(input);
		for(int i=0;i<rankedFamilies.length;i++){rankedFamilies[i].enableMemberRanking(rankingInputs,2+i,null);}
		final Path rankedSource=root.resolve("source.ranked.pddl");
		ClusterSketchCodec.write(rankedSource.toString(),false,rankedFamilies,2,64,7,false,0,0);
		fasta(root.resolve("input.faa"),input);
		write(selected,"family_id\tsource_members\n10\t1\n20\t1\n30\t1\n");
		final String originalMembers="input_index\tinput_id\tfamily_id\tstatus\n"+
			"0\tduplicate\t10\tassigned\n1\tduplicate\t20\tassigned\n2\trescue\t-1\tunassigned\n"+
			"3\tresidual\t-1\tunassigned\n4\tshort\t-1\tfiltered_length\n";
		write(baseline,originalMembers);
		final byte[] originalCodec=Files.readAllBytes(source);
		for(int threads:new int[]{1,4}){
			new DdlClustStageB(liveArgs(root,"rescue"+threads,threads)).process();
		}
		for(String suffix:new String[]{".pddl",".members.tsv",".report.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("rescue1"+suffix)),Files.readAllBytes(root.resolve("rescue4"+suffix))),
				"Live pass deterministic serial/parallel export: "+suffix);
		}
		final List<String> members=Files.readAllLines(root.resolve("rescue1.members.tsv"),StandardCharsets.US_ASCII);
		check(members.size()==input.length+1 && members.get(0).equals("input_index\tinput_id\tfamily_id\tstatus"),
			"Live export includes every original ordinal and explicit status");
		final int winner=Integer.parseInt(members.get(1).split("\t")[2]);
		check(winner>=0,"Rescue pass retains an assigned first parent");
		check(members.get(1).equals("0\tduplicate\t"+winner+"\tassigned") &&
			members.get(2).equals("1\tduplicate\t"+winner+"\tassigned") &&
			members.get(3).equals("2\trescue\t"+winner+"\tassigned"),
			"Live fusion publishes actual winners, including rescue, without collapsing duplicate input names");
		check(members.get(4).equals("3\tresidual\t-1\tunassigned") && members.get(5).equals("4\tshort\t-1\tfiltered_length"),
			"Live export distinguishes eligible residual from length-filtered input");
		final ClusterSketchCodec.Loaded loaded=ClusterSketchCodec.read(root.resolve("rescue1.pddl").toString());
		check(loaded.families.length==2,"One fusion replaces its two parents and preserves the zero-member model");
		for(ClusterSketchCodec.Entry entry:loaded.families){
			if(entry.familyId==winner){
				check(entry.members==3 && decode(entry.consensus).equals(sequence),"Published model count includes rescue and retains evaluated consensus");
			}else{
				check(entry.familyId==30 && entry.members==0 && Arrays.equals(entry.consensus,families[2].consensus()) &&
					Arrays.equals(entry.sketch,families[2].sketch()) && Arrays.equals(entry.scalars,families[2].scalars()),
					"Full current roster retains untouched zero-member model payload");
			}
		}
		liveSummary(root.resolve("rescue1.report.tsv"),"phase=fusion\tcandidates=1\tavailable=1\taccepted=1\tslot_cost=1\tthreshold=-1"+
			"\tcovered_before=2\tcovered_after=3\tresidual=1\tfiltered=1\tmodels=2\tactive_models=1\tnext_family_id=31\tnet_gain=1"+
			"\tempty_child_rejected=0");
		final String[] decision=liveDecision(root.resolve("rescue1.report.tsv"));
		check(decision[1].equals("fusion") && decision[5].equals("1") && Double.parseDouble(decision[6])==-1 &&
			decision[10].equals("1") && decision[11].equals("1") && decision[13].equals("1"),
			"Live report uses merge credit -B, measured rescue gain and actual accepted marginal gain");
		new DdlClustStageB(liveArgs(root,"fission",1,"phase=fission","localrefinerounds=1")).process();
		liveSummary(root.resolve("fission.report.tsv"),"phase=fission\tcandidates=0\tavailable=0\taccepted=0\tslot_cost=1\tthreshold=1"+
			"\tcovered_before=2\tcovered_after=2\tresidual=2\tfiltered=1\tmodels=3\tactive_models=2\tnext_family_id=31\tnet_gain=0"+
			"\tempty_child_rejected=0");
		check(new String(Files.readAllBytes(root.resolve("fission.members.tsv")),StandardCharsets.US_ASCII).equals(originalMembers),
			"Explicit fission phase with no eligible parent is a valid no-op pass");
		new DdlClustStageB(liveArgs(root,"ranked",1,"memberranking=t","codecs="+rankedSource)).process();
		new DdlClustStageB(defaultRanking(liveArgs(root,"default_ranked",1,"codecs="+rankedSource))).process();
		new DdlClustStageB(liveArgs(root,"ranked_off",1,"codecs="+rankedSource)).process();
		for(String suffix:new String[]{".pddl",".members.tsv",".report.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("ranked"+suffix)),Files.readAllBytes(root.resolve("default_ranked"+suffix))),
				"File liveprobe defaults to explicit ON on the identical ranked codec");
			check(Arrays.equals(Files.readAllBytes(root.resolve("rescue1"+suffix)),Files.readAllBytes(root.resolve("ranked_off"+suffix))),
				"Explicit OFF ignores descriptors while preserving legacy liveprobe output bytes");
		}
		boolean rankedWinner=false,rankedEmpty=false;
		for(ClusterSketchCodec.Entry entry:ClusterSketchCodec.read(root.resolve("ranked.pddl").toString()).families){
			check(entry.memberRanking!=null,"Explicit ranked live hydration writes every retained descriptor");
			if(entry.familyId==winner){
				check(entry.memberRanking.members==3,"Ranked live hydration measures the published current winners");
				rankedWinner=true;
			}else if(entry.familyId==30){
				check(entry.memberRanking.members==1 && entry.memberRanking.layer==4,
					"Ranked live hydration preserves empty last-nonempty provenance");
				rankedEmpty=true;
			}
		}
		check(rankedWinner && rankedEmpty,"Ranked live hydration keeps both the changed winner and unchanged empty template");
		// Reopen a complete current state, including its zero-member template, under the default.
		final String[] current=defaultRanking(liveArgs(root,"default_current",1,"phase=layer","currentstate=t",
			"nextfamilyid=31","selectedcount=2","selectedseeds=null","codecs="+root.resolve("ranked.pddl"),
			"baselinemembers="+root.resolve("ranked.members.tsv"),"fissionminsize=1","fissiontopfraction=0",
			"fusionbottomfraction=0","pricemode=fraction","stateout="+root.resolve("default_current.state.tsv"),
			"fusionreport="+root.resolve("default_current.fusion.tsv"),"passesreport="+root.resolve("default_current.passes.tsv"),
			"metrics="+root.resolve("default_current.metrics.tsv"),"sizehist="+root.resolve("default_current.sizes.tsv")));
		new DdlClustStageB(current).process();
		check(Arrays.equals(Files.readAllBytes(root.resolve("ranked.pddl")),Files.readAllBytes(root.resolve("default_current.pddl"))) &&
			Arrays.equals(Files.readAllBytes(root.resolve("ranked.members.tsv")),Files.readAllBytes(root.resolve("default_current.members.tsv"))),
			"Default ranked current-state no-op retains exact model, descriptor and membership bytes");
		expectFail("Legacy source lacks member ranking","Ranked hydration requires serialized member descriptors",
			liveArgs(root,"unranked_ranked",1,"memberranking=t"));
		expectFail("Default ON cannot silently fall back","Ranked hydration requires serialized member descriptors",
			defaultRanking(liveArgs(root,"unranked_default",1)));
		residentRankingModes(root,rankedSource,builder);
		check(Arrays.equals(originalCodec,Files.readAllBytes(source)) &&
			new String(Files.readAllBytes(baseline),StandardCharsets.US_ASCII).equals(originalMembers),"Live passes leave source inputs unchanged");
		liveInvalid(root);
	}

	/** A signed merge credit can accept a measured loss; those ordinals must become residuals. */
	private static void liveProbeLoss(final Path root) throws Exception{
		Files.createDirectory(root);
		final String r=repeat("A",70)+repeat("R",30),k=repeat("A",70)+repeat("K",30);
		final ProteinSequence[] input=new ProteinSequence[7];
		final StringBuilder membership=new StringBuilder("input_index\tinput_id\tfamily_id\tstatus\n");
		for(int i=0;i<input.length;i++){
			input[i]=new ProteinSequence("loss"+i,i<2 ? r : i<6 ? k : repeat("W",100));
			membership.append(i).append("\tloss").append(i).append(i<2 ? "\t100\tassigned\n" :
				i<6 ? "\t300\tassigned\n" : "\t-1\tunassigned\n");
		}
		final FamilySketchBuilder builder=new FamilySketchBuilder(64,2,7,false,0,0,null);
		final Family a=new Family(100,"r",0,input[0].enc),b=new Family(300,"k",2,input[2].enc);
		a.addMember(1);for(int i=3;i<6;i++){b.addMember(i);}
		a.refreshStats(input,builder);b.refreshStats(input,builder);
		ClusterSketchCodec.write(root.resolve("source.pddl").toString(),false,new Family[]{a,b},2,64,7,false,0,0);
		fasta(root.resolve("input.faa"),input);write(root.resolve("baseline.tsv"),membership.toString());
		write(root.resolve("selected.tsv"),"family_id\tsource_members\n100\t2\n300\t4\n");
		new DdlClustStageB(liveArgs(root,"loss",4,"selectedcount=2","minlen=100","minspan=70","minid=.8","mergeminid=.5")).process();
		final List<String> members=Files.readAllLines(root.resolve("loss.members.tsv"),StandardCharsets.US_ASCII);
		check(members.size()==8 && members.get(1).equals("0\tloss0\t-1\tunassigned") &&
			members.get(2).equals("1\tloss1\t-1\tunassigned") && members.get(7).equals("6\tloss6\t-1\tunassigned"),
			"Accepted live fusion writes both lost originals and previous residual to the next pool");
		final ClusterSketchCodec.Loaded result=ClusterSketchCodec.read(root.resolve("loss.pddl").toString());
		check(result.families.length==1 && result.families[0].members==4 && decode(result.families[0].consensus).equals(k),
			"Published membership count may shrink without replacing the evaluated consensus");
		for(int i=2;i<6;i++){
			check(members.get(i+1).equals(i+"\tloss"+i+"\t"+result.families[0].familyId+"\tassigned"),"Retained winner ordinal survives live export");
		}
		liveSummary(root.resolve("loss.report.tsv"),"phase=fusion\tcandidates=1\tavailable=1\taccepted=1\tslot_cost=4\tthreshold=-4"+
			"\tcovered_before=6\tcovered_after=4\tresidual=3\tfiltered=0\tmodels=1\tactive_models=1\tnext_family_id=301\tnet_gain=-2"+
			"\tempty_child_rejected=0");
		final String[] decision=liveDecision(root.resolve("loss.report.tsv"));
		check(decision[5].equals("4") && Double.parseDouble(decision[6])==-4 && decision[10].equals("-2") &&
			decision[11].equals("1") && decision[13].equals("-2"),"Signed -B credit accepts measured loss -2 above -4");
		new DdlClustStageB(liveArgs(root,"blocked",1,"selectedcount=2","minlen=100","minspan=70",
			"minid=.8","mergeminid=.5","slotfraction=.5")).process();
		check(new String(Files.readAllBytes(root.resolve("blocked.members.tsv")),StandardCharsets.US_ASCII).equals(membership.toString()),
			"Gain equal to the fractional signed threshold leaves the covered/residual partition unchanged");
		check(ClusterSketchCodec.read(root.resolve("blocked.pddl").toString()).families.length==2 &&
			liveDecision(root.resolve("blocked.report.tsv"))[11].equals("0"),"No-accept pass exports the full unchanged roster");
		liveSummary(root.resolve("blocked.report.tsv"),"phase=fusion\tcandidates=1\tavailable=1\taccepted=0\tslot_cost=4\tthreshold=-2"+
			"\tcovered_before=6\tcovered_after=6\tresidual=1\tfiltered=0\tmodels=2\tactive_models=2\tnext_family_id=301\tnet_gain=0"+
			"\tempty_child_rejected=0");
	}

	private static String[] liveArgs(final Path root,final String prefix,final int threads,final String...extra){
		final String[] base={"mode=liveprobe","memberranking=f","phase=fusion","gainmargin=0","in="+root.resolve("input.faa"),
			"codecs="+root.resolve("source.pddl"),"selectedseeds="+root.resolve("selected.tsv"),
			"baselinemembers="+root.resolve("baseline.tsv"),"out="+root.resolve(prefix+".pddl"),
			"members="+root.resolve(prefix+".members.tsv"),"report="+root.resolve(prefix+".report.tsv"),
			"selectedcount=3","minlen=80","minspan=40","maxalignlen=100","trimleft=f","consensustrimleft=f",
			"minid=.35","minblosumratio=.14","mergeminid=.35","mergeminratio=.14","passhits=0","maxhydrodelta=9",
			"fissionminsize=10","minchildsize=1","memberneighbors=4","memberminhits=1","evidenceneighbors=10",
			"familybudget=1","t="+threads,"frozenindexthreads="+threads,"ow=f"};
		final String[] args=Arrays.copyOf(base,base.length+extra.length);
		System.arraycopy(extra,0,args,base.length,extra.length);
		return args;
	}

	private static String[] defaultRanking(final String[] args){
		final ArrayList<String> implicit=new ArrayList<String>(Arrays.asList(args));
		check(implicit.remove("memberranking=f"),"Default witness must remove its explicit legacy control");
		return implicit.toArray(new String[implicit.size()]);
	}

	/** Reject mismatched flags/cache before a resident layer can open its reports. */
	private static void residentRankingModes(final Path root,final Path source,final FamilySketchBuilder builder)throws Exception{
		final ClusteringInput population=ClusteringInput.load(root.resolve("input.faa").toString(),-1,false);
		final MemberRankInputs cache=new MemberRankInputs(population.sequences);
		final ClusterSketchCodec.Loaded codec=ClusterSketchCodec.read(source.toString());
		final String selected=root.resolve("selected.tsv").toString(),members=root.resolve("baseline.tsv").toString();
		final LocalProbeBaseline rankedBaseline=LocalProbeBaseline.load(codec,selected,members,population.sequences,80,3,cache);
		final LiveProbeState unboundRanked=LiveProbeState.fromBaseline(rankedBaseline,population.sequences.length);
		final LiveProbeState ranked=unboundRanked.withMemberRanking(population.sequences,cache,rankedBaseline.memberRankingLayer);
		final LiveProbeState legacy=LiveProbeState.fromBaseline(LocalProbeBaseline.load(codec,selected,members,population.sequences,80,3),
			population.sequences.length);
		final ClusteringInput other=ClusteringInput.load(root.resolve("input.faa").toString(),-1,false);
		final java.lang.reflect.Method run=DdlClustStageB.class.getDeclaredMethod("runResidentLayer",String[].class,
			LiveProbeState.class,ClusteringInput.class,FamilySketchBuilder.class,int[].class);
		run.setAccessible(true);
		for(int i=0;i<4;i++){
			final String[] args={"mode=liveprobe","phase=layer","pricemode=fraction","memberranking="+(i==1 || i==2),"t=1",
				"selectedcount=3","stateout="+root.resolve("mismatch"+i+".state"),
				"passesreport="+root.resolve("mismatch"+i+".passes"),"metrics="+root.resolve("mismatch"+i+".metrics"),
				"sizehist="+root.resolve("mismatch"+i+".sizes")};
			boolean rejected=false;
			try{run.invoke(null,args,i==1 ? legacy : i==3 ? unboundRanked : ranked,i==2 ? other : population,builder,null);}
			catch(java.lang.reflect.InvocationTargetException expected){
				check(expected.getCause() instanceof IllegalArgumentException,"Resident mismatch is a clear argument failure");
				check(expected.getCause().getMessage().contains(i==2 ? "identical original population/cache" : "effective ranking mode"),
					"Resident mismatch fails for the intended cache/mode reason");
				rejected=true;
			}
			check(rejected && !Files.exists(root.resolve("mismatch"+i+".state")),"Bad resident bindings cannot publish reports");
		}
	}

	/** Summary is mandatory even for zero proposals, and contains no nondeterministic timing. */
	private static void liveSummary(final Path path,final String expected) throws Exception{
		final List<String> rows=Files.readAllLines(path,StandardCharsets.US_ASCII);
		final String[] wanted=("#live_probe_pass\t"+expected).split("\t",-1),actual=rows.get(0).split("\t",-1);
		check(wanted.length==16 && actual.length==wanted.length,"Live summary has exactly fifteen named fields");
		for(int i=0;i<wanted.length;i++){
			if(wanted[i].startsWith("threshold=")){
				check(actual[i].startsWith("threshold=") && Double.parseDouble(actual[i].substring(10))==
					Double.parseDouble(wanted[i].substring(10)),"Live summary signed threshold is exact");
			}else{check(actual[i].equals(wanted[i]),"Live summary field: "+wanted[i]);}
		}
		check(rows.get(1).equals("request_id\tphase\tparent_a\tparent_b\tchild_id\tslot_cost\tthreshold\tavailable\treason\t"+
			"strict_member_veto\tmeasured_gain\taccepted\tcommit_order\tmarginal_gain\trescued\tretained\taffected\tforced_retention\tempty_child_veto"),
			"Live decision header immediately follows its summary");
		check(rows.size()==2+Integer.parseInt(actual[2].substring("candidates=".length())),
			"Live report emits every candidate once and preserves summary/header for an empty pass");
	}

	private static String[] liveDecision(final Path path) throws Exception{
		final List<String> rows=Files.readAllLines(path,StandardCharsets.US_ASCII);
		for(String row:rows){
			if(!row.isEmpty() && Character.isDigit(row.charAt(0))){
				final String[] fields=row.split("\t",-1);
				check(fields.length==19 && fields[18].equals("0"),"Live decision schema retains measured and marginal gains separately with an explicit empty-child flag");
				return fields;
			}
		}
		throw new AssertionError("Expected one measured live proposal: "+path);
	}

	private static void liveInvalid(final Path root) throws Exception{
		final String[][] bad={{"phase=null","requires explicit phase"},{"phase=both","requires explicit phase"},
			{"metrics=unused.tsv","fusionreport, metrics and sizehist require"},
			{"gainmargin=-1","symmetric slotfraction"},{"gainmargin=NaN","symmetric slotfraction"},
			{"gainmargin=Infinity","symmetric slotfraction"},{"cohesionscale=1","Live probes do not accept"},
			{"separationscale=1","Live probes do not accept"},{"sizescale=1","Live probes do not accept"},
			{"schedule=none","Live probes do not accept"},{"temperature=1","Live probes do not accept"},
			{"probes=unused.tsv","Live probes do not accept"},{"probeseed=1","Random request counts/seed"},
			{"localrefinerounds=0","Live probes require localrefinerounds"},
			{"localrefinerounds=3","Live probes require localrefinerounds"},
			{"refinementreport=unused.tsv","refinementreport requires mode=localprobe"},
			{"members=null","mode=liveprobe requires"},{"report=null","mode=liveprobe requires"}};
		for(int i=0;i<bad.length;i++){
			final String prefix="invalid"+i;
			expectFail("Invalid live flag "+bad[i][0],bad[i][1],liveArgs(root,prefix,1,bad[i][0]));
			check(!Files.exists(root.resolve(prefix+".pddl")) && !Files.exists(root.resolve(prefix+".members.tsv")) &&
				!Files.exists(root.resolve(prefix+".report.tsv")),"Invalid live flags fail before output writers");
		}
		final ArrayList<String> missingMargin=new ArrayList<String>(Arrays.asList(liveArgs(root,"missing_margin",1)));
		check(missingMargin.remove("gainmargin=0"),"Missing-margin fixture removes its explicit gainmargin argument");
		new DdlClustStageB(missingMargin.toArray(new String[missingMargin.size()]));checks++; // Default fraction1 needs no legacy additive margin.
		final Path hard=root.resolve("source.hard.pddl"),symbolic=root.resolve("baseline.sym.tsv");
		Files.createLink(hard,root.resolve("source.pddl"));Files.createSymbolicLink(symbolic,root.resolve("baseline.tsv").toAbsolutePath());
		final byte[] codec=Files.readAllBytes(hard),baseline=Files.readAllBytes(symbolic);
		expectFail("Hardlinked live codec output","Stage-B input/output paths alias",liveArgs(root,"hard",1,"out="+hard,"ow=t"));
		expectFail("Symlinked live membership output","Stage-B input/output paths alias",liveArgs(root,"symbolic",1,"members="+symbolic,"ow=t"));
		expectFail("Live outputs alias each other","Stage-B input/output paths alias",
			liveArgs(root,"collision",1,"report="+root.resolve("collision.members.tsv")));
		liveOverwrite(root);
		check(Arrays.equals(codec,Files.readAllBytes(hard)) && Arrays.equals(baseline,Files.readAllBytes(symbolic)),
			"Live alias preflight preserves input bytes");
		check(!Files.exists(root.resolve("hard.report.tsv")) && !Files.exists(root.resolve("symbolic.pddl")) &&
			!Files.exists(root.resolve("collision.pddl")),"Live alias preflight does not publish partial output");
	}

	/** Tools.testOutputFiles asserts under -ea and returns false under -da; both must preserve outputs. */
	private static void liveOverwrite(final Path root) throws Exception{
		final Path[] outputs={root.resolve("rescue1.pddl"),root.resolve("rescue1.members.tsv"),root.resolve("rescue1.report.tsv")};
		final byte[][] before=new byte[outputs.length][];
		for(int i=0;i<outputs.length;i++){before[i]=Files.readAllBytes(outputs[i]);}
		boolean rejected=false;
		try{new DdlClustStageB(liveArgs(root,"rescue1",1)).process();}
		catch(AssertionError failure){
			check(("File "+outputs[0]+" exists and overwrite=false").equals(failure.getMessage()),
				"Only the specific Tools.testOutputFiles overwrite assertion satisfies this negative fixture");
			rejected=true;
		}catch(IllegalArgumentException failure){
			check("Live probe outputs are unavailable".equals(failure.getMessage()),
				"With assertions disabled the CLI must reject the failed output preflight");
			rejected=true;
		}
		check(rejected,"Existing live outputs must reject overwrite=false before any writer opens");
		for(int i=0;i<outputs.length;i++){
			check(Arrays.equals(before[i],Files.readAllBytes(outputs[i])),"Overwrite rejection preserves exact bytes: "+outputs[i]);
		}
	}

	/** Five identical singleton families require repeated reciprocal passes under stable ties.
	 * Their original ordinals must survive exactly once, independent of worker count.
	 */
	private static void mutualCase(Path root) throws Exception{
		Files.createDirectory(root);
		final ProteinSequence[] input=new ProteinSequence[5];
		final Family[] families=new Family[input.length];
		final FamilySketchBuilder builder=new FamilySketchBuilder(64, 2, 7, false, 0, 0, null);
		final String sequence=repeat("ACDEFGHIKLMNPQRSTVWY", 4);
		final StringBuilder map=new StringBuilder(ClusteringShardWriter.MAPPING_HEADER).append('\n');
		final StringBuilder membership=new StringBuilder("input_index\tinput_id\tfamily_id\tstatus\n");
		for(int i=0; i<input.length; i++){
			input[i]=new ProteinSequence("original"+i, sequence);
			families[i]=family(10+3*i, "source"+i, input, builder, i);
			map.append(i).append("\t0\t").append(i).append("\toriginal").append(i).append('\n');
			membership.append(i).append("\toriginal").append(i).append('\t').append(10+3*i).append("\tassigned\n");
		}
		final Path source=root.resolve("source.pddl"), fasta=root.resolve("input.faa"),
			parent=root.resolve("map.tsv"), members=root.resolve("source.members.tsv"), roster=root.resolve("roster.tsv");
		ClusterSketchCodec.write(source.toString(), false, families, 2, 64, 7, false, 0, 0);
		fasta(fasta, input); write(parent, map.toString()); write(members, membership.toString());
		new DdlClustStageB(new String[]{"mode=collect", "codecs="+source, "out="+root.resolve("collected.faa"),
			"roster="+roster, "ow=f"}).process();
		for(int threads : new int[]{1,4}){
			new DdlClustStageB(new String[]{"mode=mutual", "memberranking=f", "in="+fasta, "out="+root.resolve("merged"+threads+".pddl"),
				"roster="+roster, "partitionmap="+parent, "stageacodecs="+source, "stageamembers="+members,
				"report="+root.resolve("pairs"+threads+".tsv"), "passesreport="+root.resolve("passes"+threads+".tsv"),
				"members="+root.resolve("members"+threads+".tsv"), "minspan=5", "mergemaxsize=20",
				"shortlist=50", "t="+threads, "ow=f"}).process();
		}
		for(String pattern : new String[]{"merged#.pddl", "pairs#.tsv", "members#.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve(pattern.replace("#", "1"))),
				Files.readAllBytes(root.resolve(pattern.replace("#", "4")))), "Mutual merge worker parity: "+pattern);
		}
		final ClusterSketchCodec.Loaded result=ClusterSketchCodec.read(root.resolve("merged1.pddl").toString());
		check(result.families.length==1 && result.families[0].familyId==0 && result.families[0].members==5,
			"Repeated reciprocal passes conserve five original tails in the minimum stable-ID family");
		check(decode(result.families[0].consensus).equals(sequence), "Identical original tails retain their consensus");
		final List<String> pairs=Files.readAllLines(root.resolve("pairs1.tsv"), StandardCharsets.US_ASCII);
		check(pairs.size()==5, "Five singleton families need four disjoint pair unions in total");
		final List<String> passes=Files.readAllLines(root.resolve("passes1.tsv"), StandardCharsets.US_ASCII);
		check(passes.size()==6, "Stable ties give four reducing passes and one terminating no-pair pass");
		final List<String> parallel=Files.readAllLines(root.resolve("passes4.tsv"), StandardCharsets.US_ASCII);
		check(parallel.size()==passes.size(), "Mutual iteration has the same pass count with four workers");
		for(int i=1; i<passes.size(); i++){
			final String serial=passes.get(i), multi=parallel.get(i);
			check(serial.substring(0, serial.lastIndexOf('\t')).equals(multi.substring(0, multi.lastIndexOf('\t'))),
				"Mutual pass counters match across workers, excluding measured elapsed time");
		}
		final List<String> owners=Files.readAllLines(root.resolve("members1.tsv"), StandardCharsets.US_ASCII);
		check(owners.size()==6, "Final union writes every original assigned ordinal once");
		for(int i=0; i<input.length; i++){
			check(owners.get(i+1).equals(i+"\toriginal"+i+"\t0"), "Original ordinal and ID survive repeated unions");
		}
		frozenIndexMutual(root,source,fasta,parent,members,roster);
		structuralModes(root,source,fasta,parent,members,roster);
		probeBatch(root,source,fasta,parent,members,roster);
		new DdlClustStageB(new String[]{"mode=mutual", "memberranking=f", "in="+fasta, "out="+root.resolve("blocked.pddl"),
			"roster="+roster, "partitionmap="+parent, "stageacodecs="+source, "stageamembers="+members,
			"report="+root.resolve("blocked.pairs.tsv"), "passesreport="+root.resolve("blocked.passes.tsv"),
			"members="+root.resolve("blocked.members.tsv"), "minspan=5", "mergemaxsize=2", "t=4", "ow=f"}).process();
		final ClusterSketchCodec.Loaded blocked=ClusterSketchCodec.read(root.resolve("blocked.pddl").toString());
		check(blocked.families.length==4 && blocked.families[0].members==2,
			"A regenerated family at the threshold blocks the remaining identical small families in the next pass");
		check(Files.readAllLines(root.resolve("blocked.pairs.tsv"), StandardCharsets.US_ASCII).size()==2,
			"No second-best small pairing may bypass the regenerated large blocker");
	}

	/** Repeated rows are independent edits; rejected rows explicitly replace stale codec payloads. */
	private static void probeBatch(final Path root,final Path source,final Path fasta,final Path parent,
			final Path members,final Path roster) throws Exception{
		final Path requests=root.resolve("requests.tsv");
		final String seed=Long.toString(Long.MIN_VALUE);
		write(requests,StructuralProbeBatch.REQUEST_HEADER+"\n0\tfusion\t0\t1\t"+seed+
			"\n1\tfusion\t1\t0\t"+seed+"\n2\tfusion\t0\t2\t"+seed+"\n3\tfission\t0\t-1\t"+seed+"\n");
		for(int threads:new int[]{1,4}){
			final String prefix="probe_t"+threads;
			Files.copy(source,root.resolve(prefix+".2.pddl"));
			new DdlClustStageB(probeArgs(root,prefix,source,fasta,parent,members,roster,requests,threads)).process();
			for(int row=0; row<4; row++){
				final ClusterSketchCodec.Loaded codec=ClusterSketchCodec.read(root.resolve(prefix+"."+row+".pddl").toString());
				check(codec.families.length==(row<2 ? 1 : 0),"Feasible probes export one union; rejected probes clear stale models");
				if(row<2){check(codec.families[0].familyId==0 && codec.families[0].members==2,"Each row makes only the same original reciprocal pair edit");}
			}
			check(Arrays.equals(Files.readAllBytes(root.resolve(prefix+".0.pddl")),Files.readAllBytes(root.resolve(prefix+".1.pddl"))),
				"Repeated reciprocal probes do not merge their previous result into the initial population");
			final List<String> results=Files.readAllLines(root.resolve(prefix+".results.tsv"),StandardCharsets.US_ASCII);
			check(results.size()==5 && results.get(0).equals(StructuralProbeBatch.RESULT_HEADER),"Probe table has its complete declared schema and every request");
			for(int i=1; i<results.size(); i++){
				final String[] fields=results.get(i).split("\t",-1);
				check(fields.length==31 && fields[1].equals(seed),"Probe metadata retains the entire signed seed without wrapping");
				check(fields[3].equals(i<3 ? "1" : "0"),"Proposal feasibility remains separate from objective acceptance");
				if(i<3){check(fields[7].equals("0,1") && fields[23].equals("1.0") && fields[30].equals("2"),"Fusion lineage and slot credit use original member counts");}
				else{check(fields[19].equals("NA") && fields[24].equals("NA"),"Rejected proposals never invent numeric quality or gain");}
			}
			final List<String> lineage=Files.readAllLines(root.resolve(prefix+".members.tsv"),StandardCharsets.US_ASCII);
			check(lineage.size()==5 && lineage.get(1).equals("0\t0\toriginal0\t0") && lineage.get(4).equals("1\t1\toriginal1\t0"),
				"Only affected original tails are exported, separately for each successful probe");
			final List<String> work=Files.readAllLines(root.resolve(prefix+".work.tsv"),StandardCharsets.US_ASCII);
			check(work.size()==5 && work.get(0).equals(StructuralProbeBatch.WORK_HEADER),"Every request retains scoped native work counters");
			check(Long.parseLong(work.get(1).split("\t")[5])>0 && work.get(2).split("\t")[5].equals("0"),
				"Repeated endpoint requests reuse cached nearest choices within the batch");
		}
		for(String suffix:new String[]{".0.pddl",".1.pddl",".2.pddl",".3.pddl",".members.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("probe_t1"+suffix)),Files.readAllBytes(root.resolve("probe_t4"+suffix))),
				"Probe scientific bytes are independent of worker count: "+suffix);
		}
		final String[] invalidRows={"0\tfusion\t0\t99\t7\n","0\tfusion\t0\t1\t7\n0\tfusion\t0\t1\t7\n",
			"0\tfusion\t4294967296\t1\t7\n","0\tfusion\t0\t1\t9223372036854775808\n"};
		final String[] errors={"endpoint is absent","Probe IDs must","identifier is outside","probe numeric field"};
		for(int i=0; i<invalidRows.length; i++){
			final Path bad=root.resolve("bad_probe_"+i+".tsv");write(bad,StructuralProbeBatch.REQUEST_HEADER+"\n"+invalidRows[i]);
			expectFail("Invalid request cannot create partial probe outputs",errors[i],probeArgs(root,"bad_probe_"+i,source,fasta,parent,members,roster,bad,1));
			check(!Files.exists(root.resolve("bad_probe_"+i+".results.tsv")),"Request validation precedes output writers");
		}
		for(String flag:new String[]{"gainmargin=1","maxcycles=1","schedule=none"}){
			final String[] base=probeArgs(root,"bad_probe_flag",source,fasta,parent,members,roster,requests,1);
			final String[] bad=Arrays.copyOf(base,base.length+1);bad[base.length]=flag;
			expectFail("Raw proposals cannot use scheduling or fitted weights","Raw probes do not accept",bad);
		}
		final Path single=root.resolve("single_probe.tsv");write(single,StructuralProbeBatch.REQUEST_HEADER+"\n0\tfusion\t0\t1\t7\n");
		final String[] base=probeArgs(root,"alias_probe",source,fasta,parent,members,roster,single,1);
		final String[] alias=Arrays.copyOf(base,base.length+1);alias[base.length]="out="+source;
		final byte[] original=Files.readAllBytes(source);
		expectFail("Expanded codec path cannot alias an input","Stage-B input/output paths alias",alias);
		check(Arrays.equals(original,Files.readAllBytes(source)),"Probe preflight preserves the original source codec");
	}

	/** A legal split is exported even when its raw features would lose under the provisional objective. */
	private static void fissionProbe(final Path root) throws Exception{
		Files.createDirectory(root);
		final ProteinSequence[] input=new ProteinSequence[12];
		for(int i=0; i<input.length; i++){input[i]=new ProteinSequence("p"+i,probeSequence(i/6));}
		final FamilySketchBuilder builder=new FamilySketchBuilder(64,2,7,false,0,0,null);
		final Family original=new Family(9,"parent",0,input[0].enc),global=new Family(0,"s0_f9",0,input[0].enc);
		final StringBuilder map=new StringBuilder(ClusteringShardWriter.MAPPING_HEADER).append('\n');
		final StringBuilder members=new StringBuilder("input_index\tinput_id\tfamily_id\tstatus\n");
		for(int i=0; i<input.length; i++){
			if(i>0){original.addMember(i);global.addMember(i);}
			map.append(i).append("\t0\t").append(i).append("\tp").append(i).append('\n');
			members.append(i).append("\tp").append(i).append("\t9\tassigned\n");
		}
		original.refreshStats(input,builder);global.refreshStats(input,builder);
		final Path source=root.resolve("source.pddl"),fasta=root.resolve("input.faa"),roster=root.resolve("roster.tsv"),
			parent=root.resolve("map.tsv"),membership=root.resolve("members.tsv"),requests=root.resolve("requests.tsv"),base=root.resolve("global.pddl");
		ClusterSketchCodec.write(source.toString(),false,new Family[]{original},2,64,7,false,0,0);
		ClusterSketchCodec.write(base.toString(),false,new Family[]{global},2,64,7,false,0,0);
		fasta(fasta,input);write(parent,map.toString());write(membership,members.toString());
		write(requests,StructuralProbeBatch.REQUEST_HEADER+"\n20\tfission\t0\t-1\t7\n21\tfission\t0\t-1\t7\n");
		new DdlClustStageB(new String[]{"mode=collect","codecs="+source,"out="+root.resolve("collected.faa"),"roster="+roster,"ow=f"}).process();
		final String[] args=probeArgs(root,"split",source,fasta,parent,membership,roster,requests,4);
		final String[] split=Arrays.copyOf(args,args.length+6);
		final String[] controls={"fissionminsize=10","minchildsize=4","memberneighbors=4","familybudget=1","trimleft=f","minspan=40"};
		System.arraycopy(controls,0,split,args.length,controls.length);
		new DdlClustStageB(split).process();
		final ClusterSketchCodec.Loaded overlay=ClusterSketchCodec.read(root.resolve("split.20.pddl").toString());
		check(overlay.families.length==2 && overlay.families[0].members==6 && overlay.families[1].members==6,
			"One raw CLI split exports only the two actual six-member children");
		check(Arrays.equals(Files.readAllBytes(root.resolve("split.20.pddl")),Files.readAllBytes(root.resolve("split.21.pddl"))),
			"Repeated split rows reuse the same fresh ID and unchanged original parent");
		final ClusterSketchCodec.Loaded applied=ClusterSketchCodec.applyOverlay(ClusterSketchCodec.read(base.toString()),overlay,new int[]{0});
		check(applied.families.length==2 && applied.families[1].familyId==1,"Producer output directly satisfies the small-overlay consumer contract");
		final String[] fields=Files.readAllLines(root.resolve("split.results.tsv"),StandardCharsets.US_ASCII).get(1).split("\t",-1);
		check(fields.length==31 && fields[3].equals("1") && fields[23].equals("-12.0") && fields[30].equals("12"),
			"Raw split records its original count and consumes one budget slot");
		check(Double.parseDouble(fields[20])+Double.parseDouble(fields[21])+Double.parseDouble(fields[22])+Double.parseDouble(fields[23])<0,
			"A negative provisional gain must not censor the calibration proposal");
		check(Files.readAllLines(root.resolve("split.members.tsv"),StandardCharsets.US_ASCII).size()==25,
			"Each independent split retains exactly twelve original-member lineage rows");
		original.enableMemberRanking(new MemberRankInputs(input),5,null);
		final Path rankedSource=root.resolve("ranked_source.pddl");
		ClusterSketchCodec.write(rankedSource.toString(),false,new Family[]{original},2,64,7,false,0,0);
		for(String mode:new String[]{"off","on","default","parallel"}){
			final ArrayList<String> rankedArgs=new ArrayList<String>(Arrays.asList(probeArgs(root,"rank_"+mode,rankedSource,fasta,parent,
				membership,roster,requests,mode.equals("parallel") ? 4 : 1)));
			rankedArgs.addAll(Arrays.asList(controls));
			if(!mode.equals("off")){check(rankedArgs.remove("memberranking=f"),"Ranked raw probes must remove only their explicit legacy control");}
			if(mode.equals("on")){rankedArgs.add("memberranking=t");}
			new DdlClustStageB(rankedArgs.toArray(new String[rankedArgs.size()])).process();
			final ClusterSketchCodec.Loaded output=ClusterSketchCodec.read(root.resolve("rank_"+mode+".20.pddl").toString());
			check(output.families.length==2,"Default raw probes still publish exactly two actual children");
			for(ClusterSketchCodec.Entry child:output.families){
				check(mode.equals("off") ? child.memberRanking==null : child.memberRanking!=null && child.memberRanking.members==6 &&
					child.memberRanking.layer==6 && child.memberRanking.medianLength==80,
					"Independent raw children receive current members and one fixed hydrated layer, or no descriptors under OFF");
			}
			check(Arrays.equals(Files.readAllBytes(root.resolve("rank_"+mode+".20.pddl")),Files.readAllBytes(root.resolve("rank_"+mode+".21.pddl"))),
				"Repeated raw requests do not advance provenance or reuse published proposals");
		}
		for(String[] pair:new String[][]{{"split","rank_off"},{"rank_on","rank_default"},{"rank_on","rank_parallel"}}){
			for(String suffix:new String[]{".20.pddl",".21.pddl",".members.tsv"}){
				check(Arrays.equals(Files.readAllBytes(root.resolve(pair[0]+suffix)),Files.readAllBytes(root.resolve(pair[1]+suffix))),
					"Raw-probe codec/lineage parity on the same source payload: "+pair[0]+"/"+pair[1]);
			}
		}
		System.out.println("PROBE_MEMBER_RANK_DEFAULT_PASS children=2 members=12 layer=6");
	}
	private static String probeSequence(final int group){
		final String alphabet="ACDEFGHIKLMNPQRSTVWY";final java.util.Random random=new java.util.Random(2026103100L+group);
		final StringBuilder result=new StringBuilder(80);
		for(int i=0; i<30; i++){result.append(alphabet.charAt(random.nextInt(alphabet.length())));}
		for(int i=0; i<50; i++){result.append(alphabet.charAt(i%alphabet.length()));}
		return result.toString();
	}
	private static String[] probeArgs(final Path root,final String prefix,final Path source,final Path fasta,
			final Path parent,final Path members,final Path roster,final Path requests,final int threads){
		return new String[]{"mode=probe","memberranking=f","probes="+requests,"in="+fasta,"out="+root.resolve(prefix+".#.pddl"),
			"roster="+roster,"partitionmap="+parent,"stageacodecs="+source,"stageamembers="+members,
			"report="+root.resolve(prefix+".results.tsv"),"passesreport="+root.resolve(prefix+".work.tsv"),
			"members="+root.resolve(prefix+".members.tsv"),"minspan=5","t="+threads,"frozenindexthreads="+threads,
			"familybudget=2","temperature=0","ow=t"};
	}

	/** Uses the same sparse-ID source as legacy merging; NONE must retain each seed before any assignment. */
	private static void structuralModes(final Path root,final Path source,final Path fasta,final Path parent,
			final Path members,final Path roster) throws Exception{
		final ClusterSketchCodec.Entry[] originals=ClusterSketchCodec.read(source.toString()).families;
		for(String schedule:new String[]{"none","fission","fusion","layercake"}){
			for(int threads:new int[]{1,4}){
				final String prefix="structural_"+schedule+"_"+threads;
				new DdlClustStageB(structuralArgs(root,prefix,source,fasta,parent,members,roster,schedule,threads)).process();
				final ClusterSketchCodec.Loaded loaded=ClusterSketchCodec.read(root.resolve(prefix+".pddl").toString());
				final boolean merges=schedule.equals("fusion") || schedule.equals("layercake");
				check(loaded.families.length==(merges ? 4 : 5),"Opt-in schedule controls which structural operations execute");
				int count=0;
				for(int i=0; i<loaded.families.length; i++){
					final ClusterSketchCodec.Entry family=loaded.families[i];count+=family.members;
					if(!merges){
						check(family.familyId==i && family.name.equals("s0_f"+originals[i].familyId) &&
							family.members==originals[i].members && Arrays.equals(family.consensus,originals[i].consensus) &&
							Arrays.equals(family.sketch,originals[i].sketch) && Arrays.equals(family.monomers,originals[i].monomers),
							"No-operation export retains each original payload and assigns dense global identities");
					}
				}
				check(count==5,"Every structural mode conserves all five original tails");
				final List<String> owners=Files.readAllLines(root.resolve(prefix+".members.tsv"),StandardCharsets.US_ASCII);
				check(owners.size()==6,"Structural output includes every original ordinal once");
				for(int i=0; i<5; i++){
					check(owners.get(i+1).equals(i+"\toriginal"+i+"\t"+(merges && i==1 ? 0 : i)),
						"A single frozen fusion changes only its two endpoint memberships");
				}
				final List<String> phases=Files.readAllLines(root.resolve(prefix+".passes.tsv"),StandardCharsets.US_ASCII);
				check(phases.get(phases.size()-1).endsWith("\t"+(merges ? "CYCLE_LIMIT" : schedule.equals("none") ? "DISABLED" : "UNCHANGED")),
					"CLI reports explicit study bound separately from convergence and disabled work");
			}
			for(String suffix:new String[]{".pddl",".members.tsv",".decisions.tsv"}){
				check(Arrays.equals(Files.readAllBytes(root.resolve("structural_"+schedule+"_1"+suffix)),
					Files.readAllBytes(root.resolve("structural_"+schedule+"_4"+suffix))),"Structural CLI exact worker parity: "+schedule+suffix);
			}
		}
		for(String bad:new String[]{"cooling=1","temperature=-1","mintemperature=2","familybudget=0",
			"structuraltarget=0","maxcycles=-1","largesize=0","gainmargin=0","cohesionscale=-1",
			"separationscale=-1","sizescale=-1","fissionminsize=0","minchildsize=0","memberneighbors=0",
			"memberminhits=0","evidenceneighbors=0","fusionsizeratio=0.5","frozenindexthreads=5","mergemaxsize=20"}){
			final String[] base=structuralArgs(root,"invalid_structural",source,fasta,parent,members,roster,"fusion",4);
			final String[] invalid=Arrays.copyOf(base,base.length+1);invalid[base.length]=bad;
			boolean failed=false;
			try{new DdlClustStageB(invalid);}catch(IllegalArgumentException expected){failed=true;}
			check(failed && !Files.exists(root.resolve("invalid_structural.pddl")),"Invalid structural control fails before output: "+bad);
		}
		expectFail("Structural work needs explicit schedule","requires explicit schedule",
			new String[]{"mode=layercake"});
		expectFail("Structural controls cannot silently affect old modes","Structural schedule controls require mode=layercake",
			new String[]{"mode=mutual","schedule=fusion"});
		final byte[] originalBytes=Files.readAllBytes(source);
		final String[] aliasBase=structuralArgs(root,"alias_structural",source,fasta,parent,members,roster,"none",1);
		final String[] alias=Arrays.copyOf(aliasBase,aliasBase.length+2);
		alias[aliasBase.length]="out="+source;alias[aliasBase.length+1]="ow=t";
		expectFail("Structural output cannot overwrite source","Stage-B input/output paths alias",alias);
		check(Arrays.equals(originalBytes,Files.readAllBytes(source)),"Structural alias rejection preserves source bytes");
		structuralRanking(root,source,fasta,parent,members,roster);
	}

	/** All schedule modes hydrate original ordinals; only real publication advances their descriptor layer. */
	private static void structuralRanking(final Path root,final Path source,final Path fasta,final Path parent,
			final Path members,final Path roster)throws Exception{
		final ClusterSketchCodec.Loaded loaded=ClusterSketchCodec.read(source.toString());
		final ClusteringInput input=ClusteringInput.load(fasta.toString(),-1,false);
		final MemberRankInputs cache=new MemberRankInputs(input.sequences);
		final Family[] models=new Family[loaded.families.length];
		for(int i=0;i<models.length;i++){
			final ClusterSketchCodec.Entry entry=loaded.families[i];
			check(entry.members==1,"This literal structural fixture maps one original to each source model");
			final Family model=new Family(entry.familyId,entry.name,i,entry.consensus);
			model.setSketch(entry.sketch);model.setMonomers(entry.monomers);model.setScalars(entry.scalars,entry.scalarMask);
			model.enableMemberRanking(cache,5+i,null);models[i]=model;
		}
		final Path ranked=root.resolve("structural_ranked_source.pddl");
		ClusterSketchCodec.write(ranked.toString(),false,models,loaded.header);
		final byte[] before=Files.readAllBytes(ranked);
		for(String schedule:new String[]{"none","fission","fusion","layercake"}){
			for(String mode:new String[]{"off","on","default","parallel"}){
				final String prefix="rank_"+schedule+"_"+mode;
				final ArrayList<String> args=new ArrayList<String>(Arrays.asList(structuralArgs(root,prefix,ranked,fasta,parent,members,roster,
					schedule,mode.equals("parallel") ? 4 : 1)));
				if(!mode.equals("off")){check(args.remove("memberranking=f"),"Ranked schedule must remove its explicit OFF control");}
				if(mode.equals("on")){args.add("memberranking=t");}
				new DdlClustStageB(args.toArray(new String[args.size()])).process();
				final ClusterSketchCodec.Entry[] output=ClusterSketchCodec.read(root.resolve(prefix+".pddl").toString()).families;
				final boolean fusion=schedule.equals("fusion") || schedule.equals("layercake");
				check(output.length==(fusion ? 4 : 5),"Default schedules preserve the expected bounded family population");
				for(ClusterSketchCodec.Entry model:output){
					check(mode.equals("off") ? model.memberRanking==null : model.memberRanking!=null &&
						model.memberRanking.members==model.members && model.memberRanking.layer==(model.members==2 ? 11 : 10) &&
						model.memberRanking.medianLength==80,"Structural provenance follows imported layers and actual union membership");
				}
			}
			for(String[] pair:new String[][]{{"structural_"+schedule+"_1","rank_"+schedule+"_off"},
				{"rank_"+schedule+"_on","rank_"+schedule+"_default"},{"rank_"+schedule+"_on","rank_"+schedule+"_parallel"}}){
				for(String suffix:new String[]{".pddl",".decisions.tsv",".members.tsv"}){
					check(Arrays.equals(Files.readAllBytes(root.resolve(pair[0]+suffix)),Files.readAllBytes(root.resolve(pair[1]+suffix))),
						"Structural OFF/default/worker parity preserves complete codec, decisions and lineage");
				}
			}
		}
		check(Arrays.equals(before,Files.readAllBytes(ranked)),"Structural schedules leave source descriptors unchanged");
		// The layer setter selects the NEXT publication; build an actually measured MAX snapshot.
		final Family first=models[0];
		final Family lastLayer=new Family(first.familyId,first.id,first.memberIndex(0),first.consensus());
		lastLayer.setSketch(first.sketch());lastLayer.setMonomers(first.monomers());
		lastLayer.setScalars(first.scalars(),first.scalarMask());
		lastLayer.enableMemberRanking(cache,Integer.MAX_VALUE,null);models[0]=lastLayer;
		final Path overflow=root.resolve("structural_overflow.pddl");
		ClusterSketchCodec.write(overflow.toString(),false,models,loaded.header);
		check(ClusterSketchCodec.read(overflow.toString()).families[0].memberRanking.layer==Integer.MAX_VALUE,
			"Overflow input must serialize a measured MAX layer before testing hydration rejection");
		expectFail("Default structural hydration rejects imported overflow","Mutual hydration cannot advance beyond",
			defaultRanking(structuralArgs(root,"rank_overflow",overflow,fasta,parent,members,roster,"none",1)));
		check(!Files.exists(root.resolve("rank_overflow.pddl")),"Overflow cannot publish structural models");
		System.out.println("STRUCTURAL_MEMBER_RANK_DEFAULT_PASS initial_layer=10 union_layer=11 modes=4");
	}

	private static String[] structuralArgs(final Path root,final String prefix,final Path source,final Path fasta,
			final Path parent,final Path members,final Path roster,final String schedule,final int threads){
		return new String[]{"mode=layercake","memberranking=f","schedule="+schedule,"in="+fasta,"out="+root.resolve(prefix+".pddl"),
			"roster="+roster,"partitionmap="+parent,"stageacodecs="+source,"stageamembers="+members,
			"report="+root.resolve(prefix+".decisions.tsv"),"passesreport="+root.resolve(prefix+".passes.tsv"),
			"members="+root.resolve(prefix+".members.tsv"),"minspan=5","t="+threads,"frozenindexthreads="+threads,
			"familybudget=2","structuraltarget=.99","gainmargin=.1","maxcycles=1","temperature=0","ow=f"};
	}

	/** Hold query/rebuild threads at four and vary only frozen-index construction. */
	private static void frozenIndexMutual(Path root,Path source,Path fasta,Path parent,Path members,Path roster) throws Exception{
		final Path serialConsensus=root.resolve("frozen_omitted.faa"), serialRoster=root.resolve("frozen_omitted.roster.tsv");
		new DdlClustStageB(new String[]{"mode=collect","codecs="+root.resolve("merged4.pddl"),
			"out="+serialConsensus,"roster="+serialRoster,"t=4","ow=f"}).process();
		for(int builders:new int[]{1,2,4}){
			final String prefix="frozen"+builders;
			new DdlClustStageB(new String[]{"mode=mutual","memberranking=f","in="+fasta,"out="+root.resolve(prefix+".pddl"),
				"roster="+roster,"partitionmap="+parent,"stageacodecs="+source,"stageamembers="+members,
				"report="+root.resolve(prefix+".pairs.tsv"),"passesreport="+root.resolve(prefix+".passes.tsv"),
				"members="+root.resolve(prefix+".members.tsv"),"minspan=5","mergemaxsize=20",
				"shortlist=50","t=4","frozenindexthreads="+builders,"ow=f"}).process();
			for(String[] names:new String[][]{{"merged4.pddl",".pddl"},{"pairs4.tsv",".pairs.tsv"},{"members4.tsv",".members.tsv"}}){
				check(Arrays.equals(Files.readAllBytes(root.resolve(names[0])),Files.readAllBytes(root.resolve(prefix+names[1]))),
					"Frozen-index mutual exact byte parity, builders="+builders+": "+names[0]);
			}
			final List<String> a=Files.readAllLines(root.resolve("passes4.tsv"),StandardCharsets.US_ASCII),
				b=Files.readAllLines(root.resolve(prefix+".passes.tsv"),StandardCharsets.US_ASCII);
			check(a.size()==b.size() && a.get(0).equals(b.get(0)) && a.get(0).split("\t",-1).length==11 &&
				a.get(0).endsWith("\telapsed_ns"),"Frozen-index pass schema/count, only elapsed_ns is excluded");
			for(int row=1;row<a.size();row++){
				final String[] x=a.get(row).split("\t",-1),y=b.get(row).split("\t",-1);
				check(x.length==11 && y.length==11,"Frozen-index pass widths");
				for(int col=0;col<10;col++){check(x[col].equals(y[col]),"Frozen-index pass scientific column "+col);}
			}
			new DdlClustStageB(new String[]{"mode=collect","codecs="+root.resolve(prefix+".pddl"),
				"out="+root.resolve(prefix+".faa"),"roster="+root.resolve(prefix+".roster.tsv"),
				"t=4","frozenindexthreads=1","ow=f"}).process();
			check(Arrays.equals(Files.readAllBytes(serialConsensus),Files.readAllBytes(root.resolve(prefix+".faa"))) &&
				Arrays.equals(Files.readAllBytes(serialRoster),Files.readAllBytes(root.resolve(prefix+".roster.tsv"))),
				"Consensus/roster parity and explicit serial option remains valid in collect mode");
		}
	}

	private static void frozenIndexInvalid(Path root) throws Exception{
		final String expected="frozenindexthreads must lie in [1,t]; parallel index construction requires mode=mutual";
		for(String mode:new String[]{"mutual","collect","union"}){
			final int[] invalid=mode.equals("mutual") ? new int[]{0,-1,5} : new int[]{0,-1,5,2};
			for(int count:invalid){
				final Path out=root.resolve("bad_frozen_"+mode+"_"+count+".pddl");
				expectFail("Frozen-index bounds/mode "+mode+"/"+count,expected,
					new String[]{"mode="+mode,"out="+out,"t=4","frozenindexthreads="+count,"ow=f"});
				check(!Files.exists(out),"Rejected frozen-index configuration publishes no output");
			}
		}
	}

	/** Two independent reciprocal pairs exercise simultaneous rebuild workers, not only query parallelism. */
	private static void mutualParallelPairs(Path root) throws Exception{
		Files.createDirectory(root);
		final ProteinSequence[] input={new ProteinSequence("a", repeat("A", 80)), new ProteinSequence("w", repeat("W", 80)),
			new ProteinSequence("a2", repeat("A", 80)), new ProteinSequence("w2", repeat("W", 80))};
		final FamilySketchBuilder builder=new FamilySketchBuilder(64, 2, 7, false, 0, 0, null);
		final Family[] families=new Family[input.length];
		final StringBuilder map=new StringBuilder(ClusteringShardWriter.MAPPING_HEADER).append('\n');
		final StringBuilder members=new StringBuilder("input_index\tinput_id\tfamily_id\tstatus\n");
		for(int i=0; i<input.length; i++){
			families[i]=family(i, "f"+i, input, builder, i);
			map.append(i).append("\t0\t").append(i).append('\t').append(input[i].id).append('\n');
			members.append(i).append('\t').append(input[i].id).append('\t').append(i).append("\tassigned\n");
		}
		final Path source=root.resolve("source.pddl"), fasta=root.resolve("input.faa"), roster=root.resolve("roster.tsv");
		ClusterSketchCodec.write(source.toString(), false, families, 2, 64, 7, false, 0, 0);
		fasta(fasta, input); write(root.resolve("map.tsv"), map.toString()); write(root.resolve("source.members.tsv"), members.toString());
		new DdlClustStageB(new String[]{"mode=collect", "codecs="+source, "out="+root.resolve("collected.faa"), "roster="+roster, "ow=f"}).process();
		for(int threads : new int[]{1,4}){
			new DdlClustStageB(new String[]{"mode=mutual", "memberranking=f", "in="+fasta, "out="+root.resolve("merged"+threads+".pddl"),
				"roster="+roster, "partitionmap="+root.resolve("map.tsv"), "stageacodecs="+source,
				"stageamembers="+root.resolve("source.members.tsv"), "report="+root.resolve("pairs"+threads+".tsv"),
				"passesreport="+root.resolve("passes"+threads+".tsv"), "members="+root.resolve("members"+threads+".tsv"),
				"minspan=5", "t="+threads, "ow=f"}).process();
		}
		for(String pattern : new String[]{"merged#.pddl", "pairs#.tsv", "members#.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve(pattern.replace("#", "1"))),
				Files.readAllBytes(root.resolve(pattern.replace("#", "4")))), "Concurrent pair rebuild parity: "+pattern);
		}
		final ClusterSketchCodec.Loaded result=ClusterSketchCodec.read(root.resolve("merged4.pddl").toString());
		check(result.families.length==2 && result.families[0].members==2 && result.families[1].members==2,
			"Two distinct identical-sequence pairs yield two separate two-member unions");
		final List<String> pairs=Files.readAllLines(root.resolve("pairs4.tsv"), StandardCharsets.US_ASCII);
		check(pairs.size()==3 && pairs.get(1).startsWith("1\t0\t2\t1\t1\t0\t2\t") &&
			pairs.get(2).startsWith("1\t1\t3\t1\t1\t1\t2\t"), "Both independent unions occur in the same frozen pass");
	}

	private static Family family(int id, String name, ProteinSequence[] originals,
			FamilySketchBuilder builder, int member){
		final Family family=new Family(id, name, member, originals[member].enc);
		family.refreshStats(originals, builder);
		return family;
	}

	private static Family seedFamily(int id, String name, String consensus, int members,
			FamilySketchBuilder builder){
		final ProteinSequence[] seeds=new ProteinSequence[members];
		for(int i=0; i<members; i++){seeds[i]=new ProteinSequence("seed"+i, consensus);}
		final Family family=new Family(id, name, 0, seeds[0].enc);
		for(int i=1; i<members; i++){family.addMember(i);}
		family.refreshStats(seeds, builder);
		return family;
	}

	private static CaseResult singleUnionCase(Path root, String seed, String[] members, double minId) throws Exception{
		Files.createDirectory(root);
		final ProteinSequence[] originals=new ProteinSequence[members.length];
		final Family[] stageFamilies=new Family[members.length];
		final FamilySketchBuilder builder=new FamilySketchBuilder(64, 2, 7, false, 0, 0, null);
		final Path stage0=root.resolve("stage0.pddl"), stageB=root.resolve("stageb.pddl"),
			roster=root.resolve("roster.tsv"), originalFasta=root.resolve("originals.faa"),
			map=root.resolve("parent.tsv"), stage0Members=root.resolve("stage0.members.tsv"),
			stageBMembers=root.resolve("stageb.members.tsv"), out=root.resolve("union.pddl"),
			report=root.resolve("union.tsv"), collectFasta=root.resolve("stageb_input.faa");
		final StringBuilder mapText=new StringBuilder(ClusteringShardWriter.MAPPING_HEADER).append('\n');
		final StringBuilder stageAText=new StringBuilder("input_index\tinput_id\tfamily_id\tstatus\n");
		final StringBuilder stageBText=new StringBuilder("input_index\tinput_id\tfamily_id\tstatus\n");
		for(int i=0; i<members.length; i++){
			originals[i]=new ProteinSequence("p"+i, members[i]);
			stageFamilies[i]=family(10+i, "s0_f"+(10+i), originals, builder, i);
			mapText.append(i).append("\t0\t").append(i).append("\tp").append(i).append('\n');
			stageAText.append(i).append("\tp").append(i).append('\t').append(10+i).append("\tassigned\n");
			stageBText.append(i).append("\ts0_f").append(10+i).append("\t30\tassigned\n");
		}
		ClusterSketchCodec.write(stage0.toString(), false, stageFamilies, 2, 64, 7, false, 0, 0);
		ClusterSketchCodec.write(stageB.toString(), false,
			new Family[]{seedFamily(30, "stageb_30", seed, members.length, builder)}, 2, 64, 7, false, 0, 0);
		fasta(originalFasta, originals);
		write(map, mapText.toString()); write(stage0Members, stageAText.toString());
		write(stageBMembers, stageBText.toString());
		new DdlClustStageB(new String[]{"mode=collect", "codecs="+stage0,
			"out="+collectFasta, "roster="+roster, "ow=f"}).process();
		new DdlClustStageB(new String[]{"mode=union", "in="+originalFasta, "out="+out,
			"roster="+roster, "partitionmap="+map, "stageacodecs="+stage0,
			"stageamembers="+stage0Members, "stagebcodec="+stageB,
			"stagebmembers="+stageBMembers, "report="+report, "minspan=5",
			"minid="+minId, "t=4", "ow=f"}).process();
		final ClusterSketchCodec.Entry entry=ClusterSketchCodec.read(out.toString()).families[0];
		final String reportLine=Files.readAllLines(report, StandardCharsets.US_ASCII).get(1);
		check(entry.members==members.length, "Union keeps all source members for "+root.getFileName());
		return new CaseResult(reportLine, decode(entry.consensus));
	}

	private static String[] unionArgs(Path originalFasta, Path out, Path roster, Path map,
			Path stage0, Path stage1, Path stage0Members, Path stage1Members, Path stageB,
			Path stageBMembers, Path report){
		return new String[]{"mode=union", "in="+originalFasta, "out="+out, "roster="+roster,
			"partitionmap="+map, "stageacodecs="+stage0+","+stage1,
			"stageamembers="+stage0Members+","+stage1Members, "stagebcodec="+stageB,
			"stagebmembers="+stageBMembers, "report="+report, "minspan=5", "ow=t"};
	}

	private static void expectFail(String why, String expected, String[] args) throws Exception{
		checks++;
		try{new DdlClustStageB(args).process();}
		catch(RuntimeException failure){
			checkFailure(failure, expected);
			return;
		}
		throw new AssertionError("Expected failure: "+why);
	}

	private static void checkFailure(Throwable failure,String expected){
		for(Throwable current=failure; current!=null; current=current.getCause()){
			final String message=current.getMessage();
			if(message!=null && message.contains(expected)){return;}
		}
		throw new AssertionError("Wrong failure, expected '"+expected+"'", failure);
	}

	private static void fasta(Path path, ProteinSequence[] input) throws Exception{
		final StringBuilder sb=new StringBuilder();
		for(ProteinSequence seq : input){
			sb.append('>').append(seq.id).append('\n').append(decode(seq.enc)).append('\n');
		}
		write(path, sb.toString());
	}

	private static String decode(byte[] enc){
		return new String(ClusteringInput.decode(enc), StandardCharsets.US_ASCII);
	}

	private static String repeat(String s,int count){
		final char[] chars=new char[s.length()*count];
		for(int i=0; i<count; i++){s.getChars(0, s.length(), chars, i*s.length());}
		return new String(chars);
	}

	private static void write(Path path,String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}

	private static void check(boolean ok,String why){checks++; if(!ok){throw new AssertionError(why);}}
	private static int checks;

	private static final class CaseResult{
		CaseResult(String report_, String consensus_){
			report=report_; consensus=consensus_;
		}
		final String report, consensus;
	}

	private DdlClustStageBTest(){}
}
