package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import ddl.DDLFamilyHoldout;
import fileIO.ByteFile;
import prot.FamilyCalibrate;
import prot.FamilyCalibration.VerificationMode;
import prot.FamilyPanelRunner.Panel;
import shared.Shared;

/** Bounded resident preparation/calibration/held-out integration with synthetic inherited identities.
 * No fixture hashes are computed, including in the default-mode constructor rejection cases.
 * Legacy geometry isolates the saved negative gate; a synthetic complete dimer map also
 * exercises metric-aware preparation. Neither is an accuracy or serialized-reload parity claim.
 * @author Qiqi
 */
public final class FamilyResidentCalibrationTest{

	/** Exercise mandatory12K6 preparation, held-out scoring and serialized structural reload. */
	public static void official12(Path root)throws Exception{
		Files.createDirectory(root);final Fixture fixture=new Fixture(root.resolve("inputs"));
		write(fixture.manifest,text(fixture.manifest).replace("k=5\n","k=6\n").replace("width=1024\n","width=512\n")
			.replace("ddlreject=1024\n","ddlreject=512\n")+"aabits="+ddl.ProteinBlacklistPolicy.ALPHABET+"\n");
		final Path out=root.resolve("complete");run(fixture.manifest,out,fixture.positive,3);
		check(Files.isRegularFile(out.resolve("PASS")) && Files.isRegularFile(out.resolve("heldout/PASS")),
			"Official12K6 resident preparation and held-out scoring complete");
		final Path config=out.resolve("classifier.config");
		final String recipe=text(config);
		for(String field:new String[]{"blacklistsha80="+ddl.ProteinBlacklistPolicy.SHA80,"blacklistk=6",
			"blacklistmode=reduced_literal","blacklistwords=4096","blacklistprovenance="+ddl.ProteinBlacklistPolicy.PROVENANCE}){
			check(recipe.contains(field+"\n"),"Effective official policy serialized: "+field);
		}
		final prot.FamilyScoringPackage pack=new prot.FamilyScoringPackage(config,VerificationMode.STRUCTURAL);
		final prot.FamilyTailCascade.LoadedClassifier loaded=pack.load(1),parallel=pack.load(3);
		check(!loaded.blacklist().identityMeasured() && loaded.blacklist()==ddl.ProteinBlacklistPolicy.resolve(ddl.ProteinBlacklistPolicy.GEOMETRY,loaded.blacklist()) &&
			parallel.blacklist()==loaded.blacklist(),"Reloaded and resident filters share official immutable policy");
		for(char aa:new char[]{'A','F','C'}){
			final byte[] query=sequence(aa).getBytes(StandardCharsets.US_ASCII);
			check(!loaded.newWorker().assign("q",query).assigned() && !parallel.newWorker().assign("q",query).assigned(),
				"Frozen high-reject fixture gates agree after structural reload across thread counts");
		}
		pack.verify();
		final Path wrongPolicy=root.resolve("wrong_policy.config"),missingPolicy=root.resolve("missing_policy.config");
		write(wrongPolicy,recipe.replace("blacklistwords=4096\n","blacklistwords=8192\n"));
		write(missingPolicy,recipe.replace("blacklistsha80="+ddl.ProteinBlacklistPolicy.SHA80+"\n",""));
		for(Path invalid:new Path[]{wrongPolicy,missingPolicy}){
			reject(()->new prot.FamilyScoringPackage(invalid,VerificationMode.STRUCTURAL),"official reduced-shadow4096");
		}
		System.out.println("OFFICIAL12K6_RESIDENT_RELOAD_PASS");
	}

	public static void main(String[] args) throws Exception{
		if((args.length==3 || args.length==4) && args[0].equals("custom")){
			if(args.length==4){dna.Data.DORI_RESOURCES=Paths.get(args[3]).toAbsolutePath().toString()+java.io.File.separator;}
			verifyCustom(Paths.get(args[1]),text(Paths.get(args[2])).trim(),args.length==4);
			System.out.println("CUSTOM_STRUCTURAL_PREPARATION_PASS checks="+checks);return;
		}
		check((args.length==1 || args.length==2) && args[0].startsWith("fixture="),"Require fixture=NEW_ABSOLUTE_DIR [SELECTED_RESOURCE_DIR]");
		final Path root=Paths.get(args[0].substring(8)).normalize();
		check(root.isAbsolute() && !Files.exists(root),"Fixture directory must be absolute and fresh");
		Files.createDirectory(root);Shared.AMINO_IN=true;
		final Fixture f=new Fixture(root.resolve("inputs"));
		final Path success=root.resolve("complete");
		run(f.manifest,success,f.positive,3);
		check(Files.isRegularFile(success.resolve("PASS")) && Files.isRegularFile(success.resolve("FROZEN")),
			"Successful resident calibration must complete and freeze");
		check(!Files.exists(success.resolve("initial")),"Resident preparation must not create a provisional reload package");
		check(value(success.resolve("calibration.tsv"),"status").equals("converged"),"All negatives must converge");
		verifySummary(success.resolve("summary.tsv"),NAMES,new long[]{2,2,2,2},"negative",false);
		verifySummary(success.resolve("heldout/summary.tsv"),new String[]{"positive"},new long[]{3},"positive",false);
		verifyPopulation(success);
		verifyProvenance(success,f);
		check(value(success.resolve("heldout_gate.tsv"),"allowed").equals("true") &&
			value(success.resolve("heldout_status.tsv"),"status").equals("completed"),"Positive scoring requires the saved passing gate");
		check(value(success.resolve("heldout/scoring_time.tsv"),"queries").equals("3"),"Timing denominator must be the complete held-out population");
		check(Files.isRegularFile(success.resolve("heldout/PASS")),"Successful held-out scoring must publish its own PASS");

		// Both EOF undercount and extra records must fail, without a final evaluation success marker.
		for(int expected:new int[]{2,4}){
			final Path out=root.resolve("heldout_count_"+expected);
			reject(()->run(f.manifest,out,f.positive,expected),"count differs from its declared population");
			check(value(out.resolve("heldout_gate.tsv"),"allowed").equals("true"),"Count failure must occur after the negative gate passed");
			check(Files.isRegularFile(out.resolve("FROZEN")),"A failed positive evaluation does not undo completed negative calibration");
			check(!Files.exists(out.resolve("PASS")) && !Files.exists(out.resolve("heldout/PASS")),"Count mismatch must not publish evaluation completion");
		}

		// Distinct fixtures avoid modifying dependencies held by any previous preparation session.
		for(String mutation:new String[]{"negative_short","negative_extra","duplicate_id","member_total"}){
			final Fixture bad=new Fixture(root.resolve(mutation));
			final String reason;
			if(mutation.equals("negative_short")){
				bad.writePanel(0,1,false);reason="Wrong negative query count";
			}else if(mutation.equals("negative_extra")){
				bad.writePanel(0,3,false);reason="Too many negative controls";
			}else if(mutation.equals("duplicate_id")){
				bad.writePanel(0,2,true);reason="Duplicate or empty negative query identity";
			}else{
				bad.writeManifest(bad.manifest,false,5);reason="Family count/assigned conservation failed";
			}
			final Path out=root.resolve(mutation+"_result");
			reject(()->run(bad.manifest,out,bad.positive,3),reason);
			check(!Files.exists(out.resolve("PASS")) && !Files.exists(out.resolve("FROZEN")) &&
				!Files.exists(out.resolve("heldout_gate.tsv")),"Incomplete negative/member populations cannot reach the positive gate");
		}

		// A regular-file parent makes this positive path unusable even for privileged test users.
		final Path blocker=root.resolve("not_a_directory");write(blocker,"positive-open sentinel\n");
		final Path unreadable=blocker.resolve("positive.faa");
		reject(()->new Panel("positive",unreadable.toString(),INHERITED,true,VerificationMode.STRUCTURAL,3),"Missing sequence panel");
		final Path gateManifest=f.root.resolve("gate_failure.manifest"),skipped=root.resolve("negative_gate_failure");
		f.writeManifest(gateManifest,true,4);
		run(gateManifest,skipped,unreadable,3);
		check(value(skipped.resolve("calibration.tsv"),"status").equals("ddl_floor_exceeds_target"),
			"Exact member negatives must expose the unmodifiable DDL shortcut floor");
		check(value(skipped.resolve("heldout_gate.tsv"),"allowed").equals("false") &&
			value(skipped.resolve("heldout_gate.tsv"),"reason").equals("calibration_not_converged") &&
			value(skipped.resolve("heldout_status.tsv"),"status").equals("skipped"),"Failed negatives must persist an explicit skip");
		check(Files.isRegularFile(skipped.resolve("PASS")) && !Files.exists(skipped.resolve("FROZEN")) &&
			!Files.exists(skipped.resolve("heldout")),"Skipped positives must leave no evaluated positive result or frozen classifier");
		check(text(blocker).equals("positive-open sentinel\n"),"The positive sentinel must remain unchanged");
		verifySummary(skipped.resolve("summary.tsv"),NAMES,new long[]{2,2,2,2},"negative",true);
		verifyPooledGate(root,unreadable);

		// Constructor validation precedes preparation and therefore precedes default hash verification.
		// This checks the unchanged opt-in boundary, not the default digest algorithm itself.
		for(String mode:new String[]{null,"sha80"}){
			final Path out=root.resolve(mode==null ? "default_rejected" : "sha80_rejected");
			final String[] argv=arguments(f.manifest,out,f.positive,3,mode);
			reject(()->new FamilyCalibrate(argv),"Resident held-out scoring requires structural preparation");
			check(!Files.exists(out),"Default verification must not silently opt into resident held-out scoring");
		}
		final Path backed=root.resolve("file_backed_rejected");
		reject(()->new FamilyCalibrate(new String[]{"in="+success.resolve("classifier.config"),
			"neg="+f.panels[0],"out="+backed,"verification=structural"}),
			"Structural existing-package calibration requires inherited panel identities and exact panel size");
		check(!Files.exists(backed),"Structural package reuse requires explicit panel identities and populations");
		verifyMetricResident(root);
		verifyUnmeasuredResident(root);
		verifyFilteredResident(root);
		verifyBlacklistRejections(root);
		if(args.length==2){
			dna.Data.DORI_RESOURCES=Paths.get(args[1]).toAbsolutePath().normalize().toString()+"/";
			verifySelectedResident(root);
		}
		System.out.println("FAMILY_RESIDENT_CALIBRATION_TEST_PASS checks="+checks);
	}

	/** Real experimental maps retain full-word filters through calibration and disk reload. */
	private static void verifyCustom(Path root,String dimerGroups,boolean controls)throws Exception{
		check(root.isAbsolute() && !Files.exists(root),"Fresh custom calibration root");Files.createDirectory(root);Shared.AMINO_IN=true;
		long measuredCalls=0;
		final String[] alphabets={"all:A/C/DE/FY/G/IV/LM/NH/P/RQK/ST/W",
			"all:A/C/D/E/FY/G/H/IV/LM/N/P/Q/RK/S/T/W","dimer:1,3,5:"+dimerGroups,
			ddl.ProteinBlacklistPolicy.ALPHABET,"all:A/C/D/E/FWY/G/H/IV/K/LM/N/P/Q/R/S/T",ddl.ProteinBlacklistPolicy.DIMER120_ALPHABET};
		for(int geometry=0;geometry<(controls ? 6 : 3);geometry++){
			for(String view:new String[]{"long","short"}){
				final boolean shortView=view.equals("short");final int k=geometry%3==1 ? 5 : 6;
				final Fixture f=new Fixture(root.resolve("custom"+geometry+"_"+view+"_inputs"));f.writeMetricManifest(syntheticDimerMap());
				final ddl.ProteinKmerGeometry map=ddl.ProteinKmerGeometry.make(k,alphabets[geometry]);
				final ddl.ProteinBlacklistPolicy.Descriptor policy=ddl.ProteinBlacklistPolicy.policy(map);
				check((policy==null)==(geometry<3),"Only explicit historical controls select official defaults");
				final String mode=geometry%3==1 ? "literal" : "reduced_literal";
				final Path filter=policy==null ? f.root.resolve("filter.faa") : policy.resource();
				if(policy==null){write(filter,">excluded\n"+(k==5 ? "CCCCC" : "CCCCCC")+"\n");}
				final String identity=policy==null ? "unmeasured" : policy.sha80;
				final int words=policy==null ? 1 : policy.words;
				final String provenance=policy==null ? "custom_cutoff_fixture" : policy.provenance;
				String recipe=text(f.manifest).replace("aabits=dimer:1,2,3,4:"+syntheticDimerMap(),"aabits="+alphabets[geometry]);
				if(k==6){recipe=recipe.replace("k=5\n","k=6\n");}
				for(String key:new String[]{"families","chunks","consensus","lengthmedians","familyhydro"}){
					recipe=recipe.replace(key+"sha80="+INHERITED+"\n",key+"sha80=unmeasured\n");
				}
				String chunks=text(f.root.resolve("chunks.tsv")).replace(INHERITED,"unmeasured");
				if(shortView){
					for(String kind:new String[]{"member","negative","query"}){
						recipe=recipe.replace("min"+kind+"aa=100\n","min"+kind+"aa=40\n").replace("max"+kind+"aa=100\n","max"+kind+"aa=99\n");
					}
					recipe=recipe.replace("alignmentminspan=70\n","alignmentminspan=28\n");
					for(String name:new String[]{"consensus.faa","members.tsv","positive.faa","k5.faa","shadow.faa","reversed.faa","shuffled.faa"}){
						String content=text(f.root.resolve(name));for(char aa:new char[]{'A','F','C'}){content=content.replace(sequence(aa),sequence(aa).substring(0,60));}
						write(f.root.resolve(name),content);
					}
					write(f.root.resolve("medians.tsv"),text(f.root.resolve("medians.tsv")).replace("\t100","\t60"));
					chunks=chunks.replace("\t4\t400\n","\t4\t240\n");
				}
				recipe+="blacklist="+filter+"\nblacklistsha80="+identity+"\nblacklistk="+k+"\nblacklistwords="+words+"\nblacklistmode="+mode+"\nblacklistprovenance="+provenance+"\n";
				write(f.manifest,recipe);write(f.root.resolve("chunks.tsv"),chunks);
				final Path out=root.resolve("custom"+geometry+"_"+view+"_complete");run(f.manifest,out,f.positive,3);
				check(Files.isRegularFile(out.resolve("FROZEN")) && Files.isRegularFile(out.resolve("heldout/PASS")),"Custom preparation and held-out scoring complete");
				final prot.FamilyScoringPackage pack=new prot.FamilyScoringPackage(out.resolve("classifier.config"),VerificationMode.STRUCTURAL);
				final prot.FamilyTailCascade.LoadedClassifier one=pack.load(1),three=pack.load(3);
				check(!one.blacklist().identityMeasured() && !three.blacklist().identityMeasured(),"Structural reload never promotes a declared identity to a measured digest");
				check(one.blacklist().sha80().equals(identity) && one.blacklist().wordLength()==k &&
					one.blacklist().matchMode().label().equals(mode),"Disk reload retains exact custom word geometry and unmeasured provenance");
				one.blacklist().validate(map);three.blacklist().validate(map);
				for(char aa:new char[]{'A','F','C'}){
					final byte[] query=sequence(aa).substring(0,shortView ? 60 : 100).getBytes(StandardCharsets.US_ASCII);
					check(!one.newWorker().assign("q",query).assigned() && !three.newWorker().assign("q",query).assigned(),"High-reject control survives file reload and thread counts");
				}
				pack.verify();
				measuredCalls+=verifyAlignmentWork(pack,f,shortView);
			}
		}
		check(measuredCalls>0,"Alignment-count fixtures must exercise real alignment calls");
	}

	/** Compare ordinary multithreaded panel totals to independent per-candidate trace events. */
	private static long verifyAlignmentWork(prot.FamilyScoringPackage pack,Fixture f,boolean shortView)throws Exception{
		final prot.FamilyTailCascade.LoadedClassifier loaded=pack.loadWithOverrides(1,50,1);
		final WorkTrace positive=new WorkTrace(),negative=new WorkTrace();
		final prot.FamilyTailCascade.Worker trace=loaded.newWorker();
		for(char aa:new char[]{'A','F','C'}){trace.assignTrace("p",sequence(aa).substring(0,shortView ? 60 : 100).getBytes(StandardCharsets.US_ASCII),positive);}
		for(char aa:new char[]{'A','F'}){trace.assignTrace("n",sequence(aa).substring(0,shortView ? 60 : 100).getBytes(StandardCharsets.US_ASCII),negative);}
		final Panel[] panels={new Panel("positive",f.positive.toString(),"unmeasured",true,VerificationMode.STRUCTURAL,3),
			new Panel("shadow",f.root.resolve("shadow.faa").toString(),"unmeasured",false,VerificationMode.STRUCTURAL,2)};
		final WorkTrace[] expected={positive,negative};
		for(int threads:new int[]{1,3}){
			final prot.FamilyPanelRunner runner=new prot.FamilyPanelRunner(loaded,panels,threads,1,System.err,"alignment-work",shortView ? 40 : 100,shortView ? 99 : 100,false);
			runner.process();
			for(int panel=0;panel<panels.length;panel++){
				final WorkTrace e=expected[panel];
				check(runner.alignmentCalls(panel)==e.calls && runner.alignmentUnavailable(panel)==e.unavailable && runner.alignmentFailures(panel)==e.failed,
					"Panel work counts match trace events independently of workers and panel boundaries");
			}
		}
		return positive.calls+negative.calls;
	}
	private static final class WorkTrace implements prot.FamilyTailCascade.CandidateObserver{
		@Override public void accept(int length,String terminal,int family,int ddl,float identity,float blosum){}
		@Override public void alignedCandidate(int family,int attempt,float identity,float blosum,double upper,double blosumUpper,float lower){calls++;}
		@Override public void unalignedCandidate(int family,int attempt,String reason,prot.FamilyAlignmentMetric.State state){
			if(reason.equals("alignment_unavailable")){unavailable++;}else{calls++;failed++;}
		}
		long calls,unavailable,failed;
	}

	/** Tiny selected geometry/view inputs also support fresh-process launcher validation. */
	private static void verifySelectedResident(Path root)throws Exception{
		final String[] alphabets={ddl.ProteinBlacklistPolicy.ALPHABET,"all:A/C/D/E/FWY/G/H/IV/K/LM/N/P/Q/R/S/T",
			ddl.ProteinBlacklistPolicy.DIMER120_ALPHABET};
		for(int geometry=0;geometry<alphabets.length;geometry++){
			for(String view:new String[]{"short","long"}){
				final boolean shortView=view.equals("short");
				final Fixture f=new Fixture(root.resolve("selected"+geometry+"_"+view+"_inputs"));
				f.writeMetricManifest(syntheticDimerMap());
				String recipe=text(f.manifest).replace("aabits=dimer:1,2,3,4:"+syntheticDimerMap(),"aabits="+alphabets[geometry]);
				if(geometry!=1){recipe=recipe.replace("k=5\n","k=6\n");}
				for(String key:new String[]{"families","chunks","consensus","lengthmedians","familyhydro"}){
					recipe=recipe.replace(key+"sha80="+INHERITED+"\n",key+"sha80=unmeasured\n");
				}
				String chunks=text(f.root.resolve("chunks.tsv")).replace(INHERITED,"unmeasured");
				if(shortView){
					for(String kind:new String[]{"member","negative","query"}){
						recipe=recipe.replace("min"+kind+"aa=100\n","min"+kind+"aa=40\n")
							.replace("max"+kind+"aa=100\n","max"+kind+"aa=99\n");
					}
					recipe=recipe.replace("alignmentminspan=70\n","alignmentminspan=28\n");
					for(String name:new String[]{"consensus.faa","members.tsv","positive.faa","k5.faa","shadow.faa","reversed.faa","shuffled.faa"}){
						String content=text(f.root.resolve(name));
						for(char aa:new char[]{'A','F','C'}){content=content.replace(sequence(aa),sequence(aa).substring(0,60));}
						write(f.root.resolve(name),content);
					}
					write(f.root.resolve("medians.tsv"),text(f.root.resolve("medians.tsv")).replace("\t100","\t60"));
					chunks=chunks.replace("\t4\t400\n","\t4\t240\n");
					// Negative lengths may be narrower than the query contract; positives include both endpoints.
					recipe=recipe.replace("minnegativeaa=40\n","minnegativeaa=60\n").replace("maxnegativeaa=99\n","maxnegativeaa=60\n");
					write(f.positive,">minimum\n"+sequence('A').substring(0,40)+"\n>middle\n"+sequence('F').substring(0,60)+
						"\n>maximum\n"+sequence('C').substring(0,99)+"\n");
				}
				write(f.manifest,recipe);write(f.root.resolve("chunks.tsv"),chunks);
				final Path out=root.resolve("selected"+geometry+"_"+view+"_complete");run(f.manifest,out,f.positive,3);
				check(Files.isRegularFile(out.resolve("PASS")) && Files.isRegularFile(out.resolve("FROZEN")),"Selected preparation must converge on its tiny high-reject control");
				check(Files.isRegularFile(out.resolve("heldout/PASS")) && value(out.resolve("heldout/scoring_time.tsv"),"queries").equals("3"),
					"Selected held-out scoring must complete every query in its own SHORT or LONG length range");
				verifySummary(out.resolve("heldout/summary.tsv"),new String[]{"positive"},new long[]{3},"positive",false);
				check(value(out.resolve("verification.tsv"),"input_identity").equals("inherited_with_unmeasured_inputs"),"Selected fresh metadata identity must remain explicit");
				final ddl.ProteinKmerGeometry selected=ddl.ProteinKmerGeometry.make(geometry==1 ? 5 : 6,alphabets[geometry]);
				final ddl.ProteinBlacklistPolicy.Descriptor policy=ddl.ProteinBlacklistPolicy.policy(selected);
				final String produced=text(out.resolve("classifier.config"));
				check(produced.contains("blacklistsha80="+policy.sha80+"\n") && produced.contains("aabits="+alphabets[geometry]+"\n"),"Selected geometry and filter must survive preparation");
				final prot.FamilyScoringPackage pack=new prot.FamilyScoringPackage(out.resolve("classifier.config"),VerificationMode.STRUCTURAL);
				final ddl.DDLPentamerBlacklist loaded=pack.load(1).blacklist();
				check(!loaded.identityMeasured() && loaded==ddl.ProteinBlacklistPolicy.resolve(selected,loaded),"Reload must retain the selected structural immutable filter");pack.verify();
				if(geometry==0){
					for(int length:shortView ? new int[]{39,100} : new int[]{99,101}){
						final Path invalid=f.root.resolve("heldout_length_"+length+".faa"),failed=root.resolve("selected0_"+view+"_length_"+length);
						final char[] aa=new char[length];Arrays.fill(aa,'A');write(invalid,">outside\n"+new String(aa)+"\n");
						final String reason="Heuristic query length is outside its released interval: "+length+" not in "+(shortView ? "40..99" : "100..100");
						boolean rejected=false;
						try{run(f.manifest,failed,invalid,1);}catch(IllegalStateException expected){
							final Throwable cause=expected.getCause();
							check(cause instanceof IllegalArgumentException && cause.getMessage()!=null && cause.getMessage().contains(reason),
								"Held-out failure must identify its configured query range: "+expected);rejected=true;
						}
						check(rejected && Files.isRegularFile(failed.resolve("FROZEN")) &&
							!Files.exists(failed.resolve("PASS")) && !Files.exists(failed.resolve("heldout/PASS")),
							"Out-of-range held-out queries must fail after calibration without publishing evaluation success");
					}
				}
			}
		}
		System.out.println("SELECTED_STRUCTURAL_PREPARATION_PASS geometries=3 views=2");
	}

	/** Filtering changes sketches, never the original window/member/panel denominators. */
	private static void verifyFilteredResident(Path root) throws Exception{
		for(int scenario=0;scenario<3;scenario++){
			final int k=scenario==0 ? 5 : 6,blacklistK=scenario==2 ? 6 : 5;
			final String mode=scenario==1 ? "prefix" : "literal",label="filtered_"+scenario;
			final Fixture f=new Fixture(root.resolve(label+"_inputs"));
			configureBlacklist(f,k,blacklistK,mode);
			// Keep a separate immutable manifest for the unfiltered gate control.
			final Path control=f.root.resolve("unfiltered.manifest"),controlOut=root.resolve(label+"_control");
			final StringBuilder unfiltered=new StringBuilder();
			for(String line:text(f.manifest).split("\n")){if(!line.startsWith("blacklist")){unfiltered.append(line).append('\n');}}
			write(control,unfiltered.toString());run(control,controlOut,f.positive,3);
			check(value(controlOut.resolve("calibration.tsv"),"status").equals("ddl_floor_exceeds_target") &&
				value(controlOut.resolve("heldout_gate.tsv"),"allowed").equals("false") &&
				!Files.exists(controlOut.resolve("heldout")),"Unfiltered A negatives must fail before positive evaluation");
			verifySummary(controlOut.resolve("summary.tsv"),NAMES,new long[]{2,2,2,2},"negative",true);

			final Path out=root.resolve(label+"_complete");run(f.manifest,out,f.positive,3);
			check(Files.isRegularFile(out.resolve("PASS")) && Files.isRegularFile(out.resolve("FROZEN")) &&
				Files.isRegularFile(out.resolve("heldout/PASS")),"Filtered resident preparation, calibration and held-out scoring must complete");
			check(!Files.exists(out.resolve("initial")),"Filtered preparation must keep its resident payload without a provisional reload");
			verifyPopulation(out);verifyProvenance(out,f,true);
			verifySummary(out.resolve("summary.tsv"),NAMES,new long[]{2,2,2,2},"negative",false);
			check(value(out.resolve("heldout_gate.tsv"),"allowed").equals("true") &&
				value(out.resolve("heldout_status.tsv"),"status").equals("completed"),"Filtered negatives must pass before held-out scoring");
			final List<String> heldout=lines(out.resolve("heldout/summary.tsv"));
			check(heldout.size()==2,"One complete held-out panel must be retained");
			final String[] row=heldout.get(1).split("\t",-1),header=heldout.get(0).split("\t",-1);
			check(row.length==header.length && row[0].equals("positive") && row[1].equals("positive") &&
				row[2].equals("3") && row[3].equals("1") && row[4].equals("2") && row[9].equals("2") && row[10].equals("1"),
				"A/F/C held-out population must retain one F shortcut and two DDL rejections");
			for(int stage=11;stage<row.length;stage++){check(row[stage].equals("0"),"Alignment must not mask the filtered DDL control");}
			check(value(out.resolve("heldout/scoring_time.tsv"),"queries").equals("3"),"Filtering must not shrink the held-out timing denominator");

			final List<String> profiles=lines(out.resolve("profiles/profiles.tsv.gz"));
			check(profiles.size()==3,"Both live families must survive, including the fully excluded A sketch");
			for(int family=0;family<2;family++){
				final String[] profile=profiles.get(family+1).split("\t",-1);
				check(profile.length==10 && profile[1].equals("f"+family) && profile[3].equals("2") &&
					profile[4].equals("200") && Long.parseLong(profile[5])==2L*(101-k),
					"Excluded windows still count toward the original member/residue/window population");
				final String[] sketch=profile[7].split(",",-1);check(sketch.length==1024,"Retain every DDL bucket");
				boolean nonzero=false;for(String bucket:sketch){nonzero|=Integer.parseInt(bucket)!=0;}
				check(nonzero==(family==1),"Only A sketch contents must be removed; the F sketch must remain populated");
			}
			final List<String> metadata=lines(out.resolve("profiles/summary.tsv"));
			final String[] keys=metadata.get(0).split("\t",-1),values=metadata.get(1).split("\t",-1);
			check(keys.length==values.length,"Filtered profile metadata schema must align");
			final HashMap<String,String> fields=new HashMap<String,String>();
			for(int i=0;i<keys.length;i++){check(fields.put(keys[i],values[i])==null,"Duplicate profile metadata key");}
			check(Long.toString(4L*(101-k)).equals(fields.get("valid_words")) && "unmeasured".equals(fields.get("blacklist_sha80")) &&
				"1".equals(fields.get("blacklist_words")) && Integer.toString(blacklistK).equals(fields.get("blacklist_k")) &&
				mode.equals(fields.get("blacklist_mode")) && "resident_fixture_A".equals(fields.get("blacklist_provenance")),
				"Exact filter policy and unchanged raw-window denominator must survive preparation metadata");
			final String recipe=text(out.resolve("classifier.config"));
			for(String binding:new String[]{"blacklist="+f.root.resolve("blacklist.faa"),"blacklistsha80=unmeasured",
				"blacklistk="+blacklistK,"blacklistmode="+mode,"blacklistwords=1","blacklistprovenance=resident_fixture_A"}){
				check(("\n"+recipe).contains("\n"+binding+"\n"),"Resident classifier lost filter policy: "+binding);
			}
			check("unmeasured".equals(pairs(out.resolve("preparation_inputs.tsv")).get(f.root.resolve("blacklist.faa").toString())),
				"A structural blacklist must retain explicit unmeasured input provenance");
		}
	}

	/** Reject malformed policies and actual list contents through the real resident entry point. */
	private static void verifyBlacklistRejections(Path root) throws Exception{
		final String[] mutations={"missing_identity","missing_policy","partial_policy","orphan_policy","missing_file",
			"wrong_count","wrong_length","duplicate_word","wrong_geometry","prefix6"};
		final String[] reasons={"supplied together","complete explicit blacklist policy","specified together","missing list",
			"Preparation requires a regular file and valid declared identity","word count differs","exactly 6 amino acids","Duplicate blacklist word",
			"does not match raw geometry span","exclusion policies"};
		for(int scenario=0;scenario<mutations.length;scenario++){
			final String mutation=mutations[scenario];final Fixture f=new Fixture(root.resolve("blacklist_bad_"+mutation));
			configureBlacklist(f,6,6,"literal");String recipe=text(f.manifest);
			if(mutation.equals("missing_identity")){recipe=recipe.replace("blacklistsha80=unmeasured\n","");}
			else if(mutation.equals("missing_policy")){
				for(String field:new String[]{"blacklistk=6","blacklistmode=literal","blacklistwords=1","blacklistprovenance=resident_fixture_A"}){
					recipe=recipe.replace(field+"\n","");
				}
			}else if(mutation.equals("partial_policy")){recipe=recipe.replace("blacklistmode=literal\n","");}
			else if(mutation.equals("orphan_policy")){
				recipe=recipe.replace("blacklist="+f.root.resolve("blacklist.faa")+"\n","").replace("blacklistsha80=unmeasured\n","");
			}else if(mutation.equals("missing_file")){recipe=recipe.replace("blacklist="+f.root.resolve("blacklist.faa"),"blacklist="+f.root.resolve("missing.faa"));}
			else if(mutation.equals("wrong_count")){recipe=recipe.replace("blacklistwords=1\n","blacklistwords=2\n");}
			else if(mutation.equals("wrong_length")){write(f.root.resolve("blacklist.faa"),">bad\nAAAAA\n");}
			else if(mutation.equals("duplicate_word")){write(f.root.resolve("blacklist.faa"),">a\nAAAAAA\n>b\nAAAAAA\n");}
			else if(mutation.equals("wrong_geometry")){recipe=recipe.replace("\nk=6\n","\nk=5\n");}
			else if(mutation.equals("prefix6")){recipe=recipe.replace("blacklistmode=literal\n","blacklistmode=prefix\n");}
			write(f.manifest,recipe);
			final Path out=root.resolve("blacklist_rejected_"+mutation);
			reject(()->run(f.manifest,out,f.positive,3),reasons[scenario]);
			check(!Files.exists(out),"Invalid blacklist policy/list must reject before publishing preparation outputs: "+mutation);
		}
	}

	private static void configureBlacklist(Fixture f,int k,int blacklistK,String mode) throws Exception{
		f.writeManifest(f.manifest,true,4);
		final String recipe=text(f.manifest).replace("\nk=5\n","\nk="+k+"\n")
			.replace("\nddlreject=0\n","\nddlreject=1\n").replace("\nddlaccept=1\n","\nddlaccept=2\n");
		for(int panel=0;panel<NAMES.length;panel++){
			write(f.panels[panel],">fixture_"+NAMES[panel]+"_0\n"+sequence('A')+"\n>fixture_"+NAMES[panel]+"_1\n"+sequence('A')+"\n");
		}
		write(f.root.resolve("blacklist.faa"),">excluded_A\n"+(blacklistK==5 ? "AAAAA" : "AAAAAA")+"\n");
		write(f.manifest,recipe+"blacklist="+f.root.resolve("blacklist.faa")+"\nblacklistsha80=unmeasured\nblacklistk="+blacklistK+
			"\nblacklistmode="+mode+"\nblacklistwords=1\nblacklistprovenance=resident_fixture_A\n");
	}

	/** Same observations: pooled100/M converges, while the scientific per-class criterion fails at400/M. */
	private static void verifyPooledGate(Path root,Path unreadable) throws Exception{
		final Fixture f=new Fixture(root.resolve("pooled_inputs"));
		f.writeManifest(f.manifest,true,4);
		// Two matching positional words suffice for the A shortcut; zero-hit C tails reject immediately.
		write(f.manifest,text(f.manifest).replace("ddlreject=0\n","ddlreject=1\n").replace("ddlaccept=1\n","ddlaccept=2\n"));
		final String hit=sequence('A'),miss=sequence('C');
		final StringBuilder table=new StringBuilder("panel\tpath\tsha80\tqueries\n");
		for(int p=0;p<NAMES.length;p++){
			final StringBuilder fasta=new StringBuilder();
			for(int q=0;q<2500;q++){
				fasta.append(">fixture_pooled_").append(NAMES[p]).append('_').append(q).append('\n')
					.append(p==0 && q==0 ? hit : miss).append('\n');
			}
			write(f.panels[p],fasta.toString());
			table.append(NAMES[p]).append('\t').append(f.panels[p]).append('\t').append(INHERITED).append("\t2500\n");
		}
		write(f.root.resolve("negatives.tsv"),table.toString());
		final Path out=root.resolve("pooled_converged_class_failed");
		FamilyCalibrate.main(expectedArguments(f.manifest,out,unreadable,false,"1,0,0,0",2500));
		check(value(out.resolve("calibration.tsv"),"status").equals("converged") &&
			value(out.resolve("calibration.tsv"),"target_scope").equals("pooled"),"The pooled criterion must converge at exactly100/M");
		check(Files.isRegularFile(out.resolve("FROZEN")) && Files.isRegularFile(out.resolve("PASS")),
			"A completed pooled calibration must retain its frozen negative result despite skipping positives");
		final List<String> rounds=lines(out.resolve("rounds.tsv")),summary=lines(out.resolve("summary.tsv"));
		check(rounds.size()==2 && summary.size()==5,"One evaluation must consume exactly four complete panels");
		final String[] round=rounds.get(1).split("\t",-1);
		check(round[1].equals("10000") && round[2].equals("1") && round[4].equals("1") && round[5].equals("0") &&
			round[6].equals("converged"),"Exactly one shortcut among10000 negatives must converge without a gate update");
		final String[] header=summary.get(0).split("\t",-1);
		check(header[9].equals("ddl_reject") && header[10].equals("ddl_accept"),"Unexpected pooled summary schema");
		for(int p=0;p<NAMES.length;p++){
			final String[] row=summary.get(p+1).split("\t",-1);final long accepts=p==0 ? 1 : 0;
			check(row.length==header.length && row[0].equals(NAMES[p]) && row[1].equals("negative") && row[2].equals("2500"),
				"Pooled calibration must retain all2500 queries in every named class");
			check(Long.parseLong(row[3])==accepts && Long.parseLong(row[4])==2500-accepts &&
				Long.parseLong(row[9])==2500-accepts && Long.parseLong(row[10])==accepts,
				"Only the one member-identical A tail may shortcut; every other tail must be a DDL reject");
			for(int stage=11;stage<row.length;stage++){check(Long.parseLong(row[stage])==0,"No alignment or other terminal outcome may mask the intended gate branch");}
			check(Double.parseDouble(row[8])==(p==0 ? 400 : 0),"The leaking class must report400 accepts per million, the others zero");
		}
		check(value(out.resolve("heldout_gate.tsv"),"allowed").equals("false") &&
			value(out.resolve("heldout_gate.tsv"),"reason").equals("negative_class_exceeds_ceiling") &&
			value(out.resolve("heldout_gate.tsv"),"calibration_status").equals("converged"),
			"The separate per-class gate must override pooled convergence before positive-file access");
		check(value(out.resolve("heldout_status.tsv"),"status").equals("skipped") &&
			value(out.resolve("heldout_status.tsv"),"reason").equals("negative_class_exceeds_ceiling") &&
			!Files.exists(out.resolve("heldout")),"Failed class ceiling must skip the same unreadable positive path without an evaluation directory");
		check(text(unreadable.getParent()).equals("positive-open sentinel\n"),"Pooled gate failure must preserve the positive sentinel");

		// The scientific recipe uses perclass=t; retain pooled coverage as a separate gate regression.
		final Path perClass=root.resolve("per_class_floor_failed");
		FamilyCalibrate.main(arguments(f.manifest,perClass,unreadable,3,"structural",true));
		check(value(perClass.resolve("calibration.tsv"),"target_scope").equals("per_class") &&
			value(perClass.resolve("calibration.tsv"),"status").equals("ddl_floor_exceeds_target"),
			"The400/M class must fail scientific calibration despite the pooled100/M rate");
		check(lines(perClass.resolve("summary.tsv")).equals(summary),
			"Changing target scope must preserve all four2500-query populations and the same single DDL shortcut");
		final List<String> perClassRounds=lines(perClass.resolve("rounds.tsv"));
		check(perClassRounds.size()==2,"An unmodifiable per-class shortcut floor must stop after one evaluation");
		final String[] failedRound=perClassRounds.get(1).split("\t",-1);
		check(failedRound[1].equals("10000") && failedRound[2].equals("1") && failedRound[4].equals("1") &&
			failedRound[5].equals("0") && failedRound[6].equals("ddl_floor_exceeds_target"),
			"Per-class failure must retain the full measured population without applying an ineffective update");
		check(Files.isRegularFile(perClass.resolve("PASS")) && !Files.exists(perClass.resolve("FROZEN")) &&
			!Files.exists(perClass.resolve("heldout")),"Completed but failed per-class calibration must neither freeze nor evaluate positives");
		check(value(perClass.resolve("heldout_gate.tsv"),"allowed").equals("false") &&
			value(perClass.resolve("heldout_gate.tsv"),"reason").equals("calibration_not_converged") &&
			value(perClass.resolve("heldout_gate.tsv"),"calibration_status").equals("ddl_floor_exceeds_target") &&
			value(perClass.resolve("heldout_status.tsv"),"status").equals("skipped"),
			"Scientific per-class failure must persist its skip before accessing the same unreadable positive path");
		check(text(unreadable.getParent()).equals("positive-open sentinel\n"),"Per-class gate failure must preserve the positive sentinel");
		verifyExpectedCounts(root,f,unreadable);
	}

	/** A reordered first-round vector or wrong complete-panel population must fail before updates. */
	private static void verifyExpectedCounts(Path root,Fixture f,Path unreadable) throws Exception{
		for(int scenario=0;scenario<3;scenario++){
			final Path out=root.resolve("expected_counts_failed_"+scenario),log=root.resolve("expected_counts_failed_"+scenario+".log");
			final long population=scenario==0 ? 2500 : scenario==1 ? 2499 : 2501;
			final String[] base=expectedArguments(f.manifest,out,unreadable,true,scenario==0 ? "0,1,0,0" : "1,0,0,0",population);
			final String[] argv=Arrays.copyOf(base,base.length+1);argv[base.length]="outstream="+log;
			final String reason=scenario==0 ? "Initial calibration accepts changed before gate update: panel=k5 expected=0 observed=1" :
				"Calibration panel population changed before gate update: round=1 panel=k5 expected="+population+" observed=2500";
			reject(()->FamilyCalibrate.main(argv),reason);
			final String transcript=text(log);
			check(!transcript.contains("Calibration round 2") && !transcript.contains("; raise."),
				"An initial-count mismatch must fail before applying a gate update or entering round2");
			for(String marker:new String[]{"PASS","FROZEN","rounds.tsv","cohorts.tsv","calibration.tsv","heldout_gate.tsv","heldout_status.tsv","heldout"}){
				check(!Files.exists(out.resolve(marker)),"Initial-count failure must not publish calibration/held-out completion: "+marker);
			}
			check(text(unreadable.getParent()).equals("positive-open sentinel\n"),"Initial-count failure must not access the positive panel");
		}
		final Path malformed=root.resolve("expected_vector_length_failed");
		reject(()->new FamilyCalibrate(expectedArguments(f.manifest,malformed,unreadable,true,"1,0,0",2500)),
			"expectedinitialaccepts must have one count per negative panel, in panel order");
		check(!Files.exists(malformed),"Wrong vector length must fail during construction, before preparation");
	}

	private static String[] expectedArguments(Path manifest,Path out,Path positive,boolean perClass,String accepts,long population){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(arguments(manifest,out,positive,3,"structural",perClass)));
		final int limit=args.indexOf("rounds=1");check(limit>=0,"Expected one explicit fixture round limit");args.set(limit,"rounds=4");
		args.add("expectedinitialaccepts="+accepts);args.add("expectedpanelqueries="+population);
		return args.toArray(new String[args.size()]);
	}

	/** Synthetic 100-class map covers all400 dimers; it does not stand in for the empirical maps. */
	private static void verifyMetricResident(Path root) throws Exception{
		final Fixture f=new Fixture(root.resolve("metric_inputs"));
		final String map=syntheticDimerMap();f.writeMetricManifest(map);
		final Path out=root.resolve("metric_complete");run(f.manifest,out,f.positive,3);
		check(Files.isRegularFile(out.resolve("PASS")) && Files.isRegularFile(out.resolve("heldout/PASS")),
			"Metric-aware resident calibration and held-out evaluation must complete");
		verifyPopulation(out);verifyProvenance(out,f);
		verifySummary(out.resolve("summary.tsv"),NAMES,new long[]{2,2,2,2},"negative",false);
		verifySummary(out.resolve("heldout/summary.tsv"),new String[]{"positive"},new long[]{3},"positive",false);
		final String recipe=text(out.resolve("classifier.config"));
		check(recipe.contains("aabits=dimer:1,2,3,4:"+map+"\n") && recipe.contains("width=512\n") &&
			recipe.contains("alignmentmode=pin_right_trim_left_aligned_query\n") && recipe.contains("alignmentminspan=70\n") &&
			recipe.contains("heuristic=hitfourth_h020_c150\n"),"Metric, dimer geometry and common-pivot ranking must survive resident preparation");
		final List<String> raw=lines(out.resolve("raw_calibration.tsv.gz"));
		for(int family=0;family<2;family++){
			final String[] row=raw.get(family+1).split("\t",-1);
			check(row.length==46 && row[12].equals("pin_right_trim_left_aligned_query") && row[13].equals("70") &&
				row[14].equals("2") && row[15].equals("0"),"Both members must contribute measured alignment evidence in the selected metric");
			for(int p=0;p<NAMES.length;p++){
				check(row[18+7*p].equals("2") && row[19+7*p].equals("0") && row[20+7*p].equals("0"),
					"Each metric panel must scan both negatives and admit none at the512-hit floor");
			}
		}
		final Path mismatch=root.resolve("metric_count_mismatch");
		reject(()->run(f.manifest,mismatch,f.positive,4),"count differs from its declared population");
		check(value(mismatch.resolve("heldout_gate.tsv"),"allowed").equals("true") &&
			!Files.exists(mismatch.resolve("PASS")) && !Files.exists(mismatch.resolve("heldout/PASS")),
			"Metric-aware held-out count mismatch must fail after the gate without evaluation completion");
	}

	/** Compare fresh collection/heuristic inputs to the identical inherited fixture calculation. */
	private static void verifyUnmeasuredResident(Path root)throws Exception{
		final Fixture f=new Fixture(root.resolve("unmeasured_inputs"));f.writeMetricManifest(syntheticDimerMap());
		String recipe=text(f.manifest);
		for(String key:new String[]{"families","chunks","consensus","lengthmedians","familyhydro"}){
			recipe=recipe.replace(key+"sha80="+INHERITED+"\n",key+"sha80=unmeasured\n");
		}
		write(f.manifest,recipe);
		write(f.root.resolve("chunks.tsv"),text(f.root.resolve("chunks.tsv")).replace(INHERITED,"unmeasured"));
		final Path out=root.resolve("unmeasured_complete");run(f.manifest,out,f.positive,3);
		check(Files.isRegularFile(out.resolve("PASS")) && Files.isRegularFile(out.resolve("FROZEN")),"Fresh unmeasured collection must finish resident calibration");
		check(value(out.resolve("verification.tsv"),"input_identity").equals("inherited_with_unmeasured_inputs"),
			"Fresh unmeasured collection and heuristic tables must never be labeled all inherited");
		verifyPopulation(out);
		for(String name:new String[]{"summary.tsv","references.tsv.gz","raw_calibration.tsv.gz","profiles/profiles.tsv.gz","heldout/summary.tsv"}){
			check(lines(out.resolve(name)).equals(lines(root.resolve("metric_complete").resolve(name))),
				"Identity labels must not change scientific output bytes: "+name);
		}
		final HashMap<String,String> inputs=pairs(out.resolve("preparation_inputs.tsv"));
		for(String name:new String[]{"families.tsv","chunks.tsv","consensus.faa","medians.tsv","hydro.tsv"}){
			check("unmeasured".equals(inputs.get(f.root.resolve(name).toString())),"Fresh input provenance must remain unmeasured: "+name);
		}
		final prot.FamilyScoringPackage fresh=new prot.FamilyScoringPackage(out.resolve("classifier.config"),VerificationMode.STRUCTURAL);
		final prot.FamilyScoringPackage inherited=new prot.FamilyScoringPackage(root.resolve("metric_complete/classifier.config"),VerificationMode.STRUCTURAL);
		final prot.FamilyTailCascade.LoadedClassifier freshClassifier=fresh.load(1),inheritedClassifier=inherited.load(1);
		for(char aa:new char[]{'A','F','C'}){
			final byte[] query=sequence(aa).getBytes(StandardCharsets.US_ASCII);
			check(freshClassifier.newWorker().assign("q",query).assigned()==inheritedClassifier.newWorker().assign("q",query).assigned(),
				"Serialized structural package loading must preserve ordinary assignments with unmeasured metadata");
		}
		fresh.verify();inherited.verify();
		reject(()->new prot.FamilyScoringPackage(out.resolve("classifier.config")),"explicit structural loading");
		final Path invalidPackage=f.root.resolve("invalid_metadata.config");
		write(invalidPackage,text(out.resolve("classifier.config")).replace("lengthmedianssha80=unmeasured","lengthmedianssha80=unknown"));
		reject(()->new prot.FamilyScoringPackage(invalidPackage,VerificationMode.STRUCTURAL),"valid declared identity");
		final Path malformed=f.root.resolve("malformed.manifest");
		write(malformed,recipe.replace("familiessha80=unmeasured","familiessha80=unknown"));
		final Path rejected=root.resolve("unmeasured_malformed");
		reject(()->run(malformed,rejected,f.positive,3),"valid declared identity");
		check(!Files.exists(rejected.resolve("PASS")),"Malformed provenance must not publish a successful package");
	}

	private static String syntheticDimerMap(){
		final String alphabet="ARNDCEQGHILKMFPSTWYV";final StringBuilder groups=new StringBuilder();
		for(int pair=0;pair<400;pair++){
			if(pair>0){groups.append(pair%4==0 ? '/' : ',');}
			groups.append(alphabet.charAt(pair/20)).append(alphabet.charAt(pair%20));
		}
		return groups.toString();
	}

	private static void verifyPopulation(Path out) throws Exception{
		final List<String> rounds=lines(out.resolve("rounds.tsv")),cohorts=lines(out.resolve("cohorts.tsv"));
		check(rounds.size()==2 && cohorts.size()==5,"Exactly one complete evaluation must cover all four panels");
		final String[] round=rounds.get(1).split("\t",-1);
		check(round[1].equals("8") && round[2].equals("0") && round[6].equals("converged"),"Round denominator must sum all negative classes");
		for(int p=0;p<NAMES.length;p++){
			final String[] row=cohorts.get(p+1).split("\t",-1);
			check(row[0].equals("1") && row[1].equals(NAMES[p]) && row[2].equals("2") && row[3].equals("0"),
				"Every negative class must appear once with its full population");
		}
		final List<String> prep=lines(out.resolve("preparation.tsv"));
		check(prep.size()==7,"Preparation must record members, four negative classes and total");
		for(int p=0;p<NAMES.length;p++){
			final String[] row=prep.get(p+2).split("\t",-1);
			check(row[0].equals("negative") && row[1].equals(NAMES[p]) && row[2].equals("2") && row[3].equals("2"),
				"Preparation must conserve each panel and both live families");
		}
		final List<String> raw=lines(out.resolve("raw_calibration.tsv.gz"));
		check(raw.size()==4,"Raw preparation must retain all three physical family slots including the empty one");
		for(int i=0;i<3;i++){
			final String[] row=raw.get(i+1).split("\t",-1);
			check(row[0].equals(Integer.toString(i)) && row[1].equals(i==2 ? "empty" : "f"+i) &&
				row[2].equals(i==2 ? "0" : "2"),"Preparation lost an assigned member or changed the family roster");
		}
		final List<String> summary=lines(out.resolve("profiles/summary.tsv"));
		check(summary.size()==2,"Exactly one complete profile population must be retained");
		final String[] population=summary.get(1).split("\t",-1);
		check(population[5].equals("3") && population[6].equals("3") && population[7].equals("4") &&
			population[8].equals("400"),"Profile metadata must conserve the physical roster, four members and 400 residues");
	}

	private static void verifySummary(Path path,String[] names,long[] expected,String kind,boolean shortcuts){
		final List<String> rows=lines(path);
		check(rows.size()==names.length+1,"Summary must contain exactly the declared panels");
		final String[] header=rows.get(0).split("\t",-1);
		check(header.length>10 && header[9].equals("ddl_reject") && header[10].equals("ddl_accept"),"Unexpected terminal-stage summary schema");
		for(int p=0;p<names.length;p++){
			final String[] row=rows.get(p+1).split("\t",-1);
			check(row.length==header.length && row[0].equals(names[p]) && row[1].equals(kind),"Summary panel identity/order changed");
			final long queries=Long.parseLong(row[2]),accepted=Long.parseLong(row[3]),rejected=Long.parseLong(row[4]);
			check(queries==expected[p] && accepted+rejected==queries,"Panel population must be completely partitioned");
			long terminal=0;for(int s=9;s<row.length;s++){terminal=Math.addExact(terminal,Long.parseLong(row[s]));}
			check(terminal==queries,"Exactly one terminal stage is required for every query");
			check(shortcuts ? accepted>0 && Long.parseLong(row[10])==accepted : accepted==0 && Long.parseLong(row[9])==queries,
				"Fixture must exercise its intended DDL gate, not an incidental alignment outcome");
		}
	}

	private static void verifyProvenance(Path out,Fixture f) throws Exception{
		verifyProvenance(out,f,false);
	}
	private static void verifyProvenance(Path out,Fixture f,boolean unmeasuredBlacklist) throws Exception{
		final Path verification=out.resolve("verification.tsv");
		final String inputIdentity=unmeasuredBlacklist ? "inherited_with_unmeasured_blacklist" : "inherited";
		check(value(verification,"verification").equals("structural") && value(verification,"input_identity").equals(inputIdentity) &&
			value(verification,"output_byte_identity").equals("unmeasured") && value(verification,"serialized_reload_validated").equals("false"),
			"Provenance must distinguish inherited inputs, unmeasured output bytes and unvalidated reloads");
		if(unmeasuredBlacklist){
			check(value(verification,"blacklist_input_identity").equals("unmeasured"),"New structural blacklist identity must not be reported as inherited");
		}else{check(!pairs(verification).containsKey("blacklist_input_identity"),"Unfiltered verification must retain its original metadata shape");}
		check(value(verification,"query_identity_checks").contains("no_independent_ID_sidecar_join"),"Reader order must not be described as an independent ordered-ID proof");
		final String recipe=text(out.resolve("classifier.config"));
		check(recipe.contains("verification=structural\n") && recipe.contains("familyreferencessha80=unmeasured\n") &&
			recipe.contains("ddlprofilessha80=unmeasured\n"),"Derived classifier artifacts must not borrow inherited input digests");
		check(text(out.resolve("FROZEN")).equals(out.resolve("classifier.config")+"\tunmeasured\n"),"FROZEN must not fabricate an output fingerprint");
		check(text(out.resolve("ddlprofiles.tsv")).equals("profiles/profiles.tsv.gz\tunmeasured\tunmeasured\n"),"Derived profiles and summary must remain explicitly unmeasured");
		final HashMap<String,String> audit=pairs(out.resolve("preparation_inputs.tsv"));
		for(String relative:new String[]{"references.tsv.gz","profiles/profiles.tsv.gz","profiles/summary.tsv","raw_calibration.tsv.gz","ddlprofiles.tsv","classifier.config"}){
			check("unmeasured".equals(audit.get(out.resolve(relative).toString())),"Output provenance must be unmeasured: "+relative);
		}
		check("unmeasured".equals(audit.get(f.manifest.toString())) && INHERITED.equals(audit.get(f.families.toString())),
			"Recipe identity and inherited roster identity must not be conflated");
		final List<String> inputs=lines(out.resolve("inputs.tsv"));
		check(inputs.size()==6 && inputs.get(1).equals("package\t"+out.resolve("classifier.config")+"\tunmeasured"),"Calibration package provenance must remain explicit");
		for(int p=0;p<NAMES.length;p++){
			check(inputs.get(p+2).equals(NAMES[p]+"\t"+f.panels[p]+"\t"+INHERITED),"Each negative must retain its inherited identity");
		}
		check(text(out.resolve("heldout/inputs.tsv")).equals("panel\tpath\tinherited_sha80\npositive\t"+f.positive+"\t"+INHERITED+"\n"),
			"Held-out provenance must label its identity as inherited");
		check(value(out.resolve("heldout/evaluation_policy.tsv"),"serialized_reload_validated").equals("false") &&
			value(out.resolve("heldout/evaluation_policy.tsv"),"positive_guided_selection").equals("false"),"Held-out reporting must not claim reload parity or positive-guided selection");
	}

	private static final class Fixture{
		Fixture(Path root_) throws Exception{
			root=root_;Files.createDirectory(root);
			families=root.resolve("families.tsv");manifest=root.resolve("prepare.manifest");positive=root.resolve("positive.faa");
			write(root.resolve("SYNTHETIC.txt"),"All identities are synthetic inherited placeholders, not computed fingerprints.\n"+
				"Panel labels model four independent input classes; deterministic tiny peptides are not biological controls.\n");
			write(families,DDLFamilyHoldout.FAMILY_HEADER+"\n0\tf0\t2\n1\tf1\t2\n2\tempty\t0\n");
			write(root.resolve("consensus.faa"),">f0\n"+sequence('A')+"\n>f1\n"+sequence('F')+"\n>empty\nX\n");
			final Path part=root.resolve("members.tsv");
			write(part,DDLFamilyHoldout.MEMBER_HEADER+"\n"+member(0,0,'A')+member(1,1,'F')+member(0,2,'A')+member(1,3,'F'));
			write(root.resolve("chunks.tsv"),DDLFamilyHoldout.CHUNK_HEADER+"\nb00\t0\t"+part+"\t"+INHERITED+"\t4\t400\n");
			final StringBuilder panelTable=new StringBuilder("panel\tpath\tsha80\tqueries\n");
			for(int p=0;p<NAMES.length;p++){
				panels[p]=root.resolve(NAMES[p]+".faa");writePanel(p,2,false);
				panelTable.append(NAMES[p]).append('\t').append(panels[p]).append('\t').append(INHERITED).append("\t2\n");
			}
			write(root.resolve("negatives.tsv"),panelTable.toString());
			write(positive,">fixture_positive_0\n"+sequence('A')+"\n>fixture_positive_1\n"+sequence('F')+"\n>fixture_positive_2\n"+sequence('C')+"\n");
			final StringBuilder bounds=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
			for(String scalar:new String[]{"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"}){
				// The production scalar schema requires this historical metadata census.
				// It is synthetic schema data here, not the fixture member/family population.
				bounds.append(scalar).append("\t79565\t-20\t20\tf0\tf1\n");
			}
			write(root.resolve("scalar_bounds.tsv"),bounds.toString());writeManifest(manifest,false,4);
		}
		void writePanel(int panel,int count,boolean duplicate) throws Exception{
			final StringBuilder text=new StringBuilder();
			for(int i=0;i<count;i++){text.append(">fixture_").append(NAMES[panel]).append('_').append(duplicate ? 0 : i).append('\n').append(sequence(i%2==0 ? 'A' : 'F')).append('\n');}
			write(panels[panel],text.toString());
		}
		void writeManifest(Path path,boolean shortcut,long members) throws Exception{
			final StringBuilder recipe=new StringBuilder("# Synthetic inherited input identities; no measured file digests.\nschema=family_preparation_v1\n");
			for(String[] input:new String[][]{{"families","families.tsv"},{"chunks","chunks.tsv"},{"consensus","consensus.faa"},
				{"negativepanels","negatives.tsv"},{"scalarbounds","scalar_bounds.tsv"}}){
				recipe.append(input[0]).append('=').append(root.resolve(input[1])).append('\n').append(input[0]).append("sha80=").append(INHERITED).append('\n');
			}
			recipe.append("cosinepanel=k5\nexpectedfamilies=3\nexpectedmembers=").append(members).append("\nbuckets=1\nbatches=1\n")
				.append("maxreferenceaa=1000\nmaxbucketmembers=16\nmaxbucketbytes=16384\nmaxpreparedbytes=1048576\nmaxoutputbytes=1048576\n")
				.append("retention=metric_specific\nuppertop=2\nminmemberaa=100\nmaxmemberaa=100\nminnegativeaa=100\nmaxnegativeaa=100\n")
				.append("k=5\nseed=0\nwidth=1024\npositional=t\npartitions=4\noverlap=2\n")
				.append("scalarscale=0.1\nblosumfloor=0\nscalarz=0\nmonomer=t\ntop=50\nddlreject=").append(shortcut ? 0 : 1024)
				.append("\nddlaccept=").append(shortcut ? 1 : 1025).append("\ncosine=ranking\nrankingfeatures=combined25\nscalargates=f\n")
				.append("identityacceptfloor=30\nblosumacceptfloor=0\nindexearlyreject=f\n")
				.append("matrix=blosum62\n");
			write(path,recipe.toString());
		}
		void writeMetricManifest(String map) throws Exception{
			final Path medians=root.resolve("medians.tsv"),hydro=root.resolve("hydro.tsv");
			write(medians,"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n"+
				"f0\t2\t100\t100\t100\t2\nf1\t2\t100\t100\t100\t2\n");
			write(hydro,"family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tf0\t2\t1.8\n1\tf1\t2\t2.8\n2\tempty\t0\tNA\n");
			final String recipe=text(manifest).replace("width=1024\n","width=512\n")
				.replace("ddlreject=1024\n","ddlreject=512\n").replace("matrix=blosum62\n","");
			write(manifest,recipe+"aabits=dimer:1,2,3,4:"+map+"\nalignmentmode=pin_right_trim_left_aligned_query\n"+
				"alignmentminspan=70\nheuristic=hitfourth_h020_c150\nminqueryaa=100\nmaxqueryaa=100\nlengthmedians="+medians+"\nlengthmedianssha80="+INHERITED+
				"\nfamilyhydro="+hydro+"\nfamilyhydrosha80="+INHERITED+"\n");
		}
		private String member(int family,int ordinal,char aa){return family+"\tf"+family+"\tg|q"+ordinal+"\tg\t"+ordinal+"\t"+sequence(aa)+"\n";}
		final Path root,families,manifest,positive;
		final Path[] panels=new Path[NAMES.length];
	}

	private static String[] arguments(Path manifest,Path out,Path positive,long count,String mode){
		return arguments(manifest,out,positive,count,mode,true);
	}
	private static String[] arguments(Path manifest,Path out,Path positive,long count,String mode,boolean perClass){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("prepare="+manifest,"out="+out,"t=1","tin=1",
			"rounds=1","targetfpr=.0001","perclass="+perClass,"heldout="+positive,"heldoutsha80="+INHERITED,"heldoutqueries="+count));
		if(mode!=null){args.add("verification="+mode);}return args.toArray(new String[args.size()]);
	}
	private static void run(Path manifest,Path out,Path positive,long count) throws Exception{
		FamilyCalibrate.main(arguments(manifest,out,positive,count,"structural"));
	}
	private static List<String> lines(Path path){
		final ArrayList<String> rows=new ArrayList<String>();final ByteFile input=ByteFile.makeByteFile(path.toString(),false);
		try{for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){rows.add(new String(line,StandardCharsets.US_ASCII));}}
		finally{check(!input.close(),"Fixture output read failed: "+path);}return rows;
	}
	private static HashMap<String,String> pairs(Path path){
		final List<String> rows=lines(path);check(!rows.isEmpty(),"Missing table header: "+path);
		final HashMap<String,String> values=new HashMap<String,String>();
		for(int i=1;i<rows.size();i++){
			final String[] row=rows.get(i).split("\t",-1);
			check(row.length==2 && values.put(row[0],row[1])==null,"Malformed/duplicate key in "+path);
		}return values;
	}
	private static String value(Path path,String key){final String result=pairs(path).get(key);check(result!=null,"Missing "+key+" in "+path);return result;}
	private static String text(Path path) throws Exception{return new String(Files.readAllBytes(path),StandardCharsets.US_ASCII);}
	private static String sequence(char aa){final char[] text=new char[100];Arrays.fill(text,aa);return new String(text);}
	private static void write(Path path,String text) throws Exception{Files.write(path,text.getBytes(StandardCharsets.US_ASCII));}
	private static void reject(Action action,String reason) throws Exception{
		try{action.run();}catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(reason),"Wrong rejection, expected "+reason+": "+expected);return;
		}throw new AssertionError("Expected rejection: "+reason);
	}
	private static void check(boolean valid,String why){checks++;if(!valid){throw new AssertionError(why);}}
	private interface Action{void run() throws Exception;}
	private static final String INHERITED="0123456789abcdef0123";
	private static final String[] NAMES={"k5","shadow","reversed","shuffled"};
	private static int checks;
}
