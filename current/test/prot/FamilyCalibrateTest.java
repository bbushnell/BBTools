package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ddl.DDLFamilyProfiles;
import fileIO.ByteFile;
import prot.DigestSuffix;
import prot.FamilyCalibrate;
import prot.FamilyEval;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyReferenceExport;
import prot.FamilyReferenceLoader;
import prot.ProteinSequence;
import prot.ProteinSpectrum;
import structures.ByteBuilder;

/** Two identical families force a second-round winner to emerge after the first is rejected.
 * Tests full-package reuse, exact thread parity, target scope and honest nonconvergence.
 * Uses only synthetic populations; assertion counts belong to one isolated JVM.
 * @author Brian, Keqing
 */
public final class FamilyCalibrateTest{

	/*--------------------------------------------------------------*/
	/*----------------       Calibration checks     ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a fresh synthetic package and checks updater, reload and failure contracts.
	 * The optional legacy=t argument omits only the combined pooled/class checks.
	 * @param args fixture=absolute_directory, optionally followed by legacy=t
	 * @throws Exception if fixture creation, calibration or retained-output reading fails
	 */
	public static void main(String[] args) throws Exception{
		check((args.length==1 || args.length==2 && args[1].equals("legacy=t")) && args[0].startsWith("fixture="),
			"Require fixture=new_absolute_directory [legacy=t]");
		final Path root=Paths.get(args[0].substring(8));
		check(root.isAbsolute() && !Files.exists(root), "Calibration fixture must be fresh");
		Files.createDirectories(root);
		final String hit=repeat("ACDEFGHIKLMNPQRSTVWY", 10), miss=repeat("A", 200);
		final Path references=root.resolve("references.tsv");
		final ByteBuilder records=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA, 5, 1024, 0, 3, 2, 1,
			FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING), records);
		// Physical order is deliberately different from roster rank; the loader permits both.
		FamilyReferenceCodec.appendFiveLines(family("f2", 2, hit), records);
		FamilyReferenceCodec.appendFiveLines(FamilyReference.empty("empty", 1), records);
		FamilyReferenceCodec.appendFiveLines(family("f0", 0, hit), records);
		write(references, records.toString());
		final byte[] originalBytes=Files.readAllBytes(references);
		final Path config=root.resolve("classifier.config");
		final String recipe="familyreferences=references.tsv\nfamilyreferencessha80=unmeasured\n"+
			"top=2\nk=5\nseed=0\nwidth=1024\nddlreject=8\nddlaccept=1025\ncosine=ranking\n"+
			"rankingfeatures=monomer\nscalargates=f\npositional=f\nidentityacceptfloor=30\nblosumacceptfloor=0\nindexearlyreject=f\n";
		write(config, recipe);
		final Path legacyConfig=root.resolve("legacy.config");
		write(legacyConfig, recipe.replace("familyreferencessha80=unmeasured\n", ""));
		final Path sha80Config=root.resolve("sha80.config");
		write(sha80Config, recipe.replace("familyreferencessha80=unmeasured", "familyreferencessha80="+hash(references)));
		final Path challenge=root.resolve("challenge.faa"), clean=root.resolve("clean.faa");
		final StringBuilder a=new StringBuilder(), b=new StringBuilder();
		for(int i=0; i<600; i++){
			fasta(a, "hit"+i, hit);
			fasta(a, "miss"+i, miss);
			fasta(b, "clean"+i, miss);
		}
		write(challenge, a.toString());
		write(clean, b.toString());
		for(int threads:new int[]{1, 4}){
			final Path out=root.resolve("threads"+threads);
			run(config, challenge, clean, out, threads, "targetfpr=0", "rounds=3");
			check(status(out).equals("converged") && Files.exists(out.resolve("FROZEN")), "Third evaluation must converge and freeze");
			final List<String> rounds=Files.readAllLines(out.resolve("rounds.tsv"));
			check(rounds.size()==4, "Two different accepted families require two updates and three evaluations");
			for(int i=1; i<=3; i++){
				final String[] row=rounds.get(i).split("\t");
				check(Long.parseLong(row[1])==1800 && Long.parseLong(row[2])==(i<3 ? 600 : 0),
					"Rejecting the earlier winner must reveal the second family before all 600 false accepts disappear");
				check(Integer.parseInt(row[5])==(i<3 ? 1 : 0), "Each update raises exactly one newly exposed family");
			}
			final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(out.resolve("references.tsv.gz").toString(), 1, 10);
			check(loaded.all.length==3 && loaded.live.length==2 && loaded.all[1].members==0 && loaded.all[1].rank==1,
				"Final package preserves the empty physical roster entry and original ranks");
			check(loaded.all[0].rank==2 && loaded.all[2].rank==0,
				"Publishing updated gates must preserve physical order even when it differs from roster rank");
			for(FamilyReference f:loaded.live){
				check(f.identityUpper==101 && f.blosumFinal==1.1f && f.identityLower==50 && f.blosumLower==.5f,
					"Both identity gates rise above the exact self score; BLOSUM and member lower gates stay unchanged");
			}
			final Path replay=root.resolve("published_replay"+threads);
			FamilyEval.main(new String[]{"in="+out, "pos="+challenge, "neg="+clean, "out="+replay, "t="+threads,
				"verification=structural", "possha80=unmeasured", "negsha80=unmeasured", "poscount=1200", "negcounts=600"});
			check(Files.readAllLines(replay.resolve("summary.tsv")).get(1).split("\t")[3].equals("0"),
				"The published classifier recipe must reload and reproduce the calibrated result");
			check(Files.exists(out.resolve("PASS")), "A completed calibration publishes PASS");
		}
		for(String name:new String[]{"references.tsv.gz", "summary.tsv", "family_false_accepts.tsv.gz", "cohorts.tsv"}){
			check(lines(root.resolve("threads1/"+name)).equals(lines(root.resolve("threads4/"+name))),
				"Thread count must not change gates or observed results: "+name);
		}
		check(Arrays.equals(originalBytes, Files.readAllBytes(references)), "Calibration never rewrites its original package");
		final Path rawDefault=root.resolve("raw_default");
		FamilyCalibrate.main(new String[]{"in="+legacyConfig, "neg="+challenge+","+clean, "negnames=challenge,clean",
			"out="+rawDefault, "t=1", "targetfpr=0", "rounds=3"});
		check(status(rawDefault).equals("converged") &&
			Files.readAllLines(rawDefault.resolve("inputs.tsv"), StandardCharsets.US_ASCII).equals(Arrays.asList(
				"panel\tpath\tsha80", "package\t"+rawDefault.resolve("classifier.config")+"\tunmeasured",
				"challenge\t"+challenge+"\tunmeasured",
				"clean\t"+clean+"\tunmeasured")),
			"Raw default structural calibration must accept omitted panel pins as exact unmeasured provenance");
		final Path sha80=root.resolve("sha80");
		FamilyCalibrate.main(new String[]{"in="+sha80Config, "neg="+challenge+","+clean, "negnames=challenge,clean",
			"negsha80="+hash(challenge)+","+hash(clean), "out="+sha80, "t=1", "targetfpr=0", "rounds=3",
			"verification=sha80"});
		check(status(sha80).equals("converged"), "Explicit SHA80 calibration must still hash concrete inputs");
		final Path pooled=root.resolve("pooled");
		run(config, challenge, clean, pooled, 2, "targetfpr=.4", "rounds=3");
		check(status(pooled).equals("converged") && Files.readAllLines(pooled.resolve("rounds.tsv")).size()==2 &&
			!Files.exists(pooled.resolve("references.tsv.gz")), "Pooled 600/1800 meets .4 without changing or copying references");
		final Path perClass=root.resolve("per_class");
		run(config, challenge, clean, perClass, 2, "targetfpr=.4", "perclass=t", "rounds=3");
		check(status(perClass).equals("converged") && Files.readAllLines(perClass.resolve("rounds.tsv")).size()==4,
			"The challenge class starts at 600/1200, so a per-class .4 target requires both updates");
		final Path limit=root.resolve("limit");
		run(config, challenge, clean, limit, 2, "targetfpr=0", "rounds=2");
		check(status(limit).equals("round_limit") && !Files.exists(limit.resolve("FROZEN")), "An above-target final evaluation cannot be frozen");
		check(Files.readAllLines(limit.resolve("summary.tsv")).get(1).split("\t")[3].equals("600"),
			"The round limit retains the last evaluated package rather than applying an untested final update");
		final Path shortcut=root.resolve("shortcut.config");
		write(shortcut, recipe.replace("ddlreject=8", "ddlreject=0").replace("ddlaccept=1025", "ddlaccept=8"));
		final Path floor=root.resolve("shortcut_floor");
		run(shortcut, challenge, clean, floor, 2, "targetfpr=0", "rounds=3");
		check(status(floor).equals("ddl_floor_exceeds_target") && !Files.exists(floor.resolve("FROZEN")) &&
			!Files.exists(floor.resolve("references.tsv.gz")), "Unmodifiable shortcut errors must stop honestly without changing aligner gates");
		for(String invalid:new String[]{"targetfpr=NaN", "targetfpr=Infinity", "targetfpr=-.1",
			"rounds=0", "rounds=33", "rounds=4294967297", "pos="+challenge}){
			reject(config, challenge, clean, root.resolve("invalid"+checks), invalid);
		}
		final Path empty=root.resolve("empty.faa"), bad=root.resolve("bad.faa");
		write(empty, "");
		write(bad, ">bad\nACD*EFG\n");
		reject(config, empty, clean, root.resolve("empty_result"), "targetfpr=0");
		reject(config, bad, clean, root.resolve("bad_result"), "targetfpr=0");
		marginRepair(root, recipe, references, challenge, hit);
		classCountRepair(root, config, challenge, clean);
		if(args.length==1){pooledClassRepair(root, config, challenge, clean);}
		for(String invalid:new String[]{"marginscale=0", "marginscale=NaN", "targetaccepts=1",
			"expectedpanelqueries=1200,600,1"}){
			reject(config, challenge, clean, root.resolve("invalid_margin"+checks), invalid);
		}
		System.out.println("FAMILY_CALIBRATE_TEST_PASS checks="+checks);
	}

	/** Different integer ceilings follow named panel order, retain scalar parity,
	 * and stop on observed ceilings rather than inventing an attainable equality. */
	private static void classCountRepair(Path root, Path config, Path challenge, Path clean)throws Exception{
		final Path boundary=root.resolve("class_count_boundary");
		run(config, challenge, clean, boundary, 2, "targetacceptsbyclass=600,0", "perclass=t", "expectedinitialaccepts=600,0");
		check(status(boundary).equals("converged") && Files.readAllLines(boundary.resolve("rounds.tsv")).size()==2 &&
			!Files.exists(boundary.resolve("references.tsv.gz")), "A distinct per-class ceiling admits equality without an update");
		check(Files.readAllLines(boundary.resolve("targets.tsv")).equals(Arrays.asList(
			"panel\tqueries\tmaximum_accepts\tobserved_accepts", "challenge\t1200\t600\t600", "clean\t600\t0\t0")),
			"Saved exact targets must bind each named class to its complete population and observed count");
		check(Files.readAllLines(boundary.resolve("calibration.tsv")).contains("target_scope\tper_class_counts") &&
			Files.readAllLines(boundary.resolve("calibration.tsv")).contains("target_accepts_by_class\t600,0"),
			"Count-vector provenance must not be mislabeled as a shared rate");
		for(int threads:new int[]{1, 4}){
			final Path out=root.resolve("class_count_zero_t"+threads);
			run(config, challenge, clean, out, threads, "targetacceptsbyclass=0,0", "perclass=t", "rounds=3");
			check(status(out).equals("converged"), "Joint zero ceilings must converge after later winners are repaired");
			for(String file:new String[]{"references.tsv.gz", "summary.tsv", "family_false_accepts.tsv.gz", "cohorts.tsv"}){
				check(lines(out.resolve(file)).equals(lines(root.resolve("threads"+threads+"/"+file))),
					"Count-vector repair must preserve legacy zero-rate results and worker parity: "+file);
			}
		}
		final Path swapped=root.resolve("class_count_swapped");
		run(config, challenge, clean, swapped, 2, "targetacceptsbyclass=0,600", "perclass=t", "rounds=3", "marginscale=.985");
		check(status(swapped).equals("converged") && Files.readAllLines(swapped.resolve("rounds.tsv")).size()==4,
			"Swapping the same total budget must repair the challenge class rather than pool counts");
		check(Files.readAllLines(swapped.resolve("targets.tsv")).get(2).equals("clean\t600\t600\t0"),
			"An attainable count below its ceiling is retained without manufacturing equality");
		final Path replay=root.resolve("class_count_replay");
		FamilyEval.main(new String[]{"in="+swapped, "neg="+challenge+","+clean, "negnames=challenge,clean",
			"negsha80=unmeasured,unmeasured", "negcounts=1200,600", "verification=structural", "out="+replay,
			"marginscale=.985", "t=2"});
		check(lines(replay.resolve("summary.tsv")).equals(lines(swapped.resolve("summary.tsv"))),
			"Serialized count-vector repair must reproduce both classes at the declared deployed margin");
		final Path limit=root.resolve("class_count_limit");
		run(config, challenge, clean, limit, 2, "targetacceptsbyclass=0,0", "perclass=t", "rounds=2");
		check(status(limit).equals("round_limit") && !Files.exists(limit.resolve("FROZEN")),
			"Count-vector repair must retain honest nonconvergence at its evaluation limit");
		final Path shortcut=root.resolve("class_count_shortcut");
		run(root.resolve("shortcut.config"), challenge, clean, shortcut, 2, "targetacceptsbyclass=0,600", "perclass=t");
		check(status(shortcut).equals("ddl_floor_exceeds_target") && !Files.exists(shortcut.resolve("FROZEN")),
			"Unmodifiable DDL accepts must use the failing class ceiling, not a pooled allowance");
		reject(config, challenge, clean, root.resolve("counts_without_perclass"), "targetacceptsbyclass=0,0");
		final Path prepared=root.resolve("invalid_count_preparation");
		boolean preparationRejected=false;
		try{
			new FamilyCalibrate(new String[]{"prepare="+root.resolve("absent_preparation.config"), "out="+prepared,
				"targetacceptsbyclass=0,0", "perclass=t"});
		}catch(IllegalArgumentException expected){
			preparationRejected="Per-class acceptance counts replace targetfpr/targetaccepts and require existing-package per-class calibration".equals(expected.getMessage());
		}
		check(preparationRejected && !Files.exists(prepared), "Count-vector mode must reject preparation before opening its manifest");
		for(String[] invalid:new String[][]{
			{"targetacceptsbyclass=0"}, {"targetacceptsbyclass=0,0,0"}, {"targetacceptsbyclass=0,"},
			{"targetacceptsbyclass=-1,0"}, {"targetacceptsbyclass=1201,0"}, {"targetacceptsbyclass=0,601"},
			{"targetacceptsbyclass=9223372036854775808,0"}, {"targetacceptsbyclass=0.5,0"},
			{"targetacceptsbyclass=0,0", "targetfpr=0"}, {"targetacceptsbyclass=0,0", "targetaccepts=0"},
			{"targetacceptsbyclass=0,0", "perclass=f"}, {"targetacceptsbyclass=0,0", "expectedpanelqueries=1200"},
			{"targetacceptsbyclass=0,0", "negnames="}, {"targetacceptsbyclass=0,0", "negnames=challenge,challenge"},
			{"targetacceptsbyclass=0,0", "targetacceptsbyclass=0,0"},
			{"targetacceptsbyclass=0,0", "expectedinitialaccepts=599,0"}
		}){
			final ArrayList<String> args=new ArrayList<String>();
			args.add("perclass=t");
			args.addAll(Arrays.asList(invalid));
			reject(config, challenge, clean, root.resolve("invalid_counts"+checks), args.toArray(new String[args.size()]));
		}
	}

	/** Independent ceilings exercise the actual updater, terminal package and shortcut floor. */
	private static void pooledClassRepair(Path root, Path config, Path challenge, Path clean)throws Exception{
		final Path boundary=root.resolve("pooled_class_boundary");
		run(config, challenge, clean, boundary, 2, "targetacceptsbyclass=600,0", "targetacceptspooled=600", "perclass=t");
		check(status(boundary).equals("converged") && Files.readAllLines(boundary.resolve("rounds.tsv")).size()==2 &&
			!Files.exists(boundary.resolve("references.tsv.gz")), "Both exact ceilings admit equality without rewriting references");
		check(lines(boundary.resolve("summary.tsv")).equals(lines(root.resolve("class_count_boundary/summary.tsv"))) &&
			Arrays.equals(Files.readAllBytes(boundary.resolve("classifier.config")), Files.readAllBytes(root.resolve("class_count_boundary/classifier.config"))),
			"A nonbinding pooled ceiling preserves the legacy class-count model and observations");
		check(Files.readAllLines(boundary.resolve("calibration.tsv")).contains("target_scope\tpooled_and_per_class_counts"),
			"Combined target scope must be explicit");
		check(Files.readAllLines(boundary.resolve("gate_change_summary.tsv")).contains("unique_changed_families\t0"),
			"An already-passing pooled/class parent reports zero unique changed families");
		final String[] pooled=Files.readAllLines(boundary.resolve("pooled_target.tsv")).get(1).split("\t");
		check(pooled.length==5 && pooled[0].equals("1800") && pooled[1].equals("600") && pooled[2].equals("600") &&
			Double.parseDouble(pooled[3])==600.0/1800 && Double.parseDouble(pooled[4])==600.0/1800,
			"Pooled output records the ceiling, achieved count and both rates against the complete population");
		final String[] counts=Files.readAllLines(boundary.resolve("targets.tsv")).get(1).split("\t");
		check(counts.length==6 && counts[3].equals("600") && Double.parseDouble(counts[4])==.5 && Double.parseDouble(counts[5])==.5,
			"Combined mode reports class ceilings and achieved rates using each class population");
		for(int threads:new int[]{1, 4}){
			final Path out=root.resolve("pooled_only_failure_t"+threads);
			run(config, challenge, clean, out, threads, "targetacceptsbyclass=600,0", "targetacceptspooled=599", "perclass=t", "rounds=3");
			check(status(out).equals("converged") && Files.readAllLines(out.resolve("rounds.tsv")).size()==4 &&
				Files.exists(out.resolve("references.tsv.gz")), "Pooled-only failure must trigger both later-winner updates even while classes pass");
			final List<String> changes=Files.readAllLines(out.resolve("gate_change_summary.tsv"));
			check(changes.contains("unique_changed_families\t2") && changes.contains("identity_changed_families\t2") &&
				changes.contains("blosum_changed_families\t0") && changes.contains("identity_stored_delta_sum\t4.0") &&
				changes.contains("identity_stored_delta_max\t2.0") && changes.contains("identity_effective_delta_sum\t4.0") &&
				changes.contains("member_lowers_unchanged\ttrue"),
				"Final gate summary counts the unique changed roster and actual stored/effective deltas");
			for(String file:new String[]{"references.tsv.gz", "summary.tsv", "family_false_accepts.tsv.gz", "cohorts.tsv"}){
				check(lines(out.resolve(file)).equals(lines(root.resolve("threads"+threads+"/"+file))),
					"Pooled-triggered updates retain legacy updater and thread parity: "+file);
			}
			final Path replay=root.resolve("pooled_only_replay_t"+threads);
			FamilyEval.main(new String[]{"in="+out, "neg="+challenge+","+clean, "negnames=challenge,clean",
				"negsha80=unmeasured,unmeasured", "negcounts=1200,600", "verification=structural", "out="+replay, "t="+threads});
			check(lines(replay.resolve("summary.tsv")).equals(lines(out.resolve("summary.tsv"))),
				"Combined-mode serialized gates reproduce both class counts after reload");
		}
		check(lines(root.resolve("pooled_only_failure_t1/gate_change_summary.tsv")).equals(
			lines(root.resolve("pooled_only_failure_t4/gate_change_summary.tsv"))), "Unique-family/gate reporting preserves worker parity");
		final Path classOnly=root.resolve("class_only_failure");
		run(config, challenge, clean, classOnly, 2, "targetacceptsbyclass=599,0", "targetacceptspooled=600", "perclass=t", "rounds=3");
		check(status(classOnly).equals("converged") && Files.readAllLines(classOnly.resolve("rounds.tsv")).size()==4,
			"Class-only failure still triggers updates when the pooled ceiling passes");
		final Path limit=root.resolve("pooled_class_round_limit");
		run(config, challenge, clean, limit, 2, "targetacceptsbyclass=600,0", "targetacceptspooled=599", "perclass=t", "rounds=2");
		check(status(limit).equals("round_limit") && !Files.exists(limit.resolve("FROZEN")),
			"Pooled-only failure cannot freeze an above-target last evaluated package");
		for(int mode=0; mode<3; mode++){
			final Path out=root.resolve("pooled_class_shortcut"+mode);
			run(root.resolve("shortcut.config"), challenge, clean, out, 2, "targetacceptsbyclass="+(mode==1 ? "599,0" : "600,0"),
				"targetacceptspooled="+(mode==0 ? "599" : "600"), "perclass=t");
			check(status(out).equals(mode==2 ? "converged" : "ddl_floor_exceeds_target") &&
				Files.exists(out.resolve("FROZEN"))==(mode==2) && !Files.exists(out.resolve("references.tsv.gz")),
				"Either independent shortcut-floor violation must stop; joint equality converges without changes");
		}
		for(String[] invalid:new String[][]{
			{"targetacceptspooled=0"}, {"targetacceptspooled=-1"}, {"targetacceptspooled=0.5"},
			{"targetacceptspooled=9223372036854775808"}, {"targetacceptspooled=1801"},
			{"targetacceptspooled=0", "targetacceptspooled=0"}, {"targetacceptspooled=0", "targetfpr=0"},
			{"targetacceptspooled=0", "targetaccepts=0"}, {"targetacceptspooled=0", "perclass=f"},
			{"targetacceptspooled=0", "expectedpanelqueries=1200"}
		}){
			final ArrayList<String> args=new ArrayList<String>();
			args.add("perclass=t");
			if(invalid.length!=1 || !invalid[0].equals("targetacceptspooled=0")){args.add("targetacceptsbyclass=600,0");}
			args.addAll(Arrays.asList(invalid));
			reject(config, challenge, clean, root.resolve("invalid_pooled"+checks), args.toArray(new String[args.size()]));
		}
		final byte[] parent=Files.readAllBytes(root.resolve("references.tsv"));
		final Path bound=root.resolve("inverse_margin_representation_bound");
		String diagnostic=null;
		try{
			run(config, challenge, clean, bound, 2, "targetacceptsbyclass=0,0", "targetacceptspooled=0", "perclass=t", "marginscale=.6", "rounds=3");
		}catch(IllegalArgumentException expected){diagnostic=expected.getMessage();}
		check(diagnostic!=null && diagnostic.startsWith("Invalid upper/final alignment gates: "),
			"Low positive margins must expose the existing stored-upper representation bound, not a generic failure");
		check(Arrays.equals(parent, Files.readAllBytes(root.resolve("references.tsv"))) && !Files.exists(bound.resolve("FROZEN")) &&
			!Files.exists(bound.resolve("PASS")), "An unrepresentable inverse-margin update preserves its parent without claiming convergence");
		write(root.resolve("inverse_margin_representation_bound.tsv"), "metric\tvalue\nmarginscale\t0.6\nexpected_diagnostic\t"+diagnostic+"\nparent_unchanged\ttrue\n");
	}

	/** Actual driver/export/reload checks catch a repair that meets the target only
	 * at margin1. Two families expose a different winner after the first update. */
	private static void marginRepair(Path root, String recipe, Path original, Path challenge, String hit)throws Exception{
		for(boolean blosum:new boolean[]{false, true}){
			final String route=blosum ? "blosum" : "identity";
			final Path refs=root.resolve("margin_"+route+".tsv");
			if(!blosum){Files.copy(original, refs);}else{
				final ByteBuilder rows=new ByteBuilder();
				FamilyReferenceCodec.appendHeader(new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA, 5, 1024, 0, 3, 2, 1,
					FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING), rows);
				FamilyReferenceCodec.appendFiveLines(family("f2", 2, hit, true), rows);
				FamilyReferenceCodec.appendFiveLines(FamilyReference.empty("empty", 1), rows);
				FamilyReferenceCodec.appendFiveLines(family("f0", 0, hit, true), rows);
				write(refs, rows.toString());
			}
			final Path config=root.resolve("margin_"+route+".config");
			write(config, recipe.replace("familyreferences=references.tsv", "familyreferences="+refs)
				.replaceAll("familyreferencessha80=[^\\n]+", "familyreferencessha80=unmeasured"));
			for(int threads:new int[]{1, 4}){
				final Path out=root.resolve("margin_"+route+"_t"+threads);
				FamilyCalibrate.main(new String[]{"in="+config, "neg="+challenge, "negnames=reversed", "negsha80=unmeasured",
					"out="+out, "t="+threads, "verification=structural", "expectedpanelqueries=1200",
					"expectedinitialaccepts=600", "marginscale=.985", "targetaccepts=0", "perclass=t", "rounds=3"});
				check(status(out).equals("converged") && Files.exists(out.resolve("FROZEN")), "Scaled repair must freeze after evaluated convergence");
				final List<String> rounds=Files.readAllLines(out.resolve("rounds.tsv"));
				check(rounds.size()==4 && rounds.get(1).split("\t")[2].equals("600") &&
					rounds.get(2).split("\t")[2].equals("600") && rounds.get(3).split("\t")[2].equals("0"),
					"Scaled repair must re-evaluate later candidates, not replay only earlier winners");
				final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(out.resolve("references.tsv.gz").toString(), 1, 10);
				for(FamilyReference f:loaded.live){
					check(f.identityLower==50 && f.blosumLower==.5f, "Margin repair must retain member lower gates");
					if(blosum){
						check(f.identityUpper==101 && prot.FamilyTailCascade.evaluationGate(.5f, f.blosumFinal, 0, .985)>=1.005,
							"Published BLOSUM gate must retain the requested increase at the deployed margin");
					}else{
						check(f.blosumFinal==1.1f && prot.FamilyTailCascade.evaluationGate(50, f.identityUpper, 30, .985)>=101,
							"Published identity gate must retain the requested increase at the deployed margin");
					}
				}
				final Path replay=root.resolve("margin_replay_"+route+"_t"+threads);
				FamilyEval.main(new String[]{"in="+out, "neg="+challenge, "negnames=reversed", "negsha80=unmeasured",
					"out="+replay, "t="+threads, "verification=structural", "negcounts=1200", "marginscale=.985"});
				check(Files.readAllLines(replay.resolve("summary.tsv")).get(1).split("\t")[3].equals("0"),
					"Reloaded frozen package must reproduce the tighter-margin result");
			}
			for(String file:new String[]{"references.tsv.gz", "summary.tsv", "family_false_accepts.tsv.gz", "cohorts.tsv"}){
				check(lines(root.resolve("margin_"+route+"_t1/"+file)).equals(lines(root.resolve("margin_"+route+"_t4/"+file))),
					"Worker count must preserve the margin repair: "+route+" "+file);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fixture utilities      ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs existing-package calibration with two named panels and explicit structural counts.
	 * Extra arguments are appended unchanged; duplicate-policy rejection remains the driver's job.
	 */
	private static void run(Path config, Path challenge, Path clean, Path out, int threads, String... extra) throws Exception{
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+config, "neg="+challenge+","+clean,
			"negnames=challenge,clean", "out="+out, "t="+threads, "verification=structural",
			"negsha80=unmeasured,unmeasured", "expectedpanelqueries=1200,600"));
		args.addAll(Arrays.asList(extra));
		FamilyCalibrate.main(args.toArray(new String[args.size()]));
	}
	/** Reads the explicit status row; absence is a fixture failure, never inferred convergence. */
	private static String status(Path out) throws Exception{
		for(String row:Files.readAllLines(out.resolve("calibration.tsv"))){
			if(row.startsWith("status\t")){return row.substring(7);}
		}
		throw new AssertionError("Missing explicit calibration status");
	}
	/** Returns the digest suffix needed by explicit fixture binding checks. */
	private static String hash(Path p){return DigestSuffix.file(p.toString());}
	/** Requires a policy/state rejection and confirms the failed run did not publish PASS. */
	private static void reject(Path config, Path challenge, Path clean, Path out, String... extra) throws Exception{
		boolean failed=false;
		try{
			run(config, challenge, clean, out, 4, extra);
		}catch(IllegalArgumentException|IllegalStateException expected){failed=true;}
		check(failed && !Files.exists(out.resolve("PASS")), "Malformed calibration input must fail before completion: "+Arrays.toString(extra));
	}
	/** Constructs an identity-limited family from the supplied synthetic consensus. */
	private static FamilyReference family(String id, int rank, String sequence){
		return family(id, rank, sequence, false);
	}
	/** Builds fresh sketch/spectrum state with an identity- or BLOSUM-limited upper gate.
	 * Rank two deliberately receives a different spectrum to avoid a ranking tie.
	 */
	private static FamilyReference family(String id, int rank, String sequence, boolean blosum){
		final byte[] raw=sequence.getBytes(StandardCharsets.US_ASCII);
		final char[] ddl=new char[1024];
		DDLFamilyProfiles.sketch(raw, 5, 0, ddl);
		final ProteinSpectrum.Spectrum counts=new ProteinSpectrum.Spectrum(2);
		counts.add(raw);
		final float[] spectrum=FamilyReferenceExport.frequencies(counts);
		// Both consensuses align perfectly, but f2 ranks lower so the oracle does not depend on tie breaking.
		if(rank==2){
			Arrays.fill(spectrum, 0);
			spectrum[0]=1;
		}
		return new FamilyReference(id, rank, 1, raw.length, raw.length-4, new ProteinSequence(id, raw), ddl,
			spectrum, 0, 50, .5f, 1.1f, blosum ? 101 : 99, blosum ? .99f : 1.1f);
	}
	/** Reads decoded ASCII lines so plain and compressed outputs compare by content. */
	private static List<String> lines(Path path){
		final ArrayList<String> out=new ArrayList<String>();
		final ByteFile input=ByteFile.makeByteFile(path.toString(), false);
		try{
			for(byte[] row=input.nextLine(); row!=null; row=input.nextLine()){
				out.add(new String(row, StandardCharsets.US_ASCII));
			}
		}finally{check(!input.close(), "Fixture output must be fully readable");}
		return out;
	}
	/** Returns the requested number of copies of a literal sequence fragment. */
	private static String repeat(String text, int times){
		final StringBuilder out=new StringBuilder();
		for(int i=0; i<times; i++){out.append(text);}
		return out.toString();
	}
	/** Appends one two-line FASTA record to a caller-owned builder. */
	private static void fasta(StringBuilder out, String id, String text){
		out.append('>').append(id).append('\n').append(text).append('\n');
	}
	/** Writes fixture ASCII bytes, replacing an existing file when the caller requests it. */
	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}
	/** Fails loudly on a broken oracle and counts successful checks in this JVM. */
	private static void check(boolean valid, String why){
		if(!valid){throw new AssertionError(why);}
		checks++;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Successful assertion count; the standalone fixture is not reentrant. */
	private static int checks;
}
