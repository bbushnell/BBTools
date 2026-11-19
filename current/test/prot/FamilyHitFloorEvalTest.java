package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import prot.DigestSuffix;
import prot.FamilyEval;
import prot.FamilyHitFloorBinding;

/** Saved-floor evaluation must use the same literal classifier decisions without
 * fitting again or emitting a parent-only recipe as the complete classifier.
 * @author Brian, Keqing
 */
public final class FamilyHitFloorEvalTest{

	/** Checks ordinary and fitted evaluation, worker parity and invalid bindings in a fresh directory. */
	public static void main(final String[] args)throws Exception{
		check(args.length==1, "Supply a fresh Dori fixture directory");
		final Path root=Paths.get(args[0]).toAbsolutePath();
		check(!Files.exists(root), "Saved-floor fixture output must be fresh");
		final Inputs inputs=writeDefaultHashInputs(root);
		final Path parent=inputs.parent, fitting=inputs.fitting, floors=inputs.floors, compressed=inputs.compressed;
		final Path positive=inputs.positive, negative=inputs.negative;
		final String originalPin=pin(floors);
		final String[] ordinary={"in="+parent, "pos="+positive, "possha80="+pin(positive), "neg="+negative,
			"negsha80="+pin(negative), "negnames=k5", "marginscale=1", "threadsin=1"};
		for(int threads:new int[]{1, 2}){
			final Path output=root.resolve("ordinary"+threads);
			invoke(ordinary, output, "t="+threads);
			checkCounts(output, 0);
			check(Files.exists(output.resolve("classifier.config")) && !Files.exists(output.resolve("parent_classifier.config")),
				"Ordinary evaluation retains its existing output contract");
		}
		for(String file:new String[]{"summary.tsv", "classifier.config", "evaluation_policy.tsv", "inputs.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("ordinary1/"+file)),
				Files.readAllBytes(root.resolve("ordinary2/"+file))),
				"Ordinary semantic reports remain independent of worker count: "+file);
		}
		for(Path table:new Path[]{floors, compressed}){
			final String[] bound=bound(ordinary, parent, fitting, table);
			for(int threads:new int[]{1, 2}){
				final Path output=root.resolve(table.getFileName()+"-threads"+threads);
				invoke(bound, output, "t="+threads);
				checkCounts(output, 2);
				check(!Files.exists(output.resolve("classifier.config")) && Files.exists(output.resolve("parent_classifier.config")),
					"The parent recipe must not masquerade as the full fitted classifier");
				final String policy=new String(Files.readAllBytes(output.resolve("evaluation_policy.tsv")), StandardCharsets.UTF_8);
				check(policy.contains("admission\tfamily_fitted\n") && policy.contains("raw_hit_admission\tddl_ge_pinned_family_floor\n") &&
					policy.contains("floors_sha80\t"+pin(table)+"\n") && policy.contains("initial_support\tfrozen_from_fitting;never_reestimated_from_evaluation\n"),
					"Evaluation reports the actual saved vector and never claims a refitted support mask");
				check(Files.exists(output.resolve("admission_summary.tsv")), "Actual candidate/DP work is retained for fitted evaluation");
			}
			final String[] structural=structural(bound, root.resolve(table.getFileName()+"-structural"));
			invoke(structural, root.resolve(table.getFileName()+"-structural"), "t=1");
			checkCounts(root.resolve(table.getFileName()+"-structural"), 2);
			reject(replace(structural, root.resolve("wrong_parent"), "insha80="+pin(fitting)),
				root.resolve("wrong_parent"), "t=1");
			final String[] changed=structural(bound, root.resolve("changed_floor"));
			final FamilyEval eval=new FamilyEval(append(changed, "t=1"));
			Files.setLastModifiedTime(table,
				java.nio.file.attribute.FileTime.fromMillis(Files.getLastModifiedTime(table).toMillis()+5000));
			boolean rejected=false;
			try{
				eval.process();
			}catch(IllegalArgumentException | IllegalStateException expected){rejected=true;}
			check(rejected, "Structural saved-floor evaluation must retain the original floor-table snapshot");
			if(table.equals(floors)){
				for(String[] mutation:new String[][]{
						{"unsupported_floor5", "0", "4", "5"},
						{"initial_zero_nonzero", "0", "5", "1"},
						{"pooled_overflow", "0", "6", "1"}}){
					final Path bad=mutate(root.resolve(mutation[0]+".tsv"), floors, 6+Integer.parseInt(mutation[1]),
						Integer.parseInt(mutation[2]), mutation[3]);
					reject(replace(structural, root.resolve(mutation[0]), "floors="+bad),
						root.resolve(mutation[0]), "t=1");
				}
				reject(replace(structural, root.resolve("denominator_zero"), "floors="+
					mutate(root.resolve("denominator_zero.tsv"), floors, 4, 1, "0")),
					root.resolve("denominator_zero"), "t=1");
				reject(replace(structural, root.resolve("roster_mismatch"), "floors="+
					mutate(root.resolve("roster_mismatch.tsv"), floors, 6, 1, "other")),
					root.resolve("roster_mismatch"), "t=1");
			}
			for(String extra:new String[]{"top=49", "ddlreject=6", "admission=family_mean", "heuristic=hitfourth_h020_c150",
				"paired=t", "benchmarkrepeats=2", "minqueryaa=39", "marginscale=2", "insha80=00000000000000000000"}){
				reject(bound, root.resolve("rejected"+(rejections++)), extra);
			}
			for(String omitted:new String[]{"insha80=", "fitting=", "fittingsha80=", "floorssha80=", "marginscale=", "possha80="}){
				final ArrayList<String> incomplete=new ArrayList<String>();
				for(String arg:bound){if(!arg.startsWith(omitted)){incomplete.add(arg);}}
				reject(incomplete.toArray(new String[0]), root.resolve("rejected"+(rejections++)), "t=1");
			}
		}
		reject(ordinary, root.resolve("bare_floor"), "floors");
		check(originalPin.equals(pin(floors)), "Evaluation must never rewrite selected floors");
		System.out.println("FAMILY_HIT_FLOOR_EVAL_TEST_PASS checks="+checks+" scientific_runs=0");
	}

	/** Creates shared synthetic parent, frozen-floor and evaluation inputs without running evaluation. */
	static Inputs writeDefaultHashInputs(final Path root)throws Exception{
		Files.createDirectory(root);
		final byte[] query=FamilyHitFloorIntegrationTest.sequence();
		FamilyHitFloorIntegrationTest.fixture(root.resolve("parent"), query);
		final Path parent=root.resolve("parent/classifier.config");
		// The evaluator binds this provenance file; it must not reopen fitting sequences.
		final Path fitting=write(root.resolve("fitting.tsv"), "synthetic frozen provenance; no readable sequence paths\n");
		final FamilyHitFloorBinding binding=new FamilyHitFloorBinding(parent, pin(parent), 1, fitting, pin(fitting), 1);
		final long[] initial=new long[21];
		initial[12]=1;
		final Path floors=root.resolve("floors.tsv");
		binding.write(floors, binding.probe().withInitialSupport(1, initial), 1, initial);
		final Path compressed=root.resolve("floors.tsv.gz");
		try(java.util.zip.GZIPOutputStream out=new java.util.zip.GZIPOutputStream(Files.newOutputStream(compressed))){
			out.write(Files.readAllBytes(floors));
		}
		final Path positive=write(root.resolve("evaluation.faa"), ">one\n"+new String(query, StandardCharsets.US_ASCII)+
			"\n>two\n"+new String(query, StandardCharsets.US_ASCII)+"\n");
		final char[] unknown=new char[64];
		Arrays.fill(unknown, 'X');
		final Path negative=write(root.resolve("negative.faa"), ">unknown\n"+new String(unknown)+"\n");
		return new Inputs(parent, fitting, floors, compressed, positive, negative);
	}

	/** Paths owned by the fixture root and shared with the input-only producer. */
	static final class Inputs{
		Inputs(Path parent, Path fitting, Path floors, Path compressed, Path positive, Path negative){
			this.parent=parent;
			this.fitting=fitting;
			this.floors=floors;
			this.compressed=compressed;
			this.positive=positive;
			this.negative=negative;
		}
		final Path parent, fitting, floors, compressed, positive, negative;
	}

	/** Adds the explicit parent, fitting and floor identities to the ordinary evaluator arguments. */
	private static String[] bound(String[] ordinary, Path parent, Path fitting, Path floors){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(ordinary));
		args.add("insha80="+pin(parent));
		args.add("fitting="+fitting);
		args.add("fittingsha80="+pin(fitting));
		args.add("floors="+floors);
		args.add("floorssha80="+pin(floors));
		return args.toArray(new String[0]);
	}

	/** Switches evaluation inputs to structural verification while retaining parent/fitting bindings. */
	private static String[] structural(String[] base, Path output){
		final ArrayList<String> args=new ArrayList<String>();
		for(String arg:base){
			if(arg.startsWith("floorssha80=")){args.add("floorssha80=unmeasured");}
			else if(arg.startsWith("possha80=")){args.add("possha80=unmeasured");}
			else if(arg.startsWith("negsha80=")){args.add("negsha80=unmeasured");}
			else{args.add(arg);}
		}
		args.add("verification=structural");
		args.add("out="+output);
		return args.toArray(new String[0]);
	}

	/** Copies arguments while replacing existing output and named-option values. */
	private static String[] replace(String[] base, Path output, String replacement){
		final ArrayList<String> args=new ArrayList<String>();
		final String key=replacement.substring(0, replacement.indexOf('=')+1);
		for(String arg:base){args.add(arg.startsWith("out=") ? "out="+output : arg.startsWith(key) ? replacement : arg);}
		return args.toArray(new String[0]);
	}

	/** Returns a new argument array with the supplied suffix in order. */
	private static String[] append(String[] base, String... extra){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(base));
		args.addAll(Arrays.asList(extra));
		return args.toArray(new String[0]);
	}

	/** Copies a floor table and replaces one bounded field to create an independent malformed input. */
	static Path mutate(Path output, Path floors, int rowIndex, int column, String value)throws Exception{
		final List<String> lines=Files.readAllLines(floors, StandardCharsets.US_ASCII);
		final int row=rowIndex;
		check(row>=0 && row<lines.size(), "Floor mutation must target a floor-table line");
		final String[] fields=lines.get(row).split("\t", -1);
		check(column>=0 && column<fields.length, "Floor mutation must target an existing floor-table field");
		fields[column]=value;
		lines.set(row, String.join("\t", fields));
		return Files.write(output, lines, StandardCharsets.US_ASCII, java.nio.file.StandardOpenOption.CREATE_NEW);
	}

	/** Adds an output argument when absent, then invokes evaluation with the requested extra option. */
	private static void invoke(String[] base, Path output, String extra)throws Exception{
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(base));
		boolean hasOutput=false;
		for(String arg:base){
			if(arg.startsWith("out=")){
				hasOutput=true;
				break;
			}
		}
		if(!hasOutput){args.add("out="+output);}
		args.add(extra);
		FamilyEval.main(args.toArray(new String[0]));
	}

	/** Checks the two complete synthetic panels against literal positive and negative counts. */
	private static void checkCounts(Path output, int accepted)throws Exception{
		final List<String> rows=Files.readAllLines(output.resolve("summary.tsv"), StandardCharsets.US_ASCII);
		check(rows.size()==3, "Both complete synthetic panels must have one report row");
		final String[] positive=rows.get(1).split("\t"), negative=rows.get(2).split("\t");
		check(positive[2].equals("2") && positive[3].equals(Integer.toString(accepted)) &&
			negative[2].equals("1") && negative[3].equals("0"),
			"Literal positive/negative counts must match the independently specified outcome");
	}

	/** Requires invalid bindings to fail before an evaluation output directory is created. */
	private static void reject(String[] base, Path output, String extra)throws Exception{
		boolean rejected=false;
		try{
			invoke(base, output, extra);
		}catch(IllegalArgumentException | IllegalStateException expected){rejected=true;}
		check(rejected && !Files.exists(output), "Invalid binding/override must fail without completed evaluation output: "+extra);
	}

	/** Writes a fresh fixture file; an existing path is an error. */
	private static Path write(Path path, String text)throws Exception{
		return Files.write(path, text.getBytes(StandardCharsets.US_ASCII), java.nio.file.StandardOpenOption.CREATE_NEW);
	}
	/** Supplies the identity required by the fixture's explicit hash-verification cases. */
	private static String pin(Path path){return DigestSuffix.file(path.toString());}
	/** Counts the checked expectation and fails with its diagnostic explanation. */
	private static void check(boolean valid, String why){
		checks++;
		if(!valid){throw new AssertionError(why);}
	}
	private static int checks, rejections;
}
