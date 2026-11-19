package test.benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import benchmark.CallGenesSavedMetadata;

import shared.Shared;

/** Literal saved-table ABI, order and rejection tests; sequence paths deliberately do not exist.
 * @author Sayu
 */
public final class CallGenesSavedMetadataTest{
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh absolute fixture directory");
		final Path root=Paths.get(args[0]);
		check(root.isAbsolute() && !Files.exists(root), "Fresh absolute fixture directory required");
		Files.createDirectory(root);
		Shared.setThreads(1);
		for(int n:new int[]{
			1, 3
		}){
			final Path input=dataset(root.resolve("n"+n), n);
			final Path out=root.resolve("out"+n);
			run(input, out, "cached", n);
			String expectedInputs="policy\tgenome\tfasta\treference\tcalls\tgrades\tgrade_arm\n";
			for(int g=0; g<n; g++){
				expectedInputs+="cached_native\tg"+g+"\t/missing/f"+g+"\t/missing/r"+g+"\t/missing/n"+g+"\t"+out+"/grades.tsv\tnative\n";
				expectedInputs+="cached_formula\tg"+g+"\t/missing/f"+g+"\t/missing/r"+g+"\t/missing/b"+g+"\t"+out+"/grades.tsv\tbonus400\n";
			}
			check(read(out.resolve("inputs.tsv")).equals(expectedInputs), "Accepted order, source fields, aliases and relocated grades path");
			String expectedGrades="genome\tarm\ttp_start_reference\ttp_stop_reference\tfn_start\tfn_stop\ttp_start_query\ttp_stop_query\tfp_start\tfp_stop\n";
			for(int g=n-1; g>=0; g--){
				expectedGrades+="g"+g+"\tbonus400\t1\t2\t3\t4\t5\t6\t7\t8\n";
			}
			for(int g=n-1; g>=0; g--){
				expectedGrades+="g"+g+"\tnative\t1\t2\t3\t4\t5\t6\t7\t8\n";
			}
			check(read(out.resolve("grades.tsv")).equals(expectedGrades), "Saved raw ordering and exact eight grade counters");
			check(read(out.resolve("expected.tsv")).equals("policy\tgenomes\treference_keys\tcorrect\tmissed\textras\tstop_errors\tstart_errors\n"+
				"cached_formula\t"+n+"\t10\t11\t12\t13\t14\t15\n"+"cached_native\t"+n+"\t10\t11\t12\t13\t14\t15\n"), "Saved summary order and six literal controls");
			check(Files.isRegularFile(out.resolve("PASS")), "Complete batch publishes PASS");
			reject(()->run(input, out, "cached", n), out, "Existing output root");
		}
		bad(root, "header", "accepted.tsv", s->s.replace("genome\tfasta", "wrong\tfasta"));
		bad(root, "duplicate", "accepted.tsv", s->s+s.substring(s.indexOf('\n')+1));
		bad(root, "foreign", "retained.tsv", s->s.replace("g0\t", "other\t"));
		bad(root, "relative", "retained.tsv", s->s.replace("/missing/n0", "relative.gff"));
		bad(root, "missing-arm", "raw.tsv", s->s.substring(0, s.lastIndexOf("native\tg0")));
		bad(root, "bad-grade", "raw.tsv", s->s.replace("\t1\t2\t3\t4\t5\t6\t7\t8\n", "\t-1\t2\t3\t4\t5\t6\t7\t8\n"));
		bad(root, "overflow", "summary.tsv", s->s.replace("\t10\t11", "\t9223372036854775808\t11"));
		final Path unicode=dataset(root.resolve("unicode"), 1);
		final Path unicodeOut=root.resolve("unicode-out");
		for(String file:new String[]{
				"accepted.tsv", "retained.tsv", "raw.tsv", "summary.tsv"
		}){
			write(unicode.resolve(file), read(unicode.resolve(file)).replace("g0\t", "g-\u00e9\t"));
		}
		run(unicode, unicodeOut, "cached", 1);
		check(read(unicodeOut.resolve("grades.tsv")).contains("g-\u00e9\tbonus400\t"), "UTF-8 genome bytes survive metadata extraction");
		check(read(unicodeOut.resolve("inputs.tsv")).contains("cached_native\tg-\u00e9\t"), "UTF-8 identifiers survive input-table writing");
		final Path spelling=dataset(root.resolve("spelling"), 1);
		final Path spellingOut=root.resolve("spelling-out");
		write(spelling.resolve("raw.tsv"), read(spelling.resolve("raw.tsv")).replace("\t1\t2\t3\t", "\t0001\t02\t3\t"));
		write(spelling.resolve("summary.tsv"), read(spelling.resolve("summary.tsv")).replace("\t10\t11\t", "\t0010\t11\t"));
		run(spelling, spellingOut, "cached", 1);
		check(read(spellingOut.resolve("grades.tsv")).contains("g0\tbonus400\t0001\t02\t3\t"), "Grade byte spans retain decimal spelling");
		check(read(spellingOut.resolve("expected.tsv")).contains("cached_formula\t1\t0010\t11\t"), "Summary byte spans retain decimal spelling");
		final Path input=dataset(root.resolve("labels"), 1);
		reject(()->run(input, root.resolve("unsafe"), "bad/name", 1), root.resolve("unsafe"), "Unsafe policy");
		reject(()->run(input, root.resolve("wrong-n"), "cached", 2), root.resolve("wrong-n"), "Wrong expected cohort");
		System.out.println("SAVED_CALLER_METADATA_TEST_PASS checks="+checks+" sequence_reads=0");
	}

	/** Own a fresh table fixture with independently ordered policies and absent sequence payloads. */
	private static Path dataset(Path root, int n)throws Exception{
		Files.createDirectory(root);
		String accepted="genome\tfasta\tfasta_sha80\treference\treference_sha80\tnative\tnative_sha80\tbaseline\tbaseline_sha80\tannotations\tannotations_sha80\tgrades\tgrades_sha80\n";
		String retained="genome\tnative_gff\tformula_gff\n";
		String raw=SCORE;
		String summary=SCORE;
		for(int g=0; g<n; g++){
			accepted+="g"+g+"\t/missing/f"+g+"\tunmeasured\t/missing/r"+g+"\tunmeasured\tunused\tunmeasured\tunused\tunmeasured\tunused\tunmeasured\tunused\tunmeasured\n";
			retained+="g"+g+"\t/missing/n"+g+"\t/missing/b"+g+"\n";
		}
		// Reverse both policy and genome order independently from the accepted roster.
		for(String policy:new String[]{
				"formula", "native"
		}){
			for(int g=n-1; g>=0; g--){
				raw+=policy+"\tg"+g+"\tALL\t10\t11\t12\t13\t14\t15\t16\t17\t1\t2\t3\t4\t5\t6\t7\t8\n";
			}
			summary+=policy+"\tALL\tALL\t10\t11\t12\t13\t14\t15\t16\t17\t1\t2\t3\t4\t5\t6\t7\t8\n";
		}
		raw+="other\tg0\tALL\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		write(root.resolve("accepted.tsv"), accepted);
		write(root.resolve("retained.tsv"), retained);
		write(root.resolve("raw.tsv"), raw);
		write(root.resolve("summary.tsv"), summary);
		return root;
	}

	private interface Change{
		String apply(String text);
	}

	private static void bad(Path root, String name, String file, Change change)throws Exception{
		final Path input=dataset(root.resolve(name), 1);
		final Path out=root.resolve(name+"-out");
		write(input.resolve(file), change.apply(read(input.resolve(file))));
		reject(()->run(input, out, "cached", 1), out, name);
	}

	private static void run(Path input, Path out, String policy, int n)throws Exception{
		CallGenesSavedMetadata.extract(input.resolve("accepted.tsv"), input.resolve("raw.tsv"), input.resolve("retained.tsv"), input.resolve("summary.tsv"), out, policy, n);
	}

	private interface Action{
		void run()throws Exception;
	}

	/** An invalid fresh run must not publish PASS; an already completed root remains untouched. */
	private static void reject(Action action, Path out, String why)throws Exception{
		final boolean already=Files.exists(out.resolve("PASS"));
		try{
			action.run();
		}catch(IllegalArgumentException | IllegalStateException e){
			check(already || !Files.exists(out.resolve("PASS")), "Failure published PASS: "+why);
			return;
		}
		throw new AssertionError("Expected rejection: "+why);
	}

	private static void check(boolean valid, String why){
		if(!valid){
			throw new AssertionError(why);
		}
		checks++;
	}

	private static String read(Path path)throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);
	}

	private static void write(Path path, String text)throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8));
	}

	private static int checks;
	private static final String SCORE="policy\tgenome\tannotation_or_called_length_band\treference_keys\tcorrect\tmissed\textras\tstop_errors\tstart_errors\tunique_matched_keys\tduplicate_matching_occurrences\ttp_start_reference\ttp_stop_reference\tfn_start\tfn_stop\ttp_start_query\ttp_stop_query\tfp_start\tfp_stop\n";
}
