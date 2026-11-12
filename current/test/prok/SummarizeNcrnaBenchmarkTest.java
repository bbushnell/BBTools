package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaBenchmark;
import shared.Shared;
import structures.ByteBuilder;

/** Finite native-table oracles: median-of-summed-CPU, nonnested equal-size TP sets,
 * last/first detections, repeat/worker instability, masks, duplicate queries,
 * zero-query inputs and UTF-8 identities. Runs no caller or aligner.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaBenchmarkTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, SummarizeNcrnaBenchmarkTest.class, false);
		final Parser parser=new Parser();
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
			if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		require(parser.out1!=null, "Supply out=fresh_fixture_directory");
		final File root=new File(parser.out1);
		require(!root.exists(), "Fixture output must be fresh");
		Files.createDirectory(root.toPath());
		try{test(new Fixture(root));}
		finally{Shared.closeStream(pp.outstream);}
		System.err.println("NCRNA_BENCHMARK_TEST_PASS");
	}

	private static void test(final Fixture fixture){
		fixture.write();
		SummarizeNcrnaBenchmark.main(fixture.args("valid", fixture.times, fixture.raw, fixture.summary, fixture.manifest, fixture.profile));
		final ArrayList<String> timing=lines(fixture.path("valid.timing.tsv"));
		// ByteBuilder deliberately emits exact integers without trailing decimal zeroes.
		// Literal numeric oracles retain the statistical checks independently of spelling.
		numericRow(timing, "current\t1\t", 2, 3, 20, 10, 30, 100, 10, 100, 110, 30, 130, 9, 8, 10, 0, 0);
		numericRow(timing, TRIAL+"\t1\t", 2, 3, 20, 10, 40, 1, 1, 1, 21, 11, 41, 8, 7, 9, 80.909091, 11.111111);
		numericRow(timing, TRIAL+"\t16\t", 2, 3, 50, 45, 55, 1, 1, 1, 51, 46, 56, 5, 4, 6, -2, -25);
		require(has(lines(fixture.path("valid.log")), "NCRNA_BENCHMARK_SUMMARY_COMPLETE: differing repeat comparisons=1; differing worker comparisons=1"),
			"Completion must expose differing repeat and worker comparisons, counted once per cell comparison");
		final ArrayList<String> comparisons=lines(fixture.path("valid.comparisons.tsv"));
		require(has(comparisons, "POLICY_VS_BASELINE\tcurrent\t1\t1\t"+TRIAL+"\t1\t1\t"+RNA+"\t2\t2\t2\t2\t1\t1"),
			"Equal TP counts must retain two lost/two gained keys and first/last input changes");
		require(has(comparisons, "REPEAT_VS_FIRST\tcurrent\t16\t1\tcurrent\t16\t3\t"+RNA+"\t2\t2\t1\t1\t1\t1"),
			"Changed repeated detection sets must be reported independently of their equal counts");
		require(has(comparisons, "WORKERS_VS_FIRST\tcurrent\t1\t3\tcurrent\t16\t3\t"+RNA+"\t2\t2\t1\t1\t1\t1"),
			"Changed worker-count detection sets must be reported");
		require(has(lines(fixture.path("valid.changes.tsv")), "\tLOST\t"+GENOME+"\t"+RNA+"\tref_β"),
			"Exact reference keys must preserve original UTF-8 bytes");
		numericRow(lines(fixture.path("valid.cells.tsv")), "current\t1\t1\t", 3, 10, 100, 110, 8, 5, 10, 15, 5);
		numericRow(lines(fixture.path("valid.profile.tsv")), "COVERED_RNA_TOTAL\tCOVERED_NOT_EXHAUSTIVE_ATTRIBUTION\t", 2,
			8000000000L, 8, 12, 66.666667);
		final String globalOnly=fixture.path("global-profile.tsv");
		fixture.writeProfile(globalOnly, false);
		SummarizeNcrnaBenchmark.main(fixture.args("global", fixture.times, fixture.raw, fixture.summary, fixture.manifest, globalOnly));
		numericRow(lines(fixture.path("global.profile.tsv")), "COVERED_RNA_TOTAL\tCOVERED_NOT_EXHAUSTIVE_ATTRIBUTION\t", 2,
			5000000000L, 5, 12, 41.666667);
		testTaxonomy(fixture);
		NcrnaSubtypeComparisonTest.run(new File(fixture.root, "subtype-comparison"));

		final String missing=fixture.path("missing-times.tsv");
		writeLines(missing, withoutLast(lines(fixture.times)));
		reject(fixture.args("missing", missing, fixture.raw, fixture.summary, fixture.manifest, fixture.profile), "missing");
		final String duplicate=fixture.path("duplicate-times.tsv");
		final ArrayList<String> duplicated=lines(fixture.times);
		duplicated.add(duplicated.get(1));
		writeLines(duplicate, duplicated);
		reject(fixture.args("duplicate", duplicate, fixture.raw, fixture.summary, fixture.manifest, fixture.profile), "exactly once");
		final String badSummary=fixture.path("bad-summary.tsv");
		final ArrayList<String> changed=lines(fixture.summary);
		changed.set(1, changed.get(1).replace("\t2\t1\t2\t4\t", "\t3\t1\t1\t4\t")
			.replace("\t0.666667\t0.500000\t", "\t0.750000\t0.750000\t"));
		writeLines(badSummary, changed);
		reject(fixture.args("summary", fixture.times, fixture.raw, badSummary, fixture.manifest, fixture.profile), "counts disagree");
		final String badMatch=fixture.path("bad-match.tsv");
		final ArrayList<String> wrongMatch=lines(fixture.raw);
		for(int i=1; i<wrongMatch.size(); i++){
			if(wrongMatch.get(i).contains("\tQUERY\t"+RNA) && wrongMatch.get(i).contains("\tref_β\t")){
				wrongMatch.set(i, wrongMatch.get(i).replace("\tref_β\t", "\tunknown_parent\t"));
				break;
			}
		}
		writeLines(badMatch, wrongMatch);
		reject(fixture.args("match", fixture.times, badMatch, fixture.summary, fixture.manifest, fixture.profile), "same one-to-one");
		final String missingInput=fixture.path("missing-input.tsv");
		final ArrayList<String> incomplete=lines(fixture.manifest);
		for(int i=1; i<incomplete.size(); i++){
			if(incomplete.get(i).startsWith("current\t1\t1\tempty\t")){
				incomplete.remove(i);
				break;
			}
		}
		writeLines(missingInput, incomplete);
		reject(fixture.args("input", fixture.times, fixture.raw, fixture.summary, missingInput, fixture.profile), "each declared input");
		final String[] timingArgs={"in="+fixture.times, "out="+fixture.path("timing-only"), "timingonly=t",
			"policies=current,"+TRIAL, "baseline=current", "workers=1,16", "repeats=3", "outstream="+fixture.path("timing-only.log")};
		SummarizeNcrnaBenchmark.main(timingArgs);
		final ArrayList<String> timingOnly=lines(fixture.path("timing-only.timing.tsv"));
		require(timingOnly.size()==timing.size(), "Timing-only qualification must retain the declared report population");
		for(int i=0; i<timing.size(); i++){
			final String qualification=i==0 ? "qualification" : "Caller outputs were not compared";
			require(timingOnly.get(i).equals(timing.get(i)+"\t"+qualification),
				"Every timing-only row must qualify unchanged medians, ranges and savings inside the table");
		}
		for(final String suffix:new String[]{"grades.tsv", "comparisons.tsv", "changes.tsv"}){
			require(!new File(fixture.path("timing-only."+suffix)).exists(), "Timing-only mode must not manufacture grades or stability evidence");
		}
		final ArrayList<String> timingLog=lines(fixture.path("timing-only.log"));
		require(has(timingLog, "NCRNA_BENCHMARK_TIMING_ONLY_COMPLETE") && has(timingLog, "Accuracy and detection stability were not evaluated"),
			"Timing-only completion must explicitly withhold accuracy and stability claims");
		require(!has(timingLog, "differing repeat comparisons="), "Unevaluated stability cannot be reported as zero");
		final String[] conflicting=Arrays.copyOf(timingArgs, timingArgs.length+1);
		conflicting[1]="out="+fixture.path("timing-conflict");
		conflicting[timingArgs.length-1]="outstream="+fixture.path("timing-conflict.log");
		conflicting[conflicting.length-1]="raw="+fixture.raw;
		reject(conflicting, "does not accept grade or count inputs");
		final String[] profiled=Arrays.copyOf(timingArgs, timingArgs.length+2);
		profiled[1]="out="+fixture.path("timing-profile");
		profiled[timingArgs.length-1]="outstream="+fixture.path("timing-profile.log");
		profiled[timingArgs.length]="profile="+fixture.profile;
		profiled[timingArgs.length+1]="profiletime="+fixture.profileTime;
		SummarizeNcrnaBenchmark.main(profiled);
		require(lines(fixture.path("timing-profile.profile.tsv")).equals(lines(fixture.path("valid.profile.tsv"))),
			"Timing-only profile attribution must retain the independently checked profile clocks and denominator");
	}

	/** Hand-counted input classes distinguish copy changes from first/last detection. */
	private static void testTaxonomy(final Fixture fixture){
		final String taxonomy=fixture.path("taxonomy.tsv");
		writeTaxonomy(taxonomy, GENOME);
		final String[] args=withOption(fixture.args("lineages", fixture.times, fixture.raw, fixture.summary,
			fixture.manifest, fixture.profile), "taxonomy="+taxonomy);
		SummarizeNcrnaBenchmark.main(args);
		for(final String suffix:new String[]{"cells.tsv", "timing.tsv", "grades.tsv", "comparisons.tsv", "changes.tsv", "profile.tsv"}){
			require(lines(fixture.path("lineages."+suffix)).equals(lines(fixture.path("valid."+suffix))),
				"Optional taxonomy must not change any existing result: "+suffix);
		}
		final String comparison="POLICY_VS_BASELINE\tcurrent\t1\t1\t"+TRIAL+"\t1\t1\t"+RNA+"\t";
		final ArrayList<String> classes=lines(fixture.path("lineages.lineages.tsv"));
		require(has(classes, comparison+"Bacteria\tP_α\trare\tknown\t2\t2\t1\t2\t1\t2\t1\t0\t1"),
			"Two bacterial inputs lose two references and gain one, with one last loss and no first detection");
		require(has(classes, comparison+"Archaea\tP_C\tbreadth\tknown\t1\t0\t1\t0\t1\t0\t1\t1\t0"),
			"The previously undetected archaeal input supplies one gained reference and one first detection");
		require(has(classes, comparison+"unknown\tunknown\tunknown\tunclassified\t2\t0\t0\t0\t0\t0\t0\t0\t0"),
			"Unknown zero-query and excluded-only inputs remain in the class denominator with zero detections");
		final ArrayList<String> changes=lines(fixture.path("lineages.input_changes.tsv"));
		require(has(changes, comparison+"LAST_LOSS\tB\tBacteria\tP_α\trare\tknown\t2\t1\t0"),
			"The last-loss report must identify input B and its supplied taxonomy");
		require(has(changes, comparison+"FIRST_DETECTION\tC\tArchaea\tP_C\tbreadth\tknown\t1\t0\t1"),
			"The first-detection report must identify input C");
		require(!has(changes, comparison+"FIRST_DETECTION\t"+GENOME) && !has(changes, comparison+"LAST_LOSS\t"+GENOME),
			"A changed reference copy in an already detected input is not a first or last detection");
		final ArrayList<String> duplicate=lines(taxonomy);
		duplicate.add(duplicate.get(1));
		final String duplicatePath=fixture.path("duplicate-taxonomy.tsv");
		writeLines(duplicatePath, duplicate);
		reject(withOption(fixture.args("duplicate-taxonomy", fixture.times, fixture.raw, fixture.summary,
			fixture.manifest, fixture.profile), "taxonomy="+duplicatePath), "exactly one taxonomy row");
		final ArrayList<String> missing=lines(taxonomy);
		missing.remove(4);// The zero-query input still needs a taxonomy row.
		final String missingPath=fixture.path("missing-taxonomy.tsv");
		writeLines(missingPath, missing);
		reject(withOption(fixture.args("missing-taxonomy", fixture.times, fixture.raw, fixture.summary,
			fixture.manifest, fixture.profile), "taxonomy="+missingPath), "Every selected input needs declared taxonomy");
		final ArrayList<String> wrong=lines(taxonomy);
		wrong.set(1, wrong.get(1).replace("/"+GENOME+".fa", "/wrong.fa"));
		final String wrongPath=fixture.path("wrong-taxonomy.tsv");
		writeLines(wrongPath, wrong);
		reject(withOption(fixture.args("wrong-taxonomy", fixture.times, fixture.raw, fixture.summary,
			fixture.manifest, fixture.profile), "taxonomy="+wrongPath), "exact FASTA path");

		// Native grade fields escape a literal percent; taxonomy inputs are literal.
		final String encoded="percent%25_α", literal="percent%_α";
		final String raw=replace(fixture.raw, fixture.path("percent-raw.tsv"), GENOME, encoded);
		final String manifest=replace(fixture.manifest, fixture.path("percent-inputs.tsv"), GENOME, encoded);
		final String counts=replace(fixture.counts, fixture.path("percent-counts.tsv"), GENOME, encoded);
		final String percentTaxonomy=fixture.path("percent-taxonomy.tsv");
		writeTaxonomy(percentTaxonomy, literal);
		final String[] percent=withOption(fixture.args("percent", fixture.times, raw, fixture.summary, manifest,
			fixture.profile), "taxonomy="+percentTaxonomy);
		for(int i=0; i<percent.length; i++){if(percent[i].startsWith("counts=")){percent[i]="counts="+counts;}}
		SummarizeNcrnaBenchmark.main(percent);
		require(lines(fixture.path("percent.lineages.tsv")).equals(classes),
			"Literal percent and UTF-8 input/path metadata must bind the escaped grade keys without changing class counts");
		final String[] timing={"in="+fixture.times, "out="+fixture.path("taxonomy-timing-conflict"), "timingonly=t",
			"policies=current,"+TRIAL, "baseline=current", "taxonomy="+taxonomy};
		reject(timing, "does not accept grade or count inputs or taxonomy");
	}

	private static void writeTaxonomy(final String path, final String first){
		try(Output out=new Output(path)){
			text(out, "input_id\tfasta\tdomain\tphylum\tcorpus_inputs\tcategory\ttaxonomy_status");
			text(out, first+"\t/"+first+".fa\tBacteria\tP_α\t2\trare\tknown");
			text(out, "B\t/B.fa\tBacteria\tP_α\t2\trare\tknown");
			text(out, "C\t/C.fa\tArchaea\tP_C\t1\tbreadth\tknown");
			text(out, "empty\t/empty.fa\tunknown\tunknown\t2\tunknown\tunclassified");
			text(out, "excluded_only\t/excluded_only.fa\tunknown\tunknown\t2\tunknown\tunclassified");
			text(out, "outside\t/outside.fa\tBacteria\tOther\t1\tother\tknown");
		}
	}

	private static String[] withOption(final String[] args, final String option){
		final String[] extended=Arrays.copyOf(args, args.length+1);
		extended[args.length]=option;
		return extended;
	}

	private static String replace(final String input, final String output, final String from, final String to){
		final ArrayList<String> rows=lines(input);
		for(int i=0; i<rows.size(); i++){rows.set(i, rows.get(i).replace(from, to));}
		writeLines(output, rows);
		return output;
	}

	private static void reject(final String[] args, final String diagnostic){
		try{SummarizeNcrnaBenchmark.main(args);}
		catch(IllegalArgumentException failure){
			require(failure.getMessage().contains(diagnostic), "Wrong failure: "+failure.getMessage());
			return;
		}
		throw new AssertionError("Invalid native reports were accepted: "+diagnostic);
	}

	private static final class Fixture{
		Fixture(final File root_){root=root_;}
		String path(final String name){return new File(root, name).toString();}
		String[] args(final String name, final String time, final String raw_, final String summary_, final String inputs_, final String profile_){
			return new String[]{"in="+time, "raw="+raw_, "grades="+summary_, "gradeinputs="+inputs_, "counts="+counts,
				"profile="+profile_, "profiletime="+profileTime, "out="+path(name), "policies=current,"+TRIAL,
				"baseline=current", "workers=1,16", "repeats=3", "t=1", "outstream="+path(name+".log")};
		}
		void write(){
			times=path("times=a.tsv");
			raw=path("raw.tsv");
			summary=path("summary.tsv");
			manifest=path("inputs.tsv");
			counts=path("counts.tsv");
			profile=path("profile.tsv");
			profileTime=path("profile.time");
			try(Output time=new Output(times); Output rows=new Output(raw); Output grades=new Output(summary);
				Output inputs=new Output(manifest); Output original=new Output(counts)){
				text(time, "policy\tthreads\trepeat\tuser_seconds\tsystem_seconds\twall_seconds");
				text(rows, CONTEXT+GradeNcrnaReference.RAW_HEADER);
				text(grades, CONTEXT+GradeNcrnaReference.SUMMARY_HEADER);
				text(inputs, CONTEXT+GradeNcrnaReference.INPUT_HEADER);
				text(original, CONTEXT+"input_id\tsource_features\tCDS\ttRNA\trRNA\tprojected_RNA\tquery_features\tqualification");
				for(final String policy:new String[]{"current", TRIAL}){
					for(final int workers:new int[]{1, 16}){
						for(int repeat=1; repeat<=3; repeat++){
							final String prefix=policy+"\t"+workers+"\t"+repeat+"\t";
							final boolean trial=policy.equals(TRIAL), unstable=!trial && workers==16 && repeat==3;
							final int[] queryCounts={2, trial || unstable ? 0 : 1, trial || unstable ? 2 : 1, 0, 1};
							final int[][] observations;
							if(workers==1){observations=trial ? SINGLE_TRIAL : SINGLE_CURRENT;}
							else{observations=trial ? PARALLEL_TRIAL : PARALLEL_CURRENT;}
							final int[] measured=observations[repeat-1];
							text(time, prefix+measured[0]+"\t"+measured[1]+"\t"+measured[2]);
							text(grades, prefix+RNA+"\t5\t2\t1\t2\t4\t1\t1\t1\t0\t2\t2\t2\t0.666667\t0.500000\t0.500000\tFIXTURE");
							text(grades, prefix+"Excluded\t5\t0\t0\t0\t0\t0\t0\t0\t1\t0\t0\t0\tNA\tNA\t0.500000\tFIXTURE");
							for(int i=0; i<IDS.length; i++){
								text(inputs, prefix+IDS[i]+"\t/"+IDS[i]+".fa\t/query_"+i+"\t/reference_"+i+"\t"+queryCounts[i]+"\t0\t0\t0\t"+(i==4 ? 0 : queryCounts[i])+"\t"+(i==4 ? 1 : 0));
								text(original, prefix+IDS[i]+"\t"+(6+queryCounts[i])+"\t1\t2\t3\t"+queryCounts[i]+"\t"+queryCounts[i]+"\tFIXTURE");
							}
							ref(rows, prefix, GENOME, "ref_β", !trial ? "TP" : "FN");
							ref(rows, prefix, GENOME, "ref_B", trial ? "TP" : "FN");
							ref(rows, prefix, "B", "ref_C", !trial && !unstable ? "TP" : "FN");
							ref(rows, prefix, "C", "ref_D", trial || unstable ? "TP" : "FN");
							ref(rows, prefix, "C", "mask", "DONT_CARE");
							query(rows, prefix, GENOME, RNA, "TP", trial ? "ref_B" : "ref_β", "false");
							query(rows, prefix, trial || unstable ? "C" : "B", RNA, "TP", trial || unstable ? "ref_D" : "ref_C", "false");
							query(rows, prefix, GENOME, RNA, "FP", ".", "true");
							query(rows, prefix, "C", RNA, "MASKED", ".", "false");
							query(rows, prefix, "excluded_only", "Excluded", "EXCLUDED", ".", ".");
						}
					}
				}
			}
			writeProfile(profile, true);
			try(Output out=new Output(profileTime)){text(out, "10\t2\t5");}
		}
		void writeProfile(final String path, final boolean named){
			try(Output out=new Output(path)){
				text(out, "# Native diagnostic profile fixture");
				text(out, "family\tstage\tinvocations\tcpu_ns\talignments\taligned_bases\thbm_calls\thbm_bases\tseed_hits\twindows");
				if(named){
					span(out, RNA, "INIT", 1000000000L);
					span(out, RNA, "SCAVENGE", 2000000000L);
				}
				span(out, "GLOBAL", "SHARED_SEED_SCAN", 3000000000L);
				span(out, "GLOBAL", "SSU_APPLY", 500000000L);
				span(out, "GLOBAL", "SSU_FINISH", 0);
				span(out, "GLOBAL", "FIND_PATH", 1000000000L);
				span(out, "GLOBAL", "TRNA", 1500000000L);
			}
		}
		final File root;
		String times, raw, summary, manifest, counts, profile, profileTime;
	}

	private static void ref(final Output out, final String prefix, final String input, final String reference, final String outcome){
		text(out, prefix+input+"\tREFERENCE\t"+RNA+"\tcontig\t1\t20\t+\t"+outcome+"\t"+reference+"\t.\t.\t.\t.\t.\t.\tmodel\t.\t.");
	}
	private static void query(final Output out, final String prefix, final String input, final String family, final String outcome, final String reference, final String duplicate){
		text(out, prefix+input+"\tQUERY\t"+family+"\tcontig\t1\t20\t+\t"+outcome+"\t.\t2\t"+reference+"\t.\t.\t"+duplicate+"\t.\tmodel\t"+(outcome.equals("TP") ? "0\t0" : ".\t."));
	}
	private static void span(final Output out, final String family, final String stage, final long cpu){
		text(out, family+"\t"+stage+"\t1\t"+cpu+"\tNA\tNA\tNA\tNA\tNA\tNA");
	}
	private static void text(final Output out, final String line){out.writer.println(line.getBytes(StandardCharsets.UTF_8));}
	private static ArrayList<String> lines(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){result.add(new String(row, StandardCharsets.UTF_8));}
		}
		return result;
	}
	private static boolean has(final ArrayList<String> lines, final String prefix){
		for(final String line:lines){if(line.startsWith(prefix) || line.contains(prefix)){return true;}}
		return false;
	}
	/** Checks literal numeric oracles to the report's six-place rounding, not its optional zeroes. */
	private static void numericRow(final ArrayList<String> rows, final String prefix, final int first, final double... expected){
		final LineParser1 parser=new LineParser1('\t');
		for(final String row:rows){
			if(!row.startsWith(prefix)){continue;}
			parser.set(row.getBytes(StandardCharsets.UTF_8));
			require(parser.terms()>=first+expected.length, "Missing numeric columns: "+row);
			for(int i=0; i<expected.length; i++){
				final double actual=parser.parseDouble(first+i);
				require(Math.abs(actual-expected[i])<=0.000000500001,
					"Numeric oracle mismatch for "+prefix+" field="+(first+i)+" expected="+expected[i]+" actual="+actual);
			}
			return;
		}
		throw new AssertionError("Missing numeric report row: "+prefix);
	}
	private static ArrayList<String> withoutLast(final ArrayList<String> source){
		final ArrayList<String> result=new ArrayList<String>(source);
		result.remove(result.size()-1);
		return result;
	}
	private static void writeLines(final String path, final ArrayList<String> lines){
		try(Output out=new Output(path)){for(final String line:lines){text(out, line);}}
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private static final String CONTEXT="policy\tthreads\trepeat\t";
	private static final String RNA="RNA_δ", GENOME="genome_α", TRIAL="changed=γ";
	private static final String[] IDS={GENOME, "B", "C", "empty", "excluded_only"};
	/** Literal user/system/wall observations; these are independent statistical oracles. */
	private static final int[][] SINGLE_CURRENT={{10,100,8}, {20,10,9}, {30,100,10}};
	private static final int[][] SINGLE_TRIAL={{20,1,9}, {10,1,8}, {40,1,7}};
	private static final int[][] PARALLEL_CURRENT={{40,0,5}, {50,0,4}, {60,0,3}};
	private static final int[][] PARALLEL_TRIAL={{50,1,4}, {45,1,5}, {55,1,6}};
}
