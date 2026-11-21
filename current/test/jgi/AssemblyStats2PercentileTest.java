package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import parse.Parse;

/**
 * Ten fresh stats.sh child-process cases for AssemblyStats2 N50/L50 and
 * N90/L90 across the small-length histogram and explicit long-record list.
 * The baseline oracle is source-derived; the corrected oracle records the
 * descending minimum-prefix percentile that should hold after repair.
 * @author Jean
 */
public final class AssemblyStats2PercentileTest{

	/** Requires a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown AssemblyStats2 percentile option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblystats2percentiletest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final AssemblyStats2PercentileTest test=new AssemblyStats2PercentileTest(out, repo,
			phase.equals("baseline"));
		test.writeInputs();
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==10, "Expected ten AssemblyStats2 percentile cases: "+test.cases);
		test.check(test.children==10, "Expected one fresh stats.sh child per case: "+test.children);
		test.check(test.sourceMatches==(test.baseline ? 10 : 4),
			"Unexpected source-derived matches: "+test.sourceMatches);
		test.check(test.correctedMatches==(test.baseline ? 4 : 10),
			"Unexpected corrected-oracle matches: "+test.correctedMatches);
		final String summary="ASSEMBLY_STATS2_PERCENTILE_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" sourceMatches="+test.sourceMatches+" correctedMatches="+test.correctedMatches+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the script under test. */
	private AssemblyStats2PercentileTest(final Path out_, final Path repo_,
			final boolean baseline_) throws IOException{
		out=out_; repo=repo_; inputs=out.resolve("inputs"); baseline=baseline_;
		launcher=repo.resolve("stats.sh");
		require(Files.isRegularFile(launcher), "Missing stats.sh at "+launcher);
		Files.createDirectory(inputs);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\tinput\tcutoff\texit"+
			"\tstdoutBytes\tstderrBytes\tparseError"+
			"\tsourceN50\tsourceL50\tsourceN90\tsourceL90"+
			"\tcorrectedN50\tcorrectedL50\tcorrectedN90\tcorrectedL90"+
			"\tobservedScafN50\tobservedScafL50\tobservedScafN90\tobservedScafL90"+
			"\tobservedCtgN50\tobservedCtgL50\tobservedCtgN90\tobservedCtgL90"+
			"\tsourceMatches\tcorrectedMatches\ttotalsMatch\tstdoutEmpty\tstderrClean"+
			"\tinputsUnchanged\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
	}

	/** Writes the four retained no-N FASTA fixtures with literal lengths. */
	private void writeInputs() throws IOException{
		write("equal10.fa", 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10);
		write("one9.fa", 1, 9);
		write("mixed.fa", 3, 7, 11);
		write("exact50.fa", 7, 5, 2);
	}

	/** Writes one all-A FASTA with exactly one record of each requested length. */
	private void write(final String name, final int... lengths) throws IOException{
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<lengths.length; i++){
			final char[] bases=new char[lengths[i]];
			Arrays.fill(bases, 'A');
			sb.append('>').append(name).append('_').append(i).append('\n');
			sb.append(bases).append('\n');
		}
		final Path path=inputs.resolve(name);
		final byte[] data=sb.toString().getBytes(StandardCharsets.US_ASCII);
		Files.write(path, data, StandardOpenOption.CREATE_NEW);
		original.put(path, data);
	}

	/** Launches one stats.sh process and checks source and corrected percentile oracles. */
	private void run(final Row row) throws Exception{
		cases++;
		children++;
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		final Path input=inputs.resolve(row.inputName);
		final Path output=caseDir.resolve("stats.tsv");
		writeExpected(caseDir, row);

		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("format=3");
		command.add("swapnl=t");
		command.add("minscaf=0");
		command.add("n90=t");
		command.add("n_=f");
		command.add("simd=f");
		command.add("t=1");
		command.add("cutoff="+row.cutoff);
		command.add("in="+input);
		command.add("out="+output);
		Files.write(caseDir.resolve("argv.txt"), join(command).getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);

		final Path stdoutPath=caseDir.resolve("stdout.txt");
		final Path stderrPath=caseDir.resolve("stderr.txt");
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(stdoutPath.toFile()).redirectError(stderrPath.toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int exit=builder.start().waitFor();
		final byte[] stdout=Files.readAllBytes(stdoutPath);
		final byte[] stderr=Files.readAllBytes(stderrPath);
		final String errText=new String(stderr, StandardCharsets.UTF_8);
		Files.write(caseDir.resolve("exit.txt"), (Integer.toString(exit)+"\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);

		final boolean outputPresent=Files.isRegularFile(output);
		Map<String, String> stats=new HashMap<String, String>();
		String parseError="";
		if(outputPresent){
			try {
				stats=parseStats(output);
			} catch(final RuntimeException e) {
				parseError=e.toString();
			} catch(final IOException e) {
				parseError=e.toString();
			}
			if(!parseError.isEmpty()){
				Files.write(caseDir.resolve("parse_error.txt"),
					(parseError+"\n").getBytes(StandardCharsets.UTF_8),
					StandardOpenOption.CREATE_NEW);
			}
		}
		final long scafN50=getLength(stats, "scaf_N50");
		final long scafL50=getLong(stats, "scaf_L50");
		final long scafN90=getLength(stats, "scaf_N90");
		final long scafL90=getLong(stats, "scaf_L90");
		final long ctgN50=getLength(stats, "ctg_N50");
		final long ctgL50=getLong(stats, "ctg_L50");
		final long ctgN90=getLength(stats, "ctg_N90");
		final long ctgL90=getLong(stats, "ctg_L90");

		final boolean sourceMatch=scafN50==row.sourceN50 && scafL50==row.sourceL50 &&
			scafN90==row.sourceN90 && scafL90==row.sourceL90 &&
			ctgN50==row.sourceN50 && ctgL50==row.sourceL50 &&
			ctgN90==row.sourceN90 && ctgL90==row.sourceL90;
		final boolean correctedMatch=scafN50==row.correctedN50 && scafL50==row.correctedL50 &&
			scafN90==row.correctedN90 && scafL90==row.correctedL90 &&
			ctgN50==row.correctedN50 && ctgL50==row.correctedL50 &&
			ctgN90==row.correctedN90 && ctgL90==row.correctedL90;
		final boolean totalsMatch=totalsMatch(row, stats);
		final boolean stdoutEmpty=stdout.length==0;
		final boolean stderrClean=errText.contains("Picked up JAVA_TOOL_OPTIONS") &&
			!errText.contains("Unknown parameter") && !errText.contains("Exception") &&
			!errText.contains("AssertionError") && !errText.contains("Error reading");
		final boolean unchanged=inputsUnchanged();
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}

		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+row.inputName+"\t"+row.cutoff+
			"\t"+exit+"\t"+stdout.length+"\t"+stderr.length+"\t"+parseError+
			"\t"+row.sourceN50+"\t"+row.sourceL50+"\t"+row.sourceN90+"\t"+row.sourceL90+
			"\t"+row.correctedN50+"\t"+row.correctedL50+
			"\t"+row.correctedN90+"\t"+row.correctedL90+
			"\t"+scafN50+"\t"+scafL50+"\t"+scafN90+"\t"+scafL90+
			"\t"+ctgN50+"\t"+ctgL50+"\t"+ctgN90+"\t"+ctgL90+
			"\t"+sourceMatch+"\t"+correctedMatch+"\t"+totalsMatch+
			"\t"+stdoutEmpty+"\t"+stderrClean+"\t"+unchanged+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);

		check(exit==0, row.name+": child exit was "+exit);
		check(outputPresent, row.name+": expected stats output missing: "+output);
		check(parseError.isEmpty(), row.name+": could not parse format=3 output: "+parseError);
		check(baseline ? sourceMatch : correctedMatch,
			row.name+": "+(baseline ? "source-derived" : "corrected")+" N/L oracle mismatch at "+output);
		check(totalsMatch, row.name+": non-percentile totals mismatch at "+output);
		check(stdoutEmpty, row.name+": stats.sh wrote primary table to stdout despite out=");
		check(stderrClean, row.name+": stderr contains an unexpected diagnostic: "+errText);
		check(command.contains("-ea") && command.contains("format=3") &&
			command.contains("swapnl=t") && command.contains("minscaf=0") &&
			command.contains("n90=t") && command.contains("n_=f") &&
			command.contains("simd=f") && command.contains("t=1"),
			row.name+": child argv lost a frozen CLI pin");
		check(unchanged, row.name+": a child modified literal input bytes");
	}

	/** Writes the source-derived baseline and corrected descending oracle before launch. */
	private static void writeExpected(final Path caseDir, final Row row) throws IOException{
		final String text="stat\tsource\tcorrected\n"+
			"N50\t"+row.sourceN50+"\t"+row.correctedN50+"\n"+
			"L50\t"+row.sourceL50+"\t"+row.correctedL50+"\n"+
			"N90\t"+row.sourceN90+"\t"+row.correctedN90+"\n"+
			"L90\t"+row.sourceL90+"\t"+row.correctedL90+"\n";
		Files.write(caseDir.resolve("expected.tsv"), text.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
	}

	/** Parses the one-row format=3 TSV through a full header map. */
	private static Map<String, String> parseStats(final Path output) throws IOException{
		final List<String> lines=Files.readAllLines(output, StandardCharsets.UTF_8);
		require(lines.size()==2, "Expected header plus one data row in "+output+": "+lines.size());
		final String[] header=lines.get(0).split("\t", -1);
		final String[] values=lines.get(1).split("\t", -1);
		require(header.length==values.length,
			"Header/value width mismatch in "+output+": "+header.length+" vs "+values.length);
		final HashMap<String, String> map=new HashMap<String, String>();
		for(int i=0; i<header.length; i++){
			require(!map.containsKey(header[i]), "Duplicate output column in "+output+": "+header[i]);
			map.put(header[i], values[i]);
		}
		for(final String key : REQUIRED_COLUMNS){
			require(map.containsKey(key), "Missing required output column in "+output+": "+key);
		}
		return map;
	}

	/** Checks the non-percentile columns that should be invariant for these no-N inputs. */
	private static boolean totalsMatch(final Row row, final Map<String, String> stats){
		return getLong(stats, "scaffolds")==row.records &&
			getLong(stats, "contigs")==row.records &&
			getLong(stats, "scaf_bp")==row.bases &&
			getLong(stats, "contig_bp")==row.bases &&
			getLength(stats, "scaf_max")==row.max &&
			getLength(stats, "ctg_max")==row.max &&
			getLong(stats, "scaf_n_gt50K")==0 &&
			getDouble(stats, "scaf_pct_gt50K")==0 &&
			getDouble(stats, "gap_pct")==0 &&
			getDouble(stats, "gc_avg")==0 &&
			getDouble(stats, "gc_std")==0;
	}

	/** Fetches one whole-number length from a named output column. */
	private static long getLength(final Map<String, String> stats, final String key){
		final String value=stats.get(key);
		if(value==null || value.isEmpty()){return Long.MIN_VALUE;}
		try {
			return Long.parseLong(value.replace(",", ""));
		} catch(final NumberFormatException e) {
			return Long.MIN_VALUE;
		}
	}

	/** Fetches one integer from a named output column. */
	private static long getLong(final Map<String, String> stats, final String key){
		final String value=stats.get(key);
		if(value==null || value.isEmpty()){return Long.MIN_VALUE;}
		try {
			return Long.parseLong(value.replace(",", ""));
		} catch(final NumberFormatException e) {
			return Long.MIN_VALUE;
		}
	}

	/** Fetches one floating-point value from a named output column. */
	private static double getDouble(final Map<String, String> stats, final String key){
		final String value=stats.get(key);
		if(value==null || value.isEmpty()){return Double.NaN;}
		try {
			return Double.parseDouble(value.replace(",", ""));
		} catch(final NumberFormatException e) {
			return Double.NaN;
		}
	}

	/** Checks every retained input against the literal bytes written before children ran. */
	private boolean inputsUnchanged() throws IOException{
		for(final Map.Entry<Path, byte[]> entry : original.entrySet()){
			if(!Arrays.equals(entry.getValue(), Files.readAllBytes(entry.getKey()))){return false;}
		}
		return true;
	}

	/** Joins argv for stable evidence. */
	private static String join(final List<String> command){
		final StringBuilder sb=new StringBuilder();
		for(final String s : command){sb.append(s).append('\n');}
		return sb.toString();
	}

	/** Accumulates mismatches so every child leaves stdout, stderr and exit evidence. */
	private void check(final boolean correct, final String message){
		checks++;
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}

	/** Throws on fixture misuse. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** One child-process case with source-derived and corrected percentile oracles. */
	private static final class Row{

		Row(final String name_, final String inputName_, final int cutoff_,
				final long sourceN50_, final long sourceL50_,
				final long sourceN90_, final long sourceL90_,
				final long correctedN50_, final long correctedL50_,
				final long correctedN90_, final long correctedL90_,
				final long records_, final long bases_, final long max_){
			name=name_;
			inputName=inputName_;
			cutoff=cutoff_;
			sourceN50=sourceN50_;
			sourceL50=sourceL50_;
			sourceN90=sourceN90_;
			sourceL90=sourceL90_;
			correctedN50=correctedN50_;
			correctedL50=correctedL50_;
			correctedN90=correctedN90_;
			correctedL90=correctedL90_;
			records=records_;
			bases=bases_;
			max=max_;
		}

		final String name, inputName;
		final int cutoff;
		final long sourceN50, sourceL50, sourceN90, sourceL90;
		final long correctedN50, correctedL50, correctedN90, correctedL90;
		final long records, bases, max;
	}

	private final Path out, repo, inputs, launcher;
	private final Map<Path, byte[]> original=new HashMap<Path, byte[]>();
	private int cases, children, sourceMatches, correctedMatches, checks, failures;
	private final boolean baseline;

	private static final String[] REQUIRED_COLUMNS={
		"scaffolds", "contigs", "scaf_bp", "contig_bp", "gap_pct",
		"scaf_N50", "scaf_L50", "ctg_N50", "ctg_L50",
		"scaf_N90", "scaf_L90", "ctg_N90", "ctg_L90",
		"scaf_max", "ctg_max", "scaf_n_gt50K", "scaf_pct_gt50K", "gc_avg", "gc_std"
	};

	private static final Row[] ROWS={
		new Row("equal10-cutoff9", "equal10.fa", 9, 10, 6, 10, 10, 10, 6, 10, 10, 11, 110, 10),
		new Row("equal10-cutoff10", "equal10.fa", 10, 10, 6, 10, 10, 10, 6, 10, 10, 11, 110, 10),
		new Row("equal10-cutoff11", "equal10.fa", 11, 10, 11, 10, 11, 10, 6, 10, 10, 11, 110, 10),
		new Row("one9-cutoff2", "one9.fa", 2, 9, 1, 1, 2, 9, 1, 9, 1, 2, 10, 9),
		new Row("one9-cutoff9", "one9.fa", 9, 9, 1, 1, 2, 9, 1, 9, 1, 2, 10, 9),
		new Row("one9-cutoff10", "one9.fa", 10, 9, 1, 1, 2, 9, 1, 9, 1, 2, 10, 9),
		new Row("mixed-cutoff8", "mixed.fa", 8, 11, 1, 3, 3, 11, 1, 3, 3, 3, 21, 11),
		new Row("mixed-cutoff12", "mixed.fa", 12, 11, 1, 3, 3, 11, 1, 3, 3, 3, 21, 11),
		new Row("exact50-cutoff6", "exact50.fa", 6, 5, 2, 2, 3, 7, 1, 2, 3, 3, 14, 7),
		new Row("exact50-cutoff8", "exact50.fa", 8, 5, 2, 2, 3, 7, 1, 2, 3, 3, 14, 7),
	};
}
