package test.jgi;

import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import fileIO.ReadWrite;
import jgi.BBMerge;
import parse.Parse;
import shared.Tools;

/** Bounded helper and writer cases for BBMerge insert percentiles.
 * Direct helper cases compare only literal finite and large p50 outcomes.
 * Small writer children use BBMerge's actual 2000-entry histogram, enable
 * histogram statistics, write only nonzero rows and invoke the public writer.
 * Baseline helper calls use the current shared floor-rank helper; fixed helper
 * calls resolve BBMerge.insertPercentile(long[], int) reflectively only in fixed
 * mode so this fixture compiles before the product repair. Percentiles are
 * limited to integer p in {10,25,50,75,90}. No large-count writer case, parser
 * copy, failure/lifecycle probe, computed rank oracle or shared Tools change is
 * part of this fixture.
 * @author Jean
 */
public final class BBMergePercentileTest{

	/** mode=suite|child, phase=baseline|fixed, root, new outdir; child adds case=name. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "root", "outdir", "reference", "case").contains(pair[0]),
					"Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate fixture option: "+pair[0]);
		}
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergepercentiletest.sh with assertions enabled.");
		final String mode=options.containsKey("mode") ? options.get("mode") : "suite";
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed for the percentile repair.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		if(mode.equals("child")){
			require(!options.containsKey("reference") && Files.isDirectory(out), "Child uses an existing case directory and no reference.");
			child(out, small(options.get("case")));
			return;
		}
		require(mode.equals("suite") && !options.containsKey("case"), "Use suite mode or an internal child case.");
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require("baseline".equals(phase) ? reference==null : reference!=null,
				"Only fixed suite needs reference=<successful baseline fixture outdir>.");
		if(reference!=null){
			require(read(reference.resolve("summary.txt")).startsWith("BBMERGE_PERCENTILE_PASS phase=baseline "),
					"Reference is not a successful baseline percentile fixture directory.");
		}
		final BBMergePercentileTest test=new BBMergePercentileTest(root, out, reference, phase);
		for(final SmallSpec spec : SMALL){test.runHelpers(spec); test.runWriter(spec);}
		for(final LargeSpec spec : LARGE){test.runHelpers(spec);}
		test.check(test.helperOutcomes==33, "The fixture must retain exactly 33 helper outcomes.");
		test.check(test.children==6 && test.changed==2 && test.controls==4, "Six small writer children must contain two changed medians and four controls.");
		final String summary="BBMERGE_PERCENTILE_"+(test.failures==0 ? "PASS" : "FAIL")+" phase="+phase+
				" helperOutcomes="+test.helperOutcomes+" writerChildren="+test.children+" changed="+test.changed+
				" controls="+test.controls+" checks="+test.checks+" failures="+test.failures+"\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
		if(test.failures>0){throw new AssertionError(summary+"Inspect retained helper and writer outputs; do not fit literal oracles to output.");}
	}

	/** Creates a fresh evidence directory and resolves the proposed helper only after repair. */
	private BBMergePercentileTest(final Path root_, final Path out_, final Path reference_, final String phase_) throws Exception{
		root=root_; out=out_; reference=reference_; phase=phase_;
		if("fixed".equals(phase)){
			fixedMethod=BBMerge.class.getDeclaredMethod("insertPercentile", long[].class, int.class);
			fixedMethod.setAccessible(true);
		}else{fixedMethod=null;}
		Files.createDirectory(out);
		write(out.resolve("helper.tsv"), "case\tpercentile\tbaselineExpected\tfixedExpected\tactual\tphaseMatches\n");
		write(out.resolve("writer.tsv"), "case\tchildExit\tmedian\tphaseMedian\trawRows\tphaseOutput\tinputUnchanged\n");
	}

	/** Compares actual helper behavior to literal old/new vectors; no rank oracle is computed. */
	private void runHelpers(final Spec spec) throws Exception{
		final int[] percentiles=spec.percentiles();
		for(int i=0; i<percentiles.length; i++){
			final int p=percentiles[i];
			final long[] histogram=spec.histogram();
			final long[] before=histogram.clone();
			final int actual="baseline".equals(phase) ? Tools.percentileHistogram(histogram, p/100.0) :
				((Integer)fixedMethod.invoke(null, histogram, p)).intValue();
			final int expected="baseline".equals(phase) ? spec.baseline[i] : spec.fixed[i];
			helperOutcomes++;
			check(actual==expected, spec.name+": p"+p+" expected "+expected+" but saw "+actual);
			check(Arrays.equals(before, histogram), spec.name+": p"+p+" helper mutated the input histogram.");
			append(out.resolve("helper.tsv"), spec.name+"\t"+p+"\t"+spec.baseline[i]+"\t"+spec.fixed[i]+"\t"+actual+"\t"+(actual==expected)+"\n");
		}
	}

	/** Reads output only after Process.waitFor confirms complete child-JVM termination. */
	private void runWriter(final SmallSpec spec) throws Exception{
		if(spec.baseline[2]==spec.fixed[2]){controls++;}else{changed++;}
		final Path work=Files.createDirectory(out.resolve(spec.name));
		write(work.resolve("expected-baseline-median.txt"), "#Median\t"+spec.baseline[2]+"\n");
		write(work.resolve("expected-fixed-median.txt"), "#Median\t"+spec.fixed[2]+"\n");
		write(work.resolve("expected-raw-rows.tsv"), spec.rawRows);
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve("dev/bbmergepercentiletest.sh").toString(),
				"mode=child", "phase="+phase, "root="+root, "outdir="+work, "case="+spec.name));
		final StringBuilder saved=new StringBuilder();
		for(int i=0; i<command.size(); i++){saved.append(i).append('\t').append(command.get(i)).append('\n');}
		write(work.resolve("child.command.tsv"), saved.toString());
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile());
		builder.redirectOutput(work.resolve("child.stdout").toFile());
		builder.redirectError(work.resolve("child.stderr").toFile());
		builder.environment().remove("_JAVA_OPTIONS");
		builder.environment().remove("JDK_JAVA_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final Process process=builder.start(); children++;
		final int code=process.waitFor();
		write(work.resolve("child.exit"), code+"\n");
		check(code==0, spec.name+": child or counted writer failed.");
		final Path histogram=work.resolve("histogram.tsv");
		check(Files.isRegularFile(histogram), spec.name+": completed child did not leave a histogram.");
		if(code!=0 || !Files.isRegularFile(histogram)){
			append(out.resolve("writer.tsv"), spec.name+"\t"+code+"\tfalse\tfalse\tfalse\tfalse\tfalse\n");
			return;
		}
		final String actual=read(histogram);
		final int medians=countLines(actual, "#Median\t");
		final String phaseMedian="#Median\t"+("baseline".equals(phase) ? spec.baseline[2] : spec.fixed[2])+"\n";
		final boolean medianOk=medians==1 && actual.contains(phaseMedian);
		final boolean rawOk=rawRows(actual).equals(spec.rawRows);
		check(medianOk, spec.name+": #Median line differs from the literal phase expectation.");
		check(rawOk, spec.name+": raw histogram rows differ from the literal expected rows.");
		final boolean inputUnchanged=read(work.resolve("child-result.tsv")).equals("capacity\t2000\nhistogramUnchanged\ttrue\ncountedWritesDrained\ttrue\n");
		check(inputUnchanged, spec.name+": histogram input or counted-write completion evidence differs.");
		boolean phaseOutput=medianOk && rawOk;
		if(reference!=null){
			final Path before=reference.resolve(spec.name);
			for(final String file : new String[]{"setup.tsv", "inputs.tsv", "expected-baseline-median.txt", "expected-fixed-median.txt", "expected-raw-rows.tsv"}){
				check(Arrays.equals(Files.readAllBytes(work.resolve(file)), Files.readAllBytes(before.resolve(file))),
						spec.name+": case settings or literal oracles changed between phases: "+file);
			}
			final String normalized=replaceMedian(read(before.resolve("histogram.tsv")), spec.fixed[2]);
			phaseOutput=actual.equals(normalized);
			check(phaseOutput, spec.name+": fixed writer changed bytes other than the accepted #Median line.");
		}
		append(out.resolve("writer.tsv"), spec.name+"\t"+code+"\t"+medianOk+"\t"+(phase.equals("baseline") ? spec.baseline[2] : spec.fixed[2])+"\t"+
				rawOk+"\t"+phaseOutput+"\t"+inputUnchanged+"\n");
	}

	/** Changes test inputs only, never the real histogram array's capacity or writer code. */
	private static void child(final Path work, final SmallSpec spec) throws Exception{
		final long[] histogram=(long[])field("histTotal").get(null);
		require(histogram.length==2000 && field("histlen").getInt(null)==2000,
				"Use the actual histogram capacity; a compact replacement is not a production truth case.");
		Arrays.fill(histogram, 0);
		for(int i=0; i<spec.indices.length; i++){histogram[spec.indices[i]]=spec.counts[i];}
		field("bin").setInt(null, 1);
		field("insertMaxTotal").setInt(null, spec.maximum);
		field("showHistStats").setBoolean(null, true);
		field("NONZERO_ONLY").setBoolean(null, true);
		require(field("bin").getInt(null)==1 && field("insertMaxTotal").getInt(null)==spec.maximum &&
				field("showHistStats").getBoolean(null) && field("NONZERO_ONLY").getBoolean(null),
				"Actual writer settings differ from the isolated percentile case.");
		write(work.resolve("setup.tsv"), "capacity\t2000\nbin\t1\nmaximumInsert\t"+spec.maximum+"\nshowHistStats\ttrue\nnonzeroOnly\ttrue\n");
		final StringBuilder inputs=new StringBuilder("index\tcount\n");
		for(int i=0; i<histogram.length; i++){if(histogram[i]!=0){inputs.append(i).append('\t').append(histogram[i]).append('\n');}}
		write(work.resolve("inputs.tsv"), inputs.toString());
		final long[] before=histogram.clone();
		require(ReadWrite.countActiveThreads()==0, "Fresh child must start with no counted writes.");
		BBMerge.writeHistogram(work.resolve("histogram.tsv").toString(), 0);
		ReadWrite.waitForWritingToFinish();
		require(ReadWrite.countActiveThreads()==0 && Arrays.equals(before, histogram),
				"Writer must drain its counted write and preserve histogram counts.");
		write(work.resolve("child-result.tsv"), "capacity\t2000\nhistogramUnchanged\ttrue\ncountedWritesDrained\ttrue\n");
	}

	private static abstract class Spec{
		Spec(final String name_, final int[] baseline_, final int[] fixed_){
			name=name_; baseline=baseline_; fixed=fixed_;
		}
		abstract long[] histogram();
		abstract int[] percentiles();
		final String name;
		final int[] baseline, fixed;
	}

	private static final class SmallSpec extends Spec{
		SmallSpec(final String name_, final int maximum_, final int[] indices_, final long[] counts_,
				final int[] baseline_, final int[] fixed_, final String rawRows_){
			super(name_, baseline_, fixed_);
			maximum=maximum_; indices=indices_; counts=counts_; rawRows=rawRows_;
			require(indices.length==counts.length && maximum<2000 && baseline.length==5 && fixed.length==5,
					"Literal small cases must remain bounded and complete.");
		}
		long[] histogram(){
			final long[] hist=new long[101];
			for(int i=0; i<indices.length; i++){hist[indices[i]]=counts[i];}
			return hist;
		}
		int[] percentiles(){return PERCENTILES;}
		final String rawRows;
		final int maximum;
		final int[] indices;
		final long[] counts;
	}

	private static final class LargeSpec extends Spec{
		LargeSpec(final String name_, final int[] indices_, final long[] counts_, final int baseline_, final int fixed_){
			super(name_, new int[]{baseline_}, new int[]{fixed_});
			indices=indices_; counts=counts_;
			require(indices.length==counts.length, "Literal large cases must have matching index/count arrays.");
		}
		long[] histogram(){
			final long[] hist=new long[31];
			for(int i=0; i<indices.length; i++){hist[indices[i]]=counts[i];}
			return hist;
		}
		int[] percentiles(){return P50;}
		final int[] indices;
		final long[] counts;
	}

	private static SmallSpec small(final String name){
		for(final SmallSpec spec : SMALL){if(spec.name.equals(name)){return spec;}}
		throw new AssertionError("Unknown bounded percentile writer case: "+name);
	}

	private static Field field(final String name) throws Exception{
		final Field f=BBMerge.class.getDeclaredField(name);
		f.setAccessible(true);
		return f;
	}

	private static int countLines(final String text, final String prefix){
		int count=0;
		for(final String line : text.split("\n", -1)){if(line.startsWith(prefix)){count++;}}
		return count;
	}

	private static String rawRows(final String text){
		final String header="#InsertSize\tCount\n";
		final int index=text.indexOf(header);
		require(index>=0, "Histogram output is missing the raw insert-size header.");
		return text.substring(index+header.length());
	}

	private static String replaceMedian(final String text, final int median){
		final String[] lines=text.split("\n", -1);
		int count=0;
		for(int i=0; i<lines.length; i++){
			if(lines[i].startsWith("#Median\t")){
				lines[i]="#Median\t"+median;
				count++;
			}
		}
		require(count==1, "Reference histogram must contain exactly one #Median line.");
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<lines.length; i++){
			if(i>0){sb.append('\n');}
			sb.append(lines[i]);
		}
		return sb.toString();
	}

	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}

	private static void require(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private static Path path(final String value){require(value!=null && !value.isEmpty(), "Missing required path."); return Paths.get(value).toAbsolutePath().normalize();}
	private static void write(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);}
	private static String read(final Path path) throws IOException{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}

	private final Path root, out, reference;
	private final String phase;
	private final Method fixedMethod;
	private int checks, failures, helperOutcomes, children, changed, controls;

	private static final int[] PERCENTILES={10, 25, 50, 75, 90};
	private static final int[] P50={50};
	private static final SmallSpec[] SMALL={
		new SmallSpec("empty", 0, new int[0], new long[0],
				new int[]{0, 0, 0, 0, 0}, new int[]{0, 0, 0, 0, 0}, ""),
		new SmallSpec("single40", 40, new int[]{40}, new long[]{1},
				new int[]{0, 0, 0, 0, 0}, new int[]{40, 40, 40, 40, 40}, "40\t1\n"),
		new SmallSpec("odd10_20_30", 30, new int[]{10, 20, 30}, new long[]{1, 1, 1},
				new int[]{0, 0, 10, 20, 20}, new int[]{10, 10, 20, 30, 30}, "10\t1\n20\t1\n30\t1\n"),
		new SmallSpec("even10_20_30_40", 40, new int[]{10, 20, 30, 40}, new long[]{1, 1, 1, 1},
				new int[]{0, 10, 20, 30, 30}, new int[]{10, 10, 20, 30, 40}, "10\t1\n20\t1\n30\t1\n40\t1\n"),
		new SmallSpec("weighted10_30", 30, new int[]{10, 30}, new long[]{2, 1},
				new int[]{0, 0, 10, 10, 10}, new int[]{10, 10, 10, 30, 30}, "10\t2\n30\t1\n"),
		new SmallSpec("deciles10_100", 100, new int[]{10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, new long[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
				new int[]{10, 20, 50, 70, 90}, new int[]{10, 30, 50, 80, 90},
				"10\t1\n20\t1\n30\t1\n40\t1\n50\t1\n60\t1\n70\t1\n80\t1\n90\t1\n100\t1\n")
	};
	private static final LargeSpec[] LARGE={
		new LargeSpec("large_odd_stable", new int[]{10, 20, 30},
				new long[]{4611686018427387903L, 1L, 4611686018427387903L}, 20, 20),
		new LargeSpec("large_even_rounds_high", new int[]{10, 20},
				new long[]{4611686018427387903L, 4611686018427387903L}, 20, 10),
		new LargeSpec("large_beyond_long", new int[]{10, 20},
				new long[]{4611686018427387904L, 4611686018427387905L}, 10, 20)
	};
}
