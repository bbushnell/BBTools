package test.jgi;

import java.io.IOException;
import java.lang.reflect.Field;
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

/** Seven bounded direct histogram-writer cases for JGI088/JGI091.
 * Each child uses BBMerge's actual 2000-entry array, disables histogram statistics,
 * enables nonzero-only rows and invokes the public writer. Literal baseline wrong
 * bytes and corrected bytes are separate oracles; no billions of reads or alternate
 * implementation produce the expected text. Full bins retain ceiling-average
 * semantics; the partial tail divides by its actual entry count. Cases never need
 * a bin sum larger than Long.MAX_VALUE. Negative counts, oversized sums, histogram
 * collection, statistics and merger accuracy are outside this writer-only test.
 * Counted write completion is not a thread join. The parent waits for the entire
 * child JVM to exit before reading histogram bytes.
 * @author Sucrose
 */
public final class BBMergeHistogramTest{

	/** mode=suite|child, phase=baseline|fixed, root, new outdir; child adds case=name.
	 * Fixed suite reference is the successful baseline FIXTURE directory. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "root", "outdir", "reference", "case").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate fixture option: "+pair[0]);
		}
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergehistogramtest.sh with assertions enabled.");
		final String mode=options.containsKey("mode") ? options.get("mode") : "suite", phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed for the histogram writer repair.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		if(mode.equals("child")){
			require(!options.containsKey("reference") && Files.isDirectory(out), "Child uses an existing case directory and no reference.");
			child(out, new Spec(options.get("case")));
			return;
		}
		require(mode.equals("suite") && !options.containsKey("case"), "Use suite mode or an internal child case.");
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed suite needs reference=<successful baseline fixture outdir>.");
		if(reference!=null){require(read(reference.resolve("summary.txt")).startsWith("BBMERGE_HISTOGRAM_PASS phase=baseline cases=7 children=7 "), "Reference is not a successful seven-child baseline fixture directory.");}
		final BBMergeHistogramTest test=new BBMergeHistogramTest(root, out, reference, phase);
		for(final String name : NAMES){test.run(new Spec(name));}
		test.check(test.children==7 && test.changed==5 && test.controls==2, "Seven specified cases must retain five changed results and two unchanged controls.");
		test.check(test.sourceMatches==(phase.equals("baseline") ? 7 : 2), "Baseline-output match count differs from the distinct wrong/control oracles.");
		test.check(test.correctedMatches==(phase.equals("baseline") ? 2 : 7), "Corrected-output match count differs from the seven literal truth cases.");
		final String summary="BBMERGE_HISTOGRAM_"+(test.failures==0 ? "PASS" : "FAIL")+" phase="+phase+
			" cases=7 children="+test.children+" changed="+test.changed+" controls="+test.controls+" sourceMatches="+test.sourceMatches+
			" correctedMatches="+test.correctedMatches+" checks="+test.checks+" failures="+test.failures+"\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
		if(test.failures>0){throw new AssertionError(summary+"Inspect actual histogram bytes; do not fit the literal oracle to output.");}
	}

	/** Creates only new evidence and retains all child outcomes. */
	private BBMergeHistogramTest(final Path root_, final Path out_, final Path reference_, final String phase_) throws IOException{
		root=root_; out=out_; reference=reference_; phase=phase_;
		Files.createDirectory(out);
		write(out.resolve("cases.tsv"), "case\tchildExit\tbin\tmaximumInsert\tbaselineMatches\tcorrectedMatches\tphaseMatches\n");
	}

	/** Reads output only after Process.waitFor confirms complete child-JVM termination. */
	private void run(final Spec spec) throws Exception{
		if(spec.baseline.equals(spec.corrected)){controls++;}else{changed++;}
		final Path work=Files.createDirectory(out.resolve(spec.name));
		write(work.resolve("expected-baseline.tsv"), spec.baseline);
		write(work.resolve("expected-fixed.tsv"), spec.corrected);
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve("dev/bbmergehistogramtest.sh").toString(),
			"mode=child", "phase="+phase, "root="+root, "outdir="+work, "case="+spec.name));
		final StringBuilder saved=new StringBuilder();
		for(int i=0; i<command.size(); i++){saved.append(i).append('\t').append(command.get(i)).append('\n');}
		write(work.resolve("child.command.tsv"), saved.toString());
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile());
		builder.redirectOutput(work.resolve("child.stdout").toFile()); builder.redirectError(work.resolve("child.stderr").toFile());
		builder.environment().remove("_JAVA_OPTIONS"); builder.environment().remove("JDK_JAVA_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final Process process=builder.start(); children++;
		final int code=process.waitFor();
		write(work.resolve("child.exit"), code+"\n");
		check(code==0, spec.name+": child or its counted writer failed.");
		final Path histogram=work.resolve("histogram.tsv");
		check(Files.isRegularFile(histogram), spec.name+": completed child did not leave a histogram.");
		if(code!=0 || !Files.isRegularFile(histogram)){
			append(out.resolve("cases.tsv"), spec.name+"\t"+code+"\t"+spec.bin+"\t"+spec.maximum+"\tfalse\tfalse\tfalse\n");
			return;
		}
		final byte[] actual=Files.readAllBytes(histogram);
		final boolean source=Arrays.equals(actual, spec.baseline.getBytes(StandardCharsets.US_ASCII));
		final boolean corrected=Arrays.equals(actual, spec.corrected.getBytes(StandardCharsets.US_ASCII));
		if(source){sourceMatches++;} if(corrected){correctedMatches++;}
		final boolean matches=phase.equals("baseline") ? source : corrected;
		check(matches, spec.name+": complete header/row/index/count/newline bytes differ from the phase oracle.");
		check(read(work.resolve("child-result.tsv")).equals("capacity\t2000\nhistogramUnchanged\ttrue\ncountedWritesDrained\ttrue\n"), spec.name+": histogram input or counted-write completion evidence differs.");
		if(reference!=null){
			final Path before=reference.resolve(spec.name);
			for(final String file : new String[]{"setup.tsv", "inputs.tsv", "expected-baseline.tsv", "expected-fixed.tsv"}){
				check(Arrays.equals(Files.readAllBytes(work.resolve(file)), Files.readAllBytes(before.resolve(file))), spec.name+": case settings or literal oracles changed between phases: "+file);
			}
			if(spec.baseline.equals(spec.corrected)){check(Arrays.equals(actual, Files.readAllBytes(before.resolve("histogram.tsv"))), spec.name+": unchanged control bytes changed.");}
		}
		append(out.resolve("cases.tsv"), spec.name+"\t"+code+"\t"+spec.bin+"\t"+spec.maximum+"\t"+source+"\t"+corrected+"\t"+matches+"\n");
	}

	/** Changes test inputs only, never the real histogram array's capacity or writer code. */
	private static void child(final Path work, final Spec spec) throws Exception{
		final long[] histogram=(long[])field("histTotal").get(null);
		require(histogram.length==2000 && field("histlen").getInt(null)==2000, "Use the actual histogram capacity; a two-entry replacement is not a production truth case.");
		Arrays.fill(histogram, 0);
		for(int i=0; i<spec.indices.length; i++){histogram[spec.indices[i]]=spec.counts[i];}
		field("bin").setInt(null, spec.bin); field("insertMaxTotal").setInt(null, spec.maximum);
		field("showHistStats").setBoolean(null, false); field("NONZERO_ONLY").setBoolean(null, true);
		require(field("bin").getInt(null)==spec.bin && field("insertMaxTotal").getInt(null)==spec.maximum &&
			!field("showHistStats").getBoolean(null) && field("NONZERO_ONLY").getBoolean(null), "Actual writer settings differ from the isolated case.");
		write(work.resolve("setup.tsv"), "capacity\t2000\nbin\t"+spec.bin+"\nmaximumInsert\t"+spec.maximum+"\nshowHistStats\tfalse\nnonzeroOnly\ttrue\n");
		final StringBuilder inputs=new StringBuilder("index\tcount\n");
		for(int i=0; i<histogram.length; i++){if(histogram[i]!=0){inputs.append(i).append('\t').append(histogram[i]).append('\n');}}
		write(work.resolve("inputs.tsv"), inputs.toString());
		final long[] before=histogram.clone();
		require(ReadWrite.countActiveThreads()==0, "Fresh child must start with no counted writes; do not confuse prior work with this histogram.");
		BBMerge.writeHistogram(work.resolve("histogram.tsv").toString(), 0);
		ReadWrite.waitForWritingToFinish();
		require(ReadWrite.countActiveThreads()==0 && Arrays.equals(before, histogram), "Writer must drain its counted write and preserve histogram counts.");
		write(work.resolve("child-result.tsv"), "capacity\t2000\nhistogramUnchanged\ttrue\ncountedWritesDrained\ttrue\n");
		// Returning from main does not claim a join. Parent waits for child exit before reading.
	}

	/** Literal inputs and wrong/correct byte strings; no averaging implementation is reused. */
	private static final class Spec{
		Spec(final String name_){
			name=name_;
			require(name!=null && Arrays.asList(NAMES).contains(name), "Unknown bounded histogram case: "+name);
			maximum=name.startsWith("small-") ? 1 : 1999;
			switch(name){
				case "wide-count":
					bin=1; indices=new int[]{1}; counts=new long[]{2147483648L};
					baseline=HEADER; corrected=HEADER+"1\t2147483648\n"; break;
				case "round-overflow":
					bin=2; indices=new int[]{0}; counts=new long[]{2147483647L};
					baseline=HEADER; corrected=HEADER+"0\t1073741824\n"; break;
				case "zero-partial":
					bin=4001; indices=new int[0]; counts=new long[0];
					baseline=HEADER+"0\t2\n"; corrected=HEADER; break;
				case "tail-average":
					bin=3; indices=new int[]{1998}; counts=new long[]{2};
					baseline=HEADER+"1998\t2\n"; corrected=HEADER+"1998\t1\n"; break;
				case "small-single":
					bin=1; indices=new int[]{0, 1}; counts=new long[]{3, 5};
					baseline=HEADER+"0\t3\n1\t5\n"; corrected=baseline; break;
				case "small-double":
					bin=2; indices=new int[]{0, 1}; counts=new long[]{3, 5};
					baseline=HEADER+"0\t4\n"; corrected=baseline; break;
				case "long-ceiling":
					bin=2; indices=new int[]{0}; counts=new long[]{Long.MAX_VALUE};
					baseline=HEADER; corrected=HEADER+"0\t4611686018427387904\n"; break;
				default:throw new AssertionError("Unreachable case: "+name);
			}
			require(indices.length==counts.length && bin>0 && maximum<2000, "Literal cases must remain bounded in the actual histogram.");
		}
		final String name, baseline, corrected;
		final int bin, maximum;
		final int[] indices;
		final long[] counts;
	}

	private static Field field(final String name) throws Exception{final Field f=BBMerge.class.getDeclaredField(name); f.setAccessible(true); return f;}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}}
	private static void require(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private static Path path(final String value){require(value!=null && !value.isEmpty(), "Missing required path."); return Paths.get(value).toAbsolutePath().normalize();}
	private static void write(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);}
	private static String read(final Path path) throws IOException{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}

	private final Path root, out, reference;
	private final String phase;
	private int checks, failures, children, changed, controls, sourceMatches, correctedMatches;
	private static final String HEADER="#InsertSize\tCount\n";
	private static final String[] NAMES={"wide-count", "round-overflow", "zero-partial", "tail-average", "small-single", "small-double", "long-ceiling"};
}
