package test.bloom;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import bloom.PolyFilter;
import shared.Timer;

/** Checks PolyFilter total removed counts against polymer-specific removal counts.
 * Uses real temporary FASTQ files and the full PolyFilter processing loop with
 * tiny inputs. The polymer-specific last-run fields are read reflectively so
 * their absence is reported explicitly on the current baseline.
 * Run in a fresh JVM because construction changes process-wide I/O settings.
 * @author Brian Bushnell, Jean */
public final class PolyFilterAccountingTest{

	/** Takes one new scratch-directory path and reports all cases before failing. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one new scratch-directory path.");}
		final Path root=Files.createDirectory(Paths.get(args[0]));
		int failures=0;
		failures+=entropyOnly(root);
		failures+=truePolymer(root);
		failures+=pairedMixedCause(root);
		failures+=trimAccounting(root);
		failures+=trimTooShort(root);
		failures+=shortNoTrim(root);
		failures+=pairedTrimTooShort(root);
		failures+=pairedOriginallyShort(root);
		if(failures!=0){throw new AssertionError("PolyFilter accounting failed "+failures+" of 8 cases.");}
		System.out.println("POLYFILTER_ACCOUNTING_TEST_PASS cases=8");
	}

	/** A high-quality non-G/C low-entropy read should not count as a polymer discard. */
	private static int entropyOnly(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("entropy-only"));
		final Path in=dir.resolve("single.fq");
		final Path out=dir.resolve("kept.fq"), bad=dir.resolve("bad.fq");
		Files.write(in, FASTQ_ENTROPY);
		final Result result=run(dir, in, null, out, null, bad, null, new String[]{
				"polymers=GC", "trimpolymers=G", "trimleft=0", "trimright=0",
				"minlen=1", "minpolymer=999", "minpolymer2=999", "entropy2=1",
				"quality2=0", "ldf2=2", "merge=f"});
		final boolean totals=result.readsOut==0 && result.readsRemoved==1 &&
				result.basesOut==0 && result.basesRemoved==12;
		final boolean routing=empty(out) && equals(bad, FASTQ_ENTROPY);
		final boolean polymer=result.polymerFieldsSupported &&
				result.polymerReadsRemoved==0 && result.polymerBasesRemoved==0;
		final boolean pass=totals && routing && polymer;
		print("entropy-only", result, 1, 12, 0, 0, totals, routing, polymer, pass,
				"merge=f trim=0 minlen=1");
		return pass ? 0 : 1;
	}

	/** A true G polymer reject should count as both total and polymer-specific removal. */
	private static int truePolymer(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("true-polymer"));
		final Path in=dir.resolve("single.fq");
		final Path out=dir.resolve("kept.fq"), bad=dir.resolve("bad.fq");
		Files.write(in, FASTQ_POLY_G);
		final Result result=run(dir, in, null, out, null, bad, null, new String[]{
				"polymers=G", "trimpolymers=G", "trimleft=0", "trimright=0",
				"minlen=1", "minpolymer=999", "minpolymer2=6", "entropy2=0",
				"quality2=0", "ldf2=2", "merge=f"});
		final boolean totals=result.readsOut==0 && result.readsRemoved==1 &&
				result.basesOut==0 && result.basesRemoved==12;
		final boolean routing=empty(out) && equals(bad, FASTQ_POLY_G);
		final boolean polymer=result.polymerFieldsSupported &&
				result.polymerReadsRemoved==1 && result.polymerBasesRemoved==12;
		final boolean pass=totals && routing && polymer;
		print("true-polymer", result, 1, 12, 1, 12, totals, routing, polymer, pass,
				"merge=f trim=0 minlen=1");
		return pass ? 0 : 1;
	}

	/** A pair discarded by mixed causes should count the whole pair as polymer-caused. */
	private static int pairedMixedCause(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("paired-mixed-cause"));
		final Path in1=dir.resolve("r1.fq"), in2=dir.resolve("r2.fq");
		final Path out1=dir.resolve("kept1.fq"), out2=dir.resolve("kept2.fq");
		final Path bad1=dir.resolve("bad1.fq"), bad2=dir.resolve("bad2.fq");
		Files.write(in1, FASTQ_PAIR_ENTROPY);
		Files.write(in2, FASTQ_PAIR_POLY_G);
		final Result result=run(dir, in1, in2, out1, out2, bad1, bad2, new String[]{
				"polymers=G", "trimpolymers=G", "trimleft=0", "trimright=0",
				"minlen=1", "minpolymer=999", "minpolymer2=6", "entropy2=1",
				"quality2=0", "ldf2=2", "merge=f"});
		final boolean totals=result.readsOut==0 && result.readsRemoved==2 &&
				result.basesOut==0 && result.basesRemoved==24;
		final boolean routing=empty(out1) && empty(out2) &&
				equals(bad1, FASTQ_PAIR_ENTROPY) && equals(bad2, FASTQ_PAIR_POLY_G);
		final boolean polymer=result.polymerFieldsSupported &&
				result.polymerReadsRemoved==2 && result.polymerBasesRemoved==24;
		final boolean pass=totals && routing && polymer;
		print("paired-mixed-cause", result, 2, 24, 2, 24, totals, routing, polymer, pass,
				"merge=f trim=0 minlen=1");
		return pass ? 0 : 1;
	}

	/** Polymer bases trimmed from a kept read should count as removed polymer bases. */
	private static int trimAccounting(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("trim-accounting"));
		final Path in=dir.resolve("single.fq");
		final Path out=dir.resolve("kept.fq"), bad=dir.resolve("bad.fq");
		Files.write(in, FASTQ_TRIM_G);
		final Result result=run(dir, in, null, out, null, bad, null, new String[]{
				"polymers=G", "trimpolymers=G", "trimleft=0", "trimright=6",
				"minlen=6", "minpolymer=999", "minpolymer2=999", "entropy2=0",
				"quality2=0", "ldf2=2", "merge=f"});
		final boolean totals=result.readsOut==1 && result.readsRemoved==0 &&
				result.basesOut==6 && result.basesRemoved==6;
		final boolean routing=equals(out, FASTQ_TRIM_G_TRIMMED) && empty(bad);
		final boolean polymer=result.polymerFieldsSupported &&
				result.polymerReadsRemoved==0 && result.polymerBasesRemoved==6;
		final boolean pass=totals && routing && polymer;
		print("trim-accounting", result, 0, 6, 0, 6, totals, routing, polymer, pass,
				"merge=f trimright=6 minlen=6");
		return pass ? 0 : 1;
	}

	/** A polymer trim that makes a long-enough read too short should count the whole read. */
	private static int trimTooShort(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("trim-too-short"));
		final Path in=dir.resolve("single.fq");
		final Path out=dir.resolve("kept.fq"), bad=dir.resolve("bad.fq");
		Files.write(in, FASTQ_TRIM_SHORT);
		final Result result=run(dir, in, null, out, null, bad, null, new String[]{
				"polymers=G", "trimpolymers=G", "trimleft=0", "trimright=6",
				"minlen=6", "minpolymer=999", "minpolymer2=999", "entropy2=0",
				"quality2=0", "ldf2=2", "merge=f"});
		final boolean totals=result.readsOut==0 && result.readsRemoved==1 &&
				result.basesOut==0 && result.basesRemoved==11;
		final boolean routing=empty(out) && equals(bad, FASTQ_TRIM_SHORT_TRIMMED);
		final boolean polymer=result.polymerFieldsSupported &&
				result.polymerReadsRemoved==1 && result.polymerBasesRemoved==11;
		final boolean pass=totals && routing && polymer;
		print("trim-too-short", result, 1, 11, 1, 11, totals, routing, polymer, pass,
				"merge=f trimright=6 minlen=6 maxnonpoly=0");
		return pass ? 0 : 1;
	}

	/** An originally short read rejected by minlen should not count as a polymer discard. */
	private static int shortNoTrim(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("short-no-trim"));
		final Path in=dir.resolve("single.fq");
		final Path out=dir.resolve("kept.fq"), bad=dir.resolve("bad.fq");
		Files.write(in, FASTQ_SHORT);
		final Result result=run(dir, in, null, out, null, bad, null, new String[]{
				"polymers=G", "trimpolymers=G", "trimleft=0", "trimright=6",
				"minlen=6", "minpolymer=999", "minpolymer2=999", "entropy2=0",
				"quality2=0", "ldf2=2", "merge=f"});
		final boolean totals=result.readsOut==0 && result.readsRemoved==1 &&
				result.basesOut==0 && result.basesRemoved==5;
		final boolean routing=empty(out) && equals(bad, FASTQ_SHORT);
		final boolean polymer=result.polymerFieldsSupported &&
				result.polymerReadsRemoved==0 && result.polymerBasesRemoved==0;
		final boolean pass=totals && routing && polymer;
		print("short-no-trim", result, 1, 5, 0, 0, totals, routing, polymer, pass,
				"merge=f trimright=6 minlen=6 no-terminal-polymer");
		return pass ? 0 : 1;
	}

	/** A trim-caused short mate should count the whole originally long-enough pair. */
	private static int pairedTrimTooShort(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("paired-trim-too-short"));
		final Path in1=dir.resolve("r1.fq"), in2=dir.resolve("r2.fq");
		final Path out1=dir.resolve("kept1.fq"), out2=dir.resolve("kept2.fq");
		final Path bad1=dir.resolve("bad1.fq"), bad2=dir.resolve("bad2.fq");
		Files.write(in1, FASTQ_PAIR_TRIM_SHORT_1);
		Files.write(in2, FASTQ_PAIR_TRIM_KEEP);
		final Result result=run(dir, in1, in2, out1, out2, bad1, bad2, new String[]{
				"polymers=G", "trimpolymers=G", "trimleft=0", "trimright=6",
				"minlen=6", "minpolymer=999", "minpolymer2=999", "entropy2=0",
				"quality2=0", "ldf2=2", "merge=f"});
		final boolean totals=result.readsOut==0 && result.readsRemoved==2 &&
				result.basesOut==0 && result.basesRemoved==17;
		final boolean routing=empty(out1) && empty(out2) &&
				equals(bad1, FASTQ_PAIR_TRIM_SHORT_1_TRIMMED) &&
				equals(bad2, FASTQ_PAIR_TRIM_KEEP);
		final boolean polymer=result.polymerFieldsSupported &&
				result.polymerReadsRemoved==2 && result.polymerBasesRemoved==17;
		final boolean pass=totals && routing && polymer;
		print("paired-trim-too-short", result, 2, 17, 2, 17, totals, routing, polymer,
				pass, "merge=f trimright=6 minlen=6 both-originally-long");
		return pass ? 0 : 1;
	}

	/** A paired reject with one originally short mate should not count as polymer-caused. */
	private static int pairedOriginallyShort(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("paired-originally-short"));
		final Path in1=dir.resolve("r1.fq"), in2=dir.resolve("r2.fq");
		final Path out1=dir.resolve("kept1.fq"), out2=dir.resolve("kept2.fq");
		final Path bad1=dir.resolve("bad1.fq"), bad2=dir.resolve("bad2.fq");
		Files.write(in1, FASTQ_PAIR_SHORT);
		Files.write(in2, FASTQ_PAIR_TRIM_SHORT_2);
		final Result result=run(dir, in1, in2, out1, out2, bad1, bad2, new String[]{
				"polymers=G", "trimpolymers=G", "trimleft=0", "trimright=6",
				"minlen=6", "minpolymer=999", "minpolymer2=999", "entropy2=0",
				"quality2=0", "ldf2=2", "merge=f"});
		final boolean totals=result.readsOut==0 && result.readsRemoved==2 &&
				result.basesOut==0 && result.basesRemoved==16;
		final boolean routing=empty(out1) && empty(out2) && equals(bad1, FASTQ_PAIR_SHORT) &&
				equals(bad2, FASTQ_PAIR_TRIM_SHORT_2_TRIMMED);
		final boolean polymer=result.polymerFieldsSupported &&
				result.polymerReadsRemoved==0 && result.polymerBasesRemoved==6;
		final boolean pass=totals && routing && polymer;
		print("paired-originally-short", result, 2, 16, 0, 6, totals, routing, polymer,
				pass, "merge=f trimright=6 minlen=6 one-originally-short");
		return pass ? 0 : 1;
	}

	/** Runs PolyFilter with common guard flags and captures the last-run counters. */
	private static Result run(final Path dir, final Path in1, final Path in2, final Path out1,
			final Path out2, final Path bad1, final Path bad2, final String[] caseFlags)
			throws Exception{
		final ArrayList<String> flags=new ArrayList<String>();
		for(String flag : new String[]{
				"in="+in1, "out="+out1, "outb="+bad1, "ow=t", "ordered=t",
				"entropy=0", "quality=0", "ldf=2", "minprob=0", "kpoly=0",
				"maxnonpoly=0", "t=1", "interleaved=f", "pigz=f", "unpigz=f",
				"qin=33"}){
			flags.add(flag);
		}
		flags.addAll(Arrays.asList(caseFlags));
		if(in2!=null){flags.add("in2="+in2);}
		if(out2!=null){flags.add("out2="+out2);}
		if(bad2!=null){flags.add("outb2="+bad2);}
		Files.write(dir.resolve("flags.txt"), flags, StandardCharsets.UTF_8);

		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream old=System.err;
		try{
			System.setErr(new PrintStream(buffer, true, "UTF-8"));
			new PolyFilter(flags.toArray(new String[flags.size()])).process(new Timer());
		}finally{
			System.setErr(old);
			Files.write(dir.resolve("stderr.txt"),
					buffer.toString("UTF-8").getBytes(StandardCharsets.UTF_8));
		}
		return new Result(buffer.toString("UTF-8"));
	}

	/** Prints one result with totals and polymer-specific support side by side. */
	private static void print(final String name, final Result result, final long expectedReadsRemoved,
			final long expectedBasesRemoved, final long expectedPolymerReads,
			final long expectedPolymerBases, final boolean totals, final boolean routing,
			final boolean polymer, final boolean pass, final String branch){
		System.out.println("case="+name+" expectedRemoved="+expectedReadsRemoved+
				"/"+expectedBasesRemoved+" actualRemoved="+result.readsRemoved+
				"/"+result.basesRemoved+" expectedPolymer="+expectedPolymerReads+
				"/"+expectedPolymerBases+" actualPolymer="+result.polymerReadsRemoved+
				"/"+result.polymerBasesRemoved+" polymerSupported="+
				result.polymerFieldsSupported+" missingPolymerField="+
				result.missingPolymerField+" totals="+totals+" routing="+routing+
				" polymer="+polymer+" pass="+pass+" "+branch);
	}

	/** True if the file is absent or has no bytes. */
	private static boolean empty(final Path path) throws Exception{
		return !Files.exists(path) || Files.size(path)==0;
	}

	/** Compares a real output file to the literal expected FASTQ bytes. */
	private static boolean equals(final Path path, final byte[] expected) throws Exception{
		return Files.exists(path) && Arrays.equals(Files.readAllBytes(path), expected);
	}

	/** Reads a public static long from PolyFilter by field name. */
	private static long staticLong(final String name) throws Exception{
		final Field field=PolyFilter.class.getField(name);
		return field.getLong(null);
	}

	/** PolyFilter counters captured immediately after one processing run. */
	private static final class Result{
		Result(final String stderr_) throws Exception{
			stderr=stderr_;
			readsOut=PolyFilter.lastReadsOut;
			readsRemoved=PolyFilter.lastReadsRemoved;
			basesOut=PolyFilter.lastBasesOut;
			basesRemoved=PolyFilter.lastBasesRemoved;
			boolean supported=true;
			String missing="";
			long polymerReads=-1;
			long polymerBases=-1;
			try{
				polymerReads=staticLong("lastReadsPolymerRemoved");
				polymerBases=staticLong("lastBasesPolymerRemoved");
			}catch(NoSuchFieldException e){
				supported=false;
				missing=e.getMessage();
			}
			polymerFieldsSupported=supported;
			missingPolymerField=missing;
			polymerReadsRemoved=polymerReads;
			polymerBasesRemoved=polymerBases;
		}
		final String stderr;
		final long readsOut, readsRemoved, basesOut, basesRemoved;
		final boolean polymerFieldsSupported;
		final String missingPolymerField;
		final long polymerReadsRemoved, polymerBasesRemoved;
	}

	private static final byte[] FASTQ_ENTROPY=
			("@entropy\nATATATATATAT\n+\nIIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_POLY_G=
			("@poly-g\nGGGGGGGGGGGG\n+\nIIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_PAIR_ENTROPY=
			("@pair/1\nATATATATATAT\n+\nIIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_PAIR_POLY_G=
			("@pair/2\nGGGGGGGGGGGG\n+\nIIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_TRIM_G=
			("@trim\nAAAAAAGGGGGG\n+\nIIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_TRIM_G_TRIMMED=
			("@trim\nAAAAAA\n+\nIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_TRIM_SHORT=
			("@trim-short\nAAAAAGGGGGG\n+\nIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_TRIM_SHORT_TRIMMED=
			("@trim-short\nAAAAA\n+\nIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_SHORT=
			("@short\nAAAAA\n+\nIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_PAIR_TRIM_SHORT_1=
			("@pair-trim/1\nAAAAAGGGGGG\n+\nIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_PAIR_TRIM_SHORT_1_TRIMMED=
			("@pair-trim/1\nAAAAA\n+\nIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_PAIR_TRIM_KEEP=
			("@pair-trim/2\nAAAAAA\n+\nIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_PAIR_SHORT=
			("@pair-short/1\nAAAAA\n+\nIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_PAIR_TRIM_SHORT_2=
			("@pair-short/2\nAAAAAGGGGGG\n+\nIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_PAIR_TRIM_SHORT_2_TRIMMED=
			("@pair-short/2\nAAAAA\n+\nIIIII\n").getBytes(StandardCharsets.US_ASCII);
}
