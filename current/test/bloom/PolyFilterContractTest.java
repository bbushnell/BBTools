package test.bloom;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import bloom.PolyFilter;
import shared.Timer;

/** Checks finite PolyFilter empty-input reporting and post-trim minlen boundaries.
 * Uses real temporary FASTQ files and the full PolyFilter processing loop with
 * tiny inputs. Does not cover RQCFilter3 PolyG cause accounting.
 * Run in a fresh JVM because construction changes process-wide I/O settings.
 * @author Brian Bushnell, Jean */
public final class PolyFilterContractTest{

	/** Takes one new scratch-directory path and reports all cases before failing. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one new scratch-directory path.");}
		final Path root=Files.createDirectory(Paths.get(args[0]));
		int failures=0;
		failures+=emptyInput(root);
		failures+=singleMinLength(root);
		failures+=singleBelowMinLength(root);
		failures+=pairedMinLength(root);
		if(failures!=0){throw new AssertionError("PolyFilter contracts failed "+failures+" of 4 cases.");}
		System.out.println("POLYFILTER_CONTRACT_TEST_PASS cases=4");
	}

	/** Empty input must keep every percentage row finite at zero. */
	private static int emptyInput(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("empty"));
		final Path in=dir.resolve("empty.fq");
		final Path out=dir.resolve("kept.fq"), bad=dir.resolve("bad.fq");
		Files.write(in, new byte[0]);
		final Result result=run(dir, in, null, out, null, bad, null, true, 6);
		final boolean report=hasZeroReportLine(result.stderr, "Merged Reads:") &&
				hasZeroReportLine(result.stderr, "Poly-G Reads:") &&
				hasZeroReportLine(result.stderr, "Poly-C Reads:") &&
				hasZeroReportLine(result.stderr, "Low Entropy Reads:") &&
				hasZeroReportLine(result.stderr, "Low Quality Reads:") &&
				hasZeroReportLine(result.stderr, "Trimmed Reads:") &&
				hasZeroReportLine(result.stderr, "Total Discards:");
		final boolean finite=finite(result.stderr);
		final boolean pass=finite && report && result.readsOut==0 && result.readsRemoved==0 &&
				result.basesOut==0 && result.basesRemoved==0 && empty(out) && empty(bad);
		System.out.println("case=empty readsOut="+result.readsOut+" readsRemoved="+result.readsRemoved+
				" basesOut="+result.basesOut+" basesRemoved="+result.basesRemoved+
				" finiteReport="+finite+" zeroReportLines="+report+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** A single read trimmed exactly to minlen should be kept. */
	private static int singleMinLength(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("single-minlen"));
		final Path in=dir.resolve("single.fq");
		final Path out=dir.resolve("kept.fq"), bad=dir.resolve("bad.fq");
		Files.write(in, FASTQ_SINGLE);
		final Result result=run(dir, in, null, out, null, bad, null, false, 6);
		final boolean pass=result.readsOut==1 && result.readsRemoved==0 &&
				result.basesOut==6 && result.basesRemoved==6 &&
				equals(out, FASTQ_SINGLE_TRIMMED) && empty(bad);
		System.out.println("case=single-minlen expectedKept=1 readsOut="+result.readsOut+
				" readsRemoved="+result.readsRemoved+" basesOut="+result.basesOut+
				" basesRemoved="+result.basesRemoved+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** A single read trimmed below minlen should still be rejected. */
	private static int singleBelowMinLength(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("single-below-minlen"));
		final Path in=dir.resolve("single.fq");
		final Path out=dir.resolve("kept.fq"), bad=dir.resolve("bad.fq");
		Files.write(in, FASTQ_BELOW);
		final Result result=run(dir, in, null, out, null, bad, null, false, 6);
		final boolean pass=result.readsOut==0 && result.readsRemoved==1 &&
				result.basesOut==0 && result.basesRemoved==11 &&
				empty(out) && equals(bad, FASTQ_BELOW_TRIMMED);
		System.out.println("case=single-below-minlen expectedRejected=1 readsOut="+result.readsOut+
				" readsRemoved="+result.readsRemoved+" basesOut="+result.basesOut+
				" basesRemoved="+result.basesRemoved+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** A pair with R1 above minlen and R2 exactly at minlen should be kept. */
	private static int pairedMinLength(final Path root) throws Exception{
		final Path dir=Files.createDirectory(root.resolve("paired-minlen"));
		final Path in1=dir.resolve("r1.fq"), in2=dir.resolve("r2.fq");
		final Path out1=dir.resolve("kept1.fq"), out2=dir.resolve("kept2.fq");
		final Path bad1=dir.resolve("bad1.fq"), bad2=dir.resolve("bad2.fq");
		Files.write(in1, FASTQ_R1);
		Files.write(in2, FASTQ_R2);
		final Result result=run(dir, in1, in2, out1, out2, bad1, bad2, false, 6);
		final boolean pass=result.readsOut==2 && result.readsRemoved==0 &&
				result.basesOut==13 && result.basesRemoved==12 &&
				equals(out1, FASTQ_R1_TRIMMED) && equals(out2, FASTQ_R2_TRIMMED) &&
				empty(bad1) && empty(bad2);
		System.out.println("case=paired-minlen expectedKept=2 readsOut="+result.readsOut+
				" readsRemoved="+result.readsRemoved+" basesOut="+result.basesOut+
				" basesRemoved="+result.basesRemoved+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Runs PolyFilter and captures stderr so the finite report can be inspected. */
	private static Result run(final Path dir, final Path in1, final Path in2, final Path out1, final Path out2,
			final Path bad1, final Path bad2, final boolean emptyReport, final int minlen) throws Exception{
		final ArrayList<String> flags=new ArrayList<String>();
		for(String flag : new String[]{"in="+in1, "out="+out1, "outb="+bad1, "ow=t",
				"polymers="+(emptyReport ? "GC" : "G"), "trimpolymers=G", "trimleft=0",
				"trimright=6", "minlen="+minlen, "ldf=2", "ldf2=2",
				"entropy2="+(emptyReport ? "0.01" : "0"),
				"quality2="+(emptyReport ? "1" : "0"), "minprob=0", "t=1",
				"interleaved=f", "pigz=f", "unpigz=f", "qin=33"}){
			flags.add(flag);
		}
		if(in2!=null){flags.add("in2="+in2);}
		if(out2!=null){flags.add("out2="+out2);}
		if(bad2!=null){flags.add("outb2="+bad2);}
		if(emptyReport){flags.add("merge=t");}
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream old=System.err;
		try{
			System.setErr(new PrintStream(buffer, true, "UTF-8"));
			new PolyFilter(flags.toArray(new String[flags.size()])).process(new Timer());
		}finally{
			System.setErr(old);
			Files.write(dir.resolve("stderr.txt"), buffer.toString("UTF-8").getBytes(StandardCharsets.UTF_8));
		}
		final String stderr=buffer.toString("UTF-8");
		return new Result(stderr);
	}

	/** True when a captured report contains the exact fixed zero-count percentage row. */
	private static boolean hasZeroReportLine(final String stderr, final String label){
		final StringBuilder expected=new StringBuilder(label);
		while(expected.length()<20){expected.append(' ');}
		expected.append("       0 \t0.0000%");
		for(String line : stderr.split("\\r?\\n")){
			if(line.equals(expected.toString())){return true;}
		}
		return false;
	}

	/** True when all captured percentages and rates are finite. */
	private static boolean finite(final String stderr){
		return !stderr.contains("NaN") && !stderr.contains("Infinity");
	}

	/** True if the file is absent or has no bytes. */
	private static boolean empty(final Path path) throws Exception{
		return !Files.exists(path) || Files.size(path)==0;
	}

	/** Compares a real output file to the literal expected FASTQ bytes. */
	private static boolean equals(final Path path, final byte[] expected) throws Exception{
		return Files.exists(path) && java.util.Arrays.equals(Files.readAllBytes(path), expected);
	}

	/** PolyFilter counters captured immediately after one processing run. */
	private static final class Result{
		Result(final String stderr_){
			stderr=stderr_;
			readsOut=PolyFilter.lastReadsOut;
			readsRemoved=PolyFilter.lastReadsRemoved;
			basesOut=PolyFilter.lastBasesOut;
			basesRemoved=PolyFilter.lastBasesRemoved;
		}
		final String stderr;
		final long readsOut, readsRemoved, basesOut, basesRemoved;
	}

	private static final byte[] FASTQ_SINGLE=
			("@single\nAAAAAAGGGGGG\n+\nIIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_SINGLE_TRIMMED=
			("@single\nAAAAAA\n+\nIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_BELOW=
			("@below\nAAAAAGGGGGG\n+\nIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_BELOW_TRIMMED=
			("@below\nAAAAA\n+\nIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_R1=
			("@pair/1\nAAAAAAAGGGGGG\n+\nIIIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_R1_TRIMMED=
			("@pair/1\nAAAAAAA\n+\nIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_R2=
			("@pair/2\nAAAAAAGGGGGG\n+\nIIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_R2_TRIMMED=
			("@pair/2\nAAAAAA\n+\nIIIIII\n").getBytes(StandardCharsets.US_ASCII);
}
