package test.bloom;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import bloom.BloomFilterCorrectorWrapper;

/** Finite correction/marking reports for empty input and fully supported reads (BLM216).
 * Invokes the real wrapper with a tiny reference table, retaining each invocation's
 * flags and diagnostics. Exact read/quality bytes and original inputs must survive.
 * Run in a fresh JVM because constructors change process-wide parser/counting state.
 * This checks normal reports, not exceptional cleanup, correction accuracy or speed.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorReportTest{

	/** Captures four empty/nonempty correction/marking cases in a new directory. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path root=Files.createDirectory(Paths.get(args[0].substring(7)));
		int failures=0;
		failures+=exercise(root, "correction-empty", EMPTY, false);
		failures+=exercise(root, "correction-supported", SUPPORTED, false);
		failures+=exercise(root, "marking-empty", EMPTY, true);
		failures+=exercise(root, "marking-supported", SUPPORTED, true);
		if(failures>0){throw new AssertionError("Correction wrapper reports: "+failures+" failures of 4 cases.");}
		System.out.println("BLOOM_CORRECTOR_REPORT_TEST_PASS cases=4");
	}

	/** Checks report labels/counts/percentages and literal retained/rejected output bytes. */
	private static int exercise(final Path root, final String name, final byte[] input,
			final boolean marking) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name));
		final Path in=dir.resolve("input.fq"), ref=dir.resolve("reference.fa");
		final Path out=dir.resolve("output.fq"), rejected=dir.resolve("rejected.fq");
		Files.write(in, input);
		Files.write(ref, REFERENCE);
		final String[] flags={"in="+in, "ref="+ref, "out="+out, "outbad="+rejected,
				"ow=t", "k=15", "ksmall=15", "bits=4", "hashes=1", "cells=10007", "seed=1", "minprob=0",
				"t=1", "interleaved=f", "pigz=f", "unpigz=f", "qin=33", "qout=33", "loglog=f",
				"ecc="+!marking, "markerrors="+marking, "merge="+!marking, "ecco=f", "mincount=0", "tossjunk=f"};
		Files.write(dir.resolve("wrapper-args.txt"), Arrays.asList(flags), StandardCharsets.UTF_8);
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream prior=System.err;
		boolean completed=false;
		try(PrintStream capture=new PrintStream(buffer, true, "UTF-8")){
			System.setErr(capture);
			try{
				BloomFilterCorrectorWrapper.main(flags);
				completed=true;
			}catch(RuntimeException | AssertionError failure){failure.printStackTrace(capture);}
		}finally{
			System.setErr(prior);
			Files.write(dir.resolve("stderr.txt"), buffer.toByteArray());
		}
		final String report=buffer.toString("UTF-8");
		final boolean finite=!report.contains("NaN") && !report.contains("Infinity");
		final boolean counts;
		if(marking){
			counts=hasRow(report, "Reads marked:", "0 (0.00%)") &&
					hasRow(report, "Bases marked:", "0 (0.00%)");
		}else{
			counts=hasRow(report, "Errors detected:", "0") &&
					hasRow(report, "Errors corrected:", "0 (0 reassemble, 0 overlap)") &&
					hasRow(report, "Reads merged:", "0 (0.00%)") &&
					hasRow(report, "Reads with errors detected:", "0 (0.00%)") &&
					hasRow(report, "Reads fully corrected:", "0 (0.00% of detected)") &&
					hasRow(report, "Reads partly corrected:", "0 (0.00% of detected)") &&
					hasRow(report, "Rollbacks:", "0 (0.00% of detected)");
		}
		final boolean bytes=equals(out, input) && equals(rejected, EMPTY) && equals(in, input) && equals(ref, REFERENCE);
		final boolean pass=completed && finite && counts && bytes;
		System.out.println("case="+name+" completed="+completed+" finite="+finite+
				" expectedCounts="+counts+" exactBytes="+bytes+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Requires exactly one matching report row; ignores only spacing, not labels or values. */
	private static boolean hasRow(final String report, final String label, final String expected){
		int matches=0;
		for(String line : report.split("\\r?\\n")){
			if(line.startsWith(label)){
				if(!line.substring(label.length()).trim().replaceAll("\\s+", " ").equals(expected)){return false;}
				matches++;
			}
		}
		return matches==1;
	}

	/** Missing output is a failure, even when the expected output is empty. */
	private static boolean equals(final Path path, final byte[] expected) throws Exception{
		return Files.isRegularFile(path) && Arrays.equals(expected, Files.readAllBytes(path));
	}

	private static final byte[] EMPTY=new byte[0];
	private static final byte[] SUPPORTED="@supported\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] REFERENCE=">reference\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n".getBytes(StandardCharsets.US_ASCII);
}
