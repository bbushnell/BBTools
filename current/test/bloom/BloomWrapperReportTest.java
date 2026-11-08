package test.bloom;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import bloom.BloomFilterWrapper;

/** Checks finite match reports and exact routing for empty and single-read inputs.
 * Runs the real wrapper with a tiny reference table and saves diagnostics for each case.
 * Use a fresh JVM because the wrapper changes shared parsing and counting settings.
 * This is a reporting/routing test, not an exceptional-cleanup or performance test.
 * @author Brian Bushnell, Shinobu */
public final class BloomWrapperReportTest{

	/** Takes one new scratch directory and records all three case outcomes. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new scratch directory>");
		}
		final Path root=Files.createDirectory(Paths.get(args[0].substring(7)));
		int failures=0;
		failures+=exercise(root, "empty", EMPTY, EMPTY, EMPTY, 0, 0, "0.00%");
		failures+=exercise(root, "unmatched", UNMATCHED, UNMATCHED, EMPTY, 0, 0, "0.00%");
		failures+=exercise(root, "matched", MATCHED, EMPTY, MATCHED, 1, 8, "100.00%");
		if(failures>0){throw new AssertionError("Bloom wrapper reports failed "+failures+" of 3 cases.");}
		System.out.println("BLOOM_WRAPPER_REPORT_TEST_PASS cases=3");
	}

	/** Checks literal output bytes, unchanged inputs and both matched-count report rows. */
	private static int exercise(final Path root, final String name, final byte[] input,
			final byte[] expectedUnmatched, final byte[] expectedMatched, final long matchedReads,
			final long matchedBases, final String percent) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name));
		final Path in=dir.resolve("input.fq"), ref=dir.resolve("reference.fa");
		final Path out=dir.resolve("unmatched.fq"), outm=dir.resolve("matched.fq");
		Files.write(in, input);
		Files.write(ref, REFERENCE);
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream prior=System.err;
		boolean completed=false;
		try(PrintStream capture=new PrintStream(buffer, true, "UTF-8")){
			System.setErr(capture);
			try{
				BloomFilterWrapper.main(new String[]{"in="+in, "ref="+ref, "out="+out, "outm="+outm,
						"ow=t", "k=3", "minhits=1", "bits=4", "hashes=1", "cells=10007", "seed=1",
						"minprob=0", "t=1", "interleaved=f", "pigz=f", "unpigz=f", "qin=33", "qout=33",
						"merge=f", "ecco=f"});
				completed=true;
			}catch(RuntimeException | AssertionError failure){failure.printStackTrace(capture);}
		}finally{
			System.setErr(prior);
			Files.write(dir.resolve("stderr.txt"), buffer.toByteArray());
		}
		final String report=buffer.toString("UTF-8");
		final boolean finite=!report.contains("NaN") && !report.contains("Infinity");
		final boolean counts=hasRow(report, "Reads Matched:", matchedReads, percent) &&
				hasRow(report, "Bases Matched:", matchedBases, percent);
		final boolean bytes=equals(out, expectedUnmatched) && equals(outm, expectedMatched) &&
				equals(in, input) && equals(ref, REFERENCE);
		final boolean pass=completed && finite && counts && bytes;
		System.out.println("case="+name+" completed="+completed+" finite="+finite+
				" expectedCounts="+counts+" exactBytes="+bytes+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Requires exactly one count and percentage after the requested report label. */
	private static boolean hasRow(final String report, final String label, final long count, final String percent){
		for(String line : report.split("\\r?\\n")){
			if(line.startsWith(label)){
				final String[] fields=line.substring(label.length()).trim().split("\\s+");
				return fields.length==2 && fields[0].equals(Long.toString(count)) && fields[1].equals(percent);
			}
		}
		return false;
	}

	/** Missing files fail rather than being treated as valid empty output. */
	private static boolean equals(final Path path, final byte[] expected) throws Exception{
		return Files.isRegularFile(path) && Arrays.equals(expected, Files.readAllBytes(path));
	}

	private static final byte[] EMPTY=new byte[0];
	private static final byte[] MATCHED="@matched\nCCCCCCCC\n+\nIIIIIIII\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] UNMATCHED="@unmatched\nAAAAAAAA\n+\nIIIIIIII\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] REFERENCE=">reference\nCCCCCCCCCCCCCCCC\n".getBytes(StandardCharsets.US_ASCII);
}
