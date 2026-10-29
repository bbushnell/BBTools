package test.tax;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.zip.GZIPOutputStream;

/** Finite child-process oracle for AnalyzeAccession worker failures.
 * @author Jean */
public final class AnalyzeAccessionWorkerFailureTest{

	/** Accession2taxid input with one header row and one valid accession row. */
	private static final String HEALTHY=
			"accession\taccession.version\ttaxid\tgi\n"+
			"ABC123\tABC123.1\t7\t0\n";
	/** Exact pattern table produced by {@link #HEALTHY}. */
	private static final String EXPECTED=
			"#Pattern\tCount\tCombos\tBits\n"+
			"LLLDDD\t1\t17576000\t24.07\n";
	/** Exact pattern table produced by two copies of {@link #HEALTHY}. */
	private static final String EXPECTED_DOUBLE=
			"#Pattern\tCount\tCombos\tBits\n"+
			"LLLDDD\t2\t17576000\t24.07\n";

	/** Runs baseline or fixed expectations for failed accession workers.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture construction, child launch or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("analyze-accession-worker-");
		try{
			for(boolean assertions : new boolean[] {true, false}){
				for(boolean perFile : new boolean[] {false, true}){
					checkHealthy(root, assertions, perFile);
					checkPoison(root, assertions, perFile, baseline);
					checkTruncatedGzip(root, assertions, perFile, baseline);
				}
				checkPerFileTwoGzip(root, assertions);
				checkPerFileTruncatedGzip(root, assertions, baseline);
			}
		}finally{
			deleteTempTree(root);
		}
		System.out.println("ANALYZE_ACCESSION_WORKER_FAILURE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private AnalyzeAccessionWorkerFailureTest(){}

	/** Runs one all-healthy child and checks exact pattern output.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param perFile Value for AnalyzeAccession's perfile flag
	 * @throws Exception If the child launch or file IO fails */
	private static void checkHealthy(Path root, boolean assertions,
			boolean perFile) throws Exception{
		final String label=label("healthy", assertions, perFile);
		final Path input=root.resolve(label+".tsv");
		final Path output=root.resolve(label+".out");
		Files.write(input, HEALTHY.getBytes(StandardCharsets.US_ASCII));

		RunResult result=runAnalyzeAccession(root, label, input, output,
				assertions, perFile);
		check(result.exit==0, label+": healthy child exited "+result.exit+
				" with stderr\n"+result.stderr);
		check(EXPECTED.equals(readAscii(output)), label+
				": unexpected healthy output\n"+readAscii(output));
	}

	/** Runs one child with a non-ASCII poison byte after a valid row.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param perFile Value for AnalyzeAccession's perfile flag
	 * @param baseline True to expect the current lost-worker failure
	 * @throws Exception If the child launch or file IO fails */
	private static void checkPoison(Path root, boolean assertions,
			boolean perFile, boolean baseline) throws Exception{
		final String label=label("poison", assertions, perFile);
		final Path input=root.resolve(label+".tsv");
		final Path output=root.resolve(label+".out");
		Files.write(input, poisonInput());

		RunResult result=runAnalyzeAccession(root, label, input, output,
				assertions, perFile);
		check(result.stderr.contains("ArrayIndexOutOfBoundsException"), label+
				": poison byte did not produce the expected worker diagnostic\n"+
				result.stderr);
		check(result.stderr.contains("tax.AnalyzeAccession$ProcessThread.increment"),
				label+": poison byte did not fail in ProcessThread.increment\n"+
				result.stderr);
		if(baseline){
			check(result.exit==0, label+": baseline should lose the worker failure, exit="+
					result.exit+"\nstderr\n"+result.stderr);
			check(EXPECTED.equals(readAscii(output)), label+
					": baseline should merge the valid partial count\n"+
					readAscii(output));
		}else{
			check(result.exit!=0, label+": fixed worker failure should propagate");
		}
	}

	/** Runs valid and truncated gzip cases through scalar BF1 decoding.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param perFile Value for AnalyzeAccession's perfile flag
	 * @param baseline True to expect the current lost close-error signal
	 * @throws Exception If fixture construction, child launch or file IO fails */
	private static void checkTruncatedGzip(Path root, boolean assertions,
			boolean perFile, boolean baseline) throws Exception{
		final byte[] complete=gzip(HEALTHY.getBytes(StandardCharsets.US_ASCII));
		final byte[] truncated=Arrays.copyOf(complete, complete.length-8);

		final String healthyLabel=label("gzip", assertions, perFile);
		final Path healthyInput=root.resolve(healthyLabel+".tsv.gz");
		final Path healthyOutput=root.resolve(healthyLabel+".out");
		Files.write(healthyInput, complete);

		RunResult healthy=runAnalyzeAccession(root, healthyLabel, healthyInput,
				healthyOutput, assertions, perFile, scalarBf1Flags());
		check(healthy.exit==0, healthyLabel+": complete gzip child exited "+
				healthy.exit+" with stderr\n"+healthy.stderr);
		check(EXPECTED.equals(readAscii(healthyOutput)), healthyLabel+
				": unexpected complete gzip output\n"+readAscii(healthyOutput));

		final String poisonLabel=label("truncated-gzip", assertions, perFile);
		final Path poisonInput=root.resolve(poisonLabel+".tsv.gz");
		final Path poisonOutput=root.resolve(poisonLabel+".out");
		Files.write(poisonInput, truncated);

		RunResult poison=runAnalyzeAccession(root, poisonLabel, poisonInput,
				poisonOutput, assertions, perFile, scalarBf1Flags());
		check(poison.stderr.contains("Error: truncated or corrupt input"),
				poisonLabel+": truncated gzip did not report the BF1 read error\n"+
				poison.stderr);
		if(baseline){
			check(poison.exit==0, poisonLabel+
					": baseline should lose the BF1 close error, exit="+
					poison.exit+"\nstderr\n"+poison.stderr);
			check(EXPECTED.equals(readAscii(poisonOutput)), poisonLabel+
					": baseline should merge valid output before EOF\n"+
					readAscii(poisonOutput));
		}else{
			check(poison.exit!=0, poisonLabel+
					": fixed BF1 close error should propagate");
		}
	}

	/** Runs two complete gzip inputs through the multi-file per-file path.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @throws Exception If fixture construction, child launch or file IO fails */
	private static void checkPerFileTwoGzip(Path root, boolean assertions)
			throws Exception{
		final byte[] complete=gzip(HEALTHY.getBytes(StandardCharsets.US_ASCII));

		final String label=label("perfile-two-gzip", assertions, true);
		final Path input1=root.resolve(label+"-1.tsv.gz");
		final Path input2=root.resolve(label+"-2.tsv.gz");
		final Path output=root.resolve(label+".out");
		Files.write(input1, complete);
		Files.write(input2, complete);

		RunResult result=runAnalyzeAccession(root, label,
				input1.toString()+","+input2.toString(), output,
				assertions, true, scalarBf1Flags());
		check(result.exit==0, label+": two-gzip child exited "+
				result.exit+" with stderr\n"+result.stderr);
		check(EXPECTED_DOUBLE.equals(readAscii(output)), label+
				": unexpected two-gzip output\n"+readAscii(output));
	}

	/** Runs one good per-file input before a truncated per-file gzip input.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param baseline True to expect the current lost close-error signal
	 * @throws Exception If fixture construction, child launch or file IO fails */
	private static void checkPerFileTruncatedGzip(Path root, boolean assertions,
			boolean baseline) throws Exception{
		final byte[] complete=gzip(HEALTHY.getBytes(StandardCharsets.US_ASCII));
		final byte[] truncated=Arrays.copyOf(complete, complete.length-8);

		final String label=label("perfile-truncated-gzip", assertions, true);
		final Path goodInput=root.resolve(label+"-good.tsv.gz");
		final Path poisonInput=root.resolve(label+"-truncated.tsv.gz");
		final Path output=root.resolve(label+".out");
		Files.write(goodInput, complete);
		Files.write(poisonInput, truncated);

		RunResult result=runAnalyzeAccession(root, label,
				goodInput.toString()+","+poisonInput.toString(), output,
				assertions, true, scalarBf1Flags());
		check(result.stderr.contains("Error: truncated or corrupt input"),
				label+": truncated gzip did not report the BF1 read error\n"+
				result.stderr);
		if(baseline){
			check(result.exit==0, label+
					": baseline should lose the second BF1 close error, exit="+
					result.exit+"\nstderr\n"+result.stderr);
			check(EXPECTED_DOUBLE.equals(readAscii(output)), label+
					": baseline should merge complete plus partial output\n"+
					readAscii(output));
		}else{
			check(result.exit!=0, label+
					": fixed second-reader close error should propagate");
			check(!Files.exists(output), label+
					": fixed second-reader close error should block writer construction\n"+
					readAscii(output));
		}
	}

	/** Runs AnalyzeAccession in a child JVM.
	 * @param root Temporary artifact directory
	 * @param label Case label for captured stdout/stderr
	 * @param input Input TSV file
	 * @param output Requested pattern summary path
	 * @param assertions True for -ea, false for -da
	 * @param perFile Value for AnalyzeAccession's perfile flag
	 * @param extraArgs Additional AnalyzeAccession arguments
	 * @return Child exit code and captured stderr
	 * @throws Exception If process launch, wait or file IO fails */
	private static RunResult runAnalyzeAccession(Path root, String label,
			Path input, Path output, boolean assertions,
			boolean perFile, String... extraArgs) throws Exception{
		return runAnalyzeAccession(root, label, input.toString(), output,
				assertions, perFile, extraArgs);
	}

	/** Runs AnalyzeAccession in a child JVM.
	 * @param root Temporary artifact directory
	 * @param label Case label for captured stdout/stderr
	 * @param inputNames Exact AnalyzeAccession input argument value
	 * @param output Requested pattern summary path
	 * @param assertions True for -ea, false for -da
	 * @param perFile Value for AnalyzeAccession's perfile flag
	 * @param extraArgs Additional AnalyzeAccession arguments
	 * @return Child exit code and captured stderr
	 * @throws Exception If process launch, wait or file IO fails */
	private static RunResult runAnalyzeAccession(Path root, String label,
			String inputNames, Path output, boolean assertions,
			boolean perFile, String... extraArgs) throws Exception{
		final ArrayList<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+"/bin/java");
		command.add("-Xmx128m");
		command.add("-XX:ActiveProcessorCount=2");
		command.add(assertions ? "-ea" : "-da");
		if(!"1.8".equals(System.getProperty("java.specification.version"))){
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add("tax.AnalyzeAccession");
		command.addAll(Arrays.asList("in="+inputNames,
				"out="+output.toString(), "overwrite=t", "perfile="+perFile, "t=2"));
		command.addAll(Arrays.asList(extraArgs));

		final Path stdout=root.resolve(label+".stdout");
		final Path stderr=root.resolve(label+".stderr");
		Process child=new ProcessBuilder(command)
				.redirectOutput(stdout.toFile())
				.redirectError(stderr.toFile())
				.start();
		child.getOutputStream().close();
		return new RunResult(child.waitFor(), readAscii(stderr));
	}

	/** Compresses bytes into one complete gzip member.
	 * @param payload Uncompressed payload
	 * @return Complete gzip bytes
	 * @throws Exception If gzip generation fails */
	private static byte[] gzip(byte[] payload) throws Exception{
		ByteArrayOutputStream baos=new ByteArrayOutputStream();
		try(GZIPOutputStream gzip=new GZIPOutputStream(baos)){
			gzip.write(payload);
		}
		return baos.toByteArray();
	}

	/** Builds parser flags forcing Java BF1 gzip decoding.
	 * @return AnalyzeAccession arguments for scalar in-process BF1 decoding */
	private static String[] scalarBf1Flags(){
		return new String[] {"bf1=t", "unpigz=f", "gunzip=f", "unbgzip=f", "simd=f"};
	}

	/** Builds an accession2taxid input with a valid row before a poison byte.
	 * @return Input bytes whose poison row fails AnalyzeAccession workers */
	private static byte[] poisonInput(){
		byte[] prefix=HEALTHY.getBytes(StandardCharsets.US_ASCII);
		byte[] suffix="\tBAD.1\t8\t0\n".getBytes(StandardCharsets.US_ASCII);
		byte[] result=new byte[prefix.length+1+suffix.length];
		System.arraycopy(prefix, 0, result, 0, prefix.length);
		result[prefix.length]=(byte)0x80;
		System.arraycopy(suffix, 0, result, prefix.length+1, suffix.length);
		return result;
	}

	/** Constructs a compact artifact label.
	 * @param prefix Healthy or poison case name
	 * @param assertions True for -ea, false for -da
	 * @param perFile Value for AnalyzeAccession's perfile flag
	 * @return Filename-safe label */
	private static String label(String prefix, boolean assertions, boolean perFile){
		return prefix+(assertions ? "-ea" : "-da")+(perFile ? "-pf" : "-shared");
	}

	/** Reads one US-ASCII file, returning empty text when absent.
	 * @param path File to read
	 * @return Decoded file contents
	 * @throws Exception If the file exists but cannot be read */
	private static String readAscii(Path path) throws Exception{
		return Files.exists(path) ?
				new String(Files.readAllBytes(path), StandardCharsets.US_ASCII) : "";
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTempTree(Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Child result needed for AnalyzeAccession failure oracles. */
	private static final class RunResult{

		/** Captures process status and stderr.
		 * @param exit_ Child process exit code
		 * @param stderr_ Captured standard error */
		RunResult(int exit_, String stderr_){
			exit=exit_;
			stderr=stderr_;
		}

		/** Child process exit code. */
		final int exit;
		/** Captured standard error. */
		final String stderr;
	}
}
