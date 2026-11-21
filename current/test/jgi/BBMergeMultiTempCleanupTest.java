package test.jgi;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;

import jgi.BBMergeMulti;
import parse.Parse;

/** Finite filesystem failures in the real iterative merge wrapper (JGI-085).
 * Reflection separates construction from fault injection; no production hook is added.
 * The allocation failure occurs before child launch. Concatenation cases copy raw bytes,
 * not biological records, and do not exercise a merger process.
 * @author Shinobu
 */
public final class BBMergeMultiTempCleanupTest{

	/** Requires phase=baseline|fixed and a new outdir; only orphan-count expectations differ. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown temporary-cleanup option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(BBMergeMultiTempCleanupTest.class.desiredAssertionStatus(), "Run the dev launcher with assertions enabled.");
		final BBMergeMultiTempCleanupTest test=new BBMergeMultiTempCleanupTest(Paths.get(directory).toAbsolutePath().normalize(), phase.equals("fixed"));
		test.residualAllocationFailure();
		test.aggregateInputFailure();
		test.aggregateControl();
		test.check(test.cases==3, "Expected two failure cases and one concatenation control: "+test.cases);
		test.check(test.constructions==3, "Expected one real constructor per case: "+test.constructions);
		test.check(test.invocations==3, "Expected one real method invocation per case: "+test.invocations);
		final String summary="BBMERGE_MULTI_TEMP_CLEANUP_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" constructors="+test.constructions+
			" methods="+test.invocations+" checks="+test.checks+" failures="+test.failures+"\n";
		writeNew(test.out.resolve("summary.txt"), summary.getBytes(StandardCharsets.US_ASCII));
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+test.out);}
		System.out.print(summary);
	}

	/** Resolves only the actual constructor and methods used by these three cases. */
	private BBMergeMultiTempCleanupTest(final Path out_, final boolean fixed_) throws Exception{
		out=out_;
		fixed=fixed_;
		constructor=BBMergeMulti.class.getDeclaredConstructor(String[].class);
		process=BBMergeMulti.class.getDeclaredMethod("process");
		concatenate=BBMergeMulti.class.getDeclaredMethod("concatenate", ArrayList.class);
		constructor.setAccessible(true); process.setAccessible(true); concatenate.setAccessible(true);
		require(process.getReturnType()==void.class && concatenate.getReturnType()==void.class,
			"These oracles expect the real void process/concatenate contracts.");
		Files.createDirectory(out);
		writeNew(out.resolve("cases.tsv"), ("case\tfailure\tcause\tresidualTemps\taggregateTemps"+
			"\toutputMatches\tinputUnchanged\n").getBytes(StandardCharsets.US_ASCII));
	}

	/** The second allocation fails; the baseline has not yet registered the first temp. */
	private void residualAllocationFailure() throws Exception{
		cases++;
		final Path dir=Files.createDirectory(out.resolve("residual-allocation"));
		final Path input=createInput(dir);
		final Path mergedParent=Files.createDirectory(dir.resolve("merged-parent"));
		final Path residualParent=Files.createDirectory(dir.resolve("unmerged-parent"));
		final Path merged=mergedParent.resolve("merged.fq"), residual=residualParent.resolve("unmerged.fq");
		final BBMergeMulti multi=create(input, merged, residual);
		Files.delete(mergedParent);
		writeNew(mergedParent, BLOCKER);
		require(Files.isRegularFile(mergedParent) && Files.isDirectory(residualParent),
			"Only the owned merged parent must be blocked, before process can launch a child.");
		final Throwable failure=invoke(process, multi);
		final ArrayList<Path> leaked=matching(residualParent, "bbmerge_residual_pass1_*.fq");
		final boolean outputsAbsent=!Files.exists(merged) && !Files.exists(residual);
		final boolean blockerUnchanged=sameBytes(mergedParent, BLOCKER), inputUnchanged=sameBytes(input, INPUT);
		row("residual-allocation", failure, leaked.size(), 0, outputsAbsent && blockerUnchanged, inputUnchanged);
		check(failure instanceof RuntimeException && failure.getCause() instanceof IOException &&
			"Could not launch an iterative BBMerge pass.".equals(failure.getMessage()),
			"Expected the original allocation IOException wrapped by process: "+failure);
		check(leaked.size()==(fixed ? 0 : 1), "Unexpected residual count for this phase: "+leaked);
		check(emptyRegularFiles(leaked), "Leaked residual must be an unwritten ordinary temp file.");
		check(outputsAbsent, "Neither final output may be created before merged-temp allocation succeeds.");
		check(blockerUnchanged, "The user-output parent blocker must not be changed.");
		check(inputUnchanged, "Allocation failure must not change the input.");
	}

	/** Missing first input preserves final output; the baseline leaves an empty aggregate. */
	private void aggregateInputFailure() throws Exception{
		cases++;
		final Path dir=Files.createDirectory(out.resolve("aggregate-input-failure"));
		final Path input=createInput(dir), output=dir.resolve("merged.fq"), missing=dir.resolve("missing.fq");
		final BBMergeMulti multi=create(input, output, dir.resolve("unmerged.fq"));
		writeNew(output, SENTINEL);
		final ArrayList<File> inputs=new ArrayList<File>(); inputs.add(missing.toFile());
		final Throwable failure=invoke(concatenate, multi, inputs);
		final ArrayList<Path> leaked=matching(dir, "bbmerge_aggregate_*.fq");
		final boolean outputUnchanged=sameBytes(output, SENTINEL), inputUnchanged=sameBytes(input, INPUT);
		row("aggregate-input-failure", failure, 0, leaked.size(), outputUnchanged, inputUnchanged);
		check(failure instanceof FileNotFoundException, "Expected the real missing-input failure: "+failure);
		check(leaked.size()==(fixed ? 0 : 1), "Unexpected aggregate count for this phase: "+leaked);
		check(emptyRegularFiles(leaked), "Missing first input must leave an empty ordinary aggregate.");
		check(outputUnchanged, "Failed concatenation must preserve the final-output sentinel.");
		check(!Files.exists(missing), "The missing input must not be created by concatenation.");
		check(inputUnchanged, "Failure must not change the constructor's input.");
	}

	/** Independent literal byte oracle for the ordinary two-input concatenation path. */
	private void aggregateControl() throws Exception{
		cases++;
		final Path dir=Files.createDirectory(out.resolve("aggregate-control"));
		final Path input=createInput(dir), output=dir.resolve("merged.fq");
		final Path first=dir.resolve("first.fq"), second=dir.resolve("second.fq");
		final BBMergeMulti multi=create(input, output, dir.resolve("unmerged.fq"));
		writeNew(output, SENTINEL); writeNew(first, FIRST); writeNew(second, SECOND);
		final ArrayList<File> inputs=new ArrayList<File>(); inputs.add(first.toFile()); inputs.add(second.toFile());
		final Throwable failure=invoke(concatenate, multi, inputs);
		final ArrayList<Path> leaked=matching(dir, "bbmerge_aggregate_*.fq");
		final boolean outputMatches=sameBytes(output, COMBINED), inputUnchanged=sameBytes(input, INPUT);
		row("aggregate-control", failure, 0, leaked.size(), outputMatches, inputUnchanged);
		check(failure==null, "Ordinary concatenation must return normally: "+failure);
		check(leaked.isEmpty(), "Successful move must leave no aggregate temp: "+leaked);
		check(outputMatches, "Concatenation must replace the sentinel with first\\nsecond\\n exactly.");
		check(sameBytes(first, FIRST), "First concatenation input changed.");
		check(sameBytes(second, SECOND), "Second concatenation input changed.");
		check(inputUnchanged, "Ordinary concatenation changed the constructor's input.");
	}

	/** Uses distinct named outputs and no tmpdir override, preserving separate parent paths. */
	private BBMergeMulti create(final Path input, final Path output, final Path residual) throws Exception{
		final String[] args={"k=0,1", "in="+input, "out="+output, "outu="+residual, "ow=t"};
		constructions++;
		try{return constructor.newInstance((Object)args);}
		catch(InvocationTargetException e){throw new AssertionError("Valid owned fixture setup was rejected.", e.getCause());}
	}

	/** Returns the actual target failure; reflection failures remain setup failures. */
	private Throwable invoke(final Method method, final BBMergeMulti target, final Object...args) throws Exception{
		invocations++;
		try{method.invoke(target, args); return null;}
		catch(InvocationTargetException e){return e.getCause();}
	}

	private static Path createInput(final Path dir) throws IOException{
		final Path input=dir.resolve("input.fq"); writeNew(input, INPUT); return input;
	}

	/** Enumerates only the owned case directory and the relevant temp-name family. */
	private static ArrayList<Path> matching(final Path dir, final String glob) throws IOException{
		final ArrayList<Path> paths=new ArrayList<Path>();
		try(DirectoryStream<Path> stream=Files.newDirectoryStream(dir, glob)){
			for(final Path path : stream){paths.add(path);}
		}
		return paths;
	}

	private static boolean emptyRegularFiles(final ArrayList<Path> paths) throws IOException{
		for(final Path path : paths){if(!Files.isRegularFile(path) || Files.size(path)!=0){return false;}}
		return true;
	}

	private static boolean sameBytes(final Path path, final byte[] expected) throws IOException{
		return Files.isRegularFile(path) && Arrays.equals(expected, Files.readAllBytes(path));
	}

	private static void writeNew(final Path path, final byte[] bytes) throws IOException{
		Files.write(path, bytes, StandardOpenOption.CREATE_NEW);
	}

	private void row(final String name, final Throwable failure, final int residuals, final int aggregates,
			final boolean outputMatches, final boolean inputUnchanged) throws IOException{
		final String text=name+"\t"+failureText(failure)+"\t"+failureText(failure==null ? null : failure.getCause())+
			"\t"+residuals+"\t"+aggregates+"\t"+outputMatches+"\t"+inputUnchanged+"\n";
		Files.write(out.resolve("cases.tsv"), text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	private static String failureText(final Throwable failure){
		return failure==null ? "none" : failure.toString().replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
	}

	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}

	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private final Path out;
	private final boolean fixed;
	private final Constructor<BBMergeMulti> constructor;
	private final Method process, concatenate;
	private int cases, constructions, invocations, checks, failures;
	private static final byte[] INPUT=("@pair/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"+
		"@pair/2\nACGTACGTACGT\n+\nIIIIIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] BLOCKER="blocked-parent\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] SENTINEL="prior-output\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FIRST="first\n".getBytes(StandardCharsets.US_ASCII), SECOND="second\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] COMBINED="first\nsecond\n".getBytes(StandardCharsets.US_ASCII);
}
