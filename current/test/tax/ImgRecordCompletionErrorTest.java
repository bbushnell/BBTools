package test.tax;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import tax.ImgRecord;

/** Finite /dev/full baseline for ImgRecord text-output completion errors.
 * @author Sucrose */
public final class ImgRecordCompletionErrorTest{

	/** Runs a valid-output control and a /dev/full text-output failure check.
	 * @param args Optional child-mode or baseline/fixed expectation arguments
	 * @throws Exception if setup, child execution, output reading, or cleanup fails */
	public static void main(String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			try{
				ImgRecord.main(new String[] {args[1], args[2]});
			}catch(RuntimeException e){
				System.out.println(CHILD_CAUGHT_MARKER+e.getMessage());
			}
			return;
		}

		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("imgrecord-completion-");
		try{
			final Path in=root.resolve("img.tsv");
			Files.write(in, IMG_LINES, StandardCharsets.UTF_8);

			final Path good=root.resolve("img.out.tsv");
			ImgRecord.main(new String[] {in.toString(), good.toString()});
			check(Files.readAllLines(good, StandardCharsets.UTF_8).equals(EXPECTED_LINES),
					"Valid output should contain exactly "+EXPECTED_LINES+".");

			final Path full=Paths.get("/dev/full");
			if(!Files.exists(full)){
				System.err.println("Skipping /dev/full writer case; /dev/full is absent.");
			}else{
				final ChildResult ea=runChild(root.resolve("full-ea.log"), true, in, full);
				final ChildResult da=runChild(root.resolve("full-da.log"), false, in, full);
				if(baseline){
					checkBaseline(ea, "-ea");
					checkBaseline(da, "-da");
				}else{
					check(da.code==0,
							"Fixed -da child should catch and report the failed completion, "
							+"but exited "+da.code+" with: "+da.out);
					check(da.out.contains(CHILD_CAUGHT_MARKER),
							"Fixed -da child did not report the propagated RuntimeException: "+da.out);
					check(da.out.contains("Error writing IMG records") && da.out.contains("/dev/full"),
							"Fixed -da child lacked the expected write-error diagnostic: "+da.out);
					check(ea.code!=0,
							"Fixed -ea child should halt on failed text-output completion.");
					check(!ea.out.contains(CHILD_CAUGHT_MARKER),
							"Fixed -ea child reached the catch marker after assertDie: "+ea.out);
					check(ea.out.contains("Error writing IMG records") && ea.out.contains("/dev/full"),
							"Fixed -ea child lacked the expected write-error diagnostic: "+ea.out);
				}
			}
		}finally{
			deleteTree(root);
		}
		System.out.println("IMGRECORD_COMPLETION_ERROR_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private ImgRecordCompletionErrorTest(){}

	/** Runs one bounded child JVM and returns its merged stdout/stderr.
	 * @param processLog merged stdout/stderr log for the child
	 * @param assertions true to enable child assertions, false to disable them
	 * @param in valid IMG input path
	 * @param out child output target path
	 * @return child exit code and merged stdout/stderr
	 * @throws Exception if the child cannot be run or its output cannot be read */
	private static ChildResult runChild(final Path processLog, final boolean assertions,
			final Path in, final Path out) throws Exception{
		final ProcessBuilder pb=javaProcess(assertions,
				ImgRecordCompletionErrorTest.class.getName(),
				"child", in.toString(), out.toString());
		pb.redirectErrorStream(true);
		pb.redirectOutput(processLog.toFile());
		final Process proc=pb.start();
		final int code=proc.waitFor();
		final String output=new String(Files.readAllBytes(processLog), StandardCharsets.UTF_8);
		return new ChildResult(code, output);
	}

	/** Creates a bounded child JVM with the requested assertion state.
	 * @param assertions true to enable child assertions, false to disable them
	 * @param className fully qualified class name to run
	 * @param args arguments passed after the class name
	 * @return configured child process builder */
	private static ProcessBuilder javaProcess(final boolean assertions, final String className,
			final String... args){
		final ArrayList<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+"/bin/java");
		command.add("-Xmx128m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add(assertions ? "-ea" : "-da");
		if(!"1.8".equals(System.getProperty("java.specification.version"))){
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(className);
		command.addAll(Arrays.asList(args));
		return new ProcessBuilder(command);
	}

	/** Requires the old implementation to ignore the failed child completion.
	 * @param result child result to inspect
	 * @param mode child assertion mode for diagnostics */
	private static void checkBaseline(final ChildResult result, final String mode){
		check(result.code==0,
				"Baseline "+mode+" should expose the ignored TextStreamWriter completion error, "
				+"but the child exited "+result.code+" with: "+result.out);
		check(!result.out.contains(CHILD_CAUGHT_MARKER),
				"Baseline "+mode+" should not throw through the child wrapper: "+result.out);
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root owned temporary directory root
	 * @throws Exception if the tree cannot be walked or deleted */
	private static void deleteTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Test checks stay active under -da.
	 * @param condition condition required for test success
	 * @param message diagnostic message for failure */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Exit code and merged stdout/stderr from one child JVM. */
	private static final class ChildResult{
		/** Records one completed child JVM result.
		 * @param code_ child exit code
		 * @param out_ merged child stdout/stderr */
		private ChildResult(final int code_, final String out_){
			code=code_;
			out=out_;
		}

		/** Child process exit code. */
		private final int code;

		/** Merged child stdout/stderr. */
		private final String out;
	}

	/** Marker proving the child JVM reached the RuntimeException recovery path. */
	private static final String CHILD_CAUGHT_MARKER="IMGRECORD_CHILD_CAUGHT_RUNTIME_EXCEPTION ";

	/** One exact IMG map entry: HashMap key plus ImgRecord text body. */
	private static final java.util.List<String> EXPECTED_LINES=Arrays.asList(
			"123=123\t456\tGenome name");

	/** One valid IMG dump row with the high-quality column populated. */
	private static final java.util.List<String> IMG_LINES=Arrays.asList(
			"123\tGenome name\t456\tyes\tno\tisolate\tunused\tyes");
}
