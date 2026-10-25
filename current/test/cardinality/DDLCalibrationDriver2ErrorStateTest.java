package test.cardinality;

import java.io.File;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import cardinality.DDLCalibrationDriver2;

/** Checks that DDLCalibrationDriver2 reports writer and worker failures.
 * @author Jean
 */
public final class DDLCalibrationDriver2ErrorStateTest{

	/** Runs parent or child fixture mode. */
	public static void main(String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			runChild(args);
			return;
		}
		check(args.length>=1 && args.length<=2,
			"Expected expect=baseline or expect=fixed [output root]");
		final boolean fixed=args[0].equals("expect=fixed");
		check(fixed || args[0].equals("expect=baseline"), "Unknown expectation: "+args[0]);

		final Path root=args.length==2 ? Paths.get(args[1]) :
			Files.createTempDirectory("bbtools-ddlcal2-errors-");
		Files.createDirectories(root);
		final boolean hasFull=Files.exists(Paths.get("/dev/full"));
		if(hasFull){
			checkWriterFailure(root, "out1_full", "-ea", fixed, true);
			checkWriterFailure(root, "out3_full", "-ea", fixed, true);
			checkWriterFailure(root, "out4_full", "-ea", fixed, true);
			checkWriterFailure(root, "out1_full", "-da", fixed, false);
			checkWriterFailure(root, "out3_full", "-da", fixed, false);
			checkWriterFailure(root, "out4_full", "-da", fixed, false);
		}else{
			System.err.println("Skipping /dev/full writer cases; /dev/full is absent");
		}
		checkWriterFailure(root, "out3_bad_parent", "-ea", fixed, true);
		checkWriterFailure(root, "out4_bad_parent", "-ea", fixed, true);
		checkWriterFailure(root, "out3_bad_parent", "-da", fixed, false);
		checkWriterFailure(root, "out4_bad_parent", "-da", fixed, false);
		checkPartialAggregate(root, "-ea", fixed);
		checkPartialAggregate(root, "-da", fixed);
		checkHealthy(root, "-ea");
		checkHealthy(root, "-da");
		if(fixed){
			checkInterruption(root, "-ea");
			checkInterruption(root, "-da");
		}
		System.out.println("PASS\tfixed="+fixed);
	}

	/** Runs one writer-failure subprocess and checks exit, marker and errorState. */
	private static void checkWriterFailure(final Path root, final String name,
			final String assertionFlag, final boolean fixed, final boolean gated) throws Exception{
		final ProbeResult result=runProbe(root, name+"_"+assertionFlag, assertionFlag, name);
		if(fixed && gated){
			check(result.exit!=0, name+" "+assertionFlag+" should halt loudly");
			check(!result.marker, name+" "+assertionFlag+" should halt before the after marker");
		}else{
			check(result.exit==0, name+" "+assertionFlag+" should survive; exit="+result.exit);
			check(result.marker, name+" "+assertionFlag+" should write its after marker");
			check(result.errorState==fixed,
				name+" "+assertionFlag+" errorState expected "+fixed+
				" observed "+result.errorState+" log="+result.log);
		}
	}

	/** Runs the partial-worker aggregate probe; fixed source must reject it normally. */
	private static void checkPartialAggregate(final Path root, final String assertionFlag,
			final boolean fixed) throws Exception{
		final ProbeResult result=runProbe(root, "partial_"+assertionFlag,
			assertionFlag, "partial_worker");
		if(fixed){
			check(result.exit!=0, "Fixed partial-worker aggregate should fail");
			check(!result.marker, "Fixed partial-worker aggregate should not mark success");
			check(!Files.exists(result.root.resolve("out1.tsv")),
				"Fixed partial-worker aggregate should not publish out1");
			check(!Files.exists(result.root.resolve("out3.cf")),
				"Fixed partial-worker aggregate should not publish out3");
			check(logContains(result.log, "ArrayIndexOutOfBoundsException"),
				"Fixed partial-worker aggregate should preserve its worker cause");
		}else{
			check(result.exit==0, "Baseline partial-worker aggregate should be swallowed");
			check(result.marker, "Baseline partial-worker aggregate should write its marker");
		}
	}

	/** Runs a healthy subprocess and checks it writes all three outputs. */
	private static void checkHealthy(final Path root, final String assertionFlag) throws Exception{
		final ProbeResult result=runProbe(root, "healthy_"+assertionFlag,
			assertionFlag, "healthy");
		check(result.exit==0, "Healthy "+assertionFlag+" failed; exit="+result.exit+
			" log="+result.log);
		check(result.marker, "Healthy "+assertionFlag+" did not write its marker");
		check(!result.errorState, "Healthy "+assertionFlag+" should not set errorState");
		check(Files.size(result.root.resolve("out1.tsv"))>0,
			"Healthy "+assertionFlag+" should write out1");
		check(Files.size(result.root.resolve("out3.cf"))>0,
			"Healthy "+assertionFlag+" should write out3");
		check(Files.size(result.root.resolve("out4.tsv"))>0,
			"Healthy "+assertionFlag+" should write out4");
	}

	/** Runs a caller-interruption probe; fixed source must fail without stranding a worker. */
	private static void checkInterruption(final Path root, final String assertionFlag)
			throws Exception{
		final ProbeResult result=runProbe(root, "interrupted_"+assertionFlag,
			assertionFlag, "interrupted_caller");
		check(result.exit==0, "Interrupted caller "+assertionFlag+
			" failed; exit="+result.exit+" log="+result.log);
		check(!result.marker, "Interrupted caller "+assertionFlag+
			" should not mark a normal process return");
		check(logContains(result.log,
			"interruption\tfailure=true\tinterruptionCause=true"+
			"\tinterrupted=true\tliveCalThread=false\tpublished=false"),
			"Interrupted caller "+assertionFlag+
			" did not report the fixed failure contract; log="+result.log);
	}

	/** Launches a child JVM with a case-specific watchdog and captured output. */
	private static ProbeResult runProbe(final Path root, final String dirName,
			final String assertionFlag, final String childCase) throws Exception{
		final Path caseRoot=root.resolve(dirName);
		Files.createDirectories(caseRoot);
		final String java=System.getProperty("java.home")+File.separator+
			"bin"+File.separator+"java";
		final ArrayList<String> command=new ArrayList<String>();
		command.add(java);
		command.add("-Xmx256m");
		command.add("-Xms64m");
		command.add("-XX:+ExitOnOutOfMemoryError");
		if(!System.getProperty("java.specification.version", "8").equals("1.8")){
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}
		command.add(assertionFlag);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(DDLCalibrationDriver2ErrorStateTest.class.getName());
		command.add("child");
		command.add(childCase);
		command.add(caseRoot.toString());

		final Path log=caseRoot.resolve("child.log");
		Files.deleteIfExists(log);
		Files.deleteIfExists(caseRoot.resolve("errorState.txt"));
		Files.deleteIfExists(caseRoot.resolve("after.marker"));
		Files.deleteIfExists(caseRoot.resolve("out1.tsv"));
		Files.deleteIfExists(caseRoot.resolve("out3.cf"));
		Files.deleteIfExists(caseRoot.resolve("out4.tsv"));
		final ProcessBuilder pb=new ProcessBuilder(command);
		pb.redirectErrorStream(true);
		pb.redirectOutput(log.toFile());
		final Process child=pb.start();
		if(!child.waitFor(30, TimeUnit.SECONDS)){
			child.destroyForcibly();
			child.waitFor();
			throw new AssertionError(childCase+" "+assertionFlag+
				" timed out; log="+log);
		}
		final Path errorFile=caseRoot.resolve("errorState.txt");
		final boolean marker=Files.exists(caseRoot.resolve("after.marker"));
		final boolean errorState=Files.exists(errorFile) &&
			new String(Files.readAllBytes(errorFile), StandardCharsets.UTF_8).contains("true");
		final ProbeResult result=new ProbeResult(caseRoot, log, child.exitValue(),
			marker, errorState);
		System.out.println(dirName+"\texit="+result.exit+"\tmarker="+result.marker+
			"\terrorState="+result.errorState);
		return result;
	}

	/** Runs a child case; the after marker is only written after the tested path returns. */
	private static void runChild(final String[] args) throws Exception{
		check(args.length==3, "Expected child <case> <root>, got "+args.length);
		final String childCase=args[1];
		final Path root=Paths.get(args[2]);
		Files.createDirectories(root);
		if("interrupted_caller".equals(childCase)){
			runInterruptedCaller(root);
			return;
		}
		final DDLCalibrationDriver2 driver=newDriver(root,
			"partial_worker".equals(childCase) ? "healthy" : childCase);
		if("partial_worker".equals(childCase)){
			setNumExtendedTypes(0);
		}
		invokeProcess(driver);
		writeErrorState(root.resolve("errorState.txt"), driver);
		Files.write(root.resolve("after.marker"), "after\n".getBytes(StandardCharsets.UTF_8));
	}

	/** Interrupts the caller while its CalThread is definitely blocked before work starts. */
	private static void runInterruptedCaller(final Path root) throws Exception{
		final DDLCalibrationDriver2 driver=newDriver(root, "healthy");
		final Throwable[] failure=new Throwable[1];
		final boolean[] interrupted=new boolean[1];
		final Thread caller=new Thread(new Runnable(){
			@Override
			public void run(){
				try{
					invokeProcess(driver);
					Files.write(root.resolve("after.marker"),
						"after\n".getBytes(StandardCharsets.UTF_8));
				}catch(Throwable t){
					failure[0]=t;
				}finally{
					interrupted[0]=Thread.currentThread().isInterrupted();
				}
			}
		}, "DDLCalibrationDriver2-interrupt-caller");

		synchronized(DDLCalibrationDriver2.class){
			caller.start();
			waitForState(caller, Thread.State.WAITING, 10000);
			caller.interrupt();
			Thread.sleep(100);
			check(caller.isAlive(),
				"Caller returned before the blocked CalThread could start");
		}
		joinOrFail(caller, 10000, "interrupted caller");
		if(failure[0]!=null){failure[0].printStackTrace();}
		final boolean interruptionCause=causedByInterruption(failure[0]);
		final boolean live=hasLiveCalThread();
		final boolean published=Files.exists(root.resolve("out1.tsv")) ||
			Files.exists(root.resolve("out3.cf"));
		System.out.println("interruption\tfailure="+(failure[0]!=null)+
			"\tinterruptionCause="+interruptionCause+
			"\tinterrupted="+interrupted[0]+"\tliveCalThread="+live+
			"\tpublished="+published);
		check(failure[0]!=null, "Interrupted caller should fail");
		check(interruptionCause,
			"Interrupted caller should report interruption as its failure cause");
		check(interrupted[0], "Interrupted caller should restore the interrupt flag");
		check(!live, "Interrupted caller stranded a live CalThread");
		check(!published, "Interrupted caller published a parent aggregate table");
	}

	/** Builds a bounded one-DDL command line for each subprocess probe. */
	private static String[] driverArgs(final Path root, final String childCase){
		String out1=root.resolve("out1.tsv").toString();
		String out3=root.resolve("out3.cf").toString();
		String out4=root.resolve("out4.tsv").toString();
		if("out1_full".equals(childCase)){out1="/dev/full"; out3=null; out4=null;}
		else if("out3_full".equals(childCase)){out3="/dev/full"; out4=null;}
		else if("out4_full".equals(childCase)){out3=null; out4="/dev/full";}
		else if("out3_bad_parent".equals(childCase)){
			out3=root.resolve("missing").resolve("out3.cf").toString();
			out4=null;
		}else if("out4_bad_parent".equals(childCase)){
			out3=null;
			out4=root.resolve("missing").resolve("out4.tsv").toString();
		}else{
			check("healthy".equals(childCase), "Unknown child case: "+childCase);
		}

		final ArrayList<String> args=new ArrayList<String>();
		args.add("card=1");
		args.add("ddls=1");
		args.add("buckets=64");
		args.add("k=1");
		args.add("t=1");
		args.add("reportfrac=1");
		args.add("type=dll4");
		args.add("out="+out1);
		if(out3!=null){args.add("out3="+out3);}
		if(out4!=null){args.add("out4="+out4);}
		return args.toArray(new String[args.size()]);
	}

	/** Forces a valid worker to fail after partial metric accumulation. */
	private static void setNumExtendedTypes(final int value) throws Exception{
		final Field f=DDLCalibrationDriver2.class.getDeclaredField("NUM_EXTENDED_TYPES");
		f.setAccessible(true);
		f.setInt(null, value);
	}

	/** Constructs one child driver without running its process. */
	private static DDLCalibrationDriver2 newDriver(final Path root, final String childCase)
			throws Exception{
		final Constructor<DDLCalibrationDriver2> cons=
			DDLCalibrationDriver2.class.getDeclaredConstructor(String[].class);
		cons.setAccessible(true);
		return cons.newInstance((Object)driverArgs(root, childCase));
	}

	/** Invokes process() and preserves the production exception cause. */
	private static void invokeProcess(final DDLCalibrationDriver2 driver) throws Exception{
		final Method m=DDLCalibrationDriver2.class.getDeclaredMethod("process", long.class);
		m.setAccessible(true);
		try{m.invoke(driver, System.nanoTime());}
		catch(InvocationTargetException ite){throwCause(ite.getCause());}
	}

	/** Records DDLCalibrationDriver2.errorState(), when present in fixed source. */
	private static void writeErrorState(final Path path, final DDLCalibrationDriver2 driver)
			throws Exception{
		boolean errorState=false;
		try{
			final Method m=DDLCalibrationDriver2.class.getDeclaredMethod("errorState");
			m.setAccessible(true);
			errorState=Boolean.TRUE.equals(m.invoke(driver));
		}catch(NoSuchMethodException nsme){
			errorState=false;
		}
		Files.write(path, ("errorState="+errorState+"\n").getBytes(StandardCharsets.UTF_8));
	}

	/** Waits until a thread reaches a deterministic rendezvous state. */
	private static void waitForState(final Thread thread, final Thread.State state,
			final long timeoutMillis) throws Exception{
		final long deadline=System.nanoTime()+TimeUnit.MILLISECONDS.toNanos(timeoutMillis);
		while(System.nanoTime()<deadline){
			if(thread.getState()==state){return;}
			check(thread.isAlive(), thread.getName()+" exited before reaching "+state);
			Thread.sleep(10);
		}
		throw new AssertionError(thread.getName()+" timed out waiting for "+state+
			"; observed "+thread.getState());
	}

	/** Joins a bounded helper thread and fails loudly on timeout. */
	private static void joinOrFail(final Thread thread, final long timeoutMillis,
			final String description) throws Exception{
		thread.join(timeoutMillis);
		check(!thread.isAlive(), description+" did not finish in "+timeoutMillis+" ms");
	}

	/** Finds any live calibration worker left after process() should have joined them. */
	private static boolean hasLiveCalThread(){
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread.isAlive() && "cardinality.DDLCalibrationDriver2$CalThread".equals(
					thread.getClass().getName())){
				return true;
			}
		}
		return false;
	}

	/** Checks whether a failure chain identifies the caller interruption. */
	private static boolean causedByInterruption(Throwable failure){
		for(Throwable t=failure; t!=null; t=t.getCause()){
			if(t instanceof InterruptedException){return true;}
			if(String.valueOf(t.getMessage()).toLowerCase().contains("interrupt")){return true;}
		}
		return false;
	}

	/** Rethrows a reflected production failure without wrapping its diagnostic type. */
	private static void throwCause(final Throwable cause) throws Exception{
		if(cause instanceof Exception){throw (Exception)cause;}
		if(cause instanceof Error){throw (Error)cause;}
		throw new RuntimeException(cause);
	}

	/** Tests captured child output for an expected diagnostic. */
	private static boolean logContains(final Path log, final String pattern) throws Exception{
		return new String(Files.readAllBytes(log), StandardCharsets.UTF_8).contains(pattern);
	}

	/** Keeps oracle checks active under both assertion settings. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** A subprocess result with explicit exit, completion marker and error-state flag. */
	private static final class ProbeResult{
		ProbeResult(Path root_, Path log_, int exit_, boolean marker_, boolean errorState_){
			root=root_;
			log=log_;
			exit=exit_;
			marker=marker_;
			errorState=errorState_;
		}

		final Path root, log;
		final int exit;
		final boolean marker, errorState;
	}
}
