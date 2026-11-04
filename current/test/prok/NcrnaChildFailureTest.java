package test.prok;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.PosixFilePermission;
import java.util.HashSet;
import java.util.Set;
import java.util.UUID;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.PreParser;
import prok.NcrnaChildProcess;
import prok.NcrnaTextIO.Output;
import fileIO.FileFormat;
import prok.RrnaResourceIO;
import shared.Shared;
import structures.ByteBuilder;

/** Independent failure-chain oracles and a real full-consumption success control.
 * Completion-log-only controls retired with their writer; real failure identities,
 * source-close errors, interruption and owned-group cleanup remain the oracles.
 * Private reflection covers no-process, throwing-source-stream and recorded-error controls;
 * ordinary exit, cancellation and consumption cases use real owned child groups.
 * Preserved controls for NcrnaChildProcess, which remains necessary for the
 * standalone fixed-library regeneration mode. New public tools use NcrnaProcess;
 * its tests do not cover this older owner. Keep these controls while the
 * standalone mode retains that owner, including any future direct-library use.
 * @author Raiden */
public final class NcrnaChildFailureTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaChildFailureTest.class, false);
		String root=null, installation=null;
		for(String arg:pp.args){
			if(arg.startsWith("outdir=")){root=arg.substring(7);}
			else if(arg.startsWith("installation=")){installation=arg.substring(13);}
			else{throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		check(root!=null && installation!=null,
			"Declare fresh outdir and the sole installation");
		check(new File(root).mkdir(), "Fixture output must be fresh: "+root);
		Shared.setThreads(1);
		ReadWrite.setZipThreads(1);
		new NcrnaChildFailureTest(root, installation).run();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaChildFailureTest(String root_, String installation_){
		root=root_;
		runner=new NcrnaChildProcess(new File(installation_));
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void run() throws Exception{
		try(Output output=new Output(FileFormat.testOutput(root+"/observations.tsv", FileFormat.TEXT, null, true, false, false, false))){
			report=output.writer;
			report.println("case\tprimary_type\tprimary_retained\tinterrupt_flag\tsource_close_failure_visible");
			ordinaryFailure();
			ordinarySuccess();
			cancelBeforeRelease();
			callbackFailure("callback_io", new IOException("original callback IO"));
			callbackFailure("callback_runtime", new IllegalStateException("original callback runtime"));
			callbackFailure("callback_error", new AssertionError("original callback error"));
			noProcessMarker();
			closeFailure(false);
			closeFailure(true);
			bodyCleanup(false);bodyCleanup(true);
		}
		check(rows==11, "All independently defined cases must finish; rows="+rows);
		System.out.println("NCRNA_CHILD_FAILURE_CONTROLS_PASS cases="+rows);
	}

	/** Exercises run itself: a known stage error must precede its throwing body cleanup.
	 * The streaming case injects a recorded failure into a finished real pump; it
	 * does not claim a physical input read failed. Only this fresh case directory
	 * loses write permission, which is restored before any inspection or return. */
	private void bodyCleanup(boolean streaming) throws Exception{
		final String label=streaming ? "stream_and_owned_cleanup" : "exit_and_owned_cleanup";
		final File directory=new File(root, label);
		check(directory.mkdir(), "Owned cleanup fixture needs a fresh directory");
		final Set<PosixFilePermission> permissions=Files.getPosixFilePermissions(directory.toPath());
		final String prefix=new File(directory, "child").getPath(), control=root+"/stream_gate";
		final IOException original=new IOException("recorded streaming failure before child exit");
		final Observation coordinatorResult=new Observation();
		final Thread parent=Thread.currentThread();
		final Thread coordinator=new Thread(()->{
			try{
				final Thread pump=awaitPump(prefix);
				write(control+".read", "read\n");
				pump.join(10000);
				check(!pump.isAlive(), "Real input producer must finish before recorded-error injection");
				awaitFile(control+".consumed");
				final Field failure=field(pump.getClass(), "failure");
				check(failure.get(pump)==null, "Injected recorded error must not replace a real streaming failure");
				failure.set(pump, original);
				final Set<PosixFilePermission> denied=new HashSet<PosixFilePermission>(permissions);
				denied.remove(PosixFilePermission.OWNER_WRITE);
				denied.remove(PosixFilePermission.GROUP_WRITE);
				denied.remove(PosixFilePermission.OTHERS_WRITE);
				Files.setPosixFilePermissions(directory.toPath(), denied);
				write(control+".exit", "exit\n");
			}catch(Throwable problem){coordinatorResult.problem=problem;parent.interrupt();}
		}, "recorded-stream-failure-control");
		coordinator.setDaemon(true);
		Throwable problem=null;
		try{
			String[] shards=null;
			if(streaming){
				try(Output input=new Output(FileFormat.testOutput(root+"/stream_input.txt", FileFormat.TEXT, null, true, false, false, false))){
					final byte[] block=new byte[1<<16];java.util.Arrays.fill(block, (byte)'A');
					for(int i=0; i<16; i++){input.writer.print(block);}
				}
				shards=new String[]{root+"/stream_input.txt"};
			}
			final String[] args=streaming ? new String[]{"mode=readwait", "control="+control}
				: new String[]{"mode=failcleanup", "denydir="+directory.getPath()};
			try{
				runner.run(FIXTURE, args, shards, new File(prefix), 1, 64, ()->{
					// time may open after chmod; an existing writable file avoids changing exit7.
					write(prefix+".time.tsv", "");
					if(streaming){coordinator.start();}
				});
			}catch(Throwable caught){problem=caught;}
		}finally{
			try{
				if(streaming){
					coordinator.join(1000);
					if(coordinator.isAlive()){coordinator.interrupt();coordinator.join(15000);}
				}
			}
			finally{Files.setPosixFilePermissions(directory.toPath(), permissions);}
		}
		check(!coordinator.isAlive() && coordinatorResult.problem==null, "Recorded-error coordinator failed: "+coordinatorResult.problem);
		final boolean retained=streaming ? problem==original
			: problem instanceof IOException && problem.getMessage().contains("exited 7");
		check(retained, "Stage failure precedence disagrees with oracle: "+label+" "+problem);
		check(ownedCleanupFailure(problem), "Permission denial must exercise owned process cleanup after the stage failure");
		boolean suppressed=false;
		for(Throwable later:problem.getSuppressed()){suppressed|=ownedCleanupFailure(later);}
		check(suppressed, "Owned cleanup failure must be suppressed beneath the stage error");
		noLiveGroup(prefix);
		observe(label, problem, retained, false, false);
	}

	private void ordinaryFailure() throws Exception{
		final String prefix=root+"/ordinary_failure";
		Throwable problem=null;
		try{runner.run(FIXTURE, new String[]{"mode=fail"}, null, new File(prefix), 1, 64);}
		catch(Throwable caught){problem=caught;}
		check(problem instanceof IOException && problem.getMessage().contains("exited 7"), "Unobstructed child failure retains exit7");
		check(!new File(prefix+".exit.tsv").exists(), "A failed child must not write a removed completion log");
		noLiveGroup(prefix);
		observe("ordinary_failure", problem, true, false, false);
	}

	/** Two encoded members plus two compressed shards and one raw shard must all arrive. */
	private void ordinarySuccess() throws Exception{
		write(root+"/a.txt.gz", "one\n");
		write(root+"/b.txt.gz", "two\n");
		write(root+"/raw.txt", "end\n");
		// Copy opaque encoded bytes: decoding here would invalidate the multi-member control.
		try(OutputStream out=Files.newOutputStream(new File(root+"/members.gz").toPath(), StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE)){
			Files.copy(new File(root+"/a.txt.gz").toPath(), out);
			Files.copy(new File(root+"/b.txt.gz").toPath(), out);
		}
		final String prefix=root+"/ordinary_success", body="one\ntwo\none\ntwo\nend\n";
		final NcrnaChildProcess.Result result=runner.run(FIXTURE, new String[]{"mode=echo", "payload=all input"},
			new String[]{root+"/members.gz", root+"/a.txt.gz", root+"/b.txt.gz", root+"/raw.txt"}, new File(prefix), 1, 64);
		check(result.exitCode==0 && result.stdinBytes==20 && read(prefix+".stdout").equals("all input\n"+body),
			"Successful input consumption must retain all20ASCIIbytes, member/shard order and exit0");
		check(!new File(prefix+".exit.tsv").exists(), "A successful child must not write a removed completion log");
		noLiveGroup(prefix);
		observe("ordinary_success", null, true, false, false);
	}

	/** The callback deliberately restores the incoming flag before rethrowing its known exception. */
	private void cancelBeforeRelease() throws Exception{
		final String prefix=root+"/cancel_before_release";
		final CountDownLatch entered=new CountDownLatch(1), hold=new CountDownLatch(1);
		final InterruptedException original=new InterruptedException("original interrupted callback");
		final Observation result=new Observation();
		final Thread task=new Thread(()->{
			try{
				runner.run(FIXTURE, new String[]{"mode=echo", "payload=must not run"}, null, new File(prefix), 1, 64, ()->{
					awaitIdentity(prefix);
					entered.countDown();
					try{hold.await();}
					catch(InterruptedException interruption){Thread.currentThread().interrupt();throw original;}
				});
			}catch(Throwable problem){result.problem=problem;}
			finally{result.interrupted=Thread.currentThread().isInterrupted();}
		}, "failure-chain-cancel");
		task.setDaemon(true);
		task.start();
		try{
			check(entered.await(30, TimeUnit.SECONDS), "Cancellation callback did not become ready");
			task.interrupt();
			task.join(15000);
			check(!task.isAlive(), "Cancellation must finish owned cleanup");
			check(result.problem==original && result.interrupted, "Retain the original interruption and captured flag independently");
			check(read(prefix+".stdout").isEmpty(), "Cancellation before release must not execute the workload");
			noLiveGroup(prefix);
			observe("cancel_before_release", result.problem, result.problem==original, result.interrupted, false);
		}finally{
			task.interrupt();
			hold.countDown();
			task.join(15000);
		}
	}

	private void callbackFailure(String label, Throwable original) throws Exception{
		final String prefix=root+"/"+label;
		Throwable problem=null;
		try{
			runner.run(FIXTURE, new String[]{"mode=echo", "payload=must not run"}, null, new File(prefix), 1, 64, ()->{
				awaitIdentity(prefix);
				if(original instanceof IOException){throw (IOException)original;}
				if(original instanceof RuntimeException){throw (RuntimeException)original;}
				throw (Error)original;
			});
		}catch(Throwable caught){problem=caught;}
		check(problem==original, "Preserve the exact stage throwable through owned cleanup: "+label);
		noLiveGroup(prefix);
		observe(label, problem, problem==original, false, false);
	}

	/** Missing parent gives Files.write a deterministic error without permission changes. */
	private void noProcessMarker() throws Exception{
		final Object owned=owned(root+"/absent/child.pid");
		final Method stop=owned.getClass().getDeclaredMethod("stop");
		stop.setAccessible(true);
		final Throwable problem=invokeFailure(stop, owned);
		check(problem instanceof IllegalStateException && problem.getCause() instanceof IOException,
			"No-process cancellation-marker failure must remain visible: "+problem);
		check((Boolean)field(owned.getClass(), "stopped").get(owned), "A failed marker still leaves startup cancelled");
		observe("no_process_marker", problem, problem!=null, false, false);
	}

	/** Tests the real finalizer without a public injection API or real child. */
	private void closeFailure(boolean withPrimary) throws Exception{
		final String label=withPrimary ? "primary_and_source_close" : "source_close_only", prefix=root+"/"+label;
		final Object owned=owned(prefix+".pid");
		final Class<?> pumpClass=Class.forName("prok.NcrnaChildProcess$InputPump");
		final Constructor<?> constructor=pumpClass.getDeclaredConstructor(owned.getClass(), String[].class);
		constructor.setAccessible(true);
		final Object pump=constructor.newInstance(owned, new String[0]);
		final IOException original=withPrimary ? new IOException("original streaming failure") : null;
		final IOException close=new IOException("source close failure");
		final ThrowingInput input=new ThrowingInput(close);
		field(pumpClass, "current").set(pump, input);
		field(pumpClass, "failure").set(pump, original);
		final Thread hook=new Thread(()->{}, "throwing-source-control");
		Runtime.getRuntime().addShutdownHook(hook);
		final Throwable problem;
		try{
			final Method finish=NcrnaChildProcess.class.getDeclaredMethod("finish", owned.getClass(), pumpClass, Thread.class, Throwable.class);
			finish.setAccessible(true);
			problem=(Throwable)finish.invoke(null, owned, pump, hook, original);
			check((Boolean)field(owned.getClass(), "stopped").get(owned), "Finalization must still cancel owned startup");
			check(!Runtime.getRuntime().removeShutdownHook(hook), "Source close failure must not skip hook removal");
		}finally{Runtime.getRuntime().removeShutdownHook(hook);}
		check(problem==(withPrimary ? original : close) && contains(problem, close, 0), "Preserve an earlier primary or promote the sole close failure");
		check(input.closes>0 && field(pumpClass, "failure").get(pump)==original, "Source close must not overwrite the pump's streaming failure");
		observe(label, problem, problem==(withPrimary ? original : close), false, contains(problem, close, 0));
	}

	private void awaitIdentity(String prefix) throws InterruptedException{
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(30);
		while(!new File(prefix+".pid").isFile() && System.nanoTime()<deadline){Thread.sleep(10);}
		check(new File(prefix+".pid").isFile(), "Owned identity must exist before the controlled callback failure");
	}

	private void noLiveGroup(String prefix) throws Exception{
		final LineParser1 parser=new LineParser1('\t');
		parser.set(read(prefix+".pid").trim().getBytes(StandardCharsets.US_ASCII));
		check(parser.terms()==2, "Owned identity must retain nonce and group");
		final int group=parser.parseInt(1);
		check(group>1, "Only inspect the recorded positive owned group");
		final Process ps=new ProcessBuilder("ps", "-g", Integer.toString(group), "-o", "stat=")
			.redirectOutput(new File(prefix+".remaining")).redirectError(new File(prefix+".ps.stderr")).start();
		try{check(ps.waitFor(5, TimeUnit.SECONDS), "Owned-group inspection did not finish");}
		finally{if(ps.isAlive()){ps.destroyForcibly();ps.waitFor(5, TimeUnit.SECONDS);}}
		final ByteFile input=RrnaResourceIO.open(prefix+".remaining");
		try{
			for(byte[] line; (line=input.nextLine())!=null;){
				final String state=new String(line, StandardCharsets.US_ASCII).trim();
				check(state.isEmpty() || state.startsWith("Z"), "A failed invocation left a live owned process: "+state);
			}
		}finally{check(!input.close(), "Owned-group status read failed");}
	}

	private void observe(String label, Throwable problem, boolean retained, boolean interrupted, boolean sourceClose){
		report.print(new ByteBuilder().append(label).tab().append(problem==null ? "none" : problem.getClass().getName())
			.tab().append(retained).tab().append(interrupted).tab().append(sourceClose).nl());
		rows++;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Finds only this invocation's live pump while its pipe is blocked by the read gate. */
	private static Thread awaitPump(String prefix) throws Exception{
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(10);
		while(System.nanoTime()<deadline){
			for(Thread thread:Thread.getAllStackTraces().keySet()){
				if(!thread.getClass().getName().equals("prok.NcrnaChildProcess$InputPump")){continue;}
				final Object child=field(thread.getClass(), "child").get(thread);
				if(field(child.getClass(), "identity").get(child).equals(prefix+".pid")){return thread;}
			}
			Thread.sleep(10);
		}
		throw new IOException("This invocation's input pump did not become visible: "+prefix);
	}
	private static void awaitFile(String path) throws Exception{
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(10);
		while(!new File(path).isFile() && System.nanoTime()<deadline){Thread.sleep(10);}
		check(new File(path).isFile(), "Child did not reach its bounded control gate: "+path);
	}
	private static boolean ownedCleanupFailure(Throwable problem){
		if(problem==null){return false;}
		if(problem instanceof IllegalStateException && "Owned child cleanup failed".equals(problem.getMessage())){return true;}
		for(Throwable suppressed:problem.getSuppressed()){
			if(ownedCleanupFailure(suppressed)){return true;}
		}
		return false;
	}

	private static Object owned(String path) throws Exception{
		final Class<?> type=Class.forName("prok.NcrnaChildProcess$Owned");
		final Constructor<?> constructor=type.getDeclaredConstructor(String.class, String.class);
		constructor.setAccessible(true);
		return constructor.newInstance(path, UUID.randomUUID().toString());
	}
	private static Field field(Class<?> type, String name) throws Exception{
		final Field field=type.getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}
	private static Throwable invokeFailure(Method method, Object target) throws Exception{
		try{method.invoke(target);return null;}
		catch(InvocationTargetException e){return e.getCause();}
	}
	private static boolean contains(Throwable problem, Throwable target, int depth){
		if(problem==target){return true;}
		if(problem==null || depth>=16){return false;}
		if(contains(problem.getCause(), target, depth+1)){return true;}
		for(Throwable suppressed:problem.getSuppressed()){
			if(contains(suppressed, target, depth+1)){return true;}
		}
		return false;
	}
	private static String read(String path){
		final ByteFile input=RrnaResourceIO.open(path);
		final ByteBuilder data=new ByteBuilder();
		try{for(byte[] line; (line=input.nextLine())!=null;){data.append(line).nl();}}
		finally{check(!input.close(), "Fixture evidence read failed: "+path);}
		return new String(data.toBytes(), StandardCharsets.UTF_8);
	}
	private static void write(String path, String text){
		try(Output out=new Output(FileFormat.testOutput(path, FileFormat.TEXT, null, true, false, false, false))){out.writer.print(text);}
	}
	private static void check(boolean condition, String reason){if(!condition){throw new AssertionError(reason);}}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Observation{
		Throwable problem;
		boolean interrupted;
	}
	private static final class ThrowingInput extends InputStream{
		ThrowingInput(IOException failure_){failure=failure_;}
		@Override public int read(){return -1;}
		@Override public void close() throws IOException{closes++;throw failure;}
		final IOException failure;
		int closes;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String root;
	private final NcrnaChildProcess runner;
	private ByteStreamWriter report;
	private int rows;
	private static final String FIXTURE="dev/ncrnachildfixture.sh";
}
