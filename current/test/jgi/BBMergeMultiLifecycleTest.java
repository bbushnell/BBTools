package test.jgi;

import java.io.IOException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Optional;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.locks.LockSupport;

import jgi.BBMergeMulti;
import parse.Parse;

/**
 * BBMergeMulti child-lifetime and cleanup failure checks.
 * @author Jean
 */
public final class BBMergeMultiLifecycleTest{

	/** Requires phase=baseline or phase=fixed and a new outdir. */
	public static void main(final String[] args) throws Throwable{
		String phase=null, directory=null, root=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){root=pair[1];}
			else{throw new IllegalArgumentException("Unknown BBMergeMulti lifecycle option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline or phase=fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(root!=null && !root.isEmpty(), "Set bbtools=<repo root>.");
		final BBMergeMultiLifecycleTest test=new BBMergeMultiLifecycleTest(phase,
			Paths.get(directory).toAbsolutePath().normalize(),
			Paths.get(root).toAbsolutePath().normalize());
		test.runInterruptedChild();
		test.runMaskedCleanupFailure();
		test.check(test.cases==2, "Expected two lifecycle cases: "+test.cases);
		final String summary="BBMERGE_MULTI_LIFECYCLE_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" checks="+test.checks+
			" failures="+test.failures+"\n";
		Files.write(test.out.resolve("summary.txt"), summary.getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+test.out);}
		System.out.print(summary);
	}

	/** Binds evidence paths and the repository root. */
	private BBMergeMultiLifecycleTest(final String phase_, final Path out_, final Path repo_) throws Exception{
		fixed=phase_.equals("fixed");
		out=out_;
		repo=repo_;
		require(Files.isDirectory(repo.resolve("current")), "Missing BBTools current/ at "+repo);
		require(Files.isRegularFile(Paths.get(System.getProperty("java.home"), "bin", "java")),
			"Missing current JVM java executable");
		require(BBMergeMultiLifecycleTest.class.desiredAssertionStatus(),
			"Run with parent assertions enabled via -ea.");
		require(hasProcessHandle(),
			"Run on JDK 9+; this fixture uses reflective ProcessHandle for child cleanup.");
		Files.createDirectory(out);
		Files.write(out.resolve("cases.tsv"), ("case\tthrown\tinterrupted\tchildAlive"+
			"\tchildObserved\tprimary\tresidualDirectory\tmergedExists\tmatch\n")
			.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Captures JGI-081: waitFor interruption returns while the child is still live. */
	private void runInterruptedChild() throws Throwable{
		cases++;
		final Path work=Files.createDirectory(out.resolve("interrupted-child"));
		final Path release=work.resolve("release.txt");
		final Path probe=work.resolve("probe.txt");
		final Path observed=work.resolve("observed.txt");
		final Path ready=work.resolve("ready.txt");
		final Path exit=work.resolve("exit.txt");
		final String[] childTokens=childTokens("wait", ready, exit, release, probe, observed);
		final Path fakeHome=fakeJavaHome(work, "wait", ready, exit, release, probe, observed);
		Owner owner=null;
		Receipt receipt=null;
		ChildProcess child=null;
		Throwable failure=null;
		try{
			setJavaHome(fakeHome);
			owner=new Owner(argsFor(work));
			owner.setDaemon(true);
			owner.start();
			receipt=awaitReceipt(ready, owner, "ready");
			child=processHandle(receipt.pid, childTokens);
			owner.interrupt();
			join(owner, "owner interrupted by fixture");
			Files.write(probe, "parent-returned\n".getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE_NEW);
			if(!fixed){awaitReceipt(observed, null, "post-return observation");}
			final boolean alive=isAlive(child);
			final boolean observedAfterReturn=Files.isRegularFile(observed);
			final boolean expectAlive=!fixed;
			final boolean expectObserved=!fixed;
			final boolean interrupted=owner.interruptedAfterReturn;
			final boolean thrown=hasMessage(owner.thrown, "Interrupted while waiting") &&
				owner.thrown.getCause() instanceof InterruptedException;
			final boolean match=thrown && interrupted && alive==expectAlive &&
				observedAfterReturn==expectObserved;
			Files.write(out.resolve("cases.tsv"), ("interrupted-child\t"+thrown+"\t"+
				interrupted+"\t"+alive+"\t"+observedAfterReturn+
				"\t\t\t\t"+match+"\n").getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.APPEND);
			check(thrown, "JGI-081 source did not wrap InterruptedException: "+owner.thrown);
			check(interrupted, "JGI-081 source did not restore owner interrupt status");
			check(alive==expectAlive, "JGI-081 child liveness did not match source oracle: "+alive);
			check(observedAfterReturn==expectObserved,
				"JGI-081 post-owner probe observation did not match source oracle: "+observedAfterReturn);
			check(match, "JGI-081 fixed="+fixed+" did not match source oracle");
		}catch(Throwable t){
			failure=t;
		}finally{
			failure=releaseChild(release, ready, exit, receipt, child, owner, childTokens, failure);
			if(owner==null || !owner.isAlive()){restoreJavaHome();}
		}
		if(failure!=null){throw failure;}
	}

	/** Captures JGI-082: first cleanup delete masks pass failure and skips later temps. */
	private void runMaskedCleanupFailure() throws Throwable{
		cases++;
		final Path work=Files.createDirectory(out.resolve("masked-cleanup"));
		final Path release=work.resolve("release.txt");
		final Path probe=work.resolve("probe.txt");
		final Path observed=work.resolve("observed.txt");
		final Path ready=work.resolve("ready.txt");
		final Path exit=work.resolve("exit.txt");
		final String[] childTokens=childTokens("exit7", ready, exit, release, probe, observed);
		final Path fakeHome=fakeJavaHome(work, "exit7", ready, exit, release, probe, observed);
		Owner owner=null;
		Receipt receipt=null;
		ChildProcess child=null;
		Throwable failure=null;
		try{
			setJavaHome(fakeHome);
			owner=new Owner(argsFor(work));
			owner.setDaemon(true);
			owner.start();
			receipt=awaitReceipt(exit, owner, "exit");
			join(owner, "owner after exit7 child");
			child=processHandleIfPresent(receipt.pid, childTokens);
			final boolean deletePrimary=hasMessage(owner.thrown, "Could not delete temporary file");
			final boolean passPrimary=hasMessage(owner.thrown,
				"BBMerge iterative pass 1 failed with exit status 7.");
			final boolean residualDirectory=Files.isDirectory(Paths.get(receipt.residual));
			final boolean mergedExists=Files.exists(Paths.get(receipt.merged));
			final boolean noSuppressed=owner.thrown!=null && owner.thrown.getSuppressed().length==0;
			final boolean suppressedDelete=hasSuppressedMessage(owner.thrown,
				"Could not delete temporary file");
			final boolean primaryMatches=fixed ? passPrimary : deletePrimary;
			final boolean suppressedMatches=fixed ? suppressedDelete : noSuppressed;
			final boolean mergedMatches=fixed ? !mergedExists : mergedExists;
			final boolean match=primaryMatches && suppressedMatches && residualDirectory &&
				mergedMatches;
			Files.write(out.resolve("cases.tsv"), ("masked-cleanup\t"+
				(owner.thrown!=null)+"\t"+owner.interruptedAfterReturn+"\tfalse\tfalse\t"+
				message(owner.thrown)+"\t"+residualDirectory+"\t"+mergedExists+"\t"+
				match+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
			check(primaryMatches, "JGI-082 primary failure did not match source oracle: "+
				owner.thrown);
			check(suppressedMatches, "JGI-082 suppressed cleanup failures did not match source oracle");
			check(residualDirectory, "JGI-082 residual blocker was not left as a directory");
			check(mergedMatches, "JGI-082 merged temp existence did not match source oracle: "+
				mergedExists);
			check(match, "JGI-082 fixed="+fixed+" did not match source oracle");
		}catch(Throwable t){
			failure=t;
		}finally{
			failure=releaseChild(release, ready, exit, receipt, child, owner, childTokens, failure);
			if(owner==null || !owner.isAlive()){restoreJavaHome();}
		}
		if(failure!=null){throw failure;}
	}

	/** Produces a fake java.home whose bin/java launches the bounded Java child. */
	private Path fakeJavaHome(final Path work, final String mode, final Path ready,
			final Path exit, final Path release, final Path probe, final Path observed)
			throws Exception{
		final Path home=Files.createDirectory(work.resolve("fake-home"));
		final Path bin=Files.createDirectory(home.resolve("bin"));
		final Path lib=Files.createDirectory(home.resolve("lib"));
		final Path java=bin.resolve("java");
		final String realHome=System.getProperty("java.home");
		final Path realJava=Paths.get(realHome, "bin", "java");
		final Path realModules=Paths.get(realHome, "lib", "modules");
		require(Files.isRegularFile(realModules), "Missing current JVM modules image: "+realModules);
		Files.createSymbolicLink(lib.resolve("modules"), realModules);
		final String script="#!/usr/bin/env bash\nexec "+quote(realJava.toString())+
			" -Dfile.encoding=UTF-8 -ea -Xmx128m -Xms32m"+
			" -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError"+
			" -cp "+quote(System.getProperty("java.class.path"))+" "+
			quote(Child.class.getName())+" mode="+quote(mode)+
			" ready="+quote(ready.toString())+" exit="+quote(exit.toString())+
			" release="+quote(release.toString())+" probe="+quote(probe.toString())+
			" observed="+quote(observed.toString())+" \"$@\"\n";
		Files.write(java, script.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		require(java.toFile().setExecutable(true), "Could not make fake java executable: "+java);
		return home;
	}

	/** Command tokens that prove a recovered child belongs to this specific case. */
	private static String[] childTokens(final String mode, final Path ready, final Path exit,
			final Path release, final Path probe, final Path observed){
		return new String[]{Child.class.getName(), "mode="+mode, "ready="+ready.toString(),
			"exit="+exit.toString(), "release="+release.toString(), "probe="+probe.toString(),
			"observed="+observed.toString()};
	}

	/** Installs fake java.home after all pre-mutation fixture setup has succeeded. */
	private void setJavaHome(final Path home){
		final String old=System.setProperty("java.home", home.toString());
		previousJavaHomes.add(old);
	}

	/** Restores java.home after a fake-home case. */
	private void restoreJavaHome(){
		if(previousJavaHomes.isEmpty()){return;}
		final String previous=previousJavaHomes.remove(previousJavaHomes.size()-1);
		if(previous==null){System.clearProperty("java.home");}
		else{System.setProperty("java.home", previous);}
	}

	/** Returns real BBMergeMulti arguments that spawn one child on pass 0. */
	private static String[] argsFor(final Path work) throws Exception{
		final Path in=work.resolve("in.fq");
		Files.write(in, "@r\nACGT\n+\nIIII\n".getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
		return new String[]{"k=0,1", "in="+in, "out="+work.resolve("merged.fq"),
			"outu="+work.resolve("unmerged.fq"), "tmpdir="+work.resolve("tmp"), "ow=t"};
	}

	/** Waits for a receipt written by the external child. */
	private static Receipt awaitReceipt(final Path receipt, final Owner owner, final String name)
			throws Exception{
		final long deadline=System.nanoTime()+TIMEOUT_NANOS;
		while(System.nanoTime()<deadline){
			if(Files.isRegularFile(receipt)){return new Receipt(receipt);}
			if(owner!=null && owner.done()){
				if(Files.isRegularFile(receipt)){return new Receipt(receipt);}
				throw ownerReturnedBeforeReceipt(name, owner.thrown);
			}
			LockSupport.parkNanos(PARK_NANOS);
		}
		throw new AssertionError("Timed out waiting for child "+name+" receipt: "+receipt);
	}

	/** Releases a waiting child and forcibly terminates it if it fails to exit. */
	private static void releaseChild(final Path release, final Path ready, final Path exit,
			Receipt receipt, ChildProcess child, final Owner owner, final String[] expectedTokens)
			throws Exception{
		if(owner==null && receipt==null && child==null){
			Files.write(release, "release\n".getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
			return;
		}
		if(owner!=null && owner.isAlive()){owner.interrupt();}
		Files.write(release, "release\n".getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
		final long deadline=System.nanoTime()+TIMEOUT_NANOS;
		while(System.nanoTime()<deadline){
			if(receipt==null){receipt=readReceiptIfPresent(ready);}
			if(receipt==null){receipt=readReceiptIfPresent(exit);}
			if(child==null && receipt!=null){child=processHandleIfPresent(receipt.pid, expectedTokens);}
			if(child!=null && !isAlive(child)){join(owner, "owner after child exit"); return;}
			if(child==null && Files.isRegularFile(exit)){join(owner, "owner after child exit"); return;}
			LockSupport.parkNanos(PARK_NANOS);
		}
		require(child!=null, "Could not observe child process before cleanup timeout.");
		if(isAlive(child)){destroyForcibly(child);}
		final long killDeadline=System.nanoTime()+TIMEOUT_NANOS;
		while(isAlive(child) && System.nanoTime()<killDeadline){
			LockSupport.parkNanos(PARK_NANOS);
		}
		if(isAlive(child)){throw new AssertionError("Could not terminate child process");}
		join(owner, "owner after forcible child termination");
	}

	/** Releases a child without discarding any earlier fixture failure. */
	private static Throwable releaseChild(final Path release, final Path ready, final Path exit,
			final Receipt receipt, final ChildProcess child, final Owner owner,
			final String[] expectedTokens, final Throwable failure){
		try{
			releaseChild(release, ready, exit, receipt, child, owner, expectedTokens);
			return failure;
		}catch(Throwable t){
			if(failure!=null){failure.addSuppressed(t); return failure;}
			return t;
		}
	}

	/** Preserves the owner failure that prevented a receipt from arriving. */
	private static AssertionError ownerReturnedBeforeReceipt(final String name,
			final Throwable ownerFailure){
		final AssertionError error=new AssertionError("Owner returned before child "+name+
			" receipt: "+ownerFailure);
		if(ownerFailure!=null){error.initCause(ownerFailure);}
		return error;
	}

	/** Reads a receipt if a complete atomic receipt is available. */
	private static Receipt readReceiptIfPresent(final Path path) throws Exception{
		return Files.isRegularFile(path) ? new Receipt(path) : null;
	}

	/** Reflective JDK9+ ProcessHandle liveness check; source stays Java-8-safe. */
	private static boolean isAlive(final ChildProcess process) throws Exception{
		return process.isAlive();
	}

	/** Reflective JDK9+ ProcessHandle forcible termination after an ownership recheck. */
	private static void destroyForcibly(final ChildProcess process) throws Exception{
		process.verify();
		Class.forName("java.lang.ProcessHandle").getMethod("destroyForcibly")
			.invoke(process.handle);
	}

	/** Looks up a ProcessHandle reflectively to preserve Java 8 source/API use. */
	private static ChildProcess processHandle(final long pid, final String[] expectedTokens)
			throws Exception{
		final ChildProcess handle=processHandleIfPresent(pid, expectedTokens);
		if(handle==null){throw new AssertionError("Could not observe child process "+pid);}
		return handle;
	}

	/** Looks up a ProcessHandle reflectively, returning null for an already-dead child. */
	private static ChildProcess processHandleIfPresent(final long pid,
			final String[] expectedTokens) throws Exception{
		final Class<?> clazz=Class.forName("java.lang.ProcessHandle");
		final Method of=clazz.getMethod("of", long.class);
		final Optional<?> optional=(Optional<?>)of.invoke(null, pid);
		if(!optional.isPresent()){return null;}
		return new ChildProcess(optional.get(), expectedTokens);
	}

	/** True when this VM has JDK9+ ProcessHandle for reflective child cleanup. */
	private static boolean hasProcessHandle(){
		try{
			Class.forName("java.lang.ProcessHandle");
			return true;
		}catch(ClassNotFoundException e){
			return false;
		}
	}

	/** Joins a thread with a finite timeout. */
	private static void join(final Thread thread, final String name) throws Exception{
		if(thread==null){return;}
		final long deadline=System.nanoTime()+TIMEOUT_NANOS;
		while(thread.isAlive() && System.nanoTime()<deadline){thread.join(10);}
		if(thread.isAlive()){throw new AssertionError("Timed out joining "+name);}
	}

	/** Records one evidence check without hiding later evidence after a failure. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}

	/** Hard setup requirement. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Single-quotes a string for the generated bash stub. */
	private static String quote(final String s){
		return "'"+s.replace("'", "'\\''")+"'";
	}

	/** True when a throwable carries this text. */
	private static boolean hasMessage(final Throwable t, final String text){
		return t!=null && t.getMessage()!=null && t.getMessage().contains(text);
	}

	/** True when a direct suppressed throwable carries this text. */
	private static boolean hasSuppressedMessage(final Throwable t, final String text){
		if(t==null){return false;}
		for(final Throwable suppressed : t.getSuppressed()){
			if(hasMessage(suppressed, text)){return true;}
		}
		return false;
	}

	/** Makes a throwable message safe for the TSV evidence file. */
	private static String message(final Throwable t){
		return t==null || t.getMessage()==null ? "" : t.getMessage().replace('\t', ' ');
	}

	/** Thread that calls the real BBMergeMulti entry point. */
	private static final class Owner extends Thread{
		Owner(final String[] args_){args=args_;}

		@Override
		public void run(){
			try{
				BBMergeMulti.main(args);
			}catch(Throwable t){
				thrown=t;
			}finally{
				interruptedAfterReturn=Thread.currentThread().isInterrupted();
				done.countDown();
			}
		}

		boolean done(){return done.getCount()==0;}

		private final CountDownLatch done=new CountDownLatch(1);
		private final String[] args;
		private volatile Throwable thrown;
		private volatile boolean interruptedAfterReturn;
	}

	/** Verified direct child ProcessHandle with this case's command arguments. */
	private static final class ChildProcess{
		ChildProcess(final Object handle_, final String[] expectedTokens_) throws Exception{
			handle=handle_;
			expectedTokens=expectedTokens_;
			pid=pid(handle);
			verify();
		}

		/** Fails closed if the handle no longer proves it owns the fixture child. */
		void verify() throws Exception{
			require(pid(handle)==pid, "ProcessHandle PID changed before cleanup signal.");
			final Object self=Class.forName("java.lang.ProcessHandle").getMethod("current")
				.invoke(null);
			final Optional<?> parent=parent();
			if(!parent.isPresent() || pid(parent.get())!=pid(self)){
				require(!isAlive(),
					"Refusing to signal process without direct fixture-JVM parentage: "+pid);
				return;
			}
			final Optional<?> commandLine=commandLine();
			if(!commandLine.isPresent()){
				require(!isAlive(), "Refusing to signal process with hidden command line: "+pid);
				return;
			}
			final String command=commandLine.get().toString();
			for(final String token : expectedTokens){
				if(command.indexOf(token)<0){
					require(!isAlive(),
						"Refusing to signal process missing command token '"+token+"': "+
						command);
					return;
				}
			}
		}

		/** Reflectively returns this handle's direct parent. */
		private Optional<?> parent() throws Exception{
			final Method method=Class.forName("java.lang.ProcessHandle").getMethod("parent");
			return (Optional<?>)method.invoke(handle);
		}

		/** Reflectively returns the JDK21 command line. */
		private Optional<?> commandLine() throws Exception{
			final Method infoMethod=Class.forName("java.lang.ProcessHandle").getMethod("info");
			final Object info=infoMethod.invoke(handle);
			final Method commandLine=Class.forName("java.lang.ProcessHandle$Info")
				.getMethod("commandLine");
			return (Optional<?>)commandLine.invoke(info);
		}

		/** Reflectively checks whether this exact cached handle is still live. */
		private boolean isAlive() throws Exception{
			final Method method=Class.forName("java.lang.ProcessHandle").getMethod("isAlive");
			return ((Boolean)method.invoke(handle)).booleanValue();
		}

		/** Reflectively returns a ProcessHandle pid. */
		private static long pid(final Object processHandle) throws Exception{
			final Class<?> clazz=Class.forName("java.lang.ProcessHandle");
			return ((Long)clazz.getMethod("pid").invoke(processHandle)).longValue();
		}

		private final Object handle;
		private final String[] expectedTokens;
		private final long pid;
	}

	/** Parsed child receipt with pid and BBMergeMulti-generated temp paths. */
	private static final class Receipt{
		Receipt(final Path path) throws Exception{
			String pidString=null, merged_=null, residual_=null;
			for(final String line : Files.readAllLines(path, StandardCharsets.UTF_8)){
				final String[] pair=Parse.splitOnFirst(line, '=');
				if(pair.length<2){continue;}
				if(pair[0].equals("pid")){pidString=pair[1];}
				else if(pair[0].equals("out")){merged_=pair[1];}
				else if(pair[0].equals("outu")){residual_=pair[1];}
			}
			require(pidString!=null, "Receipt missing pid: "+path);
			require(merged_!=null, "Receipt missing out: "+path);
			require(residual_!=null, "Receipt missing outu: "+path);
			pid=Long.parseLong(pidString);
			merged=merged_;
			residual=residual_;
		}

		final long pid;
		final String merged, residual;
	}

	/** Bounded fake child launched through BBMergeMulti's generated java command. */
	public static final class Child{

		/** Executes wait-for-release or cleanup-failure behavior. */
		public static void main(final String[] args) throws Exception{
			String mode=null, ready=null, exit=null, release=null, probe=null, observed=null;
			String merged=null, residual=null;
			for(final String arg : args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length<2){continue;}
				if(pair[0].equals("mode")){mode=pair[1];}
				else if(pair[0].equals("ready")){ready=pair[1];}
				else if(pair[0].equals("exit")){exit=pair[1];}
				else if(pair[0].equals("release")){release=pair[1];}
				else if(pair[0].equals("probe")){probe=pair[1];}
				else if(pair[0].equals("observed")){observed=pair[1];}
				else if(pair[0].equals("out")){merged=pair[1];}
				else if(pair[0].equals("outu")){residual=pair[1];}
			}
			require(mode!=null && ready!=null && exit!=null && release!=null,
				"Missing child control arguments");
			if("exit7".equals(mode)){
				blockResidual(residual);
				touch(merged);
				writeReceipt(Paths.get(exit), merged, residual);
				System.exit(7);
			}
			writeReceipt(Paths.get(ready), merged, residual);
			final Path releasePath=Paths.get(release);
			final Path probePath=Paths.get(probe);
			final Path observedPath=Paths.get(observed);
			final long deadline=System.nanoTime()+TIMEOUT_NANOS;
			while(System.nanoTime()<deadline && !Files.exists(releasePath)){
				if(Files.exists(probePath) && !Files.exists(observedPath)){
					writeReceipt(observedPath, merged, residual);
				}
				LockSupport.parkNanos(PARK_NANOS);
			}
			writeReceipt(Paths.get(exit), merged, residual);
		}

		/** Replaces the residual temp file with a nonempty directory. */
		private static void blockResidual(final String residual) throws IOException{
			final Path path=Paths.get(residual);
			Files.deleteIfExists(path);
			Files.createDirectory(path);
			Files.write(path.resolve("block.txt"), "block\n".getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE_NEW);
		}

		/** Creates or truncates the merged temp file. */
		private static void touch(final String merged) throws IOException{
			Files.write(Paths.get(merged), "merged\n".getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
		}

		/** Writes pid and generated paths for the parent fixture. */
		private static void writeReceipt(final Path path, final String merged,
				final String residual) throws Exception{
			require(merged!=null && residual!=null,
				"Missing generated child out/outu arguments for "+path);
			final ArrayList<String> lines=new ArrayList<String>();
			lines.add("pid="+pid());
			if(merged!=null){lines.add("out="+merged);}
			if(residual!=null){lines.add("outu="+residual);}
			final Path tmp=Files.createTempFile(path.toAbsolutePath().getParent(),
				path.getFileName().toString(), ".tmp");
			try{
				Files.write(tmp, lines, StandardCharsets.UTF_8, StandardOpenOption.CREATE,
					StandardOpenOption.TRUNCATE_EXISTING);
				Files.move(tmp, path, StandardCopyOption.ATOMIC_MOVE,
					StandardCopyOption.REPLACE_EXISTING);
			}finally{
				Files.deleteIfExists(tmp);
			}
		}

		/** Reflective JDK9+ current pid; no direct ProcessHandle source dependency. */
		private static long pid() throws Exception{
			final Class<?> clazz=Class.forName("java.lang.ProcessHandle");
			final Object handle=clazz.getMethod("current").invoke(null);
			return ((Long)clazz.getMethod("pid").invoke(handle)).longValue();
		}
	}

	private final boolean fixed;
	private final Path out, repo;
	private final ArrayList<String> previousJavaHomes=new ArrayList<String>();
	private int cases, checks, failures;
	private static final long TIMEOUT_NANOS=10L*1000L*1000L*1000L;
	private static final long PARK_NANOS=10L*1000L*1000L;
}
