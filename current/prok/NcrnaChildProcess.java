package prok;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.IdentityHashMap;
import java.util.UUID;
import java.util.concurrent.TimeUnit;
import fileIO.ByteFile;
import fileIO.ReadWrite;

/** Runs a native CLI stage in its own Unix process group within the current job.
 * Logs go to files; optional shards stream through native decompression without a
 * concatenated bulk copy. Cleanup targets only this invocation's adopted group.
 * Immutable installation/CWD settings permit concurrent runs with distinct stable
 * evidence prefixes; same-prefix reservation is not atomic. Requires bash, setsid,
 * ps, kill and /usr/bin/time. Source APIs remain Java-8-compatible.
 * Created for the earlier pipeline's logged native stages. The shipping model
 * tools now use NcrnaProcess; this owner remains necessary support for the
 * preserved standalone regeneration mode. Keep it while that mode needs its
 * file-based stage logs and launch contract, pending the mode's disposition.
 * @author Raiden */
public final class NcrnaChildProcess{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Binds an installation and captures user.dir as the child working directory.
	 * @param installation Non-null existing BBTools directory containing the child launcher. */
	public NcrnaChildProcess(File installation){
		this(installation, new File(System.getProperty("user.dir")));
	}
	/** Binds existing installation and working directories without taking file ownership.
	 * Relative resource paths used by children resolve in directory, not installation.
	 * @param installation Non-null existing BBTools directory with dev/ncrnachild.sh.
	 * @param directory Non-null existing child working directory. */
	public NcrnaChildProcess(File installation, File directory){
		root=installation.getAbsoluteFile();
		launcher=new File(root, "dev/ncrnachild.sh");
		workingDirectory=directory.getAbsoluteFile();
		require(root.isDirectory() && launcher.isFile(), "Child stages require one complete BBTools installation with dev/ncrnachild.sh");
		require(workingDirectory.isDirectory(), "Child working directory must already exist");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs one stage without a startup callback; see the callback overload for its contract.
	 * @return Exit code, elapsed wall seconds and decompressed bytes accepted by stdin writes.
	 * @throws IOException When launch, child execution or streaming fails.
	 * @throws InterruptedException When the calling thread is interrupted during execution. */
	public Result run(String wrapper, String[] args, String[] inputShards, File prefix, int threads, int heapMiB)
	throws IOException, InterruptedException{
		return run(wrapper, args, inputShards, prefix, threads, heapMiB, null);
	}
	/** Runs one stage, then completes owned cleanup before returning its result.
	 * Registers a shutdown hook before launch, invokes startup before nonce-bound
	 * identity adoption/release, and terminates the owned group before joining a
	 * potentially blocked input producer. Nonzero exit, streaming failure or failed
	 * cleanup is failure. Partial evidence remains available after an exception.
	 * The original stage throwable remains primary; subsequent cleanup
	 * failures are suppressed beneath it. Source-input close failures are reported,
	 * including when closing is the only failed operation. A failed cleanup never
	 * returns a successful Result.
	 * @param wrapper Existing script resolving within the bound installation.
	 * @param args Non-null vector of intact non-null scalar arguments without TSV/newline/NUL delimiters.
	 * @param inputShards Ordered caller-owned input paths, or null to close child stdin immediately.
	 * @param prefix Non-null fresh evidence prefix beneath an existing directory.
	 * @param threads Positive compute-thread budget passed to the child and JVM.
	 * @param heapMiB Positive maximum child heap in MiB.
	 * @param startup Optional callback before ownership adoption/workload release.
	 * @return A successful result only after the finally block completes.
	 * @throws IOException When launch, child execution, identity or streaming fails.
	 * @throws InterruptedException When prelaunch checking detects cancellation or execution is
	 * interrupted. Prelaunch checking preserves the flag; cleanup restores interrupts it captures.
	 * @throws IllegalArgumentException When input validation or a cleanup invariant fails.
	 * @throws IllegalStateException When owned-process cleanup fails. */
	public Result run(String wrapper, String[] args, String[] inputShards, File prefix, int threads, int heapMiB, StartupAction startup)
	throws IOException, InterruptedException{
		require(threads>0 && heapMiB>0 && args!=null, "Each child needs positive CPU/heap limits and an argument vector; threads="+threads+", heapMiB="+heapMiB+", argsPresent="+(args!=null));
		final File script=new File(root, wrapper).getCanonicalFile();
		require(script.toPath().startsWith(root.getCanonicalFile().toPath()) && script.isFile(), "Only a wrapper in the bound installation may run: "+wrapper);
		final String base=prefix.getAbsolutePath();
		RrnaResourceIO.fresh(base, ".pid", ".pid.pending", ".pid.ready", ".pid.ready.pending", ".pid.cancelled", ".stdout", ".stderr", ".time.tsv");
		require(prefix.getAbsoluteFile().getParentFile().isDirectory(), "Child evidence directory must already exist");
		if(inputShards!=null){
			for(String input:inputShards){
				require(input!=null && new File(input).isFile(), "Missing child input shard: "+input);
			}
		}
		final String nonce=UUID.randomUUID().toString();
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", launcher.getPath(), base+".pid", nonce,
			"/usr/bin/time", "-f", "%U\t%S\t%e", "-o", base+".time.tsv", "bash", script.getPath()));
		for(String arg:args){
			require(arg!=null && arg.indexOf('\0')<0 && arg.indexOf('\n')<0 && arg.indexOf('\r')<0 && arg.indexOf('\t')<0,
				"Child arguments must be intact scalar values without control delimiters");
			command.add(arg);
		}
		command.add("t="+threads);
		command.add("-Xmx"+heapMiB+"m");
		command.add("-Xms"+Math.min(64, heapMiB)+"m");
		command.add("-eoom");
		checkCancellation(base);
		final ProcessBuilder pb=new ProcessBuilder(command).directory(workingDirectory);
		pb.redirectOutput(new File(base+".stdout"));
		pb.redirectError(new File(base+".stderr"));
		pb.environment().remove("_JAVA_OPTIONS");
		pb.environment().remove("JDK_JAVA_OPTIONS");
		pb.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount="+threads);
		final Owned child=new Owned(base+".pid", nonce);
		final Thread shutdown=new Thread(child::stop, "ncRNA-child-cleanup");
		InputPump pump=null;
		int exit=-1;
		final long start=System.nanoTime();
		Result result=null;
		Throwable failure=null;
		Runtime.getRuntime().addShutdownHook(shutdown);
		try{
			final Process process=child.start(pb);
			if(startup!=null){
				startup.started();
			}
			child.readIdentity();
			child.release();
			if(inputShards==null){
				process.getOutputStream().close();
			}else{
				pump=new InputPump(child, inputShards);
				pump.start();
			}
			exit=process.waitFor();
			// Retain known stage failures before owned cleanup can throw.
			if(exit!=0){
				failure=new IOException("Child stage exited "+exit+": "+wrapper+"; inspect "+base+".stderr");
			}
			if(pump!=null){failure=appendFailure(failure, pump.failure);}
			// A descendant may keep stdin open after the leader exits. Close the
			// owned group before joining a producer that could be blocked in write.
			child.stop();
			if(pump!=null){
				pump.closeInput();
				pump.join(5000);
				require(!pump.isAlive(), "Child input producer survived owned-group cleanup");
			}
			if(pump!=null){failure=appendFailure(failure, pump.failure);}
			if(failure==null){
				result=new Result(exit, (System.nanoTime()-start)*1e-9, pump==null ? 0 : pump.bytes);
			}
		}catch(Throwable problem){
			failure=appendFailure(failure, problem);
		}finally{
			failure=finish(child, pump, shutdown, failure);
		}
		if(failure!=null){throwFailure(failure);}
		assert(result!=null) : "A failure-free stage sets Result before mandatory finalization; exit="+exit;
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Attempts every independent cleanup phase and retains failures in phase order.
	 * Group termination precedes pump joining. The producer alone writes its streaming
	 * failure; parent/producer cleanup failures are collected separately after the join.
	 * Captured interrupts are restored even when a cleanup operation fails. */
	private static Throwable finish(Owned child, InputPump pump, Thread shutdown, Throwable failure){
		final Interrupts interrupts=new Interrupts();
		try{
			try{child.stop();}
			catch(Throwable problem){failure=appendFailure(failure, problem);}
			finally{interrupts.capture();}
			if(pump!=null){
				try{pump.closeInput();}
				catch(Throwable problem){failure=appendFailure(failure, problem);}
				try{pump.interrupt();}
				catch(Throwable problem){failure=appendFailure(failure, problem);}
				try{
					interrupts.join(pump, 5000);
					require(!pump.isAlive(), "Child input reader did not stop after its owned process terminated");
				}catch(Throwable problem){failure=appendFailure(failure, problem);}
				try{failure=pump.appendFailures(failure);}
				catch(Throwable problem){failure=appendFailure(failure, problem);}
			}
			try{Runtime.getRuntime().removeShutdownHook(shutdown);}
			catch(IllegalStateException shuttingDown){}
			catch(Throwable problem){failure=appendFailure(failure, problem);}
			interrupts.capture();
			return failure;
		}finally{interrupts.restore();}
	}

	/** Retains the first throwable; later distinct failures stay visible without self-suppression. */
	private static Throwable appendFailure(Throwable first, Throwable next){
		if(next==null){return first;}
		if(first==null){return next;}
		if(!containsFailure(first, next)){first.addSuppressed(next);}
		return first;
	}

	/** Identity traversal is error-path-only; cycles and already nested reports do not duplicate a failure. */
	private static boolean containsFailure(Throwable root, Throwable target){
		if(root==target){return true;}
		final ArrayList<Throwable> pending=new ArrayList<Throwable>();
		final IdentityHashMap<Throwable, Boolean> seen=new IdentityHashMap<Throwable, Boolean>();
		pending.add(root);
		for(int i=0; i<pending.size(); i++){
			final Throwable current=pending.get(i);
			if(current==target){return true;}
			if(seen.put(current, Boolean.TRUE)!=null){continue;}
			if(current.getCause()!=null){pending.add(current.getCause());}
			for(Throwable suppressed:current.getSuppressed()){pending.add(suppressed);}
		}
		return false;
	}

	/** Propagates supported stage/cleanup throwable types without replacing their identity. */
	private static void throwFailure(Throwable failure) throws IOException, InterruptedException{
		assert(failure!=null) : "Only captured stage/cleanup failures reach the terminal rethrow";
		if(failure instanceof IOException){throw (IOException)failure;}
		if(failure instanceof InterruptedException){throw (InterruptedException)failure;}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}
		throw new IllegalStateException("Unexpected checked child-stage failure", failure);
	}

	/** Detects pending startup cancellation without clearing the caller's interrupt flag. */
	private static void checkCancellation(String path) throws InterruptedException{
		if(Thread.currentThread().isInterrupted()){
			throw new InterruptedException("Child startup cancelled before launch: "+path);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable successful-stage observations; neither wall duration is CPU time. */
	public static final class Result{
		/** Stores the successful body sample; final cleanup occurs after this sample. */
		Result(int code, double seconds, long bytes){
			exitCode=code;
			elapsedSeconds=seconds;
			stdinBytes=bytes;
		}
		/** Native wrapper exit code, zero for a returned successful result. */
		public final int exitCode;
		/** Wall seconds from just before launch through body cleanup/pump join.
		 * Excludes final cleanup, which must still succeed before run returns. */
		public final double elapsedSeconds;
		/** Decompressed bytes accepted by successful buffered stdin writes.
		 * Not compressed size, unique bases, flushed bytes or proof of child consumption. */
		public final long stdinBytes;
	}
	/** Runs before identity adoption/release; failures still require owned cleanup. */
	public interface StartupAction{void started() throws IOException, InterruptedException;}
	/** Owns one nonce-bound child group. Synchronized lifecycle methods serialize
	 * adoption, release and cleanup; a concurrent stop waits for the first cleanup. */
	private static final class Owned{
		/** Binds this invocation's fresh identity pathname and expected nonce. */
		Owned(String identity_, String nonce_){
			identity=identity_;
			nonce=nonce_;
		}
		/** Starts this invocation once unless cancelled; returns its owned leader. */
		synchronized Process start(ProcessBuilder builder) throws IOException{
			if(stopped){
				throw new IOException("Child startup was cancelled before launch");
			}
			assert(process==null) : "run creates one Owned instance per launch; replacing its process would lose cleanup ownership; identity="+identity
				+", group="+group+", stopped="+stopped+", previousProcessAlive="+(process!=null && process.isAlive());
			process=builder.start();
			return process;
		}
		/** Waits for the launcher record while its leader lives, then validates/adopts it. */
		void readIdentity() throws IOException, InterruptedException{
			final File file=new File(identity);
			while(!file.isFile() && process.isAlive()){
				Thread.sleep(10);
			}
			if(!file.isFile()){
				throw new IOException("Child failed before recording its owned process group: "+identity);
			}
			adoptIdentity();
		}
		/** Reads exactly one nonce/positive-group row and closes the native reader.
		 * An already adopted group is retained; malformed or unreadable identity fails. */
		private synchronized void adoptIdentity() throws IOException{
			if(group>1){return;}
			final ByteFile in=RrnaResourceIO.open(identity);
			final byte[] row;
			try{
				row=in.nextLine();
				if(row==null || in.nextLine()!=null){
					throw new IOException("Invalid child identity framing: "+identity);
				}
			}finally{
				if(in.close()){
					//TODO: Probable bug - a close failure here masks an earlier read or
					//identity-framing exception. Preserve the primary error and suppress cleanup.
					throw new IOException("Could not read child process identity");
				}
			}
			final String value=new String(row, StandardCharsets.UTF_8);
			final String prefix=nonce+"\t";
			if(!value.startsWith(prefix)){
				throw new IOException("Child identity nonce differs: "+identity);
			}
			final String numeric=value.substring(prefix.length());
			if(!numeric.matches("[1-9][0-9]{0,9}")){
				throw new IOException("Child process group must be one positive decimal PID");
			}
			final long parsed=Long.parseLong(numeric);
			if(parsed<=1 || parsed>Integer.MAX_VALUE){
				throw new IOException("Unsafe child process group: "+numeric);
			}
			group=(int)parsed;
		}
		/** Atomically publishes this nonce's ready marker only after group adoption. */
		synchronized void release() throws IOException{
			if(stopped || group<=1 || process==null){
				throw new IOException("Child workload cannot start before cleanup ownership is established");
			}
			final java.nio.file.Path pending=new File(identity+".ready.pending").toPath();
			Files.write(pending, (nonce+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
			Files.move(pending, new File(identity+".ready").toPath(), StandardCopyOption.ATOMIC_MOVE);
		}
		/** Cancels startup, sends TERM then bounded KILL as needed, and reaps the leader.
		 * Completed stops are idempotent; concurrent stops wait on this monitor.
		 * Captured interrupts are restored. Saved marker failures remain errors even
		 * if no process was started; stopped still prevents any later launch. */
		synchronized void stop(){
			if(stopped){return;}
			final Interrupts interrupts=new Interrupts();
			try{
				IOException markerFailure=null;
				try{
					Files.write(new File(identity+".cancelled").toPath(), new byte[0], StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
				}catch(java.nio.file.FileAlreadyExistsException alreadyCancelled){
				}catch(IOException e){
					markerFailure=e;
				}
				if(process==null){
					stopped=true;
					if(markerFailure!=null){
						throw new IOException("Startup was cancelled but its cancellation marker could not be written", markerFailure);
					}
					return;
				}
				if(group<=1 && new File(identity).isFile()){
					adoptIdentity();
				}
				if(group>1){
					if(signal("0", interrupts)==0){
						signal("TERM", interrupts);
						final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(2);
						while(signal("0", interrupts)==0 && System.nanoTime()<deadline){
							interrupts.sleep(20);
						}
						// The leader may exit before a descendant that ignores TERM.
						// Group cleanup must not be conditional on leader liveness.
						if(signal("0", interrupts)==0){
							signal("KILL", interrupts);
						}
					}
				}else{
					process.destroy();
				}
				if(!interrupts.waitFor(process, 2000)){
					process.destroyForcibly();
					require(interrupts.waitFor(process, 2000), "Owned child leader did not terminate");
				}
				stopped=true;
				if(markerFailure!=null){
					throw new IOException("Child was terminated but its cancellation marker could not be written", markerFailure);
				}
			}catch(IOException e){
				if(process!=null){
					process.destroyForcibly();
				}
				throw new IllegalStateException("Owned child cleanup failed", e);
			}finally{
				interrupts.restore();
			}
		}
		/** Signals only the positive adopted group; returns the kill utility's exit code.
		 * Signal 0 probes existence. The helper process itself has a bounded wait. */
		private int signal(String name, Interrupts interrupts) throws IOException{
			assert(group>1) : "Signals may target only the positive owned session written by the nonce-bound launcher; signal="+name+", group="+group+", identity="+identity;
			final Process kill=new ProcessBuilder("kill", "-"+name, "--", "-"+group).redirectOutput(new File("/dev/null")).redirectError(new File("/dev/null")).start();
			if(!interrupts.waitFor(kill, 2000)){
				kill.destroyForcibly();
				throw new IOException("Owned-group signal command did not complete");
			}
			return kill.exitValue();
		}
		/** Owned launcher leader, published to the pump after start. */
		volatile Process process;
		/** Fresh identity-record path and invocation-specific ownership token. */
		final String identity, nonce;
		/** Positive adopted Unix group; zero until its nonce has been checked. */
		volatile int group;
		/** Completed stop state, accessed under this object's lifecycle monitor. */
		boolean stopped;
	}
	/** Cleanup consumes repeated interrupts while finishing ownership work, then restores the flag. */
	private static final class Interrupts{
		/** Captures and clears an interrupt so bounded cleanup can continue. */
		boolean interrupted=Thread.interrupted();
		/** Accumulates another interrupt without losing an earlier one. */
		void capture(){interrupted|=Thread.interrupted();}
		/** Waits up to millis milliseconds; returns false at deadline, consuming interrupts. */
		boolean waitFor(Process process, long millis){
			final long deadline=System.nanoTime()+TimeUnit.MILLISECONDS.toNanos(millis);
			while(process.isAlive()){
				final long left=deadline-System.nanoTime();
				if(left<=0){return false;}
				try{
					if(process.waitFor(left, TimeUnit.NANOSECONDS)){return true;}
				}catch(InterruptedException e){
					interrupted=true;
				}
			}
			return true;
		}
		/** Joins for at most millis milliseconds; caller must still inspect isAlive(). */
		void join(Thread thread, long millis){
			final long deadline=System.nanoTime()+TimeUnit.MILLISECONDS.toNanos(millis);
			while(thread.isAlive()){
				final long left=deadline-System.nanoTime();
				if(left<=0){return;}
				try{
					thread.join(Math.max(1, TimeUnit.NANOSECONDS.toMillis(left)));
				}catch(InterruptedException e){
					interrupted=true;
				}
			}
		}
		/** Sleeps in milliseconds; an interrupt is remembered and may shorten the sleep. */
		void sleep(long millis){
			try{
				Thread.sleep(millis);
			}catch(InterruptedException e){
				interrupted=true;
			}
		}
		/** Restores a captured or newly arrived interrupt on the calling thread. */
		void restore(){
			if(interrupted || Thread.interrupted()){
				Thread.currentThread().interrupt();
			}
		}
	}
	/** Daemon producer with its own buffer. Reads ordered borrowed shards and owns
	 * their open streams plus child stdin; caller files are never modified/deleted.
	 * Counts bytes after buffered writes, even if a later flush fails. */
	private static final class InputPump extends Thread{
		/** Clones the shard-path array; the underlying input files remain caller-owned. */
		InputPump(Owned child_, String[] files_){
			super("ncRNA-child-input");
			child=child_;
			files=files_.clone();
			setDaemon(true);
		}
		/** Streams decompressed shards and records failure for the parent; closes input on exit. */
		@Override public void run(){
			final byte[] buffer=new byte[1<<16];
			try(OutputStream out=new BufferedOutputStream(child.process.getOutputStream(), buffer.length)){
				for(String file:files){
					try(InputStream in=ReadWrite.getInputStream(file, true, false)){
						current=in;
						for(int n; (n=in.read(buffer))>=0;){
							if(n>0){
								out.write(buffer, 0, n);
								bytes+=n;
							}
						}
						current=null;
					}
				}
			}catch(Throwable problem){
				failure=problem;
				try{child.stop();}
				catch(Throwable cleanup){recordCleanup(cleanup);}
			}finally{
				closeInput();
			}
		}
		/** Closes a snapshot of the source input, recording failures without replacing
		 * the producer's streaming error. Never holds this Thread's monitor during I/O. */
		void closeInput(){
			final InputStream in=current;
			if(in!=null){
				try{
					in.close();
				}catch(Throwable problem){recordCleanup(problem);}
			}
		}
		/** Parent and producer may report concurrently; the normal path allocates no list. */
		private synchronized void recordCleanup(Throwable problem){
			assert(problem!=null) : "Only caught source-close/producer-cleanup failures enter the shared error list";
			if(cleanupFailures==null){cleanupFailures=new ArrayList<Throwable>();}
			cleanupFailures.add(problem);
		}
		/** Called after the parent's bounded join; draining is synchronized even on failed joining. */
		private synchronized Throwable appendFailures(Throwable first){
			first=appendFailure(first, failure);
			if(cleanupFailures!=null){
				for(Throwable next:cleanupFailures){first=appendFailure(first, next);}
				cleanupFailures=null;
			}
			return first;
		}
		/** Borrowed ownership controller for this invocation only. */
		final Owned child;
		/** Private ordered copy of caller-owned input path names. */
		final String[] files;
		/** Open input visible to cancellation cleanup, or null between successful shards. */
		volatile InputStream current;
		/** Streaming throwable observed by the parent after producer completion. */
		volatile Throwable failure;
		/** Cleanup-only failures; all access uses this object's monitor, never overwrites failure. */
		private ArrayList<Throwable> cleanupFailures;
		/** Decompressed bytes accepted by successful buffered writes, not consumption. */
		volatile long bytes;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Absolute installation, child launcher and execution directory captured at construction. */
	private final File root, launcher, workingDirectory;
	/** Rejects a false condition with IllegalArgumentException, including under -da. */
	private static void require(boolean ok, String why){
		if(!ok){
			throw new IllegalArgumentException(why);
		}
	}
}
