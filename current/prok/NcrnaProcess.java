package prok;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import dna.Data;

/** Owns one ncRNA caller or trainer in a separate Unix process group.
 * The small launcher reports its group in a temporary file before accepting permission
 * to start work. Cancellation signals only that owned group, then joins the child
 * and its stream workers. Public tool launchers need no changes. Requires setsid,
 * bash and kill; Java APIs remain Java 8 compatible. Stdout/stderr are inherited,
 * so cancelling a child blocked on output never depends on joining a Java relay.
 * The temporary identity is deleted before work starts; results come from files.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaProcess{

	/*--------------------------------------------------------------*/
	/*----------------           Execution          ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs one command until child and optional input feeder finish.
	 * @param command ProcessBuilder arguments, unchanged until this method returns
	 * @param input Optional writer of child stdin; null closes unused stdin
	 * @param description Diagnostic name of the operation
	 */
	public static void run(final ArrayList<String> command, final InputWriter input, final String description){
		if(command==null || command.isEmpty() || description==null){
			throw new IllegalArgumentException("Child execution requires a command and diagnostic description");
		}
		if(Thread.currentThread().isInterrupted()){
			throw new IllegalStateException(description+" interrupted before launch");
		}
		try(Session session=new Session(command, input)){
			final int exit=session.child.waitFor();
			// A waiting shell can exit while a descendant still holds a pipe open.
			// Stop the owned group before joining its input worker.
			session.stop();
			if(session.feeder!=null){session.feeder.join();}
			final Throwable feedFailure=session.feeder==null ? null : session.feeder.failure;
			if(exit!=0){
				throw new IllegalStateException(description+" exited "+exit+
					(feedFailure==null ? "" : "; input stream failed"), feedFailure);
			}
			if(feedFailure!=null){throw new IllegalStateException(description+" input stream failed", feedFailure);}
		}catch(IOException failure){
			throw new IllegalStateException(description+" process I/O failed", failure);
		}catch(InterruptedException failure){
			Thread.currentThread().interrupt();
			throw new IllegalStateException(description+" interrupted", failure);
		}
	}

	/** Appends explicit heap/assertion settings consumed by native script parsers.
	 * Session passes the parent's explicit stack setting through the child environment.
	 */
	public static void addJvmOptions(final ArrayList<String> command, final long heapMiB){
		if(command==null || heapMiB<1){throw new IllegalArgumentException("Child JVM options require a command and positive MiB heap");}
		command.add("-Xmx"+heapMiB+"m");
		command.add("-Xms"+Math.min(heapMiB, 64)+"m");
		command.add("-ea");
		command.add("eoom");
	}

	/** Supplies bytes to borrowed child stdin; the process owner closes the stream. */
	@FunctionalInterface
	public interface InputWriter{
		/** Writes the complete declared input, propagating any read/write failure. */
		void write(OutputStream output) throws IOException;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Child Ownership      ----------------*/
	/*--------------------------------------------------------------*/

	/** One try-with-resources scope preserves a body failure if cleanup also fails. */
	private static final class Session implements AutoCloseable{
		Session(final ArrayList<String> command, final InputWriter input) throws IOException, InterruptedException{
			final String launcher=Data.findPath("?dev/ncrnaprocessgroup.sh", false);
			if(launcher==null){throw new IOException("Cannot locate dev/ncrnaprocessgroup.sh in this BBTools installation");}
			try{
				try(GroupIdentity identity=new GroupIdentity()){
					final ArrayList<String> owned=new ArrayList<String>();
					owned.add("setsid");
					owned.add("--wait");
					owned.add("bash");
					owned.add(launcher);
					owned.add(identity.path.toString());
					owned.addAll(command);
					final ProcessBuilder builder=new ProcessBuilder(owned)
						.redirectOutput(ProcessBuilder.Redirect.INHERIT).redirectError(ProcessBuilder.Redirect.INHERIT);
					String stack=null;
					for(String option:ManagementFactory.getRuntimeMXBean().getInputArguments()){
						if(option.startsWith("-Xss")){stack=option;}
					}
					if(stack!=null){
						final String prior=builder.environment().get("JAVA_TOOL_OPTIONS");
						builder.environment().put("JAVA_TOOL_OPTIONS", prior==null || prior.isEmpty() ? stack : prior+" "+stack);
					}
					child=builder.start();
					shutdown=new Thread(this::stop, "ncrna-child-shutdown");
					Runtime.getRuntime().addShutdownHook(shutdown);
					registered=true;
					final int announced=readGroup(identity.path, child);
					synchronized(this){
						if(stopped){throw new IOException("Child startup cancelled before group adoption");}
						group=announced;
					}
				}// Identity deletion must succeed before the workload is released.
				synchronized(this){
					if(stopped || Thread.currentThread().isInterrupted()){
						throw new IOException("Child startup cancelled before workload release");
					}
					final OutputStream sink=child.getOutputStream();
					sink.write(START);
					sink.flush();
					if(input==null){sink.close();}
					else{
						feeder=new Feeder(this, input);
						feeder.start();
					}
				}
			}catch(IOException | InterruptedException | RuntimeException | Error failure){
				try{
					close();
				}catch(RuntimeException | Error cleanup){
					if(cleanup!=failure){failure.addSuppressed(cleanup);}
				}
				throw failure;
			}
		}

		/** Terminates the owned group, including descendants left by an exited launcher. */
		synchronized void stop(){
			if(stopped){return;}
			stopped=true;
			try{
				if(group>1){
					final Process signal=new ProcessBuilder("kill", "-KILL", "--", "-"+group)
						.redirectError(new File("/dev/null")).start();
					final int exit=awaitExit(signal);
					// kill returns 1 when the completed group's last process has already exited.
					if(exit!=0 && exit!=1){throw new IOException("Owned-group signal exited "+exit);}
				}
			}catch(IOException failure){
				throw new IllegalStateException("Could not terminate the owned ncRNA process group", failure);
			}finally{
				if(child!=null){
					if(child.isAlive()){child.destroyForcibly();}
					awaitExit(child);
				}
			}
		}

		/** Joins owned workers even when group cleanup fails; borrowed stdout stays open. */
		@Override
		public void close(){
			try{stop();}
			finally{
				join(feeder);
				if(registered){
					try{Runtime.getRuntime().removeShutdownHook(shutdown);}
					catch(IllegalStateException shuttingDown){
						// The hook owns the same group while the VM shuts down.
					}
				}
			}
		}

		/** Direct launcher process, which may be waiting for the owned session leader. */
		Process child;
		/** Positive group announced before work starts; zero means no workload was released. */
		private int group;
		/** Serializes startup release against shutdown and makes group cleanup idempotent. */
		private boolean stopped;
		/** Hook that terminates this same group if the parent VM exits. */
		private Thread shutdown;
		/** Whether hook registration succeeded and cleanup must attempt removal. */
		private boolean registered;
		/** Optional input writer, joined before the session returns or closes. */
		Feeder feeder;
	}

	/** Separates pipe backpressure from the interruptible wait that owns cancellation. */
	private static final class Feeder extends Thread{
		Feeder(final Session owner_, final InputWriter input_){
			super("ncrna-child-input");
			owner=owner_;
			input=input_;
		}

		@Override
		public void run(){
			try(OutputStream sink=owner.child.getOutputStream()){
				input.write(sink);
			}catch(Throwable caught){
				failure=caught;
				try{owner.stop();}
				catch(RuntimeException cleanup){caught.addSuppressed(cleanup);}
			}
		}

		/** Borrowed process whose stdin is owned by this feeder's output scope. */
		private final Session owner;
		/** Borrowed callback that supplies the complete declared child input. */
		private final InputWriter input;
		/** Published to the owner by Thread.join, never read while this thread runs. */
		Throwable failure;
	}

	/** Owns the single small startup file; no markers survive workload release. */
	private static final class GroupIdentity implements AutoCloseable{
		GroupIdentity() throws IOException{path=Files.createTempFile("ncrna-group-", ".pid");}
		@Override
		public void close() throws IOException{Files.delete(path);}
		final Path path;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Startup / Cleanup     ----------------*/
	/*--------------------------------------------------------------*/

	/** Startup-only wait for a complete decimal line, capped at thirty seconds.
	 * The temporary file is a tiny IPC record, not a biological data input. A partial
	 * write is never accepted as a different numeric group. No polling occurs after release.
	 */
	private static int readGroup(final Path identity, final Process child) throws IOException, InterruptedException{
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(30);
		while(child.isAlive() && System.nanoTime()<deadline){
			final byte[] line=Files.readAllBytes(identity);
			if(line.length>11){throw new IOException("Invalid ncRNA process-group startup line");}
			if(line.length>0 && line[line.length-1]=='\n'){
				long group=0;
				for(int i=0; i<line.length-1; i++){
					final int b=line[i];
					if(b<'0' || b>'9'){throw new IOException("Invalid ncRNA process-group startup line");}
					group=group*10+b-'0';
				}
				if(group<=1 || group>Integer.MAX_VALUE){throw new IOException("Invalid owned process group: "+group);}
				return (int)group;
			}
			Thread.sleep(10);
		}
		throw new IOException(child.isAlive() ? "Timed out waiting for the owned process group" : "Child exited before reporting its process group");
	}

	/** Waits through repeated interruption so cleanup never abandons an owned process. */
	private static int awaitExit(final Process child){
		boolean interrupted=false;
		try{
			for(;;){
				try{return child.waitFor();}
				catch(InterruptedException again){interrupted=true;}
			}
		}finally{
			if(interrupted){Thread.currentThread().interrupt();}
		}
	}

	private static void join(final Thread thread){
		if(thread==null){return;}
		boolean interrupted=false;
		while(thread.isAlive()){
			try{thread.join();}
			catch(InterruptedException again){interrupted=true;}
		}
		if(interrupted){Thread.currentThread().interrupt();}
	}

	/** Private stdin token consumed by the launcher before it execs the public command. */
	private static final byte[] START={'s', 't', 'a', 'r', 't', '\n'};
	private NcrnaProcess(){}
}
