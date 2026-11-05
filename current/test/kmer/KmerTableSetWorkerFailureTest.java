package test.kmer;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import kmer.AbstractKmerTable;
import kmer.AbstractKmerTableSet;
import kmer.HashArray1D;
import kmer.HashArray2D;
import kmer.KmerTableSet;
import shared.Shared;
import stream.FastaStreamer;
import stream.StreamerFactory;

/** Checks propagation from actual TableSet workers, including their final buffer flush.
 * One real AAA FASTA record and one explicitly bounded shard avoid automatic table
 * allocation. A rejecting backend supplies the failure; no production hook is used.
 * Run one case per fresh JVM under either assertion mode. Daemon containment bounds
 * a failing baseline without treating an uncaught child error as caller propagation.
 * @author Sucrose, Shinobu
 */
public final class KmerTableSetWorkerFailureTest{

	/** Runs a flush failure, a pre-interrupted caller, or the normal count control. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected case=final-flush|ingestion-flush|interrupted|control.");
		}
		final String mode=args[0].substring(5);
		final boolean ingestion=mode.equals("ingestion-flush"), control=mode.equals("control");
		final boolean interrupted=mode.equals("interrupted");
		if(!ingestion && !control && !interrupted && !mode.equals("final-flush")){
			throw new IllegalArgumentException("Unknown TableSet worker case: "+mode);
		}
		Shared.setThreads(1);
		Shared.AMINO_IN=false;
		AbstractKmerTableSet.MASK_CORE=false;
		AbstractKmerTableSet.MASK_MIDDLE=false;
		AbstractKmerTableSet.DISPLAY_PROGRESS=false;
		ByteFile.FORCE_MODE_BF1=true;
		ByteFile.FORCE_MODE_BF2=false;
		FastaStreamer.DEFAULT_THREADS=0;
		StreamerFactory.FASTA_STREAMER_2=false;
		final int buflen=ingestion ? 1 : 4;
		final KmerTableSet tables=new KmerTableSet(new String[] {
			"k=3", "ways=1", "initialsize=17", "threads=1", "buflen="+buflen,
			"prealloc=f", "prefilter=0", "onepass=f", "rcomp=f", "minprob=0",
			"minavgquality=0", "qtrim=f", "ecco=f", "merge=f"}, 12);
		check(Shared.threads()==1 && tables.k==3 && tables.ways==1 && tables.buflen==buflen
			&& !tables.amino && !tables.rcomp() && !tables.prefilter && !tables.onePass
			&& tables.tables()==null && tables.prefilterArray==null, mode,
			"Fixture requires one DNA worker and no automatically allocated tables or prefilter");
		final AbstractKmerTable backend=(control || interrupted) ? new HashArray1D(new int[] {17}, -1L) :
			new HashArray2D(new int[] {17}, -1L);
		final AbstractKmerTable[] shards={backend};
		install(tables, shards);
		check(tables.tables()==shards && backend.arrayLength()==77 && backend.size()==0,
			mode, "Fixture must use one empty, explicitly bounded shard");
		final Path input=Files.createTempFile("tableset-worker-failure-", ".fa");
		input.toFile().deleteOnExit();
		Files.write(input, ">r\nAAA\n".getBytes(StandardCharsets.US_ASCII));
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final AtomicBoolean interruptRestored=new AtomicBoolean();
		final AtomicLong added=new AtomicLong(-1);
		final FailureGroup group=new FailureGroup();
		final Thread coordinator=new Thread(group, new Runnable(){
			@Override
			public void run(){
				if(interrupted){Thread.currentThread().interrupt();}
				try{added.set(tables.loadKmers(input.toString(), null));}
				catch(Throwable t){failure.set(t);}
				finally{interruptRestored.set(Thread.currentThread().isInterrupted());}
			}
		}, "KmerTableSetWorkerFailureTest-coordinator");
		coordinator.setDaemon(true);
		coordinator.start();
		coordinator.join(10000);
		check(!coordinator.isAlive(), mode, "TableSet loader did not finish within the bounded test");
		if(control){
			check(failure.get()==null && group.failure.get()==null && !tables.errorState,
				mode, "Normal TableSet loading reported a failure");
			check(added.get()==1 && backend.size()==1 && backend.getValue(0)==1,
				mode, "AAA must create exactly one raw key0 with count1");
			check(tables.readsIn==1 && tables.basesIn==3 && tables.kmersIn==1,
				mode, "Normal loading must report one read, three bases and one window");
		}else if(interrupted){
			check(hasInterruption(failure.get()), mode, "Caller interruption did not reach the TableSet caller");
			check(interruptRestored.get(), mode, "Caller interrupt status was not restored");
			check(tables.errorState, mode, "TableSet did not retain its interruption error state");
		}else{
			final Throwable caller=unsupported(failure.get());
			final Throwable observed=caller!=null ? caller : unsupported(group.failure.get());
			check(observed!=null, mode, "The backend did not produce the intended Unsupported failure");
			check(hasFrame(observed, ingestion ? "incrementAndReturnNumCreated" : "flush"),
				mode, "The intended HashBuffer ingestion/final-flush path was not reached");
			check(caller!=null, mode, "Worker failure did not reach the TableSet caller");
			check(tables.errorState, mode, "TableSet did not retain its worker error state");
		}
		final Thread[] threads=new Thread[16];
		final int listed=group.enumerate(threads, true);
		check(listed<threads.length, mode, "Unexpected number of fixture threads");
		for(int i=0; i<listed; i++){
			check(threads[i]==null || !threads[i].isAlive(), mode, "TableSet left a fixture worker alive");
		}
		Files.deleteIfExists(input);
		System.out.println("PASS tableset-worker case="+mode);
	}

	/** Injects only explicit storage; actual loadKmers, input and worker code remain in use. */
	private static void install(final KmerTableSet tables, final AbstractKmerTable[] shards) throws ReflectiveOperationException{
		final Field field=KmerTableSet.class.getDeclaredField("tables");
		field.setAccessible(true);
		field.set(tables, shards);
	}

	/** Finds the real backend rejection even when the caller wraps it with its original cause. */
	private static Throwable unsupported(Throwable t){
		for(; t!=null; t=t.getCause()){
			if(t instanceof RuntimeException && "Unsupported.".equals(t.getMessage())){return t;}
		}
		return null;
	}

	/** Recognizes the original interruption even when the coordinator wraps checked failures. */
	private static boolean hasInterruption(Throwable t){
		for(; t!=null; t=t.getCause()){if(t instanceof InterruptedException){return true;}}
		return false;
	}

	/** Verifies the failure phase without depending on line numbers or private worker method names. */
	private static boolean hasFrame(final Throwable failure, final String method){
		for(StackTraceElement frame : failure.getStackTrace()){
			if(frame.getClassName().equals("kmer.HashBuffer") && frame.getMethodName().equals(method)){return true;}
		}
		return false;
	}

	/** Records uncaught worker failures separately from what loadKmers reports to its caller. */
	private static final class FailureGroup extends ThreadGroup{
		/** Workers inherit this group and the coordinator's daemon status. */
		FailureGroup(){super("KmerTableSetWorkerFailureTest-workers");}
		@Override
		public void uncaughtException(Thread thread, Throwable error){failure.compareAndSet(null, error);}
		/** First uncaught child failure; does not itself establish caller propagation. */
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	}

	/** Keeps fixture and contract checks active with assertions disabled. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR011 case="+mode+": "+message);}
	}
}
