package test.kmer;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.HashArray2D;
import kmer.TableLoaderLockFree;
import stream.FastaStreamer;
import stream.FastaStreamerZT;
import stream.StreamerFactory;
import structures.ByteBuilder;

/**
 * Checks caller-visible loader failures and worker cleanup with tiny literal/file inputs.
 * A real unsupported backend operation triggers a worker error; a null literal
 * triggers a producer error after workers start. A daemon coordinator bounds the
 * failing baseline: its inherited daemon workers cannot keep the test JVM alive.
 * File cases use 80 tiny FASTA records and one-record reader batches to exceed
 * the loader queue capacity. Each case runs in a fresh JVM with bounded memory.
 * @author Shinobu
 */
public final class TableLoaderFailureTest{

	/** Runs one error path or a normal single-key control under either assertion mode. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected case=worker-failure|producer-failure|full-queue|interrupted|control|file-control.");
		}
		final String mode=args[0].substring(5);
		final boolean full=mode.equals("full-queue"), interrupt=mode.equals("interrupted");
		final boolean file=full || mode.equals("file-control");
		final boolean worker=mode.equals("worker-failure") || full, producer=mode.equals("producer-failure");
		if(!worker && !producer && !interrupt && !file && !mode.equals("control")){
			throw new IllegalArgumentException("Unknown loader case: "+mode);
		}
		final AbstractKmerTable[] tables=new AbstractKmerTable[7];
		for(int i=0; i<tables.length; i++){
			tables[i]=(worker && i==1 ? new HashArray2D(new int[] {17}, -1L) :
				new HashArray1D(new int[] {17}, -1L));
		}
		final TableLoaderLockFree loader=new TableLoaderLockFree(tables, 3, 0, 0, 0, 0, false, false);
		loader.storeMode(TableLoaderLockFree.INCREMENT);
		TableLoaderLockFree.DISPLAY_PROGRESS=false;
		TableLoaderLockFree.REPLICATE_AMBIGUOUS=false;
		final Path input;
		if(file){
			ByteFile.FORCE_MODE_BF1=true;
			FastaStreamer.DEFAULT_THREADS=0;
			StreamerFactory.FASTA_STREAMER_2=false;
			FastaStreamerZT.TARGET_LIST_SIZE=1;
			input=Files.createTempFile("table-loader-failure-", ".fa");
			input.toFile().deleteOnExit();
			final ByteBuilder fasta=new ByteBuilder();
			for(int i=0; i<80; i++){fasta.append(">r").append(i).append("\nAAA\n");}
			Files.write(input, fasta.toBytes());
		}else{input=null;}
		final String[] ref=file ? new String[] {input.toString()} : null;
		final String[] literal=file ? null : producer ? new String[] {null} : new String[] {"AAA"};
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final AtomicBoolean interruptRestored=new AtomicBoolean();
		final AtomicLong added=new AtomicLong(-1);
		final FailureGroup group=new FailureGroup();
		final Thread coordinator=new Thread(group, new Runnable(){
			@Override
			public void run(){
				if(interrupt){Thread.currentThread().interrupt();}
				try{added.set(loader.processData(ref, literal, false, false, false));}
				catch(Throwable t){failure.set(t);}
				finally{interruptRestored.set(Thread.currentThread().isInterrupted());}
			}
		}, "TableLoaderFailureTest-coordinator");
		coordinator.setDaemon(true);
		coordinator.start();
		coordinator.join(10000);
		if(input!=null){Files.deleteIfExists(input);}
		if(full){
			check(unsupported(failure.get()) || unsupported(group.failure.get()), mode,
				"The rejecting backend did not produce the intended Unsupported failure");
		}
		check(!coordinator.isAlive(), mode, "Loader call did not finish within the bounded test");

		if(worker){
			check(unsupported(failure.get()) || unsupported(group.failure.get()), mode,
				"The rejecting backend did not produce the intended Unsupported failure");
			check(unsupported(failure.get()), mode, "Worker failure did not reach the caller");
		}else if(producer){
			check(hasCause(failure.get(), NullPointerException.class), mode,
				"The null literal did not produce the intended producer failure");
		}else if(interrupt){
			check(hasCause(failure.get(), InterruptedException.class), mode,
				"The caller interrupt did not reach the caller as a failure");
			check(interruptRestored.get(), mode, "Caller interrupt status was not restored");
		}else{
			check(failure.get()==null && group.failure.get()==null, mode, "Normal loading reported a failure");
			// Both inputs create one key; the file has eighty occurrences of that same key.
			check(added.get()==1 && tables[1].getValue(64)==(file ? 80 : 1) && !loader.errorState,
				mode, "Normal loading must preserve the expected AAA count and return without an error");
		}
		final Thread[] threads=new Thread[32];
		final int listed=group.enumerate(threads, true);
		check(listed<threads.length, mode, "Unexpected number of fixture threads");
		int live=0;
		for(int i=0; i<listed; i++){if(threads[i]!=null && threads[i].isAlive()){live++;}}
		check(live==0, mode, "Loader left "+live+" worker threads alive");
		if(worker || producer || interrupt){check(loader.errorState, mode, "Loader did not retain its error state");}
		System.out.println("PASS loader case="+mode);
	}

	/** Recognizes the real multi-value backend's unsupported scalar-count operation. */
	private static boolean unsupported(Throwable t){
		for(; t!=null; t=t.getCause()){
			if(t instanceof RuntimeException && "Unsupported.".equals(t.getMessage())){return true;}
		}
		return false;
	}

	/** Allows the coordinator to wrap a failure while preserving its original cause. */
	private static boolean hasCause(Throwable t, final Class<? extends Throwable> type){
		for(; t!=null; t=t.getCause()){if(type.isInstance(t)){return true;}}
		return false;
	}

	/** Observes uncaught child errors without mistaking them for caller propagation. */
	private static final class FailureGroup extends ThreadGroup{
		FailureGroup(){super("TableLoaderFailureTest-workers");}
		@Override
		public void uncaughtException(Thread thread, Throwable error){failure.compareAndSet(null, error);}
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	}

	/** Remains active under both -ea and -da. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR116 case="+mode+": "+message);}
	}
}
