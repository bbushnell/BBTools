package stream;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Unordered single/interleaved FASTQ writer for multiple producer threads.
 * Each producer formats a complete batch before queue submission, retaining
 * mate adjacency while avoiding ordered-list head-of-line blocking.
 * Numeric batch IDs do not determine output order. The output thread writes queued
 * byte snapshots; keep submitted reads and shared FASTQ settings stable during formatting.
 * Construction opens output; normal use starts the writer, submits batches, then calls
 * poisonAndWait and checks its error result. Coordinate startup/finalization with
 * producer use; delegated output and cleanup may block.
 *
 * @author Collei
 * @date September 20, 2026
 */
public final class UnorderedFastqWriter implements Writer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains the output descriptor, sizes the queue and opens buffered output immediately.
	 * Preserves the descriptor's append and subprocess settings; does not start a thread.
	 * @param ffout_ Nonnull descriptor for unordered FASTQ output
	 * @param producerThreads Queue-capacity hint, clamped to at least one; not an output-thread count
	 * @throws IllegalArgumentException If the descriptor is null, ordered or not FASTQ
	 */
	public UnorderedFastqWriter(FileFormat ffout_, int producerThreads){
		if(ffout_==null || !ffout_.fastq() || ffout_.ordered()){
			throw new IllegalArgumentException("UnorderedFastqWriter requires unordered FASTQ output: "+ffout_);
		}
		ffout=ffout_;
		fname=ffout.name();
		final int capacity=Tools.max(4, (3*Tools.max(1, producerThreads))/2+4);
		queue=new ArrayBlockingQueue<OutputJob>(capacity);
		outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates and starts the sole output thread while holding the state lock.
	 * @throws IllegalStateException If start has already been called
	 */
	@Override
	public void start(){
		synchronized(stateLock){
			if(started){throw new IllegalStateException("Writer already started: "+fname);}
			started=true;
			outputThread=new Thread(new Runnable(){
				/** Runs this writer's queue-consumption method. */
				@Override
				public void run(){writeLoop();}
			}, "UnorderedFastqWriter-Output");
			outputThread.start();
		}
	}

	/** Wraps and submits a read list; the ID does not order this writer's output.
	 * ListNum construction may assign Read.rand under its global deterministic-number option.
	 * @param reads Nonnull list to format synchronously; null entries are skipped
	 * @param id Wrapper ID, ignored for output ordering
	 */
	@Override
	public void add(ArrayList<Read> reads, long id){addReads(new ListNum<Read>(reads, id));}

	/** Formats a nonnull batch on the calling thread and submits its copied bytes.
	 * Reads, mates and serialization settings must remain stable during formatting.
	 * Null wrappers are ignored even before start; nonnull wrappers need a nonnull list.
	 * Formatting/submission failures are recorded and rethrown through asRuntime.
	 * @param reads Batch whose ID is ignored; may be null
	 * @throws IllegalStateException If a nonnull batch arrives before start or after acceptance ends
	 */
	@Override
	public void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		beginAdd();
		try{
			put(format(reads.list));
		}catch(Throwable t){
			fail(t, false);
			throw asRuntime(t);
		}finally{
			endAdd();
		}
	}

	/** Converts nonnull SAM entries to independent Read wrappers and submits them as FASTQ.
	 * Uses sequence, quality and query name; does not transfer SAM flags or mate links.
	 * Null wrappers and entries are skipped. Payload arrays are borrowed during formatting.
	 * @param lines SAM records with a nonnull list, or null to do nothing
	 */
	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		final ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines){
			if(sl!=null){reads.add(new Read(sl.seq, sl.qual, sl.qname, -1, false));}
		}
		addReads(new ListNum<Read>(reads, lines.id));
	}

	/** Serializes a batch using this caller's reusable builder, then copies its bytes.
	 * Pair-number-zero entries emit themselves and an attached mate; pair-number-one
	 * entries emit only themselves. Null entries are skipped. A successful copy is followed
	 * by clearing the builder's length, retaining its capacity for this thread's next batch.
	 * @param reads Nonnull source list kept stable by the caller
	 * @return Owned byte snapshot with individual-read and base counts
	 */
	private OutputJob format(ArrayList<Read> reads){
		final ByteBuilder bb=builders.get();
		assert(bb.length()==0) : "Thread-local FASTQ builder was not cleared after its prior batch";
		long bases=0;
		int count=0;
		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(r1!=null){
				r1.toFastq(bb).nl();
				count++;
				bases+=r1.length();
			}
			if(r2!=null){
				r2.toFastq(bb).nl();
				count++;
				bases+=r2.length();
			}
		}
		final byte[] bytes=bb.toBytes();
		bb.clear();
		return new OutputJob(bytes, count, bases, false, false);
	}

	/** Offers a job with periodic error-state checks while the bounded queue is full.
	 * Interrupted submission restores interruption and throws; this is not a timeout limit.
	 * @param job Nonnull formatted batch or terminal marker
	 */
	private void put(OutputJob job){
		try{
			while(!queue.offer(job, 100, TimeUnit.MILLISECONDS)){
				if(errorState){throw new RuntimeException("Writer failed while submitting output: "+fname);}
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new RuntimeException("Interrupted while submitting output: "+fname, e);
		}
		if(errorState){throw new RuntimeException("Writer failed while submitting output: "+fname);}
	}

	/** Checks start/acceptance state and registers a producer under the state lock. */
	private void beginAdd(){
		synchronized(stateLock){
			if(!started){throw new IllegalStateException("Writer has not started: "+fname);}
			if(!accepting || errorState){throw new IllegalStateException("Writer is closed: "+fname);}
			activeAdds++;
		}
	}

	/** Releases a registered producer and notifies waiters when no submissions remain active. */
	private void endAdd(){
		synchronized(stateLock){
			activeAdds--;
			assert(activeAdds>=0) : "FASTQ producer count became negative; beginAdd/endAdd are unbalanced";
			if(activeAdds==0){stateLock.notifyAll();}
		}
	}

	/** Stops accepting batches and waits for registered producers unless an error is observed.
	 * Records the poison decision and queues its marker when no error is recorded.
	 * Interruption while waiting for producers restores interruption, records failure and returns.
	 * Marker submission follows put's throwing behavior; does not wait for final output.
	 */
	@Override
	public void poison(){
		synchronized(stateLock){
			if(poisoned){return;}
			accepting=false;
			while(activeAdds>0 && !errorState){
				try{stateLock.wait();}catch(InterruptedException e){
					Thread.currentThread().interrupt();
					fail(e, false);
					return;
				}
			}
			poisoned=true;
		}
		if(!errorState){put(OutputJob.POISON);}
	}

	/** Waits for the finished flag; interruption records an error and ends this wait.
	 * Does not submit a terminal marker or start the output thread.
	 * @return Current error flag, including any error recorded after interruption
	 */
	@Override
	public boolean waitForFinish(){
		synchronized(stateLock){
			while(!finished){
				try{stateLock.wait();}catch(InterruptedException e){
					Thread.currentThread().interrupt();
					fail(e, false);
					break;
				}
			}
		}
		return errorState;
	}

	/** Requests normal termination, then waits according to waitForFinish's contract.
	 * @return Current error flag
	 */
	@Override
	public boolean poisonAndWait(){poison(); return waitForFinish();}

	/** Records an external failure, abandons queued jobs and requests output-thread interruption.
	 * Does not wait for output finalization.
	 */
	@Override
	public void finishError(){fail(null, true);}

	/** Records failure, disables acceptance, notifies waiters and replaces queued work with an abort offer.
	 * The offer is not checked; this method alone does not establish completion.
	 * @param t Optional diagnostic, printed only when called by the output thread
	 * @param interruptWriter Whether to interrupt a different output thread when present
	 */
	private void fail(Throwable t, boolean interruptWriter){
		final Thread thread;
		synchronized(stateLock){
			errorState=true;
			accepting=false;
			thread=outputThread;
			stateLock.notifyAll();
		}
		queue.clear();
		queue.offer(OutputJob.ABORT);
		if(interruptWriter && thread!=null && thread!=Thread.currentThread()){thread.interrupt();}
		if(t!=null && thread==Thread.currentThread()){t.printStackTrace();}
	}

	/** Retains RuntimeExceptions and wraps other throwables for the submitting caller.
	 * @param t Failure to expose through the submission API
	 * @return Existing or newly wrapped runtime exception
	 */
	private RuntimeException asRuntime(Throwable t){
		return t instanceof RuntimeException ? (RuntimeException)t : new RuntimeException(t);
	}

	/** Consumes queued snapshots until a terminal marker or throwable is encountered.
	 * Advances counters after OutputStream.write returns for a batch, before finalization.
	 * The finally block delegates stream/process cleanup, combines its error result and
	 * marks finished when that cleanup returns. Counter values do not certify durable output.
	 */
	private void writeLoop(){
		boolean normal=false;
		try{
			while(true){
				final OutputJob job=queue.take();
				if(job.abort){break;}
				if(job.poison){normal=true; break;}
				outstream.write(job.bytes);
				readsWritten+=job.reads;
				basesWritten+=job.bases;
			}
		}catch(Throwable t){
			if(!(t instanceof InterruptedException && errorState)){t.printStackTrace();}
			fail(t, false);
		}finally{
			final boolean closeError=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
			synchronized(stateLock){
				errorState|=closeError || !normal;
				finished=true;
				stateLock.notifyAll();
			}
		}
	}

	/** Returns the observed individual-read count accepted by completed batch write calls. */
	@Override
	public long readsWritten(){return readsWritten;}
	/** Returns the observed base count accepted by completed batch write calls. */
	@Override
	public long basesWritten(){return basesWritten;}
	/** Returns the retained output name. */
	@Override
	public String fname(){return fname;}
	/** Returns the current error flag without waiting for completion. */
	@Override
	public boolean errorState(){return errorState;}
	/** Reports the current finished/error flags; does not wait or finalize output. */
	@Override
	public boolean finishedSuccessfully(){return finished && !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Nested Types         ----------------*/
	/*--------------------------------------------------------------*/

	/** Formatted byte snapshot and counts, or a shared terminal marker. */
	private static final class OutputJob{
		/** Retains supplied payload/counts and marker flags without copying the byte array. */
		OutputJob(byte[] bytes_, int reads_, long bases_, boolean poison_, boolean abort_){
			bytes=bytes_; reads=reads_; bases=bases_; poison=poison_; abort=abort_;
		}
		/** Owned batch bytes; null for markers. */
		final byte[] bytes;
		/** Individual records represented in bytes. */
		final int reads;
		/** Sum of serialized read lengths. */
		final long bases;
		/** Marks normal termination. */
		final boolean poison;
		/** Marks abandoned output. */
		final boolean abort;
		/** Shared normal terminal with no payload. */
		static final OutputJob POISON=new OutputJob(null, 0, 0, true, false);
		/** Shared abort terminal with no payload. */
		static final OutputJob ABORT=new OutputJob(null, 0, 0, false, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reusable formatting buffer per submitting thread; queued jobs hold copies. */
	private final ThreadLocal<ByteBuilder> builders=new ThreadLocal<ByteBuilder>(){
		/** Creates the initial empty buffer for one submitting thread. */
		@Override
		protected ByteBuilder initialValue(){return new ByteBuilder();}
	};
	/** Guards lifecycle coordination and the active-producer count. */
	private final Object stateLock=new Object();
	/** Bounded FIFO of formatted byte snapshots and terminal markers. */
	private final ArrayBlockingQueue<OutputJob> queue;
	/** Retained output descriptor, also consulted during finalization. */
	private final FileFormat ffout;
	/** Retained output name for diagnostics and delegated cleanup. */
	private final String fname;
	/** Stream opened during construction and used by the output thread. */
	private final OutputStream outstream;
	/** Published error state; never cleared by this implementation. */
	private volatile boolean errorState=false;
	/** Published completion flag set after delegated cleanup returns. */
	private volatile boolean finished=false;
	/** Whether beginAdd may accept more submissions, guarded by stateLock. */
	private boolean accepting=true;
	/** Whether start has been called, guarded by stateLock. */
	private boolean started=false;
	/** Whether poison has recorded its terminal decision, guarded by stateLock. */
	private boolean poisoned=false;
	/** Number of registered submissions, guarded by stateLock. */
	private int activeAdds=0;
	/** Sole queue consumer created by start. */
	private Thread outputThread;
	/** Published count of individual reads after successful batch write calls. */
	private volatile long readsWritten=0;
	/** Published base count after successful batch write calls. */
	private volatile long basesWritten=0;
}
