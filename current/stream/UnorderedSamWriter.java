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
 * Unordered single-file SAM writer for multiple producer threads.
 * Producers format complete batches and submit byte blocks to one bounded FIFO
 * output thread. During normal coordinated termination, poison stops acceptance
 * and waits for registered producers before offering its marker, preserving their
 * data-before-terminal order. Batch IDs do not select output order.
 * Construction opens output; the output thread writes the selected header before
 * record batches. Keep header references, input objects and shared SAM settings stable
 * while they are used. Coordinate startup/finalization with producers; delegated I/O
 * and shared-header lookup may wait.
 *
 * @author Collei
 * @date September 20, 2026
 */
public final class UnorderedSamWriter implements Writer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Captures output/header policy, sizes the queue and opens buffered output immediately.
	 * Retains supplied header references without copying; does not start a thread.
	 * Header suppression captures current globals and append-to-existing-file status.
	 * @param ffout_ Nonnull unordered SAM descriptor; append/subprocess settings are respected
	 * @param header_ Optional header used if shared-header selection yields null
	 * @param useSharedHeader_ Whether to request the JVM-wide header before other sources
	 * @param producerThreads Queue-capacity hint, clamped to at least one; not a thread count
	 * @throws IllegalArgumentException If the descriptor is null, ordered or not SAM
	 */
	public UnorderedSamWriter(FileFormat ffout_, ArrayList<byte[]> header_,
			boolean useSharedHeader_, int producerThreads){
		if(ffout_==null || !ffout_.sam() || ffout_.ordered()){
			throw new IllegalArgumentException("UnorderedSamWriter requires unordered SAM output: "+ffout_);
		}
		ffout=ffout_;
		fname=ffout.name();
		header=header_;
		useSharedHeader=useSharedHeader_;
		supressHeader=(ReadStreamWriter.NO_HEADER || (ffout.append() && ffout.exists()));
		supressHeaderSequences=(ReadStreamWriter.NO_HEADER_SEQUENCES || supressHeader);
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
				/** Runs header emission and queued output for this writer. */
				@Override
				public void run(){writeLoop();}
			}, "UnorderedSamWriter-Output");
			outputThread.start();
		}
	}

	/** Wraps and submits reads; the ID does not order this writer's output.
	 * ListNum construction may assign Read.rand under its global deterministic-number option.
	 * @param reads Nonnull source list; conversion skips null entries
	 * @param id Wrapper ID, ignored for output ordering
	 */
	@Override
	public void add(ArrayList<Read> reads, long id){addReads(new ListNum<Read>(reads, id));}

	/** Converts and formats a batch on the submitting thread, then queues copied bytes.
	 * SamWriter.toSamLines includes list entries, attached mates and configured secondary
	 * alignments. It may reuse attached SamLines and change a mate's query name. Keep input
	 * references and shared conversion settings stable during the call.
	 * @param reads Wrapper with a nonnull list, or null to do nothing even before start
	 * @throws IllegalStateException If a nonnull batch arrives before start or after acceptance ends
	 */
	@Override
	public void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		beginAdd();
		try{
			put(format(SamWriter.toSamLines(reads.list)));
		}catch(Throwable t){
			fail(t, false);
			throw asRuntime(t);
		}finally{
			endAdd();
		}
	}

	/** Formats supplied SAM records directly on the submitting thread and queues copied bytes.
	 * Does not add mates or secondary alignments beyond the supplied list. Keep its records
	 * stable during formatting; null entries are skipped and the wrapper ID is ignored.
	 * @param lines Wrapper with a nonnull list, or null to do nothing even before start
	 * @throws IllegalStateException If a nonnull batch arrives before start or after acceptance ends
	 */
	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		beginAdd();
		try{
			put(format(lines.list));
		}catch(Throwable t){
			fail(t, false);
			throw asRuntime(t);
		}finally{
			endAdd();
		}
	}

	/** Serializes nonnull SAM records into this caller's reusable builder, then copies its bytes.
	 * Counts every serialized record and sums SamLine.length, including its CIGAR fallback.
	 * Clears builder length after a successful copy, retaining capacity for the next batch.
	 * @param lines Nonnull list of stable records
	 * @return Owned byte snapshot with record/base counts
	 */
	private OutputJob format(ArrayList<SamLine> lines){
		final ByteBuilder bb=builders.get();
		assert(bb.length()==0);
		long bases=0;
		int reads=0;
		for(SamLine sl : lines){
			if(sl==null){continue;}
			sl.toBytes(bb);
			bb.nl();
			reads++;
			bases+=sl.length();
		}
		final byte[] bytes=bb.toBytes();
		bb.clear();
		return new OutputJob(bytes, reads, bases, false, false);
	}

	/** Offers a job with periodic error checks while the bounded queue is full.
	 * Interruption restores interruption and throws; the offer interval is not a total timeout.
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

	/** Checks startup/acceptance and registers a producer under the state lock. */
	private void beginAdd(){
		synchronized(stateLock){
			if(!started){throw new IllegalStateException("Writer has not started: "+fname);}
			if(!accepting || errorState){throw new IllegalStateException("Writer is closed: "+fname);}
			activeAdds++;
		}
	}

	/** Releases a registered producer and notifies waiters when none remain active. */
	private void endAdd(){
		synchronized(stateLock){
			activeAdds--;
			assert(activeAdds>=0);
			if(activeAdds==0){stateLock.notifyAll();}
		}
	}

	/** Stops accepting batches and waits for registered producers unless an error is observed.
	 * Records the poison decision and queues a marker when no error is recorded.
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

	/** Waits for the finished flag; interruption records failure and ends this wait.
	 * Does not start the output thread or submit a terminal marker.
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

	/** Requests normal termination and waits according to waitForFinish's contract.
	 * @return Current error flag
	 */
	@Override
	public boolean poisonAndWait(){poison(); return waitForFinish();}

	/** Records external failure, abandons queued jobs and requests output-thread interruption.
	 * Does not wait for finalization.
	 */
	@Override
	public void finishError(){fail(null, true);}

	/** Disables acceptance, records failure, notifies waiters and offers abort after clearing the queue.
	 * The abort offer is unchecked; this method alone does not establish completion.
	 * @param t Optional diagnostic printed only when called by the output thread
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

	/** Emits the header, then consumes byte snapshots until a terminal marker or throwable.
	 * Counters advance after each batch write returns; header output is excluded.
	 * The finally block delegates stream/process cleanup, combines its error result and
	 * marks finished when cleanup returns. Counters do not certify durable output.
	 */
	private void writeLoop(){
		boolean normal=false;
		try{
			writeHeader();
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

	/** Selects and writes a header unless the captured all-header suppression flag is set.
	 * Shared-header lookup has priority and may wait for a publisher; otherwise uses a
	 * supplied header, then generates one using current chromosome bounds. Sequence-header
	 * suppression affects generation only, not an existing shared or supplied list.
	 * Appends a newline per header line and writes accumulated bytes at a 16384-byte threshold.
	 * @throws Exception If delegated header selection, generation or output fails
	 */
	private void writeHeader() throws Exception{
		if(supressHeader){return;}
		ArrayList<byte[]> lines=null;
		if(useSharedHeader){lines=SamReadInputStream.getSharedHeader(true);}
		if(lines==null && header!=null){lines=header;}
		if(lines==null){
			lines=SamHeader.makeHeaderList(supressHeaderSequences,
					ReadStreamWriter.MINCHROM, ReadStreamWriter.MAXCHROM);
		}
		if(lines==null){lines=new ArrayList<byte[]>();}
		final ByteBuilder bb=new ByteBuilder();
		for(byte[] line : lines){
			bb.append(line).nl();
			if(bb.length()>=16384){outstream.write(bb.toBytes()); bb.clear();}
		}
		if(bb.length()>0){outstream.write(bb.toBytes());}
	}

	/** Returns the observed number of SAM records accepted by completed batch write calls. */
	@Override
	public long readsWritten(){return readsWritten;}
	/** Returns the observed sum of SamLine lengths accepted by completed batch write calls. */
	@Override
	public long basesWritten(){return basesWritten;}
	/** Returns the retained output name. */
	@Override
	public String fname(){return fname;}
	/** Returns current error state without waiting for finalization. */
	@Override
	public boolean errorState(){return errorState;}
	/** Reports current finished/error flags without requesting output completion. */
	@Override
	public boolean finishedSuccessfully(){return finished && !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Nested Types         ----------------*/
	/*--------------------------------------------------------------*/

	/** Formatted byte snapshot and counts, or a shared terminal marker. */
	private static final class OutputJob{
		/** Retains the supplied payload/counts and marker flags without copying bytes. */
		OutputJob(byte[] bytes_, int reads_, long bases_, boolean poison_, boolean abort_){
			bytes=bytes_; reads=reads_; bases=bases_; poison=poison_; abort=abort_;
		}
		/** Owned batch bytes, null for markers. */
		final byte[] bytes;
		/** Number of serialized SAM records. */
		final int reads;
		/** Sum of serialized SamLine lengths. */
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

	/** Reusable formatting buffer per producer; queued jobs contain byte copies. */
	private final ThreadLocal<ByteBuilder> builders=new ThreadLocal<ByteBuilder>(){
		/** Creates one submitting thread's initial empty buffer. */
		@Override
		protected ByteBuilder initialValue(){return new ByteBuilder();}
	};
	/** Guards lifecycle coordination and active-producer accounting. */
	private final Object stateLock=new Object();
	/** Bounded FIFO of byte snapshots and terminal markers. */
	private final ArrayBlockingQueue<OutputJob> queue;
	/** Retained output descriptor, also consulted during cleanup. */
	private final FileFormat ffout;
	/** Output name used for diagnostics and cleanup. */
	private final String fname;
	/** Optional borrowed header list and byte arrays; no defensive copy at construction. */
	private final ArrayList<byte[]> header;
	/** Whether shared-header selection precedes the provided/generated header. */
	private final boolean useSharedHeader;
	/** Captured all-header suppression, including append to an existing output. */
	private final boolean supressHeader;
	/** Captured sequence-header suppression applied only to generated headers. */
	private final boolean supressHeaderSequences;
	/** Stream opened at construction and used by the output thread. */
	private final OutputStream outstream;
	/** Published error flag, never cleared by this implementation. */
	private volatile boolean errorState=false;
	/** Published completion flag set after delegated cleanup returns. */
	private volatile boolean finished=false;
	/** Whether beginAdd accepts submissions, guarded by stateLock. */
	private boolean accepting=true;
	/** Whether start has been called, guarded by stateLock. */
	private boolean started=false;
	/** Whether poison has recorded its terminal decision, guarded by stateLock. */
	private boolean poisoned=false;
	/** Registered submissions, guarded by stateLock. */
	private int activeAdds=0;
	/** Sole output consumer created by start. */
	private Thread outputThread;
	/** Published SAM-record count after successful batch writes. */
	private volatile long readsWritten=0;
	/** Published SamLine-length sum after successful batch writes. */
	private volatile long basesWritten=0;
}
