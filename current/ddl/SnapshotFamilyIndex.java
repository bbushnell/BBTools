package ddl;

import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.LockSupport;

import ddl.MutableFamilyIndex.Cosine;
import ddl.MutableFamilyIndex.FrozenView;
import ddl.MutableFamilyIndex.SketchSource;
import ddl.MutableFamilyIndex.Workspace;
import template.ThreadWaiter;

/** Nonblocking greedy ranking over an immutable indexed prefix and a directly scanned tail.
 * One background thread builds replacement snapshots; workers never modify an index.
 * Source rows must be immutable and IDs must form a monotonically published dense prefix.
 * Each worker owns its cosine callback, query/output arrays and Workspace.
 * The owner starts the index before query workers and closes it after joining them.
 * The source is borrowed and must remain usable until close returns. Replacing a
 * published row is allowed, but arrays retained by older snapshots remain immutable.
 * @author Brian, Fischl
 */
public final class SnapshotFamilyIndex implements AutoCloseable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** A fraction of .02 requests rebuilding once the live tail exceeds 2% of the prefix.
	 * No worker starts here: the owner calls start() before launching query workers.
	 * @param width_ Supported power-of-two register width, validated by snapshot construction
	 * @param rebuildFraction_ Positive finite tail/prefix growth threshold
	 * @param source_ Borrowed source publishing immutable rows and a monotonic visible count
	 */
	public SnapshotFamilyIndex(final int width_, final double rebuildFraction_, final SketchSource source_){
		if(source_==null || !Double.isFinite(rebuildFraction_) || rebuildFraction_<=0){
			throw new IllegalArgumentException("Snapshot ranking requires a source and positive finite growth fraction");
		}
		width=width_;
		rebuildFraction=rebuildFraction_;
		source=source_;
		current=new AtomicReference<FrozenView>(MutableFamilyIndex.snapshot(width, 0, source));
		builder=new BuilderThread();
	}

	/** Starts the builder exactly once and requests an initial snapshot when the source is nonempty.
	 * The owner serializes startup/close; Thread.start publishes configuration to the builder.
	 * A failed Thread.start leaves this index closed and unavailable for queries.
	 */
	public void start(){
		if(started || closing){throw new IllegalStateException("Snapshot builder must be started exactly once before close");}
		final int count=source.size();
		if(count<0){throw new IllegalStateException("Published family count cannot be negative");}
		started=true;
		try{builder.start();}
		catch(RuntimeException | Error problem){
			started=false;
			closing=true;
			throw problem;
		}
		requestIfNeeded(count, 0);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Queries           ----------------*/
	/*--------------------------------------------------------------*/

	/** Captures the prefix first, then the acquire-visible family count, once per query.
	 * New snapshot publication cannot change this query's prefix/tail boundary.
	 * Cosine and subsequent alignment may use newer immutable family versions.
	 * Output and ranking semantics are those of FrozenView.rankWithTail; scratch is
	 * exclusive to this query and must not be shared or reentered by the callback.
	 * @return Number of ranked families written to the caller-owned output arrays
	 */
	public int rank(final char[] query, final int limit, final int[] outIds, final double[] outScores,
			final Cosine cosine, final Workspace work){
		if(!started || closing){throw new IllegalStateException("Queries require a started, open snapshot index");}
		checkFailure();
		final FrozenView view=current.get();
		final int count=source.size();
		requestIfNeeded(count, view.size());
		final int found=view.rankWithTail(query, count, source, limit, outIds, outScores, cosine, work);
		checkFailure();
		return found;
	}

	/** Coalesces rebuild requests; only the false-to-true transition wakes the background thread. */
	private void requestIfNeeded(final int count, final int prefix){
		if(count<prefix){throw new IllegalStateException("Published family count shrank below the snapshot prefix");}
		if(needsRebuild(count, prefix) && requested.compareAndSet(false, true)){LockSupport.unpark(builder);}
	}

	/** Integer counts with a fractional threshold; any nonempty population can replace an empty prefix. */
	private boolean needsRebuild(final int count, final int prefix){
		assert(count>=prefix && prefix>=0) : "Snapshot prefixes can only grow during one greedy pass";
		return count-prefix>prefix*rebuildFraction;
	}

	/** Throws with the retained builder failure as its cause; does not restart or change the builder.
	 * rank() checks before and after scoring, and close() checks after joining.
	 */
	public void checkFailure(){
		final Throwable problem=failure;
		if(problem!=null){throw new IllegalStateException("Background family snapshot construction failed", problem);}
	}

	/** Owner calls after every query worker has joined, including exceptional exits.
	 * Requests shutdown, wakes the builder and waits for an in-flight build to finish;
	 * it does not cancel a snapshot construction or close the borrowed source.
	 * Safe before start or after a prior close. ThreadWaiter retries interrupted joins;
	 * a retained builder failure is reported after the wait.
	 */
	@Override
	public void close(){
		closing=true;
		LockSupport.unpark(builder);
		if(started){ThreadWaiter.waitForThreadsToFinish(new Thread[]{builder});}
		checkFailure();
	}

	/** Published prefix size, for progress telemetry rather than a query boundary. */
	public int frozenFamilies(){return current.get().size();}
	/** Completed background snapshots, excluding the initial empty view; final after close. */
	public long rebuilds(){return rebuilds;}
	/** Elapsed nanoseconds for completed snapshot construction/publication, excluding waits and queries.
	 * Telemetry getters are separate reads during a run; read after close for final totals.
	 */
	public long buildNanos(){return buildNanos;}

	/*--------------------------------------------------------------*/
	/*----------------       Background Builder     ----------------*/
	/*--------------------------------------------------------------*/

	/** Sole snapshot writer. Atomic publication retains older arrays until their readers finish. */
	private final class BuilderThread extends Thread{
		/** Creates the named worker without starting it. */
		private BuilderThread(){
			super("family-snapshot-builder");
		}

		/** Retains every builder failure for the owner/query-side checks. */
		@Override
		public void run(){
			try{buildRequested();}
			catch(Throwable problem){failure=problem;}
		}

		/** Wait without a monitor or spin loop, then build entirely outside the query path. */
		private void buildRequested(){
			while(!closing){
				if(!requested.getAndSet(false)){
					LockSupport.park(SnapshotFamilyIndex.this);
					if(Thread.interrupted() && !closing){throw new IllegalStateException("Snapshot builder interrupted before close");}
					continue;
				}
				final int count=source.size(), prefix=current.get().size();
				if(count<prefix){throw new IllegalStateException("Published family count shrank during snapshot construction");}
				if(!needsRebuild(count, prefix)){continue;}
				final long start=System.nanoTime();
				final FrozenView next=MutableFamilyIndex.snapshot(width, count, source);
				current.set(next);
				buildNanos+=System.nanoTime()-start;
				rebuilds++;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact number of registers in every source row. */
	private final int width;
	/** Rebuild when tail size is strictly greater than this fraction of the frozen prefix. */
	private final double rebuildFraction;
	/** Borrowed immutable-row source; publication/lifetime remains the caller's responsibility. */
	private final SketchSource source;
	/** Atomically published immutable prefix; each query retains its own captured view. */
	private final AtomicReference<FrozenView> current;
	/** Coalesced wake request; the builder rechecks current growth before rebuilding. */
	private final AtomicBoolean requested=new AtomicBoolean();
	/** Sole writer of snapshots and completed-build counters. */
	private final BuilderThread builder;
	/** Owner-controlled lifecycle flags, visible to query workers and the builder. */
	private volatile boolean started, closing;
	/** Builder failure published for later owner/query checks. */
	private volatile Throwable failure;
	/** Completed-build count and accumulated elapsed nanoseconds; independently visible telemetry. */
	private volatile long rebuilds, buildNanos;
}
