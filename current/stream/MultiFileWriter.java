package stream;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import cardinality.CardinalityTracker;
import map.LongObjectMap;
import structures.ByteBuilder;
import structures.HeapLoc;
import structures.SetLoc;

/**
 * Routes Read batches to ordinary Writers under a physical-output limit.
 * Targets own file formats and reopening; this container never serializes a read.
 * Small destination sets remain open. Larger sets retire the least recently used
 * reopenable target, waiting for its Writer before releasing the output slots.
 * No dispatcher thread is created by this container.
 *
 * Submission methods serialize access. Unordered batches preserve serialized
 * submission order to each leaf; physical write scheduling belongs to that leaf.
 * Ordered submission requires one complete ArrayListSet per dense input ID,
 * including empty batches, beginning at the configured first ID. Future batches
 * may wait for byte/batch capacity; the caller must allow the missing earlier ID
 * to arrive from another producer. Leaf IDs are dense from zero on every reopen.
 *
 * Read payloads are borrowed and must remain stable until their output completes.
 * Routing buffers use Read.countPairBytes estimates; these are not measured heap
 * bounds. Ordered mode has a separate reorder budget of the same size, plus a
 * 64-batch cap. Caller-held submissions, leaf queues/builders and target metadata
 * are outside those estimates. A single oversized read is processed as one unit.
 * Below-threshold data cannot be discarded to relieve pressure: increase the
 * budget or reduce the threshold if no eligible buffer can be flushed.
 *
 * Finish producers before close. Close flushes eligible destinations and waits
 * for open Writers; below-threshold reads remain for an explicit residual drain.
 * Routing counts include attached mates and repeated assignments. Output getters
 * instead sum the concrete Writers' reported counts, including their format-specific
 * treatment of headers/alignments. Neither counter is a durability guarantee.
 * Optional cardinality tracking adds a per-destination estimate to report rows.
 * Trackers persist across output generations and are outside the routing budget.
 * Configure global CardinalityTracker defaults before routing and keep them stable
 * while new destinations can be created, as with legacy multi-output tracking.
 *
 * @author Shinobu
 * @date October 1, 2026
 */
public final class MultiFileWriter{

	/*--------------------------------------------------------------*/
	/*----------------        Target Contracts      ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves names without opening outputs. Distinct names must use disjoint outputs. */
	public interface TargetFactory{
		/** @param name Nonnull destination name
		 * @return Stable, nonnull target for this logical destination */
		Target create(String name);
	}

	/** Owns serialization choices and the lifecycle of one logical destination. */
	public interface Target{
		/** @return Positive, stable number of simultaneously open physical outputs */
		int outputCount();
		/** @return Whether this target may be retired and later reopened without data loss */
		boolean canReopen();
		/** Creates an unstarted ST/ZT Writer with no hidden MT backend.
		 * @param reopen False for the first open; true requires append without another file header
		 * @return New Writer for this generation; the target owns startup cleanup if opening fails */
		Writer open(boolean reopen);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an unordered container with 256-entry, at most 1MB estimated batches.
	 * @param targets_ Destination resolver
	 * @param maxOpenOutputs_ Physical output cap, including mates/QUAL
	 * @param bufferBytes_ Estimated routing-buffer budget in bytes */
	public MultiFileWriter(final TargetFactory targets_, final int maxOpenOutputs_, final long bufferBytes_){
		this(targets_, maxOpenOutputs_, bufferBytes_, 256, Math.min(1_000_000L, bufferBytes_), 0, -1);
	}

	/** Creates a container with explicit buffering, threshold and ordering policies.
	 * @param targets_ Destination resolver
	 * @param maxOpenOutputs_ Positive physical-output cap
	 * @param bufferBytes_ Positive routing budget; ordered mode has an additional equal reorder budget
	 * @param batchReads_ Positive maximum buffered list entries before submission
	 * @param batchBytes_ Positive estimated bytes per batch before submission
	 * @param minReads_ Minimum cumulative routed individual reads before creating a destination
	 * @param firstListId Negative for unordered mode; otherwise the first complete input-batch ID */
	public MultiFileWriter(final TargetFactory targets_, final int maxOpenOutputs_, final long bufferBytes_,
			final int batchReads_, final long batchBytes_, final long minReads_, final long firstListId){
		this(targets_, maxOpenOutputs_, bufferBytes_, batchReads_, batchBytes_, minReads_, firstListId, false);
	}

	/** Creates a container with optional per-destination cardinality tracking.
	 * Tracking is fixed for this container; tracker defaults are read when each destination is first observed.
	 * @param targets_ Destination resolver
	 * @param maxOpenOutputs_ Positive physical-output cap
	 * @param bufferBytes_ Positive routing budget; tracker storage is additional
	 * @param batchReads_ Positive maximum buffered list entries before submission
	 * @param batchBytes_ Positive estimated bytes per batch before submission
	 * @param minReads_ Minimum cumulative routed individual reads before creating a destination
	 * @param firstListId Negative for unordered mode; otherwise the first complete input-batch ID
	 * @param trackCardinality_ Track paired input sequences using the configured CardinalityTracker defaults */
	public MultiFileWriter(final TargetFactory targets_, final int maxOpenOutputs_, final long bufferBytes_,
			final int batchReads_, final long batchBytes_, final long minReads_, final long firstListId,
			final boolean trackCardinality_){
		if(targets_==null || maxOpenOutputs_<1 || bufferBytes_<1 || batchReads_<1 || batchBytes_<1 || minReads_<0){
			throw new IllegalArgumentException("Targets and positive output/buffer limits are required; minReads must be nonnegative");
		}
		targets=targets_;
		maxOpenOutputs=maxOpenOutputs_;
		bufferBytes=bufferBytes_;
		batchReads=batchReads_;
		batchBytes=Math.min(batchBytes_, bufferBytes_);
		minReads=minReads_;
		trackCardinality=trackCardinality_;
		ordered=(firstListId>=0);
		nextListId=firstListId;
		pending=(ordered ? new LongObjectMap<Batch>(16, Batch.class) : null);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Submission          ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks that submission is still permitted; this container starts no threads. */
	public synchronized void start(){ensureOpen();}

	/** Routes one root and its linked mate to a name in unordered mode.
	 * @param read Borrowed read; null is ignored
	 * @param name Destination; null is ignored */
	public synchronized void add(final Read read, final String name){
		ensureUnordered();
		try{route(read, name);}catch(RuntimeException|Error e){fail(); throw e;}
	}

	/** Routes a named batch in unordered mode; copies references into destination buffers.
	 * @param reads Borrowed payloads; the input list may be reused after return
	 * @param name Destination; null is ignored */
	public synchronized void add(final ArrayList<Read> reads, final String name){
		ensureUnordered();
		try{route(reads, name);}catch(RuntimeException|Error e){fail(); throw e;}
	}

	/** Routes an unordered batch by String-valued Read.obj; null names/reads are ignored.
	 * @param reads Borrowed payloads; the input list may be reused after return */
	public synchronized void add(final ArrayList<Read> reads){
		ensureUnordered();
		try{
			if(reads!=null){for(Read r : reads){if(r!=null && r.obj!=null){route(r, (String)r.obj);}}}
		}catch(RuntimeException|Error e){fail(); throw e;}
	}

	/**
	 * Detaches and submits one complete grouped batch. Ordered IDs must be unique
	 * and dense from the configured first ID; unordered mode ignores the ID.
	 * The caller may reuse ArrayListSet after return, but must keep Read payloads stable.
	 * @param groups Named lists to detach; null represents an empty complete batch
	 * @param id Input batch ID, including empty batches in ordered mode
	 */
	public synchronized void add(final ArrayListSet groups, final long id){
		ensureOpen();
		try{
			if(!ordered){route(groups); return;}
			checkId(id);
			if(id==nextListId){route(groups); advance(); return;}
			final Batch batch=new Batch(groups);
			while(id!=nextListId && (pending.size()>=MAX_PENDING_BATCHES || batch.bytes>bufferBytes-pendingBytes)){
				try{wait();}catch(InterruptedException e){Thread.currentThread().interrupt(); throw new RuntimeException("Interrupted waiting for ordered output capacity", e);}
				ensureOpen();
				checkId(id);
			}
			if(id==nextListId){route(batch); advance();}
			else{pending.put(id, batch); pendingBytes+=batch.bytes;}
		}catch(RuntimeException|Error e){fail(); throw e;}
	}

	/** Advances the dense input sequence and drains available consecutive batches. */
	private void advance(){
		nextListId=Math.addExact(nextListId, 1);
		for(Batch batch=pending.remove(nextListId); batch!=null; batch=pending.remove(nextListId)){
			pendingBytes-=batch.bytes;
			route(batch);
			nextListId=Math.addExact(nextListId, 1);
		}
		notifyAll();
	}

	/** Rejects duplicate or already-consumed ordered batch IDs. */
	private void checkId(final long id){
		if(id<nextListId || pending.get(id)!=null){throw new IllegalArgumentException("Duplicate/old ordered batch ID "+id+", next="+nextListId);}
	}

	/** Routes detached complete groups without modifying Read.obj. */
	private void route(final Batch batch){for(Group group : batch.groups){route(group.reads, group.name);}}

	/** Drains named groups directly when no reordering storage is needed. */
	private void route(final ArrayListSet groups){
		if(groups==null){return;}
		for(String name : groups.getNames()){if(name!=null){route(groups.getAndClear(name), name);}}
	}

	/** Routes a list without retaining its list object. */
	private void route(final ArrayList<Read> reads, final String name){
		if(reads!=null && name!=null){for(Read r : reads){route(r, name);}}
	}

	/** Accounts for one assignment and enforces batch/routing-buffer thresholds. */
	private void route(final Read read, final String name){
		if(read==null || name==null){return;}
		Buffer b=buffers.get(name);
		if(b==null){
			final Target target=targets.create(name);
			if(target==null || target.outputCount()<1 || target.outputCount()>maxOpenOutputs){
				throw new IllegalArgumentException("Destination cannot fit the physical-output cap: "+name);
			}
			b=new Buffer(name, target, trackCardinality ? CardinalityTracker.makeTracker() : null);
			buffers.put(name, b);
		}
		if(b.reads==null){b.reads=new ArrayList<Read>(Math.min(16, batchReads));}
		b.reads.add(read);
		final long bytes=read.countPairBytes();
		assert(bytes>=0) : "Read.countPairBytes must be nonnegative for routing-budget accounting";
		b.bytes+=bytes;
		bufferedBytes+=bytes;
		b.readsIn+=read.pairCount();
		b.basesIn+=(long)read.length()+read.mateLength();
		if(b.cardinality!=null){b.cardinality.hash(read);}//The tracker already includes the mate.
		if(b.readsIn>=minReads){
			if(b.loc()<0){
				if(!heap.hasRoom()){heap=heap.resizeNew(Math.addExact(Math.multiplyExact(heap.size(), 2), 1));}
				heap.add(b);
			}else{heap.jiggle(b);}
			if(b.reads.size()>=batchReads || b.bytes>=batchBytes){flush(b);}
		}
		while(bufferedBytes>bufferBytes){
			final Buffer largest=heap.peek();
			if(largest==null || largest.bytes==0){
				throw new IllegalStateException("Routing budget filled by destinations below minReads; increase the budget or reduce minReads");
			}
			flush(largest);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Writer Residency       ----------------*/
	/*--------------------------------------------------------------*/

	/** Submits a nonempty eligible buffer with a dense generation-local leaf ID. */
	private void flush(final Buffer b){
		if(b.bytes==0){return;}
		assert(b.readsIn>=minReads && b.reads!=null && !b.reads.isEmpty()) : "Only eligible, nonempty buffers may reach a Writer";
		final Writer writer=writer(b);
		writer.add(b.reads, b.nextId);
		b.nextId++;
		bufferedBytes-=b.bytes;
		b.reads=null;
		b.bytes=0;
		heap.jiggle(b);
	}

	/** Opens a target after retiring enough reopenable outputs, or touches an existing one. */
	private Writer writer(final Buffer b){
		if(b.writer!=null){open.get(b.name); return b.writer;}
		while(b.outputs>maxOpenOutputs-openOutputs){retireOne();}
		final Writer writer=b.target.open(b.opened);
		if(writer==null){throw new IllegalStateException("Target returned no Writer: "+b.name);}
		b.writer=writer;
		b.opened=true;
		b.nextId=0;
		open.put(b.name, b);
		openOutputs+=b.outputs;
		peakOpenOutputs=Math.max(peakOpenOutputs, openOutputs);
		writersOpened++;
		assert(openOutputs<=maxOpenOutputs) : "Writer creation must respect the physical-output cap";
		writer.start();
		return writer;
	}

	/** Retires the least recently used target that can later reopen without data loss. */
	private void retireOne(){
		for(Iterator<Map.Entry<String, Buffer>> it=open.entrySet().iterator(); it.hasNext();){
			final Buffer b=it.next().getValue();
			if(b.target.canReopen()){
				finish(b);
				it.remove();
				rotations++;
				if(error){throw new IllegalStateException("Writer reported an error during retirement: "+b.name);}
				return;
			}
		}
		throw new IllegalStateException("Output cap reached with only nonreopenable targets; increase the cap");
	}

	/** Waits for one leaf and retains its final counters before releasing output slots. */
	private void finish(final Buffer b){
		assert(b.writer!=null) : "Only open Writers may release physical-output slots";
		error|=b.writer.poisonAndWait();
		b.readsOut+=b.writer.readsWritten();
		b.basesOut+=b.writer.basesWritten();
		b.writer=null;
		openOutputs-=b.outputs;
		assert(openOutputs>=0) : "Each opened physical output must be released once";
	}

	/*--------------------------------------------------------------*/
	/*----------------       Finish and Residuals    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Finishes eligible output and waits for every open leaf. Call after producers
	 * finish. Below-threshold lists remain for dumpResidual. Repeated calls return
	 * the stored error indication; a missing ordered batch ID is an explicit error.
	 * @return True if a leaf reported an error
	 */
	public synchronized boolean close(){
		if(closed){return error;}
		Throwable problem=null;
		try{
			if(ordered && !pending.isEmpty()){throw new IllegalStateException("Missing ordered input batch "+nextListId+" before close");}
			if(error){throw new IllegalStateException("Cannot finish normal output after a submission error");}
			for(Buffer b=heap.peek(); b!=null && b.bytes>0; b=heap.peek()){flush(b);}
		}catch(RuntimeException|Error e){fail(); problem=e;}
		closed=true;
		for(Iterator<Buffer> it=open.values().iterator(); it.hasNext();){
			final Buffer b=it.next();
			try{finish(b); it.remove();}
			catch(RuntimeException|Error e){
				fail();
				if(problem==null){problem=e;}else if(problem!=e){problem.addSuppressed(e);}
			}
		}
		notifyAll();
		if(problem instanceof Error){throw (Error)problem;}
		if(problem!=null){throw (RuntimeException)problem;}
		return error;
	}

	/** @return The same completion/error result as close */
	public boolean poisonAndWait(){return close();}

	/**
	 * Transfers below-threshold lists to a caller-owned residual Writer after close,
	 * or discards them if sink is null. Does not start/close the sink. The supplied
	 * first ID must continue that sink's dense sequence if it requires ordering.
	 * Repeated draining does nothing. Counters measure routed reads, including mates.
	 * @param sink Residual recipient, or null to discard deliberately
	 * @param firstId First residual batch ID accepted by the sink
	 * @return Number of individual reads drained by this call
	 */
	public synchronized long dumpResidual(final Writer sink, final long firstId){
		if(!closed || error || firstId<0){throw new IllegalStateException("Residual draining requires successful close and a nonnegative first ID");}
		long count=0, id=firstId;
		try{
			for(Buffer b : buffers.values()){
				if(b.reads==null){continue;}
				assert(b.readsIn<minReads && !b.opened) : "Only never-emitted below-threshold destinations may be residuals";
				if(sink!=null){sink.add(b.reads, id++);}
				count+=b.readsIn;
				residualReads+=b.readsIn;
				residualBases+=b.basesIn;
				bufferedBytes-=b.bytes;
				b.bytes=0;
				b.reads=null;
			}
		}catch(RuntimeException|Error e){fail(); throw e;}
		return count;
	}

	/** @return Snapshot of all observed nonnull destination names, including residuals */
	public synchronized Set<String> getKeys(){return new LinkedHashSet<String>(buffers.keySet());}
	/** @return Current physical-output count */
	public synchronized int openOutputs(){return openOutputs;}
	/** @return Maximum physical-output count observed so far */
	public synchronized int peakOpenOutputs(){return peakOpenOutputs;}
	/** @return Number of leaf instances opened, including reopenings */
	public synchronized long writersOpened(){return writersOpened;}
	/** @return Number of retirements to free slots, excluding terminal close */
	public synchronized long rotations(){return rotations;}
	/** @return Estimated bytes in destination read buffers, excluding leaf queues */
	public synchronized long bufferedBytes(){return bufferedBytes;}
	/** @return Estimated bytes in accepted out-of-order input batches */
	public synchronized long pendingBytes(){return pendingBytes;}
	/** @return Individual routed reads drained as residuals */
	public synchronized long residualReads(){return residualReads;}
	/** @return Routed bases drained as residuals */
	public synchronized long residualBases(){return residualBases;}
	/** @return Sticky submission or leaf-completion error indication */
	public synchronized boolean errorState(){return error;}
	/** @return Whether normal output is closed without a reported error */
	public synchronized boolean finishedSuccessfully(){return closed && !error;}

	/** @return Sum of concrete Writer-reported read counts, including live snapshots */
	public synchronized long readsWritten(){
		long sum=0;
		for(Buffer b : buffers.values()){sum+=b.readsOut+(b.writer==null ? 0 : b.writer.readsWritten());}
		return sum;
	}

	/** @return Sum of concrete Writer-reported base counts, including live snapshots */
	public synchronized long basesWritten(){
		long sum=0;
		for(Buffer b : buffers.values()){sum+=b.basesOut+(b.writer==null ? 0 : b.writer.basesWritten());}
		return sum;
	}

	/** Returns routed input counts for destinations that were actually opened.
	 * Estimates, when enabled, include all routed sequences for that destination across rotations.
	 * Below-threshold destinations and caller-owned residual summary rows are not included.
	 * @return Tab-separated Name, Reads, Bases and optional Cardinality rows, without a header */
	public synchronized ByteBuilder report(){
		final ByteBuilder bb=new ByteBuilder();
		for(Buffer b : buffers.values()){
			if(b.opened){
				bb.append(b.name).tab().append(b.readsIn).tab().append(b.basesIn);
				if(b.cardinality!=null){bb.tab().append(b.cardinality.cardinality());}
				bb.nl();
			}
		}
		return bb;
	}

	/** Rejects submission after closure/error. */
	private void ensureOpen(){
		if(closed || error){throw new IllegalStateException("Multi-file writer is closed or failed");}
	}
	/** Prevents bypassing the complete-batch ordering contract. */
	private void ensureUnordered(){
		ensureOpen();
		if(ordered){throw new IllegalStateException("Ordered mode requires complete ArrayListSet batches and IDs");}
	}
	/** Publishes an error to ordered producers waiting for capacity. */
	private void fail(){error=true; notifyAll();}

	/*--------------------------------------------------------------*/
	/*----------------        Private State         ----------------*/
	/*--------------------------------------------------------------*/

	/** Destination bookkeeping, kept after retirement without an eager empty list. */
	private static final class Buffer implements SetLoc<Buffer>{
		/** Captures stable target metadata without opening a Writer. */
		Buffer(final String name_, final Target target_, final CardinalityTracker cardinality_){
			name=name_; target=target_; outputs=target.outputCount(); cardinality=cardinality_;
		}
		/** Largest pending byte count sorts first in HeapLoc's min-heap. */
		@Override
		public int compareTo(final Buffer other){return Long.compare(other.bytes, bytes);}
		@Override
		public void setLoc(final int value){location=value;}
		@Override
		public int loc(){return location;}
		/** Logical name and stable writer-opening recipe. */
		final String name;
		final Target target;
		/** Optional input-sequence estimate, retained independently of writer generations. */
		final CardinalityTracker cardinality;
		/** Physical output slots consumed by an open generation. */
		final int outputs;
		/** Pending borrowed read references; null after submission or residual draining. */
		ArrayList<Read> reads;
		/** Current generation, or null while closed/retired. */
		Writer writer;
		/** Pending bytes, routed totals, completed-generation output totals and next leaf ID. */
		long bytes, readsIn, basesIn, readsOut, basesOut, nextId;
		/** Whether any generation has opened, determining append on the next one. */
		boolean opened;
		/** HeapLoc position, or -1 before the destination reaches its threshold. */
		int location=-1;
	}

	/** Detached named list belonging to one complete input batch. */
	private static final class Group{
		/** Stores an already detached list without copying its payloads. */
		Group(final String name_, final ArrayList<Read> reads_){name=name_; reads=reads_;}
		/** Destination of this detached list. */
		final String name;
		/** Borrowed immutable payloads retained until this batch can be routed. */
		final ArrayList<Read> reads;
	}

	/** Retains only populated string groups, not the caller's reusable ArrayListSet. */
	private static final class Batch{
		/** Detaches populated named lists and estimates retained bytes, counting repeated assignments. */
		Batch(final ArrayListSet input){
			long size=0;
			if(input!=null){
				for(String name : input.getNames()){
					if(name==null){continue;}
					final ArrayList<Read> reads=input.getAndClear(name);
					if(reads==null || reads.isEmpty()){continue;}
					groups.add(new Group(name, reads));
					for(Read r : reads){if(r!=null){size=Math.addExact(size, r.countPairBytes());}}
				}
			}
			bytes=size;
		}
		/** Populated string groups in caller order. */
		final ArrayList<Group> groups=new ArrayList<Group>();
		/** Estimated retained Read bytes, counting repeated assignments. */
		final long bytes;
	}

	/** Immutable destination resolver; called only while holding the container monitor. */
	private final TargetFactory targets;
	/** Physical-output cap and per-batch entry threshold. */
	private final int maxOpenOutputs, batchReads;
	/** Routing byte budget, per-batch byte threshold and cumulative individual-read threshold. */
	private final long bufferBytes, batchBytes, minReads;
	/** Whether complete input batches must be restored to their dense ID order. */
	private final boolean ordered;
	/** Whether newly observed destinations receive a cardinality tracker. */
	private final boolean trackCardinality;
	/** All observed destinations, retained in first-observation order for reports/residuals. */
	private final LinkedHashMap<String, Buffer> buffers=new LinkedHashMap<String, Buffer>();
	/** Open generations in least-recently-used submission order. */
	private final LinkedHashMap<String, Buffer> open=new LinkedHashMap<String, Buffer>(16, 0.75f, true);
	/** Qualified destinations; largest pending byte count sorts first, including zero-byte entries. */
	private HeapLoc<Buffer> heap=new HeapLoc<Buffer>(15, false);
	/** Complete future batches, guarded by this container's monitor; null when unordered. */
	private final LongObjectMap<Batch> pending;
	/** Routing/reorder occupancy, ordering cursor, generation counts and terminal residual totals. */
	private long bufferedBytes, pendingBytes, nextListId, writersOpened, rotations, residualReads, residualBases;
	/** Current and peak physical-output counts. */
	private int openOutputs, peakOpenOutputs;
	/** Completion flag and sticky error indication. */
	private boolean closed, error;
	/** Bounds metadata even for arbitrarily many empty future input batches. */
	private static final int MAX_PENDING_BATCHES=64;
}
