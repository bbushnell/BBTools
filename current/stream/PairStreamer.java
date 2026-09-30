package stream;

import java.util.ArrayList;

import structures.ListNum;

/**
 * Pairs corresponding retained reads from separate R1/R2 streams, reconciling
 * unequal batch sizes with per-side carry buffers. Children must supply compatible
 * ordered records and sampling decisions; pairing asserts matching numeric IDs,
 * not matching read-name strings. Read objects are reused and their mate links set.
 * Configure children before one start. nextList serializes pairing, but this does
 * not synchronize other lifecycle/configuration calls or establish child completion.
 *
 * @author Isla, Shinobu
 * @date October 31, 2025
 */
public class PairStreamer implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains the two children without starting them or copying their state.
	 * Supplying unstarted children is a caller precondition, not a lifecycle check.
	 * Assertions require pair markers 0/1 and non-interleaved child modes.
	 * @param s1_ Nonnull R1 reader with pairnum 0 and paired false
	 * @param s2_ Nonnull R2 reader with pairnum 1 and paired false
	 */
	public PairStreamer(final Streamer s1_, final Streamer s2_){
		s1=s1_;
		s2=s2_;
		assert(s1.pairnum()==0) : "First stream must be R1 (pairnum 0)";
		assert(s2.pairnum()==1) : "Second stream must be R2 (pairnum 1)";
		assert(!s1.paired()) : "First stream should not be interleaved";
		assert(!s2.paired()) : "Second stream should not be interleaved";
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Calls each child's start in R1/R2 order. Call once; local carries, terminal
	 * flags, output IDs and mismatch status are not reset for reuse.
	 */
	@Override
	public void start(){
		s1.start();
		s2.start();
	}

	/** Calls each child's close in R1/R2 order. Adds no join or local-state reset;
	 * completion and interruption behavior remain defined by the children.
	 */
	@Override
	public void close(){
		s1.close();
		s2.close();
	}

	/** @return Current child names joined with a comma */
	@Override
	public String fname(){return s1.fname()+","+s2.fname();}

	/** Combines child availability hints and buffered carry contents. This adapter
	 * neither calls nextList nor removes carries here; delegated behavior is child-defined.
	 * This is not a terminal or worker-completion guarantee; consume nextList to null.
	 * @return Whether either child hints at more input or either carry is nonempty
	 */
	@Override
	public boolean hasMore(){
		//Must reflect the carry-over too, not just s1: a caller driving a loop on
		//hasMore() must not stop while carried-over (already-read, not-yet-returned)
		//reads are still waiting to come out of nextList().
		return s1.hasMore() || s2.hasMore() || !carry1.isEmpty() || !carry2.isEmpty();
	}

	/** @return true; each returned R1 Read is linked to its R2 mate */
	@Override
	public boolean paired(){return true;}

	/** @return 0, the marker for returned root reads (R1) */
	@Override
	public int pairnum(){return 0;}

	/** @return Sum of currently reported child counts, in their units; no finality barrier */
	@Override
	public long readsProcessed(){return s1.readsProcessed()+s2.readsProcessed();}

	/** @return Sum of currently reported child base counts, without waiting for aggregation */
	@Override
	public long basesProcessed(){return s1.basesProcessed()+s2.basesProcessed();}

	/** Forwards sampling to both children before start; this adapter does not resample
	 * paired output. A negative seed chooses one nonnegative seed shared by both calls.
	 * Matching subsets still require compatible child sampling implementations.
	 * @param rate Rate forwarded unchanged, with no validation here
	 * @param seed Shared seed, or a negative value to choose one here
	 */
	@Override
	public void setSampleRate(final float rate, long seed){
		if(seed<0){seed=(long)(Long.MAX_VALUE*Math.random());}
		s1.setSampleRate(rate, seed);//Sample in each child, rather than resampling paired output.
		s2.setSampleRate(rate, seed);
	}

	/**
	 * Returns the next batch of paired reads, reconciling batch-size differences between
	 * s1 and s2 via per-side carry-over buffers (2026-08-21, Nepgear+Amber).
	 *
	 * BACKGROUND: s1/s2 are independent Streamers, each free to return however many
	 * complete records fit in whatever their underlying source handed back on its last
	 * physical read. For plain text this is nearly always aligned between two similarly-
	 * sized files; for gzip it commonly is NOT, because each side is an independently
	 * deflate-compressed stream and GZIPInputStream.read() returns however many bytes
	 * one inflate call happens to produce -- unrelated to record boundaries and unrelated
	 * between s1 and s2. A real (non-adversarial) paired-gzip-fastq run was observed to
	 * desync by one or two records per batch on an otherwise perfectly-valid file. The old
	 * code asserted the two batch sizes were always equal, which is a real, reproducible,
	 * false-positive crash on ordinary gzip input (see bloom/ReadCounter.java's CountThread
	 * for how that manifested).
	 *
	 * FIX: never assume equal batch sizes. Pull into carry1/carry2 only when that side
	 * is empty, retaining at most one current nonempty child batch per side. Its size
	 * is controlled by the child, not a fixed memory bound in this adapter. Take the
	 * shared prefix of length min(carry1.size(), carry2.size()), pair and return that, and
	 * leave any surplus in the carry-over for the next call. Carry reconciliation copies
	 * references; it does not reread or decompress records. Refilling is delegated and
	 * may block while a child produces its next batch.
	 *
	 * Empty nonterminal child batches are skipped until records or an actual null.
	 * A GENUINE mismatch (actually different retained counts between R1/R2 -- a real data
	 * problem, not a batching artifact) still surfaces loud: it shows up as one side
	 * reaching end-of-stream while the other still has an un-drainable carry-over, which
	 * sets mismatchError and asserts; with assertions off that path returns null.
	 * Numeric-ID mismatch is assertion-only: with assertions off it still links the
	 * reads without setting mismatchError. No read-name comparison is performed.
	 * Returned batch/list containers are new, Read objects are reused, and batch IDs
	 * increase locally rather than preserving either child's batch boundaries.
	 * @return Nonempty batch of R1 reads linked to mates, or null at exhaustion or a count mismatch
	 */
	@Override
	public synchronized ListNum<Read> nextList(){
		//STR-035: sampled empty batches are not EOF; refill until records or an actual null.
		while(carry1.isEmpty() && !finished1){
			final ListNum<Read> ln=s1.nextList();
			if(ln==null){finished1=true;}else{carry1.addAll(ln.list);}
		}
		while(carry2.isEmpty() && !finished2){
			final ListNum<Read> ln=s2.nextList();
			if(ln==null){finished2=true;}else{carry2.addAll(ln.list);}
		}

		if(carry1.isEmpty() && carry2.isEmpty()){
			if(finished1!=finished2){mismatchError=true;}
			assert(finished1==finished2) : "Paired files have different read counts! "+fname();
			return null;
		}

		final int n=Math.min(carry1.size(), carry2.size());
		if(n<1){
			//Unreachable except as a genuine mismatch: if carry1 (or carry2) is empty here,
			//the refill block above already tried to top it up and failed (finished1/2==true),
			//while the OTHER side still has leftover reads. Must return null here (not an
			//empty-but-non-null ListNum) or a caller looping on "while(nextList()!=null)" spins
			//forever once assertions are compiled out (-da) and the assert below is a no-op.
			mismatchError=true;
			assert(false) : "Paired files have different read counts (unpaired remnant): "+
				carry1.size()+" vs "+carry2.size()+" in "+fname();
			return null;
		}

		final ArrayList<Read> reads1=new ArrayList<Read>(carry1.subList(0, n));
		final ArrayList<Read> reads2=new ArrayList<Read>(carry2.subList(0, n));
		carry1.subList(0, n).clear();
		carry2.subList(0, n).clear();

		// Mate the reads
		for(int i=0; i<n; i++){
			final Read r1=reads1.get(i);
			final Read r2=reads2.get(i);
			assert(r1.numericID==r2.numericID) : r1.numericID+"!="+r2.numericID+"\n"+r1.id+"\n"+r2.id+"\n";
			r1.mate=r2;
			r2.mate=r1;
		}

		return new ListNum<Read>(reads1, nextListID++);
	}

	/** Rejects the SamLine representation, which this pairing adapter does not supply.
	 * @return Never returns normally
	 * @throws UnsupportedOperationException Always
	 */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("PairStreamer does not support SamLine");
	}

	/** @return Currently reported child errors or local count mismatch; does not await completion */
	@Override
	public boolean errorState(){return s1.errorState() || s2.errorState() || mismatchError;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained R1 reader, configured to report pair marker zero. */
	private final Streamer s1;
	/** Retained R2 reader, configured to report pair marker one. */
	private final Streamer s2;

	/** Unpaired retained R1 references awaiting a matching R2 prefix. */
	private final ArrayList<Read> carry1=new ArrayList<Read>();
	/** Unpaired retained R2 references awaiting a matching R1 prefix. */
	private final ArrayList<Read> carry2=new ArrayList<Read>();
	/** True once s1/s2 has returned null (no more batches). */
	private boolean finished1=false, finished2=false;
	/** Set (in addition to the assert, which may be disabled under -da) on a genuine R1/R2 count
	 * mismatch, so errorState() reports it even if assertions are off. */
	private boolean mismatchError=false;
	/** Sequential id for the ListNum batches this class returns; independent of s1/s2's own ids
	 * since a returned batch's reads may be assembled from more than one underlying batch. */
	private long nextListID=0;

}
