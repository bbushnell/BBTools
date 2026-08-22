package stream;

import java.util.ArrayList;

import structures.ListNum;

/**
 * Pairs reads from two separate Streamers (typically R1 and R2 files).
 * Ensures mate references are set correctly.
 * 
 * @author Isla
 * @date October 31, 2025
 */
public class PairStreamer implements Streamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public PairStreamer(Streamer s1_, Streamer s2_){
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
	
	@Override
	public void start(){
		s1.start();
		s2.start();
	}
	
	@Override
	public void close(){
		s1.close();
		s2.close();
	}
	
	@Override
	public String fname() {return s1.fname()+","+s2.fname();}
	
	@Override
	public boolean hasMore(){
		//Must reflect the carry-over too, not just s1: a caller driving a loop on
		//hasMore() must not stop while carried-over (already-read, not-yet-returned)
		//reads are still waiting to come out of nextList().
		return s1.hasMore() || s2.hasMore() || !carry1.isEmpty() || !carry2.isEmpty();
	}
	
	@Override
	public boolean paired(){return true;}
	
	@Override
	public int pairnum(){return 0;} // Paired data returns as R1
	
	@Override
	public long readsProcessed(){
		return s1.readsProcessed()+s2.readsProcessed();
	}
	
	@Override
	public long basesProcessed(){
		return s1.basesProcessed()+s2.basesProcessed();
	}
	
	@Override
	public void setSampleRate(float rate, long seed){
		if(seed<0) {seed=(long)(Long.MAX_VALUE*Math.random());}
		s1.setSampleRate(rate, seed);//Much faster to pass through than process in this stream
		s2.setSampleRate(rate, seed);
	}
	
	/**
	 * Returns the next batch of paired reads, reconciling batch-size mismatches between
	 * s1 and s2 via a small per-side carry-over buffer (2026-08-21, Nepgear+Amber).
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
	 * FIX: never assume equal batch sizes. Pull into carry1/carry2 (only refilling a side
	 * whose carry-over is currently empty, so we never grow it unboundedly), take the
	 * shared prefix of length min(carry1.size(), carry2.size()), pair and return that, and
	 * leave any surplus in the carry-over for the next call. This is pure list bookkeeping
	 * on data already in memory -- no extra I/O, no extra decompression -- so it does not
	 * change how much is read, only how it's chunked before pairing.
	 *
	 * A GENUINE mismatch (actually different read counts between R1/R2 -- a real data
	 * problem, not a batching artifact) still surfaces loud: it shows up as one side
	 * reaching end-of-stream while the other still has an un-drainable carry-over, which
	 * hits the assert below instead of being silently dropped.
	 */
	@Override
	public synchronized ListNum<Read> nextList(){
		if(carry1.isEmpty() && !finished1){
			ListNum<Read> ln=s1.nextList();
			if(ln==null){finished1=true;}else{carry1.addAll(ln.list);}
		}
		if(carry2.isEmpty() && !finished2){
			ListNum<Read> ln=s2.nextList();
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

		ArrayList<Read> reads1=new ArrayList<Read>(carry1.subList(0, n));
		ArrayList<Read> reads2=new ArrayList<Read>(carry2.subList(0, n));
		carry1.subList(0, n).clear();
		carry2.subList(0, n).clear();

		// Mate the reads
		for(int i=0; i<n; i++){
			Read r1=reads1.get(i);
			Read r2=reads2.get(i);
			assert(r1.numericID==r2.numericID) : r1.numericID+"!="+r2.numericID+"\n"+r1.id+"\n"+r2.id+"\n";
			r1.mate=r2;
			r2.mate=r1;
		}

//		// Apply subsampling if needed
//		if(samplerate<1f && randy!=null){
//			int nulled=0;
//			for(int i=0; i<reads1.size(); i++){
//				if(randy.nextFloat()>=samplerate){
//					reads1.set(i, null);
//					nulled++;
//				}
//			}
//			if(nulled>0) {Tools.condenseStrict(reads1);}
//		}

		return new ListNum<Read>(reads1, nextListID++);
	}
	
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("PairStreamer does not support SamLine");
	}
	
	@Override
	public boolean errorState(){return s1.errorState() || s2.errorState() || mismatchError;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Streamer s1; // R1
	private final Streamer s2; // R2

	/** Reads pulled from s1/s2 but not yet returned by nextList() -- the batch-size reconciliation buffer. */
	private final ArrayList<Read> carry1=new ArrayList<Read>();
	private final ArrayList<Read> carry2=new ArrayList<Read>();
	/** True once s1/s2 has returned null (no more batches). */
	private boolean finished1=false, finished2=false;
	/** Set (in addition to the assert, which may be disabled under -da) on a genuine R1/R2 count
	 * mismatch, so errorState() reports it even if assertions are off. */
	private boolean mismatchError=false;
	/** Sequential id for the ListNum batches this class returns; independent of s1/s2's own ids
	 * since a returned batch's reads may be assembled from more than one underlying batch. */
	private long nextListID=0;
//	private float samplerate=1f;
//	private shared.Random randy=null;
	
}