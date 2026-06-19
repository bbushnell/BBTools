package stream;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import dna.Data;
import shared.Shared;
import structures.ListNum;

/**
 * A ConcurrentReadInputStream backed by one or two in-memory List<Read> sources instead of a file.
 * Used to re-stream reads already loaded in RAM: Dedupe/Dedupe2/DedupeProtein, Clumpify (ClumpTools),
 * and IdentityMatrix feed their loaded read list back through the cris API. A single producer thread
 * copies reads from producer1 (paired with producer2 by index, or interleaved with mates already
 * attached) into a ConcurrentDepot, terminated by a poison (empty) list - the same contract as the
 * file-backed cris variants. Live callers always pass (list, null, -1): producer2 is null and maxReads
 * is unlimited, so the two-list index-pairing path is not exercised in practice.
 *
 * NOTE: ConcurrentReadInputStream.start() deliberately launches run() on a FRESH thread ("prevents a
 * strange deadlock in ConcurrentCollectionReadInputStream", ~L289) - the close()/shutdown lifecycle
 * here is delicate; see #001.
 *
 * @author Brian Bushnell
 */
public class ConcurrentCollectionReadInputStream extends ConcurrentReadInputStream {
	
	public ConcurrentCollectionReadInputStream(List<Read> source1, List<Read> source2, long maxReadsToGenerate){
		super("list");
		assert(source1!=source2);
		producer1=source1;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		producer2=source2;
		maxReads=maxReadsToGenerate>=0 ? maxReadsToGenerate : Long.MAX_VALUE;
		if(maxReads==0){
			System.err.println("Warning - created a read stream for 0 reads.");
			assert(false);
		}
		
	}
	
	@Override
	public synchronized ListNum<Read> nextList() {
		ArrayList<Read> list=null;
		if(verbose){System.err.println("**************** nextList() was called; shutdown="+shutdown+", depot.full="+depot.full.size());}
		while(list==null){
			if(shutdown){
				if(verbose){System.err.println("**************** nextList() returning null; shutdown="+shutdown+", depot.full="+depot.full.size());}
				return null;
			}
			try {
				list=depot.full.take();
				assert(list!=null);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		if(verbose){System.err.println("**************** nextList() returning list of size "+list.size()+"; shutdown="+shutdown+", depot.full="+depot.full.size());}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}
	
	@Override
	public void returnList(long listNumber, boolean poison){
		if(poison){
			if(verbose){System.err.println("crisC:    A: Adding empty list to full.");}
			depot.full.add(new ArrayList<Read>(0));
		}else{
			if(verbose){System.err.println("crisC:    A: Adding empty list to empty.");}
			depot.empty.add(new ArrayList<Read>(BUF_LEN));
		}
	}
	
	/**
	 * Main execution method that processes reads from the source collections.
	 * Reads singles/pairs from collections, fills depot buffers, and adds poison
	 * pills when complete. Handles remaining empty buffers during shutdown.
	 */
	@Override
	public void run() {
//		producer.start();
		threads=new Thread[] {Thread.currentThread()};
		if(verbose){System.err.println("crisC started, thread="+threads[0]);}

//		readLists();
		readSingles();

		addPoison();
		
		//End thread

		while(!depot.empty.isEmpty() && !shutdown){
//			System.out.println("Ending");
			if(verbose){System.err.println("B: Adding empty lists to full.");}
			depot.full.add(depot.empty.poll());
		}
//		System.err.println("cris thread terminated. Final depot size: "+depot.full.size()+", "+depot.empty.size());
	}
	
	private final void addPoison(){
		//System.err.println("Adding poison.");
		//Add poison pills
		if(verbose){System.err.println("C: Adding poison to full.");}
		depot.full.add(new ArrayList<Read>());
		for(int i=1; i<depot.bufferCount; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
//					System.err.println("Do not be alarmed by the following error message:");
//					e.printStackTrace();
					if(shutdown){
						i=depot.bufferCount;
						break;
					}
				}
			}
			if(list!=null){
				if(verbose){System.err.println("D: Adding list("+list.size()+") to full "+depot.full.size()+"/"+depot.bufferCount);}
				depot.full.add(list);
			}
		}
		//System.err.println("Added poison.");
	}
	
	private final void readSingles(){

		for(int i=0; !shutdown && i<producer1.size() && generated<maxReads; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=depot.empty.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					if(shutdown){break;}
				}
			}
			if(shutdown || list==null){break;}
			
			long bases=0;
			final long lim=producer1.size();
			while(list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA && generated<lim){
				Read a=producer1.get((int)generated);
				Read b=(producer2==null ? null : producer2.get((int)generated));
				if(a==null){break;}
				readsIn++;
				basesIn+=a.length();
				if(b!=null){
					readsIn++;
					basesIn+=b.length();
				}
				if(randy==null || randy.nextFloat()<samplerate){//Subsampled-IN reads are added; skipped reads still count toward readsIn/basesIn/generated above (they WERE read), only excluded from output.
					list.add(a);//Interleaved-pair convention: only 'a' enters the list; its mate 'b' rides along as a.mate (set just below). b is NEVER added to the list directly.
					if(b!=null){
						assert(a.numericID==b.numericID) : "\n"+a.numericID+", "+b.numericID+"\n"+a.toText(false)+"\n"+b.toText(false)+"\n";
						a.mate=b;
						b.mate=a;

						assert(a.pairnum()==0);
						b.setPairnum(1);
						bases+=(b.bases==null ? 0 : b.length());
					}
					bases+=(a.bases==null ? 0 : a.length());
				}
				incrementGenerated(1);
			}

			if(verbose){System.err.println("E: Adding list("+list.size()+") to full "+depot.full.size()+"/"+depot.bufferCount);}
			//This list is non-empty for every real read; it becomes EMPTY only after generated>=lim/maxReads (all reads already delivered), which the consumer correctly reads as the end-of-stream poison - so no early empty list, no data loss. The post-completion "push empty lists" iterations are bounded by the outer loop's !shutdown guard (the consumer's close() sets shutdown), NOT by producer1.size(): a handful of iterations, never O(reads). [traced - not a spin or data-loss bug]
			depot.full.add(list);
		}
	}
	
	private boolean shutdown=false;
	
	/** Initiates shutdown of the read stream.
	 * Sets shutdown flag and interrupts any running threads. */
	@Override
	public void shutdown(){
		if(verbose){System.out.println("Called shutdown.");}
		shutdown=true;
		//TODO: Possible bug [stream/ConcurrentCollectionReadInputStream#001] - dead branch: shutdown was just set true, so !shutdown is always false and the thread-interrupt body below NEVER runs. Latent LOW: termination still works because close() recycles buffers (unblocking a producer parked in depot.empty.take()) and readSingles breaks on the shutdown flag (~L129) - the interrupt is redundant here. Same template pattern as ConcurrentReadInputStreamD#001. NOT auto-fixed: adding a real interrupt risks the delicate close() lifecycle (cf. the start()-on-fresh-thread deadlock workaround); defer to a deliberate systemic decision.
		if(!shutdown){
			if(verbose){System.out.println("shutdown 2.");}
			for(Thread t : threads){
				if(verbose){System.out.println("shutdown 3.");}
				if(t!=null && t.isAlive()){
					if(verbose){System.out.println("shutdown 4.");}
					t.interrupt();
					if(verbose){System.out.println("shutdown 5.");}
				}
			}
		}
		if(verbose){System.out.println("shutdown 6.");}
	}
	
	/** Resets the stream to initial state for reuse.
	 * Clears shutdown flag, recreates depot, and resets counters. */
	@Override
	public synchronized void restart(){
		shutdown=false;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		generated=0;
		basesIn=0;
		readsIn=0;
		nextProgress=PROGRESS_INCR;
	}
	
	@Override
	public synchronized void close(){
		if(verbose){System.out.println("Thread "+Thread.currentThread().getId()+" called close.");}
		shutdown();
//		producer1.close();
//		if(producer2!=null){producer2.close();}
//		System.out.println("A");
		if(threads!=null && threads[0]!=null && threads[0].isAlive()){
			if(verbose){System.out.println("close 1.");}
			
			while(threads[0].isAlive()){
				if(verbose){System.out.println("close 2: Thread "+Thread.currentThread().getId()+" closing thread "+threads[0].getId()+" "+threads[0].getState());}
//				System.out.println("B");
				ArrayList<Read> list=null;
				for(int i=0; i<1 && list==null && threads[0].isAlive(); i++){
					if(verbose){System.out.println("close 3.");}
					try {
						if(verbose){System.out.println("close 4.");}
						list=depot.full.poll(100, TimeUnit.MILLISECONDS);
						if(verbose){System.out.println("close 5; list.size()="+depot.full.size()+", list="+(list==null ? "null" : list.size()+""));}
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						System.err.println("Do not be alarmed by the following error message:");
						e.printStackTrace();
						break;
					}
				}
				
				if(list!=null){
					list.clear();
					depot.empty.add(list);
				}
				if(verbose){System.out.println("close 6.");}
				
//				System.out.println("isAlive? "+threads[0].isAlive());
			}
			if(verbose){System.out.println("close 7.");}
			
		}
		if(verbose){System.out.println("close 8.");}
		
		if(threads!=null){
			if(verbose){System.out.println("close 9.");}
			for(int i=1; i<threads.length; i++){
				if(verbose){System.out.println("close 10.");}
				while(threads[i]!=null && threads[i].isAlive()){
					if(verbose){System.out.println("close 11.");}
					try {
						if(verbose){System.out.println("close 12.");}
						threads[i].join();
						if(verbose){System.out.println("close 13.");}
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		if(verbose){System.out.println("close 14.");}
		
	}

	@Override
	public boolean paired() {//Paired if a second list was given; otherwise infer from producer1: empty->unpaired, else interleaved iff its first read has a mate attached.
		return producer2!=null ? true : (producer1==null || producer1.isEmpty() ? false : producer1.get(0).mate!=null);
	}
	
	@Override
	public boolean verbose(){return verbose;}
	
	private void incrementGenerated(long amt){
		generated+=amt;
		if(SHOW_PROGRESS && generated>=nextProgress){
			Data.sysout.print('.');
			nextProgress+=PROGRESS_INCR;
		}
	}
	
	/**
	 * Sets the sampling rate for read selection.
	 * Creates random number generator for subsampling if rate is less than 1.0.
	 * @param rate Fraction of reads to keep (0.0 to 1.0)
	 * @param seed Random seed for reproducible sampling, or negative for random seed
	 */
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		if(rate>=1f){
			randy=null;
		}else{
			randy=Shared.threadLocalRandom(seed);
		}
	}
	
	/** Returns total number of bases processed from input.
	 * @return Total bases read from source collections */
	@Override
	public long basesIn(){return basesIn;}
	/** Returns total number of reads processed from input.
	 * @return Total reads processed from source collections */
	@Override
	public long readsIn(){return readsIn;}
	
	/** Returns current error state of the stream.
	 * @return true if an error has occurred, false otherwise */
	@Override
	public boolean errorState(){return errorState;}
	/** Flag tracking error state of the stream */
	private boolean errorState=false;
	
	private float samplerate=1f;
	private shared.Random randy=null;
	
	private Thread[] threads;
	
	@Override
	public Object[] producers(){return new Object[] {producer1, producer2};}
	
	public final List<Read> producer1;
	public final List<Read> producer2;
	private ConcurrentDepot<Read> depot;
	
	private long basesIn=0;
	private long readsIn=0;
	
	private long maxReads;
	private long generated=0;
	private long listnum=0;
	private long nextProgress=PROGRESS_INCR;
	
	public static boolean verbose=false;
	
	private static final ArrayList<Read> poison=new ArrayList<Read>(0);
	
}
