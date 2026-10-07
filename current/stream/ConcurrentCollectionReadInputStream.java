package stream;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import dna.Data;
import shared.Shared;
import structures.ListNum;

/**
 * Streams borrowed in-memory Read references through the ConcurrentReadInputStream API.
 * A primary list supplies fragments with optional attached mates. An optional second
 * list supplies mates by matching index; selected entries have their mate fields linked
 * and the second read's pair number set to one. Read objects are not copied.
 * Keep source-list structure and shared configuration stable while streaming.
 * Current construction sites in Dedupe/Dedupe2/DedupeProtein, ClumpTools and IdentityMatrix
 * pass (list, null, -1); this observation does not establish two-list behavior by testing.
 *
 * Use inherited start to launch a separate producer thread. Its historical fresh-thread
 * rationale is retained in ConcurrentReadInputStream.start; lifecycle note #001 remains
 * below. Consumers recognize empty lists as terminal and return each batch appropriately.
 * Coordinate lifecycle calls; restart retains some state and does not start another run.
 * @author Brian Bushnell
 */
public class ConcurrentCollectionReadInputStream extends ConcurrentReadInputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains source lists and creates a depot using inherited captured buffer settings.
	 * Source lists must be distinct; primary entries must be nonnull. A secondary list
	 * must cover each accessed index; present mates must have matching numeric IDs.
	 * Selected primaries with an explicit secondary must have pairnum zero.
	 * @param source1 Nonnull primary list, with optional mates already attached
	 * @param source2 Optional second list of mates, paired by index when selected
	 * @param maxReadsToGenerate Fragment limit before sampling; negative means unlimited,
	 * zero prints a warning and is rejected with assertions enabled
	 */
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

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Takes the next batch and assigns an increasing list ID, including empty terminals.
	 * Uses NORMAL ListNum wrappers; terminal status is conveyed by an empty list rather
	 * than its poison type flag. Wrapping may assign Read.rand under ListNum's global option.
	 * Interrupted waits print a diagnostic and retry; shutdown is checked before taking.
	 * @return Numbered borrowed Read references, or null when shutdown is observed
	 */
	@Override
	public synchronized ListNum<Read> nextList(){
		ArrayList<Read> list=null;
		if(verbose){System.err.println("**************** nextList() was called; shutdown="+shutdown+", depot.full="+depot.full.size());}
		while(list==null){
			if(shutdown){
				if(verbose){System.err.println("**************** nextList() returning null; shutdown="+shutdown+", depot.full="+depot.full.size());}
				return null;
			}
			try{
				list=depot.full.take();
				assert(list!=null);
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		if(verbose){System.err.println("**************** nextList() returning list of size "+list.size()+"; shutdown="+shutdown+", depot.full="+depot.full.size());}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}

	/** Adds a replacement list rather than recycling or clearing the caller's batch.
	 * The inherited ListNum overload forwards only ID and emptiness to this method.
	 * @param listNumber Accepted but ignored
	 * @param poison Add an empty terminal to full when true, otherwise a fresh buffer to empty
	 */
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

	/** Records the executing thread, submits selected references, then terminal lists.
	 * Completion depends on consumer participation in the buffer-return protocol.
	 * This method does not latch thrown failures into the local error flag.
	 */
	@Override
	public void run(){
//		producer.start();
		threads=new Thread[]{Thread.currentThread()};
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

	/** Submits a new empty terminal and obtains further lists from the available queue.
	 * Retains the existing polling/interruption behavior; no prompt-return guarantee.
	 */
	private final void addPoison(){
		//System.err.println("Adding poison.");
		//Add poison pills
		if(verbose){System.err.println("C: Adding poison to full.");}
		depot.full.add(new ArrayList<Read>());
		for(int i=1; i<depot.bufferCount; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try{
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				}catch(InterruptedException e){
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

	/** Fills batches by primary index, applying the fragment limit before sampling.
	 * Counters include each primary and its explicit secondary or attached mate before selection.
	 * Batch data includes both members of selected fragments; a whole fragment may cross MAX_DATA.
	 * Two-list mate links are changed only for selected entries; source references are retained.
	 */
	private final void readSingles(){

		for(int i=0; !shutdown && i<producer1.size() && generated<maxReads; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try{
					list=depot.empty.take();
				}catch(InterruptedException e){
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
				final Read countedMate=(b==null ? a.mate : b);
				readsIn++;
				basesIn+=a.length();
				if(countedMate!=null){
					readsIn++;
					basesIn+=countedMate.length();
				}
				if(randy==null || randy.nextFloat()<samplerate){//Subsampled-IN reads are added; skipped reads still count toward readsIn/basesIn/generated above (they WERE read), only excluded from output.
					list.add(a);//Interleaved-pair convention: only 'a' enters the list; its mate 'b' rides along as a.mate (set just below). b is NEVER added to the list directly.
					if(b!=null){
						assert(a.numericID==b.numericID) : "\n"+a.numericID+", "+b.numericID+"\n"+a.toText(false)+"\n"+b.toText(false)+"\n";
						a.mate=b;
						b.mate=a;

						assert(a.pairnum()==0);
						b.setPairnum(1);
					}
					//Resolved #003: one-list callers retain attached mates, so their bases belong in the batch budget too.
					if(countedMate!=null){bases+=(countedMate.bases==null ? 0 : countedMate.length());}
					bases+=(a.bases==null ? 0 : a.length());
				}
				incrementGenerated(1);
			}

			if(verbose){System.err.println("E: Adding list("+list.size()+") to full "+depot.full.size()+"/"+depot.bufferCount);}
			//An empty selected batch can follow processing all remaining entries or reaching maxReads;
			//sampled-out input counts as processed, not delivered. Consumers treat empty lists as terminal.
			//Historical analysis attributed post-completion empty submissions to the outer shutdown guard
			//and consumer close. Preserve that rationale without claiming a universal iteration bound,
			//no-loss guarantee or lifecycle validation from this documentation review.
			depot.full.add(list);
		}
	}

	/** Sets the shutdown flag. Under coordinated lifecycle use, the subsequent interrupt
	 * guard is false (#001); this call alone does not promise to wake a waiting producer.
	 */
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

	/** Clears shutdown, replaces the depot and resets generated/input/progress counters.
	 * Finish prior use before calling. Retains list IDs, sampling RNG/state, source lists,
	 * inherited started flag and local error flag; does not launch another producer.
	 */
	@Override
	public synchronized void restart(){
		shutdown=false;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		generated=0;
		basesIn=0;
		readsIn=0;
		nextProgress=PROGRESS_INCR;
	}

	/** Requests shutdown and drains queued lists back to the available queue while waiting.
	 * Does not close or clear the source collections, latch errors, or reset counters.
	 * Caller must coordinate consumers; this implementation does not promise prompt return.
	 */
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
					try{
						if(verbose){System.out.println("close 4.");}
						list=depot.full.poll(100, TimeUnit.MILLISECONDS);
						if(verbose){System.out.println("close 5; list.size()="+depot.full.size()+", list="+(list==null ? "null" : list.size()+""));}
					}catch(InterruptedException e){
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
					try{
						if(verbose){System.out.println("close 12.");}
						threads[i].join();
						if(verbose){System.out.println("close 13.");}
					}catch(InterruptedException e){
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		if(verbose){System.out.println("close 14.");}

	}

	/** Reports paired when a second list exists, otherwise checks the first primary's mate.
	 * @return false for a null/empty primary without a second source; otherwise the inferred mode
	 */
	@Override
	public boolean paired(){//Paired if a second list was given; otherwise infer from producer1: empty->unpaired, else interleaved iff its first read has a mate attached.
		return producer2!=null ? true : (producer1==null || producer1.isEmpty() ? false : producer1.get(0).mate!=null);
	}

	/** Returns the shared verbose flag. */
	@Override
	public boolean verbose(){return verbose;}

	/** Advances the source-fragment count and prints at most one progress dot per call.
	 * @param amt Number of source entries consumed, including sampled-out entries
	 */
	private void incrementGenerated(long amt){
		generated+=amt;
		if(SHOW_PROGRESS && generated>=nextProgress){
			Data.sysout.print('.');
			nextProgress+=PROGRESS_INCR;
		}
	}

	/** Stores the sampling rate and configures the selection RNG.
	 * Rates at least one disable sampling; lower rates create a new RNG. Configure before
	 * starting; restart does not rewind this RNG. The nominal range is not checked here.
	 * @param rate Requested fraction in [0,1], applied per primary entry with its mate
	 * @param seed Seed forwarded to Shared.threadLocalRandom when sampling is enabled
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

	/** Returns observed primary and mate bases before sampling.
	 * A nonnull explicit secondary takes precedence over an already attached mate.
	 * @return Current counter snapshot, not a completion barrier
	 */
	@Override
	public long basesIn(){return basesIn;}
	/** Returns observed individual reads, including mates, before sampling.
	 * A nonnull explicit secondary takes precedence over an already attached mate.
	 * @return Current counter snapshot, not a completion barrier
	 */
	@Override
	public long readsIn(){return readsIn;}

	/** Returns the local flag, initialized false and not updated by this implementation.
	 * @return Cached flag; processing exceptions are not automatically recorded here
	 */
	@Override
	public boolean errorState(){return errorState;}
	/** Returns a new two-element array holding the retained primary/secondary list references. */
	@Override
	public Object[] producers(){return new Object[]{producer1, producer2};}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Shutdown request flag; lifecycle coordination belongs to the caller. */
	private boolean shutdown=false;
	/** Local status flag, initialized false and not subsequently assigned here. */
	private boolean errorState=false;
	/** Stored sampling fraction, retained across restart. */
	private float samplerate=1f;
	/** Selection RNG, null when sampling is disabled; retained across restart. */
	private shared.Random randy=null;
	/** Executing producer thread array recorded by run; not cleared by restart. */
	private Thread[] threads;

	/** Borrowed primary source list. */
	public final List<Read> producer1;
	/** Optional borrowed secondary list paired by primary index. */
	public final List<Read> producer2;
	/** Current buffer depot, replaced by restart. */
	private ConcurrentDepot<Read> depot;

	/** Input bases, including mate bases, counted before sampling. */
	private long basesIn=0;
	/** Individual input reads, including mates, counted before sampling. */
	private long readsIn=0;

	/** Primary-entry limit, with negative constructor input normalized to Long.MAX_VALUE. */
	private long maxReads;
	/** Primary entries consumed, including sampled-out entries; reset by restart. */
	private long generated=0;
	/** Next delivered batch ID, retained across restart. */
	private long listnum=0;
	/** Next progress-dot threshold, reset by restart. */
	private long nextProgress=PROGRESS_INCR;

	/** Shared verbose diagnostic setting. */
	public static boolean verbose=false;

	/** Unused retained marker; terminal paths currently allocate their own empty lists. */
	private static final ArrayList<Read> poison=new ArrayList<Read>(0);

}
