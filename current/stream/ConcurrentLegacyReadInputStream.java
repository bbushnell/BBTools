package stream;

import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import structures.ListNum;

/**
 * Wraps one ReadInputStream and transfers borrowed Read references through a depot.
 * Empty output lists are terminal signals; normal list returns replenish empty
 * buffers rather than returning the same list object. The factory uses this wrapper
 * for bread and sequential inputs; direct mapper, variant and pacbio callers remain.
 * Input counters include mates and are updated before sampling, while generated
 * counts root entries. Sampling and lifecycle changes require coordinated use.
 * This legacy implementation is intended for replacement by the generic wrapper;
 * close only closes its source and is not a producer-thread completion barrier.
 *
 * @author Brian Bushnell
 */
public class ConcurrentLegacyReadInputStream extends ConcurrentReadInputStream{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Captures one source and allocates the depot without starting the wrapper thread.
	 * @param source Nonnull source supporting batch reads
	 * @param maxReadsToGenerate Requested root-entry limit before sampling; negative means unlimited
	 * @throws AssertionError If a zero limit is supplied while assertions are enabled */
	public ConcurrentLegacyReadInputStream(ReadInputStream source, long maxReadsToGenerate){
		super(source.fname());
		producer=source;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		maxReads=maxReadsToGenerate>=0 ? maxReadsToGenerate : Long.MAX_VALUE;
		if(maxReads==0){
			System.err.println("Warning - created a read stream for 0 reads.");
			assert(false);
		}
//		if(maxReads<Long.MAX_VALUE){System.err.println("maxReads="+maxReads);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Submission          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Takes the next ready list and assigns a monotonically increasing wrapper ID.
	 * Interrupted waits print a trace and retry unless shutdown is observed.
	 * @return Borrowed list, including an empty terminal list; null only after an
	 * interrupted wait observes shutdown. Return each delivered list exactly once. */
	@Override
	public synchronized ListNum<Read> nextList(){
		ArrayList<Read> list=null;
		while(list==null){
			try{
				list=depot.full.take();
			}catch(InterruptedException e){
				//Report the interruption, then retry unless shutdown was requested.
				e.printStackTrace();
				if(shutdown){return null;}
			}
		}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}
	
	/** Replenishes one depot slot or publishes another terminal signal.
	 * This method allocates a fresh empty list; it does not match or recycle the caller's list.
	 * @param listNumber Unused wrapper ID
	 * @param poison Whether to add an empty terminal to full instead of an available buffer to empty */
	@Override
	public void returnList(long listNumber, boolean poison){
		if(poison){
			if(verbose){System.err.println("cris_:    A: Adding empty list to full.");}
			depot.full.add(new ArrayList<Read>(0));
		}else{
			if(verbose){System.err.println("cris_:    A: Adding empty list to empty.");}
			depot.empty.add(new ArrayList<Read>(BUF_LEN));
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Producer           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Records the current worker, fills output lists and publishes terminal lists on normal return.
	 * Does not separately start the source or catch every failure from readLists(). */
	@Override
	public void run(){
//		producer.start();
		threads=new Thread[]{Thread.currentThread()};
		
		readLists();
		
		addPoison();
		
		//End thread
		
		while(!depot.empty.isEmpty()){
			depot.full.add(depot.empty.poll());
		}
//		System.err.println(depot.full.size()+", "+depot.empty.size());
	}
	
	/** Publishes an empty terminal, then transfers available depot buffers to ready output.
	 * Uses timed polls until buffers arrive; only an interrupted poll checks shutdown. */
	private final void addPoison(){
		//System.err.println("Adding poison.");
		//Add poison pills
		depot.full.add(new ArrayList<Read>());
		for(int i=1; i<depot.bufferCount; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try{
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				}catch(InterruptedException e){
					//An interrupt ends this wait only when shutdown was requested.
//					System.err.println("Do not be alarmed by the following error message:");
//					e.printStackTrace();
					if(shutdown){
						i=depot.bufferCount;
						break;
					}
				}
			}
			if(list!=null){depot.full.add(list);}
		}
		//System.err.println("Added poison.");
	}
	
	/**
	 * Consumes source batches through bulk-copy and individual-entry paths.
	 * Generated counts entries examined before sampling; input counts include mates.
	 * The per-entry path bounds retained list size and tests retained base totals;
	 * the bulk path tests source-list capacity and generation limits before copying.
	 */
	private final void readLists(){
		
		ArrayList<Read> buffer=null;
		ArrayList<Read> list=null;
		int next=0;
		//TODO: Probable bug - STR255: a retained source tail bypasses the outer maxReads check.
		//At the limit the inner loop cannot consume it, so empty lists can be emitted repeatedly.
		while(buffer!=null || (!shutdown && producer.hasMore() && generated<maxReads)){
			while(list==null){
				//System.err.println("Fetching a list: generated="+generated+"/"+maxReads);
				try{
					list=depot.empty.take();
				}catch(InterruptedException e){
					//Report the interruption, then retry unless shutdown was requested.
					e.printStackTrace();
					if(shutdown){break;}
				}
				//System.err.println("Fetched");
			}
			if(shutdown || list==null){
				//System.err.println("Shutdown triggered; breaking.");
				break;
			}
			
			long bases=0;
			while(list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
				if(buffer==null || next>=buffer.size()){
					buffer=producer.nextList();
					next=0;
				}
				if(buffer==null){break;}
				assert(buffer.size()<=BUF_LEN);//Although this is not really necessary.
				
				if(next==0 && buffer.size()<=(BUF_LEN-list.size()) && (buffer.size()+generated)<maxReads && randy==null){
					//Historical bulk-path rationale: copy a fitting unsampled source batch without per-entry insertion.
					//STR254: addAll starts at index zero, so only an untouched source batch may use it.
					//A retained tail must resume at next through the per-entry path to avoid repeating its prefix.
					list.addAll(buffer);
					for(Read a : buffer){
						readsIn++;
						basesIn+=a.length();
						bases+=a.length();
						if(a.mate!=null){
							readsIn++;
							basesIn+=a.mateLength();
							bases+=a.mateLength();
						}
					}
					generated+=buffer.size();
					next=0;
					buffer=null;
				}else{
					while(next<buffer.size() && list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
						Read r=buffer.get(next);
						readsIn++;
						basesIn+=r.length();
						if(r.mate!=null){
							readsIn++;
							basesIn+=r.mateLength();
						}
						if(randy==null || randy.nextFloat()<samplerate){
							list.add(r);
							bases+=r.length();
							bases+=(r.mate==null || r.mate.bases==null ? 0 : r.mateLength());
						}
						generated++;
//						if(generated>1 && (generated%1000000)==0){System.err.println("Generated read #"+generated);}
						next++;
					}
					
					if(next>=buffer.size()){
						buffer=null;
						next=0;
					}
				}
			}
			//System.err.println("Adding list to full depot.");
			depot.full.add(list);
			//System.err.println("Added.");
			list=null;
		}

	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Lifecycle           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Requests shutdown and interrupts the recorded worker when it is present and alive.
	 * Does not close the source or join the worker. */
	@Override
	public void shutdown(){
		shutdown=true;
		//#001-fix [stream/ConcurrentLegacyReadInputStream#001]: run() installs threads.
		//Guard the pre-run null reference; the shutdown flag also gates initial source reads.
		if(threads!=null && threads[0]!=null && threads[0].isAlive()){
			threads[0].interrupt();
		}
	}
	
	/** Restarts the source and replaces the depot after prior use has become quiescent.
	 * Resets shutdown and input/generated counters, but retains list IDs, sampling,
	 * local error state and worker references. Does not start a new wrapper thread. */
	@Override
	public synchronized void restart(){
		shutdown=false;
		producer.restart();
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		generated=0;
		basesIn=0;
		readsIn=0;
	}
	
	/** Closes only the underlying source, without joining or interrupting the wrapper worker.
	 * The source's boolean close result is not retained; errorState() consults its stored flag. */
	@Override
	public synchronized void close(){
//		System.err.println("Closing cris: "+maxReads+", "+generated);
//		if(threads!=null){
//			for(int i=0; i<threads.length; i++){
//				if(threads[i]!=null){System.err.println(i+": "+threads[i].isAlive());}
//			}
//		}
		//TODO: Possible bug [stream/ConcurrentLegacyReadInputStream#002, STR256] - closing
		//the source does not release a wrapper worker waiting for a depot buffer after early
		//consumer abandonment. This method neither calls shutdown() nor replenishes buffers.
		//A worker started from a non-daemon context may then keep the process alive.
		//Retain the structural lifecycle follow-up separately from this documentation pass;
		//normal completion also needs the retained-tail limit concern in STR255 considered.
		producer.close();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Configuration/Status   ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns whether the underlying stream contains paired-end reads.
	 * @return true if stream contains paired reads */
	@Override
	public boolean paired(){return producer.paired();}
	
	/** @return Current global legacy-wrapper diagnostic flag */
	@Override
	public boolean verbose(){return verbose;}
	
	/** Configures sampling before processing; this method does not validate or clamp the rate.
	 * Input/generated counters are updated before a sampled entry is retained or rejected.
	 * @param rate Rates at least one disable sampling; lower values use a random comparison
	 * @param seed Nonnegative seed passed to the random factory; negative selects its default factory */
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		if(rate>=1f){
			randy=null;
		}else if(seed>-1){
			randy=shared.Shared.random(seed);
		}else{
			randy=shared.Shared.random();
		}
	}
	
	/** @return Observed bases examined before sampling, including mates; not a completion wait */
	@Override
	public long basesIn(){return basesIn;}
	/** @return Observed individual reads examined before sampling, including mates */
	@Override
	public long readsIn(){return readsIn;}
	
	/** Combines currently observed stored error flags without closing or waiting.
	 * @return true if the local flag or the source's errorState() reports an error */
	@Override
	public boolean errorState(){return errorState || (producer!=null && producer.errorState());}
	/** @return New one-element array containing the retained source, without copying it */
	@Override
	public Object[] producers(){return new Object[]{producer};}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Plain shutdown request flag; lifecycle callers must coordinate access. */
	private boolean shutdown=false;
	/** Retained local error flag; this class does not currently set it true. */
	private boolean errorState=false;
	/** Sampling comparison threshold, configured before processing. */
	private float samplerate=1f;
	/** Null when sampling is disabled; otherwise the configured random source. */
	private shared.Random randy=null;
	/** Worker reference installed by run(), consulted by shutdown(). */
	private Thread[] threads;

	/** Retained source; callers must coordinate any direct access with this wrapper. */
	public final ReadInputStream producer;
	/** Available and ready lists for the current run; replaced by restart(). */
	private ConcurrentDepot<Read> depot;
	
	/** Global diagnostic flag shared by legacy-wrapper instances. */
	public static boolean verbose=false;
	
	/** Pair-aware examined-base count, reset by restart(). */
	private long basesIn=0;
	/** Individual examined-read count, reset by restart(). */
	private long readsIn=0;
	
	/** Requested root-entry limit, with negative constructor values mapped to Long.MAX_VALUE. */
	private long maxReads;
	/** Root entries counted before sampling, reset by restart(). */
	private long generated=0;
	/** Next delivered-list ID; not reset by restart(). */
	private long listnum=0;
	

}
