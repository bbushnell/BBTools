package stream;

import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import shared.Shared;
import structures.ListNum;

/**
 * Experimental depot wrapper for an unfinished MPI transport.
 * The master path publishes borrowed source lists and attempts rank-based delivery;
 * the receiving path depends on listening stubs. The concrete class has no working
 * MPI transport: remote unicast, general broadcast and listening paths throw TODO
 * exceptions. Unicast to the current rank returns without sending, and broadcast
 * can return through that path. Pairing/keepAll broadcasts are silent no-ops.
 *
 * The factory selects this class only with MPI enabled and USE_CRISMPI false.
 * USE_CRISMPI defaults true and selects a separate disabled factory branch.
 * The factory's null source for a receiving rank is dereferenced in this constructor
 * before it can reach the listening methods. Lifecycle and ownership plumbing here
 * is unfinished scaffolding, not a supported end-to-end distributed reader.
 * Local read/base counters are not updated by either transfer loop.
 *
 * @author Brian Bushnell
 * @date Oct 7, 2014
 */
public class ConcurrentReadInputStreamD extends ConcurrentReadInputStream{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Captures rank settings and allocates the depot without starting this wrapper.
	 * The master reads pairing from its source; status broadcasts do not transport it.
	 * @param cris_ Master source; a null receiving-rank source currently fails before role handling
	 * @param master_ Whether this instance owns the local source
	 * @param keepAll_ Whether the master requests delivery to every rank
	 * @throws NullPointerException If cris_ is null, while obtaining its filename
	 * @throws AssertionError If master_ disagrees with source presence under enabled assertions */
	public ConcurrentReadInputStreamD(ConcurrentReadInputStream cris_, boolean master_, boolean keepAll_){
		//TODO: Probable bug - STR258: the factory passes null for a receiving rank,
		//but this filename dereference occurs before the role assertion or listening stubs.
		super(cris_.fname());
		source=cris_;
		master=master_;
		rank=Shared.MPI_RANK;
		ranks=Shared.MPI_NUM_RANKS;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		assert(master==(cris_!=null));
		
		if(master){
			paired=source.paired();
			broadcastPaired(paired);
			keepAll=keepAll_;
			broadcastKeepall(keepAll);
		}else{//Concrete listening methods throw TODO; a null factory source fails earlier in super().
			paired=listenPaired();
			keepAll=listenKeepall();
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Takes a ready list and assigns a wrapper-local sequential ID.
	 * Shutdown is checked before each take; setting the flag alone does not wake a
	 * caller already waiting on the queue. Interrupted takes print a trace and retry.
	 * @return Borrowed list, possibly empty, or null when shutdown is observed before taking
	 */
	@Override
	public synchronized ListNum<Read> nextList(){
		ArrayList<Read> list=null;
		if(verbose){System.err.println("crisD:    **************** nextList() was called; shutdown="+shutdown+", depot.full="+depot.full.size());}
		while(list==null){
			if(shutdown){
				if(verbose){System.err.println("crisD:    **************** nextList() returning null; shutdown="+shutdown+", depot.full="+depot.full.size());}
				return null;
			}
			try{
				list=depot.full.take();
				assert(list!=null);
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
		
		if(verbose){System.err.println("crisD:    **************** nextList() returning list of size "+list.size()+"; shutdown="+shutdown+", depot.full="+depot.full.size());}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}
	
	/**
	 * Adds a fresh empty queue token without matching or recycling the caller's list.
	 * Normal returns replenish the available queue; terminal returns replenish the
	 * ready queue. Return each delivered list once to preserve queue accounting.
	 * @param listNumber Unused wrapper ID
	 * @param poison Whether to publish another empty terminal signal
	 */
	@Override
	public void returnList(long listNumber, boolean poison){
		if(poison){
			if(verbose){System.err.println("crisD:    A: Adding empty list to full.");}
			depot.full.add(new ArrayList<Read>(0));
		}else{
			if(verbose){System.err.println("crisD:    A: Adding empty list to empty.");}
			depot.empty.add(new ArrayList<Read>(BUF_LEN));//The token could have zero capacity; its entries are not used.
		}
	}
	
	/** Records the worker and dispatches to the role-specific transfer loop.
	 * Only normal return publishes terminal lists and clears the running marker;
	 * transport exceptions are not caught here. The master source must already be started. */
	@Override
	public void run(){
		synchronized(running){
			assert(!running[0]) : "This cris was started by multiple threads.";
			running[0]=true;
		}
		if(verbose){System.err.println("crisD:    cris started.");}
		threads=new Thread[]{Thread.currentThread()};
		
		if(master){readLists_master();}else{readLists_slave();}

		addPoison();
		
		//End thread

		while(!depot.empty.isEmpty() && !shutdown){
//			System.out.println("crisD:    Ending");
			if(verbose){System.err.println("crisD:    B: Adding empty lists to full.");}
			depot.full.add(depot.empty.poll());
		}
		if(verbose){System.err.println("crisD:    cris thread syncing before shutdown.");}
		
		synchronized(running){//Use the running monitor: close() holds this while waiting for the worker.
			//Taking this here would prevent termination while close() waits. Preserve the historical lock-order rationale.
			assert(running[0]);
			running[0]=false;
		}
		if(verbose){System.err.println("crisD:    cris thread terminated. Final depot size: "+depot.full.size()+", "+depot.empty.size());}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Publishes one empty terminal, then moves available queue tokens to ready output.
	 * Timed polls retry until a token arrives; an interrupted poll checks shutdown. */
	private final void addPoison(){
		//System.err.println("crisD:    Adding poison.");
		//Add poison pills
		if(verbose){System.err.println("crisD:    C: Adding poison to full.");}
		depot.full.add(new ArrayList<Read>());
		for(int i=1; i<depot.bufferCount; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try{
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				}catch(InterruptedException e){
//					System.err.println("crisD:    Do not be alarmed by the following error message:");
//					e.printStackTrace();
					if(shutdown){
						i=depot.bufferCount;
						break;
					}
				}
			}
			if(list!=null){
				if(verbose){System.err.println("crisD:    D: Adding list("+list.size()+") to full.");}
				depot.full.add(list);
			}
		}
		if(verbose){System.err.println("crisD:    Added poison.");}
	}
	
	/** Publishes selected borrowed source lists and attempts transport before returning source IDs.
	 * Empty source lists terminate the normal loop and are retained for final unicasts.
	 * No source-list copy is made and local input counters are not updated. */
	private final void readLists_master(){

		if(verbose){System.err.println("crisD:    Entered readLists_master().");}
		ListNum<Read> lnForUnicastShutdown=null;
		//TODO: Possible bug [stream/ConcurrentReadInputStreamD#002, STR260] - ln.list is read
		//without checking ln for null; a source returning null cannot be consumed here.
		//The receiving loop checks ln!=null. Deferred with the unfinished transport;
		//a current-rank unicast can return normally, so a second iteration is not ruled out.
		for(ListNum<Read> ln=source.nextList(); !shutdown && ln.list!=null; ln=source.nextList()){
			final ArrayList<Read> reads=ln.list;
			final int count=(reads==null ? 0 : reads.size());
			
			if(verbose){System.err.println("crisD:    Master fetched "+count+" reads.");}
			
			if(keepAll || count==0 || (ln.id%ranks)==rank){//Decide whether to process this list
				
				{
					ArrayList<Read> dummy=null;
					while(dummy==null && !shutdown){
						try{
							dummy=depot.empty.take();
						}catch(InterruptedException e){
							e.printStackTrace();
							if(shutdown){break;}
						}
					}
//					if(shutdown){break;}
				}
				
				try{
					depot.full.put(reads);
					if(verbose){System.err.println("crisD:    Master added reads to depot.");}
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			broadcast(ln);
			lnForUnicastShutdown=ln;
			if(verbose){System.err.println("crisD:    Master broadcasted.");}
			source.returnList(ln.id, count<1);
			if(verbose){System.err.println("crisD:    Master returned a list.");}
			if(count<1){break;}
		}
		if(!keepAll){//Shutdown all slaves if unicasting
			for(int i=1; i<ranks; i++){
				unicast(lnForUnicastShutdown, i);
			}
		}
		if(verbose){System.err.println("crisD:    Finished readLists_master().");}
	}
	
	/** Intended receiving loop, selecting lists by rank and publishing their references.
	 * The concrete listen() stub throws; no local input counters are updated. */
	private final void readLists_slave(){
		
		if(verbose){System.err.println("crisD:    Entered readLists_slave().");}
		for(ListNum<Read> ln=listen(); !shutdown && ln!=null; ln=listen()){
			
			final ArrayList<Read> reads=ln.list;
			final int count=(reads==null ? 0 : reads.size());
			
			if(verbose){System.err.println("crisD:    Slave fetched "+count+" reads.");}

			if(keepAll || count==0 || (ln.id%ranks)==rank){//Decide whether to process this list
				{
					ArrayList<Read> dummy=null;
					while(dummy==null && !shutdown){
						try{
							dummy=depot.empty.take();
						}catch(InterruptedException e){
							e.printStackTrace();
							if(shutdown){break;}
						}
					}
//					if(shutdown){break;}
				}
				
				try{
					depot.full.put(reads);
					if(verbose){System.err.println("crisD:    Slave added reads to depot.");}
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			if(count<1){break;}
		}
		if(verbose){System.err.println("crisD:    Finished readLists_slave().");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Concurrency Methods     ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Routes a nonempty list to its rank when keepAll is false; otherwise reaches a TODO stub.
	 * Current-rank unicast returns without sending; all other concrete paths throw.
	 * @param ln Nonnull source list whose ID selects its destination rank */
	protected void broadcast(ListNum<Read> ln){
		if(!keepAll && ln.size()>0){//Decide how to send this list
			final int toRank=(int)(ln.id%ranks);
			unicast(ln, toRank);
			return;
		}
		
		if(verbose){System.err.println("crisD "+(master ? "master" : "slave ")+":    Broadcasting reads.");}
		
		boolean success=false;
		while(!success && !shutdown){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
		
	/** Returns for the current rank; remote delivery is unimplemented and throws.
	 * @param ln Intended payload, unused by the concrete stub
	 * @param toRank Destination rank; no range validation is performed */
	protected void unicast(ListNum<Read> ln, final int toRank){
		if(toRank==rank){return;}
		if(verbose){System.err.println("crisD "+(master ? "master" : "slave ")+":    Unicasting reads to "+toRank+".");}

		boolean success=false;
		while(!success && !shutdown){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	/** Pairing-status broadcast placeholder; performs no transport and ignores the value.
	 * @param b Intended pairing status */
	protected void broadcastPaired(boolean b){
		if(verbose){System.err.println("crisD "+(master ? "master" : "slave ")+":    Broadcasting pairing status.");}
		boolean success=false;
		while(!success && !shutdown){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
//		throw new RuntimeException("TODO");
	}
		
	/** Distribution-mode broadcast placeholder; performs no transport and ignores the value.
	 * @param b Intended keepAll status */
	protected void broadcastKeepall(boolean b){
		if(verbose){System.err.println("crisD "+(master ? "master" : "slave ")+":    Broadcasting keepAll status.");}
		boolean success=false;
		while(!success && !shutdown){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
//		throw new RuntimeException("TODO");
	}

	/** Intended list receiver; the concrete implementation always throws a TODO exception.
	 * @return No value from the current implementation */
	protected ListNum<Read> listen(){
		if(verbose){System.err.println("crisD "+(master ? "master" : "slave ")+":    Listening to "+0+" for reads.");}
		boolean success=false;
		while(!success && !shutdown){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	/** Intended pairing receiver; the concrete implementation always throws a TODO exception.
	 * @return No value from the current implementation */
	protected boolean listenPaired(){
		if(verbose){System.err.println("crisD "+(master ? "master" : "slave ")+":    Listening to "+0+" for pairing status.");}
		boolean success=false;
		while(!success && !shutdown){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	/** Intended distribution-mode receiver; always throws a TODO exception.
	 * @return No value from the current implementation */
	protected boolean listenKeepall(){
		if(verbose){System.err.println("crisD "+(master ? "master" : "slave ")+":    Listening to "+0+" for keepAll status.");}
		boolean success=false;
		while(!success && !shutdown){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Termination          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Sets the shutdown flag; the existing propagation branch below is unreachable
	 * in coordinated use because it tests the flag just set to true.
	 * Does not currently interrupt the worker or propagate to the master source. */
	@Override
	public void shutdown(){
		if(verbose){System.out.println("crisD:    Called shutdown.");}
		
		shutdown=true;
		//TODO: Possible bug [stream/ConcurrentReadInputStreamD#001, STR259] - the branch
		//tests the new shutdown state rather than its prior value, skipping source shutdown
		//and worker interruption in coordinated use. Deferred with the unfinished transport.
		if(!shutdown){//The flag was just set above; see STR259.
			
			if(master){source.shutdown();}
			for(Thread t : threads){
				if(t!=null && t.isAlive()){t.interrupt();}
			}
		}
	}
	
	/** Replaces the depot and resets shutdown, local counters and list IDs; restarts the master source.
	 * Requires prior use to be quiescent. Retains local error state, running marker,
	 * worker references and captured rank/pairing settings; does not start a worker. */
	@Override
	public synchronized void restart(){
		shutdown=false;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		basesIn=0;
		readsIn=0;
		listnum=0;//Added Oct 9, 2014
		if(master){source.restart();}
	}
	
	/** Sets shutdown and closes the master source, then waits for recorded workers.
	 * While the first worker is alive, polls ready lists, clears them and replenishes
	 * available tokens. This is not a guarantee of completion for the unfinished transport. */
	@Override
	public synchronized void close(){
		shutdown();
		
		if(master){
			source.close();
		}else{
			
		}
		
		if(threads!=null && threads[0]!=null && threads[0].isAlive()){
			
			while(threads[0].isAlive()){
//				System.out.println("crisD:    B");
				ArrayList<Read> list=null;
				for(int i=0; i<1000 && list==null && threads[0].isAlive(); i++){
					try{
						list=depot.full.poll(200, TimeUnit.MILLISECONDS);
					}catch(InterruptedException e){
						System.err.println("crisD:    Do not be alarmed by the following error message:");
						e.printStackTrace();
						break;
					}
				}
				
				if(list!=null){
					list.clear();
					depot.empty.add(list);
				}
				
//				System.out.println("crisD:    isAlive? "+threads[0].isAlive());
			}
			
		}
		
		if(threads!=null){
			for(int i=1; i<threads.length; i++){
				while(threads[i]!=null && threads[i].isAlive()){
					try{
						threads[i].join();
					}catch(InterruptedException e){
						e.printStackTrace();
					}
				}
			}
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns whether input contains paired-end reads. */
	@Override
	public boolean paired(){return paired;}
	
	/** Returns the current verbose logging setting. */
	@Override
	public boolean verbose(){return verbose;}
	
	/**
	 * Forwards sampling configuration to the master source; receiving ranks do nothing.
	 * This wrapper does not validate the rate or manage a random generator.
	 * @param rate Requested rate interpreted by the source
	 * @param seed Requested seed interpreted by the source
	 */
	@Override
	public void setSampleRate(float rate, long seed){
		if(master){source.setSampleRate(rate, seed);}
	}
	
	/** @return Local base counter, which current transfer paths never increment */
	@Override
	public long basesIn(){return basesIn;}
	/** @return Local read counter, which current transfer paths never increment */
	@Override
	public long readsIn(){return readsIn;}
	
	/**
	 * Reads stored error flags without closing or waiting for completion.
	 * This class does not set its local flag; masters also consult the source.
	 * @return Current local flag, combined with the source flag for masters
	 */
	@Override
	public boolean errorState(){
		if(master){return errorState|source.errorState();}
		return errorState;
	}
	
	/** Returns array of producer objects for masters; slaves return null since they don't have direct producers.
	 * @return Producer array for masters, null for slaves */
	@Override
	public Object[] producers(){
		if(master){return source.producers();}
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Retained master source; the intended null receiving-rank source currently fails construction. */
	private ConcurrentReadInputStream source;
	/** Captured role controlling source access and transfer-loop selection. */
	private final boolean master;
	/** Captured distribution policy; receiving-rank initialization depends on a stub. */
	protected final boolean keepAll;
	/** Captured current rank and rank count; this class does not validate their range. */
	protected final int rank, ranks;
	
	/** Local error flag, never set true by this concrete implementation. */
	private boolean errorState=false;
	
	/** Worker activity marker and its dedicated monitor, separate from the close monitor. */
	private boolean[] running=new boolean[]{false};
	
	/** Plain request flag; lifecycle operations require coordinated access. */
	private boolean shutdown=false;

	/** Ready borrowed lists and available tokens; replaced on restart. */
	private ConcurrentDepot<Read> depot;

	/** Worker references installed by run(), inspected by lifecycle methods. */
	private Thread[] threads;
	
	/** Local base count, reset but never incremented here. */
	private long basesIn=0;
	/** Local individual-read count, reset but never incremented here. */
	private long readsIn=0;
	
	/** Next wrapper-local list ID, reset by restart. */
	private long listnum=0;
	
	/** Pairing from the master source; the concrete receiving-rank listener throws before assignment. */
	private final boolean paired;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Global diagnostic flag shared by instances. */
	public static boolean verbose=false;
	
}
