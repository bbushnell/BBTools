package stream;

import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import shared.Shared;
import structures.ListNum;

/**
 * This class is designed for distributed environments.
 * The 'master' reads from the filesystem, creates reads, and broadcasts them.
 * The 'slaves' listen for broadcasts.
 * @author Brian Bushnell
 * @date Oct 7, 2014
 *
 */
public class ConcurrentReadInputStreamD extends ConcurrentReadInputStream {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a distributed read input stream.
	 * Master nodes wrap an existing stream source, slaves listen for broadcasts.
	 * Initializes MPI parameters and concurrent depot for buffering.
	 *
	 * @param cris_ Source stream for master nodes, null for slaves
	 * @param master_ True if this is the master node that reads from filesystem
	 * @param keepAll_ True to process all reads, false for distributed processing
	 */
	public ConcurrentReadInputStreamD(ConcurrentReadInputStream cris_, boolean master_, boolean keepAll_){
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
		}else{
			paired=listenPaired();
			keepAll=listenKeepall();
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public synchronized ListNum<Read> nextList() {
		ArrayList<Read> list=null;
		if(verbose){System.err.println("crisD:    **************** nextList() was called; shutdown="+shutdown+", depot.full="+depot.full.size());}
		while(list==null){
			if(shutdown){
				if(verbose){System.err.println("crisD:    **************** nextList() returning null; shutdown="+shutdown+", depot.full="+depot.full.size());}
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
		
		if(verbose){System.err.println("crisD:    **************** nextList() returning list of size "+list.size()+"; shutdown="+shutdown+", depot.full="+depot.full.size());}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}
	
	@Override
	public void returnList(long listNumber, boolean poison){
		if(poison){
			if(verbose){System.err.println("crisD:    A: Adding empty list to full.");}
			depot.full.add(new ArrayList<Read>(0));
		}else{
			if(verbose){System.err.println("crisD:    A: Adding empty list to empty.");}
			depot.empty.add(new ArrayList<Read>(BUF_LEN)); //Technically this could be a length-0 list since it is never used.
		}
	}
	
	@Override
	public void run() {
		synchronized(running){
			assert(!running[0]) : "This cris was started by multiple threads.";
			running[0]=true;
		}
		if(verbose){System.err.println("crisD:    cris started.");}
		threads=new Thread[] {Thread.currentThread()};
		
		if(master){
			readLists_master();
		}else{
			readLists_slave();
		}

		addPoison();
		
		//End thread

		while(!depot.empty.isEmpty() && !shutdown){
//			System.out.println("crisD:    Ending");
			if(verbose){System.err.println("crisD:    B: Adding empty lists to full.");}
			depot.full.add(depot.empty.poll());
		}
		if(verbose){System.err.println("crisD:    cris thread syncing before shutdown.");}
		
		synchronized(running){//TODO Note: for some reason syncing on 'this' instead of 'running' causes a hang.  Something else must be syncing improperly on this.
			assert(running[0]);
			running[0]=false;
		}
		if(verbose){System.err.println("crisD:    cris thread terminated. Final depot size: "+depot.full.size()+", "+depot.empty.size());}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds poison pills to signal end of processing.
	 * Creates empty lists and adds them to full depot to unblock consumers.
	 * Handles interruption gracefully during shutdown.
	 */
	private final void addPoison(){
		//System.err.println("crisD:    Adding poison.");
		//Add poison pills
		if(verbose){System.err.println("crisD:    C: Adding poison to full.");}
		depot.full.add(new ArrayList<Read>());
		for(int i=1; i<depot.bufferCount; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
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
	
	/**
	 * Master node processing loop that reads from source stream.
	 * Distributes reads based on rank assignment and broadcasts to slaves.
	 * Handles unicast shutdown signaling when not keeping all reads.
	 */
	private final void readLists_master(){

		if(verbose){System.err.println("crisD:    Entered readLists_master().");}
		ListNum<Read> lnForUnicastShutdown=null;
		for(ListNum<Read> ln=source.nextList(); !shutdown && ln.list!=null; ln=source.nextList()){
			final ArrayList<Read> reads=ln.list;
			final int count=(reads==null ? 0 : reads.size());
			
			if(verbose){System.err.println("crisD:    Master fetched "+count+" reads.");}
			
			if(keepAll || count==0 || (ln.id%ranks)==rank){//Decide whether to process this list
				
				{
					ArrayList<Read> dummy=null;
					while(dummy==null && !shutdown){
						try {
							dummy=depot.empty.take();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
							if(shutdown){break;}
						}
					}
//					if(shutdown){break;}
				}
				
				try {
					depot.full.put(reads);
					if(verbose){System.err.println("crisD:    Master added reads to depot.");}
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
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
	
	/**
	 * Slave node processing loop that listens for broadcasted reads.
	 * Processes reads assigned to this rank and adds them to local depot.
	 * Continues until master signals completion.
	 */
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
						try {
							dummy=depot.empty.take();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
							if(shutdown){break;}
						}
					}
//					if(shutdown){break;}
				}
				
				try {
					depot.full.put(reads);
					if(verbose){System.err.println("crisD:    Slave added reads to depot.");}
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
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
	
	/**
	 * Broadcasts read list to all slave processes.
	 * Uses unicast for targeted distribution when not keeping all reads.
	 * Currently incomplete - throws RuntimeException as placeholder for MPI implementation.
	 * @param ln ListNum containing reads to broadcast
	 */
	protected void broadcast(ListNum<Read> ln){
		if(!keepAll && ln.size()>0){//Decide how to send this list
			final int toRank=(int)(ln.id%ranks);
			unicast(ln, toRank);
			return;
		}
		
		if(verbose){System.err.println("crisD "+(master?"master":"slave ")+":    Broadcasting reads.");}
		
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
		
	/**
	 * Sends read list to specific slave process rank.
	 * Skips sending if target rank matches current rank.
	 * Currently incomplete - throws RuntimeException as placeholder for MPI implementation.
	 *
	 * @param ln ListNum containing reads to send
	 * @param toRank Target MPI rank to receive the reads
	 */
	protected void unicast(ListNum<Read> ln, final int toRank){
		if(toRank==rank){return;}
		if(verbose){System.err.println("crisD "+(master?"master":"slave ")+":    Unicasting reads to "+toRank+".");}

		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	/**
	 * Broadcasts paired-end read status to all slaves.
	 * Ensures all processes know whether input contains paired reads.
	 * Currently incomplete - placeholder for MPI implementation.
	 * @param b True if reads are paired-end
	 */
	protected void broadcastPaired(boolean b){
		if(verbose){System.err.println("crisD "+(master?"master":"slave ")+":    Broadcasting pairing status.");}
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
//		throw new RuntimeException("TODO");
	}
		
	/**
	 * Broadcasts keepAll processing mode to all slaves.
	 * Determines whether all processes handle all reads or distribute by rank.
	 * Currently incomplete - placeholder for MPI implementation.
	 * @param b True to process all reads on all nodes
	 */
	protected void broadcastKeepall(boolean b){
		if(verbose){System.err.println("crisD "+(master?"master":"slave ")+":    Broadcasting keepAll status.");}
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
//		throw new RuntimeException("TODO");
	}

	/**
	 * Listens for broadcasted read lists from master process.
	 * Blocks until data is received or shutdown occurs.
	 * Currently incomplete - throws RuntimeException as placeholder for MPI implementation.
	 * @return ListNum containing received reads, null on shutdown
	 */
	protected ListNum<Read> listen(){
		if(verbose){System.err.println("crisD "+(master?"master":"slave ")+":    Listening to "+0+" for reads.");}
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	/**
	 * Listens for paired-end status broadcast from master.
	 * Receives initial configuration about whether reads are paired.
	 * Currently incomplete - throws RuntimeException as placeholder for MPI implementation.
	 * @return True if reads are paired-end
	 */
	protected boolean listenPaired(){
		if(verbose){System.err.println("crisD "+(master?"master":"slave ")+":    Listening to "+0+" for pairing status.");}
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	/**
	 * Listens for keepAll mode broadcast from master.
	 * Receives configuration about read distribution strategy.
	 * Currently incomplete - throws RuntimeException as placeholder for MPI implementation.
	 * @return True if all nodes should process all reads
	 */
	protected boolean listenKeepall(){
		if(verbose){System.err.println("crisD "+(master?"master":"slave ")+":    Listening to "+0+" for keepAll status.");}
		boolean success=false;
		while(!success && !shutdown){
			try {
				//Do some MPI stuff
				success=true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Termination          ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void shutdown(){
		if(verbose){System.out.println("crisD:    Called shutdown.");}
		
		shutdown=true;
		if(!shutdown){ //Possible bug: This condition will never be true since shutdown was just set to true above
			
			if(master){
				source.shutdown();
			}
			for(Thread t : threads){
				if(t!=null && t.isAlive()){
					t.interrupt();
				}
			}
		}
	}
	
	@Override
	public synchronized void restart(){
		shutdown=false;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		basesIn=0;
		readsIn=0;
		listnum=0; //Added Oct 9, 2014
		if(master){
			source.restart();
		}
	}
	
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
					try {
						list=depot.full.poll(200, TimeUnit.MILLISECONDS);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
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
					try {
						threads[i].join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public boolean paired() {return paired;}
	
	@Override
	public boolean verbose(){return verbose;}
	
	@Override
	public void setSampleRate(float rate, long seed){
		if(master){source.setSampleRate(rate, seed);}
	}
	
	@Override
	public long basesIn(){return basesIn;}
	@Override
	public long readsIn(){return readsIn;}
	
	@Override
	public boolean errorState(){
		if(master){return errorState|source.errorState();}
		return errorState;
	}
	
	@Override
	public Object[] producers(){
		if(master){return source.producers();}
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Wrapped source of reads.  Null for slaves. */
	private ConcurrentReadInputStream source;
	/** True if this is the master process that reads from filesystem */
	private final boolean master;
	/** True to process all reads on all nodes, false for distributed processing */
	protected final boolean keepAll;
	protected final int rank, ranks;
	
	/** Flag indicating whether stream has encountered errors */
	private boolean errorState=false;
	
	/** Array tracking whether processing thread is running */
	private boolean[] running=new boolean[] {false};
	
	/** Flag to signal shutdown of processing threads */
	private boolean shutdown=false;

	/** Concurrent buffer depot for thread-safe read list storage */
	private ConcurrentDepot<Read> depot;

	/** Array of worker threads for concurrent processing */
	private Thread[] threads;
	
	/** Counter for total bases processed by this stream */
	private long basesIn=0;
	/** Counter for total reads processed by this stream */
	private long readsIn=0;
	
	/** Sequential counter for read list numbering */
	private long listnum=0;
	
	/** This should be set in the first broadcast */
	private final boolean paired;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Global flag controlling verbose debug output for all instances */
	public static boolean verbose=false;
	
}
