package stream;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.Shared;
import structures.ListNum;

/**
 * Distributed version of ConcurrentReadOutputStream for MPI-based parallel processing.
 * Provides master-slave architecture where the master node collects reads from slave nodes
 * and writes them to the underlying output stream.
 *
 * DORMANT / UNIMPLEMENTED SCAFFOLDING (the javadoc's old "incomplete implementation" note, made precise):
 * every MPI transport method (unicast/listen/listenFinishedSuccessfully/broadcastFinishedSuccessfully/
 * broadcastJoin/listenForJoin) is a stub that throws RuntimeException("TODO"). Unreachable by default:
 * the factory (ConcurrentReadOutputStream.getStream) only builds this when mpi==true AND crismpi=f, while
 * the default crismpi=true hits a deliberate assert(false) fence first (L~92). Even if forced, a slave
 * throws on its first add()->unicast(), and a master's ListenThreads throw on listen() (see #001). Treat
 * as never-exercised-in-production.
 *
 * @author Brian Bushnell
 * @date Jan 26, 2015
 */
public class ConcurrentReadOutputStreamD extends ConcurrentReadOutputStream{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public ConcurrentReadOutputStreamD(ConcurrentReadOutputStream cros_, boolean master_){
		super(cros_==null ? null : cros_.ff1, cros_==null ? null : cros_.ff2);
		dest=cros_;
		master=master_;
		rank=Shared.MPI_RANK;
		ranks=Shared.MPI_NUM_RANKS;
		assert(master==(cros_!=null));
	}
	
	@Override
	public synchronized void start(){
		if(started){
			System.err.println("Resetting output stream.");
			throw new RuntimeException();
		}
		
		started=true;
		if(master){
			terminatedCount.set(0);
			dest.start();
			startThreads();
		}
	}
	
	private void startThreads(){
		assert(master);
		for(int i=0; i<ranks; i++){
			if(i!=rank){
				ListenThread lt=new ListenThread(i);
				lt.start();
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Adds reads to the output stream using distributed processing.
	 * Master nodes write directly to destination stream, slaves send to master via MPI.
	 * @param list List of reads to output
	 * @param listnum Sequential list identifier for ordering
	 */
	@Override
	public synchronized void add(ArrayList<Read> list, long listnum){
		if(master){
			dest.add(list, listnum);
		}else{
			unicast(list, listnum, 0);
		}
	}

	@Override
	public void close(){
		if(master){
			int count=terminatedCount.incrementAndGet();
			//This loop blocks until every ListenThread increments terminatedCount; a ListenThread that dies WITHOUT incrementing (see #001) hangs it forever.
			while(count<ranks){
				synchronized(terminatedCount){
					count=terminatedCount.intValue();
					if(count<ranks){
						try {
							terminatedCount.wait(1000);
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
			}
			dest.close();
		}else{
			unicast(new ListNum<Read>(new ArrayList<Read>(1), Long.MAX_VALUE), 0);
		}
	}

	/**
	 * Waits for all processing to complete across all nodes.
	 * Master joins the destination stream and broadcasts completion to slaves.
	 * Slaves listen for master's join completion signal.
	 */
	@Override
	public void join(){
		if(master){
			dest.join();
			broadcastJoin(true);
		}else{
			boolean b=listenForJoin();
			assert(b);
		}
	}

	/** Resets the list ID counter and processing state on master node.
	 * Clears termination count and success flag for stream reuse. */
	@Override
	public synchronized void resetNextListID(){
		if(master){
			dest.resetNextListID();
			terminatedCount.set(0);
			finishedSuccessfully=false;
		}
	}
	
	/** Returns the output filename from the first file formatter.
	 * @return Output filename string */
	@Override
	public String fname(){
		return ff1.name();
	}
	
	/**
	 * Checks if the output stream is in an error state.
	 * Master nodes check both local and destination stream error states.
	 * @return true if any component is in error state
	 */
	@Override
	public boolean errorState(){
		if(master){
			return errorState || dest.errorState();
		}else{
			return errorState;
		}
	}

	/**
	 * Determines if all processing completed successfully across all nodes.
	 * Master queries destination stream and broadcasts result to slaves.
	 * Slaves listen for master's completion status.
	 * @return true if all processing finished without errors
	 */
	@Override
	public boolean finishedSuccessfully(){
		if(finishedSuccessfully){return true;}//Benign DCL: finishedSuccessfully is monotonic (false->true only), so an unsynchronized stale-false read merely re-enters the sync block below - it can never observe true->false, so it never returns a wrong result.

		synchronized(this){
			if(finishedSuccessfully){return true;}
			if(master){
				finishedSuccessfully=dest.finishedSuccessfully();
				broadcastFinishedSuccessfully(finishedSuccessfully);
			}else{
				finishedSuccessfully=listenFinishedSuccessfully();
			}
		}
		return finishedSuccessfully;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private void unicast(ArrayList<Read> list, long listnum, int i) {
		unicast(new ListNum<Read>(list, listnum), i);
	}
	
	protected void unicast(ListNum<Read> ln, int i) {//UNIMPLEMENTED: no-op MPI body, unconditionally throws RuntimeException("TODO"). Same for listen/listenFinishedSuccessfully/broadcastFinishedSuccessfully/broadcastJoin/listenForJoin.
		if(verbose){System.err.println("crosD "+(master?"master":"slave ")+":    Unicasting reads to "+i+".");}
		assert(!master);
		
		boolean success=false;
		while(!success){
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
	
	protected ListNum<Read> listen(int i){
		if(verbose){System.err.println("crosD "+(master?"master":"slave ")+":    Listening for reads from "+i+".");}
		assert(master);
		
		boolean success=false;
		while(!success){
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
	 * Slaves listen to master's finishedSuccessfully status.
	 * Currently incomplete with TODO placeholder for MPI implementation.
	 * @return Master's success status
	 * @throws RuntimeException Always thrown due to incomplete implementation
	 */
	protected boolean listenFinishedSuccessfully() {
		if(verbose){System.err.println("crosD "+(master?"master":"slave ")+":    listenFinishedSuccessfully.");}
		assert(!master);
		
		boolean success=false;
		while(!success){
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
	 * Master reports finishedSuccessfully status to slaves.
	 * Currently incomplete with TODO placeholder for MPI implementation.
	 * @param b Success status to broadcast
	 * @throws RuntimeException Always thrown due to incomplete implementation
	 */
	protected void broadcastFinishedSuccessfully(boolean b) {
		if(verbose){System.err.println("crosD "+(master?"master":"slave ")+":    broadcastFinishedSuccessfully.");}
		assert(master);
		
		boolean success=false;
		while(!success){
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
	
	protected void broadcastJoin(boolean b) {
		if(verbose){System.err.println("crosD "+(master?"master":"slave ")+":    broadcastJoin.");}
		assert(master);
		
		boolean success=false;
		while(!success){
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

	protected boolean listenForJoin() {
		if(verbose){System.err.println("crosD "+(master?"master":"slave ")+":    listenForJoin.");}
		assert(!master);
		
		boolean success=false;
		while(!success){
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
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ListenThread extends Thread{
		
		ListenThread(int sourceNum_){
			sourceNum=sourceNum_;
			assert(sourceNum_!=rank);
			assert(sourceNum>=0 && sourceNum<ranks);
		}
		
		@Override
		public void run(){
			assert(master);
			//TODO: Possible bug [stream/ConcurrentReadOutputStreamD#001] - hang on worker death: listen() currently throws RuntimeException("TODO") (unimplemented), and more generally ANY uncaught exception in this thread skips the terminatedCount.addAndGet(1) below, so master close() waits forever for terminatedCount to reach ranks. Latent LOW (dormant). Crash-loudly fix if ever implemented: guarantee the increment in a finally, or assertDie so a worker failure exits loud instead of hanging the master. See [[assertdie-idiom]].
			ListNum<Read> ln=listen(sourceNum);
			while(ln!=null && ln.id>=0){
				dest.add(ln.list, ln.id);
				ln=listen(sourceNum);
			}
			final int count=terminatedCount.addAndGet(1);
			if(count>=ranks){
				synchronized(terminatedCount){
					terminatedCount.notify();
				}
			}
		}
		
		final int sourceNum;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets the first read stream writer from destination (master only).
	 * @return First stream writer or null for slave nodes */
	@Override
	public ReadStreamWriter getRS1(){return master ? dest.getRS1() : null;}
	/** Gets the second read stream writer from destination (master only).
	 * @return Second stream writer or null for slave nodes */
	@Override
	public ReadStreamWriter getRS2(){return master ? dest.getRS2() : null;}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	protected final AtomicInteger terminatedCount=new AtomicInteger(0);
	protected final ConcurrentReadOutputStreamD thisPointer=this;
	
	/** Wrapped destination of reads. Null for slaves. */
	protected ConcurrentReadOutputStream dest;
	protected final boolean master;
	/** Total number of MPI ranks in the computation */
	protected final int rank, ranks;
	
}
