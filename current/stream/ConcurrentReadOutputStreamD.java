package stream;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.Shared;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Jan 26, 2015
 *
 */
public class ConcurrentReadOutputStreamD extends ConcurrentReadOutputStream{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a distributed read output stream wrapper.
	 * Master nodes wrap an existing output stream, slave nodes have null destination.
	 * @param cros_ The underlying output stream for master nodes, null for slaves
	 * @param master_ true if this is the master node that will write to files
	 */
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
	
	/** Launches listener threads on master node to receive data from each slave rank.
	 * Creates one ListenThread per slave rank to handle incoming read data. */
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

	@Override
	public synchronized void resetNextListID(){
		if(master){
			dest.resetNextListID();
			terminatedCount.set(0);
			finishedSuccessfully=false;
		}
	}
	
	@Override
	public String fname(){
		return ff1.name();
	}
	
	@Override
	public boolean errorState(){
		if(master){
			return errorState || dest.errorState();
		}else{
			return errorState;
		}
	}

	@Override
	public boolean finishedSuccessfully(){
		if(finishedSuccessfully){return true;}
		
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
	
	/**
	 * Sends read list to specific rank via unicast MPI communication.
	 * Convenience wrapper that creates ListNum object.
	 *
	 * @param list Read data to send
	 * @param listnum List identifier
	 * @param i Target rank number
	 */
	private void unicast(ArrayList<Read> list, long listnum, int i) {
		unicast(new ListNum<Read>(list, listnum), i);
	}
	
	/**
	 * Sends packaged read data to specified rank using MPI unicast.
	 * Currently incomplete with TODO placeholder for MPI implementation.
	 *
	 * @param ln Packaged read list with identifier
	 * @param i Target rank to send data to
	 * @throws RuntimeException Always thrown due to incomplete implementation
	 */
	protected void unicast(ListNum<Read> ln, int i) {
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
	
	/**
	 * Listens for incoming read data from specified rank using MPI.
	 * Master nodes use this to receive data from slave ranks.
	 * Currently incomplete with TODO placeholder for MPI implementation.
	 *
	 * @param i Source rank to listen for data from
	 * @return Received read data package
	 * @throws RuntimeException Always thrown due to incomplete implementation
	 */
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
	
	/** Master tells slaves that 'join' was successful. */
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

	/** Slave listens to see that master 'join' was successful. */
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
	
	/**
	 * Worker thread that listens for read data from a specific slave rank.
	 * Master node creates one thread per slave to handle incoming data streams.
	 * Processes received reads until termination signal, then notifies completion.
	 */
	private class ListenThread extends Thread{
		
		/** Creates listener thread for specified source rank.
		 * @param sourceNum_ The rank number to listen for data from */
		ListenThread(int sourceNum_){
			sourceNum=sourceNum_;
			assert(sourceNum_!=rank);
			assert(sourceNum>=0 && sourceNum<ranks);
		}
		
		@Override
		public void run(){
			assert(master);
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
		
		/** The rank number this thread listens to for incoming data */
		final int sourceNum;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public ReadStreamWriter getRS1(){return master ? dest.getRS1() : null;}
	@Override
	public ReadStreamWriter getRS2(){return master ? dest.getRS2() : null;}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	/** Thread-safe counter tracking number of slave nodes that have terminated */
	protected final AtomicInteger terminatedCount=new AtomicInteger(0);
	/** Self-reference for inner class access */
	protected final ConcurrentReadOutputStreamD thisPointer=this;
	
	/** Wrapped destination of reads.  Null for slaves. */
	protected ConcurrentReadOutputStream dest;
	/** True if this node is the master that writes to files */
	protected final boolean master;
	/** Total number of MPI ranks in the computation */
	/** This node's MPI rank identifier */
	protected final int rank, ranks;
	
}
