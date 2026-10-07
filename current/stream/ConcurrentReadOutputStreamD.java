package stream;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.Shared;
import structures.ListNum;

/**
 * Experimental output wrapper with unfinished MPI transport.
 * Master submissions delegate directly to the local destination; remote submissions,
 * listener threads and cross-rank completion depend on concrete methods that throw
 * TODO exceptions after their role assertions. Master listeners are created only
 * for other ranks. This is dormant scaffolding, not operational distributed output.
 * The factory selects it with MPI enabled and USE_CRISMPI false; the default true
 * setting selects a separate disabled factory branch.
 *
 * This wrapper does not copy submitted lists or Read payloads. Ordering and local
 * serialization belong to the destination; close, join and ordering resets require
 * caller coordination after submissions finish. Normal close/join and the terminal convention
 * still need transport/lifecycle work; the existing abort path has a monitor limitation.
 *
 * @author Brian Bushnell
 * @date Jan 26, 2015
 */
public class ConcurrentReadOutputStreamD extends ConcurrentReadOutputStream{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Retains the destination and captures rank settings without starting writers or listeners.
	 * Factory-created non-master ranks have null destination/formats; the base permits this.
	 * @param cros_ Local destination for a master, null for a non-master rank
	 * @param master_ Whether this instance owns the local destination
	 * @throws AssertionError If role and destination presence disagree while assertions are enabled */
	public ConcurrentReadOutputStreamD(ConcurrentReadOutputStream cros_, boolean master_){
		super(cros_==null ? null : cros_.ff1, cros_==null ? null : cros_.ff2);
		dest=cros_;
		master=master_;
		rank=Shared.MPI_RANK;
		ranks=Shared.MPI_NUM_RANKS;
		assert(master==(cros_!=null));
	}
	
	/** Marks this wrapper started; masters also start the destination and per-rank listeners.
	 * One-shot operation: a repeated call throws without resetting existing state.
	 * Non-master ranks only set the flag; concrete listener transport remains unimplemented. */
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
	
	/** Starts one listener for each other rank without retaining the thread references.
	 * Intended only for a master; close uses completion counts rather than these references. */
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
	 * Forwards borrowed references to the master destination, or attempts the non-master stub.
	 * Throws after abort; otherwise the wrapper monitor is held across the delegated call.
	 * This method does not copy payloads or translate list IDs.
	 * @param list Caller-owned list, subject to the destination's ownership contract
	 * @param listnum Destination ordering ID, forwarded unchanged
	 */
	@Override
	public synchronized void add(ArrayList<Read> list, long listnum){
		if(aborted){throw new RuntimeException("Cannot add list "+listnum+" to an aborted distributed output stream.");}
		if(master){dest.add(list, listnum);}else{unicast(list, listnum, 0);}
	}

	/** Returns after abort; otherwise attempts role-specific orderly closure.
	 * Masters count their own close, wait for all listener completions, then close the destination.
	 * Other ranks attempt an empty Long.MAX_VALUE-ID message through an unimplemented stub.
	 * This protocol is unfinished and does not guarantee completion. */
	@Override
	public void close(){
		if(aborted){return;}
		if(master){
			int count=terminatedCount.incrementAndGet();
			//This loop blocks until every ListenThread increments terminatedCount; a ListenThread that dies WITHOUT incrementing (see #001) hangs it forever.
			while(count<ranks){
				synchronized(terminatedCount){
					count=terminatedCount.intValue();
					if(count<ranks){
						try{
							terminatedCount.wait(1000);
						}catch(InterruptedException e){
							e.printStackTrace();
						}
					}
				}
			}
			dest.close();
		}else{
			//TODO: Probable bug - STR262: this positive close ID does not match the listener's
			//negative-ID terminal predicate if a future transport forwards it unchanged.
			//Concrete transport is still unimplemented; reconcile this source-only protocol concern first.
			unicast(new ListNum<Read>(new ArrayList<Read>(1), Long.MAX_VALUE), 0);
		}
	}

	/**
	 * Aborts the wrapped destination when this wrapper monitor is available.
	 * TODO STR263: add() and abort() are both synchronized on this wrapper; add() can
	 * hold the monitor while dest.add() waits, so abort() cannot then acquire it
	 * and does not guarantee wakeup of that blocked producer. A general MPI
	 * implementation or wrapper-lock redesign is outside this bounded port.
	 */
	@Override
	public synchronized void abort(){
		if(aborted){return;}
		errorState=true;
		finishedSuccessfully=false;
		aborted=true;
		if(master && dest!=null){dest.abort();}
	}

	/**
	 * Attempts completion coordination; normal transport paths currently throw TODO exceptions.
	 * Masters join the destination before calling the broadcast stub; non-master ranks
	 * call the listen stub. After abort, masters only join the destination and other
	 * ranks return, without transport. This method does not request normal close.
	 */
	@Override
	public void join(){
		if(aborted){
			if(master && dest!=null){dest.join();}
			return;
		}
		if(master){
			dest.join();
			broadcastJoin(true);
		}else{
			boolean b=listenForJoin();
			assert(b);
		}
	}

	/** Forwards the master's ordering reset, then clears its completion count and success flag.
	 * Non-master ranks do nothing. Does not reset started/aborted or restart listeners;
	 * callers must coordinate this operation with all previous submissions and listeners. */
	@Override
	public synchronized void resetNextListID(){
		if(master){
			dest.resetNextListID();
			terminatedCount.set(0);
			finishedSuccessfully=false;
		}
	}
	
	/** Reads the primary descriptor without a null check.
	 * @return Primary output name when a descriptor is present
	 * @throws NullPointerException If no primary descriptor exists, as on factory-created non-master ranks */
	@Override
	public String fname(){return ff1.name();}
	
	/**
	 * Reads current stored error flags without waiting for completion.
	 * Masters also consult their destination; this is not cross-rank error collection.
	 * @return true after abort, or when the observed local/master-destination flag is set
	 */
	@Override
	public boolean errorState(){
		if(aborted){return true;}
		if(master){return errorState || dest.errorState();}else{return errorState;}
	}

	/**
	 * Attempts a role-specific status exchange while holding this wrapper's monitor.
	 * Masters first read destination status. Both concrete transport paths then throw;
	 * an aborted wrapper instead returns false without a status exchange.
	 * @return false after abort; the concrete non-aborted transport path does not return
	 */
	@Override
	public boolean finishedSuccessfully(){
		synchronized(this){
			if(aborted){return false;}
			if(master){
				finishedSuccessfully=dest.finishedSuccessfully();
				broadcastFinishedSuccessfully(finishedSuccessfully);
			}else{
				finishedSuccessfully=listenFinishedSuccessfully();
			}
		}
		return !aborted && finishedSuccessfully;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Wraps a borrowed list and forwards its ID to the unimplemented transport method.
	 * @param list Borrowed payload list
	 * @param listnum Unchanged destination ordering ID
	 * @param i Intended destination rank */
	private void unicast(ArrayList<Read> list, long listnum, int i){unicast(new ListNum<Read>(list, listnum), i);}
	
	/** Non-master send placeholder; performs no transport and ends in a TODO exception.
	 * @param ln Intended list/ID payload, unused by this stub
	 * @param i Intended destination rank, used only for diagnostics */
	protected void unicast(ListNum<Read> ln, int i){
		if(verbose){System.err.println("crosD "+(master ? "master" : "slave ")+":    Unicasting reads to "+i+".");}
		assert(!master);
		
		boolean success=false;
		while(!success){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	/** Master receive placeholder; asserts role then ends in a TODO exception.
	 * @param i Intended source rank, used only for diagnostics
	 * @return No value from the concrete implementation */
	protected ListNum<Read> listen(int i){
		if(verbose){System.err.println("crosD "+(master ? "master" : "slave ")+":    Listening for reads from "+i+".");}
		assert(master);
		
		boolean success=false;
		while(!success){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	
	/**
	 * Non-master status placeholder with no implemented transport.
	 * @return No value from the concrete implementation
	 * @throws RuntimeException TODO exception after the non-master-role assertion
	 */
	protected boolean listenFinishedSuccessfully(){
		if(verbose){System.err.println("crosD "+(master ? "master" : "slave ")+":    listenFinishedSuccessfully.");}
		assert(!master);
		
		boolean success=false;
		while(!success){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}

	/**
	 * Master status-broadcast placeholder with no implemented transport.
	 * @param b Intended success status, unused by this stub
	 * @throws RuntimeException TODO exception after the master-role assertion
	 */
	protected void broadcastFinishedSuccessfully(boolean b){
		if(verbose){System.err.println("crosD "+(master ? "master" : "slave ")+":    broadcastFinishedSuccessfully.");}
		assert(master);
		
		boolean success=false;
		while(!success){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}
	
	/** Master completion-broadcast placeholder, ending in a TODO exception after its role assertion.
	 * @param b Intended completion signal, unused by this stub */
	protected void broadcastJoin(boolean b){
		if(verbose){System.err.println("crosD "+(master ? "master" : "slave ")+":    broadcastJoin.");}
		assert(master);
		
		boolean success=false;
		while(!success){
			try{
				//Do some MPI stuff
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		throw new RuntimeException("TODO");
	}

	/** Non-master completion placeholder, ending in a TODO exception after its role assertion.
	 * @return No value from the concrete implementation */
	protected boolean listenForJoin(){
		if(verbose){System.err.println("crosD "+(master ? "master" : "slave ")+":    listenForJoin.");}
		assert(!master);
		
		boolean success=false;
		while(!success){
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
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Master-side receiver for one other rank; the concrete receive operation is a stub. */
	private class ListenThread extends Thread{
		
		/** Captures one other rank and checks its range with assertions.
		 * @param sourceNum_ Rank whose submitted lists would be forwarded to the destination */
		ListenThread(int sourceNum_){
			sourceNum=sourceNum_;
			assert(sourceNum_!=rank);
			assert(sourceNum>=0 && sourceNum<ranks);
		}
		
		/** Forwards nonnegative-ID lists, then counts completion after null or a negative ID.
		 * The concrete listen stub throws; completion accounting is not protected by finally. */
		@Override
		public void run(){
			assert(master);
			//TODO: Possible bug [stream/ConcurrentReadOutputStreamD#001, STR264] - a listener
			//exception skips completion accounting, leaving master close waiting for this rank.
			//Concrete listen() is still a throwing stub. Preserve a structural completion/error
			//handoff (or a loud unrecoverable failure) when transport is implemented; deferred here.
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
		
		/** Captured nonlocal source rank. */
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

	/** Counts the master's close and normally finished listeners; also used as their wait monitor. */
	protected final AtomicInteger terminatedCount=new AtomicInteger(0);
	/** Retained self-reference; unused by this concrete implementation. */
	protected final ConcurrentReadOutputStreamD thisPointer=this;
	
	/** Wrapped destination of reads. Null for slaves. */
	protected ConcurrentReadOutputStream dest;
	/** Captured role controlling destination access and transport attempts. */
	protected final boolean master;
	/** Captured current rank and total rank count; constructor does not validate their range. */
	protected final int rank, ranks;
	/** Latched abort status; ordering reset does not clear it. */
	private volatile boolean aborted=false;
	
}
