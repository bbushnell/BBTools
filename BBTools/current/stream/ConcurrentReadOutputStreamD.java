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
			unicast(new ListNum<Read>(new ArrayList<Read>(1), -1), 0);
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
	
	private void unicast(ArrayList<Read> list, long listnum, int i) {
		unicast(new ListNum<Read>(list, listnum), i);
	}
	
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
	
	private class ListenThread extends Thread{
		
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

	protected final AtomicInteger terminatedCount=new AtomicInteger(0);
	protected final ConcurrentReadOutputStreamD thisPointer=this;
	
	/** Wrapped destination of reads.  Null for slaves. */
	protected ConcurrentReadOutputStream dest;
	protected final boolean master;
	protected final int rank, ranks;
	
}
