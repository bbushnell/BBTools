package fileIO;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Jan 2, 2013
 *
 */
public class LoadThread<X> extends Thread{
	
	public static <Y> LoadThread<Y> load(String fname, Class<Y> c){
		LoadThread<Y> lt=new LoadThread<Y>(fname, c);
		lt.start();
		return lt;
	}
	
	private LoadThread(String fname_, Class<X> c_){
		fname=fname_;
		c=c_;
		addThread(1);
	}
	
	@Override
	public void run(){
		addRunningThread(1);
		output=ReadWrite.read(c, fname, false);
		addRunningThread(-1);
		synchronized(this){this.notify();}
	}
	
	
	private static final int addThread(int x){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				activeThreads[0]+=x;
				activeThreads[1]+=x;
			}else{
				addRunningThread(x);
			}
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
					
			return activeThreads[0];
		}
	}
	
	private static final int addRunningThread(int x){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				assert(activeThreads[1]>=x);
				while(activeThreads[2]>=lim){
					try {
						activeThreads.wait();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				activeThreads[1]-=x; //Remove from waiting
			}else{
				activeThreads[0]+=x; //Remove from active
			}
			activeThreads[2]+=x; //Change number running
			
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
			
			if(activeThreads[2]==0 || (activeThreads[2]<lim && activeThreads[1]>0)){activeThreads.notify();}
//			System.err.println(activeThreads[2]);
//			try {
//				activeThreads.wait(5000);
//			} catch (InterruptedException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
			return activeThreads[2];
		}
	}
	
	public static final int countActiveThreads(){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
			return activeThreads[0];
		}
	}
	
	public static final void waitForReadingToFinish(){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			while(activeThreads[0]>0){
				assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
						activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
				try {
					activeThreads.wait(8000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(activeThreads[2]==0 || (activeThreads[2]<lim && activeThreads[1]>0)){activeThreads.notify();}
			}
		}
	}
	
	public final void waitForThisToFinish(){
		if(output==null){
			while(this.getState()!=State.TERMINATED){
				try {
					this.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	/** {active, waiting, running} <br>
	 * Active means running or waiting.
	 */
	public static int[] activeThreads={0, 0, 0};
	
	private final String fname;
	private final Class<X> c;
	public X output=null;
	
	private static final int[] RUNNING=new int[1];
	public static int LIMIT=Tools.min(12, Tools.max(Shared.threads(), 1));
	
}
