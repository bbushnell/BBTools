package test.sort;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;

import sort.Shuffle;
import sort.Shuffle.ShuffleThread;

/** Checks that worker completion and failed starts preserve Shuffle's thread accounting.
 * @author Shinobu
 */
public final class ShuffleCompletionTest{

	/** Runs without files; a latch holds one worker below the configured capacity. */
	public static void main(String[] args) throws InterruptedException{
		Shuffle.setMaxThreads(2);
		final CountDownLatch ready=new CountDownLatch(1), release=new CountDownLatch(1);
		final CountDownLatch entered=new CountDownLatch(1), returned=new CountDownLatch(1);
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final ShuffleThread worker=new ShuffleThread(null, null, null, null, Shuffle.SHUFFLE, true){
			@Override
			public void run(){
				ready.countDown();
				try{release.await();}
				catch(InterruptedException e){failure.set(e);}
				finally{Shuffle.addThread(-1);}
			}
		};
		final Thread waiter=new Thread(){
			@Override
			public void run(){
				entered.countDown();
				try{Shuffle.waitForFinish();}
				catch(Throwable t){failure.set(t);}
				finally{returned.countDown();}
			}
		};
		worker.setDaemon(true);
		waiter.setDaemon(true);
		worker.start();
		ready.await();
		boolean early;
		try{
			waiter.start();
			entered.await();
			early=returned.await(250, TimeUnit.MILLISECONDS);
		}finally{release.countDown();}
		worker.join(5000);
		waiter.join(5000);
		if(early || worker.isAlive() || waiter.isAlive() || failure.get()!=null){
			throw new AssertionError("Completion wait returned early, stalled, or failed: "+failure.get());
		}
		boolean rejected=false;
		try{worker.start();}
		catch(IllegalThreadStateException expected){rejected=true;}
		if(!rejected || Shuffle.addThread(0)!=0){
			throw new AssertionError("Failed thread restart leaked an active-worker slot");
		}
		boolean invalidLimit=false;
		try{Shuffle.setMaxThreads(0);}
		catch(IllegalArgumentException expected){invalidLimit=true;}
		if(!invalidLimit){throw new AssertionError("Zero worker limit would prevent all startup");}
		System.out.println("PASS\tcompletion-below-capacity\tfailed-start-accounting\tinvalid-limit");
	}
}
