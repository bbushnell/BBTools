package test.assemble;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import parse.Parser;

/** Two isolated checks of the gate fixture's worker lifetime and exception precedence.
 * Reflects into the private helper without adding a production test hook. No model,
 * FMA setting, sequence file or timing measurement is involved.
 * @author Brian Bushnell, Sucrose
 */
public final class LiveBranchGateLifecycleTest{

	/** Runs interrupted-join and partial-startup cases; recognizes forwarded JVM flags. */
	public static void main(final String[] args) throws Exception{
		for(String arg : args){
			if(!Parser.isJavaFlag(arg)){throw new IllegalArgumentException("Unexpected argument: "+arg);}
		}
		final Method method=LiveBranchGateTest.class.getDeclaredMethod("startAndWait", ArrayList.class, AtomicReference.class);
		method.setAccessible(true);
		test(method, false);
		test(method, true);
		System.out.println("LIVE_BRANCH_GATE_LIFECYCLE_TEST_PASS cases=2");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Holds one started worker until the interrupted caller is observed waiting again.
	 * Startup mode adds a failing start and a later NEW worker; ordinary mode records
	 * a worker failure only after release, so suppression requires completed joining. */
	private static void test(final Method method, final boolean startup) throws InterruptedException{
		final CountDownLatch ready=new CountDownLatch(1), release=new CountDownLatch(1);
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final AssertionError workerFailure=new AssertionError("worker sentinel");
		final RuntimeException startupFailure=new IllegalStateException("startup sentinel");
		final ArrayList<Thread> workers=new ArrayList<Thread>();
		final Thread worker=new Thread(new Runnable(){
			@Override public void run(){
				ready.countDown();
				try{
					release.await();
					if(!startup){throw workerFailure;}
				}catch(Throwable problem){failure.compareAndSet(null, problem);}
			}
		}, "gate-lifetime-worker");
		workers.add(worker);
		if(startup){
			workers.add(new Thread("gate-lifetime-failed-start"){
				@Override public synchronized void start(){
					Thread.currentThread().interrupt();
					throw startupFailure;
				}
			});
			workers.add(new Thread("gate-lifetime-not-started"));
		}
		final Caller caller=new Caller(method, workers, failure, ready, !startup);
		boolean waited=false;
		caller.start();
		try{
			ready.await();
			waited=awaitJoin(caller);
		}finally{
			release.countDown();
			caller.join();
			worker.join();
		}
		check(waited, "Caller escaped while its started worker was held by the release latch.");
		check(caller.allTerminated, "Helper exited before every started worker terminated.");
		check(worker.getState()==Thread.State.TERMINATED, "Started worker must be joined before this case completes.");
		check(caller.interrupted, "Caller must regain interrupt status after waiting.");
		if(startup){
			check(caller.problem==startupFailure, "Startup failure must remain the exact primary exception.");
			final Throwable[] suppressed=startupFailure.getSuppressed();
			check(suppressed.length==1 && suppressed[0] instanceof InterruptedException,
					"Interrupted cleanup must be suppressed on the startup failure.");
			check(failure.get()==null, "Startup case must not manufacture a worker failure.");
			check(workers.get(1).getState()==Thread.State.NEW && workers.get(2).getState()==Thread.State.NEW,
					"Failed and later startup entries must remain unstarted.");
		}else{
			check(caller.problem instanceof InterruptedException, "Interrupted join must remain the primary exception.");
			final Throwable[] suppressed=caller.problem.getSuppressed();
			check(suppressed.length==1 && suppressed[0]==workerFailure && failure.get()==workerFailure,
					"Interrupted join must preserve the exact worker failure captured before termination.");
		}
		System.out.println("PASS: "+(startup ? "partial-startup-interruption" : "interruption-worker-failure"));
	}

	/** Observes a join after the initial interrupt was consumed, or detects early exit.
	 * The caller has no other latch wait; only the worker waits on the release latch.
	 * Yielding avoids sleeps and timeout-based scheduling assumptions. */
	private static boolean awaitJoin(final Thread caller){
		while(caller.isAlive()){
			if(caller.getState()==Thread.State.WAITING && !caller.isInterrupted()){return true;}
			Thread.yield();
		}
		return false;
	}

	/** Keeps checks active regardless of JVM assertion mode. */
	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Invokes the real helper and records its exit state before harness cleanup joins.
	 * Fields are read only after this thread is joined, which publishes their values. */
	private static final class Caller extends Thread{
		/** Retains this case's invocation state and signals early exit as well as worker startup. */
		Caller(final Method method_, final ArrayList<Thread> workers_,
				final AtomicReference<Throwable> failure_, final CountDownLatch ready_, final boolean interruptBefore_){
			super("gate-lifetime-caller");
			method=method_; workers=workers_; failure=failure_; ready=ready_; interruptBefore=interruptBefore_;
		}
		@Override public void run(){
			if(interruptBefore){interrupt();}
			try{
				method.invoke(null, workers, failure);
			}catch(InvocationTargetException invocation){
				problem=invocation.getCause();
			}catch(Throwable unexpected){
				problem=unexpected;
			}finally{
				interrupted=isInterrupted();
				allTerminated=true;
				for(Thread worker : workers){allTerminated&=!worker.isAlive();}
				ready.countDown();
			}
		}
		/*--------------------------------------------------------------*/
		/*----------------            Fields            ----------------*/
		/*--------------------------------------------------------------*/

		/** Borrowed stable invocation inputs, confined to this case. */
		final Method method;
		final ArrayList<Thread> workers;
		final AtomicReference<Throwable> failure;
		final CountDownLatch ready;
		final boolean interruptBefore;
		/** Actual helper outcome and status observed before this caller terminates. */
		Throwable problem;
		boolean interrupted, allTerminated;
	}
}
