package test.structures;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import structures.CoverageArray2A;

/** Concurrent same-word lane oracle for STU-010 CoverageArray2A CAS repair.
 * @author Sucrose
 */
public final class CoverageArray2AConcurrentCase{
	
	private static final int ITERATIONS=10000;
	private static final int UPPER_LOC=1;
	private static final int LOWER_LOC=2;
	
	/** Runs the two-lane concurrent increment oracle. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		
		sameWordLaneIncrements(expect);
	}
	
	/** Verifies concurrent lower and upper lane updates do not lose retries. */
	private static void sameWordLaneIncrements(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		final CountDownLatch ready=new CountDownLatch(2);
		final CountDownLatch start=new CountDownLatch(1);
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		
		final Thread upper=worker("coverage-upper", ca, UPPER_LOC, ready, start, failure);
		final Thread lower=worker("coverage-lower", ca, LOWER_LOC, ready, start, failure);
		upper.start();
		lower.start();
		await("ready", ready);
		start.countDown();
		join(upper);
		join(lower);
		
		final Throwable thrown=failure.get();
		if(thrown!=null){
			throw new AssertionError("Worker failed", thrown);
		}
		
		checkEquals("upper-lane", ITERATIONS, ca.get(UPPER_LOC));
		checkEquals("lower-lane", ITERATIONS, ca.get(LOWER_LOC));
		pass("same-word-lane-increments",
				ca.get(UPPER_LOC)+","+ca.get(LOWER_LOC), expect);
	}
	
	/** Creates an increment worker that waits for the shared start gate. */
	private static Thread worker(final String name, final CoverageArray2A ca,
			final int loc, final CountDownLatch ready, final CountDownLatch start,
			final AtomicReference<Throwable> failure){
		return new Thread(new Runnable(){
			@Override
			public void run(){
				try{
					ready.countDown();
					start.await();
					for(int i=0; i<ITERATIONS; i++){
						ca.increment(loc);
					}
				}catch(Throwable t){
					failure.compareAndSet(null, t);
				}
			}
		}, name);
	}
	
	/** Awaits a latch, converting interruption to AssertionError. */
	private static void await(final String name, final CountDownLatch latch){
		try{
			latch.await();
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new AssertionError("Interrupted while waiting for "+name, e);
		}
	}
	
	/** Joins a worker, converting interruption to AssertionError. */
	private static void join(final Thread worker){
		try{
			worker.join();
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new AssertionError("Interrupted while joining "+worker.getName(), e);
		}
	}
	
	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}
	
	/** Prints one successful case row. */
	private static void pass(final String name, final String observed, final String expect){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}
	
	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
