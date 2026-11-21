package test.jgi;

import java.lang.management.LockInfo;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadInfo;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.LockSupport;

import jgi.BBMask;
import parse.Parse;

/**
 * Finite BBMask worker-join interruption checks for private awaitWorkers.
 * @author Jean
 */
public final class BBMaskWorkerInterruptionTest{

	/** Requires phase=baseline|fixed and a new outdir. */
	public static void main(final String[] args) throws Throwable{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown BBMask worker interruption option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/bbmaskworkerinterruptiontest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		require(Files.isDirectory(repo.resolve("current")), "Missing BBTools current/ at "+repo);
		Files.createDirectory(out);
		final BBMaskWorkerInterruptionTest test=new BBMaskWorkerInterruptionTest(out,
			phase.equals("baseline"));
		for(final Scenario scenario : SCENARIOS){test.run(scenario);}
		test.check(test.cases==3, "Expected three awaitWorkers cases: "+test.cases);
		final String summary="BBMASK_WORKER_INTERRUPTION_"+
			(test.failures==0 ? "PASS" : "FAIL")+" phase="+phase+
			" cases="+test.cases+" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Binds the reflected helper and writes the evidence table. */
	private BBMaskWorkerInterruptionTest(final Path out_, final boolean baseline_)
			throws Exception{
		out=out_;
		baseline=baseline_;
		awaitWorkers=BBMask.class.getDeclaredMethod("awaitWorkers", Thread[].class);
		awaitWorkers.setAccessible(true);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\tworkers\tinterrupts"+
			"\treturned\texpectedReturned\tsuppressed\texpectedSuppressed"+
			"\tsuppressedTypes\townerInterrupted\tworkersTerminatedAtReturn"+
			"\tworker1Released\tworker2Released\tmatch\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Runs one live-worker owner interruption case with bounded cleanup. */
	private void run(final Scenario scenario) throws Throwable{
		cases++;
		final LatchWorker[] workers=new LatchWorker[scenario.workers];
		for(int i=0; i<workers.length; i++){workers[i]=new LatchWorker();}
		final Owner owner=new Owner(workers);
		owner.setDaemon(true);
		for(final LatchWorker worker : workers){worker.setDaemon(true);}
		boolean worker1Released=false;
		boolean worker2Released=false;
		Throwable failure=null;
		try{
			for(final LatchWorker worker : workers){worker.start();}
			for(final LatchWorker worker : workers){await(worker.entered, "worker entered");}
			owner.start();

			waitForOwnerWaitingOn(owner, workers[0]);
			check(!owner.done(), scenario.name+": owner returned before worker1 release");
			if(scenario.interruptFirst){
				owner.interrupt();
				waitForOwnerInterruptConsumed(owner, workers[0]);
				check(workers[0].isAlive(), scenario.name+": worker1 died before release");
			}
			workers[0].release.countDown();
			worker1Released=true;

			if(workers.length>1){
				waitForOwnerWaitingOn(owner, workers[1]);
				check(!owner.done(), scenario.name+": owner returned before worker2 release");
				if(scenario.interruptSecond){
					owner.interrupt();
					waitForOwnerInterruptConsumed(owner, workers[1]);
					check(workers[1].isAlive(), scenario.name+": worker2 died before release");
				}
				workers[1].release.countDown();
				worker2Released=true;
			}
		}catch(Throwable t){
			failure=t;
		}finally{
			for(final LatchWorker worker : workers){worker.release.countDown();}
			try{
				join(owner, "owner");
				for(int i=0; i<workers.length; i++){join(workers[i], "worker"+i);}
			}catch(Throwable t){
				if(failure==null){failure=t;}
				else{failure.addSuppressed(t);}
			}
		}
		if(failure!=null){throw failure;}

		for(final LatchWorker worker : workers){
			check(worker.failure==null, scenario.name+": worker failure "+worker.failure);
			check(worker.getState()==Thread.State.TERMINATED,
				scenario.name+": worker not terminated: "+worker.getState());
		}
		check(owner.thrown==null, scenario.name+": awaitWorkers threw "+owner.thrown);
		check(owner.workersTerminatedAtReturn,
			scenario.name+": awaitWorkers returned before all workers terminated");

		final boolean expectedReturned=!baseline && scenario.interruptions()>0;
		final int expectedSuppressed=expectedReturned ? scenario.interruptions()-1 : 0;
		final int suppressed=owner.result==null ? 0 : owner.result.getSuppressed().length;
		final boolean returned=owner.result!=null;
		final boolean suppressedTypes=suppressedAreInterrupted(owner.result);
		final boolean match=returned==expectedReturned && suppressed==expectedSuppressed &&
			suppressedTypes && !owner.interruptedAfterReturn && owner.done() &&
			owner.workersTerminatedAtReturn && worker1Released &&
			(workers.length<2 || worker2Released);
		Files.write(out.resolve("cases.tsv"), (scenario.name+"\t"+
			(baseline ? "baseline" : "fixed")+"\t"+workers.length+"\t"+
			scenario.interruptions()+"\t"+returned+"\t"+expectedReturned+"\t"+
			suppressed+"\t"+expectedSuppressed+"\t"+suppressedTypes+"\t"+
			owner.interruptedAfterReturn+"\t"+owner.workersTerminatedAtReturn+"\t"+
			worker1Released+"\t"+worker2Released+"\t"+match+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check(returned==expectedReturned, scenario.name+": returned "+returned+
			" != "+expectedReturned);
		check(suppressed==expectedSuppressed, scenario.name+": suppressed "+suppressed+
			" != "+expectedSuppressed);
		check(suppressedTypes, scenario.name+": a suppressed entry was not InterruptedException");
		check(!owner.interruptedAfterReturn,
			scenario.name+": awaitWorkers restored interrupt status internally");
		check(owner.workersTerminatedAtReturn,
			scenario.name+": awaitWorkers did not join every owned worker");
		check(match, scenario.name+": case did not match oracle");
	}

	/** Invokes the private static BBMask.awaitWorkers(Thread[]) helper. */
	private InterruptedException invokeAwaitWorkers(final Thread[] workers) throws Throwable{
		try{
			return (InterruptedException)awaitWorkers.invoke(null, new Object[]{workers});
		}catch(InvocationTargetException e){
			throw e.getCause();
		}
	}

	/** Waits until the owner is waiting in Thread.join on this exact worker monitor. */
	private static void waitForOwnerWaitingOn(final Thread owner, final Thread worker)
			throws Exception{
		final ThreadMXBean bean=ManagementFactory.getThreadMXBean();
		final int target=System.identityHashCode(worker);
		final long deadline=System.nanoTime()+TIMEOUT_NANOS;
		while(System.nanoTime()<deadline){
			final ThreadInfo info=bean.getThreadInfo(owner.getId());
			if(info==null){throw new AssertionError("Owner terminated before joining "+worker);}
			final LockInfo lock=info.getLockInfo();
			if(info.getThreadState()==Thread.State.WAITING && lock!=null &&
					lock.getIdentityHashCode()==target){
				return;
			}
			LockSupport.parkNanos(PARK_NANOS);
		}
		throw new AssertionError("Timed out waiting for owner to join "+worker+
			"; state="+owner.getState());
	}

	/** Waits until join has consumed the owner's interrupt status. */
	private static void waitForOwnerInterruptConsumed(final Thread owner, final Thread worker)
			throws Exception{
		final long deadline=System.nanoTime()+TIMEOUT_NANOS;
		while(System.nanoTime()<deadline){
			if(!owner.isAlive()){
				throw new AssertionError("Owner terminated before consuming interrupt on "+worker);
			}
			if(!owner.isInterrupted()){return;}
			LockSupport.parkNanos(PARK_NANOS);
		}
		throw new AssertionError("Timed out waiting for owner interrupt consumption on "+
			worker+"; state="+owner.getState());
	}

	/** Waits for a latch with a finite timeout. */
	private static void await(final CountDownLatch latch, final String name) throws Exception{
		if(!latch.await(TIMEOUT_NANOS, TimeUnit.NANOSECONDS)){
			throw new AssertionError("Timed out waiting for "+name);
		}
	}

	/** Joins a thread with a finite timeout. */
	private static void join(final Thread thread, final String name) throws Exception{
		final long deadline=System.nanoTime()+TIMEOUT_NANOS;
		while(thread.isAlive() && System.nanoTime()<deadline){thread.join(10);}
		if(thread.isAlive()){throw new AssertionError("Timed out joining "+name);}
	}

	/** Records one assertion without hiding later evidence after a failure. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}

	/** Hard construction-time requirement. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Requires every suppressed join result to be an InterruptedException. */
	private static boolean suppressedAreInterrupted(final InterruptedException result){
		if(result==null){return true;}
		for(final Throwable t : result.getSuppressed()){
			if(!(t instanceof InterruptedException)){return false;}
		}
		return true;
	}

	/** Owner thread that calls BBMask.awaitWorkers. */
	private final class Owner extends Thread{
		Owner(final Thread[] workers_){workers=workers_;}

		@Override
		public void run(){
			try{
				result=invokeAwaitWorkers(workers);
				workersTerminatedAtReturn=allTerminated(workers);
			}catch(Throwable t){
				thrown=t;
			}finally{
				interruptedAfterReturn=Thread.currentThread().isInterrupted();
				done.countDown();
			}
		}

		boolean done(){return done.getCount()==0;}

		final Thread[] workers;
		final CountDownLatch done=new CountDownLatch(1);
		volatile InterruptedException result;
		volatile Throwable thrown;
		volatile boolean interruptedAfterReturn, workersTerminatedAtReturn;
	}

	/** Worker that stays live until the fixture releases it. */
	private static final class LatchWorker extends Thread{
		@Override
		public void run(){
			try{
				entered.countDown();
				release.await();
			}catch(Throwable t){
				failure=t;
			}
		}

		final CountDownLatch entered=new CountDownLatch(1);
		final CountDownLatch release=new CountDownLatch(1);
		volatile Throwable failure;
	}

	/** One phase-independent awaitWorkers owner-interruption case. */
	private static final class Scenario{
		Scenario(final String name_, final int workers_, final boolean interruptFirst_,
				final boolean interruptSecond_){
			name=name_;
			workers=workers_;
			interruptFirst=interruptFirst_;
			interruptSecond=interruptSecond_;
		}

		int interruptions(){return (interruptFirst ? 1 : 0)+(interruptSecond ? 1 : 0);}

		final String name;
		final int workers;
		final boolean interruptFirst, interruptSecond;
	}

	private static final Scenario[] SCENARIOS={
		new Scenario("no-interrupt", 1, false, false),
		new Scenario("one-interrupt", 1, true, false),
		new Scenario("two-interrupts-two-workers", 2, true, true)
	};
	private static final long TIMEOUT_NANOS=10L*1000L*1000L*1000L;
	private static final long PARK_NANOS=1L*1000L*1000L;

	private static boolean allTerminated(final Thread[] threads){
		for(final Thread thread : threads){
			if(thread.getState()!=Thread.State.TERMINATED){return false;}
		}
		return true;
	}

	private final Path out;
	private final boolean baseline;
	private final Method awaitWorkers;
	private int cases, checks, failures;
}
