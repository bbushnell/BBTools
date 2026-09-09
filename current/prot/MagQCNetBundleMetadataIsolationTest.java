package prot;

import ml.CellNet;
import java.nio.file.*;
import java.util.Arrays;
import java.util.concurrent.*;

/** THREAD_SAFE_BUNDLE_INFERENCE_20260909.md item 3 ("Exercise copied metadata isolation") plus,
 * folded in per Yoimiya's scoping (2026-09-09), an empirical regression check for the
 * executionDims fix (root's read-only checkpoint review, same day): Subnet.dims is public+final,
 * but the array it points to is still mutable in place, and netForCurrentThread() used to pass
 * that same array by reference into every cold InferenceNet construction -- so an external
 * mutation of the public field could have silently corrupted a FUTURE cold construction on a
 * different thread. The fix (a separate private executionDims field, cloned independently at
 * Subnet construction) is proven here by actually mutating the public dims field from one thread
 * and then verifying a cold construction on a DIFFERENT thread is still correctly shaped and
 * scores correctly.
 *
 * <p>Every blocking wait below is BOUNDED (Yoimiya's review, 2026-09-09): an unbounded
 * barrier/latch await hangs forever, never reports failure, if a worker throws before reaching
 * the synchronization point -- the crash-loud-never-hang violation this package's own review
 * culture exists to catch. Every await has a timeout; a timeout is itself a loud, diagnosed
 * AssertionError, and workers still alive after a timeout are interrupted and given a bounded
 * join before the failure is reported.
 *
 * <p>Test 1 (tags/metadata isolation): CellNet.setFrom() copies tags via
 * `tags=new LinkedHashMap<>(cn.tags)` (ml/CellNet.java:1357) -- a fresh copy, not a shared
 * reference -- verified here empirically via the public setTag()/getTag() pair rather than
 * resting on the source read alone: two threads each obtain a CellNet for the SAME Subnet
 * (guaranteed distinct via ThreadLocal), one sets a tag, the other's read must be unaffected.
 *
 * <p>Test 2 (dims-mutation cold-thread isolation): thread A gets a cold copy and its baseline
 * ALL-outputs score (via scoreAll(), not score() -- score() throws for any subnet with
 * expectedOutputs!=4, so this must work uniformly for legacy AND four-output fixtures, per
 * Yoimiya's review); the PUBLIC subnet.dims array is then mutated in place; thread B (a different
 * physical thread, so a genuinely cold ThreadLocal miss) gets its own cold copy AFTER the
 * mutation and must still be shaped and score correctly (every output bit-identical) -- proving
 * construction used the isolated executionDims, not the corrupted public field. */
public final class MagQCNetBundleMetadataIsolationTest {
	static final long TIMEOUT_S=10L;

	public static void main(String[] args) throws Exception {
		String bundlePath=null;
		for(String a:args) {
			int e=a.indexOf('='); if(e<1) throw new IllegalArgumentException("Expected key=value: "+a);
			String k=a.substring(0,e).toLowerCase(), v=a.substring(e+1);
			if(k.equals("bundle")) bundlePath=v;
		}
		if(bundlePath==null) throw new IllegalArgumentException("Required: bundle= "+
			"(an already-packed bundle path, any schema, legacy or four-output)");

		testTagIsolation(bundlePath);
		testDimsMutationColdThreadIsolation(bundlePath);
		System.out.println("MAGQC_BUNDLE_METADATA_ISOLATION PASS");
	}

	static void testTagIsolation(String bundlePath) throws Exception {
		MagQCNetBundle b=MagQCNetBundle.loadMultiOutput(Paths.get(bundlePath));
		final MagQCNetBundle.Subnet subnet=b.subnet(0);
		final String probeKey="isolation_probe_key_20260909";
		final String probeValue="thread_A_only_value";

		final CyclicBarrier gotOwnInstance=new CyclicBarrier(2);
		final CyclicBarrier aHasSetTag=new CyclicBarrier(2);
		final CountDownLatch done=new CountDownLatch(2);
		final Object[] results=new Object[6];
		final Throwable[] failure=new Throwable[1];

		Thread threadA=new Thread(new Runnable(){ public void run(){
			try {
				CellNet a=subnet.netForCurrentThread();
				results[0]=System.identityHashCode(a);
				results[2]=a.getTag(probeKey);
				gotOwnInstance.await(TIMEOUT_S,TimeUnit.SECONDS);
				a.setTag(probeKey,probeValue);
				aHasSetTag.await(TIMEOUT_S,TimeUnit.SECONDS);
				results[4]=a.getTag(probeKey);
			} catch(Throwable t) { synchronized(failure){ if(failure[0]==null) failure[0]=t; } }
			finally { done.countDown(); }
		}},"probe-thread-A");

		Thread threadB=new Thread(new Runnable(){ public void run(){
			try {
				CellNet bNet=subnet.netForCurrentThread();
				results[1]=System.identityHashCode(bNet);
				results[3]=bNet.getTag(probeKey);
				gotOwnInstance.await(TIMEOUT_S,TimeUnit.SECONDS);
				aHasSetTag.await(TIMEOUT_S,TimeUnit.SECONDS);
				results[5]=bNet.getTag(probeKey);
			} catch(Throwable t) { synchronized(failure){ if(failure[0]==null) failure[0]=t; } }
			finally { done.countDown(); }
		}},"probe-thread-B");

		threadA.start(); threadB.start();
		awaitBounded(done,TIMEOUT_S,new Thread[]{threadA,threadB},"tag-isolation threads to finish");
		joinBounded(threadA,threadB);
		if(failure[0]!=null) throw new RuntimeException("tag isolation worker failed",failure[0]);

		check(!results[0].equals(results[1]),"expected two DISTINCT CellNet identities (one per thread), "+
			"got the same object for both threads -- ThreadLocal caching is broken");
		check(results[2]==null,"A's getTag() before any setTag() call should be null, got "+results[2]);
		check(results[3]==null,"B's getTag() before any setTag() call should be null, got "+results[3]);
		check(probeValue.equals(results[4]),"A's own instance should see the tag it set itself, expected "+
			probeValue+" got "+results[4]);
		check(results[5]==null,"METADATA ISOLATION FAILURE: B's getTag() after A's setTag() should still "+
			"be null (distinct tags maps, per setFrom()'s defensive copy), got "+results[5]);

		System.out.println("tag isolation: two InferenceNet instances for the SAME Subnet from two "+
			"different threads (distinct identities "+results[0]+" vs "+results[1]+") have fully "+
			"independent tags maps -- setTag() on A had zero effect on B");
	}

	static void testDimsMutationColdThreadIsolation(String bundlePath) throws Exception {
		MagQCNetBundle b=MagQCNetBundle.loadMultiOutput(Paths.get(bundlePath));
		final MagQCNetBundle.Subnet subnet=b.subnet(0);
		final int origExpectedInputs=subnet.expectedInputs;
		final int origExpectedOutputs=subnet.expectedOutputs;
		final int[] origDimsSnapshot=subnet.dims.clone(); // for reporting only, never re-used for scoring

		final float[] probe=new float[origExpectedInputs];
		for(int i=0;i<probe.length;i++) probe[i]=(float)((i*5)%11)/13f;

		final float[][] baselineScoreHolder=new float[1][];
		final Throwable[] failure=new Throwable[1];
		final CountDownLatch aDone=new CountDownLatch(1);
		final CountDownLatch mutateDone=new CountDownLatch(1);
		final CountDownLatch bDone=new CountDownLatch(1);

		// Thread A: cold copy + baseline ALL-outputs score, BEFORE any mutation of the public dims
		// field. scoreAll() (not score()) so this works uniformly for legacy (1-output) AND
		// four-output subnets -- score() throws IllegalStateException for expectedOutputs!=1.
		Thread threadA=new Thread(new Runnable(){ public void run(){
			try {
				CellNet a=subnet.netForCurrentThread(); // cold on this (new) thread
				check(a.numInputs()==origExpectedInputs,"thread A cold copy has wrong input width "+
					"BEFORE mutation: got="+a.numInputs()+" expected="+origExpectedInputs);
				check(a.numOutputs()==origExpectedOutputs,"thread A cold copy has wrong output width "+
					"BEFORE mutation: got="+a.numOutputs()+" expected="+origExpectedOutputs);
				baselineScoreHolder[0]=subnet.scoreAll(probe.clone());
			} catch(Throwable t) { synchronized(failure){ if(failure[0]==null) failure[0]=t; } }
			finally { aDone.countDown(); }
		}},"dims-probe-thread-A");
		threadA.start();
		awaitBounded(aDone,TIMEOUT_S,new Thread[]{threadA},"thread A's baseline cold construction/score");
		if(failure[0]!=null) throw new RuntimeException("thread A failed",failure[0]);

		// Mutate the PUBLIC dims array IN PLACE -- exactly the vulnerability executionDims fixes.
		// A real caller would never legitimately do this; the test exists to prove that if some
		// caller (buggy or malicious) did, it could not corrupt a FUTURE cold construction.
		subnet.dims[0]=999999;
		check(subnet.dims[0]==999999,"mutation of the public dims field did not take effect -- test setup invalid");
		mutateDone.countDown();

		// Thread B: a DIFFERENT physical thread -> genuinely cold ThreadLocal miss -> exercises
		// netForCurrentThread()'s construction path AFTER the public field was corrupted.
		Thread threadB=new Thread(new Runnable(){ public void run(){
			try {
				boolean gotSignal=mutateDone.await(TIMEOUT_S,TimeUnit.SECONDS);
				if(!gotSignal) throw new AssertionError("thread B timed out waiting for the mutation signal");
				CellNet bNet=subnet.netForCurrentThread(); // cold on this (different, new) thread
				check(bNet.numInputs()==origExpectedInputs,"EXECUTIONDIMS ISOLATION FAILURE: thread B's "+
					"cold copy has wrong input width AFTER public-dims mutation: got="+bNet.numInputs()+
					" expected="+origExpectedInputs+" (mutated public dims[0]=999999) -- cold construction "+
					"used the corrupted public field instead of the isolated executionDims copy");
				check(bNet.numOutputs()==origExpectedOutputs,"EXECUTIONDIMS ISOLATION FAILURE: thread B's "+
					"cold copy has wrong output width AFTER public-dims mutation: got="+bNet.numOutputs()+
					" expected="+origExpectedOutputs);
				float[] bScore=subnet.scoreAll(probe.clone());
				check(bitsEqual(bScore,baselineScoreHolder[0]),"EXECUTIONDIMS ISOLATION FAILURE: thread B's "+
					"post-mutation scoreAll() ("+Arrays.toString(bScore)+") does not bit-match thread A's "+
					"pre-mutation baseline ("+Arrays.toString(baselineScoreHolder[0])+") on the same input "+
					"-- cold construction was corrupted by the public-dims mutation");
			} catch(Throwable t) { synchronized(failure){ if(failure[0]==null) failure[0]=t; } }
			finally { bDone.countDown(); }
		}},"dims-probe-thread-B");
		threadB.start();
		awaitBounded(bDone,TIMEOUT_S,new Thread[]{threadB},"thread B's post-mutation cold construction/score");
		joinBounded(threadB);
		if(failure[0]!=null) throw new RuntimeException("thread B failed",failure[0]);

		System.out.println("dims-mutation cold-thread isolation: mutated the PUBLIC Subnet.dims field "+
			"in place (dims[0] "+origDimsSnapshot[0]+" -> 999999) between thread A's and thread B's cold "+
			"netForCurrentThread() calls; thread B's cold construction (on a DIFFERENT physical thread, "+
			"AFTER the mutation) still had the correct input/output width ("+origExpectedInputs+"/"+
			origExpectedOutputs+") and scored bit-identically ("+Arrays.toString(baselineScoreHolder[0])+
			") to thread A's pre-mutation baseline across all "+origExpectedOutputs+" output(s) -- confirms "+
			"executionDims genuinely isolates live construction from the public field's mutability");
	}

	/** Bounded latch wait: on timeout, interrupts the given workers, gives them a short bounded
	 * join, then throws a loud, diagnosed AssertionError naming which (if any) are still stuck --
	 * converts an otherwise-infinite hang (a worker threw before reaching a synchronization point)
	 * into a bounded, explained failure. */
	static void awaitBounded(CountDownLatch latch,long timeoutS,Thread[] workers,String what) throws Exception {
		boolean ok=latch.await(timeoutS,TimeUnit.SECONDS);
		if(!ok) {
			for(Thread t:workers) t.interrupt();
			for(Thread t:workers) t.join(2000);
			StringBuilder stuck=new StringBuilder();
			for(Thread t:workers) if(t.isAlive()) stuck.append(t.getName()).append(' ');
			throw new AssertionError("TIMEOUT waiting for "+what+" after "+timeoutS+"s -- likely a hung "+
				"barrier/latch (a worker threw before reaching a synchronization point, or is deadlocked). "+
				"Interrupted workers; still alive after interrupt+2s join: "+
				(stuck.length()==0 ? "(none -- they exited after interrupt)" : stuck));
		}
	}
	static void joinBounded(Thread... workers) throws InterruptedException {
		for(Thread t:workers) t.join(5000);
		for(Thread t:workers) check(!t.isAlive(),"worker thread "+t.getName()+" still alive 5s after its "+
			"latch/barrier signalled completion -- unexpected stuck state");
	}
	static boolean bitsEqual(float[] a,float[] b) {
		if(a==null||b==null||a.length!=b.length) return false;
		for(int i=0;i<a.length;i++) if(Float.floatToIntBits(a[i])!=Float.floatToIntBits(b[i])) return false;
		return true;
	}
	static void check(boolean ok,String what){ if(!ok) throw new AssertionError("FAIL: "+what); }
}
