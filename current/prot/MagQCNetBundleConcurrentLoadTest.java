package prot;

import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/** THREAD_SAFE_BUNDLE_INFERENCE_20260909.md item 2: "Check this component's concurrent
 * dense/sparse bundle loads separately without an arbitrary uncoordinated external parser." This
 * targets construction-time concurrency (parseNet()/CellNetParser master-net construction, which
 * still mutates the shared ml.CellNet.DENSE global as a side effect of parsing) specifically --
 * a genuinely different code path from scoring-time concurrency (see
 * MagQCNetBundleConcurrentScoringTest). Multiple threads each call
 * MagQCNetBundle.loadMultiOutput() on DIFFERENT bundle paths (one legacy/all-dense, one
 * mixed/real-candidate-sparse) AT THE SAME TIME, verifying every load succeeds and produces the
 * correct canonical_payload_sha256 -- confirming the synchronized(CellNet.class) lock in
 * parseNet() actually prevents concurrent master-net construction from corrupting each other.
 *
 * <p>Every blocking wait is BOUNDED (Yoimiya's review, 2026-09-09): the original version used a
 * bare barrier.await() and Future.get() with no timeout -- if one worker threw before reaching
 * the barrier on some later round, its peers would wait at that barrier forever, and the failure
 * would never be reported (crash-loud-never-hang violation). It also cast
 * `(Exception)e.getCause()`, which throws ClassCastException the moment a worker throws
 * AssertionError (an Error, not an Exception -- exactly what this test's own check()/hash-mismatch
 * failures do), masking the real failure behind an unrelated cast error. And the pool was shut
 * down only on the success path, so a failing run left the pool's worker threads alive (non-daemon,
 * blocking JVM exit). All three are fixed here: bounded timeouts on every barrier.await()/
 * Future.get(), Throwable (not Exception) propagation with no cast, and pool.shutdownNow() in a
 * finally block so any failure path -- including a deliberately wrong expected hash -- exits
 * promptly instead of hanging. */
public final class MagQCNetBundleConcurrentLoadTest {
	static final long BARRIER_TIMEOUT_S=10L;
	static final long FUTURE_TIMEOUT_S=30L;

	public static void main(String[] args) throws Exception {
		String legacyPath=null, realPath=null, legacyHash=null, realHash=null;
		for(String a:args) {
			int e=a.indexOf('='); if(e<1) throw new IllegalArgumentException("Expected key=value: "+a);
			String k=a.substring(0,e).toLowerCase(), v=a.substring(e+1);
			if(k.equals("legacybundle")) legacyPath=v;
			else if(k.equals("realbundle")) realPath=v;
			else if(k.equals("legacyhash")) legacyHash=v;
			else if(k.equals("realhash")) realHash=v;
		}
		if(legacyPath==null||realPath==null||legacyHash==null||realHash==null)
			throw new IllegalArgumentException("Required: legacybundle= realbundle= legacyhash= realhash= "+
				"(two DIFFERENT already-packed bundle paths and their expected canonical_payload_sha256 values)");

		final String fLegacyPath=legacyPath, fRealPath=realPath, fLegacyHash=legacyHash, fRealHash=realHash;
		final int THREADS_PER_PATH=4, ITERS_PER_THREAD=20;
		final int TOTAL_THREADS=THREADS_PER_PATH*2;
		final CyclicBarrier barrier=new CyclicBarrier(TOTAL_THREADS);
		final AtomicInteger legacySuccesses=new AtomicInteger(0), realSuccesses=new AtomicInteger(0);

		ExecutorService pool=Executors.newFixedThreadPool(TOTAL_THREADS);
		try {
			List<Future<?>> futures=new ArrayList<Future<?>>();
			for(int t=0;t<THREADS_PER_PATH;t++){
				futures.add(pool.submit(new Runnable(){ public void run(){
					for(int i=0;i<ITERS_PER_THREAD;i++){
						try { barrier.await(BARRIER_TIMEOUT_S,TimeUnit.SECONDS); }
						catch(Exception be) { throw new RuntimeException("legacy worker barrier wait failed at iter "+i,be); }
						MagQCNetBundle b;
						try { b=MagQCNetBundle.loadMultiOutput(Paths.get(fLegacyPath)); }
						catch(Exception le) { throw new RuntimeException("legacy load failed at iter "+i,le); }
						String got=b.metadata("canonical_payload_sha256");
						if(!fLegacyHash.equals(got))
							throw new AssertionError("legacy load hash mismatch: got="+got+" expected="+fLegacyHash);
						if(b.size()!=176) throw new AssertionError("legacy load size mismatch: got="+b.size());
						legacySuccesses.incrementAndGet();
					}
				}}));
			}
			for(int t=0;t<THREADS_PER_PATH;t++){
				futures.add(pool.submit(new Runnable(){ public void run(){
					for(int i=0;i<ITERS_PER_THREAD;i++){
						try { barrier.await(BARRIER_TIMEOUT_S,TimeUnit.SECONDS); }
						catch(Exception be) { throw new RuntimeException("real worker barrier wait failed at iter "+i,be); }
						MagQCNetBundle b;
						try { b=MagQCNetBundle.loadMultiOutput(Paths.get(fRealPath)); }
						catch(Exception le) { throw new RuntimeException("real-bundle load failed at iter "+i,le); }
						String got=b.metadata("canonical_payload_sha256");
						if(!fRealHash.equals(got))
							throw new AssertionError("real-bundle load hash mismatch: got="+got+" expected="+fRealHash);
						if(b.size()!=176) throw new AssertionError("real-bundle load size mismatch: got="+b.size());
						realSuccesses.incrementAndGet();
					}
				}}));
			}

			Throwable firstFailure=null;
			for(Future<?> f:futures) {
				try { f.get(FUTURE_TIMEOUT_S,TimeUnit.SECONDS); }
				catch(TimeoutException te) { if(firstFailure==null) firstFailure=te; }
				catch(ExecutionException ee) { if(firstFailure==null) firstFailure=ee.getCause(); }
			}
			if(firstFailure!=null) throw new RuntimeException("concurrent load worker failed",firstFailure);

			int expectedEach=THREADS_PER_PATH*ITERS_PER_THREAD;
			check(legacySuccesses.get()==expectedEach,"legacy successes="+legacySuccesses.get()+" expected="+expectedEach);
			check(realSuccesses.get()==expectedEach,"real-bundle successes="+realSuccesses.get()+" expected="+expectedEach);
		} finally {
			// Prompt shutdown on EVERY path (success, assertion failure, or timeout) -- a failed
			// run must not leave non-daemon worker threads alive blocking JVM exit.
			pool.shutdownNow();
		}

		System.out.println("PASS: "+TOTAL_THREADS+" threads ("+THREADS_PER_PATH+" legacy + "+THREADS_PER_PATH+
			" real) x "+ITERS_PER_THREAD+" barrier-synchronized concurrent loadMultiOutput() calls on TWO "+
			"DIFFERENT bundle paths; legacy="+legacySuccesses.get()+" successes, real="+realSuccesses.get()+
			" successes, every canonical_payload_sha256 matched, zero corruption from concurrent parseNet()/"+
			"CellNetParser master construction");
		System.out.println("MAGQC_BUNDLE_CONCURRENT_LOAD PASS");
	}

	static void check(boolean ok,String what){ if(!ok) throw new AssertionError("FAIL: "+what); }
}
