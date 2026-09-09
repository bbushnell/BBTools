package prot;

import ml.CellNet;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/** Explicitly separates and labels cold-copy (first netForCurrentThread() call on a Subnet
 * instance, no cached ThreadLocal entry yet) from warm-copy (cached reuse) construction and
 * scoring -- THREAD_SAFE_BUNDLE_INFERENCE_20260909.md item 1. MagQCNetBundleFourOutputTest and
 * the general concurrency coverage in this package exercise both cold and warm paths implicitly,
 * but never count or assert the split; this makes it explicit in three passes: (A) forced-cold
 * via a fresh MagQCNetBundle.loadMultiOutput() per iteration (guarantees a brand-new Subnet, and
 * therefore a brand-new empty ThreadLocal, every time -- independent of physical thread), (B)
 * forced-warm via one load + repeated same-thread reuse with an explicit identity check, (C) a
 * thread-pool pass that counts and asserts the exact cold/warm split (1 cold + (cycles-1) warm
 * per thread).
 *
 * <p>Pass (C)'s pool is guarded by try/finally (Yoimiya's review, 2026-09-09): the original
 * version called Future.get() then pool.shutdown() as two separate unguarded statements, so a
 * worker's AssertionError propagating out of f.get() skipped shutdown entirely, leaving the
 * pool's non-daemon worker threads alive and the JVM unable to exit -- a failed assertion left
 * workers alive, converting a loud failure into a hang. Fixed with a bounded Future.get() timeout
 * and pool.shutdownNow() in a finally block covering every exit path. */
public final class MagQCNetBundleColdWarmTest {
	static final long FUTURE_TIMEOUT_S=30L;

	public static void main(String[] args) throws Exception {
		String bundlePath=null, rawValPath=null, torchPath=null;
		for(String a:args) {
			int e=a.indexOf('='); if(e<1) throw new IllegalArgumentException("Expected key=value: "+a);
			String k=a.substring(0,e).toLowerCase(), v=a.substring(e+1);
			if(k.equals("bundle")) bundlePath=v;
			else if(k.equals("rawval")) rawValPath=v;
			else if(k.equals("torch")) torchPath=v;
		}
		if(bundlePath==null||rawValPath==null||torchPath==null)
			throw new IllegalArgumentException("Required: bundle= rawval= torch= "+
				"(a mixed bundle with a sparse four-output famset_0 + dense legacy famset_1, and its "+
				"raw-input/Torch-reference TSVs)");

		float[][] rows=readFloatRows(rawValPath,0);
		float[][] torch=readFloatRows(torchPath,1);

		MagQCNetBundle warmupBundle=MagQCNetBundle.loadMultiOutput(Paths.get(bundlePath));
		MagQCNetBundle.Subnet warmupSparse=warmupBundle.subnet("famset_0");
		MagQCNetBundle.Subnet warmupDense=warmupBundle.subnet("famset_1");
		validateFixture(rows,torch,warmupSparse.expectedInputs,4);
		float[] denseInput=new float[warmupDense.expectedInputs];
		for(int i=0;i<denseInput.length;i++) denseInput[i]=(float)((i*3)%7)/10f;
		float expectedDense=warmupDense.score(denseInput.clone());

		// (A) COLD PASS
		final int COLD_ITERS=8;
		int coldCount=0;
		Set<Integer> coldIdentities=new HashSet<Integer>();
		for(int i=0;i<COLD_ITERS;i++){
			MagQCNetBundle b=MagQCNetBundle.loadMultiOutput(Paths.get(bundlePath));
			MagQCNetBundle.Subnet sparse=b.subnet("famset_0");
			MagQCNetBundle.Subnet dense=b.subnet("famset_1");

			CellNet sparseNet=sparse.netForCurrentThread();
			coldIdentities.add(System.identityHashCode(sparseNet));
			float[] got=sparse.scoreAll(rows[0]);
			for(int c=0;c<4;c++) check(Float.isFinite(got[c]) && Math.abs(got[c]-torch[0][c])<=0.001f,
				"COLD sparse parity iter="+i+" col="+c+" got="+got[c]+" expected="+torch[0][c]);
			coldCount++;

			CellNet denseNet=dense.netForCurrentThread();
			coldIdentities.add(System.identityHashCode(denseNet));
			float ds=dense.score(denseInput.clone());
			check(Float.floatToIntBits(ds)==Float.floatToIntBits(expectedDense),
				"COLD dense parity iter="+i+" got="+ds+" expected="+expectedDense);
			coldCount++;
		}
		check(coldIdentities.size()==coldCount,"COLD pass expected "+coldCount+
			" distinct CellNet identities (fresh construction every time), got only "+
			coldIdentities.size()+" distinct -- construction was unexpectedly reused");

		// (B) WARM PASS
		MagQCNetBundle warmBundle=MagQCNetBundle.loadMultiOutput(Paths.get(bundlePath));
		MagQCNetBundle.Subnet warmSparse=warmBundle.subnet("famset_0");
		MagQCNetBundle.Subnet warmDense=warmBundle.subnet("famset_1");
		final int WARM_CYCLES=200;
		CellNet firstSparse=warmSparse.netForCurrentThread();
		CellNet firstDense=warmDense.netForCurrentThread();
		int warmCount=0;
		for(int cycle=0;cycle<WARM_CYCLES;cycle++){
			CellNet sn=warmSparse.netForCurrentThread();
			check(sn==firstSparse,"WARM sparse identity changed at cycle "+cycle+" -- caching broken");
			for(int r=0;r<rows.length;r++){
				float[] got=warmSparse.scoreAll(rows[r]);
				for(int c=0;c<4;c++) check(Float.isFinite(got[c]) && Math.abs(got[c]-torch[r][c])<=0.001f,
					"WARM sparse parity cycle="+cycle+" row="+r+" col="+c);
			}
			warmCount++;

			CellNet dn=warmDense.netForCurrentThread();
			check(dn==firstDense,"WARM dense identity changed at cycle "+cycle+" -- caching broken");
			float ds=warmDense.score(denseInput.clone());
			check(Float.floatToIntBits(ds)==Float.floatToIntBits(expectedDense),
				"WARM dense parity cycle="+cycle);
			warmCount++;
		}

		// (C) THREAD-POOL PASS -- pool guarded by try/finally, Future.get() bounded.
		final int THREADS=4, CYCLES_PER_THREAD=16;
		MagQCNetBundle threadBundle=MagQCNetBundle.loadMultiOutput(Paths.get(bundlePath));
		final MagQCNetBundle.Subnet threadSparse=threadBundle.subnet("famset_0");
		final AtomicInteger threadColdCount=new AtomicInteger(0);
		final AtomicInteger threadWarmCount=new AtomicInteger(0);
		final ConcurrentHashMap<Long,Integer> identityByThread=new ConcurrentHashMap<Long,Integer>();
		final float[] row0=rows[0]; final float[][] torchRef=torch;

		ExecutorService pool=Executors.newFixedThreadPool(THREADS);
		try {
			List<Future<?>> futures=new ArrayList<Future<?>>();
			for(int t=0;t<THREADS;t++){
				futures.add(pool.submit(new Runnable(){ public void run(){
					for(int cycle=0;cycle<CYCLES_PER_THREAD;cycle++){
						CellNet n=threadSparse.netForCurrentThread();
						int idHash=System.identityHashCode(n);
						Integer prior=identityByThread.putIfAbsent(Thread.currentThread().getId(),idHash);
						if(prior==null) threadColdCount.incrementAndGet();
						else {
							if(prior!=idHash) throw new AssertionError("thread "+Thread.currentThread().getId()+
								" identity changed between cycles -- caching broken");
							threadWarmCount.incrementAndGet();
						}
						float[] got=threadSparse.scoreAll(row0);
						for(int c=0;c<4;c++) if(!Float.isFinite(got[c]) || Math.abs(got[c]-torchRef[0][c])>0.001f)
							throw new AssertionError("thread-pool sparse parity cycle="+cycle+" col="+c);
					}
				}}));
			}
			Throwable firstFailure=null;
			for(Future<?> f:futures) {
				try { f.get(FUTURE_TIMEOUT_S,TimeUnit.SECONDS); }
				catch(TimeoutException te) { if(firstFailure==null) firstFailure=te; }
				catch(ExecutionException ee) { if(firstFailure==null) firstFailure=ee.getCause(); }
			}
			if(firstFailure!=null) throw new RuntimeException("thread-pool worker failed",firstFailure);

			check(threadColdCount.get()==THREADS,"expected exactly "+THREADS+" cold events across the pool "+
				"(one per physical thread), got "+threadColdCount.get());
			check(threadWarmCount.get()==THREADS*(CYCLES_PER_THREAD-1),"expected "+
				(THREADS*(CYCLES_PER_THREAD-1))+" warm events, got "+threadWarmCount.get());
		} finally {
			// Prompt shutdown on EVERY path -- a failed assertion inside a worker must not leave
			// non-daemon pool threads alive, blocking JVM exit.
			pool.shutdownNow();
		}

		System.out.println("COLD: "+coldCount+" fresh-bundle-load cold constructions ("+coldIdentities.size()+
			" distinct CellNet identities), all scored correctly");
		System.out.println("WARM: "+warmCount+" cached-reuse calls (single load, single thread), identity "+
			"preserved, all scored correctly");
		System.out.println("THREAD-POOL: "+THREADS+" threads x "+CYCLES_PER_THREAD+" cycles -> "+
			threadColdCount.get()+" cold + "+threadWarmCount.get()+" warm events, identity preserved per "+
			"thread, all scored correctly");
		System.out.println("MAGQC_BUNDLE_COLDWARM PASS");
	}

	static float[][] readFloatRows(String path,int skipCols) throws Exception {
		List<float[]> out=new ArrayList<float[]>();
		for(String line: Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8)) {
			if(line.length()==0||line.startsWith("#")) continue;
			String[] f=line.split("\t");
			float[] v=new float[f.length-skipCols];
			for(int i=0;i<v.length;i++) v[i]=Float.parseFloat(f[i+skipCols]);
			out.add(v);
		}
		return out.toArray(new float[0][]);
	}

	/** Fails loud, before any pass (cold, warm, or thread-pool) runs, on a malformed or degenerate
	 * fixture that would otherwise make a pass a silent no-op (zero rows -> the per-row parity
	 * loop never runs -> vacuous PASS) or let a corrupted reference value evade detection
	 * (Yoimiya's review, 2026-09-09: a NaN Torch value makes `Math.abs(got-NaN)>tol` evaluate to
	 * FALSE in Java -- any comparison with NaN using >,<,>=,<= is false -- so the tolerance check
	 * alone would silently accept a garbage reference row; only checking the SCORED value's
	 * finiteness, as the original parity checks did, misses this). */
	static void validateFixture(float[][] rows,float[][] torch,int expectedInputs,int expectedTorchWidth) {
		check(rows.length>0,"raw-value fixture has ZERO rows -- the per-row parity loop would never "+
			"run, making every pass a vacuous no-op PASS");
		check(rows.length==torch.length,"raw/Torch row count mismatch: raw="+rows.length+" torch="+torch.length);
		for(int r=0;r<rows.length;r++) {
			check(rows[r].length==expectedInputs,"row "+r+" raw width="+rows[r].length+" expected="+expectedInputs);
			check(torch[r].length==expectedTorchWidth,"row "+r+" torch width="+torch[r].length+" expected="+expectedTorchWidth);
			for(int c=0;c<rows[r].length;c++) check(Float.isFinite(rows[r][c]),"row "+r+" raw col "+c+" is not finite: "+rows[r][c]);
			for(int c=0;c<torch[r].length;c++) check(Float.isFinite(torch[r][c]),"row "+r+" torch col "+c+" is not finite: "+torch[r][c]);
		}
	}
	static void check(boolean ok,String what){ if(!ok) throw new AssertionError("FAIL: "+what); }
}
