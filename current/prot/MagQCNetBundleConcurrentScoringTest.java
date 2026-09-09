package prot;

import ml.CellNet;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;

/** Genuine mixed-format concurrency test through the ACTUAL MagQCNetBundle.Subnet API (not a
 * standalone CellNet subclass probe): both direct netForCurrentThread().feedForward() calls (the
 * path MagQCVectorMaker.formatAggRowCore() already uses) and the score()/scoreAll() wrappers, for
 * a sparse (real candidate) and a dense (legacy) subnet, scored CONCURRENTLY from multiple worker
 * threads, while a separate thread continuously toggles the shared, deprecated-for-inference
 * ml.CellNet.DENSE flag -- proving the InferenceNet dispatch fix (per-instance denseMode, no
 * shared mutable state at inference time) is genuinely race-free, not merely correct when scored
 * serially. */
public final class MagQCNetBundleConcurrentScoringTest {
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

		MagQCNetBundle b=MagQCNetBundle.loadMultiOutput(Paths.get(bundlePath));
		MagQCNetBundle.Subnet sparse=b.subnet("famset_0");
		MagQCNetBundle.Subnet dense=b.subnet("famset_1");

		float[][] rows=readFloatRows(rawValPath,0);
		float[][] torch=readFloatRows(torchPath,1);
		validateFixture(rows,torch,sparse.expectedInputs,4);
		float[] denseInput=new float[dense.expectedInputs];
		for(int i=0;i<denseInput.length;i++) denseInput[i]=(float)((i*3)%7)/10f;
		float expectedDense=dense.score(denseInput.clone());

		AtomicBoolean running=new AtomicBoolean(true);
		Thread toggler=new Thread(()->{
			while(running.get()){ CellNet.DENSE=true; Thread.yield(); CellNet.DENSE=false; Thread.yield(); }
		});
		ExecutorService pool=Executors.newFixedThreadPool(4);
		ArrayList<Future<?>> work=new ArrayList<>();
		toggler.start();
		final int JOBS=64, CYCLES=30;
		try {
			for(int job=0;job<JOBS;job++){
				work.add(pool.submit(()->{
					for(int cycle=0;cycle<CYCLES;cycle++){
						for(int r=0;r<rows.length;r++){
							float[] got=sparse.scoreAll(rows[r]);
							for(int i=0;i<4;i++){
								if(!Float.isFinite(got[i]) || Math.abs(got[i]-torch[r][i])>0.001f){
									throw new AssertionError("sparse scoreAll parity row="+r+" col="+i+
										" got="+got[i]+" expected="+torch[r][i]);
								}
							}
						}
						float ds=dense.score(denseInput.clone());
						if(Float.floatToIntBits(ds)!=Float.floatToIntBits(expectedDense)){
							throw new AssertionError("dense score() parity: got="+ds+" expected="+expectedDense);
						}
						CellNet sn=sparse.netForCurrentThread();
						for(int r=0;r<rows.length;r++){
							sn.applyInput(rows[r]); sn.feedForward();
							float[] got=sn.getOutput();
							for(int i=0;i<4;i++){
								if(!Float.isFinite(got[i]) || Math.abs(got[i]-torch[r][i])>0.001f){
									throw new AssertionError("sparse direct-feedForward parity row="+r+" col="+i);
								}
							}
						}
						CellNet dn=dense.netForCurrentThread();
						dn.applyInput(denseInput.clone()); dn.feedForward();
						if(Float.floatToIntBits(dn.getOutput(0))!=Float.floatToIntBits(expectedDense)){
							throw new AssertionError("dense direct-feedForward parity");
						}
					}
				}));
			}
			for(Future<?> f:work) f.get();
		} finally {
			running.set(false); toggler.join(); pool.shutdownNow();
		}
		System.out.println("PASS: "+JOBS+" tasks x "+CYCLES+" cycles, 4 workers, external DENSE toggler; "+
			"score()/scoreAll() wrappers AND direct netForCurrentThread().feedForward() both correct "+
			"for concurrently-scored sparse+dense subnets");
		System.out.println("MAGQC_BUNDLE_CONCURRENT_SCORING PASS");
	}

	static float[][] readFloatRows(String path,int skipCols) throws Exception {
		List<float[]> out=new ArrayList<>();
		for(String line: Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8)) {
			if(line.length()==0||line.startsWith("#")) continue;
			String[] f=line.split("\t");
			float[] v=new float[f.length-skipCols];
			for(int i=0;i<v.length;i++) v[i]=Float.parseFloat(f[i+skipCols]);
			out.add(v);
		}
		return out.toArray(new float[0][]);
	}

	/** Fails loud, before any worker is launched, on a malformed or degenerate fixture that would
	 * otherwise make the test a silent no-op (zero rows -> the scoring loop never runs -> vacuous
	 * PASS) or let a corrupted reference value evade detection (Yoimiya's review, 2026-09-09: a
	 * NaN Torch value makes `Math.abs(got-NaN)>tol` evaluate to FALSE in Java -- any comparison
	 * with NaN using >,<,>=,<= is false -- so the tolerance check alone would silently accept a
	 * garbage reference row; only checking the SCORED value's finiteness, as the original parity
	 * check did, misses this). */
	static void validateFixture(float[][] rows,float[][] torch,int expectedInputs,int expectedTorchWidth) {
		check(rows.length>0,"raw-value fixture has ZERO rows -- the scoring loop would never run, "+
			"making this test a vacuous no-op PASS");
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
