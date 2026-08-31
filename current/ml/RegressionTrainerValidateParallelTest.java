package ml;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.stream.Stream;

import shared.Timer;

/** Correctness and speed check for the phase-1 validate() parallelization (2026-08-31, closed
 * out 2026-08-31 per Elly's review). Trains a tiny model at production width
 * (dims=9659,2048,512,1) so predict()'s per-sample cost matches the real cdiag jobs exactly, then
 * invokes the private validate() method twice via reflection on the SAME trained weights: once
 * through the parallel worker-pool path, once forced through the serial fallback (by temporarily
 * nulling gradientPool). Since validate() is read-only with respect to model state, both calls
 * see an identical model and any difference in result is attributable only to the reduction
 * path, not to training divergence. Covers BOTH the external-validation branch (valin=) and the
 * internal vfraction split branch (validationData==null, perm[numTrain+i] indexing) -- the
 * pre-2026-08-31 version of this file only exercised the external branch. */
public class RegressionTrainerValidateParallelTest {

	private static final int INPUTS=9659;
	private static final int NUM_TRAIN=64;
	private static final int NUM_VALID=2000;//1/22 of production's 44000, same per-sample cost
	private static final int THREADS=8;

	public static void main(String[] args) throws Exception{
		testExternalValidation();
		testInternalSplit();
		System.err.println("RegressionTrainerValidateParallelTest: PASS");
	}

	/** External validation (valin=): the branch the pre-2026-08-31 version of this test already
	 * covered. Kept as its own method so the internal-split addition doesn't disturb it. */
	private static void testExternalValidation() throws Exception{
		final Path dir=Files.createTempDirectory("regressiontrainer-validate-parallel-ext-");
		try{
			final Path train=dir.resolve("train.tsv");
			final Path valid=dir.resolve("valid.tsv");
			final Random gen=new Random(7);
			writeSyntheticRows(train, NUM_TRAIN, gen);
			writeSyntheticRows(valid, NUM_VALID, gen);
			final Path out=dir.resolve("out.bbnet");

			final RegressionTrainer trainer=new RegressionTrainer(new String[]{
				"in="+train, "valin="+valid, "out="+out, "dims=9659,2048,512,1",
				"epochs=1", "batch=32", "seed=3", "final=sigmoid", "simd=t",
				"threads="+THREADS, "printevery=1"
			});
			trainer.process(new Timer());
			compareParallelVsSerial(trainer, "external (valin=)", NUM_VALID);
		}finally{
			deleteTree(dir);
		}
	}

	/** Internal vfraction split (Elly's review, 2026-08-31): the same parallel-vs-serial gate,
	 * but with NO valin= -- validationData stays null and validate() must use the
	 * perm[numTrain+i] indexing branch instead. This branch was previously only exercised
	 * serially (by the pre-existing RegressionTrainerRowTest's testLegacyFraction, threads=1). */
	private static void testInternalSplit() throws Exception{
		final Path dir=Files.createTempDirectory("regressiontrainer-validate-parallel-int-");
		try{
			final Path train=dir.resolve("train.tsv");
			final Random gen=new Random(13);
			//vfraction carves the internal split out of one file, so it needs train+valid rows
			//together; NUM_TRAIN covers the training share, NUM_VALID/(1-vfraction)-ish covers
			//the rest -- use a plain total large enough to leave >=NUM_VALID rows after the split.
			final int total=NUM_TRAIN+NUM_VALID;
			writeSyntheticRows(train, total, gen);
			final Path out=dir.resolve("out.bbnet");

			final double vfraction=NUM_VALID/(double)total;
			final RegressionTrainer trainer=new RegressionTrainer(new String[]{
				"in="+train, "out="+out, "dims=9659,2048,512,1", "vfraction="+vfraction,
				"epochs=1", "batch=32", "seed=3", "final=sigmoid", "simd=t",
				"threads="+THREADS, "printevery=1"
			});
			trainer.process(new Timer());
			final int numValid=getInt(trainer, "numValid");
			compareParallelVsSerial(trainer, "internal (vfraction="+vfraction+")", numValid);
		}finally{
			deleteTree(dir);
		}
	}

	/** Shared parallel-vs-serial comparison: rebuilds the pool the trainer shut down at the end
	 * of process(), runs validate() through both paths on the identical final model, asserts
	 * numerical equivalence, and reports timing as EVIDENCE rather than a hard pass/fail -- a
	 * relative wall-clock assertion can flake under machine load and isn't a correctness signal
	 * (Elly's review, 2026-08-31). */
	private static void compareParallelVsSerial(RegressionTrainer trainer, String label,
			int numValid) throws Exception{
		final Method setup=RegressionTrainer.class.getDeclaredMethod("setupGradientWorkers");
		setup.setAccessible(true);
		setup.invoke(trainer);

		final Object pool=getField(trainer, "gradientPool");
		final int gradientThreadCount=getInt(trainer, "gradientThreadCount");
		if(pool==null || gradientThreadCount<=1){
			throw new AssertionError("["+label+"] Test requires the parallel path to have been "
				+"set up (gradientPool=null or gradientThreadCount<=1); check threads="+THREADS
				+" was actually honored.");
		}

		final Timer tPar=new Timer();
		final double mseParallel=invokeValidate(trainer);
		tPar.stop();

		setField(trainer, "gradientPool", null);//force the serial fallback path
		final Timer tSer=new Timer();
		final double mseSerial=invokeValidate(trainer);
		tSer.stop();
		setField(trainer, "gradientPool", pool);//restore, in case anything reads it later

		final double relDelta=Math.abs(mseParallel-mseSerial)/Math.max(1e-300, Math.abs(mseSerial));
		final double speedup=tSer.timeInSeconds()/Math.max(1e-9, tPar.timeInSeconds());
		System.err.println("RegressionTrainerValidateParallelTest ["+label+"]: threads="+THREADS
			+" numValid="+numValid+" inputs="+INPUTS);
		System.err.println("  serial   mse="+mseSerial+"  time="+tSer.timeInSeconds()+"s");
		System.err.println("  parallel mse="+mseParallel+"  time="+tPar.timeInSeconds()+"s");
		System.err.println("  relative MSE delta="+relDelta+"  speedup="+speedup
			+"x (evidence only -- NOT asserted; wall-clock ratios can flake under machine load)");

		if(relDelta>1e-9){
			throw new AssertionError("["+label+"] Parallel and serial validate() disagree beyond "
				+"floating-point reordering tolerance: serial="+mseSerial+" parallel="+mseParallel
				+" relDelta="+relDelta);
		}

		shutdownIfPresent(trainer);
	}

	private static void writeSyntheticRows(Path file, int rows, Random gen) throws IOException{
		try(BufferedWriter w=Files.newBufferedWriter(file, StandardCharsets.US_ASCII)){
			w.write("#dims\t"+INPUTS+"\t1\n");
			final StringBuilder sb=new StringBuilder(INPUTS*4);
			for(int r=0; r<rows; r++){
				sb.setLength(0);
				for(int i=0; i<INPUTS; i++){
					sb.append(gen.nextInt(1000)/1000.0);
					sb.append('\t');
				}
				sb.append(gen.nextInt(2));//single binary-ish output, sigmoid-compatible
				sb.append('\n');
				w.write(sb.toString());
			}
		}
	}

	private static double invokeValidate(Object trainer) throws Exception{
		final Method m=RegressionTrainer.class.getDeclaredMethod("validate");
		m.setAccessible(true);
		return (Double)m.invoke(trainer);
	}

	private static void shutdownIfPresent(Object trainer) throws Exception{
		final Object pool=getField(trainer, "gradientPool");
		if(pool instanceof ExecutorService){((ExecutorService)pool).shutdownNow();}
	}

	private static void deleteTree(Path dir) throws IOException{
		try(Stream<Path> paths=Files.walk(dir)){
			paths.sorted(Comparator.reverseOrder()).forEach(path -> {
				try{Files.deleteIfExists(path);}
				catch(Exception e){throw new RuntimeException(e);}
			});
		}
	}

	private static Object getField(Object owner, String name) throws Exception{
		final Field field=RegressionTrainer.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.get(owner);
	}

	private static void setField(Object owner, String name, Object value) throws Exception{
		final Field field=RegressionTrainer.class.getDeclaredField(name);
		field.setAccessible(true);
		field.set(owner, value);
	}

	private static int getInt(Object owner, String name) throws Exception{
		return ((Integer)getField(owner, name)).intValue();
	}
}
