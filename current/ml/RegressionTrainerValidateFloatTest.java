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
import java.util.stream.Stream;

import shared.Timer;

/** Phase-2 gate (Elly's spec, 2026-08-31): the float/SIMD validate() path (vsimd=t) must be
 * evaluated independently of phase 1's threading change, against: max/mean per-prediction delta,
 * aggregate MSE delta, deterministic repeatability, and speed vs. phase 1 (parallel double). Uses
 * production width (dims=9659,2048,512,1) so per-sample cost matches the real cdiag jobs. */
public class RegressionTrainerValidateFloatTest {

	private static final int INPUTS=9659;
	private static final int NUM_TRAIN=64;
	private static final int NUM_VALID=2000;
	private static final int THREADS=8;

	public static void main(String[] args) throws Exception{
		final Path dir=Files.createTempDirectory("regressiontrainer-validate-float-test-");
		try{
			final Path train=dir.resolve("train.tsv");
			final Path valid=dir.resolve("valid.tsv");
			final Random gen=new Random(7);
			writeSyntheticRows(train, NUM_TRAIN, gen);
			writeSyntheticRows(valid, NUM_VALID, gen);
			final Path out=dir.resolve("out.bbnet");

			final RegressionTrainer trainer=new RegressionTrainer(new String[]{
				"in="+train, "valin="+valid, "out="+out, "dims=9659,2048,512,1",
				"epochs=1", "batch=32", "seed=3", "final=sigmoid", "simd=t", "vsimd=t",
				"threads="+THREADS, "printevery=1"
			});
			trainer.process(new Timer());

			//Rebuild the pool against the final trained state (process() shuts it down).
			final Method setup=RegressionTrainer.class.getDeclaredMethod("setupGradientWorkers");
			setup.setAccessible(true);
			setup.invoke(trainer);

			final Method validate=RegressionTrainer.class.getDeclaredMethod("validate");
			validate.setAccessible(true);

			//Phase-1 (parallel double) timing/MSE reference.
			setField(trainer, "vsimd", false);
			final Timer tDouble=new Timer();
			final double msePhase1=(Double)validate.invoke(trainer);
			tDouble.stop();

			//Phase-2 (parallel float/SIMD) timing/MSE, run twice for determinism.
			setField(trainer, "vsimd", true);
			final Timer tFloatA=new Timer();
			final double mseFloatA=(Double)validate.invoke(trainer);
			tFloatA.stop();
			final Timer tFloatB=new Timer();
			final double mseFloatB=(Double)validate.invoke(trainer);
			tFloatB.stop();

			//Per-prediction max/mean delta between the double and float forward passes on the
			//SAME trained model, over every validation sample.
			final Method predict=RegressionTrainer.class.getDeclaredMethod("predict",
				float[].class, double[][].class, double[][].class, double[].class, double[].class);
			predict.setAccessible(true);
			final Method predictF=RegressionTrainer.class.getDeclaredMethod("predictF",
				float[].class, float[][].class);
			predictF.setAccessible(true);
			final Object weights=getField(trainer, "weights");
			final Object bias=getField(trainer, "bias");
			final Object validationData=getField(trainer, "validationData");
			final Method sizeM=validationData.getClass().getDeclaredMethod("size");
			sizeM.setAccessible(true);
			final int n=(Integer)sizeM.invoke(validationData);
			final Field inputsField=validationData.getClass().getDeclaredField("inputs");
			inputsField.setAccessible(true);
			final float[][] inputs=(float[][])inputsField.get(validationData);
			final double[] scratchA=new double[10000], scratchB=new double[10000];
			final int[] dims=(int[])getField(trainer, "dims");
			final float[][] actBufF=new float[dims.length][];
			for(int i=0; i<dims.length; i++){actBufF[i]=new float[dims[i]];}

			double maxDelta=0, sumDelta=0, sumSqPredDouble=0, sumSqPredDelta=0;
			for(int i=0; i<n; i++){
				final double[] pd=(double[])predict.invoke(trainer, inputs[i], weights, bias,
					scratchA, scratchB);
				final float[] pf=(float[])predictF.invoke(trainer, inputs[i], (Object)actBufF);
				final double delta=Math.abs(pd[0]-pf[0]);
				maxDelta=Math.max(maxDelta, delta);
				sumDelta+=delta;
				sumSqPredDouble+=pd[0]*pd[0];
				sumSqPredDelta+=delta*delta;
			}
			final double meanDelta=sumDelta/n;
			//R^2-style relative magnitude: how big is the float-vs-double prediction delta
			//compared to the double predictions' own scale.
			final double relPredMagnitude=Math.sqrt(sumSqPredDelta/Math.max(1e-300, sumSqPredDouble));

			final double mseRelDelta=Math.abs(msePhase1-mseFloatA)/Math.max(1e-300, Math.abs(msePhase1));
			final double determinismDelta=Math.abs(mseFloatA-mseFloatB);

			System.err.println("RegressionTrainerValidateFloatTest results (threads="+THREADS
				+" numValid="+NUM_VALID+" inputs="+INPUTS+"):");
			System.err.println("  Phase 1 (parallel double): mse="+msePhase1
				+"  time="+tDouble.timeInSeconds()+"s");
			System.err.println("  Phase 2 (parallel float,  run A): mse="+mseFloatA
				+"  time="+tFloatA.timeInSeconds()+"s");
			System.err.println("  Phase 2 (parallel float,  run B): mse="+mseFloatB
				+"  time="+tFloatB.timeInSeconds()+"s");
			System.err.println("  MSE relative delta (phase1 vs phase2): "+mseRelDelta);
			System.err.println("  Determinism: |run A - run B| = "+determinismDelta
				+" (expect exactly 0.0 -- no randomness in the float forward pass or reduction)");
			System.err.println("  Per-prediction max delta="+maxDelta+"  mean delta="+meanDelta
				+"  relative magnitude (RMS delta / RMS double-pred)="+relPredMagnitude);
			System.err.println("  Speed: phase1="+tDouble.timeInSeconds()+"s phase2="
				+tFloatA.timeInSeconds()+"s  additional speedup="
				+(tDouble.timeInSeconds()/Math.max(1e-9, tFloatA.timeInSeconds()))+"x");

			if(determinismDelta!=0.0){
				throw new AssertionError("Float validate() is not deterministic across repeated "
					+"calls on the same model: run A="+mseFloatA+" run B="+mseFloatB);
			}
			if(mseRelDelta>0.05){
				throw new AssertionError("Float validate() MSE differs from double validate() by "
					+"more than 5%: phase1="+msePhase1+" phase2="+mseFloatA
					+" relDelta="+mseRelDelta);
			}
			System.err.println("RegressionTrainerValidateFloatTest: PASS (evidence-only -- "
				+"whether this magnitude of float/double difference is ACCEPTABLE for production "
				+"is a team decision, not asserted here as pass/fail beyond the 5% sanity bound)");

			shutdownIfPresent(trainer);
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
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
				sb.append(gen.nextInt(2));
				sb.append('\n');
				w.write(sb.toString());
			}
		}
	}

	private static void shutdownIfPresent(Object trainer) throws Exception{
		final Object pool=getField(trainer, "gradientPool");
		if(pool instanceof java.util.concurrent.ExecutorService){
			((java.util.concurrent.ExecutorService)pool).shutdownNow();
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
}
