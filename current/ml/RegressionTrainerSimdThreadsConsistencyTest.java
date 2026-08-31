package ml;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.Random;
import java.util.stream.Stream;

import shared.Timer;

/** End-to-end regression check (Brian, 2026-08-31): confirms the EXISTING SIMD and multithreading
 * paths do not change training results by more than a small, float-precision/reduction-order
 * amount, on a small net that trains fully in about a minute. Three configs, same seed/data:
 *   (1) simd=f            -- scalar double, single-threaded: the original baseline path.
 *   (2) simd=t threads=1  -- float/SIMD, but gradientPool==null so still single-threaded
 *                            (isolates the SIMD/float-precision effect alone).
 *   (3) simd=t threads=N  -- float/SIMD AND multithreaded (isolates the added threading/
 *                            parallel-reduction-order effect on top of (2)).
 * (1)-vs-(2) isolates SIMD's effect; (2)-vs-(3) isolates multithreading's effect (including this
 * session's new parallel validate(), which only activates when threads>1). */
public class RegressionTrainerSimdThreadsConsistencyTest {

	private static final int INPUTS=64, HIDDEN=32;
	private static final int NUM_TRAIN=1000, NUM_VALID=300;
	private static final int EPOCHS=25, BATCH=32, THREADS=8;
	private static final double MAX_RELATIVE_DELTA=0.05;//5%, a sanity bound, not bit-identity

	public static void main(String[] args) throws Exception{
		final Path dir=Files.createTempDirectory("regressiontrainer-simd-threads-consistency-");
		try{
			final Path train=dir.resolve("train.tsv");
			final Path valid=dir.resolve("valid.tsv");
			final Random gen=new Random(11);
			writeSyntheticRows(train, NUM_TRAIN, gen);
			writeSyntheticRows(valid, NUM_VALID, gen);

			final double scalarValid=runConfig(dir, train, valid, "scalar", false, 1);
			final double simdSerialValid=runConfig(dir, train, valid, "simd-serial", true, 1);
			final double simdParallelValid=runConfig(dir, train, valid, "simd-parallel", true, THREADS);

			final double simdDelta=relDelta(scalarValid, simdSerialValid);
			final double threadDelta=relDelta(simdSerialValid, simdParallelValid);

			System.err.println("RegressionTrainerSimdThreadsConsistencyTest results:");
			System.err.println("  scalar (simd=f, t=1)       bestValid="+scalarValid);
			System.err.println("  simd-serial (simd=t, t=1)  bestValid="+simdSerialValid);
			System.err.println("  simd-parallel (simd=t, t="+THREADS+") bestValid="+simdParallelValid);
			System.err.println("  SIMD effect (scalar vs simd-serial):        relDelta="+simdDelta);
			System.err.println("  threading effect (simd-serial vs parallel): relDelta="+threadDelta);

			if(simdDelta>MAX_RELATIVE_DELTA){
				throw new AssertionError("SIMD changed bestValid by more than "
					+(MAX_RELATIVE_DELTA*100)+"%: scalar="+scalarValid+" simd-serial="+simdSerialValid
					+" relDelta="+simdDelta);
			}
			if(threadDelta>MAX_RELATIVE_DELTA){
				throw new AssertionError("Multithreading changed bestValid by more than "
					+(MAX_RELATIVE_DELTA*100)+"%: simd-serial="+simdSerialValid
					+" simd-parallel="+simdParallelValid+" relDelta="+threadDelta);
			}
			System.err.println("RegressionTrainerSimdThreadsConsistencyTest: PASS");
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	private static double runConfig(Path dir, Path train, Path valid, String label,
			boolean simd, int threads) throws Exception{
		final Path out=dir.resolve(label+".bbnet");
		final Timer t=new Timer();
		final RegressionTrainer trainer=new RegressionTrainer(new String[]{
			"in="+train, "valin="+valid, "out="+out, "dims="+INPUTS+","+HIDDEN+",1",
			"epochs="+EPOCHS, "batch="+BATCH, "seed=5", "final=sigmoid",
			"simd="+(simd ? "t" : "f"), "threads="+threads, "printevery=0"
		});
		trainer.process(new Timer());
		t.stop();
		final double bestValid=getDouble(trainer, "bestValid");
		System.err.println("  ["+label+"] trained "+EPOCHS+" epochs in "+t.timeInSeconds()
			+"s, bestValid="+bestValid);
		return bestValid;
	}

	private static double relDelta(double a, double b){
		return Math.abs(a-b)/Math.max(1e-300, Math.abs(a));
	}

	private static void writeSyntheticRows(Path file, int rows, Random gen) throws IOException{
		try(BufferedWriter w=Files.newBufferedWriter(file, StandardCharsets.US_ASCII)){
			w.write("#dims\t"+INPUTS+"\t1\n");
			final StringBuilder sb=new StringBuilder(INPUTS*4);
			for(int r=0; r<rows; r++){
				sb.setLength(0);
				double sum=0;
				for(int i=0; i<INPUTS; i++){
					final double v=gen.nextInt(1000)/1000.0;
					sum+=v;
					sb.append(v);
					sb.append('\t');
				}
				//A real (if simple) learnable relationship, not pure noise, so training actually
				//converges and bestValid is a meaningful, comparable number across configs.
				final double target=1.0/(1.0+Math.exp(-(sum/INPUTS-0.5)*4));
				sb.append(target);
				sb.append('\n');
				w.write(sb.toString());
			}
		}
	}

	private static double getDouble(Object owner, String name) throws Exception{
		final Field field=RegressionTrainer.class.getDeclaredField(name);
		field.setAccessible(true);
		return ((Double)field.get(owner)).doubleValue();
	}
}
