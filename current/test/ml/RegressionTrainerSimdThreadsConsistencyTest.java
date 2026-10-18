package test.ml;

import ml.RegressionTrainer;

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

/**
 * End-to-end regression check (Brian, 2026-08-31) that runs three configurations
 * on the same seed/data and asserts their bestValid values stay within a 5%
 * relative sanity bound: scalar double with one thread; SIMD/float with one
 * thread; and SIMD/float with a multithreaded worker pool.
 */
public class RegressionTrainerSimdThreadsConsistencyTest{

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

			assertRelativeBestValid("SIMD", "scalar", scalarValid,
				"simd-serial", simdSerialValid);
			assertRelativeBestValid("Multithreading", "simd-serial", simdSerialValid,
				"simd-parallel", simdParallelValid);
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

	/** Runs one trainer configuration and returns its best validation loss. */
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

	/** Asserts a named finite relative bestValid bound. */
	private static void assertRelativeBestValid(String label, String aName, double a,
			String bName, double b){
		final double relDelta=relDelta(a, b);
		if(!(relDelta<=MAX_RELATIVE_DELTA)){
			throw new AssertionError(label+" changed bestValid by more than "
				+(MAX_RELATIVE_DELTA*100)+"%: "+aName+"="+a+" "+bName+"="+b
				+" relDelta="+relDelta);
		}
	}

	/** Returns relative absolute delta with a tiny denominator floor. */
	private static double relDelta(double a, double b){
		return Math.abs(a-b)/Math.max(1e-300, Math.abs(a));
	}

	/** Writes deterministic, learnable sigmoid rows for one fixture input file. */
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

	/** Reads a private RegressionTrainer double field. */
	private static double getDouble(Object owner, String name) throws Exception{
		final Field field=RegressionTrainer.class.getDeclaredField(name);
		field.setAccessible(true);
		return ((Double)field.get(owner)).doubleValue();
	}
}
