package test.icecream;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.atomic.AtomicLong;

import icecream.IceCreamMaker;
import shared.Random;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/**
 * Finite literal oracles for insertion, deletion, substitution and clean calls.
 * Uses a scripted RNG to select each error branch, without a statistical test.
 * The baseline mode records the existing substitution-only defect explicitly.
 * @author Brian Bushnell, Shinobu
 */
public final class IceCreamMakerErrorTest{

	/** Runs all cases in one fresh JVM; no streams or worker threads are started.
	 * @param args Exactly expect=baseline or expect=fixed
	 * @throws Exception if the existing private test seam cannot be invoked */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !(args[0].equals("expect=baseline") || args[0].equals("expect=fixed"))){
			throw new IllegalArgumentException("Use expect=baseline or expect=fixed");
		}
		final boolean baseline=args[0].equals("expect=baseline");
		final IceCreamMaker maker=new IceCreamMaker(new String[]{"zmws=1", "threads=1",
				"genomesize=100", "len=10", "irlen=2", "seed=1"});
		final Class<?> workerClass=Class.forName("icecream.IceCreamMaker$ProcessThread");
		final Constructor<?> constructor=workerClass.getDeclaredConstructor(IceCreamMaker.class,
				ConcurrentReadOutputStream.class, int.class, AtomicLong.class, long.class);
		constructor.setAccessible(true);
		final Object worker=constructor.newInstance(maker, null, 0, new AtomicLong(), 1L);
		final Field random=workerClass.getDeclaredField("randy");
		random.setAccessible(true);
		final Method call=workerClass.getDeclaredMethod("baseCallOnePass", byte[].class,
				float.class, int.class, int.class, int.class, long.class);
		call.setAccessible(true);
		System.out.println("case\texpected\tactual");
		check("insertion", baseline ? "CC" : "GAC", 0.5f, new float[]{0.1f, 0.9f, 0.9f}, worker, random, call);
		check("deletion", baseline ? "CC" : "C", 0.5f, new float[]{0.3f, 0.9f}, worker, random, call);
		check("substitution", "CC", 0.5f, new float[]{0.45f, 0.9f}, worker, random, call);
		check("clean", "AC", 0.5f, new float[]{0.9f, 0.9f}, worker, random, call);
		check("zero-rate", "AC", 0f, new float[]{0.1f, 0.3f}, worker, random, call);
		check("unit-rate", baseline ? "CG" : "GG", 1f, new float[]{0.2f, 0.6f, 0.9f}, worker, random, call);
	}

	/** Compares emitted bases with an independently specified short sequence. */
	private static void check(String name, String expected, float rate, float[] draws,
			Object worker, Field random, Method call) throws Exception{
		random.set(worker, new ScriptedRandom(draws));
		final Object builder=call.invoke(worker, "AC".getBytes(StandardCharsets.US_ASCII), rate, 0, 0, 10, 7L);
		final Method toRead=builder.getClass().getDeclaredMethod("toRead");
		toRead.setAccessible(true);
		final Read read=(Read)toRead.invoke(builder);
		final String actual=new String(read.bases, StandardCharsets.US_ASCII);
		System.out.println(name+"\t"+expected+"\t"+actual);
		if(!actual.equals(expected)){
			throw new AssertionError(name+": literal error-branch oracle expected "+expected+", found "+actual);
		}
	}

	/** Supplies branch draws; insertion emits G and substitution advances one base. */
	private static final class ScriptedRandom implements Random{
		ScriptedRandom(float[] draws_){draws=draws_;}
		@Override
		public long nextLong(){throw new AssertionError("Unexpected RNG method in finite base-calling fixture");}
		@Override
		public long nextLong(long bound){throw new AssertionError("Unexpected bounded long draw in base-calling fixture");}
		@Override
		public double nextGaussian(){throw new AssertionError("Unexpected Gaussian draw in base-calling fixture");}
		@Override
		public void setSeed(long seed){throw new AssertionError("Finite base-calling draws cannot be reseeded");}
		@Override
		public float nextFloat(){return index<draws.length ? draws[index++] : 0.9f;}
		@Override
		public int nextInt(){return 2;}
		@Override
		public int nextInt3(){return 0;}
		private final float[] draws;
		private int index;
	}
}
