package test.icecream;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import icecream.IceCreamMaker;
import shared.Random;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/** Literal simulator checks for the last fragment start and disabled errors.
 * No worker threads or streams are started.
 * @author Brian Bushnell, Shinobu
 */
public final class IceCreamMakerSamplingTest{

	/** Runs independent finite cases and reports every mismatch before failing.
	 * @param args Optional scope=feasibility for finite rejection probes or scope=preflight for request validation
	 * @throws Exception If the private simulator test seam changes
	 */
	public static void main(String[] args) throws Exception{
		if(args.length==1 && args[0].equals("scope=feasibility")){
			probeFeasibility();
			return;
		}
		if(args.length==1 && args[0].equals("scope=preflight")){
			checkPreflight();
			return;
		}
		if(args.length!=0){throw new IllegalArgumentException("Use no arguments, scope=feasibility or scope=preflight");}
		int failures=0;
		System.out.println("case\texpected\tobserved\tresult");
		failures+=fragment("last-start", "NNAC", "AC");
		failures+=fragment("equal-length", "AC", "AC");
		failures+=fragment("short-reference", "A", "A");
		failures+=repeat("repeat-exact-fit", "ACAA", "ACGT");
		failures+=repeat("repeat-last-start", "AAAAAA", "AAAATT");
		failures+=errors(false, "AC");
		failures+=errors(true, "CG");
		if(failures!=0){throw new AssertionError(failures+" simulator observations differ from literal expectations");}
	}

	/** Literal impossible cases and neighboring allowed bounds, without starting streams. */
	private static void checkPreflight() throws Exception{
		System.out.println("case\texpected\tobserved\tresult");
		int failures=0;
		failures+=preflight("empty", "", false, false, 2, 0, false);
		failures+=preflight("all-N", "NN", false, false, 2, 0, false);
		failures+=preflight("separated-bases", "ANA", false, false, 2, 0, false);
		failures+=preflight("last-defined-window", "NNAC", false, false, 2, 0, true);
		failures+=preflight("short-reference", "A", false, false, 2, 0, true);
		failures+=preflight("ccs-no-deletions-below", "AC", true, false, 93, 0, false);
		failures+=preflight("ccs-no-deletions-bound", "AC", true, false, 94, 0, true);
		failures+=preflight("ccs-deletions-below", "AC", true, true, 90, 0, false);
		failures+=preflight("ccs-deletions-bound", "AC", true, true, 91, 0, true);
		failures+=preflight("missing-no-deletions-below", "AC", true, false, 140, 1, false);
		failures+=preflight("missing-no-deletions-bound", "AC", true, false, 141, 1, true);
		failures+=preflight("missing-near-one", "AC", true, true, 91, 1-0x1.0p-25, false);
		failures+=preflight("missing-highest-float", "AC", true, true, 91, 1-0x1.0p-24, true);
		failures+=repeatRejected("repeat-all-N", "NNNN");
		failures+=repeatRejected("repeat-too-short", "AC");
		failures+=repeat("repeat-exact-fit", "ACAA", "ACGT");
		if(failures!=0){throw new AssertionError(failures+" simulator preflight checks differed");}
	}

	/** Invokes only the production preflight; an allowed result is not a sampling-time guarantee. */
	private static int preflight(String name, String source, boolean ccs, boolean errors,
			int movie, double missing, boolean allowed) throws Exception{
		final IceCreamMaker maker=new IceCreamMaker(new String[]{"zmws=1", "threads=1",
				"genomesize=100", "len=2", "irlen=2", "movie="+movie, "missing="+missing,
				"ccs="+ccs, "adderrors="+errors});
		final Field ref=IceCreamMaker.class.getDeclaredField("ref");
		ref.setAccessible(true);
		ref.set(maker, bytes(source));
		final Method validate=IceCreamMaker.class.getDeclaredMethod("validateSimulation");
		validate.setAccessible(true);
		String observed="allowed";
		try{validate.invoke(maker);}
		catch(InvocationTargetException e){
			if(!(e.getCause() instanceof IllegalArgumentException)){throw e;}
			observed="rejected";
		}
		return observe(name, allowed ? "allowed" : "rejected", observed);
	}

	/** Caps candidate draws so a missing repeat guard fails instead of hanging this check. */
	private static int repeatRejected(String name, String source) throws Exception{
		final IceCreamMaker maker=new IceCreamMaker(new String[]{"zmws=1", "threads=1",
				"genomesize=100", "len=2", "irlen=2", "irrate=0.5", "seed=1"});
		final Method add=IceCreamMaker.class.getDeclaredMethod("addInvertedRepeats", byte[].class, Random.class);
		add.setAccessible(true);
		String actual="allowed";
		try{add.invoke(maker, bytes(source), new EndpointRandom(false, 16));}
		catch(InvocationTargetException e){
			actual=e.getCause() instanceof IllegalArgumentException ? "rejected" :
					"throws="+e.getCause().getClass().getSimpleName();
		}
		return observe(name, "rejected", actual);
	}

	/** Observes one candidate at a time; never calls the unbounded generateList retry loop. */
	private static void probeFeasibility() throws Exception{
		System.out.println("case\texpected\tobserved\tresult");
		int failures=0;
		failures+=probe("empty-reference", "", false, "null");
		failures+=probe("undefined-reference", "NN", false, "null");
		failures+=probe("short-ccs-movie", "AC", true, "null");
		failures+=probe("normal-request", "AC", false, "CG");
		if(failures!=0){throw new AssertionError(failures+" finite candidate-rejection probes differed");}
	}

	/** Empty/all-N references and a two-unit CCS movie cannot yield a usable candidate. */
	@SuppressWarnings("unchecked")
	private static int probe(String name, String source, boolean ccs, String expected) throws Exception{
		final Object worker=worker(true, source, false, ccs);
		final Method generate=worker.getClass().getDeclaredMethod("generateZMW", long.class);
		generate.setAccessible(true);
		final ArrayList<Read> reads=(ArrayList<Read>)generate.invoke(worker, 7L);
		final String actual=reads==null ? "null" : reads.size()!=1 ? "count="+reads.size() :
				new String(reads.get(0).bases, StandardCharsets.US_ASCII);
		return observe(name, expected, actual);
	}

	/** Places one inverted repeat at the last fitting start; all candidate bases are defined. */
	private static int repeat(String name, String source, String expected) throws Exception{
		final IceCreamMaker maker=new IceCreamMaker(new String[]{"zmws=1", "threads=1",
				"genomesize=100", "len=2", "irlen=2", "irrate=0.5", "seed=1"});
		final Method add=IceCreamMaker.class.getDeclaredMethod("addInvertedRepeats", byte[].class, Random.class);
		add.setAccessible(true);
		final byte[] bases=bytes(source);
		try{
			add.invoke(maker, bases, new EndpointRandom(true));
		}catch(InvocationTargetException e){
			return observe(name, expected, "throws="+e.getCause().getClass().getSimpleName());
		}
		return observe(name, expected, new String(bases, StandardCharsets.US_ASCII));
	}

	/** Selects the last permitted RNG position; only the final two bases are defined. */
	private static int fragment(String name, String source, String expected) throws Exception{
		final Object worker=worker(true, source, true);
		final Method fetch=worker.getClass().getDeclaredMethod("fetchBases", byte[].class);
		fetch.setAccessible(true);
		final byte[] result=(byte[])fetch.invoke(worker, bytes(source));
		return observe(name, expected, result==null ? "null" : new String(result, StandardCharsets.US_ASCII));
	}

	/** er=1 with substitution-only scripted draws advances A to C and C to G. */
	@SuppressWarnings("unchecked")
	private static int errors(boolean enabled, String expected) throws Exception{
		final Object worker=worker(enabled, "AC", false);
		final Method generate=worker.getClass().getDeclaredMethod("generateZMW", long.class);
		generate.setAccessible(true);
		final ArrayList<Read> reads=(ArrayList<Read>)generate.invoke(worker, 7L);
		final String actual=reads==null ? "null" : reads.size()!=1 ? "count="+reads.size() :
				new String(reads.get(0).bases, StandardCharsets.US_ASCII);
		return observe("add-errors-"+enabled, expected, actual);
	}

	/** Builds an unstarted production worker with a tiny reference and deterministic RNG. */
	private static Object worker(boolean errors, String source, boolean last) throws Exception{
		return worker(errors, source, last, false);
	}

	/** Also selects CCS mode for the finite impossible-request probe. */
	private static Object worker(boolean errors, String source, boolean last, boolean ccs) throws Exception{
		final IceCreamMaker maker=new IceCreamMaker(new String[]{"zmws=1", "threads=1",
				"genomesize=100", "len=2", "irlen=2", "movie=2", "er=1", "seed=1",
				"missing=0", "hidden=0", "adderrors="+errors, "ccs="+ccs});
		final Field ref=IceCreamMaker.class.getDeclaredField("ref");
		ref.setAccessible(true);
		ref.set(maker, bytes(source));
		final Class<?> type=Class.forName("icecream.IceCreamMaker$ProcessThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(IceCreamMaker.class,
				ConcurrentReadOutputStream.class, int.class, AtomicLong.class, long.class);
		constructor.setAccessible(true);
		final Object worker=constructor.newInstance(maker, null, 0, new AtomicLong(), 1L);
		final Field random=type.getDeclaredField("randy");
		random.setAccessible(true);
		random.set(worker, new EndpointRandom(last));
		return worker;
	}

	/** Converts the literal ASCII DNA strings without platform-dependent encoding. */
	private static byte[] bytes(String sequence){return sequence.getBytes(StandardCharsets.US_ASCII);}

	/** Returns one on a mismatch while retaining the observation for the baseline. */
	private static int observe(String name, String expected, String actual){
		final boolean pass=expected.equals(actual);
		System.out.println(name+"\t"+expected+"\t"+actual+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Supplies endpoint choices and substitution-event draws, with no real randomness. */
	private static final class EndpointRandom implements Random{
		EndpointRandom(boolean last_){this(last_, Integer.MAX_VALUE);}
		EndpointRandom(boolean last_, int draws){last=last_; drawsLeft=draws;}
		@Override
		public int nextInt(int bound){
			if(--drawsLeft<0){throw new AssertionError("Candidate draw limit reached; expected preflight rejection");}
			if(bound<=0){throw new AssertionError("Simulator requested a nonpositive RNG bound: "+bound);}
			return last ? bound-1 : 0;
		}
		@Override
		public int nextInt(){return 0;}
		@Override
		public int nextInt3(){return 0;}
		@Override
		public boolean nextBoolean(){return false;}
		@Override
		public float nextFloat(){return 0.9f;}
		@Override
		public long nextLong(){throw new AssertionError("Unexpected long draw in finite simulator fixture");}
		@Override
		public long nextLong(long bound){throw new AssertionError("Unexpected bounded long draw in finite simulator fixture");}
		@Override
		public double nextGaussian(){throw new AssertionError("Unexpected Gaussian draw in finite simulator fixture");}
		@Override
		public void setSeed(long seed){throw new AssertionError("Finite simulator RNG cannot be reseeded");}
		private final boolean last;
		private int drawsLeft;
	}

	/** Static-only fixture. */
	private IceCreamMakerSamplingTest(){}
}
