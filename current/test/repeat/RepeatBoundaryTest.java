package test.repeat;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import repeat.Palindrome;
import repeat.PalindromeFinder;
import repeat.Repeat;
import repeat.RepeatFinder;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/** Finite repeat boundaries with literal expectations and independent nearby controls.
 * Private per-read methods are invoked on unstarted workers; no input/output
 * streams or processing threads are started. A tiny FASTA supplies only the
 * RepeatFinder constructor's file-validation requirement. This is not a CLI
 * throughput or worker-lifecycle test.
 * @author Sucrose
 */
public final class RepeatBoundaryTest{

	/** Requires one work=directory argument; preserves its newly created fixture. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("work=")){
			throw new IllegalArgumentException("Expected work=directory");
		}
		final Path work=Paths.get(args[0].substring(5));
		Files.createDirectories(work);
		final Path input=work.resolve("repeat-boundary-input.fa");
		require(!Files.exists(input), "Refusing to replace prior fixture: "+input);
		Files.write(input, ">fixture\nACGTACGT\n".getBytes(StandardCharsets.US_ASCII));
		final RepeatBoundaryTest test=new RepeatBoundaryTest();
		test.check("first-window depth peak", () -> depthPeak(true));
		test.check("later-window depth control", () -> depthPeak(false));
		test.check("terminal entropy count", () -> {
			final Object worker=worker(input, true);
			final int count=invokeCount(worker, "processReadEntropy", new Class<?>[] {Read.class}, read("AAAAAAAA"));
			require(count==8, "All eight low-entropy bases must be counted, actual="+count);
			requireSingleRange(worker, 0, 7);
		});
		test.check("high-entropy control", () -> {
			final Object worker=worker(input, true);
			final int count=invokeCount(worker, "processReadEntropy", new Class<?>[] {Read.class}, read("ACGTACGT"));
			require(count==0 && closedRepeats(worker).isEmpty(), "Balanced four-base windows must not be low entropy");
		});
		test.check("exact two-copy tandem", () -> tandem(input, "ATAT"));
		//ATAT at 0..3 and TATA at 1..4 overlap; their literal union is all five bases.
		test.check("overlapping tandem union", () -> tandem(input, "ATATA", 5));
		test.check("tandem sentinel control", () -> tandem(input, "ATATN"));
		final StringBuilder longTandem=new StringBuilder(128);
		for(int i=0; i<64; i++){longTandem.append("AT");}
		test.check("long terminal tandem", () -> tandem(input, longTandem.toString(), 128));
		test.check("long tandem sentinel control", () -> tandem(input, longTandem.toString()+"N", 128));
		test.check("zero-loop even palindrome", () -> palindrome("AT", 0, 1));
		test.check("one-loop odd palindrome control", () -> palindrome("ACT", 1, 2));
		test.check("two interior mismatches rejected at zero", () -> imperfectPalindrome("ACGATAAT", 0, false));
		test.check("two interior mismatches rejected at one", () -> imperfectPalindrome("ACGATAAT", 1, false));
		test.check("two interior mismatches accepted at two", () -> imperfectPalindrome("ACGATAAT", 2, true));
		test.check("reverse-complement mode control", () -> imperfectPalindrome("ACGTTGCA", 2, false));
		System.out.println("RepeatBoundaryTest cases="+test.cases+" failures="+test.failures);
		if(test.failures!=0){throw new AssertionError("Repeat boundary failures="+test.failures);}
	}

	/** The maximum observed depth must include either the first or the second window. */
	private static void depthPeak(boolean first){
		final Read contig=read("AAAC");
		final Repeat repeat=new Repeat(null, -1, 2, 3, 0, 4, 'R');
		require(repeat.increment(contig, 2, first ? 3 : 2)==null, "First window must open one region");
		require(repeat.increment(contig, 3, first ? 2 : 3)==null, "Adjacent window must extend the same region");
		require(repeat.start==0 && repeat.stop==3 && repeat.length()==4, "Literal AAAC region must span 0..3");
		require(repeat.maxDepth==3, "Maximum of observed depths 2 and 3 must be 3, actual="+repeat.maxDepth);
	}

	/** The exact ATAT input and its N-terminated control each cover four bases. */
	private static void tandem(Path input, String bases) throws Exception{
		tandem(input, bases, 4);
	}

	/** Checks a literal expected prefix, retaining the longer terminal regression too. */
	private static void tandem(Path input, String bases, int expected) throws Exception{
		final Object worker=worker(input, false);
		final int count=invokeCount(worker, "processReadShortTandem",
				new Class<?>[] {Read.class, int.class, int.class, int.class, int.class}, read(bases), 2, 2, 2, 4);
		require(count==expected, "Expected "+expected+" tandem bases in "+bases+", actual="+count);
		requireSingleRange(worker, 0, expected-1);
	}

	/** Literal one-base reverse-complement arms surround the requested loop. */
	private static void palindrome(String bases, int loop, int stop){
		final PalindromeFinder finder=new PalindromeFinder(1, loop, loop, 1, 0,
				0, Integer.MAX_VALUE, Integer.MAX_VALUE);
		final Palindrome p=finder.longestPalindrome(bases.getBytes(StandardCharsets.US_ASCII));
		require(p!=null, "Expected complementary A/T arms in "+bases+" with loop="+loop);
		require(p.a==0 && p.b==stop && p.matches==1 && p.mismatches==0 && p.plen()==1 && p.loop()==loop,
				"Unexpected literal palindrome: "+p);
	}

	/** Exact-region tails forbid accepting a shorter subregion instead of the literal stem.
	 * ACGATAAT pairs A/T, C/A, G/A, A/T: two complementary matches and two interior mismatches.
	 * ACGTTGCA has four equal pairs but no complementary pairs, distinguishing direct matching.
	 * The allowance of two exercises CrisprFinder's default maxPalMismatches setting. */
	private static void imperfectPalindrome(String bases, int allowance, boolean accepted){
		final PalindromeFinder finder=new PalindromeFinder(4, 0, 0, 2, allowance, 0, 0, 0);
		final Palindrome p=finder.longestPalindrome(bases.getBytes(StandardCharsets.US_ASCII));
		if(!accepted){
			require(p==null, "Expected rejection of "+bases+" with mismatch allowance="+allowance+", actual="+p);
		}else{
			require(p!=null && p.a==0 && p.b==7 && p.matches==2 && p.mismatches==2
					&& p.plen()==4 && p.loop()==0,
					"Expected exact reverse-complement stem with two interior mismatches, actual="+p);
		}
	}

	/** Creates a detector and an unstarted worker with tiny, valid entropy settings. */
	private static Object worker(Path input, boolean entropy) throws Exception{
		final RepeatFinder finder=new RepeatFinder(new String[] {"in="+input, "out=null", "depth=f",
				"entropy="+entropy, "tandem="+(!entropy), "k=3", "minrepeat=4", "ew=4", "ek=1", "mask=f"});
		final Class<?> type=Class.forName("repeat.RepeatFinder$ProcessThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(RepeatFinder.class,
				ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class, int.class);
		constructor.setAccessible(true);
		return constructor.newInstance(finder, null, null, 0);
	}

	/** Calls the existing per-read implementation and preserves its underlying failure. */
	private static int invokeCount(Object worker, String name, Class<?>[] types, Object... args) throws Exception{
		final Method method=worker.getClass().getDeclaredMethod(name, types);
		method.setAccessible(true);
		try{return ((Integer)method.invoke(worker, args)).intValue();}
		catch(InvocationTargetException failure){
			final Throwable cause=failure.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new AssertionError("Unexpected reflected failure", cause);
		}
	}

	/** Reads the worker-owned closed regions without making production members public. */
	private static List<?> closedRepeats(Object worker) throws Exception{
		final Field setField=worker.getClass().getDeclaredField("repeatSet");
		setField.setAccessible(true);
		final Object set=setField.get(worker);
		final Field listField=set.getClass().getDeclaredField("closedRepeats");
		listField.setAccessible(true);
		return (List<?>)listField.get(set);
	}

	/** Checks the literal inclusive coordinates as well as the reported base count. */
	private static void requireSingleRange(Object worker, int start, int stop) throws Exception{
		final List<?> repeats=closedRepeats(worker);
		require(repeats.size()==1, "Expected one closed repeat, actual="+repeats.size());
		final Repeat repeat=(Repeat)repeats.get(0);
		require(repeat.start==start && repeat.stop==stop, "Unexpected closed region: "+repeat.start+".."+repeat.stop);
	}

	/** Creates a separate caller-owned nucleotide array and contig for each case. */
	private static Read read(String bases){
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), null, "literal", 0, false);
	}

	/** Reports every independent case even when an earlier baseline expectation fails. */
	private void check(String name, Action action){
		cases++;
		try{action.run();System.out.println(name+" PASS");}
		catch(Exception | AssertionError failure){
			failures++;
			System.out.println(name+" FAIL: "+failure);
		}
	}

	/** Keeps literal expectations active under both JVM assertion modes. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** One bounded case, possibly using reflective access to the actual implementation. */
	private interface Action{void run() throws Exception;}
	/** Attempted and failed independent cases. */
	private int cases, failures;
}
