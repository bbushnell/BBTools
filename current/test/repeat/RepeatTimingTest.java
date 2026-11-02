package test.repeat;

import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

import repeat.Palindrome;
import repeat.PalindromeFinder;
import repeat.Repeat;
import repeat.RepeatFinder;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/** Deterministic, single-caller timing of repeat detectors on preallocated inputs.
 * Runs literal output checks before warmup and checks both accumulated checksums.
 * Entropy/tandem use unstarted workers with cached method handles; the timed scope
 * includes clearing their two output lists after each call. Production allocation
 * and statistics calculation remain included. Construction, input creation and
 * file validation are excluded. This is not full CLI or worker-pipeline throughput.
 * Use identical JVM/settings and enough rounds for at least ten seconds of compute;
 * the small defaults are a correctness/calibration run, not performance evidence.
 * @author Sucrose
 */
public final class RepeatTimingTest{

	/** Accepts mode=entropy|tandem|palindrome|imperfect, rounds=N, warmup=N, work=fresh-directory.
	 * The two palindrome modes share literal inputs/oracles; imperfect permits one mismatch. */
	public static void main(String[] args) throws Throwable{
		String mode=null;
		Path work=null;
		long rounds=10000, warmup=1000;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0, eq), value=arg.substring(eq+1);
			if(key.equals("mode")){mode=value;}
			else if(key.equals("rounds")){rounds=Long.parseLong(value);}
			else if(key.equals("warmup")){warmup=Long.parseLong(value);}
			else if(key.equals("work")){work=Paths.get(value);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		require("entropy".equals(mode) || "tandem".equals(mode) || "palindrome".equals(mode) || "imperfect".equals(mode),
				"Choose mode=entropy, tandem, palindrome or imperfect");
		require(work!=null, "work=fresh-directory is required");
		require(rounds>0 && warmup>=0, "Rounds must be positive and warmup nonnegative");
		final long expectedPerRound=(mode.equals("entropy") ? 64 : mode.equals("tandem") ? 256 : 770);
		require(rounds<=Long.MAX_VALUE/expectedPerRound && warmup<=Long.MAX_VALUE/expectedPerRound,
				"Round count would overflow the exact accumulated checksum");
		Files.createDirectory(work);
		final Path input=work.resolve("repeat-timing-input.fa");
		Files.write(input, ">fixture\nACGTACGT\n".getBytes(StandardCharsets.US_ASCII));
		final RepeatTimingTest test=new RepeatTimingTest(mode, input);
		test.verify();
		final long warmed=test.run(warmup);
		require(warmed==warmup*expectedPerRound, "Warmup checksum differs from literal expectation: "+warmed);
		final long start=System.nanoTime();
		final long checksum=test.run(rounds);
		final long elapsed=System.nanoTime()-start;
		require(checksum==rounds*expectedPerRound, "Measured checksum differs from literal expectation: "+checksum);
		System.out.println("mode="+mode+" rounds="+rounds+" warmup="+warmup+" calls="+(2*rounds)
				+" nanoseconds="+elapsed+" checksum="+checksum+" expected="+(rounds*expectedPerRound));
	}

	/** Builds one detector and two alternating contigs; no processing thread is started. */
	private RepeatTimingTest(String mode, Path input) throws Exception{
		palindromeMode=mode.equals("palindrome") || mode.equals("imperfect");
		if(palindromeMode){
			palindromeFinder=new PalindromeFinder(4, 1, 1, 4, mode.equals("imperfect") ? 1 : 0,
					0, Integer.MAX_VALUE, Integer.MAX_VALUE);
			reads=new Read[] {read(palindromeBases(0), 0), read(palindromeBases(1), 1)};
			detector=null;
			closed=recent=null;
			expectedBases=257;
		}else{
			palindromeFinder=null;
			final boolean entropy=mode.equals("entropy");
			expectedBases=(entropy ? 32 : 128);
			reads=(entropy ? new Read[] {read(entropyBases(0), 0), read(entropyBases(1), 1)} :
					new Read[] {read(tandemBases('A', 'T'), 0), read(tandemBases('C', 'G'), 1)});
			final RepeatFinder finder=new RepeatFinder(new String[] {"in="+input, "out=null", "depth=f",
					"entropy="+(entropy ? "0.1" : "f"), "tandem="+(!entropy), "k=3",
					"minrepeat="+(entropy ? 8 : 4), "ew=4", "ek=1", "mask=f"});
			final Class<?> type=Class.forName("repeat.RepeatFinder$ProcessThread");
			final Constructor<?> constructor=type.getDeclaredConstructor(RepeatFinder.class,
					ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class, int.class);
			constructor.setAccessible(true);
			final Object worker=constructor.newInstance(finder, null, null, 0);
			final Method method=(entropy ? type.getDeclaredMethod("processReadEntropy", Read.class) :
					type.getDeclaredMethod("processReadShortTandem", Read.class, int.class, int.class, int.class, int.class));
			method.setAccessible(true);
			final MethodHandle bound=MethodHandles.lookup().unreflect(method).bindTo(worker);
			detector=(entropy ? bound : MethodHandles.insertArguments(bound, 1, 2, 2, 2, 4));
			final Field setField=type.getDeclaredField("repeatSet");
			setField.setAccessible(true);
			final Object set=setField.get(worker);
			closed=list(set, "closedRepeats");
			recent=list(set, "recent");
		}
	}

	/** Checks literal full coordinates/counts on each input, outside measured work. */
	private void verify() throws Throwable{
		for(Read read : reads){
			if(palindromeMode){
				final Palindrome p=palindromeFinder.longestPalindrome(read.bases);
				require(p!=null && p.a==0 && p.b==256 && p.matches==128 && p.mismatches==0
						&& p.plen()==128 && p.loop()==1, "Expected two 128-base arms around one N");
			}else{
				final int count=(int)detector.invokeExact(read);
				require(count==expectedBases && closed.size()==1 && recent.size()==1,
						"Expected one region and "+expectedBases+" bases, actual="+count+" regions="+closed.size());
				final Repeat r=(Repeat)closed.get(0);
				require(r.contig==read && r.start==0 && r.stop==expectedBases-1,
						"Unexpected coordinates/contig for literal detector input");
				closed.clear();
				recent.clear();
			}
		}
	}

	/** Times actual detector calls plus the disclosed result-list reset; dispatch is outside loops. */
	private long run(long rounds) throws Throwable{
		long sum=0;
		if(palindromeMode){
			for(long i=0; i<rounds; i++){
				for(Read read : reads){
					final Palindrome p=palindromeFinder.longestPalindrome(read.bases);
					if(p==null){throw new AssertionError("Literal palindrome disappeared during timing");}
					sum+=p.a+p.b+p.matches+p.mismatches+p.loop();
				}
			}
		}else{
			for(long i=0; i<rounds; i++){
				for(Read read : reads){
					sum+=(int)detector.invokeExact(read);
					closed.clear();
					recent.clear();
				}
			}
		}
		return sum;
	}

	/** One 32-base uniform region, balanced middle, and a subthreshold terminal region.
	 * This makes the count equal on both sides of the separate terminal-count repair. */
	private static byte[] entropyBases(int phase){
		final byte[] bases=new byte[100];
		final byte[] alphabet="ACGT".getBytes(StandardCharsets.US_ASCII);
		Arrays.fill(bases, 0, 32, alphabet[phase]);
		for(int i=32; i<96; i++){bases[i]=alphabet[(i-32+phase+1)&3];}
		Arrays.fill(bases, 96, 100, alphabet[phase]);
		return bases;
	}

	/** Sixty-four exact two-base copies followed by N; the unmatched sentinel makes
	 * the terminal seed visible under both old and repaired scan bounds. */
	private static byte[] tandemBases(char a, char b){
		final byte[] bases=new byte[129];
		for(int i=0; i<128; i++){bases[i]=(byte)((i&1)==0 ? a : b);}
		bases[128]='N';
		return bases;
	}

	/** Two exact 128-base reverse-complement arms surrounding one unmatched center. */
	private static byte[] palindromeBases(int phase){
		final byte[] bases=new byte[257];
		final byte[] alphabet="ACGT".getBytes(StandardCharsets.US_ASCII);
		for(int i=0; i<128; i++){
			final int index=(i+phase)&3;
			bases[i]=alphabet[index];
			bases[256-i]=alphabet[3-index];
		}
		bases[128]='N';
		return bases;
	}

	/** Two alternating IDs prevent last-repeat subsumption across calls after list clearing. */
	private static Read read(byte[] bases, long id){return new Read(bases, null, "timing"+id, id, false);}

	/** Borrows a worker list once, avoiding reflective field access inside timed loops. */
	private static List<?> list(Object set, String name) throws Exception{
		final Field field=set.getClass().getDeclaredField(name);
		field.setAccessible(true);
		return (List<?>)field.get(set);
	}

	/** Keeps independent checks active with assertions disabled too. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Immutable mode and exact expected region size for the finite input checks. */
	private final boolean palindromeMode;
	private final int expectedBases;
	/** Caller-owned inputs, reused without mutation. */
	private final Read[] reads;
	/** Exactly one of the direct finder or bound worker method is used in a run. */
	private final PalindromeFinder palindromeFinder;
	private final MethodHandle detector;
	/** Worker-owned result lists cleared after each invocation to bound retained memory. */
	private final List<?> closed, recent;
}
