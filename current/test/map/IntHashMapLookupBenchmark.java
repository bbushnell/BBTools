package test.map;

import java.lang.management.ManagementFactory;
import java.util.Locale;
import map.IntHashMap;
import parse.Parse;
import shared.Primes;

/** Fixed-work primitive-key benchmark for IntHashMap get and contains lookups.
 * @author Jean
 */
public final class IntHashMapLookupBenchmark{

	/** Benchmarks an existing/missing get/contains mode or all four modes. */
	public static void main(String[] args){
		String mode="all";
		String pattern="sequential";
		int size=131072, rounds=256;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=arg.substring(0, equals).toLowerCase(Locale.ROOT), b=arg.substring(equals+1);
			if(a.equals("mode")){mode=canonicalMode(b);}
			else if(a.equals("pattern")){pattern=canonicalPattern(b);}
			else if(a.equals("size")){size=Parse.parseIntKMG(b);}
			else if(a.equals("rounds")){rounds=Parse.parseIntKMG(b);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(size<2 || size>1<<24 || Integer.bitCount(size)!=1 || rounds<1 || rounds>1000000){
			throw new IllegalArgumentException("Require power-of-two size in [2,16777216] and rounds in [1,1000000]");
		}

		final int[] existingKeys=new int[size], missingKeys=new int[size];
		final boolean mixed=pattern.equals("mixed");
		for(int i=0; i<size; i++){
			final int key=(mixed ? (i+1)*0x9E3779B9 : i+1);
			existingKeys[i]=(key<<1)-1;
			missingKeys[i]=key<<1;
		}
		if(mixed){checkMixedWorkload(existingKeys);}
		final String prefix=(mixed ? "mixed-" : "");
		if(mode.equals("all") || mode.equals("existingget")){
			bench(prefix+"existingget", "existingget", existingKeys, missingKeys, rounds);
		}
		if(mode.equals("all") || mode.equals("missingget")){
			bench(prefix+"missingget", "missingget", existingKeys, missingKeys, rounds);
		}
		if(mode.equals("all") || mode.equals("existingcontains")){
			bench(prefix+"existingcontains", "existingcontains", existingKeys, missingKeys, rounds);
		}
		if(mode.equals("all") || mode.equals("missingcontains")){
			bench(prefix+"missingcontains", "missingcontains", existingKeys, missingKeys, rounds);
		}
	}

	/** Normalizes command-line mode names. */
	private static String canonicalMode(String mode){
		mode=mode.toLowerCase(Locale.ROOT).replace("-", "").replace("_", "");
		if(mode.equals("all") || mode.equals("existingget") || mode.equals("missingget") ||
				mode.equals("existingcontains") || mode.equals("missingcontains")){
			return mode;
		}
		if(mode.equals("getexisting")){return "existingget";}
		if(mode.equals("getmissing")){return "missingget";}
		if(mode.equals("containsexisting")){return "existingcontains";}
		if(mode.equals("containsmissing")){return "missingcontains";}
		throw new IllegalArgumentException("Unknown mode: "+mode);
	}

	/** Normalizes command-line key-pattern names. */
	private static String canonicalPattern(String pattern){
		pattern=pattern.toLowerCase(Locale.ROOT).replace("-", "").replace("_", "");
		if(pattern.equals("sequential") || pattern.equals("seq")){return "sequential";}
		if(pattern.equals("mixed") || pattern.equals("mix")){return "mixed";}
		throw new IllegalArgumentException("Unknown pattern: "+pattern);
	}

	/** Runs one warmup and one measured fixed-work lookup loop. */
	private static void bench(final String label, final String mode,
			final int[] existingKeys, final int[] missingKeys,
			final int rounds){
		final IntHashMap map=makeMap(existingKeys);
		run(map, existingKeys, missingKeys, mode, 16);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long checksum=run(map, existingKeys, missingKeys, mode, rounds);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		checkChecksum(mode, existingKeys.length, rounds, checksum);
		System.out.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d%n",
			label, existingKeys.length, rounds, (long)existingKeys.length*rounds,
			elapsed*1e-9, cpu*1e-9, checksum);
	}

	/** Builds the shared read-only map outside the timed loop. */
	private static IntHashMap makeMap(final int[] keys){
		final IntHashMap map=new IntHashMap(keys.length*2);
		for(int i=0; i<keys.length; i++){map.put(keys[i], i+1);}
		if(map.size()!=keys.length){throw new AssertionError("Missing benchmark keys after construction");}
		return map;
	}

	/** Verifies that the mixed pattern reaches collision-continuation cells. */
	private static void checkMixedWorkload(final int[] existingKeys){
		final IntHashMap map=makeMap(existingKeys);
		final int modulus=(int)Primes.primeAtLeast(2L*existingKeys.length);
		int displaced=0;
		for(int key : existingKeys){
			final int initial=(key&Integer.MAX_VALUE)%modulus;
			if(map.findCell(key)!=initial){displaced++;}
		}
		if(displaced<1){throw new AssertionError("Mixed key pattern did not displace any keys");}
		System.err.println("mixed-displaced\t"+displaced);
	}

	/** Dispatches once per measured run, outside all lookup loops. */
	private static long run(final IntHashMap map, final int[] existingKeys, final int[] missingKeys,
			final String mode, final int rounds){
		assert(existingKeys.length==missingKeys.length && existingKeys.length>0 && rounds>0) :
			"Benchmark must execute a positive fixed workload";
		if(mode.equals("existingget")){return runGet(map, existingKeys, rounds);}
		if(mode.equals("missingget")){return runGet(map, missingKeys, rounds);}
		if(mode.equals("existingcontains")){return runContains(map, existingKeys, rounds);}
		if(mode.equals("missingcontains")){return runContains(map, missingKeys, rounds);}
		throw new AssertionError("Unknown mode: "+mode);
	}

	/** Consumes IntHashMap.get results for a primitive precomputed key array. */
	private static long runGet(final IntHashMap map, final int[] keys, final int rounds){
		final int mask=keys.length-1;
		long checksum=0;
		for(int r=0; r<rounds; r++){
			for(int i=0; i<keys.length; i++){checksum+=map.get(keys[(i*8191)&mask]);}
		}
		return checksum;
	}

	/** Consumes IntHashMap.contains results for a primitive precomputed key array. */
	private static long runContains(final IntHashMap map, final int[] keys, final int rounds){
		final int mask=keys.length-1;
		long checksum=0;
		for(int r=0; r<rounds; r++){
			for(int i=0; i<keys.length; i++){
				checksum+=map.contains(keys[(i*8191)&mask]) ? 1 : 0;
			}
		}
		return checksum;
	}

	/** Checks the deterministic sum expected from the chosen lookup mode. */
	private static void checkChecksum(final String mode, final int size, final int rounds,
			final long checksum){
		final long expected;
		if(mode.equals("existingget")){expected=(long)size*(size+1)/2*rounds;}
		else if(mode.equals("missingget")){expected=-(long)size*rounds;}
		else if(mode.equals("existingcontains")){expected=(long)size*rounds;}
		else if(mode.equals("missingcontains")){expected=0;}
		else{throw new AssertionError("Unknown mode: "+mode);}
		if(checksum!=expected){throw new AssertionError(mode+" checksum "+checksum+" != "+expected);}
	}

	private IntHashMapLookupBenchmark(){
		throw new AssertionError();
	}
}
