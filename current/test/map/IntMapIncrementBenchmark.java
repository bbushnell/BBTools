package test.map;

import java.lang.management.ManagementFactory;
import java.util.Locale;
import map.IntHashMap2;
import map.IntHashMap3;
import map.IntLongHashMap2;
import map.IntLongHashMap3;
import parse.Parse;

/** Fixed-work int-key increment benchmark with primitive precomputed keys.
 * @author Sucrose
 */
public final class IntMapIncrementBenchmark{

	/** Benchmarks one map type or all four confirmed stale-return siblings. */
	public static void main(String[] args){
		String type="all", mode="existing";
		int size=131072, rounds=256;
		boolean fixed=true;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=arg.substring(0, equals).toLowerCase(Locale.ROOT), b=arg.substring(equals+1);
			if(a.equals("type")){type=b.toLowerCase(Locale.ROOT);}
			else if(a.equals("mode")){mode=b.toLowerCase(Locale.ROOT);}
			else if(a.equals("size")){size=Parse.parseIntKMG(b);}
			else if(a.equals("rounds")){rounds=Parse.parseIntKMG(b);}
			else if(a.equals("fixed")){fixed=Parse.parseBoolean(b);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(size<2 || size>1<<24 || Integer.bitCount(size)!=1 || rounds<1 || rounds>1000000){
			throw new IllegalArgumentException("Require power-of-two size in [2,16777216] and rounds in [1,1000000]");
		}
		if(!mode.equals("existing") && !mode.equals("fresh") && !mode.equals("growth")){
			throw new IllegalArgumentException("Unknown mode: "+mode);
		}
		if(!type.equals("all") && !type.equals("ihm2") && !type.equals("ihm3") &&
			!type.equals("ilhm2") && !type.equals("ilhm3")){
			throw new IllegalArgumentException("Unknown type: "+type);
		}
		final int[] keys=new int[size];
		for(int i=0; i<size; i++){keys[i]=i*0x9E3779B9;}
		if(type.equals("all") || type.equals("ihm2")){bench("ihm2", keys, mode, rounds, fixed);}
		if(type.equals("all") || type.equals("ihm3")){bench("ihm3", keys, mode, rounds, fixed);}
		if(type.equals("all") || type.equals("ilhm2")){bench("ilhm2", keys, mode, rounds, fixed);}
		if(type.equals("all") || type.equals("ilhm3")){bench("ilhm3", keys, mode, rounds, fixed);}
	}

	/** Runs a warmup, a timed fixed workload and one tabular report row. */
	private static void bench(final String type, final int[] keys, final String mode, final int rounds,
		final boolean fixed){
		run(type, keys, mode, 16, fixed);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long checksum=run(type, keys, mode, rounds, fixed);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		System.out.printf(Locale.ROOT, "%s	%s	%d	%d	%d	%.6f	%.6f	%d%n",
			type, mode, keys.length, rounds, (long)keys.length*rounds, elapsed*1e-9, cpu*1e-9, checksum);
	}

	/** Dispatches to the monomorphic primitive-key loop for the requested map type. */
	private static long run(final String type, final int[] keys, final String mode, final int rounds,
		final boolean fixed){
		if(type.equals("ihm2")){return runIntHashMap2(keys, mode, rounds, fixed);}
		if(type.equals("ihm3")){return runIntHashMap3(keys, mode, rounds, fixed);}
		if(type.equals("ilhm2")){return runIntLongHashMap2(keys, mode, rounds, fixed);}
		if(type.equals("ilhm3")){return runIntLongHashMap3(keys, mode, rounds, fixed);}
		throw new AssertionError("Unknown type: "+type);
	}

	/** Consumes IntHashMap2 returned values and checks final stored counts. */
	private static long runIntHashMap2(final int[] keys, final String mode, final int rounds, final boolean fixed){
		assert(keys.length>0 && rounds>0) : "Benchmark must execute a positive fixed workload";
		final boolean existing=mode.equals("existing"), growth=mode.equals("growth");
		IntHashMap2 map=null;
		if(existing){
			map=new IntHashMap2(keys.length*2);
			for(int key : keys){map.put(key, 0);}
		}
		long checksum=0;
		for(int r=0; r<rounds; r++){
			if(!existing){map=new IntHashMap2(growth ? 256 : keys.length*2);}
			for(int i=0; i<keys.length; i++){checksum+=map.increment(keys[(i*8191)&(keys.length-1)]);}
			if(map.size()!=keys.length){throw new AssertionError("Missing IntHashMap2 benchmark keys");}
		}
		checkIntCounts(map, keys, existing ? rounds : 1);
		checkChecksum(checksum, keys.length, rounds, existing, fixed || !growth);
		return checksum;
	}

	/** Consumes IntHashMap3 returned values and checks final stored counts. */
	private static long runIntHashMap3(final int[] keys, final String mode, final int rounds, final boolean fixed){
		assert(keys.length>0 && rounds>0) : "Benchmark must execute a positive fixed workload";
		final boolean existing=mode.equals("existing"), growth=mode.equals("growth");
		IntHashMap3 map=null;
		if(existing){
			map=new IntHashMap3(keys.length*2);
			for(int key : keys){map.put(key, 0);}
		}
		long checksum=0;
		for(int r=0; r<rounds; r++){
			if(!existing){map=new IntHashMap3(growth ? 256 : keys.length*2);}
			for(int i=0; i<keys.length; i++){checksum+=map.increment(keys[(i*8191)&(keys.length-1)]);}
			if(map.size()!=keys.length){throw new AssertionError("Missing IntHashMap3 benchmark keys");}
		}
		checkIntCounts(map, keys, existing ? rounds : 1);
		checkChecksum(checksum, keys.length, rounds, existing, fixed || !growth);
		return checksum;
	}

	/** Consumes IntLongHashMap2 returned values and checks final stored counts. */
	private static long runIntLongHashMap2(final int[] keys, final String mode, final int rounds, final boolean fixed){
		assert(keys.length>0 && rounds>0) : "Benchmark must execute a positive fixed workload";
		final boolean existing=mode.equals("existing"), growth=mode.equals("growth");
		IntLongHashMap2 map=null;
		if(existing){
			map=new IntLongHashMap2(keys.length*2);
			for(int key : keys){map.put(key, 0);}
		}
		long checksum=0;
		for(int r=0; r<rounds; r++){
			if(!existing){map=new IntLongHashMap2(growth ? 256 : keys.length*2);}
			for(int i=0; i<keys.length; i++){checksum+=map.increment(keys[(i*8191)&(keys.length-1)]);}
			if(map.size()!=keys.length){throw new AssertionError("Missing IntLongHashMap2 benchmark keys");}
		}
		checkLongCounts(map, keys, existing ? rounds : 1);
		checkChecksum(checksum, keys.length, rounds, existing, fixed || !growth);
		return checksum;
	}

	/** Consumes IntLongHashMap3 returned values and checks final stored counts. */
	private static long runIntLongHashMap3(final int[] keys, final String mode, final int rounds, final boolean fixed){
		assert(keys.length>0 && rounds>0) : "Benchmark must execute a positive fixed workload";
		final boolean existing=mode.equals("existing"), growth=mode.equals("growth");
		IntLongHashMap3 map=null;
		if(existing){
			map=new IntLongHashMap3(keys.length*2);
			for(int key : keys){map.put(key, 0);}
		}
		long checksum=0;
		for(int r=0; r<rounds; r++){
			if(!existing){map=new IntLongHashMap3(growth ? 256 : keys.length*2);}
			for(int i=0; i<keys.length; i++){checksum+=map.increment(keys[(i*8191)&(keys.length-1)]);}
			if(map.size()!=keys.length){throw new AssertionError("Missing IntLongHashMap3 benchmark keys");}
		}
		checkLongCounts(map, keys, existing ? rounds : 1);
		checkChecksum(checksum, keys.length, rounds, existing, fixed || !growth);
		return checksum;
	}

	/** Checks final IntHashMap2 counts outside the measured increment loop. */
	private static void checkIntCounts(final IntHashMap2 map, final int[] keys, final int expected){
		for(int key : keys){
			if(map.get(key)!=expected){throw new AssertionError("Stored IntHashMap2 count differs from fixed-work oracle");}
		}
	}

	/** Checks final IntHashMap3 counts outside the measured increment loop. */
	private static void checkIntCounts(final IntHashMap3 map, final int[] keys, final int expected){
		for(int key : keys){
			if(map.get(key)!=expected){throw new AssertionError("Stored IntHashMap3 count differs from fixed-work oracle");}
		}
	}

	/** Checks final IntLongHashMap2 counts outside the measured increment loop. */
	private static void checkLongCounts(final IntLongHashMap2 map, final int[] keys, final int expected){
		for(int key : keys){
			if(map.get(key)!=expected){throw new AssertionError("Stored IntLongHashMap2 count differs from fixed-work oracle");}
		}
	}

	/** Checks final IntLongHashMap3 counts outside the measured increment loop. */
	private static void checkLongCounts(final IntLongHashMap3 map, final int[] keys, final int expected){
		for(int key : keys){
			if(map.get(key)!=expected){throw new AssertionError("Stored IntLongHashMap3 count differs from fixed-work oracle");}
		}
	}

	/** Checks the deterministic sum of all returned increment values when the mode is not expected to fail. */
	private static void checkChecksum(final long checksum, final int size, final int rounds,
		final boolean existing, final boolean expectCorrect){
		final long expected=existing ? (long)size*rounds*(rounds+1)/2 : (long)size*rounds;
		if(expectCorrect && checksum!=expected){
			throw new AssertionError("Increment return checksum "+checksum+" != "+expected);
		}
	}
}
