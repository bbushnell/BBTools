package test.map;

import java.lang.management.ManagementFactory;
import java.util.Locale;
import map.LongIntMap;
import parse.Parse;

/** Fixed-work primitive-key benchmark with no per-increment fixture allocation.
 * @author Shinobu
 */
public final class LongIntMapBenchmark{

	public static void main(String[] args){
		String mode="existing";
		int size=131072, rounds=256;
		boolean fixed=true;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=arg.substring(0, equals).toLowerCase(Locale.ROOT), b=arg.substring(equals+1);
			if(a.equals("mode")){mode=b;}
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
		final long[] keys=new long[size];
		for(int i=0; i<size; i++){keys[i]=i*0x9E3779B97F4A7C15L;}
		run(keys, mode, 16, fixed);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long checksum=run(keys, mode, rounds, fixed);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		System.out.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d%n",
			mode, size, rounds, (long)size*rounds, elapsed*1e-9, cpu*1e-9, checksum);
	}

	/** Consumes returned values and checks final counts outside the increment loop. */
	private static long run(final long[] keys, final String mode, final int rounds, final boolean fixed){
		assert(keys.length>0 && rounds>0) : "Benchmark must execute a positive fixed workload";
		final boolean existing=mode.equals("existing"), growth=mode.equals("growth");
		LongIntMap map=null;
		if(existing){
			map=new LongIntMap(keys.length*2);
			for(long key : keys){map.put(key, 0);}
		}
		long checksum=0;
		for(int r=0; r<rounds; r++){
			if(!existing){map=new LongIntMap(growth ? 256 : keys.length*2);}
			for(int i=0; i<keys.length; i++){checksum+=map.increment(keys[(i*8191)&(keys.length-1)]);}
			if(map.size()!=keys.length){throw new AssertionError("Missing benchmark keys");}
		}
		final int expected=existing ? rounds : 1;
		for(long key : keys){
			if(map.get(key)!=expected){throw new AssertionError("Stored count differs from fixed-work oracle");}
		}
		final long expectedSum=existing ? (long)keys.length*rounds*(rounds+1)/2 : (long)keys.length*rounds;
		if((fixed || !growth) && checksum!=expectedSum){
			throw new AssertionError("Increment return checksum "+checksum+" != "+expectedSum);
		}
		return checksum;
	}
}
