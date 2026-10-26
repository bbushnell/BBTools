package test.map;

import java.lang.management.ManagementFactory;
import java.util.Locale;
import map.IntHashMap2;
import parse.Parse;

/** Fixed-work positive-count merge benchmark for the variant used by bin.Cluster.
 * @author Shinobu
 */
public final class MapSetToMaxBenchmark{

	public static void main(String[] args){
		String mode="existing";
		int size=131072, rounds=256;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=arg.substring(0, equals).toLowerCase(Locale.ROOT), b=arg.substring(equals+1);
			if(a.equals("mode")){mode=b;}
			else if(a.equals("size")){size=Parse.parseIntKMG(b);}
			else if(a.equals("rounds")){rounds=Parse.parseIntKMG(b);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if((!mode.equals("existing") && !mode.equals("absent")) || size<2 || size>1<<24 ||
			Integer.bitCount(size)!=1 || rounds<1 || rounds>1000000){
			throw new IllegalArgumentException("Require mode=existing|absent, power-of-two size 2..16777216, rounds 1..1000000");
		}
		final IntHashMap2 source=new IntHashMap2(size*2), dest=new IntHashMap2(size*2);
		final int[] keys=new int[size];
		for(int i=0; i<size; i++){
			keys[i]=i*0x9E3779B9;
			source.put(keys[i], i+1);
			dest.put(keys[i], 0);
		}
		final boolean absent=mode.equals("absent");
		run(source, dest, keys, absent, 16);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long checksum=run(source, dest, keys, absent, rounds);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		final long remainder=rounds%size;
		final long expected=(long)(rounds/size)*size*(size+1)/2+remainder*(remainder+1)/2;
		check(checksum==expected && dest.size()==size, "Merge benchmark checksum/size");
		for(int i=0; i<size; i++){
			check(source.get(keys[i])==i+1 && dest.get(keys[i])==i+1, "Stored count at "+i);
		}
		System.out.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d%n",
			mode, size, rounds, (long)size*rounds, elapsed*1e-9, cpu*1e-9, checksum);
	}

	private static long run(final IntHashMap2 source, final IntHashMap2 dest,
		final int[] keys, final boolean absent, final int rounds){
		assert(source!=dest && keys.length>0 && rounds>0) : "Benchmark needs distinct maps and positive work";
		long checksum=0;
		for(int r=0; r<rounds; r++){
			if(absent){dest.clear();}
			dest.setToMax(source);
			checksum+=dest.get(keys[r&(keys.length-1)]);
		}
		return checksum;
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
