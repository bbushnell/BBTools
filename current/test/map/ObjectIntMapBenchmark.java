package test.map;

import java.lang.management.ManagementFactory;
import java.util.Locale;
import map.ObjectIntMap;

/** Fixed-work increment benchmark; keys are allocated before warmup and timing.
 * Existing-key counts have an identical oracle on both revisions. Fresh/growth
 * modes report the old off-by-one values explicitly instead of hiding the repair.
 * @author Shinobu
 */
public final class ObjectIntMapBenchmark{

	public static void main(String[] args){
		String mode="existing";
		int size=131072, rounds=64;
		boolean fixed=true;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("size")){size=Integer.parseInt(pair[1]);}
			else if(pair[0].equals("rounds")){rounds=Integer.parseInt(pair[1]);}
			else if(pair[0].equals("fixed")){fixed=Boolean.parseBoolean(pair[1]);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(size<2 || rounds<1 || rounds>1000000 || Integer.bitCount(size)!=1){
			throw new IllegalArgumentException("Power-of-two size and positive bounded rounds required");
		}
		if(!mode.equals("existing") && !mode.equals("fresh") && !mode.equals("growth")){
			throw new IllegalArgumentException("Unknown mode: "+mode);
		}
		final String[] keys=new String[size];
		for(int i=0; i<size; i++){keys[i]="read_"+i; keys[i].hashCode();}
		run(keys, mode, 16, fixed);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long checksum=run(keys, mode, rounds, fixed);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		System.out.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d%n",
			mode, size, rounds, (long)size*rounds, elapsed*1e-9, cpu*1e-9, checksum);
	}

	/** Measures fixed work and consumes stored values; validates outside inner loops. */
	private static long run(String[] keys, String mode, int rounds, boolean fixed){
		assert(keys.length>0 && rounds>0) : "Benchmark must perform measurable work";
		final boolean existing=mode.equals("existing"), growth=mode.equals("growth");
		ObjectIntMap<String> map=null;
		if(existing){
			map=new ObjectIntMap<String>(keys.length*2, String.class);
			for(String key : keys){map.put(key, 0);}
		}
		long checksum=0;
		for(int r=0; r<rounds; r++){
			if(!existing){map=new ObjectIntMap<String>(growth ? 256 : keys.length*2, String.class);}
			for(int i=0; i<keys.length; i++){checksum+=map.increment(keys[(i*8191)&(keys.length-1)]);}
			if(map.size()!=keys.length){throw new AssertionError("Missing benchmark keys");}
		}
		final int expected=existing ? rounds : fixed ? 1 : 0;
		for(String key : keys){
			if(map.get(key)!=expected){throw new AssertionError("Stored count disagrees with revision oracle: "+key);}
		}
		if(existing && checksum!=(long)keys.length*rounds*(rounds+1)/2){
			throw new AssertionError("Increment return values disagree with arithmetic sum");
		}
		return checksum;
	}
}
