package test.map;

import java.lang.management.ManagementFactory;
import java.util.Locale;
import map.IntHashSetList;
import parse.Parse;

/** Fixed-work sparse and dense IntHashSetList clear benchmark.
 * @author Sucrose
 */
public final class IntHashSetListClearBenchmark{

	public static void main(String[] args){
		String mode="sparse";
		int size=65536, rounds=256;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=arg.substring(0, equals).toLowerCase(Locale.ROOT), b=arg.substring(equals+1);
			if(a.equals("mode")){mode=b.toLowerCase(Locale.ROOT);}
			else if(a.equals("size")){size=Parse.parseIntKMG(b);}
			else if(a.equals("rounds")){rounds=Parse.parseIntKMG(b);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(!mode.equals("sparse") && !mode.equals("dense")){
			throw new IllegalArgumentException("Unknown mode: "+mode);
		}
		if(size<1 || size>1<<24 || rounds<1 || rounds>1000000){
			throw new IllegalArgumentException("Require size in [1,16777216] and rounds in [1,1000000]");
		}

		final int[] keys=new int[size];
		for(int i=0; i<size; i++){keys[i]=i*0x9E3779B9;}
		run(keys, mode, 16);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long checksum=run(keys, mode, rounds);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		System.out.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d%n",
			mode, size, rounds, (long)size*rounds, elapsed*1e-9, cpu*1e-9, checksum);
	}

	private static long run(final int[] keys, final String mode, final int rounds){
		assert(keys.length>0 && rounds>0) : "Benchmark must execute a positive fixed workload";
		final boolean sparse=mode.equals("sparse");
		final IntHashSetList set=new IntHashSetList(sparse ? keys.length*8 : keys.length*2);
		long checksum=0;
		for(int r=0; r<rounds; r++){
			for(int key : keys){
				if(set.add(key)){checksum++;}
			}
			if(set.size()!=keys.length){throw new AssertionError("Missing benchmark keys");}
			if(!set.contains(keys[0]) || !set.contains(keys[keys.length-1])){
				throw new AssertionError("Benchmark set lost an endpoint key");
			}
			set.clear();
			if(!set.isEmpty()){
				throw new AssertionError("clear did not empty the set");
			}
		}
		if(set.toArray().length!=0){
			throw new AssertionError("clear did not empty the set array");
		}
		if(checksum!=(long)keys.length*rounds){
			throw new AssertionError("Benchmark checksum mismatch: "+checksum);
		}
		return checksum;
	}

	private IntHashSetListClearBenchmark(){
		throw new AssertionError();
	}
}
