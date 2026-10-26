package test.map;

import java.lang.management.ManagementFactory;
import java.util.Locale;
import map.LongIntMapX;
import map.LongLongHashMap2;
import parse.Parse;

/** Fixed-work primitive-key benchmarks for LongIntMapX and LongLongHashMap2 increments.
 * @author Jean
 */
public final class LongMapHotPathBenchmark{

	public static void main(String[] args){
		String mapName="intx", mode="existing";
		int size=131072, rounds=256;
		boolean fixed=true;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=arg.substring(0, equals).toLowerCase(Locale.ROOT), b=arg.substring(equals+1);
			if(a.equals("map")){mapName=b.toLowerCase(Locale.ROOT);}
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
		if(!mapName.equals("intx") && !mapName.equals("long2")){
			throw new IllegalArgumentException("Unknown map: "+mapName);
		}

		final long[] keys=new long[size];
		for(int i=0; i<size; i++){keys[i]=i*0x9E3779B97F4A7C15L;}
		run(mapName, keys, mode, 16, fixed);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long checksum=run(mapName, keys, mode, rounds, fixed);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		System.out.printf(Locale.ROOT, "%s\t%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d%n",
			mapName, mode, size, rounds, (long)size*rounds, elapsed*1e-9, cpu*1e-9, checksum);
	}

	/** Dispatches once per measured run, outside all increment loops. */
	private static long run(final String mapName, final long[] keys, final String mode,
			final int rounds, final boolean fixed){
		assert(keys.length>0 && rounds>0) : "Benchmark must execute a positive fixed workload";
		if(mapName.equals("intx")){return runIntX(keys, mode, rounds, fixed);}
		if(mapName.equals("long2")){return runLong2(keys, mode, rounds, fixed);}
		throw new IllegalArgumentException("Unknown map: "+mapName);
	}

	/** Consumes returned int counts and checks final stored counts outside the hot loop. */
	private static long runIntX(final long[] keys, final String mode, final int rounds, final boolean fixed){
		final boolean existing=mode.equals("existing"), growth=mode.equals("growth");
		LongIntMapX map=null;
		if(existing){
			map=new LongIntMapX(keys.length*2);
			for(long key : keys){map.put(key, 1);}
		}
		long checksum=0, expectedChecksum=0;
		for(int r=0; r<rounds; r++){
			if(!existing){map=new LongIntMapX(growth ? 256 : keys.length*2);}
			for(int i=0; i<keys.length; i++){
				checksum+=map.increment(keys[(i*8191)&(keys.length-1)]);
			}
			if(map.size()!=keys.length){throw new AssertionError("LongIntMapX missing benchmark keys");}
			expectedChecksum+=(long)keys.length*(existing ? r+2 : 1);
		}
		final int expected=existing ? rounds+1 : 1;
		for(long key : keys){
			if(map.get(key)!=expected){throw new AssertionError("LongIntMapX stored count differs from oracle");}
		}
		if((fixed || !growth) && checksum!=expectedChecksum){
			throw new AssertionError("LongIntMapX checksum "+checksum+" != "+expectedChecksum);
		}
		return checksum;
	}

	/** Consumes returned long counts and checks final stored counts outside the hot loop. */
	private static long runLong2(final long[] keys, final String mode, final int rounds, final boolean fixed){
		final boolean existing=mode.equals("existing"), growth=mode.equals("growth");
		LongLongHashMap2 map=null;
		if(existing){
			map=new LongLongHashMap2(keys.length*2);
			for(long key : keys){map.put(key, 1);}
		}
		long checksum=0, expectedChecksum=0;
		for(int r=0; r<rounds; r++){
			if(!existing){map=new LongLongHashMap2(growth ? 256 : keys.length*2);}
			for(int i=0; i<keys.length; i++){
				checksum+=map.increment(keys[(i*8191)&(keys.length-1)]);
			}
			if(map.size()!=keys.length){throw new AssertionError("LongLongHashMap2 missing benchmark keys");}
			expectedChecksum+=(long)keys.length*(existing ? r+2 : 1);
		}
		final long expected=existing ? rounds+1 : 1;
		for(long key : keys){
			if(map.get(key)!=expected){throw new AssertionError("LongLongHashMap2 stored count differs from oracle");}
		}
		if((fixed || !growth) && checksum!=expectedChecksum){
			throw new AssertionError("LongLongHashMap2 checksum "+checksum+" != "+expectedChecksum);
		}
		return checksum;
	}

	private LongMapHotPathBenchmark(){
		throw new AssertionError();
	}
}
