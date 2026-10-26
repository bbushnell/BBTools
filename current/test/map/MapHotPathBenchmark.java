package test.map;

import java.lang.management.ManagementFactory;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;
import map.AbstractBitSet;
import map.ObjectDoubleMap;
import map.SeqMap;
import map.SeqPosM;
import parse.Parse;

/** Fixed-work measurements for object counts, packed counters and sequence lookup.
 * @author Shinobu
 */
public final class MapHotPathBenchmark{

	public static void main(String[] args){
		String mode="double-existing";
		int rounds=256;
		boolean fixed=true;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=arg.substring(0, equals).toLowerCase(Locale.ROOT), b=arg.substring(equals+1);
			if(a.equals("mode")){mode=b;}
			else if(a.equals("rounds")){rounds=Parse.parseIntKMG(b);}
			else if(a.equals("fixed")){fixed=Parse.parseBoolean(b);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(rounds<1 || rounds>10000000){throw new IllegalArgumentException("Rounds must be in [1,10000000]");}
		final String[] keys=new String[131072];
		if(mode.startsWith("double-")){
			for(int i=0; i<keys.length; i++){keys[i]="key"+i; keys[i].hashCode();}
		}
		run(mode, 8, fixed, keys);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final double checksum=run(mode, rounds, fixed, keys);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		System.out.printf(Locale.ROOT, "%s\t%d\t%.6f\t%.6f\t%.0f%n",
			mode, rounds, elapsed*1e-9, cpu*1e-9, checksum);
	}

	/** Dispatches once per measured run, outside all operation loops. */
	private static double run(String mode, int rounds, boolean fixed, String[] keys){
		assert(rounds>0) : "A benchmark must perform positive fixed work";
		if(mode.equals("double-existing") || mode.equals("double-fresh") || mode.equals("double-growth")){
			return count(mode, rounds, fixed, keys);
		}
		if(mode.equals("bit-increment")){return bits(rounds);}
		if(mode.equals("sequence-fetch")){return sequences(rounds);}
		throw new IllegalArgumentException("Unknown benchmark mode: "+mode);
	}

	private static double count(String mode, int rounds, boolean fixed, String[] keys){
		final boolean existing=mode.equals("double-existing"), growth=mode.equals("double-growth");
		ObjectDoubleMap<String> map=null;
		if(existing){
			map=new ObjectDoubleMap<String>(keys.length*2, String.class);
			for(String key : keys){map.put(key, 0);}
		}
		double checksum=0;
		for(int r=0; r<rounds; r++){
			if(!existing){map=new ObjectDoubleMap<String>(growth ? 256 : keys.length*2, String.class);}
			for(int i=0; i<keys.length; i++){checksum+=map.increment(keys[(i*8191)&(keys.length-1)]);}
			if(map.size()!=keys.length){throw new AssertionError("Missing counted keys");}
		}
		final double expected=existing ? rounds : fixed ? 1 : 0;
		for(String key : keys){if(map.get(key)!=expected){throw new AssertionError("Stored count oracle");}}
		final double expectedSum=existing ? (double)keys.length*rounds*(rounds+1)/2 : (double)keys.length*rounds*expected;
		if((fixed || !growth) && checksum!=expectedSum){throw new AssertionError("Returned count oracle");}
		return checksum;
	}

	/** Ordinary small increments exercise all four counter states without overflow. */
	private static long bits(int rounds){
		final AbstractBitSet set=AbstractBitSet.make(131072, 2);
		long checksum=0;
		for(int r=0; r<rounds; r++){
			set.clear();
			for(int pass=0; pass<4; pass++){
				for(int i=0; i<131072; i++){set.increment(i);}
			}
			checksum+=set.getCount(r&131071);
		}
		if(checksum!=3L*rounds || set.cardinality()!=131072){throw new AssertionError("Packed count oracle");}
		return checksum;
	}

	/** Identical valid query regions have one full-length, first-region-relative candidate. */
	private static long sequences(int rounds){
		final byte[] reference="ACGTTGCATGTCAGTCGATCGTAC".getBytes(StandardCharsets.US_ASCII);
		final byte[] query=new byte[reference.length*2];
		System.arraycopy(reference, 0, query, 0, reference.length);
		System.arraycopy(reference, 0, query, reference.length, reference.length);
		final SeqMap map=new SeqMap(7, 0, 0);
		map.add(reference, 1, 1);
		long checksum=0;
		for(int r=0; r<rounds; r++){
			for(int i=0; i<1024; i++){
				final ArrayList<SeqPosM> list=map.doubleFetch(query, 0, reference.length-1,
					reference.length, query.length-1, reference.length, 0, 0, 0, 1, false);
				if(list==null || list.size()!=1 || list.get(0).pos()!=0){throw new AssertionError("Sequence candidate oracle");}
				checksum+=list.get(0).count;
			}
		}
		if(checksum!=1024L*rounds){throw new AssertionError("Sequence count oracle");}
		return checksum;
	}
}
