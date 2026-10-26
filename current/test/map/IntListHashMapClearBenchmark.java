package test.map;

import java.lang.management.ManagementFactory;
import java.util.Locale;
import map.IntListHashMap;
import map.IntListHashMap3;
import parse.Parse;
import structures.IntList;

/** Fixed-work no-collision timing and collision validation for IntListHashMap clear.
 * @author Sucrose
 */
public final class IntListHashMapClearBenchmark{

	/** Benchmarks one map type or both IntListHashMap clear siblings. */
	public static void main(String[] args){
		String type="all", mode="nocollide";
		int size=65535, rounds=256;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=arg.substring(0, equals).toLowerCase(Locale.ROOT), b=arg.substring(equals+1);
			if(a.equals("type")){type=b.toLowerCase(Locale.ROOT);}
			else if(a.equals("mode")){mode=b.toLowerCase(Locale.ROOT);}
			else if(a.equals("size")){size=Parse.parseIntKMG(b);}
			else if(a.equals("rounds")){rounds=Parse.parseIntKMG(b);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(!type.equals("all") && !type.equals("prime") && !type.equals("pow2")){
			throw new IllegalArgumentException("Unknown type: "+type);
		}
		if(!mode.equals("nocollide") && !mode.equals("collision")){
			throw new IllegalArgumentException("Unknown mode: "+mode);
		}
		if(size<3 || size>1<<20 || rounds<1 || rounds>1000000){
			throw new IllegalArgumentException("Require size in [3,1048576] and rounds in [1,1000000]");
		}
		if(mode.equals("collision") && size%3!=0){
			throw new IllegalArgumentException("Collision mode requires size to be divisible by 3");
		}
		if(type.equals("all") || type.equals("prime")){benchPrime(mode, size, rounds);}
		if(type.equals("all") || type.equals("pow2")){benchPow2(mode, size, rounds);}
	}

	/** Runs a warmup, a timed fixed workload and one tabular report row. */
	private static void benchPrime(final String mode, final int size, final int rounds){
		final IntListHashMap map=new IntListHashMap(size*4, 0.9);
		final int[] keys=mode.equals("collision") ? primeColliders(map, size) : primeNoColliders(size);
		run(map, keys, 16);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long checksum=run(map, keys, rounds);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		System.out.printf(Locale.ROOT, "prime\t%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d%n",
			mode, keys.length, rounds, (long)keys.length*rounds, elapsed*1e-9, cpu*1e-9, checksum);
	}

	/** Runs a warmup, a timed fixed workload and one tabular report row. */
	private static void benchPow2(final String mode, final int size, final int rounds){
		final IntListHashMap3 map=new IntListHashMap3(size*4, 0.9);
		final int mask=map.keys().length-11;
		final int[] keys=mode.equals("collision") ? pow2Colliders(size, mask) : pow2NoColliders(size, mask);
		run(map, keys, 16);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long checksum=run(map, keys, rounds);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		System.out.printf(Locale.ROOT, "pow2\t%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d%n",
			mode, keys.length, rounds, (long)keys.length*rounds, elapsed*1e-9, cpu*1e-9, checksum);
	}

	/** Populates and clears the same IntListHashMap for a fixed number of rounds. */
	private static long run(final IntListHashMap map, final int[] keys, final int rounds){
		assert(keys.length>0 && rounds>0) : "Benchmark must execute a positive fixed workload";
		long checksum=0;
		for(int r=0; r<rounds; r++){
			for(int key : keys){map.put(key, r);}
			check(map.size()==keys.length, "IntListHashMap clear/reput size drift: "+map.size());
			checkValue(map.get(keys[0]), r, "IntListHashMap first key");
			checkValue(map.get(keys[keys.length-1]), r, "IntListHashMap last key");
			checksum+=map.size();
			map.clear();
			check(map.isEmpty() && map.toArray().length==0, "IntListHashMap clear metadata");
		}
		return checksum;
	}

	/** Populates and clears the same IntListHashMap3 for a fixed number of rounds. */
	private static long run(final IntListHashMap3 map, final int[] keys, final int rounds){
		assert(keys.length>0 && rounds>0) : "Benchmark must execute a positive fixed workload";
		long checksum=0;
		for(int r=0; r<rounds; r++){
			for(int key : keys){map.put(key, r);}
			check(map.size()==keys.length, "IntListHashMap3 clear/reput size drift: "+map.size());
			checkValue(map.get(keys[0]), r, "IntListHashMap3 first key");
			checkValue(map.get(keys[keys.length-1]), r, "IntListHashMap3 last key");
			checksum+=map.size();
			map.clear();
			check(map.isEmpty() && map.toArray().length==0, "IntListHashMap3 clear metadata");
		}
		return checksum;
	}

	/** Generates keys with unique initial slots in a prime-modulus map. */
	private static int[] primeNoColliders(final int size){
		final int[] keys=new int[size];
		for(int i=0; i<size; i++){keys[i]=i;}
		return keys;
	}

	/** Generates three-key same-slot groups in a prime-modulus map. */
	private static int[] primeColliders(final IntListHashMap map, final int size){
		final int modulus=map.keys().length-10, groups=size/3;
		final int[] keys=new int[size];
		for(int i=0; i<groups; i++){
			final int initial=4*i;
			keys[3*i]=initial;
			keys[3*i+1]=initial+modulus;
			keys[3*i+2]=initial+2*modulus;
		}
		return keys;
	}

	/** Generates keys with unique initial slots in an IntListHashMap3 table. */
	private static int[] pow2NoColliders(final int size, final int mask){
		final boolean[] used=new boolean[mask+1];
		final int[] keys=new int[size];
		final int limit=Math.max(1024, size*256);
		int found=0;
		for(int key=0; key<limit && found<size; key++){
			final int cell=hash(key)&mask;
			if(!used[cell]){
				used[cell]=true;
				keys[found++]=key;
			}
		}
		check(found==size, "Expected "+size+" no-collision pow2 keys under mask "+mask+"; found "+found);
		return keys;
	}

	/** Generates three-key same-slot groups in an IntListHashMap3 table. */
	private static int[] pow2Colliders(final int size, final int mask){
		final int[] counts=new int[mask+1], pending=new int[(mask+1)*3], keys=new int[size];
		final int limit=Math.max(1024, size*256);
		int found=0;
		for(int key=0; key<limit && found<size; key++){
			final int cell=hash(key)&mask, count=counts[cell];
			if((cell&3)==0 && count<3){
				pending[cell*3+count]=key;
				counts[cell]=count+1;
				if(count==2){
					final int offset=cell*3;
					keys[found++]=pending[offset];
					keys[found++]=pending[offset+1];
					keys[found++]=pending[offset+2];
				}
			}
		}
		check(found==size, "Expected "+size+" colliding pow2 keys under mask "+mask+"; found "+found);
		return keys;
	}

	/** IntListHashMap3's non-negative integer mixer. */
	private static int hash(int key){
		key=~key+(key<<15);
		key=key^(key>>>12);
		key=key+(key<<2);
		key=key^(key>>>4);
		key=key*2057;
		key=key^(key>>>16);
		key=key&(0x7FFFFFFF);
		return key<0x7FFFF800 ? key : (key-0x7FFFF800)*64;
	}

	/** Confirms that a freshly inserted list holds one expected value. */
	private static void checkValue(final IntList list, final int expected, final String label){
		check(list!=null && list.size()==1 && list.get(0)==expected, label+" value mismatch");
	}

	/** Assertions remain meaningful under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
