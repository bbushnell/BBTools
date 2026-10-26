package test.map;

import map.IntListHashMap;
import map.IntListHashMap3;
import structures.IntList;

/** Finite collision-clear and reuse oracle for IntListHashMap variants.
 * @author Sucrose
 */
public final class IntListHashMapClearTest{

	/** The baseline expectation records the orphaned-cluster bug without accepting it as a contract. */
	public static void main(String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation");
		checkPrimeMap(fixed);
		checkPow2Map(fixed);
		checkPrimeWraparound(fixed);
		checkPow2Wraparound(fixed);
		checkPrimeRemoveBeforeClear(fixed);
		checkPow2RemoveBeforeClear(fixed);
		checkPrimeSentinel();
		checkPow2Sentinel();
		if(fixed){
			checkPrimeSparseWraparound();
			checkPow2SparseWraparound();
		}
		System.out.println("PASS\tfixed="+fixed);
	}

	/** Exercises the sparse prime-modulus clear branch on a cluster of three same-initial-slot keys. */
	private static void checkPrimeMap(final boolean fixed){
		final int[] keys={1, 18, 35};
		final IntListHashMap map=new IntListHashMap(16, 0.9);
		fill(map, keys);
		checkPrimeContents(map, keys, 10);
		map.clear();
		check(map.isEmpty() && map.size()==0 && map.toArray().length==0, "IntListHashMap clear metadata");
		check(occupied(map.keys(), map.invalid())==(fixed ? 0 : 2), "IntListHashMap residual cells");
		map.put(keys[0], 100);
		final boolean stale=map.contains(keys[1]) || map.contains(keys[2]) || map.get(keys[1])!=null || map.get(keys[2])!=null;
		check(stale==!fixed, "IntListHashMap stale collision lookup after reuse: "+stale);
		map.put(keys[1], 200);
		check(map.size()==(fixed ? 2 : 1), "IntListHashMap reused stale key must be tracked exactly once");
		check(hasValue(map.get(keys[0]), 100) &&
				hasValues(map.get(keys[1]), fixed ? new int[]{200} : new int[]{11, 200}), "IntListHashMap reuse values");
		System.out.println("IntListHashMap\tstaleAfterReuse="+stale+"\tsizeAfterReput="+map.size());
	}

	/** Exercises the sparse pow2/hash clear branch on a cluster of three same-initial-slot keys. */
	private static void checkPow2Map(final boolean fixed){
		final int[] keys=pow2Colliders();
		final IntListHashMap3 map=new IntListHashMap3(16, 0.9);
		fill(map, keys);
		checkPow2Contents(map, keys, 10);
		map.clear();
		check(map.isEmpty() && map.size()==0 && map.toArray().length==0, "IntListHashMap3 clear metadata");
		check(occupied(map.keys(), map.invalid())==(fixed ? 0 : 2), "IntListHashMap3 residual cells");
		map.put(keys[0], 100);
		final boolean stale=map.contains(keys[1]) || map.contains(keys[2]) || map.get(keys[1])!=null || map.get(keys[2])!=null;
		check(stale==!fixed, "IntListHashMap3 stale collision lookup after reuse: "+stale);
		map.put(keys[1], 200);
		check(map.size()==(fixed ? 2 : 1), "IntListHashMap3 reused stale key must be tracked exactly once");
		check(hasValue(map.get(keys[0]), 100) &&
				hasValues(map.get(keys[1]), fixed ? new int[]{200} : new int[]{11, 200}), "IntListHashMap3 reuse values");
		System.out.println("IntListHashMap3\tstaleAfterReuse="+stale+"\tsizeAfterReput="+map.size());
	}

	/** Exercises the dense prime-modulus clear branch on a same-slot cluster that wraps. */
	private static void checkPrimeWraparound(final boolean fixed){
		final int[] keys=samePrimeCell(16, 12, 17);
		final IntListHashMap map=new IntListHashMap(16, 0.9);
		fill(map, keys);
		checkPrimeContents(map, keys, 10);
		map.clear();
		check(map.isEmpty() && map.toArray().length==0, "IntListHashMap wrap metadata");
		check(occupied(map.keys(), map.invalid())==(fixed ? 0 : keys.length-1), "IntListHashMap wrap residual cells");
		map.put(keys[0], 100);
		check((map.get(keys[keys.length-1])!=null)==!fixed, "IntListHashMap wrap resurrected tail");
		System.out.println("IntListHashMap\twrapResiduals="+occupied(map.keys(), map.invalid()));
	}

	/** Exercises the dense pow2 clear branch on a same-slot cluster that wraps beyond the array tail's extra cells. */
	private static void checkPow2Wraparound(final boolean fixed){
		final int[] keys=pow2Colliders(14, 15);
		final IntListHashMap3 map=new IntListHashMap3(16, 0.9);
		fill(map, keys);
		checkPow2Contents(map, keys, 10);
		map.clear();
		check(map.isEmpty() && map.toArray().length==0, "IntListHashMap3 wrap metadata");
		check(occupied(map.keys(), map.invalid())==(fixed ? 0 : keys.length-1), "IntListHashMap3 wrap residual cells");
		map.put(keys[0], 100);
		check((map.get(keys[keys.length-1])!=null)==!fixed, "IntListHashMap3 wrap resurrected tail");
		System.out.println("IntListHashMap3\twrapResiduals="+occupied(map.keys(), map.invalid()));
	}

	/** Fixed-only control for the sparse prime-modulus clear branch wrapping at prime257/home256. */
	private static void checkPrimeSparseWraparound(){
		final int[] keys=samePrimeCell(256, 12, 257);
		final IntListHashMap map=new IntListHashMap(256, 0.9);
		fill(map, keys);
		checkPrimeContents(map, keys, 10);
		map.clear();
		check(map.isEmpty() && map.size()==0 && map.toArray().length==0, "IntListHashMap sparse wrap metadata");
		check(occupied(map.keys(), map.invalid())==0, "IntListHashMap sparse wrap residual cells");
		map.put(keys[0], 100);
		check(map.get(keys[keys.length-1])==null, "IntListHashMap sparse wrap resurrected tail");
	}

	/** Fixed-only control for the sparse pow2 clear branch wrapping from mask255/home255. */
	private static void checkPow2SparseWraparound(){
		final int[] keys=pow2Colliders(12, 255, 255);
		final IntListHashMap3 map=new IntListHashMap3(256, 0.9);
		fill(map, keys);
		checkPow2Contents(map, keys, 10);
		map.clear();
		check(map.isEmpty() && map.size()==0 && map.toArray().length==0, "IntListHashMap3 sparse wrap metadata");
		check(occupied(map.keys(), map.invalid())==0, "IntListHashMap3 sparse wrap residual cells");
		map.put(keys[0], 100);
		check(map.get(keys[keys.length-1])==null, "IntListHashMap3 sparse wrap resurrected tail");
	}

	/** Exercises clear after remove has already compacted the tracked list and probe cluster. */
	private static void checkPrimeRemoveBeforeClear(final boolean fixed){
		final int[] keys={1, 18, 35};
		final IntListHashMap map=new IntListHashMap(16, 0.9);
		fill(map, keys);
		check(map.remove(keys[1]), "IntListHashMap remove middle collider");
		map.clear();
		check(map.isEmpty() && map.toArray().length==0, "IntListHashMap remove-before-clear metadata");
		check(occupied(map.keys(), map.invalid())==(fixed ? 0 : 1), "IntListHashMap remove-before-clear residual");
		map.put(keys[0], 100);
		check((map.get(keys[2])!=null)==!fixed, "IntListHashMap remove-before-clear resurrected tail");
	}

	/** Exercises clear after remove has already compacted the tracked list and probe cluster. */
	private static void checkPow2RemoveBeforeClear(final boolean fixed){
		final int[] keys=pow2Colliders();
		final IntListHashMap3 map=new IntListHashMap3(16, 0.9);
		fill(map, keys);
		check(map.remove(keys[1]), "IntListHashMap3 remove middle collider");
		map.clear();
		check(map.isEmpty() && map.toArray().length==0, "IntListHashMap3 remove-before-clear metadata");
		check(occupied(map.keys(), map.invalid())==(fixed ? 0 : 1), "IntListHashMap3 remove-before-clear residual");
		map.put(keys[0], 100);
		check((map.get(keys[2])!=null)==!fixed, "IntListHashMap3 remove-before-clear resurrected tail");
	}

	/** Exercises a real key equal to the initial invalid sentinel. */
	private static void checkPrimeSentinel(){
		final IntListHashMap map=new IntListHashMap(16, 0.9);
		final int sentinel=map.invalid();
		map.put(sentinel, 17);
		check(map.invalid()!=sentinel && hasValue(map.get(sentinel), 17), "IntListHashMap sentinel insertion");
		map.clear();
		check(map.isEmpty() && occupied(map.keys(), map.invalid())==0, "IntListHashMap sentinel clear");
		map.put(sentinel, 19);
		check(hasValue(map.get(sentinel), 19), "IntListHashMap sentinel reuse after clear");
	}

	/** Exercises a real key equal to the initial invalid sentinel. */
	private static void checkPow2Sentinel(){
		final IntListHashMap3 map=new IntListHashMap3(16, 0.9);
		final int sentinel=map.invalid();
		map.put(sentinel, 17);
		check(map.invalid()!=sentinel && hasValue(map.get(sentinel), 17), "IntListHashMap3 sentinel insertion");
		map.clear();
		check(map.isEmpty() && occupied(map.keys(), map.invalid())==0, "IntListHashMap3 sentinel clear");
		map.put(sentinel, 19);
		check(hasValue(map.get(sentinel), 19), "IntListHashMap3 sentinel reuse after clear");
	}

	/** Populates each key with one distinct value. */
	private static void fill(final IntListHashMap map, final int[] keys){
		for(int i=0; i<keys.length; i++){map.put(keys[i], 10+i);}
	}

	/** Populates each key with one distinct value. */
	private static void fill(final IntListHashMap3 map, final int[] keys){
		for(int i=0; i<keys.length; i++){map.put(keys[i], 10+i);}
	}

	/** Confirms prime-map membership before the destructive clear. */
	private static void checkPrimeContents(final IntListHashMap map, final int[] keys, final int firstValue){
		check(map.size()==keys.length && map.toArray().length==keys.length, "IntListHashMap populated size");
		for(int i=0; i<keys.length; i++){
			check(map.contains(keys[i]) && hasValue(map.get(keys[i]), firstValue+i), "IntListHashMap populated key "+i);
		}
	}

	/** Confirms pow2-map membership before the destructive clear. */
	private static void checkPow2Contents(final IntListHashMap3 map, final int[] keys, final int firstValue){
		check(map.size()==keys.length && map.toArray().length==keys.length, "IntListHashMap3 populated size");
		for(int i=0; i<keys.length; i++){
			check(map.contains(keys[i]) && hasValue(map.get(keys[i]), firstValue+i), "IntListHashMap3 populated key "+i);
		}
	}

	/** Finds three keys that collide in a fresh IntListHashMap3(16)'s mask-15 table. */
	private static int[] pow2Colliders(){
		return pow2Colliders(3, 0);
	}

	/** Finds same-slot keys that collide in a fresh IntListHashMap3(16)'s mask-15 table. */
	private static int[] pow2Colliders(final int length, final int home){
		return pow2Colliders(length, home, 15);
	}

	/** Finds same-slot keys that collide in a pow2 table with the specified mask. */
	private static int[] pow2Colliders(final int length, final int home, final int mask){
		check(home>=0 && home<=mask, "Pow2 home slot must fit mask "+mask+": "+home);
		final int[] keys=new int[length];
		int found=0;
		for(int key=0; key<65536 && found<keys.length; key++){
			if((hash(key)&mask)==home){keys[found++]=key;}
		}
		check(found==keys.length, "Expected "+length+" pow2 colliders at "+home+
				" with mask "+mask+" among first 65536 integers; found "+found);
		return keys;
	}

	/** Finds keys in the same initial slot for the prime-modulus map. */
	private static int[] samePrimeCell(final int home, final int length, final int modulus){
		final int[] keys=new int[length];
		for(int i=0; i<length; i++){keys[i]=home+i*modulus;}
		return keys;
	}

	/** Counts raw occupied cells. */
	private static int occupied(final int[] keys, final int invalid){
		int count=0;
		for(int key : keys){
			if(key!=invalid){count++;}
		}
		return count;
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

	/** Returns true when list holds exactly one expected value. */
	private static boolean hasValue(final IntList list, final int expected){
		return list!=null && list.size()==1 && list.get(0)==expected;
	}

	/** Returns true when list exactly matches all expected values. */
	private static boolean hasValues(final IntList list, final int[] expected){
		if(list==null || list.size()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){
			if(list.get(i)!=expected[i]){return false;}
		}
		return true;
	}

	/** Assertions remain meaningful under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
