package test.map;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import map.IntHashMap;

/** Regression test for IntHashMap lookups of the live empty-cell sentinel.
 * @author Jean
 */
public final class IntHashMapLookupTest{

	/** Exercises the baseline failure and the fixed behavior. */
	public static void main(String[] args){
		final boolean fixed=parseExpectation(args);
		final IntHashMap map=new IntHashMap(4);

		checkAbsent(map, 17, "ordinary absent key in empty map");
		checkCurrentSentinel(map, fixed, "empty map");

		map.put(17, 29);
		checkPresent(map, 17, 29, "ordinary inserted key");
		checkAbsent(map, 19, "ordinary absent key in populated map");
		checkCurrentSentinel(map, fixed, "populated map");

		final int insertedSentinel=map.invalid();
		map.put(insertedSentinel, 31);
		if(map.invalid()==insertedSentinel){
			throw new AssertionError("Putting the live sentinel did not choose a new sentinel");
		}
		checkPresent(map, insertedSentinel, 31, "key inserted after sentinel reset");
		checkPresent(map, 17, 29, "ordinary key after sentinel reset");
		checkAbsent(map, 19, "ordinary absent key after sentinel reset");
		checkCurrentSentinel(map, fixed, "after sentinel reset");

		map.clear();
		if(!map.isEmpty() || map.size()!=0){
			throw new AssertionError("clear() did not empty the map");
		}
		checkAbsent(map, insertedSentinel, "former sentinel key after clear");
		checkAbsent(map, 17, "ordinary key after clear");
		checkCurrentSentinel(map, fixed, "after clear");
		checkCollisionWrap(fixed);
		checkMixed(fixed);
	}

	/** Forces a probe chain across the backing-array boundary, including signed keys and removal. */
	private static void checkCollisionWrap(final boolean fixed){
		final IntHashMap map=new IntHashMap(127);
		final int[] keys=new int[32];
		for(int i=0; i<keys.length; i++){
			final int positive=126+127*i;
			keys[i]=(i&1)==0 ? positive : positive|Integer.MIN_VALUE;
			map.put(keys[i], -i-1);
		}
		//The constructor's prime modulus is127; 32 entries exceed the tail at home cell126.
		if(map.findCell(keys[31])>=126){throw new AssertionError("Fixture did not wrap its probe chain");}
		for(int i=0; i<keys.length; i++){checkPresent(map, keys[i], -i-1, "wrapped chain");}
		for(int i=0; i<keys.length; i+=3){map.remove(keys[i]);}
		for(int i=0; i<keys.length; i++){
			if(i%3==0){checkAbsent(map, keys[i], "removed from wrapped chain");}
			else{checkPresent(map, keys[i], -i-1, "rehash after wrapped removal");}
		}
		checkAbsent(map, 126+127*63, "missing key after wrapped chain");
		checkCurrentSentinel(map, fixed, "wrapped chain");
	}

	/** Checks signed keys, growth, removal and sentinel changes against an independent map. */
	private static void checkMixed(final boolean fixed){
		final IntHashMap map=new IntHashMap(1);
		final HashMap<Integer, Integer> oracle=new HashMap<Integer, Integer>();
		for(int key : new int[]{Integer.MIN_VALUE, Integer.MAX_VALUE, -1, 0, 1}){
			map.put(key, key); oracle.put(key, key);
		}
		final Random random=new Random(7231);
		for(int i=0; i<20000; i++){
			final int key=i%97==0 ? map.invalid() : (random.nextInt(1024)-512)*127;
			if((i&3)==0){
				final boolean expected=oracle.remove(key)!=null;
				if(map.remove(key)!=expected){throw new AssertionError("Removal mismatch at "+i);}
			}else{
				final int value=random.nextInt();
				map.put(key, value);
				oracle.put(key, value);
			}
			if(map.size()!=oracle.size()){throw new AssertionError("Size mismatch at "+i);}
			if(i%251==0){
				for(Map.Entry<Integer, Integer> entry : oracle.entrySet()){
					checkPresent(map, entry.getKey(), entry.getValue(), "mixed operations");
				}
				checkCurrentSentinel(map, fixed, "mixed operations");
			}
		}
		System.out.println("PASS mixed lookup operations=20000 keys="+map.size());
	}

	/** Parses expect=fixed or expect=baseline. */
	private static boolean parseExpectation(String[] args){
		Boolean fixed=null;
		for(String arg : args){
			if(arg.equals("expect=fixed")){fixed=Boolean.TRUE;}
			else if(arg.equals("expect=baseline")){fixed=Boolean.FALSE;}
			else{throw new IllegalArgumentException("Expected only expect=fixed or expect=baseline: "+arg);}
		}
		if(fixed==null){throw new IllegalArgumentException("Missing expect=fixed or expect=baseline");}
		return fixed.booleanValue();
	}

	/** Checks that a normal absent key is never found. */
	private static void checkAbsent(final IntHashMap map, final int key, final String label){
		if(map.contains(key)){throw new AssertionError(label+" unexpectedly present: "+key);}
		final int value=map.get(key);
		if(value!=-1){throw new AssertionError(label+" get() returned "+value+" instead of -1");}
	}

	/** Checks that an inserted key remains reachable with the expected value. */
	private static void checkPresent(final IntHashMap map, final int key, final int expected, final String label){
		if(!map.contains(key)){throw new AssertionError(label+" missing: "+key);}
		final int value=map.get(key);
		if(value!=expected){throw new AssertionError(label+" value "+value+" != "+expected);}
	}

	/** Checks the live sentinel query in either known-baseline or fixed mode. */
	private static void checkCurrentSentinel(final IntHashMap map, final boolean fixed, final String label){
		final int sentinel=map.invalid();
		final boolean contains=map.contains(sentinel);
		final int value=map.get(sentinel);
		if(contains==fixed){
			throw new AssertionError(label+" contains(invalid)="+contains+" for sentinel "+sentinel);
		}
		final int expected=fixed ? -1 : 0;
		if(value!=expected){
			throw new AssertionError(label+" get(invalid)="+value+" != "+expected+" for sentinel "+sentinel);
		}
	}

	private IntHashMapLookupTest(){
		throw new AssertionError();
	}
}
