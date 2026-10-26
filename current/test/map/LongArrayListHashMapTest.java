package test.map;

import java.util.ArrayList;
import map.LongArrayListHashMap;

/** Small correctness oracle for LongArrayListHashMap boxed Long values.
 * @author Jean
 */
public final class LongArrayListHashMapTest{

	/** Exercises Long values, copied input lists, removal and clear without running the legacy benchmark. */
	public static void main(String[] args){
		final LongArrayListHashMap<Long> map=new LongArrayListHashMap<Long>(4, 0.7f);
		check(map.put(7, 14L), "first Long value creates a key");
		check(!map.put(7, 21L), "second Long value appends to the same key");
		check(hasValues(map.get(7), 14L, 21L), "Long values should remain boxed as Long");
		check(!map.get(7).contains(Integer.valueOf(14)), "Integer must not match boxed Long values");

		final ArrayList<Long> input=new ArrayList<Long>();
		input.add(28L);
		input.add(35L);
		check(map.put(-5, input), "caller-owned ArrayList should create a key");
		input.clear();
		check(hasValues(map.get(-5), 28L, 35L), "put(ArrayList) should copy caller-owned values");

		check(map.remove(7), "remove should report existing key");
		check(!map.contains(7) && map.get(7)==null, "removed key should be absent");
		map.clear();
		check(map.isEmpty() && map.size()==0 && map.toArray().length==0, "clear should empty metadata");
		System.out.println("PASS");
	}

	/** Returns true when a list exactly matches all expected Long values. */
	private static boolean hasValues(final ArrayList<Long> list, final long... expected){
		if(list==null || list.size()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){
			if(!Long.valueOf(expected[i]).equals(list.get(i))){return false;}
		}
		return true;
	}

	/** Assertions remain meaningful under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private LongArrayListHashMapTest(){
		throw new AssertionError();
	}
}
