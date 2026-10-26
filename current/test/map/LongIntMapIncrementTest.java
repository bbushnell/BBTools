package test.map;

import java.util.HashMap;
import java.util.Random;
import map.LongIntMap;
import shared.Tools;

/** Independent count and return-value oracles across growth, removal and reuse.
 * @author Shinobu
 */
public final class LongIntMapIncrementTest{

	/** The baseline expectation records the known return defect without accepting it as a contract. */
	public static void main(String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation");
		checkRelocation(fixed);
		checkControls();
		checkMixed(fixed);
		System.out.println("PASS\tfixed="+fixed);
	}

	/** A one-bucket table grows on its first insert; an odd hash moves cell zero to one. */
	private static void checkRelocation(final boolean fixed){
		long key=0;
		while(key<1024 && (Tools.hash64plus(key)&1)==0){key++;}
		check(key<1024, "Fixture must find a key that relocates on first growth");
		final LongIntMap map=new LongIntMap(1);
		final int capacity=map.keys().length;
		final int returned=map.increment(key, 7);
		check(map.keys().length>capacity, "Fixture must actually grow the table");
		check(map.keys()[1]==key && map.values()[1]==7, "Inserted key must move to cell one");
		check(map.get(key)==7 && map.size()==1, "Growth must preserve the inserted count");
		check(returned==(fixed ? 7 : 0), "Increment must return the new count after relocation: "+returned);
		System.out.println("relocation\treturned="+returned+"\tstored="+map.get(key));
	}

	/** Preserves absent sentinels, arbitrary signed values, positive saturation and bulk operations. */
	private static void checkControls(){
		final LongIntMap map=new LongIntMap(64);
		check(map.get(10)==-1 && map.getOrDefault(10, 43)==43 && !map.contains(10), "Absent lookup");
		check(map.put(10, -1)==0 && map.contains(10), "Absent put returns zero, even for stored -1");
		check(map.getOrDefault(10, 43)==-1 && map.increment(10)==0, "Present -1 is a real value");
		check(map.increment(11, -3)==-3, "Absent negative increment");
		map.put(12, Integer.MAX_VALUE-1);
		check(map.increment(12, 7)==Integer.MAX_VALUE, "Positive addition saturates");
		check(map.increment(12, -2)==Integer.MAX_VALUE-2, "Negative increment from saturated count");
		final long sentinel=map.invalid();
		check(map.put(sentinel, Integer.MIN_VALUE)==0, "All long keys, including the old sentinel, are supported");
		check(map.invalid()!=sentinel && map.getOrDefault(sentinel, 0)==Integer.MIN_VALUE, "Sentinel replacement retains signed bits");
		check(map.remove(11) && !map.remove(11), "Removal reports membership");
		check(map.increment(11, 5)==5, "Removed slot starts from zero");
		final LongIntMap copy=new LongIntMap(1);
		copy.putAll(map);
		for(long key : map.toArray()){check(copy.get(key)==map.get(key), "putAll preserves signed values");}
		map.clear();
		check(map.isEmpty() && map.increment(11)==1, "Clear resets count origin");
		final LongIntMap source=new LongIntMap(2);
		source.put(11, 2); source.put(15, 4);
		map.incrementAll(source);
		check(map.get(11)==3 && map.get(15)==4, "Bulk addition preserves stored counts");
		System.out.println("controls\tabsent=-1\tabsentPut=0\tsignedValues=true\tsaturation=true");
	}

	/** A boxed reference map supplies values independently of probe positions and resizing. */
	private static void checkMixed(final boolean fixed){
		final LongIntMap map=new LongIntMap(1);
		final HashMap<Long, Integer> expected=new HashMap<Long, Integer>();
		final Random random=new Random(714);
		int badReturns=0;
		for(int i=0; i<20000; i++){
			final long key=random.nextInt(2048)-1024;
			final int operation=random.nextInt(8);
			if(operation==0){
				check(map.remove(key)==(expected.remove(key)!=null), "Removal membership at "+i);
			}else if(operation==1){
				final int value=random.nextInt(101)-50;
				final Integer previous=expected.put(key, value);
				check(map.put(key, value)==(previous==null ? 0 : previous), "Put previous value at "+i);
			}else{
				final int amount=random.nextInt(11)-5;
				final Integer previous=expected.get(key);
				final int value=(previous==null ? 0 : previous)+amount;
				expected.put(key, value);
				final int returned=map.increment(key, amount);
				if(returned!=value){badReturns++;}
				check(!fixed || returned==value, "Increment return at "+i+": "+returned+" != "+value);
				check(map.get(key)==value, "Stored value at "+i);
			}
			check(map.size()==expected.size(), "Map size after operation "+i);
		}
		for(long key=-1024; key<1024; key++){
			final Integer value=expected.get(key);
			check(map.contains(key)==(value!=null), "Final membership for "+key);
			check(map.getOrDefault(key, 123456)==(value==null ? 123456 : value), "Final stored value for "+key);
		}
		System.out.println("mixed\toperations=20000\tkeys="+map.size()+"\tincorrectReturns="+badReturns);
	}

	/** Assertions remain meaningful under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
