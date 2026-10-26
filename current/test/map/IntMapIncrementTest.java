package test.map;

import java.util.HashMap;
import java.util.Random;
import map.IntHashMap2;
import map.IntHashMap3;
import map.IntLongHashMap2;
import map.IntLongHashMap3;
import shared.Tools;

/** Independent return-value and storage oracles for int-key increment maps.
 * @author Sucrose
 */
public final class IntMapIncrementTest{

	/** The baseline expectation records the known return defect without accepting it as a contract. */
	public static void main(String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation");
		checkIntHashMap2(fixed);
		checkIntHashMap3(fixed);
		checkIntLongHashMap2(fixed);
		checkIntLongHashMap3(fixed);
		System.out.println("PASS	fixed="+fixed);
	}

	/** Exercises IntHashMap2 first-insert relocation, ordinary mutation and mixed growth. */
	private static void checkIntHashMap2(final boolean fixed){
		final int key=firstHash32plusOdd();
		final IntHashMap2 map=new IntHashMap2(1);
		final int capacity=map.keys().length;
		final int returned=map.increment(key, 7);
		check(map.keys().length>capacity, "IntHashMap2 fixture must grow the table");
		check(map.keys()[1]==key && map.values()[1]==7, "IntHashMap2 inserted key must move to cell one");
		check(map.get(key)==7 && map.size()==1, "IntHashMap2 growth must preserve the inserted count");
		check(returned==(fixed ? 7 : 0), "IntHashMap2 increment return after relocation: "+returned);

		final IntHashMap2 controls=new IntHashMap2(64);
		check(controls.get(10)==-1 && !controls.contains(10), "IntHashMap2 absent lookup");
		check(controls.put(10, -1)==0 && controls.contains(10), "IntHashMap2 absent put");
		check(controls.increment(10)==0, "IntHashMap2 present -1 is a real value");
		check(controls.increment(11, -3)==-3, "IntHashMap2 absent negative increment");
		controls.put(12, Integer.MAX_VALUE-1);
		check(controls.increment(12, 7)==Integer.MAX_VALUE, "IntHashMap2 positive saturation");
		check(controls.increment(12, -2)==Integer.MAX_VALUE-2, "IntHashMap2 negative increment after saturation");
		check(controls.remove(11) && !controls.remove(11), "IntHashMap2 remove membership");
		check(controls.increment(11, 5)==5, "IntHashMap2 removed slot starts from zero");
		final IntHashMap2 copy=new IntHashMap2(1);
		copy.putAll(controls);
		for(int storedKey : controls.toArray()){check(copy.get(storedKey)==controls.get(storedKey), "IntHashMap2 putAll");}
		controls.clear();
		check(controls.isEmpty() && controls.increment(11)==1, "IntHashMap2 clear resets count origin");
		final IntHashMap2 source=new IntHashMap2(2);
		source.put(11, 2); source.put(15, 4);
		controls.incrementAll(source);
		check(controls.get(11)==3 && controls.get(15)==4, "IntHashMap2 incrementAll");

		final int badReturns=checkIntHashMap2Mixed(fixed);
		check(fixed ? badReturns==0 : badReturns>0, "IntHashMap2 mixed bad returns: "+badReturns);
		System.out.println("IntHashMap2	relocation="+returned+"	mixedBadReturns="+badReturns);
	}

	/** Exercises IntHashMap3 first-insert relocation, ordinary mutation and mixed growth. */
	private static void checkIntHashMap3(final boolean fixed){
		final int key=firstHash32plusOdd();
		final IntHashMap3 map=new IntHashMap3(1);
		final int capacity=map.keys().length;
		final int returned=map.increment(key, 7);
		check(map.keys().length>capacity, "IntHashMap3 fixture must grow the table");
		check(map.keys()[0]!=key, "IntHashMap3 inserted key must move away from cell zero");
		final int cell=findCell(map.keys(), map.invalid(), key);
		check(cell>0 && map.values()[cell]==7, "IntHashMap3 relocated cell must store 7");
		check(map.get(key)==7 && map.size()==1, "IntHashMap3 growth must preserve the inserted count");
		check(returned==(fixed ? 7 : 0), "IntHashMap3 increment return after relocation: "+returned);

		final IntHashMap3 controls=new IntHashMap3(64);
		check(controls.get(10)==-1 && !controls.contains(10), "IntHashMap3 absent lookup");
		check(controls.put(10, -1)==0 && controls.contains(10), "IntHashMap3 absent put");
		check(controls.increment(10)==0, "IntHashMap3 present -1 is a real value");
		check(controls.increment(11, -3)==-3, "IntHashMap3 absent negative increment");
		controls.put(12, Integer.MAX_VALUE-1);
		check(controls.increment(12, 7)==Integer.MAX_VALUE, "IntHashMap3 positive saturation");
		check(controls.increment(12, -2)==Integer.MAX_VALUE-2, "IntHashMap3 negative increment after saturation");
		check(controls.remove(11) && !controls.remove(11), "IntHashMap3 remove membership");
		check(controls.increment(11, 5)==5, "IntHashMap3 removed slot starts from zero");
		final IntHashMap3 copy=new IntHashMap3(1);
		copy.putAll(controls);
		for(int storedKey : controls.toArray()){check(copy.get(storedKey)==controls.get(storedKey), "IntHashMap3 putAll");}
		controls.clear();
		check(controls.isEmpty() && controls.increment(11)==1, "IntHashMap3 clear resets count origin");
		final IntHashMap3 source=new IntHashMap3(2);
		source.put(11, 2); source.put(15, 4);
		controls.incrementAll(source);
		check(controls.get(11)==3 && controls.get(15)==4, "IntHashMap3 incrementAll");

		final int badReturns=checkIntHashMap3Mixed(fixed);
		check(fixed ? badReturns==0 : badReturns>0, "IntHashMap3 mixed bad returns: "+badReturns);
		System.out.println("IntHashMap3	relocation="+returned+"	mixedBadReturns="+badReturns);
	}

	/** Exercises IntLongHashMap2 first-insert relocation, getAndIncrement and mixed growth. */
	private static void checkIntLongHashMap2(final boolean fixed){
		final int key=firstHash32shiftOdd();
		final IntLongHashMap2 map=new IntLongHashMap2(1);
		final int capacity=map.keys().length;
		final long returned=map.increment(key, 7);
		check(map.keys().length>capacity, "IntLongHashMap2 fixture must grow the table");
		check(map.keys()[1]==key && map.values()[1]==7, "IntLongHashMap2 inserted key must move to cell one");
		check(map.get(key)==7 && map.size()==1, "IntLongHashMap2 growth must preserve the inserted count");
		check(returned==(fixed ? 7 : 0), "IntLongHashMap2 increment return after relocation: "+returned);

		final IntLongHashMap2 controls=new IntLongHashMap2(64);
		check(controls.get(10)==-1 && !controls.contains(10), "IntLongHashMap2 absent lookup");
		check(controls.put(10, -1)==0 && controls.contains(10), "IntLongHashMap2 absent put");
		check(controls.increment(10)==0, "IntLongHashMap2 present -1 is a real value");
		check(controls.increment(11, -3)==-3, "IntLongHashMap2 absent negative increment");
		check(controls.put(12, 1234567890123L)==0, "IntLongHashMap2 stores long values");
		check(controls.increment(12, 7)==1234567890130L, "IntLongHashMap2 long increment");
		check(controls.remove(11) && !controls.remove(11), "IntLongHashMap2 remove membership");
		check(controls.increment(11, 5)==5, "IntLongHashMap2 removed slot starts from zero");
		final IntLongHashMap2 grow=new IntLongHashMap2(1);
		check(grow.getAndIncrement(key)==0 && grow.get(key)==1 && grow.size()==1, "IntLongHashMap2 getAndIncrement");
		final IntLongHashMap2 copy=new IntLongHashMap2(1);
		copy.putAll(controls);
		for(int storedKey : controls.toArray()){check(copy.get(storedKey)==controls.get(storedKey), "IntLongHashMap2 putAll");}
		controls.clear();
		check(controls.isEmpty() && controls.increment(11)==1, "IntLongHashMap2 clear resets count origin");
		final IntLongHashMap2 source=new IntLongHashMap2(2);
		source.put(11, 2); source.put(15, 4);
		controls.incrementAll(source);
		check(controls.get(11)==3 && controls.get(15)==4, "IntLongHashMap2 incrementAll");

		final int badReturns=checkIntLongHashMap2Mixed(fixed);
		check(fixed ? badReturns==0 : badReturns>0, "IntLongHashMap2 mixed bad returns: "+badReturns);
		System.out.println("IntLongHashMap2	relocation="+returned+"	mixedBadReturns="+badReturns);
	}

	/** Exercises IntLongHashMap3 first-insert relocation, getAndIncrement and mixed growth. */
	private static void checkIntLongHashMap3(final boolean fixed){
		final int key=firstHash32plusOdd();
		final IntLongHashMap3 map=new IntLongHashMap3(1);
		final int capacity=map.keys().length;
		final long returned=map.increment(key, 7);
		check(map.keys().length>capacity, "IntLongHashMap3 fixture must grow the table");
		check(map.keys()[0]!=key, "IntLongHashMap3 inserted key must move away from cell zero");
		final int cell=findCell(map.keys(), map.invalid(), key);
		check(cell>0 && map.values()[cell]==7, "IntLongHashMap3 relocated cell must store 7");
		check(map.get(key)==7 && map.size()==1, "IntLongHashMap3 growth must preserve the inserted count");
		check(returned==(fixed ? 7 : 0), "IntLongHashMap3 increment return after relocation: "+returned);

		final IntLongHashMap3 controls=new IntLongHashMap3(64);
		check(controls.get(10)==-1 && !controls.contains(10), "IntLongHashMap3 absent lookup");
		check(controls.put(10, -1)==0 && controls.contains(10), "IntLongHashMap3 absent put");
		check(controls.increment(10)==0, "IntLongHashMap3 present -1 is a real value");
		check(controls.increment(11, -3)==-3, "IntLongHashMap3 absent negative increment");
		check(controls.put(12, 1234567890123L)==0, "IntLongHashMap3 stores long values");
		check(controls.increment(12, 7)==1234567890130L, "IntLongHashMap3 long increment");
		check(controls.remove(11) && !controls.remove(11), "IntLongHashMap3 remove membership");
		check(controls.increment(11, 5)==5, "IntLongHashMap3 removed slot starts from zero");
		final IntLongHashMap3 grow=new IntLongHashMap3(1);
		check(grow.getAndIncrement(key)==0 && grow.get(key)==1 && grow.size()==1, "IntLongHashMap3 getAndIncrement");
		final IntLongHashMap3 copy=new IntLongHashMap3(1);
		copy.putAll(controls);
		for(int storedKey : controls.toArray()){check(copy.get(storedKey)==controls.get(storedKey), "IntLongHashMap3 putAll");}
		controls.clear();
		check(controls.isEmpty() && controls.increment(11)==1, "IntLongHashMap3 clear resets count origin");
		final IntLongHashMap3 source=new IntLongHashMap3(2);
		source.put(11, 2); source.put(15, 4);
		controls.incrementAll(source);
		check(controls.get(11)==3 && controls.get(15)==4, "IntLongHashMap3 incrementAll");

		final int badReturns=checkIntLongHashMap3Mixed(fixed);
		check(fixed ? badReturns==0 : badReturns>0, "IntLongHashMap3 mixed bad returns: "+badReturns);
		System.out.println("IntLongHashMap3	relocation="+returned+"	mixedBadReturns="+badReturns);
	}

	/** Boxed reference model for IntHashMap2. */
	private static int checkIntHashMap2Mixed(final boolean fixed){
		final IntHashMap2 map=new IntHashMap2(1);
		final HashMap<Integer, Integer> expected=new HashMap<Integer, Integer>();
		final Random random=new Random(2351);
		int badReturns=0;
		for(int i=0; i<20000; i++){
			final int key=random.nextInt(4096)-2048;
			final int operation=random.nextInt(8);
			if(operation==0){
				check(map.remove(key)==(expected.remove(key)!=null), "IntHashMap2 remove at "+i);
			}else if(operation==1){
				final int value=random.nextInt(101)-50;
				final Integer previous=expected.put(key, value);
				check(map.put(key, value)==(previous==null ? 0 : previous), "IntHashMap2 put at "+i);
			}else{
				final int amount=random.nextInt(11)-5;
				final Integer previous=expected.get(key);
				final int value=(previous==null ? 0 : previous)+amount;
				expected.put(key, value);
				final int returned=map.increment(key, amount);
				if(returned!=value){badReturns++;}
				check(!fixed || returned==value, "IntHashMap2 return at "+i);
				check(map.get(key)==value, "IntHashMap2 stored value at "+i);
			}
			check(map.size()==expected.size(), "IntHashMap2 size at "+i);
		}
		checkIntFinal(map, expected, "IntHashMap2");
		return badReturns;
	}

	/** Boxed reference model for IntHashMap3. */
	private static int checkIntHashMap3Mixed(final boolean fixed){
		final IntHashMap3 map=new IntHashMap3(1);
		final HashMap<Integer, Integer> expected=new HashMap<Integer, Integer>();
		final Random random=new Random(2351);
		int badReturns=0;
		for(int i=0; i<20000; i++){
			final int key=random.nextInt(4096)-2048;
			final int operation=random.nextInt(8);
			if(operation==0){
				check(map.remove(key)==(expected.remove(key)!=null), "IntHashMap3 remove at "+i);
			}else if(operation==1){
				final int value=random.nextInt(101)-50;
				final Integer previous=expected.put(key, value);
				check(map.put(key, value)==(previous==null ? 0 : previous), "IntHashMap3 put at "+i);
			}else{
				final int amount=random.nextInt(11)-5;
				final Integer previous=expected.get(key);
				final int value=(previous==null ? 0 : previous)+amount;
				expected.put(key, value);
				final int returned=map.increment(key, amount);
				if(returned!=value){badReturns++;}
				check(!fixed || returned==value, "IntHashMap3 return at "+i);
				check(map.get(key)==value, "IntHashMap3 stored value at "+i);
			}
			check(map.size()==expected.size(), "IntHashMap3 size at "+i);
		}
		checkIntFinal(map, expected, "IntHashMap3");
		return badReturns;
	}

	/** Boxed reference model for IntLongHashMap2. */
	private static int checkIntLongHashMap2Mixed(final boolean fixed){
		final IntLongHashMap2 map=new IntLongHashMap2(1);
		final HashMap<Integer, Long> expected=new HashMap<Integer, Long>();
		final Random random=new Random(6971);
		int badReturns=0;
		for(int i=0; i<20000; i++){
			final int key=random.nextInt(4096)-2048;
			final int operation=random.nextInt(8);
			if(operation==0){
				check(map.remove(key)==(expected.remove(key)!=null), "IntLongHashMap2 remove at "+i);
			}else if(operation==1){
				final long value=random.nextInt(101)-50;
				final Long previous=expected.put(key, value);
				check(map.put(key, value)==(previous==null ? 0 : previous), "IntLongHashMap2 put at "+i);
			}else{
				final long amount=random.nextInt(11)-5;
				final Long previous=expected.get(key);
				final long value=(previous==null ? 0 : previous)+amount;
				expected.put(key, value);
				final long returned=map.increment(key, amount);
				if(returned!=value){badReturns++;}
				check(!fixed || returned==value, "IntLongHashMap2 return at "+i);
				check(map.get(key)==value, "IntLongHashMap2 stored value at "+i);
			}
			check(map.size()==expected.size(), "IntLongHashMap2 size at "+i);
		}
		checkLongFinal(map, expected, "IntLongHashMap2");
		return badReturns;
	}

	/** Boxed reference model for IntLongHashMap3. */
	private static int checkIntLongHashMap3Mixed(final boolean fixed){
		final IntLongHashMap3 map=new IntLongHashMap3(1);
		final HashMap<Integer, Long> expected=new HashMap<Integer, Long>();
		final Random random=new Random(6971);
		int badReturns=0;
		for(int i=0; i<20000; i++){
			final int key=random.nextInt(4096)-2048;
			final int operation=random.nextInt(8);
			if(operation==0){
				check(map.remove(key)==(expected.remove(key)!=null), "IntLongHashMap3 remove at "+i);
			}else if(operation==1){
				final long value=random.nextInt(101)-50;
				final Long previous=expected.put(key, value);
				check(map.put(key, value)==(previous==null ? 0 : previous), "IntLongHashMap3 put at "+i);
			}else{
				final long amount=random.nextInt(11)-5;
				final Long previous=expected.get(key);
				final long value=(previous==null ? 0 : previous)+amount;
				expected.put(key, value);
				final long returned=map.increment(key, amount);
				if(returned!=value){badReturns++;}
				check(!fixed || returned==value, "IntLongHashMap3 return at "+i);
				check(map.get(key)==value, "IntLongHashMap3 stored value at "+i);
			}
			check(map.size()==expected.size(), "IntLongHashMap3 size at "+i);
		}
		checkLongFinal(map, expected, "IntLongHashMap3");
		return badReturns;
	}

	/** Confirms final IntHashMap2 membership and stored values. */
	private static void checkIntFinal(final IntHashMap2 map, final HashMap<Integer, Integer> expected, final String label){
		for(int key=-2048; key<2048; key++){
			final Integer value=expected.get(key);
			check(map.contains(key)==(value!=null), label+" final membership for "+key);
			check(map.get(key)==(value==null ? -1 : value), label+" final value for "+key);
		}
	}

	/** Confirms final IntHashMap3 membership and stored values. */
	private static void checkIntFinal(final IntHashMap3 map, final HashMap<Integer, Integer> expected, final String label){
		for(int key=-2048; key<2048; key++){
			final Integer value=expected.get(key);
			check(map.contains(key)==(value!=null), label+" final membership for "+key);
			check(map.get(key)==(value==null ? -1 : value), label+" final value for "+key);
		}
	}

	/** Confirms final IntLongHashMap2 membership and stored values. */
	private static void checkLongFinal(final IntLongHashMap2 map, final HashMap<Integer, Long> expected, final String label){
		for(int key=-2048; key<2048; key++){
			final Long value=expected.get(key);
			check(map.contains(key)==(value!=null), label+" final membership for "+key);
			check(map.get(key)==(value==null ? -1 : value), label+" final value for "+key);
		}
	}

	/** Confirms final IntLongHashMap3 membership and stored values. */
	private static void checkLongFinal(final IntLongHashMap3 map, final HashMap<Integer, Long> expected, final String label){
		for(int key=-2048; key<2048; key++){
			final Long value=expected.get(key);
			check(map.contains(key)==(value!=null), label+" final membership for "+key);
			check(map.get(key)==(value==null ? -1 : value), label+" final value for "+key);
		}
	}

	/** Finds a non-sentinel key that starts in cell zero and moves after first growth. */
	private static int firstHash32plusOdd(){
		int key=0;
		while(key<1024 && (Tools.hash32plus(key)&1)==0){key++;}
		check(key<1024, "Fixture must find a hash32plus odd key");
		return key;
	}

	/** Locates a known live key in an exposed backing array. */
	private static int findCell(final int[] keys, final int invalid, final int key){
		for(int cell=0; cell<keys.length; cell++){
			if(keys[cell]!=invalid && keys[cell]==key){return cell;}
		}
		return -1;
	}

	/** Finds a non-sentinel key that maps from cell zero to cell one after 1->2 growth. */
	private static int firstHash32shiftOdd(){
		int key=0;
		while(key<1024 && (Tools.hash32shift(key)&1)==0){key++;}
		check(key<1024, "Fixture must find a hash32shift odd key");
		return key;
	}

	/** Assertions remain meaningful under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
