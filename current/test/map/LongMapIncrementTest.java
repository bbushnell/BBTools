package test.map;

import map.LongIntMapX;
import map.LongLongHashMap2;
import shared.Tools;

/** Regression fixture for LongIntMapX and LongLongHashMap2 increment return values. */
public final class LongMapIncrementTest{

	public static void main(String[] args){
		check(args.length==1, "Usage: expect=baseline|fixed");
		String[] split=args[0].split("=", 2);
		check(split.length==2 && split[0].equals("expect"), "Expected expect=baseline|fixed");
		boolean fixed;
		if(split[1].equals("baseline")){fixed=false;}
		else if(split[1].equals("fixed")){fixed=true;}
		else{throw new AssertionError("Expected expect=baseline|fixed, got "+args[0]);}

		checkLongIntMapXRelocation(fixed);
		checkLongLongHashMap2Relocation(fixed);
		checkLongLongHashMap2Overflow(fixed);
		System.err.println("PASS LongMapIncrementTest expect="+split[1]);
	}

	private static void checkLongIntMapXRelocation(boolean fixed){
		long key=oddHash64plus2Key();
		LongIntMapX map=new LongIntMapX(1);
		int oldLength=map.keys().length;
		int returned=map.increment(key, 7);
		check(map.keys().length>oldLength, "LongIntMapX should grow on first insert");
		check(map.keys()[1]==key, "LongIntMapX key should relocate from cell 0 to cell 1");
		check(map.values()[1]==7, "LongIntMapX value should remain stored after resize");
		check(map.get(key)==7 && map.size()==1, "LongIntMapX stored value/size mismatch");
		check(returned==(fixed ? 7 : 0), "LongIntMapX returned "+returned);
	}

	private static void checkLongLongHashMap2Relocation(boolean fixed){
		long key=oddHash64plusKey();
		LongLongHashMap2 map=new LongLongHashMap2(1);
		int oldLength=map.keys().length;
		long returned=map.increment(key, 7);
		check(map.keys().length>oldLength, "LongLongHashMap2 should grow on first insert");
		check(map.keys()[1]==key, "LongLongHashMap2 key should relocate from cell 0 to cell 1");
		check(map.values()[1]==7, "LongLongHashMap2 value should remain stored after resize");
		check(map.get(key)==7 && map.size()==1, "LongLongHashMap2 stored value/size mismatch");
		check(returned==(fixed ? 7 : 0), "LongLongHashMap2 returned "+returned);
	}

	private static void checkLongLongHashMap2Overflow(boolean fixed){
		LongLongHashMap2 map=new LongLongHashMap2(4);
		map.put(0, Long.MAX_VALUE-1);
		long returned=map.increment(0, 7);
		long expected=(fixed ? Long.MAX_VALUE : Long.MIN_VALUE+5);
		check(returned==expected, "LongLongHashMap2 overflow returned "+returned);
		check(map.get(0)==expected, "LongLongHashMap2 overflow stored "+map.get(0));
	}

	private static long oddHash64plusKey(){
		for(long key=0; key<1024; key++){
			if((Tools.hash64plus(key)&1L)==1){return key;}
		}
		throw new AssertionError("No odd Tools.hash64plus key found");
	}

	private static long oddHash64plus2Key(){
		for(long key=0; key<1024; key++){
			if(((Tools.hash64plus2(key)>>>32)&1L)==1){return key;}
		}
		throw new AssertionError("No odd Tools.hash64plus2 key found");
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	private LongMapIncrementTest(){
		throw new AssertionError();
	}
}
