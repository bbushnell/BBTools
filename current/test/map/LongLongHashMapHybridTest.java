package test.map;

import map.LongLongHashMapHybrid;
import structures.LongList;

/**
 * Finite regression checks for LongLongHashMapHybrid side-list reuse.
 * @author Jean
 * @date October 25, 2026
 */
public class LongLongHashMapHybridTest{

	public static void main(String[] args){
		final boolean fixed=parseExpectation(args);
		LongLongHashMapHybrid map=new LongLongHashMapHybrid(4, 0.7f);
		LongList buffer=new LongList();

		check(map.put(7, 70), "first value should insert key 7");
		check(!map.put(7, 70), "duplicate inline value should be suppressed");
		check(map.put(7, 77), "second distinct value should promote key 7");
		check(map.put(7, 77), "duplicate side-list value should be appended");
		long firstCode=map.get(7);
		check(firstCode==-2, "first multi-list code should be -2, was "+firstCode);
		map.fill(7, buffer);
		check(buffer.size==3, "side list should contain three values, found "+buffer.size);
		check(buffer.array[0]==70 && buffer.array[1]==77 && buffer.array[2]==77,
				"side list should preserve promoted and duplicate values");

		map.clear();
		check(map.isEmpty(), "clear should remove all keys");
		check(map.get(7)==LongLongHashMapHybrid.NOTPRESENT, "cleared key should be absent");

		check(map.put(9, 90), "first value should insert key 9");
		check(map.put(9, 99), "second value should promote key 9");
		long secondCode=map.get(9);
		checkCode(secondCode, fixed ? -2 : -3, "clear with one live key should reset side-list codes");

		check(map.put(11, 110), "first value should insert key 11");
		check(map.put(11, 111), "second value should promote key 11");
		check(map.put(12, 120), "first value should insert key 12");
		check(map.put(12, 121), "second value should promote key 12");
		check(map.remove(11)<0, "removing key 11 should return a multi-list code");
		check(map.remove(12)<0, "removing key 12 should return a multi-list code");
		check(map.remove(9)<0, "removing key 9 should return a multi-list code");
		check(map.isEmpty(), "removing all multi-list keys should leave the map empty");
		map.clear();

		check(map.put(13, 130), "first value should insert key 13");
		check(map.put(13, 131), "second value should promote key 13");
		long thirdCode=map.get(13);
		checkCode(thirdCode, fixed ? -2 : -5, "clear after remove-created holes should reset side-list codes");

		System.err.println("PASS");
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

	private static void checkCode(long actual, long expected, String message){
		if(actual!=expected){throw new AssertionError(message+": "+actual+" != "+expected);}
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	private LongLongHashMapHybridTest(){}
}
