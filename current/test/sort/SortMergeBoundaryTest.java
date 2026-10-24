package test.sort;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;

import sort.ReadComparator;
import sort.ReadComparatorName;
import sort.SortByName;
import stream.Read;

/** Literal upper-bound cases for merging equal read names without opening files.
 * @author Shinobu
 */
public final class SortMergeBoundaryTest{

	/** Runs expect=baseline or expect=fixed in a fresh test JVM. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || (!args[0].equals("expect=baseline") && !args[0].equals("expect=fixed"))){
			throw new IllegalArgumentException("Expected expect=baseline or expect=fixed");
		}
		final boolean fixed=args[0].equals("expect=fixed");
		boolean assertions=false;
		assert(assertions=true);
		final Field field=SortByName.class.getDeclaredField("comparator");
		field.setAccessible(true);
		final ReadComparator old=(ReadComparator)field.get(null);
		final Method method=SortByName.class.getDeclaredMethod("indexOfLowestAbovePivot", ArrayList.class, Read.class);
		method.setAccessible(true);
		SortByName.setComparator(ReadComparatorName.ascending);
		try{
			check(method, "empty", new String[0], "b", 0, 0, false, fixed, assertions);
			check(method, "before-first", new String[]{"b", "c"}, "a", 0, 0, false, fixed, assertions);
			check(method, "after-last", new String[]{"a", "b"}, "c", 2, 2, false, fixed, assertions);
			check(method, "single-equal", new String[]{"b"}, "b", 1, 1, false, fixed, assertions);
			check(method, "one-equal", new String[]{"a", "b", "c"}, "b", 2, 2, false, fixed, assertions);
			check(method, "duplicate-middle", new String[]{"a", "b", "b", "c"}, "b", 3, 2, true, fixed, assertions);
			check(method, "duplicate-tail", new String[]{"a", "b", "b"}, "b", 3, 2, true, fixed, assertions);
			check(method, "all-equal", new String[]{"b", "b", "b", "b"}, "b", 4, 2, true, fixed, assertions);
		}finally{SortByName.setComparator(old);}
	}

	/** Compares the observed boundary with literal counts of keys at or before the pivot. */
	private static void check(Method method, String label, String[] names, String pivot,
			int expected, int baselineValue, boolean baselineAssertion, boolean fixed, boolean assertions) throws Exception{
		final ArrayList<Read> reads=new ArrayList<Read>();
		for(String name : names){reads.add(read(name));}
		Integer actual=null;
		boolean assertionFailure=false;
		try{actual=(Integer)method.invoke(null, reads, read(pivot));}
		catch(InvocationTargetException e){
			if(e.getCause() instanceof AssertionError){assertionFailure=true;}
			else{throw e;}
		}
		final boolean expectedAssertion=!fixed && baselineAssertion && assertions;
		final int expectedValue=fixed ? expected : baselineValue;
		if(assertionFailure!=expectedAssertion || (!assertionFailure && (actual==null || actual!=expectedValue))){
			throw new AssertionError(label+": boundary="+actual+", assertion="+assertionFailure+
					", expected boundary="+expectedValue+", expected assertion="+expectedAssertion);
		}
		System.out.println("PASS\t"+(fixed ? "fixed" : "baseline")+'\t'+label+'\t'+
				(assertionFailure ? "AssertionError" : actual.toString()));
	}

	/** Read names alone determine these fixtures' order; bases and pair numbers are identical. */
	private static Read read(String name){return new Read(new byte[]{'A'}, null, name, 0);}
}
