package test.structures;

import structures.AtomicStringNum;
import structures.LongM;
import structures.StringNum;

/** Literal equals-contract checks for mutable structures value wrappers.
 * @author Jean
 */
public final class MutableEqualsContractCase{
	
	/** Runs current and proposed equals-contract expectations. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation: "+args[0]);
		final String expect=(fixed ? "fixed" : "baseline");
		
		longMControls(expect);
		longMBoundaries(expect, fixed);
		stringNumControls(expect);
		stringNumBoundaries(expect, fixed);
		atomicStringNumControls(expect);
		atomicStringNumBoundaries(expect, fixed);
	}
	
	/** Checks literal LongM equality controls. */
	private static void longMControls(final String expect){
		final LongM value=new LongM(7);
		checkBoolean(expect, "longm-equal-typed", true, value.equals(new LongM(7)));
		checkBoolean(expect, "longm-equal-object", true, value.equals((Object)new LongM(7)));
		checkBoolean(expect, "longm-different-value", false, value.equals(new LongM(8)));
	}
	
	/** Checks current LongM throws versus proposed Object.equals-compatible false. */
	private static void longMBoundaries(final String expect, final boolean fixed){
		checkObservation(expect, "longm-object-null", fixed ? "false" : "NullPointerException",
				observeLongMObject(null));
		checkObservation(expect, "longm-object-unrelated", fixed ? "false" : "ClassCastException",
				observeLongMObject("not a LongM"));
		checkObservation(expect, "longm-typed-null", fixed ? "false" : "NullPointerException",
				observeLongMTypedNull());
	}
	
	/** Checks literal StringNum equality controls. */
	private static void stringNumControls(final String expect){
		final StringNum value=new StringNum("label", 7);
		checkBoolean(expect, "stringnum-equal-typed", true, value.equals(new StringNum("label", 7)));
		checkBoolean(expect, "stringnum-equal-object", true,
				value.equals((Object)new StringNum("label", 7)));
		checkBoolean(expect, "stringnum-different-value", false,
				value.equals(new StringNum("label", 8)));
		checkBoolean(expect, "stringnum-different-label", false,
				value.equals(new StringNum("other", 7)));
	}
	
	/** Checks current StringNum throws versus proposed Object.equals-compatible false. */
	private static void stringNumBoundaries(final String expect, final boolean fixed){
		checkObservation(expect, "stringnum-object-null", "false", observeStringNumObject(null));
		checkObservation(expect, "stringnum-object-unrelated", fixed ? "false" : "ClassCastException",
				observeStringNumObject("not a StringNum"));
		checkObservation(expect, "stringnum-typed-null", "false", observeStringNumTypedNull());
	}
	
	/** Checks literal AtomicStringNum equality controls. */
	private static void atomicStringNumControls(final String expect){
		final AtomicStringNum value=new AtomicStringNum("label", 7);
		checkBoolean(expect, "atomicstringnum-equal-typed", true,
				value.equals(new AtomicStringNum("label", 7)));
		checkBoolean(expect, "atomicstringnum-equal-object", true,
				value.equals((Object)new AtomicStringNum("label", 7)));
		checkBoolean(expect, "atomicstringnum-different-value", false,
				value.equals(new AtomicStringNum("label", 8)));
		checkBoolean(expect, "atomicstringnum-different-label", false,
				value.equals(new AtomicStringNum("other", 7)));
	}
	
	/** Checks current AtomicStringNum throws versus proposed Object.equals-compatible false. */
	private static void atomicStringNumBoundaries(final String expect, final boolean fixed){
		checkObservation(expect, "atomicstringnum-object-null", "false",
				observeAtomicStringNumObject(null));
		checkObservation(expect, "atomicstringnum-object-unrelated",
				fixed ? "false" : "ClassCastException",
				observeAtomicStringNumObject("not an AtomicStringNum"));
		checkObservation(expect, "atomicstringnum-typed-null", "false",
				observeAtomicStringNumTypedNull());
	}
	
	/** Calls LongM.equals(Object) and returns a value or thrown type. */
	private static String observeLongMObject(final Object other){
		try{
			return Boolean.toString(new LongM(7).equals(other));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}
	
	/** Calls LongM.equals(LongM) with a typed null and returns a value or thrown type. */
	private static String observeLongMTypedNull(){
		try{
			return Boolean.toString(new LongM(7).equals((LongM)null));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}
	
	/** Calls StringNum.equals(Object) and returns a value or thrown type. */
	private static String observeStringNumObject(final Object other){
		try{
			return Boolean.toString(new StringNum("label", 7).equals(other));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}
	
	/** Calls StringNum.equals(StringNum) with a typed null and returns a value or thrown type. */
	private static String observeStringNumTypedNull(){
		try{
			return Boolean.toString(new StringNum("label", 7).equals((StringNum)null));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}
	
	/** Calls AtomicStringNum.equals(Object) and returns a value or thrown type. */
	private static String observeAtomicStringNumObject(final Object other){
		try{
			return Boolean.toString(new AtomicStringNum("label", 7).equals(other));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}
	
	/** Calls AtomicStringNum.equals(AtomicStringNum) with a typed null and returns a value or thrown type. */
	private static String observeAtomicStringNumTypedNull(){
		try{
			return Boolean.toString(new AtomicStringNum("label", 7).equals((AtomicStringNum)null));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}
	
	/** Checks and prints a boolean outcome. */
	private static void checkBoolean(final String expect, final String name,
			final boolean expected, final boolean observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
		pass(expect, name, Boolean.toString(observed));
	}
	
	/** Checks and prints a value-or-exception outcome. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		pass(expect, name, observed);
	}
	
	/** Prints a checked outcome. */
	private static void pass(final String expect, final String name, final String observed){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}
	
	/** Fails independently of assertion flags so an oracle cannot silently disappear. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
