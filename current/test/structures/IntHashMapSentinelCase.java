package test.structures;

import structures.IntHashMap;
import structures.IntHashMapBinary;

/** Baseline and fixed oracles for STU-009 IntHashMap sentinel lookups.
 * @author Jean
 */
public final class IntHashMapSentinelCase{
	
	/** Runs the IntHashMap sentinel lookup checks. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		
		intHashMapInvalidLookup(expect);
		intHashMapBinaryInvalidLookup(expect);
		intHashMapFormerInvalidKey(expect);
		intHashMapBinaryFormerInvalidKey(expect);
	}
	
	/** Verifies STU-009: IntHashMap must not find its empty-cell sentinel. */
	private static void intHashMapInvalidLookup(final String expect){
		final IntHashMap map=new IntHashMap(3);
		final int invalid=map.invalid();
		final boolean baseline="baseline".equals(expect);
		checkEquals("int-map-contains-invalid", baseline, map.contains(invalid));
		checkEquals("int-map-contains-key-invalid", baseline, map.containsKey(invalid));
		checkEquals("int-map-get-invalid", baseline ? 0 : -1, map.get(invalid));
		pass("int-map-invalid-lookup", map.contains(invalid)+","+map.get(invalid), expect);
		
		final IntHashMap removeControl=new IntHashMap(3);
		checkEquals("int-map-remove-invalid-control", false,
				removeControl.remove(removeControl.invalid()));
		pass("int-map-remove-invalid-control", "false", expect);
	}
	
	/** Verifies the IntHashMapBinary lookup control. */
	private static void intHashMapBinaryInvalidLookup(final String expect){
		final IntHashMapBinary map=new IntHashMapBinary(4);
		final int invalid=map.invalid();
		checkEquals("binary-map-contains-invalid-control", false, map.contains(invalid));
		checkEquals("binary-map-contains-key-invalid-control", false,
				map.containsKey(invalid));
		checkEquals("binary-map-get-invalid-control", -1, map.get(invalid));
		
		final IntHashMapBinary removeControl=new IntHashMapBinary(4);
		checkEquals("binary-map-remove-invalid-control", false,
				removeControl.remove(removeControl.invalid()));
		pass("binary-map-invalid-controls", "false,-1,false", expect);
	}
	
	/** Verifies an old IntHashMap sentinel remains usable after resetInvalid(). */
	private static void intHashMapFormerInvalidKey(final String expect){
		final IntHashMap map=new IntHashMap(3);
		final int formerInvalid=map.invalid();
		final int old=map.set(formerInvalid, 7);
		checkEquals("int-map-former-invalid-old", 0, old);
		check(map.invalid()!=formerInvalid,
				"IntHashMap should choose a new invalid sentinel before storing the old one");
		checkEquals("int-map-former-invalid-contains", true, map.contains(formerInvalid));
		checkEquals("int-map-former-invalid-get", 7, map.get(formerInvalid));
		checkEquals("int-map-former-invalid-remove", true, map.remove(formerInvalid));
		checkEquals("int-map-former-invalid-removed", false, map.contains(formerInvalid));
		pass("int-map-former-invalid-key", Integer.toString(map.get(formerInvalid)), expect);
	}
	
	/** Verifies the IntHashMapBinary old-sentinel insertion control. */
	private static void intHashMapBinaryFormerInvalidKey(final String expect){
		final IntHashMapBinary map=new IntHashMapBinary(4);
		final int formerInvalid=map.invalid();
		final int old=map.set(formerInvalid, 7);
		checkEquals("binary-map-former-invalid-old", 0, old);
		check(map.invalid()!=formerInvalid,
				"IntHashMapBinary should choose a new invalid sentinel before storing the old one");
		checkEquals("binary-map-former-invalid-contains", true, map.contains(formerInvalid));
		checkEquals("binary-map-former-invalid-get", 7, map.get(formerInvalid));
		checkEquals("binary-map-former-invalid-remove", true, map.remove(formerInvalid));
		checkEquals("binary-map-former-invalid-removed", false, map.contains(formerInvalid));
		pass("binary-map-former-invalid-key", Integer.toString(map.get(formerInvalid)), expect);
	}
	
	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}
	
	/** Checks a boolean oracle. */
	private static void checkEquals(final String name, final boolean expected,
			final boolean observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}
	
	/** Prints one successful case row. */
	private static void pass(final String name, final String observed, final String expect){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}
	
	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
