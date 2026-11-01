package test.shared;

import shared.NameMapper;

/**
 * Regression tests for exact and whitespace-truncated sequence-name matching.
 */
public class NameMapperTest{

	/** Runs six finite in-memory name-resolution scenarios; arguments are unused. */
	public static void main(String[] args){
		testShortToFull();
		testFullToShort();
		testStoredAmbiguity();
		testQueryAmbiguity();
		testExactShortNamePrecedence();
		testReplacement();
		System.out.println("NameMapper tests passed.");
	}

	/** A stored short name accepts one expanded query. */
	private static void testShortToFull(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("chr1", "value");
		check("value".equals(map.get("chr1 descriptive text")), "short-to-full lookup failed");
	}

	/** A stored descriptive name supplies its unique short alias. */
	private static void testFullToShort(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("chr1 descriptive text", "value");
		check("value".equals(map.get("chr1")), "full-to-short lookup failed");
	}

	/** Colliding stored aliases preserve exact names and reject the short lookup. */
	private static void testStoredAmbiguity(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("dup first", "first");
		map.put("dup second", "second");
		check("first".equals(map.get("dup first")), "exact first lookup failed");
		check("second".equals(map.get("dup second")), "exact second lookup failed");
		expectAmbiguous(map, "dup");
	}

	/** Distinct expanded queries may not silently share one fallback token. */
	private static void testQueryAmbiguity(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("dup", "value");
		check("value".equals(map.get("dup first")), "first fallback query failed");
		expectAmbiguous(map, "dup second");
	}

	/** An explicitly stored short name wins without erasing the full exact name. */
	private static void testExactShortNamePrecedence(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("dup first", "first");
		map.put("dup", "exact");
		check("exact".equals(map.get("dup")), "exact short name did not take precedence");
		check("first".equals(map.get("dup first")), "full exact name was lost");
		expectAmbiguous(map, "dup second");
	}

	/** Replacing an exact name updates its otherwise unambiguous alias. */
	private static void testReplacement(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("chr1 descriptive text", "old");
		String old=map.put("chr1 descriptive text", "new");
		check("old".equals(old), "replacement did not return old value");
		check("new".equals(map.get("chr1")), "replacement did not update alias");
	}

	/** Requires the documented ambiguity exception rather than an arbitrary failure. */
	private static void expectAmbiguous(NameMapper<String> map, String name){
		try{
			map.get(name);
			throw new AssertionError("Expected ambiguous lookup failure for "+name);
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains("Ambiguous"), "wrong ambiguity exception");
		}
	}

	/** Checks the oracle independently of JVM assertion settings. */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
