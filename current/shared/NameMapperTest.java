package shared;

/**
 * Regression tests for exact and whitespace-truncated sequence-name matching.
 */
public class NameMapperTest {

	public static void main(String[] args){
		testShortToFull();
		testFullToShort();
		testStoredAmbiguity();
		testQueryAmbiguity();
		testExactShortNamePrecedence();
		testReplacement();
		System.out.println("NameMapper tests passed.");
	}

	private static void testShortToFull(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("chr1", "value");
		check("value".equals(map.get("chr1 descriptive text")), "short-to-full lookup failed");
	}

	private static void testFullToShort(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("chr1 descriptive text", "value");
		check("value".equals(map.get("chr1")), "full-to-short lookup failed");
	}

	private static void testStoredAmbiguity(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("dup first", "first");
		map.put("dup second", "second");
		check("first".equals(map.get("dup first")), "exact first lookup failed");
		check("second".equals(map.get("dup second")), "exact second lookup failed");
		expectAmbiguous(map, "dup");
	}

	private static void testQueryAmbiguity(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("dup", "value");
		check("value".equals(map.get("dup first")), "first fallback query failed");
		expectAmbiguous(map, "dup second");
	}

	private static void testExactShortNamePrecedence(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("dup first", "first");
		map.put("dup", "exact");
		check("exact".equals(map.get("dup")), "exact short name did not take precedence");
		check("first".equals(map.get("dup first")), "full exact name was lost");
		expectAmbiguous(map, "dup second");
	}

	private static void testReplacement(){
		NameMapper<String> map=new NameMapper<String>();
		map.put("chr1 descriptive text", "old");
		String old=map.put("chr1 descriptive text", "new");
		check("old".equals(old), "replacement did not return old value");
		check("new".equals(map.get("chr1")), "replacement did not update alias");
	}

	private static void expectAmbiguous(NameMapper<String> map, String name){
		try{
			map.get(name);
			throw new AssertionError("Expected ambiguous lookup failure for "+name);
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains("Ambiguous"), "wrong ambiguity exception");
		}
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
