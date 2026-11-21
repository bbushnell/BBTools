package test.processor;

import java.lang.reflect.Field;
import java.util.HashSet;

import map.ObjectIntMap;
import processor.UniqueReadNameAllocator;

/**
 * Independent literal name, Unicode, storage-bound and exhaustion controls.
 * Reflection installs only a terminal counter history that would otherwise need
 * billions of records; assertions exercise the public allocation result and state.
 * No production allocator method supplies expected names.
 * @author Sucrose
 */
public final class UniqueReadNameAllocatorTest{

	public static void main(final String[] args) throws Exception{
		literal(false, new String[]{"x", "x", "x_2147483648", "x_2147483648",
			"x_9223372036854775808", "x_9223372036854775808",
			"x_999999999999999999999999999999", "x_999999999999999999999999999999", "x_02", "x_02"},
			new String[]{"x", "x_2", "x_2147483648", "x_2147483648_2",
			"x_9223372036854775808", "x_9223372036854775808_2",
			"x_999999999999999999999999999999", "x_999999999999999999999999999999_2", "x_02", "x_02_2"});
		literal(false, new String[]{"Aa", "BB", "Aa", "BB"}, new String[]{"Aa", "BB", "Aa_2", "BB_2"});
		//Shinobu's source review: stripping a numeric suffix can reconstruct a parent with a different insertion boundary.
		literal(false, new String[]{"x/1", "x/1", "x/1_2"}, new String[]{"x/1", "x_2/1", "x/1_2"});
		literal(true, new String[]{"x  /1", "x  /1", "x _2 /1"},
			new String[]{"x  /1", "x_2  /1", "x _2 /1"});
		literal(true, new String[]{"p 1:", "p 1:", "p 1:_2"}, new String[]{"p 1:", "p_2 1:", "p 1:_2"});
		literal(false, new String[]{"Aa/1", "BB/1", "Aa/1", "BB/1", "Aa_2/1", "BB_2/1", "Aa/1", "BB/1"},
			new String[]{"Aa/1", "BB/1", "Aa_2/1", "BB_2/1", "Aa_2_2/1", "BB_2_2/1", "Aa_3/1", "BB_3/1"});
		literal(false, new String[]{"é description", "é\tother", "é_2 last", "é final"},
			new String[]{"é description", "é_2\tother", "é_2_2 last", "é_3 final"});
		literal(true, new String[]{"é desc /1", "é desc /1", "é desc_2 /1"},
			new String[]{"é desc /1", "é desc_2 /1", "é desc_2_2 /1"});
		literal(false, new String[]{"😀/1 note", "😀/1 other", "😀_2/1 last"},
			new String[]{"😀/1 note", "😀_2/1 other", "😀_2_2/1 last"});
		literal(true, new String[]{"p desc 1:", "p desc 1:", "p desc_2 1:"},
			new String[]{"p desc 1:", "p desc_2 1:", "p desc_2_2 1:"});
		literal(false, new String[]{"", "", "_2", " a", "_02"},
			new String[]{"", "_2", "_2_2", "_3 a", "_02"});
		literal(true, new String[]{"p /1", "p /1", "p_2 /1", "p /2", "p /2"},
			new String[]{"p /1", "p_2 /1", "p_2_2 /1", "p /2", "p_2 /2"});
		literal(false, new String[]{"p_1", "p_0", "p_02", "p", "p", "p_2"},
			new String[]{"p_1", "p_0", "p_02", "p", "p_2", "p_2_2"});
		final String fresh=new String("fresh description");
		check(new UniqueReadNameAllocator(false).allocate(fresh)==fresh, "Unused IDs must retain the original String.");
		check(!UniqueReadNameAllocator.parseScope("identifier") && UniqueReadNameAllocator.parseScope("header"),
			"Scope selection must distinguish identifier and header.");
		check(!UniqueReadNameAllocator.parseScope("IDENTIFIER") && UniqueReadNameAllocator.parseScope("HEADER"),
			"Scope selector spelling is case-insensitive.");
		for(final String value : new String[]{null, "", "other"}){
			boolean rejected=false;
			try{UniqueReadNameAllocator.parseScope(value);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Unknown scope must fail clearly: "+value);
		}
		storageBound();
		exhaustion();
		System.out.println("UNIQUE_READ_NAME_ALLOCATOR_PASS checks="+checks);
	}

	/** Literal oracles are unrelated to the allocator's counter/parent lookup mechanism. */
	private static void literal(final boolean header, final String[] input, final String[] expected){
		check(input.length==expected.length, "Literal case requires one expected name per input.");
		final UniqueReadNameAllocator allocator=new UniqueReadNameAllocator(header);
		for(int i=0; i<input.length; i++){
			final String actual=allocator.allocate(input[i]);
			check(expected[i].equals(actual), "Expected ["+expected[i]+"], observed ["+actual+"].");
		}
	}

	/** Ten thousand generated outputs retain one history; a new natural input adds exactly one. */
	private static void storageBound() throws Exception{
		final UniqueReadNameAllocator allocator=new UniqueReadNameAllocator(false);
		final HashSet<String> output=new HashSet<String>();
		String last=null;
		for(int i=0; i<10000; i++){
			last=allocator.allocate("bounded");
			check(output.add(last), "Generated output repeated: "+last);
		}
		check("bounded_10000".equals(last), "The serial ordinal must survive repeated lookups and updates.");
		check(histories(allocator).size()==1, "Generated IDs must not create per-output history entries.");
		check("bounded_5000_2".equals(allocator.allocate("bounded_5000")), "A prior generated name remains occupied.");
		check(histories(allocator).size()==2, "Only the newly encountered natural input adds a second history.");
	}

	/** Tests the last int ordinal, an occupied terminal ordinal, and failure without wrapped history. */
	private static void exhaustion() throws Exception{
		final UniqueReadNameAllocator allocator=new UniqueReadNameAllocator(false);
		check("edge".equals(allocator.allocate("edge")), "Seed key must be initially unused.");
		final ObjectIntMap<?> map=histories(allocator);
		final int cell=onlyCell(map);
		map.values()[cell]=Integer.MAX_VALUE-1;
		check("edge_2147483647".equals(allocator.allocate("edge")), "The last representable suffix must be usable.");
		check(map.values()[cell]==Integer.MAX_VALUE && map.size()==1, "The last output updates only its original history.");
		boolean rejected=false;
		try{allocator.allocate("edge");}
		catch(IllegalStateException expected){rejected=true;}
		check(rejected && map.values()[cell]==Integer.MAX_VALUE, "Exhaustion must fail without wrapping or mutating the counter.");
		check("edge_2147483648".equals(allocator.allocate("edge_2147483648")), "An oversized natural suffix is still valid after exhaustion.");

		final UniqueReadNameAllocator reserved=new UniqueReadNameAllocator(false);
		reserved.allocate("last");
		final ObjectIntMap<?> other=histories(reserved);
		final int otherCell=onlyCell(other);
		other.values()[otherCell]=Integer.MAX_VALUE-1;
		check("last_2147483647".equals(reserved.allocate("last_2147483647")), "A not-yet-generated terminal suffix can be a natural input.");
		rejected=false;
		try{reserved.allocate("last");}
		catch(IllegalStateException expected){rejected=true;}
		check(rejected && other.values()[otherCell]==Integer.MAX_VALUE-1 && other.size()==2,
			"Skipping an occupied terminal suffix must fail before altering either input history.");
	}

	/** Reads the private history only to verify its storage bound or seed a numeric boundary. */
	private static ObjectIntMap<?> histories(final UniqueReadNameAllocator allocator) throws Exception{
		final Field field=UniqueReadNameAllocator.class.getDeclaredField("histories");
		field.setAccessible(true);
		return (ObjectIntMap<?>)field.get(allocator);
	}

	/** Finds the sole occupied key slot without depending on its hash/table position. */
	private static int onlyCell(final ObjectIntMap<?> map){
		check(map.size()==1, "Counter seeding requires exactly one input history.");
		final Object[] keys=map.keys();
		for(int i=0; i<keys.length; i++){if(keys[i]!=null){return i;}}
		throw new AssertionError("The single history must have a stored key.");
	}

	/** Literal checks remain enabled regardless of the JVM assertion setting. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
}
