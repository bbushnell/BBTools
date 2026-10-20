package test.clump;

import java.lang.reflect.Field;

import clump.Clump;

/** Finite parser-state oracle for Clump maxqi handling.
 * @author Sucrose
 */
public final class ClumpMaxQITest{

	/** Runs baseline or fixed expectations for maxqi after maxci controls. */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);

		System.out.println("case\tmaxci\tmaxqi\tmaxCIncorrect\tmaxQIncorrect");
		for(int maxci : new int[] {-1, 0, 3}){
			checkCase(expect, fixed, maxci, -1);
			checkCase(expect, fixed, maxci, 0);
			checkCase(expect, fixed, maxci, 7);
		}
		reset();
	}

	/** Sets maxci through the parser, then verifies maxqi output state. */
	private static void checkCase(final String expect, final boolean fixed,
			final int maxci, final int maxqi) throws Exception{
		reset();
		check(parse("maxci", maxci), "maxci must be accepted");
		final int expectedC=(maxci<0 ? Integer.MAX_VALUE : maxci);
		check(value("maxCIncorrect")==expectedC,
				"maxci="+maxci+" expected maxCIncorrect "+expectedC+
				", got "+value("maxCIncorrect"));

		check(parse("maxqi", maxqi), "maxqi must be accepted");
		final int expectedQ=(maxqi<0 && fixed ? Integer.MAX_VALUE : maxqi);
		final int observedQ=value("maxQIncorrect");
		check(observedQ==expectedQ,
				"maxci="+maxci+", maxqi="+maxqi+" expected maxQIncorrect "+
				expectedQ+", got "+observedQ);
		System.out.println(expect+'\t'+maxci+'\t'+maxqi+'\t'+expectedC+'\t'+observedQ);
	}

	/** Calls Clump.parseStatic with a conventional arg string. */
	private static boolean parse(final String name, final int value){
		final String string=Integer.toString(value);
		return Clump.parseStatic(name+'='+string, name, string);
	}

	/** Reads an integer parser field from Clump. */
	private static int value(final String name) throws Exception{
		final Field field=Clump.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getInt(null);
	}

	/** Restores the static fields touched by these parser cases. */
	private static void reset() throws Exception{
		set("maxCIncorrect", Integer.MAX_VALUE);
		set("maxQIncorrect", Integer.MAX_VALUE);
	}

	/** Writes an integer parser field in Clump. */
	private static void set(final String name, final int value) throws Exception{
		final Field field=Clump.class.getDeclaredField(name);
		field.setAccessible(true);
		field.setInt(null, value);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private ClumpMaxQITest(){}
}
