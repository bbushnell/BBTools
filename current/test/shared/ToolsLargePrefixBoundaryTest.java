package test.shared;

import shared.Tools;

/** Tests valid-offset prefix and contains guards with one shared large zeroed array.
 * The same array is the query in the insufficient-suffix case, avoiding a second
 * large allocation. Only the last byte or the immediately following invalid
 * position can be visited by the current baseline; the fixture never scans the
 * large array itself. Run large lengths only in the designated cluster job.
 * @author Sucrose
 */
public final class ToolsLargePrefixBoundaryTest{

	/** Accepts optional length=N (default16); the overflow fixture uses1073741825. */
	public static void main(String[] args){
		if(args.length>1 || (args.length==1 && !args[0].startsWith("length="))){
			throw new IllegalArgumentException("Expected zero arguments or length=N");
		}
		final int length=args.length==0 ? 16 : Integer.parseInt(args[0].substring(7));
		if(length<1){throw new IllegalArgumentException("length must be positive: "+length);}
		final byte[] array=new byte[length];
		final int index=length-1;
		final ToolsLargePrefixBoundaryTest test=new ToolsLargePrefixBoundaryTest();
		test.check("last zero byte", () -> {
			require(Tools.startsWith(array, (byte)0, index), "Valid last-byte zero match returned false");
		});
		test.check("same-array query at last index", () -> {
			final boolean actual=Tools.startsWith(array, array, index);
			require(actual==(length==1), "Expected "+(length==1)+" with one byte remaining, actual="+actual);
		});
		test.check("same-array contains at last index", () -> {
			final boolean actual=Tools.contains(array, array, index);
			require(actual==(length==1), "Expected "+(length==1)+" with one byte remaining, actual="+actual);
		});
		test.check("tiny separate controls", () -> {
			final byte[] tiny={0, 1, 2}, suffix={1, 2};
			require(Tools.startsWith(tiny, (byte)0, 0), "Tiny first-byte match");
			require(!Tools.startsWith(tiny, (byte)0, 1), "Tiny byte mismatch");
			require(Tools.startsWith(tiny, suffix, 1), "Tiny fitting suffix match");
			require(!Tools.startsWith(tiny, suffix, 0), "Tiny prefix mismatch");
			require(!Tools.startsWith(tiny, suffix, 2), "Tiny insufficient suffix");
		});
		test.check("tiny contains controls", () -> {
			final byte[] tiny={0, 1, 2}, suffix={1, 2};
			require(Tools.contains(tiny, suffix, 1), "Tiny contains match at the requested offset");
			require(!Tools.contains(tiny, suffix, 0), "Contains must not search later offsets after a mismatch");
			require(!Tools.contains(tiny, suffix, 2), "Tiny contains insufficient suffix");
		});
		System.out.println("ToolsLargePrefixBoundaryTest length="+length+" index="+index+
				" cases="+test.cases+" failures="+test.failures);
		if(test.failures>0){throw new AssertionError("Valid-offset prefix/contains failures="+test.failures);}
	}

	/** Reports a thrown exception as a failed case, then continues the other bounded checks. */
	private void check(String name, Runnable action){
		cases++;
		try{action.run();System.out.println(name+" PASS");}
		catch(Exception | AssertionError failure){
			failures++;
			System.out.println(name+" FAIL: "+failure);
		}
	}

	/** Keeps literal expectations active with JVM assertions disabled too. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Attempted and failed groups. */
	private int cases, failures;
}
