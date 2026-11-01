package test.parse;

import parse.Parse;

/** Literal signed and fractional boundary oracles for command-line number parsing.
 * This fixture covers valid inputs, not a new overflow or malformed-input policy.
 * @author Brian Bushnell, Shinobu */
public final class NumericBoundaryTest {

	/** Executes each case independently in either assertion mode.
	 * @param args No arguments are supported */
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No arguments expected");}
		final NumericBoundaryTest test=new NumericBoundaryTest();
		test.check("int minimum", () -> Parse.parseIntKMG("-2147483648"), Integer.MIN_VALUE);
		test.check("int maximum", () -> Parse.parseIntKMG("2147483647"), Integer.MAX_VALUE);
		test.check("byte int minimum", () -> Parse.parseInt("-2147483648".getBytes(java.nio.charset.StandardCharsets.US_ASCII), 0, 11), Integer.MIN_VALUE);
		test.check("string int minimum", () -> Parse.parseInt("-2147483648", 0, 11), Integer.MIN_VALUE);
		test.check("negative decimal scale", () -> Parse.parseKMG("-1k"), -1000);
		test.check("negative fractional scale", () -> Parse.parseKMG("-1.5k"), -1500);
		test.check("positive fractional count", () -> Parse.parseKMG("1.5"), 1);
		test.check("positive subunit count", () -> Parse.parseKMG("0.5"), 0);
		test.check("negative fractional count", () -> Parse.parseKMG("-1.5"), -1);
		test.check("negative subunit count", () -> Parse.parseKMG("-0.5"), 0);
		test.check("positive scaled subunit", () -> Parse.parseKMG("0.0005k"), 0);
		test.check("negative scaled subunit", () -> Parse.parseKMG("-0.0005k"), 0);
		test.check("positive decimal scale", () -> Parse.parseKMG("2k"), 2000);
		test.check("positive fractional scale", () -> Parse.parseKMG("1.5k"), 1500);
		test.check("binary negative scale", () -> Parse.parseKMGBinary("-2k"), -2048);
		test.check("binary fractional scale", () -> Parse.parseKMGBinary("1.5k"), 1536);
		System.out.println("NUMERIC_BOUNDARY_RESULT cases="+test.cases+" failures="+test.failures);
		if(test.failures!=0){throw new AssertionError("Numeric boundary failures="+test.failures);}
	}

	/** Compares a literal oracle while retaining all failures in the result. */
	private void check(String name, java.util.function.LongSupplier action, long expected){
		cases++;
		try{
			final long actual=action.getAsLong();
			if(actual!=expected){throw new AssertionError(actual+" != "+expected);}
		}catch(RuntimeException | AssertionError e){
			failures++;
			System.out.println("FAIL "+name+": "+e);
		}
	}

	/** Number of attempted checks. */
	private int cases;
	/** Number of failed checks. */
	private int failures;
}
