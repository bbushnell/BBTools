package test.structures;

import structures.ByteBuilder;

/** Literal substring and zero-decimal rounding checks for ByteBuilder.
 * @author Shinobu
 */
public final class ByteBuilderContractCase{

	/** Runs baseline or repaired expectations for valid empty spans and negative fractions. */
	public static void main(final String[] args){
		check(args.length==1 && ("expect=baseline".equals(args[0]) || "expect=fixed".equals(args[0])),
				"Expected expect=baseline or expect=fixed");
		final boolean baseline="expect=baseline".equals(args[0]);
		checkOutput("empty-whole-control", "p", new ByteBuilder("p").append(""));
		checkOutput("empty-input-slice", baseline ? "pnull" : "p",
				new ByteBuilder("p").append("", 0, 0));
		checkOutput("empty-middle-slice", baseline ? "pnull" : "p",
				new ByteBuilder("p").append("abc", 1, 1));
		checkOutput("ordinary-slice-control", "pbc", new ByteBuilder("p").append("abc", 1, 3));
		checkOutput("null-slice-control", "pnull", new ByteBuilder("p").append((String)null, 0, 0));
		checkOutput("negative-round-near", baseline ? "0" : "-1", new ByteBuilder().append(-1.2, 0));
		checkOutput("negative-round-far", baseline ? "-1" : "-2", new ByteBuilder().append(-1.6, 0));
		checkOutput("positive-round-control", "2", new ByteBuilder().append(1.6, 0));
		checkOutput("integer-control", "-1", new ByteBuilder().append(-1.0, 0));
		checkOutput("nonzero-decimals-control", "-1.2", new ByteBuilder().append(-1.2, 1));
	}

	/** Requires a literal string and emits the observed bytes as text. */
	private static void checkOutput(final String name, final String expected, final ByteBuilder builder){
		final String observed=builder.toString();
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		System.out.println("PASS\t"+name+'\t'+observed);
	}

	/** Fails independently of JVM assertions when an expected result differs. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
