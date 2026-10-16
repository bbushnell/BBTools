package test.structures;

import structures.ByteBuilder;

/** Finite growth and scientific-format checks for ByteBuilder.
 * Does not invoke bulk growth on a zero-capacity baseline builder.
 * @author Shinobu
 */
public final class ByteBuilderEdgeCase{

	/** Runs the saved baseline or repaired-output expectations with assertions enabled. */
	public static void main(final String[] args){
		check(ByteBuilder.class.desiredAssertionStatus(),
				"ByteBuilder assertions must be enabled to bound the scientific convergence check");
		check(args.length==1 && ("expect=baseline".equals(args[0]) || "expect=fixed".equals(args[0])),
				"Expected expect=baseline or expect=fixed");
		final boolean baseline="expect=baseline".equals(args[0]);
		final ByteBuilder control=new ByteBuilder(1);
		control.append((byte)'a').append((byte)'b');
		check("ab".equals(control.toString()), "Ordinary growth must retain both bytes: "+control);
		pass("ordinary-growth", "ab");
		checkEmpty("empty-wrap", new ByteBuilder(new byte[0]), baseline);
		checkEmpty("empty-copy", new ByteBuilder(new ByteBuilder()), baseline);
		checkEmpty("empty-subsequence", (ByteBuilder)new ByteBuilder().subSequence(0, 0), baseline);
		final ByteBuilder shrunk=new ByteBuilder();
		shrunk.shrinkTo(0);
		checkEmpty("empty-shrink", shrunk, baseline);
		final String seventeen=new ByteBuilder().appendSci(1.0, 17).toString();
		check("1.00000000000000000e+00".equals(seventeen),
				"Seventeen-digit control differs: "+seventeen);
		pass("scientific-seventeen", seventeen);
		checkScientificEighteen(baseline);
		if(!baseline){checkBulkGrowth();}
	}

	/** Exercises repaired bulk growth without entering the nonterminating baseline path. */
	private static void checkBulkGrowth(){
		final ByteBuilder builder=new ByteBuilder(new byte[0]);
		builder.append(new byte[] {'x', 'y'});
		check(builder.length()==2 && "xy".equals(builder.toString()),
				"Bulk growth from an empty wrapper must retain both bytes: "+builder);
		pass("empty-bulk-growth", "xy");
	}

	/** Checks finite single-byte growth from a supported empty representation. */
	private static void checkEmpty(final String name, final ByteBuilder builder, final boolean baseline){
		check(builder.length()==0, name+" must begin empty");
		RuntimeException failure=null;
		try{
			builder.append((byte)'x');
		}catch(RuntimeException e){failure=e;}
		if(baseline){
			check(failure!=null && "Array overflow: 0<=0".equals(failure.getMessage()),
					name+" expected the zero-capacity growth failure, got "+failure);
			pass(name, failure.getMessage());
		}else{
			check(failure==null, name+" should grow, got "+failure);
			check(builder.length()==1 && "x".equals(builder.toString()),
					name+" must contain exactly the appended byte: "+builder);
			pass(name, "x");
		}
	}

	/** Checks the highest documented precision using the existing finite assertion guard. */
	private static void checkScientificEighteen(final boolean baseline){
		final ByteBuilder builder=new ByteBuilder();
		AssertionError failure=null;
		try{
			builder.appendSci(1.0, 18);
		}catch(AssertionError e){failure=e;}
		if(baseline){
			check(failure!=null && failure.getMessage()!=null
					&& failure.getMessage().startsWith("appendSci did not converge: num=1.0, decimals=18,"),
					"Expected the scientific bound convergence failure, got "+failure);
			pass("scientific-eighteen", "convergence assertion");
		}else{
			check(failure==null, "Eighteen-digit formatting failed: "+failure);
			check("1.000000000000000000e+00".equals(builder.toString()),
					"Eighteen-digit representation differs: "+builder);
			pass("scientific-eighteen", builder.toString());
		}
	}

	/** Emits the observed result after its literal expectation has been checked. */
	private static void pass(final String name, final String value){
		System.out.println("PASS\t"+name+'\t'+value);
	}

	/** Fails independently of JVM assertions if a fixture expectation is violated. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
