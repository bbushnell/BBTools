package test.shared;

import java.math.BigInteger;
import java.util.Arrays;

import shared.Primes;

/** Finite large-capacity checks against the real sparse prime resource.
 * Inequality oracles use BigInteger, so their own products cannot overflow.
 * The fallback remains a product of table entries, not a claim that every result
 * is prime. Near-MAX overflow is classified separately using the legacy first
 * factor and the smallest table second factor sufficient for the requested bound.
 * @author Sucrose
 */
public final class PrimeBoundsTest{

	/** Runs ordinary table controls and six positive large requests in both directions. */
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("Expected no arguments");}
		final long[] table=Primes.primes;
		require(table.length>1 && table[0]>0, "Nonempty positive prime table required");
		for(int i=1; i<table.length; i++){require(table[i]>table[i-1], "Prime table must be strictly ascending");}
		final PrimeBoundsTest test=new PrimeBoundsTest();
		test.check("table/int-capacity controls", () -> {
			for(long request : new long[] {table[0], 1000, Integer.MAX_VALUE}){
				require(Primes.primeAtLeast(request)==ceiling(table, request), "Sparse ceiling changed for "+request);
				final int index=Arrays.binarySearch(table, request);
				final long expected=table[index>=0 ? index : -index-2];
				require(Primes.primeAtMost(request)==expected, "Sparse floor changed for "+request);
			}
		});
		final long pivot=1L<<61;
		for(final long request : new long[] {pivot-1, pivot, pivot+1, Long.MAX_VALUE-(1L<<32), Long.MAX_VALUE-1, Long.MAX_VALUE}){
			require(request>table[table.length-1], "Fixture request must exercise product fallback: "+request);
			test.check("at-least "+request, () -> atLeast(request, table));
			test.check("at-most "+request, () -> atMost(request));
		}
		System.out.println("PrimeBoundsTest cases="+test.cases+" failures="+test.failures);
		if(test.failures>0){throw new AssertionError("Prime bound failures="+test.failures);}
	}

	/** Checks representable upper results, or requires explicit overflow rejection.
	 * ceil(request/a) is computed exactly: using floor division here would allow
	 * a product below request to masquerade as a valid at-least capacity.
	 */
	private static void atLeast(long request, long[] table){
		final BigInteger x=BigInteger.valueOf(request);
		final long divisor=(long)Math.pow(request, 0.4);
		final long a=ceiling(table, request/divisor);
		final BigInteger bigA=BigInteger.valueOf(a);
		final long requiredB=x.add(bigA).subtract(BigInteger.ONE).divide(bigA).longValueExact();
		final long b=ceiling(table, requiredB);
		final BigInteger minimumProduct=bigA.multiply(BigInteger.valueOf(b));
		System.out.println("upper classification request="+request+" first="+a+" second="+b+" exactProduct="+minimumProduct);
		if(minimumProduct.compareTo(LONG_MAX)>0){
			boolean rejected=false;
			try{
				final long actual=Primes.primeAtLeast(request);
				System.out.println("unrepresentable upper returned="+actual);
			}catch(ArithmeticException expected){
				rejected=expected.getMessage()!=null && !expected.getMessage().isEmpty();
			}
			require(rejected, "Unrepresentable table-product capacity must fail explicitly, not wrap");
			return;
		}
		final long actual=Primes.primeAtLeast(request);
		final BigInteger result=BigInteger.valueOf(actual), upper=x.multiply(NINE).divide(EIGHT);
		require(actual>0 && result.compareTo(x)>=0 && result.compareTo(upper)<0,
				"Expected "+x+" <= result < "+upper+", actual="+actual);
	}

	/** Checks the lower result against the exact strict seven-eighths bound. */
	private static void atMost(long request){
		final long actual=Primes.primeAtMost(request);
		final BigInteger x=BigInteger.valueOf(request), result=BigInteger.valueOf(actual);
		final BigInteger lower=x.multiply(SEVEN).divide(EIGHT);
		require(actual>0 && result.compareTo(x)<=0 && result.compareTo(lower)>0,
				"Expected "+lower+" < result <= "+x+", actual="+actual);
	}

	/** Selects a table entry directly; never calls the production bound methods for an oracle. */
	private static long ceiling(long[] table, long request){
		final int found=Arrays.binarySearch(table, request), index=found>=0 ? found : -found-1;
		require(index<table.length, "Oracle factor exceeds the real table: "+request);
		return table[index];
	}

	/** Reports every finite case so assertion-enabled baseline failures do not mask controls. */
	private void check(String name, Runnable action){
		cases++;
		try{action.run();System.out.println(name+" PASS");}
		catch(RuntimeException | AssertionError failure){failures++;System.out.println(name+" FAIL: "+failure);}
	}

	/** Checks the independent mathematical oracle with JVM assertions either enabled or disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Attempted and failed checks. */
	private int cases, failures;
	/** Exact arithmetic constants. */
	private static final BigInteger SEVEN=BigInteger.valueOf(7), EIGHT=BigInteger.valueOf(8), NINE=BigInteger.valueOf(9);
	/** Largest representable positive long capacity. */
	private static final BigInteger LONG_MAX=BigInteger.valueOf(Long.MAX_VALUE);
}
