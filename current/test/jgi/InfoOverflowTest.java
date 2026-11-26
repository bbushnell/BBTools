package test.jgi;

import java.math.BigInteger;

import jgi.Info;

/** Checks Info's multinomial arithmetic around representable and saturating counts.
 * Expected values are validated with direct factorial arithmetic before production calls.
 * @author Brian Bushnell, Jean */
public final class InfoOverflowTest{

	/** Validates each literal oracle, then compares production helpers. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("Info overflow tests take no arguments.");}
		int failures=0;
		for(final Case test : CASES){
			validateExpected(test);
			final LongResult info=callInfo(test);
			final IntResult prefix=callPrefix(test);
			final boolean pass=info.matches(test.expectedInfo) && prefix.matches(test.expectedPrefix);
			System.out.println((pass ? "PASS " : "FAIL ")+test.name+
					" info="+info+" expected="+test.expectedInfo+
					" prefix="+prefix+" expected="+test.expectedPrefix);
			if(!pass){failures++;}
		}
		System.out.println("InfoOverflowTest cases="+CASES.length+" failures="+failures);
		if(failures>0){throw new AssertionError(failures+" failed cases");}
	}

	/** Static fixture entry point only. */
	private InfoOverflowTest(){}

	/** Validates the literal count and prefix endpoint independently of Info's recurrence. */
	private static void validateExpected(final Case test){
		final BigInteger exact=count(test.bases, test.from, test.len);
		final BigInteger max=BigInteger.valueOf(Long.MAX_VALUE);
		if(test.expectedInfo==Long.MAX_VALUE){
			require(exact.compareTo(max)>0,
					"Expected a genuinely saturating count for "+test.name+": "+exact);
		}else{
			require(exact.equals(BigInteger.valueOf(test.expectedInfo)),
					"Literal count mismatch for "+test.name+": "+exact+" != "+test.expectedInfo);
		}
		final int endpoint=expectedPrefix(test.bases, test.from, test.threshold);
		require(endpoint==test.expectedPrefix,
				"Literal prefix mismatch for "+test.name+": "+endpoint+" != "+test.expectedPrefix);
	}

	/** Computes the multinomial count for the requested slice by factorials. */
	private static BigInteger count(final byte[] bases, final int from, final int len){
		final int[] counts=new int[4];
		int used=0;
		for(int i=from, lim=from+len; i<lim; i++){
			final int number=number(bases[i]);
			if(number>=0){
				counts[number]++;
				used++;
			}
		}
		BigInteger result=factorial(used);
		for(final int c : counts){result=result.divide(factorial(c));}
		return result;
	}

	/** Finds the first absolute exclusive endpoint whose exact count reaches the threshold. */
	private static int expectedPrefix(final byte[] bases, final int from, final long threshold){
		final BigInteger target=BigInteger.valueOf(threshold);
		for(int end=from; end<=bases.length; end++){
			if(count(bases, from, end-from).compareTo(target)>=0){return end;}
		}
		return -1;
	}

	/** Returns a small factorial as a BigInteger. */
	private static BigInteger factorial(final int n){
		BigInteger result=BigInteger.ONE;
		for(int i=2; i<=n; i++){result=result.multiply(BigInteger.valueOf(i));}
		return result;
	}

	/** Maps only the bases present in these literals. */
	private static int number(final byte base){
		switch(base){
			case 'A': return 0;
			case 'C': return 1;
			case 'G': return 2;
			case 'T': return 3;
			default: return -1;
		}
	}

	/** Calls the production count helper while preserving per-case failures. */
	private static LongResult callInfo(final Case test){
		try{
			return new LongResult(Info.info(test.bases, test.from, test.len), null);
		}catch(RuntimeException e){
			return new LongResult(Long.MIN_VALUE, e);
		}
	}

	/** Calls the production prefix helper while preserving per-case failures. */
	private static IntResult callPrefix(final Case test){
		try{
			return new IntResult(Info.prefixForInfo(test.bases, test.threshold, test.from), null);
		}catch(RuntimeException e){
			return new IntResult(Integer.MIN_VALUE, e);
		}
	}

	/** Enforces expectations under any JVM assertion setting. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Test case with expected helper results. */
	private static final class Case{

		/** Stores the exact helper inputs and expected outputs. */
		Case(final String name_, final String bases_, final int from_, final int len_,
				final long threshold_, final long expectedInfo_, final int expectedPrefix_){
			name=name_;
			bases=bases_.getBytes();
			from=from_;
			len=len_;
			threshold=threshold_;
			expectedInfo=expectedInfo_;
			expectedPrefix=expectedPrefix_;
		}

		/** Stable case name. */
		final String name;
		/** Production input bases. */
		final byte[] bases;
		/** Slice start for both helpers. */
		final int from;
		/** Slice length for the info helper. */
		final int len;
		/** Threshold for the prefix helper. */
		final long threshold;
		/** Expected saturated multinomial count. */
		final long expectedInfo;
		/** Expected absolute exclusive prefix endpoint, or -1 if unreachable. */
		final int expectedPrefix;
	}

	/** Production long-result wrapper. */
	private static final class LongResult{

		/** Stores either a successful result or a runtime failure. */
		LongResult(final long value_, final RuntimeException exception_){
			value=value_;
			exception=exception_;
		}

		/** Whether this result equals the expected value without an exception. */
		boolean matches(final long expected){return exception==null && value==expected;}

		/** Stable value or exception rendering for baseline output. */
		@Override
		public String toString(){
			return exception==null ? Long.toString(value) : exception.getClass().getSimpleName();
		}

		/** Successful value, meaningful only when exception is null. */
		final long value;
		/** Production helper failure, if any. */
		final RuntimeException exception;
	}

	/** Production int-result wrapper. */
	private static final class IntResult{

		/** Stores either a successful result or a runtime failure. */
		IntResult(final int value_, final RuntimeException exception_){
			value=value_;
			exception=exception_;
		}

		/** Whether this result equals the expected value without an exception. */
		boolean matches(final int expected){return exception==null && value==expected;}

		/** Stable value or exception rendering for baseline output. */
		@Override
		public String toString(){
			return exception==null ? Integer.toString(value) : exception.getClass().getSimpleName();
		}

		/** Successful value, meaningful only when exception is null. */
		final int value;
		/** Production helper failure, if any. */
		final RuntimeException exception;
	}

	private static final String COUNT_33="AAAAAAACCCCCCCCGGGGGGGGTTTTTTTTTT";
	private static final String COUNT_34=COUNT_33+"A";
	private static final String OVERFLOW_36="ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
	private static final long COUNT_33_VALUE=292045871187813600L;
	private static final long COUNT_34_VALUE=1241194952548207800L;
	private static final Case[] CASES={
			new Case("empty", "", 0, 0, 1, 1, 0),
			new Case("AACC", "AACC", 0, 4, 6, 6, 4),
			new Case("count33", COUNT_33, 0, 33, COUNT_33_VALUE, COUNT_33_VALUE, 33),
			new Case("count34_2^60", COUNT_34, 0, 34, 1L<<60, COUNT_34_VALUE, 34),
			new Case("count34_2^61", COUNT_34, 0, 34, 1L<<61, COUNT_34_VALUE, -1),
			new Case("count34_exact", COUNT_34, 0, 34, COUNT_34_VALUE, COUNT_34_VALUE, 34),
			new Case("overflow36_MAX", OVERFLOW_36, 0, 36, Long.MAX_VALUE, Long.MAX_VALUE, 36),
			new Case("offset_count34_2^60", "NNN"+COUNT_34, 3, 34, 1L<<60, COUNT_34_VALUE, 37)
	};
}
