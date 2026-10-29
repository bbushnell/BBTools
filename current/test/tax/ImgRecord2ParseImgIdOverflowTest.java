package test.tax;

import tax.ImgRecord2;

/** Finite baseline for unchecked ImgRecord2.parseImgId long overflow.
 * @author Sucrose */
public final class ImgRecord2ParseImgIdOverflowTest{

	/** Runs valid IMG header controls and overflow checks.
	 * @param args Optional expect=baseline or expect=fixed argument */
	public static void main(String[] args){
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		check(ImgRecord2.parseImgId("read|img|1|tid|7", true)==1,
				"Small pipe-delimited IMG ID should parse exactly.");
		check(ImgRecord2.parseImgId("read~img~123456789~tid~7", true)==123456789L,
				"Small tilde-delimited IMG ID should parse exactly.");
		check(ImgRecord2.parseImgId("read|img|9223372036854775807|tid|7", true)
				==Long.MAX_VALUE,
				"Long.MAX_VALUE should remain the largest accepted IMG ID.");
		check(ImgRecord2.parseImgId("read|img|000000000000000000001|tid|7", true)==1,
				"Leading zeroes should not make a small IMG ID overflow.");
		check(ImgRecord2.parseImgId("read|img|0009223372036854775807|tid|7", true)
				==Long.MAX_VALUE,
				"Leading zeroes should not make Long.MAX_VALUE overflow.");
		check(ImgRecord2.parseImgId("read|img|19 rest of header", true)==19,
				"Parsing should stop at the first space after the IMG token.");
		check(ImgRecord2.parseImgId("read without old-style delimiter", false)==-1,
				"Headers without an old-style IMG delimiter should return -1.");
		check(ImgRecord2.parseImgId("read|tid|7", false)==-1,
				"Headers without an IMG token should return -1.");
		check(ImgRecord2.parseImgId("read|img|0|tid|7", false)==-1,
				"Zero IMG ID should return -1.");

		if(baseline){
			check(ImgRecord2.parseImgId(MIN_VALUE_WRAPS_NEGATIVE, false)==-1,
					"Baseline should silently wrap Long.MAX_VALUE+1 nonpositive and return -1.");
			check(ImgRecord2.parseImgId(ONE_AFTER_UNSIGNED_LONG_MAX, false)==1,
					"Baseline should silently wrap 2^64+1 to positive IMG ID 1.");
			check(ImgRecord2.parseImgId(ONE_AFTER_UNSIGNED_LONG_MAX, true)==1,
					"Baseline should silently wrap 2^64+1 under assertions to positive IMG ID 1.");
		}else{
			checkOverflowRejected(MIN_VALUE_WRAPS_NEGATIVE, false);
			checkOverflowRejected(ONE_AFTER_UNSIGNED_LONG_MAX, false);
			checkOverflowRejected(ONE_AFTER_UNSIGNED_LONG_MAX, true);
			checkOverflowRejected(LEADING_ZERO_OVERFLOW, false);
		}

		System.out.println("IMGRECORD2_PARSE_IMG_ID_OVERFLOW_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private ImgRecord2ParseImgIdOverflowTest(){}

	/** Requires an overflowing IMG ID to fail loudly.
	 * @param header header containing an IMG ID larger than Long.MAX_VALUE
	 * @param doAssertions whether to request the parser's missing-token assertions */
	private static void checkOverflowRejected(final String header, final boolean doAssertions){
		try{
			ImgRecord2.parseImgId(header, doAssertions);
			throw new AssertionError("Overflowing IMG ID should throw NumberFormatException: "+
					header);
		}catch(NumberFormatException expected){
			//expected
		}
	}

	/** Test checks stay active under -da.
	 * @param condition condition required for test success
	 * @param message diagnostic message for failure */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Long.MAX_VALUE+1 wraps to Long.MIN_VALUE in the baseline accumulator. */
	private static final String MIN_VALUE_WRAPS_NEGATIVE=
			"read|img|9223372036854775808|tid|7";

	/** 2^64+1 wraps to the positive value 1 in the baseline accumulator. */
	private static final String ONE_AFTER_UNSIGNED_LONG_MAX=
			"read|img|18446744073709551617|tid|7";

	/** Long.MAX_VALUE+1 with leading zeroes still overflows the valid ID range. */
	private static final String LEADING_ZERO_OVERFLOW=
			"read|img|0009223372036854775808|tid|7";
}
