package test.server;

import server.PercentEncoding;

/** Pure finite checks for server.PercentEncoding; no network or server launch. */
public class PercentEncodingTest {

	public static void main(String[] args){
		if(args.length!=1 || !(args[0].equals("baseline") || args[0].equals("fixed"))){
			throw new IllegalArgumentException("Expected baseline or fixed");
		}
		final boolean fixed=args[0].equals("fixed");

		checkAsciiControls();
		checkReservedRoundTrip();
		checkCommonRoundTrip();
		checkMalformedPercentLiterals();
		checkUnicodeCases(fixed);
	}

	private static void checkAsciiControls(){
		expect("plain-symbol-unchanged",
				"AZaz09._~-".equals(PercentEncoding.symbolToCode("AZaz09._~-")), true);
		expect("plain-common-unchanged",
				"AZaz09._~-/?:@".equals(PercentEncoding.commonSymbolToCode("AZaz09._~-/?:@")), true);
		expect("plain-decode-unchanged",
				"AZaz09._~-".equals(PercentEncoding.codeToSymbol("AZaz09._~-")), true);
		expect("ascii-decode-before-plain-hex",
				" abc+z".equals(PercentEncoding.codeToSymbol("%20abc+z")), true);
		expect("null-special-false", !PercentEncoding.containsSpecialSymbol(null), true);
		expect("null-common-false", !PercentEncoding.containsCommonSymbol(null), true);
	}

	private static void checkReservedRoundTrip(){
		final String symbols="!#$&'()*+,/:;=?@[]";
		final String codes="%21%23%24%26%27%28%29%2A%2B%2C%2F%3A%3B%3D%3F%40%5B%5D";
		expect("reserved-encode", codes.equals(PercentEncoding.symbolToCode(symbols)), true);
		expect("reserved-decode", symbols.equals(PercentEncoding.codeToSymbol(codes)), true);
		expect("lowercase-hex-decode", "/?<>\\".equals(PercentEncoding.codeToSymbol("%2f%3f%3c%3e%5c")), true);
	}

	private static void checkCommonRoundTrip(){
		final String symbols="\n \"<>\\|%";
		final String codes="%0A%20%22%3C%3E%5C%7C%25";
		expect("common-encode", codes.equals(PercentEncoding.commonSymbolToCode(symbols)), true);
		expect("common-decode", symbols.equals(PercentEncoding.codeToSymbol(codes)), true);
		expect("percent-is-escaped-on-second-pass",
				"%2520".equals(PercentEncoding.commonSymbolToCode("%20")), true);
	}

	/** Malformed percent spans are an intentional literal-preservation tolerance. */
	private static void checkMalformedPercentLiterals(){
		expect("malformed-bare-percent", "%".equals(PercentEncoding.codeToSymbol("%")), true);
		expect("malformed-short-percent", "%F".equals(PercentEncoding.codeToSymbol("%F")), true);
		expect("malformed-bad-high-nybble", "%GG".equals(PercentEncoding.codeToSymbol("%GG")), true);
		expect("malformed-bad-low-nybble", "%2G".equals(PercentEncoding.codeToSymbol("%2G")), true);
		expect("malformed-mixed-valid",
				"a%2Gb%zz c%".equals(PercentEncoding.codeToSymbol("a%2Gb%zz%20c%")), true);
	}

	private static void checkUnicodeCases(boolean fixed){
		expect("unicode-full-no-crash", !throwsRuntime(new Runnable(){
			@Override
			public void run(){PercentEncoding.symbolToCode("\u00e9?");}
		}), fixed);
		expect("unicode-common-preserved",
				"\u00e9\u6f22%20".equals(tryCommonSymbolToCode("\u00e9\u6f22 ")), fixed);
		expect("unicode-full-preserved",
				"\u00e9\u6f22%3F".equals(trySymbolToCode("\u00e9\u6f22?")), fixed);
		expect("unicode-decode-preserved",
				"\u00e9\u6f22 ".equals(PercentEncoding.codeToSymbol("\u00e9\u6f22%20")), fixed);
		expect("utf8-decode-preserved",
				"\u00e9\u6f22 ".equals(PercentEncoding.codeToSymbol("%C3%A9%E6%BC%A2%20")), fixed);
		expect("invalid-utf8-decode-replaces",
				"\uFFFD(".equals(PercentEncoding.codeToSymbol("%C3%28")), true);
	}

	private static String trySymbolToCode(String s){
		try{return PercentEncoding.symbolToCode(s);}
		catch(RuntimeException e){return null;}
	}

	private static String tryCommonSymbolToCode(String s){
		try{return PercentEncoding.commonSymbolToCode(s);}
		catch(RuntimeException e){return null;}
	}

	private static boolean throwsRuntime(Runnable r){
		try{r.run();}
		catch(RuntimeException e){return true;}
		return false;
	}

	private static void expect(String name, boolean actual, boolean expected){
		System.out.println(name+"\t"+actual+"\t"+expected);
		if(actual!=expected){throw new AssertionError(name+": expected valid="+expected+", observed="+actual);}
	}
}
