package test.parse;

import java.nio.charset.StandardCharsets;

import parse.Parse;

/** Regression tests for byte/String ASCII-insensitive comparison contracts. */
public final class ParseByteStringCaseTest{

	/** Runs literal correctness oracles for both byte/String overload orders.
	 * @param args Optional implementation=production|guarded|table|shorttable|hybrid */
	public static void main(String[] args){
		String implementation="production";
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Bad argument "+arg);}
			if(split[0].equals("implementation")){implementation=split[1];}
			else{throw new IllegalArgumentException("Unknown argument "+arg);}
		}
		final ParseByteStringCaseTest test=new ParseByteStringCaseTest(candidate(implementation));
		test.check("ASCII single-character oracle", () -> test.checkAsciiPairs());
		test.check("punctuation remains literal", () -> {
			test.expectFalse("@", "`");
			test.expectFalse("[", "{");
			test.expectFalse("\\", "|");
			test.expectFalse("^", "~");
			test.expectFalse("_", "\u007f");
		});
		test.check("null and length rules", () -> {
			test.expectFalse((byte[])null, null);
			test.expectFalse((String)null, (byte[])null);
			test.expectFalse((byte[])null, "");
			test.expectFalse("", (byte[])null);
			test.expectFalse(ascii("A"), (String)null);
			test.expectFalse((String)null, ascii("A"));
			test.expectTrue(ascii(""), "");
			test.expectTrue("", ascii(""));
			test.expectFalse(ascii("A"), "");
			test.expectFalse("", ascii("A"));
			test.expectFalse(ascii("A"), "AA");
			test.expectFalse("AA", ascii("A"));
		});
		test.check("mixed ASCII sequences", () -> {
			test.expectTrue("AbZ09!~", "aBz09!~");
			test.expectTrue("case-Fold_128", "CASE-fOLD_128");
			test.expectTrue("same-prefixAb", "same-prefixaB");
			test.expectFalse("AbZ", "aBy");
			test.expectFalse("prefixA", "prefixb");
			test.expectFalse("prefixA", "prefix\u00c1");
			test.expectFalse("a@", "A`");
		});
		test.check("non-ASCII rejected", () -> {
			test.expectFalse(new byte[] {(byte)0x80}, "\u0080");
			test.expectFalse(new byte[] {(byte)0xff}, "\u00ff");
			test.expectFalse(new byte[] {(byte)0xc1}, "A");
			test.expectFalse(ascii("A"), "\u00c1");
			test.expectFalse("A", new byte[] {(byte)0xc1});
			test.expectFalse("\u00c1", ascii("A"));
		});
		System.out.println("BYTE_STRING_CASE_RESULT implementation="+implementation+
				" cases="+test.cases+" failures="+test.failures);
		if(test.failures!=0){
			throw new AssertionError("Byte/String case failures="+test.failures);
		}
	}

	/** Constructor for a selected implementation. */
	private ParseByteStringCaseTest(Candidate candidate_){
		candidate=candidate_;
	}

	/** Checks every ASCII byte/char pair against Java's String oracle. */
	private void checkAsciiPairs(){
		int mismatches=0;
		for(int a=0; a<128; a++){
			final byte[] bytes={(byte)a};
			for(int b=0; b<128; b++){
				final String s=String.valueOf((char)b);
				final boolean expected=String.valueOf((char)a).equalsIgnoreCase(s);
				if(matches(bytes, s)!=expected){mismatches++;}
				if(matches(s, bytes)!=expected){mismatches++;}
			}
		}
		require(mismatches==0, "ASCII mismatches="+mismatches);
	}

	/** Checks both byte/String overload orders should match. */
	private void expectTrue(String bytes, String string){
		expectTrue(ascii(bytes), string);
		expectTrue(string, ascii(bytes));
	}

	/** Checks the byte-array-first form should match. */
	private void expectTrue(byte[] bytes, String string){
		require(matches(bytes, string), "expected match");
	}

	/** Checks the string-first form should match. */
	private void expectTrue(String string, byte[] bytes){
		require(matches(string, bytes), "expected match");
	}

	/** Checks both byte/String overload orders should reject a pair. */
	private void expectFalse(String bytes, String string){
		expectFalse(ascii(bytes), string);
		expectFalse(string, ascii(bytes));
	}

	/** Checks the byte-array-first form should reject a pair. */
	private void expectFalse(byte[] bytes, String string){
		require(!matches(bytes, string), "unexpected match");
	}

	/** Checks the string-first form should reject a pair. */
	private void expectFalse(String string, byte[] bytes){
		require(!matches(string, bytes), "unexpected match");
	}

	/** Invokes byte/String production or the selected prototype. */
	private boolean matches(byte[] bytes, String string){
		return candidate.matches(bytes, string);
	}

	/** Invokes real String/byte production or mirrors the prototype helper. */
	private boolean matches(String string, byte[] bytes){
		return candidate.production ? Parse.equalsIgnoreCase(string, bytes) : candidate.matches(bytes, string);
	}

	/** Selects an implementation to test. */
	private static Candidate candidate(String implementation){
		if(implementation.equals("production")){
			return new Candidate(true){
				@Override
				boolean matches(byte[] bytes, String string){return Parse.equalsIgnoreCase(bytes, string);}
			};
		}
		if(implementation.equals("guarded")){
			return new Candidate(false){
				@Override
				boolean matches(byte[] bytes, String string){
					return ParseByteStringCaseTimingTest.equalsGuarded(bytes, string);
				}
			};
		}
		if(implementation.equals("table")){
			return new Candidate(false){
				@Override
				boolean matches(byte[] bytes, String string){
					return ParseByteStringCaseTimingTest.equalsTable(bytes, string);
				}
			};
		}
		if(implementation.equals("shorttable")){
			return new Candidate(false){
				@Override
				boolean matches(byte[] bytes, String string){
					return ParseByteStringCaseTimingTest.equalsShortTable(bytes, string);
				}
			};
		}
		if(implementation.equals("hybrid")){
			return new Candidate(false){
				@Override
				boolean matches(byte[] bytes, String string){
					return ParseByteStringCaseTimingTest.equalsHybrid(bytes, string);
				}
			};
		}
		throw new IllegalArgumentException("Unknown implementation "+implementation);
	}

	/** ASCII bytes for a fixture literal. */
	private static byte[] ascii(String s){
		return s.getBytes(StandardCharsets.US_ASCII);
	}

	/** Records failures independently of whether Java assertions are enabled. */
	private void check(String name, Runnable action){
		cases++;
		try{action.run();}
		catch(RuntimeException | AssertionError e){
			failures++;
			System.out.println("FAIL "+name+": "+e);
		}
	}

	/** Fails an expectation even when the JVM runs with -da. */
	private static void require(boolean condition, String reason){
		if(!condition){throw new AssertionError(reason);}
	}

	/** Base for a selected comparison implementation. */
	private abstract static class Candidate{
		Candidate(boolean production_){production=production_;}
		abstract boolean matches(byte[] bytes, String string);
		final boolean production;
	}

	/** Selected implementation. */
	private final Candidate candidate;
	/** Number of attempted independent cases. */
	private int cases;
	/** Number of failed cases. */
	private int failures;
}
