package test.parse;

import parse.Parse;
import parse.PreParser;

/** Regression tests for parse package scalar helpers. */
public final class ParseRegressionTest {

	public static void main(String[] args){
		int failures=0;
		failures+=check("parseIntHexDecOctBin decimal 0", new TestCase(){
			@Override
			public void run(){checkHexDecOctBin("0", 0);}
		});
		failures+=check("parseIntHexDecOctBin decimal 16", new TestCase(){
			@Override
			public void run(){checkHexDecOctBin("16", 16);}
		});
		failures+=check("parseIntHexDecOctBin hex lowercase", new TestCase(){
			@Override
			public void run(){checkHexDecOctBin("0x10", 16);}
		});
		failures+=check("parseIntHexDecOctBin hex uppercase", new TestCase(){
			@Override
			public void run(){checkHexDecOctBin("0X1f", 31);}
		});
		failures+=check("parseIntHexDecOctBin binary lowercase", new TestCase(){
			@Override
			public void run(){checkHexDecOctBin("0b1010", 10);}
		});
		failures+=check("parseIntHexDecOctBin binary uppercase", new TestCase(){
			@Override
			public void run(){checkHexDecOctBin("0B11", 3);}
		});
		failures+=check("parseIntHexDecOctBin octal lowercase", new TestCase(){
			@Override
			public void run(){checkHexDecOctBin("0o10", 8);}
		});
		failures+=check("parseIntHexDecOctBin octal uppercase", new TestCase(){
			@Override
			public void run(){checkHexDecOctBin("0O7", 7);}
		});
		failures+=check("PreParser.isAmino ignores amino-prefixed flags", new TestCase(){
			@Override
			public void run(){checkIsAminoIgnoresOtherAminoPrefixedFlags();}
		});
		if(failures>0){throw new AssertionError("ParseRegressionTest failures="+failures);}
		System.out.println("PASS ParseRegressionTest");
	}

	private static int check(String name, TestCase test){
		try{
			test.run();
			return 0;
		}catch(Throwable t){
			System.err.println("FAIL "+name);
			t.printStackTrace(System.err);
			return 1;
		}
	}

	private static void checkHexDecOctBin(String text, int expected){
		final int observed=Parse.parseIntHexDecOctBin(text);
		if(observed!=expected){
			throw new AssertionError("parseIntHexDecOctBin("+text+") returned "+observed+", expected "+expected);
		}
	}

	private static void checkIsAminoIgnoresOtherAminoPrefixedFlags(){
		if(PreParser.isAmino(new String[] {"aminopath=left=right"})){
			throw new AssertionError("Unrelated amino-prefixed argument was consumed as amino=true");
		}
	}

	private static interface TestCase {
		void run();
	}
}
