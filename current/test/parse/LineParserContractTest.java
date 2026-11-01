package test.parse;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import parse.LineParser;
import parse.LineParser1;
import parse.LineParser2;
import parse.LineParser4;
import parse.LineParserS;
import parse.LineParserS1;
import parse.LineParserS2;
import parse.LineParserS3;
import parse.LineParserS4;
import parse.LineParserS4Reverse;
import parse.LineParserSimple;
import shared.Shared;
import simd.Vector;

/** Independent literal-field and character contracts for the LineParser family.
 * Run through the team's committed Dori recipe with assertions enabled and disabled.
 * This is a correctness fixture, not a timing measurement.
 * @author Brian Bushnell, Sucrose */
public final class LineParserContractTest {

	/** Runs all cases, retaining individual failures before failing the process.
	 * @param args Optional simd=t or simd=f; true requires available SIMD support */
	public static void main(String[] args){
		boolean simd=false;
		for(String arg : args){
			if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(simd && !Vector.simd128){throw new IllegalStateException("SIMD requested but unavailable");}
		final boolean saved=Shared.SIMD;
		Shared.SIMD=simd;
		final LineParserContractTest test=new LineParserContractTest();
		try{test.run();}finally{Shared.SIMD=saved;}
		System.out.println("LINE_PARSER_RESULT cases="+test.cases+" failures="+test.failures+" simd="+simd);
		if(test.failures>0){throw new AssertionError("LineParser contract failures="+test.failures);}
	}

	/** Executes cases with distinct literal expectations; production parsers are not oracles. */
	private void run(){
		final LineParser[] forward={new LineParser1(','), new LineParser2(','),
			new LineParserS1(','), new LineParserS2(','), new LineParserS3(',')};
		for(final LineParser p : forward){
			final String name=p.getClass().getSimpleName();
			check(name+" empty fields", () -> fields(p, ",a,,b,", "", "a", "", "b", ""));
			check(name+" empty line", () -> fields(p, "", ""));
			check(name+" numeric fields", () -> {
				set(p, "17,-29,1.5,2.25");
				require(p.parseInt(0)==17, "int field");
				require(p.parseLong(1)==-29, "long field");
				require(p.parseFloat(2)==1.5f, "float field");
				require(p.parseDouble(3)==2.25, "double field");
			});
			check(name+" reuse", () -> {
				fields(p, "one,two,three", "one", "two", "three");
				p.clear();
				fields(p, "x", "x");
			});
		}
		check("byte ordered delimiters", () -> fields(new LineParser4(":/"), "a:b/c:d", "a", "b", "c:d"));
		check("string ordered delimiters", () -> fields(new LineParserS4(":/"), "a:b/c:d", "a", "b", "c:d"));
		check("reverse preserves prefix", () -> fields(new LineParserS4Reverse(":/"), "a:b:c/d", "a:b", "c", "d"));
		check("reverse null empty field", () -> {
			final LineParserS4Reverse p=new LineParserS4Reverse(":/").set("a:/d");
			require(p.parseString(1)==null, "reverse parser deliberately returns null for empty field");
		});
		final LineParserS[] strings={new LineParserS1(','), new LineParserS2(','),
			new LineParserS3(','), new LineParserS4(","), new LineParserS4Reverse(",")};
		for(final LineParserS p : strings){
			for(final char c : new char[]{'A', '\u00e9', '\u0100', '\u4e2d'}){
				check(p.getClass().getSimpleName()+" char "+(int)c, () -> {
					p.set("x,"+c);
					require(p.parseChar(1, 0)==c, "parseChar must retain the original UTF-16 code unit "+(int)c);
				});
			}
		}
		check("S3 full whitespace index", () -> {
			final LineParserS3 p=new LineParserS3(':').set("a:b c:d");
			require(p.indexOfWhitespace()==1, "space terminates field 1");
		});
		check("S3 limited whitespace index", () -> {
			final LineParserS3 p=new LineParserS3(':').set("a:b c:d", 2);
			require(p.indexOfWhitespace()==1, "space terminating field 1 is inside the parsed prefix");
		});
		check("S3 non-ASCII delimiter is not whitespace", () -> {
			final LineParserS3 p=new LineParserS3('\u0100').set("a\u0100b c");
			require(p.indexOfWhitespace()==1, "primary delimiter at field 0 is not a space/tab");
		});
		check("S3 limited no whitespace control", () -> {
			final LineParserS3 p=new LineParserS3(':').set("a:b:c:d", 2);
			require(p.indexOfWhitespace()==-1, "no whitespace in parsed prefix");
		});
		check("S3 primary space before tab", () -> {
			final LineParserS3 p=new LineParserS3(' ').set("a b c\td");
			require(p.indexOfWhitespace()==2, "primary spaces do not count; tab terminates field 2");
		});
		check("S3 primary tab before space", () -> {
			final LineParserS3 p=new LineParserS3('\t').set("a\tb\tc d");
			require(p.indexOfWhitespace()==2, "primary tabs do not count; space terminates field 2");
		});
		check("S3 whitespace outside limited prefix", () -> {
			final LineParserS3 p=new LineParserS3(':').set("a:b c:d", 0);
			require(p.indexOfWhitespace()==-1, "unscanned whitespace must not be reported");
		});
		check("S3 whitespace cache resets on reuse", () -> {
			final LineParserS3 p=new LineParserS3(':').set("a b:c");
			require(p.indexOfWhitespace()==0, "first line has whitespace at field 0");
			p.set("a:b:c");
			require(p.indexOfWhitespace()==-1, "second line has only primary delimiters");
		});
		check("LP1 long delimiter positions", () -> {
			final byte[] data=new byte[257];
			Arrays.fill(data, (byte)'x');
			final int[] separators={0, 15, 16, 31, 32, 63, 64, 127, 128, 255, 256};
			for(int i : separators){data[i]=',';}
			final int[] lengths={0, 14, 0, 14, 0, 30, 0, 62, 0, 126, 0, 0};
			final LineParser1 p=new LineParser1(',').set(data);
			require(p.terms()==lengths.length, "all delimiters including final empty field");
			for(int i=0; i<lengths.length; i++){
				require(p.length(i)==lengths[i], "field length at "+i);
			}
		});
		check("simple bounded scan and reset", () -> {
			final LineParserSimple p=new LineParserSimple((byte)',', "a,,bc".getBytes(StandardCharsets.US_ASCII));
			require(p.advance()==1 && p.advance()==0 && p.advance()==2, "literal field lengths 1,0,2");
			p.reset();
			require(p.advance()==1, "reset restores first field");
		});
	}

	/** Selects the native setter so string parsers never exercise their discouraged byte adapters. */
	private static void set(LineParser p, String line){
		if(p instanceof LineParserS){((LineParserS)p).set(line);}
		else{p.set(line.getBytes(StandardCharsets.US_ASCII));}
	}

	/** Checks literal fields in increasing order, compatible with forward-only variants. */
	private static void fields(LineParser p, String line, String... expected){
		set(p, line);
		for(int i=0; i<expected.length; i++){
			final String actual=p.parseString(i);
			require(expected[i].equals(actual), "field "+i+" expected="+expected[i]+" actual="+actual);
		}
	}

	/** Records each failed contract in either JVM assertion mode and continues the remaining cases. */
	private void check(String name, Runnable action){
		cases++;
		try{action.run();}
		catch(AssertionError | RuntimeException e){
			failures++;
			System.out.println("FAIL "+name+": "+e);
		}
	}

	/** Uses an unconditional failure so disabled JVM assertions cannot hide test failures. */
	private static void require(boolean condition, String explanation){
		if(!condition){throw new AssertionError(explanation);}
	}

	/** Number of independently reported cases attempted. */
	private int cases;
	/** Number of cases that failed or threw unexpectedly. */
	private int failures;
}
