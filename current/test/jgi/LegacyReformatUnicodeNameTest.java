package test.jgi;

import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;

import jgi.ReformatReads;
import jgi.ReformatReadsOld;

/**
 * Literal Unicode and ASCII controls for both legacy String conversion helpers.
 * The baseline uses UTF-8 because the legacy code indexes encoded bytes with a
 * UTF-16 offset. Correct output replaces the final underscore without changing
 * any other characters. No-underscore calls must return the original reference.
 * Only public static String helpers are called; no driver, Read, parser or worker
 * is created. Null Strings and fixHeader policy are outside this test.
 * @author Sucrose
 */
public final class LegacyReformatUnicodeNameTest{

	/** Requires phase=baseline|fixed; the launcher explicitly selects UTF-8. */
	public static void main(final String[] args){
		if(args.length!=1 || (!args[0].equals("phase=baseline") && !args[0].equals("phase=fixed"))){
			throw new IllegalArgumentException("Set phase=baseline|fixed.");
		}
		final boolean baseline=args[0].equals("phase=baseline");
		final Charset charset=Charset.defaultCharset();
		check(charset.equals(StandardCharsets.UTF_8),
			"The literal broken outputs require UTF-8 because legacy String helpers use the default charset: "+charset);
		System.out.println("legacy_unicode_encoding\tdefault="+charset.name()+
			"\tfile.encoding="+System.getProperty("file.encoding"));
		int mismatches=0, identities=0;
		for(final boolean old : new boolean[]{false, true}){
			for(int i=0; i<CASES.length; i++){
				final String[] row=CASES[i];
				final String label=(old ? "legacyold-" : "legacy-")+row[0];
				final String input=new String(row[1]);
				final String expected=baseline ? row[2] : row[3];
				final String actual=old ? ReformatReadsOld.underscoreToSpace(input) : ReformatReads.underscoreToSpace(input);
				check(expected.equals(actual), label+": expected "+display(expected)+", observed "+display(actual));
				if(i<IDENTITY_CASES){
					check(actual==input, label+": a no-underscore input must return its original String reference.");
					identities++;
				}
				final boolean matches=row[3].equals(actual);
				if(!matches){mismatches++;}
				System.out.println("legacy_unicode_case\t"+label+"\tinput="+display(input)+
					"\tactual="+display(actual)+"\tintended="+display(row[3])+
					"\tmatchesIntended="+matches+"\tsameReference="+(actual==input));
			}
		}
		check(mismatches==(baseline ? 8 : 0),
			"Four literal Unicode offset cases per legacy helper must differ only before the repair: "+mismatches);
		if(checks!=34 || identities!=8){
			throw new AssertionError("Frozen Unicode fixture requires 34 checks/eight identity checks: "+checks+"/"+identities);
		}
		System.out.println("LEGACY_REFORMAT_UNICODE_NAME_PASS phase="+(baseline ? "baseline" : "fixed")+
			" checks="+checks+" cases=24 mismatches="+mismatches+" identityChecks="+identities);
	}

	/** ASCII-only observations preserve exact UTF-16 units, including supplementary pairs. */
	private static String display(final String value){
		if(value==null){return "<null>";}
		final StringBuilder out=new StringBuilder("[");
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(c>=32 && c<127 && c!='\\'){out.append(c);}
			else{
				out.append("\\u");
				for(int shift=12; shift>=0; shift-=4){out.append(HEX[(c>>>shift)&15]);}
			}
		}
		return out.append(']').toString();
	}

	/** Literal assertions remain active independently of the JVM assertion setting. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/** Label, input, literal UTF-8 baseline output, intended output; first four require identity. */
	private static final String[][] CASES={
		{"empty", "", "", ""},
		{"plain-ascii", "plain", "plain", "plain"},
		{"plain-bmp", "plain\u00e9", "plain\u00e9", "plain\u00e9"},
		{"plain-supplementary", "plain\ud83d\ude00", "plain\ud83d\ude00", "plain\ud83d\ude00"},
		{"ascii-one", "alpha_beta", "alpha beta", "alpha beta"},
		{"ascii-last", "alpha_beta_tail", "alpha_beta tail", "alpha_beta tail"},
		{"unicode-after", "head_\u00e9", "head \u00e9", "head \u00e9"},
		{"leading-underscore", "_\u6f22", " \u6f22", " \u6f22"},
		{"bmp-last", "\u00e9_a_b", "\u00e9_ _b", "\u00e9_a b"},
		{"bmp-prefix", "\u00e9ab_x", "\u00e9a _x", "\u00e9ab x"},
		{"cjk-prefix", "\u6f22ab_x", "\u6f22 b_x", "\u6f22ab x"},
		{"supplementary-prefix", "\ud83d\ude00ab_x", "\ud83d\ude00 b_x", "\ud83d\ude00ab x"}
	};
	private static final int IDENTITY_CASES=4;
	private static final char[] HEX="0123456789abcdef".toCharArray();
	private static int checks;
}
