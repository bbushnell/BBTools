package test.structures;

import repeat.Palindrome;
import structures.Crispr;

/** Literal annotation checks for the selected palindrome's enclosing range.
 * @author Shinobu
 */
public final class CrisprPalindromeAnnotationCase{

	/** Checks the known old second-range defect or its bounded correction. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation: "+args[0]);
		final Crispr crispr=new Crispr(0, 9, 20, 31);
		checkText("no-palindrome", "[0-9,20-31]", crispr.toString());

		// First interval 0..9 encloses palindrome 2..7 with two bases on each side.
		crispr.pa=new Palindrome(2, 7, 2, 0);
		checkText("first-only", "[0-9,20-31],(2-7,P=2,M=2,L=2,T=2+2)", crispr.toString());
		crispr.pb=new Palindrome(22, 29, 3, 0);
		checkText("first-preferred", "[0-9,20-31],(2-7,P=2,M=2,L=2,T=2+2)", crispr.toString());

		// Second interval 20..31 encloses 22..29; its tails are also two bases each.
		// The old use of 0..9 instead gives 22 and -20, printed in ascending order.
		crispr.pa=null;
		final String expected=(fixed ? "[0-9,20-31](22-29,P=3,M=3,L=2,T=2+2)"
				: "[0-9,20-31](22-29,P=3,M=3,L=2,T=-20+22)");
		checkText("second-only", expected, crispr.toString());
		final byte[] bases=new byte[32];
		for(int i=0; i<bases.length; i++){bases[i]=(byte)(i<20 ? 'A' : 'C');}
		final String expectedBases=(fixed
				? "[0-9,20-31;32](22-29,P=3,M=3,L=2,T=2+2)\nAAAAAAAAAA\nCCCCCCCCCCCC"
				: "[0-9,20-31;32](22-29,P=3,M=3,L=2,T=-20+22)\nAAAAAAAAAA\nCCCCCCCCCCCC");
		checkText("second-with-bases", expectedBases, crispr.toString(bases));
	}

	/** Checks a literal string and prints an unambiguous single-line result. */
	private static void checkText(final String name, final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		System.out.println("PASS\t"+name+'\t'+observed.replace("\n", "\\n"));
	}

	/** Fails independently of assertion flags so an oracle cannot silently disappear. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
