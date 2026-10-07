package stream;

import java.util.Arrays;

/** Historical printing diagnostic for MDWalker and MDWalker2.
 * The five fixtures are examples, not validated regression oracles. Walker results
 * or caught Throwable details go to stdout; no equivalence or success is asserted.
 * Output is intended for manual inspection. */
public class TestMDWalker{

	/** Prints five fixed comparisons; command-line arguments are ignored. */
	public static void main(String[] args){
		// 1) Simple substitution after 3 matches: MD:Z:3A2
		runCase("simple-sub", "MD:Z:3A2", "6M", "mmmmmm");

		// 2) Insertion before substitution should be skipped by match counter
		// CIGAR: 2M1I3M, MD: 2A3
		runCase("ins-before-sub", "MD:Z:2A3", "2M1I3M", "mmImmm");

		// 3) Deletion before substitution: 2M2D1M1T2M → MD: 2^CC1T2
		runCase("del-before-sub", "MD:Z:2^CC1T2", "2M2D1M1X2M", "mmDDmSmm");

		// 4) Leading clipping should be skipped safely
		runCase("leading-clip", "MD:Z:1A1", "1S1M1X1M", "C m S m".replace(" ", ""));

		// 5) Mixed I/D around matches
		runCase("mixed-indels", "MD:Z:1A1^G2", "1M1I1M1D2M", "mImDmm");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Encodes fixture text with the platform default charset, as in the historical driver. */
	private static byte[] lm(String s){return s.getBytes();}

	/** Runs the two walkers independently on copied match arrays and prints their results.
	 * Passes the MD/CIGAR strings unchanged and null query bases to fixMatch. Each
	 * Throwable is caught and printed independently; no equality or validity is asserted.
	 * Prints the original match text after both attempts. */
	private static void runCase(String name, String md, String cigar, String longmatch){
		System.out.println("== "+name+" ==");
		byte[] lm0=lm(longmatch);
		byte[] lm1=Arrays.copyOf(lm0, lm0.length);
		byte[] lm2=Arrays.copyOf(lm0, lm0.length);

		try{
			MDWalker w=new MDWalker(md, cigar, lm1, null);
			w.fixMatch(null);
			System.out.println("MDWalker : "+new String(lm1));
		}catch(Throwable t){
			System.out.println("MDWalker : threw "+t.getClass().getSimpleName()+" - "+t.getMessage());
		}

		try{
			MDWalker2 w2=new MDWalker2(md, cigar, lm2, null);
			w2.fixMatch(null);
			System.out.println("MDWalker2: "+new String(lm2));
		}catch(Throwable t){
			System.out.println("MDWalker2: threw "+t.getClass().getSimpleName()+" - "+t.getMessage());
		}

		System.out.println("Start    : "+longmatch);
		System.out.println();
	}
}
