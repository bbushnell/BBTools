package test.repeat;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import jgi.CrisprFinder;
import repeat.PalindromeFinder;
import structures.Crispr;

/** Literal checks of the direct CRISPR caller's retention and annotation behavior. */
public final class CrisprPalindromeCallerTest{
	/** Uses no files, networks or classifier models. */
	public static void main(String[] args){
		check("ATNNAT", true, true);
		check("AANNAA", true, false);
		check("AANNAA", false, true);
		System.out.println("CrisprPalindromeCallerTest PASS");
	}

	/** The two two-base repeat regions are fixed at 0..1 and 4..5. */
	private static void check(String sequence, boolean requirePalindrome, boolean retained){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		final ArrayList<Crispr> pairs=new ArrayList<Crispr>();
		final Crispr pair=new Crispr(0, 1, 4, 5);
		pairs.add(pair);
		final PalindromeFinder finder=new PalindromeFinder(1, 0, 0, 1, 0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
		CrisprFinder.findPalindromes(bases, pairs, finder, false, requirePalindrome, 1);
		require(pairs.size()==(retained ? 1 : 0), "Wrong CRISPR pair retention for "+sequence);
		if(sequence.equals("ATNNAT")){
			require(pair.pa!=null && pair.pb!=null && pair.pa.a==0 && pair.pa.b==1
					&& pair.pb.a==4 && pair.pb.b==5 && finder.tracker.found==1,
					"Both exact palindrome annotations and one pair-level tracker record are required");
		}else{require(pair.pa==null && pair.pb==null, "Noncomplementary AA arms must not be annotated");}
	}

	/** Retains independent acceptance checks under both assertion modes. */
	private static void require(boolean valid, String message){if(!valid){throw new AssertionError(message);}}
}
