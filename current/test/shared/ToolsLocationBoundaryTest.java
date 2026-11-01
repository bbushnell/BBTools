package test.shared;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import shared.Shared;
import shared.Tools;
import stream.SiteScore;

/** Literal boundary checks for mapping-location counts and delimited token searches.
 * Preserves the existing same-chromosome endpoint heuristic and score cutoff.
 * @author Shinobu
 */
public final class ToolsLocationBoundaryTest{

	/** Runs finite cases without files, workers, or global reference data. */
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("Expected no arguments");}
		final ToolsLocationBoundaryTest test=new ToolsLocationBoundaryTest();
		test.check("different chromosomes", () -> {
			require(Tools.countTopScores(sites(2, Shared.PLUS, 100, 109, 100), 0)==2,
					"Equal coordinates on different chromosomes are distinct mapping locations");
		});
		test.check("mapping controls", () -> {
			require(Tools.countTopScores(null, 0)==0, "null sites");
			require(Tools.countTopScores(new ArrayList<SiteScore>(), 0)==0, "empty sites");
			final ArrayList<SiteScore> single=new ArrayList<SiteScore>(1);
			single.add(new SiteScore(1, Shared.PLUS, 100, 109, 1, 100));
			require(Tools.countTopScores(single, 0)==1, "single best site");
			require(Tools.countTopScores(sites(1, Shared.PLUS, 100, 119, 100), 0)==1, "shared start heuristic");
			require(Tools.countTopScores(sites(1, Shared.PLUS, 90, 109, 100), 0)==1, "shared stop heuristic");
			require(Tools.countTopScores(sites(2, Shared.PLUS, 90, 109, 100), 0)==2, "shared stop on different chromosome");
			require(Tools.countTopScores(sites(1, Shared.PLUS, 100^Integer.MIN_VALUE,
					109^Integer.MIN_VALUE, 100), 0)==2, "coordinate differences with only sign bits set");
			require(Tools.countTopScores(sites(1, Shared.MINUS, 100, 109, 100), 0)==1, "existing strand heuristic");
			require(Tools.countTopScores(sites(1, Shared.PLUS, 200, 209, 95), 5)==2, "inclusive score cutoff");
			require(Tools.countTopScores(sites(2, Shared.PLUS, 100, 109, 94), 5)==1, "below cutoff");
		});
		test.check("delimiter at search start", () -> {
			require(find(";TYP=X", "TYP=", 0)==1, "leading delimiter");
			require(find("junk;TYP=X", "TYP=", 4)==5, "delimiter at nonzero search offset");
		});
		test.check("delimited search controls", () -> {
			require(find("TYP=X", "TYP=", 0)==0, "token at search start");
			require(find("XXTYP=X;TYP=Y", "TYP=", 0)==8, "skip undelimited substring");
			require(find("junk;TYP=X", "TYP=", 5)==5, "token at nonzero start");
			require(find("ABC=X;", "TYP=", 0)==-1, "missing token");
		});
		System.out.println("ToolsLocationBoundaryTest cases="+test.cases+" failures="+test.failures);
		if(test.failures!=0){throw new AssertionError("Location boundary failures="+test.failures);}
	}

	/** Creates descending-score sites with a fixed best site at chromosome1:100..109. */
	private static ArrayList<SiteScore> sites(int chrom, byte strand, int start, int stop, int score){
		final ArrayList<SiteScore> list=new ArrayList<SiteScore>(2);
		list.add(new SiteScore(1, Shared.PLUS, 100, 109, 1, 100));
		list.add(new SiteScore(chrom, strand, start, stop, 1, score));
		return list;
	}

	/** Searches an ASCII fixture with a semicolon separator. */
	private static int find(String text, String query, int start){
		return Tools.indexOfDelimited(text.getBytes(StandardCharsets.US_ASCII), query, start, (byte)';');
	}

	/** Reports each group before failing the whole run. */
	private void check(String name, Runnable action){
		cases++;
		try{action.run();System.out.println(name+" PASS");}
		catch(RuntimeException | AssertionError failure){
			failures++;
			System.out.println(name+" FAIL: "+failure);
		}
	}

	/** Keeps expectations active under both assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Attempted and failed groups. */
	private int cases, failures;
}
