package prok;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Focused regression test for the 2026-08-27 subtractClaimed fix (TrnaCaller, mirrored in
 * NcrnaScavenger): a claim landing STRICTLY INSIDE a window must yield BOTH surviving remainders
 * (left and right), not just the right one. The prior version's fallthrough "else{lo=c[1]+1;}"
 * silently discarded the left remainder's real candidate territory -- found via a direct runtime
 * trace of a tRNA shortlist tie-break that flipped which of two adjacent tied-identity candidates
 * got claimed first, orphaning the other's true region on the discarded side of this bug.
 *
 * Not a JUnit test (none of this package's other classes use one) -- a self-contained main() with
 * loud asserts, run via `java -ea`. Java 8 safe throughout (no String.repeat() or other Java 11+ API).
 *
 * MIN_TRNA=40 (TrnaCaller's own minimum), so surviving segments below that (by the EXISTING,
 * unfixed hi-lo>=MIN_TRNA convention -- not hi-lo+1) are filtered; test intervals are sized well
 * clear of that boundary so the filter never ambiguously eats a segment the test cares about.
 *
 * @author G11
 */
public class TrnaCallerSubtractClaimedTest {

	public static void main(String[] args){
		testNoOverlap();
		testFullCover();
		testLeftClip();
		testRightClip();
		testInteriorSplit();
		testMultipleClaims();
		testOrderingIndependence();
		testMultiInteriorClaimsThreeFragments();
		System.out.println("TrnaCallerSubtractClaimedTest: ALL TESTS PASSED");
	}

	static ArrayList<int[]> windows(int[]... arr){
		ArrayList<int[]> l=new ArrayList<>();
		for(int[] a : arr){l.add(a);}
		return l;
	}

	static String show(ArrayList<int[]> segs){
		StringBuilder sb=new StringBuilder("[");
		for(int i=0; i<segs.size(); i++){
			if(i>0){sb.append(",");}
			sb.append(segs.get(i)[0]).append("-").append(segs.get(i)[1]);
		}
		return sb.append("]").toString();
	}

	/** Asserts result contains exactly the expected [lo,hi] pairs, order-independent (sorted by lo). */
	static void assertSegments(String testName, ArrayList<int[]> result, int[]... expected){
		ArrayList<int[]> sortedResult=new ArrayList<>(result);
		sortedResult.sort((a, b) -> a[0]-b[0]);
		ArrayList<int[]> sortedExpected=new ArrayList<>(Arrays.asList(expected));
		sortedExpected.sort((a, b) -> a[0]-b[0]);
		assert(sortedResult.size()==sortedExpected.size()) : testName+": expected "+sortedExpected.size()
			+" segments "+show(sortedExpected)+", got "+sortedResult.size()+" "+show(sortedResult);
		for(int i=0; i<sortedResult.size(); i++){
			assert(sortedResult.get(i)[0]==sortedExpected.get(i)[0] && sortedResult.get(i)[1]==sortedExpected.get(i)[1])
				: testName+": segment "+i+" expected "+sortedExpected.get(i)[0]+"-"+sortedExpected.get(i)[1]
				+", got "+sortedResult.get(i)[0]+"-"+sortedResult.get(i)[1];
		}
		System.out.println(testName+" PASSED: "+show(sortedResult));
	}

	static void testNoOverlap(){
		ArrayList<int[]> result=TrnaCaller.subtractClaimed(windows(new int[]{10, 100}), windows(new int[]{200, 300}));
		assertSegments("testNoOverlap", result, new int[]{10, 100});
	}

	static void testFullCover(){
		ArrayList<int[]> result=TrnaCaller.subtractClaimed(windows(new int[]{10, 100}), windows(new int[]{0, 200}));
		assertSegments("testFullCover", result);//nothing survives
	}

	static void testLeftClip(){
		//claim overlaps the window's LEFT edge -- only the right remainder can survive.
		ArrayList<int[]> result=TrnaCaller.subtractClaimed(windows(new int[]{10, 100}), windows(new int[]{0, 40}));
		assertSegments("testLeftClip", result, new int[]{41, 100});
	}

	static void testRightClip(){
		//claim overlaps the window's RIGHT edge -- only the left remainder can survive.
		ArrayList<int[]> result=TrnaCaller.subtractClaimed(windows(new int[]{10, 100}), windows(new int[]{60, 200}));
		assertSegments("testRightClip", result, new int[]{10, 59});
	}

	/** THE core regression case: a claim strictly inside the window must yield BOTH remainders. */
	static void testInteriorSplit(){
		ArrayList<int[]> result=TrnaCaller.subtractClaimed(windows(new int[]{10, 200}), windows(new int[]{80, 120}));
		assertSegments("testInteriorSplit", result, new int[]{10, 79}, new int[]{121, 200});
	}

	static void testMultipleClaims(){
		//Two non-adjacent interior claims against one window -> three surviving fragments. Window starts
		//at 0 (not 10) so the first remainder [0,49] clears MIN_TRNA under the EXISTING hi-lo>=MIN_TRNA
		//convention (49-0=49>=40); a window starting at 10 would make that remainder hi-lo=39, one short
		//of 40 by the pre-existing (deliberately unfixed) off-by-one -- a test-sizing pitfall, not a
		//split-logic bug, caught by this test's own first run.
		ArrayList<int[]> result=TrnaCaller.subtractClaimed(windows(new int[]{0, 400}),
			windows(new int[]{50, 90}, new int[]{200, 250}));
		assertSegments("testMultipleClaims", result, new int[]{0, 49}, new int[]{91, 199}, new int[]{251, 400});
	}

	/** Same as testMultipleClaims but claims given in REVERSE order -- result must be identical
	 * (subtraction is order-independent: each claim is subtracted from every currently-surviving
	 * segment regardless of which claim came first). */
	static void testOrderingIndependence(){
		ArrayList<int[]> forward=TrnaCaller.subtractClaimed(windows(new int[]{0, 400}),
			windows(new int[]{50, 90}, new int[]{200, 250}));
		ArrayList<int[]> reversed=TrnaCaller.subtractClaimed(windows(new int[]{0, 400}),
			windows(new int[]{200, 250}, new int[]{50, 90}));
		assertSegments("testOrderingIndependence(reversed)", reversed, new int[]{0, 49}, new int[]{91, 199}, new int[]{251, 400});
		assert(forward.size()==reversed.size()) : "testOrderingIndependence: forward/reversed size mismatch";
		System.out.println("testOrderingIndependence PASSED: forward="+show(forward)+" reversed="+show(reversed));
	}

	/** Explicit >=3-fragment case from a single window carved by 3 separate interior claims, per
	 * Citan's requirement -- verifies the segment list correctly carries through MORE than one split
	 * (each subsequent claim must apply to ALL currently-surviving segments, not just the original pair). */
	static void testMultiInteriorClaimsThreeFragments(){
		ArrayList<int[]> result=TrnaCaller.subtractClaimed(windows(new int[]{0, 600}),
			windows(new int[]{100, 150}, new int[]{300, 350}, new int[]{500, 550}));
		assertSegments("testMultiInteriorClaimsThreeFragments", result,
			new int[]{0, 99}, new int[]{151, 299}, new int[]{351, 499}, new int[]{551, 600});
	}
}
