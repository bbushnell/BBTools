package ifa;

/**
 * Assertion-driven fixtures for the original {@link IndelFreeAligner#alignClipped}
 * pointwise scorer. Covers left/right overhangs, clipping penalties, substitutions,
 * small sequences and nonoverlap scores. It does not test candidate enumeration,
 * versions 2 through 4, SIMD dispatch or the alignment pipeline.
 * <p>
 * Run with assertions enabled ({@code -ea}): main invokes every test group only
 * inside an assert expression. With assertions disabled the calls are skipped,
 * but the final success message still prints. A failed group prints its fixture
 * diagnostics and returns false, causing the enclosing assertion to fail.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date June 6, 2025
 */
public class ClippingTest{

	/** Runs the nine fixture groups through assertions; requires -ea and ignores args. */
	public static void main(String[] args){
		// Set up basic parameters for testing
//		Query.setMode(11, 1, true); // k=11, mm=1, indexing enabled

		assert(testLeftClipping()) : "Left clipping test failed";
		assert(testRightClipping()) : "Right clipping test failed";
		assert(testBothSidesClipping()) : "Both sides clipping test failed";
		assert(testClippingLimits()) : "Clipping limits test failed";
		assert(testNoClippingNeeded()) : "No clipping test failed";
		assert(testExactMatch()) : "Exact match test failed";
		assert(testClippingWithSubstitutions()) : "Clipping with substitutions test failed";
		assert(testEdgeCases()) : "Edge cases test failed";
		assert(testNonOverlapping()) : "Non-overlapping test failed";

		System.out.println("All clipping tests passed!");
	}

	/**
	 * Tests left clipping where the query extends before reference start.
	 * Checks free overhangs and substitution cost for clipping beyond the allowance.
	 * @return true if all left clipping tests pass, false otherwise
	 */
	private static boolean testLeftClipping(){
		System.out.println("Testing left clipping...");

		// Reference: ATCGATCGATCG (12 bases)
		byte[] ref="ATCGATCGATCG".getBytes();

		// Three left clips (GGG), then an exact match to the reference.
		byte[] query="GGGATCGATCGATCG".getBytes();

		// Three free clips, zero substitutions.
		int result=IndelFreeAligner.alignClipped(query, ref, 2, 3, -3);
		if(result!=0){
			System.err.println("FAIL: Expected 0 subs with maxClips=3, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		// One free clip leaves two clips charged as substitutions.
		result=IndelFreeAligner.alignClipped(query, ref, 5, 1, -3);
		if(result!=2){
			System.err.println("FAIL: Expected 2 subs with maxClips=1, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		System.out.println("Left clipping tests passed");
		return true;
	}

	/**
	 * Tests right clipping where the query extends past reference end.
	 * @return true if all right clipping tests pass, false otherwise
	 */
	private static boolean testRightClipping(){
		System.out.println("Testing right clipping...");

		// Reference: ATCGATCGATCG (12 bases)
		byte[] ref="ATCGATCGATCG".getBytes();

		// Exact reference match followed by three right clips (GGG).
		byte[] query="ATCGATCGATCGGGG".getBytes();

		// Three free clips, zero substitutions.
		int result=IndelFreeAligner.alignClipped(query, ref, 2, 3, 0);
		if(result!=0){
			System.err.println("FAIL: Expected 0 subs with maxClips=3, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		// One free clip leaves two clips charged as substitutions.
		result=IndelFreeAligner.alignClipped(query, ref, 5, 1, 0);
		if(result!=2){
			System.err.println("FAIL: Expected 2 subs with maxClips=1, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		System.out.println("Right clipping tests passed");
		return true;
	}

	/**
	 * Checks that clipping at both ends contributes to one total allowance.
	 * @return true if all both-sides clipping tests pass, false otherwise
	 */
	private static boolean testBothSidesClipping(){
		System.out.println("Testing both sides clipping...");

		byte[] ref="ATCGATCG".getBytes();
		byte[] query="GGATCGATCGTT".getBytes();

		// Two left plus two right clips, all four free.
		int result=IndelFreeAligner.alignClipped(query, ref, 2, 4, -2);
		if(result!=0){
			System.err.println("FAIL: Expected 0 subs with maxClips=4, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		System.out.println("Both sides clipping tests passed");
		return true;
	}

	/**
	 * Checks substitution cost when total clipping exceeds maxClips.
	 * @return true if all clipping limit tests pass, false otherwise
	 */
	private static boolean testClippingLimits(){
		System.out.println("Testing clipping limits...");

		byte[] ref="ATCG".getBytes();
		byte[] query="GGGGGATCGTTTT".getBytes();//Five left plus four right clips.

		// Nine clips minus two free clips gives cost seven.
		int result=IndelFreeAligner.alignClipped(query, ref, 10, 2, -5);
		if(result!=7){
			System.err.println("FAIL: Expected 7 subs with maxClips=2, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		System.out.println("Clipping limits tests passed");
		return true;
	}

	/**
	 * Checks an exact query placement wholly inside a longer reference.
	 * @return true if all no-clipping tests pass, false otherwise
	 */
	private static boolean testNoClippingNeeded(){
		System.out.println("Testing no clipping needed...");

		byte[] ref="ATCGATCGATCG".getBytes();
		byte[] query="ATCGATCG".getBytes();

		int result=IndelFreeAligner.alignClipped(query, ref, 2, 0, 0);
		if(result!=0){
			System.err.println("FAIL: Expected 0 subs for exact internal match, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		System.out.println("No clipping tests passed");
		return true;
	}

	/**
	 * Checks an exact match between equal-length query and reference.
	 * @return true if exact match test passes, false otherwise
	 */
	private static boolean testExactMatch(){
		System.out.println("Testing exact match...");

		byte[] ref="ATCGATCG".getBytes();
		byte[] query="ATCGATCG".getBytes();

		int result=IndelFreeAligner.alignClipped(query, ref, 2, 0, 0);
		if(result!=0){
			System.err.println("FAIL: Expected 0 subs for exact match, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		System.out.println("Exact match tests passed");
		return true;
	}

	/**
	 * Checks substitution cost together with free and excess clipping.
	 * @return true if all mixed clipping/substitution tests pass, false otherwise
	 */
	private static boolean testClippingWithSubstitutions(){
		System.out.println("Testing clipping with substitutions...");

		byte[] ref="ATCGATCG".getBytes();
		byte[] query="GGATCGATCG".getBytes();
		query[6]='T';//A to T at reference position four, after two left clips.

		// Two free clips plus one substitution gives cost one.
		int result=IndelFreeAligner.alignClipped(query, ref, 2, 2, -2);
		if(result!=1){
			System.err.println("FAIL: Expected 1 sub (2 clips allowed + 1 real sub), got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		// One free clip leaves one excess clip plus one substitution.
		result=IndelFreeAligner.alignClipped(query, ref, 3, 1, -2);
		if(result!=2){
			System.err.println("FAIL: Expected 2 subs (1 excess clip + 1 real sub), got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		System.out.println("Clipping with substitutions tests passed");
		return true;
	}

	/**
	 * Checks one-base inputs and unequal lengths with free and excess clipping.
	 * @return true if all edge case tests pass, false otherwise
	 */
	private static boolean testEdgeCases(){
		System.out.println("Testing edge cases...");

		byte[] tinyRef="A".getBytes();
		byte[] normalQuery="ATCG".getBytes();

		// One matching base and three free right clips.
		int result=IndelFreeAligner.alignClipped(normalQuery, tinyRef, 5, 3, 0);
		if(result!=0){
			System.err.println("FAIL: 1bp ref test - Expected 0 subs with maxClips=3, got "+result);
			System.err.println("Query: "+new String(normalQuery));
			System.err.println("Ref:   "+new String(tinyRef));
			return false;
		}

		byte[] normalRef="ATCGATCG".getBytes();
		byte[] tinyQuery="A".getBytes();

		result=IndelFreeAligner.alignClipped(tinyQuery, normalRef, 2, 0, 0);
		if(result!=0){
			System.err.println("FAIL: 1bp query test - Expected 0 subs for exact match, got "+result);
			System.err.println("Query: "+new String(tinyQuery));
			System.err.println("Ref:   "+new String(normalRef));
			return false;
		}

		byte[] shortRef="AT".getBytes();
		byte[] longQuery="GTCGAA".getBytes();//G!=A, T=T, then four clips.

		result=IndelFreeAligner.alignClipped(longQuery, shortRef, 5, 4, 0);
		if(result!=1){
			System.err.println("FAIL: Long query test - Expected 1 sub (1 mismatch + 4 allowed clips), got "+result);
			System.err.println("Query: "+new String(longQuery));
			System.err.println("Ref:   "+new String(shortRef));
			return false;
		}

		byte[] veryShortRef="CG".getBytes();
		byte[] veryLongQuery="AAACGTTTT".getBytes();//Three left plus four right clips.

		result=IndelFreeAligner.alignClipped(veryLongQuery, veryShortRef, 10, 2, -3);
		if(result!=5){//Seven clips minus two free clips gives cost five.
			System.err.println("FAIL: Very long query test - Expected 5 subs (7 clips - 2 allowed), got "+result);
			System.err.println("Query: "+new String(veryLongQuery));
			System.err.println("Ref:   "+new String(veryShortRef));
			return false;
		}

		System.out.println("Edge case tests passed");
		return true;
	}

	/**
	 * Checks nonoverlap's query-length score and one partial-overlap score.
	 * These are pointwise scores, not claims that candidate enumeration emits them.
	 * @return true if all non-overlapping tests pass, false otherwise
	 */
	private static boolean testNonOverlapping(){
		System.out.println("Testing non-overlapping cases...");

		byte[] ref="ATCG".getBytes();
		byte[] query="GGGG".getBytes();

		// No overlap: start=-10, exclusive stop=-6. Scorer returns query.length.
		int result=IndelFreeAligner.alignClipped(query, ref, 10, 10, -10);
		if(result!=4){
			System.err.println("FAIL: Far left non-overlap - Expected 4 clipped bases, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		// No overlap: start=10 is past the reference's exclusive stop at four.
		result=IndelFreeAligner.alignClipped(query, ref, 10, 10, 10);
		if(result!=4){
			System.err.println("FAIL: Far right non-overlap - Expected 4 clipped bases, got "+result);
			System.err.println("Query: "+new String(query));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		byte[] partialQuery="GGGGATCG".getBytes();
		result=IndelFreeAligner.alignClipped(partialQuery, ref, 5, 3, -4);
		// Four left clips, four matching in-bounds bases and no right clip.
		// Three free clips leave one clip charged as a substitution.
		if(result!=1){
			System.err.println("FAIL: Partial left overlap - Expected 1 sub (4 clips - 3 allowed), got "+result);
			System.err.println("Query: "+new String(partialQuery));
			System.err.println("Ref:   "+new String(ref));
			return false;
		}

		System.out.println("Non-overlapping tests passed");
		return true;
	}
}
