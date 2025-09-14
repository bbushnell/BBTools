package aligner;

/**
 * Unit tests for clipping functionality in IndelFreeAligner.
 * Tests various scenarios including left/right clipping, clipping limits,
 * edge cases with tiny sequences, and non-overlapping alignments.
 * 
 * @author Isla SOS
 * @contributor Brian Bushnell
 * @date June 6, 2025
 */
public class ClippingTest {
    
    public static void main(String[] args) {
        // Set up basic parameters for testing
        Query.setMode(11, 1, true); // k=11, mm=1, indexing enabled
        
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
     * Tests left clipping functionality where query extends before reference start.
     * Verifies that clipped bases are handled correctly and excess clipping is penalized.
     * @return true if all left clipping tests pass, false otherwise
     */
    private static boolean testLeftClipping() {
        System.out.println("Testing left clipping...");
        
        // Reference: ATCGATCGATCG (12 bases)
        byte[] ref = "ATCGATCGATCG".getBytes();
        
        // Query extends 3 bases before reference start
        // Query: GGGATCGATCGATCG (should clip GGG, align ATCGATCGATCG)
        byte[] query = "GGGATCGATCGATCG".getBytes();
        
        // Test with maxClips=3 (should allow this - 3 clipped bases, 0 subs)
        int result = IndelFreeAligner.alignClipped(query, ref, 2, 3, -3);
        if(result != 0) {
            System.err.println("FAIL: Expected 0 subs with maxClips=3, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        // Test with maxClips=1 (should penalize 2 excess clipped bases as subs)
        result = IndelFreeAligner.alignClipped(query, ref, 5, 1, -3);
        if(result != 2) {
            System.err.println("FAIL: Expected 2 subs with maxClips=1, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        System.out.println("Left clipping tests passed");
        return true;
    }
    
    /**
     * Tests right clipping functionality where query extends past reference end.
     * Verifies that right-side clipped bases are handled correctly.
     * @return true if all right clipping tests pass, false otherwise
     */
    private static boolean testRightClipping() {
        System.out.println("Testing right clipping...");
        
        // Reference: ATCGATCGATCG (12 bases)
        byte[] ref = "ATCGATCGATCG".getBytes();
        
        // Query extends 3 bases past reference end
        // Query: ATCGATCGATCGGGG (should clip GGG, align ATCGATCGATCG)
        byte[] query = "ATCGATCGATCGGGG".getBytes();
        
        // Test with maxClips=3 (should allow this - 3 clipped bases, 0 subs)
        int result = IndelFreeAligner.alignClipped(query, ref, 2, 3, 0);
        if(result != 0) {
            System.err.println("FAIL: Expected 0 subs with maxClips=3, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        // Test with maxClips=1 (should penalize 2 excess clipped bases)
        result = IndelFreeAligner.alignClipped(query, ref, 5, 1, 0);
        if(result != 2) {
            System.err.println("FAIL: Expected 2 subs with maxClips=1, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        System.out.println("Right clipping tests passed");
        return true;
    }
    
    /**
     * Tests scenarios where query extends past both ends of the reference.
     * Verifies that total clipping from both sides is calculated correctly.
     * @return true if all both-sides clipping tests pass, false otherwise
     */
    private static boolean testBothSidesClipping() {
        System.out.println("Testing both sides clipping...");
        
        // Reference: ATCGATCG (8 bases)
        byte[] ref = "ATCGATCG".getBytes();
        
        // Query extends on both sides: GGATCGATCGTT
        byte[] query = "GGATCGATCGTT".getBytes();
        
        // Test with maxClips=4 (2 left + 2 right = 4 total clips)
        int result = IndelFreeAligner.alignClipped(query, ref, 2, 4, -2);
        if(result != 0) {
            System.err.println("FAIL: Expected 0 subs with maxClips=4, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        System.out.println("Both sides clipping tests passed");
        return true;
    }
    
    /**
     * Tests enforcement of clipping limits when total clips exceed maxClips.
     * Verifies that excess clipping is properly converted to substitution penalties.
     * @return true if all clipping limit tests pass, false otherwise
     */
    private static boolean testClippingLimits() {
        System.out.println("Testing clipping limits...");
        
        byte[] ref = "ATCG".getBytes();
        byte[] query = "GGGGGATCGTTTT".getBytes(); // 5 left + 4 right = 9 clips total
        
        // With maxClips=2, should get 7 excess clips counted as subs
        int result = IndelFreeAligner.alignClipped(query, ref, 10, 2, -5);
        if(result != 7) {
            System.err.println("FAIL: Expected 7 subs with maxClips=2, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        System.out.println("Clipping limits tests passed");
        return true;
    }
    
    /**
     * Tests alignments where no clipping is needed (query fits within reference).
     * Verifies that normal alignment works correctly when clipping is not required.
     * @return true if all no-clipping tests pass, false otherwise
     */
    private static boolean testNoClippingNeeded() {
        System.out.println("Testing no clipping needed...");
        
        byte[] ref = "ATCGATCGATCG".getBytes();
        byte[] query = "ATCGATCG".getBytes(); // Fits entirely within reference
        
        // Should have 0 clips and 0 subs for exact match
        int result = IndelFreeAligner.alignClipped(query, ref, 2, 0, 0);
        if(result != 0) {
            System.err.println("FAIL: Expected 0 subs for exact internal match, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        System.out.println("No clipping tests passed");
        return true;
    }
    
    /**
     * Tests exact matches between query and reference of the same length.
     * Basic sanity check for the alignment function.
     * @return true if exact match test passes, false otherwise
     */
    private static boolean testExactMatch() {
        System.out.println("Testing exact match...");
        
        byte[] ref = "ATCGATCG".getBytes();
        byte[] query = "ATCGATCG".getBytes();
        
        int result = IndelFreeAligner.alignClipped(query, ref, 2, 0, 0);
        if(result != 0) {
            System.err.println("FAIL: Expected 0 subs for exact match, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        System.out.println("Exact match tests passed");
        return true;
    }
    
    /**
     * Tests scenarios combining clipping with substitutions in the aligned region.
     * Verifies that both clipping penalties and substitution penalties are calculated correctly.
     * @return true if all mixed clipping/substitution tests pass, false otherwise
     */
    private static boolean testClippingWithSubstitutions() {
        System.out.println("Testing clipping with substitutions...");
        
        // Reference: ATCGATCG (8 bases)
        byte[] ref = "ATCGATCG".getBytes();
        
        // Query: GGATCGTTCG (2 left clips + 1 sub in aligned region)
        //        GG ATCG T TCG
        //           ATCG A TCG  (T->A substitution at position 4)
        byte[] query = "GGATCGATCG".getBytes();
        query[6] = 'T'; // Change A to T to create substitution
        
        // Test: 2 clips + 1 sub, maxClips=2, maxSubs=2 -> should be 1 total
        int result = IndelFreeAligner.alignClipped(query, ref, 2, 2, -2);
        if(result != 1) {
            System.err.println("FAIL: Expected 1 sub (2 clips allowed + 1 real sub), got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        // Test: Same alignment but maxClips=1 -> 1 excess clip + 1 sub = 2 total
        result = IndelFreeAligner.alignClipped(query, ref, 3, 1, -2);
        if(result != 2) {
            System.err.println("FAIL: Expected 2 subs (1 excess clip + 1 real sub), got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        System.out.println("Clipping with substitutions tests passed");
        return true;
    }
    
    /**
     * Tests edge cases with very small sequences that might cause boundary issues.
     * Includes 1bp references, 1bp queries, and extreme size mismatches.
     * @return true if all edge case tests pass, false otherwise
     */
    private static boolean testEdgeCases() {
        System.out.println("Testing edge cases...");
        
        // Test 1: 1bp reference
        byte[] tinyRef = "A".getBytes();
        byte[] normalQuery = "ATCG".getBytes();
        
        // Should clip 3 bases on the right
        int result = IndelFreeAligner.alignClipped(normalQuery, tinyRef, 5, 3, 0);
        if(result != 0) {
            System.err.println("FAIL: 1bp ref test - Expected 0 subs with maxClips=3, got " + result);
            System.err.println("Query: " + new String(normalQuery));
            System.err.println("Ref:   " + new String(tinyRef));
            return false;
        }
        
        // Test 2: 1bp query
        byte[] normalRef = "ATCGATCG".getBytes();
        byte[] tinyQuery = "A".getBytes();
        
        result = IndelFreeAligner.alignClipped(tinyQuery, normalRef, 2, 0, 0);
        if(result != 0) {
            System.err.println("FAIL: 1bp query test - Expected 0 subs for exact match, got " + result);
            System.err.println("Query: " + new String(tinyQuery));
            System.err.println("Ref:   " + new String(normalRef));
            return false;
        }
        
        // Test 3: Query longer than reference with mismatches
        byte[] shortRef = "AT".getBytes();
        byte[] longQuery = "GTCGAA".getBytes(); // G!=A, T=T, then 4 clips
        
        result = IndelFreeAligner.alignClipped(longQuery, shortRef, 5, 4, 0);
        if(result != 1) {
            System.err.println("FAIL: Long query test - Expected 1 sub (1 mismatch + 4 allowed clips), got " + result);
            System.err.println("Query: " + new String(longQuery));
            System.err.println("Ref:   " + new String(shortRef));
            return false;
        }
        
        // Test 4: Query much longer than reference, exceeding clip limits
        byte[] veryShortRef = "CG".getBytes();
        byte[] veryLongQuery = "AAACGTTTT".getBytes(); // 3 left + 4 right clips = 7 total
        
        result = IndelFreeAligner.alignClipped(veryLongQuery, veryShortRef, 10, 2, -3);
        if(result != 5) { // 7 clips - 2 allowed = 5 penalty subs
            System.err.println("FAIL: Very long query test - Expected 5 subs (7 clips - 2 allowed), got " + result);
            System.err.println("Query: " + new String(veryLongQuery));
            System.err.println("Ref:   " + new String(veryShortRef));
            return false;
        }
        
        System.out.println("Edge case tests passed");
        return true;
    }
    
    /**
     * Tests alignment scenarios where query and reference don't overlap at all.
     * Includes cases where rStart is far negative or extends far beyond reference end.
     * @return true if all non-overlapping tests pass, false otherwise
     */
    private static boolean testNonOverlapping() {
        System.out.println("Testing non-overlapping cases...");
        
        byte[] ref = "ATCG".getBytes(); // 4 bases
        byte[] query = "GGGG".getBytes(); // 4 bases
        
        // Test 1: Query starts way before reference (no overlap)
        // rStart = -10, query ends at position -6, ref starts at 0
        int result = IndelFreeAligner.alignClipped(query, ref, 10, 10, -10);
        if(result != 4) { // All 4 bases should be clipped
            System.err.println("FAIL: Far left non-overlap - Expected 4 clipped bases, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        // Test 2: Query starts way after reference ends (no overlap)
        // rStart = 10, ref ends at 3, query starts at 10
        result = IndelFreeAligner.alignClipped(query, ref, 10, 10, 10);
        if(result != 4) { // All 4 bases should be clipped
            System.err.println("FAIL: Far right non-overlap - Expected 4 clipped bases, got " + result);
            System.err.println("Query: " + new String(query));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        // Test 3: Partial overlap at left edge
        byte[] partialQuery = "GGGGATCG".getBytes(); // 8 bases
        result = IndelFreeAligner.alignClipped(partialQuery, ref, 5, 3, -4);
        // 2 left clips + 0 subs in overlapping AT + 2 right clips = 4 clips total
        // With maxClips=3, should get 1 penalty sub
        if(result != 1) {
            System.err.println("FAIL: Partial left overlap - Expected 1 sub (4 clips - 3 allowed), got " + result);
            System.err.println("Query: " + new String(partialQuery));
            System.err.println("Ref:   " + new String(ref));
            return false;
        }
        
        System.out.println("Non-overlapping tests passed");
        return true;
    }
}