package prok;

/**
 * Focused regression test for the 2026-08-27 unique-query-k-mer change to TrnaKmerIndex.shortlist()
 * (numerator dedup) and its qKmers denominator fix (2026-08-27, same day, second pass).
 *
 * Demonstrates two bounds the fix is supposed to restore:
 * (1) a model's shared count for one query must not exceed that model's OWN total index-k-mer
 *     occurrence count, even when the query repeats a k-mer many times (numerator dedup);
 * (2) the adaptive qFrac cutoff's denominator (qKmers) must be the query's UNIQUE k-mer count, not
 *     the raw positional count -- a repetitive query must not inflate qKmers and silently raise
 *     minHits past what a correct model match would clear.
 *
 * Not a JUnit test (none of this package's other classes use one) -- a self-contained main() with loud
 * asserts, run via `java -ea`, matching this codebase's assertion-first testing convention. Java 8
 * safe throughout (BBTools ships to Java 8 -- no String.repeat(), a Java 11+ method).
 *
 * @author G11
 */
public class TrnaKmerIndexBoundTest {

	public static void main(String[] args){
		testRepetitiveQueryBound();
		testUniqueKmerCollapse();
		testAdaptiveDenominatorIsUniqueNotPositional();
		System.out.println("TrnaKmerIndexBoundTest: ALL TESTS PASSED");
	}

	/** Model 0's sequence has k-mer "AAAAAAA" (k=7, all-A) at 5 distinct positions (it repeats "A" runs
	 * separated by "C" to keep the 7-mer identical each time without making the whole sequence one giant
	 * repeat that the assembler might collapse conceptually -- doesn't matter here, this is pure counting).
	 * Model 0's TOTAL k-mer occurrence count is its sequence length-k+1. A query built from many repeats
	 * of the SAME "AAAAAAA" 7-mer must not push counts[0] past model 0's actual occurrence count for that
	 * k-mer (5), regardless of how many times the query repeats it. */
	static void testRepetitiveQueryBound(){
		final int k=7;
		//Model 0: "AAAAAAAC" x5 -- gives exactly 5 occurrences of the 7-mer "AAAAAAA" (one per "AAAAAAAC"
		//block; the trailing C breaks the run so there isn't a 6th overlapping AAAAAAA at the boundary).
		final byte[] model0=toBytes(repeat("AAAAAAAC", 5));
		final byte[][] library={model0};

		//adaptive=false, fixedMinHits=1 -- deterministic minHits, not exercising the adaptive formula here.
		final TrnaKmerIndex idx=new TrnaKmerIndex(library, k, false, 0f, 0f, 0f, 1);

		//Query: the SAME 7-mer "AAAAAAA" repeated 200 times back-to-back (a single long run of 'A').
		//Old (pre-fix) code would traverse postings["AAAAAAA"] once per raw query position (~194 valid
		//positions in a 200-mer run), each traversal adding model 0's 5 occurrences -> counts[0] ~970.
		//Fixed code collapses to ONE unique query k-mer ("AAAAAAA"), so counts[0] must be exactly 5.
		final byte[] query=toBytes(repeat("A", 200));
		final int[] result=idx.shortlist(query, 10);
		final int count0=idx.lastSharedCount(0);

		assert(count0==5) : "Expected counts[0]==5 (model 0's own AAAAAAA occurrence count), got "+count0
			+" -- unique-query-k-mer dedup is not bounding the count as intended.";
		//Model 0's total index-k-mer occurrence count = sequence length - k + 1 = 40-7+1 = 34.
		final int model0TotalKmerOccurrences=model0.length-k+1;
		assert(count0<=model0TotalKmerOccurrences) : "counts[0]="+count0
			+" exceeds model 0's total k-mer occurrence count "+model0TotalKmerOccurrences+" -- bound violated.";
		assert(result.length==1 && result[0]==0) : "Expected model 0 shortlisted alone, got "
			+java.util.Arrays.toString(result);
		System.out.println("testRepetitiveQueryBound PASSED: counts[0]="+count0
			+" (model total k-mer occurrences="+model0TotalKmerOccurrences+")");
	}

	/** A query with NO repeated k-mers (every 7-mer window distinct) should give the SAME result under
	 * old and new logic -- the dedup pass only changes behavior when the query actually repeats a k-mer.
	 * Sanity check that the common case (non-repetitive real sequence) is unaffected. */
	static void testUniqueKmerCollapse(){
		final int k=7;
		//A sequence with no internal repeats at k=7 (simple non-periodic base pattern).
		final String seqStr="ACGTACGGTACGTTACGGATTCGGACTGATCGATGCATGCTAGCTAGCTAGGCATCGTAGCTAGT";
		final byte[] model0=toBytes(seqStr);
		final byte[][] library={model0};
		final TrnaKmerIndex idx=new TrnaKmerIndex(library, k, false, 0f, 0f, 0f, 1);

		final byte[] query=toBytes(seqStr);//query identical to the model -> every k-mer shared exactly once
		idx.shortlist(query, 10);
		final int count0=idx.lastSharedCount(0);
		final int expected=seqStr.length()-k+1;//every window is a distinct k-mer here, one query occurrence each
		assert(count0==expected) : "Expected counts[0]=="+expected+" for a self-match with no repeated "
			+"k-mers, got "+count0+" -- dedup pass should not alter non-repetitive-query behavior.";
		System.out.println("testUniqueKmerCollapse PASSED: counts[0]="+count0+" (expected "+expected+")");
	}

	/** Distinguishes qKmers=unique-count (correct, current) from qKmers=positional-count (the
	 * superseded intermediate state) via the ADAPTIVE cutoff formula, with >=2 matching models so a
	 * keep-top-1 fallback (triggered when zero models organically clear minHits) cannot disguise a
	 * regression as "still returns something plausible" -- a regression back to the positional
	 * denominator changes the RESULT SIZE (2 -> 1), not just an internal count.
	 *
	 * Setup: two models sharing the query's one repeated k-mer type, with occurrence counts 20 and 8.
	 * Query is 500x the same base -> 494 raw k-mer positions but exactly 1 UNIQUE k-mer type.
	 * floor=0, topFrac=0, qFrac=0.5 isolates the qFrac*qKmers term as the only source of minHits.
	 *   - Correct (qKmers=unique=1): minHits=ceil(0.5*1)=1. Both models (20>=1, 8>=1) clear ->
	 *     shortlist=[model0,model1], size 2 (count-desc order: 20 before 8).
	 *   - Regressed (qKmers=positional=494): minHits=ceil(0.5*494)=247. NEITHER model (20, 8) clears
	 *     that organically -> the keep-top-1 fallback forces size 1 ([model0] only). This is the
	 *     "fallback masks a real threshold difference" case Citan's review named explicitly -- this
	 *     test catches it via result LENGTH, which the fallback cannot make look identical. */
	static void testAdaptiveDenominatorIsUniqueNotPositional(){
		final int k=7;
		final byte[] model0=toBytes(repeat("AAAAAAAC", 20));//20 occurrences of "AAAAAAA"
		final byte[] model1=toBytes(repeat("AAAAAAAC", 8));//8 occurrences of "AAAAAAA"
		final byte[][] library={model0, model1};

		final TrnaKmerIndex idx=new TrnaKmerIndex(library, k, true, 0f, 0f, 0.5f, 1);//adaptive, qFrac=0.5

		final byte[] query=toBytes(repeat("A", 500));//494 raw positions, 1 unique k-mer type
		final int[] result=idx.shortlist(query, 10);

		assert(result.length==2) : "Expected BOTH models shortlisted (qKmers=unique=1 -> minHits=1, "
			+"both 20 and 8 clear), got "+result.length+" model(s): "+java.util.Arrays.toString(result)
			+" -- if this is 1, qKmers has regressed to the positional count (minHits would balloon to "
			+"ceil(0.5*494)=247, neither model clears organically, and the keep-top-1 fallback silently "
			+"returns just the best one instead of both).";
		assert(result[0]==0 && result[1]==1) : "Expected count-desc order [0,1] (20 before 8), got "
			+java.util.Arrays.toString(result);
		System.out.println("testAdaptiveDenominatorIsUniqueNotPositional PASSED: shortlist="
			+java.util.Arrays.toString(result)+" (both models cleared minHits=1, confirming qKmers="
			+"unique-count=1, not positional-count=494)");
	}

	static byte[] toBytes(String s){return s.getBytes(java.nio.charset.StandardCharsets.US_ASCII);}

	/** Java-8-safe replacement for String.repeat() (Java 11+) -- BBTools compiles on Java 8. */
	static String repeat(String s, int n){
		final StringBuilder sb=new StringBuilder(s.length()*n);
		for(int i=0; i<n; i++){sb.append(s);}
		return sb.toString();
	}
}
