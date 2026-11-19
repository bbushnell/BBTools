package test.prok;

import prok.*;

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
 * @author G11, Raiden
 */
public class TrnaKmerIndexBoundTest {

	public static void main(String[] args){
		testRepetitiveQueryBound();
		testUniqueKmerCollapse();
		testSingleModelStrictCutoff();
		testSingleModelDirectAlignmentControl();
		testAdaptiveDenominatorIsUniqueNotPositional();
		testCountOnlyLiteralAndReplacement();
		testCountOnlySingleModelControl();
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
		//The one-model set prefilter collapses repeats on both sides, so the shared count is exactly 1.
		final byte[] query=toBytes(repeat("A", 200));
		final int[] result=idx.shortlist(query, 10);
		final int count0=idx.lastSharedCount(0);

		assert(count0==1) : "Expected counts[0]==1 distinct shared k-mer type, got "+count0
			+" -- the one-model path must use set-membership rather than posting multiplicity.";
		//Model 0's total index-k-mer occurrence count = sequence length - k + 1 = 40-7+1 = 34.
		final int model0TotalKmerOccurrences=model0.length-k+1;
		assert(count0<=model0TotalKmerOccurrences) : "counts[0]="+count0
			+" exceeds model 0's total k-mer occurrence count "+model0TotalKmerOccurrences+" -- bound violated.";
		assert(result.length==1 && result[0]==0) : "Expected model 0 shortlisted alone, got "
			+java.util.Arrays.toString(result);
		System.out.println("testRepetitiveQueryBound PASSED: counts[0]="+count0
			+" (model total k-mer occurrences="+model0TotalKmerOccurrences+")");
	}

	/** A self-match must report the number of distinct model/query k-mer types. */
	static void testUniqueKmerCollapse(){
		final int k=7;
		//This sequence has 59 positional 7-mers and 55 distinct 7-mer types.
		final String seqStr="ACGTACGGTACGTTACGGATTCGGACTGATCGATGCATGCTAGCTAGCTAGGCATCGTAGCTAGT";
		final byte[] model0=toBytes(seqStr);
		final byte[][] library={model0};
		final TrnaKmerIndex idx=new TrnaKmerIndex(library, k, false, 0f, 0f, 0f, 1);

		final byte[] query=toBytes(seqStr);//query identical to the model -> every k-mer shared exactly once
		idx.shortlist(query, 10);
		final int count0=idx.lastSharedCount(0);
		final int expected=55;
		assert(count0==expected) : "Expected counts[0]=="+expected+" distinct shared k-mer types for the "
			+"self-match, got "+count0+".";
		System.out.println("testUniqueKmerCollapse PASSED: counts[0]="+count0+" (expected "+expected+")");
	}

	/** A one-model family must be rejectable by its configured count; the multi-model keep-best fallback
	 * must not force the sole model through after it misses the cutoff. */
	static void testSingleModelStrictCutoff(){
		final int k=7;
		final byte[] model=toBytes("ACGTACGGTACGTTACGGATTCGGACTGATCGATGCATGCTAGCTAGCTAGGCATCGTAGCTAGT");
		final TrnaKmerIndex pass=new TrnaKmerIndex(new byte[][]{model}, k, false, 0f, 0f, 0f, 1);
		final TrnaKmerIndex reject=new TrnaKmerIndex(new byte[][]{model}, k, false, 0f, 0f, 0f, model.length);
		final byte[] query=toBytes("ACGTACGG");
		assert(pass.shortlist(query, 100).length==1) : "One shared k-mer should pass fixedMinHits=1";
		assert(reject.shortlist(query, 100).length==0) : "One-model fixedMinHits must reject below-cutoff windows; "
			+"the keep-best fallback is only valid when ranking multiple models";
		System.out.println("testSingleModelStrictCutoff PASSED");
	}

	/** fixedMinHits=0 is the explicit no-prefilter control for one-model calibration. */
	static void testSingleModelDirectAlignmentControl(){
		final int k=7;
		final byte[] model=toBytes("AAAAAAAAAAAAAA");
		final TrnaKmerIndex direct=new TrnaKmerIndex(new byte[][]{model}, k, false, 0f, 0f, 0f, 0);
		final byte[] noSharedKmer=toBytes("CCCCCCCCCCCCCC");
		assert(direct.shortlist(noSharedKmer, 1).length==1) : "fixedMinHits=0 must align the sole model even "
			+"when the candidate window shares no index k-mer; this is the direct-alignment calibration control";
		System.out.println("testSingleModelDirectAlignmentControl PASSED");
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

	/** A unique query key contributes all model occurrences, irrespective of the
	 * shortlist floor/cap; empty replacement must clear every previous output. */
	static void testCountOnlyLiteralAndReplacement(){
		final byte[][] library={toBytes("AAAAAAAAA"), toBytes("CCCCC"), toBytes("AAANAAA")};
		final TrnaKmerIndex index=new TrnaKmerIndex(library, 3, false, 0f, 0f, 0f, 999);
		final char[] counts=new char[3];
		final byte[] query=toBytes(repeat("A", 200)+'N'+repeat("C", 200));
		index.countShared(query, counts);
		assert(counts[0]==7 && counts[1]==3 && counts[2]==2)
			: "Count-only must retain literal model occurrences7/3/2 once per unique query key, not apply floor999";
		assert(index.lastMaxShared()==7 && index.lastUniqueQueryKmers()==2 && index.queriesProcessed()==1)
			: "Count-only must refresh actual query diagnostics and record exactly one query";
		final int[] kept=index.shortlist(query, 1);
		assert(kept.length==1 && kept[0]==0 && index.lastSharedCount(1)==3 && index.lastSharedCount(2)==2)
			: "Legacy cap/floor/keep-best behavior and complete count diagnostics must survive count-only use";
		index.countShared(toBytes(""), counts);
		assert(counts[0]==0 && counts[1]==0 && counts[2]==0 && index.lastMaxShared()==0)
			: "An empty replacement must clear both the last native counts and caller-owned destination";
		index.countShared(toBytes("AAAAAAAAAAAAAAAA"), counts);
		assert(counts[0]==7 && counts[1]==0 && counts[2]==2)
			: "A second nonempty query cannot retain an absent model's preceding counts";
		index.setPolicy(true, 100f, 1f, 1f, 1);
		index.countShared(query, counts);
		assert(counts[0]==7 && counts[1]==3 && counts[2]==2)
			: "Adaptive shortlist settings cannot alter the count-only numerator";
		try{
			index.countShared(query, new char[2]);
			throw new AssertionError("A truncated native-count destination was accepted");
		}catch(IllegalArgumentException expected){
			assert(expected.getMessage().contains("one output slot")) : "Reject the destination shape before changing count scratch";
		}
		System.out.println("testCountOnlyLiteralAndReplacement PASSED");
	}

	/** One-model counts are strict set matches; fixed0 remains a separate legacy
	 * no-prefilter control, before and after the new API replaces its scratch. */
	static void testCountOnlySingleModelControl(){
		final TrnaKmerIndex index=new TrnaKmerIndex(new byte[][]{toBytes("AAAAAAAAA")}, 3, false, 0f, 0f, 0f, 0);
		final char[] counts=new char[1];
		final byte[] absent=toBytes("CCCCCCCCCCCC");
		assert(index.shortlist(absent, 1).length==1) : "fixed0 must retain the explicit direct-alignment control";
		index.countShared(toBytes("AAAAAAAAAAAAAAAA"), counts);
		assert(counts[0]==1 && index.lastMaxShared()==1)
			: "One-model count-only uses one distinct shared type, not the multi-model occurrence weight or fixed0 bypass";
		index.countShared(absent, counts);
		assert(counts[0]==0 && index.lastMaxShared()==0) : "No one-model match must produce actual zero support";
		assert(index.shortlist(absent, 1).length==1 && index.lastSharedCount(0)==0)
			: "Count-only must not change the following fixed0 control's policy or diagnostic zero";
		index.setFixedMinHits(99);
		index.countShared(toBytes("AAAAAAAAAAAAAAAA"), counts);
		assert(counts[0]==1) : "Count-only must not filter the one-model positive count below fixed99";
		assert(index.shortlist(toBytes("AAAAAAAAAAAAAAAA"), 1).length==0)
			: "The unchanged one-model shortlister must still reject that count below fixed99";
		System.out.println("testCountOnlySingleModelControl PASSED");
	}

	static byte[] toBytes(String s){return s.getBytes(java.nio.charset.StandardCharsets.US_ASCII);}

	/** Java-8-safe replacement for String.repeat() (Java 11+) -- BBTools compiles on Java 8. */
	static String repeat(String s, int n){
		final StringBuilder sb=new StringBuilder(s.length()*n);
		for(int i=0; i<n; i++){sb.append(s);}
		return sb.toString();
	}
}
