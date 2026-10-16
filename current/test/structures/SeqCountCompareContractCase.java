package test.structures;

import java.util.HashMap;

import sort.ReadComparatorCrispr;
import stream.Read;
import structures.SeqCount;
import structures.SeqCountM;

/** Literal compareTo contract oracle for SeqCount and SeqCountM.
 * @author Jean
 */
public final class SeqCountCompareContractCase{

	/** Runs current STU046 SeqCount/SeqCountM comparison expectations. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline|extended-baseline|fixed");
		check(args[0].startsWith("expect="),
				"Expected expect=baseline|extended-baseline|fixed, got "+args[0]);
		final String expect=args[0].substring("expect=".length());
		check("baseline".equals(expect) || "extended-baseline".equals(expect) ||
				"fixed".equals(expect),
				"Expected expect=baseline|extended-baseline|fixed, got "+args[0]);

		seqCountControls(expect);
		seqCountMControls(expect);
		mixedZeroConsistencyFailure(expect);
		inheritedReceiverFailure(expect);
		if(extended(expect)){
			countOverflowControl(expect);
			nonfiniteScoreControls(expect);
			hashMapControls(expect);
			readComparatorCrisprControls(expect);
		}
		if(fixed(expect)){
			orderingLawControls(expect);
		}
	}

	/** Checks homogeneous SeqCount ordering by length, then bases. */
	private static void seqCountControls(final String expect){
		checkCompare(expect, "seqcount_equal", 0, seqCount("AA"), seqCount("AA"));
		checkCompare(expect, "seqcount_length", -1, seqCount("A"), seqCount("AA"));
		checkCompare(expect, "seqcount_bases", -1, seqCount("AA"), seqCount("CC"));
	}

	/** Checks homogeneous SeqCountM ordering by count, length, score, then bases. */
	private static void seqCountMControls(final String expect){
		checkCompare(expect, "seqcountm_equal", 0,
				seqCountM("AA", 1, 2), seqCountM("AA", 1, 2));
		checkCompare(expect, "seqcountm_count", -1,
				seqCountM("AA", 1, 2), seqCountM("AA", 2, 2));
		checkCompare(expect, "seqcountm_length", -1,
				seqCountM("A", 1, 2), seqCountM("AA", 1, 1));
		checkCompare(expect, "seqcountm_score", 1,
				seqCountM("AA", 1, 2), seqCountM("CC", 1, 1));
		checkCompare(expect, "seqcountm_bases", -1,
				seqCountM("AA", 1, 2), seqCountM("CC", 1, 2));
	}

	/** Checks the exact SeqCount/SeqCountM zero-consistency failure. */
	private static void mixedZeroConsistencyFailure(final String expect){
		final SeqCountM mutableA=seqCountM("AA", 1, 2);
		final SeqCount baseA=seqCount("AA");
		final SeqCountM mutableC=seqCountM("CC", 1, 1);

		checkCompare(expect, "mixed_mutable_base_equal",
				fixed(expect) ? 1 : 0, mutableA, baseA);
		if(extended(expect)){
			checkCompare(expect, "mixed_base_mutable_equal",
					fixed(expect) ? -1 : 0, baseA, mutableA);
		}
		checkCompare(expect, "mixed_mutable_scored", 1, mutableA, mutableC);
		checkCompare(expect, "mixed_base_scored", -1, baseA, mutableC);
	}

	/** Checks the inherited SeqCountM receiver symmetry failure. */
	private static void inheritedReceiverFailure(final String expect){
		final SeqCountM child=new ChildSeqCountM("AA", 1, 2);
		final SeqCountM exact=seqCountM("CC", 1, 1);

		checkCompare(expect, "inherited_child_exact", 1, child, exact);
		checkCompare(expect, "inherited_exact_child",
				fixed(expect) ? -1 : 1, exact, child);
	}

	/** Checks overflow-safe count comparison after the repair. */
	private static void countOverflowControl(final String expect){
		checkCompare(expect, "count_extremes", fixed(expect) ? -1 : 1,
				seqCountM("AA", Integer.MIN_VALUE, 2),
				seqCountM("AA", Integer.MAX_VALUE, 2));
	}

	/** Checks nonfinite score ordering after the repair. */
	private static void nonfiniteScoreControls(final String expect){
		checkCompare(expect, "score_finite_nan", -1,
				seqCountM("AA", 1, 2), seqCountM("AA", 1, Float.NaN));
		checkCompare(expect, "score_nan_finite", fixed(expect) ? 1 : -1,
				seqCountM("AA", 1, Float.NaN), seqCountM("AA", 1, 2));
		checkCompare(expect, "score_nan_bases", -1,
				seqCountM("AA", 1, Float.NaN), seqCountM("CC", 1, Float.NaN));
		checkCompare(expect, "score_nan_bases_reverse", fixed(expect) ? 1 : -1,
				seqCountM("CC", 1, Float.NaN), seqCountM("AA", 1, Float.NaN));
		checkCompare(expect, "score_infinity_nan", -1,
				seqCountM("AA", 1, Float.POSITIVE_INFINITY),
				seqCountM("AA", 1, Float.NaN));
		checkCompare(expect, "score_signed_zero", 0,
				seqCountM("AA", 1, -0.0f), seqCountM("AA", 1, 0.0f));
	}

	/** Checks ordering laws over representative repaired values. */
	private static void orderingLawControls(final String expect){
		final SeqCount[] values={
			seqCountM("AA", Integer.MIN_VALUE, 0),
			seqCountM("A", 1, 0),
			seqCount("AA"),
			seqCount("CC"),
			seqCountM("AA", 1, -1),
			seqCountM("AA", 1, Float.NEGATIVE_INFINITY),
			seqCountM("AA", 1, -0.0f),
			seqCountM("AA", 1, 0.0f),
			new ChildSeqCountM("CC", 1, 1),
			seqCountM("AA", 1, Float.POSITIVE_INFINITY),
			seqCountM("AA", 1, Float.NaN),
			seqCountM("CC", 1, Float.NaN),
			seqCountM("AA", 2, 0),
			seqCountM("AA", Integer.MAX_VALUE, 0)
		};
		int pairs=0, triples=0;
		for(int i=0; i<values.length; i++){
			for(int j=0; j<values.length; j++){
				final int ij=compare(values[i], values[j]);
				final int ji=compare(values[j], values[i]);
				check(ij==-ji, "pair "+i+","+j+" has signs "+ij+","+ji);
				pairs++;
				for(int k=0; k<values.length; k++){
					final int jk=compare(values[j], values[k]);
					final int ik=compare(values[i], values[k]);
					check(!(ij<=0 && jk<=0) || ik<=0,
							"ascending triple "+i+","+j+","+k+
							" has signs "+ij+","+jk+","+ik);
					check(!(ij>=0 && jk>=0) || ik>=0,
							"descending triple "+i+","+j+","+k+
							" has signs "+ij+","+jk+","+ik);
					triples++;
				}
			}
		}
		System.out.println("PASS\t"+expect+"\tall_pairs\t"+pairs);
		System.out.println("PASS\t"+expect+"\tall_triples\t"+triples);
	}

	/** Checks that equals and hashCode stay sequence-only across base/mutable keys. */
	private static void hashMapControls(final String expect){
		final HashMap<SeqCount, String> baseMap=new HashMap<SeqCount, String>();
		baseMap.put(seqCount("ACGT"), "base");
		checkText(expect, "hash_base_mutable", "base",
				baseMap.get(seqCountM("ACGT", 7, 3)));

		final HashMap<SeqCount, String> mutableMap=new HashMap<SeqCount, String>();
		mutableMap.put(seqCountM("CC", 5, 2), "mutable");
		checkText(expect, "hash_mutable_base", "mutable",
				mutableMap.get(seqCount("CC")));
	}

	/** Checks the CRISPR read comparator with explicit metadata and no network scoring. */
	private static void readComparatorCrisprControls(final String expect){
		checkCompare(expect, "crispr_count", -1,
				read("AA", "count=1 score=9 name=A"),
				read("AA", "count=2 score=0 name=B"));
		checkCompare(expect, "crispr_score", 1,
				read("AA", "count=1 score=2 name=A"),
				read("CC", "count=1 score=1 name=B"));
		checkCompare(expect, "crispr_id", -1,
				read("AA", "count=1 score=2 name=A"),
				read("AA", "count=1 score=2 name=B"));
	}

	/** Creates a sequence key from literal A/C/G/T bases. */
	private static SeqCount seqCount(final String bases){
		return new SeqCount(bytes(bases));
	}

	/** Creates a mutable sequence key from literal bases, count and score. */
	private static SeqCountM seqCountM(final String bases, final int count,
			final float score){
		final SeqCountM value=new SeqCountM(bytes(bases));
		value.count=count;
		value.score=score;
		return value;
	}

	/** Returns US-ASCII bytes for the supplied literal bases. */
	private static byte[] bytes(final String bases){
		return bases.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
	}

	/** Creates a read with literal bases and an explicit count/score header. */
	private static Read read(final String bases, final String id){
		return new Read(bytes(bases), null, id, 1);
	}

	/** Checks a comparison's sign and prints the observed sign. */
	private static void checkCompare(final String expect, final String name,
			final int expected, final SeqCount left, final SeqCount right){
		final int observed=Integer.signum(left.compareTo(right));
		check(expected==observed, name+" expected "+expected+", got "+observed);
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Checks a CRISPR comparator result's sign and prints the observed sign. */
	private static void checkCompare(final String expect, final String name,
			final int expected, final Read left, final Read right){
		final int observed=Integer.signum(
				ReadComparatorCrispr.ascending.compare(left, right));
		check(expected==observed, name+" expected "+expected+", got "+observed);
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Checks text equality and prints the observed value. */
	private static void checkText(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Compares two sequence counters and returns only the sign. */
	private static int compare(final SeqCount left, final SeqCount right){
		return Integer.signum(left.compareTo(right));
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Returns true when running the repaired comparator expectations. */
	private static boolean fixed(final String expect){return "fixed".equals(expect);}

	/** Returns true when running expectations beyond the original 13-row suite. */
	private static boolean extended(final String expect){return !"baseline".equals(expect);}

	/** SeqCountM subclass for inherited-receiver comparison checks. */
	private static final class ChildSeqCountM extends SeqCountM{
		ChildSeqCountM(final String bases, final int count, final float score){
			super(bytes(bases));
			this.count=count;
			this.score=score;
		}
	}

}
