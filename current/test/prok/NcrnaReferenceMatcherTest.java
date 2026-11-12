package test.prok;

import java.util.Arrays;
import java.util.Collections;

import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceMatcher;
import prok.NcrnaReferenceMatcher.Result;

/** Hand-calculated mask and one-to-one expectations; no biological admission.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceMatcherTest{

	static void run(){
		final Locus core=locus("seq", "A", 10, 29, '+');
		final Locus mask=locus("seq", "A", 1, 100, '+');
		final Locus distant=locus("seq", "A", 200, 209, '+');
		final Locus[] calls={core, core, locus("seq", "A", 90, 109, '+'),
			locus("seq", "A", 90, 109, '+'), distant};
		final Result mixed=match(calls, new int[]{1, 1, 1, 1, 1},
			new Locus[]{mask, core, distant}, new int[]{1, 1, 2}, new boolean[]{true, false, true});
		counts(mixed, 1, 2, 0, 1, 2, 2, 2);
		check(Arrays.equals(mixed.callToReference, new int[]{1, -1, -1, -1, -1})
			&& Arrays.equals(mixed.referenceToCall, new int[]{-1, 0, -1}),
			"Only the first core query consumes scored truth; result indices retain the mask before it");
		check(Arrays.equals(mixed.callMaskReference, new int[]{-1, -1, 0, 0, -1}),
			"Mask-only duplicates are both unscored; a mask in another input cannot hide a false positive");
		check(mixed.scoredEligible[0] && mixed.scoredEligible[1] && !mixed.scoredEligible[2],
			"Both scored-eligible duplicates retain eligibility before one-to-one assignment");
		check(!mixed.duplicateCalls[0] && mixed.duplicateCalls[1] && mixed.duplicateCalls[3],
			"Exact duplicates are counted independently of scored or masked outcome");

		final Result isolated=match(new Locus[]{core, locus("seq", "B", 10, 29, '+'),
			locus("seq", "A", 10, 29, '-'), locus("other", "A", 10, 29, '+'), core},
			new int[]{1, 1, 1, 1, 2}, new Locus[]{mask}, new int[]{1}, new boolean[]{true});
		counts(isolated, 0, 4, 0, 0, 1, 1, 0);
		check(Arrays.equals(isolated.callMaskReference, new int[]{0, -1, -1, -1, -1}),
			"Mask qualification requires exact input, family, contig and strand");
		final Result scoredIdentity=match(new Locus[]{core, core, locus("seq", "B", 10, 29, '+')},
			new int[]{2, 1, 1}, new Locus[]{core}, new int[]{1}, new boolean[]{false});
		counts(scoredIdentity, 1, 2, 0, 1, 0, 0, 0);
		check(scoredIdentity.callToReference[1]==0, "Scored matching also enforces input and family identity");

		// C0 can match R0 or R1; C1 can only match R0. A greedy first edge loses one TP.
		final Result reroute=match(new Locus[]{locus("s", "A", 1, 10, '+'), locus("s", "A", 1, 5, '+')},
			new int[]{1, 1}, new Locus[]{locus("s", "A", 1, 8, '+'), locus("s", "A", 4, 10, '+')},
			new int[]{1, 1}, new boolean[]{false, false});
		counts(reroute, 2, 0, 0, 2, 0, 0, 0);
		check(Arrays.equals(reroute.callToReference, new int[]{1, 0}),
			"Maximum-cardinality reassignment must survive mask integration");
		final Result repeated=match(new Locus[]{core, core}, new int[]{1, 1},
			new Locus[]{core, locus("seq", "A", 11, 30, '+'), mask},
			new int[]{1, 1, 1}, new boolean[]{false, false, true});
		counts(repeated, 1, 1, 1, 2, 1, 0, 1);
		check(repeated.scoredEligible[1] && repeated.callMaskReference[1]<0,
			"An exact repeat cannot consume a second overlapping reference or escape through its mask");

		final Result noCalls=match(new Locus[0], new int[0], new Locus[]{core, mask},
			new int[]{1, 1}, new boolean[]{false, true});
		counts(noCalls, 0, 0, 1, 1, 1, 0, 0);
		counts(match(new Locus[]{core}, new int[]{1}, new Locus[0], new int[0], new boolean[0]),
			0, 1, 0, 0, 0, 0, 0);
		counts(match(new Locus[0], new int[0], new Locus[0], new int[0], new boolean[0]),
			0, 0, 0, 0, 0, 0, 0);
		final Result contact=match(new Locus[]{locus("seq", "A", 100, 120, '+'), locus("seq", "A", 101, 120, '+')},
			new int[]{1, 1}, new Locus[]{mask}, new int[]{1}, new boolean[]{true});
		counts(contact, 0, 1, 0, 0, 1, 1, 0);
		check(contact.callMaskReference[0]==0 && contact.callMaskReference[1]<0,
			"Mask overlap is positive inclusive overlap, even one base; adjacency is not overlap");
		expect(()->match(new Locus[]{core}, new int[]{0}, new Locus[0], new int[0], new boolean[0]), "positive");
		expect(()->match(new Locus[]{core}, new int[0], new Locus[0], new int[0], new boolean[0]), "exactly one");
		expect(()->match(new Locus[0], new int[0], new Locus[]{core}, new int[]{1}, new boolean[0]), "exactly one");
		expect(()->match(new Locus[]{locus("s", "A", 1, 2, '.')}, new int[]{1}, new Locus[0], new int[0], new boolean[0]), "Unknown strand");
		for(final float invalid:new float[]{0, -1, Float.NaN, Float.POSITIVE_INFINITY, 1.01f}){
			expect(()->NcrnaReferenceMatcher.match(Collections.<Locus>emptyList(), new int[0],
				Collections.<Locus>emptyList(), new int[0], new boolean[0], invalid), "overlap");
		}
		System.err.println("NCRNA_REFERENCE_MATCHER_TEST_PASS");
	}

	private static Result match(final Locus[] calls, final int[] inputs, final Locus[] references,
		final int[] refInputs, final boolean[] masks){
		return NcrnaReferenceMatcher.match(Arrays.asList(calls), inputs, Arrays.asList(references), refInputs, masks, .5f);
	}
	private static Locus locus(final String contig, final String family, final int start, final int stop, final char strand){
		return new Locus(contig, family, start, stop, strand);
	}
	private static void counts(final Result result, final int tp, final int fp, final int fn,
		final int scored, final int masks, final int unscored, final int duplicates){
		check(result.matched==tp && result.falsePositives==fp && result.falseNegatives==fn
			&& result.scoredReferences==scored && result.maskedReferences==masks
			&& result.maskedCalls==unscored && result.duplicates==duplicates,
			"Expected TP/FP/FN="+tp+"/"+fp+"/"+fn+", scored/mask refs="+scored+"/"+masks
			+", masked calls/duplicates="+unscored+"/"+duplicates+"; observed="+result.matched+"/"
			+result.falsePositives+"/"+result.falseNegatives+"; "+result.scoredReferences+"/"
			+result.maskedReferences+"; "+result.maskedCalls+"/"+result.duplicates);
	}
	private static void expect(final Runnable action, final String fragment){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(fragment), "Unexpected input failure: "+failure);
			return;
		}
		throw new AssertionError("Expected failure containing "+fragment);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
}
