package test.assemble;

import java.util.Arrays;
import java.util.Random;

import assemble.KmerWindowDepth;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import ukmer.Kmer;

/** Checks fixed-anchor span queries against independent base-level fingerprints.
 * @author Brian Bushnell, Fischl */
public final class KmerWindowDepthTest {

	/** Covers word boundaries, reverse complements, canonical ties and all eight neighbors. */
	public static void main(final String[] args){
		Kmer.PACKED=true;
		final Random random=new Random(42121);
		for(int span : new int[]{66, 72, 96, 124, 127, 128, 129, 160}){
			for(int k : new int[]{66, 96, 127, 128, 129, 160}){
				for(int trial=0; trial<20; trial++){
					final byte[] bases=new byte[span];
					for(int i=0; i<span; i++){bases[i]=AminoAcid.numberToBase[random.nextInt(4)];}
					if(trial==0){
						final byte[] arm=Arrays.copyOf(bases, 32);
						AminoAcid.reverseComplementBasesInPlace(arm);
						System.arraycopy(arm, 0, bases, span-32, 32);
					}
					exercise(bases, k);
				}
			}
		}
		boolean rejected=false;
		try{KmerWindowDepth.setHollowQuery(new Kmer(65, true), new byte[128], 0, 128);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "K65 cannot preserve both nonoverlapping33bp neighbor contexts.");
		rejected=false;
		try{KmerWindowDepth.setHollowQuery(new Kmer(128, false), new byte[128], 0, 128);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Synthetic middle bases must never be used with solid hashing.");
		check(KmerWindowDepth.rejectPathStep(10, 0, 10, 4, 3), "Supported G must veto an unsupported saved C in the insertion fixture.");
		check(!KmerWindowDepth.rejectPathStep(1, 0, 10, 4, 3), "An unsupported predecessor must not veto the path.");
		check(!KmerWindowDepth.rejectPathStep(10, 0, 2, 4, 3), "An alternative below minimum support must not veto.");
		check(!KmerWindowDepth.rejectPathStep(10, 3, 11, 4, 3), "Below-ratio support must retain the path.");
		check(KmerWindowDepth.rejectPathStep(10, 3, 12, 4, 3), "Ratio equality must veto consistently.");
		check(!KmerWindowDepth.rejectPathStep(10, 10, 0, 4, 3), "Agreement with the saved next base must retain it.");
		check(!KmerWindowDepth.rejectPathStep(10, 0, 0, 4, 3), "Missing neighbor evidence must not trigger a veto.");
		check(!KmerWindowDepth.rejectPathStep(Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 4, 3),
				"Large counts must not overflow the ratio multiplication.");
		check(KmerWindowDepth.ambiguousPathStep(185, 121, 61, 4, 3), "The insertion-path competing branch must be remembered.");
		check(KmerWindowDepth.ambiguousPathStep(58, 40, 18, 4, 3), "The deletion-path competing branch must be remembered.");
		check(!KmerWindowDepth.ambiguousPathStep(2, 10, 10, 4, 3), "Low-context evidence must not create a checkpoint.");
		check(!KmerWindowDepth.ambiguousPathStep(10, 10, 2, 4, 3), "A weak alternative must not create a checkpoint.");
		check(!KmerWindowDepth.ambiguousPathStep(12, 12, 3, 4, 3), "Chosen-side ratio equality is unambiguous.");
		check(KmerWindowDepth.ambiguousPathStep(12, 12, 4, 4, 3), "Both sufficiently supported sides within the ratio are ambiguous.");
		check(KmerWindowDepth.ambiguousPathStep(12, 4, 12, 4, 3), "Ambiguity must be symmetric in the two supported counts.");
		check(!KmerWindowDepth.ambiguousPathStep(12, 3, 12, 4, 3), "Alternative-side ratio equality is a hard veto, not ambiguity.");
		check(KmerWindowDepth.ambiguousPathStep(Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 4, 3), "Large balanced counts must not overflow.");
		check(!KmerWindowDepth.ambiguousPathStep(10, 0, 10, 4, 3), "An unsupported chosen base must not create an ambiguity checkpoint.");
		check(KmerWindowDepth.ambiguousPathStep(256, 170, 91, 4, 3), "Cohnella's alternate-span branch must be visible without a later contradiction.");
		check(!KmerWindowDepth.rejectPathStep(256, 170, 91, 4, 3), "The Cohnella branch is not a hard veto; the separate branch-stop option is required.");
		check(KmerWindowDepth.ambiguousPathStep(38, 17, 21, 4, 3), "Persicimonas's earlier supported branch must be eligible independently of truth coordinates.");
		check(KmerWindowDepth.ambiguousPathStep(19, 22, 19, 4, 3), "Persicimonas's later wrong-choice branch is also supported ambiguity.");
		nominalPredicates();
		backwardActualFloor();
		lowForwardBranch();
		solidContradiction();
		nominalContradiction();
		nominalWeakChoice();
		nominalRise();
		sourceSupport();
		riseSourceSupport();
		solidBackwardSeed();
		shortDominance();
		seedFork();
		shortDrop();
		predecessorOrientation(random);
		if(args.length>0){writeFixtures(args[0]);}
		System.out.println("KMER_WINDOW_SPAN_TEST_PASS checks="+checks);
	}

	/** Separates local dominance and depth transitions from blanket branch/depth stops. */
	private static void nominalPredicates(){
		check(KmerWindowDepth.rejectBackwardStep(54, 35, 77, 2, 3), "The observed g0048 predecessor is outweighed by its backward alternative.");
		check(KmerWindowDepth.rejectBackwardStep(30, 4, 27, 2, 3), "The observed g0014 predecessor is outweighed by its backward alternative.");
		check(!KmerWindowDepth.rejectBackwardStep(30, 27, 4, 2, 3), "A weak backward alternative must not stop a dominant actual predecessor.");
		check(!KmerWindowDepth.rejectBackwardStep(30, 4, 7, 2, 3), "Backward dominance below ratio two must retain the path.");
		check(KmerWindowDepth.rejectBackwardStep(30, 4, 8, 2, 3), "Backward dominance includes exact ratio equality.");
		check(!KmerWindowDepth.rejectBackwardStep(30, 4, 4, 2, 3), "A balanced backward junction alone is not this diagnostic's target.");
		check(!KmerWindowDepth.rejectBackwardStep(2, 4, 8, 2, 3), "Absent context support must not manufacture backward evidence.");
		check(!KmerWindowDepth.rejectBackwardStep(30, 3, 2, 2, 3), "A backward competitor still needs independent minimum support.");
		check(!KmerWindowDepth.rejectBackwardStep(15, 1, 14, 2, 3), "The known position121 seed error at origin27500 must not veto its exact suffix.");
		check(!KmerWindowDepth.rejectBackwardStep(17, 1, 16, 2, 3), "The known position121 seed error at origin45250 must not veto its exact suffix.");
		check(!KmerWindowDepth.rejectBackwardStep(13, 2, 11, 2, 3), "Actual predecessor depth just below the support floor must abstain.");
		check(!KmerWindowDepth.rejectBackwardStep(6, 0, 6, 2, 3), "A missing actual predecessor alone does not establish two supported incoming paths.");
		check(KmerWindowDepth.rejectBackwardStep(6, 3, 6, 2, 3), "Actual support equality must retain the intended backward veto.");
		check(!KmerWindowDepth.rejectBackwardStep(Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 2, 3),
				"Large backward counts must not overflow the ratio multiplication.");
		check(KmerWindowDepth.dropPathStep(80, 2, 16, 32, 2), "The observed g0058 depth80-to2 collapse must be detected.");
		check(KmerWindowDepth.dropPathStep(32, 2, 16, 32, 2), "Both drop threshold equalities must be inclusive.");
		check(!KmerWindowDepth.dropPathStep(31, 1, 16, 32, 2), "A transition below the context floor must retain the path.");
		check(!KmerWindowDepth.dropPathStep(32, 2, 17, 32, 2), "A transition below the configured ratio must retain the path.");
		check(!KmerWindowDepth.dropPathStep(80, 3, 16, 32, 2), "Depth above the low-depth ceiling must retain the path.");
		check(!KmerWindowDepth.dropPathStep(2, 2, 16, 32, 2), "A persistently shallow path must not be globally prohibited.");
		check(KmerWindowDepth.dropPathStep(32, 0, 16, 32, 2), "A supported context followed by no chosen support is a collapse.");
		check(!KmerWindowDepth.dropPathStep(0, 0, 16, 32, 2), "Two absent keys do not provide evidence for a collapse.");
		check(!KmerWindowDepth.dropPathStep(Integer.MAX_VALUE, Integer.MAX_VALUE, 2, 1, Integer.MAX_VALUE),
				"Depth-drop multiplication must not overflow when thresholds are large.");
		check(KmerWindowDepth.risePathStep(6, 238, 16, 3, 32), "The g0857 transition into fourteen genomic endpoint copies must be detectable without reference truth.");
		check(KmerWindowDepth.risePathStep(3, 48, 16, 3, 32), "Support and rise-ratio equality must be inclusive.");
		check(KmerWindowDepth.risePathStep(4, 32, 8, 3, 32), "The next-depth floor equality must be inclusive.");
		check(!KmerWindowDepth.risePathStep(3, 47, 16, 3, 32), "A rise below the ratio must retain the path.");
		check(!KmerWindowDepth.risePathStep(3, 31, 8, 3, 32), "A rise below the absolute depth floor must retain the path.");
		check(!KmerWindowDepth.risePathStep(2, 238, 16, 3, 32), "A weak seed-error context must not establish a supported depth rise.");
		check(!KmerWindowDepth.risePathStep(0, 238, 16, 3, 32), "Missing context must not become a vacuous ratio comparison.");
		check(!KmerWindowDepth.risePathStep(238, 238, 16, 3, 32), "An already deep path must not be rejected for high depth alone.");
		check(!KmerWindowDepth.risePathStep(Integer.MAX_VALUE, Integer.MAX_VALUE, 16, 3, 32), "Rise multiplication must not overflow.");
	}

	/** Keeps weak-actual experiments separate from the established context/competitor floor. */
	private static void backwardActualFloor(){
		check(KmerWindowDepth.rejectBackwardStep(11, 1, 10, 2, 3, 1), "The diagnosed g0385 predecessor becomes eligible only with an explicit weak-actual floor.");
		check(!KmerWindowDepth.rejectBackwardStep(11, 1, 10, 2, 3), "The original API must retain its coupled support requirement.");
		check(!KmerWindowDepth.rejectBackwardStep(2, 1, 10, 2, 3, 1), "Lowering actual support must not lower context support.");
		check(!KmerWindowDepth.rejectBackwardStep(11, 1, 2, 2, 3, 1), "Lowering actual support must not lower competitor support.");
		check(!KmerWindowDepth.rejectBackwardStep(11, 0, 10, 2, 3, 1), "A missing actual key must remain ineligible.");
		check(!KmerWindowDepth.rejectBackwardStep(11, 1, 10, 2, 3, 2), "Actual support below its independently configured floor must abstain.");
		check(KmerWindowDepth.rejectBackwardStep(3, 2, 4, 2, 3, 2), "Context, actual-floor and ratio equalities must be inclusive.");
		check(!KmerWindowDepth.rejectBackwardStep(3, 2, 3, 2, 3, 2), "Support without ratio dominance must retain the path.");
		check(KmerWindowDepth.rejectBackwardStep(15, 1, 14, 2, 3, 1), "The aggressive diagnostic also stops the known benign seed error; biology must measure that cost.");
		check(!KmerWindowDepth.rejectBackwardStep(Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 2, 3, 1),
				"Independent actual support must not change overflow-safe ratio arithmetic.");
		for(int floor=1; floor<=4; floor++){
			for(int actual=0; actual<=5; actual++){
				check(KmerWindowDepth.rejectBackwardStep(6, actual, 8, 2, floor)==
						KmerWindowDepth.rejectBackwardStep(6, actual, 8, 2, floor, floor),
						"Inheriting the support floor must preserve every original API decision.");
			}
		}
	}

	/** Separates count2/3 competitors from count1 noise, weak contexts and deeper branches. */
	private static void lowForwardBranch(){
		check(KmerWindowDepth.rejectLowBranchStep(37, 34, 2, 32, 3), "The g0680 wrong34/true2 branch must be eligible.");
		check(!KmerWindowDepth.rejectLowBranchStep(37, 34, 1, 32, 3), "Depth1 tolerance must remain unchanged.");
		check(!KmerWindowDepth.rejectLowBranchStep(37, 34, 0, 32, 3), "No competitor provides no branch evidence.");
		check(KmerWindowDepth.rejectLowBranchStep(37, 34, 3, 32, 3), "The upper minor-count bound is inclusive.");
		check(!KmerWindowDepth.rejectLowBranchStep(37, 34, 4, 32, 3), "Deeper competitors are outside this narrow experiment.");
		check(!KmerWindowDepth.rejectLowBranchStep(64, 64, 2, 32, 3), "Count2 ratio equality retains the saved choice.");
		check(KmerWindowDepth.rejectLowBranchStep(63, 63, 2, 32, 3), "One count below ratio equality must stop.");
		check(!KmerWindowDepth.rejectLowBranchStep(96, 96, 3, 32, 3), "Count3 ratio equality also retains the saved choice.");
		check(KmerWindowDepth.rejectLowBranchStep(95, 95, 3, 32, 3), "Count3 below-ratio evidence must stop.");
		check(!KmerWindowDepth.rejectLowBranchStep(2, 34, 2, 32, 3), "Weak context cannot trigger this veto.");
		check(!KmerWindowDepth.rejectLowBranchStep(37, 2, 2, 32, 3), "The saved choice must meet its independent support floor.");
		check(KmerWindowDepth.rejectLowBranchStep(37, 3, 2, 32, 3), "Chosen support equality is eligible.");
		check(KmerWindowDepth.rejectLowBranchStep(Integer.MAX_VALUE, Integer.MAX_VALUE, 2, Integer.MAX_VALUE, 3), "Ratio multiplication must not overflow.");
		check(!KmerWindowDepth.rejectLowBranchStep(Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 32, 3), "Very deep alternatives remain outside the explicit count range.");
	}

	/** Checks that solid support can veto a contradiction but never impose coverage. */
	private static void solidContradiction(){
		check(KmerWindowDepth.rejectSolidStep(2, 0, 2, 2), "The selected g0680 contradiction meets both support equalities.");
		check(!KmerWindowDepth.rejectSolidStep(0, 0, 2, 2), "Absent context cannot establish a solid contradiction.");
		check(!KmerWindowDepth.rejectSolidStep(1, 0, 2, 2), "One context observation is insufficient.");
		check(!KmerWindowDepth.rejectSolidStep(2, 0, 1, 2), "One alternative observation is insufficient.");
		check(!KmerWindowDepth.rejectSolidStep(2, 0, 0, 2), "A dead end alone must retain the saved path.");
		check(!KmerWindowDepth.rejectSolidStep(0, 0, 0, 2), "Missing solid coverage must not erase Hollow's benefit.");
		check(!KmerWindowDepth.rejectSolidStep(2, 1, 20, 2), "Any saved-next support makes this conservative veto abstain.");
		check(!KmerWindowDepth.rejectSolidStep(10, 10, 0, 2), "The first correct same-Hollow-triple case must be retained.");
		check(!KmerWindowDepth.rejectSolidStep(11, 9, 1, 2), "The second correct same-Hollow-triple case must be retained.");
		check(!KmerWindowDepth.rejectSolidStep(2, 0, 2, 3), "A higher floor must be honored on both counts.");
		check(KmerWindowDepth.rejectSolidStep(Integer.MAX_VALUE, 0, Integer.MAX_VALUE, 2), "Large counts must remain valid without overflow.");
	}

	/** Tests the read-supported Hollow160 contradiction independently of truth coordinates. */
	private static void nominalContradiction(){
		check(KmerWindowDepth.rejectUnsupportedStep(18, 0, 19, 3), "The diagnosed first wrong g0643 choice has independent support for an alternative.");
		check(!KmerWindowDepth.rejectUnsupportedStep(16, 16, 0, 3), "The same saved choice remains supported at Hollow128; span matters.");
		check(KmerWindowDepth.rejectUnsupportedStep(3, 0, 3, 3), "Both nominal support equalities are inclusive.");
		check(!KmerWindowDepth.rejectUnsupportedStep(2, 0, 3, 3), "Weak nominal context must abstain.");
		check(!KmerWindowDepth.rejectUnsupportedStep(3, 0, 2, 3), "Weak alternative support must abstain.");
		check(!KmerWindowDepth.rejectUnsupportedStep(3, 1, 100, 3), "Even one observation of the saved choice makes this diagnostic abstain.");
		check(!KmerWindowDepth.rejectUnsupportedStep(0, 0, 0, 3), "A coverage gap is not a contradiction.");
		check(!KmerWindowDepth.rejectUnsupportedStep(18, 0, 0, 3), "A supported dead end alone must not discard the saved suffix.");
	}

	/** Restricts the new nominal experiment to one noisy observation and a dominant alternative. */
	private static void nominalWeakChoice(){
		check(KmerWindowDepth.rejectNominalWeakStep(8, 1, 7, 4, 3), "The remaining g0600 wrong choice has count1 opposed by seven observations.");
		check(KmerWindowDepth.rejectNominalWeakStep(3, 1, 4, 4, 3), "Context-floor and ratio equalities must be inclusive.");
		check(!KmerWindowDepth.rejectNominalWeakStep(3, 1, 3, 4, 3), "Supported alternatives below fourfold must abstain.");
		check(!KmerWindowDepth.rejectNominalWeakStep(2, 1, 4, 4, 3), "Weak context must abstain despite alternative dominance.");
		check(!KmerWindowDepth.rejectNominalWeakStep(3, 1, 2, 2, 3), "Ratio equality cannot bypass the independent alternative floor.");
		check(!KmerWindowDepth.rejectNominalWeakStep(3, 0, 4, 4, 3), "Zero choices retain the existing contradiction predicate and reason.");
		check(!KmerWindowDepth.rejectNominalWeakStep(3, 2, 100, 4, 3), "Count2 choices remain outside this narrow experiment.");
		check(!KmerWindowDepth.rejectNominalWeakStep(0, 1, 100, 4, 3), "Absent context is not evidence of a supported transition.");
		check(!KmerWindowDepth.rejectNominalWeakStep(3, 1, 0, 4, 3), "No alternative must retain the saved choice.");
		check(!KmerWindowDepth.rejectNominalWeakStep(3, 1, 6, 7, 3), "The configured ratio must be honored above its default.");
		check(KmerWindowDepth.rejectNominalWeakStep(3, 1, 7, 7, 3), "Equality at a nondefault ratio remains inclusive.");
		check(KmerWindowDepth.rejectNominalWeakStep(Integer.MAX_VALUE, 1, Integer.MAX_VALUE, Integer.MAX_VALUE, 3), "Large ratio equality must retain overflow-safe arithmetic.");
	}

	/** Reuses the established rise thresholds for a separately tested auxiliary geometry. */
	private static void nominalRise(){
		check(KmerWindowDepth.risePathStep(12, 231, 16, 3, 32), "The g0859 auxiliary span enters nineteen genomic endpoint copies at its wrong C.");
		check(KmerWindowDepth.risePathStep(12, 192, 16, 3, 32), "Auxiliary rise equality must stop before the new base.");
		check(!KmerWindowDepth.risePathStep(12, 191, 16, 3, 32), "One count below the auxiliary rise ratio must abstain.");
		check(!KmerWindowDepth.risePathStep(2, 231, 16, 3, 32), "Weak context may be seed-error recovery and must abstain.");
		check(!KmerWindowDepth.risePathStep(12, 31, 2, 3, 32), "The next-depth floor remains independent of the ratio.");
		check(!KmerWindowDepth.risePathStep(12, 12, 16, 3, 32), "Flat auxiliary coverage is not a rise.");
		check(KmerWindowDepth.risePathStep(3, 48, 16, 3, 32), "Nominal scanning must retain the existing support and ratio equalities.");
		check(!KmerWindowDepth.risePathStep(Integer.MAX_VALUE, Integer.MAX_VALUE, 16, 3, 32), "Large counts must not wrap the auxiliary comparison.");
	}

	/** Protects tied or stronger source choices without changing auxiliary support floors. */
	private static void sourceSupport(){
		check(!KmerWindowDepth.sourceChoiceOutvoted(8, 17), "The false g0680 veto must not overrule stronger source support.");
		check(KmerWindowDepth.sourceChoiceOutvoted(19, 16), "The g0643 alternative strictly outweighs the saved source choice.");
		check(KmerWindowDepth.sourceChoiceOutvoted(11, 2), "The g0600 shifted path is weak at its source span.");
		check(!KmerWindowDepth.sourceChoiceOutvoted(12, 16), "The second g0600 veto also abstains under this more conservative comparison.");
		check(!KmerWindowDepth.sourceChoiceOutvoted(3, 3), "Equal counts across distinct spans must not overrule the source path.");
		check(!KmerWindowDepth.sourceChoiceOutvoted(3, 4), "One additional source observation must prevent this veto.");
		check(!KmerWindowDepth.sourceChoiceOutvoted(3, 0), "An absent source key provides no supported choice to compare.");
		check(!KmerWindowDepth.sourceChoiceOutvoted(Integer.MAX_VALUE, Integer.MAX_VALUE), "Large equal counts still abstain without arithmetic overflow.");
		check(!KmerWindowDepth.sourceChoiceOutvoted(4, 4), "The g0490 mixed-error false stop compares equally supported different repeat copies.");
		check(KmerWindowDepth.sourceChoiceOutvoted(4, 3), "One additional alternative observation meets strict dominance.");
		check(!KmerWindowDepth.sourceChoiceOutvoted(0, 1), "No alternative evidence cannot outweigh a supported source choice.");
		check(KmerWindowDepth.sourceChoiceOutvoted(Integer.MAX_VALUE, Integer.MAX_VALUE-1), "Strict dominance must remain valid at the largest counts.");
	}

	/** Checks source-choice boundaries independently of the auxiliary rise predicate. */
	private static void riseSourceSupport(){
		check(!KmerWindowDepth.weakSourceChoice(0, 3), "Missing source evidence must not establish a weak supported choice.");
		check(KmerWindowDepth.weakSourceChoice(1, 3), "One source observation is positive but below the support floor.");
		check(KmerWindowDepth.weakSourceChoice(2, 3), "The diagnosed wrong g0859 choice remains weak at its source span.");
		check(!KmerWindowDepth.weakSourceChoice(3, 3), "Support equality protects legitimate chosen-key recovery.");
		check(!KmerWindowDepth.weakSourceChoice(579, 3), "A strong choice is protected even when its source context was weak.");
		check(!KmerWindowDepth.weakSourceChoice(Integer.MAX_VALUE, 3), "High counts must not become weak through overflow.");
		check(KmerWindowDepth.weakSourceChoice(1, 2), "The configured support floor is used rather than a hardcoded three.");
		check(!KmerWindowDepth.weakSourceChoice(2, 2), "Equality at a different support floor still abstains.");
	}

	/** Keeps the fixed seed-start conjunction distinct from blanket low-Solid rejection. */
	private static void solidBackwardSeed(){
		check(KmerWindowDepth.rejectSolidBackwardSeed(11, 1, 10, 1), "The selected bad g0385 seed lacks corroborating full Solid support.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(15, 1, 14, 5), "The correct g0058 origin27500 seed has Solid depth5.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(17, 1, 16, 10), "The correct g0058 origin45250 seed has Solid depth10.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(13, 2, 11, 6), "The correct g0058 origin89250 seed has Solid depth6.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(11, 1, 10, 2), "Solid support equality retains the saved extension.");
		check(KmerWindowDepth.rejectSolidBackwardSeed(3, 2, 4, 0), "Zero Solid support is eligible only with the independent Hollow conjunction.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(2, 1, 10, 0), "Weak Hollow context must abstain despite absent Solid support.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(11, 1, 2, 0), "A weak competitor must abstain.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(11, 0, 10, 0), "An absent actual predecessor is outside this experiment.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(11, 3, 10, 0), "An already-supported predecessor is not a weak-actual candidate.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(11, 2, 3, 0), "Below-ratio evidence must retain the saved path.");
		check(!KmerWindowDepth.rejectSolidBackwardSeed(Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 0), "Large counts must neither overflow nor enter the weak-actual range.");
	}

	/** Requires agreement of weak long support and strong opposing short support. */
	private static void shortDominance(){
		check(KmerWindowDepth.rejectShortStep(2, 42, 102), "The diagnosed wrong repeat choice is weak-long and short-dominated.");
		check(!KmerWindowDepth.rejectShortStep(2, 172, 3), "The correct shallow origin534000 path dominates its short alternative.");
		check(!KmerWindowDepth.rejectShortStep(2, 209, 3), "The two correct bases before the diagnosed error must remain.");
		check(!KmerWindowDepth.rejectShortStep(2, 20, 24), "A balanced short split below twofold does not meet this fixed conjunction.");
		check(!KmerWindowDepth.rejectShortStep(3, 3, 100), "Adequate Hollow support abstains despite opposing short evidence.");
		check(KmerWindowDepth.rejectShortStep(2, 3, 6), "Both short support and ratio equalities are inclusive.");
		check(!KmerWindowDepth.rejectShortStep(2, 3, 5), "A competitor below twofold must retain the path.");
		check(!KmerWindowDepth.rejectShortStep(2, 2, 100), "A weak short chosen key does not establish this repeat competition.");
		check(!KmerWindowDepth.rejectShortStep(0, 0, 100), "Absent short evidence abstains; this is not a coverage requirement.");
		check(KmerWindowDepth.rejectShortStep(0, 3, 6), "Absent Hollow support still requires both supported short paths.");
		check(!KmerWindowDepth.rejectShortStep(2, Integer.MAX_VALUE, Integer.MAX_VALUE), "The twofold comparison must not overflow.");
		check(!KmerWindowDepth.rejectShortStep(2, 3, 0), "No short alternative means no contradiction.");
	}

	/** Checks short-fork boundaries independently of the seed-candidate gate. */
	private static void seedFork(){
		check(KmerWindowDepth.rejectSeedFork(24, 28, 23), "The selected bad path has two supported near-balanced short choices.");
		check(!KmerWindowDepth.rejectSeedFork(12, 11, 0), "The selected correct path has no supported alternative.");
		check(KmerWindowDepth.rejectSeedFork(3, 6, 3), "Support floors and exact twofold equality must stop.");
		check(!KmerWindowDepth.rejectSeedFork(3, 7, 3), "A weaker than half-depth alternative must abstain.");
		check(!KmerWindowDepth.rejectSeedFork(2, 6, 3), "Weak context must abstain.");
		check(!KmerWindowDepth.rejectSeedFork(3, 2, 100), "Weak chosen support must abstain.");
		check(!KmerWindowDepth.rejectSeedFork(3, 3, 2), "Weak alternative support must abstain.");
		check(!KmerWindowDepth.rejectSeedFork(0, 0, 0), "Absence alone must never truncate.");
		check(KmerWindowDepth.rejectSeedFork(3, 3, 100), "Stronger alternatives remain eligible.");
		check(KmerWindowDepth.rejectSeedFork(Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE), "The comparison must not overflow.");
		check(KmerWindowDepth.seedForkHorizon(250, 350, 128)==314, "The bound is the ignored interior width, not all saved bases.");
		check(KmerWindowDepth.seedForkHorizon(250, 281, 128)==281, "A short saved path bounds the inspection.");
		check(KmerWindowDepth.seedForkHorizon(250, 250, 128)==250, "An unextended seed has no inspected choices.");
		check(KmerWindowDepth.seedForkHorizon(Integer.MAX_VALUE-10, Integer.MAX_VALUE, 128)==Integer.MAX_VALUE, "Seed plus horizon must not overflow.");
	}

	/** Separates a supported transition from the false short-key stops seen in the panel. */
	private static void shortDrop(){
		check(KmerWindowDepth.rejectShortDrop(148, 13, 44, 0, 44), "The diagnosed wrong g0772 choice has both independent signals.");
		for(int[] pair : new int[][]{{5,6},{18,17},{24,24},{19,19},{14,15},{15,15},{23,23},
				{14,14},{15,13},{19,17},{17,17},{95,94},{14,14},{21,20}}){
			check(!KmerWindowDepth.rejectShortDrop(pair[0], pair[1], 44, 0, 44),
					"Every observed false stopping transition must abstain even with strong short contradiction.");
		}
		check(KmerWindowDepth.rejectShortDrop(32, 8, 2, 0, 2), "All fixed threshold equalities are inclusive.");
		check(!KmerWindowDepth.rejectShortDrop(31, 1, 2, 0, 2), "Low Hollow context must abstain.");
		check(!KmerWindowDepth.rejectShortDrop(32, 9, 2, 0, 2), "A drop below fourfold must abstain.");
		check(!KmerWindowDepth.rejectShortDrop(32, 8, 1, 0, 2), "One short context observation is insufficient.");
		check(!KmerWindowDepth.rejectShortDrop(32, 8, 2, 1, 20), "Any short chosen support must abstain.");
		check(!KmerWindowDepth.rejectShortDrop(32, 8, 2, 0, 1), "One alternative observation is insufficient.");
		check(!KmerWindowDepth.rejectShortDrop(32, 0, 0, 0, 0), "Missing short coverage must not discard a Hollow path.");
		check(KmerWindowDepth.rejectShortDrop(32, 0, 2, 0, 2), "Supported-to-zero Hollow collapse still needs a short contradiction.");
		check(!KmerWindowDepth.rejectShortDrop(Integer.MAX_VALUE, Integer.MAX_VALUE, 2, 0, 2), "The ratio multiplication must not overflow.");
	}

	/** Reconstructs the actual predecessor in both orientations, independent of key canonicalization. */
	private static void predecessorOrientation(final Random random){
		assert(random!=null) : "A fixed test seed makes the orientation cases reproducible.";
		for(int k : new int[]{66, 96, 124, 128, 160}){
			for(int previous=0; previous<4; previous++){
				final byte[] bases=new byte[k+1];
				for(int i=0; i<bases.length; i++){bases[i]=AminoAcid.numberToBase[random.nextInt(4)];}
				bases[0]=AminoAcid.numberToBase[previous];
				for(int orientation=0; orientation<2; orientation++){
					final Kmer current=new Kmer(k, true), expected=new Kmer(k, true);
					for(int i=1; i<=k; i++){current.addRight(bases[i]);}
					for(int i=0; i<k; i++){expected.addRight(bases[i]);}
					current.addLeft(bases[0]);
					check(current.xor()==expected.xor() && current.xor2()==expected.xor2(),
							"The actual predecessor must use the query's preceding base, not a canonical-orientation base.");
					AminoAcid.reverseComplementBasesInPlace(bases);
				}
			}
		}
	}

	/** Writes exact-count CLI fixtures; each evidence record contains just one nominal key. */
	private static void writeFixtures(final String directory){
		assert(directory!=null && !directory.isEmpty()) : "The test launcher supplies a fresh fixture directory, never an implicit working path.";
		final Random random=new Random(421914);
		final byte[] query=new byte[133];
		for(int i=0; i<query.length; i++){query[i]=AminoAcid.numberToBase[random.nextInt(4)];}
		query[1]='C';
		final byte[] seed=Arrays.copyOf(query, 130);
		final byte[] context=Arrays.copyOfRange(query, 2, 130);
		final byte[] next=Arrays.copyOfRange(query, 3, 131);
		final byte[] actual=Arrays.copyOfRange(query, 1, 129), other=actual.clone();
		other[0]='A';
		writeFasta(directory+"/seed.fa", new byte[][]{seed}, new int[]{1}, true);
		writeFasta(directory+"/query.fa", new byte[][]{query}, new int[]{1}, true);
		writeFasta(directory+"/window_queries.fa", new byte[][]{query, Arrays.copyOf(query, 131), seed}, new int[]{1, 1, 1}, false);
		final byte[] ambiguous=query.clone();
		ambiguous[131]='N';
		writeFasta(directory+"/window_ambiguous.fa", new byte[][]{ambiguous}, new int[]{1}, true);
		final ByteStreamWriter malformed=new ByteStreamWriter(directory+"/window_tab_id.fa", false, false, false);
		malformed.start();
		try{malformed.print(new ByteBuilder().append(">bad\tidentifier\n").append(query).nl());}
		finally{
			if(malformed.poisonAndWait()){throw new RuntimeException("Malformed-ID fixture output failed.");}
		}
		writeFasta(directory+"/back.fa", new byte[][]{context, actual, other, next}, new int[]{6, 3, 6, 2}, false);
		writeFasta(directory+"/weak1.fa", new byte[][]{context, actual, other, next}, new int[]{6, 1, 6, 2}, false);
		writeFasta(directory+"/weak2.fa", new byte[][]{context, actual, other, next}, new int[]{6, 2, 6, 2}, false);
		writeFasta(directory+"/missing.fa", new byte[][]{context, other, next}, new int[]{6, 6, 2}, false);
		writeFasta(directory+"/drop.fa", new byte[][]{context, next}, new int[]{32, 2}, false);
		writeFasta(directory+"/shallow.fa", new byte[][]{context, next}, new int[]{2, 2}, false);
		writeFasta(directory+"/below.fa", new byte[][]{context, next}, new int[]{31, 1}, false);
		writeFasta(directory+"/above.fa", new byte[][]{context, next}, new int[]{80, 3}, false);
		writeFasta(directory+"/rise.fa", new byte[][]{context, next}, new int[]{3, 48}, false);
		writeFasta(directory+"/rise_below.fa", new byte[][]{context, next}, new int[]{3, 47}, false);
		writeFasta(directory+"/rise_weak.fa", new byte[][]{context, next}, new int[]{2, 48}, false);
		final byte[] otherNext=next.clone();
		otherNext[127]=next[127]=='A' ? (byte)'C' : (byte)'A';
		for(int minor=1; minor<=4; minor++){
			writeFasta(directory+"/lowbranch"+minor+".fa", new byte[][]{context, next, otherNext}, new int[]{37, 34, minor}, false);
		}
		writeFasta(directory+"/lowbranch_equal.fa", new byte[][]{context, next, otherNext}, new int[]{64, 64, 2}, false);
		writeFasta(directory+"/lowbranch_weak.fa", new byte[][]{context, next, otherNext}, new int[]{2, 34, 2}, false);
		final byte[] secondNext=Arrays.copyOfRange(query, 4, 132), otherSecond=secondNext.clone();
		otherSecond[127]=secondNext[127]=='A' ? (byte)'C' : (byte)'A';
		writeFasta(directory+"/lowbranch_later.fa", new byte[][]{context, next, otherNext, secondNext, otherSecond}, new int[]{37, 34, 1, 34, 2}, false);
		writeFasta(directory+"/one_step.fa", new byte[][]{Arrays.copyOf(query, 131)}, new int[]{1}, true);
		writeFasta(directory+"/solid_conflict.fa", new byte[][]{context, otherNext}, new int[]{2, 2}, false);
		writeFasta(directory+"/solid_supported.fa", new byte[][]{context, next, otherNext}, new int[]{2, 1, 2}, false);
		writeFasta(directory+"/solid_weak_context.fa", new byte[][]{context, otherNext}, new int[]{1, 2}, false);
		final byte[] hollowTwin=context.clone();
		hollowTwin[64]=hollowTwin[64]=='A' ? (byte)'C' : (byte)'A';
		writeFasta(directory+"/solid_hollow_twin.fa", new byte[][]{hollowTwin}, new int[]{10}, false);
		writeFasta(directory+"/solid_weak_other.fa", new byte[][]{context, otherNext}, new int[]{2, 1}, false);
		writeFasta(directory+"/solid_later.fa", new byte[][]{context, next, otherSecond}, new int[]{2, 2, 2}, false);
		writeFasta(directory+"/solid_regain.fa", new byte[][]{next, otherSecond}, new int[]{2, 2}, false);
		final byte[] otherThird=Arrays.copyOfRange(query, 5, 133);
		otherThird[127]=otherThird[127]=='A' ? (byte)'C' : (byte)'A';
		writeFasta(directory+"/solid_restore.fa", new byte[][]{context, otherNext, secondNext, otherThird}, new int[]{2, 1, 2, 2}, false);
		writeFasta(directory+"/two_steps.fa", new byte[][]{Arrays.copyOf(query, 132)}, new int[]{1}, true);
		final byte[] shortContext=Arrays.copyOf(query, 128), shortOther=Arrays.copyOfRange(query, 1, 129);
		shortOther[127]=shortOther[127]=='A' ? (byte)'C' : (byte)'A';
		writeFasta(directory+"/solid_short.fa", new byte[][]{shortContext, shortOther}, new int[]{2, 2}, false);
		writeFasta(directory+"/short_seed.fa", new byte[][]{shortContext}, new int[]{1}, true);
		final byte[] unrelated=query.clone();
		for(int i=0; i<unrelated.length; i++){unrelated[i]=AminoAcid.numberToBase[random.nextInt(4)];}
		writeFasta(directory+"/absent.fa", new byte[][]{Arrays.copyOf(unrelated, 128)}, new int[]{1}, false);
		writeShortFixtures(directory, context, next, secondNext);
		writeCandidates(directory);
		writeSeedForkFixtures(directory);
		writeShortDropFixtures(directory);
		writeNominalFixtures(directory);
		writeNominalWeakFixtures(directory);
		writeSpacingFixtures(directory);
	}

	/** Isolates spacing protection from nominal contradictions with exact192bp records.
	 * The two observed spacings are encoded independently with33bp arms; every
	 * evidence copy contributes one table key and no extended query contributes counts. */
	private static void writeSpacingFixtures(final String directory){
		assert(directory!=null) : "Spacing fixtures must use the launcher's fresh test directory.";
		final byte[] query=new byte[280];
		final Random random=new Random(1007);
		for(int i=0; i<query.length; i++){query[i]=AminoAcid.numberToBase[random.nextInt(4)];}
		writeFasta(directory+"/spacing_query.fa", new byte[][]{query}, new int[]{1}, true);
		writeFasta(directory+"/spacing_seed.fa", new byte[][]{Arrays.copyOf(query, 200)}, new int[]{1}, true);
		writeFasta(directory+"/spacing_batch.fa", new byte[][]{query, Arrays.copyOf(query, 200), query}, new int[]{1, 1, 1}, false);
		final byte[] ambiguous=query.clone();
		ambiguous[7]='N';//Only the k+1 context includes this base at the first choice.
		writeFasta(directory+"/spacing_ambiguous.fa", new byte[][]{ambiguous}, new int[]{1}, true);
		for(int span : new int[]{191, 193}){
			for(int end : new int[]{200, 202, 279}){
				final byte[] context=Arrays.copyOfRange(query, end-192, end);
				final byte[] other=Arrays.copyOfRange(query, end-191, end+1);
				other[191]=other[191]=='A' ? (byte)'C' : (byte)'A';
				final byte[] spaced=spacingRecord(query, end-span, span);
				final byte[] choice=spacingRecord(query, end-span+1, span);
				writeFasta(directory+"/spacing_"+span+"_"+end+".fa", new byte[][]{context, other, spaced, choice}, new int[]{3, 3, 3, 3}, false);
				if(end==200){
					writeFasta(directory+"/spacing_"+span+"_weak_context.fa", new byte[][]{context, other, spaced, choice}, new int[]{3, 3, 2, 3}, false);
					writeFasta(directory+"/spacing_"+span+"_weak_choice.fa", new byte[][]{context, other, spaced, choice}, new int[]{3, 3, 3, 2}, false);
					writeFasta(directory+"/spacing_"+span+"_absent.fa", new byte[][]{context, other}, new int[]{3, 3}, false);
					final byte[] later=Arrays.copyOfRange(query, 87, 279);
					final byte[] laterOther=Arrays.copyOfRange(query, 88, 280);
					laterOther[191]=laterOther[191]=='A' ? (byte)'C' : (byte)'A';
					writeFasta(directory+"/spacing_"+span+"_later.fa", new byte[][]{context, other, spaced, choice, later, laterOther}, new int[]{3, 3, 3, 3, 3, 3}, false);
					final byte[] source=Arrays.copyOfRange(query, 73, 201);
					writeFasta(directory+"/spacing_"+span+"_source.fa", new byte[][]{context, other, spaced, choice, source}, new int[]{3, 3, 3, 3, 2}, false);
					AminoAcid.reverseComplementBasesInPlace(context);
					AminoAcid.reverseComplementBasesInPlace(other);
					AminoAcid.reverseComplementBasesInPlace(spaced);
					AminoAcid.reverseComplementBasesInPlace(choice);
					writeFasta(directory+"/spacing_"+span+"_reverse.fa", new byte[][]{context, other, spaced, choice}, new int[]{3, 3, 3, 3}, false);
				}
			}
		}
	}

	/** Builds independent fixed-span evidence from the observed arms, including neighbor context. */
	private static byte[] spacingRecord(final byte[] query, final int start, final int span){
		assert(span>=66 && start>=0 && start+span<=query.length) :
				"The fixture copies two complete33bp contexts so neither neighbor samples synthetic middle bases.";
		final byte[] record=new byte[192];
		Arrays.fill(record, (byte)'A');
		System.arraycopy(query, start, record, 0, 33);
		System.arraycopy(query, start+span-33, record, 159, 33);
		return record;
	}

	/** Exact192bp records isolate count1 evidence at the first, later and final saved choices. */
	private static void writeNominalWeakFixtures(final String directory){
		assert(directory!=null) : "Weak-choice fixtures use the same fresh output directory as the prior nominal tests.";
		final byte[] query=new byte[280];
		final Random random=new Random(995);
		for(int i=0; i<query.length; i++){query[i]=AminoAcid.numberToBase[random.nextInt(4)];}
		writeFasta(directory+"/weak_query.fa", new byte[][]{query}, new int[]{1}, true);
		writeFasta(directory+"/weak_seed.fa", new byte[][]{Arrays.copyOf(query, 200)}, new int[]{1}, true);
		writeFasta(directory+"/weak_batch.fa", new byte[][]{query, Arrays.copyOf(query, 200)}, new int[]{1, 1}, false);
		for(int end : new int[]{200, 202, 279}){
			final byte[] context=Arrays.copyOfRange(query, end-192, end);
			final byte[] chosen=Arrays.copyOfRange(query, end-191, end+1), other=chosen.clone();
			context[96]=context[96]=='A' ? (byte)'C' : (byte)'A';
			other[191]=other[191]=='A' ? (byte)'C' : (byte)'A';
			writeFasta(directory+"/weak_"+end+".fa", new byte[][]{context, chosen, other}, new int[]{3, 1, 4}, false);
			if(end==200){
				writeFasta(directory+"/weak_below.fa", new byte[][]{context, chosen, other}, new int[]{3, 1, 3}, false);
				writeFasta(directory+"/weak_context.fa", new byte[][]{context, chosen, other}, new int[]{2, 1, 4}, false);
				writeFasta(directory+"/weak_two.fa", new byte[][]{context, chosen, other}, new int[]{3, 2, 100}, false);
				writeFasta(directory+"/weak_zero.fa", new byte[][]{context, other}, new int[]{3, 3}, false);
				final byte[] source=Arrays.copyOfRange(query, 73, 201);
				writeFasta(directory+"/weak_source.fa", new byte[][]{context, chosen, other, source}, new int[]{3, 1, 4, 8}, false);
				final byte[] later=Arrays.copyOfRange(query, 10, 202);
				final byte[] laterChoice=Arrays.copyOfRange(query, 11, 203), laterOther=laterChoice.clone();
				laterOther[191]=laterOther[191]=='A' ? (byte)'C' : (byte)'A';
				writeFasta(directory+"/weak_later.fa", new byte[][]{context, chosen, other, later, laterChoice, laterOther},
						new int[]{3, 1, 3, 3, 1, 4}, false);
				writeFasta(directory+"/weak_then_zero.fa", new byte[][]{context, chosen, other, later, laterOther},
						new int[]{3, 1, 3, 3, 3}, false);
				AminoAcid.reverseComplementBasesInPlace(context);
				AminoAcid.reverseComplementBasesInPlace(chosen);
				AminoAcid.reverseComplementBasesInPlace(other);
				writeFasta(directory+"/weak_reverse.fa", new byte[][]{context, chosen, other}, new int[]{3, 1, 4}, false);
			}
		}
	}

	/** Isolates each Hollow160 next-base choice using exact-span evidence records.
	 * An interior mutation in the context confirms that this is Hollow, not Solid,
	 * evidence. Unsupported prefixes also exercise recovery before a later veto. */
	private static void writeNominalFixtures(final String directory){
		assert(directory!=null) : "Nominal fixtures share the fresh directory supplied to the test launcher.";
		final Random random=new Random(971);
		final byte[] query=new byte[240];
		for(int i=0; i<query.length; i++){query[i]=AminoAcid.numberToBase[random.nextInt(4)];}
		writeFasta(directory+"/nominal_query.fa", new byte[][]{query}, new int[]{1}, true);
		writeFasta(directory+"/nominal_seed.fa", new byte[][]{Arrays.copyOf(query, 180)}, new int[]{1}, true);
		writeFasta(directory+"/nominal_seed160.fa", new byte[][]{Arrays.copyOf(query, 160)}, new int[]{1}, true);
		for(int end : new int[]{180, 182, 239}){
			final byte[] context=Arrays.copyOfRange(query, end-160, end);
			final byte[] chosen=Arrays.copyOfRange(query, end-159, end+1), other=chosen.clone();
			context[80]=context[80]=='A' ? (byte)'C' : (byte)'A';
			other[159]=other[159]=='A' ? (byte)'C' : (byte)'A';
			writeFasta(directory+"/nominal_"+end+".fa", new byte[][]{context, other}, new int[]{3, 3}, false);
			final byte[] sourceChosen=Arrays.copyOfRange(query, end-127, end+1);
			writeFasta(directory+"/source_"+end+".fa", new byte[][]{context, other, sourceChosen}, new int[]{3, 3, 2}, false);
			writeFasta(directory+"/source_equal_"+end+".fa", new byte[][]{context, other, sourceChosen}, new int[]{3, 3, 3}, false);
			if(end==180){
				writeFasta(directory+"/source_stronger.fa", new byte[][]{context, other, sourceChosen}, new int[]{3, 3, 4}, false);
				writeFasta(directory+"/source_equal.fa", new byte[][]{context, other, sourceChosen}, new int[]{3, 3, 3}, false);
				writeFasta(directory+"/nominal_supported.fa", new byte[][]{context, other, chosen}, new int[]{3, 30, 1}, false);
				writeFasta(directory+"/nominal_weak_context.fa", new byte[][]{context, other}, new int[]{2, 3}, false);
				writeFasta(directory+"/nominal_weak_other.fa", new byte[][]{context, other}, new int[]{3, 2}, false);
				writeFasta(directory+"/nominal_deadend.fa", new byte[][]{context}, new int[]{3}, false);
				final byte[] later=Arrays.copyOfRange(query, 22, 182), laterOther=Arrays.copyOfRange(query, 23, 183);
				laterOther[159]=laterOther[159]=='A' ? (byte)'C' : (byte)'A';
				writeFasta(directory+"/nominal_restore.fa", new byte[][]{context, other, later, laterOther}, new int[]{3, 2, 3, 3}, false);
				final byte[] laterSource=Arrays.copyOfRange(query, 55, 183);
				writeFasta(directory+"/source_later.fa", new byte[][]{context, other, sourceChosen, later, laterOther, laterSource},
						new int[]{3, 3, 4, 3, 3, 2}, false);
				// Separate the later evidence far enough that it cannot add source counts at the first tie.
				final byte[] farContext=Arrays.copyOfRange(query, 79, 239), farOther=Arrays.copyOfRange(query, 80, 240);
				final byte[] farSource=Arrays.copyOfRange(query, 112, 240);
				farOther[159]=farOther[159]=='A' ? (byte)'C' : (byte)'A';
				writeFasta(directory+"/source_equal_later.fa", new byte[][]{context, other, sourceChosen, farContext, farOther, farSource},
						new int[]{3, 3, 3, 3, 3, 2}, false);
				AminoAcid.reverseComplementBasesInPlace(context);
				AminoAcid.reverseComplementBasesInPlace(other);
				writeFasta(directory+"/nominal_reverse.fa", new byte[][]{context, other}, new int[]{3, 3}, false);
				AminoAcid.reverseComplementBasesInPlace(sourceChosen);
				writeFasta(directory+"/source_reverse.fa", new byte[][]{context, other, sourceChosen}, new int[]{3, 3, 2}, false);
				writeFasta(directory+"/source_equal_reverse.fa", new byte[][]{context, other, sourceChosen}, new int[]{3, 3, 3}, false);
			}
		}
		writeNominalRiseFixtures(directory, query);
		writeRiseSourceFixtures(directory, query);
	}

	/** Separates H160 rises from H128 source support with ignored-middle mutations.
	 * Position40 lies inside H128's left arm but outside H160's two sampled arms.
	 * Exact128bp records then control source chosen counts without changing H160. */
	private static void writeRiseSourceFixtures(final String directory, final byte[] query){
		assert(query.length==240) : "First/later/final source-choice fixtures share the nominal240bp query.";
		for(int end : new int[]{180, 182, 239}){
			final byte[] context=Arrays.copyOfRange(query, end-160, end);
			final byte[] chosen=Arrays.copyOfRange(query, end-159, end+1);
			context[40]=context[40]=='A' ? (byte)'C' : (byte)'A';
			chosen[40]=chosen[40]=='A' ? (byte)'C' : (byte)'A';
			final byte[] source=Arrays.copyOfRange(query, end-127, end+1);
			writeFasta(directory+"/rise_source_"+end+".fa", new byte[][]{context, chosen, source}, new int[]{3, 48, 2}, false);
			if(end==180){
				for(int count : new int[]{0, 1, 3, 4}){
					writeFasta(directory+"/rise_source_count"+count+".fa", new byte[][]{context, chosen, source}, new int[]{3, 48, count}, false);
				}
				final byte[] laterContext=Arrays.copyOfRange(query, 22, 182);
				final byte[] laterChoice=Arrays.copyOfRange(query, 23, 183);
				laterContext[40]=laterContext[40]=='A' ? (byte)'C' : (byte)'A';
				laterChoice[40]=laterChoice[40]=='A' ? (byte)'C' : (byte)'A';
				final byte[] laterSource=Arrays.copyOfRange(query, 55, 183);
				writeFasta(directory+"/rise_source_later.fa", new byte[][]{context, chosen, source, laterContext, laterChoice, laterSource},
						new int[]{3, 48, 3, 3, 48, 2}, false);
				final byte[] other=laterChoice.clone();
				other[159]=other[159]=='A' ? (byte)'C' : (byte)'A';
				writeFasta(directory+"/rise_source_then_conflict.fa", new byte[][]{context, chosen, source, laterContext, other, laterSource},
						new int[]{3, 48, 3, 3, 3, 2}, false);
				AminoAcid.reverseComplementBasesInPlace(context);
				AminoAcid.reverseComplementBasesInPlace(chosen);
				AminoAcid.reverseComplementBasesInPlace(source);
				writeFasta(directory+"/rise_source_reverse.fa", new byte[][]{context, chosen, source}, new int[]{3, 48, 2}, false);
			}
		}
	}

	/** Isolates rises from zero-choice contradictions and checks continuation after abstention.
	 * Evidence records each span exactly160 bases, so every copy contributes one key. */
	private static void writeNominalRiseFixtures(final String directory, final byte[] query){
		assert(query.length==240) : "The CLI fixtures share a180bp seed and test first, later and final saved choices.";
		for(int end : new int[]{180, 182, 239}){
			final byte[] context=Arrays.copyOfRange(query, end-160, end);
			final byte[] chosen=Arrays.copyOfRange(query, end-159, end+1);
			context[80]=context[80]=='A' ? (byte)'C' : (byte)'A';
			writeFasta(directory+"/nominal_rise_"+end+".fa", new byte[][]{context, chosen}, new int[]{3, 48}, false);
			if(end==180){
				writeFasta(directory+"/nominal_rise_below.fa", new byte[][]{context, chosen}, new int[]{3, 47}, false);
				writeFasta(directory+"/nominal_rise_weak.fa", new byte[][]{context, chosen}, new int[]{2, 48}, false);
				writeFasta(directory+"/nominal_rise_absent.fa", new byte[][]{chosen}, new int[]{48}, false);
				writeFasta(directory+"/nominal_rise_floor.fa", new byte[][]{context, chosen}, new int[]{3, 31}, false);
				writeFasta(directory+"/nominal_rise_minimum.fa", new byte[][]{context, chosen}, new int[]{4, 32}, false);
				writeFasta(directory+"/nominal_rise_flat.fa", new byte[][]{context, chosen}, new int[]{48, 48}, false);
				final byte[] other=chosen.clone();
				other[159]=other[159]=='A' ? (byte)'C' : (byte)'A';
				final byte[] source=Arrays.copyOfRange(query, 53, 181);
				final byte[] laterContext=Arrays.copyOfRange(query, 22, 182);
				final byte[] laterChoice=Arrays.copyOfRange(query, 23, 183);
				writeFasta(directory+"/nominal_rise_after_abstain.fa", new byte[][]{context, other, source, laterContext, laterChoice},
						new int[]{3, 3, 4, 3, 48}, false);
				AminoAcid.reverseComplementBasesInPlace(context);
				AminoAcid.reverseComplementBasesInPlace(chosen);
				writeFasta(directory+"/nominal_rise_reverse.fa", new byte[][]{context, chosen}, new int[]{3, 48}, false);
			}
		}
		writeFasta(directory+"/nominal_rise_batch.fa", new byte[][]{query, Arrays.copyOf(query, 180)}, new int[]{1, 1}, false);
	}

	/** Makes positive Hollow next counts without any matching Solid64 next key.
	 * The changed interior base is ignored by Hollow128 but lies in the Solid64
	 * suffix. Exact-length evidence isolates the intended stopping position. */
	private static void writeShortDropFixtures(final String directory){
		assert(directory!=null) : "All diagnostic fixtures use the supplied fresh directory.";
		final byte[] query=new byte[200];
		final Random random=new Random(965);
		for(int i=0; i<query.length; i++){query[i]=AminoAcid.numberToBase[random.nextInt(4)];}
		writeFasta(directory+"/shortdrop_query.fa", new byte[][]{query}, new int[]{1}, true);
		writeFasta(directory+"/shortdrop_seed.fa", new byte[][]{Arrays.copyOf(query, 130)}, new int[]{1}, true);
		for(int end : new int[]{130, 132, 199}){
			final byte[] context=Arrays.copyOfRange(query, end-128, end);
			final byte[] next=Arrays.copyOfRange(query, end-127, end+1), twin=next.clone();
			twin[80]=twin[80]=='A' ? (byte)'C' : (byte)'A';
			final byte[] other=Arrays.copyOfRange(next, 64, 128);
			other[63]=other[63]=='A' ? (byte)'C' : (byte)'A';
			writeFasta(directory+"/shortdrop_"+end+".fa", new byte[][]{context, twin, other}, new int[]{32, 8, 2}, false);
			if(end==130){
				writeFasta(directory+"/shortdrop_flat.fa", new byte[][]{context, twin, other}, new int[]{32, 32, 2}, false);
				writeFasta(directory+"/shortdrop_supported.fa", new byte[][]{context, twin, other, next}, new int[]{64, 8, 2, 1}, false);
				writeFasta(directory+"/shortdrop_weak_other.fa", new byte[][]{context, twin, other}, new int[]{32, 8, 1}, false);
				AminoAcid.reverseComplementBasesInPlace(context);
				AminoAcid.reverseComplementBasesInPlace(twin);
				AminoAcid.reverseComplementBasesInPlace(other);
				writeFasta(directory+"/shortdrop_reverse_evidence.fa", new byte[][]{context, twin, other}, new int[]{32, 8, 2}, false);
			}
		}
	}

	/** Isolates short forks from full-word recovery with exact-length count records. */
	private static void writeSeedForkFixtures(final String directory){
		assert(directory!=null) : "Fork fixtures must live beside the candidate and output-preservation fixtures.";
		final byte[] query=new byte[200];
		final Random random=new Random(953);
		for(int i=0; i<query.length; i++){query[i]=AminoAcid.numberToBase[random.nextInt(4)];}
		final byte[] seedWord=Arrays.copyOfRange(query, 2, 130);
		writeFasta(directory+"/fork_query.fa", new byte[][]{query}, new int[]{1}, true);
		writeFasta(directory+"/fork_seed.fa", new byte[][]{Arrays.copyOf(query, 130)}, new int[]{1}, true);
		writeFasta(directory+"/fork_absent.fa", new byte[][]{seedWord}, new int[]{1}, false);
		for(int end : new int[]{130, 132, 193, 194}){
			final byte[] context=Arrays.copyOfRange(query, end-64, end);
			final byte[] chosen=Arrays.copyOfRange(query, end-63, end+1), other=chosen.clone();
			other[63]=chosen[63]=='A' ? (byte)'C' : (byte)'A';
			writeFasta(directory+"/fork_"+end+".fa", new byte[][]{seedWord, context, chosen, other}, new int[]{1, 3, 6, 3}, false);
			if(end==132){
				final byte[] recovered=Arrays.copyOfRange(query, 3, 131);
				writeFasta(directory+"/fork_recovery.fa", new byte[][]{seedWord, context, chosen, other, recovered}, new int[]{1, 3, 6, 3, 2}, false);
				writeFasta(directory+"/fork_strong_seed.fa", new byte[][]{seedWord, context, chosen, other}, new int[]{2, 3, 6, 3}, false);
			}
		}
		final byte[] qualities=new byte[query.length];
		Arrays.fill(qualities, (byte)'I');
		final ByteStreamWriter writer=new ByteStreamWriter(directory+"/fork_query.fq", false, false, false);
		writer.start();
		try{writer.print(new ByteBuilder().append("@seed\n").append(query).append("\n+\n").append(qualities).nl());}
		finally{
			if(writer.poisonAndWait()){throw new RuntimeException("Fork quality fixture output failed.");}
		}
	}

	/** Adds short-only evidence without increasing the Hollow128 counts. */
	private static void writeShortFixtures(final String directory, final byte[] context, final byte[] next, final byte[] second){
		assert(context.length==128 && next.length==128 && second.length==128) :
				"Long fixture records each count one Hollow key; added64bp records affect only the Solid table.";
		final byte[] chosen=Arrays.copyOfRange(next, 64, 128), other=chosen.clone();
		other[63]=chosen[63]=='A' ? (byte)'C' : (byte)'A';
		final byte[] chosen2=Arrays.copyOfRange(second, 64, 128), other2=chosen2.clone();
		other2[63]=chosen2[63]=='A' ? (byte)'C' : (byte)'A';
		writeFasta(directory+"/short_conflict.fa", new byte[][]{context, next, chosen, other}, new int[]{2, 2, 3, 10}, false);
		writeFasta(directory+"/short_below.fa", new byte[][]{context, next, chosen, other}, new int[]{2, 2, 3, 9}, false);
		writeFasta(directory+"/short_strong_hollow.fa", new byte[][]{context, next, chosen, other}, new int[]{2, 3, 3, 100}, false);
		writeFasta(directory+"/short_weak_solid.fa", new byte[][]{context, next, other}, new int[]{2, 2, 100}, false);
		writeFasta(directory+"/short_zero_hollow.fa", new byte[][]{context, chosen, other}, new int[]{2, 3, 6}, false);
		writeFasta(directory+"/short_no_choice.fa", new byte[][]{context, other}, new int[]{2, 100}, false);
		writeFasta(directory+"/short_later.fa", new byte[][]{context, next, second, chosen2, other2}, new int[]{2, 3, 2, 3, 10}, false);
		writeFasta(directory+"/short_restore.fa", new byte[][]{context, next, chosen, other, second, chosen2, other2}, new int[]{2, 2, 3, 9, 2, 3, 10}, false);
	}

	/** Writes native-schema candidates, including malformed joins and unrelated stop positions. */
	private static void writeCandidates(final String directory){
		assert(directory!=null) : "Candidate fixtures share the fresh sequence-fixture directory.";
		candidate(directory, "seed", "seed\t130\t350\t130\t128\t128\t11\t1\t10\tA\t-1\t-1\tbackward_dominance", 1);
		candidate(directory, "later", "seed\t130\t350\t131\t128\t128\t11\t1\t10\tA\t-1\t-1\tbackward_dominance", 1);
		candidate(directory, "other", "seed\t130\t350\t130\t128\t129\t11\t1\t10\tA\t-1\t-1\talternate_span", 1);
		candidate(directory, "none", "seed\t130\t350\t350\t128\t-1\t-1\t-1\t-1\t.\t-1\t-1\tnone", 1);
		candidate(directory, "bad_id", "wrong\t130\t350\t130\t128\t128\t11\t1\t10\tA\t-1\t-1\tbackward_dominance", 1);
		candidate(directory, "bad_seed", "seed\t129\t350\t130\t128\t128\t11\t1\t10\tA\t-1\t-1\tbackward_dominance", 1);
		candidate(directory, "bad_k", "seed\t130\t350\t130\t127\t128\t11\t1\t10\tA\t-1\t-1\tbackward_dominance", 1);
		candidate(directory, "bad_span", "seed\t130\t350\t130\t128\t129\t11\t1\t10\tA\t-1\t-1\tbackward_dominance", 1);
		candidate(directory, "bad_bounds", "seed\t130\t132\t130\t128\t128\t11\t1\t10\tA\t-1\t-1\tbackward_dominance", 1);
		candidate(directory, "bad_count", "seed\t130\t350\t130\t128\t128\t11\t-2\t10\tA\t-1\t-1\tbackward_dominance", 1);
		candidate(directory, "bad_reason", "seed\t130\t350\t130\t128\t128\t11\t1\t10\tA\t-1\t-1\tunknown", 1);
		candidate(directory, "extra", "seed\t130\t350\t130\t128\t128\t11\t1\t10\tA\t-1\t-1\tbackward_dominance", 2);
		candidate(directory, "empty", "", 0);
	}

	/** Uses native output for tiny tabular fixtures; no shell-generated scientific evidence. */
	private static void candidate(final String directory, final String name, final String row, final int copies){
		assert(copies>=0 && copies<=2) : "Zero/two copies exercise missing/extra rows; ordinary fixtures have one.";
		final ByteStreamWriter writer=new ByteStreamWriter(directory+"/candidate_"+name+".tsv", false, false, false);
		writer.start();
		try{
			writer.print("id\tseed_length\toriginal_length\tretained_length\tk\tquery_span\tcontext_depth\tchosen_count\tother_count\tother_base\tself_count\texternal_count\tstop_reason\n");
			for(int copy=0; copy<copies; copy++){writer.print(row+"\n");}
		}finally{
			if(writer.poisonAndWait()){throw new RuntimeException("Candidate fixture output failed: "+name);}
		}
	}

	/** Uses native output and distinct evidence IDs, with a shared seed/query ID for prefix comparisons. */
	private static void writeFasta(final String path, final byte[][] sequences, final int[] copies, final boolean seed){
		assert(sequences.length==copies.length && (!seed || sequences.length==1)) :
				"Fixture copy counts correspond one-to-one to words; seed/query files each contain one record.";
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		final ByteBuilder row=new ByteBuilder();
		writer.start();
		try{
			for(int i=0; i<sequences.length; i++){
				for(int copy=0; copy<copies[i]; copy++){
					row.clear().append('>');
					if(seed){row.append("seed");}
					else{row.append("evidence_").append(i).append('_').append(copy);}
					row.nl().append(sequences[i]).nl();
					writer.print(row);
				}
			}
		}finally{
			if(writer.poisonAndWait()){throw new RuntimeException("Fixture output failed: "+path);}
		}
	}

	/** The independent oracle uses original query endpoints and the requested table-span salt. */
	private static void exercise(final byte[] bases, final int k){
		assert(bases.length>=66 && k>=66) : "These fixtures exercise nonoverlapping anchor and neighbor contexts.";
		final Kmer word=new Kmer(k, true);
		// Nonzero offsets catch accidental absolute indexing of the observed query.
		final byte[] padded=new byte[bases.length+9];
		System.arraycopy(bases, 0, padded, 4, bases.length);
		KmerWindowDepth.setHollowQuery(word, padded, 4, bases.length);
		checkHash(word, bases, k);
		final byte[] reverse=bases.clone();
		AminoAcid.reverseComplementBasesInPlace(reverse);
		final Kmer rc=new Kmer(k, true);
		KmerWindowDepth.setHollowQuery(rc, reverse, 0, reverse.length);
		check(word.xor()==rc.xor() && word.xor2()==rc.xor2(), "Reverse-complement query keys differ.");
		for(int base=0; base<4; base++){
			final Kmer right=word.clone(), left=word.clone();
			right.addRightNumeric(base); left.addLeftNumeric(base);
			final byte[] shiftedRight=new byte[bases.length], shiftedLeft=new byte[bases.length];
			System.arraycopy(bases, 1, shiftedRight, 0, bases.length-1);
			shiftedRight[bases.length-1]=AminoAcid.numberToBase[base];
			System.arraycopy(bases, 0, shiftedLeft, 1, bases.length-1);
			shiftedLeft[0]=AminoAcid.numberToBase[base];
			checkHash(right, shiftedRight, k);
			checkHash(left, shiftedLeft, k);
		}
	}

	/** Encodes endpoints directly, without using the synthetic-word construction. */
	private static void checkHash(final Kmer word, final byte[] query, final int k){
		assert(query.length>=64) : "The independent oracle reads two32bp endpoints.";
		final byte[] rc=query.clone();
		AminoAcid.reverseComplementBasesInPlace(rc);
		final long f1=arm(query, 0), f2=arm(query, query.length-32);
		final long r1=arm(rc, 0), r2=arm(rc, rc.length-32);
		final int comparison=Long.compareUnsigned(f1, r1);
		final boolean forward=comparison>0 || (comparison==0 && Long.compareUnsigned(f2, r2)>=0);
		final long a=forward ? f1 : r1, b=forward ? f2 : r2;
		final long first=(Long.rotateLeft(Tools.hash64shift(a), 25)^Tools.hash64shift(b)^Tools.hash64shift(k))&Long.MAX_VALUE;
		final long second=Tools.splitMix64(Long.rotateRight(Tools.splitMix64(a), 23)^b^Tools.splitMix64(k))&Long.MAX_VALUE;
		check(word.xor()==first && word.xor2()==second, "Observed anchors or one-step neighbor context changed at table K="+k+", queryspan="+query.length);
	}

	/** Packs one anchor base by base independently of Kmer's packed-word boundaries. */
	private static long arm(final byte[] bases, final int start){
		assert(start>=0 && start<=bases.length-32) : "An oracle anchor must fit the original query.";
		long value=0;
		for(int i=start; i<start+32; i++){value=(value<<2)|AminoAcid.baseToNumber[bases[i]];}
		return value;
	}

	/** Keeps test failures fatal even when the caller disables assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
}
