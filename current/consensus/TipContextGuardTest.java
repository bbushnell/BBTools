package consensus;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Random;

/** Exact evidence-count and overlap-geometry tests, separate from biological safety. @author Fischl */
public final class TipContextGuardTest {

	/** Runs deterministic cases without changing any files or ordinary consensus defaults. */
	public static void main(final String[] args){
		counts();
		geometry();
		kmerEvidence();
		System.out.println("TIP_CONTEXT_GUARD_TEST_PASS checks="+checks);
	}

	/** A reference pseudocount must not turn one real observation into two. */
	private static void counts(){
		final BaseGraph graph=new BaseGraph("target", "ACGT".getBytes(StandardCharsets.US_ASCII), null, 0, 0);
		final BaseNode node=graph.ref[0];
		expect(TipContextGuard.observations(node, (byte)'A')==0, "Reference alone is not read support.");
		node.add((byte)'A', 31);
		expect(TipContextGuard.observations(node, (byte)'A')==1, "Subtract exactly one reference observation.");
		node.add((byte)'A', 31);
		expect(TipContextGuard.observations(node, (byte)'A')==2, "Two real observations must protect the allele.");
		node.add((byte)'C', 31);
		expect(TipContextGuard.observations(node, (byte)'C')==1, "Alternate observations have no pseudocount.");
	}

	/** Tests placement multiplicity, both target sides, donor strands and strict boundaries. */
	private static void geometry(){
		final byte[] target=new byte[160];
		final Random random=new Random(246L);
		for(int i=0; i<target.length; i++){target[i]=(byte)"ACGT".charAt(random.nextInt(4));}
		final byte[] donor=new byte[160], context=new byte[240];
		System.arraycopy(target, 80, donor, 0, 80);
		for(int i=80; i<donor.length; i++){donor[i]=(byte)"ACGT".charAt(random.nextInt(4));}
		final int pos=120;
		donor[pos-80]=target[pos]=='A' ? (byte)'C' : (byte)'A';
		System.arraycopy(target, 0, context, 0, 160);
		System.arraycopy(donor, 80, context, 160, 80);
		final TipContextGuard.CandidateSet set=scan(target, donor);
		expect(set.count==1 && set.only.overlap==80, "The nonrepetitive example has one 80bp placement.");
		expect(TipContextGuard.matchesContext(set.only, context, 0, 160, false), "Added right context must match donor.");
		expect(TipContextGuard.matchesContext(set.only, TipContextGuard.reverse(context), 80, 160, true),
				"Left context in reverse orientation must describe the same physical bases.");
		final byte to=donor[pos-80];
		for(int original=0; original<=3; original++){
			for(int additional=0; additional<=6; additional++){
				final String expected=original>=2 ? "preserve_original_support"
						: additional<2 ? "insufficient_alternate_support" : "accept";
				expect(TipContextGuard.decision(original, additional, pos, to, 160, set, true).equals(expected),
						"Support boundary changed: original="+original+" alternate="+additional);
			}
		}
		expect(TipContextGuard.decision(1, 5, 79, to, 160, set, true).equals("outside_overlap"), "No extrapolation outside overlap.");
		expect(TipContextGuard.decision(1, 5, pos, target[pos], 160, set, true).equals("donor_disagrees"), "Donor must agree with proposed base.");
		expect(TipContextGuard.decision(1, 5, pos, to, 160, set, false).equals("context_not_from_unique_donor"), "Arbitrary context is not evidence.");
		final byte[] wrong=context.clone();
		wrong[200]=wrong[200]=='A' ? (byte)'C' : (byte)'A';
		expect(!TipContextGuard.matchesContext(set.only, wrong, 0, 160, false), "Wrong borrowed context must be rejected.");
		expect(!TipContextGuard.matchesContext(set.only, target, 0, 160, false), "No extension means no additional context.");
		expect(!TipContextGuard.matchesContext(set.only, Arrays.copyOf(context, 241), 0, 160, false), "Do not read beyond donor.");
		final byte[] competing=donor.clone();
		competing[pos-80]=target[pos];
		final TipContextGuard.CandidateSet conflict=scan(target, donor);
		TipContextGuard.scan(target, competing, 60, 100, conflict);
		expect(conflict.count==2, "A clean competing allele is a second placement, not grounds to pick the mismatch.");
		expect(TipContextGuard.decision(1, 5, pos, to, 160, conflict, true).equals("nonunique_or_missing_overlap"), "Conflicting donors veto.");
		TipContextGuard.scan(target, donor, 60, 100, set);
		expect(set.count==2, "Identical repeated copies also count as distinct candidates.");
		expect(TipContextGuard.decision(1, 5, pos, to, 160, set, true).equals("nonunique_or_missing_overlap"), "Duplicate copies veto.");
		final byte[] two=donor.clone();
		two[20]=two[20]=='A' ? (byte)'C' : (byte)'A';
		expect(scan(target, two).count==0, "Two mismatches exceed the bounded experiment.");
		final byte[] ambiguous=donor.clone();
		ambiguous[0]='N';
		expect(scan(target, ambiguous).count==0, "Unknown bases cannot prove an overlap.");
		for(int boundary : new int[]{0, 15, 16, 63, 64, 79}){
			final byte[] test=competing.clone();
			test[boundary]=test[boundary]=='A' ? (byte)'C' : (byte)'A';
			expect(scan(target, test).count==((boundary>=16 && boundary<=63) ? 1 : 0), "16bp mismatch anchors, offset="+boundary);
		}
		// The same one-read evidence cannot be given different decisions by changing its truth label.
		final TipContextGuard.CandidateSet single=scan(target, donor);
		expect(TipContextGuard.decision(1, 6, pos, to, 160, single, true).equals("accept"),
				"One-read true variants remain an explicit adverse counterexample, not guaranteed protected.");
	}

	/** Mirrors the command's two-strand donor search without mutating either input. */
	private static TipContextGuard.CandidateSet scan(final byte[] target, final byte[] donor){
		assert(target!=null && donor!=null) : "Fixture arrays must exist before searching their end overlaps.";
		final TipContextGuard.CandidateSet set=new TipContextGuard.CandidateSet();
		TipContextGuard.scan(target, donor, 60, 100, set);
		TipContextGuard.scan(target, TipContextGuard.reverse(donor), 60, 100, set);
		return set;
	}

	/** Missing queries differ from observed zero counts; either end uses original, not oriented, coordinates. */
	private static void kmerEvidence(){
		final TipContextGuard.KmerEvidence evidence=new TipContextGuard.KmerEvidence(96);
		expect(evidence.support(32)==-1, "Missing independent evidence must abstain.");
		Arrays.fill(evidence.depths, 1);
		expect(evidence.support(32)==1, "Depth-one words remain eligible for the other guards.");
		evidence.depths[63]=2;
		expect(evidence.support(32)==2, "Last overlapping word includes the proposed base and protects it.");
		expect(evidence.support(31)==1, "The next word must not protect a base it does not contain.");
		evidence.depths[32]=-1;
		expect(evidence.support(32)==2, "Known repeated support already suffices to veto, even if another query is missing.");
		evidence.depths[63]=1;
		expect(evidence.support(32)==-1, "Weak support needs every overlapping query before allowing a change.");
		evidence.depths[32]=0;
		expect(evidence.support(32)==1, "An actual zero observation is not a missing query.");
		evidence.depths[31]=4;
		expect(evidence.support(0)==4, "First base has exactly one complete overlapping32mer.");
		evidence.depths[95]=3;
		expect(evidence.support(95)==3, "Last base has exactly one complete overlapping32mer.");
	}

	/** Reports the violated scientific contract even without -ea. */
	private static void expect(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
}
