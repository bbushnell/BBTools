package test.prot;

import prot.FamilyShortlistScores;
import prot.HbmPositionModel;
import prot.HbmTailGate;

import java.util.Arrays;

/** Literal tail paths distinguish paired stop evidence from mere endpoint consumption. */
public final class HbmTailGateTest{
	public static void main(String[] args){
		try{test();}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}
	private static void test(){
		caps();
		final HbmTailGate.Metrics m=new HbmTailGate.Metrics();
		final byte[] query=new byte[120], ref=new byte[80], path=ops(80, 'm');
		check(evaluate(query, ref, 40, 0, path, 0, 79, m)==HbmTailGate.SUPPORTED, "Perfect short protein must survive a freely clipped query prefix");
		check(m.paired==80 && m.pairedCore==80 && m.terminalPaired==20 && m.querySuffix==0 && m.referenceSuffix==0,
			"Literal eighty-column perfect tail has exact paired/terminal counts");
		check(evaluate(new byte[511], new byte[1200], 11, 700, ops(500, 'm'), 0, 1199, m)==HbmTailGate.SUPPORTED
			&& m.coreLength==500 && m.coreCoverage==1, "Long reference coverage denominator must be capped to its tail");
		check(evaluate(new byte[111], new byte[100], 0, 0, concat(ops(100, 'm'), ops(11, 'I')), 0, 99, m)==HbmTailGate.UNSUPPORTED_END
			&& m.querySuffix==11 && m.referenceSuffix==0, "Long inserted query suffix cannot prove its stop");
		check(evaluate(new byte[100], new byte[111], 0, 0, concat(ops(100, 'm'), ops(11, 'D')), 0, 110, m)==HbmTailGate.UNSUPPORTED_END
			&& m.querySuffix==0 && m.referenceSuffix==11, "Internal model match cannot prove terminality by deleting the remaining tail");
		check(evaluate(new byte[110], new byte[100], 0, 0, concat(ops(100, 'm'), ops(10, 'I')), 0, 99, m)==HbmTailGate.SUPPORTED
			&& m.terminalPaired==10, "Experimental ten-residue suffix and ten terminal pairs are inclusive boundaries");
		check(evaluate(new byte[29], new byte[29], 0, 0, ops(29, 'm'), 0, 28, m)==HbmTailGate.SHORT_SUPPORT,
			"Twenty-nine perfect residues fail the explicit thirty-residue experimental floor");
		check(evaluate(new byte[30], new byte[30], 0, 0, ops(30, 'm'), 0, 29, m)==HbmTailGate.SUPPORTED,
			"Thirty perfect residues meet the experimental support floor");
		final byte[] lowIdentity=new byte[100]; Arrays.fill(lowIdentity, 40, 100, (byte)1);
		check(evaluate(lowIdentity, new byte[100], 0, 0, ops(100, 'm'), 0, 99, m)==HbmTailGate.LOW_IDENTITY
			&& m.identity==40, "Identity is computed from literal matches, not profile score");
		check(evaluate(new byte[79], new byte[100], 0, 21, ops(79, 'm'), 0, 99, m)==HbmTailGate.LOW_COVERAGE,
			"Free model-prefix clipping does not erase the capped core coverage requirement");
		check(evaluate(new byte[80], new byte[100], 0, 20, ops(80, 'm'), 0, 99, m)==HbmTailGate.SUPPORTED,
			"Eighty percent capped core coverage is inclusive");
		check(evaluate(new byte[100], new byte[100], 0, 0, ops(100, 'm'), 0, 19, m)==HbmTailGate.SHORT_SUPPORT,
			"Only twenty paired core columns are insufficient despite a long total path");
		check(HbmTailGate.evaluate(query, ref, null, 0, 79, HbmTailGate.EXPERIMENTAL, m)==HbmTailGate.NO_ALIGNMENT
			&& m.paired==0 && Float.isNaN(m.identity), "No-hit evaluation must clear earlier winner measurements");
		boolean rejected=false;
		try{evaluate(query, ref, 40, 0, ops(79, 'm'), 0, 79, m);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A malformed trace that fails to consume both ends must fail loudly");
		final int[][] hits={{4,4}};
		final float[][] composition={{0.5f,0.5f}}, scores={{0,0}};
		final double[][] lengths={{500,1200}}, spreads={{0.05,0.05}};
		FamilyShortlistScores.score(hits, composition, lengths, spreads, 500, scores);
		check(scores[0][0]>scores[0][1], "Ordinary shortlist retains its full-protein length penalty");
		FamilyShortlistScores.scoreTail(hits, composition, lengths, spreads, 500, scores);
		check(scores[0][0]==scores[0][1], "Identical tail evidence must not penalize the longer full protein");
		final float[][] splitScores={{0},{0}};
		FamilyShortlistScores.scoreTail(new int[][]{{4},{4}}, new float[][]{{0.5f},{0.5f}},
			new double[][]{{500},{1200}}, new double[][]{{0.05},{0.05}}, 500, splitScores);
		check(scores[0][0]==splitScores[0][0] && scores[0][1]==splitScores[1][0], "Tail normalization is invariant to pack partition");
		System.err.println("HBM_TAIL_GATE_TEST_PASS literal_boundaries=true capped_core=true terminal_gaps=true identity=true no_hit_reset=true");
	}
	/** Window changes must not silently recalibrate either named evidence policy. */
	private static void caps(){
		check(HbmTailGate.forName("geometry")==HbmTailGate.EXPERIMENTAL
			&& HbmTailGate.forName("profile")==HbmTailGate.PROFILE_JOINT, "Legacy names retain default singleton policies");
		for(String name : new String[]{"geometry", "PROFILE"}){
			final HbmTailGate.Policy base=HbmTailGate.forName(name);
			check(HbmTailGate.forName(name, 500)==base, "Explicit default cap preserves the named policy object");
			for(int cap : new int[]{30, 200, 300, 400, 500, 600}){
				final HbmTailGate.Policy policy=HbmTailGate.forName(name, cap);
				check(policy.maxLength==cap && policy.minPaired==30 && policy.maxTerminalGap==10
					&& policy.terminalWindow==20 && policy.minTerminalPaired==10
					&& Float.floatToIntBits(policy.minIdentity)==Float.floatToIntBits(40.612846f)
					&& Float.floatToIntBits(policy.minCoverage)==Float.floatToIntBits(0.8f)
					&& policy.profileEvidence==base.profileEvidence,
					"Only maxLength may change; literal geometry thresholds and profile evidence flag are fixed");
				final byte[] query=new byte[cap+1], reference=new byte[cap+7];
				final HbmTailGate.Metrics metrics=new HbmTailGate.Metrics();
				final HbmPositionModel.TailResult exact=new HbmPositionModel.TailResult(100, 1, cap, 7, cap+6, ops(cap, 'm'));
				check(HbmTailGate.evaluate(query, reference, exact, 0, cap+6, policy, metrics)==HbmTailGate.SUPPORTED
					&& metrics.paired==cap && metrics.coreLength==cap && metrics.coreCoverage==1,
					"Inclusive cap boundary must preserve exact paired/core counts for both policies");
				boolean rejected=false;
				try{
					HbmTailGate.evaluate(query, reference,
						new HbmPositionModel.TailResult(100, 0, cap, 6, cap+6, ops(cap+1, 'm')), 0, cap+6, policy, metrics);
				}catch(IllegalArgumentException expected){rejected=true;}
				check(rejected, "One residue beyond either capped window must fail rather than change its denominator");
			}
			for(int cap : new int[]{-1, 0, 29}){
				boolean rejected=false;
				try{HbmTailGate.forName(name, cap);}catch(IllegalArgumentException expected){rejected=true;}
				check(rejected, "A cap below the unchanged thirty-pair minimum cannot form a valid policy");
			}
		}
		boolean rejected=false;
		try{HbmTailGate.forName("unknown", 200);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Custom caps must not bypass named-policy validation");
		System.err.println("HBM_TAIL_CAP_TEST_PASS caps_200_600=true default_identity=true profile_preserved=true boundaries=true");
	}
	private static int evaluate(byte[] q, byte[] r, int qs, int rs, byte[] path, int first, int last, HbmTailGate.Metrics m){
		final HbmPositionModel.TailResult a=new HbmPositionModel.TailResult(100, qs, q.length-1, rs, r.length-1, path);
		return HbmTailGate.evaluate(q, r, a, first, last, HbmTailGate.EXPERIMENTAL, m);
	}
	private static byte[] ops(int length, char op){final byte[] out=new byte[length]; Arrays.fill(out, (byte)op); return out;}
	private static byte[] concat(byte[] a, byte[] b){final byte[] out=Arrays.copyOf(a, a.length+b.length); System.arraycopy(b, 0, out, a.length, b.length); return out;}
	private static void check(boolean value, String reason){if(!value){throw new AssertionError(reason);}}
}
