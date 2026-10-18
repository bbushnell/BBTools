package test.prot;

import prot.AAGraph;
import prot.HbmCompetitiveGate;
import prot.HbmPositionModel;
import prot.HbmTailGate;

import java.util.Arrays;

/** Independent literal paths exercise paired identity and terminal evidence boundaries. @author Keqing */
public final class HbmTailProfileTest{
	public static void main(String[] args){
		final HbmTailGate.Metrics m=new HbmTailGate.Metrics();
		final double[] background=new double[20]; Arrays.fill(background, 0.05);
		final byte[] ref=new byte[50], query=new byte[50];
		final HbmCompetitiveGate gate=gate(ref, background);
		final HbmPositionModel.TailResult full=path(50, 50, ops(50, 'm'));
		check(gate.evaluateTail(query, 0, full, HbmTailGate.PROFILE_JOINT, m)==HbmTailGate.SUPPORTED
			&& m.terminalScore64>0, "A conserved terminal window must retain positive profile evidence");
		Arrays.fill(query, 30, 50, (byte)1);
		check(gate.evaluateTail(query, 0, full, HbmTailGate.EXPERIMENTAL, m)==HbmTailGate.SUPPORTED,
			"Geometry baseline accepts this sixty-percent-identity literal path");
		check(gate.evaluateTail(query, 0, full, HbmTailGate.PROFILE_JOINT, m)==HbmTailGate.UNSUPPORTED_END
			&& m.terminalScore64<0, "Negative terminal score must reject despite adequate whole-path identity and geometry");

		// The dominant A background rounds an A-only tail to zero; rare R prefix still carries signal.
		final double[] skewed=new double[20]; Arrays.fill(skewed, .001/19); skewed[0]=.999;
		final byte[] zeroEnd=new byte[50]; Arrays.fill(zeroEnd, 0, 30, (byte)1);
		final HbmCompetitiveGate zeroGate=gate(zeroEnd, skewed);
		check(zeroGate.tailWindowScore(zeroEnd, 0, full, 20)==0, "Literal background fixture must have exactly zero terminal score");
		check(zeroGate.evaluateTail(zeroEnd, 0, full, HbmTailGate.PROFILE_JOINT, m)==HbmTailGate.UNSUPPORTED_END,
			"Candidate requires strictly positive, not nonnegative, terminal evidence");

		final byte[] gappedQuery=new byte[50], gappedRef=new byte[60]; Arrays.fill(gappedQuery, 0, 29, (byte)1);
		final HbmPositionModel.TailResult gaps=path(50, 60, concat(ops(40, 'm'), ops(10, 'D'), ops(10, 'm')));
		check(HbmTailGate.evaluate(gappedQuery, gappedRef, gaps, 0, 59, HbmTailGate.EXPERIMENTAL, m)==HbmTailGate.LOW_IDENTITY
			&& m.identity==35f, "Geometry identity includes ten deletion columns:21/60");
		check(HbmTailGate.evaluate(gappedQuery, gappedRef, gaps, 0, 59, HbmTailGate.PROFILE_JOINT, m)==HbmTailGate.SUPPORTED
			&& m.identity==42f && m.terminalPaired==10, "Profile geometry uses21/50 paired identity and preserves inclusive joint-ten pairing");
		final HbmPositionModel.TailResult shifted=path(69, 82, concat(ops(60, 'm'), ops(11, 'D'), ops(9, 'm'), ops(2, 'D')));
		check(HbmTailGate.evaluate(new byte[69], new byte[82], shifted, 0, 81, HbmTailGate.PROFILE_JOINT, m)==HbmTailGate.UNSUPPORTED_END
			&& m.queryTerminalPaired==20 && m.terminalPaired==9, "Query-terminal pairing cannot substitute for joint stop-proximal pairing");
		check(gate.evaluateTail(query, 0, null, HbmTailGate.PROFILE_JOINT, m)==HbmTailGate.NO_ALIGNMENT
			&& m.terminalScore64==Integer.MIN_VALUE && m.queryTerminalPaired==0, "No-hit clears prior profile evidence");
		check(HbmTailGate.forName("geometry")==HbmTailGate.EXPERIMENTAL && HbmTailGate.forName("profile")==HbmTailGate.PROFILE_JOINT,
			"CLI policy names must select the intended immutable settings");
		boolean rejected=false;
		try{HbmTailGate.forName("unknown");}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Unknown policy must fail before caller output");
		System.err.println("HBM_TAIL_PROFILE_TEST_PASS positive_negative_zero=true paired_identity=true joint_pairing=true no_hit_reset=true");
	}
	private static HbmCompetitiveGate gate(byte[] reference, double[] background){
		final HbmPositionModel profile=HbmPositionModel.deriveOne(new AAGraph(reference, 0), "logodds", .01, true, background, -4);
		return new HbmCompetitiveGate(new String[]{"fixture"}, new byte[][]{reference}, new HbmPositionModel[]{profile}, new int[]{0}, new int[]{reference.length-1});
	}
	private static HbmPositionModel.TailResult path(int q, int r, byte[] ops){return new HbmPositionModel.TailResult(1, 0, q-1, 0, r-1, ops);}
	private static byte[] ops(int n, char op){final byte[] a=new byte[n]; Arrays.fill(a, (byte)op); return a;}
	private static byte[] concat(byte[]... parts){
		int n=0; for(byte[] part : parts){n+=part.length;}
		final byte[] all=new byte[n]; int at=0;
		for(byte[] part : parts){System.arraycopy(part, 0, all, at, part.length); at+=part.length;}
		return all;
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private HbmTailProfileTest(){}
}
