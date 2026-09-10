package bin;

import java.util.ArrayList;

/**
 * Coverage log-likelihood ratio for bin-pair triage (the covlr gate/heuristic
 * in Oracle).  Answers "could these two depth profiles come from ONE genome?"
 * with a calibrated, length-aware statistic instead of a hard ratio cutoff.
 *
 * Model, per sample s with depths cA, cB and bin sizes LA, LB (bp):
 *   H0 (same genome):  d = cB-cA ~ Normal(0, U),
 *       U = K*mid*(1/LA + 1/LB) + s2*mid^2,  mid = max((cA+cB)/2, 0.05)
 *       (first term: read-sampling noise, shrinks with length;
 *        second term: biological/compositional dispersion, s2 = 2*sigma^2)
 *   H1 (different genomes):  d ~ Laplace(scale max(mid, 1.5))
 *   A sample where both bins are absent (max depth < DET) contributes the
 *   capped co-absence term coabs[s] = min(-ln(absence rate of sample s), 2)
 *   instead: shared absence is weak evidence where most genomes are absent
 *   (e.g. sparse soil libraries) and stronger where few are.
 * The sum over samples is the LR in nats; positive favors same-genome.
 *
 * CITATION: this math is reimplemented from warpbin by the NeLLi team
 * (Frederik Schulz's group, LBNL): github.com/NeLLi-team/warpbin, method.md
 * and binner.py (Model.lr_cov).  The warpbin CODE is under an LBNL
 * non-commercial license and none of it is copied; formulas and default
 * constants (K=350, DET=0.5, sigma=0.12, co-absence cap 2 nats) come from
 * its documentation and were re-derived here in Java.
 *
 * DEPARTURE from warpbin (design: Amber's oracle_lr_triage_design.md,
 * 2026-09-09): warpbin sums per-sample evidence as if samples were
 * independent, which the correlation-trap benchmark showed is catastrophic
 * on correlated libraries (8 nominal = 2 logical x 4 near-duplicates:
 * warpbin 976 -> 108 Total).  Here the statistic is intended to run on
 * CONDENSED logical samples (see autocondense / SampleCondenser); when
 * samples are left uncondensed and appear correlated, DataLoader prints a
 * warning rather than silently over-counting evidence.
 *
 * @author Amber
 * @date September 9, 2026
 */
public class CovLR {

	/**
	 * Coverage log-likelihood ratio between two bins, in nats.
	 * Positive favors H0 (same genome); negative favors H1 (different).
	 * O(numDepths) with two logs per active sample.
	 */
	public static float lrCov(Bin a, Bin b){
		final int n=a.numDepths();
		assert(n==b.numDepths()) : n+", "+b.numDepths();
		assert(coabs!=null && coabs.length>=n) :
			"CovLR.calibrate was not run for "+n+" samples (see DataLoader.loadDepth)";
		final double invLA=1.0/Math.max(a.size(), 1), invLB=1.0/Math.max(b.size(), 1);
		final double s2=2*sigma*sigma;
		double sum=0;
		for(int s=0; s<n; s++){
			final float ca=a.depth(s), cb=b.depth(s);
			if(ca<DET && cb<DET){sum+=coabs[s]; continue;}//Co-absence carries capped evidence
			final double mid=Math.max(0.5*(ca+cb), 0.05);
			final double U=K*mid*(invLA+invLB)+s2*mid*mid;
			final double d=cb-ca;
			final double l0=-0.5*(LOG2PI+Math.log(U))-d*d/(2*U);
			final double lb=Math.max(mid, 1.5);
			final double l1=-Math.log(2*lb)-Math.abs(d)/lb;
			sum+=l0-l1;
		}
		return (float)sum;
	}

	/**
	 * Computes per-sample absence rates (bp-weighted fraction of the assembly
	 * below DET) and the capped co-absence evidence terms.  Call after depth
	 * loading AND after any condensation, so the terms match the live columns.
	 */
	public static synchronized void calibrate(ArrayList<Contig> contigs, int numDepths){
		final double[] absentBp=new double[numDepths];
		double totalBp=0;
		for(Contig c : contigs){
			final long size=c.size();
			totalBp+=size;
			for(int s=0; s<numDepths; s++){
				if(c.depth(s)<DET){absentBp[s]+=size;}
			}
		}
		coabs=new float[numDepths];
		StringBuilder sb=new StringBuilder("CovLR co-absence evidence (nats):");
		for(int s=0; s<numDepths; s++){
			final double rate=absentBp[s]/Math.max(totalBp, 1);
			coabs[s]=(float)Math.min(-Math.log(Math.max(rate, 1e-9)), COABS_MAX);
			sb.append(String.format(" %.2f", coabs[s]));
		}
		System.err.println(sb.toString());
	}

	static final double LOG2PI=Math.log(2*Math.PI);
	/** Read-sampling noise constant; warpbin default, flag covlrk. */
	public static float K=350f;
	/** Depth below which a sample counts as absent; flag covlrdet. */
	public static float DET=0.5f;
	/** Relative biological dispersion (s2=2*sigma^2); warpbin SIG_DEF, flag covlrsigma. */
	public static float sigma=0.12f;
	/** Cap on per-sample co-absence evidence, nats. */
	public static float COABS_MAX=2.0f;
	/** Per-sample co-absence terms; filled by calibrate(). */
	private static float[] coabs;
}
