package bin;

import java.util.ArrayList;

import shared.Tools;

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

	/**
	 * Combined depth-compatibility score in log-odds (nats) for depthoracle=fused.
	 * Returns z = BIAS[n] + W_LNRATIO[n]*ln(depthRatio) + W_COV[n]*covariance
	 * + W_COVLR[n]*covLR, the pre-sigmoid activation of a logistic regression fit
	 * per sample count on the labeled depthvectors corpus (scripts/43_fit_depth.sh,
	 * negatives 4x oversampled for contamination priority): P(same genome)=sigmoid(z),
	 * so z>=0 favors merging.  ONE threshold in log-odds space replaces the three
	 * separate legacy depth cutoffs; stringency enters additively at the Oracle call
	 * site.  Constants are the fitted first-layer weights (RegressionTrainer folds
	 * input standardization into layer-1 at export, so they consume RAW inputs);
	 * verified to reproduce ml.BBNetApply to <1e-6 (scripts/44_verify_mapping.sh,
	 * 2026-09-10).  Signs match the data: bigger depthRatio and bigger covariance
	 * difference LOWER z (less same); bigger covLR RAISES it.  Booleanizing before
	 * combining would discard vote confidence (Brian: a 0.4999 reject must not outvote
	 * a 0.943 accept) - magnitude decides, thresholded ONCE.
	 */
	public static float depthScore(final float depthRatio, final float covariance,
			final float covLR, final int numDepths, final long minSize){
		final int n=Tools.max(1, Tools.min(numDepths, TABLE_MAX));
		final float lnRatio=(float)Math.log(Math.max(depthRatio, 1e-6f));
		return BIAS[n]+mLnRatio*W_LNRATIO[n]*lnRatio+mCov*W_COV[n]*covariance
				+mLr*W_COVLR[n]*covLR;
	}

	/** Global term multipliers for the fused score, default 1 (=raw fitted fusion).  These
	 * are the FIRST descent knobs (flags covlrmratio/covlrmcov/covlrmlr): the logistic fit
	 * minimized P(same|depth) MSE, not binning Total/Contam, so these reweight each term
	 * toward the real objective via the sensitivity descent (Brian: tune constants, follow
	 * the best gradient on score AND contamination to breakeven, then switch). */
	public static float mLnRatio=1f, mCov=1f, mLr=1f;

	/** Highest sample count with its own fitted coefficient row; larger n clamps here. */
	public static final int TABLE_MAX=9;
	/**
	 * Fitted logistic-regression coefficients indexed by sample count (raw-input space:
	 * bias, then weights on ln(depthRatio), covariance, covLR).  Rows 1,2,3,4,6,9 are
	 * fitted (scripts/43_fit_depth.sh on the 2026-09-10 corpus); rows 5,7,8 are linear
	 * interpolations of their fitted neighbors; row 0 mirrors row 1 (n clamps to [1,9]).
	 * These REPLACE the earlier SEED constants; regenerate by refitting the corpus and
	 * re-running scripts/44_verify_mapping.sh (must reproduce ml.BBNetApply to <1e-6). */
	public static float[] BIAS=     {-3.5375f, -3.5375f, -0.824739f, 0.462639f, 1.446205f, 1.950495f, 2.454785f, 2.595168f, 2.735552f, 2.875935f};
	public static float[] W_LNRATIO={-3.019216f, -3.019216f, -6.83011f, -4.566097f, -4.481306f, -4.022081f, -3.562856f, -3.231043f, -2.899229f, -2.567416f};
	public static float[] W_COV=    {0f, 0f, -17.713409f, -7.092163f, -9.09537f, -10.957198f, -12.819026f, -15.505455f, -18.191884f, -20.878313f};
	public static float[] W_COVLR=  {0.371015f, 0.371015f, 0.080259f, 0.101604f, 0.088055f, 0.072849f, 0.057642f, 0.051325f, 0.045008f, 0.038691f};

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
