package prot;

import java.util.Arrays;

/**
 * Bounded, I/O-independent metric core for the Step-1 family-artifact assay
 * (plans/STEP1_ASSAY_HARNESS_DESIGN_v1.md #4.1/4.3/4.4). Every method operates purely on
 * already-extracted score arrays -- no manifest parsing, no artifact loading, no file I/O, no
 * cross-family/cross-cell orchestration. Margin arithmetic itself is NOT reimplemented here --
 * callers combine {@link #maxOtherScore(double[], int)} with {@link FamilyArtifactScore#margin}
 * / {@link FamilyArtifactScore#clearsMargin} directly.
 *
 * @author Eru
 */
public final class Step1AssayMetrics {

	private Step1AssayMetrics(){}

	// ==================== 4.1: exhaustive top-1 + confusion ====================

	/** Every held-out positive query falls into exactly one of these, exhaustive by
	 * construction (design doc #4.1 v3 -- widened from the v2 own/other-tie-only definition
	 * after UMP45 found a top-tie-among-two-wrong-families case that satisfied neither v2
	 * bucket). */
	public enum Top1Bucket { CORRECT, CONFUSED, TIED }

	/**
	 * Classifies one positive query's top-1 outcome against its full ChallengeSet(X).
	 * @param scores one score per ChallengeSet(X) artifact, in any fixed order; must be finite.
	 * @param ownIndex index of X's own artifact within scores.
	 * @return CORRECT if X's own artifact is the UNIQUE maximum; CONFUSED if a single
	 *         different artifact is the unique maximum; TIED if the maximum is achieved by
	 *         MORE THAN ONE artifact -- regardless of whether X's own artifact is among the
	 *         tied set (scope line 165: a tie "does NOT receive top-1 credit", full stop).
	 *         Ties are exact double equality, matching real quantized scorer outputs.
	 */
	public static Top1Bucket classifyTop1(double[] scores, int ownIndex){
		requireFinite(scores, "scores");
		if(scores.length==0){throw new IllegalArgumentException("scores must be non-empty.");}
		if(ownIndex<0 || ownIndex>=scores.length){
			throw new IllegalArgumentException("ownIndex "+ownIndex+" out of range [0,"+scores.length+").");
		}
		final double max=maxOf(scores);
		int tieCount=0;
		boolean ownIsMax=false;
		for(int i=0; i<scores.length; i++){
			if(scores[i]==max){
				tieCount++;
				if(i==ownIndex){ownIsMax=true;}
			}
		}
		assert(tieCount>=1) : "At least the maximum's own index must match itself.";
		if(tieCount>1){return Top1Bucket.TIED;}
		return ownIsMax ? Top1Bucket.CORRECT : Top1Bucket.CONFUSED;
	}

	// ==================== 4.3: margin support ====================

	/**
	 * Max score among ChallengeSet(X) artifacts OTHER than ownIndex -- the "strongest other"
	 * term {@link FamilyArtifactScore#margin} and {@code clearsMargin} expect as their second
	 * argument.
	 * @throws IllegalArgumentException if fewer than 2 artifacts are present (no "other" to
	 *         compare against) or ownIndex is out of range.
	 */
	public static double maxOtherScore(double[] scores, int ownIndex){
		requireFinite(scores, "scores");
		if(scores.length<2){
			throw new IllegalArgumentException("Need at least 2 ChallengeSet artifacts (own + "
				+">=1 other) to compute a margin; got "+scores.length+".");
		}
		if(ownIndex<0 || ownIndex>=scores.length){
			throw new IllegalArgumentException("ownIndex "+ownIndex+" out of range [0,"+scores.length+").");
		}
		double max=Double.NEGATIVE_INFINITY;
		for(int i=0; i<scores.length; i++){
			if(i!=ownIndex && scores[i]>max){max=scores[i];}
		}
		return max;
	}

	// ==================== 4.4: AUROC (Mann-Whitney / rank-sum, midrank ties) ====================

	/**
	 * Empirical Mann-Whitney U / rank-sum AUROC estimator, exactly equivalent to trapezoidal
	 * ROC integration with no binning step (design doc #4.4). Positive and negative values are
	 * ranked together ascending; tied values (across the pooled set, positive+negative alike)
	 * receive the MIDRANK -- the average of the ranks the tie group would occupy.
	 * {@code AUROC = (sumRankPos - nPos*(nPos+1)/2) / (nPos*nNeg)}, computed in WIDENED
	 * (double) arithmetic throughout -- {@code nPos*(nPos+1)} as plain int multiplication
	 * overflows {@code Integer.MAX_VALUE} once {@code nPos>46,340} (Elly's review catch), a
	 * scale a real family/query pool can plausibly reach.
	 * <p>Algebraic identity (used as a self-test, holds exactly regardless of ties):
	 * {@code auroc(pos,neg) + auroc(neg,pos) == 1.0}.
	 */
	public static double auroc(double[] positive, double[] negative){
		requireFinite(positive, "positive");
		requireFinite(negative, "negative");
		final int nPos=positive.length, nNeg=negative.length;
		if(nPos==0 || nNeg==0){
			throw new IllegalArgumentException("AUROC needs at least one positive and one negative "
				+"value; got nPos="+nPos+" nNeg="+nNeg+".");
		}
		final int n=nPos+nNeg;
		final double[] allVals=new double[n];
		final boolean[] isPos=new boolean[n];
		System.arraycopy(positive, 0, allVals, 0, nPos);
		Arrays.fill(isPos, 0, nPos, true);
		System.arraycopy(negative, 0, allVals, nPos, nNeg);

		//TODO: perf - primitive sort to avoid Integer boxing. Non-blocking (UMP45/Elly review):
		//this runs at sample scale (~60 families), not a hot path, so boxing cost is negligible
		//here -- flagged for a future cleanup pass, not this correctness fix.
		final Integer[] idx=new Integer[n];
		for(int i=0; i<n; i++){idx[i]=i;}
		Arrays.sort(idx, (a, b)->Double.compare(allVals[a], allVals[b]));

		final double[] rank=new double[n];
		int i=0;
		while(i<n){
			int j=i;
			while(j+1<n && allVals[idx[j+1]]==allVals[idx[i]]){j++;}
			//1-based ranks i+1..j+1, midrank is their average.
			final double midrank=((i+1)+(j+1))/2.0;
			for(int k=i; k<=j; k++){rank[idx[k]]=midrank;}
			i=j+1;
		}

		double sumRankPos=0;
		for(int k=0; k<n; k++){ if(isPos[k]){sumRankPos+=rank[k];} }

		//nPos*(nPos+1.0) forces double arithmetic for the WHOLE product (the +1.0 promotes
		//nPos before the multiply) -- plain int nPos*(nPos+1) silently overflows
		//Integer.MAX_VALUE past nPos=46,340 and this project's real query pools can plausibly
		//exceed that (Elly's review catch).
		final double result=(sumRankPos-nPos*(nPos+1.0)/2.0)/((double)nPos*nNeg);
		assert(result>=-1e-9 && result<=1.0+1e-9) : "AUROC is mathematically bounded to [0,1] "
			+"(a Mann-Whitney rank-sum statistic can never exceed its own pair count); result="
			+result+" nPos="+nPos+" nNeg="+nNeg+" sumRankPos="+sumRankPos+" -- this exact bound "
			+"is what the nPos*(nPos+1) int-overflow bug violated before it was fixed.";
		return result;
	}

	// ==================== 4.4: conservative TPR@FPR<=maxFpr ====================

	/** One chosen operating point: the threshold, its TPR, and its REALIZED FPR (guaranteed
	 * {@code <= maxFpr} by construction). */
	public static final class RocPoint {
		public final double threshold, tpr, fpr;
		RocPoint(double threshold, double tpr, double fpr){
			this.threshold=threshold; this.tpr=tpr; this.fpr=fpr;
		}
		@Override
		public String toString(){
			return "RocPoint(threshold="+threshold+", tpr="+tpr+", fpr="+fpr+")";
		}
	}

	/**
	 * Conservative empirical ROC operating point (design doc #4.4, round-2 fix replacing an
	 * unsound nearest-rank threshold pick that could exceed the FPR bound under ties; round-3
	 * fix, Elly, correcting the candidate SET itself). Candidate thresholds T = every distinct
	 * value across BOTH positive and negative scores, plus +Infinity (which always yields
	 * FPR=0, guaranteeing the valid set is never empty). A {@code >=t} threshold rule's
	 * classification can only change AT an observed score -- between two observed values
	 * nothing crosses the boundary -- so the complete candidate set must draw from both classes,
	 * not negatives alone. **Round-3 counterexample (Elly) that the negatives-only set missed:**
	 * negatives all {@code =1}, positives all {@code =2}, {@code maxFpr=0.05}. The
	 * negatives-only set was {1, +Inf}: t=1 gives FPR=1 (invalid), t=+Inf gives TPR=0 -- so the
	 * old algorithm reported TPR=0 even though the real threshold t=2 (a POSITIVE-only value,
	 * absent from the negatives-only set) achieves FPR=0 AND TPR=1. For each {@code t in T}:
	 * {@code FPR(t) = count(negative >= t) / nNeg}; restricted to {@code FPR(t) <= maxFpr};
	 * among those, {@code TPR(t) = count(positive >= t) / nPos} is maximized; ties for the max
	 * TPR are broken by taking the LARGEST threshold (fewest negatives admitted -- the smallest
	 * realized FPR among equally-good options). No monotonicity is assumed -- every distinct
	 * candidate is evaluated independently, matching the design text's "enumerate candidate
	 * thresholds" literally.
	 */
	public static RocPoint conservativeTprAtFpr(double[] positive, double[] negative, double maxFpr){
		requireFinite(positive, "positive");
		requireFinite(negative, "negative");
		final int nPos=positive.length, nNeg=negative.length;
		if(nPos==0 || nNeg==0){
			throw new IllegalArgumentException("TPR@FPR needs at least one positive and one negative "
				+"value; got nPos="+nPos+" nNeg="+nNeg+".");
		}
		if(!Double.isFinite(maxFpr) || maxFpr<0 || maxFpr>1){
			throw new IllegalArgumentException("maxFpr must be finite in [0,1]; got "+maxFpr);
		}

		//Distinct values from BOTH classes, descending, plus +Infinity (round-3 fix, Elly: a
		//negatives-only candidate set misses achievable zero-FPR operating points whenever the
		//real boundary sits at a positive-only value strictly above every negative).
		final java.util.TreeSet<Double> candidates=new java.util.TreeSet<Double>(java.util.Collections.reverseOrder());
		candidates.add(Double.POSITIVE_INFINITY);
		for(double v : negative){candidates.add(v);}
		for(double v : positive){candidates.add(v);}

		double bestTpr=-1, bestThreshold=Double.NaN, bestFpr=Double.NaN;
		for(double t : candidates){
			final int fp=countGE(negative, t);
			final double fprAtT=fp/(double)nNeg;
			if(fprAtT>maxFpr){continue;}
			final int tp=countGE(positive, t);
			final double tprAtT=tp/(double)nPos;
			//Descending iteration order means the FIRST t to reach a given TPR is already the
			//LARGEST such t -- a strict '>' here never lets a smaller-t tie overwrite it.
			if(tprAtT>bestTpr){
				bestTpr=tprAtT; bestThreshold=t; bestFpr=fprAtT;
			}
		}
		assert(bestTpr>=0) : "+Infinity always yields FPR=0<=maxFpr, so the valid set can never be empty.";
		return new RocPoint(bestThreshold, bestTpr, bestFpr);
	}

	// ==================== shared helpers ====================

	static int countGE(double[] values, double threshold){
		int c=0;
		for(double v : values){ if(v>=threshold){c++;} }
		return c;
	}

	static double maxOf(double[] values){
		double max=Double.NEGATIVE_INFINITY;
		for(double v : values){ if(v>max){max=v;} }
		return max;
	}

	static void requireFinite(double[] values, String name){
		if(values==null){throw new IllegalArgumentException(name+" must not be null.");}
		for(int i=0; i<values.length; i++){
			if(!Double.isFinite(values[i])){
				throw new IllegalArgumentException(name+"["+i+"]="+values[i]+" is not finite -- a "
					+"missing/failed score must never reach this pure metric core (see harness design "
					+"doc #2: a missing assay row is a failed assay, thrown loud upstream, never a "
					+"silent 0 or NaN passed down here).");
			}
		}
	}

	// ==================== self-test ====================

	public static void main(String[] args){
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){selftest(); return;}
		throw new RuntimeException("Usage: java -ea prot.Step1AssayMetrics selftest");
	}

	static void selftest(){
		testClassifyTop1();
		testMaxOtherScore();
		testAuroc();
		testConservativeTprAtFpr();
		System.err.println("Step1AssayMetrics SELFTEST PASS.");
	}

	static void testClassifyTop1(){
		//Clear win.
		expect(Top1Bucket.CORRECT, classifyTop1(new double[]{10,3,7}, 0));
		//Clear loss, unique other max.
		expect(Top1Bucket.CONFUSED, classifyTop1(new double[]{3,10,7}, 0));
		//Own/other two-way tie (own IS in the tied set).
		expect(Top1Bucket.TIED, classifyTop1(new double[]{10,10,3}, 0));
		//THE v3 fix case (UMP45's finding): a tie between two DIFFERENT WRONG families, own
		//absent from the tied set entirely -- must be TIED, not uncategorized.
		expect(Top1Bucket.TIED, classifyTop1(new double[]{3,10,10}, 0));
		//Three-way tie including own.
		expect(Top1Bucket.TIED, classifyTop1(new double[]{10,10,10}, 1));
		//Degenerate single-artifact challenge set: own is trivially the unique max.
		expect(Top1Bucket.CORRECT, classifyTop1(new double[]{5}, 0));
		//Exhaustiveness bound: correct+confused+tied partition must be well-defined for a large
		//randomized batch, including deliberately-injected ties (never an exception, never an
		//unmapped case).
		final java.util.Random rnd=new java.util.Random(42);
		for(int trial=0; trial<2000; trial++){
			final int n=1+rnd.nextInt(6);
			final double[] scores=new double[n];
			//Small integer range forces frequent ties.
			for(int i=0; i<n; i++){scores[i]=rnd.nextInt(4);}
			final int own=rnd.nextInt(n);
			final Top1Bucket b=classifyTop1(scores, own);
			if(b==null){throw new RuntimeException("SELFTEST FAILED: classifyTop1 returned null.");}
		}
		//Invalid inputs.
		expectThrows(()->classifyTop1(new double[]{1,2}, 5), "ownIndex out of range");
		expectThrows(()->classifyTop1(new double[]{1,Double.NaN}, 0), "NaN score");
		expectThrows(()->classifyTop1(new double[]{}, 0), "empty scores");
		System.err.println("  testClassifyTop1: PASS (own-win, other-win, own-tie, own-absent-tie "
			+"[v3 fix], 3-way tie, singleton, 2000 randomized exhaustiveness trials, 3 invalid-input cases)");
	}

	static void testMaxOtherScore(){
		if(maxOtherScore(new double[]{5,1,9,3}, 0)!=9){throw new RuntimeException("SELFTEST FAILED: maxOtherScore basic.");}
		if(maxOtherScore(new double[]{5,1,9,3}, 2)!=5){throw new RuntimeException("SELFTEST FAILED: maxOtherScore excludes own index, not just own value.");}
		//Cross-check against FamilyArtifactScore.margin using this helper's output.
		final double[] s={7,2,9};
		final double margin=FamilyArtifactScore.margin(s[0], maxOtherScore(s, 0));
		if(margin!=7-9){throw new RuntimeException("SELFTEST FAILED: margin composition mismatch.");}
		expectThrows(()->maxOtherScore(new double[]{5}, 0), "fewer than 2 artifacts");
		expectThrows(()->maxOtherScore(new double[]{5,1}, 9), "ownIndex out of range");
		System.err.println("  testMaxOtherScore: PASS (basic, own-index-excluded-not-value, "
			+"margin composition, 2 invalid-input cases)");
	}

	static void testAuroc(){
		//Perfect separation.
		if(auroc(new double[]{10,11,12}, new double[]{1,2,3})!=1.0){throw new RuntimeException("SELFTEST FAILED: perfect-separation AUROC != 1.0.");}
		//Perfect anti-separation.
		if(auroc(new double[]{1,2,3}, new double[]{10,11,12})!=0.0){throw new RuntimeException("SELFTEST FAILED: perfect-anti-separation AUROC != 0.0.");}
		//Total overlap, single tied value both sides -> exactly 0.5.
		if(auroc(new double[]{5,5}, new double[]{5,5})!=0.5){throw new RuntimeException("SELFTEST FAILED: fully-tied AUROC != 0.5.");}
		//Hand-computed worked example: pos=[1,2,3], neg=[2,2,4].
		//Pooled ascending: 1(pos),2(pos),2(neg),2(neg),3(pos),4(neg) -> ranks 1,2,3,4,5,6.
		//The three 2's occupy ranks 2,3,4 -> midrank 3 each. Final ranks: 1(pos)=1, 2(pos)=3,
		//2(neg)=3, 2(neg)=3, 3(pos)=5, 4(neg)=6. sumRankPos = 1+3+5 = 9. nPos=3,nNeg=3.
		//AUROC = (9 - 3*4/2) / (3*3) = (9-6)/9 = 3/9 = 1/3.
		final double handComputed=auroc(new double[]{1,2,3}, new double[]{2,2,4});
		if(Math.abs(handComputed-(1.0/3.0))>1e-12){
			throw new RuntimeException("SELFTEST FAILED: hand-computed midrank AUROC example gave "
				+handComputed+", expected 0.3333333333333333.");
		}
		//Algebraic identity: swapping labels must sum to exactly 1.0, regardless of ties. Checked
		//across several randomized configurations, not just one.
		final java.util.Random rnd=new java.util.Random(7);
		for(int trial=0; trial<200; trial++){
			final int nPos=1+rnd.nextInt(8), nNeg=1+rnd.nextInt(8);
			final double[] pos=new double[nPos], neg=new double[nNeg];
			for(int i=0; i<nPos; i++){pos[i]=rnd.nextInt(5);}//small range forces ties
			for(int i=0; i<nNeg; i++){neg[i]=rnd.nextInt(5);}
			final double a=auroc(pos, neg), b=auroc(neg, pos);
			if(Math.abs((a+b)-1.0)>1e-9){
				throw new RuntimeException("SELFTEST FAILED: auroc(pos,neg)+auroc(neg,pos) = "+(a+b)
					+" != 1.0 for pos="+Arrays.toString(pos)+" neg="+Arrays.toString(neg));
			}
		}
		expectThrows(()->auroc(new double[]{}, new double[]{1}), "empty positive");
		expectThrows(()->auroc(new double[]{1}, new double[]{}), "empty negative");
		expectThrows(()->auroc(new double[]{Double.NaN}, new double[]{1}), "NaN positive");

		//THE overflow regression (Elly's review catch): nPos*(nPos+1) as plain int arithmetic
		//overflows Integer.MAX_VALUE (2,147,483,647) once nPos>46,340 -- 50,000 crosses it
		//(50000*50001=2,500,050,000). Fully tied (every value identical, both classes) must
		//still give exactly 0.5 by symmetry; the pre-fix bug instead wrapped to a large
		//negative int, producing an AUROC wildly outside [0,1] (the assert in auroc() now also
		//catches this class directly, but this regression pins the exact scale that triggers it).
		final double[] bigTiedPos=new double[50000];
		Arrays.fill(bigTiedPos, 3.0);
		final double[] bigTiedNeg=new double[3];
		Arrays.fill(bigTiedNeg, 3.0);
		final double bigTiedAuroc=auroc(bigTiedPos, bigTiedNeg);
		if(Math.abs(bigTiedAuroc-0.5)>1e-9){
			throw new RuntimeException("SELFTEST FAILED: 50,000-positive fully-tied regression gave "
				+bigTiedAuroc+", expected exactly 0.5 -- the nPos*(nPos+1) int-overflow bug is back.");
		}
		System.err.println("  testAuroc: PASS (perfect separation, perfect anti-separation, full-tie=0.5, "
			+"hand-computed midrank worked example, 200-trial swap-identity property test, "
			+"50,000-positive fully-tied overflow regression, 3 invalid-input cases)");
	}

	static void testConservativeTprAtFpr(){
		//THE round-3 fix case, Elly's EXACT minimal counterexample: negatives all =1, positives
		//all =2, maxFpr=0.05. The negatives-only candidate set {1,+Inf} could never find t=2 --
		//it would report TPR=0 despite a real threshold (t=2, a POSITIVE-only value) achieving
		//FPR=0 AND TPR=1. Candidates now = {+inf,2,1}. t=2: FP=0(no neg>=2),FPR=0<=0.05 valid,
		//TP=2(both pos>=2),TPR=1 -- strictly better than +inf's TPR=0, so it must win.
		RocPoint pMin=conservativeTprAtFpr(new double[]{2,2}, new double[]{1,1}, 0.05);
		if(pMin.threshold!=2 || pMin.tpr!=1.0 || pMin.fpr!=0.0){
			throw new RuntimeException("SELFTEST FAILED: Elly's minimal positive-only-threshold "
				+"counterexample gave "+pMin+", expected threshold=2 tpr=1 fpr=0.");
		}

		//Richer positive-only-threshold case, hand-verified: positives {10,8,3}, negatives
		//{9,5,1}, maxFpr=0.34 (FP<=1 of 3, i.e. FPR<=1/3). Full candidate set (union, desc):
		//+inf,10,9,8,5,3,1.
		//t=+inf: FP=0,TPR=0.  t=10: FP=0,TP=1(10),TPR=1/3.  t=9: FP=1(9),FPR=1/3 valid,TP=1,TPR=1/3.
		//t=8: FP=1(9>=8),FPR=1/3 valid, TP=2(10,8>=8),TPR=2/3 -- a POSITIVE-only value (8, not a
		//negative score) gives the best valid point. t=5/3/1: FP=2 or 3, FPR>0.34, invalid.
		RocPoint p=conservativeTprAtFpr(new double[]{10,8,3}, new double[]{9,5,1}, 0.34);
		if(p.threshold!=8 || Math.abs(p.tpr-2.0/3.0)>1e-12 || Math.abs(p.fpr-1.0/3.0)>1e-12){
			throw new RuntimeException("SELFTEST FAILED: positive-only-threshold case gave "+p
				+", expected threshold=8 tpr=2/3 fpr=1/3 (the negatives-only candidate set would "
				+"have missed threshold=8 entirely).");
		}

		//THE round-2 fix case: a cluster of TIED negatives right at the boundary, with NO escape
		//route (positives sit AT the same tied value, not above it, so no positive-only threshold
		//can rescue this one -- forces the algorithm to actually confront the tie). Negatives
		//{5,5,5,5,1}, positives {5,5}, maxFpr=0.05. Union candidates: +inf,5,1. A naive
		//"nearest-rank at ceil(5%*5)=1" pick would choose value 5, then a >= comparison sweeps in
		//all four 5's at once (FP 0->4, FPR 0->0.8), blowing far past 0.05. t=5 here: FP=4/5=0.8
		//>0.05, correctly rejected; t=1: FP=5,FPR=1, rejected. Only +inf survives.
		final double[] tiedNeg={5,5,5,5,1};
		final double[] tiedPos={5,5};
		RocPoint p2=conservativeTprAtFpr(tiedPos, tiedNeg, 0.05);
		if(p2.fpr>0.05+1e-12){
			throw new RuntimeException("SELFTEST FAILED: realized FPR "+p2.fpr+" exceeds the 0.05 bound "
				+"-- this is exactly the tied-boundary bug the conservative algorithm exists to prevent.");
		}
		if(p2.threshold!=Double.POSITIVE_INFINITY || p2.tpr!=0.0 || p2.fpr!=0.0){
			throw new RuntimeException("SELFTEST FAILED: tied-boundary (no-escape) case gave "+p2
				+", expected threshold=+Infinity tpr=0 fpr=0 (no finite threshold keeps FPR<=0.05 at nNeg=5).");
		}

		//Tie-break-for-max-TPR case with a NON-zero tied TPR: positives {10}, negatives
		//{8,6,1}, maxFpr=0.4 (FP<=1 of 3). Union candidates desc: +inf,10,8,6,1.
		//t=+inf: FP0,TPR0.  t=10: FP=0(no neg>=10),FPR=0 valid,TP=1,TPR=1.
		//t=8: FP=1(8>=8),FPR=1/3<=0.4 valid,TP=1(10>=8),TPR=1 -- TIES with t=10 at TPR=1.
		//t=6: FP=2(8,6),FPR=2/3 invalid.  t=1: FP=3,FPR=1 invalid.
		//Two valid thresholds (10 and 8) both reach the max TPR=1; the LARGER one (10, giving the
		//smaller realized FPR=0 instead of 8's FPR=1/3) must be chosen.
		RocPoint p3=conservativeTprAtFpr(new double[]{10}, new double[]{8,6,1}, 0.4);
		if(p3.threshold!=10 || p3.tpr!=1.0 || p3.fpr!=0.0){
			throw new RuntimeException("SELFTEST FAILED: TPR-tie-break case gave "+p3
				+", expected the largest tied threshold (10, tpr=1, fpr=0), not the smaller one (8, fpr=1/3).");
		}

		//Randomized property test: realized FPR must NEVER exceed maxFpr, across many random
		//configurations including heavy tie clustering (small integer range).
		final java.util.Random rnd=new java.util.Random(99);
		for(int trial=0; trial<500; trial++){
			final int nPos=1+rnd.nextInt(6), nNeg=1+rnd.nextInt(6);
			final double[] posR=new double[nPos], negR=new double[nNeg];
			for(int i=0; i<nPos; i++){posR[i]=rnd.nextInt(5);}
			for(int i=0; i<nNeg; i++){negR[i]=rnd.nextInt(5);}
			final double maxFpr=rnd.nextDouble();
			final RocPoint rp=conservativeTprAtFpr(posR, negR, maxFpr);
			if(rp.fpr>maxFpr+1e-12){
				throw new RuntimeException("SELFTEST FAILED: randomized trial realized FPR "+rp.fpr
					+" > maxFpr "+maxFpr+" for pos="+Arrays.toString(posR)+" neg="+Arrays.toString(negR));
			}
		}

		expectThrows(()->conservativeTprAtFpr(new double[]{}, new double[]{1}, 0.05), "empty positive");
		expectThrows(()->conservativeTprAtFpr(new double[]{1}, new double[]{}, 0.05), "empty negative");
		expectThrows(()->conservativeTprAtFpr(new double[]{1}, new double[]{1}, 1.5), "maxFpr out of range");
		expectThrows(()->conservativeTprAtFpr(new double[]{1}, new double[]{1}, -0.1), "negative maxFpr");
		System.err.println("  testConservativeTprAtFpr: PASS (Elly's minimal positive-only-threshold "
			+"counterexample [round-3 fix], a richer positive-only-threshold case, the no-escape "
			+"tied-boundary round-2-fix case [never exceeds 5% FPR], non-zero-TPR largest-threshold "
			+"tie-break, 500-trial randomized FPR-bound property test, 4 invalid-input cases)");
	}

	interface ThrowingRunnable { void run(); }

	static void expectThrows(ThrowingRunnable r, String label){
		try{
			r.run();
			throw new RuntimeException("SELFTEST FAILED: expected an exception for '"+label+"' but none was thrown.");
		}catch(IllegalArgumentException expected){
			//correct
		}
	}

	static void expect(Object want, Object got){
		if(!want.equals(got)){throw new RuntimeException("SELFTEST FAILED: expected "+want+" but got "+got+".");}
	}
}
