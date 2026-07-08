package cardinality;

/**
 * TailMLE: analytic maximum-likelihood extractor for ExpandedTwinTailLogLog5.
 * <p>
 * ETTLL5's update rule is set-semantic (max-anchored window, monotone OR-shift),
 * so its state is a pure function of the set of hashes observed.  That makes the
 * exact likelihood factorize with no calibration tables:
 * <p>
 *   L(bank | mu) = exp(-2*mu*Q(M)) * PROD over tracked cells Bernoulli(1-e^-lambda)
 * <p>
 * where M is the bank's max expanded tier, Q(M) is the hash mass above tier M,
 * mu = n/(2*numBuckets) is the per-tail Poisson rate, and lambda = mu*q(tier).
 * The tier masses q form an exact geometric ladder of ratio sqrt(2) because the
 * mantissa threshold is (2-sqrt(2)) (the geometric midpoint of the octave).
 * <p>
 * All bank terms aggregate into per-tier set/unset counts plus one scalar
 * (total mass-above), so a full Newton solve on ln(mu) costs O(128) per
 * iteration after a single O(numBuckets) scan.  The same likelihood subsumes
 * linear counting: untouched banks contribute exp(-2*mu), which is exactly the
 * LC empty-bucket term.
 *
 * @author Amber
 * @date July 2026
 */
public final class TailMLE {

	/*--------------------------------------------------------------*/
	/*----------------        Tier Mass Ladder       ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of expanded tiers (2 per NLZ octave, 64 octaves). */
	static final int MAX_TIER=128;
	/** q[t] = probability a uniform hash lands in expanded tier t. */
	static final double[] Q_MASS=new double[MAX_TIER];
	/** ABOVE[t] = sum of Q_MASS[t..MAX_TIER-1]; Q(m) = ABOVE[m+1]. */
	static final double[] ABOVE=new double[MAX_TIER+1];

	static{
		final double lo=2.0-Math.sqrt(2.0);   // P(mantissa=0) within an octave
		final double hi=Math.sqrt(2.0)-1.0;   // P(mantissa=1) within an octave
		for(int nlz=0; nlz<64; nlz++){
			final double oct=Math.pow(2.0, -(nlz+1)); // P(rawNLZ==nlz)
			Q_MASS[2*nlz]=oct*lo;
			Q_MASS[2*nlz+1]=oct*hi;
		}
		double s=0;
		ABOVE[MAX_TIER]=0;
		for(int t=MAX_TIER-1; t>=0; t--){
			s+=Q_MASS[t];
			ABOVE[t]=s;
		}
	}

	/** Mass strictly above tier m (m may be -1 for "no max yet"). */
	private static double massAbove(int m){
		if(m<-1){m=-1;}
		if(m>=MAX_TIER-1){return 0;}
		return ABOVE[m+1];
	}

	/*--------------------------------------------------------------*/
	/*----------------          Estimation           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Maximum-likelihood cardinality estimate from raw ETTLL5 state.
	 * @param regs packed register array (4 x 15-bit banks per long)
	 * @param globalExp shared floor in expanded-tier units
	 * @param numBuckets number of banks (each holding 2 tails)
	 * @return estimated number of distinct elements added
	 */
	public static double estimate(final long[] regs, final int globalExp, final int numBuckets){
		final int HIST_LEN=ExpandedTwinTailLogLog5.HIST_LEN;      // 5
		final int NUM_TAILS=ExpandedTwinTailLogLog5.NUM_TAILS;    // 2
		final int EXP_SHIFT=ExpandedTwinTailLogLog5.EXP_SHIFT;    // 10
		final int EXP_MASK=ExpandedTwinTailLogLog5.EXP_MASK;      // 31
		final int REG_BITS=ExpandedTwinTailLogLog5.REG_BITS;      // 15
		final int RPW=ExpandedTwinTailLogLog5.REGS_PER_WORD;      // 4
		final long RMASK=ExpandedTwinTailLogLog5.REG_MASK_L;

		final long[] setCount=new long[MAX_TIER];
		final long[] unsetCount=new long[MAX_TIER];
		double massAboveSum=0;     // sum over banks of Q(M_bank); untouched banks contribute 1.0
		long totalSet=0;

		for(int i=0; i<numBuckets; i++){
			final int reg=(int)((regs[i/RPW]>>>((i%RPW)*REG_BITS))&RMASK);
			if(reg==0 && globalExp==0){
				massAboveSum+=1.0;   // untouched: all mass is "above" a nonexistent max
				continue;
			}
			final int localExp=(reg>>>EXP_SHIFT)&EXP_MASK;
			final int m=globalExp+localExp;
			massAboveSum+=massAbove(m);
			for(int d=0; d<HIST_LEN; d++){
				final int tier=m-d;
				if(tier<0){break;}
				if(tier>=MAX_TIER){continue;}
				final int bitPos=(HIST_LEN-1)-d;
				for(int t=0; t<NUM_TAILS; t++){
					final int bit=(reg>>>(t*HIST_LEN+bitPos))&1;
					if(bit!=0){setCount[tier]++; totalSet++;}
					else{unsetCount[tier]++;}
				}
			}
		}

		if(totalSet==0){return 0;}

		final double x=solve(setCount, unsetCount, massAboveSum);
		final double mu=Math.exp(x);
		return 2.0*numBuckets*mu;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Newton Solver         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Solves dL/dx = 0 for x = ln(mu).  The log-likelihood is strictly concave
	 * in x, so a bracketed Newton iteration is globally convergent.
	 */
	private static double solve(final long[] S, final long[] U, final double A){
		// Bracket: dL/dx > 0 at very small mu (approaches sum(S) > 0),
		// dL/dx -> -infinity at large mu.
		double xLo=-34, xHi=44;   // mu from ~2e-15 to ~1e19
		double x=0;
		// Ensure bracket validity
		if(deriv(xLo, S, U, A)<=0){return xLo;}   // degenerate: essentially zero
		for(int iter=0; iter<200; iter++){
			final double d1=deriv(x, S, U, A);
			if(d1>0){xLo=x;}else{xHi=x;}
			final double d2=deriv2(x, S, U, A);
			double xNew=(d2<0) ? x-d1/d2 : Double.NaN;
			if(!(xNew>xLo && xNew<xHi)){xNew=0.5*(xLo+xHi);}   // bisection safeguard
			if(Math.abs(xNew-x)<1e-12){return xNew;}
			x=xNew;
			if(xHi-xLo<1e-12){return 0.5*(xLo+xHi);}
		}
		return x;
	}

	/** h(lambda) = lambda / (e^lambda - 1); the per-set-cell score term. */
	private static double h(final double lam){
		if(lam<1e-9){return 1.0-0.5*lam;}
		if(lam>45){return 0;}
		return lam/Math.expm1(lam);
	}

	/** dL/dx at x=ln(mu): -2*mu*A + sum_t [ S*h(lam) - U*lam ]. */
	private static double deriv(final double x, final long[] S, final long[] U, final double A){
		final double mu=Math.exp(x);
		double sum=-2.0*mu*A;
		for(int t=0; t<MAX_TIER; t++){
			final long s=S[t], u=U[t];
			if(s==0 && u==0){continue;}
			final double lam=mu*Q_MASS[t];
			if(s>0){sum+=s*h(lam);}
			if(u>0){sum-=u*lam;}
		}
		return sum;
	}

	/** d2L/dx2: -2*mu*A + sum_t [ S*h*(1-lam-h) - U*lam ].  Always negative. */
	private static double deriv2(final double x, final long[] S, final long[] U, final double A){
		final double mu=Math.exp(x);
		double sum=-2.0*mu*A;
		for(int t=0; t<MAX_TIER; t++){
			final long s=S[t], u=U[t];
			if(s==0 && u==0){continue;}
			final double lam=mu*Q_MASS[t];
			if(s>0){
				final double hh=h(lam);
				sum+=s*hh*(1.0-lam-hh);
			}
			if(u>0){sum-=u*lam;}
		}
		return sum;
	}
}
