package cardinality;

/**
 * AvllMLE: analytic maximum-likelihood extractor for ArithmeticVariableLogLog.
 * <p>
 * An AVLL register is a 3-cell tail anchored to its own max: the exponent
 * records the highest NLZ observed (exactly), and the two history bits are
 * presence cells at ranks max-1 and max-2.  The update rule is a monotone
 * OR-shift, so the state is a set-function and the likelihood factorizes:
 * <p>
 *   state 0:            exp(-mu * A(g))                    [nothing at rank >= g]
 *   exact max at r:     exp(-lam) * (1 - exp(-lam)),  lam = mu * q(r)
 *                       [none above r] * [>=1 at r];  A(r+1) == q(r) == 2^-(r+1)
 *   history cell at r': Bernoulli(1 - exp(-mu * q(r')))
 *   clamped (exp==18):  1 - exp(-mu * A(r))                [max >= r, censored]
 * <p>
 * where mu = n / modBuckets, q(r) = 2^-(r+1) is the NLZ mass, A(r) = 2^-r is
 * the mass at rank >= r, and g = globalNLZ.  History-cell validity follows the
 * state map: both cells for stored exp 2..13, the upper cell only for exp 1,
 * none for exp >= 14.  Everything aggregates into per-rank set/unset counts
 * plus one exposure scalar, so a Newton solve costs O(64) per iteration.
 * No calibration tables.
 *
 * @author Amber
 * @date July 2026
 */
public final class AvllMLE {

	/** q[r] = P(NLZ == r) = 2^-(r+1). */
	static final double[] Q=new double[64];
	static{
		for(int r=0; r<64; r++){Q[r]=Math.pow(2.0, -(r+1));}
	}

	/** Mass at rank >= r (r may be negative). */
	private static double massAtOrAbove(int r){
		if(r<=0){return 1.0;}
		if(r>=64){return 0;}
		return Math.pow(2.0, -r);
	}

	/**
	 * Maximum-likelihood cardinality estimate from raw AVLL state.
	 * @param avll the sketch
	 * @return estimated number of distinct elements added
	 */
	public static double estimate(final ArithmeticVariableLogLog avll){
		final int B=avll.actualBuckets();
		final int g=avll.getGlobalNLZ();

		final long[] S=new long[64];
		final long[] U=new long[64];
		double exposure=0;   // sum of "none-above" masses; term is -mu*exposure
		long touched=0;

		for(int i=0; i<B; i++){
			final int reg=avll.getReg(i);
			if(reg==0){
				exposure+=massAtOrAbove(g);   // no element with absNlz >= g
				continue;
			}
			touched++;
			final int es=ArithmeticVariableLogLog.stateToExp(reg);
			final int hist=ArithmeticVariableLogLog.stateToHist(reg);
			final int top=es+g-1;             // absolute NLZ of the register max
			if(top<0 || top>=64){exposure+=massAtOrAbove(g); continue;} // degenerate
			if(es>=18){
				// Clamped: max >= top (censored).  Model as a set cell of mass
				// A(top) = 2^-top = q(top-1).
				if(top>=1){S[top-1]++;}
				continue;
			}
			// Exact max: none above top (mass q(top)) and >=1 at top.
			exposure+=Q[top];                  // A(top+1) == q(top)
			S[top]++;
			// History cells
			if(es>=2 && es<=13){
				final int r1=top-1, r0=top-2;
				if(r1>=0){if(((hist>>>1)&1)!=0){S[r1]++;}else{U[r1]++;}}
				if(r0>=0){if((hist&1)!=0){S[r0]++;}else{U[r0]++;}}
			}else if(es==1){
				final int r1=top-1;
				if(r1>=0){if(((hist>>>1)&1)!=0){S[r1]++;}else{U[r1]++;}}
			}
		}

		if(touched==0){return 0;}

		final double x=solve(S, U, exposure);
		return B*Math.exp(x);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Newton Solver         ----------------*/
	/*--------------------------------------------------------------*/

	private static double solve(final long[] S, final long[] U, final double A){
		double xLo=-34, xHi=48;
		double x=0;
		if(deriv(xLo, S, U, A)<=0){return xLo;}
		for(int iter=0; iter<200; iter++){
			final double d1=deriv(x, S, U, A);
			if(d1>0){xLo=x;}else{xHi=x;}
			final double d2=deriv2(x, S, U, A);
			double xNew=(d2<0) ? x-d1/d2 : Double.NaN;
			if(!(xNew>xLo && xNew<xHi)){xNew=0.5*(xLo+xHi);}
			if(Math.abs(xNew-x)<1e-12){return xNew;}
			x=xNew;
			if(xHi-xLo<1e-12){return 0.5*(xLo+xHi);}
		}
		return x;
	}

	private static double h(final double lam){
		if(lam<1e-9){return 1.0-0.5*lam;}
		if(lam>45){return 0;}
		return lam/Math.expm1(lam);
	}

	private static double deriv(final double x, final long[] S, final long[] U, final double A){
		final double mu=Math.exp(x);
		double sum=-mu*A;
		for(int r=0; r<64; r++){
			final long s=S[r], u=U[r];
			if(s==0 && u==0){continue;}
			final double lam=mu*Q[r];
			if(s>0){sum+=s*h(lam);}
			if(u>0){sum-=u*lam;}
		}
		return sum;
	}

	private static double deriv2(final double x, final long[] S, final long[] U, final double A){
		final double mu=Math.exp(x);
		double sum=-mu*A;
		for(int r=0; r<64; r++){
			final long s=S[r], u=U[r];
			if(s==0 && u==0){continue;}
			final double lam=mu*Q[r];
			if(s>0){
				final double hh=h(lam);
				sum+=s*hh*(1.0-lam-hh);
			}
			if(u>0){sum-=u*lam;}
		}
		return sum;
	}
}
