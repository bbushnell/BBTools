package ml;

/** Standard hidden activations used for matched Java/PyTorch comparisons.
 * LeakyReLU uses slope .01, ELU alpha 1, and Softplus beta 1/threshold 20.
 * GELU uses the Gaussian CDF; GELUTANH is the distinct tanh approximation.
 * All methods allocate nothing. Nonfinite inputs propagate except where a
 * limiting value is explicit. Existing Functions arithmetic is unchanged.
 * @author Brian Bushnell, Yoimiya */
public final class StandardFunctions{

	public static double relu(final double x){return Math.max(0, x);}
	/** PyTorch's zero-point subgradient is zero. */
	public static double reluDerivative(final double x){return x>0 ? 1 : x<=0 ? 0 : Double.NaN;}
	public static double leakyRelu(final double x){return x>0 ? x : .01*x;}
	/** PyTorch selects the negative branch at zero. */
	public static double leakyReluDerivative(final double x){return x>0 ? 1 : x<=0 ? .01 : Double.NaN;}
	public static double elu(final double x){return x>0 ? x : Math.expm1(x);}
	public static double eluDerivative(final double x){return x>0 ? 1 : Math.exp(x);}
	public static double softplus(final double x){return x>20 ? x : Math.log1p(Math.exp(x));}
	public static double softplusDerivative(final double x){return x>20 ? 1 : Functions.sigmoid(x);}

	public static double gelu(final double x){
		if(x==Double.POSITIVE_INFINITY){return x;}
		if(x==Double.NEGATIVE_INFINITY){return -0.0;}
		return .5*x*(1+erf(x*INV_SQRT_TWO));
	}
	public static double geluDerivative(final double x){
		if(x==Double.POSITIVE_INFINITY){return 1;}
		if(x==Double.NEGATIVE_INFINITY){return 0;}
		return .5*(1+erf(x*INV_SQRT_TWO))+x*INV_SQRT_TWO_PI*Math.exp(-.5*x*x);
	}
	public static double geluTanh(final double x){
		// Beyond this point tanh has rounded to its endpoint in double arithmetic.
		if(x>20){return x;}
		if(x<-20){return -0.0;}
		return .5*x*(1+Math.tanh(SQRT_TWO_OVER_PI*x*(1+.044715*x*x)));
	}
	public static double geluTanhDerivative(final double x){
		if(x>20){return 1;}
		if(x<-20){return 0;}
		final double x2=x*x, t=Math.tanh(SQRT_TWO_OVER_PI*x*(1+.044715*x2));
		return .5*(1+t)+.5*x*(1-t*t)*SQRT_TWO_OVER_PI*(1+.134145*x2);
	}
	public static double mish(final double x){
		if(x==Double.NEGATIVE_INFINITY){return -0.0;}
		return x*Math.tanh(softplus(x));
	}
	public static double mishDerivative(final double x){
		if(x>20){return 1;}
		if(x==Double.NEGATIVE_INFINITY){return 0;}
		final double t=Math.tanh(softplus(x));
		return t+x*Functions.sigmoid(x)*(1-t*t);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Error function         ----------------*/
	/*--------------------------------------------------------------*/

	/* Rational coefficients and range reduction adapted from fdlibm s_erf.c:
	 * https://netlib.org/fdlibm/s_erf.c
	 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
	 * Developed at SunSoft, a Sun Microsystems, Inc. business.
	 * Permission to use, copy, modify, and distribute this software is freely
	 * granted, provided that this notice is preserved.
	 */
	/** Java 8-compatible erf for the exact Gaussian-CDF GELU, not a tanh surrogate.
	 * The rational pieces avoid a data-dependent iterative integral in each neuron. */
	public static double erf(final double x){
		final double a=Math.abs(x);
		if(Double.isNaN(x)){return x;}
		if(a<0.84375){
			if(a<0x1.0p-28){return x+0.128379167095512586316*x;}
			final double z=x*x;
			return x+x*polynomial(z, ERF_P)/polynomial(z, ERF_Q);
		}
		if(a<1.25){
			final double z=a-1, y=0.845062911510467529297+polynomial(z, ERF_PA)/polynomial(z, ERF_QA);
			return Math.copySign(y, x);
		}
		if(a>=6){return Math.copySign(1, x);}
		final double s=1/(a*a);
		final double ratio=a<2.857142857142857 ? polynomial(s, ERF_RA)/polynomial(s, ERF_SA)
			: polynomial(s, ERF_RB)/polynomial(s, ERF_SB);
		final double z=Double.longBitsToDouble(Double.doubleToRawLongBits(a)&0xffffffff00000000L);
		final double tail=Math.exp(-z*z-.5625)*Math.exp((z-a)*(z+a)+ratio)/a;
		return x>=0 ? 1-tail : tail-1;
	}
	private static double polynomial(final double x, final double[] coefficients){
		assert(coefficients.length>0) : "The erf rational pieces require a nonempty coefficient polynomial";
		double value=coefficients[coefficients.length-1];
		for(int i=coefficients.length-2; i>=0; i--){value=coefficients[i]+x*value;}
		return value;
	}
	private StandardFunctions(){}
	private static final double INV_SQRT_TWO=0.70710678118654752440;
	private static final double INV_SQRT_TWO_PI=0.39894228040143267794;
	private static final double SQRT_TWO_OVER_PI=0.79788456080286535588;
	private static final double[] ERF_P={.128379167095512558561, -.325042107247001499370,
		-.0284817495755985104766, -.00577027029648944159157, -.0000237630166566501626084};
	private static final double[] ERF_Q={1, .397917223959155352819, .0650222499887672944485,
		.00508130628187576562776, .000132494738004321644526, -.00000396022827877536812320};
	private static final double[] ERF_PA={-.00236211856075265944077, .414856118683748331666,
		-.372207876035701323847, .318346619901161753674, -.110894694282396677476,
		.0354783043256182359371, -.00216637559486879084300};
	private static final double[] ERF_QA={1, .106420880400844228286, .540397917702171048937,
		.0718286544141962662868, .126171219808761642112, .0136370839120290507362, .0119844998467991074170};
	private static final double[] ERF_RA={-.00986494403484714822705, -.693858572707181764372,
		-10.5586262253232909814, -62.3753324503260060396, -162.396669462573470355,
		-184.605092906711035994, -81.2874355063065934246, -9.81432934416914548592};
	private static final double[] ERF_SA={1, 19.6512716674392571292, 137.657754143519042600,
		434.565877475229228821, 645.387271733267880336, 429.008140027567833386,
		108.635005541779435134, 6.57024977031928170135, -.0604244152148580987438};
	private static final double[] ERF_RB={-.00986494292470009928597, -.799283237680523006574,
		-17.7579549177547519889, -160.636384855821916062, -637.566443368389627722,
		-1025.09513161107724954, -483.519191608651397019};
	private static final double[] ERF_SB={1, 30.3380607434824582924, 325.792512996573918826,
		1536.72958608443695994, 3199.85821950859553908, 2553.05040643316442583,
		474.528541206955367215, -22.4409524465858183362};
}
