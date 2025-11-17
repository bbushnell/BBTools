package ml;

/**
 * Collection of activation functions and their derivatives for machine learning.
 * Includes sigmoid, tanh, swish, Gaussian, and specialized variants with optimized implementations.
 * All functions provide both forward evaluation and derivative calculation methods.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class Functions {
	
	/*--------------------------------------------------------------*/
	/*----------------           Sigmoid            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Standard sigmoid activation function.
	 * Computes 1/(1+e^(-x)) with output range (0,1).
	 * @param x Input value
	 * @return Sigmoid of x
	 */
	public static double sigmoid(double x) {return 1.0/(1.0+Math.exp(-x));}
	
	/**
	 * Computes the derivative of sigmoid function at x.
	 * Uses sigmoid(x) * (1 - sigmoid(x)) formula.
	 * @param x Input value
	 * @return Derivative of sigmoid at x
	 */
	static final double sigmoidDerivativeX(double x) {return sigmoidDerivativeFX(sigmoid(x));}
	
	/**
	 * Computes sigmoid derivative given sigmoid output.
	 * More efficient when sigmoid value is already available.
	 * @param fx Already computed sigmoid(x) value
	 * @return Derivative of sigmoid
	 */
	static final double sigmoidDerivativeFX(double fx) {return fx*(1-fx);}
	
	/**
	 * Computes sigmoid derivative given both input and sigmoid output.
	 * Uses the fx value for efficiency, ignoring x parameter.
	 *
	 * @param x Input value (unused)
	 * @param fx Already computed sigmoid(x) value
	 * @return Derivative of sigmoid
	 */
	static final double sigmoidDerivativeXFX(double x, double fx) {return fx*(1-fx);}
	
	/*--------------------------------------------------------------*/
	/*----------------       Extended Sigmoid       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Extended sigmoid activation function.
	 * Computes 2*sigmoid(x)-1 with output range (-1,1) instead of (0,1).
	 * @param x Input value
	 * @return Extended sigmoid of x
	 */
	public static double eSigmoid(double x) {
		return 2*sigmoid(x)-1;
	}
	
	/**
	 * Computes the derivative of extended sigmoid at x.
	 * Simply 2 times the standard sigmoid derivative.
	 * @param x Input value
	 * @return Derivative of extended sigmoid at x
	 */
	static final double eSigmoidDerivativeX(double x) {
		return 2*sigmoidDerivativeX(x);
	}
	
	/**
	 * Computes extended sigmoid derivative given function output.
	 * Uses optimized formula 0.5*(fx+1)*(2-fx-1) for efficiency.
	 * @param fx Already computed eSigmoid(x) value
	 * @return Derivative of extended sigmoid
	 */
	static final double eSigmoidDerivativeFX(double fx) {
		final double fx2=fx+1;
//		final double sfx=fx2*0.5;
//		final double d=2*sigmoidDerivativeFX(sfx); //This should be correct
		final double d2=0.5*fx2*(2-fx2); //This now matches the correct one and should be faster.
//		assert(d==d2) : fx+", "+fx2+", "+sfx+", "+d+", "+d2;
		return d2;
	}
	
	/**
	 * Computes extended sigmoid derivative given both input and output.
	 * Uses the fx value for efficiency, ignoring x parameter.
	 *
	 * @param x Input value (unused)
	 * @param fx Already computed eSigmoid(x) value
	 * @return Derivative of extended sigmoid
	 */
	static final double eSigmoidDerivativeXFX(double x, double fx) {
		return eSigmoidDerivativeFX(fx);//TODO: Change to x if x is faster
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             TanH             ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Hyperbolic tangent activation function with optimized implementation.
	 * Computes (e^x - e^(-x))/(e^x + e^(-x)) with clamping for extreme values.
	 * Returns -1 for x < -20 and 1 for x > 20 to prevent overflow.
	 *
	 * @param x Input value
	 * @return Hyperbolic tangent of x as float
	 */
	public static float tanh(double x) {
//		return (float)Math.tanh(x); //This is slower and gives the same results
		if(x<-20) {return -1;}
		if(x>20) {return 1;}
		
		double ex=Math.exp(x);
		double emx=Math.exp(-x);
		return (float)((ex-emx)/(ex+emx));
		
//		float ex=(float)Math.exp(x);This gives totally different results
//		float emx=(float)Math.exp(-x);
//		return ((ex-emx)/(ex+emx));
	}
	
	/**
	 * Computes the derivative of tanh function at x.
	 * Uses 1 - tanh^2(x) formula.
	 * @param x Input value
	 * @return Derivative of tanh at x
	 */
	static final double tanhDerivativeX(double x) {
		return tanhDerivativeFX(tanh(x));
	}
	
	/**
	 * Computes tanh derivative given tanh output.
	 * More efficient when tanh value is already available.
	 * @param fx Already computed tanh(x) value
	 * @return Derivative of tanh
	 */
	static final double tanhDerivativeFX(double fx) {
		return 1-fx*fx;
	}
	
	/**
	 * Computes tanh derivative given both input and tanh output.
	 * Uses the fx value for efficiency, ignoring x parameter.
	 *
	 * @param x Input value (unused)
	 * @param fx Already computed tanh(x) value
	 * @return Derivative of tanh
	 */
	static final double tanhDerivativeXFX(double x, double fx) {
		return 1-fx*fx;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Swish             ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Swish activation function.
	 * Computes x * sigmoid(x) = x / (1 + e^(-x)) for smooth, non-monotonic activation.
	 * @param x Input value
	 * @return Swish activation of x
	 */
	public static double swish(double x) {
//		double b=1;
//		return x*sigmoidD(b*x);
//		return x*sigmoidD(x);
		final double y=x/(1+Math.exp(-x)); //Maybe faster
//		return y<100000 ? y : Math.min(y, 100000+Math.log(y));
		return y;
	}
	
	/**
	 * Computes the derivative of swish function at x.
	 * Uses sigmoid value and swish output for calculation.
	 * @param x Input value
	 * @return Derivative of swish at x
	 */
	public static double swishDerivativeX(double x) {
		double sigx=sigmoid(x);
		double fx=x*sigx;
		//double fx=swish(x);
		return swishDerivativeFXSIGX(fx, sigx);
	}
	
	/**
	 * Computes swish derivative given only swish output.
	 * Currently unimplemented - throws RuntimeException.
	 *
	 * @param fx Already computed swish(x) value
	 * @return Derivative of swish
	 * @throws RuntimeException Always thrown as this method is unimplemented
	 */
	static final double swishDerivativeFX(double fx) {
		throw new RuntimeException("Unimplemented.");
	}
	
	/**
	 * Computes swish derivative given both input and swish output.
	 * Calculates sigmoid(x) internally and uses both values.
	 *
	 * @param x Input value
	 * @param fx Already computed swish(x) value
	 * @return Derivative of swish
	 */
	public static double swishDerivativeXFX(double x, double fx) {
		double sigx=sigmoid(x);
		return swishDerivativeFXSIGX(fx, sigx);
	}
	
	/**
	 * Computes swish derivative given swish output and sigmoid input.
	 * Uses formula fx + sigmoid(x) * (1 - fx) for efficient calculation.
	 *
	 * @param fx Already computed swish(x) value
	 * @param sigx Already computed sigmoid(x) value
	 * @return Derivative of swish
	 */
	public static double swishDerivativeFXSIGX(double fx, double sigx) {
		return fx+sigx*(1-fx);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            RSLog             ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Reverse-symmetric logarithm activation function.
	 * Computes log(|x|+1) while preserving sign: -rslog(-x) for negative x.
	 * Provides smooth activation that handles negative inputs gracefully.
	 *
	 * @param x Input value
	 * @return Reverse-symmetric log of x
	 */
	public static double rslog(double x) {
		if(x<0) {return -rslog(-x);}
		return Math.log(x+1);
	}
	
	/**
	 * Computes the derivative of rslog function at x.
	 * Uses 1/(|x|+1) formula, handling sign appropriately.
	 * @param x Input value
	 * @return Derivative of rslog at x
	 */
	static final double rslogDerivativeX(double x) {
		if(x<0) {x=-x;}
		return 1/(x+1);
	}
	
	/**
	 * Computes rslog derivative given rslog output.
	 * Contains assertion failure - likely unfinished implementation.
	 * Attempts to recover x from fx using inverse log operation.
	 *
	 * @param fx Already computed rslog(x) value
	 * @return Derivative of rslog
	 */
	static final double rslogDerivativeFX(double fx) {
		assert(false);
		if(fx<0) {fx=-fx;}
		//if fx=10 then log(x+1)=10 so e^10=x+1
		double x=Math.exp(fx)-1;
		return rslogDerivativeX(x);
	}
	
	/**
	 * Computes rslog derivative given both input and rslog output.
	 * Uses the x value for efficiency, ignoring fx parameter.
	 *
	 * @param x Input value
	 * @param fx Already computed rslog(x) value (unused)
	 * @return Derivative of rslog
	 */
	static final double rslogDerivativeXFX(double x, double fx) {
		return rslogDerivativeX(x);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Mirrored Sigmoid       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * X-axis offset for mirrored sigmoid function, controls peak height and width
	 */
	private static final double MSIG_X_OFFSET=5; //Bigger makes the peak higher, wider, and less sharp
	/**
	 * X-axis multiplier for mirrored sigmoid function, useful for 0-1 input ranges
	 */
	private static final double MSIG_X_MULT=2; //Useful if you design around a 0-1 range
	/**
	 * Y-axis multiplier for mirrored sigmoid, normalizes output to be slightly above 1
	 */
	private static final double MSIG_Y_MULT=1.0/sigmoid(MSIG_X_OFFSET);//Should be very slightly higher than 1.
	
	/**
	 * Mirrored sigmoid activation function.
	 * Creates a symmetric bell-shaped curve using two sigmoid functions.
	 * For x < 0: uses sigmoid(2x + 5), for x >= 0: uses sigmoid(-2x + 5).
	 * Provides peak activation at x=0 with symmetric decay.
	 *
	 * @param x Input value
	 * @return Mirrored sigmoid of x
	 */
	public static double mSig(double x) {
		final double offset=MSIG_X_OFFSET;
		final double xmult=MSIG_X_MULT;
		final double ymult=MSIG_Y_MULT;
		
		//sigmoid:  1.0/(1.0+Math.exp(-x));
		if(x<0) {
			//=1/(1+EXP(-(2*C7+$E$5)))
			double y=1.0/(1.0+Math.exp(-(xmult*x+offset)));
//			double y2=sigmoid(mult*x+offset);
//			assert(y2==y) : y+", "+y2+", "+x;
			return ymult*y;
		}else {
			double y=1.0/(1.0+Math.exp(xmult*x-offset));
//			double y2=sigmoid(-(mult*x-offset));
//			assert(y2==y) : y+", "+y2+", "+x;
			return ymult*y;
		}
	}
	
	/**
	 * Computes the derivative of mirrored sigmoid at x.
	 * Calculates mSig(x) and uses both values for derivative.
	 * @param x Input value
	 * @return Derivative of mirrored sigmoid at x
	 */
	static final double mSigDerivativeX(double x) {
		double fx=mSig(x);
		return mSigDerivativeXFX(x, fx);
	}
	
	/**
	 * Computes mirrored sigmoid derivative given only function output.
	 * Cannot be calculated without x value - throws RuntimeException.
	 *
	 * @param fx Already computed mSig(x) value
	 * @return Derivative of mirrored sigmoid
	 * @throws RuntimeException Always thrown as calculation requires x value
	 */
	static final double mSigDerivativeFX(double fx) {
		throw new RuntimeException("Cannot be calculated.");
	}
	
	/**
	 * Computes mirrored sigmoid derivative given both input and output.
	 * Uses x to determine sign and applies appropriate multiplier to sigmoid derivative.
	 *
	 * @param x Input value
	 * @param fx Already computed mSig(x) value
	 * @return Derivative of mirrored sigmoid
	 */
	static final double mSigDerivativeXFX(double x, double fx) {
		final double xmult=MSIG_X_MULT;
		double d=sigmoidDerivativeFX(fx);
		if(x<0){return xmult*d;}
		return -xmult*d;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------  Extended Mirrored Sigmoid   ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Extended mirrored sigmoid activation function.
	 * Computes 2*mSig(x)-1 to shift output range from (0,1) to (-1,1).
	 * @param x Input value
	 * @return Extended mirrored sigmoid of x
	 */
	public static double emSig(double x) {
		return 2*mSig(x)-1;
	}
	
	/**
	 * Computes the derivative of extended mirrored sigmoid at x.
	 * Simply 2 times the mirrored sigmoid derivative.
	 * @param x Input value
	 * @return Derivative of extended mirrored sigmoid at x
	 */
	static final double emSigDerivativeX(double x) {
		return 2*mSigDerivativeX(x);
	}
	
	/**
	 * Computes extended mirrored sigmoid derivative given function output.
	 * Cannot be calculated without x value - throws RuntimeException.
	 *
	 * @param fx Already computed emSig(x) value
	 * @return Derivative of extended mirrored sigmoid
	 * @throws RuntimeException Always thrown as calculation requires x value
	 */
	static final double emSigDerivativeFX(double fx) {
		throw new RuntimeException("Cannot be calculated.");
	}
	
	/**
	 * Computes extended mirrored sigmoid derivative given input and output.
	 * Uses x value to recalculate mirrored sigmoid derivative.
	 *
	 * @param x Input value
	 * @param fx Already computed emSig(x) value (unused)
	 * @return Derivative of extended mirrored sigmoid
	 */
	static final double emSigDerivativeXFX(double x, double fx) {
		return 2*mSigDerivativeX(x); //Possibly slow, but simple.
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Gaussian           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Gaussian bell curve activation function.
	 * Computes e^(-x^2) providing symmetric activation centered at x=0.
	 * Output range is (0,1] with maximum at x=0.
	 *
	 * @param x Input value
	 * @return Gaussian bell curve of x
	 */
	public static double bell(double x) {
		return Math.exp(-(x*x));
	}
	
	/**
	 * Computes the derivative of Gaussian bell function at x.
	 * Uses -2x * e^(-x^2) formula, equivalent to -2x * fx.
	 * @param x Input value
	 * @return Derivative of Gaussian bell at x
	 */
	static final double bellDerivativeX(double x) {
		return -2*x*Math.exp(-x*x); //i.e. -2*x*fx
	}
	
	/**
	 * Computes Gaussian bell derivative given only function output.
	 * Cannot be calculated without x value - throws RuntimeException.
	 *
	 * @param fx Already computed bell(x) value
	 * @return Derivative of Gaussian bell
	 * @throws RuntimeException Always thrown as calculation requires x value
	 */
	static final double bellDerivativeFX(double fx) {
		throw new RuntimeException("Cannot be calculated.");
	}
	
	/**
	 * Computes Gaussian bell derivative given both input and output.
	 * Uses efficient -2x * fx formula where fx is the bell function value.
	 *
	 * @param x Input value
	 * @param fx Already computed bell(x) value
	 * @return Derivative of Gaussian bell
	 */
	static final double bellDerivativeXFX(double x, double fx) {
		return -2*x*fx;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Other             ----------------*/
	/*--------------------------------------------------------------*/
	
	//For multiple tests of a single output neuron
	/**
	 * Computes mean squared error between target and actual values.
	 * Used for evaluating multiple predictions against expected outputs.
	 * Calculates sum of squared differences divided by array length.
	 *
	 * @param target Expected values array
	 * @param actual Predicted values array
	 * @return Mean squared error between the arrays
	 */
	public static double mse(float[] target, float[] actual) {
		assert(target.length==actual.length);
		double sum=0;
		final int len=target.length;
		for(int i=0; i<len; i++) {
			float t=target[i], a=actual[i];
			float dif=t-a;
			sum+=(dif*dif);
		}
		return sum/len;
	}

}
