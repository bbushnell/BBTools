package ml;

public class Functions {
	
	/*--------------------------------------------------------------*/
	/*----------------           Sigmoid            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static double sigmoid(double x) {return 1.0/(1.0+Math.exp(-x));}
	
	static final double sigmoidDerivativeX(double x) {return sigmoidDerivativeFX(sigmoid(x));}
	
	static final double sigmoidDerivativeFX(double fx) {return fx*(1-fx);}
	
	static final double sigmoidDerivativeXFX(double x, double fx) {return fx*(1-fx);}
	
	/*--------------------------------------------------------------*/
	/*----------------       Extended Sigmoid       ----------------*/
	/*--------------------------------------------------------------*/
	
	public static double eSigmoid(double x) {
		return 2*sigmoid(x)-1;
	}
	
	static final double eSigmoidDerivativeX(double x) {
		return 2*sigmoidDerivativeX(x);
	}
	
	static final double eSigmoidDerivativeFX(double fx) {
		final double fx2=fx+1;
//		final double sfx=fx2*0.5;
//		final double d=2*sigmoidDerivativeFX(sfx); //This should be correct
		final double d2=0.5*fx2*(2-fx2); //This now matches the correct one and should be faster.
//		assert(d==d2) : fx+", "+fx2+", "+sfx+", "+d+", "+d2;
		return d2;
	}
	
	static final double eSigmoidDerivativeXFX(double x, double fx) {
		return eSigmoidDerivativeFX(fx);//TODO: Change to x if x is faster
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             TanH             ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	static final double tanhDerivativeX(double x) {
		return tanhDerivativeFX(tanh(x));
	}
	
	static final double tanhDerivativeFX(double fx) {
		return 1-fx*fx;
	}
	
	static final double tanhDerivativeXFX(double x, double fx) {
		return 1-fx*fx;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Swish             ----------------*/
	/*--------------------------------------------------------------*/
	
	public static double swish(double x) {
//		double b=1;
//		return x*sigmoidD(b*x);
//		return x*sigmoidD(x);
		final double y=x/(1+Math.exp(-x)); //Maybe faster
//		return y<100000 ? y : Math.min(y, 100000+Math.log(y));
		return y;
	}
	
	public static double swishDerivativeX(double x) {
		double sigx=sigmoid(x);
		double fx=x*sigx;
		//double fx=swish(x);
		return swishDerivativeFXSIGX(fx, sigx);
	}
	
	static final double swishDerivativeFX(double fx) {
		throw new RuntimeException("Unimplemented.");
	}
	
	public static double swishDerivativeXFX(double x, double fx) {
		double sigx=sigmoid(x);
		return swishDerivativeFXSIGX(fx, sigx);
	}
	
	public static double swishDerivativeFXSIGX(double fx, double sigx) {
		return fx+sigx*(1-fx);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            RSLog             ----------------*/
	/*--------------------------------------------------------------*/
	
	public static double rslog(double x) {
		if(x<0) {return -rslog(-x);}
		return Math.log(x+1);
	}
	
	static final double rslogDerivativeX(double x) {
		if(x<0) {x=-x;}
		return 1/(x+1);
	}
	
	static final double rslogDerivativeFX(double fx) {
		assert(false);
		if(fx<0) {fx=-fx;}
		//if fx=10 then log(x+1)=10 so e^10=x+1
		double x=Math.exp(fx)-1;
		return rslogDerivativeX(x);
	}
	
	static final double rslogDerivativeXFX(double x, double fx) {
		return rslogDerivativeX(x);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Mirrored Sigmoid       ----------------*/
	/*--------------------------------------------------------------*/

	private static final double MSIG_X_OFFSET=5; //Bigger makes the peak higher, wider, and less sharp
	private static final double MSIG_X_MULT=2; //Useful if you design around a 0-1 range
	private static final double MSIG_Y_MULT=1.0/sigmoid(MSIG_X_OFFSET);//Should be very slightly higher than 1.
	
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
	
	static final double mSigDerivativeX(double x) {
		double fx=mSig(x);
		return mSigDerivativeXFX(x, fx);
	}
	
	static final double mSigDerivativeFX(double fx) {
		throw new RuntimeException("Cannot be calculated.");
	}
	
	static final double mSigDerivativeXFX(double x, double fx) {
		final double xmult=MSIG_X_MULT;
		double d=sigmoidDerivativeFX(fx);
		if(x<0){return xmult*d;}
		return -xmult*d;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------  Extended Mirrored Sigmoid   ----------------*/
	/*--------------------------------------------------------------*/
	
	public static double emSig(double x) {
		return 2*mSig(x)-1;
	}
	
	static final double emSigDerivativeX(double x) {
		return 2*mSigDerivativeX(x);
	}
	
	static final double emSigDerivativeFX(double fx) {
		throw new RuntimeException("Cannot be calculated.");
	}
	
	static final double emSigDerivativeXFX(double x, double fx) {
		return 2*mSigDerivativeX(x); //Possibly slow, but simple.
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Gaussian           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static double bell(double x) {
		return Math.exp(-(x*x));
	}
	
	static final double bellDerivativeX(double x) {
		return -2*x*Math.exp(-x*x); //i.e. -2*x*fx
	}
	
	static final double bellDerivativeFX(double fx) {
		throw new RuntimeException("Cannot be calculated.");
	}
	
	static final double bellDerivativeXFX(double x, double fx) {
		return -2*x*fx;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Other             ----------------*/
	/*--------------------------------------------------------------*/
	
	//For multiple tests of a single output neuron
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
