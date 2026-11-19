package ml;

/** Tanh with a linear leak and an independent, output-range restoring gradient.
 * @author Brian Bushnell, Yoimiya */
public final class BTanh extends Function{

	private BTanh(){epsilon=EPSILON; serializedName="BTANH";}

	/** Immutable per-model leak; a parameterized token cannot load as legacy .01. */
	public BTanh(final double epsilon_){
		if(!Double.isFinite(epsilon_) || epsilon_<0){throw new IllegalArgumentException("BTANH epsilon must be finite and nonnegative: "+epsilon_);}
		epsilon=epsilon_;
		serializedName="BTANH:"+Double.toString(epsilon);
	}

	@Override
	public double activate(final double x){
		assert(Double.isFinite(x)) : "Bungee activation requires finite pre-activation to preserve finite training state: "+x;
		return epsilon==0 ? Tanh.instance.activate(x) : Math.tanh(x)+epsilon*x;
	}

	@Override
	public double derivativeX(final double x){
		assert(Double.isFinite(x)) : "Bungee derivative requires finite pre-activation; output-only tanh algebra is invalid after the leak: "+x;
		if(epsilon==0){return Tanh.instance.derivativeX(x);}
		final double base=Math.tanh(x);
		return (1-base*base)+epsilon;
	}

	@Override
	public double derivativeFX(final double fx){
		if(epsilon==0){return Tanh.instance.derivativeFX(fx);}
		throw new UnsupportedOperationException("BTANH derivative requires pre-activation; leaked output is not tanh(x)");
	}

	@Override
	public double derivativeXFX(final double x, final double fx){
		return epsilon==0 ? Tanh.instance.derivativeXFX(x, fx) : derivativeX(x);
	}

	@Override
	public double restoringGradient(final double x){
		if(epsilon==0){return 0;}
		final double out=activate(x);
		return out>1 ? out-1 : out< -1 ? out+1 : 0;
	}

	@Override
	public int type(){return Function.BTANH;}
	@Override
	public String name(){return serializedName;}

	public final double epsilon;
	private final String serializedName;

	/** Brian's BF16 training choice: epsilon*x near saturation is several ULPs
	 * at the unit rail (spacing above 1 is 2^-7), giving a robust visible leak.
	 * Do not silently reduce this to .001; that choice is reserved for a separately
	 * specified pure-FP32 run. Sub-ULP increments do not universally round to zero. */
	public static final double EPSILON=0.01;
	static final BTanh instance=new BTanh();
}
