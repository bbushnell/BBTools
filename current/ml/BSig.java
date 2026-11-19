package ml;

/** Sigmoid with a linear leak and restoring force below zero or above one.
 * @author Brian Bushnell, Yoimiya */
public final class BSig extends Function{

	private BSig(){epsilon=EPSILON; serializedName="BSIG";}

	/** Immutable per-model leak; a parameterized token cannot load as legacy .01. */
	public BSig(final double epsilon_){
		if(!Double.isFinite(epsilon_) || epsilon_<0){throw new IllegalArgumentException("BSIG epsilon must be finite and nonnegative: "+epsilon_);}
		epsilon=epsilon_;
		serializedName="BSIG:"+Double.toString(epsilon);
	}

	@Override
	public double activate(final double x){
		assert(Double.isFinite(x)) : "Bungee activation requires finite pre-activation to preserve finite training state: "+x;
		return epsilon==0 ? Sigmoid.instance.activate(x) : 1/(1+Math.exp(-x))+epsilon*x;
	}

	@Override
	public double derivativeX(final double x){
		assert(Double.isFinite(x)) : "Bungee derivative requires finite pre-activation; leaked output is not a sigmoid probability: "+x;
		if(epsilon==0){return Sigmoid.instance.derivativeX(x);}
		final double base=1/(1+Math.exp(-x));
		return base*(1-base)+epsilon;
	}

	@Override
	public double derivativeFX(final double fx){
		if(epsilon==0){return Sigmoid.instance.derivativeFX(fx);}
		throw new UnsupportedOperationException("BSIG derivative requires pre-activation; leaked output is not sigmoid(x)");
	}

	@Override
	public double derivativeXFX(final double x, final double fx){
		return epsilon==0 ? Sigmoid.instance.derivativeXFX(x, fx) : derivativeX(x);
	}

	@Override
	public double restoringGradient(final double x){
		if(epsilon==0){return 0;}
		final double out=activate(x);
		return out>1 ? out-1 : out<0 ? out : 0;
	}

	@Override
	public int type(){return Function.BSIG;}
	@Override
	public String name(){return serializedName;}

	public final double epsilon;
	private final String serializedName;

	/** Same fixed BF16-compatible leak as BTanh; see its precision rationale. */
	public static final double EPSILON=BTanh.EPSILON;
	static final BSig instance=new BSig();
}
