package ml;

/** Bounded affine sigmoid/tanh with absolute overshoot and backprop-only restoring force.
 * @author Brian Bushnell, Nilou */
public final class BungeeAffine extends Function{

	public BungeeAffine(final boolean tanh_, final double overshoot_){
		if(!Double.isFinite(overshoot_) || overshoot_<0 || overshoot_>0.01){
			throw new IllegalArgumentException("Affine overshoot must be finite and in [0,0.01]: "+overshoot_);
		}
		tanh=tanh_;
		overshoot=overshoot_;
		scale=tanh ? 1+overshoot : 1+2*overshoot;
		crossing=overshoot==0 ? Double.POSITIVE_INFINITY :
			tanh ? 0.5*(Math.log(2+overshoot)-Math.log(overshoot)) :
				Math.log1p(overshoot)-Math.log(overshoot);
		serializedName=(tanh ? "BTANH" : "BSIG")+":AFFINE:"+Double.toString(overshoot);
	}

	/** Parsing occurs once per exact token in CellNetParser; each token owns its overshoot. */
	public static BungeeAffine parse(final String token){
		if(token==null){throw new IllegalArgumentException("Missing affine activation token");}
		final String[] fields=token.split(":", -1);
		if(fields.length!=3 || (!fields[0].equals("BTANH") && !fields[0].equals("BSIG")) ||
				!fields[1].equals("AFFINE")){
			throw new IllegalArgumentException("Malformed affine activation: "+token);
		}
		return new BungeeAffine(fields[0].equals("BTANH"), Double.parseDouble(fields[2]));
	}

	@Override
	public double activate(final double x){
		assert(Double.isFinite(x)) : "Affine bungee activation requires finite pre-activation: "+x;
		if(overshoot==0){return base().activate(x);}
		return tanh ? scale*Math.tanh(x) : scale*sigmoid(x)-overshoot;
	}

	@Override
	public double derivativeX(final double x){
		assert(Double.isFinite(x)) : "Affine bungee derivative requires finite pre-activation: "+x;
		if(overshoot==0){return base().derivativeX(x);}
		final double b=tanh ? Math.tanh(x) : sigmoid(x);
		return scale*(tanh ? 1-b*b : b*(1-b));
	}

	@Override
	public double derivativeFX(final double fx){
		if(overshoot==0){return base().derivativeFX(fx);}
		throw new UnsupportedOperationException("Affine bungee derivative requires pre-activation; output is scaled and shifted");
	}

	@Override
	public double derivativeXFX(final double x, final double fx){
		return overshoot==0 ? base().derivativeXFX(x, fx) : derivativeX(x);
	}

	@Override
	public double restoringGradient(final double x){
		assert(Double.isFinite(x)) : "Affine restoring force requires finite pre-activation: "+x;
		if(overshoot==0){return 0;}
		final double z=Math.abs(x);
		if(z<=crossing){return 0;}
		final double k=tanh ? 2 : 1;
		final double q=overshoot/(k*scale);
		final double u=k*(z-crossing);
		final double change=Math.expm1(-u);
		final double overflow=overshoot*(1-q)*(-change)/(1+q*change);
		return Math.copySign(overflow, x);
	}

	private Function base(){return Function.getFunction(tanh ? Function.TANH : Function.SIG);}

	private static double sigmoid(final double x){
		if(x>=0){return 1/(1+Math.exp(-x));}
		final double e=Math.exp(x);return e/(1+e);
	}

	@Override
	public int type(){return tanh ? Function.BTANH : Function.BSIG;}
	@Override
	public String name(){return serializedName;}

	public final boolean tanh;
	public final double overshoot;
	private final double scale, crossing;
	private final String serializedName;
}
