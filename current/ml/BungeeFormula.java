package ml;

/** Immutable model-owned correction with its matching output-overflow potential.
 * Critical FORMULA tokens make old readers reject unsupported mathematics.
 * @author Brian Bushnell, Yoimiya */
public final class BungeeFormula extends Function{

	private BungeeFormula(final boolean tanh_, final Shape shape_, final double strength_,
			final double first_, final double second_){
		if(!Double.isFinite(strength_) || strength_<0){throw new IllegalArgumentException("Correction strength must be finite and nonnegative");}
		if(!Double.isFinite(first_) || !Double.isFinite(second_)){
			throw new IllegalArgumentException("Correction parameters must be finite");
		}
		if(shape_==Shape.RATIONAL && !scale(first_)){
			throw new IllegalArgumentException("Rational k must be in the supported positive FP32 scale range");
		}
		if(shape_==Shape.GATED && (first_<0 || first_>1e18 || !scale(second_))){
			throw new IllegalArgumentException("Gated c must be nonnegative and s a supported positive FP32 scale");
		}
		tanh=tanh_; shape=shape_; strength=strength_; first=first_; second=second_;
		base=Function.getFunction(tanh ? Function.TANH : Function.SIG);
		crossing=strength==0 ? Double.POSITIVE_INFINITY : findCrossing();
		gatedBaseline=shape==Shape.GATED && strength>0 ? gatedRemainder(crossing) : 0;
		serializedName=(tanh ? "BTANH" : "BSIG")+":FORMULA:"+shape.token+":"+Double.toString(strength)+
			(shape==Shape.RATIONAL || shape==Shape.GATED ? ":"+Double.toString(first) : "")+
			(shape==Shape.GATED ? ":"+Double.toString(second) : "");
		assert(strength==0 || Double.isFinite(crossing) && crossing>0) : "A positive correction needs a finite positive rail crossing for its restoring potential";
	}

	/** Parsing occurs once on model loading; inference never allocates tokens. */
	public static BungeeFormula parse(final String token){
		if(token==null){throw new IllegalArgumentException("Missing formula activation token");}
		final String[] fields=token.split(":", -1);
		if(fields.length<4 || (!fields[0].equals("BTANH") && !fields[0].equals("BSIG")) || !fields[1].equals("FORMULA")){
			throw new IllegalArgumentException("Malformed formula activation: "+token);
		}
		final Shape shape=Shape.parse(fields[2]);
		final int count=shape==Shape.GATED ? 6 : shape==Shape.RATIONAL ? 5 : 4;
		if(fields.length!=count){throw new IllegalArgumentException("Wrong number of correction parameters: "+token);}
		return new BungeeFormula(fields[0].equals("BTANH"), shape, Double.parseDouble(fields[3]),
			count>4 ? Double.parseDouble(fields[4]) : 0, count>5 ? Double.parseDouble(fields[5]) : 0);
	}

	@Override
	public double activate(final double x){
		assert(Double.isFinite(x)) : "Formula activation requires finite pre-activation before network propagation: "+x;
		if(strength==0){return base.activate(x);}
		return (tanh ? Math.tanh(x) : sigmoid(x))+strength*shape.value(x, first, second);
	}

	@Override
	public double derivativeX(final double x){
		assert(Double.isFinite(x)) : "Formula derivatives use pre-activation, not leaked output: "+x;
		if(strength==0){return base.derivativeX(x);}
		final double value=tanh ? Math.tanh(x) : sigmoid(x);
		return (tanh ? 1-value*value : value*(1-value))+strength*shape.derivative(x, first, second);
	}

	@Override
	public double derivativeFX(final double fx){
		if(strength==0){return base.derivativeFX(fx);}
		throw new UnsupportedOperationException("A formula correction needs pre-activation for its derivative");
	}

	@Override
	public double derivativeXFX(final double x, final double fx){
		return strength==0 ? base.derivativeXFX(x, fx) : derivativeX(x);
	}

	@Override
	public double restoringGradient(final double x){
		assert(Double.isFinite(x)) : "Restoring force requires finite pre-activation: "+x;
		final double z=Math.abs(x);
		if(strength==0 || z<=crossing){return 0;}
		final double k=tanh ? 2 : 1, e=Math.exp(-k*z);
		return Math.copySign(1, x)*(strength*shape.value(z, first, second)-k*e/(1+e));
	}

	/** Integral of restoringGradient, zero within the original activation rails. */
	public double restoringPotential(final double x){
		assert(Double.isFinite(x)) : "Restoring potential requires a finite integration endpoint: "+x;
		final double z=Math.abs(x), t=crossing;
		if(strength==0 || z<=t){return 0;}
		final double d=z-t, k=tanh ? 2 : 1, e=Math.exp(-k*t), q=e/(1+e);
		final double correction;
		switch(shape){
			case LINEAR: correction=.5*d*(d+2*t); break;
			case RATIONAL:{
				final double difference=d*(d+2*t), denominator=first*first+t*t, v=difference/denominator;
				correction=.5*((t*t/denominator)*difference+first*first*vMinusLog1p(v)); break;
			}
			case ASINH:{
				final double ht=Math.hypot(t, 1), hz=Math.hypot(z, 1), rootDifference=d*(d+2*t)/(hz+ht);
				correction=d*asinh(z)+t*Math.log1p((d+rootDifference)/(t+ht))-rootDifference; break;
			}
			case GATED:{
				final double quadratic=t>=first ? .5*d*(z+t) : .5*Math.max(0, z-first)*(z+first);
				correction=quadratic+gatedRemainder(z)-gatedBaseline; break;
			}
			default: throw new IllegalStateException("Unregistered correction primitive");
		}
		return Math.log1p(q*Math.expm1(-k*d))+strength*correction;
	}

	private double findCrossing(){
		final double logStrength=Math.log(strength), k=tanh ? 2 : 1;
		double lo=0, hi=1;
		while(logStrength+shape.logValue(hi, first, second)<Math.log(k)-k*hi-Math.log1p(Math.exp(-k*hi))){
			lo=hi; hi*=2;
			if(!Double.isFinite(hi)){throw new IllegalArgumentException("No finite correction rail bracket");}
		}
		if(lo==0){
			lo=hi/2;
			while(lo>0 && logStrength+shape.logValue(lo, first, second)>=Math.log(k)-k*lo-Math.log1p(Math.exp(-k*lo))){
				hi=lo;lo/=2;
			}
		}
		for(int i=0; i<80; i++){
			final double mid=lo+(hi-lo)/2;
			if(mid==0){return hi;}
			if(logStrength+shape.logValue(mid, first, second)>=Math.log(k)-k*mid-Math.log1p(Math.exp(-k*mid))){hi=mid;}
			else{lo=mid;}
		}
		return lo+(hi-lo)/2;
	}

	private double gatedRemainder(final double z){
		assert(shape==Shape.GATED) : "Dilogarithm primitive belongs only to the sigmoid-gated correction";
		final double u=(z-first)/second, r=Math.exp(-Math.abs(u)), log=Math.log1p(r), dilog=li2NegativeUnit(r);
		return second*z*log+second*second*(u>=0 ? -Math.PI*Math.PI/6-dilog : dilog);
	}

	private static double li2NegativeUnit(final double r){
		assert(r>=0 && r<=1) : "Reflected dilogarithm series requires an argument in [0,1]: "+r;
		final double y=r/(1+r), log=Math.log1p(r);
		double polynomial=1./(48*48);
		for(int n=47; n>=1; n--){polynomial=polynomial*y+1./(n*n);}
		return -y*polynomial-.5*log*log;
	}

	private static double vMinusLog1p(final double v){
		assert(v>=0) : "Positive-tail primitive requires an endpoint at or above its rail crossing: "+v;
		if(v>=.01){return v-Math.log1p(v);}
		double polynomial=-1./9;
		for(int n=8; n>=2; n--){polynomial=polynomial*v+(n%2==0 ? 1. : -1.)/n;}
		return v*v*polynomial;
	}

	private static double sigmoid(final double x){
		if(x>=0){return 1/(1+Math.exp(-x));}
		final double e=Math.exp(x);return e/(1+e);
	}

	private static double asinh(final double x){
		final double z=Math.abs(x);
		return Math.copySign(z>1e150 ? Math.log(z)+Math.log(2) : Math.log1p(z+z*z/(Math.hypot(z, 1)+1)), x);
	}

	private static boolean scale(final double value){return value>=1e-18 && value<=1e18;}
	@Override
	public int type(){return tanh ? Function.BTANH : Function.BSIG;}
	@Override
	public String name(){return serializedName;}
	public double railCrossing(){return crossing;}

	/** Each registration owns its forward, derivative and log-domain bracket term. */
	private enum Shape{
		LINEAR("linear"), RATIONAL("rational"), ASINH("asinh"), GATED("gated");
		Shape(final String token_){token=token_;}
		static Shape parse(final String token){
			for(final Shape value:values()){if(value.token.equals(token)){return value;}}
			throw new IllegalArgumentException("Unknown correction formula: "+token);
		}
		double value(final double x, final double a, final double b){
			switch(this){
				case LINEAR:return x;
				case RATIONAL:{final double r=x/Math.hypot(x, a);return x*r*r;}
				case ASINH:return asinh(x);
				case GATED:return x*sigmoid((Math.abs(x)-a)/b);
				default:throw new IllegalStateException("Missing forward correction");
			}
		}
		double derivative(final double x, final double a, final double b){
			switch(this){
				case LINEAR:return 1;
				case RATIONAL:{final double ratio=x/Math.hypot(x, a), r=ratio*ratio;return r*(3-2*r);}
				case ASINH:return 1/Math.hypot(x, 1);
				case GATED:{final double z=Math.abs(x), q=sigmoid((z-a)/b);return q==0 || q==1 ? q : q+(z/b)*q*(1-q);}
				default:throw new IllegalStateException("Missing correction derivative");
			}
		}
		double logValue(final double z, final double a, final double b){
			assert(z>0) : "Log-domain rail search requires positive pre-activation";
			switch(this){
				case LINEAR:return Math.log(z);
				case RATIONAL:return Math.log(z)+2*(Math.log(z)-Math.log(Math.hypot(z, a)));
				case ASINH:return Math.log(asinh(z));
				case GATED:{final double u=(z-a)/b;return Math.log(z)-Math.max(-u, 0)-Math.log1p(Math.exp(-Math.abs(u)));}
				default:throw new IllegalStateException("Missing correction rail bracket");
			}
		}
		final String token;
	}

	public final double strength;
	private final boolean tanh;
	private final Shape shape;
	private final double first, second, crossing, gatedBaseline;
	private final Function base;
	private final String serializedName;
}
