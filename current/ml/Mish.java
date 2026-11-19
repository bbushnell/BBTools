package ml;

/** Mish: x times tanh of softplus(x).
 * @author Brian Bushnell, Yoimiya */
public final class Mish extends Function{

	private Mish(){}
	@Override
	public double activate(final double x){return StandardFunctions.mish(x);}
	@Override
	public double derivativeX(final double x){return StandardFunctions.mishDerivative(x);}
	@Override
	public double derivativeXFX(final double x, final double fx){return StandardFunctions.mishDerivative(x);}
	@Override
	public double derivativeFX(final double fx){throw new UnsupportedOperationException("Mish is not invertible; differentiation requires x");}
	@Override
	public int type(){return Function.MISH;}
	@Override
	public String name(){return "MISH";}

	static final Mish instance=new Mish();
}
