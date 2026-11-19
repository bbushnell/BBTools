package ml;

/** Gaussian error linear unit using the erf-based Gaussian CDF.
 * @author Brian Bushnell, Yoimiya */
public final class GELU extends Function{

	private GELU(){}
	@Override
	public double activate(final double x){return StandardFunctions.gelu(x);}
	@Override
	public double derivativeX(final double x){return StandardFunctions.geluDerivative(x);}
	@Override
	public double derivativeXFX(final double x, final double fx){return StandardFunctions.geluDerivative(x);}
	@Override
	public double derivativeFX(final double fx){throw new UnsupportedOperationException("GELU is not invertible; differentiation requires x");}
	@Override
	public int type(){return Function.GELU;}
	@Override
	public String name(){return "GELU";}

	static final GELU instance=new GELU();
}
