package ml;

/** Tanh approximation to GELU, serialized separately from the erf-based form.
 * @author Brian Bushnell, Yoimiya */
public final class GELUTanh extends Function{

	private GELUTanh(){}
	@Override
	public double activate(final double x){return StandardFunctions.geluTanh(x);}
	@Override
	public double derivativeX(final double x){return StandardFunctions.geluTanhDerivative(x);}
	@Override
	public double derivativeXFX(final double x, final double fx){return StandardFunctions.geluTanhDerivative(x);}
	@Override
	public double derivativeFX(final double fx){throw new UnsupportedOperationException("GELUTANH is not invertible; differentiation requires x");}
	@Override
	public int type(){return Function.GELUTANH;}
	@Override
	public String name(){return "GELUTANH";}

	static final GELUTanh instance=new GELUTanh();
}
