package ml;

/** Softplus with beta 1 and the PyTorch linear threshold 20.
 * @author Brian Bushnell, Yoimiya */
public final class Softplus extends Function{

	private Softplus(){}
	@Override
	public double activate(final double x){return StandardFunctions.softplus(x);}
	@Override
	public double derivativeX(final double x){return StandardFunctions.softplusDerivative(x);}
	@Override
	public double derivativeXFX(final double x, final double fx){return StandardFunctions.softplusDerivative(x);}
	@Override
	public double derivativeFX(final double fx){return fx>20 ? 1 : -Math.expm1(-fx);}
	@Override
	public int type(){return Function.SOFTPLUS;}
	@Override
	public String name(){return "SOFTPLUS";}

	static final Softplus instance=new Softplus();
}
