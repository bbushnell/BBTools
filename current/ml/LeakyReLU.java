package ml;

/** Leaky rectified linear unit with negative slope .01, including at zero.
 * @author Brian Bushnell, Yoimiya */
public final class LeakyReLU extends Function{

	private LeakyReLU(){}
	@Override
	public double activate(final double x){return StandardFunctions.leakyRelu(x);}
	@Override
	public double derivativeX(final double x){return StandardFunctions.leakyReluDerivative(x);}
	@Override
	public double derivativeXFX(final double x, final double fx){return StandardFunctions.leakyReluDerivative(x);}
	@Override
	public double derivativeFX(final double fx){return StandardFunctions.leakyReluDerivative(fx);}
	@Override
	public int type(){return Function.LEAKYRELU;}
	@Override
	public String name(){return "LEAKYRELU";}

	static final LeakyReLU instance=new LeakyReLU();
}
