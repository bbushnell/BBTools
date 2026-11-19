package ml;

/** Rectified linear unit; derivative at zero is zero.
 * @author Brian Bushnell, Yoimiya */
public final class ReLU extends Function{

	private ReLU(){}
	@Override
	public double activate(final double x){return StandardFunctions.relu(x);}
	@Override
	public double derivativeX(final double x){return StandardFunctions.reluDerivative(x);}
	@Override
	public double derivativeXFX(final double x, final double fx){return StandardFunctions.reluDerivative(x);}
	@Override
	public double derivativeFX(final double fx){return StandardFunctions.reluDerivative(fx);}
	@Override
	public int type(){return Function.RELU;}
	@Override
	public String name(){return "RELU";}

	static final ReLU instance=new ReLU();
}
