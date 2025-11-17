package ml;

/**
 * Sigmoid activation function implementation for machine learning applications.
 * Provides the logistic sigmoid function f(x) = 1/(1+e^(-x)) and its derivatives.
 * Used as an activation function in neural networks with output range (0,1).
 *
 * @author Brian Bushnell
 * @date 2014
 */
public class Sigmoid extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private Sigmoid() {}

	@Override
	public double activate(double x) {return Functions.sigmoid(x);}

	@Override
	public double derivativeX(double x) {return Functions.sigmoidDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.sigmoidDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.sigmoidDerivativeFX(fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** String identifier for the sigmoid function type */
	static final String name="SIG";
	/** Integer type identifier for the sigmoid function */
	static final int type=Function.toType(name, true);
	/** Singleton instance of the sigmoid function implementation */
	static final Sigmoid instance=new Sigmoid();

}
