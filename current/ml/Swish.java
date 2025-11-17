package ml;

/**
 * Swish activation function implementation for neural network architectures.
 * Implements the swish(x) = x * sigmoid(x) activation function, providing smoother
 * gradients and improved performance compared to traditional activation functions like ReLU.
 *
 * @author Brian Bushnell
 * @date 2024
 */
public class Swish extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Private constructor preventing external instantiation - use static instance
	 */
	private Swish() {}

	@Override
	public double activate(double x) {return Functions.swish(x);}

	@Override
	public double derivativeX(double x) {return Functions.swishDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.swishDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.swishDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** String identifier for the Swish activation function */
	static final String name="SWISH";
	/** Numeric type identifier computed from the function name */
	static final int type=Function.toType(name, true);
	/** Singleton instance of the Swish activation function */
	static final Swish instance=new Swish();

}
