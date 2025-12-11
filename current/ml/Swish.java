package ml;

/**
 * Swish activation function implementation, swish(x) = x * sigmoid(x).
 * Provides smoother gradients than ReLU and is used in newer neural architectures.
 * Singleton implementation; use static instance rather than constructing directly.
 *
 * @author Brian Bushnell
 * @date 2024
 */
public class Swish extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private Swish() {}

	/**
	 * Computes the Swish activation for the given input.
	 * @param x Input value
	 * @return swish(x) = x * sigmoid(x)
	 */
	@Override
	public double activate(double x) {return Functions.swish(x);}

	/**
	 * Computes the derivative of Swish with respect to x.
	 * @param x Input value
	 * @return d/dx swish(x)
	 */
	@Override
	public double derivativeX(double x) {return Functions.swishDerivativeX(x);}

	/**
	 * Calculates the derivative of Swish given the precomputed activation value.
	 * Use when swish(x) has already been computed to avoid recomputation.
	 * @param fx Pre-computed Swish function output
	 * @return Derivative of Swish function
	 */
	@Override
	public double derivativeFX(double fx) {return Functions.swishDerivativeFX(fx);}

	/**
	 * Calculates the derivative using both input and precomputed output values.
	 * Convenient when both x and swish(x) are available.
	 *
	 * @param x Original input value
	 * @param fx Pre-computed Swish function output
	 * @return Derivative of Swish function
	 */
	@Override
	public double derivativeXFX(double x, double fx) {return Functions.swishDerivativeXFX(x, fx);}

	/** Returns the numeric type identifier for this activation function.
	 * @return Numeric type constant for Swish function */
	@Override
	public int type() {return type;}

	/** Returns the string name identifier for this activation function.
	 * @return String name "SWISH" */
	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static final String name="SWISH";
	static final int type=Function.toType(name, true);
	static final Swish instance=new Swish();

}
