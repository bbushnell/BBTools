package ml;

/**
 * Swish activation function implementation, swish(x) = x * sigmoid(x).
 * Obtain the singleton through {@link Function#getFunction(int)}.
 *
 * @author Brian Bushnell
 * @date 2024
 */
public class Swish extends Function{
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates the shared stateless activation instance. */
	private Swish(){}

	/**
	 * Computes the Swish activation for the given input.
	 * @param x Input value
	 * @return swish(x) = x * sigmoid(x)
	 */
	@Override
	public double activate(double x){return Functions.swish(x);}

	/**
	 * Computes the derivative of Swish with respect to x.
	 * @param x Input value
	 * @return d/dx swish(x)
	 */
	@Override
	public double derivativeX(double x){return Functions.swishDerivativeX(x);}

	/**
	 * Rejects output-only differentiation; this overload is unimplemented.
	 * @param fx Pre-computed Swish function output
	 * @return Never returns normally
	 * @throws RuntimeException Always thrown; use an input-aware derivative method
	 */
	@Override
	public double derivativeFX(double fx){return Functions.swishDerivativeFX(fx);}

	/**
	 * Calculates the derivative using both input and precomputed output values.
	 * Convenient when both x and swish(x) are available.
	 *
	 * @param x Original input value
	 * @param fx Pre-computed Swish function output
	 * @return Derivative of Swish function
	 */
	@Override
	public double derivativeXFX(double x, double fx){return Functions.swishDerivativeXFX(x, fx);}

	/** Returns the numeric type identifier for this activation function.
	 * @return Numeric type constant for Swish function */
	@Override
	public int type(){return type;}

	/** Returns the string name identifier for this activation function.
	 * @return String name "SWISH" */
	@Override
	public String name(){return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Serialized activation name. */
	static final String name="SWISH";
	/** Registered index for the activation name. */
	static final int type=Function.toType(name, true);
	/** Shared stateless activation instance. */
	static final Swish instance=new Swish();

}
