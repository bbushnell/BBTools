package ml;

/** Logistic activation and derivatives, obtained through {@link Function#getFunction(int)}. */
public class Sigmoid extends Function{
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates the shared stateless activation instance. */
	private Sigmoid(){}

	/**
	 * Applies the sigmoid activation function to the input value.
	 * Computes f(x) = 1/(1+e^(-x)); floating-point evaluation can saturate at 0 or 1.
	 * @param x Input value to the activation function
	 * @return Sigmoid activation of finite x, in range [0,1]
	 */
	@Override
	public double activate(double x){return Functions.sigmoid(x);}

	/**
	 * Computes the derivative of the sigmoid function with respect to input x.
	 * Uses the formula f'(x) = f(x) * (1 - f(x)) where f(x) is the sigmoid function.
	 * @param x Input value where derivative should be computed
	 * @return Derivative of sigmoid function at point x
	 */
	@Override
	public double derivativeX(double x){return Functions.sigmoidDerivativeX(x);}

	/**
	 * Computes the derivative of the sigmoid function given the function output.
	 * Uses the supplied output without recomputing the activation.
	 * Uses the formula f'(x) = fx * (1 - fx) where fx = f(x).
	 *
	 * @param fx Pre-computed sigmoid function output f(x)
	 * @return Derivative of sigmoid function
	 */
	@Override
	public double derivativeFX(double fx){return Functions.sigmoidDerivativeFX(fx);}

	/**
	 * Computes the derivative of the sigmoid function given both input and output.
	 * Uses the supplied output; the input argument is ignored.
	 *
	 * @param x Input value (unused in this implementation)
	 * @param fx Pre-computed sigmoid function output f(x)
	 * @return Derivative of sigmoid function
	 */
	@Override
	public double derivativeXFX(double x, double fx){return Functions.sigmoidDerivativeFX(fx);}

	/** Returns the integer type identifier for the sigmoid function.
	 * @return Function type constant from parent class */
	@Override
	public int type(){return type;}

	/** Returns the string name identifier for the sigmoid function.
	 * @return Function name "SIG" */
	@Override
	public String name(){return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Serialized activation name. */
	static final String name="SIG";
	/** Registered index for the activation name. */
	static final int type=Function.toType(name, true);
	/** Shared stateless activation instance. */
	static final Sigmoid instance=new Sigmoid();

}
