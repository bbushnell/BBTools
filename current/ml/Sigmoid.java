package ml;

public class Sigmoid extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private Sigmoid() {}

	/**
	 * Applies the sigmoid activation function to the input value.
	 * Computes f(x) = 1/(1+e^(-x)) which maps any real number to the range (0,1).
	 * @param x Input value to the activation function
	 * @return Sigmoid activation of x, in range (0,1)
	 */
	@Override
	public double activate(double x) {return Functions.sigmoid(x);}

	/**
	 * Computes the derivative of the sigmoid function with respect to input x.
	 * Uses the formula f'(x) = f(x) * (1 - f(x)) where f(x) is the sigmoid function.
	 * @param x Input value where derivative should be computed
	 * @return Derivative of sigmoid function at point x
	 */
	@Override
	public double derivativeX(double x) {return Functions.sigmoidDerivativeX(x);}

	/**
	 * Computes the derivative of the sigmoid function given the function output.
	 * More efficient than derivativeX when the sigmoid output is already known.
	 * Uses the formula f'(x) = fx * (1 - fx) where fx = f(x).
	 *
	 * @param fx Pre-computed sigmoid function output f(x)
	 * @return Derivative of sigmoid function
	 */
	@Override
	public double derivativeFX(double fx) {return Functions.sigmoidDerivativeFX(fx);}

	/**
	 * Computes the derivative of the sigmoid function given both input and output.
	 * Implementation uses the pre-computed function output for efficiency.
	 *
	 * @param x Input value (unused in this implementation)
	 * @param fx Pre-computed sigmoid function output f(x)
	 * @return Derivative of sigmoid function
	 */
	@Override
	public double derivativeXFX(double x, double fx) {return Functions.sigmoidDerivativeFX(fx);}

	/** Returns the integer type identifier for the sigmoid function.
	 * @return Function type constant from parent class */
	@Override
	public int type() {return type;}

	/** Returns the string name identifier for the sigmoid function.
	 * @return Function name "SIG" */
	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static final String name="SIG";
	static final int type=Function.toType(name, true);
	static final Sigmoid instance=new Sigmoid();

}
