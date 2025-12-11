package ml;

public class Tanh extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private Tanh() {}

	/**
	 * Computes the hyperbolic tangent activation of the input value.
	 * @param x Input value to transform
	 * @return Hyperbolic tangent of x, ranging from -1 to 1
	 */
	@Override
	public double activate(double x) {return Functions.tanh(x);}

	/**
	 * Computes the derivative of tanh with respect to the input value x.
	 * @param x Input value at which to compute the derivative
	 * @return Derivative of tanh(x) with respect to x
	 */
	@Override
	public double derivativeX(double x) {return Functions.tanhDerivativeX(x);}

	/**
	 * Computes the derivative of tanh using the already-computed function value.
	 * More efficient than derivativeX when the function value is already known.
	 * @param fx Pre-computed tanh(x) value
	 * @return Derivative of tanh at the point where tanh(x) = fx
	 */
	@Override
	public double derivativeFX(double fx) {return Functions.tanhDerivativeFX(fx);}

	/**
	 * Computes the derivative of tanh using both input and function values.
	 * Provides flexibility for optimization based on available computed values.
	 *
	 * @param x Input value
	 * @param fx Pre-computed tanh(x) value
	 * @return Derivative of tanh at x
	 */
	@Override
	public double derivativeXFX(double x, double fx) {return Functions.tanhDerivativeXFX(x, fx);}

	/** Returns the numeric type identifier for the tanh function.
	 * @return Type constant for tanh activation function */
	@Override
	public int type() {return type;}

	/** Returns the string name identifier for the tanh function.
	 * @return "TANH" string identifier */
	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static final String name="TANH";
	static final int type=Function.toType(name, true);
	static final Tanh instance=new Tanh();

}
