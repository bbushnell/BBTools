package ml;

/**
 * Identity activation function: f(x)=x, f'(x)=1.
 * Useful as a final-layer activation for regression networks whose target is
 * unbounded (or already in a convenient space), where a squashing or
 * compressive function would distort the output.  Zero cost, exact gradients.
 *
 * @author Amber
 * @date July 2026
 */
public class Linear extends Function {

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private Linear() {}

	/**
	 * Identity: returns the input unchanged.
	 * @param x Input value
	 * @return x
	 */
	@Override
	public double activate(double x) {return x;}

	/**
	 * Derivative of the identity function is 1 everywhere.
	 * @param x Input value (unused)
	 * @return 1
	 */
	@Override
	public double derivativeX(double x) {return 1;}

	/**
	 * Derivative of the identity function is 1 everywhere.
	 * @param fx Pre-computed function value (unused)
	 * @return 1
	 */
	@Override
	public double derivativeFX(double fx) {return 1;}

	/**
	 * Derivative of the identity function is 1 everywhere.
	 * @param x Input value (unused)
	 * @param fx Pre-computed function value (unused)
	 * @return 1
	 */
	@Override
	public double derivativeXFX(double x, double fx) {return 1;}

	/** Returns the numeric type identifier for the linear function.
	 * @return Type constant for linear activation function */
	@Override
	public int type() {return type;}

	/** Returns the string name identifier for the linear function.
	 * @return "LINEAR" string identifier */
	@Override
	public String name() {return name;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	static final String name="LINEAR";
	static final int type=Function.toType(name, true);
	static final Linear instance=new Linear();

}
