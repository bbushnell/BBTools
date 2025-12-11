package ml;

/**
 * Bell curve (Gaussian) activation function f(x)=e^(-x^2) with outputs in [0,1] and peak at x=0.
 * Singleton implementation used in neural nets requiring smooth, normalized activations.
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class Bell extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Private constructor to enforce singleton usage via the static instance. */
	private Bell() {}

	/**
	 * Computes bell curve activation: e^(-x^2).
	 * @param x Input value
	 * @return Activation output in range (0,1]
	 */
	@Override
	public double activate(double x) {return Functions.bell(x);}

	/**
	 * Computes derivative with respect to x: -2x * e^(-x^2).
	 * @param x Input value
	 * @return Derivative of bell activation at x
	 */
	@Override
	public double derivativeX(double x) {return Functions.bellDerivativeX(x);}

	/**
	 * Computes derivative using function output.
	 * @param fx Output of bell activation
	 * @return Derivative corresponding to fx
	 */
	@Override
	public double derivativeFX(double fx) {return Functions.bellDerivativeFX(fx);}

	/**
	 * Computes derivative using both x and function output: -2x * fx.
	 * @param x Input value
	 * @param fx Output of bell activation
	 * @return Derivative of bell activation at x
	 */
	@Override
	public double derivativeXFX(double x, double fx) {return Functions.bellDerivativeXFX(x, fx);}

	/** Returns the numeric type identifier for this function.
	 * @return Type constant for Bell */
	@Override
	public int type() {return type;}

	/** Returns the string name identifier for this function.
	 * @return \"BELL\" */
	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static final String name="BELL";
	static final int type=Function.toType(name, true);
	static final Bell instance=new Bell();

}
