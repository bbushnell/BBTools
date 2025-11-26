package ml;

/**
 * Bell curve (Gaussian) activation function for neural networks.
 * Implements the mathematical function f(x) = e^(-x²) which produces
 * a symmetric bell-shaped curve with maximum value of 1 at x=0.
 * Output values range from 0 to 1, making it suitable for normalized
 * activation in neural network architectures.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class Bell extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Private constructor prevents direct instantiation. Use Bell.instance instead.
	 */
	private Bell() {}

	@Override
	public double activate(double x) {return Functions.bell(x);}

	@Override
	public double derivativeX(double x) {return Functions.bellDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.bellDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.bellDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** String identifier for the Bell activation function. */
	static final String name="BELL";
	/** Numeric type identifier derived from the string name. */
	static final int type=Function.toType(name, true);
	/** Singleton instance of the Bell activation function. */
	static final Bell instance=new Bell();

}
