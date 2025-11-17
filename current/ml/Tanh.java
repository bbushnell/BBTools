package ml;

/**
 * Hyperbolic tangent activation function for neural networks.
 * Provides symmetric output range from -1 to 1, making it suitable for
 * centered data and neural network layers requiring zero-centered outputs.
 * Delegates mathematical calculations to the Functions utility class.
 *
 * @author Brian Bushnell
 */
public class Tanh extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private Tanh() {}

	@Override
	public double activate(double x) {return Functions.tanh(x);}

	@Override
	public double derivativeX(double x) {return Functions.tanhDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.tanhDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.tanhDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** String identifier for the tanh activation function */
	static final String name="TANH";
	/** Numeric type identifier for the tanh activation function */
	static final int type=Function.toType(name, true);
	/** Singleton instance of the Tanh activation function */
	static final Tanh instance=new Tanh();

}
