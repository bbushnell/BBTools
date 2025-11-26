package ml;

/**
 * Extended sigmoid activation function for neural networks with symmetric output range [-1, 1].
 * Implements the transformation eSigmoid(x) = 2*sigmoid(x) - 1, shifting standard sigmoid
 * from [0,1] range to [-1,1] range for improved neural network training dynamics.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class ExtendedSigmoid extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Private constructor for singleton pattern implementation */
	private ExtendedSigmoid() {}

	@Override
	public double activate(double x) {return Functions.eSigmoid(x);}

	@Override
	public double derivativeX(double x) {return Functions.eSigmoidDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.eSigmoidDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.eSigmoidDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** String identifier for extended sigmoid function type */
	static final String name="ESIG";
	/** Numeric type identifier derived from the function name. */
	static final int type=Function.toType(name, true);
	/** Singleton instance for extended sigmoid activation function */
	static final ExtendedSigmoid instance=new ExtendedSigmoid();

}
