package ml;

/**
 * Extended mirrored sigmoid activation function for neural networks.
 * Implements an extended version of mirrored sigmoid that maps input to range [-1,1].
 * The function is symmetric around the origin and has a bell-like shape with
 * maximum values at approximately ±2.5 and minimum at 0.
 *
 * @author Brian Bushnell
 */
public class ExtendedMSig extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Private constructor prevents external instantiation. */
	private ExtendedMSig() {}

	@Override
	public double activate(double x) {return Functions.emSig(x);}

	@Override
	public double derivativeX(double x) {return Functions.emSigDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {throw new RuntimeException("Cannot be calculated.");}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.emSigDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** String identifier for extended mirrored sigmoid function. */
	static final String name="EMSIG";
	/** Numeric type identifier computed from function name. */
	static final int type=Function.toType(name, true);
	/** Singleton instance of the extended mirrored sigmoid function. */
	static final ExtendedMSig instance=new ExtendedMSig();

}
