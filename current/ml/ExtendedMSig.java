package ml;

/**
 * Extended mirrored sigmoid activation mapping inputs to [-1,1].
 * Even, bell-shaped function with a maximum of 1 at zero, approaching -1
 * as the input magnitude grows. The derivative API uses zero at the cusp.
 * Singleton implementation for neural network layers requiring centered outputs.
 * @author Brian Bushnell
 */
public class ExtendedMSig extends Function{
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates the shared stateless extended mirrored activation. */
	private ExtendedMSig(){}

	/**
	 * Applies the extended mirrored sigmoid activation to input x.
	 * @param x Input value
	 * @return emSig(x) in range [-1,1]
	 */
	@Override
	public double activate(double x){return Functions.emSig(x);}

	/**
	 * Computes the derivative of extended mirrored sigmoid with respect to x.
	 * @param x Input value
	 * @return Derivative at x
	 */
	@Override
	public double derivativeX(double x){return Functions.emSigDerivativeX(x);}

	/**
	 * Attempts to compute derivative from function output.
	 * Extended mirrored sigmoid derivative cannot be calculated from output alone.
	 * @param fx Function output value
	 * @return Never returns normally
	 * @throws RuntimeException Always thrown as calculation is impossible
	 */
	@Override
	public double derivativeFX(double fx){throw new RuntimeException("Cannot be calculated.");}

	/**
	 * Recomputes the mirrored derivative from x and doubles it; fx is ignored.
	 * @param x Input value
	 * @param fx Function output value (unused)
	 * @return Derivative at input x
	 */
	@Override
	public double derivativeXFX(double x, double fx){return Functions.emSigDerivativeXFX(x, fx);}

	/** Returns the numeric type identifier for this function.
	 * @return Function type constant */
	@Override
	public int type(){return type;}

	/** Returns the string name of this activation function.
	 * @return Function name "EMSIG" */
	@Override
	public String name(){return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Serialized activation name. */
	static final String name="EMSIG";
	/** Registered index for the activation name. */
	static final int type=Function.toType(name, true);
	/** Shared stateless activation instance. */
	static final ExtendedMSig instance=new ExtendedMSig();

}
