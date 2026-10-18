package ml;

/**
 * Mirrored sigmoid activation function for neural networks.
 * Provides a symmetric sigmoid activation function centered around zero.
 * Produces finite-input outputs in [0,1], with a normalized peak at zero.
 * The derivative API returns negative zero at the nondifferentiable peak.
 *
 * @author Brian Bushnell
 * @date December 12, 2024
 */
public class MSig extends Function{
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Private constructor for singleton pattern. */
	private MSig(){}

	/** Evaluates the normalized mirrored sigmoid at x. */
	@Override
	public double activate(double x){return Functions.mSig(x);}

	/** Recomputes the activation and evaluates the input-aware derivative. */
	@Override
	public double derivativeX(double x){return Functions.mSigDerivativeX(x);}

	/** Rejects output-only differentiation because the output does not determine the input sign.
	 * @param fx Activation output
	 * @return Never returns normally
	 * @throws RuntimeException Always thrown */
	@Override
	public double derivativeFX(double fx){throw new RuntimeException("Cannot be calculated.");}

	/** Uses x for the derivative sign and fx for the normalized derivative magnitude. */
	@Override
	public double derivativeXFX(double x, double fx){return Functions.mSigDerivativeXFX(x, fx);}

	/** Returns the registered MSIG type index. */
	@Override
	public int type(){return type;}

	/** Returns the serialized activation name MSIG. */
	@Override
	public String name(){return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** String identifier for the mirrored sigmoid activation function. */
	static final String name="MSIG";
	/** Numeric type identifier derived from the function name. */
	static final int type=Function.toType(name, true);
	/** Singleton instance of the mirrored sigmoid activation function. */
	static final MSig instance=new MSig();

}
