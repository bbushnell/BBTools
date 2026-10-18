package ml;

/**
 * Rotationally symmetric logarithmic activation function.
 * Implements rslog(x) = log(|x|+1) * sign(x), providing a smooth logarithmic
 * transformation that preserves sign and handles both positive and negative inputs.
 * Part of the machine learning function library for neural network activations.
 *
 * @author Brian Bushnell
 */
public class RSLog extends Function{
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates the shared stateless activation instance. */
	private RSLog(){}

	/** Evaluates the sign-symmetric logarithm through Functions.rslog. */
	@Override
	public double activate(double x){return Functions.rslog(x);}

	/** Evaluates 1/(abs(x)+1) from the input. */
	@Override
	public double derivativeX(double x){return Functions.rslogDerivativeX(x);}

	/** Legacy output-only path: asserts under {@code -ea}, otherwise uses inverse-log arithmetic.
	 * @param fx Precomputed activation value
	 * @return Derivative reconstructed from the activation when assertions are disabled
	 * @throws AssertionError Always thrown when assertions are enabled */
	@Override
	public double derivativeFX(double fx){return Functions.rslogDerivativeFX(fx);} //Possible bug: Functions.rslogDerivativeFX contains assert(false)

	/** Uses the input-aware derivative; the supplied activation is ignored. */
	@Override
	public double derivativeXFX(double x, double fx){return Functions.rslogDerivativeXFX(x, fx);}

	/** Returns the registered RSLOG type index. */
	@Override
	public int type(){return type;}

	/** Returns the serialized activation name RSLOG. */
	@Override
	public String name(){return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Function name identifier for this activation function */
	static final String name="RSLOG";
	/** Function type constant derived from the name for lookup operations */
	static final int type=Function.toType(name, true);
	/** Singleton instance of the RSLog activation function */
	static final RSLog instance=new RSLog();

}
