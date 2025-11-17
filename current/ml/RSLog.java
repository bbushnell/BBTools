package ml;

/**
 * Rotationally symmetric logarithmic activation function.
 * Implements rslog(x) = log(|x|+1) * sign(x), providing a smooth logarithmic
 * transformation that preserves sign and handles both positive and negative inputs.
 * Part of the machine learning function library for neural network activations.
 *
 * @author Brian Bushnell
 */
public class RSLog extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private RSLog() {}

	@Override
	public double activate(double x) {return Functions.rslog(x);}

	@Override
	public double derivativeX(double x) {return Functions.rslogDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.rslogDerivativeFX(fx);} //Possible bug: Functions.rslogDerivativeFX contains assert(false)

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.rslogDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
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
