package ml;

/**
 * Mirrored sigmoid activation function for neural networks.
 * Provides a symmetric sigmoid activation function centered around zero.
 * The function applies sigmoid to both positive and negative inputs with mirroring,
 * producing outputs in the range [0, 1] with symmetric behavior around x=0.
 *
 * @author Brian Bushnell
 * @date December 12, 2024
 */
public class MSig extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Private constructor for singleton pattern. */
	private MSig() {}

	@Override
	public double activate(double x) {return Functions.mSig(x);}

	@Override
	public double derivativeX(double x) {return Functions.mSigDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {throw new RuntimeException("Cannot be calculated.");}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.mSigDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
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
