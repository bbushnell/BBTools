package ml;

/**
 * Bell curve (Gaussian) activation f(x)=e^(-x^2), with finite-input outputs in [0,1]
 * and a peak at x=0. NaN inputs propagate through the underlying arithmetic.
 * Obtain the singleton through {@link Function#getFunction(int)}.
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class Bell extends Function{
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Private constructor to enforce singleton usage via the static instance. */
	private Bell(){}

	/**
	 * Computes bell curve activation: e^(-x^2).
	 * @param x Input value
	 * @return Activation output in [0,1], including zero after floating-point underflow
	 */
	@Override
	public double activate(double x){return Functions.bell(x);}

	/**
	 * Computes derivative with respect to x: -2x * e^(-x^2).
	 * @param x Input value
	 * @return Derivative of bell activation at x
	 */
	@Override
	public double derivativeX(double x){return Functions.bellDerivativeX(x);}

	/**
	 * Rejects output-only differentiation because the output does not determine x's sign.
	 * @param fx Output of bell activation
	 * @return Never returns normally
	 * @throws RuntimeException Always thrown; use an input-aware derivative method
	 */
	@Override
	public double derivativeFX(double fx){return Functions.bellDerivativeFX(fx);}

	/**
	 * Computes derivative using both x and function output: -2x * fx.
	 * @param x Input value
	 * @param fx Output of bell activation
	 * @return Derivative of bell activation at x
	 */
	@Override
	public double derivativeXFX(double x, double fx){return Functions.bellDerivativeXFX(x, fx);}

	/** Returns the numeric type identifier for this function.
	 * @return Type constant for Bell */
	@Override
	public int type(){return type;}

	/** Returns the string name identifier for this function.
	 * @return "BELL" */
	@Override
	public String name(){return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Serialized activation name. */
	static final String name="BELL";
	/** Registered index for the activation name. */
	static final int type=Function.toType(name, true);
	/** Shared stateless activation instance. */
	static final Bell instance=new Bell();

}
