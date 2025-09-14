package ml;

public class ExtendedSigmoid extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	static final String name="ESIG";
	static final int type=Function.toType(name, true);
	static final ExtendedSigmoid instance=new ExtendedSigmoid();

}
