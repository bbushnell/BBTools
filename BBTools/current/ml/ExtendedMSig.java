package ml;

public class ExtendedMSig extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	static final String name="EMSIG";
	static final int type=Function.toType(name, true);
	static final ExtendedMSig instance=new ExtendedMSig();

}
