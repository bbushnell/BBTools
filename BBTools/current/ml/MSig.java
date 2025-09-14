package ml;

public class MSig extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	static final String name="MSIG";
	static final int type=Function.toType(name, true);
	static final MSig instance=new MSig();

}
