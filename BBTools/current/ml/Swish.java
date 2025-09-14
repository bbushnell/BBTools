package ml;

public class Swish extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private Swish() {}

	@Override
	public double activate(double x) {return Functions.swish(x);}

	@Override
	public double derivativeX(double x) {return Functions.swishDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.swishDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.swishDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static final String name="SWISH";
	static final int type=Function.toType(name, true);
	static final Swish instance=new Swish();

}
