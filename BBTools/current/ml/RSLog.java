package ml;

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
	public double derivativeFX(double fx) {return Functions.rslogDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.rslogDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static final String name="RSLOG";
	static final int type=Function.toType(name, true);
	static final RSLog instance=new RSLog();

}
