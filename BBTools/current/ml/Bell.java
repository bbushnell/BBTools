package ml;

public class Bell extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private Bell() {}

	@Override
	public double activate(double x) {return Functions.bell(x);}

	@Override
	public double derivativeX(double x) {return Functions.bellDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.bellDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.bellDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static final String name="BELL";
	static final int type=Function.toType(name, true);
	static final Bell instance=new Bell();

}
