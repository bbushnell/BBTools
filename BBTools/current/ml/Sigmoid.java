package ml;

public class Sigmoid extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private Sigmoid() {}

	@Override
	public double activate(double x) {return Functions.sigmoid(x);}

	@Override
	public double derivativeX(double x) {return Functions.sigmoidDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.sigmoidDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.sigmoidDerivativeFX(fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static final String name="SIG";
	static final int type=Function.toType(name, true);
	static final Sigmoid instance=new Sigmoid();

}
