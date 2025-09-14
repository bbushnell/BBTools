package ml;

public class Tanh extends Function {
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	private Tanh() {}

	@Override
	public double activate(double x) {return Functions.tanh(x);}

	@Override
	public double derivativeX(double x) {return Functions.tanhDerivativeX(x);}

	@Override
	public double derivativeFX(double fx) {return Functions.tanhDerivativeFX(fx);}

	@Override
	public double derivativeXFX(double x, double fx) {return Functions.tanhDerivativeXFX(x, fx);}

	@Override
	public int type() {return type;}

	@Override
	public String name() {return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	static final String name="TANH";
	static final int type=Function.toType(name, true);
	static final Tanh instance=new Tanh();

}
