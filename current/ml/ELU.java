package ml;

/** Exponential linear unit with alpha 1.
 * @author Brian Bushnell, Yoimiya */
public final class ELU extends Function{

	private ELU(){}
	@Override
	public double activate(final double x){return StandardFunctions.elu(x);}
	@Override
	public double derivativeX(final double x){return StandardFunctions.eluDerivative(x);}
	@Override
	public double derivativeXFX(final double x, final double fx){return StandardFunctions.eluDerivative(x);}
	@Override
	public double derivativeFX(final double fx){return fx>0 ? 1 : fx+1;}
	@Override
	public int type(){return Function.ELU;}
	@Override
	public String name(){return "ELU";}

	static final ELU instance=new ELU();
}
