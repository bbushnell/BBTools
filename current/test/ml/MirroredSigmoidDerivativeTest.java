package test.ml;

import ml.Function;

/** Independent analytic and finite-difference checks for mirrored sigmoid derivatives.
 * @author Brian Bushnell, Shinobu
 */
public final class MirroredSigmoidDerivativeTest{

	/** Runs one isolated case, or all cases when no selector is supplied. */
	public static void main(String[] args){
		if(args.length>1){throw new IllegalArgumentException("Expected at most one case name");}
		if(args.length==0){
			run("msig_nonzero");
			run("emsig_nonzero");
			run("peak_controls");
		}else{run(args[0]);}
	}

	/** Keeps the existing zero-cusp convention separate from differentiable points. */
	private static void run(String name){
		//Initialize the base registry before touching subclass singletons.
		if(name.equals("peak_controls")){
			for(Function function : new Function[]{Function.getFunction(Function.MSIG),
				Function.getFunction(Function.EMSIG)}){
				near(function.activate(0), 1, 1e-14, function.name()+" peak");
				near(function.activate(-1), function.activate(1), 1e-14, function.name()+" symmetry");
				near(function.derivativeX(0), 0, 0, function.name()+" existing zero-cusp convention");
				near(function.derivativeXFX(0, 1), 0, 0, function.name()+" cached zero-cusp convention");
			}
		}else if(name.equals("msig_nonzero")){
			checkNonzero(Function.getFunction(Function.MSIG), 1);
		}else if(name.equals("emsig_nonzero")){
			checkNonzero(Function.getFunction(Function.EMSIG), 2);
		}else{throw new IllegalArgumentException("Unknown case: "+name);}
		System.err.println("MirroredSigmoidDerivativeTest PASS "+name);
	}

	/** Checks both derivative APIs against an exponential oracle and a numerical slope. */
	private static void checkNonzero(Function function, int scale){
		final double[] points={-5, -2.5, -0.5, 0.5, 2.5, 5};
		final double h=1e-5;
		boolean correct=true;
		for(double x : points){
			//For x!=0, f(x)=M/(1+exp(2*abs(x)-5)), M=1+exp(-5).
			//Differentiate that expression directly; do not derive the oracle from fx.
			final double exponential=Math.exp(2*Math.abs(x)-5);
			final double normalizer=1+Math.exp(-5);
			final double denominator=1+exponential;
			final double expected=(x<0 ? 1 : -1)*scale*2*normalizer*exponential/
				(denominator*denominator);
			final double value=scale*normalizer/denominator-(scale-1);
			final double fx=function.activate(x);
			near(fx, value, 1e-13, function.name()+" activation at "+x);
			final double numerical=(function.activate(x+h)-function.activate(x-h))/(2*h);
			near(numerical, expected, 2e-8, function.name()+" independent numerical slope at "+x);
			final double fromX=function.derivativeX(x);
			final double fromXFX=function.derivativeXFX(x, fx);
			System.out.println(function.name()+"\t"+x+"\t"+expected+"\t"+numerical+
				"\t"+fromX+"\t"+fromXFX);
			correct&=Double.isFinite(fromX) && Math.abs(fromX-expected)<=1e-12;
			correct&=Double.isFinite(fromXFX) && Math.abs(fromXFX-expected)<=1e-12;
		}
		if(!correct){throw new AssertionError(function.name()+" derivative disagrees with independent oracle");}
	}

	/** Rejects non-finite values and discrepancies without depending on JVM assertion flags. */
	private static void near(double actual, double expected, double tolerance, String context){
		if(!Double.isFinite(actual) || Math.abs(actual-expected)>tolerance){
			throw new AssertionError(context+": actual="+actual+", expected="+expected);
		}
	}
}
