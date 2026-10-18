package test.ml;

import ml.BBNetExporter;
import ml.Cell;
import ml.CellNet;
import ml.Function;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Same-JVM state and literal inference fixture for {@link BBNetExporter}.
 * Each selected case restores Function's process-global activation rates after checking export
 * behavior or ordinary {@link CellNet#randomize()} activation selection.
 * @author Jean
 */
public final class BBNetExporterTest{

	/** Runs one assertion-backed case selected by name. */
	public static void main(String[] args){
		boolean enabled=false; assert(enabled=true);
		if(!enabled || args.length!=1){throw new IllegalArgumentException("Run this fixture with -ea and one case name");}
		final float[] savedRates=Function.TYPE_RATES.clone();
		final float[] savedCumulative=Function.TYPE_RATES_CUM==null ? null : Function.TYPE_RATES_CUM.clone();
		try{
			if(args[0].equals("fresh-default")){
				resetDefaultRates();
				checkExport(args[0]);
			}else if(args[0].equals("already-normalized")){
				resetDefaultRates();
				Function.normalizeTypeRates();
				checkExport(args[0]);
			}else if(args[0].equals("randomize-normalized")){
				checkRandomize(args[0]);
			}else{
				throw new IllegalArgumentException("Unknown case: "+args[0]);
			}
		}finally{
			System.arraycopy(savedRates, 0, Function.TYPE_RATES, 0, savedRates.length);
			Function.TYPE_RATES_CUM=savedCumulative;
		}
		System.out.println("BBNetExporterTest PASS "+args[0]);
	}

	/** Exercises export while requiring unchanged activation-selection globals. */
	private static void checkExport(String label){
		final float[] rates=Function.TYPE_RATES.clone();
		final float[] cumulative=Function.TYPE_RATES_CUM==null ? null : Function.TYPE_RATES_CUM.clone();
		final double[][] w1={{0.5, -1.5, 2.0}, {-0.25, 0.75, 1.25}};
		final double[] b1={0.125, -0.5};
		final double[] w2={1.5, -0.25};
		final double b2=0.75;
		final double[] mean={1.0, -2.0, 0.5};
		final double[] sd={2.0, 4.0, 0.25};
		final CellNet net=BBNetExporter.export(w1, b1, w2, b2, mean, sd, Function.LINEAR);
		assertUnchanged(label, rates, cumulative);
		final float[] raw={3.0f, 6.0f, 0.75f};
		final double h0=Math.tanh(b1[0]+w1[0][0]*(raw[0]-mean[0])/sd[0]
				+w1[0][1]*(raw[1]-mean[1])/sd[1]+w1[0][2]*(raw[2]-mean[2])/sd[2]);
		final double h1=Math.tanh(b1[1]+w1[1][0]*(raw[0]-mean[0])/sd[0]
				+w1[1][1]*(raw[1]-mean[1])/sd[1]+w1[1][2]*(raw[2]-mean[2])/sd[2]);
		final float expected=(float)(b2+w2[0]*h0+w2[1]*h1);
		net.applyInput(raw);
		final float actual=net.feedForward();
		assert(Math.abs(actual-expected)<1e-6f) : label+": literal oracle mismatch "+actual+" != "+expected;
		assert(net.net[1][0].function.type()==Function.TANH && net.net[2][0].function.type()==Function.LINEAR) :
			label+": exporter must set deterministic activation functions after topology allocation";
	}

	/** Verifies that normal randomize() still allocates topology and selects hidden activations. */
	private static void checkRandomize(String label){
		Arrays.fill(Function.TYPE_RATES, 0);
		Function.TYPE_RATES[Function.LINEAR]=1;
		Function.TYPE_RATES_CUM=null;
		Function.normalizeTypeRates();
		final float[] rates=Function.TYPE_RATES.clone();
		final float[] cumulative=Function.TYPE_RATES_CUM.clone();
		final CellNet net=new CellNet(new int[]{2, 3, 1}, 5, 1f, 0f, 1, new ArrayList<String>());
		net.randomize();
		assertUnchanged(label, rates, cumulative);
		for(Cell c : net.net[1]){
			assert(c.function.type()==Function.LINEAR) : label+": randomize() did not select hidden activations";
			assert(c.weights!=null && c.weights.length==2) : label+": missing hidden dense weights";
		}
		final Cell out=net.net[2][0];
		assert(out.weights!=null && out.weights.length==3) : label+": missing output dense weights";
	}

	/** Restores Function's class-load defaults without normalizing the cumulative table. */
	private static void resetDefaultRates(){
		Arrays.fill(Function.TYPE_RATES, 0);
		Function.TYPE_RATES[Function.TANH]=0.4f;
		Function.TYPE_RATES[Function.SIG]=0.6f;
		Function.TYPE_RATES[Function.MSIG]=0.02f;
		Function.TYPE_RATES[Function.RSLOG]=0.02f;
		Function.TYPE_RATES_CUM=null;
	}

	/** Checks that export acted as a reusable helper, not a process-global rate mutator. */
	private static void assertUnchanged(String label, float[] rates, float[] cumulative){
		assert(Arrays.equals(Function.TYPE_RATES, rates)) : label+": TYPE_RATES changed";
		assert(Arrays.equals(Function.TYPE_RATES_CUM, cumulative)) : label+": TYPE_RATES_CUM changed";
	}

	/** Utility class only. */
	private BBNetExporterTest(){}
}
