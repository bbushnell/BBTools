package test.ml;

import ml.Cell;
import ml.CellNet;
import ml.CellNetParser;
import ml.Function;
import ml.RegressionTrainer;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.Stream;

import shared.Timer;

/**
 * Finite same-JVM checks for RegressionTrainer activation-rate isolation.
 * Each invocation trains one small scalar model with valid input, checks its
 * explicit activation types and finite predictions, and restores Function globals.
 * Run one named case per assertion-enabled JVM; unrelated process globals are not restored.
 * @author Shinobu
 */
public final class RegressionTrainerStateTest{

	/** Runs fresh-default, already-normalized, or zero-rates-control. */
	public static void main(String[] args) throws Exception{
		boolean enabled=false; assert(enabled=true);
		if(!enabled || args.length!=1){throw new IllegalArgumentException("Run with -ea and one case name");}
		final String name=args[0];
		final float[] savedRates=Function.TYPE_RATES.clone();
		final float[] savedCumulative=Function.TYPE_RATES_CUM;
		final Path dir=Files.createTempDirectory("regression-trainer-state-");
		try{
			Arrays.fill(Function.TYPE_RATES, 0f);
			Function.TYPE_RATES_CUM=null;
			if(name.equals("fresh-default") || name.equals("already-normalized")){
				Function.TYPE_RATES[Function.TANH]=0.4f;
				Function.TYPE_RATES[Function.SIG]=0.6f;
				Function.TYPE_RATES[Function.MSIG]=0.02f;
				Function.TYPE_RATES[Function.RSLOG]=0.02f;
				if(name.equals("already-normalized")){Function.normalizeTypeRates();}
			}else if(!name.equals("zero-rates-control")){
				throw new IllegalArgumentException("Unknown case: "+name);
			}
			final float[] rates=Function.TYPE_RATES.clone();
			final float[] cumulative=Function.TYPE_RATES_CUM;
			final float[] cumulativeValues=cumulative==null ? null : cumulative.clone();
			final Path input=dir.resolve("rows.tsv"), output=dir.resolve("trained.bbnet");
			Files.write(input, ROWS.getBytes(StandardCharsets.US_ASCII));
			final RegressionTrainer trainer=new RegressionTrainer(new String[]{
				"in="+input, "out="+output, "dims=1,2,1", "epochs=1", "batch=2",
				"seed=3", "final=linear", "vfraction=0.25", "simd=f", "threads=1"
			});
			trainer.process(new Timer());
			assert(!trainer.errorState) : "Finite scalar training and round-trip checks must complete without an error latch";
			checkModel(output);
			assert(Arrays.equals(rates, Function.TYPE_RATES)) : name+": TYPE_RATES changed";
			assert(cumulative==Function.TYPE_RATES_CUM) : name+": TYPE_RATES_CUM reference changed";
			assert(Arrays.equals(cumulativeValues, Function.TYPE_RATES_CUM)) : name+": TYPE_RATES_CUM values changed";
			System.out.println("RegressionTrainerStateTest PASS "+name);
		}finally{
			System.arraycopy(savedRates, 0, Function.TYPE_RATES, 0, savedRates.length);
			Function.TYPE_RATES_CUM=savedCumulative;
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path->{
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	/** Checks explicit model shape/functions and emits predictions for cross-case parity. */
	private static void checkModel(Path path){
		final CellNet net=CellNetParser.load(path.toString());
		assert(net!=null && net.numInputs()==1 && net.numOutputs()==1) : "Training must write the requested one-input/one-output model";
		for(Cell cell : net.net[1]){
			assert(cell.function.type()==Function.TANH) : "Default hidden activation must remain tanh regardless of global rates";
		}
		assert(net.net[2][0].function.type()==Function.LINEAR) : "final=linear must determine the exported output activation";
		final float[] input=new float[1];
		for(int i=0; i<4; i++){
			input[0]=i*2;
			net.applyInput(input);
			final float prediction=net.feedForward();
			assert(Float.isFinite(prediction)) : "Finite fixture inputs must produce finite predictions: "+prediction;
			System.out.println("PRED\t"+i+"\t"+Float.toHexString(prediction));
		}
	}

	/** Eight exact binary-fraction targets; the last two shuffled rows form validation. */
	private static final String ROWS="#dims\t1\t1\n"
		+"0\t0\n1\t0.125\n2\t0.25\n3\t0.375\n"
		+"4\t0.5\n5\t0.625\n6\t0.75\n7\t0.875\n";

	private RegressionTrainerStateTest(){}
}
