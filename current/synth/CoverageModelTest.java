package synth;

import shared.Random;

/** Finite analytic oracles for the single-wave coverage failure and nearby controls.
 * @author Brian Bushnell, Jean
 */
public final class CoverageModelTest{

	/** Checks twelve model/decision cases without entering the generator's retry loop.
	 * @param args One phase=baseline or phase=fixed argument
	 */
	public static void main(String[] args){
		require(args.length==1 && (args[0].equals("phase=baseline") || args[0].equals("phase=fixed")),
				"Specify phase=baseline or phase=fixed for the expected one-wave behavior");
		require(CoverageModel.class.desiredAssertionStatus(), "Production assertions must be enabled");
		final boolean baseline=args[0].equals("phase=baseline");
		int cases=0, mismatches=0;
		System.out.println("waves\tcontig_length\tposition\texpected_probability\tactual_probability\texpected_keep\tactual_keep\tconstructor_floats\tconstructor_doubles\tmatch");
		for(int waves : new int[]{1, 2, 4}){
			for(int length : new int[]{1999, 2000}){
				for(int position : new int[]{0, 500}){
					final FixedRandom setup=new FixedRandom(0);
					final CoverageModel model=new CoverageModel(waves, 1, 0, 0.25, setup, 2000, 2000);
					// Short contigs bypass the model. Otherwise sin(0)=0 and sin(pi/2)=1:
					// at zero the eighth-power valley is1/256, yielding1-3/1024.
					final float expected=length<2000 ? 1f : position==0 ? 1021f/1024f : 0.25f;
					final float actual=model.getProbability(position, length);
					final FixedRandom decision=new FixedRandom(0.5);
					final boolean expectedKeep=0.5<=expected;
					final boolean actualKeep=model.shouldGenerateReadAt(position, length, decision);
					require(setup.floatCalls==waves+1 && setup.doubleCalls==waves,
							"Constructor must consume one origin draw and one period/phase pair per wave");
					require(decision.doubleCalls==1 && decision.floatCalls==0,
							"The read decision must consume exactly one double draw");
					final boolean match=Float.isFinite(actual) && Math.abs(actual-expected)<=0.000001f
							&& actualKeep==expectedKeep;
					final boolean knownFailure=baseline && waves==1 && length==2000;
					if(knownFailure){
						require(Float.isNaN(actual) && !actualKeep, "Baseline one-wave case must expose NaN/rejection");
					}else{require(match, "Analytic probability/decision mismatch for waves="+waves+", length="+length+", position="+position);}
					if(!match){mismatches++;}
					cases++;
					System.out.println(waves+"\t"+length+"\t"+position+"\t"+expected+"\t"+actual
							+"\t"+expectedKeep+"\t"+actualKeep+"\t"+setup.floatCalls+"\t"+setup.doubleCalls+"\t"+match);
				}
			}
		}
		require(cases==12 && mismatches==(baseline ? 2 : 0), "Unexpected case or mismatch count");
	}

	/** Supplies zero period draws and a chosen phase/decision draw; counts every use. */
	private static final class FixedRandom implements Random{
		FixedRandom(double value){doubleValue=value;}
		@Override
		public float nextFloat(){floatCalls++; return 0;}
		@Override
		public double nextDouble(){doubleCalls++; return doubleValue;}
		@Override
		public long nextLong(){throw new AssertionError("Unexpected long draw in the coverage oracle");}
		@Override
		public long nextLong(long bound){throw new AssertionError("Unexpected bounded long draw in the coverage oracle");}
		@Override
		public double nextGaussian(){throw new AssertionError("Unexpected Gaussian draw in the coverage oracle");}
		@Override
		public void setSeed(long seed){throw new AssertionError("Unexpected reseeding in the coverage oracle");}
		final double doubleValue;
		int floatCalls, doubleCalls;
	}

	/** Keeps fixture checks active independently of assertion mode. */
	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Static-only test entry point. */
	private CoverageModelTest(){}
}
