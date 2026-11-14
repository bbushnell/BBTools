package test.synth;

import synth.CoverageModel;

import shared.Random;

/** Finite oracle for zero wave-amplitude coverage-model construction. */
public final class CoverageModelZeroTest{

	/** Verifies current zero-amplitude behavior and a nearby valid control.
	 * @param args One phase=baseline or phase=fixed argument
	 */
	public static void main(String[] args){
		require(args.length==1 && (args[0].equals("phase=baseline") || args[0].equals("phase=fixed")),
				"Specify phase=baseline or phase=fixed for the expected zero-amplitude behavior");
		final boolean baseline=args[0].equals("phase=baseline");
		int cases=0, mismatches=0;
		System.out.println("case\texpected\tactual\tconstructor_floats\tconstructor_doubles\tmatch");

		final FixedRandom zeroSetup=new FixedRandom(0);
		boolean zeroMatch=false;
		try {
			final CoverageModel model=new CoverageModel(1, 0, 0, 0.25, zeroSetup, 2000, 2000);
			final float actual=model.getProbability(0, 2000);
			zeroMatch=baseline && Float.isNaN(actual);
			require(zeroMatch, "Fixed phase must reject zero wave amplitude before NaN probability");
			System.out.println("zero_wave_amplitude\tIllegalArgumentException\t"+actual+"\t"
					+zeroSetup.floatCalls+"\t"+zeroSetup.doubleCalls+"\t"+zeroMatch);
		}catch(IllegalArgumentException e){
			zeroMatch=!baseline && e.getMessage()!=null && e.getMessage().indexOf("maxAmplitude")>=0
					&& zeroSetup.floatCalls==0 && zeroSetup.doubleCalls==0;
			require(zeroMatch, "Baseline should expose NaN; fixed should reject maxAmplitude before draws");
			System.out.println("zero_wave_amplitude\tIllegalArgumentException\t"+e.getClass().getSimpleName()
					+"\t"+zeroSetup.floatCalls+"\t"+zeroSetup.doubleCalls+"\t"+zeroMatch);
		}
		cases++;
		if(!zeroMatch){mismatches++;}

		final FixedRandom controlSetup=new FixedRandom(0);
		final CoverageModel control=new CoverageModel(1, 1, 0, 0.25, controlSetup, 2000, 2000);
		final float actual=control.getProbability(0, 2000);
		final float expected=1021f/1024f;
		final boolean controlMatch=Math.abs(actual-expected)<=0.000001f
				&& controlSetup.floatCalls==2 && controlSetup.doubleCalls==1;
		require(controlMatch, "Positive wave amplitude control mismatch");
		cases++;
		System.out.println("positive_wave_amplitude\t"+expected+"\t"+actual+"\t"
				+controlSetup.floatCalls+"\t"+controlSetup.doubleCalls+"\t"+controlMatch);

		require(cases==2 && mismatches==0, "Unexpected zero-amplitude case count or mismatch count");
	}

	/** Supplies zero period and phase draws; counts every use. */
	private static final class FixedRandom implements Random{
		FixedRandom(double value){doubleValue=value;}
		@Override
		public float nextFloat(){floatCalls++; return 0;}
		@Override
		public double nextDouble(){doubleCalls++; return doubleValue;}
		@Override
		public long nextLong(){throw new AssertionError("Unexpected long draw in the zero-amplitude oracle");}
		@Override
		public long nextLong(long bound){throw new AssertionError("Unexpected bounded long draw");}
		@Override
		public double nextGaussian(){throw new AssertionError("Unexpected Gaussian draw");}
		@Override
		public void setSeed(long seed){throw new AssertionError("Unexpected reseeding");}
		final double doubleValue;
		int floatCalls, doubleCalls;
	}

	/** Keeps fixture checks active independently of assertion mode. */
	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Static-only test entry point. */
	private CoverageModelZeroTest(){}
}
