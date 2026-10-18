package test.ml;

import ml.Matrix;
import ml.Sample;
import ml.SampleSet;
import ml.Seed;
import ml.Trainer;

import java.util.ArrayList;
import java.util.HashSet;

/**
 * Focused residual checks for the ML025, ML048, and ML049 deferred fixes.
 * @author Brian Bushnell, Sucrose
 */
public final class ResidualDeferralTest{

	/** Runs one residual-fix case. */
	public static void main(String[] args) throws Exception{
		boolean enabled=false; assert(enabled=true);
		if(!enabled || args.length!=2){
			throw new IllegalArgumentException("Run with -ea, baseline|fixed, and one case name");
		}
		final boolean fixed=parsePhase(args[0]);
		if(args[1].equals("fully-connected-edges")){
			testFullyConnectedEdges(fixed);
		}else if(args[1].equals("seed-equality")){
			testSeedEquality(fixed);
		}else if(args[1].equals("roc-overrun")){
			testRocOverrunAssertion(fixed);
		}else{
			throw new IllegalArgumentException("Unknown case: "+args[1]);
		}
		System.out.println("ResidualDeferralTest PASS "+args[0]+" "+args[1]);
	}

	/** Parses the expected source phase for baseline/fixed paired checks. */
	private static boolean parsePhase(String phase){
		if(phase.equals("baseline")){return false;}
		if(phase.equals("fixed")){return true;}
		throw new IllegalArgumentException("Unknown phase: "+phase);
	}

	/** Confirms fullyConnectedEdges widens operands before multiplication. */
	private static void testFullyConnectedEdges(final boolean fixed){
		final long expected=fixed ? 4295098368L : 131072L;
		final long observed=Trainer.fullyConnectedEdges(new int[]{65536, 65536, 2});
		assert(observed==expected) : "fullyConnectedEdges ML049 "+phaseName(fixed)+
			" expected "+expected+", observed "+observed;
		assert(Trainer.fullyConnectedEdges(new int[]{3, 4, 5})==32) :
			"fullyConnectedEdges unaffected small case should return 3*4+4*5";
	}

	/** Confirms seed-only equality and hash codes are null/type safe. */
	private static void testSeedEquality(final boolean fixed) throws Exception{
		final Seed low=new Seed(7, 0.1f);
		final Seed high=new Seed(7, 0.9f);
		final Seed other=new Seed(8, 0.1f);
		final Seed highSeedA=new Seed(1L<<40, 0.1f);
		final Seed highSeedB=new Seed(1L<<40, 0.9f);
		final Seed lowSeedSamePivot=new Seed(0, 0.1f);
		assert(low.equals(high)) : "Seed equality is seed-only (ML025)";
		assert(!low.equals(other)) : "Different low32 netSeed values must not compare equal";
		assert(highSeedA.equals(highSeedB)) : "Upper32 netSeed equality is seed-only";
		assert(!highSeedA.equals(lowSeedSamePivot)) :
			"Full long netSeed bits must participate in equality";
		assert(low.compareTo(high)<0) :
			"Seed compareTo should still order equal seeds by pivot";
		if(fixed){
			final Object nullObject=null;
			final Seed nullSeed=null;
			assert(!low.equals(nullSeed)) : "Seed.equals(Seed) must be null-safe";
			assert(!low.equals(nullObject)) : "Seed.equals(Object) must be null-safe";
			assert(!low.equals("7")) : "Seed.equals(Object) must be type-safe";
			assert(highSeedA.hashCode()==highSeedB.hashCode()) :
				"Seed hashCode must match upper32 seed-only equality (ML025)";
			final HashSet<Seed> set=new HashSet<Seed>();
			set.add(highSeedA);
			set.add(highSeedB);
			assert(set.size()==1) : "HashSet should deduplicate seed-only equals";
		}else{
			assertThrows(NullPointerException.class, () -> low.equals((Seed)null),
				"baseline Seed.equals(Seed) should throw for null");
			assertThrows(NullPointerException.class, () -> low.equals((Object)null),
				"baseline Seed.equals(Object) should throw for null");
			assertThrows(ClassCastException.class, () -> low.equals("7"),
				"baseline Seed.equals(Object) should throw for non-Seed");
			assert(Seed.class.getMethod("hashCode").getDeclaringClass()==Object.class) :
				"baseline Seed should inherit Object.hashCode";
		}
	}

	/** Confirms repeated rounded ROC targets now trip a diagnostic assertion. */
	private static void testRocOverrunAssertion(final boolean fixed){
		final SampleSet set=makeRocOverrunSet();
		if(fixed){
			final AssertionError expected=assertThrows(AssertionError.class,
				() -> set.calcROC(3), "fixed calcROC should assert before overrun");
			assert(expected.getMessage().contains("ML048")) : expected.getMessage();
		}else{
			assertThrows(ArrayIndexOutOfBoundsException.class, () -> set.calcROC(3),
				"baseline calcROC should index past repeated rounded ROC targets");
		}
		final float[] roc=makeRocControlSet().calcROC(2);
		assert(roc.length==2 && roc[0]==0.5f && roc[1]==1f) :
			"Unrepaired calcROC two-bin control changed: "+roc[0]+", "+roc[1];
	}

	/** Builds two negative rows before one positive row in sorted result order. */
	private static SampleSet makeRocOverrunSet(){
		resetMatrixStatics();
		final Matrix m=new Matrix();
		m.columns=new ArrayList<String>();
		m.columns.add("input");
		m.columns.add("target");
		m.numInputs=1;
		m.numOutputs=1;
		m.inputs=new float[][]{{0}, {1}, {2}};
		m.outputs=new float[][]{{0}, {0}, {1}};
		m.weights=new float[][]{{1}, {1}, {1}};
		m.data=new float[][][]{m.inputs, m.outputs, m.weights};
		m.validLines=m.inputs.length;
		m.initializeRange();
		final SampleSet set=new SampleSet(m);
		set.makeSamples();
		for(Sample s : set.samples){
			s.result[0]=s.goal[0];
		}
		set.sortByValue();
		return set;
	}

	/** Builds an unaffected two-bin ROC case with no repeated rounded target. */
	private static SampleSet makeRocControlSet(){
		final SampleSet set=makeSampleSet(new float[]{0, 1.1f, 0.2f, 1.3f, 0.4f},
			new float[]{0, 1, 2, 3, 4});
		set.sortByValue();
		return set;
	}

	/** Builds a sample set with caller-supplied goals and result ordering. */
	private static SampleSet makeSampleSet(final float[] goals, final float[] results){
		resetMatrixStatics();
		final Matrix m=new Matrix();
		m.columns=new ArrayList<String>();
		m.columns.add("input");
		m.columns.add("target");
		m.numInputs=1;
		m.numOutputs=1;
		m.inputs=new float[goals.length][1];
		m.outputs=new float[goals.length][1];
		m.weights=new float[goals.length][1];
		for(int i=0; i<goals.length; i++){
			m.inputs[i][0]=i;
			m.outputs[i][0]=goals[i];
			m.weights[i][0]=1;
		}
		m.data=new float[][][]{m.inputs, m.outputs, m.weights};
		m.validLines=m.inputs.length;
		m.initializeRange();
		final SampleSet set=new SampleSet(m);
		set.makeSamples();
		for(int i=0; i<set.samples.length; i++){
			final Sample s=set.samples[i];
			final float goal=s.goal[0];
			final int row=Math.round(goal>=1 ? (goal-1)*10 : goal*10);
			s.result[0]=results[row];
		}
		return set;
	}

	/** Asserts that a runnable throws exactly the expected class. */
	private static <T extends Throwable> T assertThrows(final Class<T> clazz,
			final ThrowingRunnable runnable, final String message){
		try{
			runnable.run();
		}catch(Throwable t){
			if(t.getClass()==clazz){return clazz.cast(t);}
			throw new AssertionError(message+"; expected "+clazz.getName()+
				", observed "+t.getClass().getName(), t);
		}
		throw new AssertionError(message+"; expected "+clazz.getName());
	}

	/** Returns the short name for the expected source phase. */
	private static String phaseName(final boolean fixed){return fixed ? "fixed" : "baseline";}

	/** Resets Matrix output remapping so this fixture uses literal 0/1 targets. */
	private static void resetMatrixStatics(){
		Matrix.convertTo01=false;
		Matrix.targetOutputRangeMin=0;
		Matrix.targetOutputRangeMax=0;
		Matrix.setTargetOutputRangeMin=false;
		Matrix.setTargetOutputRangeMax=false;
	}

	private ResidualDeferralTest(){}

	/** Runnable that can throw checked or unchecked exceptions. */
	private interface ThrowingRunnable{
		void run() throws Throwable;
	}
}
