package test.synth;

import synth.CrossContaminate;

/** Validates CrossContaminate probability-bound rejection before I/O opens. */
public class CrossContaminateValidationTest {

	public static void main(String[] args){
		int checks=0;
		checks+=rejectsZeroAndNonfinite("minprob=0", "maxprob=0.025");
		checks+=rejectsZeroAndNonfinite("minprob=0.001", "maxprob=0");
		checks+=rejectsZeroAndNonfinite("minprob=NaN", "maxprob=0.025");
		checks+=rejectsZeroAndNonfinite("minprob=0.5", "maxprob=0.25");
		checks+=rejectsZeroAndNonfinite("minprob=0.25", "maxprob=1.5");
		checks+=acceptsPositiveBoundsThenReachesInputValidation();
		System.out.println("CrossContaminateValidationTest PASS checks="+checks);
	}

	private static int rejectsZeroAndNonfinite(String minArg, String maxArg){
		try{
			make(minArg, maxArg);
		}catch(IllegalArgumentException e){
			if(!e.getMessage().contains("positive probabilities")){throw e;}
			return 1;
		}
		throw new AssertionError("Expected probability rejection for "+minArg+", "+maxArg);
	}

	private static int acceptsPositiveBoundsThenReachesInputValidation(){
		try{
			make("minprob=0.001", "maxprob=0.025");
		}catch(RuntimeException e){
			if(e instanceof IllegalArgumentException){throw e;}
			if(!e.getMessage().contains("Can't read file")){throw e;}
			return 1;
		}
		throw new AssertionError("Expected missing-input rejection after probability validation");
	}

	private static void make(String minArg, String maxArg){
		new CrossContaminate(new String[] {
				"in=__bbtools_missing_cross_1.fq,__bbtools_missing_cross_2.fq",
				"out=__bbtools_missing_cross_1.out,__bbtools_missing_cross_2.out",
				"minsinks=1", "maxsinks=1", minArg, maxArg});
	}

}
