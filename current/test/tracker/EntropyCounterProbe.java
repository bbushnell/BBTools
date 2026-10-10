package test.tracker;

import tracker.EntropyTracker;

/** Finite JT011 DNA counter observations; accepted all-A windows should stay homogeneous.
 * This does not decide which window sizes the API should support.
 * @author Brian Bushnell, Jean */
public final class EntropyCounterProbe {

	public static void main(String[] args){
		assert(args.length==0) : "JT011 counter fixtures use only DNA k1/k2 and three fixed small windows.";
		System.out.println("k\twindow\taccepted\tfilled_entropy\tfilled_monomer\tslid_entropy\tslid_monomer\tfailure\tmetric_agrees");
		for(int k=1; k<=2; k++){
			for(int window=32766; window<=32768; window++){observe(k, window);}
		}
	}

	private static void observe(int k, int window){
		assert(k>=1 && k<=2 && window>=32766 && window<=32768) :
			"Keep this counter probe finite: at most16 k-mer bins and32768 base slots.";
		boolean accepted=false;
		float filledEntropy=Float.NaN, filledMonomer=Float.NaN;
		float slidEntropy=Float.NaN, slidMonomer=Float.NaN;
		String failure="none", stage="constructor";
		try{
			final EntropyTracker tracker=new EntropyTracker(k, window, false);
			accepted=true;
			stage="fill";
			for(int i=0; i<window; i++){tracker.add((byte)'A');}
			filledEntropy=tracker.calcEntropy();
			filledMonomer=tracker.calcMaxMonomerFraction();
			stage="slide";
			tracker.add((byte)'A');
			slidEntropy=tracker.calcEntropy();
			slidMonomer=tracker.calcMaxMonomerFraction();
		}catch(RuntimeException | AssertionError e){
			failure=stage+":"+e.getClass().getName();
			System.err.println("JT011 k="+k+" window="+window+" stage="+stage);
			e.printStackTrace(System.err);
		}
		//Every completed all-A window has one k-mer species and one monomer species.
		final boolean agrees=accepted && failure.equals("none") &&
			Math.abs(filledEntropy)<1e-6 && Math.abs(slidEntropy)<1e-6 &&
			Math.abs(filledMonomer-1)<1e-6 && Math.abs(slidMonomer-1)<1e-6;
		System.out.println(k+"\t"+window+"\t"+accepted+"\t"+filledEntropy+"\t"+filledMonomer+
			"\t"+slidEntropy+"\t"+slidMonomer+"\t"+failure+"\t"+agrees);
	}
}
