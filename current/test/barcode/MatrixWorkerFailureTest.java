package test.barcode;

import java.util.ArrayList;
import java.util.Locale;

import barcode.Barcode;
import barcode.PCRMatrixHDist;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Exercises an actual Hamming-matrix worker failure and a matching valid control.
 * Two workers receive alternating entries; a bad first entry kills worker zero
 * in the baseline, leaving exactly 256 observations from worker one.
 * @author Brian Bushnell, Shinobu
 */
public final class MatrixWorkerFailureTest{

	/** Requires mode=baseline|fixed; assertions must be enabled. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, MatrixWorkerFailureTest.class, false);
		try{
			String mode=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			if(!"baseline".equals(mode) && !"fixed".equals(mode)){throw new IllegalArgumentException("Set mode=baseline or mode=fixed.");}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("Worker-failure reproduction requires -ea.");}
			Shared.setThreads(2);
			final Matrix control=new Matrix();
			control.populateCounts(queries(false), 1);
			check(control.counted()==512 && !control.errorState, "Valid two-worker control must count all 512 records.");
			final Matrix broken=new Matrix();
			boolean threw=false;
			try{broken.populateCounts(queries(true), 1);}
			catch(final RuntimeException e){threw=true;}
			if(mode.equals("baseline")){
				check(!threw && !broken.errorState && broken.counted()==256,
						"Baseline must demonstrate silent partial success: threw="+threw+", error="+broken.errorState+", count="+broken.counted());
			}else{
				check(threw && broken.errorState, "Worker failure must propagate with errorState set.");
			}
			System.out.println("MATRIX_WORKER mode="+mode+" control="+control.counted()+" threw="+threw+
					" error="+broken.errorState+" partial="+broken.counted());
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; one process owns one finite reproduction. */
	private MatrixWorkerFailureTest(){}

	/** Makes 512 valid observations or replaces only worker zero's first observation. */
	private static ArrayList<Barcode> queries(final boolean bad){
		final ArrayList<Barcode> list=new ArrayList<Barcode>(512);
		for(int i=0; i<512; i++){list.add(new Barcode(bad && i==0 ? "A" : "AAAAAAAA", 1));}
		return list;
	}

	/** Exposes only the protected total for the independent finite-count oracle. */
	private static final class Matrix extends PCRMatrixHDist{
		/** Two hundred distinct eight-base references force the two-worker workload path. */
		Matrix(){
			super(8, 0, 0, false);
			final ArrayList<String> expected=new ArrayList<String>(200);
			final char[] bases={'A', 'C', 'G', 'T'};
			for(int code=0; code<200; code++){
				final char[] sequence=new char[8];
				int value=code;
				for(int pos=7; pos>=0; pos--){sequence[pos]=bases[value&3]; value>>>=2;}
				expected.add(new String(sequence));
			}
			populateExpected(expected);
			populateSplitCodes();
		}

		/** Number of observations actually processed, including rejected assignments. */
		long counted(){return totalCounted;}
	}

	/** Always enforce oracle checks, independent of assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
