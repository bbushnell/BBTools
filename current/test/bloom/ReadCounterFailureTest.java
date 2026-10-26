package test.bloom;

import java.io.IOException;
import java.util.Arrays;

import bloom.KCountArray;
import bloom.ReadCounter;
import stream.FastaReadInputStream;

/** Ensures failed input passes cannot masquerade as complete prefilter counts.
 * @author Brian Bushnell, Fischl */
public final class ReadCounterFailureTest {

	/** Covers every primary/extra pass, checked/runtime exceptions and Errors. */
	public static void main(final String[] args){
		for(int passes=1; passes<=3; passes++){
			for(int failure=1; failure<=2*passes; failure++){
				for(int kind=0; kind<3; kind++){exercise(passes, failure, kind);}
			}
			exercise(passes, 0, 0);
		}
		System.out.println("READ_COUNTER_FAILURE_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Test Methods          ----------------*/
	/*--------------------------------------------------------------*/

	/** Verifies no partial result escapes and all temporary global settings recover. */
	private static void exercise(final int passes, final int failure, final int kind){
		assert(passes>0 && failure>=0 && failure<=2*passes) : "Each pass has one primary and one extra input.";
		final byte oldq=ReadCounter.minQuality;
		final long oldmax=ReadCounter.maxReads;
		final boolean oldsplit=FastaReadInputStream.SPLIT_READS;
		ReadCounter.minQuality=7;
		ReadCounter.maxReads=123;
		FastaReadInputStream.SPLIT_READS=false;
		final InjectedCounter counter=new InjectedCounter(failure, kind);
		Throwable caught=null;
		KCountArray result=null;
		try{
			try{
				result=counter.makeKca("primary.fq", null, Arrays.asList("extra.fq"),
						2, 4096, 2, 0, 19, passes, 1, 1, null, 0);
			}catch(Throwable e){caught=e;}
			check(ReadCounter.minQuality==7 && ReadCounter.maxReads==123
					&& !FastaReadInputStream.SPLIT_READS,
					"Success and every failure kind must restore the caller's counting settings.");
			if(failure>0){
				check(result==null && caught!=null, "A failed input must not return a partial count array.");
				check(counter.calls==failure, "Later input passes must not run after failure.");
				while(caught.getCause()!=null){caught=caught.getCause();}
				check(caught==counter.failure, "The original failure must remain available to diagnose incomplete counts.");
			}else{
				check(caught==null && result!=null, "Successful passes must still return their array.");
				check(counter.calls==2*passes, "Every primary and extra pass must execute exactly once.");
				check(result.read(1234)==2, "The final array must contain both successful input contributions.");
			}
		}finally{
			ReadCounter.minQuality=oldq;
			ReadCounter.maxReads=oldmax;
			FastaReadInputStream.SPLIT_READS=oldsplit;
		}
	}

	/** Unit failures remain fatal without relying on the launcher's assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Test Counter          ----------------*/
	/*--------------------------------------------------------------*/

	/** Isolates makeKca's failure contract from filesystem and streamer behavior. */
	private static final class InjectedCounter extends ReadCounter {

		/** Selects an input-call ordinal and the exception family to inject there. */
		InjectedCounter(final int failAt_, final int kind){
			super(64);
			failAt=failAt_;
			failure=(kind==0 ? new IOException("Injected input failure") :
					kind==1 ? new IllegalStateException("Injected incomplete counts") :
					new AssertionError("Injected counting invariant failure"));
		}

		@Override
		public KCountArray count(final String first, final String second, final KCountArray counts) throws Exception{
			return contribute(first, counts);
		}

		@Override
		public KCountArray count(final String first, final String second, final int bits,
				final KCountArray counts, final KCountArray trusted, final long limit,
				final int threshold, final boolean conservative) throws Exception{
			return contribute(first, counts);
		}

		/** Adds one observable count, then fails only on the selected input call. */
		private KCountArray contribute(final String first, final KCountArray counts) throws Exception{
			calls++;
			check(first.equals((calls&1)==1 ? "primary.fq" : "extra.fq"), "Input order changed across passes.");
			check(minQuality==0, "The requested temporary quality setting must apply during counting.");
			counts.increment(1234);
			if(calls==failAt){
				//The wrapper must restore this global even for an Error, not just Exceptions.
				FastaReadInputStream.SPLIT_READS=true;
				if(failure instanceof Error){throw (Error)failure;}
				throw (Exception)failure;
			}
			return counts;
		}

		private final int failAt;
		private final Throwable failure;
		private int calls;
	}

	private static int checks;
}
