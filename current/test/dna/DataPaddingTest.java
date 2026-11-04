package test.dna;

import dna.Data;

/** Finite padding controls and a bounded empty-symbol no-progress oracle.
 * Each case requires a fresh JVM. The faulty call runs only on a daemon worker,
 * so a failed two-second join cannot keep the test JVM alive after main exits.
 * This fixture does not claim a production cancellation or threading contract.
 * @author Sucrose
 */
public final class DataPaddingTest{

	/** Runs case=controls|empty-symbol with assertions enabled or disabled. */
	public static void main(String[] args) throws InterruptedException{
		require(args.length==1 && args[0].startsWith("case="), "Require one case= argument");
		final String mode=args[0].substring(5);
		if("controls".equals(mode)){controls();}
		else if("empty-symbol".equals(mode)){emptySymbol();}
		else{throw new IllegalArgumentException(mode);}
		System.out.println("PASS "+mode);
	}

	/** Literal outputs preserve whole-symbol concatenation and Java's null-string conversion. */
	private static void controls(){
		equal("single-symbol", Data.padFront("7", 3, "0"), "007");
		equal("whole-symbol-overshoot", Data.padFront("x", 4, "ab"), "ababx");
		equal("empty-text", Data.padFront("", 3, "ab"), "abab");
		equal("null-symbol", Data.padFront("x", 2, null), "nullx");
		equal("repeated-null-symbol", Data.padFront("x", 7, null), "nullnullx");
		final String original=new String("long");
		int unchanged=0;
		for(int width : new int[]{4, 3, 0, -1}){
			for(String symbol : new String[]{"", null, "ab"}){
				require(Data.padFront(original, width, symbol)==original,
						"No-padding call must return the original object at width "+width+", symbol="+symbol);
				unchanged++;
			}
		}
		final String empty=new String("");
		require(Data.padFront(empty, 0, "")==empty && Data.padFront(empty, 0, null)==empty,
				"Empty text at zero width needs no padding or symbol validation");
		int nullFailures=0;
		for(int width : new int[]{-1, 0, 1}){
			for(String symbol : new String[]{null, "", "x"}){
				RuntimeException failure=null;
				try{Data.padFront(null, width, symbol);}
				catch(RuntimeException problem){failure=problem;}
				require(failure instanceof NullPointerException,
						"Null num must retain NPE at width "+width+", symbol="+symbol+", failure="+failure);
				nullFailures++;
			}
		}
		System.out.println("noPaddingIdentity="+unchanged+" emptyIdentity=2 nullNumNPE="+nullFailures);
	}

	/** Prewarms Data before starting the only potentially nonterminating call. */
	private static void emptySymbol() throws InterruptedException{
		final String warm=new String("warm");
		require(Data.padFront(warm, 4, "")==warm, "Prewarm must finish without entering the padding loop");
		final PaddingCall call=new PaddingCall();
		final Thread worker=new Thread(call, "empty-padding-probe");
		worker.setDaemon(true);
		worker.start();
		worker.join(JOIN_MILLIS);
		final boolean alive=worker.isAlive();
		final Throwable failure=call.failure;
		final String result=call.result;
		System.out.println("terminated="+!alive+" result="+result+" exception="+
				(failure==null ? "none" : failure.getClass().getSimpleName()));
		require(!alive, "Empty symbol did not terminate within the bounded daemon-worker probe");
		require(failure instanceof IllegalArgumentException,
				"Empty symbol requiring padding must reject with IAE; result="+result+", failure="+failure);
	}

	/** Safely publishes observations even if the timed join returns before worker termination. */
	private static final class PaddingCall implements Runnable{
		@Override
		public void run(){
			try{result=Data.padFront("x", 2, "");}
			catch(Throwable problem){failure=problem;}
		}
		/** Normal return, if any; null until published. */
		volatile String result;
		/** Any thrown worker failure; only IAE satisfies the oracle. */
		volatile Throwable failure;
	}

	/** Prints and checks an independently specified finite result. */
	private static void equal(String label, String actual, String expected){
		System.out.println(label+" actual="+actual+" expected="+expected);
		require(expected.equals(actual), "Unexpected padding output for "+label);
	}

	/** Keeps the test oracle active under -da. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Maximum join duration in milliseconds; no unbounded wait or non-daemon probe follows. */
	private static final long JOIN_MILLIS=2000;
}
