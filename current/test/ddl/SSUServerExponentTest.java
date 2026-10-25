package test.ddl;

import cardinality.DynamicDemiLog;
import ddl.SSUServer;
import stream.Read;

/** Checks server option parsing without loading references or starting a listener.
 * @author Shinobu
 */
public final class SSUServerExponentTest{

	/** Runs the captured baseline or the corrected argument contract. */
	public static void main(String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed=args[0].equals("expect=fixed");
		check(fixed || args[0].equals("expect=baseline"), "Unknown expectation: "+args[0]);
		final int initialExponent=DynamicDemiLog.exponentBits();
		final boolean initialUtoT=Read.U_TO_T;
		try{
			checkCase(fixed, "default", 4);
			checkCase(fixed, "explicit-default", 4, "exponent=4");
			checkCase(fixed, "exponent", 5, "exponent=5");
			checkCase(fixed, "alias", 6, "ebits=6");
			checkCase(fixed, "case-insensitive", 5, "EXPONENT=5");
			checkCase(fixed, "last-option-wins", 6, "exponent=5", "ebits=6");
			checkCase(fixed, "minimum", 1, "exponent=1");
			checkCase(fixed, "maximum", 8, "ebits=8");
			checkCase(fixed, "new-default-after-custom", 4);
		}finally{
			DynamicDemiLog.setExponent(initialExponent);
			Read.U_TO_T=initialUtoT;
		}
	}

	/** Checks the process-global setting consumed by subsequent sketch construction. */
	private static void checkCase(boolean fixed, String name, int requested, String... options){
		new SSUServer(options);
		final int observed=DynamicDemiLog.exponentBits();
		final int expected=fixed ? requested : 4;
		check(observed==expected, name+": expected exponent "+expected+", observed "+observed);
		System.out.println("PASS\t"+name+"\trequested="+requested+"\tobserved="+observed);
	}

	/** Keeps oracle checks active under both assertion settings. */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
