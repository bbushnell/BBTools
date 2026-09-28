package prok;

/**
 * Focused tests for P1's per-family ncRNA window-pad CLI overrides (Brian via Citan,
 * 2026-08-28): rnaseppad=/srpsmallpad=/srplargepad=. Tests the two methods factored out
 * of parse()/loadNcrnaResources() for direct testability (parsePadOverride, resolvePad) --
 * no real staged resource files or genome data needed.
 *
 * <p>Restored from G11's C3 fixture; also checks current targeted-sweep precedence
 * and isolation. Run with -ea. Helper defaults below are synthetic examples.
 * @author G11, Raiden
 */
public class CallGenesPadOverrideTest {

	public static void main(String[] args){
		int failures=0;
		failures+=testDefaultSentinelUsesShippedDefault() ? 0 : 1;
		failures+=testOverrideValueIsUsed() ? 0 : 1;
		failures+=testZeroOverrideIsNotTreatedAsUnset() ? 0 : 1;
		failures+=testNegativeValueRejectedLoud() ? 0 : 1;
		failures+=testNonNumericValueRejectedLoud() ? 0 : 1;
		failures+=testAllThreeFamiliesIndependentSentinels() ? 0 : 1;
		testTargetedSweepPrecedence();
		if(failures>0){
			throw new AssertionError(failures+" padding checks failed");
		}
		System.err.println("CallGenesPadOverrideTest: ALL TESTS PASSED");
	}

	/** The -1 sentinel (flag never passed) must resolve to the family's shipped default,
	 * byte-identical to pre-P1 behavior. */
	private static boolean testDefaultSentinelUsesShippedDefault(){
		final boolean ok=CallGenes.resolvePad(-1, 450)==450
			&& CallGenes.resolvePad(-1, 100)==100
			&& CallGenes.resolvePad(-1, 350)==350;
		System.err.println("testDefaultSentinelUsesShippedDefault "+(ok ? "PASSED" : "FAILED"));
		return ok;
	}

	/** A non-negative override must win over the shipped default. */
	private static boolean testOverrideValueIsUsed(){
		final boolean ok=CallGenes.resolvePad(600, 450)==600 && CallGenes.resolvePad(1, 100)==1;
		System.err.println("testOverrideValueIsUsed "+(ok ? "PASSED" : "FAILED"));
		return ok;
	}

	/** 0 is a legitimate pad value, not a second "unset" sentinel -- must NOT fall back to
	 * the shipped default. */
	private static boolean testZeroOverrideIsNotTreatedAsUnset(){
		final boolean ok=CallGenes.resolvePad(0, 450)==0;
		System.err.println("testZeroOverrideIsNotTreatedAsUnset "+(ok ? "PASSED" : "FAILED (returned "+CallGenes.resolvePad(0,450)+", expected 0)"));
		return ok;
	}

	/** A negative CLI value must fail loud with a diagnostic naming the flag, not silently
	 * clamp or accept it. */
	private static boolean testNegativeValueRejectedLoud(){
		try{
			CallGenes.parsePadOverride("rnaseppad", "-5");
			System.err.println("testNegativeValueRejectedLoud FAILED: expected IllegalArgumentException, none thrown");
			return false;
		}catch(IllegalArgumentException e){
			final boolean ok=e.getMessage()!=null && e.getMessage().contains("rnaseppad") && e.getMessage().contains(">=0");
			System.err.println("testNegativeValueRejectedLoud "+(ok ? "PASSED" : "FAILED (wrong message): "+e.getMessage()));
			return ok;
		}
	}

	/** A non-numeric CLI value must fail loud (NumberFormatException from the underlying
	 * Integer.parseInt), not silently default to 0 or -1. */
	private static boolean testNonNumericValueRejectedLoud(){
		try{
			CallGenes.parsePadOverride("srpsmallpad", "abc");
			System.err.println("testNonNumericValueRejectedLoud FAILED: expected NumberFormatException, none thrown");
			return false;
		}catch(NumberFormatException e){
			System.err.println("testNonNumericValueRejectedLoud PASSED");
			return true;
		}
	}

	/** Setting one family's override must not affect the sentinel value CallGenes.parse()
	 * would use for the other two -- each flag maps to its own independent static field. */
	private static boolean testAllThreeFamiliesIndependentSentinels(){
		final int rnasep=CallGenes.parsePadOverride("rnaseppad", "999");
		final boolean ok=rnasep==999
			&& CallGenes.resolvePad(-1, 100)==100   //srp_small default unaffected
			&& CallGenes.resolvePad(-1, 350)==350;  //srp_large default unaffected
		System.err.println("testAllThreeFamiliesIndependentSentinels "+(ok ? "PASSED" : "FAILED"));
		return ok;
	}

	/** A targeted pad wins only for its named family; unset preserves per-family pads. */
	private static void testTargetedSweepPrecedence(){
		final String oldFamily=CallGenes.NCRNA_FAMILY_FILTER;
		final int oldPad=CallGenes.NCRNA_WINDOW_PAD_OVERRIDE;
		try{
			CallGenes.NCRNA_FAMILY_FILTER="rnasep";
			CallGenes.NCRNA_WINDOW_PAD_OVERRIDE=0;
			if(CallGenes.resolveSweepPad("rnasep", 600, 320)!=0
					|| CallGenes.resolveSweepPad("srp_small", 80, 75)!=80
					|| CallGenes.resolveSweepPad("srp_large", -1, 250)!=250){
				throw new AssertionError("Targeted padding leaked to another family or lost zero override");
			}
			CallGenes.NCRNA_WINDOW_PAD_OVERRIDE=-1;
			if(CallGenes.resolveSweepPad("rnasep", 600, 320)!=600
					|| CallGenes.resolveSweepPad("rnasep", -1, 320)!=320){
				throw new AssertionError("Unset sweep pad must preserve family override or shipped fallback");
			}
		}finally{
			CallGenes.NCRNA_FAMILY_FILTER=oldFamily;
			CallGenes.NCRNA_WINDOW_PAD_OVERRIDE=oldPad;
		}
	}
}
