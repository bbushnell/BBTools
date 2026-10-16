package test.structures;

import structures.SeqCount;
import structures.SeqCountM;

/** Literal equals-contract checks for SeqCount and SeqCountM sequence keys.
 * @author Jean
 */
public final class SeqCountEqualsContractCase{

	/** Runs current and proposed SeqCount.equals(Object) expectations. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation: "+args[0]);
		final String expect=(fixed ? "fixed" : "baseline");

		seqCountControls(expect);
		seqCountBoundaries(expect, fixed);
		seqCountMControls(expect);
	}

	/** Checks literal SeqCount equality controls through typed and Object calls. */
	private static void seqCountControls(final String expect){
		final SeqCount value=seqCount("ACGT");
		checkBoolean(expect, "seqcount-equal-typed", true, value.equals(seqCount("ACGT")));
		checkBoolean(expect, "seqcount-equal-object", true,
				value.equals((Object)seqCount("ACGT")));
		checkBoolean(expect, "seqcount-different-typed", false,
				value.equals(seqCount("ACGA")));
		checkBoolean(expect, "seqcount-different-object", false,
				value.equals((Object)seqCount("ACGA")));
	}

	/** Checks current SeqCount throws versus proposed Object.equals-compatible false. */
	private static void seqCountBoundaries(final String expect, final boolean fixed){
		checkObservation(expect, "seqcount-object-null", "false", observeSeqCountObject(null));
		checkObservation(expect, "seqcount-object-unrelated",
				fixed ? "false" : "ClassCastException",
				observeSeqCountObject(new Object()));
		checkObservation(expect, "seqcount-typed-null", "false", observeSeqCountTypedNull());
	}

	/** Checks inherited SeqCountM equality compatibility remains sequence-only. */
	private static void seqCountMControls(final String expect){
		final SeqCount base=seqCount("ACGT");
		final SeqCountM mutable=seqCountM("ACGT");
		checkBoolean(expect, "seqcountm-receiver-equal-typed", true,
				mutable.equals(seqCount("ACGT")));
		checkBoolean(expect, "seqcountm-receiver-equal-object", true,
				mutable.equals((Object)seqCount("ACGT")));
		checkBoolean(expect, "seqcountm-argument-equal-typed", true, base.equals(mutable));
		checkBoolean(expect, "seqcountm-different-typed", false,
				mutable.equals(seqCount("ACGA")));
		checkBoolean(expect, "seqcountm-different-object", false,
				mutable.equals((Object)seqCount("ACGA")));
	}

	/** Calls SeqCount.equals(Object) and returns a value or thrown type. */
	private static String observeSeqCountObject(final Object other){
		try{
			return Boolean.toString(seqCount("ACGT").equals(other));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Calls SeqCount.equals(SeqCount) with a typed null and returns a value or thrown type. */
	private static String observeSeqCountTypedNull(){
		try{
			return Boolean.toString(seqCount("ACGT").equals((SeqCount)null));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Creates a sequence key from literal A/C/G/T bases. */
	private static SeqCount seqCount(final String bases){
		return new SeqCount(bytes(bases));
	}

	/** Creates a mutable sequence key from literal A/C/G/T bases. */
	private static SeqCountM seqCountM(final String bases){
		return new SeqCountM(seqCount(bases));
	}

	/** Returns US-ASCII bytes for the supplied literal bases. */
	private static byte[] bytes(final String bases){
		return bases.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
	}

	/** Checks and prints a boolean outcome. */
	private static void checkBoolean(final String expect, final String name,
			final boolean expected, final boolean observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
		pass(expect, name, Boolean.toString(observed));
	}

	/** Checks and prints a value-or-exception outcome. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		pass(expect, name, observed);
	}

	/** Prints a checked outcome. */
	private static void pass(final String expect, final String name, final String observed){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Fails independently of assertion flags so an oracle cannot silently disappear. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
