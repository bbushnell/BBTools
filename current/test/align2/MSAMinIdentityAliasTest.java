package test.align2;

import align2.MSA;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Checks six literal minimum-identity ratio contracts without constructing aligners.
 * Baseline mode requires the two original explicit-JNI dispatch assertions;
 * fixed mode requires those aliases to return the ordinary 11ts ratio.
 * Four controls preserve fractional, percent, perfect and minimum-ratio behavior.
 * No native implementation is initialized and no sequence alignment is performed.
 * Run in an isolated JVM with assertions enabled and an explicit mode.
 * @author Brian Bushnell, Sucrose
 */
public final class MSAMinIdentityAliasTest{

	/** Accepts mode=baseline or mode=fixed and normal launcher JVM flags. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, MSAMinIdentityAliasTest.class, false);
		try{
			String mode=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2 || !pair[0].equalsIgnoreCase("mode") ||
						(!pair[1].equalsIgnoreCase("baseline") && !pair[1].equalsIgnoreCase("fixed"))){
					throw new IllegalArgumentException("Expected mode=baseline or mode=fixed: "+arg);
				}
				mode=pair[1];
			}
			if(mode==null){throw new IllegalArgumentException("Specify mode=baseline or mode=fixed explicitly.");}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("Alias regression checks require -ea to observe baseline dispatch failures.");}
			runChecks(mode.equalsIgnoreCase("baseline"));
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Static fixture entry point only. */
	private MSAMinIdentityAliasTest(){}

	/** Uses literal expectations derived independently from the existing 11ts constants. */
	private static void runChecks(final boolean baseline){
		final ByteBuilder out=new ByteBuilder(512);
		out.append("case\texpected\tobserved\tpassed\n");
		int checks=0, assertions=0;
		checks+=checkRatio(out, "perfect_control", "MultiStateAligner11ts", 1.0, 1.0f);
		checks+=checkRatio(out, "fractional_control", "MultiStateAligner11ts", 0.75, 0.5405f);
		checks+=checkRatio(out, "percent_control", "MultiStateAligner11ts", 75.0, 0.5405f);
		checks+=checkRatio(out, "floor_control", "MultiStateAligner11ts", 0.5, 0.1f);
		if(baseline){
			assertions+=checkUnhandled(out, "explicit_jni_alias", "MultiStateAligner11tsJNI", 0.75);
			assertions+=checkUnhandled(out, "mixed_jni_percent", "mUlTiStAtEaLiGnEr11TsJnI", 75.0);
			checks+=assertions;
		}else{
			checks+=checkRatio(out, "explicit_jni_alias", "MultiStateAligner11tsJNI", 0.75, 0.5405f);
			checks+=checkRatio(out, "mixed_jni_percent", "mUlTiStAtEaLiGnEr11TsJnI", 75.0, 0.5405f);
		}
		assert(checks==6) : "The saved regression contract requires four controls and two alias cases; checks="+checks;
		System.out.print(out);
		System.out.println("MSA_MIN_IDENTITY_ALIAS mode="+(baseline ? "baseline" : "fixed")+
				" checks="+checks+" assertions="+assertions+" failed=0");
	}

	/** Appends one exact float-bit comparison; the expected value is a caller literal. */
	private static int checkRatio(final ByteBuilder out, final String name, final String classname,
			final double identity, final float expected){
		assert(out!=null && name!=null && classname!=null && identity>0 && Float.isFinite(expected)) :
			"Ratio fixtures require a named valid input and finite literal expectation: "+name;
		final float observed=MSA.minIdToMinRatio(identity, classname);
		final boolean passed=Float.floatToIntBits(observed)==Float.floatToIntBits(expected);
		out.append(name).tab().append(Float.toString(expected)).tab().append(Float.toString(observed))
				.tab().append(passed ? 1 : 0).nl();
		if(!passed){
			System.out.print(out);
			throw new AssertionError("Minimum-identity ratio differs for "+name+": expected="+expected+" observed="+observed);
		}
		return 1;
	}

	/** Accepts only the original dispatch assertion, never an unrelated assertion or error. */
	private static int checkUnhandled(final ByteBuilder out, final String name, final String classname,
			final double identity){
		assert(out!=null && name!=null && classname!=null && identity>0) :
			"Baseline alias fixtures require a named valid input: "+name;
		final float observed;
		try{
			observed=MSA.minIdToMinRatio(identity, classname);
		}catch(final AssertionError error){
			if(!("Unhandled MSA type: "+classname).equals(error.getMessage())){throw error;}
			out.append(name).append("\tAssertionError\tAssertionError\t1\n");
			return 1;
		}
		out.append(name).append("\tAssertionError\t").append(Float.toString(observed)).append("\t0\n");
		System.out.print(out);
		throw new AssertionError("Original JNI alias unexpectedly returned a ratio: "+name+" observed="+observed);
	}
}
