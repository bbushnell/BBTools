package test.jgi;

import java.util.Arrays;

import jgi.ReformatReads;
import jgi.ReformatReadsOld;
import stream.Read;

/**
 * Direct nullable-ID controls for the two legacy Read conversion overloads.
 * Driver callers convert before their existing numeric-name fallback; a missing
 * ID must therefore remain available for that later caller-owned fallback.
 * This fixture calls only public static Read helpers, never a driver constructor,
 * main method, input parser, worker or CLI. String-null and Unicode policy are
 * outside this bounded nullable-Read test; named controls use literal ASCII IDs.
 * @author Sucrose
 */
public final class LegacyReformatNullableNameTest{

	/** Requires phase=baseline|fixed; both phases retain the same literal IDs and payload oracles. */
	public static void main(final String[] args){
		if(args.length!=1 || (!args[0].equals("phase=baseline") && !args[0].equals("phase=fixed"))){
			throw new IllegalArgumentException("Set phase=baseline|fixed.");
		}
		final boolean baseline=args[0].equals("phase=baseline");
		for(final boolean old : new boolean[]{false, true}){
			final String prefix=old ? "legacyold" : "legacy";
			one(old, baseline, prefix+"-null-read", true, null, null);
			one(old, baseline, prefix+"-null-id", false, null, null);
			one(old, baseline, prefix+"-empty-id", false, "", "");
			one(old, baseline, prefix+"-named", false, "alpha_beta_tail", "alpha_beta tail");
			pair(old, baseline, prefix+"-mixed-pair");
		}
		if(checks!=28){throw new AssertionError("Frozen legacy nullable-ID check count: "+checks+" != 28");}
		System.out.println("LEGACY_REFORMAT_NULLABLE_NAME_PASS phase="+(baseline ? "baseline" : "fixed")+" checks="+checks);
	}

	/** A missing Read is already a no-op; a valid Read can independently have a null ID. */
	private static void one(final boolean old, final boolean baseline, final String label,
		final boolean absent, final String input, final String expected){
		final Read read=absent ? null : read(input);
		NullPointerException failure=null;
		try{convert(old, read);}
		catch(NullPointerException e){failure=e;}
		expectFailure(old, label, failure, baseline && !absent && input==null);
		if(read!=null){
			check(same(expected, read.id), label+": unexpected resulting ID.");
			check(payload(read), label+": conversion changed bases, qualities or numeric ID.");
		}
		observe(label, failure, read, null);
	}

	/** Models just the two public helper calls, leaving numeric fallback explicitly to the driver. */
	private static void pair(final boolean old, final boolean baseline, final String label){
		final Read left=read("alpha_beta"), right=read(null);
		left.mate=right; right.mate=left;
		left.setPairnum(0); right.setPairnum(1);
		NullPointerException failure=null;
		try{convert(old, left); convert(old, right);}
		catch(NullPointerException e){failure=e;}
		expectFailure(old, label, failure, baseline);
		check("alpha beta".equals(left.id), label+": the first helper must convert the named mate.");
		check(payload(left), label+": left payload/numeric ID changed.");
		check(right.id==null && payload(right) && left.mate==right && right.mate==left &&
			left.pairnum()==0 && right.pairnum()==1, label+": missing ID or mate/payload state changed before fallback.");
		observe(label, failure, left, right);
	}

	/** Uses the actual public Read overloads, with no reflected or uninitialized driver instance. */
	private static void convert(final boolean old, final Read read){
		if(old){ReformatReadsOld.underscoreToSpace(read);}
		else{ReformatReads.underscoreToSpace(read);}
	}

	/** Baseline failures must originate in the selected legacy conversion helper. */
	private static void expectFailure(final boolean old, final String label,
		final NullPointerException failure, final boolean expected){
		final String type=old ? "jgi.ReformatReadsOld" : "jgi.ReformatReads";
		final boolean location=failure!=null && failure.getStackTrace().length>0 &&
			failure.getStackTrace()[0].getClassName().equals(type) &&
			failure.getStackTrace()[0].getMethodName().equals("underscoreToSpace");
		check(expected ? location : failure==null, label+": unexpected nullable-ID failure or origin.");
	}

	/** Fresh literal payloads and numeric ID42 are independent of all header conversion code. */
	private static Read read(final String id){return new Read(BASES.clone(), QUALITIES.clone(), id, 42);}

	/** Checks exact payload/numeric state after each direct call. */
	private static boolean payload(final Read read){
		return read.numericID==42 && Arrays.equals(BASES, read.bases) && Arrays.equals(QUALITIES, read.quality);
	}

	/** A null ID and an existing empty String are different inputs. */
	private static boolean same(final String expected, final String actual){return expected==null ? actual==null : expected.equals(actual);}

	/** Prints direct-call outcomes without converting missing IDs into fallback names. */
	private static void observe(final String label, final NullPointerException failure, final Read left, final Read right){
		System.out.println("legacy_nullable_case\t"+label+"\t"+(failure==null ? "ok" : "conversion-null-id")+
			"\tleft="+id(left)+"\tright="+id(right));
	}

	/** Makes absent Reads, missing IDs and empty IDs distinguishable in saved observations. */
	private static String id(final Read read){return read==null ? "<absent>" : read.id==null ? "<null>" : "["+read.id+"]";}

	/** Literal assertions remain active when JVM assertions are disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static final byte[] BASES={'A', 'C', 'G', 'T'}, QUALITIES={40, 40, 40, 40};
	private static int checks;
}
