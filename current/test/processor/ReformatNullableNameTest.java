package test.processor;

import java.util.Arrays;

import parse.Parser;
import processor.ReformatProcessor;
import stream.Read;

/**
 * Direct-call controls for nullable IDs before last-underscore conversion.
 * Null IDs are retained by Read and can result from trimmed empty FASTQ headers.
 * Baseline expects null-ID conversion to throw at underscoreToSpace;
 * fixed mode expects the Read overload to leave a missing ID for existing fallback.
 * No driver, reader thread or CLI process is started.
 * @author Sucrose
 */
public final class ReformatNullableNameTest{

	/** Requires phase=baseline|fixed; observations and the final check count go to stdout. */
	public static void main(final String[] args){
		if(args.length!=1 || (!args[0].equals("phase=baseline") && !args[0].equals("phase=fixed"))){
			throw new IllegalArgumentException("Set phase=baseline|fixed.");
		}
		final boolean baseline=args[0].equals("phase=baseline");
		boolean failed=false;
		try{ReformatProcessor.underscoreToSpace((Read)null);}
		catch(NullPointerException e){failed=true;}
		check(!failed, "A missing Read is already a supported no-op.");
		System.out.println("nullable_case\tnull-read\tok");
		direct(baseline, "direct-null-id", null, null, true);
		direct(baseline, "direct-empty-id", "", "", false);
		direct(baseline, "direct-named", "alpha_beta_tail", "alpha_beta tail", false);
		process(baseline, "disabled-null", null, null, false, false, false, null, null);
		process(baseline, "disabled-numeric-fallback", null, null, false, false, true, "42", null);
		process(baseline, "enabled-null", null, null, false, true, false, null, null);
		process(baseline, "enabled-numeric-fallback", null, null, false, true, true, baseline ? null : "42", null);
		process(baseline, "enabled-empty", "", null, false, true, false, "", null);
		process(baseline, "enabled-mixed-pair", "alpha_beta", null, true, true, true,
			"alpha beta", baseline ? null : "alpha beta");
		final int expected=baseline ? 32 : 35;
		if(checks!=expected){throw new AssertionError("Frozen nullable-ID check count: "+checks+" != "+expected);}
		System.out.println("REFORMAT_NULLABLE_NAME_PASS phase="+(baseline ? "baseline" : "fixed")+" checks="+checks);
	}

	/** Calls only the Read overload; the String overload retains its documented nonnull precondition. */
	private static void direct(final boolean baseline, final String label, final String input,
		final String expected, final boolean defective){
		final Read read=read(input);
		NullPointerException failure=null;
		try{ReformatProcessor.underscoreToSpace(read);}
		catch(NullPointerException e){failure=e;}
		expectFailure(label, failure, baseline && defective);
		check(same(expected, read.id), label+": unexpected ID "+read.id);
		check(payload(read), label+": conversion changed bases, qualities or numeric ID.");
		observe(label, failure, read, null);
	}

	/** Calls the processor synchronously, retaining normal disabled-option/fallback controls. */
	private static void process(final boolean baseline, final String label, final String first,
		final String second, final boolean paired, final boolean convert, final boolean unique,
		final String expectedFirst, final String expectedSecond){
		final ReformatProcessor proc=new ReformatProcessor();
		proc.setFromParser(new Parser());
		proc.minReadLength=0; proc.maxReadLength=0;
		proc.lastUnderscoreToSpace=convert; proc.uniqueNames=unique;
		proc.postParse();
		final Read left=read(first), right=paired ? read(second) : null;
		if(paired){left.mate=right; right.mate=left; left.setPairnum(0); right.setPairnum(1);}
		int selected=-1;
		NullPointerException failure=null;
		try{selected=proc.processReadPair(left, right);}
		catch(NullPointerException e){failure=e;}
		final boolean defective=convert && (first==null || (paired && second==null));
		expectFailure(label, failure, baseline && defective);
		check(same(expectedFirst, left.id) && (!paired || same(expectedSecond, right.id)), label+": unexpected resulting ID.");
		check(payload(left), label+": left payload/numeric ID changed.");
		if(paired){
			check(payload(right) && left.mate==right && right.mate==left && left.pairnum()==0 && right.pairnum()==1,
				label+": right payload/numeric ID or reciprocal mate state changed.");
		}
		if(failure==null){check(selected==(paired ? 3 : 1), label+": valid reads must survive unchanged length filters.");}
		observe(label, failure, left, right);
	}

	/** Baseline failures must originate in this specific conversion helper, not another processor path. */
	private static void expectFailure(final String label, final NullPointerException failure, final boolean expected){
		final boolean location=failure!=null && failure.getStackTrace().length>0 &&
			failure.getStackTrace()[0].getClassName().equals("processor.ReformatProcessor") &&
			failure.getStackTrace()[0].getMethodName().equals("underscoreToSpace");
		check(expected ? location : failure==null, label+": unexpected null-ID failure outcome.");
	}

	/** Every fresh read uses a valid four-base payload and an independent numeric fallback42. */
	private static Read read(final String id){return new Read(BASES.clone(), QUALITIES.clone(), id, 42);}

	/** The production conversion must preserve these independent literal payload/numeric properties. */
	private static boolean payload(final Read read){
		return read.numericID==42 && Arrays.equals(BASES, read.bases) && Arrays.equals(QUALITIES, read.quality);
	}

	/** Null is distinct from an existing empty header. */
	private static boolean same(final String expected, final String actual){return expected==null ? actual==null : expected.equals(actual);}

	/** Emits the direct call's outcome and resulting fixture IDs without starting a worker. */
	private static void observe(final String label, final NullPointerException failure, final Read left, final Read right){
		System.out.println("nullable_case\t"+label+"\t"+(failure==null ? "ok" : "conversion-null-id")+
			"\tleft="+(left.id==null ? "<null>" : "["+left.id+"]")+
			"\tright="+(right==null ? "<absent>" : right.id==null ? "<null>" : "["+right.id+"]"));
	}

	/** Literal expectations remain active with or without JVM assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static final byte[] BASES={'A', 'C', 'G', 'T'}, QUALITIES={40, 40, 40, 40};
	private static int checks;
}
