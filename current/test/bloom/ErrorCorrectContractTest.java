package test.bloom;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.BitSet;

import bloom.ErrorCorrect;
import bloom.KCountArray;
import bloom.KCountArray7MTA;
import shared.Shared;
import stream.Read;

/** Finite legacy detection/error-count and optional-quality contracts.
 * Uses a real table with literal CCC/GGG keys, populated on both strands so
 * candidate orientation cannot confound the quality checks. Run in a fresh JVM.
 * Does not test CLI streams, burst limits or short-mate handling.
 * @author Brian Bushnell, Sucrose */
public final class ErrorCorrectContractTest{

	/** Reports all eleven cases before failing; checks remain active with assertions disabled. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("This finite test accepts no arguments.");}
		final int oldThreads=Shared.threads(), oldAdvance=ErrorCorrect.MIN_ADVANCE;
		final int oldBurst=ErrorCorrect.MAX_ERROR_BURST;
		final boolean oldInfo=ErrorCorrect.OUTPUT_INFO, oldVerbose=ErrorCorrect.verbose;
		int failures=0;
		try{
			Shared.setThreads(1);
			ErrorCorrect.MIN_ADVANCE=2;
			ErrorCorrect.MAX_ERROR_BURST=3;
			ErrorCorrect.OUTPUT_INFO=false;
			ErrorCorrect.verbose=false;
			KCountArray7MTA.setSeed(20261105);
			final KCountArray table=KCountArray.makeNew(2000003, 4, 3);
			table.increment(21L, 5);//CCC
			table.increment(42L, 5);//GGG, reverse complement of CCC
			for(int key=0; key<64; key++){
				check(table.read(key)==(key==21 || key==42 ? 5 : 0),
						"Every three-base key must match the literal table; collisions invalidate this fixture.");
			}
			failures+=detection(table, "CCC", 3, 0);
			failures+=detection(table, "CCCA", 3, 1);
			failures+=detection(table, "AAA", 0, 3);
			for(boolean quality : new boolean[]{true, false}){
				failures+=correction(table, quality, "CCCA", "CCCC", 4, 0, 4);
				failures+=correction(table, quality, "CCCAA", "CCCNA", 3, 2, 0);
			}
			for(boolean quality : new boolean[]{true, false}){
				failures+=fromRight(table, quality, "ACCC", "CCCC", 0, 4);
				failures+=fromRight(table, quality, "AACCC", "ANCCC", 2, 0);
			}
		}finally{
			Shared.setThreads(oldThreads);
			ErrorCorrect.MIN_ADVANCE=oldAdvance;
			ErrorCorrect.MAX_ERROR_BURST=oldBurst;
			ErrorCorrect.OUTPUT_INFO=oldInfo;
			ErrorCorrect.verbose=oldVerbose;
		}
		if(failures!=0){throw new AssertionError("Legacy error count/optional-quality contract failed "+failures+" of 11 cases.");}
		System.out.println("ERROR_CORRECT_CONTRACT_TEST_PASS cases=11");
	}

	/** Checks literal base coverage and a nonnegative count of unsupported bases. */
	private static int detection(final KCountArray table, final String sequence, final int supported, final int errors){
		final Read read=new Read(sequence.getBytes(StandardCharsets.US_ASCII), null, sequence, 0, false);
		final BitSet actual=ErrorCorrect.detectErrorsBulk(read, table, 3, 2, 1);
		final BitSet expected=new BitSet();
		expected.set(0, supported);
		final boolean pass=expected.equals(actual) && read.errors==errors &&
				Arrays.equals(read.bases, sequence.getBytes(StandardCharsets.US_ASCII));
		System.out.println("detection input="+sequence+" expectedCoverage="+expected+" actualCoverage="+actual+
				" expectedErrors="+errors+" actualErrors="+read.errors+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Compares identical correction inputs with and without a quality array; catches only the targeted null failure. */
	private static int correction(final KCountArray table, final boolean quality, final String sequence,
			final String expectedSequence, final int supported, final int errors, final int editedQuality){
		final byte[] quals=(quality ? new byte[sequence.length()] : null);
		if(quals!=null){Arrays.fill(quals, (byte)30);}
		final Read read=new Read(sequence.getBytes(StandardCharsets.US_ASCII), quals, sequence, 0, false);
		final BitSet initial=ErrorCorrect.detectErrorsBulk(read, table, 3, 2, 1);
		final BitSet initialExpected=new BitSet();
		initialExpected.set(0, 3);
		check(initialExpected.equals(initial), "Both correction inputs must start with exactly the supported CCC prefix.");
		BitSet actual=null;
		String exception="none";
		try{actual=ErrorCorrect.correctErrors(read, table, 3, 2, initial, 6, 3);}
		catch(NullPointerException e){exception=e.getClass().getSimpleName();}
		final BitSet expected=new BitSet();
		expected.set(0, supported);
		final byte[] expectedQuals=(quality ? new byte[sequence.length()] : null);
		if(expectedQuals!=null){Arrays.fill(expectedQuals, (byte)30); expectedQuals[3]=(byte)editedQuality;}
		final boolean pass=exception.equals("none") && expected.equals(actual) && read.errors==errors &&
				Arrays.equals(read.bases, expectedSequence.getBytes(StandardCharsets.US_ASCII)) && Arrays.equals(read.quality, expectedQuals);
		System.out.println("correction input="+sequence+" qualities="+quality+" expected="+expectedSequence+
				" actual="+new String(read.bases, StandardCharsets.US_ASCII)+" expectedCoverage="+expected+
				" actualCoverage="+actual+" expectedErrors="+errors+" actualErrors="+read.errors+
				" expectedQualities="+Arrays.toString(expectedQuals)+" actualQualities="+Arrays.toString(read.quality)+
				" exception="+exception+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Mirrors the two prefix-anchored cases using a supported CCC suffix; original seven cases stay unchanged. */
	private static int fromRight(final KCountArray table, final boolean quality, final String sequence,
			final String expectedSequence, final int errors, final int editedQuality){
		final int length=sequence.length(), edited=length-4;
		final byte[] quals=(quality ? new byte[length] : null);
		if(quals!=null){Arrays.fill(quals, (byte)30);}
		final Read read=new Read(sequence.getBytes(StandardCharsets.US_ASCII), quals, sequence, 0, false);
		final BitSet initial=ErrorCorrect.detectErrorsBulk(read, table, 3, 2, 1);
		final BitSet initialExpected=new BitSet();
		initialExpected.set(length-3, length);
		check(initialExpected.equals(initial), "Right-origin inputs must start with exactly the supported CCC suffix.");
		BitSet actual=null;
		String exception="none";
		try{actual=ErrorCorrect.correctErrors(read, table, 3, 2, initial, 6, 3);}
		catch(NullPointerException e){exception=e.getClass().getSimpleName();}
		final BitSet expected=new BitSet();
		expected.set(errors, length);
		final byte[] expectedQuals=(quality ? new byte[length] : null);
		if(expectedQuals!=null){Arrays.fill(expectedQuals, (byte)30); expectedQuals[edited]=(byte)editedQuality;}
		final boolean pass=exception.equals("none") && expected.equals(actual) && read.errors==errors &&
				Arrays.equals(read.bases, expectedSequence.getBytes(StandardCharsets.US_ASCII)) && Arrays.equals(read.quality, expectedQuals);
		System.out.println("rightOrigin input="+sequence+" qualities="+quality+" expected="+expectedSequence+
				" actual="+new String(read.bases, StandardCharsets.US_ASCII)+" expectedCoverage="+expected+
				" actualCoverage="+actual+" expectedErrors="+errors+" actualErrors="+read.errors+
				" expectedQualities="+Arrays.toString(expectedQuals)+" actualQualities="+Arrays.toString(read.quality)+
				" exception="+exception+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Makes invalid fixture setup fail under both assertion modes. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
