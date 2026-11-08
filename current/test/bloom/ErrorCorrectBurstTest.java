package test.bloom;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.BitSet;

import bloom.ErrorCorrect;
import bloom.KCountArray;
import bloom.KCountArray7MTA;
import shared.Shared;
import stream.Read;

/** Direct API contract for maxBurst; current production callers pass the global value.
 * This deliberately supplies a different argument and makes no CLI-divergence claim.
 * @author Brian Bushnell, Sucrose */
public final class ErrorCorrectBurstTest{

	/** Compares burst limits one and three against the same two-adjacent-error input. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("This finite test accepts no arguments.");}
		final int oldThreads=Shared.threads(), oldAdvance=ErrorCorrect.MIN_ADVANCE;
		final int oldBurst=ErrorCorrect.MAX_ERROR_BURST, oldDistance=ErrorCorrect.BURST_THRESH;
		final boolean oldInfo=ErrorCorrect.OUTPUT_INFO, oldVerbose=ErrorCorrect.verbose;
		int failures=0;
		try{
			Shared.setThreads(1);
			ErrorCorrect.MIN_ADVANCE=1; ErrorCorrect.MAX_ERROR_BURST=3; ErrorCorrect.BURST_THRESH=2;
			ErrorCorrect.OUTPUT_INFO=ErrorCorrect.verbose=false;
			KCountArray7MTA.setSeed(20261105);
			final KCountArray table=KCountArray.makeNew(2000003, 4, 3);
			table.increment(341, 5); table.increment(682, 5);//CCCCC and GGGGG
			for(int key=0; key<1024; key++){
				require(table.read(key)==(key==341 || key==682 ? 5 : 0), "All five-base keys must match the literal table.");
			}
			for(int limit : new int[]{1, 3}){
				final byte[] quals={30, 30, 30, 30, 30, 30, 30};
				final Read read=new Read("AACCCCC".getBytes(StandardCharsets.US_ASCII), quals, "burst", 0, false);
				final BitSet initial=ErrorCorrect.detectErrorsBulk(read, table, 5, 2, 1);
				final BitSet expectedInitial=new BitSet(); expectedInitial.set(2, 7);
				require(expectedInitial.equals(initial), "Only the terminal CCCCC window supports the initial read.");
				final BitSet actual=ErrorCorrect.correctErrors(read, table, 5, 2, initial, 6, limit);
				final boolean limited=limit==1;
				final String sequence=(limited ? "NCCCCCC" : "CCCCCCC");
				final byte[] expectedQuals={(byte)(limited ? 0 : 4), 4, 30, 30, 30, 30, 30};
				final BitSet expected=new BitSet(); expected.set(limited ? 1 : 0, 7);
				final boolean pass=expected.equals(actual) && read.errors==(limited ? 1 : 0) &&
						Arrays.equals(read.bases, sequence.getBytes(StandardCharsets.US_ASCII)) && Arrays.equals(read.quality, expectedQuals);
				System.out.println("apiOnly maxBurst="+limit+" globalBurst=3 expected="+sequence+
						" actual="+new String(read.bases, StandardCharsets.US_ASCII)+" expectedCoverage="+expected+
						" actualCoverage="+actual+" errors="+read.errors+" expectedQualities="+Arrays.toString(expectedQuals)+
						" actualQualities="+Arrays.toString(read.quality)+" pass="+pass);
				if(!pass){failures++;}
			}
		}finally{
			Shared.setThreads(oldThreads);
			ErrorCorrect.MIN_ADVANCE=oldAdvance; ErrorCorrect.MAX_ERROR_BURST=oldBurst; ErrorCorrect.BURST_THRESH=oldDistance;
			ErrorCorrect.OUTPUT_INFO=oldInfo; ErrorCorrect.verbose=oldVerbose;
		}
		if(failures>0){throw new AssertionError("Passed maxBurst contract failed "+failures+" of 2 API-only cases.");}
		System.out.println("ERROR_CORRECT_BURST_TEST_PASS cases=2 apiOnly=true");
	}

	/** Rejects a faulty fixture independently of the JVM assertion setting. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
