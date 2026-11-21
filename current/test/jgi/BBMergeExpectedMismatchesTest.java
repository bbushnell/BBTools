package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;

import jgi.BBMergeOverlapper;
import parse.Parse;
import stream.Read;

/** Six direct-helper controls for the insert-coordinate quality sum (JGI-083).
 * Coordinates and decimal sums are literal oracles, not calculated by the helper.
 * Reads already have the orientation expected after reverse-complementing mate B.
 * No merger search, CLI pipeline, model inference or feature mutation is exercised.
 * @author Sucrose
 */
public final class BBMergeExpectedMismatchesTest{

	/** Requires phase=baseline|fixed and a new evidence directory. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown mismatch fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergeexpectedmismatchestest.sh with assertions enabled.");
		final BBMergeExpectedMismatchesTest test=new BBMergeExpectedMismatchesTest(
			Paths.get(directory).toAbsolutePath().normalize(), phase.equals("baseline"));
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==6, "Expected six literal rows: "+test.cases);
		test.check(test.calls==6, "Expected six direct helper calls: "+test.calls);
		test.check(test.changed==2, "Expected two affected unequal-length rows: "+test.changed);
		test.check(test.controls==4, "Expected four unchanged controls: "+test.controls);
		test.check(test.sourceMatches==(test.baseline ? 6 : 4), "Unexpected source matches: "+test.sourceMatches);
		test.check(test.correctedMatches==(test.baseline ? 4 : 6), "Unexpected corrected matches: "+test.correctedMatches);
		final String summary="BBMERGE_EXPECTED_MISMATCHES_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" direct="+test.calls+" changed="+test.changed+
			" controls="+test.controls+" sourceMatches="+test.sourceMatches+
			" correctedMatches="+test.correctedMatches+" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(test.out.resolve("summary.txt"), summary.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+test.out);}
		System.out.print(summary);
	}

	/** Creates only a new output directory; existing evidence is never overwritten. */
	private BBMergeExpectedMismatchesTest(final Path out_, final boolean baseline_) throws IOException{
		out=out_; baseline=baseline_;
		Files.createDirectory(out);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\tinsert\taBases\tbBases\taQualities\tbQualities"+
			"\tsourcePairsOracle\tcorrectedPairsOracle\tsourceSumOracle\tcorrectedSumOracle"+
			"\tobservedFloat\tobservedRawBitsHex\tsourceMatches\tcorrectedMatches\tinputsUnchanged\n")
			.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Calls only the real public helper, retaining floats and both independent oracles. */
	private void run(final Row row) throws IOException{
		cases++;
		if(row.changed){changed++;}else{controls++;}
		final byte[] ab=row.aBases.getBytes(StandardCharsets.US_ASCII), bb=row.bBases.getBytes(StandardCharsets.US_ASCII);
		final byte[] aq=row.aQualities.clone(), bq=row.bQualities.clone();
		final byte[] abBefore=ab.clone(), bbBefore=bb.clone(), aqBefore=aq.clone(), bqBefore=bq.clone();
		final Read a=new Read(ab, aq, row.name+"/A", 0, false);
		final Read b=new Read(bb, bq, row.name+"/B-rc", 1, false);
		calls++;
		final float observed=BBMergeOverlapper.expectedMismatches(a, b, row.insert);
		final boolean finite=Float.isFinite(observed) && observed>=0;
		final boolean sourceMatch=finite && Math.abs((double)observed-row.sourceSum)<=TOLERANCE;
		final boolean correctedMatch=finite && Math.abs((double)observed-row.correctedSum)<=TOLERANCE;
		final boolean unchanged=a.bases==ab && b.bases==bb && a.quality==aq && b.quality==bq &&
			Arrays.equals(abBefore, ab) && Arrays.equals(bbBefore, bb) &&
			Arrays.equals(aqBefore, aq) && Arrays.equals(bqBefore, bq);
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+"\t"+row.insert+
			"\t"+row.aBases+"\t"+row.bBases+"\t"+Arrays.toString(aqBefore)+"\t"+Arrays.toString(bqBefore)+
			"\t"+row.sourcePairs+"\t"+row.correctedPairs+"\t"+row.sourceSum+"\t"+row.correctedSum+
			"\t"+Float.toString(observed)+"\t"+Integer.toHexString(Float.floatToRawIntBits(observed))+
			"\t"+sourceMatch+"\t"+correctedMatch+"\t"+unchanged+"\n")
			.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
		check(finite, row.name+": expected finite nonnegative error sum, observed="+observed);
		check(sourceMatch==(!row.changed || baseline), row.name+": source-oracle match="+sourceMatch+", observed="+observed);
		check(correctedMatch==(!row.changed || !baseline), row.name+": corrected-oracle match="+correctedMatch+", observed="+observed);
		check(unchanged, row.name+": helper replaced or changed base/quality arrays");
	}

	/** Keeps all six observed rows even when an earlier value fails its oracle. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}

	/** Rejects fixture misuse before helper execution. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Literal indices use A:B pairs; they are descriptive expected coordinates, not observations. */
	private static final class Row{
		Row(final String name_, final String aBases_, final String bBases_, final byte[] aq_, final byte[] bq_,
				final int insert_, final String sourcePairs_, final String correctedPairs_,
				final double sourceSum_, final double correctedSum_, final boolean changed_){
			name=name_; aBases=aBases_; bBases=bBases_; aQualities=aq_; bQualities=bq_; insert=insert_;
			sourcePairs=sourcePairs_; correctedPairs=correctedPairs_;
			sourceSum=sourceSum_; correctedSum=correctedSum_; changed=changed_;
			require(aBases.length()==aQualities.length && bBases.length()==bQualities.length,
				"Each literal base must have its own quality: "+name);
		}
		final String name, aBases, bBases, sourcePairs, correctedPairs;
		final byte[] aQualities, bQualities;
		final int insert;
		final double sourceSum, correctedSum;
		final boolean changed;
	}

	private final Path out;
	private final boolean baseline;
	private int cases, calls, changed, controls, sourceMatches, correctedMatches, checks, failures;
	private static final double TOLERANCE=0.000001;
	private static final byte[] Q4={10, 20, 30, 40}, Q6={10, 20, 30, 40, 10, 20};

	/* On the fragment axis, intersect A=[0,alen), B=[insert-blen,insert), fragment=[0,insert).
	 * Map a shared x to A index x and B index x-(insert-blen); pairs below are explicit.
	 * The helper's probCorrect4 table gives .9,.99,.999,.9999 at Q10,Q20,Q30,Q40.
	 * Decimal sums are sums of 1-C(qa)*C(qb), not phred reconstructions or helper outputs.
	 * No N bases/null qualities/overlap search/model decisions are covered here. */
	private static final Row[] ROWS={
		new Row("equal-short-insert", "AAAA", "AAAA", Q4, Q4, 3,
			"0:1,1:2,2:3", "0:1,1:2,2:3", .1210899, .1210899, false),
		new Row("equal-long-insert", "AAAA", "AAAA", Q4, Q4, 6,
			"2:0,3:1", "2:0,3:1", .110999, .110999, false),
		new Row("short-a-short-insert", "AAAA", "AAAAAA", Q4, Q6, 4,
			"0:0,1:1,2:2,3:3", "0:2,1:3,2:4,3:5", .21209899, .221998, true),
		new Row("short-b-short-insert", "AAAAAA", "AAAA", Q6, Q4, 4,
			"0:2,1:3", "0:0,1:1,2:2,3:3", .110999, .21209899, true),
		new Row("short-a-long-insert", "AAAA", "AAAAAA", Q4, Q6, 7,
			"1:0,2:1,3:2", "1:0,2:1,3:2", .1210899, .1210899, false),
		new Row("short-b-long-insert", "AAAAAA", "AAAA", Q6, Q4, 7,
			"3:0,4:1,5:2", "3:0,4:1,5:2", .22008, .22008, false)
	};
}
