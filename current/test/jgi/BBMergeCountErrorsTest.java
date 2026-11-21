package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;

import jgi.BBMerge;
import parse.Parse;
import stream.Read;

/** Literal raw-orientation error-count cases for JGI090.
 * Both public raw overloads see fresh arrays for each call. Expected coordinates and
 * error counts are hand-enumerated; no tested overlap formula generates an oracle.
 * Five unequal-length cases distinguish wrong counts, false zeroes and a skipped
 * overlap. Ten controls cover equal lengths, short inserts, nonoverlap and null bases.
 * Read objects are nonnull, including the null-base controls. Null Read objects and
 * the joined-Read overload are outside this contract. No merge search or throughput
 * is measured, and known internal callers do not prove these public APIs unused.
 * @author Sucrose
 */
public final class BBMergeCountErrorsTest{

	/** Requires phase=baseline|fixed and a new evidence directory. */
	public static void main(final String[] args) throws IOException{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected documented flag=value: "+arg);
			if(pair[0].equals("phase")){require(phase==null, "Duplicate phase."); phase=pair[1];}
			else if(pair[0].equals("outdir")){require(directory==null, "Duplicate outdir."); directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown raw-error fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergecounterrorstest.sh with assertions enabled.");
		final BBMergeCountErrorsTest test=new BBMergeCountErrorsTest(Paths.get(directory).toAbsolutePath().normalize(), phase.equals("baseline"));
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==15, "Expected fifteen literal input rows: "+test.cases);
		test.check(test.calls==30, "Both public raw overloads must run for every row: "+test.calls);
		test.check(test.changed==5, "Expected five affected unequal-length rows: "+test.changed);
		test.check(test.controls==10, "Expected ten unchanged controls: "+test.controls);
		test.check(test.sourceMatches==(test.baseline ? 30 : 20), "Unexpected source-oracle call matches: "+test.sourceMatches);
		test.check(test.correctedMatches==(test.baseline ? 20 : 30), "Unexpected corrected-oracle call matches: "+test.correctedMatches);
		final String summary="BBMERGE_COUNT_ERRORS_"+(test.failures==0 ? "PASS" : "FAIL")+" phase="+phase+
			" cases="+test.cases+" calls="+test.calls+" changed="+test.changed+" controls="+test.controls+
			" sourceMatches="+test.sourceMatches+" correctedMatches="+test.correctedMatches+" checks="+test.checks+" failures="+test.failures+"\n";
		write(test.out.resolve("summary.txt"), summary); System.out.print(summary);
		if(test.failures>0){throw new AssertionError(summary+"Inspect saved literal inputs, pair oracles and actual calls; do not fit expected counts to output.");}
	}

	private BBMergeCountErrorsTest(final Path out_, final boolean baseline_) throws IOException{
		out=out_; baseline=baseline_;
		Files.createDirectory(out);
		write(out.resolve("inputs.tsv"), "case\tinsert\taBases\tbRawBases\taQualities\tbQualities\tsourcePairsOracle\tcorrectedPairsOracle\tsourceErrorsOracle\tcorrectedErrorsOracle\n");
		write(out.resolve("cases.tsv"), "case\tapi\tobservedErrors\tsourceMatches\tcorrectedMatches\tinputsUnchanged\n");
	}

	/** Retains one literal input row; each API invocation constructs independent arrays. */
	private void run(final Row row) throws IOException{
		cases++;
		if(row.sourceErrors!=row.correctedErrors){changed++;}else{controls++;}
		append(out.resolve("inputs.tsv"), row.name+"\t"+row.insert+"\t"+text(row.a)+"\t"+text(row.b)+
			"\t"+Arrays.toString(qualities(row.a))+"\t"+Arrays.toString(qualities(row.b))+"\t"+row.sourcePairs+"\t"+row.correctedPairs+
			"\t"+row.sourceErrors+"\t"+row.correctedErrors+"\n");
		call(row, false);
		call(row, true);
	}

	/** Calls only the chosen real public raw API and checks arrays after that call. */
	private void call(final Row row, final boolean readApi) throws IOException{
		final byte[] ab=bytes(row.a), bb=bytes(row.b), aq=qualities(row.a), bq=qualities(row.b);
		final byte[] abBefore=copy(ab), bbBefore=copy(bb), aqBefore=copy(aq), bqBefore=copy(bq);
		// The retaining constructor avoids changing literal inputs during normalization.
		final Read a=new Read(ab, aq, row.name+"/A", 0, false), b=new Read(bb, bq, row.name+"/B-raw", 1, false);
		final int observed=readApi ? BBMerge.countErrors(a, b, row.insert) : BBMerge.countErrors(ab, bb, row.insert);
		calls++;
		final boolean sourceMatch=observed==row.sourceErrors, correctedMatch=observed==row.correctedErrors;
		final boolean unchanged=a.bases==ab && b.bases==bb && a.quality==aq && b.quality==bq &&
			Arrays.equals(abBefore, ab) && Arrays.equals(bbBefore, bb) && Arrays.equals(aqBefore, aq) && Arrays.equals(bqBefore, bq);
		if(sourceMatch){sourceMatches++;} if(correctedMatch){correctedMatches++;}
		final String api=readApi ? "Read,Read,int" : "byte[],byte[],int";
		append(out.resolve("cases.tsv"), row.name+"\t"+api+"\t"+observed+"\t"+sourceMatch+"\t"+correctedMatch+"\t"+unchanged+"\n");
		check(sourceMatch==(row.sourceErrors==row.correctedErrors || baseline), row.name+"/"+api+": source-oracle mismatch; observed="+observed);
		check(correctedMatch==(row.sourceErrors==row.correctedErrors || !baseline), row.name+"/"+api+": corrected-oracle mismatch; observed="+observed);
		check(unchanged, row.name+"/"+api+": raw error counting must not replace/reorient/change any base or quality array");
	}

	/** Deliberately stores expected index pairs and integer counts as literals. */
	private static final class Row{
		Row(final String name_, final String a_, final String b_, final int insert_, final String sourcePairs_,
				final String correctedPairs_, final int sourceErrors_, final int correctedErrors_){
			name=name_; a=a_; b=b_; insert=insert_; sourcePairs=sourcePairs_; correctedPairs=correctedPairs_;
			sourceErrors=sourceErrors_; correctedErrors=correctedErrors_;
			require(sourceErrors>=0 && correctedErrors>=0, "Literal error counts must be nonnegative: "+name);
		}
		final String name, a, b, sourcePairs, correctedPairs;
		final int insert, sourceErrors, correctedErrors;
	}

	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private static String text(final String s){return s==null ? "NULL" : s;}
	private static byte[] bytes(final String s){return s==null ? null : s.getBytes(StandardCharsets.US_ASCII);}
	private static byte[] copy(final byte[] a){return a==null ? null : a.clone();}
	/** Distinct quality bytes reveal accidental reversal, although neither API uses them. */
	private static byte[] qualities(final String s){if(s==null){return null;} final byte[] q=new byte[s.length()]; for(int i=0; i<q.length; i++){q[i]=(byte)(10+i);} return q;}
	private static void write(final Path path, final String s) throws IOException{Files.write(path, s.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String s) throws IOException{Files.write(path, s.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);}

	private final Path out;
	private final boolean baseline;
	private int cases, calls, changed, controls, sourceMatches, correctedMatches, checks, failures;

	/* Fragment coordinates: A[x] lies at x; raw B[j] lies at insert-1-j.
	 * A shared coordinate pairs x with j=insert-1-x and compares A to complement(B).
	 * The following pairs/counts were hand-enumerated, not produced with either
	 * source/corrected overlap-start expression. No joined sequence is supplied. */
	private static final Row[] ROWS={
		new Row("short-a-long-insert", "ACGT", "CGTACA", 7, "3:5", "1:5,2:4,3:3", 0, 1),
		new Row("short-b-long-insert", "ACGTAC", "CGTA", 7, "1:3,2:2,3:1,4:0", "3:3,4:2,5:1", 4, 0),
		new Row("short-a-middle-insert", "ACGT", "CACGTG", 5, "1:4,2:3,3:2", "0:4,1:3,2:2,3:1", 3, 0),
		new Row("short-b-middle-insert", "ACGTAC", "TACG", 5, "0:3,1:2,2:1,3:0", "1:3,2:2,3:1,4:0", 4, 0),
		new Row("skipped-a-overlap", "ACGT", "ACGTAAAA", 8, "none", "0:7,1:6,2:5,3:4", 0, 3),
		new Row("unequal-short-insert-control", "ACGT", "ACGTTT", 3, "0:2,1:1,2:0", "0:2,1:1,2:0", 3, 3),
		new Row("equal-short-one-error", "ACGT", "CGAA", 3, "0:2,1:1,2:0", "0:2,1:1,2:0", 1, 1),
		new Row("equal-long-perfect", "ACGT", "AAAC", 6, "2:3,3:2", "2:3,3:2", 0, 0),
		new Row("nonoverlap-boundary", "ACGT", "ACGTTT", 10, "none", "none", 0, 0),
		new Row("nonoverlap-gap", "ACGT", "ACGTTT", 11, "none", "none", 0, 0),
		new Row("zero-insert", "ACGT", "ACGTTT", 0, "none", "none", 0, 0),
		new Row("negative-insert", "ACGT", "ACGTTT", -1, "none", "none", 0, 0),
		new Row("null-a-bases", null, "ACGT", 3, "none", "none", 0, 0),
		new Row("null-b-bases", "ACGT", null, 3, "none", "none", 0, 0),
		new Row("null-both-bases", null, null, 3, "none", "none", 0, 0)
	};
}
