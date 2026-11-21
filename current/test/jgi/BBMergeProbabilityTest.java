package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jgi.BBMergeOverlapper;
import parse.Parse;
import stream.Read;

/** Finite numerical checks for the public overlap likelihood ratio.
 * Mates are already aligned; no overlap search, model inference or timing occurs.
 * Oracles cancel common factors analytically, using pinned rounded quality-table
 * values and explicit counts of adverse match/mismatch groups. They never accumulate
 * the legacy float products. The original twelve Q4 cases retain their input arrays.
 * Additional cases cover the normal/subnormal boundary, numerator-only underflow,
 * legitimate final float zero, Q0's adverse matches, skipped N and one literal
 * unequal-length covered range. That control does not establish general geometry.
 * Baseline records current discrepancies. Fixed requires every oracle and compares
 * ordinary-region raw bits to a successful expanded baseline with identical inputs.
 * Neither mode establishes caller reachability or biological/model accuracy.
 * @author Shinobu, Sucrose
 */
public final class BBMergeProbabilityTest{

	/** Requires phase=baseline|fixed, a new outdir, and reference=baseline/cases for fixed. */
	public static void main(final String[] args) throws IOException{
		String phase=null, directory=null, reference=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){require(phase==null, "Duplicate phase"); phase=pair[1];}
			else if(pair[0].equals("outdir")){require(directory==null, "Duplicate outdir"); directory=pair[1];}
			else if(pair[0].equals("reference")){require(reference==null, "Duplicate reference"); reference=pair[1];}
			else{throw new IllegalArgumentException("Unknown probability fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline or phase=fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		final boolean fixed="fixed".equals(phase);
		require(fixed ? reference!=null && !reference.isEmpty() : reference==null,
			"Only phase=fixed requires reference=PATH to the expanded baseline's fixture output directory.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergeprobabilitytest.sh with assertions enabled.");
		final BBMergeProbabilityTest test=new BBMergeProbabilityTest(fixed,
			Paths.get(directory).toAbsolutePath().normalize());
		test.process(reference==null ? null : Paths.get(reference).toAbsolutePath().normalize());
	}

	private BBMergeProbabilityTest(final boolean fixed_, final Path out_){fixed=fixed_; out=out_;}

	/** Builds all inputs before a call; saves the same manifest in either phase. */
	private void process(final Path reference) throws IOException{
		final ArrayList<Case> inputs=makeCases();
		final byte[] manifest=manifest(inputs).getBytes(StandardCharsets.US_ASCII);
		final List<String> baseline=fixed ? readReference(reference, manifest) : null;
		Files.createDirectory(out);
		Files.write(out.resolve("inputs.tsv"), manifest);
		rows.append(HEADER).append('\n');
		final float[] observed=new float[inputs.size()];
		for(int i=0; i<inputs.size(); i++){
			observed[i]=run(inputs.get(i), baseline==null ? null : baseline.get(i+1));
		}
		check(cases==23 && normalCases==7 && earlyCases==2 && subnormalCases==14,
			"Expected23 cases:7 normal,2 absent-quality and14 subnormal; got "+cases+","+normalCases+","+earlyCases+","+subnormalCases);
		check(originalDiscrepancies==(fixed ? 0 : 4), "Original twelve cases must retain baseline's four discrepancies or have none after repair: "+originalDiscrepancies);
		check(fixed ? close(observed[4], inputs.get(4).expected) && close(observed[5], inputs.get(5).expected) : Math.abs(observed[4]-observed[5])>0.1f,
			"Length256 first/last mismatch: expected baseline position dependence or both fixed oracle matches.");
		check(fixed ? close(observed[6], inputs.get(6).expected) && close(observed[7], inputs.get(7).expected) : Math.abs(observed[6]-observed[7])>0.1f,
			"Length1024 first/last mismatch: expected baseline position dependence or both fixed oracle matches.");
		check(bitComparisons==(fixed ? 9 : 0), "Fixed must compare all seven normal and two early-return scores to saved bits: "+bitComparisons);
		Files.write(out.resolve("cases.tsv"), rows.toString().getBytes(StandardCharsets.US_ASCII));
		final String summary="BBMERGE_PROBABILITY_"+(failures==0 ? "PASS" : "FAIL")+
			" phase="+(fixed ? "fixed" : "baseline")+" version=2 cases="+cases+
			" normalCases="+normalCases+" earlyCases="+earlyCases+" subnormalCases="+subnormalCases+
			" discrepancies="+discrepancies+" originalDiscrepancies="+originalDiscrepancies+
			" bitComparisons="+bitComparisons+" checks="+checks+" failures="+failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.US_ASCII));
		System.out.print(summary);
		if(failures!=0){throw new AssertionError(summary+"Inspect actual values and the oracle; do not fit expectations to output.");}
	}

	/** Keeps the original twelve cases and adds eleven explicit numerical controls. */
	private static ArrayList<Case> makeCases(){
		final ArrayList<Case> list=new ArrayList<Case>();
		list.add(original("short-perfect", 16, -1, 16, false, false, false, false, NORMAL));
		list.add(original("short-first", 16, 0, 16, false, false, false, false, NORMAL));
		list.add(original("short-last", 16, 15, 16, false, false, false, false, NORMAL));
		list.add(original("long-perfect", 256, -1, 256, false, false, false, false, SUBNORMAL));
		list.add(original("long-first", 256, 0, 256, false, false, false, true, SUBNORMAL));
		list.add(original("long-last", 256, 255, 256, false, false, false, true, SUBNORMAL));
		list.add(original("longer-first", 1024, 0, 1024, false, false, false, true, SUBNORMAL));
		list.add(original("longer-last", 1024, 1023, 1024, false, false, false, true, SUBNORMAL));
		list.add(original("all-N", 256, -1, 256, true, false, false, false, NORMAL));
		list.add(original("missing-a-quality", 256, 255, 256, false, true, false, false, EARLY));
		list.add(original("missing-b-quality", 256, 255, 256, false, false, true, false, EARLY));
		list.add(original("nonoverlap", 256, 255, 512, false, false, false, false, NORMAL));
		// At Q0, M=1/4 and E=3/4. All matches are adverse, NOT neutral.
		// (1/4)^63 is exactly Float.MIN_NORMAL; length64 is immediately below.
		list.add(uniform("q0-normal-boundary-63", 63, 0, 0, Math.pow(1.0/3, 31.5), NORMAL, "63 Q0 matches"));
		list.add(uniform("q0-subnormal-boundary-64", 64, 0, 0, Math.pow(1.0/3, 32), SUBNORMAL, "64 Q0 matches"));
		list.add(q0AdverseMatch("q0-long-first-match", 0));
		list.add(q0AdverseMatch("q0-long-last-match", 1023));
		// Taking half the exponent directly avoids underflow of a double product
		// whose square root is still representable as a double (100 mismatches).
		list.add(q40("q40-normal-8-mismatches", 16, 8, NORMAL));
		list.add(q40("q40-subnormal-10-mismatches", 16, 10, SUBNORMAL));
		list.add(q40("q40-zero-numerator-12-mismatches", 16, 12, SUBNORMAL));
		list.add(q40("q40-float-zero-24-mismatches", 32, 24, SUBNORMAL));
		list.add(q40("q40-double-product-zero-100-mismatches", 128, 100, SUBNORMAL));
		list.add(mixed());
		list.add(asymmetric());
		return list;
	}

	/** Reconstructs the original arrays exactly, including a possible C in all-N B. */
	private static Case original(final String name, final int length, final int mismatch, final int insert,
			final boolean allN, final boolean missingA, final boolean missingB, final boolean discrepancy, final int region){
		require(length>0 && mismatch<length, "Original fixture positions must fit their arrays.");
		final byte[] ab=bases(length, allN ? 'N' : 'A'), bb=bases(length, allN ? 'N' : 'A');
		if(mismatch>=0){bb[mismatch]='C';}
		final boolean neutral=mismatch<0 || allN || missingA || missingB || insert>=2*length;
		return new Case(name, ab, bb, missingA ? null : quality(length, 4), missingB ? null : quality(length, 4),
			insert, neutral ? 1 : Math.sqrt(ratio(4, 4, false)), region, true, discrepancy,
			"Q4;mismatch="+mismatch+";allN="+allN, neutral ? "no adverse covered positions" : "1 Q4 mismatch");
	}

	/** Uniform arrays with the first mismatchCount positions of B changed to C. */
	private static Case uniform(final String name, final int length, final int q, final int mismatchCount,
			final double expected, final int region, final String groups){
		require(mismatchCount>=0 && mismatchCount<=length, "Uniform mismatch count exceeds array length.");
		final byte[] bb=bases(length, 'A');
		Arrays.fill(bb, 0, mismatchCount, (byte)'C');
		return new Case(name, bases(length, 'A'), bb, quality(length, q), quality(length, q), length,
			expected, region, false, false, "Q"+q+";prefixMismatches="+mismatchCount, groups);
	}

	private static Case q40(final String name, final int length, final int mismatches, final int region){
		return uniform(name, length, 40, mismatches, Math.pow(ratio(40, 40, false), mismatches*0.5),
			region, mismatches+" Q40 mismatches; remaining matches cancel");
	}

	/** All 1023 mismatches cancel at Q0; only the one match penalizes the ratio. */
	private static Case q0AdverseMatch(final String name, final int matchIndex){
		require(matchIndex==0 || matchIndex==1023, "Q0 position control must select an endpoint.");
		final byte[] bb=bases(1024, 'C'); bb[matchIndex]='A';
		return new Case(name, bases(1024, 'A'), bb, quality(1024, 0), quality(1024, 0), 1024,
			Math.sqrt(1.0/3), SUBNORMAL, false, false, "Q0;onlyMatch="+matchIndex,
			"1 Q0 match;1023 neutral Q0 mismatches");
	}

	/** Four adverse groups with separate pinned qualities; two N positions are skipped. */
	private static Case mixed(){
		final byte[] ab=bases(512, 'A'), bb=bases(512, 'A'), aq=quality(512, 4), bq=quality(512, 4);
		bb[0]='C'; aq[0]=bq[0]=40;
		bb[255]='C';
		aq[511]=bq[511]=0;
		bb[400]='C'; bq[400]=40;
		bb[100]='N'; aq[100]=bq[100]=40;
		ab[200]='N'; bb[200]='C'; aq[200]=bq[200]=40;
		final double expected=Math.sqrt(ratio(40, 40, false)*ratio(4, 4, false)*ratio(0, 0, true)*ratio(4, 40, false));
		return new Case("mixed-quality-N", ab, bb, aq, bq, 512, expected, SUBNORMAL, false, false,
			"Q40 mismatch0;Q4 mismatch255;Q0 match511;Q4/40 mismatch400;N100,200",
			"1 Q40 mismatch;1 Q4 mismatch;1 Q0 match;1 Q4/40 mismatch;2 skipped N");
	}

	/** Literal geometry: A[0..255] pairs with B[44..299]; B[0] is outside. */
	private static Case asymmetric(){
		final byte[] bb=bases(300, 'A'); bb[0]='C'; bb[44]='C'; bb[100]='N';
		return new Case("asymmetric-covered-range", bases(256, 'A'), bb, quality(256, 4), quality(300, 4),
			256, Math.sqrt(ratio(4, 4, false)), SUBNORMAL, false, false,
			"A256;B300;insert256;B0C excluded;B44C covered;B100N skipped",
			"1 covered Q4 mismatch;254 covered Q4 matches cancel;1 covered N skipped");
	}

	/** Calls the real helper; source-classified normal/early cases use saved baseline bits. */
	private float run(final Case c, final String baselineRow){
		final byte[] abBefore=c.ab.clone(), bbBefore=c.bb.clone();
		final byte[] aqBefore=c.aq==null ? null : c.aq.clone(), bqBefore=c.bq==null ? null : c.bq.clone();
		final Read a=new Read(c.ab, c.aq, c.name+"/A", 0, false), b=new Read(c.bb, c.bq, c.name+"/B-aligned", 1, false);
		final float observed=BBMergeOverlapper.probability(a, b, c.insert);
		final boolean matches=close(observed, c.expected);
		final boolean unchanged=a.bases==c.ab && b.bases==c.bb && a.quality==c.aq && b.quality==c.bq &&
			Arrays.equals(c.ab, abBefore) && Arrays.equals(c.bb, bbBefore) && Arrays.equals(c.aq, aqBefore) && Arrays.equals(c.bq, bqBefore);
		cases++;
		if(c.region==NORMAL){normalCases++;}else if(c.region==EARLY){earlyCases++;}else{subnormalCases++;}
		if(!matches){discrepancies++; if(c.original){originalDiscrepancies++;}}
		check(Float.isFinite(observed) && observed>=0 && observed<=1, c.name+": expected finite score in[0,1], got "+observed);
		check(unchanged, c.name+": probability must not replace or mutate input arrays.");
		// New subnormal-region baseline outcomes are recorded, not fitted in advance.
		// Legitimate final-zero controls must already return +0 in either phase.
		if(fixed || c.original || c.region==NORMAL || (float)c.expected==0){
			check(!fixed && c.original && c.baselineDiscrepancy ? !matches : matches,
				c.name+": unexpected oracle match="+matches+", observed="+observed+", oracle="+c.expected);
		}
		String referenceMatch="NA";
		if(fixed && c.region!=SUBNORMAL){
			final String[] fields=baselineRow.split("\t", -1);
			final boolean same=bits(observed).equals(fields[9]);
			bitComparisons++; referenceMatch=Boolean.toString(same);
			check(same, c.name+": ordinary/early score bits differ from saved expanded baseline: "+fields[9]+" versus "+bits(observed));
		}
		rows.append(c.name).append('\t').append(c.ab.length).append('\t').append(c.bb.length).append('\t').append(c.insert)
			.append('\t').append(c.pattern).append('\t').append(c.groups).append('\t').append(c.expected)
			.append('\t').append(bits((float)c.expected)).append('\t').append(observed).append('\t').append(bits(observed))
			.append('\t').append(Math.abs(observed-c.expected)).append('\t').append(tolerance(c.expected))
			.append('\t').append(matches).append('\t').append(unchanged).append('\t').append(REGIONS[c.region])
			.append('\t').append(referenceMatch).append('\n');
		return observed;
	}

	/** Rejects old/incomplete baselines and mismatched inputs before creating output. */
	private List<String> readReference(final Path reference, final byte[] manifest) throws IOException{
		final String summary=new String(Files.readAllBytes(reference.resolve("summary.txt")), StandardCharsets.US_ASCII);
		require(summary.matches("BBMERGE_PROBABILITY_PASS phase=baseline version=2 cases=23 normalCases=7 earlyCases=2 subnormalCases=14 discrepancies=[0-9]+ originalDiscrepancies=4 bitComparisons=0 checks=67 failures=0\\n"),
			"Reference must be a successful version2 expanded baseline fixture, not the old twelve-case run: "+reference);
		require(Arrays.equals(manifest, Files.readAllBytes(reference.resolve("inputs.tsv"))),
			"Fixed/reference inputs and grouped oracles must be byte-identical.");
		final List<String> lines=Files.readAllLines(reference.resolve("cases.tsv"), StandardCharsets.US_ASCII);
		require(lines.size()==24 && HEADER.equals(lines.get(0)), "Expanded reference requires the version2 header and23 rows.");
		for(int i=1; i<lines.size(); i++){
			final String[] fields=lines.get(i).split("\t", -1);
			require(fields.length==16 && fields[9].matches("[0-9a-f]{8}") && "true".equals(fields[13]),
				"Malformed baseline row or changed baseline input arrays at row "+i);
		}
		// Matching the complete case manifest also fixes row order and expected values;
		// require row names explicitly so a reordered/truncated result cannot supply bits.
		final String[] inputLines=new String(manifest, StandardCharsets.US_ASCII).split("\n");
		for(int i=1; i<lines.size(); i++){
			require(lines.get(i).substring(0, lines.get(i).indexOf('\t')).equals(inputLines[i].substring(0, inputLines[i].indexOf('\t'))),
				"Reference result/input case order differs at row "+i);
		}
		return lines;
	}

	/** Includes literal bases, numeric qualities, coverage, groups and source region. */
	private static String manifest(final List<Case> inputs){
		final StringBuilder text=new StringBuilder("case\tinsert\tA\tB-aligned\tA-quality\tB-quality\toracle\tregion\tgroups\tpattern\n");
		for(final Case c : inputs){
			text.append(c.name).append('\t').append(c.insert).append('\t').append(new String(c.ab, StandardCharsets.US_ASCII))
				.append('\t').append(new String(c.bb, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(c.aq))
				.append('\t').append(Arrays.toString(c.bq)).append('\t').append(c.expected).append('\t').append(REGIONS[c.region])
				.append('\t').append(c.groups).append('\t').append(c.pattern).append('\n');
		}
		return text.toString();
	}

	/** Single-pair float rounding is contractual; no accumulated-product oracle. */
	private static double ratio(final int qa, final int qb, final boolean matched){
		final float correct=correct(qa)*correct(qb), match=correct+(1-correct)*0.25f, error=1-match;
		return (double)(matched ? match : error)/Math.max(match, error);
	}
	private static float correct(final int q){
		if(q==0){return 0;}if(q==4){return 0.6019f;}if(q==40){return 0.9999f;}
		throw new IllegalArgumentException("Unpinned fixture quality: "+q);
	}
	/** Nonzero expected floats cannot pass with zero, regardless of the ULP allowance. */
	private static boolean close(final float observed, final double expected){
		if((float)expected==0){return Float.floatToRawIntBits(observed)==0;}
		return observed>0 && Math.abs(observed-expected)<=tolerance(expected);
	}
	private static double tolerance(final double expected){return Math.max(8.0*Math.ulp((float)expected), 0.000002*Math.abs(expected));}
	private static byte[] bases(final int length, final char base){final byte[] b=new byte[length]; Arrays.fill(b, (byte)base); return b;}
	private static byte[] quality(final int length, final int q){final byte[] b=new byte[length]; Arrays.fill(b, (byte)q); return b;}
	private static String bits(final float value){return String.format(java.util.Locale.ROOT, "%08x", Float.floatToRawIntBits(value));}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}}

	/** Owned finite inputs plus an independently specified grouped oracle. */
	private static final class Case{
		Case(final String name_, final byte[] ab_, final byte[] bb_, final byte[] aq_, final byte[] bq_, final int insert_,
				final double expected_, final int region_, final boolean original_, final boolean discrepancy_, final String pattern_, final String groups_){
			require(expected_>=0 && expected_<=1 && Double.isFinite(expected_), "Grouped oracle must be a finite score: "+name_);
			name=name_; ab=ab_; bb=bb_; aq=aq_; bq=bq_; insert=insert_; expected=expected_; region=region_;
			original=original_; baselineDiscrepancy=discrepancy_; pattern=pattern_; groups=groups_;
		}
		final String name, pattern, groups;
		final byte[] ab, bb, aq, bq;
		final int insert, region;
		final double expected;
		final boolean original, baselineDiscrepancy;
	}

	private final boolean fixed;
	private final Path out;
	private final StringBuilder rows=new StringBuilder();
	private int cases, normalCases, earlyCases, subnormalCases, discrepancies, originalDiscrepancies, bitComparisons, checks, failures;
	// Regions classify the legacy final products from the stated finite inputs;
	// they are not instrumented observations of the production helper's branch.
	private static final int NORMAL=0, EARLY=1, SUBNORMAL=2;
	private static final String[] REGIONS={"normal", "absent-quality", "subnormal"};
	private static final String HEADER="case\taLength\tbLength\tinsert\tpattern\toracleGroups\toracle\texpectedFloatBits\tobserved\tobservedBits\tabsError\ttolerance\toracleMatches\tinputsUnchanged\tlegacyRegion\treferenceBitsMatch";
}
