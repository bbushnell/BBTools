package test.align2;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import align2.MakeQualityHistogram;
import align2.MultiStateAligner9PacBio;
import align2.PacBioScoreParameters;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks the empty quality-bin print policy and PacBio DEFAULT-equivalent
 * imperfect-score helper policy with finite, independently-derived expectations.
 * @author Brian Bushnell, Jean
 */
public final class QualityPacBioPolicyTest{

	/** Accepts explicit mode=baseline|fixed, case=histogram|helper|bench and iterations=N arguments. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, QualityPacBioPolicyTest.class, false);
		try{
			String mode=null, name=null;
			long iterations=-1;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected key=value argument: "+arg);}
				if(pair[0].equalsIgnoreCase("mode")){mode=pair[1];}
				else if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
				else if(pair[0].equalsIgnoreCase("iterations")){iterations=Parse.parseKMG(pair[1]);}
				else if(Parser.parseStatic(arg, pair[0].toLowerCase(), pair[1])){
					// Accept normal javasetup flags such as simd=f before fixture keys.
				}else{throw new IllegalArgumentException("Unexpected argument: "+arg);}
			}
			if(!"baseline".equalsIgnoreCase(mode) && !"fixed".equalsIgnoreCase(mode)){
				throw new IllegalArgumentException("Expected explicit mode=baseline or mode=fixed.");
			}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("Quality/PacBio checks require -ea.");}
			final boolean baseline=mode.equalsIgnoreCase("baseline");
			final Case kase=Case.forName(name);
			if(kase==Case.HISTOGRAM){runHistogramSuite(baseline);}
			else if(kase==Case.HELPER){runHelperSuite(baseline);}
			else if(kase==Case.BENCH){runBench(baseline, iterations);}
			else{throw new IllegalStateException("Unhandled case: "+kase);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Static fixture entry point only. */
	private QualityPacBioPolicyTest(){}

	/** Checks both quality histogram printers on empty, ordinary and large bins. */
	private static void runHistogramSuite(final boolean baseline){
		final String empty=baseline ? "NaN" : "0.000";
		final int large=1<<29;
		final int[][] mapped={{0, 1, 3*large, Integer.MAX_VALUE},
				{0, 3, large, Integer.MAX_VALUE}};
		final int[][] paired={{0, 2, large, Integer.MAX_VALUE},
				{0, 3, 3*large, Integer.MAX_VALUE}};
		final String mappedExpected=
				"#Error Quality Histogram\n"+
				"Quality\tMapped\tUnmapped\tPercent Mapped\n"+
				"0\t0\t0\t"+empty+"\n"+
				"1\t1\t3\t25.000\n"+
				"2\t1610612736\t536870912\t75.000\n"+
				"3\t2147483647\t2147483647\t50.000\n";
		final String pairedExpected=
				"#Error Quality Histogram\n"+
				"Quality\tPaired\tSingle\tPercent Paired\n"+
				"0\t0\t0\t"+empty+"\n"+
				"1\t2\t3\t40.000\n"+
				"2\t536870912\t1610612736\t25.000\n"+
				"3\t2147483647\t2147483647\t50.000\n";
		checkText("mapped", mappedExpected, captureMapped(mapped));
		checkText("paired", pairedExpected, capturePaired(paired));
		System.out.println("QUALITY_PACBIO case=histogram mode="+
				(baseline ? "baseline" : "fixed")+" checks=2 failed=0");
	}

	/** Checks DEFAULT-equivalent helpers, unsupported costs, and experimental alignment availability. */
	private static void runHelperSuite(final boolean baseline){
		final PacBioScoreParameters equivalent=PacBioScoreParameters.experimental(
				MultiStateAligner9PacBio.POINTS_DEL4,
				MultiStateAligner9PacBio.POINTS_DEL5,
				MultiStateAligner9PacBio.POINTS_INS4,
				MultiStateAligner9PacBio.POINTS_INS,
				MultiStateAligner9PacBio.POINTS_DEL, 1f);
		if(equivalent==PacBioScoreParameters.DEFAULT ||
				!equivalent.sameValues(PacBioScoreParameters.DEFAULT)){
			throw new AssertionError("Expected a distinct DEFAULT-equivalent cost object.");
		}
		checkAccepted("default", PacBioScoreParameters.DEFAULT);
		checkAccepted("equivalent", equivalent);
		final PacBioScoreParameters del4=PacBioScoreParameters.experimental(
				MultiStateAligner9PacBio.POINTS_DEL4+1,
				MultiStateAligner9PacBio.POINTS_DEL5,
				MultiStateAligner9PacBio.POINTS_INS4,
				MultiStateAligner9PacBio.POINTS_INS,
				MultiStateAligner9PacBio.POINTS_DEL, 1f);
		final PacBioScoreParameters zero=PacBioScoreParameters.experimental(0, 0, 0, 0, 0, 1f);
		checkCondition(!del4.sameValues(PacBioScoreParameters.DEFAULT),
				"del4-only costs should differ from DEFAULT");
		checkConditional("del4", del4, baseline);
		checkUnsupported("zero", zero, baseline);
		checkExperimentalAlignment(zero);
		System.out.println("QUALITY_PACBIO case=helper mode="+
				(baseline ? "baseline" : "fixed")+" checks=9 failed=0");
	}

	/** Verifies both helper overloads on one allowed cost object. */
	private static void checkAccepted(final String label, final PacBioScoreParameters costs){
		final MultiStateAligner9PacBio msa=new MultiStateAligner9PacBio(4, 4, costs);
		final int intScore=msa.maxImperfectScore(10);
		final int byteScore=msa.maxImperfectScore(new byte[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
		// DEFAULT: maxQuality(10) + min(POINTS_DEL, POINTS_INS - POINTS_MATCH2) = 990 - 305.
		if(intScore!=685 || byteScore!=740){
			throw new AssertionError(label+" expected 685/740 but observed "+
					intScore+"/"+byteScore);
		}
	}

	/** Verifies both helper overloads accept or reject an extension-only non-DEFAULT cost object. */
	private static void checkConditional(final String label, final PacBioScoreParameters costs,
			final boolean baseline){
		if(baseline){checkAccepted(label, costs);}
		else{checkDisallowed(label+"_int", costs, false, IllegalStateException.class);
			checkDisallowed(label+"_bytes", costs, true, IllegalStateException.class);}
	}

	/** Verifies both helper overloads reject zero-gap costs by the expected mechanism. */
	private static void checkUnsupported(final String label, final PacBioScoreParameters costs,
			final boolean baseline){
		final Class<?> expected=baseline ? AssertionError.class : IllegalStateException.class;
		checkDisallowed(label+"_int", costs, false, expected);
		checkDisallowed(label+"_bytes", costs, true, expected);
	}

	/** Verifies one helper overload rejects one unsupported cost object. */
	private static void checkDisallowed(final String label, final PacBioScoreParameters costs,
			final boolean bytes, final Class<?> expected){
		final MultiStateAligner9PacBio msa=new MultiStateAligner9PacBio(4, 4, costs);
		try{
			if(bytes){msa.maxImperfectScore(new byte[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10});}
			else{msa.maxImperfectScore(10);}
		}catch(final Throwable t){
			if(expected.isInstance(t)){
				System.out.println("unsupported\t"+label+"\t"+expected.getSimpleName());
				return;
			}
			throw new AssertionError(label+" expected "+expected.getSimpleName()+
					" but observed "+t, t);
		}
		throw new AssertionError(label+" unexpectedly accepted unsupported costs");
	}

	/** Verifies experimental alignment still fills, scores and traces when the helper is out of domain. */
	private static void checkExperimentalAlignment(final PacBioScoreParameters costs){
		final byte[] read={'A', 'C', 'G', 'T'};
		final byte[] ref={'A', 'C', 'G', 'T'};
		final MultiStateAligner9PacBio msa=new MultiStateAligner9PacBio(read.length, ref.length, costs);
		final int[] max=msa.fillUnlimited(read, ref, 0, ref.length-1, null);
		checkCondition(max!=null && max[0]==4 && max[1]==4 && max[3]==390,
				"Expected perfect fill row=4 col=4 score=390 but observed "+arrayString(max));
		final int[] score=msa.score(read, ref, 0, ref.length-1, max[0], max[1], max[2], false);
		final byte[] match=msa.traceback(read, ref, 0, ref.length-1, max[0], max[1], max[2], false);
		checkCondition(score[0]==390 && score[1]==0 && score[2]==3,
				"Expected score 390,0,3 but observed "+arrayString(score));
		checkCondition("mmmm".equals(new String(match)),
				"Expected perfect match mmmm but observed "+new String(match));
	}

	/** Times default helper calls with reused inputs and a rolling checksum. */
	private static void runBench(final boolean baseline, final long iterations){
		if(iterations<1){throw new IllegalArgumentException("case=bench requires iterations=N with N>0");}
		if(iterations>Long.MAX_VALUE/2){
			throw new IllegalArgumentException("iterations is too large to report helper_calls safely: "+
					iterations);
		}
		final MultiStateAligner9PacBio msa=new MultiStateAligner9PacBio(8, 8, PacBioScoreParameters.DEFAULT);
		final long warmup=Math.min(iterations, 1000000L);
		long checksum=0;
		for(long i=0; i<warmup; i++){checksum=roll(checksum, benchOnce(msa, i));}
		final long start=System.nanoTime();
		for(long i=0; i<iterations; i++){checksum=roll(checksum, benchOnce(msa, i));}
		final long elapsed=System.nanoTime()-start;
		System.out.println("QUALITY_PACBIO case=bench mode="+
				(baseline ? "baseline" : "fixed")+
				" iterations="+iterations+" warmup="+warmup+
				" helper_calls="+(2*iterations)+
				" checksum="+checksum+" elapsed_nanos="+elapsed);
	}

	/** Makes one paired int/byte helper call for the benchmark. */
	private static int benchOnce(final MultiStateAligner9PacBio msa, final long round){
		final int a=msa.maxImperfectScore(INT_LENGTHS[(int)(round&3)]);
		final int b=msa.maxImperfectScore(BYTE_SCORES[(int)((round>>>2)&3)]);
		return a*31+b;
	}

	/** Updates the rolling benchmark checksum. */
	private static long roll(final long checksum, final int value){
		return checksum*1000003L+value;
	}

	/** Captures one mapped histogram. */
	private static String captureMapped(final int[][] counts){
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream old=System.out;
		System.setOut(new PrintStream(buffer));
		try{
			MakeQualityHistogram.printMappedHistogram(counts);
		}finally{
			System.out.flush();
			System.setOut(old);
		}
		return buffer.toString();
	}

	/** Captures one paired histogram. */
	private static String capturePaired(final int[][] counts){
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream old=System.out;
		System.setOut(new PrintStream(buffer));
		try{
			MakeQualityHistogram.printPairedHistogram(counts);
		}finally{
			System.out.flush();
			System.setOut(old);
		}
		return buffer.toString();
	}

	/** Requires exact bounded text. */
	private static void checkText(final String label, final String expected, final String observed){
		if(!expected.equals(observed)){
			throw new AssertionError(label+" expected:\n"+expected+"\nbut observed:\n"+observed);
		}
	}

	/** Keeps test failures loud without relying on the fixture's own -ea gate. */
	private static void checkCondition(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
	}

	/** Formats a nullable int array for diagnostics without importing Arrays only for failure text. */
	private static String arrayString(final int[] array){
		if(array==null){return "null";}
		final StringBuilder sb=new StringBuilder("[");
		for(int i=0; i<array.length; i++){
			if(i>0){sb.append(", ");}
			sb.append(array[i]);
		}
		return sb.append(']').toString();
	}

	/** One top-level fixture case, each run in a separate JVM by the outer recipe. */
	private enum Case{
		HISTOGRAM("histogram"),
		HELPER("helper"),
		BENCH("bench");

		Case(final String name_){name=name_;}

		static Case forName(final String name){
			for(final Case c : values()){
				if(c.name.equalsIgnoreCase(name)){return c;}
			}
			throw new IllegalArgumentException("Unknown case: "+name);
		}

		final String name;
	}

	/** Small helper lengths for the no-I/O timed loop. */
	private static final int[] INT_LENGTHS={5, 6, 7, 8};

	/** Reused byte-score arrays for the no-I/O timed loop. */
	private static final byte[][] BYTE_SCORES={
		{1, 2, 4, 8},
		{2, 3, 5, 7, 11},
		{1, 1, 2, 3, 5, 8},
		{13, 8, 5, 3, 2, 1, 1}
	};
}
