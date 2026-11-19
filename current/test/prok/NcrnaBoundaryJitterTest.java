package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Locale;
import java.util.Random;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaBoundaryJitter;
import prok.NcrnaSignedEndpointErrors;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Hand-counted conditioning and sampling oracle, without biological inputs.
 * Every retained integer weight is visited once using a controlled RNG; this
 * detects edge clamping without a flaky statistical-frequency assertion.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBoundaryJitterTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaBoundaryJitterTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown finite fixture option: "+arg);
			}
			check(parser.out1!=null && new File(parser.out1).mkdir(), "Use one fresh out= directory");
			final File root=new File(parser.out1);
			final String histogram=HEADER
				+"RF00174\t5prime\tsigned\t-20\t10\nRF00174\t5prime\tsigned\t-6\t1\n"
				+"RF00174\t5prime\tsigned\t-1\t2\nRF00174\t5prime\tsigned\t0\t3\n"
				+"RF00174\t5prime\tsigned\t5\t5\nRF00174\t5prime\tsigned\t6\t4\n"
				+"RF00174\t5prime\tsigned\t20\t20\nRF00174\t5prime\tabsolute\t0\t3\n"
				+"RF00174\t3prime\tsigned\t-30\t20\nRF00174\t3prime\tsigned\t-6\t5\n"
				+"RF00174\t3prime\tsigned\t0\t6\nRF00174\t3prime\tsigned\t6\t4\n"
				+"RF00174\t3prime\tsigned\t30\t10\n";
			final String path=write(root, "histogram.tsv", histogram);
			final NcrnaBoundaryJitter jitter=NcrnaBoundaryJitter.load(path, "RF00174", 6);
			check(jitter.observed(0)==45 && jitter.retained(0)==15 && jitter.excludedLow(0)==10 && jitter.excludedHigh(0)==20,
				"Five-prime conditioning retains15of45, excluding10below and20above, with no edge clamping");
			check(jitter.observed(1)==45 && jitter.retained(1)==15 && jitter.excludedLow(1)==20 && jitter.excludedHigh(1)==10,
				"Three-prime counts preserve their own measured tails and matched denominator");
			final int[] five={-6, -1, -1, 0, 0, 0, 5, 5, 5, 5, 5, 6, 6, 6, 6};
			final int[] three={-6, -6, -6, -6, -6, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6};
			final ControlledRandom random=new ControlledRandom();
			for(int i=0; i<15; i++){
				random.value=i;
				check(jitter.draw(0, random)==five[i] && jitter.draw(1, random)==three[i],
					"Every retained integer weight must map to its literal observed signed-error bin: "+i);
			}
			final NcrnaBoundaryJitter narrow=NcrnaBoundaryJitter.load(path, "RF00174", 1);
			check(narrow.retained(0)==5 && narrow.excludedLow(0)==11 && narrow.excludedHigh(0)==29
				&& narrow.retained(1)==6 && narrow.excludedLow(1)==25 && narrow.excludedHigh(1)==14,
				"Changing radius conditions the original histogram, not a previously clipped distribution");
			final String empty=HEADER+"RF00174\t5prime\tsigned\t-20\t2\nRF00174\t5prime\tsigned\t20\t3\n"
				+"RF00174\t3prime\tsigned\t-20\t2\nRF00174\t3prime\tsigned\t20\t3\n";
			final NcrnaBoundaryJitter missing=NcrnaBoundaryJitter.load(write(root, "empty.tsv", empty), "RF00174", 6);
			check(!missing.available(0) && !missing.available(1) && missing.observed(0)==5 && missing.retained(0)==0,
				"No measured mass in the window must remain explicitly unavailable");
			fails(()->missing.draw(0, random), "Cannot draw");
			missing.write(new File(root, "empty_summary.tsv").getPath());
			final String summary=new String(Files.readAllBytes(new File(root, "empty_summary.tsv").toPath()), StandardCharsets.UTF_8);
			check(summary.contains("RF00174\t5prime\t6\t5\t0\t2\t3\tEMPTY_CONDITIONING_WINDOW\n")
				&& summary.contains("RF00174\t3prime\t6\t5\t0\t2\t3\tEMPTY_CONDITIONING_WINDOW\n"),
				"Empty windows must still emit exact tail accounting, not invented vectors");
			fails(()->NcrnaBoundaryJitter.load(path, "RF01482", 6), "selected emitted member");
			final String duplicate=write(root, "duplicate.tsv", histogram+"RF00174\t5prime\tsigned\t-1\t2\n");
			fails(()->NcrnaBoundaryJitter.load(duplicate, "RF00174", 6), "strictly increasing");
			final String unequal=write(root, "unequal.tsv", histogram.replace("3prime\tsigned\t30\t10", "3prime\tsigned\t30\t9"));
			fails(()->NcrnaBoundaryJitter.load(unequal, "RF00174", 6), "matched-locus denominator");
			final String wrong=write(root, "wrong_measure.tsv", histogram.replace("\tsigned\t", "\tmagnitude\t"));
			fails(()->NcrnaBoundaryJitter.load(wrong, "RF00174", 6), "Unknown histogram measure");
			checkFullErrors(root, path, duplicate, unequal, wrong);
			pp.outstream.println("NCRNA_BOUNDARY_JITTER_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Visits all 45 integer masses, including both distant tails and the zero bin. */
	private static void checkFullErrors(final File root, final String path, final String duplicate,
			final String unequal, final String wrong) throws Exception{
		final NcrnaSignedEndpointErrors errors=NcrnaSignedEndpointErrors.load(path, "RF00174");
		final NcrnaSignedEndpointErrors.Distribution five=errors.distribution(0), three=errors.distribution(1);
		check(five.total==45 && five.zero==3 && five.incorrect()==42
			&& three.total==45 && three.zero==6 && three.incorrect()==39,
			"Full signed distributions must retain their separate measured zero fractions and every distant count");
		check(five.backgroundWeight(3)==0 && five.backgroundWeight(0)==10.0/42
			&& three.backgroundWeight(0)==20.0/39,
			"Wrong-end weights exclude zero while retaining each end's signed tail and nonzero denominator");
		final ControlledRandom random=new ControlledRandom();
		for(int mass=0; mass<45; mass++){
			random.value=mass;
			final int expectedFive=mass<10 ? -20 : mass<11 ? -6 : mass<13 ? -1 : mass<16 ? 0 : mass<21 ? 5 : mass<25 ? 6 : 20;
			final int expectedThree=mass<20 ? -30 : mass<25 ? -6 : mass<31 ? 0 : mass<35 ? 6 : 30;
			check(five.draw(random)==expectedFive && three.draw(random)==expectedThree,
				"Every original integer mass must select its measured signed error without conditioning or clamping: "+mass);
		}
		final String bandsPath=new File(root, "full_bands.tsv").getPath();
		errors.writeBands(bandsPath);
		final String bands=new String(Files.readAllBytes(new File(bandsPath).toPath()), StandardCharsets.UTF_8);
		check(bands.contains("RF00174\t5prime\t6\t10\t3\t3\t9\t20\t45\n")
			&& bands.contains("RF00174\t3prime\t6\t20\t5\t6\t4\t10\t45\n")
			&& bands.contains("RF00174\t5prime\t20\t0\t13\t3\t29\t0\t45\n"),
			"Directional radius rows must partition the full signed counts with zero separate");
		final String cliBands=new File(root, "full_bands_cli.tsv").getPath();
		NcrnaSignedEndpointErrors.main(new String[]{"in="+path, "out="+cliBands, "family=RF00174"});
		check(java.util.Arrays.equals(Files.readAllBytes(new File(bandsPath).toPath()), Files.readAllBytes(new File(cliBands).toPath())),
			"The small native report CLI must preserve the independently checked directional table bytes");
		final int[] offsets={-7, 0, 11};
		final long[] counts={2, 3, 1};
		final NcrnaSignedEndpointErrors.Distribution copied=new NcrnaSignedEndpointErrors.Distribution(offsets, counts);
		offsets[0]=999;
		counts[0]=999;
		check(copied.offset(0)==-7 && copied.count(0)==2 && copied.total==6,
			"Immutable distributions must own their primitive support and counts");
		final NcrnaSignedEndpointErrors.Distribution correct=new NcrnaSignedEndpointErrors.Distribution(new int[]{0}, new long[]{5});
		check(correct.draw(random)==0 && correct.zero==5, "An all-correct measured distribution draws only zero");
		fails(()->correct.backgroundWeight(0), "no measured wrong-end background");
		fails(()->new NcrnaSignedEndpointErrors.Distribution(new int[]{0, 0}, new long[]{1, 2}), "unique increasing");
		fails(()->new NcrnaSignedEndpointErrors.Distribution(new int[]{0}, new long[]{0}), "positive counts");
		fails(()->NcrnaSignedEndpointErrors.load(path, "RF01482"), "selected member");
		fails(()->NcrnaSignedEndpointErrors.load(duplicate, "RF00174"), "unique increasing");
		fails(()->NcrnaSignedEndpointErrors.load(unequal, "RF00174"), "same matched-locus population");
		fails(()->NcrnaSignedEndpointErrors.load(wrong, "RF00174"), "Unknown measure");
	}

	private static String write(final File root, final String name, final String contents){
		final String path=new File(root, name).getPath();
		try(Output output=new Output(path)){output.writer.print(contents);}
		return path;
	}
	private static void fails(final Runnable action, final String expected){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(expected), "Wrong failure instead of the expected input limitation: "+failure);
			return;
		}
		throw new AssertionError("Invalid conditioning input was accepted: "+expected);
	}
	private static final class ControlledRandom extends Random{
		@Override
		public long nextLong(){return ((long)value)<<1;}
		int value;
		private static final long serialVersionUID=1L;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String HEADER="family\tend\tmeasure\toffset_bp\tloci\n";
}
