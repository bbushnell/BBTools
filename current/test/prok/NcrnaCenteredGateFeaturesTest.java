package test.prok;

import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaCenteredGateFeatures;
import prok.NcrnaEndpointGateFeatures;
import shared.Shared;

/** Literal whole-move inputs and targets; no learned payload, fit or I/O.
 * @author Brian, Ganyu
 */
public final class NcrnaCenteredGateFeaturesTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCenteredGateFeaturesTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown feature-test option");
			}
			check(parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null,
				"Pure feature controls consume no files; results go to the captured test log");
			run();
			pp.outstream.println("NCRNA_CENTERED_GATE_35_FEATURES_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(){
		final float[] values=new float[69];
		values[63]=1;
		values[64]=.4f;
		values[65]=.8f;
		values[66]=.1f;
		values[67]=.2f;
		values[68]=.3f;
		final float[] scores=new float[22];
		final float[] scratch=new float[32];
		final float[] result=new float[35];
		final float[] legacy=new float[32];
		scores[12]=1;
		long proposal=NcrnaCenteredGateFeatures.fill(scores, values, .45f, .5f, 110,
			100, 209, 100, 209, 0, scratch, result);
		NcrnaEndpointGateFeatures.fill(scores, values, .45f, .5f, 108, 110, legacy);
		for(int i=0; i<32; i++){
			check(Float.floatToRawIntBits(legacy[i])==Float.floatToRawIntBits(result[i]),
				"First32fields retain exact legacy arithmetic when raw and center coincide");
		}
		check(proposal==102 && result[32]==0 && result[33]==0 && result[34]==1,
			"An uncentered neural move has zero applied median and an explicit movement flag");
		Arrays.fill(scores, 0);
		scores[10]=1;
		scores[21]=2;
		values[63]=108f/110;
		proposal=NcrnaCenteredGateFeatures.fill(scores, values, .45f, .5f, 110,
			100, 209, 105, 212, 0, scratch, result);
		check(proposal==105 && result[32]==5 && result[33]==5 && result[34]==0,
			"A current-best neural output still leaves a real offset-only whole-move candidate");
		scores[9]=2;
		scores[11]=2;
		check(NcrnaCenteredGateFeatures.adjustment(scores)==0, "A reachable tie contributes no neural displacement");
		Arrays.fill(scores, 0);
		scores[5]=1;
		proposal=NcrnaCenteredGateFeatures.fill(scores, values, .45f, .5f, 110,
			100, 209, 105, 212, 0, scratch, result);
		check(proposal==100 && result[32]==5 && result[33]==5 && result[34]==1,
			"A neural reversal can cancel the median completely; candidate-equals-raw is a deterministic hold");
		Arrays.fill(scores, 0);
		scores[10]=1;
		values[63]=112f/110;
		proposal=NcrnaCenteredGateFeatures.fill(scores, values, .45f, .5f, 110,
			100, 209, 95, 206, 1, scratch, result);
		check(proposal==206 && result[32]==3 && result[33]==-3 && result[34]==0,
			"Signed applied stop displacement retains gene-sense direction separately from its magnitude");
		scores[12]=2;
		scores[21]=1000;
		check(NcrnaCenteredGateFeatures.adjustment(scores)==2, "The out-of-range head is evidence, not a mandatory proposal veto");
		check(!NcrnaCenteredGateFeatures.improves(100, 104, 100)
			&& NcrnaCenteredGateFeatures.improves(100, 110, 115)
			&& !NcrnaCenteredGateFeatures.improves(100, 110, 105),
			"Targets compare the entire proposal with raw: already-correct harm is negative, strict improvement positive, equal distance negative");
		scores[0]=Float.NaN;
		boolean rejected=false;
		try{NcrnaCenteredGateFeatures.adjustment(scores);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Nonfinite first-stage evidence must fail, never become a fabricated gate feature");
		scores[0]=0;
		values[63]=1;
		rejected=false;
		try{NcrnaCenteredGateFeatures.fill(scores, values, .45f, .5f, 110,
			100, 209, 105, 212, 0, scratch, result);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A raw-frame length feature must not silently enter a centered-frame gate");
	}

	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
