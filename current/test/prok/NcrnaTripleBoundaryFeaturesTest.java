package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import idaligner.AlignmentStats;
import idaligner.ScrabbleAligner;
import map.LongHashSet;
import ml.CellNet;
import parse.PreParser;
import prok.NcrnaBoundaryRefiner;
import prok.NcrnaFamily;
import prok.NcrnaPositionalBoundaryFeatures;
import prok.NcrnaPositionalBoundaryRefiner;
import prok.Orf;
import prok.ProkObject;
import prok.RrnaPositionalKmerTable;
import prok.RrnaPositionalKmerTable.Table;
import prok.TrnaBoundaryFeatures;
import shared.Shared;

/** Optional K9/R10 input-order, omitted-locus train/serve and shared-HBM oracle.
 * Literal rebuilt tables are the LOO oracle; the unchanged single-end sliding
 * helper supplies independent coverage values for the new paired helper.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaTripleBoundaryFeaturesTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaTripleBoundaryFeaturesTest.class, false);
		try{
			if(pp.args.length!=0){throw new IllegalArgumentException("This finite oracle accepts no parameters: "+pp.args[0]);}
			run();
			pp.outstream.println("NCRNA_TRIPLE_BOUNDARY_FEATURES_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(){
		final byte[] consensus="ACGTTGCAAGTCGATCGTACCTGAGACTGCTAGTCAGGCTACTGACGATTCGATGCTAGC".getBytes(StandardCharsets.US_ASCII);
		final int start=80, stop=start+consensus.length-1;
		final byte[] bases=new byte[240];
		Arrays.fill(bases, (byte)'A');
		System.arraycopy(consensus, 0, bases, start, consensus.length);
		final BaseGraph graph=model(consensus);
		final Table[] five=tables(bases, "5prime", -2, start, 2), three=tables(bases, "3prime", -7, stop, 2);
		final NcrnaPositionalBoundaryFeatures ordinary=new NcrnaPositionalBoundaryFeatures(five, three, false);
		final NcrnaPositionalBoundaryFeatures training=new NcrnaPositionalBoundaryFeatures(five, three, true);
		final NcrnaPositionalBoundaryFeatures rebuilt=new NcrnaPositionalBoundaryFeatures(
			tables(bases, "5prime", -2, start, 1), tables(bases, "3prime", -7, stop, 1), false);
		training.setLeaveOneOut(bases, start, stop);
		for(final int left:new int[]{-10, 0, 10}){
			for(final int right:new int[]{-10, 0, 10}){
				final int s=start+left, e=stop+right;
				check(training.fill(bases, s, e, consensus, graph, 10, .37f, 0)==3
					&& rebuilt.fill(bases, s, e, consensus, graph, 10, .37f, 0)==3,
					"Both jitter extremes must retain all three positional probes and actual model features");
				for(int end=0; end<2; end++){
					check(Arrays.equals(training.values[end], rebuilt.values[end]),
						"Every optional training input must equal inference on a literal omitted-locus table");
					final byte[] candidate=Arrays.copyOfRange(bases, s, e+1);
					final float[] expected=TrnaBoundaryFeatures.tipFuzzinessFeature(candidate, graph, end==0);
					for(int i=0; i<3; i++){
						check(training.values[end][66+i]==expected[i],
							"Fuzziness must use the candidate-aligned HBM end and the unchanged sliding coverage definition");
					}
				}
			}
		}
		check(ordinary.fill(bases, start, stop, consensus, graph, 10, .37f, 0)==3,
			"The original complete span must have two usable optional vectors");
		for(int end=0; end<2; end++){
			check(ordinary.values[end].length==69 && ordinary.values[end][63]==1
				&& ordinary.values[end][64]==.37f && ordinary.values[end][65]==1,
				"The69inputs end with selected-consensus length ratio, full-contig GC, Quantum identity, then HBM coverage");
		}
		check(ordinary.values[0][66]==0 && ordinary.values[0][67]==.1f && ordinary.values[0][68]==.2f
			&& ordinary.values[1][66]==.2f && ordinary.values[1][67]==.1f && ordinary.values[1][68]==0,
			"Literal endpoint coverage is0/1/2and2/1/0, divided by the model's maximum10");
		final float[][] expected={ordinary.values[0].clone(), ordinary.values[1].clone()};
		final float[][] captured=new float[2][];
		final LongHashSet seeds=new LongHashSet();
		seeds.add(0);
		final NcrnaFamily family=new NcrnaFamily("cobalamin_rf00174", new byte[][]{consensus}, new BaseGraph[]{graph},
			new String[]{"m0 n=2"}, seeds, 17, 16, 64, 6, 2, false, 1, 0, 0, 1, 0, 1, .6f, .6f);
		final NcrnaPositionalBoundaryRefiner.Resources resources=new NcrnaPositionalBoundaryRefiner.Resources(family,
			five, three, new CaptureNet(8, 0, captured), new CaptureNet(13, 1, captured), 0);
		for(int strand=0; strand<2; strand++){
			final NcrnaBoundaryRefiner refiner=resources.create(family);
			final Orf orf=new Orf("fixture", start, stop, strand, 0, bases, false, ProkObject.RNA);
			orf.ncrnaFamily=family.name;
			orf.trnaModel=family.modelNames[0];
			check(refiner.refine(orf, bases, 0, .37f, .8f, true, 0, 16, 1000)
				&& orf.start==start-2 && orf.stop==stop+3,
				"Optional inference must apply both literal classes in sense coordinates on either genomic strand");
			check(Arrays.equals(expected[0], captured[0]) && Arrays.equals(expected[1], captured[1]),
				"Both networks must see the same UNMOVED span, not a stop vector recomputed after moving the start");
		}
		check(ordinary.fill(bases, start, stop, consensus, graph, 10, .37f, 1)==2,
			"Only the unfixed end is available for inference");
		for(final float value:ordinary.values[0]){check(Float.isNaN(value), "A fixed optional end must clear every stale feature");}
		final byte[] ambiguous=bases.clone();
		ambiguous[start-9]='N';
		check(ordinary.fill(ambiguous, start, stop, consensus, graph, 10, .37f, 0)==2,
			"An ambiguous exclusive-left probe disables that end instead of fabricating a score");
		final int word=RrnaPositionalKmerTable.encode(bases, start+five[0].offset-five[0].radius, five[0].k);
		final double savedScore=five[0].scores[word];
		for(final double bad:new double[]{Double.NaN, Double.POSITIVE_INFINITY, Double.MAX_VALUE}){
			boolean rejected=false;
			five[0].scores[word]=bad;
			try{ordinary.fill(bases, start, stop, consensus, graph, 10, .37f, 0);}
			catch(IllegalArgumentException expectedFailure){rejected=true;}
			finally{five[0].scores[word]=savedScore;}
			check(rejected, "A computed nonfinite or float-overflowed score must throw, never masquerade as a missing word");
		}
		final float[] missingFive=new float[3], missingThree=new float[3];
		TrnaBoundaryFeatures.tipFuzzinessFeatures(consensus, null, 0, new AlignmentStats(true), missingFive, 0, missingThree, 0);
		for(int i=0; i<3; i++){check(missingFive[i]==-1 && missingThree[i]==-1, "Missing models retain the shared minus-one sentinel");}
		final byte[] noAlignment=new byte[consensus.length];
		Arrays.fill(noAlignment, (byte)'N');
		final AlignmentStats diagnostic=new AlignmentStats(true);
		final float returnedIdentity=ScrabbleAligner.alignAndTraceStatic(noAlignment, consensus, diagnostic);
		System.err.println("Legacy all-N identity diagnostic: returned="+returnedIdentity+", traced="+diagnostic.identity);
		TrnaBoundaryFeatures.tipFuzzinessFeatures(noAlignment, graph, 10, new AlignmentStats(true), missingFive, 0, missingThree, 0);
		check(Arrays.equals(missingFive, TrnaBoundaryFeatures.tipFuzzinessFeature(noAlignment, graph, true))
			&& Arrays.equals(missingThree, TrnaBoundaryFeatures.tipFuzzinessFeature(noAlignment, graph, false)),
			"Ambiguous candidates must preserve the actual legacy helper, not an assumed rejection rule");
		// Ns have neutral Scrabble score and do not guarantee rejection by its
		// pre-trace confidence estimate. C versus an all-A model is a literal
		// zero-match fixture, independent of ambiguous-base score conventions.
		final byte[] allA=new byte[consensus.length];
		Arrays.fill(allA, (byte)'A');
		Arrays.fill(noAlignment, (byte)'C');
		final BaseGraph zeroMatch=model(allA);
		check(TrnaBoundaryFeatures.alignToModelFrame(noAlignment, zeroMatch)==null,
			"The unambiguous C/A fixture must actually take the no-alignment path before testing its sentinel");
		TrnaBoundaryFeatures.tipFuzzinessFeatures(noAlignment, zeroMatch, 10, new AlignmentStats(true), missingFive, 0, missingThree, 0);
		for(int i=0; i<3; i++){
			check(missingFive[i]==-2 && missingThree[i]==-2,
				"Zero-match C/A alignment must retain minus-two; actual="+missingFive[i]+","+missingThree[i]);
		}
	}

	private static BaseGraph model(final byte[] original){
		final BaseGraph model=new BaseGraph("m0", original, null, 0, 0);
		for(int position=0; position<original.length; position++){
			final int depth=position<2 ? position+1 : position>=original.length-2 ? original.length-position : 10;
			final byte base=original[position];
			BaseGraphHelper.setCounts(model, position, base=='A' ? depth : 0, base=='C' ? depth : 0,
				base=='G' ? depth : 0, base=='T' ? depth : 0, 0);
		}
		BaseGraphHelper.initForScoring(model);
		return model;
	}

	private static Table[] tables(final byte[] bases, final String end, final int first, final int truth, final int copies){
		final int[] offsets={first, -9, 1};
		final Table[] tables=new Table[3];
		for(int probe=0; probe<3; probe++){
			tables[probe]=new Table("RF00174", end, 9, offsets[probe], 1.0/(1<<18), 10);
			for(int copy=0; copy<copies; copy++){tables[probe].add(bases, truth);}
			tables[probe].finish();
		}
		return tables;
	}

	private static final class CaptureNet extends CellNet{
		CaptureNet(final int peak_, final int end_, final float[][] captured_){
			peak=peak_;
			end=end_;
			captured=captured_;
			output[peak]=1;
		}
		@Override
		public int numInputs(){return 69;}
		@Override
		public int numOutputs(){return 21;}
		@Override
		public CellNet copy(final boolean deltas){return new CaptureNet(peak, end, captured);}
		@Override
		public CellNet applyInput(final float[] values){
			captured[end]=values.clone();
			return this;
		}
		@Override
		public float feedForward(){return output[peak];}
		@Override
		public float[] getOutput(){return output;}
		private final int peak, end;
		private final float[][] captured;
		private final float[] output=new float[21];
	}

	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
