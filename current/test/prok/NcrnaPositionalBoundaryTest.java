package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import map.LongHashSet;
import ml.CellNet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaBoundaryRefiner;
import prok.NcrnaFamily;
import prok.NcrnaPositionalBoundaryFeatures;
import prok.NcrnaPositionalBoundaryRefiner;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.ProkObject;
import prok.RrnaPositionalKmerTable;
import shared.Shared;

/** Literal positional feature, worker ownership and post-trim inference oracle.
 * Homopolymer counts yield zero log enrichment at every position. Mock nets
 * choose explicit endpoint classes independently of the implementation.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPositionalBoundaryTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPositionalBoundaryTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown positional-fixture option: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.maxReads<0 && parser.samplerate==1
				&& !parser.overwrite && !parser.append && new File(parser.out1).mkdir(),
				"Supply out=fresh_fixture_directory");
			run(new File(parser.out1));
			pp.outstream.println("NCRNA_POSITIONAL_BOUNDARY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File directory){
		final byte[] bases=new byte[240];
		Arrays.fill(bases, (byte)'A');
		final byte[] original=bases.clone();
		final NcrnaFamily family=family();
		final RrnaPositionalKmerTable.Table five=table(bases, "5prime", 0, 80), three=table(bases, "3prime", -8, 139);
		final String fivePath=new File(directory, "five.tsv.gz").getPath(), threePath=new File(directory, "three.tsv.gz").getPath();
		five.write(fivePath);
		three.write(threePath);
		final RrnaPositionalKmerTable.Table loadedFive=RrnaPositionalKmerTable.loadScores(fivePath, "RF00174");
		final RrnaPositionalKmerTable.Table loadedThree=RrnaPositionalKmerTable.loadScores(threePath, "RF00174");
		final NcrnaPositionalBoundaryFeatures training=new NcrnaPositionalBoundaryFeatures(five, three);
		final NcrnaPositionalBoundaryFeatures inference=new NcrnaPositionalBoundaryFeatures(loadedFive, loadedThree);
		check(training.fill(bases, 80, 139, family.library[0], .37f, 0)==3
			&& inference.fill(bases, 80, 139, family.library[0], .37f, 0)==3,
			"Both ends of the complete candidate must have usable native positional features");
		for(int end=0; end<2; end++){
			check(Arrays.equals(training.values[end], inference.values[end]), "Training/full-table and inference/loaded-score features must agree bit-for-bit");
			for(int site=0; site<13; site++){
				check(Math.abs(training.values[end][site])<1e-12,
					"Uniform words have analytic log enrichment zero; summing twelve double frequencies permits rounding: "+training.values[end][site]);
			}
			check(training.values[end][13]==1 && training.values[end][14]==.37f && training.values[end][15]==1,
				"Extras must be literal60/60 selected-consensus length, measured whole-contig GC and candidate-span Quantum identity");
		}
		training.fill(bases, 80, 139, family.library[1], .37f, 0);
		check(training.values[0][13]==60f/90, "The same member table must normalize by the selected consensus, never member median");
		final byte[] missing=bases.clone();
		missing[80]='N';
		check(training.fill(missing, 80, 139, family.library[0], .37f, 0)==2,
			"A missing five-prime word cannot produce a fabricated usable row; the other end remains available");
		check(training.fill(bases, 80, 139, family.library[0], .37f, 1)==2 && Float.isNaN(training.values[0][0]),
			"A fixed end must clear stale scratch and never become an inference row");
		final int[] calls=new int[2];
		final NcrnaPositionalBoundaryRefiner.Resources resources=new NcrnaPositionalBoundaryRefiner.Resources(family,
			loadedFive, loadedThree, new PositionNet(4, false, calls), new PositionNet(9, false, calls), .02f);
		final NcrnaBoundaryRefiner first=resources.create(family), second=resources.create(family);
		check(calls[1]==4 && first!=second && first.afterAlignmentTrim(), "Every worker owns two network copies and explicitly uses the post-trim seam");
		for(int strand=0; strand<2; strand++){
			final Orf candidate=candidate(bases, strand, 80, 139);
			check((strand==0 ? first : second).refine(candidate, bases, 0, .37f, .8f, true, 0, 16, 1000)
				&& candidate.start==78 && candidate.stop==142, "Classes4/9 at radius6 mean raw start-2/stop+3 on either sense-oriented strand");
			preserved(candidate);
			final TrimSpy caller=new TrimSpy(family);
			caller.setBoundaryRefiner(resources.create(family));
			final Orf raw=candidate(bases, strand, 75, 145);
			check(caller.finishAcceptedBoundary(raw, bases, 0, 0, bases.length-1, 0, .8f, 71, .8f, true, null, null)
				&& caller.trims==1 && raw.start==78 && raw.stop==142, "Actual dispatch must refine the literal80..139 trim result, not the original75..145 window");
			preserved(raw);
		}
		for(int mask=1; mask<=3; mask++){
			final Orf candidate=candidate(bases, 0, 80, 139);
			final int before=calls[0];
			first.refine(candidate, bases, 0, .37f, .8f, true, mask, 16, 1000);
			check(candidate.start==((mask&1)!=0 ? 80 : 78) && candidate.stop==((mask&2)!=0 ? 139 : 142),
				"Each fixed-end bit must retain its original coordinate");
			check(calls[0]-before==(mask==3 ? 0 : 1), "Only the unfixed usable end evaluates its network");
		}
		final Orf tooLong=candidate(bases, 0, 80, 139);
		check(!first.refine(tooLong, bases, 0, .37f, .8f, true, 0, 16, 64)
			&& tooLong.start==80 && tooLong.stop==139, "An invalid joint65-base prediction must retain the accepted60-base span");
		final NcrnaBoundaryRefiner uncertain=new NcrnaPositionalBoundaryRefiner.Resources(family, loadedFive, loadedThree,
			new PositionNet(4, true, calls), new PositionNet(9, true, calls), .005f).create(family);
		final Orf held=candidate(bases, 0, 80, 139);
		check(!uncertain.refine(held, bases, 0, .37f, .8f, true, 0, 16, 1000), "Equal best alternatives must retain both caller endpoints under a positive margin");
		fails(()->resources.create(family()));
		fails(()->new NcrnaPositionalBoundaryFeatures(loadedThree, loadedFive));
		final RrnaPositionalKmerTable.Table wrongFive=table(bases, "RF01482", "5prime", 0, 80);
		final RrnaPositionalKmerTable.Table wrongThree=table(bases, "RF01482", "3prime", -8, 139);
		fails(()->new NcrnaPositionalBoundaryRefiner.Resources(family, wrongFive, wrongThree,
			new PositionNet(4, false, calls), new PositionNet(9, false, calls), 0));
		fails(()->new NcrnaPositionalBoundaryRefiner.Resources(family, loadedFive, loadedThree,
			new PositionNet(4, false, calls){@Override public int numInputs(){return 17;}}, new PositionNet(9, false, calls), 0));
		check(Arrays.equals(original, bases), "Feature computation and refinement must not mutate genomic bases");
	}

	private static RrnaPositionalKmerTable.Table table(final byte[] bases, final String end, final int anchor, final int truth){
		return table(bases, "RF00174", end, anchor, truth);
	}
	private static RrnaPositionalKmerTable.Table table(final byte[] bases, final String member,
			final String end, final int anchor, final int truth){
		final RrnaPositionalKmerTable.Table table=new RrnaPositionalKmerTable.Table(member, end, 9, anchor, .5, 6);
		table.add(bases, truth);
		table.finish();
		return table;
	}
	private static NcrnaFamily family(){
		final byte[][] library={new byte[60], new byte[90]};
		for(final byte[] sequence:library){Arrays.fill(sequence, (byte)'A');}
		final LongHashSet seeds=new LongHashSet();
		seeds.add(0);
		return new NcrnaFamily("cobalamin_rf00174", library, null, NAMES, seeds, 17, 16, 64, 6, 2, false, 1, 0, 0, 1, 0, 1, .6f, .6f);
	}
	private static Orf candidate(final byte[] bases, final int strand, final int start, final int stop){
		final Orf result=new Orf("fixture", start, stop, strand, 0, bases, false, ProkObject.RNA);
		result.ncrnaFamily="cobalamin_rf00174";
		result.trnaModel=NAMES[0];
		result.orfScore=42;
		return result;
	}
	private static void preserved(final Orf candidate){
		check(candidate.orfScore==42 && NAMES[0].equals(candidate.trnaModel)
			&& "cobalamin_rf00174".equals(candidate.ncrnaFamily), "Refinement preserves detection score and exact descriptive model/family labels");
	}
	private static final class TrimSpy extends NcrnaScavenger{
		TrimSpy(final NcrnaFamily family){super(family.library, null, family.modelNames, family.kmerSet, 17, 16, 64);}
		@Override
		public boolean trimToAlignmentExtent(final Orf candidate, final byte[] bases, final int model,
				final int start, final int stop, final float identity){
			trims++;
			candidate.start=80;
			candidate.stop=139;
			return true;
		}
		int trims;
	}
	private static class PositionNet extends CellNet{
		PositionNet(final int peak_, final boolean tie_, final int[] calls_){
			peak=peak_;
			tie=tie_;
			calls=calls_;
			Arrays.fill(output, .125f);
			output[peak]=.875f;
			if(tie){output[peak+1]=.875f;}
		}
		@Override
		public int numInputs(){return 16;}
		@Override
		public int numOutputs(){return 13;}
		@Override
		public CellNet copy(final boolean deltas){
			calls[1]++;
			return new PositionNet(peak, tie, calls);
		}
		@Override
		public CellNet applyInput(final float[] values){
			NcrnaPositionalBoundaryTest.check(values.length==16, "Mock receives the native radius6 feature width");
			return this;
		}
		@Override
		public float feedForward(){
			calls[0]++;
			return output[peak];
		}
		@Override
		public float[] getOutput(){return output;}
		final int peak;
		final boolean tie;
		final int[] calls;
		final float[] output=new float[13];
	}
	private static void fails(final Runnable call){
		boolean failed=false;
		try{call.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Wrong positional binding or schema must fail before inference");
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private static final String[] NAMES={"selected_model n=1", "longer_model n=1"};
}
