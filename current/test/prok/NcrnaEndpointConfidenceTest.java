package test.prok;

import java.util.Arrays;
import java.util.Locale;

import ml.CellNet;
import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEndpointConfidence;
import prok.NcrnaBoundaryScorer;
import prok.TrnaBoundaryFeatures;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.ProkObject;
import shared.Shared;

/** Literal score-window oracle for the shared optional confidence rule.
 * Legacy scalar ties keep current; legacy positional ties keep first argmax.
 * These historical choices stay distinct at zero margin. Full-error inference
 * instead holds every tie, including ties with its non-positional veto head.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointConfidenceTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointConfidenceTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown confidence-fixture option: "+arg);
			}
			check(parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null
				&& parser.maxReads<0 && parser.samplerate==1, "Confidence fixture has fixed literal data; input/output or sampling overrides are unsupported");
			run();
			pp.outstream.println("NCRNA_ENDPOINT_CONFIDENCE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private static void run(){
		final float[] tied={.8f, .8f, .2f};
		check(NcrnaEndpointConfidence.select(tied, 1, 1, 0)==1, "Zero margin preserves scalar no-move on a current-position tie");
		check(NcrnaEndpointConfidence.select(tied, 1, 0, 0)==0, "Zero margin preserves positional first-argmax ties");
		check(NcrnaEndpointConfidence.select(tied, 1, 0, .005f)==1, "A best/current tie cannot support a positive confidence move");
		final float[] ambiguous={.9f, .8995f, .2f};
		check(NcrnaEndpointConfidence.select(ambiguous, 2, 0, .005f)==2,
			"A large improvement over current is inconclusive when runner-up is almost as good");
		final float[] smallGain={.9f, .2f, .8995f};
		check(NcrnaEndpointConfidence.select(smallGain, 2, 0, .005f)==2,
			"A best position must convincingly beat current even with a distant alternative");
		final float[] clear={.875f, .5f, .625f, 0};
		final boolean[] clearAvailable={true, true, true, false};
		final float[] saved=clear.clone();
		check(NcrnaEndpointConfidence.select(clear, 2, 0, .02f, clearAvailable)==0,
			"A unique best exceeding current and runner-up by more than margin may move");
		check(NcrnaEndpointConfidence.select(clear, 2, 0, .25f, clearAvailable)==2,
			"A gap exactly equal to the margin is inconclusive, not a strict improvement");
		check(Arrays.equals(clear, saved), "Confidence selection borrows scores without mutating them");
		final float[] runnerTie={.875f, .875f, .25f};
		check(NcrnaEndpointConfidence.select(runnerTie, 2, 0, .005f)==2, "Tied best alternatives must hold even when current is worse");
		final float[] negative={-.125f, -.5f, -.375f};
		check(NcrnaEndpointConfidence.select(negative, 2, 0, .02f)==0,
			"Margins are score differences and do not assume outputs are positive probabilities");
		final float[] unavailable={0, .125f, .625f};
		check(NcrnaEndpointConfidence.select(unavailable, 1, 2, .1f, new boolean[]{false, true, true})==2,
			"Unavailable positions do not become fictitious runners-up or neutral scores");
		for(final float margin:new float[]{0, .005f, .02f, .1f}){
			check(NcrnaEndpointConfidence.select(clear, 2, 2, margin, clearAvailable)==2,
				"An ungated no-move choice cannot become a new location when confidence increases");
		}
		fails(()->NcrnaEndpointConfidence.select(null, 0, 0, 0));
		fails(()->NcrnaEndpointConfidence.select(clear, -1, 0, .1f, clearAvailable));
		fails(()->NcrnaEndpointConfidence.select(clear, 2, 4, .1f, clearAvailable));
		fails(()->NcrnaEndpointConfidence.select(clear, 3, 0, .1f, clearAvailable));
		fails(()->NcrnaEndpointConfidence.select(clear, 2, 3, .1f, clearAvailable));
		fails(()->NcrnaEndpointConfidence.select(clear, 2, 0, Float.NaN, clearAvailable));
		fails(()->NcrnaEndpointConfidence.select(clear, 2, 0, -.1f, clearAvailable));
		fails(()->NcrnaEndpointConfidence.select(new float[]{.8f, Float.POSITIVE_INFINITY, .2f}, 2, 0, .1f));
		fails(()->NcrnaEndpointConfidence.select(new float[]{Float.NaN, .5f, .25f}, 1, 1, 0));
		fails(()->NcrnaEndpointConfidence.select(new float[]{.8f, .5f, Float.POSITIVE_INFINITY}, 1, 0, 0));
		check(NcrnaEndpointConfidence.select(new float[]{-Float.MAX_VALUE, 0, Float.MAX_VALUE}, 1, 2, 0)==2,
			"Every finite signed score remains a measurement, including both representable extremes");
		clippedSweep();
		scalarSweeps();
		familyDispatch();
		fullErrors();
		signedHeadMargins();
	}

	/** Dyadic gaps distinguish a signed head allowance from relaxed positional ties. */
	private static void signedHeadMargins(){
		final float[] scores={.75f, .25f, .5f, 1};
		final float[] saved=scores.clone();
		check(NcrnaEndpointConfidence.selectFull(scores, 1, 0, -.25f)==1,
			"A head gap exactly equal to the signed threshold must hold");
		check(NcrnaEndpointConfidence.selectFull(scores, 1, 0, Math.nextDown(-.25f))==0,
			"A strictly looser negative head margin permits the unique reachable winner");
		check(NcrnaEndpointConfidence.selectFull(scores, 1, .25f, 0, false)==1,
			"Ignoring the head must not relax the positional runner-up margin");
		check(NcrnaEndpointConfidence.selectFull(scores, 1, 0, 0, false)==0,
			"Ignore-head returns the unique best reachable offset even below the extra head");
		check(NcrnaEndpointConfidence.selectFull(scores, 1, 0, -1)==0,
			"Negative one is a real signed margin, not absence or ignore state");
		check(NcrnaEndpointConfidence.selectFull(scores, 1, 0, 0)==1,
			"An explicit numeric zero still compares the head and differs from ignoring it");
		check(Arrays.equals(scores, saved), "Signed-head selection cannot modify the borrowed score buffer");
		scores[2]=.75f;
		check(NcrnaEndpointConfidence.selectFull(scores, 1, 0, 0, false)==1,
			"Two reachable winners hold even when the out-of-range head is ignored");
		scores[1]=1;
		check(NcrnaEndpointConfidence.selectFull(scores, 1, 0, -100)==1, "Current-best must hold at any signed head allowance");
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 1, -.1f, -.1f));
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 1, 0, Float.NaN));
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 1, 0, Float.POSITIVE_INFINITY));
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 1, 0, Float.NEGATIVE_INFINITY));
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 1, 0, Float.NaN, false));
		scores[3]=Float.NaN;
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 1, 0, 0, false));
	}

	/** Literal dyadic scores distinguish strict gaps from rounded float subtraction. */
	private static void fullErrors(){
		final float[] scores=new float[22];
		Arrays.fill(scores, -.5f);
		scores[10]=.25f;
		scores[4]=.75f;
		scores[21]=.5f;
		final float[] saved=scores.clone();
		check(NcrnaEndpointConfidence.selectFull(scores, 10, 0)==4, "A strict reachable winner moves at zero margin");
		check(NcrnaEndpointConfidence.selectFull(scores, 10, .25f)==10, "A gap equal to the out-of-range margin holds");
		check(Arrays.equals(scores, saved), "Full-error selection must not modify borrowed scores");
		scores[21]=.75f;
		check(NcrnaEndpointConfidence.selectFull(scores, 10, 0)==10, "An out-of-range tie vetoes even at zero margin");
		scores[21]=1;
		check(NcrnaEndpointConfidence.selectFull(scores, 10, 0)==10, "The out-of-range head must never become a coordinate");
		scores[21]=0;
		scores[7]=.75f;
		check(NcrnaEndpointConfidence.selectFull(scores, 10, 0)==10, "A tied pair of reachable alternatives holds at zero margin");
		scores[7]=.625f;
		check(NcrnaEndpointConfidence.selectFull(scores, 10, .125f)==10, "A runner gap equal to the margin holds");
		check(NcrnaEndpointConfidence.selectFull(scores, 10, Math.nextDown(.125f))==4, "A strictly smaller margin permits the reachable winner");
		scores[10]=.75f;
		check(NcrnaEndpointConfidence.selectFull(scores, 10, 0)==10, "A current-position tie cannot move an already best end");
		fails(()->NcrnaEndpointConfidence.selectFull(null, 10, 0));
		fails(()->NcrnaEndpointConfidence.selectFull(new float[21], 10, 0));
		fails(()->NcrnaEndpointConfidence.selectFull(scores, -1, 0));
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 10, -.1f));
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 10, Float.NaN));
		scores[21]=Float.NaN;
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 10, 0));
		scores[21]=Float.POSITIVE_INFINITY;
		fails(()->NcrnaEndpointConfidence.selectFull(scores, 10, 0));
	}

	/** Checks the actual family-control handoff, not just policy storage and scorer helpers. */
	private static void familyDispatch(){
		final byte[] bases=new byte[240], consensus=new byte[60];
		Arrays.fill(bases, (byte)'A');
		Arrays.fill(consensus, (byte)'A');
		for(final boolean shared:new boolean[]{false, true}){
			final TrnaBoundaryFeatures.NinemerTable table=(seq, pos)->0;
			final LongHashSet seeds=new LongHashSet();
			seeds.add(0);
			final int[] offsets={-1, 0, 1};
			final NcrnaFamily family=new NcrnaFamily("sam1", new byte[][]{consensus}, null,
				new String[]{"sam1_fixture"}, seeds, 17, 16, 64, 6, 1, false, 1, 0, 0, 1,
				0, 1, .6f, .6f, 1.01f, .9f, new ScoreNet(.9f, .2f, .8995f, true),
				new ScoreNet(1, 1, 1, true), table, table, 9, 0, 9, 0, 1000, offsets, offsets, .005f, .005f);
			family.boundarySharedConfidence=shared;
			final FamilyTrimSpy caller=new FamilyTrimSpy(family);
			GeneCaller.applyNcrnaFamilyControls(caller, family);
			final Orf candidate=new Orf("fixture", 75, 145, 0, 0, bases, false, ProkObject.RNA);
			candidate.ncrnaFamily=family.name;
			candidate.trnaModel=family.modelNames[0];
			candidate.orfScore=42;
			check(caller.finishAcceptedBoundary(candidate, bases, 0, 0, 239, 0, .8f, 60, .8f, true, null, null)
				&& caller.trims==1 && candidate.start==(shared ? 80 : 79) && candidate.stop==139 && candidate.orfScore==42,
				"Real family controls must connect the opt-in scalar gate: OFF moves an ambiguous end, ON holds after actual trim dispatch");
		}
	}

	private static final class FamilyTrimSpy extends NcrnaScavenger{
		FamilyTrimSpy(final NcrnaFamily family){
			super(family.library, family.models, family.modelNames, family.kmerSet, family.kLong, family.minLen, family.windowPad,
				family.indexK, family.indexTopN, family.adaptive, family.adaptFloor, family.adaptTopFrac, family.adaptQFrac,
				family.fixedMinHits, family.scoreA, family.scoreB, family.idPass, family.idBorderline,
				family.boundary5NetTemplate, family.boundary3NetTemplate, family.boundaryStartTable, family.boundaryStopTable,
				family.boundaryStartInside, family.boundaryStartOutside, family.boundaryStopInside, family.boundaryStopOutside,
				family.boundaryMeanLen, family.boundaryStartOffsets, family.boundaryStopOffsets, family.boundaryMarginStart, family.boundaryMarginStop);
		}
		@Override
		public boolean trimToAlignmentExtent(final Orf candidate, final byte[] bases, final int model,
				final int start, final int stop, final float identity){
			trims++;
			candidate.start=80;
			candidate.stop=139;
			return true;
		}
		private int trims;
	}

	/** The off-contig start offset is absent, while both legal sites retain real scores. */
	private static void clippedSweep(){
		final byte[] bases=new byte[240];
		Arrays.fill(bases, (byte)'A');
		final TrnaBoundaryFeatures.NinemerTable table=(seq, pos)->0;
		final int[] offsets={-1, 0, 1}, captures={0};
		final NcrnaBoundaryScorer.ScoreSink sink=(end, context, start, stop, values, available)->{
			check(end==0 && start==0 && stop==59 && available!=null && !available[0] && available[1] && available[2]
				&& values[1]==.2f && values[2]==.5f, "Actual scalar scoring must mark the off-contig site absent and keep exact current/alternate scores");
			captures[0]++;
		};
		final NcrnaBoundaryScorer.RefinementResult result=NcrnaBoundaryScorer.refineBoundariesWithCounts(
			new ScoreNet(.1f, .2f, .5f), new ScoreNet(1, 1, 1), bases, 0, 59, new byte[60], null,
			table, table, 9, 0, 9, 0, .5f, 1000f, offsets, offsets, 0, 0, true, false,
			NcrnaBoundaryScorer.FEATURES_V3, .8f, null, null, null, null, 2, Integer.MAX_VALUE, true, sink);
		check(captures[0]==1 && result.startSitesScored==2 && result.stopSitesScored==0
			&& result.startOffset==1 && result.stopOffset==0 && result.stopChosenScore==0,
			"Fixed ends have no score; legal start candidates preserve the literal scalar choice");
	}

	/** Literal length-dependent outputs exercise the real sequential scalar sweep. */
	private static void scalarSweeps(){
		for(final boolean stopFirst:new boolean[]{false, true}){
			final NcrnaBoundaryScorer.RefinementResult old=scalar(stopFirst, .9f, .8995f, .005f, false, 0);
			final NcrnaBoundaryScorer.RefinementResult zero=scalar(stopFirst, .9f, .8995f, 0, true, 0);
			final NcrnaBoundaryScorer.RefinementResult held=scalar(stopFirst, .9f, .8995f, .005f, true, 0);
			check(stopFirst ? old.stopOffset==1 && old.startOffset==0 : old.startOffset==-1 && old.stopOffset==0,
				"Legacy scalar margins compare only current; the closer alternative previously did not prevent this move");
			check(zero.startOffset==old.startOffset && zero.stopOffset==old.stopOffset,
				"Explicit shared confidence at zero preserves the actual scalar selector in both sequential orders");
			check(held.startOffset==0 && held.stopOffset==0
				&& (stopFirst ? held.stopChosenScore : held.startChosenScore)==.2f,
				"Positive shared confidence holds an ambiguous scalar end and reports its retained current score");
			final NcrnaBoundaryScorer.RefinementResult clear=scalar(stopFirst, .875f, .5f, .02f, true, 0);
			check(clear.startOffset==old.startOffset && clear.stopOffset==old.stopOffset,
				"A clear scalar winner may move in either sequential order");
			final NcrnaBoundaryScorer.RefinementResult tie=scalar(stopFirst, .875f, .875f, .005f, true, 0);
			final NcrnaBoundaryScorer.RefinementResult exact=scalar(stopFirst, .875f, .625f, .25f, true, 0);
			check(tie.startOffset==0 && tie.stopOffset==0 && exact.startOffset==0 && exact.stopOffset==0,
				"Tied alternatives and an exactly equal runner-up gap must not establish confidence");
			final NcrnaBoundaryScorer.RefinementResult fixed=scalar(stopFirst, .9f, .8995f, .005f, true, stopFirst ? 1 : 2);
			check(fixed.startOffset==0 && fixed.stopOffset==0,
				"A fixed physical end stays fixed while its ambiguous opposite end holds");
		}
	}

	private static NcrnaBoundaryScorer.RefinementResult scalar(final boolean stopFirst, final float best,
			final float runnerUp, final float margin, final boolean shared, final int fixed){
		final byte[] bases=new byte[240];
		Arrays.fill(bases, (byte)'A');
		final ScoreNet moving=new ScoreNet(best, .2f, runnerUp), constant=new ScoreNet(1, 1, 1);
		final TrnaBoundaryFeatures.NinemerTable table=(seq, pos)->0;
		final int[] offsets={-1, 0, 1};
		final NcrnaBoundaryScorer.RefinementResult result=NcrnaBoundaryScorer.refineBoundariesWithCounts(
			stopFirst ? constant : moving, stopFirst ? moving : constant, bases, 80, 139, new byte[60], null,
			table, table, 9, 0, 9, 0, .5f, 1000f, offsets, offsets, margin, margin,
			false, false, NcrnaBoundaryScorer.FEATURES_V3, .8f, null, null, null, null, fixed, Integer.MAX_VALUE, shared);
		check(moving.calls+constant.calls==(fixed==0 ? 8 : 3),
			"Runner-up confidence must not add a network evaluation to either ordinary or fixed-end sweeps");
		if(fixed==0){
			final boolean moved=result.startOffset!=0 || result.stopOffset!=0;
			final int center=moved ? 61 : 60;
			check(constant.calls==4 && constant.lengths[0]==60
				&& constant.lengths[1]==center+(stopFirst ? 1 : -1)
				&& constant.lengths[2]==center && constant.lengths[3]==center+(stopFirst ? -1 : 1),
				"Second-end inputs must use the actually held or moved first end in both sequential orders");
		}
		return result;
	}

	/** The V3 length feature identifies three literal candidate scores, independent of the gate. */
	private static final class ScoreNet extends CellNet{
		ScoreNet(final float longer_, final float current_, final float shorter_){
			this(longer_, current_, shorter_, false);
		}
		ScoreNet(final float longer_, final float current_, final float shorter_, final boolean v1_){
			longer=longer_;
			current=current_;
			shorter=shorter_;
			v1=v1_;
		}
		@Override
		public CellNet copy(final boolean deltas){return new ScoreNet(longer, current, shorter, v1);}
		@Override
		public CellNet applyInput(final float[] values){
			NcrnaEndpointConfidenceTest.check(values.length==10, "Scalar oracle uses the real ten-feature V3 contract");
			length=Math.round(values[v1 ? 8 : 4]*1000);
			return this;
		}
		@Override
		public float feedForward(){
			lengths[calls++]=length;
			return length>60 ? longer : length<60 ? shorter : current;
		}
		private final float longer, current, shorter;
		private final boolean v1;
		private final int[] lengths=new int[4];
		private int length, calls;
	}
	private static void fails(final Runnable call){
		boolean failed=false;
		try{call.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Malformed confidence inputs must fail rather than invent a move");
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
}
