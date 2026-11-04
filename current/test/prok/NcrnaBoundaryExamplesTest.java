package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import consensus.BaseGraph;
import map.LongHashSet;
import ml.CellNet;
import prok.NcrnaBoundaryExamples;
import prok.NcrnaBoundaryScorer;
import prok.NcrnaEndpointTable;
import prok.NcrnaFamily;
import prok.NcrnaLearningContext;
import prok.TrnaBoundaryFeatures;
import prok.TrnaBoundaryFeatures.BoundaryType;
import stream.Read;

/** Compares candidate features with the actual inference entry point, before network evaluation.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBoundaryExamplesTest{

	public static void main(String[] args){
		assert(args.length==0) : "Feature fixture uses literal small sequences and models";
		final byte[] gene=bytes("ACGTTGCATGTCAGTACGATCGTAGCTAGCATGACTACGA");
		final byte[][] library={gene, bytes("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")};
		final BaseGraph[] models={new BaseGraph("exact", library[0], null, 0, 0), new BaseGraph("other", library[1], null, 1, 0)};
		final NcrnaFamily family=new NcrnaFamily("fixture", library, models, new String[]{"exact", "other"},
			new LongHashSet(), 17, 16, 8, 7, 2, false, 0, 0, 0, 1);
		final byte[] window=bytes("GGGG"+new String(gene, StandardCharsets.US_ASCII)+"CCCC");
		final NcrnaLearningContext context=new NcrnaLearningContext();
		context.read=new Read(window, null, "fixture", 0, false);
		context.start=4;
		context.stop=3+gene.length;
		context.gc=0.5f;
		context.status=NcrnaLearningContext.AVAILABLE;
		final NcrnaEndpointTable start=new NcrnaEndpointTable(BoundaryType.START, 3, 0);
		final NcrnaEndpointTable stop=new NcrnaEndpointTable(BoundaryType.STOP, 2, 1);
		start.add(window, context.start);
		stop.add(window, context.stop);
		final NcrnaBoundaryExamples examples=new NcrnaBoundaryExamples(family, start, stop);
		assert(examples.upperMedianLength==library[1].length) : "Two-model normalization uses the upper median, not their arithmetic mean";
		assert(examples.prepare(context, false) && examples.model()==0) : "Exact accepted consensus must win mapping shortlist plus Quantum selection";
		final Recorder recorder=new Recorder();
		for(boolean isStop:new boolean[]{false, true}){
			for(int offset=-3; offset<=3; offset++){
				assert(examples.fill(isStop, offset)==NcrnaBoundaryExamples.VALID) : "All fourteen candidate shifts fit this context";
				final NcrnaEndpointTable table=isStop ? stop : start;
				NcrnaBoundaryScorer.score(recorder, window, context.start+(isStop ? 0 : offset),
					context.stop+(isStop ? offset : 0), isStop, library[0], models[0], table,
					table.inside, table.outside, context.gc, examples.upperMedianLength);
				assert(Arrays.equals(examples.values, recorder.input))
					: "Every candidate's ten training values must equal the actual caller input in the same order";
				assert(examples.label==(offset==0 ? 1 : 0) && examples.values[4]==(isStop ? 1 : 0)
					&& examples.values[8]==(gene.length+(isStop ? offset : -offset))/examples.upperMedianLength
					&& examples.values[9]==0.5f) : "Labels, endpoint direction, candidate length and whole-contig GC have literal expectations";
			}
		}
		assert(examples.prepare(context, true)) : "The same eligible training record remains alignable";
		assert(examples.fill(false, 0)==NcrnaBoundaryExamples.VALID && examples.label==1
			&& examples.values[1]==0 && examples.values[2]==1 && examples.values[3]==1)
			: "Removing the sole true word leaves zero observed counts and the inference pseudocount profile [0,1,1]";
		assert(start.rawCount(window, context.start)==1) : "Training view must not mutate the table used for validation";
		NcrnaBoundaryScorer.score(recorder, window, context.start, context.stop, false,
			library[0], null, start, start.inside, start.outside, context.gc, examples.upperMedianLength);
		assert(recorder.input[0]==1 && recorder.input[5]==TrnaBoundaryFeatures.PENDING_DORI
			&& recorder.input[6]==TrnaBoundaryFeatures.PENDING_DORI && recorder.input[7]==TrnaBoundaryFeatures.PENDING_DORI)
			: "Missing model preserves exact ANI and the distinct unavailable-fuzziness sentinel";
		checkEnrichment();
		assert(examples.fill(false, -5)==NcrnaBoundaryExamples.BEFORE_CONTEXT) : "A start before the context is not emitted";
		assert(examples.fill(true, 5)==NcrnaBoundaryExamples.AFTER_CONTEXT) : "A stop after the context is not emitted";
		assert(examples.fill(false, gene.length-15)==NcrnaBoundaryExamples.TOO_SHORT) : "Fifteen-base candidates are below the sixteen-base minimum";
		final byte[] unknown=new byte[window.length];
		Arrays.fill(unknown, (byte)'N');
		context.read=new Read(unknown, null, "ambiguous", 1, false);
		assert(!examples.prepare(context, false) && examples.model()==-1) : "No shared mapping word yields an explicit whole-locus skip";
		System.out.println("Boundary examples passed: accepted-model selection, fourteen candidate inputs equal to inference, labels, upper median, leave-one-out and candidate limits");
	}

	/** Five literal local counts define the expected center and rival ratios. */
	private static void checkEnrichment(){
		final int[] counts={0, 1, 2, 3, 4};
		final TrnaBoundaryFeatures.NinemerTable table=(sequence, position)->
			position<0 || position>=counts.length ? 0 : counts[position];
		final float magnitude=(float)(Math.log(5)/Math.log(2)/8);
		final float[] inference=TrnaBoundaryFeatures.enrichmentProfile(null, 2, BoundaryType.START, 1, 0, table);
		assert(Arrays.equals(inference, new float[]{magnitude, .6f, 1}))
			: "Counts[0,1,2,3,4] give F5, center3/5 and rival5/5";
		final float[] training=TrnaBoundaryFeatures.enrichmentProfile(null, 2, 2, BoundaryType.START, 1, 0, table);
		assert(Arrays.equals(training, new float[]{magnitude, .4f, 1}))
			: "Excluding one true-center observation changes only center3/5 to2/5";
		final float[] missing=TrnaBoundaryFeatures.enrichmentProfile(null, 2, BoundaryType.START, 1, 0, null);
		assert(Arrays.equals(missing, new float[]{-1, -1, -1})) : "Absent tables have three unavailable-feature sentinels";
		final float[] scratch={99, 99, 99, 1234};
		assert(TrnaBoundaryFeatures.enrichmentProfile(null, 2, BoundaryType.START, 1, 0, table, scratch)==scratch)
			: "The buffer overload must return caller-owned storage";
		assert(scratch[0]==magnitude && scratch[1]==.6f && scratch[2]==1 && scratch[3]==1234)
			: "The inference profile overwrites exactly its three output slots";
		TrnaBoundaryFeatures.enrichmentProfile(null, 2, BoundaryType.START, 1, 0, null, scratch);
		assert(scratch[0]==-1 && scratch[1]==-1 && scratch[2]==-1 && scratch[3]==1234)
			: "Absent-table sentinels replace prior values without touching adjacent storage";
		TrnaBoundaryFeatures.enrichmentProfile(null, 2, 2, BoundaryType.START, 1, 0, table, scratch);
		assert(scratch[0]==magnitude && scratch[1]==.4f && scratch[2]==1 && scratch[3]==1234)
			: "Training reuse must replace all prior sentinels and retain the literal exclusion result";
	}

	/** Records exactly what the real scorer hands to the network, without imposing a learned decision. */
	private static final class Recorder extends CellNet{
		@Override
		public CellNet applyInput(float[] values){
			input=values.clone();
			return this;
		}
		@Override
		public float feedForward(){return 0.5f;}
		float[] input;
	}

	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}
}
