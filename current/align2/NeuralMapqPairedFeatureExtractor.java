package align2;

import java.util.Arrays;

import dna.Data;
import stream.Read;
import stream.SiteScore;

/**
 * Extracts 108 raw fields: anchor 42, mate 42, and 24 pair fields.
 * Reads must be reciprocal mates. Each mapped end needs the single-end
 * extractor's primary/site/traceback state; an unmapped mate's raw block is zero.
 * Scratch and output arrays belong to one caller thread and are reused.
 * No truth labels or neural predictions enter these features.
 *
 * @author Collei
 */
public final class NeuralMapqPairedFeatureExtractor{

	private NeuralMapqPairedFeatureExtractor(){}

	/** Overwrites all fields for one mapped primary anchor; the mate may be unmapped. */
	public static void fill(final Read anchor, final Read mate,
			final int averagePairDistance, final boolean requireCorrectStrands,
			final boolean sameStrandPairs, final float[] vector, final Scratch scratch){
		fill(anchor, mate, averagePairDistance, requireCorrectStrands, sameStrandPairs, vector, scratch, false);
	}

	/** Fills V1 inference features; composition fields 37..41 in each end block are zero.
	 * Use fill for full raw exports or transforms that consume composition. */
	static void fillRuntime(final Read anchor, final Read mate,
			final int averagePairDistance, final boolean requireCorrectStrands,
			final boolean sameStrandPairs, final float[] vector, final Scratch scratch){
		fill(anchor, mate, averagePairDistance, requireCorrectStrands, sameStrandPairs, vector, scratch, true);
	}

	private static void fill(final Read anchor, final Read mate,
			final int averagePairDistance, final boolean requireCorrectStrands,
			final boolean sameStrandPairs, final float[] vector, final Scratch scratch,
			final boolean runtime){
		validatePair(anchor, mate, averagePairDistance, vector, scratch);
		fillEnd(anchor, scratch.anchorRaw, scratch.anchorScratch, runtime);
		if(mate.mapped()){
			fillEnd(mate, scratch.mateRaw, scratch.mateScratch, runtime);
		}else{Arrays.fill(scratch.mateRaw, 0);}
		fillPrepared(anchor, mate, averagePairDistance, requireCorrectStrands, sameStrandPairs,
				scratch.anchorRaw, scratch.mateRaw, vector);
	}

	/** Extracts each mapped end once and fills its anchor orientation.
	 * At least one end must be mapped. Both output buffers must have width 108
	 * and must be distinct when both ends are mapped. An unmapped end's output
	 * is left untouched, not zeroed: callers must check mapped() before using it. */
	public static void fillBoth(final Read first, final Read second,
			final int averagePairDistance, final boolean requireCorrectStrands,
			final boolean sameStrandPairs, final float[] firstVector,
			final float[] secondVector, final Scratch scratch){
		fillBoth(first, second, averagePairDistance, requireCorrectStrands, sameStrandPairs,
				firstVector, secondVector, scratch, false);
	}

	/** V1 inference variant of fillBoth; zeros omitted composition fields in both end blocks.
	 * Retains fillBoth's validation and leaves an unmapped end's output untouched. */
	static void fillBothRuntime(final Read first, final Read second,
			final int averagePairDistance, final boolean requireCorrectStrands,
			final boolean sameStrandPairs, final float[] firstVector,
			final float[] secondVector, final Scratch scratch){
		fillBoth(first, second, averagePairDistance, requireCorrectStrands, sameStrandPairs,
				firstVector, secondVector, scratch, true);
	}

	private static void fillBoth(final Read first, final Read second,
			final int averagePairDistance, final boolean requireCorrectStrands,
			final boolean sameStrandPairs, final float[] firstVector,
			final float[] secondVector, final Scratch scratch, final boolean runtime){
		if(first==null || second==null || first.mate!=second || second.mate!=first){
			throw new IllegalArgumentException("Paired neural MAPQ requires reciprocal mates");
		}
		if(!first.mapped() && !second.mapped()){
			throw new IllegalArgumentException("Paired neural MAPQ requires at least one mapped mate");
		}
		if(averagePairDistance<0){throw new IllegalArgumentException("Average pair distance must be nonnegative: "+averagePairDistance);}
		validateVector(firstVector, scratch); validateVector(secondVector, scratch);
		if(first.mapped() && second.mapped() && firstVector==secondVector){
			throw new IllegalArgumentException("Mapped mates need distinct neural MAPQ output vectors; the second orientation would overwrite the first");
		}
		if(first.mapped()){
			if(!first.primary()){throw new IllegalArgumentException("Paired neural MAPQ first mate must be primary");}
			fillEnd(first, scratch.anchorRaw, scratch.anchorScratch, runtime);
		}else{Arrays.fill(scratch.anchorRaw, 0);}
		if(second.mapped()){
			if(!second.primary()){throw new IllegalArgumentException("Paired neural MAPQ second mate must be primary");}
			fillEnd(second, scratch.mateRaw, scratch.mateScratch, runtime);
		}else{Arrays.fill(scratch.mateRaw, 0);}
		if(first.mapped()){
			fillPrepared(first, second, averagePairDistance, requireCorrectStrands, sameStrandPairs,
					scratch.anchorRaw, scratch.mateRaw, firstVector);
		}
		if(second.mapped()){
			fillPrepared(second, first, averagePairDistance, requireCorrectStrands, sameStrandPairs,
					scratch.mateRaw, scratch.anchorRaw, secondVector);
		}
	}

	/** Full exports retain composition; the frozen V1 inference path does not use it. */
	private static void fillEnd(final Read read, final float[] raw,
			final NeuralMapqFeatureExtractor.Scratch scratch, final boolean runtime){
		if(runtime){NeuralMapqFeatureExtractor.fillPairedEndRuntime(read, raw, scratch);}
		else{NeuralMapqFeatureExtractor.fillPairedEnd(read, raw, scratch);}
	}

	private static void validatePair(final Read anchor, final Read mate,
			final int averagePairDistance, final float[] vector, final Scratch scratch){
		if(anchor==null || mate==null || anchor.mate!=mate || mate.mate!=anchor){
			throw new IllegalArgumentException("Paired neural MAPQ requires reciprocal mates");
		}
		if(!anchor.mapped() || !anchor.primary()){
			throw new IllegalArgumentException("Paired neural MAPQ anchor must be mapped and primary");
		}
		if(averagePairDistance<0){
			throw new IllegalArgumentException("Average pair distance must be nonnegative: "+averagePairDistance);
		}
		validateVector(vector, scratch);
	}

	private static void validateVector(final float[] vector, final Scratch scratch){
		if(vector==null || vector.length!=NeuralMapqPairedFeatureSchema.WIDTH || scratch==null){
			throw new IllegalArgumentException("Paired neural MAPQ vector/scratch differs from the pilot schema");
		}
	}

	/** Copies prepared end blocks and appends pair fields in schema order. */
	private static void fillPrepared(final Read anchor, final Read mate,
			final int averagePairDistance, final boolean requireCorrectStrands,
			final boolean sameStrandPairs, final float[] anchorRaw,
			final float[] mateRaw, final float[] vector){
		System.arraycopy(anchorRaw, 0, vector, 0, anchorRaw.length);
		System.arraycopy(mateRaw, 0, vector, anchorRaw.length, mateRaw.length);
		final boolean mateMapped=mate.mapped();
		final boolean sameChrom=mateMapped && anchor.chrom==mate.chrom;
		final int left=mateMapped ? Math.min(anchor.start, mate.start) : 0;
		final int right=mateMapped ? Math.max(anchor.stop, mate.stop) : 0;
		final boolean sameScaffold=sameChrom && Data.isSingleScaffold(anchor.chrom, left, right);
		final boolean sameStrand=mateMapped && anchor.strand()==mate.strand();
		// This feature tests strand relation only, not left/right order or concordance.
		final boolean expectedOrientation=mateMapped && (sameStrand==sameStrandPairs);
		final int observedInsert=observedInsert(anchor, mate, requireCorrectStrands, sameStrandPairs);
		final boolean insertMissing=observedInsert<1;
		final int storedInsert=anchor.insert();
		//Widen before arithmetic: the public API accepts any nonnegative int distance.
		final long expectedFragment=(long)averagePairDistance+anchor.length()+mate.length();
		final int innerDistance=(sameChrom ? innerDistance(anchor, mate, requireCorrectStrands) : 0);
		final long signedDeviation=(sameChrom ? (long)innerDistance-averagePairDistance : 0);
		final SiteScore topA=anchor.topSite();
		final SiteScore topM=mateMapped ? mate.topSite() : null;
		final int pairedA=topA==null ? 0 : topA.pairedScore;
		final int pairedM=topM==null ? 0 : topM.pairedScore;
		final int slowA=topA==null ? 0 : topA.slowScore;
		final int slowM=topM==null ? 0 : topM.slowScore;

		int i=NeuralMapqFeatureSchema.WIDTH*2;
		vector[i++]=anchor.pairnum();
		vector[i++]=mateMapped ? 1 : 0;
		vector[i++]=anchor.paired() ? 1 : 0;
		vector[i++]=mate.paired() ? 1 : 0;
		vector[i++]=sameChrom ? 1 : 0;
		vector[i++]=sameScaffold ? 1 : 0;
		vector[i++]=sameStrand ? 1 : 0;
		vector[i++]=expectedOrientation ? 1 : 0;
		vector[i++]=anchor.insertvalid() ? 1 : 0;
		vector[i++]=insertMissing ? 1 : 0;
		vector[i++]=insertMissing ? 0 : observedInsert;
		vector[i++]=storedInsert<0 ? 0 : storedInsert;
		vector[i++]=averagePairDistance;
		vector[i++]=expectedFragment;
		vector[i++]=innerDistance;
		vector[i++]=signedDeviation;
		vector[i++]=Math.abs(signedDeviation);
		vector[i++]=pairedA;
		vector[i++]=pairedM;
		vector[i++]=Math.max(0L, (long)pairedA-slowA);
		vector[i++]=Math.max(0L, (long)pairedM-slowM);
		vector[i++]=(topA==null ? 0L : topA.score)+(topM==null ? 0L : topM.score);
		vector[i++]=(long)anchor.length()+mate.length();
		vector[i++]=(anchor.perfect() && mate.perfect()) ? 1 : 0;
		if(i!=vector.length){throw new AssertionError("Paired extractor wrote "+i+" of "+vector.length);}
		NeuralMapqPairedFeatureSchema.validateVector(vector);
	}

	/** Uses Read's insert-size convention, including its cross-chromosome/missing zero. */
	static int observedInsert(final Read anchor, final Read mate,
			final boolean requireCorrectStrands, final boolean sameStrandPairs){
		if(anchor==null || mate==null || !anchor.mapped() || !mate.mapped()){return 0;}
		return Read.insertSizeMapped(anchor, mate, sameStrandPairs || !requireCorrectStrands);
	}

	/** Signed start-minus-stop distance; adjacent ends yield 1, overlap yields zero or less.
	 * This is not the count of intervening bases (which subtracts another one).
	 * Opposite-strand required pairs use plus-to-minus order; others use left-to-right. */
	static int innerDistance(final Read a, final Read b, final boolean requireCorrectStrands){
		if(a==null || b==null || !a.mapped() || !b.mapped() || a.chrom!=b.chrom){return 0;}
		if(requireCorrectStrands && a.strand()!=b.strand()){
			return a.strand()==0 ? b.start-a.stop : a.start-b.stop;
		}
		return a.start<=b.start ? b.start-a.stop : a.start-b.stop;
	}

	/** Worker-local raw end blocks and counters, reused by either extraction entry point. */
	public static final class Scratch{
		final float[] anchorRaw=new float[NeuralMapqFeatureSchema.WIDTH];
		final float[] mateRaw=new float[NeuralMapqFeatureSchema.WIDTH];
		final NeuralMapqFeatureExtractor.Scratch anchorScratch=new NeuralMapqFeatureExtractor.Scratch();
		final NeuralMapqFeatureExtractor.Scratch mateScratch=new NeuralMapqFeatureExtractor.Scratch();
	}
}
