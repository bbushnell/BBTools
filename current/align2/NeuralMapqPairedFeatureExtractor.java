package align2;

import java.util.Arrays;

import dna.Data;
import stream.Read;
import stream.SiteScore;

/** Allocation-free raw feature extraction for one mapped anchor and its mate. */
public final class NeuralMapqPairedFeatureExtractor {

	private NeuralMapqPairedFeatureExtractor(){}

	public static void fill(final Read anchor, final Read mate,
			final int averagePairDistance, final boolean requireCorrectStrands,
			final boolean sameStrandPairs, final float[] vector, final Scratch scratch){
		validatePair(anchor,mate,averagePairDistance,vector,scratch);
		NeuralMapqFeatureExtractor.fillPairedEnd(anchor,scratch.anchorRaw,scratch.anchorScratch);
		if(mate.mapped()){
			NeuralMapqFeatureExtractor.fillPairedEnd(mate,scratch.mateRaw,scratch.mateScratch);
		}else{Arrays.fill(scratch.mateRaw,0);}
		fillPrepared(anchor,mate,averagePairDistance,requireCorrectStrands,sameStrandPairs,
				scratch.anchorRaw,scratch.mateRaw,vector);
	}

	/** Fills both anchor orientations while extracting each end's raw features once. */
	public static void fillBoth(final Read first,final Read second,
			final int averagePairDistance,final boolean requireCorrectStrands,
			final boolean sameStrandPairs,final float[] firstVector,
			final float[] secondVector,final Scratch scratch){
		if(first==null||second==null||first.mate!=second||second.mate!=first){
			throw new IllegalArgumentException("Paired neural MAPQ requires reciprocal mates");
		}
		if(!first.mapped()&&!second.mapped()){
			throw new IllegalArgumentException("Paired neural MAPQ requires at least one mapped mate");
		}
		if(averagePairDistance<0){throw new IllegalArgumentException("Average pair distance must be nonnegative: "+averagePairDistance);}
		validateVector(firstVector,scratch);validateVector(secondVector,scratch);
		if(first.mapped()){
			if(!first.primary()){throw new IllegalArgumentException("Paired neural MAPQ first mate must be primary");}
			NeuralMapqFeatureExtractor.fillPairedEnd(first,scratch.anchorRaw,scratch.anchorScratch);
		}else{Arrays.fill(scratch.anchorRaw,0);}
		if(second.mapped()){
			if(!second.primary()){throw new IllegalArgumentException("Paired neural MAPQ second mate must be primary");}
			NeuralMapqFeatureExtractor.fillPairedEnd(second,scratch.mateRaw,scratch.mateScratch);
		}else{Arrays.fill(scratch.mateRaw,0);}
		if(first.mapped()){
			fillPrepared(first,second,averagePairDistance,requireCorrectStrands,sameStrandPairs,
					scratch.anchorRaw,scratch.mateRaw,firstVector);
		}
		if(second.mapped()){
			fillPrepared(second,first,averagePairDistance,requireCorrectStrands,sameStrandPairs,
					scratch.mateRaw,scratch.anchorRaw,secondVector);
		}
	}

	private static void validatePair(final Read anchor,final Read mate,
			final int averagePairDistance,final float[] vector,final Scratch scratch){
		if(anchor==null || mate==null || anchor.mate!=mate || mate.mate!=anchor){
			throw new IllegalArgumentException("Paired neural MAPQ requires reciprocal mates");
		}
		if(!anchor.mapped() || !anchor.primary()){
			throw new IllegalArgumentException("Paired neural MAPQ anchor must be mapped and primary");
		}
		if(averagePairDistance<0){
			throw new IllegalArgumentException("Average pair distance must be nonnegative: "+averagePairDistance);
		}
		validateVector(vector,scratch);
	}

	private static void validateVector(final float[] vector,final Scratch scratch){
		if(vector==null || vector.length!=NeuralMapqPairedFeatureSchema.WIDTH || scratch==null){
			throw new IllegalArgumentException("Paired neural MAPQ vector/scratch differs from the pilot schema");
		}
	}

	private static void fillPrepared(final Read anchor,final Read mate,
			final int averagePairDistance,final boolean requireCorrectStrands,
			final boolean sameStrandPairs,final float[] anchorRaw,
			final float[] mateRaw,final float[] vector){
		System.arraycopy(anchorRaw,0,vector,0,anchorRaw.length);
		System.arraycopy(mateRaw,0,vector,anchorRaw.length,mateRaw.length);
		final boolean mateMapped=mate.mapped();
		final boolean sameChrom=mateMapped && anchor.chrom==mate.chrom;
		final int left=mateMapped ? Math.min(anchor.start,mate.start) : 0;
		final int right=mateMapped ? Math.max(anchor.stop,mate.stop) : 0;
		final boolean sameScaffold=sameChrom && Data.isSingleScaffold(anchor.chrom,left,right);
		final boolean sameStrand=mateMapped && anchor.strand()==mate.strand();
		final boolean expectedOrientation=mateMapped && (sameStrand==sameStrandPairs);
		final int observedInsert=observedInsert(anchor,mate,requireCorrectStrands,sameStrandPairs);
		final boolean insertMissing=observedInsert<1;
		final int storedInsert=anchor.insert();
		final int expectedFragment=averagePairDistance+anchor.length()+mate.length();
		final int innerDistance=(sameChrom ? innerDistance(anchor,mate,requireCorrectStrands) : 0);
		final int signedDeviation=(sameChrom ? innerDistance-averagePairDistance : 0);
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
		vector[i++]=Math.abs((long)signedDeviation);
		vector[i++]=pairedA;
		vector[i++]=pairedM;
		vector[i++]=Math.max(0,pairedA-slowA);
		vector[i++]=Math.max(0,pairedM-slowM);
		vector[i++]=(topA==null ? 0 : topA.score)+(topM==null ? 0 : topM.score);
		vector[i++]=anchor.length()+mate.length();
		vector[i++]=(anchor.perfect() && mate.perfect()) ? 1 : 0;
		if(i!=vector.length){throw new AssertionError("Paired extractor wrote "+i+" of "+vector.length);}
		NeuralMapqPairedFeatureSchema.validateVector(vector);
	}

	static int observedInsert(final Read anchor, final Read mate,
			final boolean requireCorrectStrands, final boolean sameStrandPairs){
		if(anchor==null || mate==null || !anchor.mapped() || !mate.mapped()){return 0;}
		return Read.insertSizeMapped(anchor,mate,sameStrandPairs || !requireCorrectStrands);
	}

	static int innerDistance(final Read a, final Read b, final boolean requireCorrectStrands){
		if(a==null || b==null || !a.mapped() || !b.mapped() || a.chrom!=b.chrom){return 0;}
		if(requireCorrectStrands && a.strand()!=b.strand()){
			return a.strand()==0 ? b.start-a.stop : a.start-b.stop;
		}
		return a.start<=b.start ? b.start-a.stop : a.start-b.stop;
	}

	public static final class Scratch {
		final float[] anchorRaw=new float[NeuralMapqFeatureSchema.WIDTH];
		final float[] mateRaw=new float[NeuralMapqFeatureSchema.WIDTH];
		final NeuralMapqFeatureExtractor.Scratch anchorScratch=new NeuralMapqFeatureExtractor.Scratch();
		final NeuralMapqFeatureExtractor.Scratch mateScratch=new NeuralMapqFeatureExtractor.Scratch();
	}
}
