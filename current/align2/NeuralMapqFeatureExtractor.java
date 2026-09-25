package align2;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import stream.Read;
import stream.SiteScore;

/** Fills the ordered single-end neural-MAPQ pilot vector. */
public final class NeuralMapqFeatureExtractor {

	private NeuralMapqFeatureExtractor(){}

	/**
	 * Fills {@code vector} without allocating.  The caller owns one reusable
	 * Scratch and vector per mapping thread.
	 */
	public static void fill(final Read read, final float[] vector, final Scratch scratch){
		fill(read,vector,scratch,true,true);
	}

	/** Fills only the accepted 37-input runtime boundary; dropped composition fields become zero. */
	public static void fillRuntime(final Read read, final float[] vector, final Scratch scratch){
		fill(read,vector,scratch,false,true);
	}

	/** Raw single-end feature block for one end of a paired example. */
	static void fillPairedEnd(final Read read, final float[] vector, final Scratch scratch){
		fill(read,vector,scratch,true,false);
	}

	private static void fill(final Read read, final float[] vector, final Scratch scratch,
			final boolean includeDroppedComposition, final boolean requireUnpaired){
		if(read==null || vector==null || scratch==null){
			throw new IllegalArgumentException("Neural MAPQ extraction requires read, vector, and scratch");
		}
		if(requireUnpaired && read.mate!=null){
			throw new IllegalArgumentException("Single-end neural MAPQ V1 cannot extract paired reads");
		}
		if(!read.mapped() || !read.primary() || read.bases==null || read.bases.length<1){
			throw new IllegalArgumentException("Neural MAPQ extraction requires a mapped primary read with bases");
		}
		if(read.match==null){
			throw new IllegalArgumentException("Neural MAPQ extraction requires the existing primary traceback");
		}
		final ArrayList<SiteScore> sites=read.sites;
		if(sites==null || sites.isEmpty()){
			throw new IllegalArgumentException("Neural MAPQ extraction requires the retained primary SiteScore");
		}
		if(vector.length!=NeuralMapqFeatureSchema.WIDTH){
			throw new IllegalArgumentException("Neural MAPQ vector width "+vector.length+
					" differs from schema width "+NeuralMapqFeatureSchema.WIDTH);
		}

		final SiteScore top=sites.get(0);
		final SiteScore second=sites.size()>1 ? sites.get(1) : null;
		final SiteScore third=sites.size()>2 ? sites.get(2) : null;
		final int length=read.length();
		final float inverseLength=1f/length;
		final int topScore=top.score;
		int equalTopCount=1;
		while(equalTopCount<sites.size() && sites.get(equalTopCount).score==topScore){equalTopCount++;}

		scratch.clear();
		parseMatch(read.match, scratch);
		if(includeDroppedComposition){parseRead(read.bases, scratch);}

		int i=0;
		vector[i++]=length;
		vector[i++]=read.rescued() ? 1 : 0;
		vector[i++]=read.perfect() ? 1 : 0;
		vector[i++]=top.semiperfect ? 1 : 0;
		vector[i++]=read.ambiguous() ? 1 : 0;
		vector[i++]=read.mapScore;
		vector[i++]=read.mapScore*inverseLength;
		vector[i++]=top.quickScore;
		vector[i++]=top.slowScore;
		vector[i++]=topScore;
		vector[i++]=topScore/(100f*length);
		vector[i++]=sites.size();
		vector[i++]=equalTopCount;
		vector[i++]=second==null ? 1 : 0;
		vector[i++]=third==null ? 1 : 0;
		vector[i++]=second==null ? 0 : second.score;
		vector[i++]=third==null ? 0 : third.score;
		vector[i++]=second==null ? 0 : topScore-second.score;
		vector[i++]=third==null ? 0 : topScore-third.score;
		vector[i++]=second==null ? 0 : ratio(second.score,topScore);
		vector[i++]=third==null ? 0 : ratio(third.score,topScore);
		vector[i++]=top.hits;
		vector[i++]=Read.identity(read.match);
		vector[i++]=scratch.substitutionEvents;
		vector[i++]=scratch.substitutedBases;
		vector[i++]=scratch.insertionEvents;
		vector[i++]=scratch.insertedBases;
		vector[i++]=scratch.deletionEvents;
		vector[i++]=scratch.deletedBases;
		vector[i++]=scratch.substitutedBases+scratch.insertedBases+scratch.deletedBases;
		vector[i++]=scratch.longestIndel;
		vector[i++]=scratch.insertionEvents+scratch.deletionEvents;
		vector[i++]=scratch.clippedBases;
		vector[i++]=scratch.alignmentNBases;
		vector[i++]=(float)read.avgQualityByScoreDouble(length);
		vector[i++]=read.minQuality();
		vector[i++]=read.expectedErrors(true,length);
		if(includeDroppedComposition){
			vector[i++]=scratch.readNBases*inverseLength;
			vector[i++]=scratch.definedBases<1 ? 0 :
					(scratch.monomers[1]+scratch.monomers[2])/(float)scratch.definedBases;
			vector[i++]=read.longestHomopolymer()*inverseLength;
			vector[i++]=entropy(scratch.monomers,scratch.definedBases,2.0);
			vector[i++]=entropy(scratch.dimers,scratch.definedDimers,4.0);
		}else{
			while(i<NeuralMapqFeatureSchema.WIDTH){vector[i++]=0;}
		}
		assert(i==NeuralMapqFeatureSchema.WIDTH) : "Extractor wrote "+i+
				" values for schema width "+NeuralMapqFeatureSchema.WIDTH;
		NeuralMapqFeatureSchema.validateVector(vector);
	}

	private static float ratio(final int numerator, final int denominator){
		return denominator==0 ? 0 : numerator/(float)denominator;
	}

	private static float entropy(final int[] counts, final int total, final double maximumBits){
		if(total<2){return 0;}
		double entropy=0;
		for(final int count : counts){
			if(count>0){
				final double probability=count/(double)total;
				entropy-=probability*(Math.log(probability)/Math.log(2));
			}
		}
		return (float)(entropy/maximumBits);
	}

	private static void parseRead(final byte[] bases, final Scratch scratch){
		int previous=-1;
		for(final byte base : bases){
			final int number=AminoAcid.baseToNumber[base];
			if(number<0){
				scratch.readNBases++;
				previous=-1;
			}else{
				scratch.monomers[number]++;
				scratch.definedBases++;
				if(previous>=0){
					scratch.dimers[(previous<<2)|number]++;
					scratch.definedDimers++;
				}
				previous=number;
			}
		}
	}

	private static void parseMatch(final byte[] match, final Scratch scratch){
		byte mode=0;
		int count=0;
		for(final byte symbol : match){
			if(symbol>='0' && symbol<='9'){
				count=count*10+symbol-'0';
			}else{
				if(mode!=0){scratch.accept(mode,Math.max(1,count));}
				mode=symbol;
				count=0;
			}
		}
		if(mode!=0){scratch.accept(mode,Math.max(1,count));}
	}

	/** Reusable primitive counters; one instance belongs to each mapper thread. */
	public static final class Scratch {
		public void clear(){
			Arrays.fill(monomers,0);
			Arrays.fill(dimers,0);
			substitutionEvents=substitutedBases=insertionEvents=insertedBases=0;
			deletionEvents=deletedBases=longestIndel=clippedBases=alignmentNBases=0;
			readNBases=definedBases=definedDimers=0;
			previousClass=0;
			currentIndelLength=0;
		}

		private void accept(final byte symbol, final int length){
			if(length<1){throw new IllegalArgumentException("Nonpositive match run: "+length);}
			final byte eventClass=eventClass(symbol);
			if(eventClass==SUBSTITUTION){
				substitutedBases+=length;
				if(previousClass!=eventClass){substitutionEvents++;}
			}else if(eventClass==INSERTION){
				insertedBases+=length;
				if(previousClass!=eventClass){insertionEvents++;currentIndelLength=0;}
				currentIndelLength+=length;
				longestIndel=Math.max(longestIndel,currentIndelLength);
			}else if(eventClass==DELETION){
				deletedBases+=length;
				if(previousClass!=eventClass){deletionEvents++;currentIndelLength=0;}
				currentIndelLength+=length;
				longestIndel=Math.max(longestIndel,currentIndelLength);
			}else if(eventClass==CLIP){clippedBases+=length;currentIndelLength=0;}
			else if(eventClass==NOCALL){alignmentNBases+=length;currentIndelLength=0;}
			else{currentIndelLength=0;}
			previousClass=eventClass;
		}

		private static byte eventClass(final byte symbol){
			switch(symbol){
				case 'm': case 'M': return MATCH;
				case 'S': case 's': return SUBSTITUTION;
				case 'I': case 'i': case 'X': case 'Y': return INSERTION;
				case 'D': case 'd': return DELETION;
				case 'C': case 'V': return CLIP;
				case 'N': case 'R': case 'B': return NOCALL;
				default: throw new IllegalArgumentException("Unsupported BBMap match symbol: "+(char)symbol);
			}
		}

		final int[] monomers=new int[4];
		final int[] dimers=new int[16];
		int substitutionEvents,substitutedBases,insertionEvents,insertedBases;
		int deletionEvents,deletedBases,longestIndel,clippedBases,alignmentNBases;
		int readNBases,definedBases,definedDimers;
		private byte previousClass;
		private int currentIndelLength;
	}

	private static final byte MATCH=1,SUBSTITUTION=2,INSERTION=3,DELETION=4,CLIP=5,NOCALL=6;
}
