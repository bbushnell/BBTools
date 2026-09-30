package align2;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import stream.Read;
import stream.SiteScore;

/**
 * Extracts the 42 raw neural-MAPQ fields shared by training and inference.
 * The caller supplies a mapped primary read, its existing traceback, and retained
 * sites with the selected primary first and competitors in descending score order.
 * This class neither aligns the read nor sorts or modifies its candidate sites.
 * Network input scaling and field selection belong to NeuralMapqFeatureTransform.
 *
 * @author Collei
 */
public final class NeuralMapqFeatureExtractor{

	private NeuralMapqFeatureExtractor(){}

	/**
	 * Fills all 42 fields for an unpaired read, including composition statistics.
	 * The caller owns one reusable Scratch and vector per mapping thread; successful
	 * extraction uses primitive counters without allocating a feature object.
	 * No truth/oracle labels are read. Every vector position is overwritten.
	 */
	public static void fill(final Read read, final float[] vector, final Scratch scratch){
		fill(read, vector, scratch, true, true);
	}

	/** Fills a 42-field raw vector for an unpaired read. Fields 37..41 become zero;
	 * the separate transform produces the accepted 37 network inputs. */
	public static void fillRuntime(final Read read, final float[] vector, final Scratch scratch){
		fill(read, vector, scratch, false, true);
	}

	/** Fills all 42 fields for one mapped primary end, allowing a mate pointer.
	 * Reciprocal-pair validation is the paired extractor's responsibility. */
	static void fillPairedEnd(final Read read, final float[] vector, final Scratch scratch){
		fill(read, vector, scratch, true, false);
	}

	/** Paired inference omits composition fields 37..41, which V1 does not consume. */
	static void fillPairedEndRuntime(final Read read, final float[] vector, final Scratch scratch){
		fill(read, vector, scratch, false, false);
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
		// The retained list defines the competitor population, not every discovered site.
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
		vector[i++]=second==null ? 0 : (long)topScore-second.score;
		vector[i++]=third==null ? 0 : (long)topScore-third.score;
		vector[i++]=second==null ? 0 : ratio(second.score, topScore);
		vector[i++]=third==null ? 0 : ratio(third.score, topScore);
		vector[i++]=top.hits;
		// Preserve Read's configured flat/skewed identity definition for training parity.
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
		// Read's missing-quality defaults are mean=40, minimum=41, expected errors=0.
		vector[i++]=(float)read.avgQualityByScoreDouble(length);
		vector[i++]=read.minQuality();
		vector[i++]=read.expectedErrors(true, length);
		if(includeDroppedComposition){
			vector[i++]=scratch.readNBases*inverseLength;
			vector[i++]=scratch.definedBases<1 ? 0 :
					(scratch.monomers[1]+scratch.monomers[2])/(float)scratch.definedBases;
			vector[i++]=read.longestHomopolymer()*inverseLength;
			vector[i++]=entropy(scratch.monomers, scratch.definedBases, 2.0);
			vector[i++]=entropy(scratch.dimers, scratch.definedDimers, 4.0);
		}else{
			while(i<NeuralMapqFeatureSchema.WIDTH){vector[i++]=0;}
		}
		assert(i==NeuralMapqFeatureSchema.WIDTH) : "Extractor wrote "+i+
				" values for schema width "+NeuralMapqFeatureSchema.WIDTH;
		NeuralMapqFeatureSchema.validateVector(vector);
	}

	/** Missing zero-score denominator yields zero; signed competitor scores are retained. */
	private static float ratio(final int numerator, final int denominator){
		return denominator==0 ? 0 : numerator/(float)denominator;
	}

	/** Shannon entropy normalized by the alphabet's maximum bits (2 for bases, 4 for dimers). */
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

	/** Counts defined bases and adjacent defined dimers; ambiguous bases break dimer runs. */
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

	/** Reads BBMap long matches or symbol-then-count RLE (m4S2), never SAM CIGAR.
	 * Adjacent tokens of the same event class are merged by Scratch.accept.
	 * Runs must be positive and expanded length must fit the int counters used
	 * here and by Read.identity. Read.toShortMatchString emits this format. */
	private static void parseMatch(final byte[] match, final Scratch scratch){
		if(match.length==0){throw new IllegalArgumentException("Neural MAPQ requires a nonempty match string");}
		byte mode=0;
		int count=0;
		boolean haveMode=false, haveCount=false;
		for(final byte symbol : match){
			if(symbol>='0' && symbol<='9'){
				if(!haveMode){throw new IllegalArgumentException("MAPQ match run length precedes its symbol");}
				final int digit=symbol-'0';
				if(count>(Integer.MAX_VALUE-digit)/10){
					throw new IllegalArgumentException("MAPQ match run exceeds int range for symbol "+(char)mode);
				}
				count=count*10+digit;
				haveCount=true;
			}else{
				if(haveMode){scratch.accept(mode, haveCount ? count : 1);}
				mode=symbol;
				count=0;
				haveMode=true;
				haveCount=false;
			}
		}
		if(haveMode){scratch.accept(mode, haveCount ? count : 1);}
	}

	/** Reusable primitive counters; one instance belongs to each mapper thread. */
	public static final class Scratch{
		/** Resets all counts and run continuity before extracting another read. */
		public void clear(){
			Arrays.fill(monomers, 0);
			Arrays.fill(dimers, 0);
			substitutionEvents=substitutedBases=insertionEvents=insertedBases=0;
			deletionEvents=deletedBases=longestIndel=clippedBases=alignmentNBases=0;
			readNBases=definedBases=definedDimers=0;
			previousClass=0;
			currentIndelLength=0;
			expandedLength=0;
		}

		/** Adds a run, counting contiguous event classes rather than encoded tokens. */
		private void accept(final byte symbol, final int length){
			if(length<1){throw new IllegalArgumentException("Nonpositive match run: "+length);}
			final byte eventClass=eventClass(symbol);
			//All event totals and Read.identity's counts are bounded by expanded length.
			if(length>Integer.MAX_VALUE-expandedLength){
				throw new IllegalArgumentException("Expanded MAPQ match exceeds int counter range: prior="+
						expandedLength+", next run="+length);
			}
			expandedLength+=length;
			if(eventClass==SUBSTITUTION){
				substitutedBases+=length;
				if(previousClass!=eventClass){substitutionEvents++;}
			}else if(eventClass==INSERTION){
				insertedBases+=length;
				if(previousClass!=eventClass){insertionEvents++; currentIndelLength=0;}
				currentIndelLength+=length;
				longestIndel=Math.max(longestIndel, currentIndelLength);
			}else if(eventClass==DELETION){
				deletedBases+=length;
				if(previousClass!=eventClass){deletionEvents++; currentIndelLength=0;}
				currentIndelLength+=length;
				longestIndel=Math.max(longestIndel, currentIndelLength);
			}else if(eventClass==CLIP){clippedBases+=length; currentIndelLength=0;}
			else if(eventClass==NOCALL){alignmentNBases+=length; currentIndelLength=0;}
			else{currentIndelLength=0;}
			previousClass=eventClass;
		}

		/** Collapses BBMap traceback symbols into the pilot's six event classes. */
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
		int substitutionEvents, substitutedBases, insertionEvents, insertedBases;
		int deletionEvents, deletedBases, longestIndel, clippedBases, alignmentNBases;
		int readNBases, definedBases, definedDimers;
		private byte previousClass;
		private int currentIndelLength;
		private int expandedLength;
	}

	private static final byte MATCH=1, SUBSTITUTION=2, INSERTION=3, DELETION=4, CLIP=5, NOCALL=6;
}
