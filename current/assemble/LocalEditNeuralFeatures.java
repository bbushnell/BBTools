package assemble;

import java.util.Arrays;

import structures.ByteBuilder;
import structures.IntList;
import ukmer.Kmer;

/** Shared 28-feature legacy and 39-feature read-depth neural contracts.
 * Candidate order is the three substitutions in A/C/G/T order excluding the
 * observed base, deletion, then insertions A/C/G/T.
 * @author Fischl */
public final class LocalEditNeuralFeatures {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private LocalEditNeuralFeatures(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Fill unrounded feature values. Returns false when evidence is unavailable. */
	public static boolean fill(final double[] out, final int currentDepth, final int[] candidateDepths,
			final int candidate, final int originalBase, final int leftFlankDepth, final int rightFlankDepth,
			final int runLength, final boolean lowComplexity, final int readLength, final int coordinate){
		if(out==null || out.length<FEATURE_COUNT || candidateDepths==null || candidateDepths.length!=CANDIDATE_COUNT ||
				candidate<0 || candidate>=CANDIDATE_COUNT || originalBase<0 || originalBase>3 || currentDepth<0 ||
				(leftFlankDepth<0 && rightFlankDepth<0) || runLength<0 || readLength<1 || coordinate<0 || coordinate>=readLength){
			return false;
		}
		for(final int depth:candidateDepths){if(depth<0){return false;}}
		int next=0;
		out[next++]=Math.log1p(currentDepth);
		for(final int depth:candidateDepths){out[next++]=Math.log1p(depth);}
		final double normal=localDepth(leftFlankDepth, rightFlankDepth);
		out[next++]=Math.log1p(normal);
		final int depth=candidateDepths[candidate];
		out[next++]=Math.log((depth+1.0)/(currentDepth+1.0));
		out[next++]=Math.log((depth+1.0)/(normal+1.0));
		out[next++]=Math.log((depth+1.0)/(bestOther(candidateDepths, candidate)+1.0));
		final int operation=operation(candidate), base=candidateBase(candidate, originalBase);
		for(int op=SUBSTITUTION; op<=INSERTION; op++){out[next++]=operation==op ? 1 : 0;}
		for(int b=0; b<4; b++){out[next++]=base==b ? 1 : 0;}
		for(int b=0; b<4; b++){out[next++]=originalBase==b ? 1 : 0;}
		out[next++]=Math.log1p(runLength);
		out[next++]=lowComplexity ? 1 : 0;
		out[next++]=Math.log1p(readLength)/10.0;
		out[next++]=coordinate/(double)Math.max(1, readLength-1);
		assert(next==FEATURE_COUNT) : "Feature contract must emit exactly "+FEATURE_COUNT+" values.";
		return true;
	}

	/** Append P0,P10,...,P100 from the original read, scaled by log2(depth)/8.
	 * The first 28 values retain their legacy order and scaling. */
	public static boolean appendReadDepths(final double[] out, final int[] percentiles){
		if(out==null || out.length!=DEPTH_FEATURE_COUNT || percentiles==null || percentiles.length!=PERCENTILE_COUNT){return false;}
		for(int i=0; i<PERCENTILE_COUNT; i++){
			if(percentiles[i]<0 || (i>0 && percentiles[i]<percentiles[i-1])){return false;}
			out[FEATURE_COUNT+i]=Math.log(Math.max(1, percentiles[i]))*DEPTH_SCALE;
		}
		return true;
	}

	/** Count every valid original-read kmer once; omit N-containing windows, not
	 * measured zero depths. Sort reusable scratch and select nearest-rank deciles.
	 * Returns the number of count lookups; zero means no usable read context. */
	static int readDepthPercentiles(final byte[] bases, final Kmer key,
			final HomopolymerIndelProposal.CountLookup lookup, final IntList scratch, final int[] out){
		if(bases==null || key==null || lookup==null || scratch==null || out==null || out.length!=PERCENTILE_COUNT){
			throw new IllegalArgumentException("Read-depth percentiles require bases, count lookup, scratch and eleven output cells.");
		}
		scratch.clear();
		key.clearFast();
		for(final byte base:bases){
			final int numeric=baseIndex(base);
			if(numeric<0){key.clearFast(); continue;}
			key.addRightNumeric(numeric);
			if(key.len()<key.kbig){continue;}
			final int depth=lookup.count(key);
			if(depth<-1){throw new IllegalStateException("Invalid read-depth count: "+depth);}
			scratch.add(Math.max(0, depth));
		}
		selectDepthPercentiles(scratch, out);
		return scratch.size;
	}

	/** Nearest rank: P0=min; P(10*i)=sorted[ceil(i*n/10)-1] for i>0.
	 * Empty input yields zeros, but callers must abstain when no valid kmers exist. */
	static void selectDepthPercentiles(final IntList depths, final int[] out){
		assert(depths!=null && out!=null && out.length==PERCENTILE_COUNT) :
				"The neural contract appends exactly eleven read-depth deciles.";
		if(depths.size==0){Arrays.fill(out, 0); return;}
		Arrays.sort(depths.array, 0, depths.size);
		if(depths.array[0]<0){throw new IllegalArgumentException("Only measured nonnegative depths enter read percentiles.");}
		out[0]=depths.array[0];
		for(int i=1; i<PERCENTILE_COUNT; i++){
			final int index=(int)(((long)i*depths.size+9)/10)-1;
			out[i]=depths.array[index];
		}
		assert(out[10]==depths.array[depths.size-1]) : "P100 must retain the maximum observed depth.";
	}

	/** Serialize the same fixed-eight-decimal representation used to train the frozen model. */
	public static void append(final ByteBuilder bb, final double[] features, final int label){
		if(bb==null || features==null || !supportedCount(features.length) || (label!=0 && label!=1)){
			throw new IllegalArgumentException("Feature serialization requires 28 or 39 values and a binary label.");
		}
		for(int i=0; i<features.length; i++){
			if(!Double.isFinite(features[i])){throw new IllegalArgumentException("Nonfinite feature at index "+i+'.');}
			bb.append(features[i], DECIMALS).tab();
		}
		bb.append(label).nl();
	}

	/** Convert raw values exactly as BBTools' eight-decimal writer and parser do. */
	public static boolean toModelInput(final double[] features, final float[] input){
		if(features==null || !supportedCount(features.length) || input==null || input.length!=features.length){return false;}
		for(int i=0; i<features.length; i++){
			final double value=features[i];
			if(!Double.isFinite(value)){return false;}
			input[i]=quantizedFloat(value);
			if(!Float.isFinite(input[i])){return false;}
		}
		return true;
	}

	/** Only the legacy and depth-extended layouts have defined feature semantics. */
	public static boolean supportedCount(final int count){
		return count==FEATURE_COUNT || count==DEPTH_FEATURE_COUNT;
	}

	/** Map an operation/base pair to the frozen candidate order. */
	public static int candidateIndex(final int originalBase, final int operation, final int base){
		if(originalBase<0 || originalBase>3){return -1;}
		if(operation==DELETION){return base<0 ? 3 : -1;}
		if(operation==INSERTION){return base>=0 && base<4 ? 4+base : -1;}
		if(operation!=SUBSTITUTION || base<0 || base>3 || base==originalBase){return -1;}
		int index=0;
		for(int b=0; b<4; b++){if(b!=originalBase){if(b==base){return index;}index++;}}
		return -1;
	}

	/** Returns the operation encoded by a frozen candidate index, or -1 when invalid. */
	public static int operation(final int candidate){
		return candidate<0 || candidate>=CANDIDATE_COUNT ? -1 : candidate<3 ? SUBSTITUTION : candidate==3 ? DELETION : INSERTION;
	}

	/** Returns the proposed base in numeric A/C/G/T order, or -1 for deletion. */
	public static int candidateBase(final int candidate, final int originalBase){
		if(candidate<0 || candidate>=CANDIDATE_COUNT || originalBase<0 || originalBase>3){return -1;}
		if(candidate==3){return -1;}
		if(candidate>3){return candidate-4;}
		int index=0;
		for(int b=0; b<4; b++){if(b!=originalBase){if(index++==candidate){return b;}}}
		return -1;
	}

	/** Returns the homopolymer run length after applying a substitution. */
	public static int substitutionRun(final byte[] bases, final int position, final int base){
		return bases==null || position<0 || position>=bases.length || base<0 || base>3 ? -1 :
				1+runLeft(bases, position-1, base)+runRight(bases, position+1, base);
	}

	/** Returns the homopolymer run length after applying an insertion. */
	public static int insertionRun(final byte[] bases, final int position, final int base){
		return bases==null || position<0 || position>=bases.length || base<0 || base>3 ? -1 :
				1+runLeft(bases, position-1, base)+runRight(bases, position, base);
	}

	/** Returns the joined neighboring run length after applying a deletion. */
	public static int deletionRun(final byte[] bases, final int position){
		if(bases==null || position<1 || position+1>=bases.length){return 0;}
		final int left=baseIndex(bases[position-1]), right=baseIndex(bases[position+1]);
		return left>=0 && left==right ? runLeft(bases, position-1, left)+runRight(bases, position+1, left) : 0;
	}

	/** True when the observed window is periodic with period 1, 2, or 3 away from the proposed locus. */
	public static boolean lowComplexity(final byte[] bases, final int start, final int k, final int position){
		if(bases==null || start<0 || k<1 || (long)start+k>bases.length || position<start || position>=start+k){return false;}
		for(int period=1; period<=3; period++){
			boolean periodic=true;
			for(int i=start+period; i<start+k; i++){
				if(i==position || i-period==position){continue;}
				final int a=baseIndex(bases[i]), b=baseIndex(bases[i-period]);
				if(a<0 || b<0 || a!=b){periodic=false; break;}
			}
			if(periodic){return true;}
		}
		return false;
	}

	/** Maps an ASCII nucleotide to numeric A/C/G/T order, or -1 if unsupported. */
	public static int baseIndex(final byte base){
		final byte b=base>='a' && base<='z' ? (byte)(base-32) : base;
		return b=='A' ? 0 : b=='C' ? 1 : b=='G' ? 2 : b=='T' ? 3 : -1;
	}

	private static int runLeft(final byte[] bases, int position, final int base){
		int length=0;
		while(position>=0 && baseIndex(bases[position])==base){length++; position--;}
		return length;
	}

	private static int runRight(final byte[] bases, int position, final int base){
		int length=0;
		while(position<bases.length && baseIndex(bases[position])==base){length++; position++;}
		return length;
	}

	private static double localDepth(final int left, final int right){
		if(left<0){return right;}
		if(right<0){return left;}
		return ((double)left+right)*0.5;
	}

	private static int bestOther(final int[] depths, final int excluded){
		int best=0;
		for(int i=0; i<depths.length; i++){if(i!=excluded){best=Math.max(best, depths[i]);}}
		return best;
	}

	private static float quantizedFloat(final double value){
		if(value==(long)value){return (float)value;}
		final boolean negative=value<0;
		double magnitude=negative ? -value : value;
		magnitude+=0.5*DECIMAL_SCALE_INVERSE;
		final long upper=(long)magnitude;
		final long lower=(long)((magnitude-upper)*DECIMAL_SCALE);
		final double parsed=upper+lower*DECIMAL_SCALE_INVERSE;
		return (float)(negative ? -parsed : parsed);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	public static final int FEATURE_COUNT=28, CANDIDATE_COUNT=8;
	public static final int PERCENTILE_COUNT=11, DEPTH_FEATURE_COUNT=FEATURE_COUNT+PERCENTILE_COUNT;
	private static final double DEPTH_SCALE=0.125/Math.log(2);
	public static final int SUBSTITUTION=1, DELETION=2, INSERTION=3;
	private static final int DECIMALS=8;
	private static final double DECIMAL_SCALE=100000000.0, DECIMAL_SCALE_INVERSE=0.00000001;
}
