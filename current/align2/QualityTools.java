package align2;

import java.util.Arrays;
import shared.Random;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;



/**
 * 
 *  @author Brian Bushnell
 *  @date Jul 17, 2011 12:04:06 PM
 */
public class QualityTools {
	
	/*-------------------- Main --------------------*/

	/** Program entry point for testing quality matrix operations.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		
		for(int i=0; i<MATRIX_SIZE; i++){
			for(int j=0; j<MATRIX_SIZE; j++){
				System.err.print((int)qualsToPhredSafe((byte)i, (byte)j)+",");
			}
			System.err.println();
		}
		
//		byte[] quals=new byte[] {15, 12, 20, 9, 10, 16, 14, 7, 11, 10, 10, 10, 10, 4, 4, 30, 30, 30, 30};
//		float[] probs=makeKeyProbs(quals, 4);
//		float[] probs2=makeKeyProbs(quals, 4);
//
//		int[] scores=makeKeyScores(quals, 4, 50, 50, null);
//
//		System.out.println(Arrays.toString(probs)+"\n");
//		System.out.println(Arrays.toString(probs2)+"\n");
//		System.out.println(Arrays.toString(scores)+"\n");
//
//		bench(50, 20000000);
//		bench2(50, 20000000);
//		bench(50, 20000000);
//		bench2(50, 20000000);
//		bench(50, 20000000);
//		bench2(50, 20000000);
		
//		System.out.println(1-((1-.1f)*(1-.1f)*(1-.1f)));
//		System.out.println("\n"+Arrays.toString(PROB));
//		System.out.println("\n"+Arrays.toString(INVERSE));
//		System.out.println("\n"+Arrays.toString(SUB_PROB));
//		System.out.println("\n"+Arrays.toString(SUB_INVERSE));
		
//		initializeq102matrix(null);
//		for(int a=0; a<42; a++){
//			for(int b=0; b<42; b++){
//				for(int c=0; c<42; c++){
//					System.out.println(a+"\t"+b+"\t"+c+"\t"+q3ProbMatrix[a][b][c]);
//				}
//			}
//		}
		
	}
	
	/*-------------------- Constructors --------------------*/

	/** Default constructor for QualityTools */
	public QualityTools(){}
	
	/*-------------------- Methods --------------------*/

	/*-------------------- Overridden Methods --------------------*/

	/*-------------------- Abstract Methods --------------------*/

	/*-------------------- Static Methods --------------------*/
	
	/**
	 * Benchmarks performance of makeKeyProbs method.
	 * @param length Length of quality array to generate
	 * @param rounds Number of benchmark iterations
	 */
	public static void bench(int length, int rounds){
		
		long time=System.nanoTime();
		
		byte[] qual=new byte[length];
		for(int i=0; i<qual.length; i++){
			qual[i]=(byte)(Math.random()*30+5);
		}
		for(int i=0; i<rounds; i++){
			float[] r=makeKeyProbs(qual, null, 8, false);
			if(r[r.length-1]>1 || r[r.length-1]<0){
				System.err.println("Ooops! "+Arrays.toString(r));
			}
		}
		
		time=System.nanoTime()-time;
		float seconds=(float)(time/1000000000d);
		System.out.println("Bench Time: "+Tools.format("%.3f",seconds)+" s");
	}
	
	/**
	 * Benchmarks performance of makeKeyProbs2 method.
	 * @param length Length of quality array to generate
	 * @param rounds Number of benchmark iterations
	 */
	public static void bench2(int length, int rounds){
		
		long time=System.nanoTime();
		
		byte[] qual=new byte[length];
		for(int i=0; i<qual.length; i++){
			qual[i]=(byte)(Math.random()*30+5);
		}
		for(int i=0; i<rounds; i++){
			float[] r=makeKeyProbs2(qual, 8);
			if(r[r.length-1]>1 || r[r.length-1]<0){
				System.err.println("Ooops! "+Arrays.toString(r));
			}
		}
		
		time=System.nanoTime()-time;
		float seconds=(float)(time/1000000000d);
		System.out.println("Bench2 Time: "+Tools.format("%.3f",seconds)+" s");
	}
	
	/**
	 * Creates key quality scores from quality and base arrays.
	 * Combines probability calculations with scoring transformation.
	 *
	 * @param qual Quality scores array
	 * @param bases Base sequence array
	 * @param keylen Length of each key
	 * @param range Score range for transformation
	 * @param baseScore Base score for calculations
	 * @param out Output array (created if null)
	 * @param useModulo Whether to use modulo filtering for large references
	 * @return Array of quality scores for each key position
	 */
	public static int[] makeKeyScores(byte[] qual, byte[] bases, int keylen, int range, int baseScore, int[] out, boolean useModulo){
		float[] probs=makeKeyProbs(qual, bases, keylen, useModulo);
		return makeKeyScores(probs, (qual.length-keylen+1), range, baseScore, out);
	}
	
	/**
	 * Converts probability array to score array using linear transformation.
	 *
	 * @param probs Error probability array
	 * @param numProbs Number of probabilities to process
	 * @param range Score range for transformation
	 * @param baseScore Base score for calculations
	 * @param out Output array (created if null)
	 * @return Array of transformed scores
	 */
	public static int[] makeKeyScores(float[] probs, int numProbs, int range, int baseScore, int[] out){
		if(out==null){out=new int[numProbs];}
//		assert(out.length==probs.length);
		assert(out.length>=numProbs);
		for(int i=0; i<numProbs; i++){
			out[i]=baseScore+(int)Math.round(range*(1-(probs[i])));
		}
		return out;
	}
	
	/**
	 * Creates integer score array from quality scores.
	 * Scales correct probabilities to maxScore range.
	 *
	 * @param qual Quality scores array
	 * @param maxScore Maximum score value
	 * @param out Output array (created if null)
	 * @return Integer score array
	 */
	public static int[] makeIntScoreArray(byte[] qual, int maxScore, int[] out){
		if(out==null){out=new int[qual.length];}
		assert(out.length==qual.length);
		for(int i=0; i<qual.length; i++){
			float probM=PROB_CORRECT[qual[i]];
			out[i]=(int)Math.round(maxScore*probM);
		}
		return out;
	}
	
	/**
	 * Creates byte score array from quality scores.
	 *
	 * @param qual Quality scores array (if null, uses fixed high quality)
	 * @param maxScore Maximum score value
	 * @param out Output array (created if null)
	 * @param negative If true, scores are negative (offset from maxScore)
	 * @return Byte score array
	 */
	public static byte[] makeByteScoreArray(byte[] qual, int maxScore, byte[] out, boolean negative){
		if(qual==null){return makeByteScoreArray(maxScore, out, negative);}
		if(out==null){out=new byte[qual.length];}
		assert(out.length==qual.length);
		for(int i=0; i<qual.length; i++){
			float probM=PROB_CORRECT[qual[i]];
			int x=(int)Math.round(maxScore*probM);
			assert(x>=Byte.MIN_VALUE && x<=Byte.MAX_VALUE);
			if(negative){
				x=x-maxScore;
				assert(x<=0);
			}else{
				assert(x>=0 && x<=maxScore);
			}
			out[i]=(byte)x;
		}
		return out;
	}
	
	/**
	 * Creates byte score array filled with zeros (high quality assumption).
	 *
	 * @param maxScore Maximum score value (unused in current implementation)
	 * @param out Output array to fill
	 * @param negative Negative scoring flag (unused in current implementation)
	 * @return Zero-filled byte array
	 */
	public static byte[] makeByteScoreArray(int maxScore, byte[] out, boolean negative){
		assert(out!=null);
//		for(int i=0; i<out.length; i++){
//			float probM=SUB_PROB[30];
//			int x=(int)Math.round(maxScore*probM);
//			assert(x>=Byte.MIN_VALUE && x<=Byte.MAX_VALUE);
//			if(negative){
//				x=x-maxScore;
//				assert(x<=0);
//			}else{
//				assert(x>=0 && x<=maxScore);
//			}
//			out[i]=(byte)x;
//		}
		Arrays.fill(out, (byte)0);
		return out;
	}
	
	/** Returns prob of error for each key */
	public static float[] makeKeyProbs(byte[] quality, byte[] bases, int keylen, boolean useModulo){
		return makeKeyProbs(quality, bases, keylen, null, useModulo);
	}
	
	/** Returns prob of error for each key */
	public static float[] makeKeyProbs(byte[] quality, byte[] bases, int keylen, float[] out, boolean useModulo){
		if(quality==null){return makeKeyProbs(bases, keylen, out, useModulo);}
		if(out==null){out=new float[quality.length-keylen+1];}
		assert(out.length>=quality.length-keylen+1) : quality.length+", "+keylen+", "+out.length;
//		assert(out.length==quality.length-keylen+1);
		float key1=1;
		
		int timeSinceZero=0;
		for(int i=0; i<keylen; i++){
//			byte q=(bases==null || bases[i]!='N' ? quality[i] : 0);
			byte q=quality[i];
			if(q>0){timeSinceZero++;}else{timeSinceZero=0;} //Tracks location of N's
			assert(q<PROB_CORRECT.length) : Arrays.toString(quality);
			float f=PROB_CORRECT[q];
			key1*=f;
		}
		out[0]=1-key1;
		if(timeSinceZero<keylen){out[0]=1;}
		
		for(int a=0, b=keylen; b<quality.length; a++, b++){
//			byte qa=(bases==null || bases[a]!='N' ? quality[a] : 0);
//			byte qb=(bases==null || bases[b]!='N' ? quality[b] : 0);
			byte qa=quality[a];
			byte qb=quality[b];
			if(qb>0){timeSinceZero++;}else{timeSinceZero=0;}
			float ipa=PROB_CORRECT_INVERSE[qa];
			float pb=PROB_CORRECT[qb];
			key1=key1*ipa*pb;
			out[a+1]=1-key1;
			if(timeSinceZero<keylen){out[a+1]=1;}
		}
		
		if(bases!=null){
			if(useModulo){//Rare case for large references
				final int shift=2*keylen;
				final int shift2=shift-2;
				final int mask=~((-1)<<shift);
				int kmer=0, rkmer=0;
				
				int len=0;
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final int x=AminoAcid.baseToNumber[b];
					final int x2=AminoAcid.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
					
					if(x<0){len=0; rkmer=0;}else{len++;}
					if(len>=keylen){
						if(kmer%IndexMaker4.MODULO!=0 && rkmer%IndexMaker4.MODULO!=0){
							out[i-keylen+1]=1f;
//							assert(false) : kmer;
						}
					}
				}
			}
		}
		
		return out;
	}
	
	/** Returns prob of error for each key */
	public static float[] makeKeyProbs(byte[] bases, int keylen, float[] out, boolean useModulo){
		assert(out!=null) : "Must provide array if no quality vector";
		Arrays.fill(out, 0);
		
		if(bases!=null){
			if(useModulo){//Rare case for large references
				final int shift=2*keylen;
				final int shift2=shift-2;
				final int mask=~((-1)<<shift);
				int kmer=0, rkmer=0;
				
				int len=0;
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final int x=AminoAcid.baseToNumber[b];
					final int x2=AminoAcid.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
					
					if(x<0){len=0; rkmer=0;}else{len++;}
					if(len>=keylen){
						if(kmer%IndexMaker4.MODULO!=0 && rkmer%IndexMaker4.MODULO!=0){
							out[i-keylen+1]=1f;
//							assert(false) : kmer;
						}
					}
				}
			}
		}
		return out;
	}
	
	/**
	 * Alternative k-mer probability calculation using two-pointer approach.
	 * Processes from both ends of array simultaneously for benchmarking.
	 *
	 * @param quality Quality scores array
	 * @param keylen Length of each k-mer key
	 * @return Error probability for each k-mer position
	 */
	public static float[] makeKeyProbs2(byte[] quality, int keylen){
		float[] out=new float[quality.length-keylen+1];
		
		final int mid=out.length/2;
		
		float key1=1;
		float key2=1;
		for(int i=0, j=mid; i<keylen; i++, j++){
			byte q1=quality[i];
			float f1=PROB_CORRECT[q1];
			key1*=f1;
			byte q2=quality[j];
			float f2=PROB_CORRECT[q2];
			key2*=f2;
		}
		out[0]=1-key1;
		out[mid]=1-key2;
		
		for(int a=0, b=keylen, c=mid, d=mid+keylen; d<quality.length;
				a++, b++, c++, d++){
			byte qa=quality[a];
			byte qb=quality[b];
			byte qc=quality[c];
			byte qd=quality[d];
			float ipa=PROB_CORRECT_INVERSE[qa];
			float ipc=PROB_CORRECT_INVERSE[qc];
			float pb=PROB_CORRECT[qb];
			float pd=PROB_CORRECT[qd];
			key1=key1*ipa*pb;
			key2=key2*ipc*pd;
			out[a+1]=1-key1;
			out[c+1]=1-key2;
		}
		return out;
	}

	/**
	 * Generates synthetic quality array with realistic quality distribution.
	 * Creates degrading quality toward sequence ends and applies random variance.
	 *
	 * @param length Length of quality array to generate
	 * @param randyQual Random number generator
	 * @param minQual Minimum quality value
	 * @param maxQual Maximum quality value
	 * @param baseQuality Starting quality value
	 * @param slant Quality degradation rate along sequence
	 * @param variance Random variance to apply
	 * @return Synthetic quality array
	 */
	public static byte[] makeQualityArray(int length, Random randyQual,
			int minQual, int maxQual, byte baseQuality, byte slant, int variance) {
		byte[] out=new byte[length];
		
		for(int i=0; i<length; i++){
			byte q=(byte)(baseQuality-(slant*i)/length);
			
			int hilo=randyQual.nextInt();
			
//			if((hilo&7)>0){
//				int range=Tools.max(1, maxQual-q+1);
//				int delta=Tools.min(randyQual.nextInt(range), randyQual.nextInt(range));
//				q=(byte)(q+delta);
//			}else{
//				int range=Tools.max(1, q-minQual+1);
//				int delta=Tools.min(randyQual.nextInt(range), randyQual.nextInt(range), randyQual.nextInt(range));
//				q=(byte)(q-delta);
//			}
			
			if((hilo&15)>0){
				int range=Tools.max(1, maxQual-q+1);
				int delta=(randyQual.nextInt(range)+randyQual.nextInt(range+1))/2;
				q=(byte)(q+delta);
			}else{
				int range=Tools.max(1, q-minQual+1);
				int delta=Tools.min(randyQual.nextInt(range), randyQual.nextInt(range));
				q=(byte)(q-delta);
			}
			q=(byte)Tools.min(Tools.max(q, minQual), maxQual);
			out[i]=q;
		}
		
		if(length>50){
			final int x=length/10;
			for(int i=0; i<x; i++){
				int y=x-i;
				out[i]=(byte)Tools.max(out[i]-(y+randyQual.nextInt(y+1))/2, minQual);
				out[length-i-1]=(byte)Tools.max(out[length-i-1]-(y+randyQual.nextInt(y+1))/2, minQual);
			}
		}
		
		int delta=0;
		if(variance>0){
			delta=(byte)(randyQual.nextInt(variance+1)+randyQual.nextInt(variance+1)-variance);
		}
		for(int i=0; i<out.length; i++){
			int x=Tools.mid(2, out[i]+delta, 41);
			out[i]=(byte)x;
		}
		
		return out;
	}
	
	/**
	 * Modifies offset array by removing positions with very high error probability.
	 * Removes middle elements with probability ≥0.98 and adjusts adjacent positions.
	 *
	 * @param offsets Array of offset positions
	 * @param keyProbs Error probabilities for each position
	 * @return Modified offset array with problematic positions removed
	 */
	public static int[] modifyOffsets(int[] offsets, float[] keyProbs) {
		if(offsets==null || offsets.length<3){return offsets;}

		int index=0;
		float max=keyProbs[offsets[0]];
		final int maxOffset=offsets[offsets.length-1];
		
		for(int i=1; i<offsets.length; i++){
			float f=keyProbs[offsets[i]];
			if(f>max){
				max=f;
				index=i;
			}
		}
		
		if(index==0 || index==offsets.length-1){return offsets;}
		if(max<.98f){return offsets;}
		
		final int removed=offsets[index];
		{
			int[] offsets2=new int[offsets.length-1];
			for(int i=0; i<index; i++){offsets2[i]=offsets[i];}
			for(int i=index; i<offsets2.length; i++){offsets2[i]=offsets[i+1];}
			offsets=offsets2;
			offsets2=null;
		}
		
		if(index==0){
			assert(false);
//			int i=offsets[0];
//			assert(i>removed && removed>=0);
//			while(i>removed && keyProbs[i-1]>=keyProbs[i]){i--;}
//			offsets[0]=i;
		}else if(index==offsets.length){
			assert(false);
//			int i=offsets[offsets.length-1];
//			assert(i<removed && removed==maxOffset);
//			while(i<removed && keyProbs[i+1]>=keyProbs[i]){i++;}
//			offsets[offsets.length-1]=i;
		}else if(offsets.length>2){
			if(index==offsets.length-1){
				assert(index>1);
				int i=offsets[index-1]; //5, 7, 9, 5, 6
				assert(i<removed && removed<maxOffset) : i+", "+removed+", "+maxOffset+", "+index+", "+offsets.length;
				while(i<removed-1 && keyProbs[i+1]>=keyProbs[i]){i++;}
				offsets[index-1]=i;
			}else{
				assert(index<offsets.length-1 && index>0);
				int i=offsets[index];
				assert(i>removed && removed>=0);
				while(i>removed+1 && keyProbs[i-1]>=keyProbs[i]){i--;}
				offsets[index]=i;
			}
		}
		
		return offsets;
	}
	
	/** Requires qualities under MATRIX_SIZE */
	public static byte qualsToPhred(byte qa, byte qb){
		return PHRED_MATRIX[qa][qb];
	}
	
	/** Safe version for qualities >=MATRIX_SIZE */
	public static byte qualsToPhredSafe(byte qa, byte qb){
		qa=Tools.max((byte)0, Tools.min(qa, MATRIX_SIZE));
		qb=Tools.max((byte)0, Tools.min(qb, MATRIX_SIZE));
		return (qa<=qb) ? PHRED_MATRIX[qa][qb] : PHRED_MATRIX[qb][qa];
	}
	
	/**
	 * Computes combined error probability from two quality scores.
	 * @param qa First quality score
	 * @param qb Second quality score
	 * @return Combined error probability
	 */
	public static float qualsToProbError(byte qa, byte qb){
		return ERROR_MATRIX[qa][qb];
	}
	
	/**
	 * Computes combined correct probability from two quality scores.
	 * @param qa First quality score
	 * @param qb Second quality score
	 * @return Combined correct probability (1 - error probability)
	 */
	public static float qualsToProbCorrect(byte qa, byte qb){
		return 1-qualsToProbError(qa, qb);
	}
	
	/**
	 * Safe version of qualsToProbError that handles quality values ≥MATRIX_SIZE.
	 * @param qa First quality score
	 * @param qb Second quality score
	 * @return Combined error probability
	 */
	public static float qualsToProbErrorSafe(byte qa, byte qb){
		qa=Tools.max((byte)0, Tools.min(qa, MATRIX_SIZE));
		qb=Tools.max((byte)0, Tools.min(qb, MATRIX_SIZE));
		return (qa<=qb) ? ERROR_MATRIX[qa][qb] : ERROR_MATRIX[qb][qa];
	}
	
	/**
	 * Safe version of qualsToProbCorrect that handles quality values ≥MATRIX_SIZE.
	 * @param qa First quality score
	 * @param qb Second quality score
	 * @return Combined correct probability
	 */
	public static float qualsToProbCorrectSafe(byte qa, byte qb){
		return 1-qualsToProbErrorSafe(qa, qb);
	}

	/**
	 * Creates fake quality array filled with single quality value.
	 * @param q Quality value to fill array (0-127)
	 * @param len Length of array to create
	 * @return Quality array filled with specified value
	 */
	public static byte[] fakeQuality(int q, int len){
		assert(q>=0 && q<=127);
		byte[] r=new byte[len];
		Arrays.fill(r, (byte)q);
		return r;
	}
	
	/*-------------------- Fields --------------------*/

	/*-------------------- Final Fields --------------------*/

	/*-------------------- Static Fields --------------------*/
	
	/** Maximum quality value for matrix operations */
	public static final byte MATRIX_SIZE=50;
	
	/** Probability that this base is an error */
	public static final float[] PROB_ERROR=makeQualityToFloat(128);
	/** 1/PROB */
	public static final float[] PROB_ERROR_INVERSE=makeInverse(PROB_ERROR);
	
	/** Lookup array for converting Phred scores to correct probabilities */
	public static final float[] PROB_CORRECT=oneMinus(PROB_ERROR);
	/** Inverse of PROB_CORRECT array for efficient calculations */
	public static final float[] PROB_CORRECT_INVERSE=makeInverse(PROB_CORRECT);
	
	/** Probability that at least one base will be incorrect, given two quality scores */
	public static final float[][] ERROR_MATRIX=makeErrorMatrix(PROB_ERROR, MATRIX_SIZE);
	
	/** Combined phred score given two quality scores */
	public static final byte[][] PHRED_MATRIX=makePhredMatrix(ERROR_MATRIX);
	
	/*-------------------- Constants --------------------*/

	/*-------------------- Initializers --------------------*/

	/**
	 * Converts array of Phred scores to error probabilities.
	 * @param trimq Array of Phred scores (may be null)
	 * @return Array of error probabilities, or null if input is null
	 */
	public static float[] phredToProbError(float[] trimq){
		if(trimq==null){return null;}
		float[] trimE=trimq.clone();
		for(int i=0; i<trimE.length; i++){
			trimE[i]=(float)QualityTools.phredToProbError(trimE[i]);
		}
		return trimE;
	}
	
	/**
	 * Converts correct probability to Phred score.
	 * @param prob Correct probability (0.0-1.0)
	 * @return Phred score
	 */
	public static byte probCorrectToPhred(double prob){
		return probErrorToPhred(1-prob);
	}
	
	/**
	 * Converts error probability to Phred score with rounding.
	 * @param prob Error probability (0.0-1.0)
	 * @return Phred score
	 */
	public static byte probErrorToPhred(double prob){
		return probErrorToPhred(prob, true);
	}
	
	/**
	 * Converts Phred score to error probability using standard formula.
	 * Handles special cases for very low quality scores (≤1).
	 * @param q Phred score
	 * @return Error probability
	 */
	public static double phredToProbError(double q){
		if(q<=0){return 0.75;}
		if(q<=1){return 0.75-q*0.05;}
		return Tools.min(0.7, Math.pow(10, -0.1*q));
	}
	
	/**
	 * Converts error probability to Phred score with optional rounding.
	 * @param prob Error probability (0.0-1.0)
	 * @param round Whether to round the result
	 * @return Phred score clamped to valid range
	 */
	public static byte probErrorToPhred(double prob, boolean round){
		double phred=probErrorToPhredDouble(prob);
		final int q=round ? (int)Math.round(phred) : (int)phred;
		return  (byte)Tools.mid(0, q, Read.MAX_CALLED_QUALITY());
	}
	
	/**
	 * Converts error probability to Phred score as double precision.
	 * Uses standard formula: -10 * log10(probability).
	 * @param prob Error probability (0.0-1.0)
	 * @return Phred score as double
	 */
	public static double probErrorToPhredDouble(double prob){
		if(prob>=1){return 0;}
		if(prob<=0.000001){return 60;}
		
		double phred=-10*Math.log10(prob);
		return phred;
	}
	
	/**
	 * Creates lookup table for converting Phred scores to error probabilities.
	 * Sets special values for quality scores 0 and 1.
	 * @param n Size of lookup table
	 * @return Array mapping Phred scores to error probabilities
	 */
	private static final float[] makeQualityToFloat(int n){
		float[] r=new float[n];
		for(int i=0; i<n; i++){
			float x=(float)Math.pow(10, 0-.1*i);
			r[i]=x;
		}
		r[0]=.75f;
		r[1]=.7f;
//		assert(false) : Arrays.toString(r);
		return r;
	}
	
	/**
	 * Creates inverse probability array for efficient calculations.
	 * @param prob Probability array
	 * @return Array of inverse probabilities (1/prob[i])
	 */
	private static final float[] makeInverse(float[] prob){
		float[] r=new float[prob.length];
		for(int i=0; i<r.length; i++){r[i]=1/prob[i];}
		return r;
	}
	
	/**
	 * Creates complement probability array.
	 * @param prob Probability array
	 * @return Array of complement probabilities (1-prob[i])
	 */
	private static final float[] oneMinus(float[] prob){
		float[] r=new float[prob.length];
		for(int i=0; i<r.length; i++){r[i]=1-prob[i];}
		return r;
	}
	
	/**
	 * Creates matrix for combining two quality scores into error probability.
	 * Uses formula: 1-((1-a)*(1-b)) for independent error events.
	 *
	 * @param prob Error probability lookup array
	 * @param maxq Maximum quality value for matrix
	 * @return 2D matrix of combined error probabilities
	 */
	private static final float[][] makeErrorMatrix(float[] prob, byte maxq){
		maxq++;
		float[][] matrix=new float[maxq][maxq];
		for(int i=0; i<maxq; i++){
			for(int j=0; j<maxq; j++){
				float a=prob[i], b=prob[j];
				matrix[i][j]=1-((1-a)*(1-b));
			}
		}
		return matrix;
	}
	
	/**
	 * Creates Phred score matrix from error probability matrix.
	 * Converts error probabilities back to Phred scores.
	 * @param error 2D error probability matrix
	 * @return 2D matrix of combined Phred scores
	 */
	private static final byte[][] makePhredMatrix(float[][] error){
		final int maxq=error.length;
		byte[][] matrix=new byte[maxq][maxq];
		for(int i=0; i<maxq; i++){
			for(int j=0; j<maxq; j++){
				matrix[i][j]=probCorrectToPhred(1-error[i][j]);
			}
		}
		return matrix;
	}

	/*-------------------- Notes --------------------*/

}
