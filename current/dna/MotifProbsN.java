package dna;

import shared.KillSwitch;
import shared.Tools;

/**
 * Implements N-base motif probability calculations for biological sequence analysis.
 * Extends the base Motif class to handle variable-length N-tuples with probabilistic
 * matching and position-specific scoring matrices. Supports exact matching and
 * probabilistic strength calculations with importance weighting.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class MotifProbsN extends Motif {
	
	/**
	 * Test method demonstrating MotifProbsN functionality with example sequences.
	 * Creates an "Exon Stops MP3" motif and tests exact matching and strength
	 * calculation against a sample DNA sequence.
	 * @param args Command-line arguments (currently unused)
	 */
	public static void main(String args[]){
		
		String s1="ATN";
		String s2="CTATGCCCATCTGATGGCATGAGGATGAA";
		
//		if(args.length>0){s1=args[0];}
//		if(args.length>1){s2=args[1];}
		
		MotifProbsN m=makeMotif("Exon Stops MP3", 10, 3, 3);
		
		System.out.println("Made motif "+m.name);
		
		String source=s2;
		
		
		int x=m.countExact(source);
		System.out.println(x+" matches.");
		
		byte[] sbytes=source.getBytes();
		
		for(int i=0; i<s2.length(); i++){
			String sub=s2.substring(i, min(i+m.probs.length, s2.length()));
			float p=m.matchStrength(sbytes, i);
			System.out.println(sub+Tools.format(": \t%.5f ->\t%.5f", p, m.normalize(p)));
		}
		
	}
	
	/**
	 * Factory method to create MotifProbsN instances from matrix data.
	 * Loads position-specific probability matrices from the Matrix registry
	 * and optionally associates percentile data for normalization.
	 *
	 * @param name_ Name of the motif to load from Matrix registry
	 * @param length_ Length of the motif in positions
	 * @param center_ Center position for alignment calculations
	 * @param n_ Number of bases per position (N-tuple size)
	 * @return Configured MotifProbsN instance with loaded probability matrices
	 */
	public static MotifProbsN makeMotif(String name_, int length_, int center_, int n_){
		Matrix mat=Matrix.get(name_);
		assert(mat!=null) : "\nCan't find '"+name_+"' in:\n"+Matrix.keys()+"\n\n";
		float[][] sub=mat.subGrid(center_, length_);
		
//		System.out.println("Found "+name+":\n"+Arrays.toString(sub[preLen]));
		
		assert(sub[0].length==(1<<(2*n_)));
		
		MotifProbsN r=new MotifProbsN(name_, sub, center_, n_);
		
		Matrix percentMatrix=null;
		
		
		try {
			percentMatrix=Matrix.get(name_+", "+r.length+", "+r.center);
		} catch (Exception e) {
			// TODO Auto-generated catch block
//			System.out.println("\nIgnoring missing percentMatrix:\n"+e);
		}
		
		if(percentMatrix!=null){
			r.percentile=percentMatrix.grid[0];
		}
//		r.percentile=percentTable.get(name);
		
		return r;
	}
	
	/**
	 * Constructs MotifProbsN with probability matrix and parameters.
	 * Initializes probability arrays, calculates position importance weights,
	 * adjusts for base probability, and pre-computes normalization constants.
	 *
	 * @param name_ Name identifier for this motif
	 * @param p Two-dimensional probability matrix [position][n-tuple]
	 * @param cen Center position for alignment
	 * @param n Number of bases per position (N-tuple size)
	 */
	public MotifProbsN(String name_, float[][] p, int cen, int n){
		super(name_, p.length, cen);
		
		N=n;
		chunk=KillSwitch.allocByte1D(N);
		baseProb=Motif.baseProbN[N];
		
		probs=p;
		importance=positionImportance(probs);
		
		adjustForBaseProb(probs, baseProb);
		
		double pmin=1, pmax=1;
		
		double sum=0;
		for(int i=0; i<p.length; i++){
			for(int j=0; j<p[i].length; j++){
				sum+=p[i][j];
			}
		}
		matrixAvg=(float)(sum/(p.length*p[0].length));
		
		
		//Adjusts for importance
		for(int i=0; i<probs.length; i++){
			for(int j=0; j<probs[i].length; j++){
				probs[i][j]=(float)Math.pow(probs[i][j], 1+(importance[i]*.8));
			}
		}
		
		
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<probs.length; i++){
			int x=maxPos(probs[i]);
			int y=minPos(probs[i]);
			sb.append((char)numberToBase[x>>(2*(N-1))]);

//			pmax*=probs[i][x]*4; //TODO Note the .25; could be an empirical inverse probability, but that causes complications
//			pmin*=probs[i][y]*4;

			pmax*=probs[i][x];
			pmin*=probs[i][y];
			
//			pmax*=Math.pow(probs[i][x], 1+importance[i]);
//			pmin*=Math.pow(probs[i][y], 1+importance[i]);
			
//			pmax*=(probs[i][x]+(matrixAvg*importance[i]*.1f));
//			pmin*=(probs[i][y]+(matrixAvg*importance[i]*.1f));
		}
		

		maxProb=(float)pmax;
		minProb=(float)pmin;

		invProbDif=1f/(maxProb-minProb);
		invLength=1f/(length);
		
		commonLetters=sb.toString();
		
		lettersUpper=commonLetters.toUpperCase().getBytes();
		lettersLower=commonLetters.toLowerCase().getBytes();
		
		numbers=new byte[commonLetters.length()];
		numbersExtended=new byte[commonLetters.length()];
		
		for(int i=0; i<lettersUpper.length; i++){
			byte b=lettersUpper[i];
			numbers[i]=baseToNumber[b];
			numbersExtended[i]=baseToNumberExtended[b];
		}
		
	}
	
	
	/**
	 * Adjusts probability matrix values by dividing by background base probabilities.
	 * Normalizes position-specific probabilities against expected base frequencies
	 * to emphasize deviations from background.
	 *
	 * @param grid Probability matrix to adjust [position][n-tuple]
	 * @param base Background probability array for each n-tuple
	 */
	public void adjustForBaseProb(float[][] grid, float[] base){
		for(int i=0; i<grid.length; i++){
			for(int j=0; j<grid[i].length; j++){
				grid[i][j]/=base[j];
			}
		}
	}
	
	
	@Override
	public float normalize(double strength){
		double r=strength-minProb;
//		r=r/(maxProb-minProb);
//		r=Math.pow(r, 1d/length);
		r=r*invProbDif;
		r=Math.pow(r, invLength);
		return (float)r;
	}
	
	
	/**
	 * Alternative normalization using logarithmic scaling.
	 * Provides log-space normalization between minimum and maximum
	 * probability bounds for improved dynamic range.
	 *
	 * @param strength Raw probability strength to normalize
	 * @return Log-normalized strength value
	 */
	public float normalize2(double strength){
		double r=Math.log(strength)-Math.log(minProb);
		
		double r2=Math.log(maxProb)-Math.log(minProb);
		
		r=r/r2;
		return (float)r;
	}
	
	
	@Override
	public boolean matchesExactly(byte[] source, int a){
		
		a=a-center;
		if(a<0 || a+length>source.length){return false;}
		
		for(int i=0; i<lettersUpper.length; i++){
			int x=i+a;
			if(source[x]!=lettersUpper[i] && source[x]!=lettersLower[i]){
				return false;
			}
		}
		return true;
	}
	
	
	@Override
	public float matchStrength(byte[] source, int a){
		
		a=a-center;
		if(a<0 || a+length+1>source.length){return minProb;}
		
		float r=1;
		
		for(int i=0; i<probs.length; i++){
			int x=i+a;
			
			for(int c=0; c<N; c++){
				chunk[c]=source[x+c];
			}
			
			int n=AminoAcid.baseTupleToNumber(chunk);
			if(n<0 || n>baseProb.length){return minProb;}
			
//			float p1=(probs[i][n]+(matrixAvg*importance[i]*.1f));
			
//			float p1=(float)Math.pow(probs[i][n], 1+importance[i]); //Note:  Assumes (A,C,G,T) only.
			float p1=probs[i][n]; //Note:  Assumes (A,C,G,T) only.
			
//			float p2=invBaseProb2[n];
//			float p2=4; //TODO
//
//			r=r*p1*p2;
			
			r=r*p1;
		}
		return r;
	}
	
	
	/**
	 * Calculates position importance weights based on deviation from background.
	 * Measures how much each position deviates from background base probabilities,
	 * applying power transformations to emphasize highly informative positions.
	 *
	 * @param rawProbs Raw probability matrix before adjustment
	 * @return Array of importance weights normalized to 0-1 range
	 */
	public float[] positionImportance(float[][] rawProbs){
		float[] base=baseProb;
		float[] out=new float[rawProbs.length];
		
		double maxSum=0;
		
		for(int i=0; i<out.length; i++){
			float[] array=rawProbs[i];
			double sum=0;
			for(int code=0; code<array.length; code++){
				double dif=Math.abs(array[code]-base[code]);
				sum+=Math.pow(dif,1.5); //Raise to a power to increase the effect
			}
			sum=Math.pow(sum, 0.75);
			out[i]=(float)sum;
			if(sum>maxSum){
				maxSum=sum;
			}
		}
		
		for(int i=0; i<out.length; i++){
			out[i]=(float)(out[i]/maxSum);
//			out[i]=out[i]*.9f+.1f; //Weakens the effect
//			out[i]=out[i]*.5f; //makes the scale 0 to .5
		}
		
		return out;
	}

	@Override
	public int numBases() {
		return N;
	}
	
	/** Number of bases per position in the N-tuple motif */
	public final int N;
	
	/** Position-specific probability matrix [position][n-tuple] after adjustment */
	public final float[][] probs;
	/** Position importance weights based on deviation from background */
	public final float[] importance;
	/** Average probability value across all matrix positions */
	public final float matrixAvg;

	/** Uppercase consensus sequence bytes for exact matching */
	public final byte[] lettersUpper;
	/** Lowercase consensus sequence bytes for exact matching */
	public final byte[] lettersLower;
	/** Numeric representation of consensus sequence for calculations */
	public final byte[] numbers;
	/** Extended numeric representation supporting ambiguous bases */
	public final byte[] numbersExtended;
	
	/** Temporary buffer for extracting N-tuples from sequences */
	private final byte[] chunk;
	/** Background probability array for N-tuples used in normalization */
	private final float[] baseProb;
	
	/** Maximum possible probability for this motif across all positions */
	public float maxProb;
	/** Minimum possible probability for this motif across all positions */
	public float minProb;
	
	/** Inverse of probability difference (maxProb - minProb) for normalization */
	public final float invProbDif;
	/** Inverse of motif length for power scaling in normalization */
	public final float invLength;
	
}
