package dna;

/**
 * Abstract base class for sequence motif matching and analysis.
 * Provides methods for exact and extended motif matching in nucleotide sequences,
 * along with statistical analysis capabilities for motif strength evaluation.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public abstract class Motif {
	
	/**
	 * Constructs a Motif with specified parameters.
	 * @param name_ Name identifier for this motif
	 * @param length_ Total length of the motif in bases
	 * @param center_ Position of the central base (0-based indexing)
	 */
	public Motif(String name_, int length_, int center_){
		center=center_;
		length=length_;
		suffix=length-center-1;
		name=name_;
		
//		assert(center>=0 && center<length);
	}
	
	
	/**
	 * Counts exact matches of this motif in the entire string.
	 * @param s The sequence string to search
	 * @return Number of exact motif matches found
	 */
	public final int countExact(String s){
		return countExact(s, 0, s.length());
	}
	
	
	/**
	 * Counts exact matches of this motif within a substring range.
	 *
	 * @param s The sequence string to search
	 * @param a Start position (inclusive)
	 * @param b End position (exclusive)
	 * @return Number of exact motif matches found in the range
	 */
	public final int countExact(String s, int a, int b){
		return countExact(s.getBytes(), a, b);
	}
	
	
	/**
	 * Counts extended matches of this motif in the entire string.
	 * Extended matching allows for ambiguous nucleotides and fuzzy matching.
	 * @param s The sequence string to search
	 * @return Number of extended motif matches found
	 */
	public final int countExtended(String s){
		return countExtended(s, 0, s.length());
	}
	
	
	/**
	 * Counts extended matches of this motif within a substring range.
	 * Extended matching allows for ambiguous nucleotides and fuzzy matching.
	 *
	 * @param s The sequence string to search
	 * @param a Start position (inclusive)
	 * @param b End position (exclusive)
	 * @return Number of extended motif matches found in the range
	 */
	public final int countExtended(String s, int a, int b){
		return countExtended(s.getBytes(), a, b);
	}
	
	
	/**
	 * Counts exact matches of this motif within a byte array range.
	 * Scans through the specified range looking for positions where the motif
	 * matches exactly at each position.
	 *
	 * @param source The sequence byte array to search
	 * @param a Start position (inclusive)
	 * @param b End position (exclusive)
	 * @return Number of exact motif matches found in the range
	 */
	public final int countExact(byte[] source, int a, int b){
		
		int max=min(b, source.length-1)-length+1;
		
		int count=0;
		
		for(int i=a; i<=max; i++){
			if(matchesExactly(source, i)){count++;}
		}
		
		return count;
		
	}
	
	
	/**
	 * Counts extended matches of this motif within a byte array range.
	 * Extended matching allows for ambiguous nucleotides and fuzzy matching.
	 * Scans through the specified range looking for positions where the motif
	 * matches in extended mode at each position.
	 *
	 * @param source The sequence byte array to search
	 * @param a Start position (inclusive)
	 * @param b End position (exclusive)
	 * @return Number of extended motif matches found in the range
	 */
	public final int countExtended(byte[] source, int a, int b){
		
		int max=min(b, source.length-1)-length+1;
		
		int count=0;
		
		for(int i=a; i<=max; i++){
			if(matchesExtended(source, i)){count++;}
		}
		
		return count;
		
	}
	
	
	/**
	 * Tests if the motif matches exactly at a specific position.
	 * Abstract method to be implemented by subclasses for exact matching logic.
	 *
	 * @param source The sequence byte array
	 * @param a The position to test for a match
	 * @return true if motif matches exactly at position a
	 * @throws RuntimeException Always thrown in base class implementation
	 */
	public boolean matchesExactly(byte[] source, int a){
		throw new RuntimeException();
	}
	
	
	/**
	 * Tests if the motif matches in extended mode at a specific position.
	 * Abstract method to be implemented by subclasses for extended matching logic.
	 * Extended matching typically allows for ambiguous nucleotides.
	 *
	 * @param source The sequence byte array
	 * @param a The position to test for a match
	 * @return true if motif matches in extended mode at position a
	 * @throws RuntimeException Always thrown in base class implementation
	 */
	public boolean matchesExtended(byte[] source, int a){
		throw new RuntimeException();
	}
	
	/**
	 * Normalizes a strength value to a float.
	 * Default implementation simply casts to float.
	 * @param strength The strength value to normalize
	 * @return The normalized strength as a float
	 */
	public float normalize(double strength){
		return (float)strength;
	}
	
	
	/**
	 * Calculates the match strength at a specific position.
	 * Returns 1.0 for exact matches, 0.0 for no match.
	 *
	 * @param source The sequence byte array
	 * @param a The position to evaluate
	 * @return Match strength (1.0 for exact match, 0.0 for no match)
	 */
	public float matchStrength(byte[] source, int a){
		return(matchesExactly(source, a) ? 1 : 0);
	}
	
	
	/**
	 * Finds the index of the minimum value in a float array.
	 * @param array The array to search
	 * @return Index of the minimum value
	 */
	public static final int minPos(float[] array){
		int pos=0;
		for(int i=1; i<array.length; i++){
			if(array[i]<array[pos]){pos=i;}
		}
		return pos;
	}
	
	
	/**
	 * Finds the index of the maximum value in a float array.
	 * @param array The array to search
	 * @return Index of the maximum value
	 */
	public static final int maxPos(float[] array){
		int pos=0;
		for(int i=1; i<array.length; i++){
			if(array[i]>array[pos]){pos=i;}
		}
		return pos;
	}
	
	@Override
	public String toString(){
		return name+", "+length+", "+center;
	}
	

	/** Name identifier for this motif */
	public final String name;
	/** String containing common letters or characters associated with this motif */
	public String commonLetters;
	/** Position of the central base in the motif (0-based indexing) */
	public final int center;
	/** Total length of the motif in bases */
	public final int length;
	/**
	 * Length of the suffix portion after the center position (length - center - 1)
	 */
	public final int suffix;
	

	/** Returns the smaller of two integers */
	static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger of two integers */
	static final int max(int x, int y){return x>y ? x : y;}
	/** Returns the smaller of two floats */
	static final float min(float x, float y){return x<y ? x : y;}
	/** Returns the larger of two floats */
	static final float max(float x, float y){return x>y ? x : y;}
	
	/** Array mapping numeric values to nucleotide bases */
	static final byte[] numberToBase=AminoAcid.numberToBase;
	/**
	 * Array mapping numeric values to extended nucleotide bases including ambiguous codes
	 */
	static final byte[] numberToBaseExtended=AminoAcid.numberToBaseExtended;
	/** Array mapping nucleotide bases (ACGTN) to numeric values */
	static final byte[] baseToNumber=AminoAcid.baseToNumberACGTN;
	/** Array mapping extended nucleotide bases to numeric values */
	static final byte[] baseToNumberExtended=AminoAcid.baseToNumberExtended;
	
	/** Base probability distribution for 4 standard nucleotides (A, C, G, T) */
	static final float[] baseProb1={0.256614f, 0.226617f, 0.238012f, 0.278756f};
	
	//Within 200 of exon and gene ends only
	/**
	 * Base probability distribution for 16-element extended nucleotide set, calculated within 200bp of exon and gene ends
	 */
	static final float[] baseProb2={
		0.076019f, 0.046405f, 0.071754f, 0.062437f, 0.067143f, 0.066057f, 0.020333f, 0.073085f,
		0.060553f, 0.054897f, 0.068741f, 0.053822f, 0.052896f, 0.059260f, 0.077188f, 0.089412f
	};
	
	//name: Overall Frequency MP3
	/**
	 * Base probability distribution for 64-element extended nucleotide set, overall frequency MP3
	 */
	static final float[] baseProb3={
		0.027343f, 0.011857f, 0.018295f, 0.018524f, 0.015942f, 0.012337f, 0.003792f, 0.014333f,
		0.019988f, 0.015837f, 0.020411f, 0.015518f, 0.014382f, 0.011355f, 0.016466f, 0.020234f,
		0.014364f, 0.014299f, 0.022875f, 0.015605f, 0.018893f, 0.019412f, 0.006677f, 0.021076f,
		0.003629f, 0.005854f, 0.006783f, 0.004067f, 0.010491f, 0.018413f, 0.024257f, 0.019924f,
		0.018029f, 0.010640f, 0.019427f, 0.012458f, 0.015158f, 0.017025f, 0.006167f, 0.016547f,
		0.018098f, 0.016891f, 0.020042f, 0.013710f, 0.010580f, 0.010773f, 0.018026f, 0.014443f,
		0.016281f, 0.009609f, 0.011157f, 0.015849f, 0.017150f, 0.017284f, 0.003696f, 0.021130f,
		0.018839f, 0.016316f, 0.021506f, 0.020527f, 0.017442f, 0.018720f, 0.018440f, 0.034811f
	};
	
//	protected static final Hashtable<String, float[]> percentTable=makePercentTable();
//
//	private static final Hashtable<String, float[]> makePercentTable(){
//
//		String[] keys={
//				"Exon Stops MP3",
//		};
//
//		float[][] values={
//				{
//					0.00234f, 0.01071f, 0.02476f, 0.05155f, 0.08682f, 0.1453f, 0.22434f, 0.29615f, 0.36233f, 0.41034f,
//					0.46028f, 0.52224f, 0.58198f, 0.63879f, 0.68356f, 0.70622f, 0.7268f, 0.75131f, 0.77065f, 0.79546f,
//					0.82445f, 0.85279f, 0.86899f, 0.88287f, 0.89197f, 0.90166f, 0.91405f, 0.93129f, 0.94708f, 0.95521f,
//					0.96106f, 0.96293f, 0.9663f, 0.97242f, 0.97662f, 0.97866f, 0.98017f, 0.98242f, 0.98459f, 0.98703f,
//					0.98957f, 0.99064f, 0.99157f, 0.99286f, 0.9952f, 0.99721f, 0.99858f, 0.99914f, 0.99967f, 0.9999f, 0.99998f
//				},
//		};
//
//		Hashtable<String, float[]> r= new Hashtable<String, float[]>();
//		for(int i=0; i<keys.length; i++){
//			r.put(keys[i], values[i]);
//		}
//
//		return r;
//	}

	/** Inverted base probabilities for baseProb1 (reciprocals) */
	static final float[] invBaseProb1=invert(baseProb1);
	
	/** Inverted base probabilities for baseProb2 (reciprocals) */
	static final float[] invBaseProb2=invert(baseProb2);
	
	/** Inverted base probabilities for baseProb3 (reciprocals) */
	static final float[] invBaseProb3=invert(baseProb3);
	
	/** Array of base probability distributions indexed by nucleotide set size */
	static final float[][] baseProbN={
		null,
		baseProb1,
		baseProb2,
		baseProb3
	};
	
	/**
	 * Array of inverted base probability distributions indexed by nucleotide set size
	 */
	static final float[][] invBaseProbN={
		null,
		invBaseProb1,
		invBaseProb2,
		invBaseProb3
	};
	
	/**
	 * Creates an inverted array where each element is the reciprocal of the input.
	 * Used to convert probability arrays to inverse probability arrays.
	 * @param in Input probability array
	 * @return Array where each element is 1/in[i]
	 */
	private static final float[] invert(float[] in){
		float[] out=new float[in.length];
		for(int i=0; i<in.length; i++){
			out[i]=1f/in[i];
		}
		return out;
	}
	
	/** Percentile array for strength calculations specific to this motif */
	protected float[] percentile;

	/**
	 * Returns the number of bases used by this motif type.
	 * Abstract method to be implemented by subclasses.
	 * @return Number of bases in the motif's alphabet
	 */
	public abstract int numBases();
	
	/**
	 * Converts a raw strength value to a percentile using the motif's percentile array.
	 * Performs linear interpolation between array values for smooth results.
	 *
	 * @param strength Raw strength value to convert
	 * @return Percentile value (0.0 to 1.0) corresponding to the strength
	 * @throws RuntimeException if percentile array is null
	 */
	public float percentile(float strength){
//		float[] array=percentiles[numBases()];
		
		if(percentile==null){
			throw new RuntimeException("Can't find percentile array for "+this);
		}
		
		float[] array=percentile;
		
		int index=(int)(strength*array.length);
		
//		System.out.print(" *** index = "+index+" -> "+array[index]+" -> "+array[index+1]+" *** ");
		
		if(index>=array.length-1){return 1;}
		
		float a, b;
		if(index==0){
			a=0;
			b=array[0];
		}else{
			a=array[index];
			b=array[index+1];
		}
		
		float ratio=strength-(index/((float)array.length));
		
		return ratio*b+(1-ratio)*a;
		
	}
	
}
