package prot;

import dna.AminoAcid;

/**
 * BLOSUM62 substitution matrix plus affine-gap costs and Karlin-Altschul
 * statistical parameters for protein local alignment.
 *
 * <p>Residues are encoded with BBTools' {@link AminoAcid#acidToNumberExtended}
 * table: the 20 standard amino acids map to 0-19, and X (unknown, plus the
 * ambiguity codes B/Z/J that are folded into X on input) maps to
 * {@link #X_CODE}. The score matrix is indexed directly by those encoded
 * values, so callers never need to know the internal residue ordering.</p>
 *
 * <p>Gap costs and statistical parameters are the published NCBI BLAST defaults
 * for gapped BLOSUM62 with gap-open 11 / gap-extend 1. A gap of length L costs
 * {@code GAP_OPEN + L*GAP_EXTEND} (i.e. a length-1 gap costs 12), matching
 * BLAST's "Existence: 11, Extension: 1" convention.</p>
 *
 * @author Eru
 */
public final class Blosum62 {

	/** Utility class; no instances. */
	private Blosum62(){}

	/** Gap-open (existence) penalty; BLAST default for BLOSUM62. */
	public static final int GAP_OPEN=11;
	/** Gap-extension penalty; BLAST default for BLOSUM62. */
	public static final int GAP_EXTEND=1;

	/** Karlin-Altschul lambda for gapped BLOSUM62 (11,1); NCBI value. */
	public static final double LAMBDA=0.267;
	/** Karlin-Altschul K for gapped BLOSUM62 (11,1); NCBI value. */
	public static final double K=0.041;
	/** Natural log of K, precomputed for bitscore conversion. */
	public static final double LN_K=Math.log(K);
	/** Natural log of 2, precomputed for bitscore conversion. */
	public static final double LN_2=Math.log(2);

	/** Encoded value of X (unknown residue) under acidToNumberExtended. */
	public static final byte X_CODE=(byte)AminoAcid.acidToNumberExtended['X'];

	/** Sentinel used for gap columns in an alignment (never a valid residue). */
	public static final byte GAP=-1;

	/** Residue letters, in the order the raw DATA rows/columns are written. */
	private static final String LETTERS="ARNDCQEGHILKMFPSTWYVX";

	/**
	 * Raw BLOSUM62 scores restricted to the 20 standard amino acids plus X,
	 * in {@link #LETTERS} order. This is the canonical NCBI BLOSUM62 (the B/Z
	 * and stop rows are omitted; B/Z are folded into X on input).
	 */
	private static final int[][] DATA={
		{ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0, 0},//A
		{-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1},//R
		{-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,-1},//N
		{-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,-1},//D
		{ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-2},//C
		{-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,-1},//Q
		{-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,-1},//E
		{ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1},//G
		{-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,-1},//H
		{-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-1},//I
		{-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-1},//L
		{-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,-1},//K
		{-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-1},//M
		{-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-1},//F
		{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2},//P
		{ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0},//S
		{ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0, 0},//T
		{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-2},//W
		{-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-1},//Y
		{ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-1},//V
		{ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1} //X
	};

	/** Highest encoded residue value handled (X_CODE). */
	private static final int SIZE=X_CODE+1;

	/** Score matrix indexed by encoded residue value (acidToNumberExtended). */
	private static final int[][] MATRIX=buildMatrix();

	private static int[][] buildMatrix(){
		//The extended-alphabet encoding must place X at 21; crash loud if BBTools drifts.
		assert(X_CODE==21) : "acidToNumberExtended['X'] changed to "+X_CODE+"; update Blosum62.";
		assert(DATA.length==LETTERS.length());
		final int[][] m=new int[SIZE][SIZE];
		//Fill with a sentinel so any unmapped pair is caught, not silently zero.
		for(int i=0; i<SIZE; i++){
			for(int j=0; j<SIZE; j++){m[i][j]=Integer.MIN_VALUE/4;}
		}
		for(int i=0; i<LETTERS.length(); i++){
			final int ei=AminoAcid.acidToNumberExtended[LETTERS.charAt(i)];
			assert(ei>=0 && ei<SIZE) : "Bad encoding for "+LETTERS.charAt(i);
			for(int j=0; j<LETTERS.length(); j++){
				final int ej=AminoAcid.acidToNumberExtended[LETTERS.charAt(j)];
				m[ei][ej]=DATA[i][j];
			}
		}
		return m;
	}

	/**
	 * Substitution score for a pair of encoded residues.
	 * @param encA First residue, encoded via acidToNumberExtended (0-19 or X_CODE).
	 * @param encB Second residue, encoded the same way.
	 * @return BLOSUM62 score for the pair.
	 */
	public static final int score(final byte encA, final byte encB){
		assert(encA>=0 && encA<SIZE && encB>=0 && encB<SIZE) :
			"Unencoded residue in scoring: "+encA+","+encB;
		final int s=MATRIX[encA][encB];
		assert(s>Integer.MIN_VALUE/8) : "Unmapped residue pair: "+encA+","+encB;
		return s;
	}

	/**
	 * True if the encoded residue is one of the 20 standard amino acids
	 * (i.e. eligible to count as an identity). X is deliberately excluded.
	 * @param enc Encoded residue value.
	 * @return True for standard residues 0-19.
	 */
	public static final boolean isStandard(final byte enc){return enc>=0 && enc<=19;}

	/**
	 * Validates and encodes a raw amino-acid sequence per the frozen protein
	 * search contract: uppercase on input; B/Z/J folded to X; X accepted; and
	 * '*', gap characters, digits, and any other symbol rejected loudly.
	 *
	 * @param raw Raw ASCII amino-acid bytes.
	 * @param id Sequence identifier, used only in error messages.
	 * @return Encoded residue array (values 0-19 or X_CODE).
	 * @throws RuntimeException If an illegal residue is present (crash loud).
	 */
	public static final byte[] encode(final byte[] raw, final String id){
		final byte[] enc=new byte[raw.length];
		for(int i=0; i<raw.length; i++){
			final int c=raw[i]&0xFF;
			final byte e=encodeResidue(c);
			if(e<0){
				throw new RuntimeException("Illegal residue '"+(char)c+"' (0x"+
					Integer.toHexString(c)+") in sequence '"+id+"' at position "+(i+1)+".");
			}
			enc[i]=e;
		}
		return enc;
	}

	/**
	 * Encodes a single residue character, returning -1 for illegal input.
	 * @param c Residue character code (0-255).
	 * @return Encoded value (0-19 or X_CODE), or -1 if illegal.
	 */
	private static final byte encodeResidue(final int c){
		if(c<0 || c>=128){return -1;}
		//Standard residue (upper or lower via acidToNumber, which handles both).
		final byte std=AminoAcid.acidToNumber[c];
		if(std>=0 && std<=19){return std;}
		switch(Character.toUpperCase(c)){
			case 'X': case 'B': case 'Z': case 'J': return X_CODE;//ambiguity -> X
			default: return -1;//'*', '-', '.', digits, punctuation, U, O, etc.
		}
	}
}
