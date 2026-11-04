package dna;

import java.util.Arrays;
import java.util.HashMap;

import align2.QualityTools;
import shared.KillSwitch;
import shared.Tools;
import simd.Vector;
import structures.ByteBuilder;

/**
 * Amino-acid descriptors, shared alphabet tables and sequence conversion kernels.
 * Packed nucleotide codes use A=0, C=1, G=2, T/U=3 in either letter case. Most
 * character/codon methods directly index 128-entry ASCII tables; their inputs must
 * be valid indices unless the method explicitly guards signed bytes. The separate
 * 256-entry byte-conversion tables are indexed with unsigned byte values.
 * Public arrays/maps are shared mutable state; callers must not mutate them during
 * concurrent use. Descriptor codon arrays are exposed too. Translation overloads
 * retain different unknown-base behavior, documented below.
 * @author Brian Bushnell
 * @date July 1, 2010
 *
 */
public final class AminoAcid{

	/**
	 * Demonstrates translation and prints the complement table. Before printing,
	 * replaces any nonpositive complement entry with its identity mapping in place.
	 * @param args Command line arguments (not used)
	 */
	public static void main(String[] args){
		//		for(String s : stringToAA.keySet()){
		//			System.out.println(s+"\t->\t"+stringToAA.get(s));
		//		}

		String bases="atctgatTGGcgcgatatatcg";
		String acids=stringToAAs(bases);

		System.out.println(bases+" -> "+acids);
		for(int i=0; i<baseToComplementExtended.length; i++){
			if(baseToComplementExtended[i]<1){baseToComplementExtended[i]=(byte)i;}
		}
		System.out.println(new String(baseToComplementExtended));
		System.out.println(Arrays.toString(baseToComplementExtended));
	}

	/** Unused no-argument path; delegates null to the descriptor parser and fails there.
	 * Normal descriptors are created by the other private constructors during initialization. */
	private AminoAcid(){
		this(null);
		assert(false);
		System.exit(0);
	}

	/**
	 * Constructor for creating amino acid from comma-separated string.
	 * Splits on the exact comma-space separator: name, symbol, letter, codon1, ...
	 * @param line Nonnull descriptor with at least three fields and a nonempty letter field
	 */
	private AminoAcid(String line){
		String[] s2=line.split(", ");
		String[] s3=new String[s2.length-3];
		for(int i=3; i<s2.length; i++){
			s3[i-3]=s2[i];
		}

		name=s2[0];
		symbol=s2[1];
		letter=s2[2].charAt(0);
		codeStrings=s3;
	}

	/**
	 * Constructor for creating amino acid with specified properties.
	 *
	 * @param n Full name of the amino acid
	 * @param c3 Three-letter symbol code
	 * @param c1 Single-letter code
	 * @param bases Borrowed codon array, retained without copying
	 */
	private AminoAcid(String n, String c3, String c1, String[] bases){
		name=n;
		symbol=c3;
		letter=c1.charAt(0);
		codeStrings=bases;
	}

	/** Formats the name, symbol, letter and current codon-array contents. */
	@Override
	public String toString(){
		return name+", "+symbol+", "+letter+", "+Arrays.toString(codeStrings);
	}

	/**
	 * Converts a binary k-mer representation to DNA string.
	 * Each base is encoded in 2 bits: A=0, C=1, G=2, T=3.
	 *
	 * @param kmer Binary representation of the k-mer
	 * @param k Length of the k-mer
	 * @return DNA sequence string
	 */
	public static String kmerToString(long kmer, int k){
		ByteBuilder sb=new ByteBuilder(k);
		for(int i=0; i<k; i++){
			int x=(int)(kmer&3);
			sb.append((char)numberToBase[x]);
			kmer>>=2;
		}
		return sb.reverse().toString();
	}

	/**
	 * Converts a DNA string to binary k-mer representation.
	 * Each base is encoded in 2 bits: A=0, C=1, G=2, T=3.
	 * Does not reject unknown symbols or prevent shifts from discarding high bits.
	 * @param s Up to 32 defined ASCII bases; accepts lowercase and U as T
	 * @return Packed representation; a 32-base value may be negative
	 */
	public static long stringToKmer(String s){
		long kmer=0;
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			kmer=(kmer<<2)|(baseToNumber[c]);
		}
		return kmer;
	}

	/**
	 * Converts a binary amino acid k-mer to protein sequence string.
	 * Uses five-bit indices into numberToAcid, which includes 20 residues and stop.
	 * Every decoded index must be 0..20; this method does not decode extended amino codes.
	 *
	 * @param kmer Binary representation of amino acid k-mer
	 * @param k Length of the amino acid k-mer
	 * @return Protein sequence string
	 */
	public static String kmerToStringAA(long kmer, int k){
		ByteBuilder sb=new ByteBuilder(k);
		for(int i=0; i<k; i++){
			int x=(int)(kmer&31);
			sb.append((char)numberToAcid[x]);
			kmer>>=5;
		}
		return sb.reverse().toString();
	}

	/**
	 * Converts a codon number to its string representation.
	 * @param codon Integer representation of the codon (0-63)
	 * @return Three-character codon string, or "NNN" if invalid
	 */
	public static final String codonToString(int codon){
		return codon>=0 && codon<codonToString.length ? codonToString[codon] : "NNN";
	}

	/** Gets the canonical (first) codon for this amino acid.
	 * @return The primary codon string for this amino acid, or null if none */
	public String canonicalCodon(){
		return codeStrings==null || codeStrings.length<1 ? null : codeStrings[0];
	}

	/** Full name of the amino acid (e.g., "Alanine") */
	public final String name;
	/** Three-letter symbol code for the amino acid (e.g., "Ala") */
	public final String symbol;
	/** Single-letter code for the amino acid (e.g., 'A') */
	public final char letter;
	/** Exposed mutable codon strings; standard descriptors are converted to DNA spelling during initialization. */
	public final String[] codeStrings;

	//a=1
	//c=2
	//g=4
	//t=8

	//	R 	G A (puRine)
	//	Y 	T C (pYrimidine)
	//	K 	G T (Ketone)
	//	M 	A C (aMino group)
	//	S 	G C (Strong interaction)
	//	W 	A T (Weak interaction)
	//	B 	G T C (not A) (B comes after A)
	//	D 	G A T (not C) (D comes after C)
	//	H 	A C T (not G) (H comes after G)
	//	V 	G C A (not T, not U) (V comes after U)
	//	N 	A G C T (aNy)
	//	X 	masked
	//	- 	gap of indeterminate length

	/** First-codon snapshot by amino index, captured in RNA spelling before codeStrings is converted to DNA. */
	public static final String[] canonicalCodons=new String[21];

	/** Converts base number (0-4) to ASCII character: A, C, G, T, N */
	public static final byte[] numberToBase={
		'A', 'C', 'G', 'T', 'N'
	};

	/** Converts amino acid number (0-20) to single-letter amino acid code */
	public static final byte[] numberToAcid=new byte[21];

	/** Converts base number (0-4) to complement ASCII character: T, G, C, A, N */
	public static final byte[] numberToComplementaryBase={
		'T', 'G', 'C', 'A', 'N'
	};

	/** Converts base number (0-4) to complement number: 3, 2, 1, 0, 4 */
	public static final byte[] numberToComplement={
		3, 2, 1, 0, 4
	};

	/** Extended base number to ASCII mapping including IUPAC ambiguity codes */
	public static final byte[] numberToBaseExtended={
		' ', 'A', 'C', 'M', 'G', 'R', 'S', 'V', //0-7
		'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N', //8-15
		'X', ' ', ' ', ' ', ' ', ' ', ' ', ' ', //16-23
	};

	/** Has 'N' in position 0.  Mainly for translating compressed arrays containing zeroes to bases. */
	public static final byte[] numberToBaseExtended2={
		'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V', //0-7
		'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N', //8-15
		'X', ' ', ' ', ' ', ' ', ' ', ' ', ' ', //16-23
	};

	/** Mapping for IUPAC degenerate base codes representing multiple bases */
	public static final byte[] degenerateBases={
		' ', ' ', ' ', 'M', ' ', 'R', 'S', 'V', //0-7
		' ', 'W', 'Y', 'H', 'K', 'D', 'B', ' ', //8-15
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', //16-23
	};

	/** Extended complement mapping including IUPAC ambiguity codes */
	public static final byte[] numberToComplementaryBaseExtended={
		' ', 'T', 'G', 'K', 'C', 'Y', 'S', 'B', //0-7
		'A', 'W', 'R', 'D', 'M', 'H', 'V', 'N', //8-15
		'X', ' ', ' ', ' ', ' ', ' ', ' ', ' ', //16-23
	};

	/** Element i is: N-bit code for a symbol, -1 otherwise */
	public static final byte[] symbolToNumber(boolean amino){
		return amino ? acidToNumber : baseToNumber;
	}

	/** Element i is: N-bit code for a symbol, 0 otherwise */
	public static final byte[] symbolToNumber0(boolean amino){
		return amino ? acidToNumber0 : baseToNumber0;
	}

	/** Element i is: N-bit code for a symbol, -1 otherwise */
	public static final byte[] symbolToComplementNumber(boolean amino){
		return amino ? acidToNumber : baseToComplementNumber;
	}

	/** Element i is: N-bit code for a symbol, 0 otherwise */
	public static final byte[] symbolToComplementNumber0(boolean amino){
		return amino ? acidToNumber0 : baseToComplementNumber0;
	}

	/** Amino-name-order index for standard residue letters in either case and stop; otherwise -1. */
	public static final byte[] acidToNumber=new byte[128];

	/** Standard residue indices from acidToNumber, with stop and unsupported symbols mapped to -1. */
	public static final byte[] acidToNumberNoStops=new byte[128];

	/** Standard residue/stop indices from acidToNumber, with unsupported ASCII entries mapped to zero. */
	public static final byte[] acidToNumber0=new byte[128];//Rename acidToNumber0

	/** Standard indices plus a shared unknown index for X/B/J/Z/U/O/dot and a separate dash index.
	 * Letters accept either case; other unsupported ASCII entries are -1. */
	public static final byte[] acidToNumberExtended=new byte[128];

	/** Uppercase amino grouping into eight classes 0..7; unsupported/lowercase symbols are -1. */
	public static final byte[] acidToNumber8=new byte[128];

	/** A/a=0, C/c=1, G/g=2, T/t/U/u=3; all other ASCII entries are -1. */
	public static final byte[] baseToNumber=new byte[128];

	/** A/a=0, C/c=1, G/g=2, T/t/U/u=3; unknown ASCII entries are zero. */
	public static final byte[] baseToNumber0=new byte[128];

	/** A/a=0, C/c=1, G/g=2, T/t/U/u=3; unknown ASCII entries are three. */
	public static final byte[] baseToNumber3=new byte[128];

	/** A/a=0, C/c=1, G/g=2, T/t/U/u=3; unknown ASCII entries are four. */
	public static final byte[] baseToNumber4=new byte[128];

	/** A/a=3, C/c=2, G/g=1, T/t/U/u=0; all other ASCII entries are -1. */
	public static final byte[] baseToComplementNumber=new byte[128];

	/** A/a=3, C/c=2, G/g=1, T/t/U/u=0; all other ASCII entries are zero. */
	public static final byte[] baseToComplementNumber0=new byte[128];

	/** A/a=0, C/c=1, G/g=2, T/t/U/u=3, N/n=4; all other ASCII entries are -1. */
	public static final byte[] baseToNumberACGTN=new byte[128];

	/** A/a/N/n=0, C/c=1, G/g=2, T/t/U/u=3; all other ASCII entries are -1. */
	public static final byte[] baseToNumberACGTN2=new byte[128];

	/** A/a=0, C/c=1, G/g=2, T/t/U/u=3; all other ASCII entries are four. */
	public static final byte[] baseToNumberACGTother=new byte[128];

	/** Seven-bit hash codes for A/C/G/T in either case; U/u uses T's code, others zero. */
	public static final byte[] baseToHashcode=new byte[128];

	/** Letters in either case map to 1..26, star to 27 and dash to 28; all others to 31.
	 * ASCII NUL also maps to 31; zero is unused by this lookup.
	 * Allows fast alphabetical sorting of sequences. */
	public static final byte[] symbolTo5Bit=new byte[128];

	/** Uppercases A/C/G/T, maps U/u to T, and maps all other supported ASCII entries to N. */
	public static final byte[] baseToACGTN=new byte[128];

	/**
	 * Extended base complement mapping including IUPAC codes and special characters
	 */
	public static final byte[] baseToComplementExtended=new byte[128];

	/** Array mapping codon numbers (0-63) to three-character codon strings */
	public static final String[] codonToString=new String[64];

	/** Uracil to Thymine, everything else unchanged */
	public static final byte[] uToT=new byte[256];
	/** Thymine to Uracil, everything else unchanged */
	public static final byte[] tToU=new byte[256];
	/** Dot, dash, X/x and lowercase n map to N; other byte values are unchanged. */
	public static final byte[] dotDashXToNocall=new byte[256];
	/** Dot, dash and X/x map to X; other byte values are unchanged. */
	public static final byte[] dotDashXToNocallAA=new byte[256];
	/** Letters to uppercase, everything else unchanged */
	public static final byte[] toUpperCase=new byte[256];
	/** Lowercase to N, everything else unchanged */
	public static final byte[] lowerCaseToNocall=new byte[256];
	/** Lowercase to ., everything else unchanged */
	public static final byte[] lowerCaseToNocallAA=new byte[256];
	/** Recognized non-defined extended symbols map to N, including N/n and X/x.
	 * A/C/G/T/U in either case and unrecognized byte values remain unchanged. */
	public static final byte[] iupacToNocall=new byte[256];

	/** Element i is the bitwise OR of constituent IUPAC base numbers in baseToNumber.<br>
	 * For example, {@code baseToNumberExtended['M'] = (1<<baseToNumber['A']) | (1<<baseToNumber['C']) = 3}.<br>
	 * Invalid characters are -1 */
	public static final byte[] baseToNumberExtended=new byte[128];
	/** Extended base to complement number mapping for IUPAC codes */
	public static final byte[] baseToComplementNumberExtended=new byte[128];
	/** Twenty standard residues in alphabetical amino-acid name order, followed by END. */
	public static final AminoAcid[] AlphabeticalAAs=new AminoAcid[21];
	/** Standard codons 0..63, unused null slot 64 and ANY at 65. */
	public static final AminoAcid[] codeToAA=new AminoAcid[66];
	/** Standard codon letters 0..63, unused zero slot 64 and X at 65. */
	public static final char[] codeToChar=new char[66];
	/** Byte counterpart of codeToChar, including the unused zero slot 64. */
	public static final byte[] codeToByte=new byte[66];
	/** Canonical codon 0..63 by amino ASCII symbol; unresolved extensions use 65, others -1. */
	public static final byte[] aminoToCode=new byte[128];
	/**
	 * Hash map for looking up AminoAcid objects by name, symbol, letter, or codon string
	 */
	public static final HashMap<String, AminoAcid> stringToAA=new HashMap<String, AminoAcid>(512);

	/** Alanine amino acid constant with its codons GCU, GCC, GCA, GCG */
	public static final AminoAcid Alanine=new AminoAcid("Alanine, Ala, A, GCU, GCC, GCA, GCG");
	/** Arginine amino acid constant with its codons CGU, CGC, CGA, CGG, AGA, AGG */
	public static final AminoAcid Arginine=new AminoAcid("Arginine, Arg, R, CGU, CGC, CGA, CGG, AGA, AGG");
	/** Asparagine amino acid constant with its codons AAU, AAC */
	public static final AminoAcid Asparagine=new AminoAcid("Asparagine, Asn, N, AAU, AAC");
	/** Aspartic acid amino acid constant with its codons GAU, GAC */
	public static final AminoAcid AsparticAcid=new AminoAcid("AsparticAcid, Asp, D, GAU, GAC");
	/** Cysteine amino acid constant with its codons UGU, UGC */
	public static final AminoAcid Cysteine=new AminoAcid("Cysteine, Cys, C, UGU, UGC");
	/** Glutamic acid amino acid constant with its codons GAA, GAG */
	public static final AminoAcid GlutamicAcid=new AminoAcid("GlutamicAcid, Glu, E, GAA, GAG");
	/** Glutamine amino acid constant with its codons CAA, CAG */
	public static final AminoAcid Glutamine=new AminoAcid("Glutamine, Gln, Q, CAA, CAG");
	/** Glycine amino acid constant with its codons GGU, GGC, GGA, GGG */
	public static final AminoAcid Glycine=new AminoAcid("Glycine, Gly, G, GGU, GGC, GGA, GGG");
	/** Histidine amino acid constant with its codons CAU, CAC */
	public static final AminoAcid Histidine=new AminoAcid("Histidine, His, H, CAU, CAC");
	/** Isoleucine amino acid constant with its codons AUU, AUC, AUA */
	public static final AminoAcid Isoleucine=new AminoAcid("Isoleucine, Ile, I, AUU, AUC, AUA");
	/** Leucine amino acid constant with its codons UUA, UUG, CUU, CUC, CUA, CUG */
	public static final AminoAcid Leucine=new AminoAcid("Leucine, Leu, L, UUA, UUG, CUU, CUC, CUA, CUG");
	/** Lysine amino acid constant with its codons AAA, AAG */
	public static final AminoAcid Lysine=new AminoAcid("Lysine, Lys, K, AAA, AAG");
	/** Methionine amino acid constant with its codon AUG (start codon) */
	public static final AminoAcid Methionine=new AminoAcid("Methionine, Met, M, AUG");
	/** Phenylalanine amino acid constant with its codons UUU, UUC */
	public static final AminoAcid Phenylalanine=new AminoAcid("Phenylalanine, Phe, F, UUU, UUC");
	/** Proline amino acid constant with its codons CCU, CCC, CCA, CCG */
	public static final AminoAcid Proline=new AminoAcid("Proline, Pro, P, CCU, CCC, CCA, CCG");
	/** Serine amino acid constant with its codons UCU, UCC, UCA, UCG, AGU, AGC */
	public static final AminoAcid Serine=new AminoAcid("Serine, Ser, S, UCU, UCC, UCA, UCG, AGU, AGC");
	/** Threonine amino acid constant with its codons ACU, ACC, ACA, ACG */
	public static final AminoAcid Threonine=new AminoAcid("Threonine, Thr, T, ACU, ACC, ACA, ACG");
	/** Tryptophan amino acid constant with its codon UGG */
	public static final AminoAcid Tryptophan=new AminoAcid("Tryptophan, Trp, W, UGG");
	/** Tyrosine amino acid constant with its codons UAU, UAC */
	public static final AminoAcid Tyrosine=new AminoAcid("Tyrosine, Tyr, Y, UAU, UAC");
	/** Valine amino acid constant with its codons GUU, GUC, GUA, GUG */
	public static final AminoAcid Valine=new AminoAcid("Valine, Val, V, GUU, GUC, GUA, GUG");

	/**
	 * Selenocysteine amino acid constant (21st amino acid, sometimes coded by UGA)
	 */
	public static final AminoAcid Selenocysteine=new AminoAcid("Selenocysteine, Sec, U"); //UGA sometimes
	/** Pyrrolysine amino acid constant (22nd amino acid found in some archaea) */
	public static final AminoAcid Pyrrolysine=new AminoAcid("Pyrrolysine, Pyl, O");

	/** Stop codon amino acid constant with codons UAA, UGA, UAG */
	public static final AminoAcid END=new AminoAcid("End, End, *, UAA, UGA, UAG");
	/** Unknown/any amino acid constant represented by X */
	public static final AminoAcid ANY=new AminoAcid("Any, Any, X, XXX");

	/** Number of bits used to encode amino acids (5 bits for 20+ amino acids) */
	public static int AMINO_SHIFT=5;

	/**
	 * Color-space encoding matrix for SOLiD sequencing technology base transitions
	 */
	public static final byte[][] COLORS=new byte[][]{
		{0, 1, 2, 3},
		{1, 0, 3, 2},
		{2, 3, 0, 1},
		{3, 2, 1, 0}
	};

	/**
	 * Checks reverse-complement equality without regard to letter case.
	 * Both arrays must be nonnull and contain ASCII symbols supported by the complement table.
	 * @param bases First DNA sequence
	 * @param bases2 Second DNA sequence
	 * @return true if bases2 is the reverse complement of bases
	 */
	public static boolean equalsReverseComp(byte[] bases, byte[] bases2){
		if(bases.length!=bases2.length){return false;}
		for(int i=0, j=bases2.length-1; i<bases.length; i++, j--){
			byte a=Tools.toUpperCase(bases[i]), b=Tools.toUpperCase(baseToComplementExtended[bases2[j]]);
			if(a!=b){return false;}
		}
		return true;
	}

	/** Returns a new reverse-complemented array using the case-preserving ASCII table.
	 * Requires nonnull input. Unsupported ASCII symbols retain their identity mapping. */
	public static final byte[] reverseComplementBases(final byte[] in){
		byte[] out=new byte[in.length];
		final int last=in.length-1;
		for(int i=0; i<in.length; i++){
			out[i]=baseToComplementExtended[in[last-i]];
		}
		return out;
	}

	/**
	 * Reverse-complements the entire ASCII array in place; null is a no-op.
	 */
	public static final void reverseComplementBasesInPlace(final byte[] in){
		if(in!=null){reverseComplementBasesInPlace(in, in.length);}
	}
	/**
	 * Complements DNA sequence in place without reversing.
	 * Changes {@code A<->T} and {@code C<->G} in the input array.
	 * @param in DNA sequence to be complemented in place
	 */
	public static final void complementBasesInPlace(final byte[] in){
		if(in==null){return;}
		complementBasesInPlace(in, in.length);
	}
	/**
	 * Complements specified length of DNA sequence in place.
	 * Changes {@code A<->T} and {@code C<->G} for the first 'length' bases.
	 * @param in DNA sequence to be complemented
	 * @param length Number of bases to complement
	 */
	public static final void complementBasesInPlace(final byte[] in, final int length){
		if(in==null){return;}
		for(int i=0; i<length; i++){
			in[i]=baseToComplementExtended[in[i]];
		}
	}

	/**
	 * Reverse complements specified length of DNA sequence in place.
	 * Efficiently handles both reversal and complementation simultaneously.
	 * @param in DNA sequence to be reverse complemented
	 * @param length Number of bases to process
	 */
	public static final void reverseComplementBasesInPlace(final byte[] in, final int length){
		if(in==null){return;}
		final int last=length-1;
		final int max=length/2;
		for(int i=0; i<max; i++){
			byte a=in[i];
			byte b=in[last-i];
			in[i]=baseToComplementExtended[b];
			in[last-i]=baseToComplementExtended[a];
		}
		if((length&1)==1){//Odd length; process middle
			in[max]=baseToComplementExtended[in[max]];
		}
	}

	/**
	 * Creates reverse complement of a DNA string.
	 * @param in Input DNA sequence string
	 * @return Reverse complement string, or null if input is null
	 */
	public static final String reverseComplementBases(String in){
		return in==null ? null : new String(reverseComplementBases(in.getBytes()));
	}

	/**
	 * Computes reverse complement of binary k-mer representation.
	 * Uses bitwise operations for efficient k-mer manipulation.
	 *
	 * @param kmer Binary k-mer (32-bit integer)
	 * @param k Length of the k-mer
	 * @return Reverse complement as binary k-mer
	 */
	public static final int reverseComplementBinary(int kmer, int k){
		int out=0;
		kmer=~kmer;
		for(int i=0; i<k; i++){
			out=((out<<2)|(kmer&3));
			kmer>>=2;
		}
		return out;
	}

	/**
	 * Computes reverse complement of binary k-mer representation.
	 * Uses bitwise operations for efficient k-mer manipulation with 64-bit precision.
	 *
	 * @param kmer Binary k-mer (64-bit long)
	 * @param k Length of the k-mer
	 * @return Reverse complement as binary k-mer
	 */
	public static final long reverseComplementBinary(long kmer, int k){
		long out=0;
		kmer=~kmer;
		for(int i=0; i<k; i++){
			out=((out<<2)|(kmer&3L));
			kmer>>=2;
		}
		return out;
	}

	/**
	 * Fast reverse complement using lookup table for 4-base chunks.
	 * More efficient than bit-by-bit processing for longer k-mers.
	 *
	 * @param kmer Binary k-mer (32-bit integer)
	 * @param k Length of the k-mer
	 * @return Reverse complement as binary k-mer
	 */
	public static final int reverseComplementBinaryFast(int kmer, int k){
		int out=0;
		int extra=k&3;
		for(int i=0; i<extra; i++){
			out=((out<<2)|((~kmer)&3));
			kmer>>=2;
		}
		k-=extra;
		for(int i=0; i<k; i+=4){
			out=((out<<8)|(rcompBinaryTable[kmer&0xFF]));
			kmer>>=8;
		}
		return out;
	}

	/** Reverse-complements the low 2*k bits using a four-base lookup, for k in 0..32. */
	public static final long reverseComplementBinaryFastLookup(long kmer, int k){
		long out=0;
		int extra=k&3;
		for(int i=0; i<extra; i++){
			out=((out<<2)|((~kmer)&3L));
			kmer>>=2;
		}
		k-=extra;
		for(int i=0; i<k; i+=4){
			out=((out<<8)|(rcompBinaryTable[(int)(kmer&0xFFL)]));
			kmer>>=8;
		}
		return out;
	}

	//70% faster at K=31
	/**
	 * Reverses packed two-bit bases by fixed-width bit swaps after complementing.
	 * This overload uses no lookup table. Requires k in 1..32: k=0 would trigger
	 * Java's masked shift count and retain a full-width result instead of an empty k-mer.
	 *
	 * @param kmer Binary k-mer (64-bit long)
	 * @param k Length of the k-mer
	 * @return Reverse complement as binary k-mer
	 */
	public static final long reverseComplementBinaryFast(long kmer, int k){
		// Complement first
		long x=~kmer;

		// Swap adjacent 2-bit pairs, then 4-bit, then 8-bit, etc
		x=((x&0x3333333333333333L)<<2)|((x&0xCCCCCCCCCCCCCCCCL)>>>2);
		x=((x&0x0F0F0F0F0F0F0F0FL)<<4)|((x&0xF0F0F0F0F0F0F0F0L)>>>4);
		x=((x&0x00FF00FF00FF00FFL)<<8)|((x&0xFF00FF00FF00FF00L)>>>8);
		x=((x&0x0000FFFF0000FFFFL)<<16)|((x&0xFFFF0000FFFF0000L)>>>16);
		x=(x<<32)|(x>>>32);

		// Right-align for k<32
		x=x>>>(2*(32-k));

		return x;
	}

	/**
	 * Converts consecutive bases to color-space encoding.
	 * Used for SOLiD sequencing color-space representation.
	 *
	 * @param base1 First base
	 * @param base2 Second base
	 * Legacy behavior maps one undefined base to T; only two undefined bases return N.
	 * @return Numeric color 0..3, or ASCII N when both bases are undefined
	 */
	public static final byte baseToColor(byte base1, byte base2){
		byte a=baseToNumber[base1];
		byte b=baseToNumber[base2];
		if(a<0 && b<0){return 'N';}
		if(a<0){a=3;}
		if(b<0){b=3;}
		return COLORS[a][b];
	}

	/**
	 * Converts color-space encoding back to base given reference base.
	 * Reverses color-space encoding using known reference base.
	 *
	 * @param base1 Reference base
	 * @param color Color code (0-3)
	 * @return Decoded base or 'N' if invalid
	 */
	public static final byte colorToBase(byte base1, byte color){
		if(!isFullyDefined(base1) || color<0 || color>3){
			return (byte)'N';
		}
		byte a=baseToNumber[base1];

		return numberToBase[COLORS[a][color]];
	}

	//	public static final byte toNumber(String code){
	//		return toNumber(code.charAt(0), code.charAt(1), code.charAt(2));
	//	}

	/**
	 * Converts three-character codon string to amino acid.
	 * Reads the first three characters; legacy character translation maps N to A.
	 * @param code ASCII string with at least three characters
	 * @return Corresponding AminoAcid object
	 */
	public static final AminoAcid toAA(String code){
		return toAA(code.charAt(0), code.charAt(1), code.charAt(2));
	}

	/**
	 * Converts three-character codon string to amino acid single-letter code.
	 * Reads the first three characters; legacy character translation maps N to A.
	 * @param code ASCII string with at least three characters
	 * @return Single-letter amino acid code
	 */
	public static final char toChar(String code){
		return toChar(code.charAt(0), code.charAt(1), code.charAt(2));
	}

	/**
	 * Splits IUPAC ambiguous base into constituent bases.
	 * For example, 'M' (A or C) returns ['A', 'C'].
	 * @param c ASCII base with extended mask 1..15, including U as T; excludes X
	 * @return Constituent bases in A/C/G/T bit order, using the current numberToBase alphabet
	 * @throws IllegalArgumentException if c is outside the table or its mask is outside 1..15
	 */
	public static final char[] splitBase(char c){
		//DNA-043: Only the four canonical-base bits can be expanded by the loop below.
		final int b=(c<baseToNumberExtended.length ? baseToNumberExtended[c] : -1);
		if(b<1 || b>15){throw new IllegalArgumentException("Cannot split character code "+(int)c+
				": expected a nonempty four-base mask (1..15), found "+b);}
		int len=Integer.bitCount(b);
		char[] out=new char[len];

		int index=0;
		for(int i=0; i<4; i++){
			if(((1<<i)&b)!=0){
				out[index]=(char)numberToBase[i];
				index++;
			}
		}
		return out;
	}

	/**
	 * Converts integer code to n-length DNA sequence.
	 * Each base occupies 2 bits in the integer representation.
	 *
	 * @param code Integer encoding of the sequence
	 * @param n Number of low-order two-bit bases to decode, normally 0..16
	 * @return DNA sequence as byte array
	 */
	public static final byte[] numberToBases(int code, int n){

		byte[] bytes=KillSwitch.allocByte1D(n);

		for(int i=n-1; i>=0; i--){
			int temp=code&3;
			code>>=2;
			bytes[i]=numberToBase[temp];
		}

		return bytes;
	}

	/**
	 * Converts DNA sequence to integer code.
	 * Each base occupies 2 bits: A=0, C=1, G=2, T=3.
	 * @param tuple ASCII bases, at most 16 to avoid discarding high bits; accepts U as T
	 * @return Packed code or -1 for an undefined base; a valid 16-base code may also equal -1
	 */
	public static final int baseTupleToNumber(byte[] tuple){

		int r=0;
		for(int i=0; i<tuple.length; i++){
			int temp=baseToNumberACGTN[tuple[i]];
			if(temp<0 || temp>3){return -1;}
			r=((r<<2)|temp);
		}

		return r;
	}

	/**
	 * Checks the defined-base table, accepting A/C/G/T/U in either case.
	 * @param base ASCII table index
	 * @return Whether the base maps to 0..3
	 */
	public static boolean isFullyDefined(char base){
		return baseToNumber[base]>=0;
	}

	/**
	 * Checks the defined-base table, accepting A/C/G/T/U in either case.
	 * @param base Byte to check
	 * @return Whether the base maps to 0..3; negative bytes return false
	 */
	public static boolean isFullyDefined(byte base){
		return base>=0 && baseToNumber[base]>=0;
	}

	/**
	 * Checks whether a byte represents one of the 20 standard residues or stop.
	 * @param acid Amino acid byte to check
	 * @return True for a defined residue in either case or '*'; negative bytes return false
	 */
	public static boolean isFullyDefinedAA(byte acid){
		return acid>=0 && acidToNumber[acid]>=0;
	}

	/**
	 * Checks if amino acid byte represents standard amino acid excluding stops.
	 * @param acid Amino acid byte to check
	 * @return true if acid is a standard amino acid but not a stop codon
	 */
	public static boolean isFullyDefinedAANoStops(byte acid){
		return acid>=0 && acidToNumberNoStops[acid]>=0;
	}

	/**
	 * Checks for A/C/G/T/U/N in either case; empty strings pass.
	 * @param bases Nonnull ASCII string
	 * @return Whether every character has a nonnegative ACGTN-table entry
	 */
	public static boolean isACGTN(String bases){
		for(int i=0; i<bases.length(); i++){
			if(!isACGTN(bases.charAt(i))){return false;}
		}
		return true;
	}

	/**
	 * Checks for A/C/G/T/U/N in either case; empty arrays pass, negative bytes fail.
	 * @param bases Nonnull byte array
	 * @return Whether every byte is accepted by the ACGTN table
	 */
	public static boolean isACGTN(byte[] bases){
		for(int i=0; i<bases.length; i++){
			if(!isACGTN(bases[i])){return false;}
		}
		return true;
	}

	/**
	 * Checks for A/C/G/T/U/N in either case.
	 * @param base ASCII table index
	 * @return Whether the table entry is nonnegative
	 */
	public static boolean isACGTN(char base){
		return baseToNumberACGTN[base]>=0;
	}

	/**
	 * Checks for A/C/G/T/U/N in either case.
	 * @param base Byte to check
	 * @return Whether the table entry is nonnegative; negative bytes return false
	 */
	public static boolean isACGTN(byte base){
		return base>=0 && baseToNumberACGTN[base]>=0;
	}

	/**
	 * Checks for A/C/G/T/U/N in either case.
	 * Returns true for null or empty strings.
	 * @param s ASCII string to validate, or null
	 * @return True for null/empty input or only accepted nucleotide symbols
	 */
	public static boolean containsOnlyACGTN(String s){
		if(s==null || s.length()==0){return true;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(baseToNumberACGTN[c]<0){return false;}
		}
		return true;
	}

	/**
	 * Checks for A/C/G/T/U/N in either case, plus question-mark placeholders.
	 * @param s ASCII string to validate, or null
	 * @return True for null/empty input or only accepted symbols
	 */
	public static boolean containsOnlyACGTNQ(String s){
		if(s==null || s.length()==0){return true;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c!='?' && baseToNumberACGTN[c]<0){return false;}
		}
		return true;
	}

	/**
	 * Checks for A/C/G/T/U/N in either case; negative bytes fail.
	 * Returns true for null or empty arrays.
	 * @param array Byte array to validate
	 * @return True for null/empty input or only accepted nucleotide symbols
	 */
	public static boolean containsOnlyACGTN(byte[] array){
		if(array==null || array.length==0){return true;}
		for(int i=0; i<array.length; i++){
			byte b=array[i];
			if(b<0 || baseToNumberACGTN[b]<0){return false;}
		}
		return true;
	}

	/**
	 * Checks for defined A/C/G/T/U bases in either case; empty strings pass.
	 * @param s Nonnull ASCII string
	 * @return Whether every character maps to 0..3
	 */
	public static boolean isFullyDefined(String s){
		for(int i=0; i<s.length(); i++){
			if(!isFullyDefined(s.charAt(i))){return false;}
		}
		return true;
	}

	/**
	 * Checks for defined A/C/G/T/U bases in either case; negative bytes fail.
	 * @param s Nonnull byte array; empty arrays pass
	 * @return Whether every byte maps to 0..3
	 */
	public static boolean isFullyDefined(byte[] s){
		for(int i=0; i<s.length; i++){
			if(!isFullyDefined(s[i])){return false;}
		}
		return true;
	}

	/**
	 * Counts number of undefined bases in sequence.
	 * Undefined means outside A/C/G/T/U in either case, including negative bytes.
	 * @param s Nonnull sequence
	 * @return Count of undefined bases
	 */
	public static int countUndefined(byte[] s){
		int x=0;
		for(int i=0; i<s.length; i++){
			if(!isFullyDefined(s[i])){x++;}
		}
		return x;
	}

	/**
	 * Counts number of defined bases in sequence.
	 * Defined bases are A/C/G/T/U in either case; negative bytes are undefined.
	 * @param s DNA sequence to analyze
	 * @return Count of defined bases, or 0 if array is null
	 */
	public static int countDefined(byte[] s){
		if(s==null){return 0;}
		int x=0;
		for(int i=0; i<s.length; i++){
			if(isFullyDefined(s[i])){x++;}
		}
		return x;
	}

	/**
	 * Converts three-character DNA string to numeric codon code.
	 * Each base occupies 2 bits in the result.
	 * @param s Three-character DNA string
	 * @return Numeric codon code (0-63), or -1 if invalid bases
	 */
	public static final byte toNumber(String s){
		assert(s.length()==3);
		int num=0;
		for(int i=0; i<3; i++){
			char c=s.charAt(i);
			int x=baseToNumber[c];
			if(x<0){return (byte)-1;}
			num=(num<<2)|x;
		}
		return (byte)num;
	}

	/**
	 * Converts three characters to numeric codon code.
	 * Unlike the String overload, this legacy method maps N to A; U maps to T.
	 *
	 * @param c1 First base character
	 * @param c2 Second base character
	 * @param c3 Third base character
	 * @return Numeric codon code (0-63)
	 */
	public static final byte toNumber(char c1, char c2, char c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		int x=(baseToNumberACGTN2[c1]<<4)|(baseToNumberACGTN2[c2]<<2)|(baseToNumberACGTN2[c3]);
		return (byte)x;
	}

	/**
	 * Converts three base characters to corresponding amino acid.
	 * Legacy ASCII character translation treats N as A and U as T, accepting either case.
	 *
	 * @param c1 First base character
	 * @param c2 Second base character
	 * @param c3 Third base character
	 * @return AminoAcid object for the codon
	 */
	public static final AminoAcid toAA(char c1, char c2, char c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		int x=(baseToNumberACGTN2[c1]<<4)|(baseToNumberACGTN2[c2]<<2)|(baseToNumberACGTN2[c3]);
		return codeToAA[x];
	}

	/**
	 * Converts three ASCII bases to a letter, treating N as A and U as T.
	 *
	 * @param c1 First base character
	 * @param c2 Second base character
	 * @param c3 Third base character
	 * @return Single-letter amino acid code
	 */
	public static final char toChar(char c1, char c2, char c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		int x=(baseToNumberACGTN2[c1]<<4)|(baseToNumberACGTN2[c2]<<2)|(baseToNumberACGTN2[c3]);
		return codeToChar[x];
	}

	/**
	 * Converts three base bytes to amino acid byte code.
	 * Returns X for undefined ASCII bases, including N; accepts U and either case.
	 * Negative bytes are not valid table indices and are not guarded here.
	 *
	 * @param c1 First base byte
	 * @param c2 Second base byte
	 * @param c3 Third base byte
	 * @return Amino acid byte code or 'X' for invalid codons
	 */
	public static final byte toByte(byte c1, byte c2, byte c3){
		int a=baseToNumber[c1], b=baseToNumber[c2], c=baseToNumber[c3];
		if(a<0 || b<0 || c<0){return (byte)'X';}
		int x=((a<<4)|(b<<2)|c);
		return codeToByte[x];
	}

	/**
	 * Converts three base bytes to amino acid character.
	 * NOTE: N maps to A here (baseToNumberACGTN2['N']=0), so N-bearing codons translate
	 * as if N were A; the n>3 '?' guard below is unreachable for valid ACGTN input (all
	 * table values are 0-3 or -1, and -1 is caught by the assert). Currently has no
	 * callers; use toByte for N-aware translation (it returns 'X' for N codons).
	 *
	 * @param c1 First base byte
	 * @param c2 Second base byte
	 * @param c3 Third base byte
	 * @return Amino acid character (N treated as A)
	 */
	public static final char toChar(byte c1, byte c2, byte c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		byte n1=baseToNumberACGTN2[c1], n2=baseToNumberACGTN2[c2], n3=baseToNumberACGTN2[c3];
		if(n1>3 || n2>3 || n3>3){return '?';}
		int x=(n1<<4)|(n2<<2)|(n3);
		//		return (x<codeToChar.length ? codeToChar[x] : '?');
		return codeToChar[x];
	}

	/**
	 * Translates DNA string to amino acid sequence string.
	 * Processes complete triplets beginning at offset zero, ignoring trailing 1–2 bases.
	 * Uses character translation, so N is treated as A rather than producing X.
	 * @param bases DNA sequence string
	 * @return Amino acid sequence string
	 */
	public static final String stringToAAs(String bases){
		StringBuilder sb=new StringBuilder(bases.length()/3);
		for(int i=2; i<bases.length(); i+=3){
			char a=toAA(bases.charAt(i-2), bases.charAt(i-1), bases.charAt(i)).letter;
			sb.append(a);
		}
		return sb.toString();
	}

	/**
	 * Translates DNA sequence to amino acids in all six reading frames.
	 * Returns frames 0–2 for the forward strand, 3–5 for its reverse complement.
	 * Null/short input yields six null entries; frames with no complete codon are null.
	 * @param bases DNA sequence to translate
	 * @return Array of 6 amino acid sequences, one per reading frame
	 */
	public static final byte[][] toAAsSixFrames(byte[] bases){
		byte[][] out=new byte[6][];
		if(bases!=null && bases.length>2){
			for(int i=0; i<3; i++){
				out[i]=toAAs(bases, i);
			}
			byte[] rcomp=reverseComplementBases(bases);
			for(int i=0; i<3; i++){
				out[i+3]=toAAs(rcomp, i);
			}
		}
		return out;
	}

	/**
	 * Computes quality scores for six-frame amino acid translation.
	 * Converts the product of three base-correctness probabilities back to Phred.
	 * Temporarily reverses quals in place for reverse frames, restoring on normal
	 * completion; callers must exclude concurrent access. Offset addition uses byte
	 * arithmetic without clamping. Null/short input yields six null entries.
	 *
	 * @param quals Quality scores for DNA bases
	 * @param offset Quality score offset to apply to results
	 * @return Array of 6 quality arrays corresponding to six reading frames
	 */
	public static final byte[][] toQualitySixFrames(byte[] quals, int offset){
		byte[][] out=new byte[6][];
		if(quals!=null && quals.length>2){
			for(int i=0; i<3; i++){
				out[i]=toAAQuality(quals, i);
			}
			Vector.reverseInPlace(quals);
			for(int i=0; i<3; i++){
				out[i+3]=toAAQuality(quals, i);
			}
			Vector.reverseInPlace(quals);
		}

		if(offset!=0){
			for(byte[] array : out){
				if(array!=null){
					for(int i=0; i<array.length; i++){
						array[i]+=offset;
					}
				}
			}
		}

		return out;
	}

	/**
	 * Translates DNA sequence to amino acids in specified reading frame.
	 * Frame 0 starts at position 0, frame 1 at position 1, frame 2 at position 2.
	 *
	 * @param bases DNA sequence to translate
	 * @param frame Reading frame (0, 1, or 2)
	 * @return New amino array, or null for null input or no complete codon in the frame
	 */
	public static final byte[] toAAs(byte[] bases, int frame){
		assert(frame>=0 && frame<3);
		if(bases==null){return null;}
		int blen=bases.length-frame;
		if(blen<3){return null;}
		blen=blen-(blen%3);
		final int stop=frame+blen;
		final int alen=blen/3;

		byte[] out=KillSwitch.allocByte1D(alen);
		for(int i=2+frame, j=0; i<stop; i+=3, j++){
			byte a=toByte(bases[i-2], bases[i-1], bases[i]);
			out[j]=a;
		}
		return out;
	}

	/**
	 * Translates a DNA subsequence to amino acids, OMITTING the terminal codon.
	 * De-facto contract (sole caller is prok.CallGenes on ORF coordinates): start and
	 * stop are INCLUSIVE base positions spanning a full ORF (first base of the start
	 * codon .. last base of the stop codon), and the terminal stop codon is
	 * intentionally not translated. Yields (stop-start-2)/3 = orf.length()/3-1 acids.
	 *
	 * @param bases DNA sequence to translate
	 * @param start Starting position (inclusive; first base of start codon)
	 * @param stop Ending position (inclusive; last base of stop codon)
	 * @return Amino acid sequence for [start,stop] minus the terminal stop codon
	 */
	public static final byte[] toAAs(byte[] bases, int start, int stop){
		if(bases==null){return null;}
		//Note: Confirmed correct. Drops the terminal stop codon (Orf coords are inclusive per
		//Orf.java:32-33), NOT an off-by-one; alen=(stop-start-2)/3=orf.length()/3-1, verified
		//against CallGenes.translate's asserted (r.length()+1)*3==orf.length(). Do not "fix" to i<=stop.
		stop-=2;
		final int blen=stop-start;
		final int alen=blen/3;

		byte[] out=KillSwitch.allocByte1D(alen);
		for(int i=2+start, j=0; i<stop; i+=3, j++){
			byte a=toByte(bases[i-2], bases[i-1], bases[i]);
			out[j]=a;
		}
		return out;
	}

	/**
	 * Computes amino acid quality scores for specific reading frame.
	 * Multiplies three correctness probabilities and converts the product back to Phred.
	 *
	 * @param quals Nonnull scores indexing QualityTools.PROB_CORRECT
	 * @param frame Reading frame (0, 1, or 2)
	 * @return New quality array, or null if the frame contains no complete triplet
	 */
	public static final byte[] toAAQuality(byte[] quals, int frame){
		assert(frame>=0 && frame<3);
		int blen=quals.length-frame;
		if(blen<3){return null;}
		blen=blen-(blen%3);
		final int stop=frame+blen;
		final int alen=blen/3;

		byte[] out=KillSwitch.allocByte1D(alen);
		for(int i=2+frame, j=0; i<stop; i+=3, j++){
			byte qa=quals[i-2], qb=quals[i-1], qc=quals[i];
			float pa=QualityTools.PROB_CORRECT[qa], pb=QualityTools.PROB_CORRECT[qb], pc=QualityTools.PROB_CORRECT[qc];
			float p=pa*pb*pc;
			byte q=QualityTools.probCorrectToPhred(p);
			out[j]=q;

			//			System.out.println();
			//			System.out.println(qa+", "+qb+", "+qc+" -> "+q);
			//			System.out.println(pa+", "+pb+", "+pc+" -> "+p);

		}
		//		System.out.println(Arrays.toString(out));
		return out;
	}

	/**
	 * Converts amino acid sequence back to canonical nucleotide codons.
	 * Uses the canonical (first) codon for each amino acid and NNN when no canonical
	 * codon is assigned. This is a representative encoding, not recovery of the original DNA.
	 * Uses cached standard DNA spellings when numberToBase is A/C/G/T. A caller's
	 * nonstandard alphabet is read through a slower fallback, preserving table substitution.
	 * @param aminos ASCII amino acid sequence, or null
	 * @return Newly allocated DNA sequence using canonical codons, or null for null input
	 */
	public static final byte[] toNTs(final byte[] aminos){
		if(aminos==null){return null;}
		final int alen=aminos.length;
		final int blen=alen*3;

		final byte[] out=KillSwitch.allocByte1D(blen);
		if(numberToBase[0]!='A' || numberToBase[1]!='C' || numberToBase[2]!='G' || numberToBase[3]!='T'){
			fillCodonsFromCurrentAlphabet(aminos, out);
			return out;
		}
		for(int i=0, j=0; i<alen; i++, j+=3){
			final int packed=PACKED_CODON_BASES[aminoToCode[aminos[i]]&255];
			out[j]=(byte)(packed>>>16);
			out[j+1]=(byte)(packed>>>8);
			out[j+2]=(byte)packed;
		}
		return out;
	}

	/** Preserves public alphabet substitutions without putting their checks in the usual per-residue loop. */
	private static void fillCodonsFromCurrentAlphabet(byte[] aminos, byte[] out){
		assert(out.length==aminos.length*3) : "Backtranslation output must reserve three bytes per residue";
		for(int i=0, j=0; i<aminos.length; i++, j+=3){
			final int code=aminoToCode[aminos[i]];
			if(code<0 || code>63){out[j]=out[j+1]=out[j+2]='N';}
			else{
				out[j]=numberToBase[(code>>>4)&3];
				out[j+1]=numberToBase[(code>>>2)&3];
				out[j+2]=numberToBase[code&3];
			}
		}
	}

	/**
	 * Tests for a repeating period up to maxRepeat bases, not just one-base homopolymers.
	 * Requires a right-aligned nonnegative encoding with k in 1..15 for this int/shift method.
	 *
	 * @param kmer Binary k-mer representation
	 * @param k Length of k-mer
	 * @param maxRepeat Maximum repeat length to check for
	 * @return true if k-mer contains repeats up to maxRepeat length
	 */
	public static final boolean isHomopolymer(int kmer, int k, int maxRepeat){
		if(maxRepeat<1){return false;}
		final int mask=~((-1)<<(2*k));
		final int inv=(~kmer)&mask;
		boolean polymer=false;
		for(int i=1; i<=maxRepeat && !polymer; i++){
			int shift=2*i;
			polymer=(kmer==(kmer|(kmer>>shift))) && (inv==(inv|(inv>>shift)));
		}
		return polymer;
	}

	/** Lookup table for fast binary reverse complement of 4-base chunks */
	public static final short[] rcompBinaryTable=makeBinaryRcompTable(4);

	/**
	 * Creates lookup table for fast binary reverse complement operations.
	 * Pre-computes reverse complements for all possible k-length bit patterns.
	 * @param k K-mer length for the lookup table
	 * @return Array mapping binary patterns to their reverse complements
	 */
	private static final short[] makeBinaryRcompTable(int k){
		int bits=2*k;
		short[] r=new short[1<<bits];
		for(int i=0; i<r.length; i++){
			r[i]=(short)reverseComplementBinary(i, k);
		}
		return r;
	}

	static{

		for(int i=0; i<uToT.length; i++){uToT[i]=(byte)i;}
		uToT['u']='t';
		uToT['U']='T';

		for(int i=0; i<tToU.length; i++){tToU[i]=(byte)i;}
		tToU['t']='u';
		tToU['T']='U';

		for(int i=0; i<dotDashXToNocall.length; i++){
			dotDashXToNocall[i]=(byte)i;
			iupacToNocall[i]=(byte)i;
		}
		dotDashXToNocall['.']='N';
		dotDashXToNocall['-']='N';
		dotDashXToNocall['X']='N';
		dotDashXToNocall['x']='N';
		dotDashXToNocall['n']='N';

		for(int i=0; i<dotDashXToNocallAA.length; i++){dotDashXToNocallAA[i]=(byte)i;}
		dotDashXToNocallAA['.']='X';
		dotDashXToNocallAA['-']='X';
		dotDashXToNocallAA['X']='X';
		dotDashXToNocallAA['x']='X';

		for(int i=0; i<toUpperCase.length; i++){
			toUpperCase[i]=(byte)((i>='a' && i<='z') ? i-32 : i);
			lowerCaseToNocall[i]=((i>='a' && i<='z') ? (byte)'N' : (byte)i);
			lowerCaseToNocallAA[i]=((i>='a' && i<='z') ? (byte)'.' : (byte)i);
		}

		Arrays.fill(baseToACGTN, (byte)'N');
		Arrays.fill(baseToNumberExtended, (byte)-1);
		for(int i=0; i<numberToBaseExtended.length; i++){
			char x=(char)numberToBaseExtended[i];
			if(!Character.isWhitespace(x)){
				baseToNumberExtended[x]=(byte)i;
				baseToNumberExtended[Tools.toLowerCase(x)]=(byte)i;
			}
		}
		baseToNumberExtended['U']=baseToNumberExtended['u']=baseToNumberExtended['T'];

		Arrays.fill(baseToComplementNumberExtended, (byte)-1);
		for(int i=0; i<numberToComplementaryBaseExtended.length; i++){
			char x=(char)numberToComplementaryBaseExtended[i];
			if(!Character.isWhitespace(x)){
				baseToComplementNumberExtended[x]=(byte)i;
				baseToComplementNumberExtended[Tools.toLowerCase(x)]=(byte)i;
			}
		}
		baseToComplementNumberExtended['U']=baseToComplementNumberExtended['u']=baseToComplementNumberExtended['T'];

		Arrays.fill(baseToNumberACGTN, (byte)-1);
		Arrays.fill(baseToNumberACGTother, (byte)4);
		for(int i=0; i<numberToBase.length; i++){
			char x=(char)numberToBase[i];
			if(!Character.isWhitespace(x)){
				baseToNumberACGTN[x]=baseToNumberACGTother[x]=(byte)i;
				baseToNumberACGTN[Tools.toLowerCase(x)]=baseToNumberACGTother[Tools.toLowerCase(x)]=(byte)i;
				baseToACGTN[x]=baseToACGTN[Tools.toLowerCase(x)]=(byte)x;
			}
		}
		baseToNumberACGTN['U']=baseToNumberACGTN['u']=3;
		baseToNumberACGTother['U']=baseToNumberACGTother['u']=3;
		baseToACGTN['U']=baseToACGTN['u']=(byte)'T';

		for(int i=0; i<baseToNumberACGTN.length; i++){baseToNumberACGTN2[i]=baseToNumberACGTN[i];}
		baseToNumberACGTN2['N']=0;
		baseToNumberACGTN2['n']=0;

		Arrays.fill(baseToNumber, (byte)-1);
		Arrays.fill(baseToNumber0, (byte)0);
		Arrays.fill(baseToNumber3, (byte)3);
		Arrays.fill(baseToNumber4, (byte)4);
		for(int i=0; i<numberToBase.length; i++){
			char x=(char)numberToBase[i];
			if(x=='A' || x=='C' || x=='G' || x=='T'){
				int x2=Tools.toLowerCase(x);
				baseToNumber4[x]=baseToNumber3[x]=baseToNumber0[x]=baseToNumber[x]=(byte)i;
				baseToNumber4[x2]=baseToNumber3[x2]=baseToNumber0[x2]=baseToNumber[x2]=(byte)i;
			}
		}
		baseToNumber4['U']=baseToNumber3['U']=baseToNumber0['U']=baseToNumber['U']=3;
		baseToNumber4['u']=baseToNumber3['u']=baseToNumber0['u']=baseToNumber['u']=3;

		Arrays.fill(baseToComplementNumber, (byte)-1);
		baseToComplementNumber['A']=baseToComplementNumber['a']=3;
		baseToComplementNumber['C']=baseToComplementNumber['c']=2;
		baseToComplementNumber['G']=baseToComplementNumber['g']=1;
		baseToComplementNumber['T']=baseToComplementNumber['t']=0;
		baseToComplementNumber['U']=baseToComplementNumber['u']=0;

		Arrays.fill(baseToComplementNumber0, (byte)0);
		baseToComplementNumber0['A']=baseToComplementNumber0['a']=3;
		baseToComplementNumber0['C']=baseToComplementNumber0['c']=2;
		baseToComplementNumber0['G']=baseToComplementNumber0['g']=1;
		baseToComplementNumber0['T']=baseToComplementNumber0['t']=0;
		baseToComplementNumber0['U']=baseToComplementNumber0['u']=0;

		//Invalid symbols are unchanged.
		//This prevents crashes from -1 being out of bounds, and allows
		//consecutive rcomp operations to restore the original sequence.
		for(int i=0; i<baseToComplementExtended.length; i++){
			baseToComplementExtended[i]=(byte)i;
		}
		//		Arrays.fill(baseToComplementExtended, (byte)-1);
		for(int i=0; i<numberToBaseExtended.length; i++){
			char x=(char)numberToBaseExtended[i];
			char x2=(char)numberToComplementaryBaseExtended[i];
			baseToComplementExtended[x]=(byte)x2;
			baseToComplementExtended[Tools.toLowerCase(x)]=(byte)Tools.toLowerCase(x2);
		}
		baseToComplementExtended['U']=(byte)'A';
		baseToComplementExtended['u']=(byte)'a';
		baseToComplementExtended['?']=(byte)'?';
		baseToComplementExtended[' ']=(byte)' ';
		baseToComplementExtended['-']=(byte)'-';
		baseToComplementExtended['*']=(byte)'*';
		baseToComplementExtended['.']=(byte)'.';

		AlphabeticalAAs[0]=Alanine;
		AlphabeticalAAs[1]=Arginine;
		AlphabeticalAAs[2]=Asparagine;
		AlphabeticalAAs[3]=AsparticAcid;
		AlphabeticalAAs[4]=Cysteine;
		AlphabeticalAAs[5]=GlutamicAcid;
		AlphabeticalAAs[6]=Glutamine;
		AlphabeticalAAs[7]=Glycine;
		AlphabeticalAAs[8]=Histidine;
		AlphabeticalAAs[9]=Isoleucine;
		AlphabeticalAAs[10]=Leucine;
		AlphabeticalAAs[11]=Lysine;
		AlphabeticalAAs[12]=Methionine;
		AlphabeticalAAs[13]=Phenylalanine;
		AlphabeticalAAs[14]=Proline;
		AlphabeticalAAs[15]=Serine;
		AlphabeticalAAs[16]=Threonine;
		AlphabeticalAAs[17]=Tryptophan;
		AlphabeticalAAs[18]=Tyrosine;
		AlphabeticalAAs[19]=Valine;
		AlphabeticalAAs[20]=END;
		//		AlphabeticalAAs[21]=Selenocysteine;
		//		AlphabeticalAAs[22]=Pyrrolysine;
		//		AlphabeticalAAs[23]=ANY;

		Arrays.fill(aminoToCode, (byte)-1);
		Arrays.fill(acidToNumber, (byte)-1);
		Arrays.fill(acidToNumber0, (byte)0);
		Arrays.fill(acidToNumber8, (byte)-1);
		for(int i=0; i<AlphabeticalAAs.length; i++){
			AminoAcid aa=AlphabeticalAAs[i];

			acidToNumber[aa.letter]=(byte)i;
			acidToNumber[Tools.toLowerCase(aa.letter)]=(byte)i;
			acidToNumber0[aa.letter]=(byte)i;
			acidToNumber0[Tools.toLowerCase(aa.letter)]=(byte)i;
			numberToAcid[i]=(byte)aa.letter;
			canonicalCodons[i]=aa.canonicalCodon();

			stringToAA.put(aa.name, aa);
			stringToAA.put(aa.symbol, aa);
			stringToAA.put(aa.letter+"", aa);
			for(int j=0; j<aa.codeStrings.length; j++){
				String s=aa.codeStrings[j];
				stringToAA.put(s, aa);
				aa.codeStrings[j]=s.replace('U', 'T');
				stringToAA.put(aa.codeStrings[j], aa);

				int x=toNumber(s);
				//				System.out.println("x="+x+", aa="+aa);
				codeToAA[x]=aa;
				codeToChar[x]=aa.letter;
				codeToByte[x]=(byte)(aa.letter);
				if(j==0){
					aminoToCode[aa.letter]=(byte)x;
					aminoToCode[Tools.toLowerCase(aa.letter)]=(byte)x;
				}
			}
		}

		for(int i=0; i<acidToNumberNoStops.length; i++){acidToNumberNoStops[i]=acidToNumber[i];}
		acidToNumberNoStops[END.letter]=-1;

		for(int i=0; i<acidToNumber.length; i++){
			acidToNumberExtended[i]=acidToNumber[i];
		}

		{
			byte anySym=(byte)(Tools.max(acidToNumberExtended)+1);
			byte dash=(byte)(anySym+1);
			acidToNumberExtended['x']=acidToNumberExtended['X']=acidToNumberExtended['.']=anySym; //Unknown
			acidToNumberExtended['b']=acidToNumberExtended['B']=anySym;
			acidToNumberExtended['z']=acidToNumberExtended['Z']=anySym;
			acidToNumberExtended['j']=acidToNumberExtended['J']=anySym;
			acidToNumberExtended['u']=acidToNumberExtended['U']=anySym; //Selenocysteine
			acidToNumberExtended['o']=acidToNumberExtended['O']=anySym; //Pyrrolysine
			acidToNumberExtended['-']=dash; //Deletion
		}

		acidToNumber8['H']=acidToNumber8['K']=acidToNumber8['R']=0;
		acidToNumber8['D']=acidToNumber8['E']=1;
		acidToNumber8['S']=acidToNumber8['T']=acidToNumber8['N']=acidToNumber8['Q']=2;
		acidToNumber8['A']=acidToNumber8['V']=acidToNumber8['L']=acidToNumber8['I']=acidToNumber8['M']=3;
		acidToNumber8['F']=acidToNumber8['Y']=acidToNumber8['W']=4;
		acidToNumber8['P']=acidToNumber8['G']=5;
		acidToNumber8['C']=acidToNumber8['*']=6;
		acidToNumber8['B']=acidToNumber8['Z']=7;

		aminoToCode['X']=aminoToCode['x']=aminoToCode['B']=aminoToCode['b']=
			aminoToCode['Z']=aminoToCode['z']=aminoToCode['J']=aminoToCode['j']=
			aminoToCode['O']=aminoToCode['o']=aminoToCode['U']=aminoToCode['u']=65;
		codeToAA[65]=ANY;
		codeToChar[65]='X';
		codeToByte[65]='X';

		stringToAA.put("X", ANY);
		stringToAA.put("Start", Methionine);
		stringToAA.put("Begin", Methionine);
		stringToAA.put("Stop", END);
		stringToAA.put("Aspartic Acid", AsparticAcid);
		stringToAA.put("Glutamic Acid", GlutamicAcid);

		String[] temp=stringToAA.keySet().toArray(new String[0]);

		for(String s : temp){
			AminoAcid aa=stringToAA.get(s);
			assert(aa!=null);
			stringToAA.put(s.toLowerCase(), aa);
		}

		for(int i=0; i<codonToString.length; i++){
			codonToString[i]=kmerToString(i, 3);
		}

		for(int i='A'; i<='z'; i++){
			if(baseToNumber[i]<0 && baseToNumberExtended[i]>=0){
				iupacToNocall[i]='N';
			}
		}

		Arrays.fill(symbolTo5Bit, (byte)31);
		for(int i=0; i<26; i++){symbolTo5Bit['A'+i]=symbolTo5Bit['a'+i]=(byte)(i+1);}
		symbolTo5Bit['*']=27;
		symbolTo5Bit['-']=28;

		baseToHashcode['A']=baseToHashcode['a']=0b1101001;
		baseToHashcode['C']=baseToHashcode['c']=0b0110011;
		baseToHashcode['G']=baseToHashcode['g']=0b0001110;
		baseToHashcode['T']=baseToHashcode['t']=0b1010010;
		baseToHashcode['U']=baseToHashcode['u']=baseToHashcode['T'];
	}

	/** Snapshots canonical DNA spellings after the shared codon-string table is initialized.
	 * Every unsigned byte code outside 0..63 maps to NNN, including -1 encoded as 255. */
	private static int[] makePackedCodonBases(){
		final int[] packed=new int[256];
		Arrays.fill(packed, ('N'<<16)|('N'<<8)|'N');
		for(int code=0; code<64; code++){
			final String codon=codonToString[code];
			assert(codon!=null && codon.length()==3) : "Packed backtranslation requires the initialized three-base codon spelling: code="+code;
			packed[code]=(codon.charAt(0)<<16)|(codon.charAt(1)<<8)|codon.charAt(2);
		}
		return packed;
	}

	/** Immutable private snapshot of codon bases, initialized after the shared alphabet tables. */
	private static final int[] PACKED_CODON_BASES=makePackedCodonBases();

}
