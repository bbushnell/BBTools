package dna;

/**
 * Simple DNA motif implementation supporting exact and extended nucleotide matching.
 * Handles standard ACGT bases and IUPAC extended nucleotide codes for pattern searching.
 * Exact matching compares letter case variants. Extended matching requires every
 * base represented by the source symbol to be allowed by the pattern symbol;
 * ambiguity is therefore directional, not merely an intersection of possible bases.
 * Matching borrows source bytes and reads public pattern arrays; callers must keep
 * those arrays and the shared alphabet mappings unchanged during concurrent calls.
 * @author Brian Bushnell
 */
public class MotifSimple extends Motif{

	/**
	 * Test program demonstrating motif matching functionality.
	 * Creates a motif from first argument and counts matches in second argument.
	 * @param args Command-line arguments: [motif_pattern] [sequence_to_search]
	 */
	public static void main(String args[]){

		String s1="ATN";
		String s2="ATGCCCATCTGATG";

		if(args.length>0){s1=args[0];}
		if(args.length>1){s2=args[1];}

		MotifSimple m=new MotifSimple(s1, 0);
		String source=s2;

		int x=m.countExtended(source);
		System.out.println(x+" matches.");
	}

	/**
	 * Constructs a motif from a DNA sequence string with specified center position.
	 * Converts sequence to uppercase and lowercase byte arrays, determines if extended
	 * nucleotides are present, and pre-computes numeric representations for matching.
	 *
	 * @param s ASCII pattern using supported alphabet symbols; case conversion uses
	 * the default locale and byte conversion uses the default charset
	 * @param cen Center offset from pattern start, stored without range validation
	 */
	public MotifSimple(String s, int cen){
		super(s, s.length(), cen);

		commonLetters=s;
		lettersUpper=commonLetters.toUpperCase().getBytes();
		lettersLower=commonLetters.toLowerCase().getBytes();

		boolean x=false;
		for(int i=0; i<lettersUpper.length; i++){
			if(lettersUpper[i]!='A' && lettersUpper[i]!='C' && lettersUpper[i]!='G' && lettersUpper[i]!='T'){
				x=true;
			}
		}
		extended=x;

		numbers=new byte[s.length()];
		numbersExtended=new byte[s.length()];

		for(int i=0; i<lettersUpper.length; i++){
			byte b=lettersUpper[i];
			numbers[i]=baseToNumber[b];
			numbersExtended[i]=baseToNumberExtended[b];
		}
	}

	/** Tests literal case-insensitive matching at center a; out-of-array placements fail.
	 * Requires an A/C/G/T-only pattern, checked by assertion. Does not mutate source. */
	@Override
	public boolean matchesExactly(byte[] source, int a){
		assert(!extended);

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

	/** Tests subset matching of each source symbol's base mask against the pattern mask.
	 * For example, pattern N accepts source A, but pattern A rejects source N.
	 * Requires supported ASCII symbols; out-of-array placements return false. */
	@Override
	public boolean matchesExtended(byte[] source, int a){

		a=a-center;
		if(a<0 || a+length>source.length){return false;}

		for(int i=0; i<lettersUpper.length; i++){
			int x=i+a;

			byte s=source[x];
			byte n=baseToNumberExtended[s];

			if((n&numbersExtended[i])!=n){
				return false;
			}
		}
		return true;
	}

	/** Returns motif length in positions, not the alphabet size. */
	@Override
	public int numBases(){
		return numbers.length;
	}

	/** Motif sequence as uppercase ASCII bytes for case-insensitive matching */
	public final byte[] lettersUpper;
	/** Motif sequence as lowercase ASCII bytes for case-insensitive matching */
	public final byte[] lettersLower;
	/** Numeric representation of motif bases for standard ACGT matching */
	public final byte[] numbers;
	/** Extended numeric representation supporting IUPAC nucleotide codes */
	public final byte[] numbersExtended;

	/**
	 * Flag indicating whether motif contains extended nucleotide codes beyond ACGT
	 */
	public final boolean extended;

}
