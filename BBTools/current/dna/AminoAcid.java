package dna;
import java.util.Arrays;
import java.util.HashMap;

import align2.QualityTools;
import shared.KillSwitch;
import shared.Tools;
import structures.ByteBuilder;


/**
 * @author Brian Bushnell
 * @date July 1, 2010
 *
 */
public final class AminoAcid {
	
	
	public static void main(String[] args){
//		for(String s : stringToAA.keySet()){
//			System.out.println(s+"\t->\t"+stringToAA.get(s));
//		}
		
		String bases="atctgatTGGcgcgatatatcg";
		String acids=stringToAAs(bases);
		
		System.out.println(bases+" -> "+acids);
		
	}
	
	
	private AminoAcid(){
		this(null);
		assert(false);
		System.exit(0);
	}
	
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
	
	private AminoAcid(String n, String c3, String c1, String[] bases){
		name=n;
		symbol=c3;
		letter=c1.charAt(0);
		codeStrings=bases;
	}
	
	@Override
	public String toString(){
		return name+", "+symbol+", "+letter+", "+Arrays.toString(codeStrings);
	}
	
	public static String kmerToString(long kmer, int k){
		ByteBuilder sb=new ByteBuilder(k);
		for(int i=0; i<k; i++){
			int x=(int)(kmer&3);
			sb.append((char)numberToBase[x]);
			kmer>>=2;
		}
		return sb.reverse().toString();
	}
	
	public static long stringToKmer(String s){
		long kmer=0;
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			kmer=(kmer<<2)|(baseToNumber[c]);
		}
		return kmer;
	}
	
	public static String kmerToStringAA(long kmer, int k){
		ByteBuilder sb=new ByteBuilder(k);
		for(int i=0; i<k; i++){
			int x=(int)(kmer&31);
			sb.append((char)numberToAcid[x]);
			kmer>>=5;
		}
		return sb.reverse().toString();
	}
	
	public static final String codonToString(int codon){
		return codon>=0 && codon<codonToString.length ? codonToString[codon] : "NNN";
	}
	
	public String canonicalCodon(){
		return codeStrings==null || codeStrings.length<1 ? null : codeStrings[0];
	}
	
	
	public final String name;
	public final String symbol;
	public final char letter;
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
	
	public static final String[] canonicalCodons=new String[21];

	public static final byte[] numberToBase={
		'A','C','G','T','N'
	};
	
	public static final byte[] numberToAcid=new byte[21];
	
	public static final byte[] numberToComplementaryBase={
		'T','G','C','A','N'
	};
	
	public static final byte[] numberToComplement={
		3,2,1,0,4
	};
	
	public static final byte[] numberToBaseExtended={
		' ','A','C','M','G','R','S','V', //0-7
		'T','W','Y','H','K','D','B','N', //8-15
		'X',' ',' ',' ',' ',' ',' ',' ', //16-23
	};
	
	/** Has 'N' in position 0.  Mainly for translating compressed arrays containing zeroes to bases. */
	public static final byte[] numberToBaseExtended2={
		'N','A','C','M','G','R','S','V', //0-7
		'T','W','Y','H','K','D','B','N', //8-15
		'X',' ',' ',' ',' ',' ',' ',' ', //16-23
	};
	
	public static final byte[] degenerateBases={
		' ',' ',' ','M',' ','R','S','V', //0-7
		' ','W','Y','H','K','D','B',' ', //8-15
		' ',' ',' ',' ',' ',' ',' ',' ', //16-23
	};
	
	public static final byte[] numberToComplementaryBaseExtended={
		' ','T','G','K','C','Y','W','B', //0-7
		'A','S','R','D','M','H','V','N', //8-15
		'X',' ',' ',' ',' ',' ',' ',' ', //16-23
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
	
	/** Element i is: 5-bit alphabetical code for a symbol, -1 otherwise */
	public static final byte[] acidToNumber=new byte[128];
	
	/** Element i is: 5-bit alphabetical code for a symbol other than stop, -1 otherwise */
	public static final byte[] acidToNumberNoStops=new byte[128];
	
	/** Element i is: 5-bit alphabetical code for a symbol, 0 otherwise */
	public static final byte[] acidToNumber0=new byte[128];//Rename acidToNumber0
	
	/** Element i is: 5-bit alphabetical code for a symbol (plus X, B, J, Z, . and -), -1 otherwise */
	public static final byte[] acidToNumberExtended=new byte[128];
	
	/** Element i is: 5-bit alphabetical code for a symbol, -1 otherwise */
	public static final byte[] acidToNumber8=new byte[128];
	
	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', -1 otherwise */
	public static final byte[] baseToNumber=new byte[128];

	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', 0 otherwise */
	public static final byte[] baseToNumber0=new byte[128];

	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', 3 otherwise */
	public static final byte[] baseToNumber3=new byte[128];

	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', 4 otherwise */
	public static final byte[] baseToNumber4=new byte[128];
	
	/** Element i is: 3 for 'A', 2 for 'C', 1 for 'G', 0 for 'T', -1 otherwise */
	public static final byte[] baseToComplementNumber=new byte[128];
	
	/** Element i is: 3 for 'A', 2 for 'C', 1 for 'G', 0 for 'T', 0 otherwise */
	public static final byte[] baseToComplementNumber0=new byte[128];
	
	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', 4 for 'N', -1 otherwise */
	public static final byte[] baseToNumberACGTN=new byte[128];
	
	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', 0 for 'N', -1 otherwise */
	public static final byte[] baseToNumberACGTN2=new byte[128];
	
	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', 4 otherwise */
	public static final byte[] baseToNumberACGTother=new byte[128];

	/** Element i is: Some 7-bit number for ACGT, 0 otherwise */
	public static final byte[] baseToHashcode=new byte[128];
	
	/** Element i is: 0 for empty, 1-26 for A-Z, plus "*-"; otherwise, 31.
	 * Allows fast alphabetical sorting of sequences. */
	public static final byte[] symbolTo5Bit=new byte[128];
	
	/** A>A, C>C, G>G, T/U>T, other>N */
	public static final byte[] baseToACGTN=new byte[128];
	
	public static final byte[] baseToComplementExtended=new byte[128];

	public static final String[] codonToString=new String[64];
	
	/** Uracil to Thymine, everything else unchanged */
	public static final byte[] uToT=new byte[256];
	/** Thymine to Uracil, everything else unchanged */
	public static final byte[] tToU=new byte[256];
	/** . - X to N, everything else unchanged */
	public static final byte[] dotDashXToNocall=new byte[256];
	/** . - X to ., everything else unchanged */
	public static final byte[] dotDashXToNocallAA=new byte[256];
	/** Letters to uppercase, everything else unchanged */
	public static final byte[] toUpperCase=new byte[256];
	/** Lowercase to N, everything else unchanged */
	public static final byte[] lowerCaseToNocall=new byte[256];
	/** Lowercase to ., everything else unchanged */
	public static final byte[] lowerCaseToNocallAA=new byte[256];
	/** Non-acgtACGT alphabet letters to N */
	public static final byte[] iupacToNocall=new byte[256];
	
	/** Element i is the bitwise OR of constituent IUPAC base numbers in baseToNumber.<br>
	 * For example, baseToNumberExtended['M'] = ((1<<baseToNumber['A']) | (1<<baseToNumber['C'])) = (1 | 2) = 3 <br>
	 * Invalid characters are -1 */
	public static final byte[] baseToNumberExtended=new byte[128];
	public static final byte[] baseToComplementNumberExtended=new byte[128];
	public static final AminoAcid[] AlphabeticalAAs=new AminoAcid[21];
	public static final AminoAcid[] codeToAA=new AminoAcid[66];
	public static final char[] codeToChar=new char[66];
	public static final byte[] codeToByte=new byte[66];
	public static final byte[] aminoToCode=new byte[128];
	public static final HashMap<String, AminoAcid> stringToAA=new HashMap<String, AminoAcid>(512);
	
	public static final AminoAcid Alanine=new AminoAcid("Alanine, Ala, A, GCU, GCC, GCA, GCG");
	public static final AminoAcid Arginine=new AminoAcid("Arginine, Arg, R, CGU, CGC, CGA, CGG, AGA, AGG");
	public static final AminoAcid Asparagine=new AminoAcid("Asparagine, Asn, N, AAU, AAC");
	public static final AminoAcid AsparticAcid=new AminoAcid("AsparticAcid, Asp, D, GAU, GAC");
	public static final AminoAcid Cysteine=new AminoAcid("Cysteine, Cys, C, UGU, UGC");
	public static final AminoAcid GlutamicAcid=new AminoAcid("GlutamicAcid, Glu, E, GAA, GAG");
	public static final AminoAcid Glutamine=new AminoAcid("Glutamine, Gln, Q, CAA, CAG");
	public static final AminoAcid Glycine=new AminoAcid("Glycine, Gly, G, GGU, GGC, GGA, GGG");
	public static final AminoAcid Histidine=new AminoAcid("Histidine, His, H, CAU, CAC");
	public static final AminoAcid Isoleucine=new AminoAcid("Isoleucine, Ile, I, AUU, AUC, AUA");
	public static final AminoAcid Leucine=new AminoAcid("Leucine, Leu, L, UUA, UUG, CUU, CUC, CUA, CUG");
	public static final AminoAcid Lysine=new AminoAcid("Lysine, Lys, K, AAA, AAG");
	public static final AminoAcid Methionine=new AminoAcid("Methionine, Met, M, AUG");
	public static final AminoAcid Phenylalanine=new AminoAcid("Phenylalanine, Phe, F, UUU, UUC");
	public static final AminoAcid Proline=new AminoAcid("Proline, Pro, P, CCU, CCC, CCA, CCG");
	public static final AminoAcid Serine=new AminoAcid("Serine, Ser, S, UCU, UCC, UCA, UCG, AGU, AGC");
	public static final AminoAcid Threonine=new AminoAcid("Threonine, Thr, T, ACU, ACC, ACA, ACG");
	public static final AminoAcid Tryptophan=new AminoAcid("Tryptophan, Trp, W, UGG");
	public static final AminoAcid Tyrosine=new AminoAcid("Tyrosine, Tyr, Y, UAU, UAC");
	public static final AminoAcid Valine=new AminoAcid("Valine, Val, V, GUU, GUC, GUA, GUG");
	
	public static final AminoAcid Selenocysteine=new AminoAcid("Selenocysteine, Sec, U"); //UGA sometimes
	public static final AminoAcid Pyrrolysine=new AminoAcid("Pyrrolysine, Pyl, O");
	
	public static final AminoAcid END=new AminoAcid("End, End, *, UAA, UGA, UAG");
	public static final AminoAcid ANY=new AminoAcid("Any, Any, X, XXX");
	
	public static int AMINO_SHIFT=5;


	public static final byte[][] COLORS=new byte[][] {
		{0, 1, 2, 3},
		{1, 0, 3, 2},
		{2, 3, 0, 1},
		{3, 2, 1, 0}
	};

	public static boolean equalsReverseComp(byte[] bases, byte[] bases2) {
		if(bases.length!=bases2.length) {return false;}
		for(int i=0, j=bases2.length-1; i<bases.length; i++, j--) {
			byte a=Tools.toUpperCase(bases[i]), b=baseToComplementExtended[bases2[j]];
			if(a!=b) {return false;}
		}
		return true;
	}
	
	/** Returns a new reverse-complemented array in ASCII coding*/
	public static final byte[] reverseComplementBases(final byte[] in){
		byte[] out=new byte[in.length];
		final int last=in.length-1;
		for(int i=0; i<in.length; i++){
			out[i]=baseToComplementExtended[in[last-i]];
		}
		return out;
	}
	

	public static final void reverseComplementBasesInPlace(final byte[] in){
		if(in!=null){reverseComplementBasesInPlace(in, in.length);}
	}
	public static final void complementBasesInPlace(final byte[] in){
		if(in==null){return;}
		complementBasesInPlace(in, in.length);
	}
	public static final void complementBasesInPlace(final byte[] in, final int length){
		if(in==null){return;}
		for(int i=0; i<length; i++){
			in[i]=baseToComplementExtended[in[i]];
		}
	}
	
	public static final void reverseComplementBasesInPlace(final byte[] in, final int length){
		if(in==null){return;}
		final int last=length-1;
		final int max=length/2;
		for(int i=0; i<max; i++){
			byte a=in[i];
			byte b=in[last-i];
//			assert(b>0 && b<baseToComplementExtended.length) : ((int)b)+"\t"+((char)b)+"\t"+Arrays.toString(in);
//			System.out.println((char)a+", "+(char)b+", "+i+", "+last);
			in[i]=baseToComplementExtended[b];
			in[last-i]=baseToComplementExtended[a];
		}
		if((length&1)==1){//Odd length; process middle
			in[max]=baseToComplementExtended[in[max]];
		}
	}
	
	public static final String reverseComplementBases(String in){
		return in==null ? null : new String(reverseComplementBases(in.getBytes()));
	}
	
	public static final int reverseComplementBinary(int kmer, int k){
		int out=0;
		kmer=~kmer;
		for(int i=0; i<k; i++){
			out=((out<<2)|(kmer&3));
			kmer>>=2;
		}
		return out;
	}
	
	public static final long reverseComplementBinary(long kmer, int k){
		long out=0;
		kmer=~kmer;
		for(int i=0; i<k; i++){
			out=((out<<2)|(kmer&3L));
			kmer>>=2;
		}
		return out;
	}
	
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
	
	public static final long reverseComplementBinaryFast(long kmer, int k){
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
	
	public static final byte baseToColor(byte base1, byte base2){
		byte a=baseToNumber[base1];
		byte b=baseToNumber[base2];
		if(a<0 && b<0){return 'N';}
		if(a<0){a=3;}
		if(b<0){b=3;}
		return COLORS[a][b];
	}
	
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
	
	public static final AminoAcid toAA(String code){
		return toAA(code.charAt(0), code.charAt(1), code.charAt(2));
	}
	
	public static final char toChar(String code){
		return toChar(code.charAt(0), code.charAt(1), code.charAt(2));
	}
	
	public static final char[] splitBase(char c){
		byte b=baseToNumberExtended[c];
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
	

	
	
	public static final byte[] numberToBases(int code, int n){
		
		byte[] bytes=KillSwitch.allocByte1D(n);
		
		for(int i=n-1; i>=0; i--){
			int temp=code&3;
			code>>=2;
			bytes[i]=numberToBase[temp];
		}
		
		return bytes;
	}
	
	public static final int baseTupleToNumber(byte[] tuple){
		
		int r=0;
		for(int i=0; i<tuple.length; i++){
			int temp=baseToNumberACGTN[tuple[i]];
			if(temp<0 || temp>3){return -1;}
			r=((r<<2)|temp);
		}
		
		return r;
	}
	
	public static boolean isFullyDefined(char base){
		return baseToNumber[base]>=0;
	}
	
	public static boolean isFullyDefined(byte base){
		return base>=0 && baseToNumber[base]>=0;
	}
	
	public static boolean isFullyDefinedAA(byte acid){
		return acid>=0 && acidToNumber[acid]>=0;
	}
	
	public static boolean isFullyDefinedAANoStops(byte acid){
		return acid>=0 && acidToNumberNoStops[acid]>=0;
	}
	
	public static boolean isACGTN(String bases){
		for(int i=0; i<bases.length(); i++) {
			if(!isACGTN(bases.charAt(i))) {return false;}
		}
		return true;
	}
	
	public static boolean isACGTN(byte[] bases){
		for(int i=0; i<bases.length; i++) {
			if(!isACGTN(bases[i])) {return false;}
		}
		return true;
	}
	
	public static boolean isACGTN(char base){
		return baseToNumberACGTN[base]>=0;
	}
	
	public static boolean isACGTN(byte base){
		return base>=0 && baseToNumberACGTN[base]>=0;
	}
	
	public static boolean containsOnlyACGTN(String s){
		if(s==null || s.length()==0){return true;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(baseToNumberACGTN[c]<0){return false;}
		}
		return true;
	}
	
	public static boolean containsOnlyACGTNQ(String s){
		if(s==null || s.length()==0){return true;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c!='?' && baseToNumberACGTN[c]<0){return false;}
		}
		return true;
	}
	
	public static boolean containsOnlyACGTN(byte[] array){
		if(array==null || array.length==0){return true;}
		for(int i=0; i<array.length; i++){
			byte b=array[i];
			if(b<0 || baseToNumberACGTN[b]<0){return false;}
		}
		return true;
	}
	
	public static boolean isFullyDefined(String s){
		for(int i=0; i<s.length(); i++){
			if(!isFullyDefined(s.charAt(i))){return false;}
		}
		return true;
	}
	
	public static boolean isFullyDefined(byte[] s){
		for(int i=0; i<s.length; i++){
			if(!isFullyDefined(s[i])){return false;}
		}
		return true;
	}
	
	public static int countUndefined(byte[] s){
		int x=0;
		for(int i=0; i<s.length; i++){
			if(!isFullyDefined(s[i])){x++;}
		}
		return x;
	}
	
	public static int countDefined(byte[] s){
		if(s==null){return 0;}
		int x=0;
		for(int i=0; i<s.length; i++){
			if(isFullyDefined(s[i])){x++;}
		}
		return x;
	}

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
	
	public static final byte toNumber(char c1, char c2, char c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		int x=(baseToNumberACGTN2[c1]<<4)|(baseToNumberACGTN2[c2]<<2)|(baseToNumberACGTN2[c3]);
		return (byte)x;
	}
	
	public static final AminoAcid toAA(char c1, char c2, char c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		int x=(baseToNumberACGTN2[c1]<<4)|(baseToNumberACGTN2[c2]<<2)|(baseToNumberACGTN2[c3]);
		return codeToAA[x];
	}
	
	public static final char toChar(char c1, char c2, char c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		int x=(baseToNumberACGTN2[c1]<<4)|(baseToNumberACGTN2[c2]<<2)|(baseToNumberACGTN2[c3]);
		return codeToChar[x];
	}
	
	public static final byte toByte(byte c1, byte c2, byte c3){
		int a=baseToNumber[c1], b=baseToNumber[c2], c=baseToNumber[c3];
		if(a<0 || b<0 || c<0){return (byte)'X';}
		int x=((a<<4)|(b<<2)|c);
		return codeToByte[x];
	}
	
	public static final char toChar(byte c1, byte c2, byte c3){
		assert(baseToNumberACGTN2[c1]>=0 && baseToNumberACGTN2[c2]>=0 && baseToNumberACGTN2[c3]>=0);
		byte n1=baseToNumberACGTN2[c1], n2=baseToNumberACGTN2[c2], n3=baseToNumberACGTN2[c3];
		if(n1>3 || n2>3 || n3>3){return '?';}
		int x=(n1<<4)|(n2<<2)|(n3);
//		return (x<codeToChar.length ? codeToChar[x] : '?');
		return codeToChar[x];
	}
	
	public static final String stringToAAs(String bases){
		StringBuilder sb=new StringBuilder(bases.length()/3);
		for(int i=2; i<bases.length(); i+=3){
			char a=toAA(bases.charAt(i-2), bases.charAt(i-1), bases.charAt(i)).letter;
			sb.append(a);
		}
		return sb.toString();
	}
	
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
	
	public static final byte[][] toQualitySixFrames(byte[] quals, int offset){
		byte[][] out=new byte[6][];
		if(quals!=null && quals.length>2){
			for(int i=0; i<3; i++){
				out[i]=toAAQuality(quals, i);
			}
			Tools.reverseInPlace(quals);
			for(int i=0; i<3; i++){
				out[i+3]=toAAQuality(quals, i);
			}
			Tools.reverseInPlace(quals);
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
	
	public static final byte[] toAAs(byte[] bases, int start, int stop){
		if(bases==null){return null;}
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
	
	public static final byte[] toNTs(final byte[] aminos){
		if(aminos==null){return null;}
		final int alen=aminos.length;
		final int blen=alen*3;
		
		final byte[] out=KillSwitch.allocByte1D(blen);
		for(int i=0, j=0; i<alen; i++, j+=3){
			int code=aminoToCode[aminos[i]];
			out[j+2]=numberToBase[(code&3)];
			out[j+1]=numberToBase[((code>>2)&3)];
			out[j]=numberToBase[((code>>4)&3)];
		}
		return out;
	}
	
	public static final boolean isHomopolymer(int kmer, int k, int maxRepeat) {
		if(maxRepeat<1) {return false;}
		final int mask=~((-1)<<(2*k));
		final int inv=(~kmer)&mask;
		boolean polymer=false;
		for(int i=1; i<=maxRepeat && !polymer; i++) {
			int shift=2*i;
			polymer=(kmer==(kmer|(kmer>>shift))) && (inv==(inv|(inv>>shift)));
		}
		return polymer;
	}
	
	public static final short[] rcompBinaryTable=makeBinaryRcompTable(4);
	
	private static final short[] makeBinaryRcompTable(int k){
		int bits=2*k;
		short[] r=new short[1<<bits];
		for(int i=0; i<r.length; i++){
			r[i]=(short)reverseComplementBinary(i, k);
		}
		return r;
	}
	
	static {
		
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
			toUpperCase[i]=(byte) ((i>='a' && i<='z') ? i-32 : i);
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
	
}
