package test.dna;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.AminoAcid;
import stream.Read;

/** Literal alphabet/translation controls and focused candidate-defect oracles.
 * Run each case in a fresh JVM with assertions enabled or disabled.
 * @author Shinobu
 */
public final class AminoAcidContractTest{

	/** Runs one case; failures record the discrepancy before any production repair. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){throw new IllegalArgumentException("Require case=controls|reverse-case|backtranslate|backtranslate-read|backtranslate-all|split-controls|split-invalid");}
		final String mode=args[0].substring(5);
		if(mode.equals("controls")){controls();}
		else if(mode.equals("reverse-case")){reverseCase();}
		else if(mode.equals("backtranslate")){backtranslate();}
		else if(mode.equals("backtranslate-read")){backtranslateRead();}
		else if(mode.equals("backtranslate-all")){backtranslateAll();}
		else if(mode.equals("split-controls")){splitControls();}
		else if(mode.equals("split-invalid")){splitInvalid();}
		else{throw new IllegalArgumentException(mode);}
		System.out.println("PASS "+mode);
	}

	/** Preserves ordinary codons, strand conversion and binary encodings. */
	private static void controls(){
		equal(AminoAcid.toNTs(bytes("ACDMW*")), "GCTTGTGATATGTGGTAA", "Canonical codons");
		equal(AminoAcid.toNTs(bytes("acdmw*")), "GCTTGTGATATGTGGTAA", "Lowercase canonical codons");
		equal(AminoAcid.toAAs(bytes("ATGTGGGCTTAA"), 0), "MWA*", "Byte translation");
		equal(AminoAcid.toAAs(bytes("ATGTGGGCTTAA"), 0, 11), "MWA", "ORF range omits terminal stop");
		require(AminoAcid.toByte((byte)'A', (byte)'N', (byte)'G')=='X', "Unknown byte codon must remain unknown");
		equal(AminoAcid.reverseComplementBases(bytes("AaCGTryN")), "NryACGtT", "IUPAC complement and case");
		long state=7;
		for(int k=1; k<=32; k++){
			for(int i=0; i<1000; i++){
				state=state*6364136223846793005L+1442695040888963407L;
				final String text=AminoAcid.kmerToString(state, k);
				final String reverse=AminoAcid.reverseComplementBases(text);
				final long expected=AminoAcid.stringToKmer(reverse);
				require(AminoAcid.reverseComplementBinary(state, k)==expected, "Loop reverse complement mismatch at k="+k);
				require(AminoAcid.reverseComplementBinaryFast(state, k)==expected, "Fast reverse complement mismatch at k="+k);
				require(AminoAcid.reverseComplementBinaryFastLookup(state, k)==expected, "Lookup reverse complement mismatch at k="+k);
			}
		}
	}

	/** Equivalent DNA letters should compare regardless of which operand is lowercase. */
	private static void reverseCase(){
		boolean allMatched=true;
		for(String a : new String[]{"ACG", "acg", "AcG"}){
			for(String b : new String[]{"CGT", "cgt", "CgT"}){
				final boolean result=AminoAcid.equalsReverseComp(bytes(a), bytes(b));
				System.out.println(a+" versus "+b+" => "+result);
				allMatched&=result;
			}
		}
		require(!AminoAcid.equalsReverseComp(bytes("ACG"), bytes("CAT")), "Mismatch control was accepted");
		require(allMatched, "Letter case changed reverse-complement equality");
	}

	/** Unknown amino acids must not become a concrete asparagine codon. */
	private static void backtranslate(){
		equal(AminoAcid.toNTs(bytes("X")), "NNN", "Unknown amino acid");
		equal(AminoAcid.toNTs(bytes("XBZJOUx")), "NNNNNNNNNNNNNNNNNNNNN", "Unresolved canonical codons");
	}

	/** Checks every ASCII input against literal representative codons, independent of lookup construction. */
	private static void backtranslateAll(){
		final String residues="ACDEFGHIKLMNPQRSTVWY*";
		final String[] codons={"GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "TTA",
				"ATG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT", "TAA"};
		for(int symbol=0; symbol<128; symbol++){
			final int index=residues.indexOf(Character.toUpperCase((char)symbol));
			final String expected=index<0 ? "NNN" : codons[index];
			final byte[] actual=AminoAcid.toNTs(new byte[]{(byte)symbol});
			require(Arrays.equals(actual, bytes(expected)), "Backtranslation mismatch at ASCII "+symbol+
					": expected "+expected+", actual "+new String(actual, StandardCharsets.US_ASCII));
		}
		require(AminoAcid.toNTs(null)==null && AminoAcid.toNTs(new byte[0]).length==0,
				"Null/empty backtranslation controls changed");
		final byte oldT=AminoAcid.numberToBase[3];
		try{
			AminoAcid.numberToBase[3]='U';
			equal(AminoAcid.toNTs(bytes("M*X")), "AUGUAANNN", "Substituted alphabet fallback");
		}finally{AminoAcid.numberToBase[3]=oldT;}
		final String[] allCodons=new String[64];
		int next=0;
		for(char first : new char[]{'A', 'C', 'G', 'T'}){
			for(char second : new char[]{'A', 'C', 'G', 'T'}){
				for(char third : new char[]{'A', 'C', 'G', 'T'}){
					allCodons[next++]=new String(new char[]{first, second, third});
				}
			}
		}
		final byte oldCode=AminoAcid.aminoToCode['?'];
		try{
			for(int code=-128; code<128; code++){
				AminoAcid.aminoToCode['?']=(byte)code;
				final String expected=code>=0 && code<64 ? allCodons[code] : "NNN";
				require(Arrays.equals(AminoAcid.toNTs(bytes("?")), bytes(expected)),
						"Retrieved amino code "+code+" must decode as "+expected);
			}
		}finally{AminoAcid.aminoToCode['?']=oldCode;}
		System.out.println("PASS 128 ASCII inputs, 256 code substitutions and null/empty controls");
	}

	/** Exercises the live Read caller independently of the direct-helper failure. */
	private static void backtranslateRead(){
		final byte[] original=bytes("MX*"), qualities={10, 20, 30};
		final Read protein=new Read(original, qualities, "protein", 1, false);
		protein.setAminoAcid(true);
		final Read converted=protein.aminoToNucleic();
		equal(converted.bases, "ATGNNNTAA", "Read aminoToNucleic caller");
		require(Arrays.equals(protein.bases, bytes("MX*")) && protein.aminoacid(), "Caller mutated original protein");
		require(!converted.aminoacid() && converted.quality.length==9, "Caller metadata or qualities did not expand");
		require(Arrays.equals(converted.quality, new byte[]{15, 15, 15, 25, 25, 25, 35, 35, 35}),
				"Caller must expand each quality to three copies with the documented five-point increase");
		require(Arrays.equals(protein.quality, new byte[]{10, 20, 30}), "Caller mutated original qualities");
	}

	/** Preserves literal IUPAC expansions, case, U and supported public-table substitutions. */
	private static void splitControls(){
		final String letters="ACGTURYSWKMBDHVN";
		final String[] expected={"A", "C", "G", "T", "T", "AG", "CT", "CG", "AT", "GT", "AC",
				"CGT", "AGT", "ACT", "ACG", "ACGT"};
		for(int i=0; i<letters.length(); i++){
			final char letter=letters.charAt(i);
			require(new String(AminoAcid.splitBase(letter)).equals(expected[i]), "Uppercase expansion: "+letter);
			require(new String(AminoAcid.splitBase(Character.toLowerCase(letter))).equals(expected[i]),
					"Lowercase expansion: "+letter);
		}
		final byte oldT=AminoAcid.numberToBase[3], oldQuestion=AminoAcid.baseToNumberExtended['?'];
		try{
			AminoAcid.numberToBase[3]='U';
			AminoAcid.baseToNumberExtended['?']=5;
			require(new String(AminoAcid.splitBase('N')).equals("ACGU"), "Splitting must use the current output alphabet");
			require(new String(AminoAcid.splitBase('?')).equals("AG"), "A substituted valid mask must remain supported");
		}finally{
			AminoAcid.numberToBase[3]=oldT;
			AminoAcid.baseToNumberExtended['?']=oldQuestion;
		}
		System.out.println("PASS 32 literal symbol expansions and two public-table controls");
	}

	/** Unsupported masks and out-of-table characters must reject instead of emitting partial/NUL arrays. */
	private static void splitInvalid(){
		boolean ok=true;
		for(char letter : new char[]{'X', 'x', '?', '-', (char)0, (char)127, (char)128, (char)65535}){
			RuntimeException failure=null;
			int length=-1;
			try{length=AminoAcid.splitBase(letter).length;}
			catch(RuntimeException e){failure=e;}
			final boolean rejected=failure instanceof IllegalArgumentException;
			System.out.println("code="+(int)letter+" length="+length+" exception="+
					(failure==null ? "none" : failure.getClass().getSimpleName())+" rejected="+rejected);
			ok&=rejected;
		}
		require(ok, "splitBase must reject masks outside 1..15 and unsupported character indices with a clear argument error");
	}

	/** Compares against an explicit ASCII oracle, printing the observed bytes. */
	private static void equal(byte[] actual, String expected, String what){
		System.out.println(what+": "+(actual==null ? "null" : new String(actual, StandardCharsets.US_ASCII)));
		require(Arrays.equals(actual, bytes(expected)), what+" expected "+expected);
	}
	private static byte[] bytes(String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
