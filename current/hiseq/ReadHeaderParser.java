package hiseq;

import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/**
 * Abstract superclass for parsing Illumina sequencer read headers.
 * Provides common functionality for extracting metadata from various
 * Illumina platforms including MiSeq, HiSeq, NextSeq, NovaSeq, and NovaSeqX.
 * Concrete variants define supported layouts and optional parsing modes.
 * Instances retain mutable parse state; callers must confine each instance.
 *
 * @author Brian Bushnell
 */
public abstract class ReadHeaderParser{

	/*--------------------------------------------------------------*/
	/*----------------        Expected Format       ----------------*/
	/*--------------------------------------------------------------*/

	//@VP2-06:112:H7LNDMCVY:2:2437:14181:20134 (Novaseq6k)
	//2402:6:1101:6337:2237/1
	//MISEQ08:172:000000000-ABYD0:1:1101:18147:1925 1:N:0:TGGATATGCGCCAATT
	//HISEQ07:419:HBFNEADXX:1:1101:1238:2072
	//A00178:38:H5NYYDSXX:2:1101:3007:1000 1:N:0:CAACCTA+CTAGGTT
	//@LH00223:28:22GLGMLT3:1:1101:5928:1016 1:N:0:CTGCTTGGTT+CTAACGACAG (NovaseqX)

	//	@HWI-Mxxxx or @Mxxxx - MiSeq
	//	@HWUSI - GAIIx
	//	@HWI-Dxxxx - HiSeq 2000/2500
	//	@Kxxxx - HiSeq 3000(?)/4000
	//	@Nxxxx - NextSeq 500/550
	//
	//	AAXX = Genome Analyzer
	//	BCXX = HiSeq v1.5
	//	ACXX = HiSeq High-Output v3
	//	ANXX = HiSeq High-Output v4
	//	ADXX = HiSeq RR v1
	//	AMXX, BCXX =HiSeq RR v2
	//	ALXX = HiSeqX
	//	BGXX, AGXX = High-Output NextSeq
	//	AFXX = Mid-Output NextSeq
	//	5 letter/number = MiSeq

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public final void test(String s){
		if(s==null){s="LH00223:28:22GLGMLT3:1:1101:5928:1016 1:N:0:CTGCTTGGTT+CTAACGACAG";}
		parse(s);
		System.err.println("ihp="+this);
		System.err.println("id="+id());


		System.err.println("whitespaceIndex="+whitespaceIndex());


		System.err.println("machine="+machine());
		System.err.println("run="+run());
		System.err.println("flowcell="+flowcell());
		System.err.println("lane="+lane());
		System.err.println("tile="+tile());
		System.err.println("xPos="+xPos());
		System.err.println("yPos="+yPos());
		System.err.println("surface="+surface());
		System.err.println("swath="+swath());
		System.err.println("pairCode="+pairCode());
		System.err.println("pairnum="+pairnum());
		System.err.println("chastityCode="+chastityCode());
		System.err.println("chastityFail="+chastityFail());
		System.err.println("controlBits="+controlBits());
		System.err.println("barcode="+barcode());
		System.err.println("extra="+extra());
		System.err.println("index3="+index3());
		System.err.println("commentSeparator='"+commentSeparator()+"' ("+(int)commentSeparator()+")");
		System.err.println("pairnum="+pairnum());
		System.err.println("barcodeDelimiter='"+barcodeDelimiter()+"' ("+(int)barcodeDelimiter()+")");
		System.err.println("barcodeLength1="+barcodeLength1());
		System.err.println("barcodeLength2="+barcodeLength2());
		System.err.println("barcodeLength="+barcodeLength());
		System.err.println("barcodeLetters="+barcodeLetters());
		System.err.println("numBarcodes="+numBarcodes());
	}

	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces this instance's input/state and returns the parser. */
	public abstract ReadHeaderParser parse(String id);
	/** Instrument identifier, or null when the variant does not extract it. */
	public abstract String machine();
	/** Instrument run number, or the variant's unavailable sentinel. */
	public abstract int run();
	/** Flowcell identifier, or null when unavailable. */
	public abstract String flowcell();
	/** Numeric lane, according to the concrete parser's supported layout. */
	public abstract int lane();
	/** Numeric tile, according to the concrete parser's supported layout. */
	public abstract int tile();
	/** X coordinate in the header's pixel-coordinate domain. */
	public abstract int xPos();
	/** Y coordinate in the header's pixel-coordinate domain. */
	public abstract int yPos();
	/** Mate marker, normally '1' or '2'; availability depends on the concrete layout. */
	public abstract char pairCode();
	/** Filter code when available: 'Y' denotes failure and 'N' denotes pass. */
	public abstract char chastityCode();
	/** Numeric comment control field, when available. */
	public abstract int controlBits();
	/** Barcode text, possibly dual-index, or null when unavailable. */
	public abstract String barcode();
	/** Additional parsed comment text, or null when unavailable. */
	public abstract String extra();
	/** Sample identifier, or null when the concrete variant does not extract it. */
	public abstract String sample();

	/*--------------------------------------------------------------*/
	/*----------------       Concrete Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses the borrowed read's ID without modifying the read. */
	public final ReadHeaderParser parse(Read r){
		return parse(r.id);
	}

	/** Returns the retained input identifier. */
	public final String id(){return id;}

	/** Returns the first literal space/slash separator, or zero if none exists. */
	public char commentSeparator(){
		for(int i=0; i<id.length(); i++){
			char c=id.charAt(i);
			if(c==' ' || c== '/'){return c;}
		}
		return 0;
	}

	/** Returns tile()/1000, the surface prefix of a conventional four-digit tile; does not validate its layout. */
	public int surface(){
		return surface(tile());
	}

	/** Returns the hundreds digit of tile(), conventionally the swath; does not validate its layout. */
	public int swath(){
		return swath(tile());
	}

	/** Converts the parsed '1'/'2' code to zero-based mate number. */
	public int pairnum(){
		return pairCode()-(int)'1';
	}

	/** Returns true for a 'Y' filter code; asserts the expected Y/N alphabet. */
	public boolean chastityFail(){
		int c=chastityCode();
		assert(c=='N' || c=='Y') : c;
		return c=='Y';
	}

	/** First barcode character rejected by Character.isLetter, or zero for null/empty/all-letter text. */
	public char barcodeDelimiter(){
		return barcodeDelimiter(barcode());
	}

	/** Leading ASCII-letter count before the first nonletter; zero for a null/empty barcode. */
	public int barcodeLength1(){
		return barcodeLength1(barcode());
	}

	/** Trailing ASCII-letter count after the last nonletter; zero for null/empty/all-letter text. */
	public int barcodeLength2(){
		return barcodeLength2(barcode());
	}

	/** Total barcode UTF-16 length including delimiters; zero for null or empty text. */
	public int barcodeLength(){
		String bc=barcode();
		return bc==null ? 0 : bc.length();
	}

	/** Counts ASCII letters throughout the barcode using Tools.isLetter; zero for null/empty text. */
	public int barcodeLetters(){
		return barcodeLetters(barcode());
	}

	/** Counts nonempty leading/trailing ASCII-letter runs (0..2); ignores interior runs and index3(). */
	public int numBarcodes(){
		return numBarcodes(barcode());//TODO: could add index3 count here too
	}

	/** Optional extra index between coordinates and comment; absent by default. */
	public String index3(){return null;}
	/** Field index terminated by the first whitespace boundary, or -1 when unavailable. */
	public int whitespaceIndex(){return -1;}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** First literal space or slash, or zero for empty/separator-free text; id must be nonnull (tabs are ignored). */
	public static char commentSeparator(String id){
		for(int i=0; i<id.length(); i++){
			char c=id.charAt(i);
			if(c==' ' || c== '/'){return c;}
		}
		return 0;
	}

	//Illumina tile number is SSWWTT-ish: surface = thousands digit, swath = hundreds digit.
	//e.g. tile 1101 -> surface 1, swath 1. (tile%1000 strips the surface, /100 takes the swath digit.)
	/** Returns tile/1000, conventionally the surface prefix; no tile-range or sign validation. */
	public static final int surface(int tile){
		return tile/1000;
	}

	/** Returns (tile%1000)/100, conventionally the swath digit; no tile-range or sign validation. */
	public static final int swath(int tile){
		return (tile%1000)/100;
	}

	/** First nonletter under Character.isLetter (Unicode char classification), or zero for null/empty/all-letter text.
	 * This differs from the ASCII-only Tools.isLetter used by the length/count helpers. */
	public static final char barcodeDelimiter(String bc){
		if(bc==null){return 0;}
		for(int i=0; i<bc.length(); i++){
			char c=bc.charAt(i);
			if(!Character.isLetter(c)){return c;}
		}
		return 0;
	}

	/** Leading run of ASCII letters, ending at the first character rejected by Tools.isLetter.
	 * Returns zero for null/empty text and the full length for all-ASCII-letter text. */
	public static int barcodeLength1(String bc){
		if(bc==null){return 0;}
		for(int i=0; i<bc.length(); i++){
			char c=bc.charAt(i);
			if(!Tools.isLetter(c)){return i;}
		}
		return bc.length();
	}

	//Dual-index barcode "AAAA+CCCC": length1 = letters before the first non-letter (the '+' delimiter),
	//length2 = letters after the last non-letter (scanning from the end). A single-index barcode (all
	//letters) has no delimiter, so length2 returns 0 here and numBarcodes() counts it as 1.
	/** Trailing run of ASCII letters after the last character rejected by Tools.isLetter.
	 * Returns zero for null/empty text or text with no nonletter, including a single all-letter index. */
	public static int barcodeLength2(String bc){
		if(bc==null){return 0;}
		for(int i=bc.length()-1; i>=0; i--){
			char c=bc.charAt(i);
			if(!Tools.isLetter(c)){return bc.length()-1-i;}
		}
		return 0;
	}

	/** Counts ASCII letters in all positions using Tools.isLetter; ignores nonletters, returns zero for null/empty text. */
	public static int barcodeLetters(String bc){
		if(bc==null){return 0;}
		int letters=0;
		for(int i=0; i<bc.length(); i++){
			char c=bc.charAt(i);
			letters+=(Tools.isLetter(c) ? 1 : 0);
		}
		return letters;
	}

	/** Counts nonempty leading and trailing ASCII-letter runs (0..2), not arbitrary delimited fields.
	 * Null/empty text gives zero; an all-ASCII-letter single index gives one; interior runs are ignored. */
	public static int numBarcodes(String bc){
		return (barcodeLength1(bc)>0 ? 1 : 0)+(barcodeLength2(bc)>0 ? 1 : 0);
	}

	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		bb.append("lane:\t").append(lane()).nl();
		bb.append("tile:\t").append(tile()).nl();
		bb.append("x:\t").append(xPos()).nl();
		bb.append("y:\t").append(yPos()).nl();
		bb.append("pairnum:\t").append(pairCode()).nl();
		bb.append("barcode:\t").append(barcode()).nl();
		bb.append("chastity:\t").append(chastityCode()).nl();
		return bb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*----------------        Private Fields        ----------------*/
	/*--------------------------------------------------------------*/

	String id;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Configure legacy coordinate parsing before use; indexed variants may ignore this flag. */
	public static boolean PARSE_COORDINATES=true;
	/** Configure legacy comment parsing before use; concrete variants may ignore this flag. */
	public static boolean PARSE_COMMENT=false;

}
