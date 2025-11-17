package hiseq;

import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/**
 * Superclass for Illumina header parsers.
 * @author Brian Bushnell
 * @date April 3, 2024
 *
 */
public abstract class ReadHeaderParser {
	
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
	
	/**
	 * Tests header parsing functionality with a sample read ID.
	 * Prints all extracted metadata fields to standard error for debugging.
	 * @param s Read header string to test (uses default NovaSeqX example if null)
	 */
	public final void test(String s) {
		if(s==null) {s="LH00223:28:22GLGMLT3:1:1101:5928:1016 1:N:0:CTGCTTGGTT+CTAACGACAG";}
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
	
	/**
	 * Parses a read header string to extract sequencer metadata.
	 * Implementation varies by sequencer platform and header format.
	 * @param id The read header string to parse
	 * @return This parser instance for method chaining
	 */
	public abstract ReadHeaderParser parse(String id);
	/** Returns the sequencer machine identifier.
	 * @return Machine name/ID string */
	public abstract String machine();
	/** Returns the sequencing run number.
	 * @return Run number as integer */
	public abstract int run();
	/** Returns the flowcell identifier.
	 * @return Flowcell ID string */
	public abstract String flowcell();
	/** Returns the flowcell lane number.
	 * @return Lane number (typically 1-8 depending on platform) */
	public abstract int lane();
	/**
	 * Returns the tile number within the lane.
	 * Encodes surface and swath information in thousands/hundreds digits.
	 * @return Tile number
	 */
	public abstract int tile();
	/** Returns the X coordinate position on the tile.
	 * @return X position coordinate */
	public abstract int xPos();
	/** Returns the Y coordinate position on the tile.
	 * @return Y position coordinate */
	public abstract int yPos();
	/** Returns the pair code character indicating read pairing.
	 * @return '1' for first read, '2' for second read in paired-end sequencing */
	public abstract char pairCode();
	/** Returns the chastity filter code.
	 * @return 'Y' if read failed chastity filter, 'N' if passed */
	public abstract char chastityCode();
	/** Returns the control bits field.
	 * @return Control bits as integer */
	public abstract int controlBits();
	/** Returns the barcode sequence(s) for sample demultiplexing.
	 * @return Barcode string, may contain multiple barcodes separated by delimiters */
	public abstract String barcode();
	/** Returns additional header information not covered by other fields.
	 * @return Extra header data as string */
	public abstract String extra();
	/** Returns the sample identifier.
	 * @return Sample name/ID string */
	public abstract String sample();
	
	/*--------------------------------------------------------------*/
	/*----------------       Concrete Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses a Read object's header to extract sequencer metadata.
	 * @param r The Read object containing the header to parse
	 * @return This parser instance for method chaining
	 */
	public final ReadHeaderParser parse(Read r){
		return parse(r.id);
	}
	
	/** Returns the currently parsed read header string */
	public final String id() {return id;}
	
	/** Finds the character separating the main header from comment fields.
	 * @return Space (' ') or slash ('/') separator, or 0 if none found */
	public char commentSeparator() {
		for(int i=0; i<id.length(); i++){
			char c=id.charAt(i);
			if(c==' ' || c== '/'){return c;}
		}
		return 0;
	}
	
	/** Extracts the surface number from the tile number.
	 * @return Surface number (thousands digit of tile number) */
	public int surface() {
		return surface(tile());
	}
	
	/** Extracts the swath number from the tile number.
	 * @return Swath number (hundreds digit of tile number) */
	public int swath() {
		return swath(tile());
	}
	
	/** Converts pair code character to numeric value.
	 * @return 0 for first read, 1 for second read */
	public int pairnum() {
		return pairCode()-(int)'1';
	}
	
	/** Determines if the read failed the chastity filter.
	 * @return true if chastity code is 'Y' (failed), false if 'N' (passed) */
	public boolean chastityFail() {
		int c=chastityCode();
		assert(c=='N' || c=='Y') : c;
		return c=='Y';
	}
	
	/** Finds delimiter character separating multiple barcodes.
	 * @return Delimiter character, or 0 if no delimiter found */
	public char barcodeDelimiter() {
		return barcodeDelimiter(barcode());
	}
	
	/** Returns the length of the first barcode sequence.
	 * @return Length of first barcode in characters */
	public int barcodeLength1() {
		return barcodeLength1(barcode());
	}
	
	/** Returns the length of the second barcode sequence.
	 * @return Length of second barcode in characters */
	public int barcodeLength2() {
		return barcodeLength2(barcode());
	}
	
	/** Returns the total length of the barcode string.
	 * @return Total barcode string length, or 0 if no barcode */
	public int barcodeLength() {
		String bc=barcode();
		return bc==null ? 0 : bc.length();
	}
	
	/** Counts the number of letter characters in the barcode.
	 * @return Count of nucleotide letters in barcode sequence(s) */
	public int barcodeLetters() {
		return barcodeLetters(barcode());
	}
	
	/** Counts the number of distinct barcode sequences.
	 * @return Number of barcodes (0, 1, or 2) */
	public int numBarcodes() {
		return numBarcodes(barcode());//TODO: could add index3 count here too
	}
	
	/** Returns the third index/barcode sequence if present.
	 * @return Third index string, or null if not available */
	public String index3() {return null;}
	/** Returns the position of whitespace separating header sections.
	 * @return Index of whitespace character, or -1 if none found */
	public int whitespaceIndex() {return -1;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Finds comment separator character in a read ID string.
	 * @param id Read header string to examine
	 * @return Space (' ') or slash ('/') separator, or 0 if none found
	 */
	public static char commentSeparator(String id) {
		for(int i=0; i<id.length(); i++){
			char c=id.charAt(i);
			if(c==' ' || c== '/'){return c;}
		}
		return 0;
	}
	
	/**
	 * Extracts surface number from tile number.
	 * @param tile Tile number containing encoded surface information
	 * @return Surface number (thousands digit)
	 */
	public static final int surface(int tile) {
		return tile/1000;
	}
	
	/**
	 * Extracts swath number from tile number.
	 * @param tile Tile number containing encoded swath information
	 * @return Swath number (hundreds digit)
	 */
	public static final int swath(int tile) {
		return (tile%1000)/100;
	}
	
	/**
	 * Finds delimiter character separating multiple barcodes.
	 * @param bc Barcode string to examine
	 * @return Delimiter character (typically '+'), or 0 if no delimiter found
	 */
	public static final char barcodeDelimiter(String bc) {
		if(bc==null) {return 0;}
		for(int i=0; i<bc.length(); i++){
			char c=bc.charAt(i);
			if(!Character.isLetter(c)){return c;}
		}
		return 0;
	}
	
	/**
	 * Calculates the length of the first barcode sequence.
	 * @param bc Barcode string containing one or more sequences
	 * @return Length of first barcode before delimiter, or total length if no delimiter
	 */
	public static int barcodeLength1(String bc) {
		if(bc==null) {return 0;}
		for(int i=0; i<bc.length(); i++){
			char c=bc.charAt(i);
			if(!Tools.isLetter(c)){return i;}
		}
		return bc.length();
	}
	
	/**
	 * Calculates the length of the second barcode sequence.
	 * @param bc Barcode string containing one or more sequences
	 * @return Length of second barcode after delimiter, or 0 if no second barcode
	 */
	public static int barcodeLength2(String bc) {
		if(bc==null) {return 0;}
		for(int i=bc.length()-1; i>=0; i--){
			char c=bc.charAt(i);
			if(!Tools.isLetter(c)){return bc.length()-1-i;}
		}
		return 0;
	}
	
	/**
	 * Counts nucleotide letters in barcode string.
	 * @param bc Barcode string to count
	 * @return Number of letter characters (A, C, G, T, N)
	 */
	public static int barcodeLetters(String bc) {
		if(bc==null) {return 0;}
		int letters=0;
		for(int i=0; i<bc.length(); i++){
			char c=bc.charAt(i);
			letters+=(Tools.isLetter(c) ? 1 : 0);
		}
		return letters;
	}
	
	/**
	 * Counts distinct barcode sequences in barcode string.
	 * @param bc Barcode string to analyze
	 * @return Number of separate barcodes (0, 1, or 2)
	 */
	public static int numBarcodes(String bc) {
		return (barcodeLength1(bc)>0 ? 1 : 0)+(barcodeLength2(bc)>0 ? 1 : 0);
	}
	
	/** Creates formatted string representation of parsed header data.
	 * @return Multi-line string with tab-separated field names and values */
	public String toString() {
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

	/** Read header */
	String id;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Parse lane, tile, x, and y coordinates */
	public static boolean PARSE_COORDINATES=true;
	/** Parse the comment field for pair number, chastity filter, and barcode */
	public static boolean PARSE_COMMENT=false;
	
}