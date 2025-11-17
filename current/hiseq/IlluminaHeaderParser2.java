package hiseq;

import shared.LineParser;
import shared.LineParserS3;
import structures.ByteBuilder;

/**
 * Faster version of IlluminaHeaderParser using LineParser.
 * @author Brian Bushnell
 * @date April 3, 2024
 *
 */
public class IlluminaHeaderParser2 extends ReadHeaderParser {
	
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
	//	@Axxxxx - NovaSeq
	//	@Vxxxxx = NextSeq 2000
	//	@AAxxxxx - NextSeq 2000 P1/P2/P3
	//	@Hxxxxxx - NovaSeq S1/S2/S4
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
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Test method for header parsing functionality.
	 * @param args Command-line arguments, first argument used as test header */
	public static void main(String[] args) {
		IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
		ihp.test(args.length>0 ? args[0] : null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses an Illumina sequencing read header and initializes internal state.
	 * Tokenizes the header using colon delimiter and locates whitespace separator.
	 * @param id_ The complete read header string to parse
	 * @return This parser instance for method chaining
	 */
	public IlluminaHeaderParser2 parse(String id_) {
		id=id_;
		lp.set(id_);
		whitespaceIndex=lp.indexOfWhitespace();
		return this;
	}
	
	/**
	 * Determines if this header can be shortened by removing redundant information.
	 * A header can be shrunk if it looks valid but is not already in shrunk format.
	 * @return true if header can be compressed to save space
	 */
	public boolean canShrink() {
		return looksValid() && !looksShrunk();
	}
	
	/**
	 * Validates header structure by checking term count and whitespace position.
	 * Valid headers have 8+ terms with whitespace at position 6 or 7.
	 * @return true if header appears to be valid Illumina format
	 */
	public boolean looksValid() {
		return(lp.terms()>=8 && whitespaceIndex>=6 && whitespaceIndex<=7);
	}
	
	/**
	 * Checks if header is already in compressed format.
	 * Shrunk headers have >3 terms with the third term starting at position 2.
	 * @return true if header is already compressed
	 */
	public boolean looksShrunk() {
		return(lp.terms()>3 && lp.bounds().get(2)==2);
	}
	
	@Override
	public String machine() {
		return lp.terms()<=0 ? null : lp.parseString(0);
	}

	@Override
	public String sample() {
		return null;
	}

	@Override
	public int run() {
		return lp.parseInt(1);
	}

	@Override
	public String flowcell() {
		return lp.parseString(2);
	}

	@Override
	public int lane() {return lp.parseInt(3);}

	@Override
	public int tile() {return lp.parseInt(4);}

	@Override
	public int xPos() {return lp.parseInt(5);}

	@Override
	public int yPos() {return lp.parseInt(6);}

	@Override
	public char pairCode() {return lp.parseChar(whitespaceIndex+1, 0);}

	@Override
	public char chastityCode() {return lp.parseChar(whitespaceIndex+2, 0);}

	@Override
	public int controlBits() {return lp.parseInt(whitespaceIndex+3);}

	@Override
	public String barcode() {
		return lp.terms()<=whitespaceIndex+4 ? null : lp.parseString(whitespaceIndex+4);
	}

	@Override
	public String index3() {
		return whitespaceIndex<7 ? null : lp.parseString(7);
	}

	@Override
	public int whitespaceIndex() {
		return whitespaceIndex;
	}

	@Override
	public String extra() {
		return lp.terms()<=whitespaceIndex+5 ? null : lp.parseString(whitespaceIndex+5);
	}
	
	/**
	 * Appends a specific term from the parsed header to a ByteBuilder.
	 * Provides efficient access to individual header components.
	 *
	 * @param bb ByteBuilder to append to
	 * @param term Term index to append
	 * @return The modified ByteBuilder for method chaining
	 */
	public ByteBuilder appendTerm(ByteBuilder bb, int term) {
		return lp.appendTerm(bb, term);
	}
	
	/**
	 * Appends coordinate information in lane:tile:x:y format to ByteBuilder.
	 * Creates standardized coordinate string representation.
	 * @param bb ByteBuilder to append coordinates to
	 * @return The modified ByteBuilder with coordinate information
	 */
	public ByteBuilder appendCoordinates(ByteBuilder bb) {
		return bb.append(lane()).colon().append(tile()).colon()
		.append(xPos()).colon().append(yPos());
	}
	
	/**
	 * Encodes coordinate information into a single long value using bit shifting.
	 * Packs lane (upper bits), tile (17 bits), x-position (20 bits),
	 * and y-position (20 bits) for efficient storage and comparison.
	 * @return Encoded coordinate value as long integer
	 */
	public long encodeCoordinates() {
		long x=lane();
		x=(x<<17)^tile();
		x=(x<<20)^xPos();
		x=(x<<20)^yPos();
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Line parser configured with colon delimiter for header tokenization */
	private final LineParserS3 lp=new LineParserS3(':');
	/** Returns the internal line parser for direct access */
	public LineParser lp() {return lp;}
	/** Index position of whitespace separator in tokenized header, -1 if unset */
	int whitespaceIndex=-1;
	
}
