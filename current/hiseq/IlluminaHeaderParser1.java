package hiseq;

import shared.KillSwitch;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Parses Illumina sequencing read headers to extract positional and metadata
 * information from various Illumina platforms (MiSeq, HiSeq, NovaSeq).
 * This parser has been superseded by IlluminaHeaderParser2 for improved performance.
 *
 * Supports multiple header formats including:
 * - NovaSeq 6000: @VP2-06:112:H7LNDMCVY:2:2437:14181:20134
 * - MiSeq: MISEQ08:172:000000000-ABYD0:1:1101:18147:1925 1:N:0:TGGATATGCGCCAATT
 * - HiSeq: HISEQ07:419:HBFNEADXX:1:1101:1238:2072
 * - NovaSeq S: A00178:38:H5NYYDSXX:2:1101:3007:1000 1:N:0:CAACCTA+CTAGGTT
 * - NovaSeq X: @LH00223:28:22GLGMLT3:1:1101:5928:1016 1:N:0:CTGCTTGGTT+CTAACGACAG
 *
 * @author Brian Bushnell
 * @date Aug 22, 2018
 */
public class IlluminaHeaderParser1 extends ReadHeaderParser {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Test hook that enables comment parsing, constructs a parser, and parses a
	 * provided header (or a default example) to display extracted fields.
	 * @param args Optional single header string to parse
	 */
	public static void main(String[] args) {
		PARSE_COMMENT=true;
		IlluminaHeaderParser1 ihp=new IlluminaHeaderParser1();
		ihp.test(args.length>0 ? args[0] : null);
	}
	
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
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses an Illumina header string, extracting coordinates and comment fields
	 * according to static flags PARSE_COORDINATES and PARSE_COMMENT. Resets parser
	 * state before parsing and tolerates multiple platform header formats.
	 *
	 * @param id_ Header string to parse
	 * @return This parser instance for chaining
	 */
	public IlluminaHeaderParser1 parse(String id_) {
		reset(id_);
		
		try {
			commentSeparator=findCommentSeparator();
			if(PARSE_COORDINATES){parseCoordinates();}
			if(PARSE_COMMENT){parseComment();}
		} catch (Throwable e) {
			System.err.println("Trouble parsing header "+id_);
			KillSwitch.throwableKill(e);
		}
		return this;
	}
	
	/**
	 * Returns sequencing instrument identifier. Not parsed by this legacy parser
	 * and always returns null.
	 * @return null (machine name not extracted)
	 */
	@Override
	public String machine() {
		return null;
	}
	
	/** Returns sample identifier. Not parsed by this legacy parser and always null.
	 * @return null (sample not extracted) */
	@Override
	public String sample() {
		return null;
	}

	/** Returns run number; not parsed in this implementation so -1 is returned. */
	@Override
	public int run() {
		return -1;
	}

	/**
	 * Returns flowcell identifier; not parsed in this implementation so null is returned.
	 */
	@Override
	public String flowcell() {
		return null;
	}

	/** Returns the flowcell lane number */
	@Override
	public int lane() {return lane;}

	/** Returns the tile number */
	@Override
	public int tile() {return tile;}

	/** Returns the x-coordinate position within the tile (in pixels) */
	@Override
	public int xPos() {return x;}

	/** Returns the y-coordinate position within the tile (in pixels) */
	@Override
	public int yPos() {return y;}

	/** Returns the pair code ('1' for read 1, '2' for read 2) */
	@Override
	public char pairCode() {return pairCode;}

	/** Returns the chastity filter code ('Y' for fail, 'N' for pass) */
	@Override
	public char chastityCode() {return chastityCode;}

	/** Returns the control bits value from the header */
	@Override
	public int controlBits() {return controlBits;}

	/**
	 * Returns the barcode sequence from the header comment field.
	 * Requires PARSE_COMMENT flag to be enabled.
	 * @return Barcode string or null if not parsed
	 */
	@Override
	public String barcode() {
		assert(PARSE_COMMENT);
		return barcode;
	}

	/**
	 * Returns any extra information after the barcode in the comment field.
	 * Requires PARSE_COMMENT flag to be enabled.
	 * @return Extra string data or null if not present
	 */
	@Override
	public String extra() {
		assert(PARSE_COMMENT);
		return extra;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses lane, tile, x, and y coordinates from the header string.
	 * Handles various header formats including HiSeq 3000 variations.
	 * Updates lane, tile, x, and y fields.
	 */
	private void parseCoordinates(){
		pos=commentSeparator;
		goBackSeveralColons(4);
		lane=parseInt();
		if(!Tools.isDigit(id.charAt(pos))){//Hiseq 3000?  I'm not really sure what the header looks like for this block
			while(pos<limit && id.charAt(pos)!=':'){pos++;}
			pos++;
			lane=parseInt();
		}

		tile=parseInt();
		x=parseInt();
		y=parseInt();
	}
	
	/**
	 * Parses the comment field to extract pair number, chastity filter, and barcode.
	 * Extracts pairCode, chastityCode, controlBits, barcode, and extra fields.
	 * Handles both space and tab-separated extra information.
	 */
	private void parseComment(){
		pos=commentSeparator+1;
		pairCode=parseChar();
		chastityCode=parseChar();
		controlBits=parseInt();
		int idx=id.indexOf(' ', pos);
		idx=(idx>=0 ? idx : id.indexOf('\t', pos));
		if(idx<0) {
			barcode=id.substring(pos);
			extra=null;
		}else {
			barcode=id.substring(pos, idx);
			extra=id.substring(idx+1);
		}
	}
	
	/**
	 * Resets all internal fields and sets the header to parse.
	 * Initializes parsing state and clears previously extracted values.
	 * @param id_ The new header string to parse
	 */
	private void reset(String id_){
		id=id_;
		limit=(id==null ? -1 : id.length());
		pos=-1;
		commentSeparator=-1;

		lane=-1;
		tile=-1;
		x=-1;
		y=-1;
		
		pairCode='?';
		barcode=null;
		chastityCode='?';
		controlBits=-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Locates the comment separator character (' ' or '/') in the header.
	 * This separator divides the coordinate section from the comment section.
	 * @return Position of the separator, or header length if not found
	 */
	private int findCommentSeparator(){
		for(int i=0; i<limit; i++){
			char c=id.charAt(i);
			if(c==' ' || c== '/'){return i;}
		}
		return limit;
	}
	
	/**
	 * Moves the parsing position backward through the specified number of colons.
	 * Used to locate the beginning of coordinate fields in the header.
	 * @param target Number of colons to backtrack through
	 */
	private void goBackSeveralColons(int target){
		for(int colons=0; pos>=0; pos--){
			if(id.charAt(pos)==':'){
				colons++;
				if(colons==target){break;}
			}
		}
		pos++;
	}
	
	/**
	 * Parses an integer from the current position and advances past it.
	 * Reads consecutive digits and moves position to the next token.
	 * @return The parsed integer value
	 */
	private int parseInt(){
		int current=0;
		assert(Tools.isDigit(id.charAt(pos))) : id;
		while(pos<limit && Tools.isDigit(id.charAt(pos))){
			current=current*10+(id.charAt(pos)-'0');
			pos++;
		}
		pos++;
		return current;
	}
	
	/**
	 * Parses a single character from the current position and advances past delimiter.
	 * Expects character to be followed by a colon delimiter.
	 * @return The parsed character
	 */
	private char parseChar(){
		char c=id.charAt(pos);
		assert(c!=':');
		pos++;
		assert(pos>=limit || id.charAt(pos)==':');
		pos++;
		assert(pos>=limit || id.charAt(pos)!=':');
		return c;
	}
	
	public String toString() {
		ByteBuilder bb=new ByteBuilder();
		bb.append("lane:\t").append(lane).nl();
		bb.append("tile:\t").append(tile).nl();
		bb.append("x:\t").append(x).nl();
		bb.append("y:\t").append(y).nl();
		bb.append("pairnum:\t").append(pairCode).nl();
		bb.append("barcode:\t").append(barcode).nl();
		bb.append("chastity:\t").append(chastityCode).nl();
		return bb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Flowcell lane number extracted from header */
	public int lane;
	/** Tile number extracted from header */
	public int tile;
	/** X-coordinate within tile (pixels) extracted from header */
	public int x;
	/** Y-coordinate within tile (pixels) extracted from header */
	public int y;
	
	/** Pair designation: '1' for read 1, '2' for read 2 */
	public char pairCode;
	/** Chastity filter result: 'Y' for fail, 'N' for pass */
	public char chastityCode;
	/** Control bits value indicating special read properties */
	public int controlBits;
	/** Barcode sequence extracted from comment field */
	public String barcode;
	/** Additional information following the barcode in comment field */
	public String extra;
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Length of the header string being parsed */
	private int limit;
	/** Current parsing position, typically start of next token */
	private int pos;
	/** Position of separator between coordinates and comment (space or slash) */
	private int commentSeparator;
	
}
