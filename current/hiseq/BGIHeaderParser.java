package hiseq;

import parse.LineParserS4;
import structures.ByteBuilder;

/**
 * Parses header metadata for BGI (Beijing Genomics Institute) sequencing data files.
 * Extends ReadHeaderParser to extract specific metadata components from BGI sequencing
 * read headers, including parsing machine, run, flowcell, lane, tile, and other
 * sequencing parameters. Uses LineParserS4 to tokenize header strings with "_LCR/\t"
 * delimiters and provides conversion to Illumina-style format.
 *
 * Expected header format: v300056266_run28L3C001R0010057888/1
 *
 * @author Brian Bushnell
 * @date April 5, 2024
 */
public class BGIHeaderParser extends ReadHeaderParser {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args) {
		BGIHeaderParser ihp=new BGIHeaderParser();
		ihp.test(args.length>0 ? args[0] : null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Expected Format       ----------------*/
	/*--------------------------------------------------------------*/
	
	//v300056266_run28L3C001R0010057888/1
	//20A_V100002704L1C001R012000000/1
	//E200008112L1C001R00100063962/1
	
	//split: [v300056266, run28, 3, 001, 0010057888, 1]
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses a BGI header string and primes the internal tokenizer.
	 * @param id_ Header string to parse
	 * @return this parser for chaining
	 */
	public BGIHeaderParser parse(String id_) {
		id=id_;
		lp.set(id_);
		return this;
	}

	//@LH00223:28:22GLGMLT3:1:1101:5928:1016 1:N:0:CTGCTTGGTT+CTAACGACAG (NovaseqX)
	/**
	 * Converts the parsed BGI header to Illumina-style format (machine:run:flowcell:lane:tile:x:y pairCode:N:controlBits:barcode).
	 * Appends optional barcode and extra metadata if present.
	 * @param barcode Optional barcode sequence
	 * @return Illumina-formatted header string
	 */
	public String toIllumina(String barcode) {
		bb.clear();
		bb.append(machine()).colon();
		bb.append(run()).colon();
		bb.append(flowcell()).colon();
		bb.append(lane()).colon();
		bb.append(tile()).colon();
		bb.append(xPos()).colon();
		bb.append(yPos()).space();
		bb.append(pairCode()).colon();
		bb.append('N').colon();
		bb.append(controlBits()).colon();
		if(barcode!=null) {bb.append(barcode);}
		String ex=extra();
		if(ex!=null) {bb.tab().append(ex);}
		return bb.toString();
	}
	
	/** Returns the sample/machine token (first term) or null if unavailable. */
	@Override
	public String sample() {
		return lp.terms()<=0 ? null : lp.parseString(0);
	}
	
	/**
	 * Returns the machine identifier (BGI headers do not provide one; returns null).
	 */
	@Override
	public String machine() {
		return null;
	}

	/** Returns the run number (BGI headers do not provide one; returns 0). */
	@Override
	public int run() {
		return 0;
	}
	
	/**
	 * Extracts the flowcell identifier from the BGI header.
	 * Returns the second parsed token which represents the flowcell ID.
	 * @return Flowcell identifier string, or null if not available
	 */
	@Override
	public String flowcell() {
		return lp.terms()<=1 ? null : lp.parseString(1);
	}

	/** Extracts the lane number from the BGI header.
	 * @return Lane number parsed from the third token */
	@Override
	public int lane() {return lp.parseInt(2);}

	/**
	 * Extracts the tile number from the BGI header.
	 * Parses a substring from the fifth token, starting at position 3 with base 10.
	 * @return Tile number
	 */
	@Override
	public int tile() {return lp.parseInt(4, 3, 10);}

	/** Extracts the x-coordinate position from the BGI header.
	 * @return X-coordinate parsed from the fourth token */
	@Override
	public int xPos() {return lp.parseInt(3);}

	/**
	 * Extracts the y-coordinate position from the BGI header.
	 * Parses a substring from the fifth token, starting at position 0 with length 3.
	 * @return Y-coordinate
	 */
	@Override
	public int yPos() {return lp.parseInt(4, 0, 3);}

	/**
	 * Extracts the pair code from the BGI header.
	 * Returns the first character of the sixth token.
	 * @return Pair code character (typically '1' or '2')
	 */
	@Override
	public char pairCode() {return lp.parseChar(5, 0);}

	/**
	 * Gets the chastity code from the BGI header.
	 * BGI headers do not contain chastity information.
	 * @return Always returns 'N' for BGI headers
	 */
	@Override
	public char chastityCode() {return 'N';}

	/**
	 * Gets the control bits from the BGI header.
	 * BGI headers do not contain control bit information.
	 * @return Always returns 0 for BGI headers
	 */
	@Override
	public int controlBits() {return 0;}

	/**
	 * Gets the barcode sequence from the BGI header.
	 * BGI headers do not contain embedded barcode information.
	 * @return Always returns null for BGI headers
	 */
	@Override
	public String barcode() {
		return null;
	}

	/**
	 * Extracts any extra information from the BGI header.
	 * Returns the seventh token if available, typically containing additional metadata.
	 * @return Extra information string, or null if not available
	 */
	@Override
	public String extra() {
		return lp.terms()<=11 ? null : lp.parseString(6);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	private final LineParserS4 lp=new LineParserS4("_LCR/\t");
	private final ByteBuilder bb=new ByteBuilder(64);
	
}
