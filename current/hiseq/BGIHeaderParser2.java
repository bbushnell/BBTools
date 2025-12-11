package hiseq;

import shared.LineParserS4Reverse;
import structures.ByteBuilder;

/**
 * Parser for BGI sequencing platform read headers with unknown prefix format.
 * Uses reverse parsing strategy to handle variable prefix lengths in BGI headers.
 * Does not support comment parsing by default due to reverse parsing approach.
 * BGI headers follow format: prefixLlaneClusterRreadNumber/pairNumber
 *
 * @author Brian Bushnell
 * @date May 6, 2024
 */
public class BGIHeaderParser2 extends ReadHeaderParser {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args) {
		BGIHeaderParser2 ihp=new BGIHeaderParser2();
		ihp.test(args.length>0 ? args[0] : null);
		ihp.parse(args[0]);
		System.err.println("toIllumina: "+ihp.toIllumina("ACGTACGT"));
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
	
	public BGIHeaderParser2 parse(String id_) {
		id=id_;
		extra=null;
		if(PARSE_EXTRA) {//Handles comments, but it's slow
			int idx=firstWhitespace(id);
			if(idx>=0) {
				extra=id.substring(idx+1);
				id=id.substring(0, idx);
			}
		}
		lp.set(id);
		return this;
	}
	
	private static int firstWhitespace(String s) {
		for(int i=0; i<s.length(); i++) {
			if(Character.isWhitespace(s.charAt(i))){
				return i;
			}
		}
		return -1;
	}

	//@LH00223:28:22GLGMLT3:1:1101:5928:1016 1:N:0:CTGCTTGGTT+CTAACGACAG (NovaseqX)
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
	
	@Override
	public String sample() {
		String s=lp.terms()<=1 ? null : lp.parseString(1);
		return s!=null ? s : "SA";
	}
	
	@Override
	public String machine() {
		return "CG";//Complete Genomics
	}

	@Override
	public int run() {
		return 0;
	}
	
	@Override
	public String flowcell() {
		String s=lp.terms()<=0 ? null : lp.parseString(0);
//		String s=lp.terms()<=1 ? null : lp.parseString(1);
		return s!=null ? s : "FC";
	}

	/**
	 * Extracts lane number from BGI header.
	 * Parses the third field (index 2) as lane identifier.
	 * @return Lane number as integer
	 */
	@Override
	public int lane() {return lp.parseInt(2);}

	/**
	 * Extracts tile number from BGI header using complex parsing.
	 * Parses from field 4, starting at position 3, using base 10.
	 * @return Tile number as integer
	 */
	@Override
	public int tile() {return lp.parseInt(4, 3, 10);}

	/**
	 * Extracts X coordinate from BGI header.
	 * Uses the fourth field (index 3) as X position.
	 * @return X coordinate as integer
	 */
	@Override
	public int xPos() {return lp.parseInt(3);}

	/**
	 * Extracts Y coordinate from BGI header using substring parsing.
	 * Parses from field 4, positions 0-3 as Y coordinate.
	 * @return Y coordinate as integer
	 */
	@Override
	public int yPos() {return lp.parseInt(4, 0, 3);}

	/**
	 * Extracts pair number from BGI header.
	 * Returns the first character of the sixth field (index 5).
	 * @return Pair code character ('1' or '2')
	 */
	@Override
	public char pairCode() {return lp.parseChar(5, 0);}

	/**
	 * Returns chastity filter status for BGI reads.
	 * BGI format does not include chastity information, returns 'N'.
	 * @return Always returns 'N' (not filtered)
	 */
	@Override
	public char chastityCode() {return 'N';}

	/**
	 * Returns control bits value for BGI headers.
	 * BGI format does not include control bits, returns 0.
	 * @return Always returns 0
	 */
	@Override
	public int controlBits() {return 0;}

	/**
	 * Returns barcode sequence from BGI header.
	 * BGI format does not include inline barcodes, returns null.
	 * @return Always returns null
	 */
	@Override
	public String barcode() {
		return null;
	}

	/**
	 * Returns extra comment information from BGI header.
	 * Returns the comment portion extracted during parsing if PARSE_EXTRA enabled.
	 * @return Comment string or null if no comments or parsing disabled
	 */
	@Override
	public String extra() {
//		return lp.terms()<=6 ? null : lp.parseString(6);
		return extra;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	private String extra=null;
	private final LineParserS4Reverse lp=new LineParserS4Reverse("_LCR/");
	private final ByteBuilder bb=new ByteBuilder(64);
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean PARSE_EXTRA=false;
	
}
