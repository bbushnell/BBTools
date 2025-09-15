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
	
	public abstract ReadHeaderParser parse(String id);
	public abstract String machine();
	public abstract int run();
	public abstract String flowcell();
	public abstract int lane();
	public abstract int tile();
	public abstract int xPos();
	public abstract int yPos();
	public abstract char pairCode();
	public abstract char chastityCode();
	public abstract int controlBits();
	public abstract String barcode();
	public abstract String extra();
	public abstract String sample();
	
	/*--------------------------------------------------------------*/
	/*----------------       Concrete Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public final ReadHeaderParser parse(Read r){
		return parse(r.id);
	}
	
	public final String id() {return id;}
	
	public char commentSeparator() {
		for(int i=0; i<id.length(); i++){
			char c=id.charAt(i);
			if(c==' ' || c== '/'){return c;}
		}
		return 0;
	}
	
	public int surface() {
		return surface(tile());
	}
	
	public int swath() {
		return swath(tile());
	}
	
	public int pairnum() {
		return pairCode()-(int)'1';
	}
	
	public boolean chastityFail() {
		int c=chastityCode();
		assert(c=='N' || c=='Y') : c;
		return c=='Y';
	}
	
	public char barcodeDelimiter() {
		return barcodeDelimiter(barcode());
	}
	
	public int barcodeLength1() {
		return barcodeLength1(barcode());
	}
	
	public int barcodeLength2() {
		return barcodeLength2(barcode());
	}
	
	public int barcodeLength() {
		String bc=barcode();
		return bc==null ? 0 : bc.length();
	}
	
	public int barcodeLetters() {
		return barcodeLetters(barcode());
	}
	
	public int numBarcodes() {
		return numBarcodes(barcode());//TODO: could add index3 count here too
	}
	
	public String index3() {return null;}
	public int whitespaceIndex() {return -1;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static char commentSeparator(String id) {
		for(int i=0; i<id.length(); i++){
			char c=id.charAt(i);
			if(c==' ' || c== '/'){return c;}
		}
		return 0;
	}
	
	public static final int surface(int tile) {
		return tile/1000;
	}
	
	public static final int swath(int tile) {
		return (tile%1000)/100;
	}
	
	public static final char barcodeDelimiter(String bc) {
		if(bc==null) {return 0;}
		for(int i=0; i<bc.length(); i++){
			char c=bc.charAt(i);
			if(!Character.isLetter(c)){return c;}
		}
		return 0;
	}
	
	public static int barcodeLength1(String bc) {
		if(bc==null) {return 0;}
		for(int i=0; i<bc.length(); i++){
			char c=bc.charAt(i);
			if(!Tools.isLetter(c)){return i;}
		}
		return bc.length();
	}
	
	public static int barcodeLength2(String bc) {
		if(bc==null) {return 0;}
		for(int i=bc.length()-1; i>=0; i--){
			char c=bc.charAt(i);
			if(!Tools.isLetter(c)){return bc.length()-1-i;}
		}
		return 0;
	}
	
	public static int barcodeLetters(String bc) {
		if(bc==null) {return 0;}
		int letters=0;
		for(int i=0; i<bc.length(); i++){
			char c=bc.charAt(i);
			letters+=(Tools.isLetter(c) ? 1 : 0);
		}
		return letters;
	}
	
	public static int numBarcodes(String bc) {
		return (barcodeLength1(bc)>0 ? 1 : 0)+(barcodeLength2(bc)>0 ? 1 : 0);
	}
	
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