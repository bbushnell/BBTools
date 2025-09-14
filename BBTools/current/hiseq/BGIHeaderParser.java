package hiseq;

import shared.LineParserS4;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date April 5, 2024
 *
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
	
	public BGIHeaderParser parse(String id_) {
		id=id_;
		lp.set(id_);
		return this;
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
		return lp.terms()<=0 ? null : lp.parseString(0);
	}
	
	@Override
	public String machine() {
		return null;
	}

	@Override
	public int run() {
		return 0;
	}
	
	@Override
	public String flowcell() {
		return lp.terms()<=1 ? null : lp.parseString(1);
	}

	@Override
	public int lane() {return lp.parseInt(2);}

	@Override
	public int tile() {return lp.parseInt(4, 3, 10);}

	@Override
	public int xPos() {return lp.parseInt(3);}

	@Override
	public int yPos() {return lp.parseInt(4, 0, 3);}

	@Override
	public char pairCode() {return lp.parseChar(5, 0);}

	@Override
	public char chastityCode() {return 'N';}

	@Override
	public int controlBits() {return 0;}

	@Override
	public String barcode() {
		return null;
	}

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
