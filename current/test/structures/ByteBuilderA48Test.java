package test.structures;

import parse.LineParser1;
import parse.LineParser2;
import structures.ByteBuilder;

/** Raw A48 floating-point bit roundtrips.
 * @author Yelan
 */
public final class ByteBuilderA48Test{

	/** Verifies ByteBuilder A48 output against both tab-delimited line parsers. */
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No ByteBuilderA48Test arguments");}
		final long[] bits={0L, Long.MIN_VALUE, Double.doubleToRawLongBits(Double.MIN_VALUE),
			Double.doubleToRawLongBits(-Double.MIN_VALUE), Double.doubleToRawLongBits(1.25),
			Double.doubleToRawLongBits(-2.5), Double.doubleToRawLongBits(Double.POSITIVE_INFINITY),
			Double.doubleToRawLongBits(Double.NEGATIVE_INFINITY), 0x7ff8000000000001L};// Quiet NaN payload.
		final ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<bits.length; i++){
			if(i>0){bb.tab();}
			bb.appendDoubleA48(Double.longBitsToDouble(bits[i]));
		}
		final LineParser1 lp1=new LineParser1('\t');
		lp1.set(bb.toBytes());
		for(int i=0; i<bits.length; i++){
			check(Double.doubleToRawLongBits(lp1.parseDoubleA48(i))==bits[i], "LP1 raw double A48 roundtrip");
		}
		final LineParser2 lp2=new LineParser2('\t');
		lp2.set(bb.toBytes());
		for(int i=0; i<bits.length; i++){
			check(Double.doubleToRawLongBits(lp2.parseDoubleA48())==bits[i], "LP2 streaming raw double A48 roundtrip");
		}
		lp2.set(bb.toBytes());
		for(int i=0; i<bits.length; i++){
			check(Double.doubleToRawLongBits(lp2.parseDoubleA48(i))==bits[i], "LP2 indexed raw double A48 roundtrip");
		}
		System.out.println("BYTE_BUILDER_A48_TEST_PASS checks="+checks);
	}

	private static void check(boolean condition, String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	private static int checks;
}
