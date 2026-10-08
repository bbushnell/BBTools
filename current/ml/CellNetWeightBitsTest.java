package ml;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import parse.LineParser2;
import structures.ByteBuilder;

/** Independent bit-mask oracles for18/24-bit dense/sparse CellNet I/O.
 * @author Brian, Yoimiya
 */
public final class CellNetWeightBitsTest {

	/** Exercises exact legacy bits, truncation boundaries, preserved metadata and malformed payloads. */
	public static void main(String[] args){
		boolean enabled=false; assert(enabled=true);
		if(!enabled || args.length!=0){throw new IllegalArgumentException("Run this fixture with -ea and no arguments");}
		for(boolean dense:new boolean[]{true, false}){for(int precision:new int[]{24, 18}){
			final CellNet legacy=parse(fixture(dense));
			check(legacy, 32);
			final String before=legacy.toBytes().toString();
			if(before.contains("#weightbits")){throw new AssertionError("Legacy output gained a reduced-precision header");}
			check(parse(before), 32);
			legacy.setWeightBits(precision);
			final String encoded=legacy.toBytes().toString();
			final String header="#weightbits "+precision;
			if(!encoded.contains(header+"\n")){throw new AssertionError("Reduced edge payload requires its critical header");}
			for(String line:before.split("\n")){
				if(line.startsWith("#input")){assert(encoded.contains(line+"\n")) : "Input preprocessing constants must remain exact32-bit";}
			}
			for(String line:encoded.split("\n")){
				if(line.startsWith("C") || line.startsWith("W")){
					final String[] fields=line.split(" ");
					for(int i=3; i<fields.length; i++){
						assert(fields[i].length()<=precision/6) : "One fewer A48 symbol must encode18-bit edges";
					}
				}
			}
			final CellNet reduced=parse(encoded);
			assert(reduced.weightBits()==precision) : "Parser must retain the critical edge precision";
			check(reduced, precision);
			check(reduced.copy(false), precision);
			check(parse(reduced.toBytes().toString()), precision);
			check(CellNetParser.loadInferenceFromLines(lines(encoded)), precision);
			reduced.setWeightBits(32);
			check(parse(reduced.toBytes().toString()), precision);
			reject(encoded.replace(header, "#weightbits 25"));
			reject(encoded.replace(header, header+"junk"));
			reject(encoded.replace(header, "#weightbitsX "+precision));
			reject(encoded.replace(header, header+"\n"+header));
			reject(encoded.replace("#coding A48\n", ""));
			reject(encoded.replace("#coding A48", "#coding decimal"));
			reject(encoded.replace("#coding A48", "#coding A48\n#coding A48"));
			reject(encoded.replace("#coding A48", "#codingX A48"));
			legacy.net[1][0].weights[0]=Float.POSITIVE_INFINITY;
			boolean rejected=false;
			try{legacy.toBytes();}catch(IllegalArgumentException expected){rejected=true;}
			assert(rejected) : "Nonfinite source weights cannot be exported as reduced-precision data";
		}}
		checkZeroAbsentHeaders();
		checkSetFromMetadata();
		checkDecimalFloor();
		for(int precision:new int[]{24, 18}){
			final String[] invalid=precision==24 ? new String[]{"", "00000", "/", "p", "Ooo0", "Oooo", "Oh00", "oh00"} :
				new String[]{"", "0000", "/", "p", "Oh0", "oh0", "Ooo"};
			for(String token:invalid){
				final LineParser2 parser=new LineParser2(' '); parser.set(token.getBytes(StandardCharsets.US_ASCII));
				boolean rejected=false;
				try{
					if(precision==24){parser.parseFloatA48Truncated24();}else{parser.parseFloatA48Truncated(18);}
				}catch(IllegalArgumentException expected){rejected=true;}
				assert(rejected) : "Malformed/out-of-range/nonfinite reduced token accepted: "+precision+" "+token;
			}
		}
		System.out.println("CellNetWeightBitsTest PASS dense/sparse, legacy32, signed zero, subnormals, roundtrip/copy, exact bias/normalization and malformed headers/tokens; precision18_24");
	}

	/** Checkpoint replacement must copy precision and zero-edge meaning, not inherit the destination's metadata. */
	private static void checkSetFromMetadata(){
		for(boolean dense:new boolean[]{true, false}){
			for(boolean zeroAbsent:new boolean[]{false, true}){
				for(int precision:new int[]{18, 24, 32}){
					final CellNet source=parse(fixture(dense)), destination=parse(fixture(dense));
					source.setWeightBits(precision);
					source.setZeroAbsent(zeroAbsent);
					destination.setWeightBits(precision==18 ? 24 : 18);
					destination.setZeroAbsent(!zeroAbsent);
					destination.setFrom(source, false);
					assert(destination.weightBits()==precision && destination.zeroAbsent()==zeroAbsent) :
						"setFrom must preserve the source wire precision and active-zero contract";
					assert(destination.toBytes().toString().equals(source.toBytes().toString())) :
						"Replacing a same-shape checkpoint must preserve its serialized edge semantics";
				}
			}
		}
	}

	/** Checks raw bits independently of A48 parsing/formatting arithmetic. */
	private static void check(CellNet net, int precision){
		final Cell cell=net.net[1][0];
		assert(cell.weights.length==BITS.length) : "Reduced encoding must preserve every stored edge including zero";
		assert(Float.floatToRawIntBits(cell.bias())==BIAS) : "Bias is not an edge and must retain all32 bits";
		for(int i=0; i<BITS.length; i++){
			final int expected=expectedReduced(BITS[i], precision);
			assert(Float.floatToRawIntBits(cell.weights[i])==expected) : "Raw-bit oracle mismatch at edge "+i;
			if(cell.inputs!=null){assert(cell.inputs[i]==i) : "Sparse edge identity changed during encoding";}
		}
		final float[] input=new float[BITS.length];
		java.util.Arrays.fill(input, Float.intBitsToFloat(MEAN));
		net.applyInput(input);
		for(float value:net.values[0]){assert(value==0) : "Exact input mean must still normalize to zero";}
		java.util.Arrays.fill(input, Float.intBitsToFloat(MEAN)+1f);
		final float expected=(input[0]-Float.intBitsToFloat(MEAN))*Float.intBitsToFloat(0x3f801337);
		net.applyInput(input);
		for(float value:net.values[0]){
			assert(Float.floatToRawIntBits(value)==Float.floatToRawIntBits(expected)) : "Inverse deviations must remain exact float32 too";
		}
	}

	/** Checks decimal edge output on both sides of the exact half-quantum boundary. */
	private static void checkDecimalFloor(){
		final boolean oldCoding=CellNet.codingA48Out, oldDense=CellNet.DENSE;
		try{
			CellNet.codingA48Out=false;
			CellNet.DENSE=true;
			final float below=floatBelow(HALF_DECIMAL_QUANTUM);
			final float above=floatAtOrAbove(HALF_DECIMAL_QUANTUM);
			final CellNet net=new CellNet(new int[]{6, 1}, 2, 1f, 0f, 1, new ArrayList<String>());
			net.net[1][0].function=Function.getFunction(Function.LINEAR);
			net.net[1][0].bias=0;
			net.net[1][0].weights=new float[]{below, -below, above, -above, 0.000002f, -0.000002f};
			net.net[1][0].deltas=new float[net.net[1][0].weights.length];
			net.setWeightBits(32);
			final String text=net.toBytes().toString();
			String cLine=null;
			for(String line:text.split("\n")){
				if(line.startsWith("C")){cLine=line;}
			}
			if(cLine==null){throw new AssertionError("Decimal fixture did not emit a dense C line");}
			final String[] fields=cLine.split(" ");
			assert("0.000001".equals(fields[3])) : "Positive below-half edge rounded to zero: "+fields[3];
			assert("-0.000001".equals(fields[4])) : "Negative below-half edge rounded to zero: "+fields[4];
			assert("0.000001".equals(fields[5])) : "Positive half-or-above edge changed: "+fields[5];
			assert("-0.000001".equals(fields[6])) : "Negative half-or-above edge changed: "+fields[6];
			assert("0.000002".equals(fields[7])) : "Ordinary positive decimal edge changed: "+fields[7];
			assert("-0.000002".equals(fields[8])) : "Ordinary negative decimal edge changed: "+fields[8];
		}finally{
			CellNet.codingA48Out=oldCoding;
			CellNet.DENSE=oldDense;
		}
	}

	/** Checks the explicit dense zero-mask contract and live-edge count. */
	private static void checkZeroAbsentHeaders(){
		final String dense=fixture(true);
		final String header="#zeroabsent true\n#liveedges "+liveEdges()+"\n";
		final String flagged=dense.replace("#coding A48\n", header+"#coding A48\n");
		assert(!parse(dense).zeroAbsent()) : "Unflagged dense files must keep legacy zero-is-active semantics";
		final CellNet parsed=parse(flagged);
		assert(parsed.zeroAbsent()) : "Dense #zeroabsent header was not retained";
		final String roundTrip=parsed.toBytes().toString();
		assert(roundTrip.contains(header)) : "Dense #zeroabsent/#liveedges headers were not round-tripped";
		reject(flagged.replace("#liveedges "+liveEdges()+"\n", ""));
		reject(flagged.replace("#liveedges "+liveEdges(), "#liveedges "+(liveEdges()-1)));
		reject(flagged.replace("#liveedges "+liveEdges(), "#liveedges "+liveEdges()+"\n#liveedges "+liveEdges()));
		reject(flagged.replace("#liveedges "+liveEdges(), "#liveedges "));
		reject(flagged.replace("#liveedges "+liveEdges(), "#liveedges -1"));
		reject(flagged.replace("#liveedges "+liveEdges(), "#liveedges 1x"));
		reject(flagged.replace("#liveedges "+liveEdges(), "#liveedges 9223372036854775808"));
		reject(flagged.replace("#liveedges "+liveEdges(), "#liveedges -1\n#liveedges "+liveEdges()));
		reject(flagged.replace("#zeroabsent true", "#zeroabsent true\n#zeroabsent true"));
		reject(flagged.replace("#zeroabsent true\n", "#zeroabsent false\n"));
		reject(duplicateDenseRow(flagged).replace("#liveedges "+liveEdges(), "#liveedges "+(2*liveEdges())));
		reject(dense.replace("#coding A48\n", "#liveedges -1\n#coding A48\n"));
		reject(fixture(false).replace("#coding A48\n", "#zeroabsent false\n#coding A48\n"));
		reject(fixture(false).replace("#coding A48\n", header+"#coding A48\n"));
	}

	/** Number of nonzero fixture bits. */
	private static int liveEdges(){
		int count=0;
		for(int bits : BITS){count+=(bits==0 || bits==0x80000000 ? 0 : 1);}
		return count;
	}

	/** Duplicates the sole dense edge row while allowing #liveedges to be forged. */
	private static String duplicateDenseRow(final String text){
		final int start=text.indexOf("\nC")+1;
		final int end=text.indexOf('\n', start)+1;
		return text.substring(0, end)+text.substring(start, end)+text.substring(end);
	}

	/** Expected retained raw float bits after reduced-A48 flooring. */
	private static int expectedReduced(final int bits, final int precision){
		if(precision==32){return bits;}
		final int mask=precision==18 ? 0xffffc000 : 0xffffff00;
		int expected=bits&mask;
		final int payload=bits>>>(32-precision);
		final int signless=payload&((1<<(precision-1))-1);
		if(bits!=0 && bits!=0x80000000 && signless==0){expected|=1<<(32-precision);}
		return expected;
	}

	/** Largest positive float below target. */
	private static float floatBelow(final double target){
		int bits=Float.floatToRawIntBits((float)target);
		while(Float.intBitsToFloat(bits)>=target){bits--;}
		while(Float.intBitsToFloat(bits+1)<target){bits++;}
		return Float.intBitsToFloat(bits);
	}

	/** Smallest positive float at or above target. */
	private static float floatAtOrAbove(final double target){
		int bits=Float.floatToRawIntBits((float)target);
		while(Float.intBitsToFloat(bits)<target){bits++;}
		while(Float.intBitsToFloat(bits-1)>=target){bits--;}
		return Float.intBitsToFloat(bits);
	}

	/** Creates an actual native32-bit file, including positive/negative zero and subnormal boundaries. */
	private static String fixture(boolean dense){
		final ByteBuilder text=new ByteBuilder();
		text.append("##bbnet\n#version 1\n#concise\n").append(dense ? "#dense\n" : "#sparse\n")
			.append("#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims ").append(BITS.length).append(" 1\n#coding A48\n");
		text.append("#inputmean_a48");
		for(int ignored:BITS){text.space().appendFloatA48(Float.intBitsToFloat(MEAN));}
		text.nl().append("#inputinversestd_a48");
		for(int ignored:BITS){text.space().appendFloatA48(Float.intBitsToFloat(0x3f801337));}
		text.nl();
		final int id=BITS.length+1;
		if(!dense){text.append('I').append(id); for(int i=0; i<BITS.length; i++){text.space().append(i);} text.nl();}
		text.append(dense ? 'C' : 'W').append(id).append(" LINEAR ").appendFloatA48(Float.intBitsToFloat(BIAS));
		for(int bits:BITS){text.space().appendFloatA48(Float.intBitsToFloat(bits));}
		return text.nl().toString();
	}

	/** Converts tiny fixture text to the production byte-line parser input. */
	private static ArrayList<byte[]> lines(String text){
		final ArrayList<byte[]> result=new ArrayList<byte[]>();
		for(String line:text.split("\n")){result.add(line.getBytes(StandardCharsets.US_ASCII));}
		return result;
	}
	private static CellNet parse(String text){return CellNetParser.loadFromLines(lines(text));}
	/** Critical malformed headers must throw even outside assertions-only validation. */
	private static void reject(String text){
		try{parse(text);}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed reduced-precision header accepted");
	}
	private static final int[] BITS={0x3f812345, 0xbf854321, 0, 0x80000000, 1, 0x80000001,
		0x000000ff, 0x00000100, 0x800000ff, 0x80000100, 0x007fffff, 0x807fffff, 0x7f7fffff,
		0x00003fff, 0x00004000, 0x80003fff, 0x80004000, 0x3f803fff, 0x3f804000};
	private static final int BIAS=0x3f812347, MEAN=0x3f800123;
	private static final double HALF_DECIMAL_QUANTUM=0.0000005;
}
