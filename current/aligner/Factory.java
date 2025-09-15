package aligner;

import java.util.Arrays;

import shared.Tools;

public class Factory {

	public static IDAligner makeIDAligner() {
		return makeIDAligner(type);
	}

	public static IDAligner makeIDAligner(int type) {
		if(type==GLOCAL) {return new GlocalAligner();}
		if(type==BANDED) {return new BandedAligner();}
		if(type==DRIFTING) {return new DriftingAligner();}
		if(type==WOBBLE) {return new WobbleAligner();}
		if(type==QUANTUM) {return new QuantumAligner();}
		if(type==CROSSCUT) {return new CrossCutAligner();}
		if(type==SSA2) {return new SingleStateAlignerFlat2();}
		if(type==SSA3) {return new SingleStateAlignerFlat3();}
		if(type==WAVE) {return new WaveFrontAligner();}
		assert(false) : type;
		return null;
	}
	
	/** Pads to a multiple of pad, assuming that is a power of 2 */
	public static final long[] encodeLong(byte[] in, byte nCode, int pad) {
		final int len=((in.length+pad-1)&~(pad-1));
		long[] out=new long[len];
		for(int i=0; i<in.length; i++) {
			final byte code=codes[in[i]];
			out[i]=(code<=8 ? code : nCode);
		}
		for(int i=in.length; i<out.length; i++) {out[i]=nCode;}
		return out;
	}
	
	public static final long[] encodeLong(byte[] in, byte nCode) {
		long[] out=new long[in.length];
		for(int i=0; i<in.length; i++) {
			final byte code=codes[in[i]];
			out[i]=(code<=8 ? code : nCode);
		}
		return out;
	}
	
	public static final int[] encodeInt(byte[] in, byte nCode) {
		int[] out=new int[in.length];
		for(int i=0; i<in.length; i++) {
			final byte code=codes[in[i]];
			out[i]=(code<=8 ? code : nCode);
		}
		return out;
	}
	
	public static final byte[] encodeByte(byte[] in, byte nCode) {
		byte[] out=new byte[in.length];
		for(int i=0; i<in.length; i++) {
			final byte code=codes[in[i]];
			out[i]=(code<=8 ? code : nCode);
		}
		return out;
	}

	public static int setType(String b) {
		if(b==null) {return type;}
		return type=Tools.find(b.toUpperCase(), types);
	}

	public static final int GLOCAL=1, BANDED=2, DRIFTING=3, 
			WOBBLE=4, QUANTUM=5, CROSSCUT=6, SSA2=7, SSA3=8, WAVE=9;
	public static final String[] types={"NULL", "GLOCAL", "BANDED", "DRIFTING", 
			"WOBBLE", "QUANTUM", "CROSSCUT", "SSA2", "SSA3", "WAVE"};
	public static int type=QUANTUM;
	
	public static final byte[] codes=makeCodes((byte)(15+16));
	public static final byte[] makeCodes(byte nCode) {
		byte[] codes=new byte[128];
		Arrays.fill(codes, nCode);
		codes['A']=codes['a']=1;
		codes['C']=codes['c']=2;
		codes['G']=codes['g']=4;
		codes['T']=codes['t']=8;
		codes['U']=codes['u']=8;
		return codes;
	}
	
}
