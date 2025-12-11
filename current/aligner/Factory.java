package aligner;

import java.util.Arrays;

import shared.Tools;

/**
 * Factory class for creating different types of sequence aligners and encoding sequences.
 * Provides a unified interface for instantiating various alignment algorithms and
 * converting DNA sequences to numeric encodings for alignment processing.
 * @author Brian Bushnell
 */
public class Factory {

	public static IDAligner makeIDAligner() {
		return makeIDAligner(type);
	}

	/**
	 * Creates an IDAligner of the specified type.
	 * Supports multiple alignment algorithms including glocal, banded, drifting,
	 * wobble, quantum, crosscut, single-state, and wavefront aligners.
	 *
	 * @param type The aligner type constant (GLOCAL, BANDED, DRIFTING, etc.)
	 * @return A new IDAligner instance of the specified type
	 */
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
	
	/**
	 * Encodes a DNA sequence as long array with padding to specified alignment.
	 * Converts DNA bases to numeric codes and pads the result to a multiple
	 * of the specified value (assumed to be a power of 2).
	 *
	 * @param in Input DNA sequence as byte array
	 * @param nCode Code to use for ambiguous bases and padding
	 * @param pad Padding alignment value (must be power of 2)
	 * @return Long array with encoded sequence padded to multiple of pad
	 */
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
	
	/**
	 * Encodes a DNA sequence as long array without padding.
	 * Converts DNA bases to numeric codes using the standard encoding table.
	 *
	 * @param in Input DNA sequence as byte array
	 * @param nCode Code to use for ambiguous bases
	 * @return Long array with encoded sequence
	 */
	public static final long[] encodeLong(byte[] in, byte nCode) {
		long[] out=new long[in.length];
		for(int i=0; i<in.length; i++) {
			final byte code=codes[in[i]];
			out[i]=(code<=8 ? code : nCode);
		}
		return out;
	}
	
	/**
	 * Encodes a DNA sequence as int array.
	 * Converts DNA bases to numeric codes using the standard encoding table.
	 *
	 * @param in Input DNA sequence as byte array
	 * @param nCode Code to use for ambiguous bases
	 * @return Int array with encoded sequence
	 */
	public static final int[] encodeInt(byte[] in, byte nCode) {
		int[] out=new int[in.length];
		for(int i=0; i<in.length; i++) {
			final byte code=codes[in[i]];
			out[i]=(code<=8 ? code : nCode);
		}
		return out;
	}
	
	/**
	 * Encodes a DNA sequence as byte array.
	 * Converts DNA bases to numeric codes using the standard encoding table.
	 *
	 * @param in Input DNA sequence as byte array
	 * @param nCode Code to use for ambiguous bases
	 * @return Byte array with encoded sequence
	 */
	public static final byte[] encodeByte(byte[] in, byte nCode) {
		byte[] out=new byte[in.length];
		for(int i=0; i<in.length; i++) {
			final byte code=codes[in[i]];
			out[i]=(code<=8 ? code : nCode);
		}
		return out;
	}

	/**
	 * Sets the default aligner type from a string name.
	 * Accepts aligner type names like "GLOCAL", "BANDED", "QUANTUM", etc.
	 * @param b String name of the aligner type (case insensitive)
	 * @return The numeric type constant, or current type if input is null
	 */
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
	/**
	 * Creates a lookup table for DNA base encoding.
	 * Maps A=1, C=2, G=4, T/U=8, with all other characters mapped to nCode.
	 * Uses bit-based encoding where each base has a unique power-of-2 value.
	 *
	 * @param nCode Code to assign to ambiguous or invalid bases
	 * @return Byte array lookup table for ASCII characters to numeric codes
	 */
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
