package structures;

import parse.Parse;

/** Compact nonnegative ratios: mantissa times2^(exponent-bias), no sign or
 * implicit leading bit. Six exponent bits retain a broad ratio range;12/18-bit
 * containers have6/12mantissa bits. Zero is exact. Canonical codes use the
 * smallest fitting exponent and nearest-even rounding. Positive values outside
 * the declared range fail rather than underflowing to a zero log argument.
 * A decoded ratio is a double, not a potentially overflowing integer count.
 * Immutable; encode/decode/append/read allocate nothing.
 * @author Brian, Raiden
 */
public final class UnsignedFloat{

	private UnsignedFloat(final int bits_, final int bias_){
		bits=bits_;
		bias=bias_;
		mantissaBits=bits-6;
		maximumMantissa=(1<<mantissaBits)-1;
		minimum=Math.scalb(1.0, -bias);
		maximum=Math.scalb((double)maximumMantissa, 63-bias);
		minimumNormal=Math.scalb((double)(1<<(mantissaBits-1)), -bias);
		maximumNormalRelativeError=1.0/(1<<mantissaBits);
		assert((bits==12 || bits==18) && minimum>0 && Double.isFinite(maximum))
			: "Fixed native A48 ratio formats need finite positive endpoints and two or three symbols";
	}

	/** The proposed ratio presets share the same lower normalized exponent. */
	public static UnsignedFloat forBits(final int bits){
		if(bits==12){return TWELVE;}
		if(bits==18){return EIGHTEEN;}
		throw new IllegalArgumentException("Unsigned ratio floats use12or18bits: "+bits);
	}

	/** Encodes exact zero or a finite ratio in the declared positive interval.
	 * Relative error is bounded by2^-mantissaBits above minimumNormal; below that
	 * threshold the absolute error is at most minimum/2, so relative error grows.
	 */
	public int encode(final double value){
		if(!Double.isFinite(value) || value<0 || value>maximum || value>0 && value<minimum){
			throw new IllegalArgumentException("Unsigned ratio must be zero or inside ["+minimum+", "+maximum+"]: "+value);
		}
		if(value==0){return 0;}
		int exponent=Math.max(0, Math.getExponent(value)+bias-mantissaBits+1);
		long mantissa=(long)Math.rint(Math.scalb(value, bias-exponent));
		if(mantissa>maximumMantissa){
			exponent++;
			mantissa=(long)Math.rint(Math.scalb(value, bias-exponent));
		}
		assert(exponent<=63 && mantissa>0 && mantissa<=maximumMantissa
			&& (exponent==0 || mantissa>maximumMantissa/2))
			: "Nearest rounding must retain a positive canonical mantissa and the smallest fitting exponent";
		return (exponent<<mantissaBits)|(int)mantissa;
	}

	/** Returns the exact binary ratio represented by a canonical code. */
	public double decode(final int code){
		if(code<0 || code>=(1<<bits)){
			throw new IllegalArgumentException("Unsigned ratio code exceeds its fixed-width container: "+code);
		}
		final int exponent=code>>>mantissaBits, mantissa=code&maximumMantissa;
		if(exponent>0 && mantissa<=maximumMantissa/2){
			throw new IllegalArgumentException("Unsigned ratio code has a smaller equivalent exponent: "+code);
		}
		return Math.scalb((double)mantissa, exponent-bias);
	}

	/** Appends the existing native A48 code, padded to exactly two/three symbols. */
	public ByteBuilder append(final ByteBuilder output, final int code){
		decode(code);
		assert(output!=null) : "Ratio fields append into a caller-owned reused output buffer";
		final int digits=Math.max(1, (32-Integer.numberOfLeadingZeros(code)+5)/6);
		for(int i=digits; i<bits/6; i++){output.append('0');}
		return output.appendA48(code);
	}

	/** Reads exactly one native A48 field; caller checks the complete row width. */
	public int read(final byte[] row, final int from){
		final int width=bits/6;
		if(row==null || from<0 || from>row.length-width){
			throw new IllegalArgumentException("Unsigned ratio input needs its complete fixed-width A48 field");
		}
		for(int i=from; i<from+width; i++){
			if(row[i]<48 || row[i]>111){throw new IllegalArgumentException("Invalid native A48 ratio symbol");}
		}
		final int code=(int)Parse.parseLongA48(row, from, from+width);
		decode(code);
		return code;
	}

	public final int bits, bias, mantissaBits, maximumMantissa;
	public final double minimum, maximum, minimumNormal, maximumNormalRelativeError;
	private static final UnsignedFloat TWELVE=new UnsignedFloat(12, 32);
	private static final UnsignedFloat EIGHTEEN=new UnsignedFloat(18, 38);
}
