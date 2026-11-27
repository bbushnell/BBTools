package structures;

import parse.Parse;

/** Unsigned compact counts: mantissa times two to the exponent, with exact zero.
 * No sign, bias, implicit leading bit, infinity or NaN. The canonical code uses
 * the smallest exponent whose rounded mantissa fits; ties round to even.
 * Twelve bits use4exponent/8mantissa bits. Eighteen bits use6exponent/12mantissa
 * bits. Decoding requires a signed long; a product above Long.MAX_VALUE is an
 * invalid resource code and fails before shifting. This is not an IEEE format.
 * Instances are immutable and allocate nothing while encoding/decoding.
 * @author Brian Bushnell, Raiden
 */
public final class CountFloat{

	private CountFloat(final int containerBits_, final int exponentBits_, final int mantissaBits_){
		containerBits=containerBits_;
		exponentBits=exponentBits_;
		mantissaBits=mantissaBits_;
		maximumMantissa=(1<<mantissaBits)-1;
		maximumExponent=(1<<exponentBits)-1;
		maximum=((long)maximumMantissa)<<Math.min(maximumExponent, 63-mantissaBits);
		payloadMask=(1<<(mantissaBits+exponentBits))-1;
		assert(containerBits%6==0 && mantissaBits+exponentBits<=containerBits && maximum>0)
			: "Count storage must fit its fixed-width native A48 container and valid signed-long count range";
	}

	/** Returns a supported unsigned-count format, never an IEEE equivalent. */
	public static CountFloat forBits(final int bits){
		if(bits==12){return TWELVE;}
		if(bits==18){return EIGHTEEN;}
		throw new IllegalArgumentException("CountFloat supports exactly the 12-bit and 18-bit count containers: "+bits);
	}

	/** Encodes a finite nonnegative quantity; out-of-range values fail, never wrap.
	 * Integers0..maximumMantissa+1are exact. Above that range, nearest rounding
	 * has relative error at most1/maximumMantissa. Fractional quantities have
	 * absolute error at most0.5while exponent0; fixed-point scaling is external.
	 */
	public int encode(final double value){
		if(!Double.isFinite(value) || value<0 || value>maximum){
			throw new IllegalArgumentException("CountFloat quantity must be finite, nonnegative and at most "+maximum+": "+value);
		}
		int exponent=Math.max(0, Math.getExponent(value)-mantissaBits+1);
		long mantissa=(long)Math.rint(Math.scalb(value, -exponent));
		if(mantissa>maximumMantissa){
			exponent++;
			mantissa=(long)Math.rint(Math.scalb(value, -exponent));
		}
		assert(exponent<=maximumExponent && mantissa>=0 && mantissa<=maximumMantissa
			&& (exponent==0 || mantissa>maximumMantissa/2))
			: "Nearest rounding must yield the unique smallest fitting exponent and a bounded unsigned mantissa";
		return (exponent<<mantissaBits)|(int)mantissa;
	}

	/** Returns the exact integer represented by a canonical code. */
	public long decode(final int code){
		if(code<0 || (code&~payloadMask)!=0){
			throw new IllegalArgumentException("CountFloat reserved/container bits must be zero: "+code);
		}
		final int exponent=code>>>mantissaBits, mantissa=code&maximumMantissa;
		if(exponent>0 && mantissa<=maximumMantissa/2){
			throw new IllegalArgumentException("Noncanonical CountFloat code has a smaller equivalent exponent: "+code);
		}
		final boolean fits=mantissa==0 || exponent<Long.numberOfLeadingZeros(mantissa);
		assert(fits) : "Stored count mantissa times two to its exponent must fit a signed long; code="+code;
		if(!fits){
			throw new IllegalArgumentException("Stored count exceeds a signed long: "+code);
		}
		return ((long)mantissa)<<exponent;
	}

	/** Appends a code in exactly two or three existing A48 symbols, zero padded. */
	public ByteBuilder append(final ByteBuilder output, final int code){
		decode(code);
		assert(output!=null) : "A count record appends into its caller-owned reusable output buffer";
		final int digits=Math.max(1, (32-Integer.numberOfLeadingZeros(code)+5)/6);
		for(int i=digits; i<containerBits/6; i++){output.append('0');}
		return output.appendA48(code);
	}

	/** Reads one fixed-width field through the native six-bit A48 parser. */
	public int read(final byte[] row, final int from){
		final int width=containerBits/6;
		if(row==null || from<0 || from>row.length-width){
			throw new IllegalArgumentException("CountFloat input must contain its complete fixed-width A48 field");
		}
		for(int i=from; i<from+width; i++){
			if(row[i]<48 || row[i]>111){throw new IllegalArgumentException("Invalid native A48 count symbol");}
		}
		final int code=(int)Parse.parseLongA48(row, from, from+width);
		decode(code);
		return code;
	}

	public final int containerBits, exponentBits, mantissaBits, maximumMantissa, maximumExponent;
	/** Largest valid count; contiguous exact integers stop at maximumMantissa+1. */
	public final long maximum;
	private final int payloadMask;
	private static final CountFloat TWELVE=new CountFloat(12, 4, 8);
	private static final CountFloat EIGHTEEN=new CountFloat(18, 6, 12);
}
