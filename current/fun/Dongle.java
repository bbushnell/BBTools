package fun;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

import shared.Tools;

/**
 * Utility class for time-based encoding, decoding, and validation of long
 * integer values using bitwise manipulation and rotation techniques.
 * Provides specialized encoding and decoding mechanisms for timestamps with
 * built-in validation and date parsing capabilities.
 *
 * @author Brian Bushnell
 */
public class Dongle {
	
	/**
	 * Main entry point that demonstrates encoding and decoding functionality.
	 * Parses the first argument as a timestamp (numeric or date string),
	 * encodes it, prints the encoded value, original timestamp, and decoded result.
	 * @param args Command line arguments; first argument is parsed as timestamp
	 */
	public static void main(String[] args) {
		long millis=System.currentTimeMillis();
		if(args.length>0) {
			millis=parse(args[0]);
		}
		long a=encode(millis);
		System.err.println(a);
		System.err.println(millis);
		System.err.println(decode(a));
	}
	
	/**
	 * Validates if the current system time falls within a specified time range.
	 * Uses decoded values from min/max bounds or provided arguments to establish
	 * the valid time window. Checks that the range is valid and within the limit.
	 * @param args Variable arguments: [1]=start time (Long), [2]=end time (Long)
	 * @return true if current time is within the valid decoded time range
	 */
	public static final boolean check(Object...args) {
		long a=min, b=max;
		if(args!=null && args.length>1) {a=(Long)args[1];}
		if(args!=null && args.length>2) {b=(Long)args[2];}
		a=decode(a);
		b=decode(b);
		long x=System.currentTimeMillis();
		return (b>a && b-a<limit && x>=a && x<=b);
	}
	
	/**
	 * Encodes a timestamp using complex bitwise XOR and bit rotation operations.
	 * Applies three stages of bit manipulation with different masks and rotations
	 * to obfuscate the original timestamp value.
	 * @param x The timestamp to encode
	 * @return The encoded timestamp
	 */
	private static long encode(long x) {
		x^=number;
		long a=x&mask;
		long b=x&(mask<<1);
		x=a|(Long.rotateLeft(b, rot));
		a=x&mask2;
		b=x&(mask2<<2);
		x=a|Long.rotateRight(b, rot2);
		a=x&mask3;
		b=x&(mask3>>>4);
		x=a|Long.rotateLeft(b, rot3);
		return x;
	}
	
	/**
	 * Decodes a previously encoded timestamp by reversing the encoding process.
	 * Applies the inverse operations of encode() in reverse order to recover
	 * the original timestamp value.
	 * @param x The encoded timestamp to decode
	 * @return The original decoded timestamp
	 */
	private static long decode(long x) {
		long a=x&mask3;
		long b=x&(mask3>>>4);
		x=a|(Long.rotateRight(b, rot3));
		a=x&mask2;
		b=x&(mask2<<2);
		x=(a|(Long.rotateLeft(b, rot2)));
		a=x&mask;
		b=x&(mask<<1);
		x=(a|(Long.rotateRight(b, rot)))^number;
		return x;
	}
	
	/**
	 * Parses a string as either a numeric timestamp or a date in MM-dd-yyyy format.
	 * If the string is numeric, returns it as a long. Otherwise, attempts to
	 * parse it as a date using the configured SimpleDateFormat.
	 * @param s The string to parse as a timestamp or date
	 * @return The parsed timestamp in milliseconds, or -1L if parsing fails
	 */
	private static long parse(String s) {
		try {
			if(Tools.isNumeric(s)) {return Long.parseLong(s);}
			Date d=sdf.parse(s);
			System.err.println(d);
			return d.getTime();
		} catch (ParseException e) {
			return -1L;
		}
	}
	
	/** Date format pattern used for parsing date strings */
	private static final String pattern="MM-dd-yyyy";
	/** SimpleDateFormat instance configured with the MM-dd-yyyy pattern */
	private static final SimpleDateFormat sdf=new SimpleDateFormat(pattern);
	/** Minimum encoded timestamp value for range validation */
	private static final long min=-8950436867421912347L;
	/** Maximum encoded timestamp value for range validation */
	private static final long max=2408489082714442483L;
	/** XOR key used in the encoding and decoding process */
	private static final long number=4964420948893066024L;
	/** Bit mask for first stage of encoding/decoding (alternating bits pattern) */
	private static final long mask=0x5555555555555555L;
	/** Bit mask for second stage of encoding/decoding (2-bit groups pattern) */
	private static final long mask2=0x3333333333333333L;
	/** Bit mask for third stage of encoding/decoding (4-bit groups pattern) */
	private static final long mask3=0xF0F0F0F0F0F0F0F0L;
	/** Maximum allowable time range in milliseconds for validation */
	private static final long limit=346896001029L;
	/** Rotation count for first stage bit rotation operations */
	private static final int rot=26;
	/** Rotation count for second stage bit rotation operations */
	private static final int rot2=44;
	/** Rotation count for third stage bit rotation operations */
	private static final int rot3=16;
	
}
