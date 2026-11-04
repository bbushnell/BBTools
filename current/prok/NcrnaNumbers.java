package prok;

import parse.Parse;

/** Numeric parsing shared by ncRNA readers and startup configuration.
 * Unsigned range methods parse borrowed bytes/characters without allocation;
 * they accept only ASCII digits in 0..Integer.MAX_VALUE. Format readers attach
 * their own file/field context. CLI methods retain native numeric syntax while
 * adding option names to errors and rejecting absent values before parsing.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaNumbers{

	/** Parses a native float option, retaining its name in missing/malformed errors. */
	public static float floatValue(final String name, final String text){
		value(name, text);
		try{return Float.parseFloat(text);}catch(NumberFormatException failure){throw invalid(name, text, failure);}
	}

	/** Parses a double option before any deliberate narrowing for native scoring. */
	public static double doubleValue(final String name, final String text){
		value(name, text);
		try{return Double.parseDouble(text);}catch(NumberFormatException failure){throw invalid(name, text, failure);}
	}

	/** Parses a signed decimal integer option; range errors include the flag name. */
	public static int intValue(final String name, final String text){
		value(name, text);
		try{return Integer.parseInt(text);}catch(NumberFormatException failure){throw invalid(name, text, failure);}
	}

	/** Parses a signed decimal long option, such as a random seed. */
	public static long longValue(final String name, final String text){
		value(name, text);
		try{return Long.parseLong(text);}catch(NumberFormatException failure){throw invalid(name, text, failure);}
	}

	/** Uses BBTools decimal KMG suffixes and checks int range before narrowing. */
	public static int intKMG(final String name, final String text){
		value(name, text);
		final long parsed;
		try{parsed=Parse.parseKMG(text);}catch(RuntimeException failure){throw invalid(name, text, failure);}
		if(parsed<Integer.MIN_VALUE || parsed>Integer.MAX_VALUE){throw invalid(name, text, null);}
		return (int)parsed;
	}

	/** Uses BBTools binary KMG units for an explicitly supplied memory budget. */
	public static long binaryValue(final String name, final String text){
		value(name, text);
		try{return Parse.parseKMGBinary(text);}catch(RuntimeException failure){throw invalid(name, text, failure);}
	}

	/** Native numeric parsers disagree on null; numeric options always require values. */
	public static String value(final String name, final String text){
		if(text==null || text.isEmpty()){throw invalid(name, text, null);}
		return text;
	}

	private static IllegalArgumentException invalid(final String name, final String text, final RuntimeException cause){
		return new IllegalArgumentException("Invalid numeric option "+name+"="+text+"; supply a valid value in range", cause);
	}

	/** Parses the half-open byte range [from,to) as a nonnegative signed-range int. */
	public static int unsignedInt(final byte[] text, final int from, final int to){
		if(text==null || from<0 || from>=to || to>text.length){
			throw new IllegalArgumentException("Unsigned decimal requires a nonempty valid byte range: "+from+".."+to);
		}
		int value=0;
		for(int i=from; i<to; i++){value=appendDigit(value, text[i]-'0');}
		return value;
	}

	/** Parses the half-open character range [from,to), without allocating a substring. */
	public static int unsignedInt(final CharSequence text, final int from, final int to){
		if(text==null || from<0 || from>=to || to>text.length()){
			throw new IllegalArgumentException("Unsigned decimal requires a nonempty valid character range: "+from+".."+to);
		}
		int value=0;
		for(int i=from; i<to; i++){value=appendDigit(value, text.charAt(i)-'0');}
		return value;
	}

	/** Checks each appended digit before narrowing; prior values never exceed int range. */
	private static int appendDigit(final int value, final int digit){
		assert(value>=0) : "Every previous decimal prefix was checked before narrowing to int";
		if(digit<0 || digit>9){throw new IllegalArgumentException("Unsigned decimal requires ASCII digits");}
		final long next=(long)value*10+digit;
		if(next>Integer.MAX_VALUE){throw new IllegalArgumentException("Unsigned decimal exceeds Integer.MAX_VALUE");}
		return (int)next;
	}

	private NcrnaNumbers(){}
}
