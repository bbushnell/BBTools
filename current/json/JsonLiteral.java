package json;

import shared.Tools;

/**
 * Immutable unquoted text for JSON output, including explicitly formatted numbers.
 * The string constructor is a trusted escape hatch: it does not validate JSON,
 * add quotes, or escape non-ASCII text. Supply an already escaped ASCII fragment.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class JsonLiteral{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	/** Borrows a preformatted literal without validation or escaping.
	 * @param s_ ASCII JSON fragment; null is rendered as JSON null by JsonObject */
	public JsonLiteral(String s_){
		s=s_;
	}
	
	/**
	 * Creates a JsonLiteral from a double value with specified decimal precision.
	 * Formats the numeric value using Tools.format with the specified number of
	 * decimal places, using the locale-independent ROOT formatter.
	 * NaN and either infinity are represented by JSON null.
	 *
	 * @param value The numeric value to format
	 * @param decimals Nonnegative number of decimal places
	 * @throws IllegalArgumentException If the requested format is invalid
	 */
	public JsonLiteral(double value, int decimals){
		if(decimals<0){throw new IllegalArgumentException("Decimal places must be nonnegative: "+decimals);}
		s=Double.isFinite(value) ? Tools.format("%."+decimals+"f", value) : "null";
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	/** Returns the stored fragment verbatim, including null when supplied by the caller.
	 * @return Stored literal text */
	@Override
	public String toString(){return s;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	/** The formatted string representation of this literal value */
	private final String s;
	
}
