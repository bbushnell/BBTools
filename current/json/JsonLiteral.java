package json;

import shared.Tools;

/**
 * Represents unquoted literal values in JSON output for precise numeric formatting.
 * Enables creation of JSON-compatible numeric literals that can be precisely formatted
 * without string quotation, supporting both pre-formatted strings and numeric values
 * with specific decimal precision.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class JsonLiteral {
	
	/** Creates a JsonLiteral from a pre-formatted string.
	 * @param s_ The string literal value to store */
	public JsonLiteral(String s_){
		s=s_;
	}
	
	/**
	 * Creates a JsonLiteral from a double value with specified decimal precision.
	 * Formats the numeric value using Tools.format with the specified number of
	 * decimal places.
	 *
	 * @param value The numeric value to format
	 * @param decimals Number of decimal places to include in formatting
	 */
	public JsonLiteral(double value, int decimals){
		s=Tools.format("%."+decimals+"f", value);
	}
	
	@Override
	public String toString(){return s;}
	
	/** The formatted string representation of this literal value */
	private final String s;
	
}
