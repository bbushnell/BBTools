package parse;

import java.util.ArrayList;

import shared.Tools;
import structures.ByteBuilder;

/** 
 * Mutable field access for delimited lines, using byte indexes in byte parsers
 * and UTF-16 indexes in {@link LineParserS} implementations. Field indexes are
 * zero-based; bounds are half-open. Instances are not thread-safe and retain
 * their input without copying it. Do not mutate borrowed byte arrays while parsing.
 * Empty fields are retained; the reverse string variant returns null for an
 * empty field where other variants return an empty string.
 * Indexed variants retain delimiter positions for random access; variants 2
 * scan forward with bounded auxiliary memory. Their indexed calls must not
 * move backwards without reset. Numeric conversion delegates to {@link Parse}.
 * Callers must provide an initialized line, existing field indexes and offsets
 * inside the selected field. Assertions are diagnostics, not complete validation
 * of malformed indexes or numeric input when assertions are disabled.
 * 
 * @author Brian Bushnell
 * @date May 24, 2023
 *
 */
public interface LineParser{

	/**
	 * Designate the line of text to be parsed.
	 * @param line Line to parse.
	 * @return this
	 */
	public LineParser set(byte[] line);

	/**
	 * Designate the line of text to be parsed.
	 * @param line Line to parse.
	 * @param maxTerm Inclusive last term to index; it must exist for indexed variants.
	 * Forward-only variants ignore the limit; the reverse variant rejects this overload.
	 * @return this
	 */
	public LineParser set(byte[] line, int maxTerm);

	/**
	 * Clear all internal state.
	 * @return this
	 */
	public LineParser clear();

	/**
	 * Reset forward-only traversal to before the first field. Indexed variants
	 * deliberately do nothing, including retaining their current bounds.
	 * @return this
	 */
	public LineParser reset();

	/**
	 * Parse the designated field as an int.
	 * @param term Field to parse.
	 * @return int value.
	 */
	public int parseInt(int term);

	/**
	 * Parse the current field as an int.
	 * @return int value.
	 */
	public int parseIntFromCurrentField();

	/**
	 * Parse the designated field as a long.
	 * @param term Field to parse.
	 * @return long value.
	 */
	public long parseLong(int term);
	
	/**
	 * Parse the designated field as a float.
	 * @param term Field to parse.
	 * @return float value.
	 */
	public float parseFloat(int term);
	
	/**
	 * Parse the designated field as a double.
	 * @param term Field to parse.
	 * @return double value.
	 */
	public double parseDouble(int term);
	
	/**
	 * Parse a byte from the designated field.
	 * @param term Field to parse.
	 * @param offset Position in the field to read.
	 * @return byte value.
	 */
	public byte parseByte(int term, int offset);
	
	/**
	 * Parse the designated field as a byte[].
	 * @param term Field to parse.
	 * @return A byte[].
	 */
	public byte[] parseByteArray(int term);
	
	/**
	 * Parse the designated field as a String.
	 * @param term Field to parse.
	 * @return A String.
	 */
	public String parseString(int term);
	
	/**
	 * USE CAUTION WITH THIS METHOD.
	 * Parse the current field as a String.
	 * @return A String.
	 */
	public String parseStringFromCurrentField();
	
	/**
	 * USE CAUTION WITH THIS METHOD.
	 * Parse the current field as a byte[].
	 * @return A byte[].
	 */
	public byte[] parseByteArrayFromCurrentField();

	/** Append the raw term to bb without parsing, and without delimiters */
	public ByteBuilder appendTerm(ByteBuilder bb, int term);
	
	/** Convenience method; determines whether the line starts with this term. 
	 * Implemented in Tools. */
	public boolean startsWith(String s);
	
	/** Convenience method; determines whether the line starts with this character. 
	 * Implemented in Tools. */
	public boolean startsWith(char c);
	
	/** Convenience method; determines whether the line starts with this character. 
	 * Implemented in Tools. */
	public boolean startsWith(byte c);

	/** True if the term starts with s */
	public boolean termStartsWith(String s, int term);
	
	/** Tests the first byte of a nonempty field with Tools.isLetter. */
	public default boolean termStartsWithLetter(int term){
		return Tools.isLetter(parseByte(term, 0));
	}
	
	/** True if the term equals s */
	public boolean termEquals(String s, int term);
	
	/** True if the term equals c */
	public boolean termEquals(char c, int term);
	
	/** True if the term equals c */
	public boolean termEquals(byte c, int term);
	
	/** Length of this field. Indexed variants select it as the current field;
	 * forward-only variants restore their prior position after measuring. */
	public int length(int term);

	/** Length of current field. */
	public int currentFieldLength();
	
	/** Adjust the current start bound without updating field indexes or checking bounds.
	 * The caller owns the resulting range; returns its new length. */
	public int incrementA(int amt);

	/** Adjust the current exclusive end bound without updating indexes or checking bounds.
	 * The caller owns the resulting range; returns its new length. */
	public int incrementB(int amt);

	/** Variant-specific cursor test, not a total-field count. Forward-only variants
	 * use this before advancing; indexed variants reflect their most recent bound selection. */
	public boolean hasMore();

	/** Length of the line being parsed */
	public int lineLength();
	
	/** Select a field and return its length. Forward-only variants cannot move backwards. */
	int setBounds(int term);
	
	/** Line being parsed */
	public Object line();
	
	/** Position of left pointer (first character of current term) */
	public int a();

	/** Position of right pointer (just after end of current term) */
	public int b();
	
	@Override
	public String toString();
	
	/**
	 * Materialize fields as strings. Indexed variants include all indexed fields;
	 * forward-only variants consume remaining fields and require an unexhausted cursor.
	 * @return A List.
	 */
	public ArrayList<String> toList();

}
