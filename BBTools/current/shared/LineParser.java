package shared;

import java.util.ArrayList;

import structures.ByteBuilder;

/** 
 * Interface for classes that parse delimited byte[] lines.
 * 
 * @author Brian Bushnell
 * @date May 24, 2023
 *
 */
public interface LineParser {

	/**
	 * Designate the line of text to be parsed.
	 * @param line Line to parse.
	 * @return this
	 */
	public LineParser set(byte[] line);

	/**
	 * Designate the line of text to be parsed.
	 * @param line Line to parse.
	 * @param maxTerm Stop parsing after this term.
	 * @return this
	 */
	public LineParser set(byte[] line, int maxTerm);

	/**
	 * Clear all internal state.
	 * @return this
	 */
	public LineParser clear();

	/**
	 * Reset state to the start of the current line.
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
	
	/** True if the term equals s */
	public boolean termEquals(String s, int term);
	
	/** True if the term equals c */
	public boolean termEquals(char c, int term);
	
	/** True if the term equals c */
	public boolean termEquals(byte c, int term);
	
	/** Length of this field. */
	public int length(int term);

	/** Length of current field. */
	public int currentFieldLength();
	
	/** Not recommended.  Secret method. */
	public int incrementA(int amt);

	/** Not recommended. Secret method. */
	public int incrementB(int amt);

	/** True if there is a subsequent term */
	public boolean hasMore();

	/** Length of the line being parsed */
	public int lineLength();
	
	/** LP1 and LP2 both can use this, but LP2 will fail if it's advancing backwards. */
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
	 * Format the line as a list of Strings, one per field.
	 * @return A List.
	 */
	public ArrayList<String> toList();

}