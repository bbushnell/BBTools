package shared;

/** 
 * Interface for classes that parse delimited string lines.
 * 
 * @author Brian Bushnell
 * @date April 3, 2024
 *
 */
public interface LineParserS extends LineParser {

	/**
	 * Designate the line of text to be parsed.
	 * @param line Line to parse.
	 * @return this
	 */
	public LineParserS set(String line);

	/**
	 * Designate the line of text to be parsed.
	 * @param line Line to parse.
	 * @param maxTerm Stop parsing after this term.
	 * @return this
	 */
	public LineParserS set(String line, int maxTerm);

	/**
	 * Clear all internal state.
	 * @return this
	 */
	public LineParserS clear();

	/**
	 * Reset state to the start of the current line.
	 * @return this
	 */
	public LineParserS reset();
	
	/**
	 * Parse a char from the designated field.
	 * @param term Field to parse.
	 * @param offset Position in the field to read.
	 * @return char value.
	 */
	public char parseChar(int term, int offset);
	

}