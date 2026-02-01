package parse;

/**
 * Interface for parsing delimited string lines with convenience setters and helpers.
 * Extends LineParser to add fluent setters and character-level access for parsing.
 * @author Brian Bushnell
 * @date April 3, 2024
 */
public interface LineParserS extends LineParser {

	/**
	 * Designates the line of text to be parsed and resets internal state.
	 * @param line Line to parse
	 * @return this parser for method chaining
	 */
	public LineParserS set(String line);

	/**
	 * Designates the line of text to be parsed with a maximum term limit.
	 * Parsing will stop after the specified term to improve performance.
	 * @param line Line to parse
	 * @param maxTerm Stop parsing after this term number
	 * @return this parser for method chaining
	 */
	public LineParserS set(String line, int maxTerm);

	/** Clears all internal parser state and resets to initial conditions.
	 * @return this parser for method chaining */
	public LineParserS clear();

	/**
	 * Resets parser state to the beginning of the current line without clearing the line data.
	 * Allows re-parsing the same line from the start.
	 * @return this parser for method chaining
	 */
	public LineParserS reset();
	
	/**
	 * Extracts a single character from a specified field at a given position.
	 * @param term Field number to parse (0-based)
	 * @param offset Position within the field to read (0-based)
	 * @return Character at the specified position
	 */
	public char parseChar(int term, int offset);
	

}