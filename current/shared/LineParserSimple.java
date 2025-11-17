package shared;

/**
 * Simple parser for delimiter-separated data using two-pointer scanning.
 * Maintains current position pointers to efficiently traverse byte arrays
 * separated by a specified delimiter character.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class LineParserSimple {

	/**
	 * Creates a LineParserSimple with the specified delimiter.
	 * Line data must be set separately before parsing.
	 * @param delimiter_ Byte value used as field separator
	 */
	public LineParserSimple(byte delimiter_) {
		delimiter=delimiter_;
	}

	/**
	 * Creates a LineParserSimple with delimiter and line data.
	 * Ready for immediate parsing after construction.
	 * @param delimiter_ Byte value used as field separator
	 * @param line_ Byte array containing the line to parse
	 */
	public LineParserSimple(byte delimiter_, byte[] line_) {
		delimiter=delimiter_;
		line=line_;
	}
	
	/**
	 * Advances to the next delimiter-separated segment without bounds checking.
	 * Updates internal pointers a and b to span the next field.
	 * Assumes caller has verified array bounds to avoid exceptions.
	 * @return Length of the current segment in bytes
	 */
	public int advanceInner() {//Does not check array bounds
		b++;
		a=b;
		assert(b<line.length);
		while(line[b]!=delimiter){b++;}
		return b-a;
	}
	
	/**
	 * Advances to the next delimiter-separated segment with bounds checking.
	 * Updates internal pointers a and b to span the next field.
	 * Safe version that checks array bounds before accessing elements.
	 * @return Length of the current segment in bytes
	 */
	public int advance() {
		b++;
		a=b;
		while(b<line.length && line[b]!=delimiter){b++;}
		return b-a;
	}
	
	/** Resets parser state to initial position.
	 * Sets all position pointers to -1 for fresh parsing. */
	public void reset() {
		a=b=segment=-1;
	}
	
	/** Start position of current segment being parsed */
	private int a=-1;
	/** End position of current segment being parsed */
	private int b=-1;
	/** Current segment number being processed */
	private int segment=-1;
	/** Byte array containing the line data to parse */
	private byte[] line;
	
	/** Delimiter character used to separate fields */
	private final byte delimiter;
}
