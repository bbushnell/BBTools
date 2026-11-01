package parse;

/**
 * Simple parser for delimiter-separated data using two-pointer scanning.
 * Maintains current position pointers to efficiently traverse byte arrays
 * separated by a specified delimiter character.
 *
 * <p>Minimal length-only scanner with borrowed bytes and unchecked cursor state.
 * Only the two-argument constructor supplies a usable line; there is no setter.
 * The delimiter-only constructor leaves line null. reset() retains the line.
 * advanceInner() requires a terminating delimiter inside the array; advance()
 * allows an unterminated final field. Neither method checks exhaustion before
 * incrementing its cursor. Instances are not thread-safe.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class LineParserSimple{

	/**
	 * Creates a LineParserSimple with the specified delimiter.
	 * Leaves line null; this class has no setter, so parsing requires the two-argument constructor.
	 * @param delimiter_ Byte value used as field separator
	 */
	public LineParserSimple(byte delimiter_){
		delimiter=delimiter_;
	}

	/**
	 * Creates a LineParserSimple with delimiter and line data.
	 * Ready for immediate parsing after construction.
	 * @param delimiter_ Byte value used as field separator
	 * @param line_ Byte array containing the line to parse
	 */
	public LineParserSimple(byte delimiter_, byte[] line_){
		delimiter=delimiter_;
		line=line_;
	}

	/**
	 * Advances to the next delimiter-separated segment without bounds checking.
	 * Updates internal pointers a and b to span the next field.
	 * Assumes caller has verified array bounds to avoid exceptions.
	 * @return Length of the current segment in bytes
	 */
	public int advanceInner(){//Does not check array bounds
		b++;
		a=b;
		assert(b<line.length);
		while(line[b]!=delimiter){b++;}
		return b-a;
	}

	/**
	 * Advances to the next delimiter-separated segment with bounds checking.
	 * Updates internal pointers a and b to span the next field.
	 * Stops scanning at array end; caller must still avoid advancing an exhausted cursor.
	 * @return Length of the current segment in bytes
	 */
	public int advance(){
		b++;
		a=b;
		while(b<line.length && line[b]!=delimiter){b++;}
		return b-a;
	}

	/** Resets parser state to initial position.
	 * Sets all position pointers to -1 for fresh parsing. */
	public void reset(){
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
