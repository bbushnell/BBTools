package parse;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextFile;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Line parser that uses multiple ordered delimiters for field separation.
 * Each position in the parsing sequence uses a different delimiter character,
 * allowing flexible parsing of structured data with varying separators.
 *
 * @author Brian Bushnell
 * @date May 24, 2023
 */
public final class LineParserS4 implements LineParserS {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	//For testing
	//Syntax: LineParser fname/literal delimiter 
	public static void main(String[] args) {
		assert(args.length==2);
		String fname=args[0];
		String dstring=args[1];
//		assert(dstring.length()==1);
		
		final String[] lines;
		if(new File(fname).exists()){
			lines=TextFile.toStringLines(fname);
		}else{
			lines=new String[] {fname};
		}
		
		LineParserS lp=new LineParserS4(dstring);
		for(String line : lines) {
			lp.set(line);
			System.out.println(lp);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public LineParserS4(String delimiters_) {
		delimiters=delimiters_;
		maxDPos=delimiters.length()-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public LineParserS4 set(byte[] line_) {
		assert(false) : "Use byte version.";
		return set(new String(line_));
	}

	@Override
	public LineParserS4 set(byte[] line_, int maxTerm) {
		assert(false) : "Use byte version.";
		return set(new String(line_), maxTerm);
	}
	
	/**
	 * Sets the line to parse and identifies all field boundaries.
	 * Uses the ordered delimiter sequence to separate fields, advancing through
	 * each delimiter position until the line is fully parsed.
	 *
	 * @param line_ String containing the line to parse
	 * @return This parser instance for method chaining
	 */
	@Override
	public LineParserS4 set(String line_) {
		clear();
		line=line_;
		for(int len=advance(); b<line.length(); len=advance()) {
			bounds.add(b);
		}
		bounds.add(b);
		return this;
	}
	
	/**
	 * Sets the line to parse with a maximum number of terms.
	 * Parses only up to the specified number of terms using ordered delimiters.
	 *
	 * @param line_ String containing the line to parse
	 * @param maxTerm Maximum number of terms to parse
	 * @return This parser instance for method chaining
	 */
	@Override
	public LineParserS4 set(String line_, int maxTerm) {
		clear();
		line=line_;
		//TODO: test performance of presumably safer loop below
//		for(int term=0; term<=maxTerm && b<line.length(); term++) {
		for(int term=0; term<=maxTerm; term++) {
			int len=advance();
			bounds.add(b);
		}
		return this;
	}
	
	/**
	 * Clears the parser state and resets all internal variables.
	 * Resets delimiter position, line reference, bounds, and indices.
	 * @return This parser instance for method chaining
	 */
	@Override
	public LineParserS4 clear() {
		delimiterPos=0;
		line=null;
		a=b=-1;
		bounds.clear();
		return this;
	}
	
	/**
	 * Resets the parser for reuse.
	 * This implementation does nothing as the parser maintains no state to reset.
	 * @return This parser instance for method chaining
	 */
	@Override
	public LineParserS4 reset() {
		//Does nothing for this class
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public int terms() {return bounds.size();}
	
	/**
	 * Parses the specified term as an integer.
	 * @param term Zero-based term index
	 * @return Integer value of the term
	 */
	@Override
	public int parseInt(int term) {
		setBounds(term);
		return Parse.parseInt(line, a, b);
	}
	
	public int parseInt(int term, int from, int to) {
		setBounds(term);
		return Parse.parseInt(line, a+from, Tools.min(line.length(), a+to));
//		return Parse.parseInt(line, a+from, Tools.min(b, a+to));
	}
	
	/**
	 * Parses the specified term as a long integer.
	 * @param term Zero-based term index
	 * @return Long value of the term
	 */
	@Override
	public long parseLong(int term) {
		setBounds(term);
		return Parse.parseLong(line, a, b);
	}
	
	@Override
	public float parseFloat(int term) {
		setBounds(term);
		return Parse.parseFloat(line, a, b);
	}
	
	/**
	 * Parses the specified term as a double.
	 * @param term Zero-based term index
	 * @return Double value of the term
	 */
	@Override
	public double parseDouble(int term) {
		setBounds(term);
		return Parse.parseDouble(line, a, b);
	}
	
	@Override
	public byte parseByte(int term, int offset) {
		return (byte)parseChar(term, offset);
	}
	
	/**
	 * Parses a single character from the specified term at the given offset.
	 * @param term Zero-based term index
	 * @param offset Character offset within the term
	 * @return Character at the specified position
	 */
	@Override
	public char parseChar(int term, int offset) {
		setBounds(term);
		final int index=a+offset;
		assert(index<b);
		return line.charAt(index);
	}
	
	/**
	 * Parses the specified term as a byte array.
	 * Each character in the term is converted to its byte representation.
	 * @param term Zero-based term index
	 * @return Byte array representation of the term
	 */
	@Override
	public byte[] parseByteArray(int term) {
		final int len=setBounds(term);
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}
	
	/**
	 * Parses the current field as a byte array.
	 * Uses the currently set field boundaries (a and b indices).
	 * @return Byte array representation of the current field
	 */
	@Override
	public byte[] parseByteArrayFromCurrentField() {
		int len=b-a;
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}
	
	/**
	 * Parses the specified term as a String.
	 * @param term Zero-based term index
	 * @return String value of the term
	 */
	@Override
	public String parseString(int term) {
		final int len=setBounds(term);
		return line.substring(a, b);
	}

	/**
	 * Appends the specified term to a ByteBuilder.
	 * Each character in the term is appended individually.
	 *
	 * @param bb ByteBuilder to append to
	 * @param term Zero-based term index
	 * @return The ByteBuilder for method chaining
	 */
	@Override
	public ByteBuilder appendTerm(ByteBuilder bb, int term) {
		final int len=setBounds(term);
		for(int i=a; i<b; i++) {bb.append(line.charAt(i));}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses the current field as an integer.
	 * Uses the currently set field boundaries (a and b indices).
	 * @return Integer value of the current field
	 */
	@Override
	public int parseIntFromCurrentField() {
		return Parse.parseInt(line, a, b);
	}
	
	/**
	 * Parses the current field as a String.
	 * Uses the currently set field boundaries (a and b indices).
	 * @return String value of the current field
	 */
	@Override
	public String parseStringFromCurrentField() {
		return line.substring(a, b);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Query Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Tests if the current line starts with the specified string.
	 * @param s String to test for
	 * @return true if the line starts with the string, false otherwise
	 */
	@Override
	public boolean startsWith(String s) {
		return line.startsWith(s);
	}
	
	/**
	 * Tests if the current line starts with the specified character.
	 * @param c Character to test for
	 * @return true if the line starts with the character, false otherwise
	 */
	@Override
	public boolean startsWith(char c) {
		return Tools.startsWith(line, c);
	}
	
	/**
	 * Tests if the current line starts with the specified byte.
	 * @param b Byte to test for
	 * @return true if the line starts with the byte, false otherwise
	 */
	@Override
	public boolean startsWith(byte b) {
		return Tools.startsWith(line, b);
	}
	
	/**
	 * Tests if the specified term starts with the given string.
	 * @param s String to test for
	 * @param term Zero-based term index
	 * @return true if the term starts with the string, false otherwise
	 */
	@Override
	public boolean termStartsWith(String s, int term) {
		final int len=setBounds(term);
		if(len<s.length()) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(c!=line.charAt(a+i)) {return false;}
		}
		return true;
	}
	
	/**
	 * Tests if the specified term equals the given string.
	 * @param s String to compare against
	 * @param term Zero-based term index
	 * @return true if the term equals the string, false otherwise
	 */
	@Override
	public boolean termEquals(String s, int term) {
		final int len=setBounds(term);
		if(len!=s.length()) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(c!=line.charAt(a+i)) {return false;}
		}
		return true;
	}
	
	/**
	 * Tests if the specified term equals the given character.
	 * @param c Character to compare against
	 * @param term Zero-based term index
	 * @return true if the term equals the character, false otherwise
	 */
	@Override
	public boolean termEquals(char c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	/**
	 * Tests if the specified term equals the given byte.
	 * @param c Byte to compare against
	 * @param term Zero-based term index
	 * @return true if the term equals the byte, false otherwise
	 */
	@Override
	public boolean termEquals(byte c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	/**
	 * Increments the start boundary (a) by the specified amount.
	 * @param amt Amount to increment
	 * @return New length of the current field (b-a)
	 */
	@Override
	public int incrementA(int amt) {
		a+=amt;
		return b-a;
	}
	
	/**
	 * Increments the start boundary (a) by the specified amount.
	 * Note: This appears to increment 'a' rather than 'b' as the name suggests.
	 * @param amt Amount to increment
	 * @return New length of the current field (b-a)
	 */
	@Override
	public int incrementB(int amt) {
		a+=amt;
		return b-a;
	}

	/**
	 * Returns the length of the specified term.
	 * @param term Zero-based term index
	 * @return Length of the term in characters
	 */
	@Override
	public int length(int term) {
		return setBounds(term);
	}

	/** Returns the length of the current field (b-a) */
	@Override
	public int currentFieldLength() {
		return b-a;
	}

	/** Tests if there are more characters to parse in the line.
	 * @return true if the end boundary is before the line end, false otherwise */
	@Override
	public boolean hasMore() {
		return b<line.length();
	}

	/** Returns the total length of the current line */
	@Override
	public int lineLength() {
		return line.length();
	}

	/** Returns the current line being parsed */
	@Override
	public String line() {return line;}
	
	/** Returns the start boundary index of the current field */
	@Override
	public int a() {return a;}
	
	/** Returns the end boundary index of the current field */
	@Override
	public int b() {return b;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets the field boundaries (a and b) for the specified term.
	 * Updates internal state to point to the requested term's boundaries.
	 * @param term Zero-based term index
	 * @return Length of the term (b-a)
	 */
	@Override
	public int setBounds(int term){
		a=(term==0 ? 0 : bounds.get(term-1)+1);
		b=bounds.get(term);
		return b-a;
	}
	
	private int advance() {
		char delimiter=(delimiterPos<delimiters.length() ? delimiters.charAt(delimiterPos) : 0);
		delimiterPos++;
		b++;
		a=b;
		while(b<line.length() && delimiter!=line.charAt(b)){b++;}
		return b-a;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns a string representation of all parsed terms.
	 * @return String representation as a list of terms */
	@Override
	public String toString() {
		return toList().toString();
	}
	
	/** Converts all parsed terms to an ArrayList of strings.
	 * @return ArrayList containing all terms as strings */
	@Override
	public ArrayList<String> toList(){
		ArrayList<String> list=new ArrayList<String>(bounds.size);
		for(int i=0; i<bounds.size; i++){
			list.add(parseString(i));
		}
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final IntList bounds=new IntList();
	
	private int a=-1;
	private int b=-1;
	private String line;
	
	public final String delimiters;
	private final int maxDPos;
	private int delimiterPos=0;
	
}
