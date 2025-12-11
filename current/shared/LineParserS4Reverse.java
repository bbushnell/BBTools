package shared;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextFile;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Parses lines right-to-left using multiple ordered delimiters.
 * Processes lines from end to start while maintaining left-to-right delimiter interpretation.
 * Designed for situations where the end structure is known but the prefix is uncertain.
 *
 * @author Brian Bushnell
 * @date May 6, 2024
 */
public final class LineParserS4Reverse implements LineParserS {
	
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
		
		LineParserS4Reverse lp=new LineParserS4Reverse(dstring);
		for(String line : lines) {
			lp.set(line);
			System.out.println(lp.bounds);
			System.out.println(lp);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public LineParserS4Reverse(String delimiters_) {
		delimiters=new String(Tools.reverseAndCopy(delimiters_.toCharArray()));
		maxDPos=delimiters.length()-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public LineParserS4Reverse set(byte[] line_) {
		assert(false) : "Use byte version.";
		return set(new String(line_));
	}

	@Override
	public LineParserS4Reverse set(byte[] line_, int maxTerm) {
		assert(false) : "Use byte version.";
		return set(new String(line_), maxTerm);
	}
	
	/**
	 * Sets the line to parse and computes term boundaries.
	 * Parses from right to left, building bounds array in reverse order.
	 * @param line_ String containing the line to parse
	 * @return This parser instance for chaining
	 */
	@Override
	public LineParserS4Reverse set(String line_) {
		clear();
		line=line_;
		a=b=line.length();
		bounds.add(a);
		for(int len=advance(); a>=0; len=advance()) {
			bounds.add(a);
		}
		while(bounds.size()<=delimiters.length()) {bounds.add(bounds.lastElement());}
		bounds.reverse();
		return this;
	}
	
	/**
	 * Sets the line to parse with maximum term limit.
	 * Not implemented for this reverse parser.
	 *
	 * @param line_ String containing the line to parse
	 * @param maxTerm Maximum number of terms to parse
	 * @return This parser instance for chaining
	 * @throws RuntimeException Always thrown as this method is not valid for reverse parsing
	 */
	@Override
	public LineParserS4Reverse set(String line_, int maxTerm) {
		throw new RuntimeException("Not valid.");
	}
	
	/**
	 * Clears the parser state and resets all fields.
	 * Prepares the parser for a new line.
	 * @return This parser instance for chaining
	 */
	@Override
	public LineParserS4Reverse clear() {
		delimiterPos=0;
		line=null;
		a=b=-1;
		bounds.clear();
		return this;
	}
	
	/**
	 * Resets the parser to initial parsing position.
	 * Does nothing for this reverse parser implementation.
	 * @return This parser instance for chaining
	 */
	@Override
	public LineParserS4Reverse reset() {
		//Does nothing for this class
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public int terms() {return bounds.size();}
	
	/**
	 * Parses the specified term as an integer.
	 * @param term Zero-based term index to parse
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
	 * @param term Zero-based term index to parse
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
	 * @param term Zero-based term index to parse
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
	 * Parses a character from the specified term at given offset.
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
	 * @param term Zero-based term index to parse
	 * @return Byte array representation of the term
	 */
	@Override
	public byte[] parseByteArray(int term) {
		final int len=setBounds(term);
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}
	
	/** Returns the current field as a byte array */
	@Override
	public byte[] parseByteArrayFromCurrentField() {
		int len=b-a;
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}
	
	/**
	 * Parses the specified term as a string.
	 * @param term Zero-based term index to parse
	 * @return String value of the term, or null if empty
	 */
	@Override
	public String parseString(int term) {
		final int len=setBounds(term);
		return a>=b ? null : line.substring(a, b);
	}

	/**
	 * Appends the specified term to a ByteBuilder.
	 * @param bb ByteBuilder to append to
	 * @param term Zero-based term index to append
	 * @return The modified ByteBuilder for chaining
	 */
	@Override
	public ByteBuilder appendTerm(ByteBuilder bb, int term) {
		final int len=setBounds(term);
		for(int i=a; i<b; i++) {bb.append(line.charAt(i));}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Parses the current field as an integer */
	@Override
	public int parseIntFromCurrentField() {
		return Parse.parseInt(line, a, b);
	}
	
	/** Parses the current field as a string */
	@Override
	public String parseStringFromCurrentField() {
		return line.substring(a, b);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Query Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Checks if the line starts with the specified string.
	 * @param s String to check for
	 * @return true if line starts with the string
	 */
	@Override
	public boolean startsWith(String s) {
		return line.startsWith(s);
	}
	
	/**
	 * Checks if the line starts with the specified character.
	 * @param c Character to check for
	 * @return true if line starts with the character
	 */
	@Override
	public boolean startsWith(char c) {
		return Tools.startsWith(line, c);
	}
	
	/**
	 * Checks if the line starts with the specified byte.
	 * @param b Byte to check for
	 * @return true if line starts with the byte
	 */
	@Override
	public boolean startsWith(byte b) {
		return Tools.startsWith(line, b);
	}
	
	/**
	 * Checks if the specified term starts with the given string.
	 * @param s String to check for
	 * @param term Zero-based term index
	 * @return true if the term starts with the string
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
	 * Checks if the specified term equals the given string.
	 * @param s String to compare against
	 * @param term Zero-based term index
	 * @return true if the term equals the string
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
	 * Checks if the specified term equals the given character.
	 * @param c Character to compare against
	 * @param term Zero-based term index
	 * @return true if the term equals the character
	 */
	@Override
	public boolean termEquals(char c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	/**
	 * Checks if the specified term equals the given byte.
	 * @param c Byte to compare against
	 * @param term Zero-based term index
	 * @return true if the term equals the byte
	 */
	@Override
	public boolean termEquals(byte c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	/**
	 * Increments the start boundary of the current field.
	 * @param amt Amount to increment by
	 * @return New length of the current field
	 */
	@Override
	public int incrementA(int amt) {
		a+=amt;
		return b-a;
	}
	
	/**
	 * Increments the start boundary of the current field.
	 * Note: This appears to be a bug as it increments 'a' instead of 'b'.
	 * @param amt Amount to increment by
	 * @return New length of the current field
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

	/** Returns the length of the current field */
	@Override
	public int currentFieldLength() {
		return b-a;
	}

	/** Returns true if there are more characters to parse */
	@Override
	public boolean hasMore() {
		return a>=0;
	}

	/** Returns the length of the entire line */
	@Override
	public int lineLength() {
		return line.length();
	}

	/** Returns the current line being parsed */
	@Override
	public String line() {return line;}
	
	/** Returns the start position of the current field */
	@Override
	public int a() {return a;}
	
	/** Returns the end position of the current field */
	@Override
	public int b() {return b;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets the parsing boundaries for the specified term.
	 * Updates internal position markers a and b.
	 * @param term Zero-based term index
	 * @return Length of the term
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
		a--;
		b=a;
		while(a>=0 && delimiter!=line.charAt(a)){a--;}
//		System.err.println("delimiter="+delimiter+", a="+a+", b="+b+", len="+(b-a));
		return b-a;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns string representation of all parsed terms */
	@Override
	public String toString() {
		return toList().toString();
	}
	
	/** Converts all parsed terms to a list of strings.
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
	
	//This is already reversed
	public final String delimiters;
	private final int maxDPos;
	private int delimiterPos=0;
	
}
