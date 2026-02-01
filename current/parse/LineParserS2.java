package parse;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextFile;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Memory-efficient line parser for delimited text strings.
 * Similar speed to other parsers but uses bounded memory footprint for very long lines.
 * Implements the LineParserS interface using string-based parsing with character delimiters.
 *
 * @author Brian Bushnell
 * @date May 24, 2023
 */
public final class LineParserS2 implements LineParserS {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	//For testing
	//Syntax: LineParser fname/literal delimiter 
	public static void main(String[] args) {
		assert(args.length==2);
		String fname=args[0];
		String dstring=args[1];
		assert(dstring.length()==1);
		
		final String[] lines;
		if(new File(fname).exists()){
			lines=TextFile.toStringLines(fname);
		}else{
			lines=new String[] {fname};
		}
		
		LineParserS lp=new LineParserS2(dstring.charAt(0));
		for(String line : lines) {
			lp.set(line);
			System.out.println(lp);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public LineParserS2(char delimiter_) {delimiter=delimiter_;}

	public LineParserS2(int delimiter_) {
		assert(delimiter_>=0 && delimiter_<=Character.MAX_VALUE);
		delimiter=(char)delimiter_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public LineParserS2 set(byte[] line_) {
		assert(false) : "Use byte version.";
		return set(new String(line_));
	}

	@Override
	public LineParserS2 set(byte[] line_, int maxTerm) {
		assert(false) : "Use byte version.";
		return set(new String(line_), maxTerm);
	}
	
	public LineParserS2 set(String line_) {
		reset();
		line=line_;
		return this;
	}
	
	public LineParserS2 set(String line_, int maxTerm) {
		return set(line_);
	}
	
	public LineParserS2 clear() {
		line=null;
		a=b=currentTerm=-1;
		return this;
	}
	
	public LineParserS2 reset() {
		a=b=currentTerm=-1;
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public int parseInt() {
		advance();
		return Parse.parseInt(line, a, b);
	}
	
	public long parseLong() {
		advance();
		return Parse.parseLong(line, a, b);
	}
	
	public float parseFloat() {
		advance();
		return Parse.parseFloat(line, a, b);
	}
	
	public double parseDouble() {
		advance();
		return Parse.parseDouble(line, a, b);
	}
	
	public byte parseByte(int offset) {
		advance();
		int index=a+offset;
		assert(index<b);
		return (byte)line.charAt(index);
	}
	
	public String parseString() {
		int len=advance();
		assert(b>a) : currentTerm+", "+line;
		return line.substring(a, b);
	}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public int parseInt(int term) {
		advanceTo(term);
		return Parse.parseInt(line, a, b);
	}

	/**
	 * Advances to specified term and parses it as a long.
	 * @param term The term number to advance to (0-based)
	 * @return The parsed long value
	 */
	@Override
	public long parseLong(int term) {
		advanceTo(term);
		return Parse.parseLong(line, a, b);
	}

	/**
	 * Advances to specified term and parses it as a float.
	 * @param term The term number to advance to (0-based)
	 * @return The parsed float value
	 */
	@Override
	public float parseFloat(int term) {
		advanceTo(term);
		return Parse.parseFloat(line, a, b);
	}

	@Override
	public double parseDouble(int term) {
		advanceTo(term);
		return Parse.parseDouble(line, a, b);
	}

	@Override
	public byte parseByte(int term, int offset) {
		advanceTo(term);
		int index=a+offset;
		assert(index<b);
		return (byte)line.charAt(index);
	}

	@Override
	public char parseChar(int term, int offset) {
		return (char)parseByte(term, offset);
	}
	
	@Override
	public byte[] parseByteArray(int term) {
		int len=advanceTo(term);
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}
	
	@Override
	public byte[] parseByteArrayFromCurrentField() {
		int len=b-a;
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}

	@Override
	public String parseString(int term) {
		int len=advanceTo(term);
		return line.substring(a, b);
	}

	@Override
	public ByteBuilder appendTerm(ByteBuilder bb, int term) {
		final int len=advanceTo(term);
		for(int i=a; i<b; i++) {bb.append(line.charAt(i));}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public int parseIntFromCurrentField() {
		return Parse.parseInt(line, a, b);
	}
	
	@Override
	public String parseStringFromCurrentField() {
		return line.substring(a, b);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Query Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean startsWith(String s) {
		return line.startsWith(s);
	}
	
	/**
	 * Tests if the line starts with the specified character.
	 * @param c The character to test for
	 * @return true if the line starts with the specified character
	 */
	@Override
	public boolean startsWith(char c) {
		return Tools.startsWith(line, c);
	}
	
	/**
	 * Tests if the line starts with the specified byte value as character.
	 * @param b The byte value to test for
	 * @return true if the line starts with the specified byte as character
	 */
	@Override
	public boolean startsWith(byte b) {
		return Tools.startsWith(line, b);
	}
	
	/**
	 * Tests if the specified term starts with the given string.
	 * @param s The prefix to test for
	 * @param term The term number to check (0-based)
	 * @return true if the term starts with the specified string
	 */
	@Override
	public boolean termStartsWith(String s, int term) {
		final int len=advanceTo(term);
		if(len<s.length()) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(c!=line.charAt(a+i)) {return false;}
		}
		return true;
	}
	
	/**
	 * Tests if the specified term equals the given string.
	 * @param s The string to compare against
	 * @param term The term number to check (0-based)
	 * @return true if the term equals the specified string
	 */
	@Override
	public boolean termEquals(String s, int term) {
		final int len=advanceTo(term);
		if(len!=s.length()) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(c!=line.charAt(a+i)) {return false;}
		}
		return true;
	}
	
	/**
	 * Tests if the specified term equals the given character.
	 * @param c The character to compare against
	 * @param term The term number to check (0-based)
	 * @return true if the term equals the specified character
	 */
	@Override
	public boolean termEquals(char c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	/**
	 * Tests if the specified term equals the given byte value as character.
	 * @param c The byte value to compare against
	 * @param term The term number to check (0-based)
	 * @return true if the term equals the specified byte as character
	 */
	@Override
	public boolean termEquals(byte c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	/**
	 * Increments the start position of current field by specified amount.
	 * @param amt Amount to increment the start position
	 * @return The new length of the current field
	 */
	@Override
	public int incrementA(int amt) {
		a+=amt;
		return b-a;
	}
	
	/**
	 * Increments the start position of current field by specified amount.
	 * Note: This appears to increment 'a' instead of 'b' as method name suggests.
	 * @param amt Amount to increment the position
	 * @return The new length of the current field
	 */
	@Override
	public int incrementB(int amt) {
		a+=amt;
		return b-a;
	}

	/**
	 * Returns the length of the specified term without changing current position.
	 * Temporarily advances to the term, measures length, then restores position.
	 * @param term The term number to measure (0-based)
	 * @return The length of the specified term
	 */
	@Override
	public int length(int term) {
		int a0=a, b0=b, c0=currentTerm;
		int len=advanceTo(term);
		a=a0; b=b0; currentTerm=c0;
		return len;
	}

	/** Returns the length of the current field.
	 * @return The number of characters in the current field */
	@Override
	public int currentFieldLength() {
		return b-a;
	}

	/** Tests if there are more characters to parse in the line.
	 * @return true if the current position has not reached end of line */
	@Override
	public boolean hasMore() {
		return b<line.length();
	}

	/** Returns the total length of the line being parsed.
	 * @return The number of characters in the line */
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
	/*----------------        Advance Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Sets the parsing bounds to the specified term.
	 * @param term The term number to advance to (0-based)
	 * @return The length of the term
	 */
	@Override
	public int setBounds(int term) {
		return advanceTo(term);
	}
	
	public final int advance() {
		currentTerm++;
		b++;
		a=b;
		while(b<line.length() && line.charAt(b)!=delimiter){b++;}
		return b-a;
	}
	
	public void advanceBy(int terms) {
		for(; terms>0; terms--) {
			advance();
		}
	}
	
	//Advances to term before toTerm
	public void advanceToBefore(int toTerm) {
		assert(toTerm>=currentTerm) : "Can't advance backwards: "+currentTerm+">"+toTerm;
		for(toTerm--; currentTerm<toTerm;) {
			advance();
		}
	}
	
	//Advances to actual term
	private int advanceTo(int toTerm) {
		assert(toTerm>=currentTerm) : "Can't advance backwards: "+currentTerm+">"+toTerm;
		for(toTerm--; currentTerm<=toTerm;) {
			advance();
		}
		return b-a;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns string representation of all terms in the line as a list.
	 * @return String representation of parsed terms */
	@Override
	public String toString() {
		return toList().toString();
	}
	
	public ArrayList<String> toList(){
		ArrayList<String> list=new ArrayList<String>();
		do{
			list.add(parseString());
		}while(b<line.length());
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int a=-1;
	private int b=-1;
	private int currentTerm=-1;
	private String line;
	
	public final char delimiter;
	
}
