package shared;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextFile;
import structures.ByteBuilder;
import structures.IntList;

public final class LineParserS1 implements LineParserS {
	
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
		
		LineParserS lp=new LineParserS1(dstring.charAt(0));
		for(String line : lines) {
			lp.set(line);
			System.out.println(lp);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public LineParserS1(char delimiter_) {delimiter=delimiter_;}

	public LineParserS1(int delimiter_) {
		assert(delimiter_>=0 && delimiter_<=Character.MAX_VALUE);
		delimiter=(char)delimiter_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public LineParserS1 set(byte[] line_) {
		assert(false) : "Use byte version.";
		return set(new String(line_));
	}

	@Override
	public LineParserS1 set(byte[] line_, int maxTerm) {
		assert(false) : "Use byte version.";
		return set(new String(line_), maxTerm);
	}
	
	@Override
	public LineParserS1 set(String line_) {
		clear();
		line=line_;
		for(int len=advance(); b<line.length(); len=advance()) {
			bounds.add(b);
		}
		bounds.add(b);
		return this;
	}
	
	@Override
	public LineParserS1 set(String line_, int maxTerm) {
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
	
	@Override
	public LineParserS1 clear() {
		line=null;
		a=b=-1;
		bounds.clear();
		return this;
	}
	
	@Override
	public LineParserS1 reset() {
		//Does nothing for this class
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public int terms() {return bounds.size();}
	
	@Override
	public int parseInt(int term) {
		setBounds(term);
		return Parse.parseInt(line, a, b);
	}
	
	/**
	 * Parses the specified field as a long integer value.
	 * @param term Zero-based field index to parse
	 * @return Long value of the field
	 */
	@Override
	public long parseLong(int term) {
		setBounds(term);
		return Parse.parseLong(line, a, b);
	}
	
	/**
	 * Parses the specified field as a floating-point value.
	 * @param term Zero-based field index to parse
	 * @return Float value of the field
	 */
	@Override
	public float parseFloat(int term) {
		setBounds(term);
		return Parse.parseFloat(line, a, b);
	}
	
	@Override
	public double parseDouble(int term) {
		setBounds(term);
		return Parse.parseDouble(line, a, b);
	}
	
	/**
	 * Extracts a single byte from within a field at the specified offset.
	 * @param term Zero-based field index
	 * @param offset Character position within the field
	 * @return Byte value of the character at the offset position
	 */
	@Override
	public byte parseByte(int term, int offset) {
		return (byte)parseChar(term, offset);
	}
	
	/**
	 * Extracts a single character from within a field at the specified offset.
	 * @param term Zero-based field index
	 * @param offset Character position within the field (0-based)
	 * @return Character at the specified position within the field
	 */
	@Override
	public char parseChar(int term, int offset) {
		setBounds(term);
		final int index=a+offset;
		assert(index<b);
		return line.charAt(index);
	}
	
	/**
	 * Converts the specified field to a byte array.
	 * Each character in the field becomes one byte in the array.
	 * @param term Zero-based field index to convert
	 * @return Byte array representation of the field content
	 */
	@Override
	public byte[] parseByteArray(int term) {
		int len=setBounds(term);
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}
	
	/**
	 * Converts the current field (as set by setBounds) to a byte array.
	 * Uses the current field boundaries without changing them.
	 * @return Byte array representation of the current field content
	 */
	@Override
	public byte[] parseByteArrayFromCurrentField() {
		int len=b-a;
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}
	
	/**
	 * Extracts the specified field as a string.
	 * @param term Zero-based field index to extract
	 * @return String content of the field
	 */
	@Override
	public String parseString(int term) {
		final int len=setBounds(term);
		return line.substring(a, b);
	}

	/**
	 * Appends the content of a field to a ByteBuilder.
	 * Efficiently adds field characters to the builder without creating intermediate strings.
	 *
	 * @param bb The ByteBuilder to append to
	 * @param term Zero-based field index to append
	 * @return The same ByteBuilder instance for method chaining
	 */
	@Override
	public ByteBuilder appendTerm(ByteBuilder bb, int term) {
		final int len=setBounds(term);
		for(int i=a; i<b; i++) {bb.append(line.charAt(i));}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses the current field as an integer using existing field boundaries.
	 * Does not change the current field position.
	 * @return Integer value of the current field
	 */
	@Override
	public int parseIntFromCurrentField() {
		return Parse.parseInt(line, a, b);
	}
	
	/**
	 * Returns the current field as a string using existing field boundaries.
	 * Does not change the current field position.
	 * @return String content of the current field
	 */
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
	
	@Override
	public boolean startsWith(char c) {
		return Tools.startsWith(line, c);
	}
	
	@Override
	public boolean startsWith(byte b) {
		return Tools.startsWith(line, b);
	}
	
	/**
	 * Tests if the specified field starts with the given string.
	 * @param s String prefix to test for
	 * @param term Zero-based field index to check
	 * @return true if the field starts with the string, false otherwise
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
	 * Tests if the specified field exactly matches the given string.
	 * @param s String to compare against
	 * @param term Zero-based field index to check
	 * @return true if the field equals the string, false otherwise
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
	 * Tests if the specified field consists of exactly one character that matches the given character.
	 * @param c Character to compare against
	 * @param term Zero-based field index to check
	 * @return true if the field is a single character matching c, false otherwise
	 */
	@Override
	public boolean termEquals(char c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	/**
	 * Tests if the specified field consists of exactly one character that matches the given byte value.
	 * @param c Byte value to compare against (treated as character)
	 * @param term Zero-based field index to check
	 * @return true if the field is a single character matching the byte value, false otherwise
	 */
	@Override
	public boolean termEquals(byte c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	/**
	 * Moves the start boundary of the current field forward by the specified amount.
	 * @param amt Number of characters to advance the start position
	 * @return New length of the current field after adjustment
	 */
	@Override
	public int incrementA(int amt) {
		a+=amt;
		return b-a;
	}
	
	/**
	 * Moves the start boundary of the current field forward by the specified amount.
	 * Note: This method appears to have a bug - it modifies 'a' instead of 'b'.
	 * @param amt Number of characters to advance the boundary
	 * @return New length of the current field after adjustment
	 */
	@Override
	public int incrementB(int amt) {
		a+=amt;
		return b-a;
	}

	/**
	 * Returns the length of the specified field in characters.
	 * @param term Zero-based field index
	 * @return Length of the field in characters
	 */
	@Override
	public int length(int term) {
		return setBounds(term);
	}

	/** Returns the length of the current field using existing boundaries.
	 * @return Length of the current field in characters */
	@Override
	public int currentFieldLength() {
		return b-a;
	}

	/** Tests if there are more characters after the current field position.
	 * @return true if more content exists beyond the current position, false otherwise */
	@Override
	public boolean hasMore() {
		return b<line.length();
	}

	@Override
	public int lineLength() {
		return line.length();
	}

	@Override
	public String line() {return line;}
	
	@Override
	public int a() {return a;}
	
	@Override
	public int b() {return b;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public int setBounds(int term){
		a=(term==0 ? 0 : bounds.get(term-1)+1);
		b=bounds.get(term);
		return b-a;
	}
	
	private int advance() {
		b++;
		a=b;
		while(b<line.length() && line.charAt(b)!=delimiter){b++;}
		return b-a;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString() {
		return toList().toString();
	}
	
	/**
	 * Converts all parsed fields into a list of strings.
	 * Creates a new ArrayList containing each field as a separate string element.
	 * @return ArrayList containing all fields as strings
	 */
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
	
	public final char delimiter;
	
}
