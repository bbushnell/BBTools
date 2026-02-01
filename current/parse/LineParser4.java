package parse;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.TextFile;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Line parser that uses multiple ordered delimiters for field separation.
 * Each delimiter in the sequence is used for successive field boundaries,
 * enabling parsing of complex formats like ",. ,," where different separators
 * are used for different field positions.
 *
 * @author Brian Bushnell
 * @date May 24, 2023
 */
public final class LineParser4 implements LineParser {
	
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
		
		LineParser lp=new LineParser4(dstring);
		for(String line : lines) {
			lp.set(line.getBytes());
			System.out.println(lp);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public LineParser4(String delimiters_) {
		this(delimiters_.getBytes());
	}

	public LineParser4(byte[] delimiters_) {
		delimiters=delimiters_;
		maxDPos=delimiters.length-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public LineParser4 set(byte[] line_) {
		clear();
		line=line_;
		for(int len=advance(); b<line.length; len=advance()) {
			bounds.add(b);
		}
		bounds.add(b);
		return this;
	}

	@Override
	public LineParser4 set(byte[] line_, int maxTerm) {
		clear();
		line=line_;
		//TODO: test performance of presumably safer loop below
//		for(int term=0; term<=maxTerm && b<line.length; term++) {
		for(int term=0; term<=maxTerm; term++) {
			int len=advance();
			bounds.add(b);
		}
		return this;
	}
	
//	public LineParser4 set(String line_) {
//		assert(false) : "Use string version.";
//		return set(line_.getBytes());
//	}
//	
//	public LineParser4 set(String line_, int maxTerm) {
//		assert(false) : "Use string version.";
//		return set(line_.getBytes(), maxTerm);
//	}
	
	@Override
	public LineParser4 clear() {
		delimiterPos=0;
		line=null;
		a=b=-1;
		bounds.clear();
		return this;
	}
	
	/**
	 * Reset method that does nothing for this implementation.
	 * Included for interface compliance.
	 * @return This parser instance for method chaining
	 */
	@Override
	public LineParser4 reset() {
		//Does nothing for this class
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public int terms() {return bounds.size();}
	
	/**
	 * Parses the specified field as an integer.
	 * @param term Field index to parse
	 * @return Integer value of the field
	 */
	@Override
	public int parseInt(int term) {
		setBounds(term);
		return Parse.parseInt(line, a, b);
	}
	
	public int parseInt(int term, int from, int to) {
		setBounds(term);
		return Parse.parseInt(line, a+from, Tools.min(line.length, a+to));
//		return Parse.parseInt(line, a+from, Tools.min(b, a+to));
	}
	
	/**
	 * Parses the specified field as a long integer.
	 * @param term Field index to parse
	 * @return Long value of the field
	 */
	@Override
	public long parseLong(int term) {
		setBounds(term);
		return Parse.parseLong(line, a, b);
	}
	
	/**
	 * Parses the specified field as a float.
	 * @param term Field index to parse
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
	 * Extracts a single byte from specified field at given offset.
	 * @param term Field index
	 * @param offset Byte offset within the field
	 * @return Byte value at the specified position
	 */
	@Override
	public byte parseByte(int term, int offset) {
		setBounds(term);
		final int index=a+offset;
		assert(index<b);
		return line[index];
	}
	
//	@Override
//	public char parseChar(int term, int offset) {
//		return (char)parseByte(term, offset);
//	}
	
	@Override
	public byte[] parseByteArray(int term) {
		final int len=setBounds(term);
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line[a+1];}
		return ret;
	}
	
	/** Returns the current field as a byte array */
	@Override
	public byte[] parseByteArrayFromCurrentField() {
		int len=b-a;
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++) {ret[i]=(byte)line[a+1];}
		return ret;
	}
	
	/**
	 * Parses the specified field as a String.
	 * @param term Field index to parse
	 * @return String representation of the field
	 */
	@Override
	public String parseString(int term) {
		final int len=setBounds(term);
		return new String(line, a, len, StandardCharsets.US_ASCII);
	}

	/**
	 * Appends the specified field to a ByteBuilder.
	 * @param bb ByteBuilder to append to
	 * @param term Field index to append
	 * @return The modified ByteBuilder for method chaining
	 */
	@Override
	public ByteBuilder appendTerm(ByteBuilder bb, int term) {
		final int len=setBounds(term);
		for(int i=a; i<b; i++) {bb.append(line[i]);}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Parses the current field as an integer */
	@Override
	public int parseIntFromCurrentField() {
		return Parse.parseInt(line, a, b);
	}
	
	/** Returns the current field as a String */
	@Override
	public String parseStringFromCurrentField() {
		final int len=b-a;
		return new String(line, a, len, StandardCharsets.US_ASCII);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Query Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Checks if the line starts with the given string.
	 * @param s String prefix to check for
	 * @return true if line starts with the string
	 */
	@Override
	public boolean startsWith(String s) {
		return Tools.startsWith(line, s);
	}
	
	/**
	 * Checks if the line starts with the given character.
	 * @param c Character prefix to check for
	 * @return true if line starts with the character
	 */
	@Override
	public boolean startsWith(char c) {
		return Tools.startsWith(line, c);
	}
	
	/**
	 * Checks if the line starts with the given byte.
	 * @param b Byte prefix to check for
	 * @return true if line starts with the byte
	 */
	@Override
	public boolean startsWith(byte b) {
		return Tools.startsWith(line, b);
	}
	
	/**
	 * Checks if the specified field starts with given string.
	 * @param s String prefix to check for
	 * @param term Field index to check
	 * @return true if field starts with the string
	 */
	@Override
	public boolean termStartsWith(String s, int term) {
		final int len=setBounds(term);
		if(len<s.length()) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(c!=line[a+i]) {return false;}
		}
		return true;
	}
	
	/**
	 * Checks if the specified field exactly equals given string.
	 * @param s String to compare against
	 * @param term Field index to check
	 * @return true if field equals the string
	 */
	@Override
	public boolean termEquals(String s, int term) {
		final int len=setBounds(term);
		if(len!=s.length()) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(c!=line[a+i]) {return false;}
		}
		return true;
	}
	
	/**
	 * Checks if the specified field exactly equals given character.
	 * @param c Character to compare against
	 * @param term Field index to check
	 * @return true if field equals the character
	 */
	@Override
	public boolean termEquals(char c, int term) {
		final int len=setBounds(term);
		return len==1 && line[a]==c;
	}
	
	/**
	 * Checks if the specified field exactly equals given byte.
	 * @param c Byte to compare against
	 * @param term Field index to check
	 * @return true if field equals the byte
	 */
	@Override
	public boolean termEquals(byte c, int term) {
		final int len=setBounds(term);
		return len==1 && line[a]==c;
	}
	
	/**
	 * Increments the start position of current field boundary.
	 * @param amt Amount to increment by
	 * @return New length of current field
	 */
	@Override
	public int incrementA(int amt) {
		a+=amt;
		return b-a;
	}
	
	/**
	 * Increments the start position (appears to be implementation error).
	 * @param amt Amount to increment by
	 * @return New length of current field
	 */
	@Override
	public int incrementB(int amt) {
		a+=amt;
		return b-a;
	}

	/**
	 * Returns the length of the specified field.
	 * @param term Field index
	 * @return Length of the field in characters
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

	/** Returns true if there are more characters to parse in the line */
	@Override
	public boolean hasMore() {
		return b<line.length;
	}

	/** Returns the total length of the line being parsed */
	@Override
	public int lineLength() {
		return line.length;
	}

	@Override
	public byte[] line() {return line;}
	
	@Override
	public int a() {return a;}
	
	@Override
	public int b() {return b;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets the current field boundaries (a,b) to the specified term.
	 * Updates internal pointers to define the field range.
	 * @param term Field index to set boundaries for
	 * @return Length of the specified field
	 */
	@Override
	public int setBounds(int term){
		a=(term==0 ? 0 : bounds.get(term-1)+1);
		b=bounds.get(term);
		return b-a;
	}
	
	private int advance() {
		byte delimiter=(delimiterPos<delimiters.length ? delimiters[delimiterPos] : 0);
		delimiterPos++;
		b++;
		a=b;
		while(b<line.length && delimiter!=line[b]){b++;}
		return b-a;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns string representation of all parsed fields */
	@Override
	public String toString() {
		return toList().toString();
	}
	
	/** Converts all parsed fields to a list of strings.
	 * @return ArrayList containing string representations of all fields */
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
	private byte[] line;
	
	public final byte[] delimiters;
	private final int maxDPos;
	private int delimiterPos=0;
	
}
