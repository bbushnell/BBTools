package shared;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextFile;
import structures.ByteBuilder;

/** Similar speed, but less powerful.
 * Main advantage is having a bounded memory footprint for very long lines.
 * 
 * @author Brian Bushnell
 * @date May 24, 2023
 *
 */
public final class LineParserS2 implements LineParserS {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	//For testing
	//Syntax: LineParser fname/literal delimiter 
	/**
	 * Test method for LineParserS2 functionality.
	 * Accepts a filename/literal string and single character delimiter.
	 * @param args Command-line arguments: [filename/literal] [delimiter]
	 */
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

	/** Constructs a LineParserS2 with the specified character delimiter.
	 * @param delimiter_ The character used to separate fields */
	public LineParserS2(char delimiter_) {delimiter=delimiter_;}

	/** Constructs a LineParserS2 with the specified integer delimiter converted to char.
	 * @param delimiter_ The delimiter as an integer (must be valid char value) */
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
	
	/**
	 * Sets the line to parse as a string and resets parsing state.
	 * @param line_ The string to parse
	 * @return This parser instance for method chaining
	 */
	public LineParserS2 set(String line_) {
		reset();
		line=line_;
		return this;
	}
	
	/**
	 * Sets the line to parse as a string. maxTerm parameter is ignored.
	 * @param line_ The string to parse
	 * @param maxTerm Maximum number of terms (ignored in this implementation)
	 * @return This parser instance for method chaining
	 */
	public LineParserS2 set(String line_, int maxTerm) {
		return set(line_);
	}
	
	/** Clears the parser by setting line to null and resetting all position markers.
	 * @return This parser instance for method chaining */
	public LineParserS2 clear() {
		line=null;
		a=b=currentTerm=-1;
		return this;
	}
	
	/** Resets parsing position to the beginning without clearing the line.
	 * @return This parser instance for method chaining */
	public LineParserS2 reset() {
		a=b=currentTerm=-1;
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Advances to next field and parses it as an integer.
	 * @return The parsed integer value */
	public int parseInt() {
		advance();
		return Parse.parseInt(line, a, b);
	}
	
	/** Advances to next field and parses it as a long.
	 * @return The parsed long value */
	public long parseLong() {
		advance();
		return Parse.parseLong(line, a, b);
	}
	
	/** Advances to next field and parses it as a float.
	 * @return The parsed float value */
	public float parseFloat() {
		advance();
		return Parse.parseFloat(line, a, b);
	}
	
	/** Advances to next field and parses it as a double.
	 * @return The parsed double value */
	public double parseDouble() {
		advance();
		return Parse.parseDouble(line, a, b);
	}
	
	/**
	 * Advances to next field and returns the byte at specified offset within that field.
	 * @param offset Position within the current field to extract byte from
	 * @return The byte value at the specified offset
	 */
	public byte parseByte(int offset) {
		advance();
		int index=a+offset;
		assert(index<b);
		return (byte)line.charAt(index);
	}
	
	/** Advances to next field and returns it as a string.
	 * @return The current field as a string */
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

	@Override
	public long parseLong(int term) {
		advanceTo(term);
		return Parse.parseLong(line, a, b);
	}

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
	
	@Override
	public boolean startsWith(char c) {
		return Tools.startsWith(line, c);
	}
	
	@Override
	public boolean startsWith(byte b) {
		return Tools.startsWith(line, b);
	}
	
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
	
	@Override
	public boolean termEquals(char c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	@Override
	public boolean termEquals(byte c, int term) {
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}
	
	@Override
	public int incrementA(int amt) {
		a+=amt;
		return b-a;
	}
	
	@Override
	public int incrementB(int amt) {
		a+=amt;
		return b-a;
	}

	@Override
	public int length(int term) {
		int a0=a, b0=b, c0=currentTerm;
		int len=advanceTo(term);
		a=a0; b=b0; currentTerm=c0;
		return len;
	}

	@Override
	public int currentFieldLength() {
		return b-a;
	}

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
	/*----------------        Advance Methods       ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int setBounds(int term) {
		return advanceTo(term);
	}
	
	/**
	 * Advances to the next field in the line.
	 * Increments term counter and finds next delimiter boundary.
	 * @return The length of the new current field
	 */
	public final int advance() {
		currentTerm++;
		b++;
		a=b;
		while(b<line.length() && line.charAt(b)!=delimiter){b++;}
		return b-a;
	}
	
	/** Advances by the specified number of terms.
	 * @param terms Number of terms to advance */
	public void advanceBy(int terms) {
		for(; terms>0; terms--) {
			advance();
		}
	}
	
	//Advances to term before toTerm
	/** Advances to the term immediately before the specified term.
	 * @param toTerm The target term number (will advance to toTerm-1) */
	public void advanceToBefore(int toTerm) {
		assert(toTerm>=currentTerm) : "Can't advance backwards: "+currentTerm+">"+toTerm;
		for(toTerm--; currentTerm<toTerm;) {
			advance();
		}
	}
	
	//Advances to actual term
	/**
	 * Advances to the specified term number.
	 * @param toTerm The term number to advance to (0-based)
	 * @return The length of the target term
	 */
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
	
	@Override
	public String toString() {
		return toList().toString();
	}
	
	/** Parses all remaining terms in the line and returns them as a list of strings.
	 * @return ArrayList containing all terms as strings */
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
	
	/** Start position of the current field in the line */
	private int a=-1;
	/** End position of the current field in the line */
	private int b=-1;
	/** Index of the current term being processed (0-based) */
	private int currentTerm=-1;
	/** The string line being parsed */
	private String line;
	
	/** The character used to separate fields in the line */
	public final char delimiter;
	
}
