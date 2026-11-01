package parse;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextFile;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Memory-efficient line parser for delimited text strings.
 * Uses bounded auxiliary memory for arbitrarily long retained input strings.
 * Implements the LineParserS interface using string-based parsing with character delimiters.
 *
 * <p>Forward-only UTF-16 scanner. maxTerm is ignored; indexed access cannot move
 * backwards without reset(). No-argument reads advance to the next field; check
 * hasMore() before calling them. Empty and trailing fields are retained. Input
 * strings are retained and each mutable instance belongs to one thread.
 *
 * @author Brian Bushnell
 * @date May 24, 2023
 */
public final class LineParserS2 implements LineParserS{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	//For testing
	//Syntax: LineParser fname/literal delimiter 
	/** Demonstrates parsing a file or literal input with the supplied delimiter. */
	public static void main(String[] args){
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
		for(String line : lines){
			lp.set(line);
			System.out.println(lp);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a parser with the given UTF-16 delimiter code unit. */
	public LineParserS2(char delimiter_){delimiter=delimiter_;}

	/** Creates a character-delimited parser; asserts the inclusive UTF-16 range 0..65535. */
	public LineParserS2(int delimiter_){
		assert(delimiter_>=0 && delimiter_<=Character.MAX_VALUE);
		delimiter=(char)delimiter_;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public LineParserS2 set(byte[] line_){
		assert(false) : "Use byte version.";
		return set(new String(line_));
	}

	@Override
	public LineParserS2 set(byte[] line_, int maxTerm){
		assert(false) : "Use byte version.";
		return set(new String(line_), maxTerm);
	}

	public LineParserS2 set(String line_){
		reset();
		line=line_;
		return this;
	}

	public LineParserS2 set(String line_, int maxTerm){
		return set(line_);
	}

	public LineParserS2 clear(){
		line=null;
		a=b=currentTerm=-1;
		return this;
	}

	public LineParserS2 reset(){
		a=b=currentTerm=-1;
		return this;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Advances to the next existing field and delegates numeric conversion to Parse. */
	public int parseInt(){
		advance();
		return Parse.parseInt(line, a, b);
	}

	/** Advances to the next existing field and delegates numeric conversion to Parse. */
	public long parseLong(){
		advance();
		return Parse.parseLong(line, a, b);
	}

	/** Advances to the next existing field and delegates numeric conversion to Parse. */
	public float parseFloat(){
		advance();
		return Parse.parseFloat(line, a, b);
	}

	/** Advances to the next existing field and delegates numeric conversion to Parse. */
	public double parseDouble(){
		advance();
		return Parse.parseDouble(line, a, b);
	}

	/** Advances to the next existing field and narrows its UTF-16 code unit at offset to a byte. */
	public byte parseByte(int offset){
		advance();
		int index=a+offset;
		assert(index<b);
		return (byte)line.charAt(index);
	}

	public String parseString(){
		int len=advance();
		//[parse/LineParserS2#001] No assert(b>a) here: an empty field (a==b) is valid delimited
		//data and must yield "", matching the byte-twin LineParser2.parseString() and all S-siblings.
		//The real "called past end" guard is hasMore(); the removed assert wrongly crashed valid empty fields.
		return line.substring(a, b);
	}

	/*--------------------------------------------------------------*/

	@Override
	public int parseInt(int term){
		advanceTo(term);
		return Parse.parseInt(line, a, b);
	}

	/**
	 * Advances to specified term and parses it as a long.
	 * @param term The term number to advance to (0-based)
	 * @return The parsed long value
	 */
	@Override
	public long parseLong(int term){
		advanceTo(term);
		return Parse.parseLong(line, a, b);
	}

	/**
	 * Advances to specified term and parses it as a float.
	 * @param term The term number to advance to (0-based)
	 * @return The parsed float value
	 */
	@Override
	public float parseFloat(int term){
		advanceTo(term);
		return Parse.parseFloat(line, a, b);
	}

	@Override
	public double parseDouble(int term){
		advanceTo(term);
		return Parse.parseDouble(line, a, b);
	}

	@Override
	public byte parseByte(int term, int offset){
		advanceTo(term);
		int index=a+offset;
		assert(index<b);
		return (byte)line.charAt(index);
	}

	/** Selects a forward-reachable field and returns its original UTF-16 code unit at a valid offset. */
	@Override
	public char parseChar(int term, int offset){
		advanceTo(term);
		final int index=a+offset;
		assert(index<b) : "Character offset must stay inside the selected field; term="+term+", offset="+offset;
		return line.charAt(index);
	}

	@Override
	public byte[] parseByteArray(int term){
		int len=advanceTo(term);
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++){ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}

	@Override
	public byte[] parseByteArrayFromCurrentField(){
		int len=b-a;
		byte[] ret=new byte[len];
		for(int i=0; i<len; i++){ret[i]=(byte)line.charAt(a+i);}
		return ret;
	}

	@Override
	public String parseString(int term){
		int len=advanceTo(term);
		return line.substring(a, b);
	}

	@Override
	public ByteBuilder appendTerm(ByteBuilder bb, int term){
		final int len=advanceTo(term);
		for(int i=a; i<b; i++){bb.append(line.charAt(i));}
		return bb;
	}

	/*--------------------------------------------------------------*/

	@Override
	public int parseIntFromCurrentField(){
		return Parse.parseInt(line, a, b);
	}

	@Override
	public String parseStringFromCurrentField(){
		return line.substring(a, b);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Query Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public boolean startsWith(String s){
		return line.startsWith(s);
	}

	/**
	 * Tests if the line starts with the specified character.
	 * @param c The character to test for
	 * @return true if the line starts with the specified character
	 */
	@Override
	public boolean startsWith(char c){
		return Tools.startsWith(line, c);
	}

	/**
	 * Tests if the line starts with the specified byte value as character.
	 * @param b The byte value to test for
	 * @return true if the line starts with the specified byte as character
	 */
	@Override
	public boolean startsWith(byte b){
		return Tools.startsWith(line, b);
	}

	/**
	 * Tests if the specified term starts with the given string.
	 * @param s The prefix to test for
	 * @param term The term number to check (0-based)
	 * @return true if the term starts with the specified string
	 */
	@Override
	public boolean termStartsWith(String s, int term){
		final int len=advanceTo(term);
		if(len<s.length()){return false;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c!=line.charAt(a+i)){return false;}
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
	public boolean termEquals(String s, int term){
		final int len=advanceTo(term);
		if(len!=s.length()){return false;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c!=line.charAt(a+i)){return false;}
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
	public boolean termEquals(char c, int term){
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
	public boolean termEquals(byte c, int term){
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}

	/**
	 * Increments the start position of current field by specified amount.
	 * @param amt Amount to increment the start position
	 * @return The new length of the current field
	 */
	@Override
	public int incrementA(int amt){
		a+=amt;
		return b-a;
	}

	/**
	 * Increments the end position (right bound b) of the current field by the specified amount,
	 * growing (or with negative amt, shrinking) the current field on the right.
	 * @param amt Amount to increment the end position
	 * @return The new length of the current field
	 */
	@Override
	public int incrementB(int amt){
		b+=amt;
		return b-a;
	}

	/**
	 * Returns the length of the specified term without changing current position.
	 * Temporarily advances to the term, measures length, then restores position.
	 * @param term The term number to measure (0-based)
	 * @return The length of the specified term
	 */
	@Override
	public int length(int term){
		int a0=a, b0=b, c0=currentTerm;
		int len=advanceTo(term);
		a=a0; b=b0; currentTerm=c0;
		return len;
	}

	/** Returns the length of the current field.
	 * @return The number of characters in the current field */
	@Override
	public int currentFieldLength(){
		return b-a;
	}

	/** Tests if there are more characters to parse in the line.
	 * @return true if the current position has not reached end of line */
	@Override
	public boolean hasMore(){
		return b<line.length();
	}

	/** Returns the total length of the line being parsed.
	 * @return The number of characters in the line */
	@Override
	public int lineLength(){
		return line.length();
	}

	/** Returns the current line being parsed */
	@Override
	public String line(){return line;}

	/** Returns the start position of the current field */
	@Override
	public int a(){return a;}

	/** Returns the end position of the current field */
	@Override
	public int b(){return b;}

	/*--------------------------------------------------------------*/
	/*----------------        Advance Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Sets the parsing bounds to the specified term.
	 * @param term The term number to advance to (0-based)
	 * @return The length of the term
	 */
	@Override
	public int setBounds(int term){
		return advanceTo(term);
	}

	/** Moves to the next field; requires hasMore() before entry and retains empty fields. */
	public final int advance(){
		currentTerm++;
		b++;
		a=b;
		while(b<line.length() && line.charAt(b)!=delimiter){b++;}
		return b-a;
	}

	/** Advances terms times; nonpositive counts do nothing. Caller must ensure fields exist. */
	public void advanceBy(int terms){
		for(; terms>0; terms--){
			advance();
		}
	}

	//Advances to term before toTerm
	/** Advances until immediately before toTerm, without moving backwards. */
	public void advanceToBefore(int toTerm){
		assert(toTerm>=currentTerm) : "Can't advance backwards: "+currentTerm+">"+toTerm;
		for(toTerm--; currentTerm<toTerm;){
			advance();
		}
	}

	//Advances to actual term
	/** Selects a forward-reachable existing field; reselecting the current field is allowed. */
	private int advanceTo(int toTerm){
		assert(toTerm>=currentTerm) : "Can't advance backwards: "+currentTerm+">"+toTerm;
		for(toTerm--; currentTerm<=toTerm;){
			advance();
		}
		return b-a;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Consumes remaining fields into a list string; requires an unexhausted cursor.
	 * @return String representation of parsed terms */
	@Override
	public String toString(){
		return toList().toString();
	}

	/** Consumes remaining fields, requiring an unexhausted cursor; reset first to include all fields. */
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

	/** Start of the current field in native line units, inclusive; -1 before initialization. */
	private int a=-1;
	/** End of the current field in native line units, exclusive; -1 before initialization. */
	private int b=-1;
	/** Current forward-only field index; -1 before the first advance. */
	private int currentTerm=-1;
	/** Retained input string; null after clear and before set. */
	private String line;

	/** Primary UTF-16 delimiter code unit. */
	public final char delimiter;

}
