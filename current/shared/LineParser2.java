package shared;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import structures.ByteBuilder;

/** Similar speed, but less powerful.
 * Main advantage is having a bounded memory footprint for very long lines.
 * 
 * @author Brian Bushnell
 * @date May 24, 2023
 *
 */
public final class LineParser2 implements LineParser {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	//For testing
	//Syntax: LineParser fname/literal delimiter 
	/** Test method for LineParser2 functionality.
	 * @param args Command-line arguments: [filename/literal, delimiter] */
	public static void main(String[] args) {
		assert(args.length==2);
		String fname=args[0];
		String dstring=args[1];
		assert(dstring.length()==1);
		
		final ArrayList<byte[]> lines;
		if(new File(fname).exists()){
			lines=ByteFile.toLines(fname);
		}else{
			lines=new ArrayList<byte[]>(1);
			lines.add(fname.getBytes());
		}
		
		LineParser2 lp=new LineParser2(dstring.charAt(0));
		for(byte[] line : lines) {
			lp.set(line);
			System.out.println(lp);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs a LineParser2 with the specified byte delimiter.
	 * @param delimiter_ The byte value to use as field delimiter */
	public LineParser2(byte delimiter_) {delimiter=delimiter_;}

	/** Constructs a LineParser2 with the specified ASCII delimiter.
	 * @param delimiter_ ASCII value of delimiter (must be 0-127) */
	public LineParser2(int delimiter_) {
		assert(delimiter_>=0 && delimiter_<=127);
		delimiter=(byte)delimiter_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public LineParser2 set(byte[] line_) {
		reset();
		line=line_;
		return this;
	}
	
	@Override
	public LineParser2 set(byte[] line_, int maxTerm) {
		return set(line_);
	}
	
	@Override
	public LineParser2 clear() {
		line=null;
		a=b=currentTerm=-1;
		return this;
	}
	
	@Override
	public LineParser2 reset() {
		a=b=currentTerm=-1;
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Advances to the next field and parses it as an integer.
	 * @return The integer value of the next field */
	public int parseInt() {
		advance();
		return Parse.parseInt(line, a, b);
	}
	
	/** Advances to the next field and parses it as a long.
	 * @return The long value of the next field */
	public long parseLong() {
		advance();
		return Parse.parseLong(line, a, b);
	}
	
	/** Advances to the next field and parses it as a float.
	 * @return The float value of the next field */
	public float parseFloat() {
		advance();
		return Parse.parseFloat(line, a, b);
	}
	
	/** Advances to the next field and parses it as a double.
	 * @return The double value of the next field */
	public double parseDouble() {
		advance();
		return Parse.parseDouble(line, a, b);
	}
	
	/**
	 * Advances to the next field and returns a byte at the specified offset.
	 * @param offset Offset from the start of the current field
	 * @return The byte at the specified offset within the field
	 */
	public byte parseByte(int offset) {
		advance();
		int index=a+offset;
		assert(index<b);
		return line[index];
	}
	
	/** Advances to the next field and parses it as a String.
	 * @return The String value of the next field */
	public String parseString() {
		int len=advance();
		return new String(line, a, len, StandardCharsets.US_ASCII);
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
		return line[index];
	}
	
	@Override
	public byte[] parseByteArray(int term) {
		int len=advanceTo(term);
		return Arrays.copyOfRange(line, a, b);
	}
	
	@Override
	public byte[] parseByteArrayFromCurrentField() {
		return Arrays.copyOfRange(line, a, b);
	}

	@Override
	public String parseString(int term) {
		int len=advanceTo(term);
		return new String(line, a, len, StandardCharsets.US_ASCII);
	}

	@Override
	public ByteBuilder appendTerm(ByteBuilder bb, int term) {
		final int len=setBounds(term);
		for(int i=a; i<b; i++) {bb.append(line[i]);}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public int parseIntFromCurrentField() {
		return Parse.parseInt(line, a, b);
	}
	
	@Override
	public String parseStringFromCurrentField() {
		return new String(line, a, b-a, StandardCharsets.US_ASCII);
	}

	/** Returns the first byte of the current field without advancing.
	 * @return The first byte of the current field */
	public byte parseByteFromCurrentField() {
		return line[a];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Query Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean startsWith(String s) {
		return Tools.startsWith(line, s);
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
			if(c!=line[a+i]) {return false;}
		}
		return true;
	}
	
	@Override
	public boolean termEquals(String s, int term) {
		final int len=advanceTo(term);
		if(len!=s.length()) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(c!=line[a+i]) {return false;}
		}
		return true;
	}
	
	@Override
	public boolean termEquals(char c, int term) {
		final int len=setBounds(term);
		return len==1 && line[a]==c;
	}
	
	@Override
	public boolean termEquals(byte c, int term) {
		final int len=setBounds(term);
		return len==1 && line[a]==c;
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
	public byte[] line() {return line;}
	
	@Override
	public int a() {return a;}
	
	@Override
	public int b() {return b;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Advance Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public int setBounds(int term){
		return advanceTo(term);
	}
	
	/**
	 * Advances to the next field by moving past the delimiter.
	 * Sets boundaries (a,b) to encompass the next field.
	 * @return Length of the field that was advanced to
	 */
	public final int advance() {
		currentTerm++;
		b++;
		a=b;
		while(b<line.length && line[b]!=delimiter){b++;}
		return b-a;
	}
	
	/** Advances forward by the specified number of terms/fields.
	 * @param terms Number of fields to advance */
	public void advanceBy(int terms) {
		for(; terms>0; terms--) {
			advance();
		}
	}
	
	//Advances to term before toTerm
	/** Advances to the field just before the specified term number.
	 * @param toTerm Target term number (will stop at toTerm-1) */
	public void advanceToBefore(int toTerm) {
		assert(toTerm>=currentTerm) : "Can't advance backwards: "+currentTerm+">"+toTerm;
		for(toTerm--; currentTerm<toTerm;) {
			advance();
		}
	}
	
	//Advances to actual term
	/**
	 * Advances to the specified term number.
	 * @param toTerm Target term number to advance to
	 * @return Length of the field at the target term
	 */
	private int advanceTo(int toTerm) {
		assert(toTerm>=currentTerm) : "Can't advance backwards: "+currentTerm+">"+toTerm;
		for(toTerm--; currentTerm<=toTerm;) {
			advance();
		}
		return b-a;
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
	
	/**
	 * Manually sets the field boundaries.
	 * @param a_ Start position of current field
	 * @param b_ End position of current field
	 */
	public void setBounds(int a_, int b_) {
		a=a_;
		b=b_;
	}

	@Override
	public boolean hasMore() {
		return b<line.length;
	}

	@Override
	public int lineLength() {
		return line.length;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString() {
		return /*toList().toString()+"\n"*/"a="+a+", b="+b+", line.length="+line.length;
	}
	
	/** Parses all remaining fields in the line as Strings.
	 * @return ArrayList containing all remaining fields as Strings */
	public ArrayList<String> toList(){
		ArrayList<String> list=new ArrayList<String>();
		do{
			list.add(parseString());
		}while(b<line.length);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Start index of the current field */
	private int a=-1;
	/** End index of the current field (exclusive) */
	private int b=-1;
	/** Index of the current term/field being processed */
	private int currentTerm=-1;
	/** The line being parsed as a byte array */
	private byte[] line;
	
	/** The byte value used as field delimiter */
	public final byte delimiter;
	
}
