package shared;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import structures.ByteBuilder;
import structures.IntList;

public final class LineParser1 implements LineParser {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	//For testing
	//Syntax: LineParser fname/literal delimiter 
	public static void main(String[] args) {
		assert(args.length==2 || args.length==3 || args.length==4);
		String fname=args[0];
		String dstring=Parse.parseSymbol(args[1]);
		final boolean benchmark=args.length>2;
		Shared.SIMD=args.length<4 ? false : 
			(args[3].equalsIgnoreCase("simd") || args[3].equalsIgnoreCase("simd=t"));
		if(benchmark) {
			System.err.println("Benchmark - SIMD="+Shared.SIMD);
		}
		assert(dstring.length()==1);
		
		final ArrayList<byte[]> lines;
		if(new File(fname).exists()){
			lines=ByteFile.toLines(fname);
		}else{
			lines=new ArrayList<byte[]>(1);
			lines.add(fname.getBytes());
		}
		Timer t=new Timer();
		long bytes=0, terms=0;
		LineParser1 lp=new LineParser1(dstring.charAt(0));
		for(byte[] line : lines) {
			lp.set(line);
			bytes+=line.length;
			terms+=lp.terms();
			if(!benchmark) {System.out.println(lp);}
		}
		t.stop();
		System.err.println(Tools.timeLinesBytesProcessed(t, lines.size(), bytes, 8));
		System.err.println(Tools.thingsProcessed(t.elapsed, terms, 8, "Terms"));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public LineParser1(byte delimiter_) {delimiter=delimiter_;}

	public LineParser1(int delimiter_) {
		assert(delimiter_>=0 && delimiter_<=127);
		delimiter=(byte)delimiter_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public LineParser1 set(byte[] line_) {
		clear();
		line=line_;
//		for(int len=advance(); b<line.length; len=advance()) {
//			bounds.add(b);
//		}
//		bounds.add(b);
		Vector.findSymbols(line, 0, line.length, delimiter, bounds);
		bounds.add(line.length);
		b=bounds.get(0);
		return this;
	}
	
	@Override
	public LineParser set(byte[] line_, int maxTerm) {
		clear();
		line=line_;
		for(int term=0; term<=maxTerm; term++) {
			int len=advance();
			bounds.add(b);
		}
		return this;
	}
	
	@Override
	public LineParser clear() {
		line=null;
		a=b=-1;
		bounds.clear();
		return this;
	}
	
	/**
	 * Reset method - no operation for this implementation.
	 * Included for interface compatibility.
	 * @return this LineParser instance for method chaining
	 */
	@Override
	public LineParser reset() {
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
	 * @return The parsed integer value
	 */
	@Override
	public int parseInt(int term) {
		setBounds(term);
		return Parse.parseInt(line, a, b);
	}
	
	public int parseInt(int term, int offset) {
		setBounds(term);
		return Parse.parseInt(line, a+offset, b);
	}
	
	/**
	 * Parses the specified term as a long integer.
	 * @param term Zero-based term index
	 * @return The parsed long value
	 */
	@Override
	public long parseLong(int term) {
		setBounds(term);
		return Parse.parseLong(line, a, b);
	}
	
	public long parseLongA48(int term) {
		setBounds(term);
		return Parse.parseLongA48(line, a, b);
	}
	
	public long[] parseLongArray(int term) {
		long[] array=new long[terms()-term];
		return parseLongArray(term, array);
	}
	
	public long[] parseLongArray(int term, long[] array) {
		for(int i=0; i<array.length; i++) {
			array[i]=parseLong(term+i);
		}
		return array;
	}
	
	public long[] parseLongArrayA48(int term, long[] array) {
		for(int i=0; i<array.length; i++) {
			array[i]=parseLongA48(term+i);
		}
		return array;
	}
	
	/**
	 * Parses the specified term as a float.
	 * @param term Zero-based term index
	 * @return The parsed float value
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
	 * Gets the byte at the specified offset within a term.
	 * @param term Zero-based term index
	 * @param offset Character offset within the term
	 * @return The byte value at the position
	 */
	@Override
	public byte parseByte(int term, int offset) {
		setBounds(term);
		final int index=a+offset;
		assert(index<b);
		return line[index];
	}
	
	public byte parseByteFromCurrentField(int offset) {
		assert(a<b);
		return line[a];
	}
	
	@Override
	public byte[] parseByteArray(int term) {
		final int len=setBounds(term);
		return Arrays.copyOfRange(line, a, b);
	}
	
	public byte[] parseByteArray(int term, int offset) {
		final int len=setBounds(term);
		return Arrays.copyOfRange(line, a+offset, b);
	}
	
	/**
	 * Gets the current field as a byte array copy.
	 * Uses the current field bounds set by setBounds().
	 * @return Copy of current field as byte array
	 */
	@Override
	public byte[] parseByteArrayFromCurrentField() {
		return Arrays.copyOfRange(line, a, b);
	}
	
	@Override
	public String parseString(int term) {
		final int len=setBounds(term);
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
	
	/*--------------------------------------------------------------*/
	/*----------------         Query Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Tests if the line starts with the specified String.
	 * @param s The prefix string to test
	 * @return true if line starts with the string
	 */
	@Override
	public boolean startsWith(String s) {
		return Tools.startsWith(line, s);
	}
	
	/**
	 * Tests if the line starts with the specified character.
	 * @param c The prefix character to test
	 * @return true if line starts with the character
	 */
	@Override
	public boolean startsWith(char c) {
		return Tools.startsWith(line, c);
	}
	
	/**
	 * Tests if the line starts with the specified byte.
	 * @param b The prefix byte to test
	 * @return true if line starts with the byte
	 */
	@Override
	public boolean startsWith(byte b) {
		return Tools.startsWith(line, b);
	}
	
	/**
	 * Tests if the specified term starts with the given string.
	 * @param s The prefix string to test
	 * @param term Zero-based term index
	 * @return true if term starts with the string
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
	 * Tests if the specified term exactly equals the given string.
	 * @param s The string to compare
	 * @param term Zero-based term index
	 * @return true if term equals the string
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
	 * Tests if the specified term exactly equals the given character.
	 * @param c The character to compare
	 * @param term Zero-based term index
	 * @return true if term equals the character
	 */
	@Override
	public boolean termEquals(char c, int term) {
		final int len=setBounds(term);
		return len==1 && line[a]==c;
	}
	
	/**
	 * Tests if the specified term exactly equals the given byte.
	 * @param c The byte to compare
	 * @param term Zero-based term index
	 * @return true if term equals the byte
	 */
	@Override
	public boolean termEquals(byte c, int term) {
		final int len=setBounds(term);
		return len==1 && line[a]==c;
	}
	
	public boolean currentTermEquals(byte c) {
		return b-a==1 && line[a]==c;
	}

	/**
	 * Gets the length of the specified term in characters.
	 * @param term Zero-based term index
	 * @return Length of the term
	 */
	@Override
	public int length(int term) {
		return setBounds(term);
	}

	/**
	 * Gets the length of the current field.
	 * Uses current field bounds set by setBounds().
	 * @return Length of current field in characters
	 */
	@Override
	public int currentFieldLength() {
		return b-a;
	}
	
	/**
	 * Increments the start position of current field by the specified amount.
	 * @param amt Amount to increment start position
	 * @return New length of current field
	 */
	@Override
	public int incrementA(int amt) {
		a+=amt;
		return b-a;
	}
	
	/**
	 * Increments the start position by the specified amount.
	 * Note: Implementation appears to increment 'a' instead of 'b'.
	 * @param amt Amount to increment position
	 * @return New length of current field
	 */
	@Override
	public int incrementB(int amt) {
		a+=amt;
		return b-a;
	}

	/** Tests if there are more characters to process in the line.
	 * @return true if current position is before end of line */
	@Override
	public boolean hasMore() {
		return b<line.length;
	}

	/** Gets the total length of the current line.
	 * @return Length of line in characters */
	@Override
	public int lineLength() {
		return line.length;
	}

	/** Gets reference to the current line byte array.
	 * @return The current line as byte array */
	@Override
	public byte[] line() {return line;}
	
	/** Gets the current field start position.
	 * @return Start index of current field */
	@Override
	public int a() {return a;}
	
	/** Gets the current field end position.
	 * @return End index of current field */
	@Override
	public int b() {return b;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets the current field boundaries for the specified term.
	 * Updates internal start (a) and end (b) positions from pre-computed bounds.
	 * @param term Zero-based term index
	 * @return Length of the term
	 */
	@Override
	public int setBounds(int term){
		a=(term==0 ? 0 : bounds.get(term-1)+1);
		b=bounds.get(term);
		return b-a;
	}
	
	/**
	 * Advances to the next delimiter position during bounds computation.
	 * Internal method used during line parsing to build the bounds list.
	 * Should not be made public as it's for internal bounds construction only.
	 * @return Length of the found field
	 */
	private int advance() {
		b++;
		a=b;
		while(b<line.length && line[b]!=delimiter){b++;}
		return b-a;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns a string representation of all parsed terms.
	 * @return String representation of the parsed terms */
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
	private byte[] line;
	
	public final byte delimiter;
}
