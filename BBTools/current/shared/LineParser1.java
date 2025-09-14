package shared;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Finds delimiters of a text line efficiently, to allow for parsing.
 * For example:<br>
 * Integer.parseInt("a b c 22 jan".split(" ")[3])<br>
 * 
 * could be redone as:<br>
 * LineParser lp=new LineParser(' ')<br>
 * lp.set("a b c 22 jan".toBytes()).parseInt(3)<br>
 * 
 * Uses memory proportional to 4*(# delimiters per line); for constant memory, use LineParser2.
 * 
 * @author Brian Bushnell
 * @date May 24, 2023
 *
 */
public final class LineParser1 implements LineParser {
	
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
		
		final ArrayList<byte[]> lines;
		if(new File(fname).exists()){
			lines=ByteFile.toLines(fname);
		}else{
			lines=new ArrayList<byte[]>(1);
			lines.add(fname.getBytes());
		}
		
		LineParser lp=new LineParser1(dstring.charAt(0));
		for(byte[] line : lines) {
			lp.set(line);
			System.out.println(lp);
		}
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
		for(int len=advance(); b<line.length; len=advance()) {
			bounds.add(b);
		}
		bounds.add(b);
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
	
	@Override
	public LineParser reset() {
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
	
	public int parseInt(int term, int offset) {
		setBounds(term);
		return Parse.parseInt(line, a+offset, b);
	}
	
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
	
	@Override
	public byte parseByte(int term, int offset) {
		setBounds(term);
		final int index=a+offset;
		assert(index<b);
		return line[index];
	}
	
	@Override
	public byte[] parseByteArray(int term) {
		final int len=setBounds(term);
		return Arrays.copyOfRange(line, a, b);
	}
	
	@Override
	public byte[] parseByteArrayFromCurrentField() {
		return Arrays.copyOfRange(line, a, b);
	}
	
	@Override
	public String parseString(int term) {
		final int len=setBounds(term);
		return new String(line, a, len);
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
		return new String(line, a, b-a);
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
		final int len=setBounds(term);
		if(len<s.length()) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(c!=line[a+i]) {return false;}
		}
		return true;
	}
	
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
		return setBounds(term);
	}

	@Override
	public int currentFieldLength() {
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

	@Override
	public boolean hasMore() {
		return b<line.length;
	}

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
	
	@Override
	public int setBounds(int term){
		a=(term==0 ? 0 : bounds.get(term-1)+1);
		b=bounds.get(term);
		return b-a;
	}
	
	/** 
	 * Do not make public.  This is for internal use making the bounds list,
	 * not for advancing to a new term like in LineParser2.
	 * @return Length of new field.
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
	
	@Override
	public String toString() {
		return toList().toString();
	}
	
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
