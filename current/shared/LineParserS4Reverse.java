package shared;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextFile;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Uses multiple ordered delimiters, e.g. ",. ,,".
 * Parses the line right-to-left, but the delimiter string
 * and public functions are still left-to-right.
 * For situations where the end of the line is known but the
 * prefix is unknown.
 * 
 * @author Brian Bushnell
 * @date May 6, 2024
 *
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
	
	@Override
	public LineParserS4Reverse set(String line_, int maxTerm) {
		throw new RuntimeException("Not valid.");
	}
	
	@Override
	public LineParserS4Reverse clear() {
		delimiterPos=0;
		line=null;
		a=b=-1;
		bounds.clear();
		return this;
	}
	
	@Override
	public LineParserS4Reverse reset() {
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
	
	public int parseInt(int term, int from, int to) {
		setBounds(term);
		return Parse.parseInt(line, a+from, Tools.min(line.length(), a+to));
//		return Parse.parseInt(line, a+from, Tools.min(b, a+to));
	}
	
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
	
	@Override
	public double parseDouble(int term) {
		setBounds(term);
		return Parse.parseDouble(line, a, b);
	}
	
	@Override
	public byte parseByte(int term, int offset) {
		return (byte)parseChar(term, offset);
	}
	
	@Override
	public char parseChar(int term, int offset) {
		setBounds(term);
		final int index=a+offset;
		assert(index<b);
		return line.charAt(index);
	}
	
	@Override
	public byte[] parseByteArray(int term) {
		final int len=setBounds(term);
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
		final int len=setBounds(term);
		return a>=b ? null : line.substring(a, b);
	}

	@Override
	public ByteBuilder appendTerm(ByteBuilder bb, int term) {
		final int len=setBounds(term);
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
		final int len=setBounds(term);
		if(len<s.length()) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(c!=line.charAt(a+i)) {return false;}
		}
		return true;
	}
	
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
		return setBounds(term);
	}

	@Override
	public int currentFieldLength() {
		return b-a;
	}

	@Override
	public boolean hasMore() {
		return a>=0;
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
	private String line;
	
	//This is already reversed
	public final String delimiters;
	private final int maxDPos;
	private int delimiterPos=0;
	
}
