package shared;

import java.io.File;
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

	public LineParser2(byte delimiter_) {delimiter=delimiter_;}

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
		return line[index];
	}
	
	public String parseString() {
		int len=advance();
		return new String(line, a, len);
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
	
	public final int advance() {
		currentTerm++;
		b++;
		a=b;
		while(b<line.length && line[b]!=delimiter){b++;}
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
	
	private int a=-1;
	private int b=-1;
	private int currentTerm=-1;
	private byte[] line;
	
	public final byte delimiter;
	
}
