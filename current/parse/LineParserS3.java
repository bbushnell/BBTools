package parse;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextFile;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Line parser that treats tab and space as implicit delimiters in addition to
 * a specified delimiter.
 * Extends standard delimiter-based parsing by automatically recognizing whitespace boundaries,
 * making it suitable for mixed-delimiter text formats. Unlike LineParserS1, this parser considers
 * any occurrence of tab or space as additional delimiter characters regardless
 * of the primary delimiter.
 *
 * <p>Indexes every primary delimiter, space and tab separately, including adjacent
 * empty fields. Inputs are retained; auxiliary memory grows with field count.
 * maxTerm is an inclusive existing field index. reset() is a no-op. Mutable
 * instances and the exposed bounds list belong to one thread.
 *
 * @author Brian Bushnell
 * @date May 24, 2023
 */
public final class LineParserS3 implements LineParserS{

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

		LineParserS lp=new LineParserS3(dstring.charAt(0));
		for(String line : lines){
			lp.set(line);
			System.out.println(lp);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a parser with the given UTF-16 delimiter code unit. */
	public LineParserS3(char delimiter_){delimiter=delimiter_;}

	/** Creates a character-delimited parser; asserts the inclusive UTF-16 range 0..65535. */
	public LineParserS3(int delimiter_){
		assert(delimiter_>=0 && delimiter_<=Character.MAX_VALUE);
		delimiter=(char)delimiter_;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public LineParserS3 set(byte[] line_){
		assert(false) : "Use byte version.";
		return set(new String(line_));
	}

	@Override
	public LineParserS3 set(byte[] line_, int maxTerm){
		assert(false) : "Use byte version.";
		return set(new String(line_), maxTerm);
	}

	@Override
	public LineParserS3 set(String line_){
		clear();
		line=line_;
		for(int len=advance(); b<line.length(); len=advance()){
			bounds.add(b);
		}
		bounds.add(b);
		return this;
	}

	@Override
	public LineParserS3 set(String line_, int maxTerm){
		clear();
		line=line_;
		//TODO: test performance of presumably safer loop below
//		for(int term=0; term<=maxTerm && b<line.length(); term++) {
		for(int term=0; term<=maxTerm; term++){
			int len=advance();
			bounds.add(b);
		}
		return this;
	}

	@Override
	public LineParserS3 clear(){
		line=null;
		a=b=-1;
		bounds.clear();
		whitespaceIndex=-1;
		return this;
	}

	@Override
	public LineParserS3 reset(){
		//Does nothing for this class
		return this;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of indexed fields, including empty fields. */
	public int terms(){return bounds.size();}

	@Override
	public int parseInt(int term){
		setBounds(term);
		return Parse.parseInt(line, a, b);
	}

	@Override
	public long parseLong(int term){
		setBounds(term);
		return Parse.parseLong(line, a, b);
	}

	@Override
	public float parseFloat(int term){
		setBounds(term);
		return Parse.parseFloat(line, a, b);
	}

	@Override
	public double parseDouble(int term){
		setBounds(term);
		return Parse.parseDouble(line, a, b);
	}

	@Override
	public byte parseByte(int term, int offset){
		return (byte)parseChar(term, offset);
	}

	@Override
	public char parseChar(int term, int offset){
		setBounds(term);
		final int index=a+offset;
//		assert(false) : bounds()+"\n"
//				+ "line="+line+"\n"+
//				"a="+a+", b="+b+", index="+index+", term="+
//				+term+", offset="+offset+", len="+line.length();
		assert(index<b);
		return line.charAt(index);
	}

	@Override
	public byte[] parseByteArray(int term){
		final int len=setBounds(term);
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
		final int len=setBounds(term);
		return line.substring(a, b);
	}

	@Override
	public ByteBuilder appendTerm(ByteBuilder bb, int term){
		final int len=setBounds(term);
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

	@Override
	public boolean startsWith(char c){
		return Tools.startsWith(line, c);
	}

	@Override
	public boolean startsWith(byte b){
		return Tools.startsWith(line, b);
	}

	@Override
	public boolean termStartsWith(String s, int term){
		final int len=setBounds(term);
		if(len<s.length()){return false;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c!=line.charAt(a+i)){return false;}
		}
		return true;
	}

	@Override
	public boolean termEquals(String s, int term){
		final int len=setBounds(term);
		if(len!=s.length()){return false;}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c!=line.charAt(a+i)){return false;}
		}
		return true;
	}

	@Override
	public boolean termEquals(char c, int term){
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}

	@Override
	public boolean termEquals(byte c, int term){
		final int len=setBounds(term);
		return len==1 && line.charAt(a)==c;
	}

	@Override
	public int incrementA(int amt){
		a+=amt;
		return b-a;
	}

	@Override
	public int incrementB(int amt){
		b+=amt;
		return b-a;
	}

	@Override
	public int length(int term){
		return setBounds(term);
	}

	@Override
	public int currentFieldLength(){
		return b-a;
	}

	@Override
	public boolean hasMore(){
		return b<line.length();
	}

	@Override
	public int lineLength(){
		return line.length();
	}

	@Override
	public String line(){return line;}

	@Override
	public int a(){return a;}

	@Override
	public int b(){return b;}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int setBounds(int term){
		a=(term==0 ? 0 : bounds.get(term-1)+1);
		b=bounds.get(term);
		return b-a;
	}

	/** Scans one field and records its first non-primary whitespace boundary, before bounds insertion. */
	private int advance(){
		b++;
		a=b;
		while(b<line.length()){
			final char c=line.charAt(b);
			if(c==delimiter){break;}
			if(c==' ' || c=='\t'){
				// Primary space/tab delimiters have already been handled above.
				if(whitespaceIndex<0){whitespaceIndex=bounds.size;}
				break;
			}
			b++;
		}
		return b-a;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public String toString(){
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

	/** Index of the first parsed field terminated by space/tab rather than the primary
	 * delimiter, or -1 when absent. With a limited setter only scanned boundaries count.
	 * A space/tab that is itself the primary delimiter does not count as extra whitespace. */
	public int indexOfWhitespace(){return whitespaceIndex;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the mutable backing delimiter-position list; caller mutation can invalidate parsing. */
	public IntList bounds(){return bounds;}
	/** Exclusive field ends, reused across lines; the final stored end is normally the line length. */
	private final IntList bounds=new IntList();
	/** First non-primary whitespace boundary; reset for each line and tracked while scanning. */
	private int whitespaceIndex=-1;

	/** Start of the current field in native line units, inclusive; -1 before initialization. */
	private int a=-1;
	/** End of the current field in native line units, exclusive; -1 before initialization. */
	private int b=-1;
	/** Retained input string; null after clear and before set. */
	private String line;

	/** Primary UTF-16 delimiter code unit. */
	public final char delimiter;

}
