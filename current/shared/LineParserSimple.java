package shared;

public class LineParserSimple {

	public LineParserSimple(byte delimiter_) {
		delimiter=delimiter_;
	}

	public LineParserSimple(byte delimiter_, byte[] line_) {
		delimiter=delimiter_;
		line=line_;
	}
	
	public int advanceInner() {//Does not check array bounds
		b++;
		a=b;
		assert(b<line.length);
		while(line[b]!=delimiter){b++;}
		return b-a;
	}
	
	public int advance() {
		b++;
		a=b;
		while(b<line.length && line[b]!=delimiter){b++;}
		return b-a;
	}
	
	public void reset() {
		a=b=segment=-1;
	}
	
	private int a=-1;
	private int b=-1;
	private int segment=-1;
	private byte[] line;
	
	private final byte delimiter;
}
