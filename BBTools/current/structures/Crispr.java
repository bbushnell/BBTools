package structures;

import repeat.Palindrome;
import shared.Tools;

public class Crispr implements Comparable<Crispr> {
	
	public Crispr() {}
	
	public Crispr(int a1, int b1, int a2, int b2) {
		a=new Range(a1, b1);
		b=new Range(a2, b2);
	}
	
	public int gap() {
		return b.a-a.b-1;
	}
	
	public boolean containsInGap(Range r) {
		return r.a>a.b && r.b<b.a;
	}
	
	public String toString() {
		return appendTo(new ByteBuilder()).toString();
	}
	
	public String toString(byte[] bases) {
		ByteBuilder bb=new ByteBuilder();
		appendTo(bb, bases.length, bases);
		return bb.toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb) {
		return appendTo(bb, 0, null);
	}
	
	public ByteBuilder appendTo(ByteBuilder bb, int len, byte[] bases) {
		bb.append('[').append(a.a).dash().append(a.b).comma();
		bb.append(b.a).dash().append(b.b);
		if(len>0) {bb.semi().append(len);}
		bb.append(']');
		if(pa!=null) {pa.appendTo(bb.comma(), a.a, a.b);}
		else if(pb!=null) {pb.appendTo(bb, a.a, a.b);}
		if(bases!=null) {
			bb.nl();
			for(int i=a.a; i<=a.b; i++) {bb.append(bases[i]);}
			bb.nl();
			for(int i=b.a; i<=b.b; i++) {bb.append(bases[i]);}
		}
		return bb;
	}
	
	public void set(int a1, int b1, int a2, int b2) {
		a.a=a1;
		a.b=b1;
		b.a=a2;
		b.b=b2;
		
		//These could happen and will probably be dealt with later, but are good to prevent if possible.
		assert(a2>b1) : this;
		assert(a1<b1) : this;
		assert(a2<b2) : this;
	}
	
	public int edgeDist(int length) {
		return Tools.min(a.a, length-b.b-1);
	}
	
	public boolean spans(int length) {
		return a.a==0 && b.b+1==length;
	}
	
	public int minLength() {
		return Tools.min(a.length(), b.length());
	}
	
	public int maxLength() {
		return Tools.max(a.length(), b.length());
	}
	
	public int lengthDif() {
		return b.length()-a.length();
	}
	
	public void fixBounds(int length) {
		a.fixBounds(length);
		b.fixBounds(length);
	}
	
	public float maxScore() {
		return Tools.max(scoreA, scoreB);
	}
	
	public boolean sameLength() {
		return a.length()==b.length();
	}

//	public boolean internal(byte[] seq) {
//		return !touchesEdge(seq.length);
//	}
	
	public boolean internal(int seqLen) {
		return a.a>0 && b.b<seqLen-1;
	}
	
//	public boolean touchesEdge(byte[] seq) {
//		return touchesEdge(seq.length);
//	}
	
	public boolean touchesEdge(int seqLen) {
		return a.a<=0 || b.b>=seqLen-1;
	}
	
	public boolean touchesBothEnds(int seqLen) {
		return a.a<=0 && b.b>=seqLen-1;
	}
	
	@Override
	/** 
	 * So this sort of compares them but it's not really optimal.
	 * Not clear how to improve it though.
	 * @param o
	 * @return
	 */
	public int compareTo(Crispr o) {
		int x=a.compareTo(o.a);
		if(x!=0) {return x;}
		return b.compareTo(o.b);
	}
	
	public Range a;
	public Range b;
	public Palindrome pa, pb;
	public float scoreA, scoreB;
	public int matches=0, mismatches=0;
	public int trimmedConsensus=0;
	public int extendedConsensus=0;
	
}
