package barcode;

import dna.AminoAcid;
import structures.ByteBuilder;

public class Transition implements Comparable<Transition> {
	
	public Transition(int pos_, byte ref_, byte query_, long count_) {
		pos=pos_;
		ref=ref_;
		query=query_;
		count=count_;
	}
	
	public int encode() {return encode(pos, ref, query);}
	
	public static int encode(int pos, int ref, int query) {
		int x1=baseToNumber[ref];
		assert(x1>=0 && x1<4);//Only defined symbols allowed for ref
		int x2=baseToNumber[query];
		int idx=((pos<<2)|x1)*5+x2;
		return idx;
	}
	
	public static Transition decode(int idx) {
		int x2=idx%5;
		idx/=5;
		int x1=idx&3;
		idx=idx>>2;
		int pos=idx;
		byte r=numberToBase[x1];
		byte q=numberToBase[x2];
		return new Transition(pos, r, q, 0);
	}
	
	public ByteBuilder appendTo(ByteBuilder bb) {
		return bb.append(pos).tab().append(ref).tab().append(query).tab().append(count);
	}
	
	@Override
	public int compareTo(Transition b) {
		if(count!=b.count) {return count<b.count ? 1 : -1;}
		if(pos!=b.pos) {return pos-b.pos;}
		if(ref!=b.ref) {return baseToNumber[ref]-baseToNumber[b.ref];}
		return baseToNumber[query]-baseToNumber[b.query];
	}
	
	public final int pos;
	public final byte ref;
	public final byte query;
	public long count;

	private static final byte[] numberToBase=AminoAcid.numberToBase;
	private static final byte[] baseToNumber=AminoAcid.baseToNumber4;
	
}
