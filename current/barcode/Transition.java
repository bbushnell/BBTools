package barcode;

import dna.AminoAcid;
import structures.ByteBuilder;

/**
 * Represents a single nucleotide transition at a specific position, tracking the
 * reference base, query base, position, and occurrence count for barcode analysis.
 * Includes compact encode/decode helpers and comparison support.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class Transition implements Comparable<Transition> {
	
	public Transition(int pos_, byte ref_, byte query_, long count_) {
		pos=pos_;
		ref=ref_;
		query=query_;
		count=count_;
	}
	
	public int encode() {return encode(pos, ref, query);}
	
	/**
	 * Encodes a transition into a compact integer representation using
	 * `((pos << 2) | ref_number) * 5 + query_number`. The reference base must be
	 * one of A/C/G/T.
	 *
	 * @param pos Position in the sequence
	 * @param ref Reference base (must be a defined nucleotide)
	 * @param query Query base (may be ambiguous)
	 * @return Encoded integer representing the transition
	 */
	public static int encode(int pos, int ref, int query) {
		int x1=baseToNumber[ref];
		assert(x1>=0 && x1<4);//Only defined symbols allowed for ref
		int x2=baseToNumber[query];
		int idx=((pos<<2)|x1)*5+x2;
		return idx;
	}
	
	/**
	 * Decodes an encoded transition back to a Transition object by reversing the
	 * encode formula. The returned Transition has count initialized to 0.
	 * @param idx Encoded transition integer
	 * @return Decoded Transition object with count set to 0
	 */
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
	
	/**
	 * Appends transition data to a ByteBuilder in tab-separated format:
	 * position, reference base, query base, count.
	 * @param bb ByteBuilder to append to
	 * @return The same ByteBuilder for method chaining
	 */
	public ByteBuilder appendTo(ByteBuilder bb) {
		return bb.append(pos).tab().append(ref).tab().append(query).tab().append(count);
	}
	
	/**
	 * Compares transitions for sorting: count descending, then position, then
	 * reference base, then query base.
	 * @param b Transition to compare against
	 * @return Negative if this < b, positive if this > b, zero if equal
	 */
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
