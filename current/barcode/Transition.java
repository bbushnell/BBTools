package barcode;

import dna.AminoAcid;
import structures.ByteBuilder;

/**
 * Represents a single nucleotide transition at a specific position.
 * Tracks reference base, query base, position, and occurrence count for barcode analysis.
 * Provides encoding/decoding functionality for compact storage and comparison capabilities.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class Transition implements Comparable<Transition> {
	
	/**
	 * Constructs a transition with specified position, bases, and count.
	 *
	 * @param pos_ Position in the sequence where transition occurs
	 * @param ref_ Reference base at this position
	 * @param query_ Query base observed at this position
	 * @param count_ Number of times this transition was observed
	 */
	public Transition(int pos_, byte ref_, byte query_, long count_) {
		pos=pos_;
		ref=ref_;
		query=query_;
		count=count_;
	}
	
	/** Encodes this transition into a compact integer representation.
	 * @return Encoded integer representing position and base changes */
	public int encode() {return encode(pos, ref, query);}
	
	/**
	 * Encodes a transition into a compact integer representation.
	 * Uses formula: ((pos << 2) | ref_number) * 5 + query_number.
	 * Reference base must be a defined nucleotide (A, C, G, T).
	 *
	 * @param pos Position in the sequence
	 * @param ref Reference base (must be defined nucleotide)
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
	 * Decodes an encoded transition back to a Transition object.
	 * Reverses the encoding process to extract position and bases.
	 * Creates transition with count=0.
	 *
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
	 * Appends transition data to a ByteBuilder in tab-separated format.
	 * Format: position, reference base, query base, count.
	 * @param bb ByteBuilder to append to
	 * @return The same ByteBuilder for method chaining
	 */
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
	
	/** Position in the sequence where this transition occurs */
	public final int pos;
	/** Reference base at this position */
	public final byte ref;
	/** Query base observed at this position */
	public final byte query;
	/** Number of times this transition was observed */
	public long count;

	/** Lookup table for converting base numbers back to nucleotide bytes */
	private static final byte[] numberToBase=AminoAcid.numberToBase;
	/**
	 * Lookup table for converting nucleotide bytes to base numbers (0-3 for A,C,G,T)
	 */
	private static final byte[] baseToNumber=AminoAcid.baseToNumber4;
	
}
