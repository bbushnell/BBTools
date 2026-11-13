package barcode;

import dna.AminoAcid;
import structures.ByteBuilder;

/**
 * Represents a single nucleotide transition at a specific position, tracking the
 * reference base, query base, position, and occurrence count for barcode analysis.
 * Includes compact encode/decode helpers and comparison support.
 * Positions and base bytes are immutable; count is mutable and requires caller
 * synchronization if shared. Encoding normalizes reference/query letter case,
 * U to T, and other ASCII query symbols to N; it does not preserve arbitrary symbols.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class Transition implements Comparable<Transition>{
	
	/** Retains the supplied values without validation or base normalization.
	 * @param pos_ Zero-based position
	 * @param ref_ Reference ASCII base
	 * @param query_ Query ASCII base or ambiguous symbol
	 * @param count_ Number of observed transitions */
	public Transition(int pos_, byte ref_, byte query_, long count_){
		pos=pos_;
		ref=ref_;
		query=query_;
		count=count_;
	}
	
	/** Encodes position and normalized bases; the observation count is not encoded. */
	public int encode(){return encode(pos, ref, query);}
	
	/**
	 * Encodes a transition into a compact integer representation using
	 * {@code ((pos << 2) | ref_number) * 5 + query_number}. Reference A/C/G/T
	 * letters and U are accepted case-insensitively. Other ASCII query symbols map
	 * to code4 (N). Position must keep the encoded value in the nonnegative int range.
	 *
	 * @param pos Nonnegative zero-based position
	 * @param ref ASCII reference base with a defined A/C/G/T/U code
	 * @param query ASCII query symbol, normalized through baseToNumber4
	 * @return Encoded integer representing the transition
	 */
	public static int encode(int pos, int ref, int query){
		int x1=baseToNumber[ref];
		assert(x1>=0 && x1<4);//Only defined symbols allowed for ref
		int x2=baseToNumber[query];
		int idx=((pos<<2)|x1)*5+x2;
		return idx;
	}
	
	/**
	 * Decodes an encoded transition back to a Transition object by reversing the
	 * encode formula. The returned Transition has count initialized to 0 and
	 * canonical uppercase bases; ambiguous query symbols decode as N.
	 * @param idx Nonnegative transition encoding produced without integer overflow
	 * @return Decoded Transition object with count set to 0
	 */
	public static Transition decode(int idx){
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
	public ByteBuilder appendTo(ByteBuilder bb){
		return bb.append(pos).tab().append(ref).tab().append(query).tab().append(count);
	}
	
	/**
	 * Compares transitions for sorting: count descending, then position, then
	 * reference base, then query base.
	 * @param b Transition to compare against
	 * @return Negative when this sorts first, positive when b sorts first, or zero for equal sort keys
	 */
	@Override
	public int compareTo(Transition b){
		if(count!=b.count){return count<b.count ? 1 : -1;}
		if(pos!=b.pos){return pos-b.pos;}
		if(ref!=b.ref){return baseToNumber[ref]-baseToNumber[b.ref];}
		return baseToNumber[query]-baseToNumber[b.query];
	}
	
	/** Zero-based position supplied at construction. */
	public final int pos;
	/** Original reference symbol; encode/decode may normalize it. */
	public final byte ref;
	/** Original query symbol; encode/decode may normalize it. */
	public final byte query;
	/** Mutable observation count, used as the primary descending sort key. */
	public long count;

	private static final byte[] numberToBase=AminoAcid.numberToBase;
	private static final byte[] baseToNumber=AminoAcid.baseToNumber4;
	
}
