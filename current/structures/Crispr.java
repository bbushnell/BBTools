package structures;

import repeat.Palindrome;
import shared.Tools;

/**
 * Mutable pair of inclusive repeat ranges in one sequence, with optional
 * palindrome annotations, scores and caller-maintained alignment counters.
 * Coordinates normally refer to a before b; the coordinate constructor checks
 * each range independently and does not enforce separation between them.
 * All fields are public and unsynchronized. Callers initialize both ranges
 * before querying them and coordinate mutation, sorting and annotation updates.
 * Ordering compares endpoints only; object equality remains identity-based.
 * @author Brian Bushnell
 */
public class Crispr implements Comparable<Crispr>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an uninitialized pair with null ranges and annotations. */
	public Crispr(){}

	/**
	 * Allocates two ranges, each with assertion-checked ordered endpoints.
	 *
	 * @param a1 Start position of first range
	 * @param b1 End position of first range
	 * @param a2 Start position of second range
	 * @param b2 End position of second range
	 */
	public Crispr(int a1, int b1, int a2, int b2){
		a=new Range(a1, b1);
		b=new Range(a2, b2);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns b.start-a.end-1: zero for adjacency and negative for overlap. */
	public int gap(){
		return b.a-a.b-1;
	}

	/**
	 * Determines if a range is completely contained within the gap between ranges a and b.
	 * @param r Range to test for containment
	 * @return true if the range is entirely within the gap
	 */
	public boolean containsInGap(Range r){
		return r.a>a.b && r.b<b.a;
	}

	/** Formats coordinates and the preferred available palindrome annotation. */
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}

	/**
	 * Returns a string representation including the actual sequence bases.
	 * @param bases Nonnull sequence array containing both inclusive ranges
	 * @return Formatted string with ranges and sequence content
	 */
	public String toString(byte[] bases){
		ByteBuilder bb=new ByteBuilder();
		appendTo(bb, bases.length, bases);
		return bb.toString();
	}

	/** Appends coordinates and one available palindrome, without length or bases. */
	public ByteBuilder appendTo(ByteBuilder bb){
		return appendTo(bb, 0, null);
	}

	/**
	 * Appends a formatted representation of this CRISPR to a ByteBuilder.
	 * Includes coordinates, positive len, and at most one palindrome: pa if present,
	 * otherwise pb. The historical delimiter difference between these paths is retained.
	 * Optional bases are appended as two newline-prefixed inclusive sequence slices.
	 *
	 * @param bb ByteBuilder to append to
	 * @param len Displayed sequence length; nonpositive values omit it
	 * @param bases Optional sequence array; both ranges must be valid array indexes
	 * @return The ByteBuilder for method chaining
	 */
	public ByteBuilder appendTo(ByteBuilder bb, int len, byte[] bases){
		bb.append('[').append(a.a).dash().append(a.b).comma();
		bb.append(b.a).dash().append(b.b);
		if(len>0){bb.semi().append(len);}
		bb.append(']');
		if(pa!=null){pa.appendTo(bb.comma(), a.a, a.b);}
		// The second palindrome belongs to b; its tail distances use b's bounds.
		else if(pb!=null){pb.appendTo(bb, b.a, b.b);}
		if(bases!=null){
			bb.nl();
			for(int i=a.a; i<=a.b; i++){bb.append(bases[i]);}
			bb.nl();
			for(int i=b.a; i<=b.b; i++){bb.append(bases[i]);}
		}
		return bb;
	}

	/**
	 * Mutates both existing ranges, then asserts separation and at least two bases
	 * per range. Assertions do not roll back those assignments on failure.
	 * Scores, palindromes and counters are not recomputed here.
	 *
	 * @param a1 Start position of first range
	 * @param b1 End position of first range
	 * @param a2 Start position of second range
	 * @param b2 End position of second range
	 */
	public void set(int a1, int b1, int a2, int b2){
		a.a=a1;
		a.b=b1;
		b.a=a2;
		b.b=b2;

		//These could happen and will probably be dealt with later, but are good to prevent if possible.
		assert(a2>b1) : this;
		assert(a1<b1) : this;
		assert(a2<b2) : this;
	}

	/**
	 * Calculates the smaller distance from a.start or b.end to its sequence edge.
	 * The result can be negative when either endpoint lies outside the sequence.
	 * @param length Total length of the sequence
	 * @return Distance to the nearest sequence edge
	 */
	public int edgeDist(int length){
		return Tools.min(a.a, length-b.b-1);
	}

	/**
	 * Tests exact placement of a.start at zero and b.end at length-1.
	 * A spacer between the two ranges is permitted.
	 * @param length Total sequence length
	 * @return true if both outer endpoints exactly match the sequence ends
	 */
	public boolean spans(int length){
		return a.a==0 && b.b+1==length;
	}

	/** Returns the length of the shorter of the two ranges */
	public int minLength(){
		return Tools.min(a.length(), b.length());
	}

	/** Returns the length of the longer of the two ranges */
	public int maxLength(){
		return Tools.max(a.length(), b.length());
	}

	/** Returns the difference in length between range b and range a */
	public int lengthDif(){
		return b.length()-a.length();
	}

	/**
	 * Applies Range's one-sided clamping to both ranges; wholly outside ranges may
	 * become empty or reversed. Does not recompute annotations, scores or counters.
	 * @param length Sequence length passed to each range's fixBounds method
	 */
	public void fixBounds(int length){
		a.fixBounds(length);
		b.fixBounds(length);
	}

	/** Returns the higher of the two range scores */
	public float maxScore(){
		return Tools.max(scoreA, scoreB);
	}

	/** Returns true if both ranges have the same length */
	public boolean sameLength(){
		return a.length()==b.length();
	}

//	public boolean internal(byte[] seq) {
//		return !touchesEdge(seq.length);
//	}

	/**
	 * Determines if this CRISPR is completely internal to the sequence.
	 * @param seqLen Total sequence length
	 * @return true if a.start is positive and b.end is below seqLen-1
	 */
	public boolean internal(int seqLen){
		return a.a>0 && b.b<seqLen-1;
	}

//	public boolean touchesEdge(byte[] seq) {
//		return touchesEdge(seq.length);
//	}

	/**
	 * Tests whether either outer endpoint reaches or extends beyond its sequence edge.
	 * @param seqLen Total sequence length
	 * @return true if a.start is at most zero or b.end is at least seqLen-1
	 */
	public boolean touchesEdge(int seqLen){
		return a.a<=0 || b.b>=seqLen-1;
	}

	/**
	 * Tests whether both outer endpoints reach or extend beyond their sequence edges.
	 * @param seqLen Total sequence length
	 * @return true if a.start is at most zero and b.end is at least seqLen-1
	 */
	public boolean touchesBothEnds(int seqLen){
		return a.a<=0 && b.b>=seqLen-1;
	}

	//Orders by first range, then second range; not a strict total order over all fields but adequate for sorting.
	@Override
	public int compareTo(Crispr o){
		int x=a.compareTo(o.a);
		if(x!=0){return x;}
		return b.compareTo(o.b);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** First genomic range of the CRISPR structure */
	public Range a;
	/** Second genomic range of the CRISPR structure */
	public Range b;
	/** Optional palindrome annotations for a and b, maintained by the caller. */
	public Palindrome pa, pb;
	/** Caller-maintained scores for the respective ranges. */
	public float scoreA, scoreB;
	/** Caller-maintained alignment counts; matches=-1 can request recomputation. */
	public int matches=0, mismatches=0;
	/** Accumulated consensus-trimming amount, not the resulting sequence length. */
	public int trimmedConsensus=0;
	/** Accumulated consensus-extension amount, not the resulting sequence length. */
	public int extendedConsensus=0;

}
