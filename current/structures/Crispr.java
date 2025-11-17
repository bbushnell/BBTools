package structures;

import repeat.Palindrome;
import shared.Tools;

/**
 * Represents a CRISPR (Clustered Regularly Interspaced Short Palindromic Repeats) structure
 * with two distinct genomic ranges and associated palindrome information.
 * Provides methods for calculating gaps, overlaps, scoring, and structural analysis.
 * @author Brian Bushnell
 */
public class Crispr implements Comparable<Crispr> {
	
	public Crispr() {}
	
	/**
	 * Constructs a CRISPR structure with two ranges.
	 *
	 * @param a1 Start position of first range
	 * @param b1 End position of first range
	 * @param a2 Start position of second range
	 * @param b2 End position of second range
	 */
	public Crispr(int a1, int b1, int a2, int b2) {
		a=new Range(a1, b1);
		b=new Range(a2, b2);
	}
	
	/** Returns the size of the gap between the two ranges */
	public int gap() {
		return b.a-a.b-1;
	}
	
	/**
	 * Determines if a range is completely contained within the gap between ranges a and b.
	 * @param r Range to test for containment
	 * @return true if the range is entirely within the gap
	 */
	public boolean containsInGap(Range r) {
		return r.a>a.b && r.b<b.a;
	}
	
	public String toString() {
		return appendTo(new ByteBuilder()).toString();
	}
	
	/**
	 * Returns a string representation including the actual sequence bases.
	 * @param bases Sequence array to extract bases from
	 * @return Formatted string with ranges and sequence content
	 */
	public String toString(byte[] bases) {
		ByteBuilder bb=new ByteBuilder();
		appendTo(bb, bases.length, bases);
		return bb.toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb) {
		return appendTo(bb, 0, null);
	}
	
	/**
	 * Appends a formatted representation of this CRISPR to a ByteBuilder.
	 * Includes range coordinates, optional length, palindrome info, and sequences.
	 *
	 * @param bb ByteBuilder to append to
	 * @param len Optional length parameter (0 to omit)
	 * @param bases Optional sequence array for displaying actual bases
	 * @return The ByteBuilder for method chaining
	 */
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
	
	/**
	 * Sets the coordinates of both ranges with validation.
	 * Ensures proper ordering constraints are met.
	 *
	 * @param a1 Start position of first range
	 * @param b1 End position of first range
	 * @param a2 Start position of second range
	 * @param b2 End position of second range
	 */
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
	
	/**
	 * Calculates the minimum distance to either edge of the sequence.
	 * @param length Total length of the sequence
	 * @return Distance to the nearest sequence edge
	 */
	public int edgeDist(int length) {
		return Tools.min(a.a, length-b.b-1);
	}
	
	/**
	 * Determines if this CRISPR spans the entire sequence from start to end.
	 * @param length Total sequence length
	 * @return true if CRISPR covers the full sequence
	 */
	public boolean spans(int length) {
		return a.a==0 && b.b+1==length;
	}
	
	/** Returns the length of the shorter of the two ranges */
	public int minLength() {
		return Tools.min(a.length(), b.length());
	}
	
	/** Returns the length of the longer of the two ranges */
	public int maxLength() {
		return Tools.max(a.length(), b.length());
	}
	
	/** Returns the difference in length between range b and range a */
	public int lengthDif() {
		return b.length()-a.length();
	}
	
	/** Ensures both ranges are within valid sequence bounds.
	 * @param length Total sequence length for boundary checking */
	public void fixBounds(int length) {
		a.fixBounds(length);
		b.fixBounds(length);
	}
	
	/** Returns the higher of the two range scores */
	public float maxScore() {
		return Tools.max(scoreA, scoreB);
	}
	
	/** Returns true if both ranges have the same length */
	public boolean sameLength() {
		return a.length()==b.length();
	}

//	public boolean internal(byte[] seq) {
//		return !touchesEdge(seq.length);
//	}
	
	/**
	 * Determines if this CRISPR is completely internal to the sequence.
	 * @param seqLen Total sequence length
	 * @return true if CRISPR doesn't touch either sequence end
	 */
	public boolean internal(int seqLen) {
		return a.a>0 && b.b<seqLen-1;
	}
	
//	public boolean touchesEdge(byte[] seq) {
//		return touchesEdge(seq.length);
//	}
	
	/**
	 * Determines if this CRISPR touches either end of the sequence.
	 * @param seqLen Total sequence length
	 * @return true if CRISPR touches start or end of sequence
	 */
	public boolean touchesEdge(int seqLen) {
		return a.a<=0 || b.b>=seqLen-1;
	}
	
	/**
	 * Determines if this CRISPR spans from the start to the end of the sequence.
	 * @param seqLen Total sequence length
	 * @return true if CRISPR touches both sequence ends
	 */
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
	
	/** First genomic range of the CRISPR structure */
	public Range a;
	/** Second genomic range of the CRISPR structure */
	public Range b;
	public Palindrome pa, pb;
	public float scoreA, scoreB;
	public int matches=0, mismatches=0;
	/** Length of consensus sequence after trimming */
	public int trimmedConsensus=0;
	/** Length of consensus sequence after extension */
	public int extendedConsensus=0;
	
}
