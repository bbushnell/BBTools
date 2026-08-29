package assemble;

import simd.Vector;
import structures.ByteBuilder;

/**
 * Represents directed edges in the De Bruijn graph connecting two contigs
 * with orientation and overlap information for assembly operations.
 * Stores connection data including orientation encoding, overlap length,
 * coverage depth, and optional sequence data for graph traversal.
 *
 * @author Brian Bushnell
 */
public class Edge {
	
	public Edge(int origin_, int destination_, int length_, int orientation_, int depth_, byte[] bases_){
		this(origin_, destination_, length_, orientation_, depth_, bases_, 0);
	}

	/** Constructs an edge with an explicit overlap; zero uses the graph kmer length. */
	public Edge(int origin_, int destination_, int length_, int orientation_, int depth_, byte[] bases_, int overlap_){
		this(origin_, destination_, length_, orientation_, depth_, bases_, overlap_, 0, 0);
	}

	/** Constructs a trim-aware cross-k edge between anchors inside two long-k tip regions. */
	public Edge(int origin_, int destination_, int length_, int orientation_, int depth_, byte[] bases_,
			int overlap_, int sourceTrim_, int destTrim_){
		origin=origin_;
		destination=destination_;
		length=length_;
		orientation=orientation_;
		depth=depth_;
		bases=bases_;
		overlap=overlap_;
		sourceTrim=sourceTrim_;
		destTrim=destTrim_;
	}
	
	/**
	 * Returns string representation of this edge for debugging.
	 * Delegates to appendTo method for efficient string construction.
	 * @return Formatted string showing edge properties
	 */
	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}
	
	/**
	 * Appends formatted edge information to the given ByteBuilder.
	 * Format: "(destination-orientation-length-depth-bases)"
	 * @param bb ByteBuilder to append to
	 * @return The same ByteBuilder for method chaining
	 */
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append('(');
		bb.append(destination).append('-')/*.append(direction).append('-')*/.append(orientation);
		bb.append('-').append(length).append('-').append(depth).append('-').append(bases);
		if(overlap>0){bb.append("-overlap=").append(overlap);}
		if(sourceTrim>0){bb.append("-sourceTrim=").append(sourceTrim);}
		if(destTrim>0){bb.append("-destTrim=").append(destTrim);}
		bb.append(')');
		return bb;
	}
	
	/**
	 * Generates Graphviz DOT format representation of this edge.
	 * Creates directed edge notation with orientation and length labels
	 * for graph visualization and debugging.
	 * @param bb ByteBuilder to append DOT format output to
	 */
	public void toDot(ByteBuilder bb){
		bb.append(origin);
		bb.append(" -> ");
		bb.append(destination);
		bb.append(" [label=\"").append(((orientation&1)==0) ? "LEFT" : "RIGHT").append("\\nlen=").append(length);
		bb.append("\\norient=").append(orientation).append("\"]").append('\n');
	}
	
	/**
	 * Tests if destination connection is on the right side.
	 * Checks bit 1 of orientation encoding.
	 * @return true if destination connects to right side, false for left
	 */
	public boolean destRight() {
		return (orientation&2)==2;
	}
	/**
	 * Tests if source connection is on the right side.
	 * Checks bit 0 of orientation encoding.
	 * @return true if source connects to right side, false for left
	 */
	public boolean sourceRight() {
		return (orientation&1)==1;
	}
	
	/**
	 * Flips the source orientation and reverse complements sequence data.
	 * Toggles bit 0 of orientation encoding and applies reverse complement
	 * transformation to stored sequence to maintain biological accuracy.
	 */
	void flipSource(){
		if(Tadpole.verbose){System.err.print("Flipping edge source "+this+" -> ");}
		if(bases!=null){Vector.reverseComplementInPlace(bases);}
		orientation^=1;
		if(Tadpole.verbose){System.err.println(this);}
	}
	
	/** Flips the destination orientation encoding.
	 * Toggles bit 1 of orientation encoding without modifying sequence data. */
	void flipDest(){
		if(Tadpole.verbose){System.err.print("Flipping edge dest "+this+" -> ");}
		orientation^=2;
		if(Tadpole.verbose){System.err.println(this);}
	}
	
	/**
	 * Merges another edge with identical topology into this edge.
	 * Combines coverage depths and preserves data from higher-coverage edge.
	 * Requires matching origin, destination, and orientation values.
	 * @param e Edge to merge with this one
	 */
	void merge(Edge e){
		assert(e.origin==origin);
		assert(e.destination==destination);
		assert(e.orientation==orientation);
		if(e.depth>depth){
			length=e.length;
			bases=e.bases;
			orientation=e.orientation;
			overlap=e.overlap;
			sourceTrim=e.sourceTrim;
			destTrim=e.destTrim;
			depth+=e.depth;
		}else{
			depth+=e.depth;
		}
	}
	
	byte[] bases;
	int origin;
	int destination;
	int length;
	int orientation; //left source to left dest; 1 right source to left dest; 2 left source to right dest; 3 right source to right dest
//	int orientation; //0 left kmer, 1 left rkmer, 2 right kmer, 3 right rkmer (of dest)
//	final int direction; //0 forward, 1 backward //They are all forward edges now
	int depth;
	/** Explicit exact overlap for cross-k tip joins; zero means use the active graph k. */
	int overlap;
	/** Bases beyond the source anchor that are replaced by the verified edge path. */
	int sourceTrim;
	/** Bases before the destination anchor that are omitted from the merged sequence. */
	int destTrim;
	
}
