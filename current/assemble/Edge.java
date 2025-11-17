package assemble;

import shared.Vector;
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
	
	/**
	 * Constructs an edge connecting two contigs in the assembly graph.
	 *
	 * @param origin_ Source contig ID
	 * @param destination_ Target contig ID
	 * @param length_ Overlap length in bases between connected contigs
	 * @param orientation_ 2-bit encoding: bit 0 = source side, bit 1 = dest side
	 * @param depth_ Coverage depth representing support evidence from reads
	 * @param bases_ Optional sequence data for the overlap region (may be null)
	 */
	public Edge(int origin_, int destination_, int length_, int orientation_, int depth_, byte[] bases_){
		origin=origin_;
		destination=destination_;
		length=length_;
		orientation=orientation_;
		depth=depth_;
		bases=bases_;
	}
	
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
			depth+=e.depth;
		}else{
			depth+=e.depth;
		}
	}
	
	/** Optional sequence data for the overlap region between contigs */
	byte[] bases;
	/** Source contig ID in the assembly graph */
	int origin;
	/** Target contig ID in the assembly graph */
	int destination;
	/** Overlap length in bases between connected contigs */
	int length;
	/**
	 * 2-bit orientation encoding: 0=left-to-left, 1=right-to-left, 2=left-to-right, 3=right-to-right
	 */
	int orientation; //left source to left dest; 1 right source to left dest; 2 left source to right dest; 3 right source to right dest
//	int orientation; //0 left kmer, 1 left rkmer, 2 right kmer, 3 right rkmer (of dest)
//	final int direction; //0 forward, 1 backward //They are all forward edges now
	/** Coverage depth representing support evidence from read data */
	int depth;
	
}
