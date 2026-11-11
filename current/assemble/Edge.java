package assemble;

import simd.Vector;
import structures.ByteBuilder;

/**
 * Directed connection between two contigs in a De Bruijn graph.
 * Stores endpoint orientation, supporting depth, optional emitted bases and
 * cross-k trim metadata used when a verified path replaces a tip interval.
 *
 * @author Brian Bushnell
 */
public class Edge{

	/**
	 * Constructs an edge whose overlap is inferred from the active graph k.
	 *
	 * @param origin_ source contig id
	 * @param destination_ destination contig id
	 * @param length_ emitted base count for dense Tadpole paths, or kmer-step
	 *        distance for graph-only paths without a payload
	 * @param orientation_ two-bit source/destination endpoint orientation
	 * @param depth_ supporting read depth
	 * @param bases_ optional emitted payload bases; retained and later mutated in
	 *        place by {@link #flipSource()}
	 */
	public Edge(int origin_, int destination_, int length_, int orientation_, int depth_, byte[] bases_){
		this(origin_, destination_, length_, orientation_, depth_, bases_, 0);
	}

	/**
	 * Constructs an edge with an explicit overlap.
	 *
	 * @param overlap_ exact overlap in bases; zero means use the active graph k
	 */
	public Edge(int origin_, int destination_, int length_, int orientation_, int depth_, byte[] bases_, int overlap_){
		this(origin_, destination_, length_, orientation_, depth_, bases_, overlap_, 0, 0);
	}

	/**
	 * Constructs a trim-aware cross-k edge between anchors inside two long-k tip regions.
	 * Cross-k overlap links use {@code length==0} and no payload; their explicit
	 * overlap and trim fields describe the sequence to retain from each contig.
	 *
	 * @param sourceTrim_ bases past the source anchor that the payload replaces
	 * @param destTrim_ bases before the destination anchor that are trimmed before joining
	 */
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

	/** Returns a debug string containing destination, orientation, length, depth and payload fields. */
	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}

	/**
	 * Appends this edge in compact debug form:
	 * {@code (destination-orientation-length-depth-bases[-overlap=N][-sourceTrim=N][-destTrim=N])}.
	 *
	 * @param bb destination builder
	 * @return the same builder for chaining
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
	 * Appends a Graphviz DOT edge from {@link #origin} to {@link #destination}.
	 *
	 * @param bb destination builder
	 */
	public void toDot(ByteBuilder bb){
		bb.append(origin);
		bb.append(" -> ");
		bb.append(destination);
		bb.append(" [label=\"").append(((orientation&1)==0) ? "LEFT" : "RIGHT").append("\\nlen=").append(length);
		bb.append("\\norient=").append(orientation).append("\"]").append('\n');
	}

	/** Returns whether the destination endpoint is on the destination contig's right side. */
	public boolean destRight(){
		return (orientation&2)==2;
	}

	/** Returns whether the source endpoint is on the source contig's right side. */
	public boolean sourceRight(){
		return (orientation&1)==1;
	}

	/**
	 * Flips the source endpoint and reverse-complements emitted bases in place.
	 */
	public void flipSource(){
		if(Tadpole.verbose){System.err.print("Flipping edge source "+this+" -> ");}
		if(bases!=null){Vector.reverseComplementInPlace(bases);}
		orientation^=1;
		if(Tadpole.verbose){System.err.println(this);}
	}

	/**
	 * Reads an emitted base in outward traversal order.
	 * Left-source payloads are read from the opposite end and complemented to
	 * restore their biological orientation.
	 */
	final byte traversalBase(final int pos){
		assert(bases!=null && pos>=0 && pos<bases.length) :
				"Graph path consumers require an emitted base inside the edge payload: "+pos+", edge="+this;
		return sourceRight() ? bases[pos] : dna.AminoAcid.baseToComplementExtended[bases[bases.length-1-pos]];
	}

	/** Flips the destination endpoint without modifying emitted bases. */
	void flipDest(){
		if(Tadpole.verbose){System.err.print("Flipping edge dest "+this+" -> ");}
		orientation^=2;
		if(Tadpole.verbose){System.err.println(this);}
	}

	/**
	 * Merges depth from a same-topology edge.
	 * Payload, length, overlap and trims are replaced only when the incoming edge
	 * has greater supporting depth than this edge had before the merge; the path
	 * depth statistics remain untouched. Accumulated supporting depth saturates at
	 * {@link Integer#MAX_VALUE}. The chosen payload reference is retained rather
	 * than copied.
	 *
	 * @param e same-origin, same-destination, same-orientation edge to merge
	 */
	void merge(Edge e){
		assert(e.origin==origin);
		assert(e.destination==destination);
		assert(e.orientation==orientation);
		final int oldDepth=depth;
		if(e.depth>oldDepth){
			length=e.length;
			bases=e.bases;
			orientation=e.orientation;
			overlap=e.overlap;
			sourceTrim=e.sourceTrim;
			destTrim=e.destTrim;
		}
		depth=(int)Math.min(Integer.MAX_VALUE, (long)oldDepth+e.depth);
	}

	/** Optional emitted payload; ordinary path payloads have {@code bases.length==length}. */
	public byte[] bases;
	/** Source contig id. */
	public int origin;
	/** Destination contig id. */
	public int destination;
	/** Emitted bases for dense paths, kmer steps for graph-only paths, or zero for cross-k overlap links. */
	public int length;
	/** Two-bit endpoint orientation: source is bit 0, destination is bit 1. */
	public int orientation; //left source to left dest; 1 right source to left dest; 2 left source to right dest; 3 right source to right dest
//	int orientation; //0 left kmer, 1 left rkmer, 2 right kmer, 3 right rkmer (of dest)
//	final int direction; //0 forward, 1 backward //They are all forward edges now
	/**
	 * Supporting read depth accumulated across merged same-topology edges, capped
	 * at {@link Integer#MAX_VALUE}.
	 */
	public int depth;
	/** Minimum read-kmer depth along a graph-only representative path; -1 if not measured. */
	public int pathMinDepth=-1;
	/** Maximum read-kmer depth along the same graph-only representative path. */
	public int pathMaxDepth=-1;
	/** Mean read-kmer depth along that path, including both endpoint kmers. */
	public double pathMeanDepth=-1;
	/** Explicit exact overlap for cross-k tip joins; zero means use the active graph k. */
	public int overlap;
	/** Bases beyond the source anchor that are replaced by the verified edge path. */
	public int sourceTrim;
	/** Bases before the destination anchor that are omitted from the merged sequence. */
	public int destTrim;

}
