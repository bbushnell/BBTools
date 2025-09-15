package assemble;

import dna.AminoAcid;
import structures.ByteBuilder;

public class Edge {
	
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
	
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append('(');
		bb.append(destination).append('-')/*.append(direction).append('-')*/.append(orientation);
		bb.append('-').append(length).append('-').append(depth).append('-').append(bases);
		bb.append(')');
		return bb;
	}
	
	public void toDot(ByteBuilder bb){
		bb.append(origin);
		bb.append(" -> ");
		bb.append(destination);
		bb.append(" [label=\"").append(((orientation&1)==0) ? "LEFT" : "RIGHT").append("\\nlen=").append(length);
		bb.append("\\norient=").append(orientation).append("\"]").append('\n');
	}
	
	public boolean destRight() {
		return (orientation&2)==2;
	}
	public boolean sourceRight() {
		return (orientation&1)==1;
	}
	
	void flipSource(){
		if(Tadpole.verbose){System.err.print("Flipping edge source "+this+" -> ");}
		if(bases!=null){AminoAcid.reverseComplementBasesInPlace(bases);}
		orientation^=1;
		if(Tadpole.verbose){System.err.println(this);}
	}
	
	void flipDest(){
		if(Tadpole.verbose){System.err.print("Flipping edge dest "+this+" -> ");}
		orientation^=2;
		if(Tadpole.verbose){System.err.println(this);}
	}
	
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
	
	byte[] bases;
	int origin;
	int destination;
	int length;
	int orientation; //left source to left dest; 1 right source to left dest; 2 left source to right dest; 3 right source to right dest
//	int orientation; //0 left kmer, 1 left rkmer, 2 right kmer, 3 right rkmer (of dest)
//	final int direction; //0 forward, 1 backward //They are all forward edges now
	int depth;
	
}
