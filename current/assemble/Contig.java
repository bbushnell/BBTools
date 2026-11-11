package assemble;

import java.util.ArrayList;
import java.util.HashMap;

import dna.AminoAcid;
import shared.Tools;
import simd.Vector;
import structures.ByteBuilder;
import tracker.KmerTracker;
import ukmer.Kmer;

/**
 * Mutable assembled sequence with coverage, endpoint state, optional graph edges
 * and graph-quality annotations.
 * Production callers create contigs from Tadpole unitigs, packed KmerCompressor
 * fragments, imported FASTA records and path-merging passes; graph passes may
 * reverse-complement the sequence, mark collapsed nodes and remove or retarget edges.
 * Graph mutation is not internally synchronized; callers coordinate exclusive
 * access while endpoint fields and edge lists are edited.
 *
 * @author Brian Bushnell
 * @date July 12, 2018
 */
public class Contig{

	/**
	 * Wraps the supplied base array without cloning it and uses id 0.
	 *
	 * @param bases_ Base array to retain as this contig's mutable sequence
	 */
	public Contig(byte[] bases_){
		this(bases_, null, 0);
	}

	/**
	 * Wraps the supplied base array without cloning it and stores the numeric id.
	 *
	 * @param bases_ Base array to retain as this contig's mutable sequence
	 * @param id_ Numeric id used in generated headers and graph references
	 */
	public Contig(byte[] bases_, int id_){
		this(bases_, null, id_);
	}

	/**
	 * Wraps the supplied base array without cloning it and stores optional header
	 * text plus the numeric id. Non-null header text is retained; generated graph
	 * classification may still be appended when present.
	 *
	 * @param bases_ Base array to retain as this contig's mutable sequence
	 * @param name_ Header text to retain when non-null
	 * @param id_ Numeric id used in generated headers and graph references
	 */
	public Contig(byte[] bases_, String name_, int id_){
		bases=bases_;
		name=name_;
		id=id_;
	}

	/**
	 * Converts this contig to FASTA format with optional line wrapping.
	 * Appends header and sequence data to the provided ByteBuilder.
	 *
	 * @param wrap Number of bases per line (unlimited if {@code <1})
	 * @param sb ByteBuilder to append FASTA output to
	 * @return The ByteBuilder with FASTA content appended
	 */
	public ByteBuilder toFasta(int wrap, ByteBuilder sb){
		if(wrap<1){wrap=Integer.MAX_VALUE;}
		sb.append('>');
		toHeader(sb);
		if(bases!=null){
			int pos=0;
			while(pos<bases.length-wrap){
				sb.append('\n');
				sb.append(bases, pos, wrap);
				pos+=wrap;
			}
			if(pos<bases.length){
				sb.append('\n');
				sb.append(bases, pos, bases.length-pos);
			}
		}
		return sb;
	}

	/**
	 * Returns string representation with extended contig information.
	 * @return String containing contig details including flip status and usage */
	@Override
	public String toString(){
		return name2();
	}

	/** Returns the FASTA header text without the leading {@code >}. */
	public String name(){
		return toHeader(new ByteBuilder()).toString();
	}

	/**
	 * Gets extended header name with additional status information.
	 * Includes flip status, usage state, and association flags.
	 * @return Extended header string with detailed contig status
	 */
	public String name2(){//Extra info
		ByteBuilder bb=toHeader(new ByteBuilder());
		bb.append(" flipped="+flipped+", used="+used+", associate="+associate);
		return bb.toString();
	}

	/**
	 * Builds contig header with comprehensive metadata.
	 * Includes coverage statistics, GC content, branching information, and edge details.
	 * @param bb ByteBuilder to append header information to
	 * @return The ByteBuilder with header content appended
	 */
	private ByteBuilder toHeader(ByteBuilder bb){
		if(name!=null){
			bb.append(name);
			if(graphClass>=0){bb.append(",class="); appendGraphClass(bb);}
			return bb;
		}
		bb.append("contig_").append(id);
		bb.append(",len=").append(length());
		bb.append(",cov=").append(coverage, 1);
		bb.append(",gc=").append(gc(), 3);
		if(graphClass>=0){bb.append(",class="); appendGraphClass(bb);}
		if(tid>=0){bb.append(",tid_").append(tid);}
		if(VERY_SHORT_NAMES){return bb;}
		bb.append(",min=").append(minCov);
		bb.append(",max=").append(maxCov);
		bb.append(",hh=").append(hh(), 3);
		bb.append(",caga=").append(caga(), 3);
		if(SHORT_NAMES){return bb;}

		bb.append(",left=").append(Tadpole.codeStrings[leftCode]);
		if(leftBranch()){
			bb.append('_').append(leftRatio, 1);
		}
		bb.append(",right=").append(Tadpole.codeStrings[rightCode]);
		if(rightBranch()){
			bb.append('_').append(rightRatio, 1);
		}
		if(leftEdges!=null){
			bb.append(' ').append("leftEdges=");
			appendEdges(leftEdges, bb);
		}
		if(rightEdges!=null){
			bb.append(' ').append("rightEdges=");
			appendEdges(rightEdges, bb);
		}
		return bb;
	}

	/** Appends the graph-derived QC classification without allocating a String. */
	ByteBuilder appendGraphClass(final ByteBuilder bb){
		appendGraphLengthClass(bb).append('/');
		appendGraphDepthClass(bb).append('/');
		return appendGraphTopology(bb);
	}

	/** Appends the long/short graph-length axis label. */
	ByteBuilder appendGraphLengthClass(final ByteBuilder bb){
		if(graphLengthClass==LENGTH_LONG){return bb.append("long");}
		if(graphLengthClass==LENGTH_SHORT){return bb.append("short");}
		return bb.append("unclassified");
	}

	/** Appends the graph-depth bucket label. */
	ByteBuilder appendGraphDepthClass(final ByteBuilder bb){
		if(graphDepthClass==DEPTH_LOW){return bb.append("low");}
		if(graphDepthClass==DEPTH_SEMILOW){return bb.append("semilow");}
		if(graphDepthClass==DEPTH_MEDIUM){return bb.append("medium");}
		if(graphDepthClass==DEPTH_HIGH){return bb.append("high");}
		return bb.append("unclassified");
	}

	/** Appends the topology label, including edge-hop distance when applicable. */
	ByteBuilder appendGraphTopology(final ByteBuilder bb){
		if(graphClass==GRAPH_UNANCHORED){return bb.append("unanchored");}
		if(graphClass==GRAPH_TERMINAL){return bb.append("terminal-").append(graphClassHop);}
		if(graphClass==GRAPH_BRANCHED_TERMINAL){return bb.append("branched-terminal-").append(graphClassHop);}
		if(graphClass==GRAPH_LOOPBACK){return bb.append("loopback-").append(graphClassHop);}
		if(graphClass==GRAPH_CONNECTED){return bb.append("connected-").append(graphClassHop);}
		if(graphClass==GRAPH_BRANCHED_CONNECTED){return bb.append("branched-connected-").append(graphClassHop);}
		if(graphClass==GRAPH_MULTI_CONNECTED){return bb.append("multi-connected-").append(graphClassHop);}
		if(graphClass==GRAPH_SELF_LOOP){return bb.append("self-loop");}
		return bb.append("unclassified");
	}

	/** Appends non-null edges in header form and returns the number appended. */
	int appendEdges(ArrayList<Edge> edges, ByteBuilder bb){
		int added=0;
		if(edges==null || edges.isEmpty()){return 0;}

//		bb.append('_').append('[');
		for(int i=0; i<edges.size(); i++){
			Edge e=edges.get(i);
			if(e!=null){
				if(added>0){bb.append(';');}
				e.appendTo(bb);
				added++;
			}
		}
//		bb.append(']');
		return added;
	}

	/**
	 * Returns the cached GC fraction, computing scalar caches through the fast
	 * mapper if needed. Null or empty input returns 0 immediately; one-base input
	 * returns 0 via {@link #calcScalarsFast()} and leaves {@link #hh()} and
	 * {@link #caga()} at their unset sentinels.
	 * @return GC fraction for supported bases, or 0 for null/short input
	 */
	public float gc(){
		if(bases==null || bases.length<1){return 0;}
		if(gc>-1){return gc;}
		if(true){return calcScalarsFast();}
		int atSum=0, gcSum=0;
		for(byte b : bases){//Only works for ACGTUN
			int gcbit=b>>1;
			atSum+=(~gcbit&1);
			gcSum+=(gcbit&~(b>>3))&1;
		}
		if(gcSum<1){return 0;}
		return gc=(gcSum*1f/Math.max(1, atSum+gcSum));
	}

	/**
	 * Returns the cached dimer-homopolymer statistic, computing it if needed.
	 * Null, empty and one-base inputs leave this cache unset at -1.
	 */
	public float hh(){
		if(hh<0){calcScalarsFast();}
		return hh;
	}

	/**
	 * Returns the cached CAGA dimer statistic, computing it if needed.
	 * Null, empty and one-base inputs leave this cache unset at -1.
	 */
	public float caga(){
		if(hh<0){calcScalarsFast();}
		return caga;//was 'return hh' (copy-paste from hh()); caga field is computed by calcScalarsFast/calcScalars [assemble/Contig#001 FIXED]
	}

	/**
	 * Computes GC and dimer statistics with AminoAcid's A/C/G/T/U mapper.
	 * Unsupported bytes reset dimer continuity and retain their legacy contribution
	 * to the bitwise GC accumulator. Null, empty and one-base inputs leave scalar
	 * caches unset.
	 */
	public void calcScalars(){
		if(bases==null || bases.length<2){return;}
		if(hh>-1){return;}
		final int[] counts=new int[16];
		final int mask=0xF;
		int kmer=0, len=0;
		int atSum=0, gcSum=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			int gcbit=b>>1;
			atSum+=(~gcbit&1);
			gcSum+=(gcbit&~(b>>3))&1;
			final int x=AminoAcid.baseToNumber[b];
			kmer=(((kmer<<2)|x)&mask);
			if(x>=0){
				len++;
				if(len>=2){counts[kmer]++;}
			}else{len=kmer=0;}
		}
		gc=gcSum*1f/Math.max(1, atSum+gcSum);
		hh=KmerTracker.HH(counts);
		caga=KmerTracker.CAGA(counts);
	}

	/**
	 * Computes GC, HH and CAGA with a branch-light ASCII mapper for normal
	 * A/C/G/T/U/N input in either case. N/n, and any other byte with bit 0x08 set,
	 * break dimer continuity; other unsupported bytes are not rejected and are
	 * interpreted from their low bits. Empty dimer subsets make
	 * {@link KmerTracker#HH(int[])} return 0 and
	 * {@link KmerTracker#CAGA(int[])} return 0.5. Null, empty and one-base inputs
	 * return 0 without setting the scalar caches.
	 * @return GC fraction, or 0 for null/short input
	 */
	public float calcScalarsFast(){
		if(bases==null || bases.length<2){return 0;}
		if(hh>-1){return gc;}
		final int[] counts=new int[16];

		int prevBad=8;//Initialize as "N" so first dimer fails
		int prevVal=0;
		int atSum=0, gcSum=0;

		for(byte b : bases){
			// GC Logic (Dimer-based number is imprecise)
			final int gcbit=b>>1;
			atSum+=(~gcbit&1);
			gcSum+=(gcbit&~(b>>3))&1;

			// Dimer Logic (Branchless Mapping)
			// Map ASCII directly to 0-3
			int val=(b&6)>>1;
			val^=(val&2)>>1;

			int bad=(b&8); // 8 if N, 0 if ACGT

			// Only count if both Previous and Current bases are valid
			// Since N is rare, this branch is highly predictable
			if((prevBad|bad)==0){
				counts[(prevVal<<2)|val]++;
			}

			prevVal=val;
			prevBad=bad;
		}
		hh=KmerTracker.HH(counts);
		caga=KmerTracker.CAGA(counts);
		return gc=gcSum*1f/Math.max(1, atSum+gcSum);
	}

	/** Returns the number of stored bases. */
	public int length(){return bases.length;}

	/** Returns the leftmost kmer as a 2-bit code. */
	long leftKmer(int k){
		assert(k<32);
		long kmer=0;
		for(int i=0; i<k; i++){
			byte b=bases[i];
			byte x=AminoAcid.baseToNumber[b];
			assert(x>=0);
			kmer=(kmer<<2)|x;
		}
		return kmer;
	}

	/** Returns the rightmost kmer as a 2-bit code. */
	long rightKmer(int k){
		assert(k<32);
		long kmer=0;
		for(int i=bases.length-k; i<bases.length; i++){
			byte b=bases[i];
			byte x=AminoAcid.baseToNumber[b];
			assert(x>=0);
			kmer=(kmer<<2)|x;
		}
		return kmer;
	}

	/** Writes the leftmost {@code kmer.kbig} bases into a reusable wide kmer. */
	public Kmer leftKmer(Kmer kmer){
		kmer.clearFast();
		for(int i=0; i<kmer.kbig; i++){
			kmer.addRight(bases[i]);
		}
		return kmer;
	}

	/** Writes the rightmost {@code kmer.kbig} bases into a reusable wide kmer. */
	public Kmer rightKmer(Kmer kmer){
		kmer.clearFast();
		for(int i=bases.length-kmer.kbig; i<bases.length; i++){
			kmer.addRight(bases[i]);
		}
		return kmer;
	}

	/**
	 * Determines if this contig is in canonical orientation.
	 * Compares the sequence to its reverse complement to determine lexicographic order.
	 * @return true if contig is in canonical form, false otherwise
	 */
	boolean canonical(){
		//[assemble/Contig#002 FIXED 2026-07-02] Was `baseToComplementExtended[i]` (indexed by loop POSITION not by a base) => baseToComplementExtended[0]==0 < any base byte => returned true on iteration 0 for every contig => canonical() was constant-true => the rcomp() canonicalization at Tadpole1:689/Tadpole2:552 NEVER fired => output contig ORIENTATION was left build-order-dependent (nondeterministic under multithreading), which the orientation-sensitive popBubbles then amplified into a nondeterministic contig COUNT. Correct test: a contig is canonical iff its sequence <= its reverse-complement, i.e. compare bases[i] to complement(bases[n-1-i]).
		for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
			final byte a=bases[i], b=AminoAcid.baseToComplementExtended[bases[j]];
			if(a<b){return true;}
			else if(b<a){return false;}
		}
		return true;
//		final long leftKmer=leftKmer(k);
//		final long rightKmer=rightKmer(k);
//		final long leftRKmer=AminoAcid.reverseComplementBinaryFast(leftKmer, k);
//		assert((rightKmer!=leftKmer && rightKmer!=leftRKmer) || length()==k /* Palindrome case */);
//		return rightKmer>=leftRKmer;
	}

	/**
	 * Reverse complements this contig's sequence and swaps endpoint metadata.
	 * Used for unconnected contigs before edge source/destination orientation
	 * updates are required.
	 */
	public void rcomp(){
		Vector.reverseComplementInPlace(bases);
		{
			int temp=leftCode;
			leftCode=rightCode;
			rightCode=temp;
		}
		{
			float temp=leftRatio;
			leftRatio=rightRatio;
			rightRatio=temp;
		}
		{
			boolean temp=leftBridgeEndpoint;
			leftBridgeEndpoint=rightBridgeEndpoint;
			rightBridgeEndpoint=temp;
		}
		{
			final boolean temp=leftFusionEndpoint;
			leftFusionEndpoint=rightFusionEndpoint;
			rightFusionEndpoint=temp;
		}
		{
			ArrayList<Edge> temp=leftEdges;
			leftEdges=rightEdges;
			rightEdges=temp;
		}
	}

	/** Returns whether the left endpoint has any branch code. */
	boolean leftBranch(){return Tadpole.isBranchCode(leftCode);}
	/** Returns whether the right endpoint has any branch code. */
	public boolean rightBranch(){return Tadpole.isBranchCode(rightCode);}

	/** Returns whether the left endpoint branches forward. */
	boolean leftForwardBranch(){return leftCode==Tadpole.F_BRANCH;}
	/** Returns whether the right endpoint branches forward. */
	public boolean rightForwardBranch(){return rightCode==Tadpole.F_BRANCH;}

	/** Returns whether the left endpoint branches backward. */
	boolean leftBackwardBranch(){return leftCode==Tadpole.B_BRANCH;}
	/** Returns whether the right endpoint branches backward. */
	public boolean rightBackwardBranch(){return rightCode==Tadpole.B_BRANCH;}

	/** Returns the number of retained left-edge records. */
	public final int leftEdgeCount(){
		return leftEdges==null ? 0 : leftEdges.size();
	}

	/** Returns the number of retained right-edge records. */
	public final int rightEdgeCount(){
		return rightEdges==null ? 0 : rightEdges.size();
	}

	/**
	 * Adds an edge connection to the left end of this contig.
	 * Merges with existing edges to the same destination if appropriate,
	 * saturating accumulated support depth at {@link Integer#MAX_VALUE}.
	 * @param e The edge to add
	 */
	public void addLeftEdge(Edge e){
		assert(!hasLeftEdge(e.destination, e.bases, e.orientation)) : e+"\n"+leftEdges;
		Edge old=getLeftEdge(e.destination, e.orientation);
		if(old!=null){
			if(e.depth>=old.depth && (old.depth==1 || old.length==e.length)){//Merge edges
				old.bases=e.bases;
				old.length=e.length;
				old.overlap=e.overlap;
				old.sourceTrim=e.sourceTrim;
				old.destTrim=e.destTrim;
				// ASM051: duplicate-edge support can exceed int range; cap the
				// retained count instead of wrapping negative.
				old.depth=(int)Math.min(Integer.MAX_VALUE, (long)old.depth+e.depth);
				return;
			}
		}
		if(leftEdges==null){leftEdges=new ArrayList<Edge>(4);}
		leftEdges.add(e);
	}

	/**
	 * Adds an edge connection to the right end of this contig.
	 * Merges with existing edges to the same destination if appropriate,
	 * saturating accumulated support depth at {@link Integer#MAX_VALUE}.
	 * @param e The edge to add
	 */
	public void addRightEdge(Edge e){
		assert(!hasRightEdge(e.destination, e.bases, e.orientation)) : e+"\n"+rightEdges;
		Edge old=getRightEdge(e.destination, e.orientation);
		if(old!=null){
			if(e.depth>=old.depth && (old.depth==1 || old.length==e.length)){//Merge edges
				old.bases=e.bases;
				old.length=e.length;
				old.overlap=e.overlap;
				old.sourceTrim=e.sourceTrim;
				old.destTrim=e.destTrim;
				// ASM051: duplicate-edge support can exceed int range; cap the
				// retained count instead of wrapping negative.
				old.depth=(int)Math.min(Integer.MAX_VALUE, (long)old.depth+e.depth);
				return;
			}
		}
		if(rightEdges==null){rightEdges=new ArrayList<Edge>(4);}
		rightEdges.add(e);
	}

//	public boolean hasLeftEdge(int dest) {
//		return hasLeftEdge(dest, null);
//	}
//
//	public boolean hasRightEdge(int dest) {
//		return hasRightEdge(dest, null);
//	}

	/** Removes left edges pointing to one destination side. */
	public void removeLeftEdge(int dest, boolean destRight){leftEdges=removeEdge(dest, destRight, leftEdges);}
	/** Removes right edges pointing to one destination side. */
	public void removeRightEdge(int dest, boolean destRight){rightEdges=removeEdge(dest, destRight, rightEdges);}

	/** Removes edges pointing to one destination side and returns null when none remain. */
	public static ArrayList<Edge> removeEdge(int dest, boolean destRight, ArrayList<Edge> edges){
		if(edges==null){return null;}
		int removed=0;
		for(int i=0; i<edges.size(); i++){
			Edge e=edges.get(i);
			assert(e!=null);
			if(e.destination==dest && e.destRight()==destRight){
				edges.set(i, null);
				removed++;
			}
		}
		if(removed>0){Tools.condenseStrict(edges);}
		return (edges.isEmpty() ? null : edges);
	}

	/** Returns the first left edge with a matching destination and optional orientation. */
	public Edge getLeftEdge(int dest, int orientation){
		if(leftEdges==null){return null;}
		for(Edge e : leftEdges){
			if(e.destination==dest && (orientation<0 || orientation==e.orientation)){
				return e;
			}
		}
		return null;
	}

	/** Returns the first right edge with a matching destination and optional orientation. */
	public Edge getRightEdge(int dest, int orientation){
		if(rightEdges==null){return null;}
		for(Edge e : rightEdges){
			if(e.destination==dest && (orientation<0 || orientation==e.orientation)){
				return e;
			}
		}
		return null;
	}

	/** Returns whether a left edge matches destination, orientation and optional payload. */
	public boolean hasLeftEdge(int dest, byte[] path, int orientation){
		if(leftEdges!=null){
			for(Edge e : leftEdges){
				if(e.destination==dest && e.orientation==orientation){
					if(path==null || Tools.equals(path, e.bases)){
						return true;
					}
				}
			}
		}
		return false;
	}

	/** Returns whether a right edge matches destination, orientation and optional payload. */
	public boolean hasRightEdge(int dest, byte[] path, int orientation){
		if(rightEdges!=null){
			for(Edge e : rightEdges){
				if(e.destination==dest && e.orientation==orientation){
					if(path==null || Tools.equals(path, e.bases)){
						return true;
					}
				}
			}
		}
		return false;
	}

	/** Returns whether {@link #flip(ArrayList)} has reversed this contig an odd number of times. */
	public final boolean flipped(){return flipped;}

	/** Captures assembly-K dead ends before TadpoleMulti clears the original graph. */
	public void markFusionEndpoints(){
		assert(!used && !associate) : "Only surviving initial contigs may supply historical fusion endpoints.";
		leftFusionEndpoint=(leftCode==Tadpole.DEAD_END && leftEdgeCount()==0);
		rightFusionEndpoint=(rightCode==Tadpole.DEAD_END && rightEdgeCount()==0);
	}

	/**
	 * Reverse complements the sequence, swaps endpoint fields and updates local
	 * edge source orientations plus every supplied inbound edge destination.
	 * Caller must supply the complete inbound edge list for this contig if all
	 * touching edge orientations are expected to stay synchronized.
	 */
	public final void flip(ArrayList<Edge> inbound){
		if(Tadpole.verbose){System.err.println("Flipping contig "+name());}
		flipped=!flipped;
		Vector.reverseComplementInPlace(bases);
		{
			int temp=leftCode;
			leftCode=rightCode;
			rightCode=temp;
		}
		{
			float temp=leftRatio;
			leftRatio=rightRatio;
			rightRatio=temp;
		}
		{
			boolean temp=leftBridgeEndpoint;
			leftBridgeEndpoint=rightBridgeEndpoint;
			rightBridgeEndpoint=temp;
		}
		{
			final boolean temp=leftFusionEndpoint;
			leftFusionEndpoint=rightFusionEndpoint;
			rightFusionEndpoint=temp;
		}
		{
			ArrayList<Edge> temp=leftEdges;
			leftEdges=rightEdges;
			rightEdges=temp;
			if(leftEdges!=null){
				for(Edge e : leftEdges){
					e.flipSource();
				}
			}
			if(rightEdges!=null){
				for(Edge e : rightEdges){
					e.flipSource();
				}
			}
		}
		if(inbound!=null){
			for(Edge e : inbound){
				assert(e.destination==id) : id+", "+e.destination+"\n"+e+"\n"+name()+"\n";
				e.flipDest();
			}
		}
		if(Tadpole.verbose){System.err.println("Flipped contig "+name());}
	}

	/** Returns whether this contig has been consumed as a path member or merge input. */
	public final boolean used(){return used;}

	/** Looks up this contig's inbound list and delegates consumed-state edge removal. */
	final void setUsed(HashMap<Integer, ArrayList<Edge>> map, ArrayList<Contig> allContigs){setUsed(map.get(id), allContigs);}

	/**
	 * Marks this contig consumed and removes local edges plus every supplied
	 * inbound edge from its live origin. Caller must pass the complete inbound
	 * list for this contig to remove all graph edges touching it.
	 */
	final void setUsed(ArrayList<Edge> inbound, ArrayList<Contig> allContigs){
		assert(!used);
		assert(!associate);
		used=true;
		removeAllEdges(inbound, allContigs);
	}

	/** Returns whether this contig is the discarded associate path of a collapsed bubble. */
	public final boolean associate(){return associate;}

	/** Looks up this contig's inbound list and delegates associate-state edge removal. */
	final void setAssociate(HashMap<Integer, ArrayList<Edge>> map, ArrayList<Contig> allContigs){setAssociate(map.get(id), allContigs);}

	/**
	 * Marks this contig as a collapsed bubble associate and removes local edges plus
	 * every supplied inbound edge from its live origin. Caller must pass the
	 * complete inbound list for this contig to remove all graph edges touching it.
	 */
	final void setAssociate(ArrayList<Edge> inbound, ArrayList<Contig> allContigs){
		assert(!used);
		assert(!associate);
		associate=true;
		removeAllEdges(inbound, allContigs);
	}

	/**
	 * Clears local endpoint edges and removes the supplied inbound references from
	 * still-live source contigs. Caller must supply every inbound edge for this
	 * contig when all touching edge references should be removed.
	 */
	void removeAllEdges(ArrayList<Edge> inbound, ArrayList<Contig> allContigs){
		leftEdges=null;
		rightEdges=null;

		if(inbound!=null){
			for(Edge e : inbound){
				if(e.destination==id){
					Contig source=allContigs.get(e.origin);
					if(source!=this && !source.used() && !source.associate()){
						source.removeLeftEdgesTo(id);
						source.removeRightEdgesTo(id);
					}
				}
			}
		}
	}

	/** Removes local left edges to a destination node after a peer is consumed. */
	private void removeLeftEdgesTo(int destNode){
		leftEdges=removeEdgesTo(leftEdges, destNode);
	}

	/** Removes local right edges to a destination node after a peer is consumed. */
	private void removeRightEdgesTo(int destNode){
		rightEdges=removeEdgesTo(rightEdges, destNode);
	}

	/** Removes all edges to {@code destNode} and condenses or clears the list. */
	private static ArrayList<Edge> removeEdgesTo(ArrayList<Edge> edges, final int destNode){
		if(edges==null || edges.isEmpty()){return null;}
		int cnt=0;
		for(int i=0; i<edges.size(); i++){
			Edge e=edges.get(i);
			if(e.destination==destNode){
				edges.set(i, null);
				cnt++;
			}
		}
		if(cnt==edges.size()){return null;}
		Tools.condenseStrict(edges);
		return edges;
	}

	/**
	 * Retargets this contig and local endpoint-edge origins from its old id to
	 * {@code newID}. Caller must supply every inbound edge when all destination
	 * references should be retargeted.
	 */
	final void renumber(final int newID, ArrayList<Edge> inbound){
		if(id==newID){return;}
		if(leftEdges!=null){
			for(Edge e : leftEdges){
				assert(e.origin==id);
				e.origin=newID;
			}
		}
		if(rightEdges!=null){
			for(Edge e : rightEdges){
				assert(e.origin==id);
				e.origin=newID;
			}
		}
		if(inbound!=null){
			for(Edge e : inbound){
				assert(e.destination==id);
				e.destination=newID;
			}
		}

		id=newID;
	}

	/** Optional explicit FASTA header; non-null values bypass generated metrics. */
	public String name;
	/** Mutable retained sequence array; callers that change bases must invalidate scalar caches. */
	public byte[] bases;
	/** Mean valid assembly-kmer depth across the contig. */
	public float coverage;
	/** Graph-only reachability was truncated by a distance or state bound at these ends. */
	public boolean graphLeftLimited, graphRightLimited;
	/** Minimum valid assembly-kmer depth observed by graph annotation. */
	public int minCov;
	/** Maximum valid assembly-kmer depth observed by graph annotation. */
	public int maxCov;
	/** Tadpole endpoint code for the left side. */
	public int leftCode;
	/** Tadpole endpoint code for the right side. */
	public int rightCode;
	/** Branch ratio for the left endpoint code. */
	public float leftRatio;
	/** Branch ratio for the right endpoint code. */
	public float rightRatio;
	/** Whether the left side remains eligible for cross-k bridge/tip searches. */
	public boolean leftBridgeEndpoint;
	/** Whether the right side remains eligible for cross-k bridge/tip searches. */
	public boolean rightBridgeEndpoint;
	/** Original assembly-K dead-end eligibility; follows the surviving outer ends. */
	public boolean leftFusionEndpoint, rightFusionEndpoint;
	/** Contig id used as the index for edge origin/destination references. */
	public int id;
	/** Combined graph QC class; nonnegative values are GRAPH_* labels. */
	public int graphClass=GRAPH_UNCLASSIFIED;
	/** Edge-hop distance for anchored topology graph classes. */
	public int graphClassHop=0;
	/** Length QC axis; nonnegative values are LENGTH_* labels. */
	public int graphLengthClass=LENGTH_UNCLASSIFIED;
	/** Depth QC axis; nonnegative values are DEPTH_* labels. */
	public int graphDepthClass=DEPTH_UNCLASSIFIED;
	/** Optional taxonomy ID carried from the Tadpole taxID setting. */
	public int tid=-1;
	/** Lazily computed sequence scalars; -1 means unset. */
	float gc=-1, hh=-1, caga=-1;

	/** Whether this contig has been reverse-complemented by flip an odd number of times. */
	private boolean flipped=false;
	/** Whether this contig has been consumed as a path member or merge input. */
	private boolean used=false;
	/** Whether this contig is the discarded path in a collapsed bubble. */
	private boolean associate=false;

//	public Edge[] leftEdges;
//	public Edge[] rightEdges;

	/** Outbound edges anchored to the left endpoint; entries must originate from {@link #id}. */
	public ArrayList<Edge> leftEdges;
	/** Outbound edges anchored to the right endpoint; entries must originate from {@link #id}. */
	public ArrayList<Edge> rightEdges;

	/** Emits generated headers without edge lists when true. */
	public static boolean SHORT_NAMES=true;
	/** Emits generated headers without min/max, HH/CAGA or edge details when true. */
	public static boolean VERY_SHORT_NAMES=false;
	/** Length QC labels used by graph classification. */
	public static final int LENGTH_UNCLASSIFIED=-1, LENGTH_SHORT=0, LENGTH_LONG=1;
	/** Depth QC labels used by graph classification. */
	public static final int DEPTH_UNCLASSIFIED=-1, DEPTH_LOW=0, DEPTH_SEMILOW=1,
			DEPTH_MEDIUM=2, DEPTH_HIGH=3;
	/** Topology QC labels used by graph classification. */
	public static final int GRAPH_UNCLASSIFIED=-1, GRAPH_UNANCHORED=0, GRAPH_TERMINAL=1,
			GRAPH_BRANCHED_TERMINAL=2, GRAPH_LOOPBACK=3, GRAPH_CONNECTED=4,
			GRAPH_BRANCHED_CONNECTED=5, GRAPH_MULTI_CONNECTED=6, GRAPH_SELF_LOOP=7;
}
