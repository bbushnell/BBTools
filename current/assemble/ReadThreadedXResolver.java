package assemble;

import java.util.ArrayList;
import java.util.HashMap;

import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/**
 * Finds closed 2-by-2 repeat junctions and records reads that traverse both
 * boundaries.  Topology alone never authorizes a graph change: an X is
 * resolvable only when both inferred paths have separate, nonconflicting
 * read support.
 */
class ReadThreadedXResolver {

	ReadThreadedXResolver(ArrayList<Contig> contigs_, HashMap<Integer, ArrayList<Edge>> destMap_,
			int k_, int minSupport_, int maxNoise_){
		contigs=contigs_;
		destMap=destMap_;
		k=k_;
		minSupport=minSupport_;
		maxNoise=maxNoise_;
		candidates=new Candidate[contigs.size()];
	}

	/** Finds structurally closed X junctions without changing the graph. */
	int findCandidates(){
		int found=0;
		for(int i=0; i<contigs.size(); i++){
			final Contig mid=contigs.get(i);
			if(mid.id!=i){
				throw new IllegalStateException("Contig IDs must match graph indexes: index="+i+", id="+mid.id);
			}
			Candidate c=makeCandidate(mid);
			if(c!=null){candidates[mid.id]=c; found++;}
		}
		return found;
	}

	/** Records every consecutive owned-contig triple in one read path. */
	void observePath(IntList path){
		observePath(path, new LongList(4));
	}

	/** Records each candidate/path combination at most once per physical read. */
	void observePath(IntList path, LongList seen){
		for(int i=2; i<path.size; i++){
			final int a=path.array[i-2], mid=path.array[i-1], b=path.array[i];
			if(a!=mid && mid!=b && a!=b){
				final long key=observationKey(a, mid, b);
				if(key>=0 && !seen.contains(key)){
					seen.add(key);
					increment(key);
				}
			}
		}
	}

	/** Records one traversal in either read orientation. */
	void observe(final int a, final int mid, final int b){
		final long key=observationKey(a, mid, b);
		if(key>=0){increment(key);}
	}

	private long observationKey(final int a, final int mid, final int b){
		if(mid<0 || mid>=candidates.length){return -1;}
		final Candidate c=candidates[mid];
		if(c==null){return -1;}
		int li=indexOf(c.left, a), ri=indexOf(c.right, b);
		if(li<0 || ri<0){
			li=indexOf(c.left, b);
			ri=indexOf(c.right, a);
		}
		return li>=0 && ri>=0 ? (((long)mid)<<2)|((li<<1)|ri) : -1;
	}

	private void increment(final long key){
		final Candidate c=candidates[(int)(key>>>2)];
		synchronized(c){c.support[(int)(key&3)]++;}
	}

	/** Returns the unique supported pairing, or -1 when evidence is insufficient/ambiguous. */
	int pairing(final int mid){
		if(mid<0 || mid>=candidates.length || candidates[mid]==null){return -1;}
		return candidates[mid].pairing(minSupport, maxNoise);
	}

	long support(final int mid, final int leftIndex, final int rightIndex){
		return candidates[mid].support[(leftIndex<<1)|rightIndex];
	}

	/** Resolves every nonoverlapping X with a unique, nonconflicting read-supported pairing. */
	int resolveSupported(){
		int resolved=0;
		for(Candidate c : candidates){
			if(c!=null && c.pairing(minSupport, maxNoise)>=0 && resolve(c)){resolved++;}
		}
		if(resolved>0){rebuildDestMap();}
		return resolved;
	}

	Candidate candidate(final int mid){return candidates[mid];}
	int contigCount(){return candidates.length;}

	long observations(){
		long count=0;
		for(Candidate c : candidates){if(c!=null){count+=c.totalSupport();}}
		return count;
	}

	int candidatesWithEvidence(){
		int count=0;
		for(Candidate c : candidates){if(c!=null && c.totalSupport()>0){count++;}}
		return count;
	}

	int resolvableCandidates(){
		int count=0;
		for(Candidate c : candidates){
			if(c!=null && c.pairing(minSupport, maxNoise)>=0){count++;}
		}
		return count;
	}

	private Candidate makeCandidate(final Contig mid){
		if(mid==null || mid.used() || mid.associate()
				|| mid.leftEdgeCount()!=2 || mid.rightEdgeCount()!=2){return null;}
		final int[] left=neighborIDs(mid.leftEdges);
		final int[] right=neighborIDs(mid.rightEdges);
		if(left==null || right==null || left[0]==left[1] || right[0]==right[1]){return null;}
		if(left[0]==mid.id || left[1]==mid.id || right[0]==mid.id || right[1]==mid.id
				|| left[0]==right[0] || left[0]==right[1]
				|| left[1]==right[0] || left[1]==right[1]){return null;}
		if(countInbound(mid.id, false)!=2 || countInbound(mid.id, true)!=2){return null;}
		for(Edge e : mid.leftEdges){if(!validClosedNeighbor(mid, e, false)){return null;}}
		for(Edge e : mid.rightEdges){if(!validClosedNeighbor(mid, e, true)){return null;}}
		return new Candidate(mid.id, left, right);
	}

	private boolean resolve(final Candidate c){
		final Contig mid=contigs.get(c.middle);
		if(mid.used() || mid.associate()){return false;}
		final Candidate current=makeCandidate(mid);
		if(current==null || !sameNeighbors(c.left, current.left) || !sameNeighbors(c.right, current.right)){
			return false;
		}
		final int pairing=c.pairing(minSupport, maxNoise);
		if(pairing<0){return false;}
		final Contig[] left={contigs.get(c.left[0]), contigs.get(c.left[1])};
		final Contig[] right={contigs.get(c.right[0]), contigs.get(c.right[1])};
		for(int i=0; i<2; i++){
			final Edge leftEdge=mid.leftEdges.get(indexOfEdge(mid.leftEdges, left[i].id));
			final Edge rightEdge=mid.rightEdges.get(indexOfEdge(mid.rightEdges, right[i].id));
			if(!exteriorSafe(left[i], !leftEdge.destRight(), c)
					|| !exteriorSafe(right[i], !rightEdge.destRight(), c)){return false;}
		}

		final boolean[] flippedLeft=new boolean[2], flippedRight=new boolean[2];
		for(int i=0; i<2; i++){
			final Edge e=mid.leftEdges.get(indexOfEdge(mid.leftEdges, left[i].id));
			if(!e.destRight()){
				left[i].flip(destMap.get(left[i].id));
				flippedLeft[i]=true;
			}
		}
		for(int i=0; i<2; i++){
			final Edge e=mid.rightEdges.get(indexOfEdge(mid.rightEdges, right[i].id));
			if(e.destRight()){
				right[i].flip(destMap.get(right[i].id));
				flippedRight[i]=true;
			}
		}

		final Contig[] pairedRight={right[pairing], right[1-pairing]};
		final byte[][] paths=new byte[2][];
		final Edge[] entries=new Edge[2], exits=new Edge[2];
		for(int i=0; i<2; i++){
			entries[i]=left[i].getRightEdge(mid.id, 1);
			exits[i]=mid.getRightEdge(pairedRight[i].id, 1);
			if(!validForwardEdge(left[i], mid, entries[i])
					|| !validForwardEdge(mid, pairedRight[i], exits[i])){
				restoreOrientations(left, right, flippedLeft, flippedRight);
				return false;
			}
			paths[i]=makePath(left[i], entries[i], mid, exits[i], pairedRight[i]);
		}

		final long support0=c.support[pairing], support1=c.support[3-pairing];
		final float fraction0=(float)(support0/((double)support0+support1));
		final float[] coverage={pathCoverage(left[0], entries[0], mid, exits[0], pairedRight[0],
				paths[0].length, fraction0),
				pathCoverage(left[1], entries[1], mid, exits[1], pairedRight[1],
						paths[1].length, 1-fraction0)};

		for(int i=0; i<2; i++){
			final Contig l=left[i], r=pairedRight[i];
			redirectExterior(r.id, l.id, true);
			final ArrayList<Edge> exterior=r.rightEdges;
			if(exterior!=null){for(Edge e : exterior){e.origin=l.id;}}
			l.rightEdges=exterior;
			r.setUsed(destMap, contigs);
			setProduct(l, r, mid, paths[i], coverage[i]);
		}
		mid.setUsed(destMap, contigs);
		return true;
	}

	private void restoreOrientations(final Contig[] left, final Contig[] right,
			final boolean[] flippedLeft, final boolean[] flippedRight){
		for(int i=1; i>=0; i--){if(flippedRight[i]){right[i].flip(destMap.get(right[i].id));}}
		for(int i=1; i>=0; i--){if(flippedLeft[i]){left[i].flip(destMap.get(left[i].id));}}
	}

	private boolean exteriorSafe(final Contig node, final boolean exteriorRight, final Candidate x){
		final ArrayList<Edge> exterior=exteriorRight ? node.rightEdges : node.leftEdges;
		if(exterior==null){return true;}
		for(Edge e : exterior){
			if(e==null || e.origin!=node.id || e.sourceRight()!=exteriorRight
					|| inCandidate(e.destination, x)){return false;}
			final Contig target=contigs.get(e.destination);
			if(target.used() || target.associate()){return false;}
		}
		return true;
	}

	private static boolean inCandidate(final int id, final Candidate c){
		return id==c.middle || id==c.left[0] || id==c.left[1] || id==c.right[0] || id==c.right[1];
	}

	private boolean validForwardEdge(final Contig source, final Contig target, final Edge e){
		if(e==null || e.origin!=source.id || e.destination!=target.id || e.orientation!=1
				|| source.length()<k || target.length()<k || e.length<1
				|| e.bases==null || e.bases.length!=e.length){return false;}
		if(e.length<=k){
			final int shared=k-e.length;
			for(int i=0; i<shared; i++){
				if(source.bases[source.length()-shared+i]!=target.bases[i]){return false;}
			}
			for(int i=0; i<e.length; i++){
				if(e.bases[i]!=target.bases[shared+i]){return false;}
			}
		}else{
			for(int i=0; i<k; i++){
				if(e.bases[e.length-k+i]!=target.bases[i]){return false;}
			}
		}
		return true;
	}

	private byte[] makePath(final Contig left, final Edge entry, final Contig mid,
			final Edge exit, final Contig right){
		bb.clear();
		bb.append(left.bases).append(entry.bases);
		for(int i=k; i<mid.bases.length; i++){bb.append(mid.bases[i]);}
		bb.append(exit.bases);
		for(int i=k; i<right.bases.length; i++){bb.append(right.bases[i]);}
		return bb.toBytes();
	}

	private float pathCoverage(final Contig left, final Edge entry, final Contig mid,
			final Edge exit, final Contig right, final int pathLength, final float middleFraction){
		final double sum=left.coverage*(left.length()-k+1)
				+mid.coverage*(mid.length()-k+1)*middleFraction
				+right.coverage*(right.length()-k+1)
				+entry.depth*Tools.max(0, entry.length-1)+exit.depth*Tools.max(0, exit.length-1);
		return (float)(sum/(pathLength-k+1));
	}

	private static void setProduct(final Contig left, final Contig right, final Contig mid,
			final byte[] bases, final float coverage){
		left.bases=bases;
		left.coverage=coverage;
		left.minCov=Tools.min(left.minCov, mid.minCov, right.minCov);
		left.maxCov=Tools.max(left.maxCov, mid.maxCov, right.maxCov);
		left.rightCode=right.rightCode;
		left.rightRatio=right.rightRatio;
		left.rightBridgeEndpoint=right.rightBridgeEndpoint;
		left.name=null;
		left.tid=-1;
		left.gc=left.hh=left.caga=-1;
	}

	private void redirectExterior(final int from, final int to, final boolean destRight){
		final ArrayList<Edge> old=destMap.get(from);
		if(old==null){return;}
		ArrayList<Edge> retained=null;
		ArrayList<Edge> moved=destMap.get(to);
		for(Edge e : old){
			if(e.destRight()==destRight){
				e.destination=to;
				if(moved==null){moved=new ArrayList<Edge>();}
				moved.add(e);
			}else{
				if(retained==null){retained=new ArrayList<Edge>();}
				retained.add(e);
			}
		}
		if(retained==null){destMap.remove(from);}else{destMap.put(from, retained);}
		if(moved!=null){destMap.put(to, moved);}
	}

	private void rebuildDestMap(){
		destMap.clear();
		for(Contig c : contigs){
			if(c.used() || c.associate()){continue;}
			addInbound(c.leftEdges);
			addInbound(c.rightEdges);
		}
	}

	private void addInbound(final ArrayList<Edge> edges){
		if(edges==null){return;}
		for(Edge e : edges){
			ArrayList<Edge> list=destMap.get(e.destination);
			if(list==null){list=new ArrayList<Edge>(); destMap.put(e.destination, list);}
			list.add(e);
		}
	}

	private static int indexOfEdge(final ArrayList<Edge> edges, final int destination){
		for(int i=0; i<edges.size(); i++){if(edges.get(i).destination==destination){return i;}}
		return -1;
	}

	private static boolean sameNeighbors(final int[] a, final int[] b){
		return (a[0]==b[0] && a[1]==b[1]) || (a[0]==b[1] && a[1]==b[0]);
	}

	private boolean validClosedNeighbor(final Contig mid, final Edge forward, final boolean midRight){
		if(forward==null || forward.origin!=mid.id || forward.sourceRight()!=midRight
				|| forward.destination<0 || forward.destination>=contigs.size()){return false;}
		final Contig neighbor=contigs.get(forward.destination);
		if(neighbor==null || neighbor.used() || neighbor.associate()){return false;}
		final boolean joinedRight=forward.destRight();
		final ArrayList<Edge> joined=joinedRight ? neighbor.rightEdges : neighbor.leftEdges;
		if(joined==null || joined.size()!=1 || countInbound(neighbor.id, joinedRight)!=1){return false;}
		final Edge reverse=joined.get(0);
		return reverse.origin==neighbor.id && reverse.destination==mid.id
				&& reverse.sourceRight()==joinedRight && reverse.destRight()==midRight
				&& reverse.length==forward.length && reverse.depth>0 && forward.depth>0;
	}

	private int countInbound(final int id, final boolean destRight){
		final ArrayList<Edge> inbound=destMap.get(id);
		if(inbound==null){return 0;}
		int count=0;
		for(Edge e : inbound){if(e.destRight()==destRight){count++;}}
		return count;
	}

	private static int[] neighborIDs(final ArrayList<Edge> edges){
		if(edges==null || edges.size()!=2){return null;}
		return new int[]{edges.get(0).destination, edges.get(1).destination};
	}

	private static int indexOf(final int[] array, final int value){
		return array[0]==value ? 0 : array[1]==value ? 1 : -1;
	}

	static final class Candidate {
		Candidate(int middle_, int[] left_, int[] right_){
			middle=middle_;
			left=left_;
			right=right_;
		}

		/** 0 pairs 00+11; 1 pairs 01+10; -1 declines. */
		int pairing(final int minSupport, final int maxNoise){
			final boolean diagonal=support[0]>=minSupport && support[3]>=minSupport
					&& support[1]<=maxNoise && support[2]<=maxNoise
					&& support[0]>support[1] && support[0]>support[2]
					&& support[3]>support[1] && support[3]>support[2];
			final boolean crossed=support[1]>=minSupport && support[2]>=minSupport
					&& support[0]<=maxNoise && support[3]<=maxNoise
					&& support[1]>support[0] && support[1]>support[3]
					&& support[2]>support[0] && support[2]>support[3];
			return diagonal==crossed ? -1 : diagonal ? 0 : 1;
		}

		long totalSupport(){return support[0]+support[1]+support[2]+support[3];}

		final int middle;
		final int[] left, right;
		final long[] support=new long[4];
	}

	private final ArrayList<Contig> contigs;
	private final HashMap<Integer, ArrayList<Edge>> destMap;
	private final Candidate[] candidates;
	private final int k;
	private final int minSupport;
	private final int maxNoise;
	private final ByteBuilder bb=new ByteBuilder();
}
