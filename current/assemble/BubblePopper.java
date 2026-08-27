package assemble;

import java.util.ArrayList;
import java.util.HashMap;

import shared.Tools;
import structures.ByteBuilder;

/**
 * Identifies and resolves bubble structures in assembly graphs.
 * Bubbles are short alternative paths between two nodes that represent
 * sequencing errors or minor sequence variants. This class collapses
 * these redundant paths to simplify the assembly graph.
 *
 * @author Brian Bushnell
 */
public class BubblePopper {
	
	/**
	 * Constructs a BubblePopper for the given assembly graph.
	 * @param allContigs_ Complete list of contigs in the assembly
	 * @param destMap_ Mapping from contig IDs to their incoming edges
	 * @param kbig_ K-mer size used in the assembly
	 */
	BubblePopper(ArrayList<Contig> allContigs_, HashMap<Integer, ArrayList<Edge>> destMap_, int kbig_){
		this(allContigs_, destMap_, kbig_, null);
	}

	/**
	 * Constructs a bubble popper with Tadpole's configured kmer-count error classifier.
	 * A null classifier disables destructive indirect bubble cleanup while leaving
	 * lossless direct unitig merging available.
	 */
	BubblePopper(ArrayList<Contig> allContigs_, HashMap<Integer, ArrayList<Edge>> destMap_, int kbig_,
			ErrorClassifier errorClassifier_){
		allContigs=allContigs_;
		destMap=destMap_;
		kbig=kbig_;
		minLen=2*kbig-1;
		errorClassifier=errorClassifier_;
	}
	
	/**
	 * Expands a contig by finding and merging compatible adjacent paths.
	 * Attempts both direct merging of simple paths and indirect bubble popping.
	 * Processes the contig in both directions by flipping when needed.
	 *
	 * @param c The contig to expand from
	 * @return Number of expansions performed
	 */
	int expand(Contig c) {
		assert(!validateGraph || validate(c));
//		if(true) {return 0;}
		if(verbose){System.err.println("\n\n*expand: "+c.name()+"\n");}
		assert(!c.used());
		center=c;
		int count=0;
		
//		debranch(c);
		
		while(popDirect && expandRightSimple()){count++;}
		while((popIndirect || unzipBubbles) && center.rightForwardBranch() && expandRight()){
			count++;
			while(popDirect && expandRightSimple()){count++;}
		}
		
		if((popDirect && (crossKMerge || (center.leftCode!=Tadpole.LOOP && center.leftCode!=Tadpole.DEAD_END)) && center.leftEdges!=null)
				|| ((popIndirect || unzipBubbles) && center.leftForwardBranch())){
			final int countBeforeFlip=count;
			center.flip(destMap.get(center.id));
			while(popDirect && expandRightSimple()){count++;}
			while((popIndirect || unzipBubbles) && center.rightForwardBranch() && expandRight()) {
				count++;
				while(popDirect && expandRightSimple()){count++;}
			}
			if(unzipBubbles && count==countBeforeFlip){center.flip(destMap.get(center.id));}
//			center.flip(destMap.get(center.id));
		}
		assert(!validateGraph || validate(c));
		return count;
	}
	
	/**
	 * Removes dead-end branches from both sides of a contig.
	 * Only operates when the debranch flag is enabled.
	 * @param c The contig to debranch
	 */
	void debranch(Contig c){
		assert(debranch);
		debranchRight(c);
		debranchLeft(c);
	}
	
	/**
	 * Removes dead-end branches from the right side of a contig.
	 * Skips processing if there are fewer than 2 outbound edges or if it's a loop.
	 * @param c The contig to debranch
	 */
	private void debranchRight(Contig c){
		if(c.rightEdgeCount()<2 || c.rightCode==Tadpole.LOOP){return;}
		debranch(c, c.rightEdges);
	}
	
	/**
	 * Removes dead-end branches from the left side of a contig.
	 * Skips processing if there are fewer than 2 outbound edges or if it's a loop.
	 * @param c The contig to debranch
	 */
	private void debranchLeft(Contig c){
		if(c.leftEdgeCount()<2 || c.leftCode==Tadpole.LOOP){return;}
		debranch(c, c.leftEdges);
	}
	
	/**
	 * Determines if the left side of a contig is a dead end.
	 * Considers length, edge codes, and incoming edge patterns.
	 * Contigs longer than 400bp are never considered dead ends.
	 *
	 * @param c The contig to test
	 * @return true if the left side is a dead end
	 */
	private boolean isDeadEndLeft(Contig c){
		if(c.length()>400){return false;}
		
		if(c.leftCode==Tadpole.DEAD_END){return true;}
		if(c.leftEdgeCount()>0){return false;}
		
		ArrayList<Edge> inbound=destMap.get(c.id);
		if(inbound==null){return true;}
		for(Edge e : inbound){
			if(!e.destRight()){return false;}
		}
		return true;
	}
	
	/**
	 * Determines if the right side of a contig is a dead end.
	 * Considers length, edge codes, and incoming edge patterns.
	 * Contigs longer than 400bp are never considered dead ends.
	 *
	 * @param c The contig to test
	 * @return true if the right side is a dead end
	 */
	private boolean isDeadEndRight(Contig c){
		if(c.length()>400){return false;}
		
		if(c.rightCode==Tadpole.DEAD_END){return true;}
		if(c.rightEdgeCount()>0){return false;}
		
		ArrayList<Edge> inbound=destMap.get(c.id);
		if(inbound==null){return true;}
		for(Edge e : inbound){
			if(e.destRight()){return false;}
		}
		return true;
	}
	
	/**
	 * Removes dead-end branches from the specified outbound edge set.
	 * Identifies dead ends, optionally keeps the longest branch if all are dead,
	 * and truncates the selected branches.
	 *
	 * @param source The source contig
	 * @param outbound The outbound edges to evaluate
	 */
	private void debranch(Contig source, ArrayList<Edge> outbound){
		if(outbound==null || outbound.size()<2) {return;}
		int deadEnds=0;
		int longest=0;
		int deepest=0;
		float mult=0;
		for(Edge e : outbound){
			if(e.destination==e.origin){return;}
			Contig c=allContigs.get(e.destination);
			if(!c.used()){
				longest=Tools.max(longest, c.length());
				deepest=Tools.max(deepest, e.depth);
				mult=Tools.max(mult, c.length()*c.coverage);
			}
			boolean dead=c.used() || e.destRight() ? isDeadEndLeft(c) : isDeadEndRight(c);
			if(dead){deadEnds++;}
		}
		if(deadEnds==0){return;}
		
		boolean keepLongest=(deadEnds==outbound.size());

		ArrayList<Edge> toTruncate=new ArrayList<Edge>(deadEnds);
		for(Edge e : outbound){
			Contig c=allContigs.get(e.destination);
			boolean dead=c.used() || e.destRight() ? isDeadEndLeft(c) : isDeadEndRight(c);
			
			if(dead && (!keepLongest || c.length()*c.coverage<mult)){toTruncate.add(e);}
		}
		
		int removed=0;
		for(Edge e : toTruncate){
			Contig c=allContigs.get(e.destination);
			truncate(e, source, c);
			removed++;
		}
		assert(outbound.size()>=1);
	}
	
	/**
	 * Removes an edge connection between two contigs.
	 * Updates both the source and destination contig edge lists
	 * and the global destination mapping.
	 *
	 * @param e The edge to remove
	 * @param from The source contig
	 * @param to The destination contig
	 */
	void truncate(Edge e, Contig from, Contig to){
		branchesRemoved++;
		{
			if(e.sourceRight()){
				from.removeRightEdge(to.id, e.destRight());
			}else{
				from.removeLeftEdge(to.id, e.destRight());
			}
			ArrayList<Edge> inbound=destMap.get(to.id);
			inbound=Contig.removeEdge(to.id, e.destRight(), inbound);
			if(inbound==null){destMap.remove(to.id);}
		}
		{
			if(e.destRight()){
				to.removeRightEdge(from.id, e.sourceRight());
			}else{
				to.removeLeftEdge(from.id, e.sourceRight());
			}
			ArrayList<Edge> inbound=destMap.get(from.id);
			inbound=Contig.removeEdge(from.id, e.sourceRight(), inbound);
			if(inbound==null){destMap.remove(from.id);}
		}
	}
	
	/**
	 * Attempts to expand the center contig rightward through a simple path.
	 * A simple path has exactly one outbound edge leading to a node with
	 * at most one return edge to the center. Validates connectivity patterns
	 * before merging.
	 *
	 * @return true if expansion was successful
	 */
	private boolean expandRightSimple(){
		assert(!validateGraph || validate(center));
		ArrayList<Edge> outbound=center.rightEdges;
		if(outbound==null || center.rightCode==Tadpole.LOOP || outbound.size()>1){return false;}
		Edge leftEdge=outbound.get(0);
		assert(leftEdge.destination<allContigs.size()) : "\n"+leftEdge.toString()+", "+center.used()+", "+center.associate()+"\n"+center.name()+"\n"+allContigs.size();
		
//		 : "\ncenter="+center.name2()+"\ndest="+dest.name2()+"\nc="+c.name2()+"\nother="+other.name2()+"\ne="+e;
		
		dest=allContigs.get(leftEdge.destination);
		
		if(dest.used() || dest==center){return false;}
		if(crossKMerge && crossKMaxDepthRatio>0){
			final float maxFlank=Tools.max(center.coverage, dest.coverage);
			if(maxFlank>0 && leftEdge.depth>maxFlank*crossKMaxDepthRatio){return false;}
		}
		ArrayList<Edge> outboundRight=(leftEdge.destRight() ? dest.rightEdges : dest.leftEdges);
		int rightCode=leftEdge.destRight() ? dest.rightCode : dest.leftCode;
		
		if(rightCode==Tadpole.LOOP){return false;}
		
//		if(outboundRight==null || outboundRight.size()>1){return false;}
		if(crossKMerge && (outboundRight==null || outboundRight.size()!=1 || outboundRight.get(0).destination!=center.id)){
			return false;
		}else if(outboundRight==null) {
			//do nothing
		}else if(outboundRight.size()>1){
			return false;
		}else if(outboundRight.get(0).destination!=center.id){
			return false;
		}
		
		int inbound=countInbound(center.id, true);
//		if(inbound==null || inbound.size()>1){return false;}
		if(crossKMerge ? inbound!=1 : inbound>1){return false;}
		
		int inboundRight=countInbound(dest.id, leftEdge.destRight());
		if(crossKMerge ? inboundRight!=1 : inboundRight>1){return false;}
		
		if(leftEdge.destRight()){
			dest.flip(destMap.get(dest.id));
		}
		assert(!validateGraph || validate(center));
		assert(!validateGraph || validate(dest));
		
		if(crossKMerge){
			System.err.println("Cross-k merge: left="+center.name()+", right="+dest.name()+
					", bridge="+leftEdge.length+", depth="+leftEdge.depth+", orientation="+leftEdge.orientation);
		}
		return merge(center, dest, leftEdge);
	}
	
	/**
	 * Counts incoming edges to a contig from a specific side.
	 * @param id The contig ID to check
	 * @param destRight Whether to count edges incoming to the right side
	 * @return Number of matching incoming edges
	 */
	private int countInbound(int id, boolean destRight){
		int count=0;
		ArrayList<Edge> inbound=destMap.get(id);
		if(inbound==null){return 0;}
		for(Edge e : inbound){
			if(e.destRight()==destRight){
				count++;
			}
		}
		return count;
	}

	/**
	 * Gets incoming edges to a contig from a specific side.
	 * @param id The contig ID to check
	 * @param destRight Whether to get edges incoming to the right side
	 * @return List of matching incoming edges, or null if none found
	 */
	private ArrayList<Edge> getInbound(int id, boolean destRight){
		ArrayList<Edge> inbound=destMap.get(id);
		if(inbound==null){return null;}
		ArrayList<Edge> inboundSide=new ArrayList<Edge>(inbound.size());
		for(Edge e : inbound){
			if(e.destRight()==destRight){
				inboundSide.add(e);
			}
		}
		return inboundSide.isEmpty() ? null : inboundSide;
	}
	
	/**
	 * Attempts to expand the center contig rightward through bubble resolution.
	 * Identifies bubble structures with multiple intermediate nodes that
	 * converge to a mutual destination. Validates bubble topology and
	 * performs the bubble popping operation.
	 *
	 * @return true if bubble was successfully popped
	 */
	private boolean expandRight(){
		assert(!validateGraph || validate(center));
		//Reset shared state
		dest=null;
		lastMutualDest=-1;
		lastMutualDestOrientation=-1;
		
		if(verbose){System.err.println("expandRight: "+center.name());}
		
		if(!center.rightForwardBranch() || center.rightEdgeCount()<1){
			if(verbose){System.err.println("Returned because not forward branch or no edges.");}
			return false;
		}
		ArrayList<Edge> outbound=center.rightEdges;
		final Edge leftMidEdge=findRepresentativeMidEdge(outbound);
		
		if(leftMidEdge==null){
			if(verbose){System.err.println("No leftMidEdge.");}
			return false;
		}
		final Contig mid=allContigs.get(leftMidEdge.destination);
		if(mid==null || mid.length()<minLen){
			if(verbose){System.err.println("No mid, or mid too short ("+(mid==null ? 0 : mid.length())+").");}
			return false;
		}
		
		if(verbose){System.err.println("\nFinding mutualDest for center node.");}
		final int mutualDest=findMutualDest(outbound);
		final int mutualDestOrientation=lastMutualDestOrientation;
		final boolean mutualDestRight=((lastMutualDestOrientation&2)==2);
		if(verbose){System.err.println("mutualDest="+mutualDest+", mutualDestOrientation="+mutualDestOrientation+", mutualDestRight="+mutualDestRight);}
		
		if(mutualDest<0 || mutualDestOrientation<0){
			if(verbose){System.err.println("Bad mutual destination.");}
			return false;
		}
		//At this point we are fairly confident everything is safe, but still need to run more tests.
		
		dest=allContigs.get(mutualDest);
		if(dest.used()){return false;}
//		assert(!dest.used());//This happened; not sure how
		if(dest.id==center.id){
			if(verbose){System.err.println("Self mutual destination.");}
			return false;
		}
		
		if(mutualDestRight && !dest.rightForwardBranch()){
			if(verbose){System.err.println("Mutual destination does not have a right F-branch: "+dest.name());}
			return false;
		}
		if(!mutualDestRight && !dest.leftForwardBranch()){
			if(verbose){System.err.println("Mutual destination does not have a left F-branch: "+dest.name());}
			return false;
		}
		final ArrayList<Edge> destOutbound=mutualDestRight ? dest.rightEdges : dest.leftEdges;
		if(destOutbound==null){
			if(verbose){System.err.println("No dest outbound edges.");}
			return false;
		}

		if(verbose){System.err.println("\nFinding mutualDest for right node.");}
		final int mutualDest2=findMutualDest(destOutbound);
		if(mutualDest2<0 || mutualDest2!=center.id){
			if(verbose){System.err.println("MutualDest2 is not the correct origin: "+center.id+"!="+mutualDest2+"; "+dest.name());}
			return false;
		}
		
		final ArrayList<Contig> flippedMids=new ArrayList<Contig>(outbound.size());
		ArrayList<Contig> midNodes=fetchMidNodes(outbound, true, flippedMids);
		if(midNodes==null){
			restoreOrientation(flippedMids, false);
			return false;
		}
		//Now we have all intermediate nodes, which are flipped into the correct orientation.
		
		if(!midNodesConcur(midNodes)){
			if(verbose){System.err.println("Mid nodes do not concur.");}
			restoreOrientation(flippedMids, false);
			return false;
		}
		
		//At this point all tests have been run and we are confident that this is a simple bubble.
		if(mutualDestRight){
			if(verbose){System.err.println("Flipping dest node.");}
			dest.flip(destMap.get(dest.id));
		}
		
		//Now the destination is also flipped into the correct orientation.
		final Edge rightMidEdge=mid.getRightEdge(dest.id, 1);
		if(rightMidEdge==null){
			restoreOrientation(flippedMids, mutualDestRight);
			return false;
		}

		final ArrayList<Contig> errorMids=findErrorMids(outbound, midNodes, mid,
				leftMidEdge, rightMidEdge);
		if(errorMids==null){
			if(verbose){System.err.println("Bubble could not be classified safely.");}
			restoreOrientation(flippedMids, mutualDestRight);
			return false;
		}
		if(errorMids.isEmpty()){
			if(verbose){System.err.println("No alternate paths classified as errors.");}
			if(unzipBubbles && unzipTrueBubble(center, dest, mid, leftMidEdge, rightMidEdge, midNodes)){
				return true;
			}
			restoreOrientation(flippedMids, mutualDestRight);
			return false;
		}
		if(!popIndirect){
			restoreOrientation(flippedMids, mutualDestRight);
			return false;
		}
		if(errorMids.size()<midNodes.size()-1){
			if(verbose){System.err.println("Pruning "+errorMids.size()+" error arms while preserving true alternatives.");}
			return pruneErrorMids(errorMids);
		}

		if(verbose){System.err.println("Popping bubble between "+center.id+" and "+dest.id);}
		return pop(center, dest, mid, leftMidEdge, rightMidEdge, midNodes);
		
	}

	/** Restores orientation after a declined candidate, in reverse mutation order. */
	private void restoreOrientation(ArrayList<Contig> flippedMids, boolean destFlipped){
		if(!unzipBubbles){return;}
		if(destFlipped){dest.flip(destMap.get(dest.id));}
		for(int i=flippedMids.size()-1; i>=0; i--){
			Contig c=flippedMids.get(i);
			c.flip(destMap.get(c.id));
		}
	}

	/** Returns alternate mid nodes whose entry and exit branch counts both indicate error. */
	private ArrayList<Contig> findErrorMids(ArrayList<Edge> entryEdges, ArrayList<Contig> midNodes,
			Contig representative, Edge representativeEntry, Edge representativeExit){
		if(errorClassifier==null || midNodes.size()<2){return null;}
		ArrayList<Contig> errors=new ArrayList<Contig>(midNodes.size()-1);
		for(Contig c : midNodes){
			if(c==representative){continue;}
			Edge entry=findEdge(entryEdges, c.id);
			Edge exit=c.getRightEdge(dest.id, 1);
			if(entry==null || exit==null){return null;}
			if(isErrorArm(representativeEntry, representativeExit, entry, exit)){
				errors.add(c);
			}
		}
		return errors;
	}

	/** Finds the unique edge to a destination in an already topology-validated edge list. */
	private static Edge findEdge(ArrayList<Edge> edges, int destination){
		Edge found=null;
		for(Edge e : edges){
			if(e.destination==destination){
				if(found!=null){return null;}
				found=e;
			}
		}
		return found;
	}

	/** Requires the alternate to be lower and error-like at both bubble boundaries. */
	private boolean isErrorArm(Edge representativeEntry, Edge representativeExit,
			Edge alternateEntry, Edge alternateExit){
		if(representativeEntry.depth<=alternateEntry.depth || representativeExit.depth<=alternateExit.depth){
			return false;
		}
		return errorClassifier.isError(representativeEntry.depth, alternateEntry.depth)
				&& errorClassifier.isError(representativeExit.depth, alternateExit.depth);
	}

	/** Detaches only error-classified arms, retaining comparable-depth biological alternatives. */
	private boolean pruneErrorMids(ArrayList<Contig> errorMids){
		assert(!errorMids.isEmpty());
		for(Contig c : errorMids){
			assert(!c.used() && !c.associate());
			c.setAssociate(destMap, allContigs);
			assert(!validateGraph || validate(c));
		}
		expansions++;
		contigsAbsorbed+=errorMids.size();
		assert(!validateGraph || validate(center));
		assert(!validateGraph || validate(dest));
		return true;
	}

	/**
	 * Linearizes an isolated two-arm biological bubble into two terminal products.
	 * This deliberately declines bubbles with any exterior connectivity; duplicating
	 * shared flanks inside a larger graph would otherwise invent unsupported phasing.
	 */
	private boolean unzipTrueBubble(Contig left, Contig right, Contig representative,
			Edge representativeEntry, Edge representativeExit, ArrayList<Contig> midNodes){
		if(midNodes.size()!=2 || left==right || left.used() || left.associate()
				|| right.used() || right.associate()){return false;}
		if(left.leftCode!=Tadpole.DEAD_END || right.rightCode!=Tadpole.DEAD_END){return false;}
		if(left.leftEdgeCount()!=0 || right.rightEdgeCount()!=0){return false;}
		if(countInbound(left.id, false)!=0 || countInbound(right.id, true)!=0){return false;}
		if(left.rightEdgeCount()!=2 || right.leftEdgeCount()!=2){return false;}

		final int representativeIndex=midNodes.indexOf(representative);
		if(representativeIndex<0){return false;}
		final Contig alternate=midNodes.get(1-representativeIndex);
		if(alternate==representative || representative.used() || representative.associate()
				|| alternate.used() || alternate.associate()){return false;}

		final Edge alternateEntry=findEdge(left.rightEdges, alternate.id);
		final Edge alternateExit=alternate.getRightEdge(right.id, 1);
		if(alternateEntry==null || alternateExit==null){return false;}
		if(isErrorArm(alternateEntry, alternateExit, representativeEntry, representativeExit)){return false;}

		final Edge representativeBack=representative.getLeftEdge(left.id, 2);
		final Edge alternateBack=alternate.getLeftEdge(left.id, 2);
		final Edge rightToRepresentative=findEdge(right.leftEdges, representative.id);
		final Edge rightToAlternate=findEdge(right.leftEdges, alternate.id);
		if(!validBubbleArm(left, right, representative, representativeEntry, representativeBack,
				representativeExit, rightToRepresentative)
				|| !validBubbleArm(left, right, alternate, alternateEntry, alternateBack,
						alternateExit, rightToAlternate)){return false;}

		if(!inboundExactly(representative.id, representativeEntry, rightToRepresentative)
				|| !inboundExactly(alternate.id, alternateEntry, rightToAlternate)
				|| !inboundExactly(left.id, representativeBack, alternateBack)
				|| !inboundExactly(right.id, representativeExit, alternateExit)){return false;}

		final byte[] representativePath=makePath(left, representativeEntry, representative,
				representativeExit, right);
		final byte[] alternatePath=makePath(left, alternateEntry, alternate, alternateExit, right);
		if(representativePath==null || alternatePath==null){return false;}

		final long entryDepthSum=(long)representativeEntry.depth+alternateEntry.depth;
		final long exitDepthSum=(long)representativeExit.depth+alternateExit.depth;
		if(entryDepthSum<1 || exitDepthSum<1){return false;}
		final float representativeCoverage=pathCoverage(left, representativeEntry, representative,
				representativeExit, right, representativePath.length,
				representativeEntry.depth/(float)entryDepthSum, representativeExit.depth/(float)exitDepthSum);
		final float alternateCoverage=pathCoverage(left, alternateEntry, alternate,
				alternateExit, right, alternatePath.length,
				alternateEntry.depth/(float)entryDepthSum, alternateExit.depth/(float)exitDepthSum);
		final int leftMin=left.minCov, leftMax=left.maxCov;
		final int rightMin=right.minCov, rightMax=right.maxCov;

		representative.setUsed(destMap, allContigs);
		alternate.setUsed(destMap, allContigs);
		assert(left.leftEdges==null && left.rightEdges==null);
		assert(right.leftEdges==null && right.rightEdges==null);

		setTerminalProduct(left, representativePath, representativeCoverage,
				Tools.min(leftMin, rightMin, representative.minCov),
				Tools.max(leftMax, rightMax, representative.maxCov));
		setTerminalProduct(right, alternatePath, alternateCoverage,
				Tools.min(leftMin, rightMin, alternate.minCov),
				Tools.max(leftMax, rightMax, alternate.maxCov));

		expansions++;
		trueBubblesUnzipped++;
		contigsAbsorbed+=2;
		assert(!validateGraph || validate(left));
		assert(!validateGraph || validate(right));
		assert(!validateGraph || validate(representative));
		assert(!validateGraph || validate(alternate));
		return true;
	}

	/** Requires one exact reciprocal, sequence-consistent path through an arm. */
	private boolean validBubbleArm(Contig left, Contig right, Contig mid, Edge entry,
			Edge back, Edge exit, Edge reverseExit){
		if(mid==left || mid==right || mid.length()<kbig){return false;}
		if(mid.leftEdgeCount()!=1 || mid.rightEdgeCount()!=1){return false;}
		if(entry==null || back==null || exit==null || reverseExit==null){return false;}
		if(entry.origin!=left.id || entry.destination!=mid.id || entry.orientation!=1){return false;}
		if(back.origin!=mid.id || back.destination!=left.id || back.orientation!=2){return false;}
		if(exit.origin!=mid.id || exit.destination!=right.id || exit.orientation!=1){return false;}
		if(reverseExit.origin!=right.id || reverseExit.destination!=mid.id || reverseExit.orientation!=2){return false;}
		if(entry.depth<1 || back.depth<1 || exit.depth<1 || reverseExit.depth<1
				|| entry.length!=back.length || exit.length!=reverseExit.length){return false;}
		return validForwardEdge(left, mid, entry) && validForwardEdge(mid, right, exit);
	}

	/** Validates the sequence overlap represented by an oriented forward edge. */
	private boolean validForwardEdge(Contig source, Contig target, Edge e){
		if(source.length()<kbig || target.length()<kbig || e.length<1
				|| e.bases==null || e.bases.length!=e.length){return false;}
		if(e.length<=kbig){
			final int shared=kbig-e.length;
			for(int i=0; i<shared; i++){
				if(source.bases[source.length()-shared+i]!=target.bases[i]){return false;}
			}
			for(int i=0; i<e.length; i++){
				if(e.bases[i]!=target.bases[shared+i]){return false;}
			}
		}else{
			for(int i=0; i<kbig; i++){
				if(e.bases[e.length-kbig+i]!=target.bases[i]){return false;}
			}
		}
		return true;
	}

	/** The destination map must contain exactly these two distinct inbound edges. */
	private boolean inboundExactly(int id, Edge a, Edge b){
		if(a==null || b==null || a==b){return false;}
		final ArrayList<Edge> inbound=destMap.get(id);
		return inbound!=null && inbound.size()==2 && inbound.contains(a) && inbound.contains(b);
	}

	/** Builds one oriented left-arm-right product using ordinary exact-k joins. */
	private byte[] makePath(Contig left, Edge entry, Contig mid, Edge exit, Contig right){
		bb.clear();
		bb.append(left.bases);
		bb.append(entry.bases);
		for(int i=kbig; i<mid.bases.length; i++){bb.append(mid.bases[i]);}
		bb.append(exit.bases);
		for(int i=kbig; i<right.bases.length; i++){bb.append(right.bases[i]);}
		return bb.length()>=kbig ? bb.toBytes() : null;
	}

	/** Returns span-weighted kmer coverage for one unzipped path. */
	private float pathCoverage(Contig left, Edge entry, Contig mid, Edge exit, Contig right,
			int pathLength, float leftFraction, float rightFraction){
		final double sum=left.coverage*(left.length()-kbig+1)*leftFraction
				+mid.coverage*(mid.length()-kbig+1)+right.coverage*(right.length()-kbig+1)*rightFraction
				+entry.depth*Tools.max(0, entry.length-1)+exit.depth*Tools.max(0, exit.length-1);
		return (float)(sum/(pathLength-kbig+1));
	}

	/** Reinitializes a reused flank contig as an edge-free linear product. */
	private static void setTerminalProduct(Contig c, byte[] bases, float coverage, int minCov, int maxCov){
		c.bases=bases;
		c.coverage=coverage;
		c.minCov=minCov;
		c.maxCov=maxCov;
		c.name=null;
		c.tid=-1;
		c.gc=c.hh=c.caga=-1;
		c.leftCode=c.rightCode=Tadpole.DEAD_END;
		c.leftRatio=c.rightRatio=0;
		c.leftBridgeEndpoint=c.rightBridgeEndpoint=true;
		c.leftEdges=c.rightEdges=null;
	}
	
	/**
	 * Pops a bubble by merging the left and right nodes through the best path.
	 * Concatenates sequence data, redirects edges, marks intermediate nodes
	 * as used, and updates coverage statistics.
	 *
	 * @param left The left (source) contig
	 * @param right The right (destination) contig
	 * @param mid The representative middle contig
	 * @param leftMidEdge Edge from left to middle
	 * @param rightMidEdge Edge from middle to right
	 * @param midNodes All intermediate nodes in the bubble
	 * @return true if pop operation completed successfully
	 */
	private boolean pop(Contig left, Contig right, Contig mid, Edge leftMidEdge, Edge rightMidEdge, ArrayList<Contig> midNodes){
		assert(!validateGraph || validate(left));
		assert(!validateGraph || validate(right));
		assert(!validateGraph || validate(mid));
		for(Contig c : midNodes){
			assert(!c.used() && !c.associate());
			assert(c!=right && c!=left);
			assert(!validateGraph || validate(c));
		}
		
		bb.clear();
		final int originalLeftLength=left.length();
		
		//Append the three contigs as two ordinary k-exact joins.  This form is
		//also correct when either connecting edge spans more than one base.
		bb.append(left.bases);
		if(leftMidEdge.bases!=null){bb.append(leftMidEdge.bases);}
		for(int i=kbig; i<mid.bases.length; i++){bb.append(mid.bases[i]);}
		if(rightMidEdge.bases!=null){bb.append(rightMidEdge.bases);}
		for(int i=kbig; i<right.bases.length; i++){bb.append(right.bases[i]);}
		left.bases=bb.toBytes();
		
		//Cleanup
		left.rightEdges.clear();
		if(right.rightEdgeCount()>0){
			for(Edge e : right.rightEdges){
				e.origin=left.id;
				left.addRightEdge(e);
			}
		}else{left.rightEdges=null;}
		
//		ArrayList<Edge> inbound=destMap.get(right.id);
//		if(inbound!=null){
//			for(Edge e : inbound){
//				e.destination=left.id;
//			}
//		}
		redirectEdges(right.id, left.id, true);
		right.setUsed(destMap, allContigs);//Don't do this until redirection is finished!
		for(Contig c : midNodes){//Don't do this until redirection is finished!
			assert(!c.used()) : "\nleft="+left.name2()+"\nright="+right.name2()+"\nc="+c.name2()+"\n\n"+midNodes;
			if(c==mid){c.setUsed(destMap, allContigs);}
			else{c.setAssociate(destMap, allContigs);}
		}

		left.maxCov=Tools.max(left.maxCov, right.maxCov, mid.maxCov);
		//[assemble/BubblePopper#001 FIXED 2026-07-02] minCov was computed from the maxCov FIELDS (min of two maximums) and
		//omitted mid; corrected to the minCov fields incl. mid. Twin of merge() below. NOTE: Contig.minCov/maxCov are only
		//printed in the FASTA header (Contig.java:101-102) and are never initialized at build time (0 for un-merged contigs),
		//so today this is a latent/correctness fix; making min=/max= actually meaningful needs build-side coverage init (a
		//separate feature, NOT done here - awaiting Brian's call).
		left.minCov=Tools.min(left.minCov, right.minCov, mid.minCov);
		left.rightCode=right.rightCode;
		left.rightRatio=right.rightRatio;
		left.rightBridgeEndpoint=right.rightBridgeEndpoint;
		final double coverageSum;
		if(crossKMerge){
			coverageSum=left.coverage*originalLeftLength+right.coverage*right.length();
			left.coverage=(float)(coverageSum/(originalLeftLength+right.length()));
		}else{
			coverageSum=left.coverage*(originalLeftLength-kbig+1)
					+mid.coverage*(mid.length()-kbig+1)
					+right.coverage*(right.length()-kbig+1)
					+leftMidEdge.depth*Tools.max(0, leftMidEdge.length-1)
					+rightMidEdge.depth*Tools.max(0, rightMidEdge.length-1);
			left.coverage=(float)(coverageSum/(left.length()-kbig+1));
		}
		
		expansions++;
		contigsAbsorbed+=(1+midNodes.size());
		
		if(isLoop(left)){
			left.leftCode=Tadpole.LOOP;
			left.rightCode=Tadpole.LOOP;
			left.leftBridgeEndpoint=left.rightBridgeEndpoint=false;
			left.removeAllEdges(destMap.get(left.id), allContigs);
		}
		
		if(verbose){
			System.err.println("*Result: "+center.name());
		}

		assert(!validateGraph || validate(left));
		assert(!validateGraph || validate(right));
		assert(!validateGraph || validate(mid));
		for(Contig c : midNodes){
			assert(!validateGraph || validate(c));
		}
		
		return true;
	}
	
	/**
	 * Redirects incoming edges from one contig to another.
	 * Updates the destination mapping to point edges from the source
	 * contig to the target contig for the specified orientation.
	 *
	 * @param from Source contig ID
	 * @param to Target contig ID
	 * @param destRight Whether to redirect edges to the right side
	 */
	private void redirectEdges(final int from, final int to, final boolean destRight){
		if(from==to){return;}
		ArrayList<Edge> inboundFrom=destMap.get(from);
		if(inboundFrom==null){return;}
		assert(inboundFrom.size()>0);
		destMap.remove(from);
		
		ArrayList<Edge> inboundTo=destMap.get(to);
		if(inboundTo==null){inboundTo=new ArrayList<Edge>(inboundFrom.size());}
		for(Edge e : inboundFrom){
			assert(e.destination==from);
			if(e.destRight()==destRight){
				e.destination=to;
				inboundTo.add(e);
			}
		}
		if(inboundTo.isEmpty()){inboundTo=null;}
		destMap.put(to, inboundTo);
	}
	
	/**
	 * Removes edges pointing to used or associated contigs.
	 * Cleans up both left and right edge lists.
	 * @param c The contig to clean
	 */
	void removeDeadEdges(Contig c){
		c.leftEdges=removeDeadEdges(c.leftEdges);
		c.rightEdges=removeDeadEdges(c.rightEdges);
	}
	
	/**
	 * Removes edges pointing to used or associated contigs from an edge list.
	 * Nullifies dead edges and condenses the list to remove gaps.
	 * @param edges The edge list to clean
	 * @return The cleaned edge list, or null if all edges were removed
	 */
	private ArrayList<Edge> removeDeadEdges(ArrayList<Edge> edges){
		if(edges==null){return null;}
		int removed=0;
		for(int i=0; i<edges.size(); i++){
			Edge e=edges.get(i);
			Contig c=allContigs.get(e.destination);
			if(c.used() || c.associate()){
				edges.set(i, null);
				removed++;
			}
		}
		if(removed>0){
			Tools.condenseStrict(edges);
		}
		return edges.isEmpty() ? null : edges;
	}
	
	/**
	 * Merges two contigs connected by a direct edge.
	 * Concatenates sequence data, transfers edges from right to left,
	 * updates coverage statistics, and marks the right contig as used.
	 *
	 * @param left The left (target) contig that will absorb the right contig
	 * @param right The right (source) contig to be merged
	 * @param leftEdge The connecting edge between them
	 * @return true if merge completed successfully
	 */
	private boolean merge(Contig left, Contig right, Edge leftEdge){

		assert(!validateGraph || validate(left));
		assert(!validateGraph || validate(right));
		
		bb.clear();
		final int originalLeftLength=left.length();
		
		//Append path
		bb.append(left.bases);
		if(leftEdge.bases!=null){bb.append(leftEdge.bases);}
		for(int i=kbig; i<right.bases.length; i++){bb.append(right.bases[i]);}
		left.bases=bb.toBytes();
		
		//Cleanup
		left.rightEdges.clear();
		if(right.rightEdgeCount()>0){
			for(Edge e : right.rightEdges){
				e.origin=left.id;
				left.addRightEdge(e);
			}
		}else{left.rightEdges=null;}
		
//		ArrayList<Edge> inbound=destMap.get(right.id);
//		if(inbound!=null){
//			for(Edge e : inbound){
//				e.destination=left.id;
//			}
//		}
		redirectEdges(right.id, left.id, true);
		right.setUsed(destMap, allContigs); //Don't do this until redirection is finished!

		left.maxCov=Tools.max(left.maxCov, right.maxCov);
		//[assemble/BubblePopper#001 FIXED 2026-07-02] TWIN of pop(): minCov was computed from the maxCov fields; corrected to
		//the minCov fields. (Same header-only/uninitialized caveat as noted in pop().)
		left.minCov=Tools.min(left.minCov, right.minCov);
		left.rightCode=right.rightCode;
		left.rightRatio=right.rightRatio;
		left.rightBridgeEndpoint=right.rightBridgeEndpoint;
		final double coverageSum;
		if(crossKMerge){
			coverageSum=left.coverage*originalLeftLength+right.coverage*right.length();
			left.coverage=(float)(coverageSum/(originalLeftLength+right.length()));
		}else{
			coverageSum=left.coverage*(originalLeftLength-kbig+1)+right.coverage*(right.length()-kbig+1)
					+leftEdge.depth*Tools.max(0, leftEdge.length-1);
			left.coverage=(float)(coverageSum/(left.length()-kbig+1));
		}
		
		if(isLoop(left)){
			left.leftCode=Tadpole.LOOP;
			left.rightCode=Tadpole.LOOP;
			left.leftBridgeEndpoint=left.rightBridgeEndpoint=false;
			left.removeAllEdges(destMap.get(left.id), allContigs);
		}
		
		expansions++;
		contigsAbsorbed++;
		
		if(verbose){
			System.err.println("*Result: "+center.name());
		}

		assert(!validateGraph || validate(left));
		assert(!validateGraph || validate(right));
		
		return true;
	}
	
	/**
	 * Determines if a contig forms a self-loop structure.
	 * Checks for loop codes or validates that all edges point back to itself
	 * with correct orientations.
	 *
	 * @param c The contig to test
	 * @return true if the contig is a loop
	 */
	boolean isLoop(Contig c){
		if(c.leftCode==Tadpole.LOOP && c.rightCode==Tadpole.LOOP){return true;}
		if(c.leftEdgeCount()!=1 || c.rightEdgeCount()!=1){return false;}
		for(Edge e : c.leftEdges){
			if(e.destination!=c.id || !e.destRight()){return false;}
		}
		for(Edge e : c.rightEdges){
			if(e.destination!=c.id || e.destRight()){return false;}
		}
		ArrayList<Edge> inbound=destMap.get(c.id);
		for(Edge e : inbound){
			if(e.origin!=c.id){return false;}
		}
		return true;
	}
	
	/**
	 * Validates the structural integrity of a contig and its edges.
	 * Performs extensive consistency checks on edge relationships,
	 * destination mappings, and contig states. Used for debugging
	 * and ensuring graph correctness.
	 *
	 * @param c The contig to validate
	 * @return true if all validation checks pass
	 */
	boolean validate(Contig c){
		assert(c!=null);
		assert(c.id>=0 && c.id<allContigs.size()) : c.id+", "+allContigs.size();
		assert(allContigs.get(c.id)==c) : c.id;
		final ArrayList<Edge> inbound=destMap.get(c.id);
		assert(inbound!=null || !destMap.containsKey(c.id)) : c.id;
		if(c.used() || c.associate()){
			assert(c.leftEdges==null);
			assert(c.rightEdges==null);
		}else{
			if(inbound!=null){
				for(Edge e : inbound){
					assert(e!=null);
					assert(e.destination==c.id) : e+", "+c.id;
					assert(e.origin>=0 && e.origin<allContigs.size()) : e+", "+allContigs.size();
					Contig other=allContigs.get(e.origin);
					if(other.used() || other.associate()){
						//ignore
					}else{
						ArrayList<Edge> sourceEdges=e.sourceRight() ? other.rightEdges : other.leftEdges;
						assert(sourceEdges!=null && sourceEdges.contains(e)) :
							"Inbound edge missing from source: c="+c.name2()+"\nother="+other.name2()+"\ne="+e;
					}
				}
			}
			assert(validateEdges(c, c.leftEdges, false));
			assert(validateEdges(c, c.rightEdges, true));
		}
		return true;
	}

	/** Validates one side of a live contig's outbound edge list. */
	private boolean validateEdges(Contig source, ArrayList<Edge> edges, boolean sourceRight){
		if(edges==null){return true;}
		assert(!edges.isEmpty());
		for(Edge e : edges){
			assert(e!=null);
			assert(e.origin==source.id) : e+", "+source.id;
			assert(e.sourceRight()==sourceRight) : e+", sourceRight="+sourceRight;
			assert(e.destination>=0 && e.destination<allContigs.size()) : e+", "+allContigs.size();
			Contig target=allContigs.get(e.destination);
			assert(!target.used() && !target.associate()) :
				"Live edge targets retired contig: source="+source.name2()+"\ntarget="+target.name2()+"\ne="+e;
			ArrayList<Edge> inbound=destMap.get(e.destination);
			assert(inbound!=null && inbound.contains(e)) :
				"Outbound edge missing from destMap: source="+source.name2()+"\ntarget="+target.name2()+"\ne="+e;
		}
		return true;
	}
	
	/**
	 * Selects the best representative edge from a set of outbound edges.
	 * Prioritizes edges to contigs meeting minimum length requirements,
	 * then selects based on depth and contig length.
	 *
	 * @param edges The edges to choose from
	 * @return The best representative edge, or null if none suitable
	 */
	private Edge findRepresentativeMidEdge(ArrayList<Edge> edges){
		//TODO: Probable bug - preferring any >=minLen arm can select a long low-depth error arm
		//over a shorter high-depth real arm.  popIndirect then compares only other arms against
		//that low representative and can miss the error; unzipTrueBubble has a symmetric guard,
		//but changing the destructive prune selector needs separate regression evidence.
		Edge midEdge=null;
		Contig mid=null;
		for(Edge e : edges){
			Contig c=allContigs.get(e.destination);
			if(midEdge==null){
				midEdge=e;
				mid=c;
			}else{
				if(mid.length()<minLen && c.length()>=minLen){
					midEdge=e;
					mid=c;
				}else if(c.length()>=minLen) {
					if(e.depth>midEdge.depth || (e.depth==midEdge.depth && c.length()>mid.length())) {
						midEdge=e;
						mid=c;
					}
				}
			}
		}
		return midEdge;
	}
	
	/**
	 * Validates that all intermediate nodes in a bubble have consistent
	 * connectivity patterns. Ensures all nodes connect to the same
	 * left and right destinations with no self-loops or external connections.
	 *
	 * @param midNodes The intermediate nodes to validate
	 * @return true if all nodes have consistent connectivity
	 */
	private boolean midNodesConcur(ArrayList<Contig> midNodes){
		int leftDest=-1;
		int rightDest=-1;
		for(Contig c : midNodes){
			if(c.leftEdges==null){
				if(verbose){System.err.println("No midnode left edges for "+c.name());}
				return false;
			}
			if(c.rightEdges==null){
				if(verbose){System.err.println("No midnode right edges for "+c.name());}
				return false;
			}
			for(Edge e : c.leftEdges){
				if(leftDest<0){leftDest=e.destination;}
				else if(leftDest!=e.destination){
					if(verbose){System.err.println("Different left destination: "+leftDest+" vs "+e.destination);}
					return false;
				}
				if(e.origin==e.destination){
					if(verbose){System.err.println("Left midnode loop for "+c.name());}
					return false;
				}
			}
			for(Edge e : c.rightEdges){
				if(rightDest<0){rightDest=e.destination;}
				else if(rightDest!=e.destination){
					if(verbose){System.err.println("Different right destination: "+rightDest+" vs "+e.destination);}
					return false;
				}
				if(e.origin==e.destination){
					if(verbose){System.err.println("Right midnode loop for "+c.name());}
					return false;
				}
			}
			ArrayList<Edge> incoming=destMap.get(c.id);
			if(incoming==null){
				if(verbose){System.err.println("No incoming edges.");}
				assert(false);
				return false;
			}
			for(Edge e : incoming){
				if(e.origin!=center.id && e.origin!=dest.id){
					if(verbose){System.err.println("Midnode incoming loop for "+c.name());}
					return false;
				}
			}
		}
		if(leftDest>=0 && leftDest!=center.id){return false;}//workaround for actual assertion failure
		assert(leftDest<0 || leftDest==center.id) : 
			leftDest+", "+center.id; //TODO: This triggered once nondeterministially; determine why
		if(rightDest>=0 && rightDest!=dest.id){return false;}//workaround for potential assertion failure
		assert(rightDest<0 || rightDest==dest.id);
		
		if(verbose){System.err.println("Mid nodes concur.");}
		return leftDest>=0 && rightDest>=0;
	}
	
	/**
	 * Retrieves all intermediate contigs reachable through outbound edges.
	 * Optionally flips contigs to ensure consistent orientation.
	 * Returns null if any duplicate nodes or used nodes are encountered.
	 *
	 * @param outbound The edges leading to intermediate nodes
	 * @param flipAsNeeded Whether to flip nodes for consistent orientation
	 * @return List of intermediate contigs, or null if invalid structure
	 */
	private ArrayList<Contig> fetchMidNodes(ArrayList<Edge> outbound, boolean flipAsNeeded,
			ArrayList<Contig> flippedMids){
		ArrayList<Contig> midNodes=new ArrayList<Contig>(outbound.size());
		for(Edge e : outbound){
			Contig mid=allContigs.get(e.destination);
			if(midNodes.contains(mid)){return null;} //It's possible for there to be 2 edges to the same node. 
//			assert(!mid.used());
			if(mid.used()){return null;}
			if(!mid.used()){
				midNodes.add(mid);
				if(flipAsNeeded && e.destRight()){
					mid.flip(destMap.get(mid.id));
					flippedMids.add(mid);
				}
			}
		}
		return midNodes;
	}
	
	/**
	 * Finds the mutual destination that all paths through intermediate nodes
	 * converge to. Validates that all intermediate nodes have outbound edges
	 * leading to the same destination with consistent orientation.
	 *
	 * @param edges The edges leading to intermediate nodes
	 * @return The mutual destination contig ID, or negative if none found
	 */
	private int	findMutualDest(ArrayList<Edge> edges){
		if(verbose){System.err.println("findMutualDest("+edges+")");}
		lastMutualDest=-2;
		lastMutualDestOrientation=-1;
		for(Edge e : edges){
			if(verbose){System.err.println("\nConsidering inbound edge "+e);}
			Contig mid=allContigs.get(e.destination);
			if(verbose){System.err.println("Considering mid node "+mid.name());}
			if(mid==center){
				if(verbose){System.err.println("Mid node is center.");}
				return -1;
			}
			ArrayList<Edge> outbound=(e.destRight() ? mid.leftEdges : mid.rightEdges);
			if(verbose){System.err.println("e.destRight()="+e.destRight()+", using "+(outbound==mid.leftEdges ? "mid.leftEdges" : "mid.rightEdges"));}
			if(outbound!=null){
				for(Edge o : outbound){
					if(verbose){System.err.println("Considering mid node edge "+o);}
					if(lastMutualDest<0){
						lastMutualDest=o.destination;
						lastMutualDestOrientation=(o.orientation&2);
						if(verbose){System.err.println("Mutual dest is now "+lastMutualDest+", orientation "+lastMutualDestOrientation);}
					}else if(lastMutualDest!=o.destination){
						if(verbose){System.err.println("Mismatched mutual dest: "+lastMutualDest+" versus "+o.destination);}
						return -1;
					}else if(lastMutualDestOrientation!=(o.orientation&2)){
						if(verbose){System.err.println("Mismatched mutual orientation: "+lastMutualDestOrientation+" versus "+(o.orientation&2));}
						return -1;
					}
				}
			}
		}
		return lastMutualDest;//Can be -2 if there is no destination
	}
	
//	boolean testRightBubble(Contig c, ArrayList<Edge> inbound_0){
//		ArrayList<Edge> outbound=c.rightEdges;
//		int mutualDest=findMutualDest(outbound);
//		if(mutualDest<0){return false;}
//		
//		IntList midContigs=new IntList(4);
//		ArrayList<Edge> inbound=new ArrayList<Edge>(inbound_0.size());
//		for(Edge e : inbound_0){
//			if(e.destRight()){
//				inbound.add(e);
//				midContigs.add(e.destination);
//			}
//		}
//		for(Edge e : outbound){
//			midContigs.add(e.destination);
//		}
//		
//		midContigs.sort();
//		for(int i=0; i<midContigs.size(); i+=2){
//			if(midContigs.size()<i+2 || midContigs.get(i)!=midContigs.get(i+1)){
//				assert(false) : midContigs; //for testing; should be removed
//				return false;
//			}
//		}
//		midContigs.shrinkToUnique();
//		//At this point we have established that all edges are bidirectional.
//		
//		int mutualSource=findMutualSource(inbound);
//		if(mutualSource<0){return false;}
//		int mutualLength=findMutualLength(outbound, 3);
//		if(mutualLength>3*kbig){return false;}
//		if(alternativeSource(outbound)){return false;}
//		
//		if(mutualSource==c.id && mutualDest>c.)
//	}
//	
//	boolean alternativeSource(ArrayList<Edge> edges){
//		int source=-2;
//		for(Edge e : edges){
//			Contig mid=allContigs.get(e.destination);
//			ArrayList<Edge> inbound=(e.destRight() ? mid.rightEdges : mid.leftEdges);
//			if(inbound==null || inbound.size()>1){return false;}
//			if(inbound!=null){
//				for(Edge i : inbound){
//					if(source<0){source=i.destination;}
//					else if(source!=i.destination){return true;}
//				}
//			}
//		}
//		return false;
//	}
//	
//	int	findMutualDest(ArrayList<Edge> edges){
//		int dest=-2;
//		int orientation=-1;
//		for(Edge e : edges){
//			Contig mid=allContigs.get(e.destination);
//			ArrayList<Edge> outbound=(e.destRight() ? mid.leftEdges : mid.rightEdges);
////			ArrayList<Edge> inbound=(e.destRight() ? mid.rightEdges : mid.leftEdges);
////			if(inbound!=null && inbound.size()>2){return -1;}
//			if(outbound!=null){
//				for(Edge o : outbound){
//					if(dest<0){
//						dest=o.destination;
//						orientation=(o.orientation&2);
//					}else if(dest!=o.destination){return -1;}
//					else if(orientation!=(o.orientation&2)){return -1;}
//				}
//			}
//		}
//		return dest;//Can be -2 if there is no destination
//	}
//	
//	int	findMutualSource(ArrayList<Edge> edges){
//		int source=-2;
//		for(Edge e : edges){
//			Contig mid=allContigs.get(e.origin);
//			ArrayList<Edge> outbound=(e.sourceRight() ? mid.leftEdges : mid.rightEdges);
//			if(outbound!=null){
//				for(Edge o : outbound){
//					if(source<0){source=o.origin;}
//					else if(source!=o.origin){return -1;}
//				}
//			}
//		}
//		return source;//Can be -2 if there is no source
//	}
//	
//	int	findMutualLength(ArrayList<Edge> edges, int leeway){
//		int min=Integer.MAX_VALUE;
//		int max=-1;
//		for(Edge e : edges){
//			min=Tools.min(min, e.length);
//			max=Tools.max(max, e.length);
//		}
//		return (max-min>leeway ? -1 : max);
//	}
	
	/** Complete list of all contigs in the assembly */
	final ArrayList<Contig> allContigs;
	/** Mapping from contig IDs to their incoming edge lists */
	final HashMap<Integer, ArrayList<Edge>> destMap;
	/** K-mer size used in the assembly */
	final int kbig;
	/** Minimum contig length for bubble resolution (2*kbig-1) */
	final int minLen;
	/** Tadpole's configured kmer-count error classifier; null preserves all indirect alternatives. */
	final ErrorClassifier errorClassifier;
	/** Reusable buffer for sequence concatenation operations */
	final ByteBuilder bb=new ByteBuilder();
	
	/** Current center contig being processed for expansion */
	Contig center=null;
	/** Current destination contig in bubble resolution */
	Contig dest=null;
	
	/** ID of the last found mutual destination contig */
	int lastMutualDest=-1;
	/** Orientation of the last found mutual destination */
	int lastMutualDestOrientation=-1;
	
	/** Counter for successful expansion operations */
	int expansions=0;
	/** Counter for contigs merged or absorbed during operations */
	int contigsAbsorbed=0;
	/** Counter for isolated biological bubbles unzipped into two terminal paths */
	int trueBubblesUnzipped=0;
	/** Counter for branches removed during debranching operations */
	long branchesRemoved=0;

	/** Whether to print detailed debugging information */
	static boolean verbose=false;
	/** Whether direct path merging is enabled */
	static boolean popDirect=true;
	/** Whether depth-gated indirect error-bubble cleanup is enabled */
	static boolean popIndirect=false;
	/** Whether isolated two-arm biological bubbles are unzipped into two terminal paths */
	static boolean unzipBubbles=false;
	/** Whether direct merges join contigs built at a longer kmer length. */
	static boolean crossKMerge=false;
	/** Maximum bridge depth divided by the greater flank coverage; nonpositive disables. */
	static float crossKMaxDepthRatio=3;
	/** Whether debranching of dead ends is enabled */
	static boolean debranch=false;
	/** Whether assertion-enabled runs perform graph consistency checks */
	static boolean validateGraph=false;

	/** Minimal callback used to share Tadpole's tuned error-count heuristic without duplicating constants. */
	interface ErrorClassifier {
		boolean isError(int high, int low);
	}
}
