package assemble;

import java.util.ArrayList;
import java.util.HashMap;

import shared.Tools;
import structures.ByteBuilder;

public class BubblePopper {
	
	BubblePopper(ArrayList<Contig> allContigs_, HashMap<Integer, ArrayList<Edge>> destMap_, int kbig_){
		allContigs=allContigs_;
		destMap=destMap_;
		kbig=kbig_;
		minLen=2*kbig-1;
	}
	
	int expand(Contig c) {
		//assert(validate(c));
//		if(true) {return 0;}
		if(verbose){System.err.println("\n\n*expand: "+c.name()+"\n");}
		assert(!c.used());
		center=c;
		int count=0;
		
//		debranch(c);
		
		while(popDirect && expandRightSimple()){count++;}
		while(popIndirect && center.rightForwardBranch() && expandRight()){
			count++;
			while(popDirect && expandRightSimple()){count++;}
		}
		
		if((popDirect && center.leftCode!=Tadpole.LOOP && center.leftCode!=Tadpole.DEAD_END && center.leftEdges!=null)
				|| (popIndirect && center.leftForwardBranch())){
			center.flip(destMap.get(center.id));
			while(popDirect && expandRightSimple()){count++;}
			while(popIndirect && center.rightForwardBranch() && expandRight()) {
				count++;
				while(popDirect && expandRightSimple()){count++;}
			}
//			center.flip(destMap.get(center.id));
		}
		//assert(validate(c));
		return count;
	}
	
	void debranch(Contig c){
		assert(debranch);
		debranchRight(c);
		debranchLeft(c);
	}
	
	private void debranchRight(Contig c){
		if(c.rightEdgeCount()<2 || c.rightCode==Tadpole.LOOP){return;}
		debranch(c, c.rightEdges);
	}
	
	private void debranchLeft(Contig c){
		if(c.leftEdgeCount()<2 || c.leftCode==Tadpole.LOOP){return;}
		debranch(c, c.leftEdges);
	}
	
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
	
	private boolean expandRightSimple(){
		//assert(validate(center));
		ArrayList<Edge> outbound=center.rightEdges;
		if(outbound==null || center.rightCode==Tadpole.LOOP || outbound.size()>1){return false;}
		Edge leftEdge=outbound.get(0);
		assert(leftEdge.destination<allContigs.size()) : "\n"+leftEdge.toString()+", "+center.used()+", "+center.associate()+"\n"+center.name()+"\n"+allContigs.size();
		
//		 : "\ncenter="+center.name2()+"\ndest="+dest.name2()+"\nc="+c.name2()+"\nother="+other.name2()+"\ne="+e;
		
		dest=allContigs.get(leftEdge.destination);
		
		if(dest.used() || dest==center){return false;}
		ArrayList<Edge> outboundRight=(leftEdge.destRight() ? dest.rightEdges : dest.leftEdges);
		int rightCode=leftEdge.destRight() ? dest.rightCode : dest.leftCode;
		
		if(rightCode==Tadpole.LOOP){return false;}
		
//		if(outboundRight==null || outboundRight.size()>1){return false;}
		if(outboundRight==null) {
			//do nothing
		}else if(outboundRight.size()>1){
			return false;
		}else if(outboundRight.get(0).destination!=center.id){
			return false;
		}
		
		int inbound=countInbound(center.id, true);
//		if(inbound==null || inbound.size()>1){return false;}
		if(inbound>1){return false;}
		
		int inboundRight=countInbound(dest.id, leftEdge.destRight());
		if(inboundRight>1){return false;}
		
		if(leftEdge.destRight()){
			dest.flip(destMap.get(dest.id));
		}
		//assert(validate(center));
		//assert(validate(dest));
		
		return merge(center, dest, leftEdge);
	}
	
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
	
	private boolean expandRight(){
		//assert(validate(center));
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
		
		ArrayList<Contig> midNodes=fetchMidNodes(outbound, true);
		if(midNodes==null){return false;}
		//Now we have all intermediate nodes, which are flipped into the correct orientation.
		
		if(!midNodesConcur(midNodes)){
			if(verbose){System.err.println("Mid nodes do not concur.");}
			return false;
		}
		
		//At this point all tests have been run and we are confident that this is a simple bubble.
		if(mutualDestRight){
			if(verbose){System.err.println("Flipping dest node.");}
			dest.flip(destMap.get(dest.id));
		}
		
		//Now the destination is also flipped into the correct orientation.
		final Edge rightMidEdge=mid.getRightEdge(dest.id, 1);
		if(rightMidEdge==null){return false;}

		if(verbose){System.err.println("Popping bubble between "+center.id+" and "+dest.id);}
		return pop(center, dest, mid, leftMidEdge, rightMidEdge, midNodes);
		
	}
	
	private boolean pop(Contig left, Contig right, Contig mid, Edge leftMidEdge, Edge rightMidEdge, ArrayList<Contig> midNodes){
		//assert(validate(left));
		//assert(validate(right));
		//assert(validate(mid));
		for(Contig c : midNodes){
			assert(!c.used() && !c.associate());
			assert(c!=right && c!=left);
			//assert(validate(c));
		}
		
		bb.clear();
		final int originalLeftLength=left.length();
		
		//Append path
		bb.append(left.bases);
		if(leftMidEdge.bases!=null){
			for(int i=0; i<leftMidEdge.bases.length-1; i++){
				bb.append(leftMidEdge.bases[i]);
			}
		}
		for(int i=kbig-1, lim=mid.length()-kbig+1; i<lim; i++){
			bb.append(mid.bases[i]);
		}
		if(rightMidEdge.bases!=null){
			for(int i=0; i<rightMidEdge.bases.length-1; i++){
				bb.append(rightMidEdge.bases[i]);
			}
		}
		bb.append(right.bases);
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
		left.minCov=Tools.min(left.maxCov, right.maxCov);
		left.rightCode=right.rightCode;
		left.rightRatio=right.rightRatio;
		double coverageSum=left.coverage*(originalLeftLength-kbig+1)+right.coverage*(right.length()-kbig+1);
		left.coverage=(float)(coverageSum/(left.length()-kbig+1));
		
		expansions++;
		contigsAbsorbed+=(1+midNodes.size());
		
		if(isLoop(left)){
			left.leftCode=Tadpole.LOOP;
			left.rightCode=Tadpole.LOOP;
			left.removeAllEdges(destMap.get(left.id), allContigs);
		}
		
		if(verbose){
			System.err.println("*Result: "+center.name());
		}

		//assert(validate(left));
		//assert(validate(right));
		//assert(validate(mid));
		for(Contig c : midNodes){
			//assert(validate(c));
		}
		
		return true;
	}
	
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
	
	void removeDeadEdges(Contig c){
		c.leftEdges=removeDeadEdges(c.leftEdges);
		c.rightEdges=removeDeadEdges(c.rightEdges);
	}
	
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
	
	private boolean merge(Contig left, Contig right, Edge leftEdge){

		//assert(validate(left));
		//assert(validate(right));
		
		bb.clear();
		final int originalLeftLength=left.length();
		
		//Append path
		bb.append(left.bases);
		if(leftEdge.bases!=null){
			for(int i=0; i<leftEdge.bases.length-1; i++){
				bb.append(leftEdge.bases[i]);
			}
		}
		bb.append(right.bases);
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
		left.minCov=Tools.min(left.maxCov, right.maxCov);
		left.rightCode=right.rightCode;
		left.rightRatio=right.rightRatio;
		double coverageSum=left.coverage*(originalLeftLength-kbig+1)+right.coverage*(right.length()-kbig+1);
		left.coverage=(float)(coverageSum/(left.length()-kbig+1));
		
		if(isLoop(left)){
			left.leftCode=Tadpole.LOOP;
			left.rightCode=Tadpole.LOOP;
			left.removeAllEdges(destMap.get(left.id), allContigs);
		}
		
		expansions++;
		contigsAbsorbed++;
		
		if(verbose){
			System.err.println("*Result: "+center.name());
		}

//		//assert(validate(left));
//		//assert(validate(right));
		
		return true;
	}
	
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
	
	boolean validate(Contig c){
		if(true){return true;}
		if(c.used() || c.associate()){
			assert(c.leftEdges==null);
			assert(c.rightEdges==null);
		}else{
			ArrayList<Edge> inbound=destMap.get(c.id);
			if(inbound!=null){
				for(Edge e : inbound){
					Contig other=allContigs.get(e.origin);
					if(other.used() || other.associate()){
						//ignore
					}else{
						assert(e.destination==c.id) : "\nc="+c.name2()+"\nother="+other.name2()+"\ne="+e+
							"\nid="+c.id+", origin="+e.origin+", dest="+e.destination+", size="+allContigs.size()+"\n"+inbound;
						assert(e.origin<allContigs.size());
						assert(e.destination<allContigs.size());
					}
				}
			}
			if(c.leftEdges!=null){
				assert(c.leftEdges.size()>0);
				for(Edge e : c.leftEdges){
					assert(e.destination<allContigs.size()) : "\nc="+c.name2()+"\ne="+e+
						"\nid="+c.id+", origin="+e.origin+", size="+allContigs.size();
					Contig other=allContigs.get(e.destination);
					assert(!e.sourceRight());
					assert(e.origin==c.id);
					assert(e.origin<allContigs.size()) : /*"\ncenter="+center.name2()+"\ndest="+dest.name2()+*/"\nc="+c.name2()+"\nother="+other.name2()+"\ne="+e+
						"\nid="+c.id+", origin="+e.origin+", size="+allContigs.size();
					assert(!other.used()) : "\nc="+c.name2()+"\nother="+other.name2()+"\ne="+e+
						"\nid="+c.id+", origin="+e.origin+", size="+allContigs.size();;
					assert(!other.associate());
				}
			}
			if(c.rightEdges!=null){
				assert(c.rightEdges.size()>0);
				for(Edge e : c.rightEdges){
					assert(e.destination<allContigs.size());
					Contig other=allContigs.get(e.destination);
					assert(other!=null) : "\ncenter="+center.name2()+"\ndest="+dest.name2()+"\nc="+c.name2()+"\ne="+e+"\nsize="+allContigs.size();
					assert(e.sourceRight());
					assert(e.origin==c.id);
					assert(e.origin<allContigs.size());
//					assert(!other.used()) : "\ncenter="+center.name2()+"\ndest="+dest.name2()+"\nc="+c.name2()+"\nother="+other.name2()+"\ne="+e;
					assert(!other.used()) : /*"\ncenter="+center.name2()+"\ndest="+dest.name2()+*/"\nc="+c.name2()+"\nother="+other.name2()+
						"\ne="+e+"\nsize="+allContigs.size();
					assert(!other.associate());
				}
			}
		}
		return true;
	}
	
	private Edge findRepresentativeMidEdge(ArrayList<Edge> edges){
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
		if(rightDest>=0 && rightDest!=center.id){return false;}//workaround for potential assertion failure
		assert(rightDest<0 || rightDest==dest.id);
		
		if(verbose){System.err.println("Mid nodes concur.");}
		return leftDest>=0 && rightDest>=0;
	}
	
	private ArrayList<Contig> fetchMidNodes(ArrayList<Edge> outbound, boolean flipAsNeeded){
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
				}
			}
		}
		return midNodes;
	}
	
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
	
	final ArrayList<Contig> allContigs;
	final HashMap<Integer, ArrayList<Edge>> destMap;
	final int kbig;
	final int minLen;
	final ByteBuilder bb=new ByteBuilder();
	
	Contig center=null;
	Contig dest=null;
	
	int lastMutualDest=-1;
	int lastMutualDestOrientation=-1;
	
	int expansions=0;
	int contigsAbsorbed=0;
	long branchesRemoved=0;

	static boolean verbose=false;
	static boolean popDirect=true;
	static boolean popIndirect=true;
	static boolean debranch=false;
}
