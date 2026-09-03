package assemble;

import java.util.ArrayList;
import java.util.Collections;

import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Extracts maximal topology-safe walks from an oriented contig graph.
 *
 * Each oriented contig is one vertex in the contig-adjacency line graph.  Maximal
 * unitigs in that graph are simple omnitigs in the underlying compacted assembly
 * graph: shared sequence may occur in multiple output walks, but no traversal ever
 * chooses a pairing across a branch.
 *
 * @author Noelle
 */
class SimpleOmnitigExtractor {

	SimpleOmnitigExtractor(ArrayList<Contig> contigs_, int k_){
		this(contigs_, k_, 1);
	}

	SimpleOmnitigExtractor(ArrayList<Contig> contigs_, int k_, int minAnchorLen_){
		contigs=contigs_;
		k=k_;
		minAnchorLen=minAnchorLen_;
		covered=new boolean[Math.multiplyExact(2, contigs.size())];
		visitStamp=new int[covered.length];
	}

	/** Returns canonical representatives of every maximal safe walk. */
	ArrayList<Contig> extract(){
		validateInput();
		final ArrayList<Contig> output=new ArrayList<Contig>();
		classifyEdges();

		/* Isolated contigs have no adjacency arc from which to seed a walk. */
		for(Contig c : contigs){
			inputBases+=c.length();
			if(outDegree(c.id<<1)==0 && outDegree((c.id<<1)|1)==0){
				covered[c.id<<1]=covered[(c.id<<1)|1]=true;
				emitSingleton(c, output);
			}
		}

		/* Every outgoing arc of a non-1-in/1-out state seeds one maximal walk. */
		for(int state=0; state<covered.length; state++){
			final int in=inDegree(state), out=outDegree(state);
			if(out<1 || (in==1 && out==1)){continue;}
			final ArrayList<Edge> edges=outEdges(state);
			for(Edge edge : edges){if(usableEdge(state, edge)){emitWalk(state, edge, output);}}
		}

		/* An unseeded component is a pure 1-in/1-out cycle.  Preserve its original
		 * contigs rather than selecting a rotation or repeating a circular walk. */
		for(Contig c : contigs){
			if(!covered[c.id<<1] && !covered[(c.id<<1)|1]){
				preservedCycleContigs++;
				covered[c.id<<1]=covered[(c.id<<1)|1]=true;
				emitSingleton(c, output);
			}
		}

		for(int i=0; i<output.size(); i++){output.get(i).id=i;}
		return output;
	}

	/**
	 * Returns a deterministic copy-aware path cover of the usable graph.
	 * Ordinary contigs occur once.  A clean 2-by-2 X center occurs twice, once
	 * for each edge on one selected boundary; the opposite boundary remains
	 * unresolved so no left-to-right pairing is invented.
	 */
	ArrayList<Contig> extractNonredundant(){
		validateInput();
		classifyEdges();
		final int statesCount=covered.length;
		final int[] degree=new int[statesCount];
		for(int state=0; state<statesCount; state++){degree[state]=outDegree(state);}
		final boolean[] xContig=new boolean[contigs.size()];
		for(Contig c : contigs){xContig[c.id]=(degree[c.id<<1]>1 && degree[(c.id<<1)|1]>1);}
		final ArrayList<XFanPlan> fanPlans=makeXFanPlans(degree);
		final boolean[] reserved=new boolean[contigs.size()];
		final ArrayList<XFanPlan> acceptedFans=new ArrayList<XFanPlan>();
		for(XFanPlan plan : fanPlans){
			final int center=plan.state>>>1;
			final int a=plan.edges[0].destination, b=plan.edges[1].destination;
			if(reserved[center] || reserved[a] || reserved[b]
					|| xContig[a] || xContig[b] || !fanHasOutputAnchor(center, a, b)){continue;}
			reserved[center]=reserved[a]=reserved[b]=true;
			acceptedFans.add(plan);
		}

		final ArrayList<PathCandidate> candidates=new ArrayList<PathCandidate>();
		for(int state=0; state<statesCount; state++){
			final ArrayList<Edge> edges=outEdges(state);
			if(edges==null){continue;}
			for(Edge edge : edges){
				if(!usableEdge(state, edge)){continue;}
				final int destination=destinationState(edge);
				if(xContig[state>>>1] || xContig[destination>>>1]
						|| reserved[state>>>1] || reserved[destination>>>1]){continue;}
				final int reverseState=destination^1;
				if(state>reverseState || (state==reverseState && destination>=(state^1))){continue;}
				final Edge reverseEdge=reciprocalEdge(state, edge);
				candidates.add(new PathCandidate(state, edge, reverseState, reverseEdge,
						degree[state]+degree[reverseState]-2,
						Tools.min(contig(state).length(), contig(destination).length())));
			}
		}
		Collections.sort(candidates);
		pathCandidates=candidates.size();

		final Edge[] selectedOut=new Edge[statesCount];
		final DisjointSet sets=new DisjointSet(contigs.size());
		for(PathCandidate candidate : candidates){
			final int state=candidate.state;
			final int destination=destinationState(candidate.edge);
			final int reverseState=candidate.reverseState;
			if(selectedOut[state]!=null || selectedOut[reverseState]!=null){continue;}
			if(sets.find(state>>>1)==sets.find(destination>>>1)){
				pathCyclesRejected++;
				continue;
			}
			selectedOut[state]=candidate.edge;
			selectedOut[reverseState]=candidate.reverseEdge;
			sets.union(state>>>1, destination>>>1);
			selectedJoins++;
		}
		discardUnanchoredPaths(selectedOut, sets);

		final ArrayList<Contig> output=new ArrayList<Contig>();
		final boolean[] used=new boolean[contigs.size()];
		for(XFanPlan plan : acceptedFans){
			final int center=plan.state>>>1;
			final long outputBefore=outputBases;
			long originalBases=contigs.get(center).length();
			for(Edge edge : plan.edges){
				final int neighbor=edge.destination;
				if(used[neighbor]){throw new IllegalStateException("X-fan reused flank contig "+neighbor);}
				used[neighbor]=true;
				originalBases+=contigs.get(neighbor).length();
				emitFanPath(plan.state, edge, output);
			}
			used[center]=true;
			xCentersDuplicated++;
			xCopiesEmitted+=2;
			xCenterBasesDuplicated+=contigs.get(center).length();
			xNetBases+=outputBases-outputBefore-originalBases;
		}
		for(int state=0; state<statesCount; state++){
			if(selectedOut[state]!=null && selectedOut[state^1]==null && !used[state>>>1]){
				emitSelectedPath(state, selectedOut, used, output);
			}
		}
		for(Contig c : contigs){
			inputBases+=c.length();
			if(!used[c.id]){
				used[c.id]=true;
				emitSingleton(c, output);
			}
		}
		for(int i=0; i<output.size(); i++){output.get(i).id=i;}
		return output;
	}

	/** Both X-fan products must contain sequence that was independently output-worthy. */
	private boolean fanHasOutputAnchor(final int center, final int a, final int b){
		return contigs.get(center).length()>=minAnchorLen
				|| (contigs.get(a).length()>=minAnchorLen && contigs.get(b).length()>=minAnchorLen);
	}

	/** Prevents joins among short graph debris from promoting it across the output threshold. */
	private void discardUnanchoredPaths(final Edge[] selectedOut, final DisjointSet sets){
		final boolean[] anchored=new boolean[contigs.size()];
		for(Contig c : contigs){
			if(c.length()>=minAnchorLen){anchored[sets.find(c.id)]=true;}
		}
		int discardedArcs=0;
		for(int state=0; state<selectedOut.length; state++){
			if(selectedOut[state]==null || anchored[sets.find(state>>>1)]){continue;}
			selectedOut[state]=null;
			discardedArcs++;
		}
		if((discardedArcs&1)!=0){throw new IllegalStateException("Unpaired selected graph arcs: "+discardedArcs);}
		smallOnlyJoinsDiscarded=discardedArcs>>>1;
		selectedJoins-=smallOnlyJoinsDiscarded;
	}

	private ArrayList<XFanPlan> makeXFanPlans(final int[] degree){
		final ArrayList<XFanPlan> plans=new ArrayList<XFanPlan>();
		for(Contig c : contigs){
			final int forward=c.id<<1, reverse=forward|1;
			if(degree[forward]!=2 || degree[reverse]!=2){continue;}
			final XFanPlan a=makeXFanPlan(forward), b=makeXFanPlan(reverse);
			if(a==null || b==null){continue;}
			plans.add(a.compareTo(b)<=0 ? a : b);
		}
		Collections.sort(plans);
		return plans;
	}

	private XFanPlan makeXFanPlan(final int state){
		final Edge[] edges=new Edge[2];
		int size=0;
		final ArrayList<Edge> list=outEdges(state);
		if(list==null){return null;}
		for(Edge edge : list){
			if(usableEdge(state, edge)){
				if(size>=2){return null;}
				edges[size++]=edge;
			}
		}
		if(size!=2 || edges[0].destination==edges[1].destination
				|| edges[0].destination==(state>>>1) || edges[1].destination==(state>>>1)){return null;}
		return new XFanPlan(state, edges, contigs);
	}

	private void emitFanPath(final int state, final Edge edge, final ArrayList<Contig> output){
		states.clear();
		pathEdges.clear();
		validateEdge(state, edge);
		states.add(state);
		states.add(destinationState(edge));
		pathEdges.add(edge);
		final byte[] bases=makeSequence();
		if(compareToReverseComplement(bases)>0){
			reverseComplement(bases);
			reversePathInPlace();
		}
		output.add(makeContig(bases));
		walksEmitted++;
		extendedWalks++;
		outputBases+=bases.length;
	}

	private void classifyEdges(){
		for(int state=0; state<covered.length; state++){
			final ArrayList<Edge> edges=outEdges(state);
			if(edges!=null){
				for(Edge edge : edges){
					if(!joinCompatible(state, edge)){inexactEdges++;}
					else if(reciprocalEdge(state, edge)==null){nonreciprocalEdges++;}
				}
			}
		}
	}

	private void emitSelectedPath(final int start, final Edge[] selectedOut, final boolean[] used,
			final ArrayList<Contig> output){
		states.clear();
		pathEdges.clear();
		int state=start;
		while(true){
			final int id=state>>>1;
			if(used[id]){throw new IllegalStateException("Graph-cover path reused contig "+id);}
			used[id]=true;
			states.add(state);
			final Edge edge=selectedOut[state];
			if(edge==null){break;}
			validateEdge(state, edge);
			pathEdges.add(edge);
			state=destinationState(edge);
		}
		final byte[] bases=makeSequence();
		if(compareToReverseComplement(bases)>0){
			reverseComplement(bases);
			reversePathInPlace();
		}
		output.add(makeContig(bases));
		walksEmitted++;
		if(states.size>1){extendedWalks++;}
		outputBases+=bases.length;
	}

	/** Reorients path metadata after canonicalizing a nonredundant product. */
	private void reversePathInPlace(){
		for(int i=0, j=states.size-1; i<j; i++, j--){
			final int a=states.get(i)^1, b=states.get(j)^1;
			states.set(i, b);
			states.set(j, a);
		}
		if((states.size&1)==1){
			final int mid=states.size>>>1;
			states.set(mid, states.get(mid)^1);
		}
	}

	private void emitWalk(final int start, final Edge first, final ArrayList<Contig> output){
		states.clear();
		pathEdges.clear();
		states.add(start);
		int current=start;
		Edge edge=first;
		final int stamp=nextStamp();
		visitStamp[start]=stamp;

		while(true){
			validateEdge(current, edge);
			final int next=destinationState(edge);
			if(visitStamp[next]==stamp){
				truncatedCycleWalks++;
				break;
			}
			pathEdges.add(edge);
			current=next;
			states.add(current);
			covered[current]=covered[current^1]=true;
			visitStamp[current]=stamp;
			if(inDegree(current)!=1 || outDegree(current)!=1){break;}
			edge=firstUsableOutEdge(current);
		}
		covered[start]=covered[start^1]=true;

		final byte[] bases=makeSequence();
		final int canonical=compareToReverseComplement(bases);
		if(canonical>0 || (canonical==0 && comparePathToReverseTwin()>0)){return;}
		output.add(makeContig(bases));
		walksEmitted++;
		if(states.size>1){extendedWalks++;}
		outputBases+=bases.length;
	}

	private byte[] makeSequence(){
		bb.clear();
		appendContig(bb, states.get(0), 0);
		for(int i=0; i<pathEdges.size(); i++){
			final Edge edge=pathEdges.get(i);
			appendEdge(bb, edge, states.get(i), states.get(i+1));
			final int overlap=(edge.overlap>0 ? edge.overlap : k);
			appendContig(bb, states.get(i+1), overlap);
		}
		return bb.toBytes();
	}

	/**
	 * Appends an edge in walk orientation.  A dense edge longer than k is a
	 * compound field: its leading length-k bytes are novel traversal bases,
	 * while its final k bytes reproduce the destination's terminal kmer in the
	 * destination contig's stored strand.  They cannot be reverse-complemented
	 * as one ordinary sequence.  Orient only the novel bytes, then spell the
	 * represented terminal portion directly from the oriented destination.
	 */
	private void appendEdge(final ByteBuilder builder, final Edge edge,
			final int sourceState, final int destinationState){
		if(edge.bases==null){return;}
		final int represented=Tools.min(k, edge.length);
		final int novel=edge.length-represented;
		for(int i=0; i<novel; i++){
			final byte b=edge.bases[i];
			builder.append(reverse(sourceState) ? AminoAcid.baseToComplementExtended[b] : b);
		}
		final Contig target=contig(destinationState);
		for(int i=k-represented; i<k; i++){
			builder.append(orientedBase(target, destinationState, i));
		}
	}

	private Contig makeContig(final byte[] bases){
		final Contig out=new Contig(bases, -1);
		double coverageSum=0;
		long coverageLength=0;
		int minCov=Integer.MAX_VALUE, maxCov=0;
		for(int i=0; i<states.size; i++){
			final Contig c=contig(states.get(i));
			final int span=Tools.max(1, c.length()-k+1);
			coverageSum+=c.coverage*span;
			coverageLength+=span;
			minCov=Tools.min(minCov, c.minCov);
			maxCov=Tools.max(maxCov, c.maxCov);
		}
		for(Edge edge : pathEdges){
			final int span=Tools.max(0, edge.length-1);
			coverageSum+=edge.depth*span;
			coverageLength+=span;
		}
		out.coverage=(float)(coverageSum/Tools.max(1L, coverageLength));
		out.minCov=(minCov==Integer.MAX_VALUE ? 0 : minCov);
		out.maxCov=maxCov;

		final int first=states.get(0), last=states.lastElement();
		final Contig left=contig(first), right=contig(last);
		final boolean firstReverse=reverse(first), lastReverse=reverse(last);
		out.leftCode=firstReverse ? left.rightCode : left.leftCode;
		out.leftRatio=firstReverse ? left.rightRatio : left.leftRatio;
		out.leftBridgeEndpoint=firstReverse ? left.rightBridgeEndpoint : left.leftBridgeEndpoint;
		out.rightCode=lastReverse ? right.leftCode : right.rightCode;
		out.rightRatio=lastReverse ? right.leftRatio : right.rightRatio;
		out.rightBridgeEndpoint=lastReverse ? right.leftBridgeEndpoint : right.rightBridgeEndpoint;
		return out;
	}

	private void emitSingleton(final Contig source, final ArrayList<Contig> output){
		final byte[] bases=source.bases.clone();
		final Contig out;
		if(compareToReverseComplement(bases)<=0){
			out=copySingleton(source, bases, false);
		}else{
			reverseComplement(bases);
			out=copySingleton(source, bases, true);
		}
		output.add(out);
		outputBases+=out.length();
	}

	private static Contig copySingleton(final Contig source, final byte[] bases, final boolean reverse){
		final Contig out=new Contig(bases, -1);
		out.coverage=source.coverage;
		out.minCov=source.minCov;
		out.maxCov=source.maxCov;
		out.leftCode=reverse ? source.rightCode : source.leftCode;
		out.leftRatio=reverse ? source.rightRatio : source.leftRatio;
		out.leftBridgeEndpoint=reverse ? source.rightBridgeEndpoint : source.leftBridgeEndpoint;
		out.rightCode=reverse ? source.leftCode : source.rightCode;
		out.rightRatio=reverse ? source.leftRatio : source.rightRatio;
		out.rightBridgeEndpoint=reverse ? source.leftBridgeEndpoint : source.rightBridgeEndpoint;
		return out;
	}

	private void appendContig(final ByteBuilder builder, final int state, final int skip){
		final Contig c=contig(state);
		if(skip<0 || skip>c.length()){
			throw new IllegalStateException("Invalid omnitig overlap "+skip+" for contig "+c.id+
					" of length "+c.length());
		}
		if(!reverse(state)){
			builder.append(c.bases, skip, c.length()-skip);
		}else{
			for(int i=skip; i<c.length(); i++){
				builder.append(AminoAcid.baseToComplementExtended[c.bases[c.length()-1-i]]);
			}
		}
	}

	private void validateInput(){
		if(k<1){throw new IllegalArgumentException("Simple-omnitig k must be positive: "+k);}
		for(int i=0; i<contigs.size(); i++){
			final Contig c=contigs.get(i);
			if(c==null || c.id!=i || c.used() || c.associate() || c.length()<k){
				throw new IllegalStateException("Invalid omnitig graph contig at index "+i+": "+c);
			}
			validateEdges(c, c.leftEdges, false);
			validateEdges(c, c.rightEdges, true);
		}
	}

	private void validateEdges(final Contig source, final ArrayList<Edge> edges, final boolean right){
		if(edges==null){return;}
		for(Edge edge : edges){
			if(edge==null || edge.origin!=source.id || edge.sourceRight()!=right
					|| edge.destination<0 || edge.destination>=contigs.size()
					|| edge.length<0 || (edge.bases==null ? edge.length!=0 : edge.bases.length!=edge.length)
					|| edge.overlap<0){
				throw new IllegalStateException("Invalid omnitig graph edge from contig "+source.id+": "+edge);
			}
		}
	}

	private void validateEdge(final int state, final Edge edge){
		if(edge.origin!=contig(state).id || edge.sourceRight()==reverse(state)){
			throw new IllegalStateException("Edge is inconsistent with oriented omnitig state "+state+": "+edge);
		}
		final int destination=destinationState(edge);
		final int overlap=(edge.overlap>0 ? edge.overlap : k);
		final Contig source=contig(state), target=contig(destination);
		if(overlap>source.length() || overlap>target.length() || (edge.length<1 && edge.overlap<1)){
			throw new IllegalStateException("Invalid omnitig edge span from state "+state+": "+edge);
		}
		if(!usableEdge(state, edge)){
			throw new IllegalStateException("Omnitig attempted an inexact graph edge from state "+state+": "+edge);
		}
	}

	/** Returns true only when an edge spells an exact linear join in this orientation. */
	private boolean joinCompatible(final int state, final Edge edge){
		final int destination=destinationState(edge);
		final Contig source=contig(state), target=contig(destination);
		final int overlap=(edge.overlap>0 ? edge.overlap : k);
		if(overlap>source.length() || overlap>target.length() || (edge.length<1 && edge.overlap<1)){
			return false;
		}
		if(edge.overlap>0){
			for(int i=0; i<overlap; i++){
				if(orientedBase(source, state, source.length()-overlap+i)!=orientedBase(target, destination, i)){
					return false;
				}
			}
		}else{
			final int shared=Tools.max(0, k-edge.length);
			for(int i=0; i<shared; i++){
				if(orientedBase(source, state, source.length()-shared+i)!=orientedBase(target, destination, i)){
					return false;
				}
			}
		}
		return true;
	}

	/** A topology-safe arc must have an exact reverse traversal of the same physical path. */
	private boolean usableEdge(final int state, final Edge edge){
		return joinCompatible(state, edge) && reciprocalEdge(state, edge)!=null;
	}

	private Edge reciprocalEdge(final int state, final Edge edge){
		final int reverseState=destinationState(edge)^1;
		final ArrayList<Edge> reverseEdges=outEdges(reverseState);
		if(reverseEdges==null){return null;}
		for(Edge candidate : reverseEdges){
			if(candidate.destination!=(state>>>1) || destinationState(candidate)!=(state^1)
					|| candidate.length!=edge.length || candidate.overlap!=edge.overlap
					|| !joinCompatible(reverseState, candidate)){continue;}
			final int novel=Tools.max(0, edge.length-k);
			int i=0;
			while(i<novel && orientedNovelBase(edge, state, i)==
					AminoAcid.baseToComplementExtended[orientedNovelBase(candidate, reverseState, novel-1-i)]){i++;}
			if(i==novel){return candidate;}
		}
		return null;
	}

	private static byte orientedNovelBase(final Edge edge, final int state, final int pos){
		final byte b=edge.bases[pos];
		return reverse(state) ? AminoAcid.baseToComplementExtended[b] : b;
	}

	private static byte orientedBase(final Contig c, final int state, final int pos){
		return reverse(state) ? AminoAcid.baseToComplementExtended[c.bases[c.length()-1-pos]] : c.bases[pos];
	}

	private int inDegree(final int state){
		return outDegree(state^1);
	}

	private int outDegree(final int state){
		final ArrayList<Edge> edges=outEdges(state);
		if(edges==null){return 0;}
		int count=0;
		for(Edge edge : edges){if(usableEdge(state, edge)){count++;}}
		return count;
	}

	private Edge firstUsableOutEdge(final int state){
		final ArrayList<Edge> edges=outEdges(state);
		if(edges!=null){for(Edge edge : edges){if(usableEdge(state, edge)){return edge;}}}
		throw new IllegalStateException("No usable edge from 1-out omnitig state "+state);
	}

	private ArrayList<Edge> outEdges(final int state){
		final Contig c=contig(state);
		return reverse(state) ? c.leftEdges : c.rightEdges;
	}

	private Contig contig(final int state){return contigs.get(state>>>1);}
	private static boolean reverse(final int state){return (state&1)!=0;}
	private static int destinationState(final Edge edge){return (edge.destination<<1)|(edge.destRight() ? 1 : 0);}

	private int nextStamp(){
		if(++stamp==0){
			for(int i=0; i<visitStamp.length; i++){visitStamp[i]=0;}
			stamp=1;
		}
		return stamp;
	}

	private int comparePathToReverseTwin(){
		for(int i=0, j=states.size-1; i<states.size; i++, j--){
			final int a=states.get(i), b=states.get(j)^1;
			if(a!=b){return a<b ? -1 : 1;}
		}
		return 0;
	}

	private static int compareToReverseComplement(final byte[] bases){
		for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
			final int a=bases[i]&0xFF;
			final int b=AminoAcid.baseToComplementExtended[bases[j]]&0xFF;
			if(a!=b){return a<b ? -1 : 1;}
		}
		return 0;
	}

	private static void reverseComplement(final byte[] bases){
		for(int i=0, j=bases.length-1; i<=j; i++, j--){
			final byte a=AminoAcid.baseToComplementExtended[bases[j]];
			final byte b=AminoAcid.baseToComplementExtended[bases[i]];
			bases[i]=a;
			bases[j]=b;
		}
	}

	private static class PathCandidate implements Comparable<PathCandidate> {

		PathCandidate(final int state_, final Edge edge_, final int reverseState_, final Edge reverseEdge_,
				final int ambiguity_, final int minLength_){
			state=state_;
			edge=edge_;
			reverseState=reverseState_;
			reverseEdge=reverseEdge_;
			ambiguity=ambiguity_;
			minLength=minLength_;
			support=Tools.min(edge.depth, reverseEdge.depth);
		}

		@Override
		public int compareTo(final PathCandidate b){
			if(ambiguity!=b.ambiguity){return ambiguity<b.ambiguity ? -1 : 1;}
			if(minLength!=b.minLength){return minLength>b.minLength ? -1 : 1;}
			if(support!=b.support){return support>b.support ? -1 : 1;}
			if(state!=b.state){return state<b.state ? -1 : 1;}
			final int destination=destinationState(edge), bDestination=destinationState(b.edge);
			return destination==bDestination ? 0 : destination<bDestination ? -1 : 1;
		}

		final int state, reverseState, ambiguity, minLength, support;
		final Edge edge, reverseEdge;
	}

	private static class XFanPlan implements Comparable<XFanPlan> {

		XFanPlan(final int state_, final Edge[] edges_, final ArrayList<Contig> contigs){
			state=state_;
			edges=edges_;
			final Contig center=contigs.get(state>>>1);
			centerLength=center.length();
			neighborLength=(long)contigs.get(edges[0].destination).length()+contigs.get(edges[1].destination).length();
			support=(long)edges[0].depth+edges[1].depth;
		}

		@Override
		public int compareTo(final XFanPlan b){
			if(neighborLength!=b.neighborLength){return neighborLength>b.neighborLength ? -1 : 1;}
			if(support!=b.support){return support>b.support ? -1 : 1;}
			if(centerLength!=b.centerLength){return centerLength>b.centerLength ? -1 : 1;}
			return state==b.state ? 0 : state<b.state ? -1 : 1;
		}

		final int state, centerLength;
		final long neighborLength, support;
		final Edge[] edges;
	}

	private static class DisjointSet {

		DisjointSet(final int size){
			parent=new int[size];
			rank=new byte[size];
			for(int i=0; i<size; i++){parent[i]=i;}
		}

		int find(int x){
			int root=x;
			while(parent[root]!=root){root=parent[root];}
			while(parent[x]!=x){
				final int next=parent[x];
				parent[x]=root;
				x=next;
			}
			return root;
		}

		void union(final int a, final int b){
			final int x=find(a), y=find(b);
			if(x==y){return;}
			if(rank[x]<rank[y]){parent[x]=y;}
			else{
				parent[y]=x;
				if(rank[x]==rank[y]){rank[x]++;}
			}
		}

		final int[] parent;
		final byte[] rank;
	}

	long inputBases(){return inputBases;}
	long outputBases(){return outputBases;}
	int walksEmitted(){return walksEmitted;}
	int extendedWalks(){return extendedWalks;}
	int preservedCycleContigs(){return preservedCycleContigs;}
	int inexactEdges(){return inexactEdges;}
	int nonreciprocalEdges(){return nonreciprocalEdges;}
	int truncatedCycleWalks(){return truncatedCycleWalks;}
	int pathCandidates(){return pathCandidates;}
	int selectedJoins(){return selectedJoins;}
	int pathCyclesRejected(){return pathCyclesRejected;}
	int smallOnlyJoinsDiscarded(){return smallOnlyJoinsDiscarded;}
	int xCentersDuplicated(){return xCentersDuplicated;}
	int xCopiesEmitted(){return xCopiesEmitted;}
	long xCenterBasesDuplicated(){return xCenterBasesDuplicated;}
	long xNetBases(){return xNetBases;}

	private final ArrayList<Contig> contigs;
	private final int k, minAnchorLen;
	private final boolean[] covered;
	private final int[] visitStamp;
	private final IntList states=new IntList();
	private final ArrayList<Edge> pathEdges=new ArrayList<Edge>();
	private final ByteBuilder bb=new ByteBuilder();
	private int stamp=0;
	private int walksEmitted=0, extendedWalks=0, preservedCycleContigs=0;
	private int inexactEdges=0, nonreciprocalEdges=0;
	private int truncatedCycleWalks=0;
	private int pathCandidates=0, selectedJoins=0, pathCyclesRejected=0, smallOnlyJoinsDiscarded=0;
	private int xCentersDuplicated=0, xCopiesEmitted=0;
	private long inputBases=0, outputBases=0, xCenterBasesDuplicated=0, xNetBases=0;
}
