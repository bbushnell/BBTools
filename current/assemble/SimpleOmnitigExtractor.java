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
 * One instance supports one extraction pass and is not thread-safe. The graph is
 * borrowed and must remain unchanged during extraction; outputs are new contigs
 * with new base arrays. Traversal state and statistics are retained after the pass.
 *
 * @author Noelle
 */
public class SimpleOmnitigExtractor{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a single-use extractor with no nonredundant path-cover length threshold.
	 *
	 * @param contigs_ borrowed graph whose contig ids match list order
	 * @param k_ graph k-mer length and default overlap for graph edges without explicit overlap
	 */
	public SimpleOmnitigExtractor(ArrayList<Contig> contigs_, int k_){
		this(contigs_, k_, 1);
	}

	/**
	 * Creates a single-use extractor for one graph.
	 * Extraction validates positive K, nonnull unused/unassociated contigs with
	 * ids matching list indices, lengths at least K, and structurally valid edges.
	 * The constructor allocates traversal arrays but does not validate those entries.
	 *
	 * @param contigs_ borrowed graph whose contig ids match list order
	 * @param k_ graph k-mer length and default overlap for graph edges without explicit overlap
	 * @param minAnchorLen_ anchor length used only by extractNonredundant; singletons are not filtered
	 */
	public SimpleOmnitigExtractor(ArrayList<Contig> contigs_, int k_, int minAnchorLen_){
		contigs=contigs_;
		k=k_;
		minAnchorLen=minAnchorLen_;
		covered=new boolean[Math.multiplyExact(2, contigs.size())];
		visitStamp=new int[covered.length];
	}

	/*--------------------------------------------------------------*/
	/*----------------        Maximal Walks         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns canonical representatives of every maximal safe walk.
	 *
	 * This single-use method records visited states and cumulative statistics in
	 * this extractor; construct a fresh instance for another extraction.
	 */
	public ArrayList<Contig> extract(){
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

	/*--------------------------------------------------------------*/
	/*----------------          Path Cover          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns a deterministic copy-aware path cover of the usable graph.
	 * Ordinary contigs occur once.  A clean 2-by-2 X center occurs twice, once
	 * for each edge on one selected boundary; the opposite boundary remains
	 * unresolved so no left-to-right pairing is invented.  This single-use
	 * method consumes the extractor state; construct a fresh instance for another
	 * extraction.
	 */
	public ArrayList<Contig> extractNonredundant(){
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

	/*--------------------------------------------------------------*/
	/*----------------        Graph Planning        ----------------*/
	/*--------------------------------------------------------------*/

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

	/** Selects and sorts one preferred side for each eligible two-in/two-out X center. */
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

	/** Returns exactly two usable exits with distinct non-self destinations, or null. */
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

	/*--------------------------------------------------------------*/
	/*----------------        Product Emission      ----------------*/
	/*--------------------------------------------------------------*/

	/** Emits one accepted X-fan path and updates walk/base counters; the caller counts duplicated centers. */
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

	/** Counts every unusable exactness or reciprocity failure in the oriented graph. */
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

	/** Emits one selected nonredundant path from its uncovered terminal state. */
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

	/** Emits one maximal topology-safe walk from a non-1-in/1-out seed arc. */
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

	/** Reconstructs bases for the currently buffered oriented states and edges. */
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
	 * Appends an edge in walk orientation. Payloads are stored in the source
	 * strand, so left exits reverse-complement them just as Edge.flipSource does.
	 * The destination supplies its represented terminal bases; only the novel
	 * prefix needs to be read from the oriented payload.
	 */
	private void appendEdge(final ByteBuilder builder, final Edge edge,
			final int sourceState, final int destinationState){
		if(edge.bases==null){return;}
		final int represented=Tools.min(k, edge.length);
		final int novel=edge.length-represented;
		for(int i=0; i<novel; i++){
			builder.append(edge.traversalBase(i));
		}
		final Contig target=contig(destinationState);
		for(int i=k-represented; i<k; i++){
			builder.append(orientedBase(target, destinationState, i));
		}
	}

	/** Creates an output contig with coverage and terminal metadata from the buffered path. */
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
		out.leftFusionEndpoint=firstReverse ? left.rightFusionEndpoint : left.leftFusionEndpoint;
		out.rightCode=lastReverse ? right.leftCode : right.rightCode;
		out.rightRatio=lastReverse ? right.leftRatio : right.rightRatio;
		out.rightBridgeEndpoint=lastReverse ? right.leftBridgeEndpoint : right.rightBridgeEndpoint;
		out.rightFusionEndpoint=lastReverse ? right.leftFusionEndpoint : right.rightFusionEndpoint;
		return out;
	}

	/** Emits one isolated or preserved-cycle contig in canonical orientation. */
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

	/** Uses the supplied base array and copies coverage plus strand-adjusted terminal metadata. */
	private static Contig copySingleton(final Contig source, final byte[] bases, final boolean reverse){
		final Contig out=new Contig(bases, -1);
		out.coverage=source.coverage;
		out.minCov=source.minCov;
		out.maxCov=source.maxCov;
		out.leftCode=reverse ? source.rightCode : source.leftCode;
		out.leftRatio=reverse ? source.rightRatio : source.leftRatio;
		out.leftBridgeEndpoint=reverse ? source.rightBridgeEndpoint : source.leftBridgeEndpoint;
		out.leftFusionEndpoint=reverse ? source.rightFusionEndpoint : source.leftFusionEndpoint;
		out.rightCode=reverse ? source.leftCode : source.rightCode;
		out.rightRatio=reverse ? source.leftRatio : source.rightRatio;
		out.rightBridgeEndpoint=reverse ? source.leftBridgeEndpoint : source.rightBridgeEndpoint;
		out.rightFusionEndpoint=reverse ? source.leftFusionEndpoint : source.rightFusionEndpoint;
		return out;
	}

	/** Appends an oriented contig, skipping already represented leading overlap. */
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

	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/

	/** Verifies the borrowed graph can be traversed by id-indexed oriented states. */
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

	/** Verifies stored edges leave the expected source side and point inside the graph. */
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

	/** Verifies one selected edge is legal and sequence-compatible in this orientation. */
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

	/*--------------------------------------------------------------*/
	/*----------------        Oriented Graph        ----------------*/
	/*--------------------------------------------------------------*/

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

	/** Returns the exact reverse traversal of an oriented edge, or null when absent. */
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

	/** Returns one novel payload base in traversal orientation. */
	private static byte orientedNovelBase(final Edge edge, final int state, final int pos){
		return edge.traversalBase(pos);
	}

	/** Returns one contig base in the orientation encoded by a state. */
	private static byte orientedBase(final Contig c, final int state, final int pos){
		return reverse(state) ? AminoAcid.baseToComplementExtended[c.bases[c.length()-1-pos]] : c.bases[pos];
	}

	/** Returns usable incoming degree for an oriented state. */
	private int inDegree(final int state){
		return outDegree(state^1);
	}

	/** Returns usable outgoing degree for an oriented state. */
	private int outDegree(final int state){
		final ArrayList<Edge> edges=outEdges(state);
		if(edges==null){return 0;}
		int count=0;
		for(Edge edge : edges){if(usableEdge(state, edge)){count++;}}
		return count;
	}

	/** Returns the first usable exit; callers establish that the state has exactly one. */
	private Edge firstUsableOutEdge(final int state){
		final ArrayList<Edge> edges=outEdges(state);
		if(edges!=null){for(Edge edge : edges){if(usableEdge(state, edge)){return edge;}}}
		throw new IllegalStateException("No usable edge from 1-out omnitig state "+state);
	}

	/** Returns the edge list leaving an oriented state. */
	private ArrayList<Edge> outEdges(final int state){
		final Contig c=contig(state);
		return reverse(state) ? c.leftEdges : c.rightEdges;
	}

	/** Returns the contig for an oriented state. */
	private Contig contig(final int state){return contigs.get(state>>>1);}
	/** Returns true when a state traverses its contig in reverse-complement orientation. */
	private static boolean reverse(final int state){return (state&1)!=0;}
	/** Returns the oriented destination state reached by an edge. */
	private static int destinationState(final Edge edge){return (edge.destination<<1)|(edge.destRight() ? 1 : 0);}

	/*--------------------------------------------------------------*/
	/*----------------           Utilities          ----------------*/
	/*--------------------------------------------------------------*/

	/** Advances the traversal stamp, clearing old stamps when it wraps to zero. */
	private int nextStamp(){
		if(++stamp==0){
			for(int i=0; i<visitStamp.length; i++){visitStamp[i]=0;}
			stamp=1;
		}
		return stamp;
	}

	/** Orders a buffered walk against the same physical walk seeded from the reverse end. */
	private int comparePathToReverseTwin(){
		for(int i=0, j=states.size-1; i<states.size; i++, j--){
			final int a=states.get(i), b=states.get(j)^1;
			if(a!=b){return a<b ? -1 : 1;}
		}
		return 0;
	}

	/** Compares a base array lexicographically to its reverse complement. */
	private static int compareToReverseComplement(final byte[] bases){
		for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
			final int a=bases[i]&0xFF;
			final int b=AminoAcid.baseToComplementExtended[bases[j]]&0xFF;
			if(a!=b){return a<b ? -1 : 1;}
		}
		return 0;
	}

	/** Reverse-complements a base array in place. */
	private static void reverseComplement(final byte[] bases){
		for(int i=0, j=bases.length-1; i<=j; i++, j--){
			final byte a=AminoAcid.baseToComplementExtended[bases[j]];
			final byte b=AminoAcid.baseToComplementExtended[bases[i]];
			bases[i]=a;
			bases[j]=b;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Pair of reciprocal arcs considered as one undirected path-cover join. */
	private static class PathCandidate implements Comparable<PathCandidate>{

		/** Borrows reciprocal edges and snapshots the keys used for deterministic join priority. */
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

		/** Orders joins by ambiguity, then length, support and topology for deterministic selection. */
		@Override
		public int compareTo(final PathCandidate b){
			if(ambiguity!=b.ambiguity){return ambiguity<b.ambiguity ? -1 : 1;}
			if(minLength!=b.minLength){return minLength>b.minLength ? -1 : 1;}
			if(support!=b.support){return support>b.support ? -1 : 1;}
			if(state!=b.state){return state<b.state ? -1 : 1;}
			final int destination=destinationState(edge), bDestination=destinationState(b.edge);
			return destination==bDestination ? 0 : destination<bDestination ? -1 : 1;
		}

		/** Exit states, excess endpoint degree, shorter contig length and minimum reciprocal depth. */
		final int state, reverseState, ambiguity, minLength, support;
		/** Borrowed forward and reverse traversals of the same physical join. */
		final Edge edge, reverseEdge;
	}

	/** One candidate side of an X center, containing its two outgoing usable edges. */
	private static class XFanPlan implements Comparable<XFanPlan>{

		/** Borrows the two-edge array and snapshots neighbor span, depth and center length. */
		XFanPlan(final int state_, final Edge[] edges_, final ArrayList<Contig> contigs){
			state=state_;
			edges=edges_;
			final Contig center=contigs.get(state>>>1);
			centerLength=center.length();
			neighborLength=(long)contigs.get(edges[0].destination).length()+contigs.get(edges[1].destination).length();
			support=(long)edges[0].depth+edges[1].depth;
		}

		/** Orders X centers by flank span, support, center length and state. */
		@Override
		public int compareTo(final XFanPlan b){
			if(neighborLength!=b.neighborLength){return neighborLength>b.neighborLength ? -1 : 1;}
			if(support!=b.support){return support>b.support ? -1 : 1;}
			if(centerLength!=b.centerLength){return centerLength>b.centerLength ? -1 : 1;}
			return state==b.state ? 0 : state<b.state ? -1 : 1;
		}

		/** Oriented center state and its contig length. */
		final int state, centerLength;
		/** Sum of the two neighbor lengths and outgoing edge depths. */
		final long neighborLength, support;
		/** Borrowed exits from this candidate side, retaining graph-list order. */
		final Edge[] edges;
	}

	/** Union-find structure used to keep the nonredundant path cover acyclic. */
	private static class DisjointSet{

		/** Starts one component per contig id. */
		DisjointSet(final int size){
			parent=new int[size];
			rank=new byte[size];
			for(int i=0; i<size; i++){parent[i]=i;}
		}

		/** Returns the canonical component for a contig id. */
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

		/** Merges the components containing two contig ids. */
		void union(final int a, final int b){
			final int x=find(a), y=find(b);
			if(x==y){return;}
			if(rank[x]<rank[y]){parent[x]=y;}
			else{
				parent[y]=x;
				if(rank[x]==rank[y]){rank[x]++;}
			}
		}

		/** Parent links; roots point to themselves. */
		final int[] parent;
		/** Union-by-rank bounds, independent of component sizes. */
		final byte[] rank;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns input bases visited by this extraction. */
	long inputBases(){return inputBases;}
	/** Returns bases emitted by this extraction. */
	long outputBases(){return outputBases;}
	/** Returns the number of maximal walks or selected paths emitted. */
	int walksEmitted(){return walksEmitted;}
	/** Returns the number of emitted walks spanning more than one contig. */
	public int extendedWalks(){return extendedWalks;}
	/** Returns pure-cycle contigs preserved as singletons. */
	public int preservedCycleContigs(){return preservedCycleContigs;}
	/** Returns edges rejected because their sequence does not spell an exact join. */
	public int inexactEdges(){return inexactEdges;}
	/** Returns exact edges rejected because the reverse traversal is absent. */
	public int nonreciprocalEdges(){return nonreciprocalEdges;}
	/** Returns branch-seeded walks truncated before repeating a state. */
	public int truncatedCycleWalks(){return truncatedCycleWalks;}
	/** Returns reciprocal candidate joins considered by the path cover. */
	int pathCandidates(){return pathCandidates;}
	/** Returns ordinary path-cover joins retained after anchor filtering; excludes X-fan joins. */
	public int selectedJoins(){return selectedJoins;}
	/** Returns candidate joins rejected because they would close a path-cover cycle. */
	int pathCyclesRejected(){return pathCyclesRejected;}
	/** Returns selected joins discarded for lacking a long anchor. */
	public int smallOnlyJoinsDiscarded(){return smallOnlyJoinsDiscarded;}
	/** Returns clean X centers represented by two selected paths. */
	public int xCentersDuplicated(){return xCentersDuplicated;}
	/** Returns X-center copies emitted into selected paths. */
	public int xCopiesEmitted(){return xCopiesEmitted;}
	/** Returns raw X-center bases duplicated by selected paths. */
	public long xCenterBasesDuplicated(){return xCenterBasesDuplicated;}
	/** Returns net bases added by duplicating accepted X centers. */
	public long xNetBases(){return xNetBases;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed id-indexed graph; extraction reads but does not intentionally mutate it. */
	private final ArrayList<Contig> contigs;
	/** Graph k-mer length and minimum anchoring length for nonredundant products. */
	private final int k, minAnchorLen;
	/** States visited by maximal-walk processing or retained as singletons; reverse twins are marked together. */
	private final boolean[] covered;
	/** Traversal stamps for cycle detection without per-walk allocation. */
	private final int[] visitStamp;
	/** Buffered oriented states for the product currently being emitted. */
	private final IntList states=new IntList();
	/** Buffered edges connecting {@link #states} for the product currently being emitted. */
	private final ArrayList<Edge> pathEdges=new ArrayList<Edge>();
	/** Reusable base builder for output sequence reconstruction. */
	private final ByteBuilder bb=new ByteBuilder();
	/** Last assigned visit stamp; initially zero, which is skipped on wrap. */
	private int stamp=0;
	/** Maximal-walk output counters. */
	private int walksEmitted=0, extendedWalks=0, preservedCycleContigs=0;
	/** Rejected-edge counters. */
	private int inexactEdges=0, nonreciprocalEdges=0;
	/** Branch-seeded cycle truncation counter. */
	private int truncatedCycleWalks=0;
	/** Nonredundant path-cover candidate and selection counters. */
	private int pathCandidates=0, selectedJoins=0, pathCyclesRejected=0, smallOnlyJoinsDiscarded=0;
	/** Copy-aware X-center counters. */
	private int xCentersDuplicated=0, xCopiesEmitted=0;
	/** Base totals for input, output and accepted X-center duplication. */
	private long inputBases=0, outputBases=0, xCenterBasesDuplicated=0, xNetBases=0;
}
