package assemble;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;

/**
 * Classifies contigs on independent length, depth, and long-anchor topology axes.
 * @author Noelle
 */
final class ContigGraphClassifier {

	ContigGraphClassifier(final ArrayList<Contig> contigs_, final int minContig_, final int k_,
			final float lowDepthCutoff_, final float mediumDepthCutoff_, final float highDepthCutoff_){
		contigs=contigs_;
		minContig=minContig_;
		k=k_;
		lowDepthCutoff=lowDepthCutoff_;
		mediumDepthCutoff=mediumDepthCutoff_;
		highDepthCutoff=highDepthCutoff_;
		if(k<1){throw new IllegalArgumentException("Graph-classification k must be positive: "+k+".");}
		if(lowDepthCutoff<0 || mediumDepthCutoff<lowDepthCutoff || highDepthCutoff<mediumDepthCutoff){
			throw new IllegalArgumentException("Invalid graph depth thresholds: "+lowDepthCutoff+", "+
					mediumDepthCutoff+", "+highDepthCutoff+".");
		}
	}

	Result classify(){
		final int states=Math.multiplyExact(contigs.size(), 2);
		final int[] anchor1=KillSwitch.allocInt1D(states), anchor2=KillSwitch.allocInt1D(states);
		final int[] distance1=KillSwitch.allocInt1D(states), distance2=KillSwitch.allocInt1D(states);
		final int queueLength=Math.multiplyExact(states, 2);
		final int[] queueState=KillSwitch.allocInt1D(queueLength);
		final int[] queueAnchor=KillSwitch.allocInt1D(queueLength);
		final int[] queueDistance=KillSwitch.allocInt1D(queueLength);
		Arrays.fill(anchor1, -1);
		Arrays.fill(anchor2, -1);
		for(int i=0; i<contigs.size(); i++){
			final Contig c=contigs.get(i);
			if(c.id!=i){throw new IllegalStateException("Nonconsecutive contig ID "+c.id+" at index "+i+".");}
			c.graphClassHop=0;
			if(c.used() || c.associate()){
				c.graphClass=Contig.GRAPH_UNCLASSIFIED;
				c.graphLengthClass=Contig.LENGTH_UNCLASSIFIED;
				c.graphDepthClass=Contig.DEPTH_UNCLASSIFIED;
			}else{
				c.graphClass=Contig.GRAPH_UNANCHORED;
				c.graphLengthClass=(c.length()>=minContig ? Contig.LENGTH_LONG : Contig.LENGTH_SHORT);
				c.graphDepthClass=depthClass(c.coverage);
			}
		}

		int head=0, tail=0;
		for(Contig c : contigs){
			if(c.graphLengthClass!=Contig.LENGTH_LONG){continue;}
			tail=seed(c, false, c.leftEdges, c.id, anchor1, anchor2, distance1, distance2,
					queueState, queueAnchor, queueDistance, tail);
			tail=seed(c, true, c.rightEdges, c.id, anchor1, anchor2, distance1, distance2,
					queueState, queueAnchor, queueDistance, tail);
		}
		while(head<tail){
			final int state=queueState[head], anchor=queueAnchor[head], distance=queueDistance[head++];
			final Contig c=contigs.get(state>>>1);
			final boolean enteredRight=(state&1)!=0;
			final boolean exitRight=!enteredRight;
			final ArrayList<Edge> exits=(exitRight ? c.rightEdges : c.leftEdges);
			if(exits==null){continue;}
			for(Edge e : exits){
				tail=visit(c, exitRight, e, anchor, distance+1, anchor1, anchor2, distance1, distance2,
						queueState, queueAnchor, queueDistance, tail);
			}
		}

		final Result result=new Result(lowDepthCutoff, mediumDepthCutoff, highDepthCutoff);
		for(Contig c : contigs){
			if(c.graphClass==Contig.GRAPH_UNCLASSIFIED){continue;}
			classifyTopology(c, anchor1, anchor2, distance1, distance2);
			result.add(c);
		}
		return result;
	}

	private int depthClass(final float coverage){
		if(coverage<=lowDepthCutoff){return Contig.DEPTH_LOW;}
		if(coverage<mediumDepthCutoff){return Contig.DEPTH_SEMILOW;}
		if(coverage<=highDepthCutoff){return Contig.DEPTH_MEDIUM;}
		return Contig.DEPTH_HIGH;
	}

	private void classifyTopology(final Contig c, final int[] anchor1, final int[] anchor2,
			final int[] distance1, final int[] distance2){
		if(selfLoop(c)){
			c.graphClass=Contig.GRAPH_SELF_LOOP;
			return;
		}
		int leftState=c.id<<1, rightState=leftState|1;
		int leftCount=anchorCount(leftState, anchor1, anchor2);
		int rightCount=anchorCount(rightState, anchor1, anchor2);
		if(leftCount>rightCount){
			final int tempState=leftState; leftState=rightState; rightState=tempState;
			final int tempCount=leftCount; leftCount=rightCount; rightCount=tempCount;
		}
		if(rightCount==0){c.graphClass=Contig.GRAPH_UNANCHORED; return;}
		if(leftCount==0){
			c.graphClass=(rightCount==1 ? Contig.GRAPH_TERMINAL : Contig.GRAPH_BRANCHED_TERMINAL);
			c.graphClassHop=nearestDistance(rightState, anchor1, anchor2, distance1, distance2);
			return;
		}
		if(leftCount==1 && rightCount==1){
			if(anchor1[leftState]==anchor1[rightState]){
				c.graphClass=Contig.GRAPH_LOOPBACK;
			}else{c.graphClass=Contig.GRAPH_CONNECTED;}
			c.graphClassHop=Tools.max(distance1[leftState], distance1[rightState]);
			return;
		}
		c.graphClass=(leftCount==1 ? Contig.GRAPH_BRANCHED_CONNECTED : Contig.GRAPH_MULTI_CONNECTED);
		c.graphClassHop=connectedHop(leftState, rightState, anchor1, anchor2, distance1, distance2);
		if(c.graphClassHop<1){
			throw new IllegalStateException("Multi-anchor topology lacked distinct anchors on contig "+c.id+".");
		}
	}

	private int seed(final Contig source, final boolean sourceRight, final ArrayList<Edge> edges,
			final int anchor, final int[] anchor1, final int[] anchor2,
			final int[] distance1, final int[] distance2, final int[] queueState,
			final int[] queueAnchor, final int[] queueDistance, int tail){
		if(edges==null){return tail;}
		for(Edge e : edges){
			tail=visit(source, sourceRight, e, anchor, 1, anchor1, anchor2, distance1, distance2,
					queueState, queueAnchor, queueDistance, tail);
		}
		return tail;
	}

	private int visit(final Contig source, final boolean sourceRight, final Edge edge,
			final int anchor, final int nextDistance, final int[] anchor1, final int[] anchor2,
			final int[] distance1, final int[] distance2, final int[] queueState,
			final int[] queueAnchor, final int[] queueDistance, int tail){
		if(edge.origin!=source.id || edge.sourceRight()!=sourceRight){
			throw new IllegalStateException("Misplaced graph edge "+edge+" on contig "+source.id+".");
		}
		if(edge.destination<0 || edge.destination>=contigs.size()){
			throw new IllegalStateException("Graph edge has invalid destination: "+edge+".");
		}
		final Contig dest=contigs.get(edge.destination);
		if(dest.used() || dest.associate() || dest.id==anchor || !usableEdge(source, sourceRight, edge, dest)){return tail;}
		final int state=(dest.id<<1)|(edge.destRight() ? 1 : 0);
		if(anchor1[state]==anchor || anchor2[state]==anchor){return tail;}
		if(anchor1[state]<0){anchor1[state]=anchor; distance1[state]=nextDistance;}
		else if(anchor2[state]<0){anchor2[state]=anchor; distance2[state]=nextDistance;}
		else{return tail;}
		if(dest.graphLengthClass==Contig.LENGTH_LONG){return tail;}
		queueState[tail]=state;
		queueAnchor[tail]=anchor;
		queueDistance[tail++]=nextDistance;
		return tail;
	}

	private static int anchorCount(final int state, final int[] anchor1, final int[] anchor2){
		return anchor1[state]<0 ? 0 : anchor2[state]<0 ? 1 : 2;
	}

	/** Returns the shortest worst-side distance supported by two distinct long anchors. */
	private static int connectedHop(final int left, final int right, final int[] anchor1,
			final int[] anchor2, final int[] distance1, final int[] distance2){
		int best=Integer.MAX_VALUE;
		if(anchor1[left]>=0 && anchor1[right]>=0 && anchor1[left]!=anchor1[right]){
			best=Tools.min(best, Tools.max(distance1[left], distance1[right]));
		}
		if(anchor1[left]>=0 && anchor2[right]>=0 && anchor1[left]!=anchor2[right]){best=Tools.min(best, Tools.max(distance1[left], distance2[right]));}
		if(anchor2[left]>=0 && anchor1[right]>=0 && anchor2[left]!=anchor1[right]){best=Tools.min(best, Tools.max(distance2[left], distance1[right]));}
		if(anchor2[left]>=0 && anchor2[right]>=0 && anchor2[left]!=anchor2[right]){best=Tools.min(best, Tools.max(distance2[left], distance2[right]));}
		return best==Integer.MAX_VALUE ? -1 : best;
	}

	private static int nearestDistance(final int state, final int[] anchor1, final int[] anchor2,
			final int[] distance1, final int[] distance2){
		if(anchor1[state]<0){return -1;}
		return anchor2[state]<0 ? distance1[state] : Tools.min(distance1[state], distance2[state]);
	}

	private boolean selfLoop(final Contig c){
		return c.leftCode==Tadpole.LOOP || c.rightCode==Tadpole.LOOP ||
				selfEdge(c, false, c.leftEdges) || selfEdge(c, true, c.rightEdges);
	}

	private boolean selfEdge(final Contig source, final boolean sourceRight, final ArrayList<Edge> edges){
		if(edges==null){return false;}
		for(Edge e : edges){
			if(e.destination==source.id && e.origin==source.id && e.sourceRight()==sourceRight &&
					usableEdge(source, sourceRight, e, source)){return true;}
		}
		return false;
	}

	/** A classifier arc must spell an exact path and have the exact reverse traversal. */
	private boolean usableEdge(final Contig source, final boolean sourceRight, final Edge edge, final Contig dest){
		final int state=(source.id<<1)|(sourceRight ? 0 : 1);
		return joinCompatible(state, edge) && reciprocalEdge(state, edge, dest)!=null;
	}

	private Edge reciprocalEdge(final int state, final Edge edge, final Contig dest){
		final int reverseState=destinationState(edge)^1;
		final ArrayList<Edge> reverse=outEdges(reverseState);
		if(reverse==null){return null;}
		for(Edge back : reverse){
			if(back.origin!=dest.id || back.destination!=(state>>>1) || destinationState(back)!=(state^1)
					|| back.length!=edge.length || back.overlap!=edge.overlap
					|| back.sourceTrim!=edge.destTrim || back.destTrim!=edge.sourceTrim
					|| !joinCompatible(reverseState, back)){continue;}
			final int novel=Tools.max(0, edge.length-k);
			int i=0;
			while(i<novel && orientedNovelBase(edge, state, i)==
					AminoAcid.baseToComplementExtended[orientedNovelBase(back, reverseState, novel-1-i)]){i++;}
			if(i==novel){return back;}
		}
		return null;
	}

	/** Returns true only when an edge spells an exact join in its oriented traversal. */
	private boolean joinCompatible(final int state, final Edge edge){
		if(edge==null || edge.origin!=(state>>>1) || edge.sourceRight()==reverse(state)
				|| edge.destination<0 || edge.destination>=contigs.size() || edge.length<0
				|| edge.sourceTrim<0 || edge.destTrim<0 || edge.overlap<0
				|| (edge.bases==null ? edge.length!=0 : edge.bases.length!=edge.length)){return false;}
		final int destination=destinationState(edge);
		final Contig source=contig(state), target=contig(destination);
		if(edge.overlap>0){
			if(edge.sourceTrim+edge.overlap>=source.length() || edge.destTrim+edge.overlap>=target.length()){return false;}
			final int sourceStart=source.length()-edge.sourceTrim-edge.overlap;
			for(int i=0; i<edge.overlap; i++){
				if(orientedBase(source, state, sourceStart+i)!=orientedBase(target, destination, edge.destTrim+i)){return false;}
			}
			return true;
		}
		if(edge.sourceTrim!=0 || edge.destTrim!=0 || edge.length<1 || k>source.length() || k>target.length()){return false;}
		/* Bytes represented inside the destination terminal kmer are redundant and
		 * retain legacy side-dependent encoding.  Validate the exact shared contig
		 * sequence here; reciprocalEdge validates the only emitted payload, the
		 * novel prefix beyond k. */
		final int shared=Tools.max(0, k-edge.length);
		for(int i=0; i<shared; i++){
			if(orientedBase(source, state, source.length()-shared+i)!=orientedBase(target, destination, i)){return false;}
		}
		return true;
	}

	private ArrayList<Edge> outEdges(final int state){
		final Contig c=contig(state);
		return reverse(state) ? c.leftEdges : c.rightEdges;
	}

	private Contig contig(final int state){return contigs.get(state>>>1);}
	private static boolean reverse(final int state){return (state&1)!=0;}
	private static int destinationState(final Edge edge){return (edge.destination<<1)|(edge.destRight() ? 1 : 0);}

	private static byte orientedNovelBase(final Edge edge, final int state, final int pos){
		final byte b=edge.bases[pos];
		return reverse(state) ? AminoAcid.baseToComplementExtended[b] : b;
	}

	private static byte orientedBase(final Contig c, final int state, final int pos){
		return reverse(state) ? AminoAcid.baseToComplementExtended[c.bases[c.length()-1-pos]] : c.bases[pos];
	}

	static final class Result {
		@SuppressWarnings("unchecked")
		Result(final float lowDepthCutoff_, final float mediumDepthCutoff_, final float highDepthCutoff_){
			lowDepthCutoff=lowDepthCutoff_;
			mediumDepthCutoff=mediumDepthCutoff_;
			highDepthCutoff=highDepthCutoff_;
			topology=(ArrayList<Contig>[])new ArrayList[8];
			for(int i=0; i<topology.length; i++){topology[i]=new ArrayList<Contig>();}
		}

		void add(final Contig c){
			topology[c.graphClass].add(c);
			topologyBases[c.graphClass]+=c.length();
			lengthCount[c.graphLengthClass]++;
			lengthBases[c.graphLengthClass]+=c.length();
			depthCount[c.graphDepthClass]++;
			depthBases[c.graphDepthClass]+=c.length();
			combinationCount[c.graphLengthClass][c.graphDepthClass][c.graphClass]++;
			combinationBases[c.graphLengthClass][c.graphDepthClass][c.graphClass]+=c.length();
			maxHop[c.graphClass]=Tools.max(maxHop[c.graphClass], c.graphClassHop);
		}

		final ArrayList<Contig>[] topology;
		final long[] topologyBases=new long[8];
		final long[] lengthCount=new long[2], lengthBases=new long[2];
		final long[] depthCount=new long[4], depthBases=new long[4];
		final long[][][] combinationCount=new long[2][4][8], combinationBases=new long[2][4][8];
		final int[] maxHop=new int[8];
		final float lowDepthCutoff, mediumDepthCutoff, highDepthCutoff;
	}

	private final ArrayList<Contig> contigs;
	private final int minContig, k;
	private final float lowDepthCutoff, mediumDepthCutoff, highDepthCutoff;
}
