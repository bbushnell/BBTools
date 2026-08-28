package assemble;

import java.util.ArrayList;
import java.util.HashMap;

import dna.AminoAcid;

/** Finds unique reciprocal exact overlaps between selected contig ends. */
class CrossKTipOverlapper {

	CrossKTipOverlapper(ArrayList<Contig> contigs_, int minOverlap_, int maxOverlap_){
		this(contigs_, minOverlap_, maxOverlap_, false);
	}

	CrossKTipOverlapper(ArrayList<Contig> contigs_, int minOverlap_, int maxOverlap_, boolean graphKEnds_){
		if(minOverlap_<1 || maxOverlap_<minOverlap_){
			throw new IllegalArgumentException("Invalid cross-k overlap range: "+minOverlap_+"-"+maxOverlap_);
		}
		contigs=contigs_;
		minOverlap=minOverlap_;
		maxOverlap=maxOverlap_;
		graphKEnds=graphKEnds_;
	}

	/** Adds reciprocal overlap edges and returns the number of acyclic pairs added. */
	int addEdges(){
		final ArrayList<Tip> tips=makeTips();
		if(tips.size()<2){printSummary(tips.size(), 0, 0, 0); return 0;}

		final HashMap<Long, ArrayList<Tip>> terminalMap=new HashMap<Long, ArrayList<Tip>>(tips.size()*2);
		final long outgoingPower=power(HASH_MULT, minOverlap-1);
		for(Tip tip : tips){
			final boolean reverse=!tip.right;
			final long hash=hash(tip.contig, reverse, tip.contig.length()-minOverlap, minOverlap);
			ArrayList<Tip> list=terminalMap.get(hash);
			if(list==null){list=new ArrayList<Tip>(1); terminalMap.put(hash, list);}
			list.add(tip);
		}

		for(Tip dest : tips){
			final Contig c=dest.contig;
			final boolean reverse=dest.right;
			final int maxStart=Math.min(maxOverlap-minOverlap, c.length()-minOverlap-1);
			if(maxStart<0){continue;}
			long hash=hash(c, reverse, 0, minOverlap);
			for(int start=0; start<=maxStart; start++){
				final ArrayList<Tip> sources=terminalMap.get(hash);
				if(sources!=null){
					final int overlap=minOverlap+start;
					for(Tip source : sources){consider(source, dest, overlap);}
				}
				if(start<maxStart){
					final int oldBase=baseAt(c, reverse, start)&0xFF;
					final int newBase=baseAt(c, reverse, start+minOverlap)&0xFF;
					hash=(hash-(oldBase+1)*outgoingPower)*HASH_MULT+(newBase+1);
				}
			}
		}

		final ArrayList<Pair> reciprocal=new ArrayList<Pair>();
		int ambiguous=0;
		for(Tip tip : tips){
			if(tip.ambiguous){ambiguous++;}
			final Tip dest=tip.best;
			if(dest!=null && !tip.ambiguous && !dest.ambiguous && dest.best==tip
					&& dest.bestOverlap==tip.bestOverlap && tip.index<dest.index){
				reciprocal.add(new Pair(tip, dest, tip.bestOverlap));
			}
		}

		final boolean[] cyclic=cyclicComponents(reciprocal);
		int added=0, cycleRejected=0;
		for(int i=0; i<reciprocal.size(); i++){
			final Pair pair=reciprocal.get(i);
			if(cyclic[i]){cycleRejected++; continue;}
			addEdges(pair);
			added++;
		}
		printSummary(tips.size(), ambiguous, reciprocal.size(), cycleRejected);
		return added;
	}

	private ArrayList<Tip> makeTips(){
		final ArrayList<Tip> tips=new ArrayList<Tip>();
		for(int i=0; i<contigs.size(); i++){
			final Contig c=contigs.get(i);
			if(c.id!=i){throw new RuntimeException("Cross-k overlap contig index mismatch: "+c.id+" != "+i);}
			if(c.length()<minOverlap+1){continue;}
			if(graphKEnds ? c.leftCode==Tadpole.KEEP_GOING : c.leftBridgeEndpoint){
				tips.add(new Tip(c, false, tips.size()));
			}
			if(graphKEnds ? c.rightCode==Tadpole.KEEP_GOING : c.rightBridgeEndpoint){
				tips.add(new Tip(c, true, tips.size()));
			}
		}
		return tips;
	}

	private void consider(Tip source, Tip dest, int overlap){
		if(source.contig==dest.contig || overlap>=source.contig.length() || overlap>=dest.contig.length()){return;}
		if(!matches(source, dest, overlap)){return;}
		exactCandidates++;
		if(overlap>source.bestOverlap){
			source.best=dest;
			source.bestOverlap=overlap;
			source.ambiguous=false;
		}else if(overlap==source.bestOverlap && source.best!=dest){
			source.ambiguous=true;
		}
	}

	private static boolean matches(Tip source, Tip dest, int overlap){
		final boolean sourceReverse=!source.right;
		final boolean destReverse=dest.right;
		final int sourceStart=source.contig.length()-overlap;
		for(int i=0; i<overlap; i++){
			if(baseAt(source.contig, sourceReverse, sourceStart+i)!=baseAt(dest.contig, destReverse, i)){
				return false;
			}
		}
		return true;
	}

	private static void addEdges(Pair pair){
		final Tip a=pair.a, b=pair.b;
		final int depth=Math.max(1, (int)Math.min(a.contig.coverage, b.contig.coverage));
		final int orientationAB=(a.right ? 1 : 0)|(b.right ? 2 : 0);
		final int orientationBA=(b.right ? 1 : 0)|(a.right ? 2 : 0);
		final Edge ab=new Edge(a.contig.id, b.contig.id, 0, orientationAB, depth, null, pair.overlap);
		final Edge ba=new Edge(b.contig.id, a.contig.id, 0, orientationBA, depth, null, pair.overlap);
		if(a.right){a.contig.addRightEdge(ab);}else{a.contig.addLeftEdge(ab);}
		if(b.right){b.contig.addRightEdge(ba);}else{b.contig.addLeftEdge(ba);}
	}

	/** Marks every pair belonging to a cyclic contig-level overlap component. */
	private boolean[] cyclicComponents(ArrayList<Pair> pairs){
		final boolean[] answer=new boolean[pairs.size()];
		if(pairs.isEmpty()){return answer;}
		final int[] parent=new int[contigs.size()];
		final boolean[] present=new boolean[contigs.size()];
		for(int i=0; i<parent.length; i++){parent[i]=i;}
		for(Pair pair : pairs){
			final int a=pair.a.contig.id, b=pair.b.contig.id;
			present[a]=present[b]=true;
			union(parent, a, b);
		}
		final int[] vertices=new int[parent.length], edges=new int[parent.length];
		for(int i=0; i<present.length; i++){if(present[i]){vertices[find(parent, i)]++;}}
		for(Pair pair : pairs){edges[find(parent, pair.a.contig.id)]++;}
		for(int i=0; i<pairs.size(); i++){
			final int root=find(parent, pairs.get(i).a.contig.id);
			answer[i]=(edges[root]>=vertices[root]);
		}
		return answer;
	}

	private static int find(int[] parent, int x){
		while(parent[x]!=x){parent[x]=parent[parent[x]]; x=parent[x];}
		return x;
	}

	private static void union(int[] parent, int a, int b){
		a=find(parent, a); b=find(parent, b);
		if(a!=b){parent[b]=a;}
	}

	private void printSummary(int tips, int ambiguous, int reciprocal, int cycleRejected){
		System.err.println((graphKEnds ? "Graph-k" : "Cross-k")+" tip overlaps: endpoints="+tips+
				", exactCandidates="+exactCandidates+
				", ambiguous="+ambiguous+", reciprocal="+reciprocal+
				", cycleRejected="+cycleRejected+", added="+(reciprocal-cycleRejected)+".");
	}

	private static long hash(Contig c, boolean reverse, int start, int length){
		long hash=0;
		for(int i=0; i<length; i++){hash=hash*HASH_MULT+((baseAt(c, reverse, start+i)&0xFF)+1);}
		return hash;
	}

	private static long power(long base, int exponent){
		long value=1;
		for(int i=0; i<exponent; i++){value*=base;}
		return value;
	}

	private static byte baseAt(Contig c, boolean reverse, int pos){
		return reverse ? AminoAcid.baseToComplementExtended[c.bases[c.length()-1-pos]] : c.bases[pos];
	}

	static class Tip {
		Tip(Contig contig_, boolean right_, int index_){contig=contig_; right=right_; index=index_;}
		final Contig contig;
		final boolean right;
		final int index;
		Tip best;
		int bestOverlap=0;
		boolean ambiguous=false;
	}

	static class Pair {
		Pair(Tip a_, Tip b_, int overlap_){a=a_; b=b_; overlap=overlap_;}
		final Tip a, b;
		final int overlap;
	}

	private final ArrayList<Contig> contigs;
	private final int minOverlap, maxOverlap;
	private final boolean graphKEnds;
	private long exactCandidates=0;
	private static final long HASH_MULT=0x9E3779B185EBCA87L;
}
