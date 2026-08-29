package assemble;

import java.util.ArrayList;

import dna.AminoAcid;
import map.LongLongHashMap2;

/** Finds unique reciprocal exact overlaps between selected contig ends. */
class CrossKTipOverlapper {

	CrossKTipOverlapper(ArrayList<Contig> contigs_, int minOverlap_, int maxOverlap_){
		this(contigs_, minOverlap_, maxOverlap_, false, 0);
	}

	CrossKTipOverlapper(ArrayList<Contig> contigs_, int minOverlap_, int maxOverlap_, boolean graphKEnds_){
		this(contigs_, minOverlap_, maxOverlap_, graphKEnds_, 0);
	}

	CrossKTipOverlapper(ArrayList<Contig> contigs_, int minOverlap_, int maxOverlap_,
			boolean graphKEnds_, int minContig_){
		if(minOverlap_<1 || maxOverlap_<minOverlap_){
			throw new IllegalArgumentException("Invalid cross-k overlap range: "+minOverlap_+"-"+maxOverlap_);
		}
		contigs=contigs_;
		minOverlap=minOverlap_;
		maxOverlap=maxOverlap_;
		graphKEnds=graphKEnds_;
		minContig=minContig_;
	}

	/** Adds reciprocal overlap edges and returns the number of acyclic pairs added. */
	int addEdges(){
		final ArrayList<Tip> tips=makeTips();
		if(tips.size()<2){printSummary(tips.size(), 0, 0, 0); return 0;}

		final int maxTrim=maxOverlap-minOverlap+1;
		final long expected=(long)tips.size()*(maxTrim+1);
		final int initialSize=(int)Math.min(Integer.MAX_VALUE/4, Math.max(256L, expected));
		final LongLongHashMap2 terminalMap=new LongLongHashMap2(initialSize);
		final LongLongHashMap2 duplicateMap=new LongLongHashMap2(Math.max(256, tips.size()));
		final long outgoingPower=power(HASH_MULT, minOverlap-1);
		for(Tip tip : tips){
			final boolean reverse=!tip.right;
			final int trimLimit=Math.min(maxTrim, tip.contig.length()-minOverlap-1);
			long hash=hash(tip.contig, reverse, tip.contig.length()-trimLimit-minOverlap, minOverlap);
			for(int trim=trimLimit; trim>=0; trim--){
				addAnchor(terminalMap, duplicateMap, hash, (((long)tip.index)<<32)|(trim&0xFFFFFFFFL));
				if(trim>0){
					final int start=tip.contig.length()-trim-minOverlap;
					final int oldBase=baseAt(tip.contig, reverse, start)&0xFF;
					final int newBase=baseAt(tip.contig, reverse, start+minOverlap)&0xFF;
					hash=(hash-(oldBase+1)*outgoingPower)*HASH_MULT+(newBase+1);
				}
			}
		}

		for(Tip dest : tips){
			final Contig c=dest.contig;
			final boolean reverse=dest.right;
			final int trimLimit=Math.min(maxTrim, c.length()-minOverlap-1);
			long hash=hash(c, reverse, 0, minOverlap);
			for(int destTrim=0; destTrim<=trimLimit; destTrim++){
				final long code=terminalMap.get(hash);
				if(code>0 && code!=AMBIGUOUS){
					consider(tips, code-1, dest, destTrim);
					final long duplicate=duplicateMap.get(hash);
					if(duplicate>0){consider(tips, duplicate-1, dest, destTrim);}
				}
				if(destTrim<trimLimit){
					final int oldBase=baseAt(c, reverse, destTrim)&0xFF;
					final int newBase=baseAt(c, reverse, destTrim+minOverlap)&0xFF;
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
					&& dest.bestOverlap==tip.bestOverlap
					&& tip.bestSourceTrim==dest.bestDestTrim && tip.bestDestTrim==dest.bestSourceTrim
					&& tip.index<dest.index){
				reciprocal.add(new Pair(tip, dest, tip.bestOverlap, tip.bestSourceTrim, tip.bestDestTrim));
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

	/** Stores one anchor per tip and at most two tips per hash; repeated or third tips are ambiguous. */
	private static void addAnchor(LongLongHashMap2 primary, LongLongHashMap2 duplicate,
			long hash, long code){
		final long encoded=code+1;
		final long old=primary.get(hash);
		if(old<0){primary.set(hash, encoded);}
		else if(old==AMBIGUOUS || old==encoded){return;}
		else if(((old-1)>>>32)==(code>>>32)){primary.set(hash, AMBIGUOUS);}
		else{
			final long old2=duplicate.get(hash);
			if(old2<0){duplicate.set(hash, encoded);}
			else if(old2!=encoded){primary.set(hash, AMBIGUOUS);}
		}
	}

	private ArrayList<Tip> makeTips(){
		final ArrayList<Tip> tips=new ArrayList<Tip>();
		for(int i=0; i<contigs.size(); i++){
			final Contig c=contigs.get(i);
			if(c.id!=i){throw new RuntimeException("Cross-k overlap contig index mismatch: "+c.id+" != "+i);}
			if(c.length()<Math.max(minOverlap+1, minContig)){continue;}
			if(graphKEnds ? c.leftCode==Tadpole.KEEP_GOING : c.leftBridgeEndpoint){
				tips.add(new Tip(c, false, tips.size()));
			}
			if(graphKEnds ? c.rightCode==Tadpole.KEEP_GOING : c.rightBridgeEndpoint){
				tips.add(new Tip(c, true, tips.size()));
			}
		}
		return tips;
	}

	private void consider(ArrayList<Tip> tips, long sourceCode, Tip dest, int destAnchorTrim){
		final Tip source=tips.get((int)(sourceCode>>>32));
		final int sourceTrim=(int)sourceCode;
		final boolean sourceReverse=!source.right, destReverse=dest.right;
		final int sourceStart=source.contig.length()-sourceTrim-minOverlap;
		final int extraLimit=Math.min(maxOverlap-minOverlap, Math.min(sourceStart, destAnchorTrim));
		int extra=0;
		while(extra<extraLimit && baseAt(source.contig, sourceReverse, sourceStart-extra-1)
				==baseAt(dest.contig, destReverse, destAnchorTrim-extra-1)){extra++;}
		final int overlap=minOverlap+extra;
		final int destTrim=destAnchorTrim-extra;
		if(source==dest || overlap>=source.contig.length()-sourceTrim
				|| overlap>=dest.contig.length()-destTrim){return;}
		if(!matches(source, sourceTrim, dest, destTrim, overlap)){return;}
		if(source.contig==dest.contig){
			selfMatches++;
			if(overlap>source.bestOverlap){
				source.best=null;
				source.bestOverlap=overlap;
				source.bestSourceTrim=sourceTrim;
				source.bestDestTrim=destTrim;
				source.ambiguous=true;
			}else if(overlap==source.bestOverlap){source.ambiguous=true;}
			return;
		}
		exactCandidates++;
		if(overlap>source.bestOverlap){
			source.best=dest;
			source.bestOverlap=overlap;
			source.bestSourceTrim=sourceTrim;
			source.bestDestTrim=destTrim;
			source.ambiguous=false;
		}else if(overlap==source.bestOverlap){
			if(source.best!=dest){source.ambiguous=true;}
			else if(sourceTrim+destTrim<source.bestSourceTrim+source.bestDestTrim){
				source.bestSourceTrim=sourceTrim;
				source.bestDestTrim=destTrim;
			}
		}
	}

	private static boolean matches(Tip source, int sourceTrim, Tip dest, int destTrim, int overlap){
		final boolean sourceReverse=!source.right;
		final boolean destReverse=dest.right;
		final int sourceStart=source.contig.length()-sourceTrim-overlap;
		for(int i=0; i<overlap; i++){
			if(baseAt(source.contig, sourceReverse, sourceStart+i)!=baseAt(dest.contig, destReverse, destTrim+i)){
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
		final Edge ab=new Edge(a.contig.id, b.contig.id, 0, orientationAB, depth, null,
				pair.overlap, pair.aTrim, pair.bTrim);
		final Edge ba=new Edge(b.contig.id, a.contig.id, 0, orientationBA, depth, null,
				pair.overlap, pair.bTrim, pair.aTrim);
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
		if(!BubblePopper.verbose){return;}
		System.err.println((graphKEnds ? "Graph-k" : "Cross-k")+" tip overlaps: endpoints="+tips+
				", exactCandidates="+exactCandidates+
				", selfMatches="+selfMatches+
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
		int bestSourceTrim=0, bestDestTrim=0;
		boolean ambiguous=false;
	}

	static class Pair {
		Pair(Tip a_, Tip b_, int overlap_, int aTrim_, int bTrim_){
			a=a_; b=b_; overlap=overlap_; aTrim=aTrim_; bTrim=bTrim_;
		}
		final Tip a, b;
		final int overlap, aTrim, bTrim;
	}

	private final ArrayList<Contig> contigs;
	private final int minOverlap, maxOverlap, minContig;
	private final boolean graphKEnds;
	private long exactCandidates=0, selfMatches=0;
	private static final long AMBIGUOUS=Long.MAX_VALUE;
	private static final long HASH_MULT=0x9E3779B185EBCA87L;
}
