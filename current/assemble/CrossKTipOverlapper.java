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
		this(contigs_, minOverlap_, maxOverlap_, graphKEnds_, minContig_, -1);
	}

	/** A negative mismatch allowance preserves legacy unchecked terminal trimming. */
	CrossKTipOverlapper(final ArrayList<Contig> contigs_, final int minOverlap_, final int maxOverlap_,
			final boolean graphKEnds_, final int minContig_, final int maxMismatches_){
		this(contigs_, minOverlap_, maxOverlap_, graphKEnds_, minContig_, maxMismatches_, false);
	}

	/** Optionally requires historical dead-end evidence in addition to current tip eligibility. */
	CrossKTipOverlapper(final ArrayList<Contig> contigs_, final int minOverlap_, final int maxOverlap_,
			final boolean graphKEnds_, final int minContig_, final int maxMismatches_, final boolean deadEndsOnly_){
		this(contigs_, minOverlap_, maxOverlap_, graphKEnds_, minContig_, maxMismatches_, deadEndsOnly_, false);
	}

	/** Optionally vetoes a selected pair when a longer conflicting placement exists in this pass. */
	CrossKTipOverlapper(final ArrayList<Contig> contigs_, final int minOverlap_, final int maxOverlap_,
			final boolean graphKEnds_, final int minContig_, final int maxMismatches_,
			final boolean deadEndsOnly_, final boolean rejectConflicts_){
		if(minOverlap_<1 || maxOverlap_<minOverlap_){
			throw new IllegalArgumentException("Invalid cross-k overlap range: "+minOverlap_+"-"+maxOverlap_);
		}
		contigs=contigs_;
		minOverlap=minOverlap_;
		maxOverlap=maxOverlap_;
		graphKEnds=graphKEnds_;
		minContig=minContig_;
		if(maxMismatches_<-1){throw new IllegalArgumentException("Invalid fusion mismatch allowance: "+maxMismatches_);}
		maxMismatches=maxMismatches_;
		deadEndsOnly=deadEndsOnly_;
		rejectConflicts=rejectConflicts_;
		if(rejectConflicts && maxMismatches<0){
			throw new IllegalArgumentException("fuseconflicts requires fusemaxmismatches>=0.");
		}
	}

	/** Adds reciprocal overlap edges and returns the number of acyclic pairs added. */
	int addEdges(){
		validateCoverageRatio(maxCoverageRatio);
		final ArrayList<Tip> tips=makeTips();
		if(BubblePopper.verbose){
			System.err.println("FusionPass\t"+minOverlap+"\t"+maxOverlap+"\t"+graphKEnds);
		}
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

		findCandidates(tips, terminalMap, duplicateMap, outgoingPower, false);
		//Freeze best partners before checking their alternatives: enumeration order
		//must not let an earlier rejected placement disappear behind a later best.
		if(rejectConflicts){findCandidates(tips, terminalMap, duplicateMap, outgoingPower, true);}

		//Log final state only after every candidate has had a chance to replace an earlier best.
		if(BubblePopper.verbose){for(Tip tip : tips){printTipState(tip);}}
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

		//Keep the original candidate graph for cycle rejection. Removing a conflict
		//first could break a cycle and inadvertently admit other previously rejected edges.
		final boolean[] cyclic=cyclicComponents(reciprocal);
		int added=0, cycleRejected=0;
		for(int i=0; i<reciprocal.size(); i++){
			final Pair pair=reciprocal.get(i);
			if(collector!=null){
				collector.record(pair.a.contig, !pair.a.right, pair.aTrim,
						pair.b.contig, pair.b.right, pair.bTrim, pair.overlap);
			}
			if(cyclic[i]){cycleRejected++; continue;}
			if(pair.a.conflictingPlacement || pair.b.conflictingPlacement){conflictRejected++; continue;}
			if(!compatibleCoverage(pair.a.contig.coverage, pair.b.contig.coverage, maxCoverageRatio)){
				coverageRejected++;
				continue;
			}
			if(support!=null && !support.supported(pair.a.contig, !pair.a.right, pair.aTrim,
					pair.b.contig, pair.b.right, pair.bTrim, pair.overlap)){
				supportRejected++;
				continue;
			}
			if(neural!=null && !neural.accepts(pair.a.contig, !pair.a.right, pair.aTrim,
					pair.b.contig, pair.b.right, pair.bTrim, pair.overlap)){
				neuralRejected++;
				continue;
			}
			if(BubblePopper.verbose){
				System.err.println("FusionSelected\t"+minOverlap+"\t"+pair.a.index+"\t"+pair.b.index+
						"\t"+pair.overlap+"\t"+pair.aTrim+"\t"+pair.bTrim);
			}
			addEdges(pair);
			added++;
		}
		printSummary(tips.size(), ambiguous, reciprocal.size(), cycleRejected);
		return added;
	}

	/** Zero disables the experimental guard; an enabled maximum ratio cannot be below one. */
	static void validateCoverageRatio(final float ratio){
		if(!Float.isFinite(ratio) || ratio<0 || (ratio>0 && ratio<1)){
			throw new IllegalArgumentException("fusecoverageratio must be zero or finite and at least one: "+ratio);
		}
	}

	/**
	 * Compares frozen whole-contig mean depths, not a second read-count table.
	 * This deliberately broad experimental veto may remove correct one-sided
	 * repeat attachments. Missing depth is not evidence of compatible copy number.
	 */
	static boolean compatibleCoverage(final float a, final float b, final float ratio){
		assert(Float.isFinite(ratio) && (ratio==0 || ratio>=1)) :
				"The coverage veto requires a validated multiplicative limit or its zero/off sentinel.";
		if(ratio==0){return true;}
		assert(Float.isFinite(a) && Float.isFinite(b) && a>=0 && b>=0) :
				"Contig mean depths must be finite and nonnegative before comparing copy-number evidence: "+a+", "+b;
		final double low=Math.min(a, b), high=Math.max(a, b);
		return low>0 && high<=low*ratio;
	}

	/** Reuses the anchor maps for either best selection or a read-only replay of selected alternatives. */
	private void findCandidates(final ArrayList<Tip> tips, final LongLongHashMap2 terminalMap,
			final LongLongHashMap2 duplicateMap, final long outgoingPower, final boolean conflictsOnly){
		assert(!conflictsOnly || rejectConflicts) : "Conflict replay needs the explicitly enabled fusion policy.";
		final int maxTrim=maxOverlap-minOverlap+1;
		for(Tip dest : tips){
			if(conflictsOnly && (dest.best==null || dest.ambiguous || dest.best.ambiguous || dest.best.best!=dest)){continue;}
			final Contig c=dest.contig;
			final boolean reverse=dest.right;
			final int trimLimit=Math.min(maxTrim, c.length()-minOverlap-1);
			long hash=hash(c, reverse, 0, minOverlap);
			for(int destTrim=0; destTrim<=trimLimit; destTrim++){
				final long code=terminalMap.get(hash);
				if(!conflictsOnly && BubblePopper.verbose && code==AMBIGUOUS){
					System.err.println("FusionBlocked\t"+minOverlap+"\t"+dest.index+"\t"+destTrim+"\t"+hash);
				}
				if(code>0 && code!=AMBIGUOUS){
					consider(tips, code-1, dest, destTrim, conflictsOnly);
					final long duplicate=duplicateMap.get(hash);
					if(duplicate>0){consider(tips, duplicate-1, dest, destTrim, conflictsOnly);}
				}
				if(destTrim<trimLimit){
					final int oldBase=baseAt(c, reverse, destTrim)&0xFF;
					final int newBase=baseAt(c, reverse, destTrim+minOverlap)&0xFF;
					hash=(hash-(oldBase+1)*outgoingPower)*HASH_MULT+(newBase+1);
				}
			}
		}

	}

	/** Stores one anchor per tip and at most two tips per hash; repeated or third tips are ambiguous. */
	private void addAnchor(LongLongHashMap2 primary, LongLongHashMap2 duplicate,
			long hash, long code){
		if(BubblePopper.verbose){
			System.err.println("FusionAnchor\t"+minOverlap+"\t"+(code>>>32)+"\t"+(int)code+"\t"+hash);
		}
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
			if((!deadEndsOnly || c.leftFusionEndpoint)
					&& (graphKEnds ? c.leftCode==Tadpole.KEEP_GOING : c.leftBridgeEndpoint)){
				tips.add(new Tip(c, false, tips.size()));
			}
			if((!deadEndsOnly || c.rightFusionEndpoint)
					&& (graphKEnds ? c.rightCode==Tadpole.KEEP_GOING : c.rightBridgeEndpoint)){
				tips.add(new Tip(c, true, tips.size()));
			}
		}
		return tips;
	}

	private void consider(ArrayList<Tip> tips, long sourceCode, Tip dest, int destAnchorTrim, boolean conflictsOnly){
		final Tip source=tips.get((int)(sourceCode>>>32));
		if(conflictsOnly && (source.best!=dest || dest.best!=source || source.conflictingPlacement)){return;}
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
		if(conflictsOnly){
			final long span=(long)overlap+sourceTrim+destTrim;
			final long bestSpan=(long)source.bestOverlap+source.bestSourceTrim+source.bestDestTrim;
			//Uncovered context is lack of evidence, not a contradiction. Compare the
			//full placement, including trims, rather than just the exact anchor length.
			if(span>bestSpan && span<=source.contig.length() && span<=dest.contig.length()
					&& !compatibleTrimmedFlanks(source.contig, sourceReverse, sourceTrim,
							dest.contig, destReverse, destTrim, overlap, maxMismatches)){
				source.conflictingPlacement=true;
				if(BubblePopper.verbose){
					System.err.println("FusionConflict\t"+minOverlap+"\t"+source.index+"\t"+dest.index+
							"\t"+overlap+"\t"+sourceTrim+"\t"+destTrim+"\t"+bestSpan);
				}
			}
			return;
		}
		if(source.contig==dest.contig){
			if(BubblePopper.verbose){printCandidate("self", source, dest, overlap, sourceTrim, destTrim);}
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
		//Self-repeat matches above remain conservative ambiguity blockers. For an
		//external partner, validate discarded context before best-partner selection.
		//Anchor offsets are not final trims: a longer exact terminal overlap can
		//start at an internal seed and extend back to destTrim=0.
		if(!allowTrim && (sourceTrim>0 || destTrim>0)){return;}
		if(maxMismatches>=0 && !compatibleTrimmedFlanks(source.contig, sourceReverse, sourceTrim,
				dest.contig, destReverse, destTrim, overlap, maxMismatches)){
			flankRejected++;
			if(BubblePopper.verbose){printCandidate("flank", source, dest, overlap, sourceTrim, destTrim);}
			return;
		}
		if(BubblePopper.verbose){printCandidate("eligible", source, dest, overlap, sourceTrim, destTrim);}
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

	/**
	 * Checks the exact anchor and both discarded flanks against the continuation.
	 * Rechecking the anchor also protects merges after other edits. One mismatch budget
	 * covers both sides; unknown bases are disagreements, not matching evidence.
	 * Uncovered trimmed bases are rejected rather than silently left unvalidated.
	 */
	static boolean compatibleTrimmedFlanks(final Contig source, final boolean sourceReverse,
			final int sourceTrim, final Contig dest, final boolean destReverse, final int destTrim,
			final int overlap, final int allowance){
		assert(source!=null && dest!=null && allowance>=0 && sourceTrim>=0 && destTrim>=0 && overlap>0) :
				"Fusion flank validation requires contigs, a positive anchor and nonnegative trims/budget.";
		final int sourceEnd=source.length()-sourceTrim;
		final int sourceStart=sourceEnd-overlap;
		final int destEnd=destTrim+overlap;
		//BubblePopper.exactOverlap also requires non-overlap sequence on both sides.
		if(sourceStart<=0 || destEnd>=dest.length()){return false;}
		if(destTrim>sourceStart || sourceTrim>dest.length()-destEnd){return false;}
		for(int i=0; i<overlap; i++){
			if(baseAt(source, sourceReverse, sourceStart+i)!=baseAt(dest, destReverse, destTrim+i)){
				return false;
			}
		}
		int mismatches=0;
		for(int i=0; i<destTrim; i++){
			final byte a=baseAt(source, sourceReverse, sourceStart-destTrim+i);
			final byte b=baseAt(dest, destReverse, i);
			if(a!=b || AminoAcid.baseToNumber[a]<0){
				if(++mismatches>allowance){return false;}
			}
		}
		for(int i=0; i<sourceTrim; i++){
			final byte a=baseAt(source, sourceReverse, sourceEnd+i);
			final byte b=baseAt(dest, destReverse, destEnd+i);
			if(a!=b || AminoAcid.baseToNumber[a]<0){
				if(++mismatches>allowance){return false;}
			}
		}
		return true;
	}

	private static void addEdges(Pair pair){
		final Tip a=pair.a, b=pair.b;
		if(BubblePopper.verbose){printPair(pair);}
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

	/** Logs oriented terminal sequences so overlap pairs survive later contig renumbering. */
	private static void printPair(final Pair pair){
		final Tip a=pair.a, b=pair.b;
		assert(pair.overlap>0) : "An accepted exact-overlap pair must retain a positive shared sequence.";
		final int aStart=Math.max(0, a.contig.length()-pair.aTrim-pair.overlap-200);
		final int bStop=Math.min(b.contig.length(), pair.bTrim+pair.overlap+200);
		final StringBuilder source=new StringBuilder(a.contig.length()-aStart);
		final StringBuilder dest=new StringBuilder(bStop);
		for(int i=aStart; i<a.contig.length(); i++){
			source.append((char)baseAt(a.contig, !a.right, i));
		}
		for(int i=0; i<bStop; i++){
			dest.append((char)baseAt(b.contig, b.right, i));
		}
		System.err.println("OverlapPair\t"+a.contig.id+"\t"+b.contig.id+"\t"+
				a.contig.length()+"\t"+b.contig.length()+"\t"+pair.overlap+"\t"+
				pair.aTrim+"\t"+pair.bTrim+"\t"+aStart+"\t"+source+"\t"+dest);
	}

	/** Logs a candidate and its prior best state without retaining any cross-pass cache. */
	private void printCandidate(final String reason, final Tip source, final Tip dest,
			final int overlap, final int sourceTrim, final int destTrim){
		assert(overlap>=minOverlap && overlap<=maxOverlap) :
				"Fusion diagnostics must describe the candidate admitted by this pass's overlap bounds.";
		System.err.println("FusionCandidate\t"+minOverlap+"\t"+source.index+"\t"+dest.index+
				"\t"+overlap+"\t"+sourceTrim+"\t"+destTrim+"\t"+reason+
				"\t"+(source.best==null ? -1 : source.best.index)+"\t"+source.bestOverlap+"\t"+source.ambiguous);
	}

	/**
	 * Logs the final selection state and outward-oriented end sequence before merging.
	 * End strings allow cross-pass matching despite renumbering; duplicate strings must
	 * remain ambiguous in analysis. These strings are temporary diagnostic output only.
	 */
	private void printTipState(final Tip tip){
		assert(tip.contig.length()>minOverlap) : "Only nonconsumed contig ends enter overlap discovery.";
		final StringBuilder end=new StringBuilder(200);
		for(int i=Math.max(0, tip.contig.length()-200); i<tip.contig.length(); i++){
			end.append((char)baseAt(tip.contig, !tip.right, i));
		}
		System.err.println("FusionTip\t"+minOverlap+"\t"+tip.index+"\t"+tip.contig.id+
				"\t"+(tip.right ? "R" : "L")+"\t"+tip.contig.length()+
				"\t"+(tip.best==null ? -1 : tip.best.index)+"\t"+tip.bestOverlap+
				"\t"+tip.bestSourceTrim+"\t"+tip.bestDestTrim+"\t"+tip.ambiguous+"\t"+end);
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
		if(!BubblePopper.verbose && maxMismatches<0 && support==null && neural==null && maxCoverageRatio==0){return;}
		System.err.println((graphKEnds ? "Graph-k" : "Cross-k")+" tip overlaps: endpoints="+tips+
				", exactCandidates="+exactCandidates+
				", selfMatches="+selfMatches+
				", flankRejected="+flankRejected+
				", ambiguous="+ambiguous+", reciprocal="+reciprocal+
				", cycleRejected="+cycleRejected+", conflictRejected="+conflictRejected+
				", supportRejected="+supportRejected+
				(maxCoverageRatio==0 ? "" : ", coverageRejected="+coverageRejected)+
				(neural==null ? "" : ", neuralRejected="+neuralRejected)+
				", added="+(reciprocal-cycleRejected-conflictRejected-supportRejected-coverageRejected-neuralRejected)+".");
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
		/** A frozen best pair has a covered, larger placement with contradictory flanks. */
		boolean conflictingPlacement=false;
	}

	static class Pair {
		Pair(Tip a_, Tip b_, int overlap_, int aTrim_, int bTrim_){
			a=a_; b=b_; overlap=overlap_; aTrim=aTrim_; bTrim=bTrim_;
		}
		final Tip a, b;
		final int overlap, aTrim, bTrim;
	}

	private final ArrayList<Contig> contigs;
	private final int minOverlap, maxOverlap, minContig, maxMismatches;
	private final boolean graphKEnds;
	private final boolean deadEndsOnly;
	private final boolean rejectConflicts;
	private long exactCandidates=0, selfMatches=0;
	private long flankRejected=0;
	private int conflictRejected=0;
	/** Optional immutable-read-table checker, installed before serial edge discovery. */
	FusionKmerSupport support;
	/** Optional read-only census of every frozen reciprocal pair before safety vetoes. */
	FusionJoinCollector collector;
	/** Optional post-guard veto; never changes the frozen partner or cycle graph. */
	FusionNeuralGate neural;
	/** Controls external-pair eligibility, not conservative self/anchor/conflict evidence. */
	boolean allowTrim=true;
	/** Experimental whole-contig depth-discontinuity veto; zero preserves legacy behavior. */
	float maxCoverageRatio=0;
	private int supportRejected=0;
	private int coverageRejected=0;
	private int neuralRejected=0;
	private static final long AMBIGUOUS=Long.MAX_VALUE;
	private static final long HASH_MULT=0x9E3779B185EBCA87L;
}
