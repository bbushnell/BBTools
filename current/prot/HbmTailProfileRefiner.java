package prot;

import java.util.Arrays;
import java.util.HashSet;

/** Builds stop-anchored profiles without changing the full-protein refiner.
 * Inputs are borrowed; returned graphs and accounting arrays belong to the caller.
 * Rejected copies must return to the assignment pool, never remain in model counts.
 * @author Keqing
 */
public final class HbmTailProfileRefiner {

	/** BLOSUM seed placement on the longest member, followed by two profile passes. */
	public static Result build(byte[][] members, String[] ids, double[] background){
		validate(members, ids);
		int pivot=0;
		for(int i=1; i<members.length; i++){
			if(members[i].length>members[pivot].length || (members[i].length==members[pivot].length && ids[i].compareTo(ids[pivot])<0)){pivot=i;}
		}
		final AAGraph scaffold=new AAGraph(members[pivot].clone(), 0);
		final HbmPositionModel blosum=HbmPositionModel.deriveOne(scaffold, "consensus", .01, true, background, -4);
		final String[] reasons=new String[members.length];
		final Stage seed=accumulate(scaffold.pivot, blosum, members, reasons, "seed");
		if(seed.accepted<MIN_COPIES){return retired(reasons, "seed_below_min_copies");}
		return refine(seed.graph, members, background, reasons, seed.clipped, null);
	}

	/** Rebuilds changed membership from a frozen previous tail graph, not a new BLOSUM seed. */
	public static Result rebuild(AAGraph previous, byte[][] members, String[] ids, double[] background){
		return rebuild(previous, members, ids, background, null);
	}

	/** Rebuilds changed membership and optionally selects the final fit set with classifier gates. */
	public static Result rebuild(AAGraph previous, byte[][] members, String[] ids, double[] background, MembershipGate gate){
		validate(members, ids);
		require(previous!=null && previous.pad==0 && previous.pivot.length>=MIN_LENGTH && previous.pivot.length<=CAP,
			"A tail rebuild needs a frozen unpadded100..200-residue reference");
		return refine(previous, members, background, new String[members.length], 0, gate);
	}

	/** First pass selects a new consensus; the second accumulates its final frozen graph. */
	private static Result refine(AAGraph previous, byte[][] members, double[] background, String[] reasons, long seedClipped,
			MembershipGate gate){
		final HbmPositionModel profile=HbmPositionModel.deriveOne(previous, "logodds", .01, true, background, -4);
		final Stage growth=gate==null ? accumulate(previous.pivot, profile, members, reasons, "refine1") :
			accumulateInclusive(previous.pivot, profile, members, "refine1");
		if(growth.accepted<MIN_COPIES){return retired(gate==null ? reasons : growth.reasons, "refine1_below_min_copies");}
		if(gate!=null && growth.accepted<members.length){return failed(growth.reasons, "refine1_placement_failed");}
		growth.graph.trimDepthFraction=.1f;// Selection-only policy; final serialized graphs retain canonical zero trimming.
		final AAGraph.Traversal selected=growth.graph.traverseProfile();
		final byte[] untrimmed=selected.consensus();
		require(untrimmed.length>0 && selected.sourceType(untrimmed.length-1)==AAGraphNode.REF
			&& selected.sourcePosition(untrimmed.length-1)==previous.pivot.length-1,
			"Tail consensus must end at the observed terminal REF column; right trimming or extension changes the represented stop");
		if(untrimmed.length<MIN_LENGTH){return retired(reasons, "consensus_below_min_length");}
		final int trim=Math.max(0, untrimmed.length-CAP);
		final byte[] consensus=Arrays.copyOfRange(untrimmed, trim, untrimmed.length);
		final HbmPositionModel selectedProfile=HbmPositionModel.deriveTraversal(selected, background);
		final HbmPositionModel next=trim==0 ? selectedProfile : suffix(selectedProfile, trim);
		final Stage last=gate==null ? accumulate(consensus, next, members, reasons, "refine2") :
			accumulateMembership(consensus, next, members, gate);
		final String[] finalReasons=gate==null ? reasons : last.reasons;
		if(last.accepted<MIN_COPIES){return retired(finalReasons, "refine2_below_min_copies");}
		HbmBundleBuilder.validateGraph(new HbmBundleBuilder.FamilyInput("stop_tail", consensus, last.graph));
		return new Result(last.graph, finalReasons, last.accepted, last.residues, last.placed,
			last.clipped, seedClipped, growth.clipped, trim, null);
	}

	/** Copies exactly the retained selected score columns; no background is re-estimated. */
	private static HbmPositionModel suffix(HbmPositionModel model, int trim){
		require(trim>0 && trim<model.length(), "Score-frame clipping requires a nonempty retained suffix");
		final int n=model.length()-trim;
		final int[][] scores=new int[22][n];
		for(int a=0; a<22; a++){for(int j=0; j<n; j++){scores[a][j]=model.scoreAt((byte)a, j+trim);}}
		return HbmPositionModel.fromScores(scores);
	}

	/** All graph observations use suffixes whose trace consumes both right endpoints. */
	private static Stage accumulate(byte[] pivot, HbmPositionModel model, byte[][] members, String[] reasons, String phase){
		require(model.length()==pivot.length, "Tail profile scores and consensus must share one coordinate frame");
		final AAGraph graph=new AAGraph(pivot.clone(), 0);
		final Stage result=new Stage(graph);
		for(int i=0; i<members.length; i++){
			if(reasons[i]!=null){continue;}
			final byte[] query=members[i];
			final HbmPositionModel.TailResult alignment=model.alignTail(query, CAP, true);
			final String rejected=rejection(query, pivot, alignment);
			if(rejected!=null){reasons[i]=phase+"_"+rejected; continue;}
			final byte[] suffix=alignment.queryStart==0 ? query : Arrays.copyOfRange(query, alignment.queryStart, query.length);
			final int dropped=graph.addTrace(suffix, alignment.start, alignment.path);
			require(dropped==0, "Accepted tail trace lost terminal residues during graph accumulation");
			result.accepted++;
			result.residues=Math.addExact(result.residues, query.length);
			result.placed=Math.addExact(result.placed, suffix.length);
			result.clipped=Math.addExact(result.clipped, alignment.queryStart);
		}
		require(result.residues==result.placed+result.clipped, "Accepted tail residues must partition into placed and upstream-clipped residues");
		checkCounts(graph, result.placed);
		require(graph.ref[graph.ref.length-1].countSum==result.accepted+1,
			"Every accepted copy must observe the terminal REF column exactly once, in addition to its scaffold count");
		return result;
	}

	/** Places every assigned member with only right-end trace sanity; any failure is a bucket diagnostic. */
	private static Stage accumulateInclusive(byte[] pivot, HbmPositionModel model, byte[][] members, String phase){
		require(model.length()==pivot.length, "Tail profile scores and consensus must share one coordinate frame");
		final AAGraph graph=new AAGraph(pivot.clone(), 0);
		final Stage result=new Stage(graph, new String[members.length]);
		for(int i=0; i<members.length; i++){
			final byte[] query=members[i];
			final HbmPositionModel.TailResult alignment=model.alignTailTerminalPair(query, CAP, true);
			final String placed=rejection(query, pivot, alignment, TailPolicy.PLACEMENT_ONLY);
			if(placed!=null){result.reasons[i]=phase+"_"+placed; continue;}
			final byte[] suffix=alignment.queryStart==0 ? query : Arrays.copyOfRange(query, alignment.queryStart, query.length);
			final int dropped=graph.addTrace(suffix, alignment.start, alignment.path);
			require(dropped==0, "Accepted inclusive trace lost terminal residues during graph accumulation");
			result.accepted++;
			result.residues=Math.addExact(result.residues, query.length);
			result.placed=Math.addExact(result.placed, suffix.length);
			result.clipped=Math.addExact(result.clipped, alignment.queryStart);
		}
		require(result.residues==result.placed+result.clipped, "Inclusive residues must partition into placed and upstream-clipped residues");
		checkCounts(graph, result.placed);
		require(graph.ref[graph.ref.length-1].countSum==result.accepted+1,
			"Every placed inclusive copy must observe the terminal REF column exactly once, in addition to its scaffold count");
		return result;
	}

	/** Final fit-set check: classifier right-anchored BLOSUM gate, terminal-paired graph placement. */
	private static Stage accumulateMembership(byte[] pivot, HbmPositionModel model, byte[][] members, MembershipGate gate){
		require(model.length()==pivot.length, "Tail profile scores and consensus must share one coordinate frame");
		final AAGraph graph=new AAGraph(pivot.clone(), 0);
		final Stage result=new Stage(graph, new String[members.length]);
		final HbmPositionModel blosum=HbmPositionModel.fromConsensus(pivot);
		for(int i=0; i<members.length; i++){
			final byte[] query=members[i];
			final String rejected=gate.rejection(query, pivot, blosum.alignTail(query, CAP, true));
			if(rejected!=null){result.reasons[i]="membership_"+rejected; continue;}
			final HbmPositionModel.TailResult alignment=model.alignTailTerminalPair(query, CAP, true);
			final String placed=rejection(query, pivot, alignment, TailPolicy.PLACEMENT_ONLY);
			if(placed!=null){result.reasons[i]="placement_"+placed; continue;}
			final byte[] suffix=alignment.queryStart==0 ? query : Arrays.copyOfRange(query, alignment.queryStart, query.length);
			final int dropped=graph.addTrace(suffix, alignment.start, alignment.path);
			require(dropped==0, "Accepted membership trace lost terminal residues during graph accumulation");
			result.accepted++;
			result.residues=Math.addExact(result.residues, query.length);
			result.placed=Math.addExact(result.placed, suffix.length);
			result.clipped=Math.addExact(result.clipped, alignment.queryStart);
		}
		require(result.residues==result.placed+result.clipped, "Membership residues must partition into placed and upstream-clipped residues");
		checkCounts(graph, result.placed);
		require(graph.ref[graph.ref.length-1].countSum==result.accepted+1,
			"Every membership copy must observe the terminal REF column exactly once, in addition to its scaffold count");
		return result;
	}

	/** Construction gates from the v1 plan; these are not calibrated detection cutoffs. */
	static String rejection(byte[] query, byte[] pivot, HbmPositionModel.TailResult alignment){
		return rejection(query, pivot, alignment, TailPolicy.LEGACY_CONSTRUCTION);
	}

	/** Named tail gates keep final whole-pool assignment distinct from seed/rebuild construction. */
	public static String rejection(byte[] query, byte[] pivot, HbmPositionModel.TailResult alignment, TailPolicy policy){
		require(query!=null && pivot!=null && policy!=null, "Tail rejection needs sequences and a named policy");
		if(alignment==null){return "no_alignment";}
		require(alignment.path!=null && alignment.path.length>0 && alignment.score>0
			&& alignment.queryStart>=0 && alignment.start>=0
			&& alignment.queryEnd==query.length-1 && alignment.end==pivot.length-1,
			"Tail construction requires a positive recorded trace ending at both sequence endpoints");
		int q=alignment.queryStart, r=alignment.start, paired=0, matches=0, terminal=0;
		for(byte op : alignment.path){
			if(op=='m'){
				require(q<query.length && r<pivot.length, "Tail paired column exceeds its source sequences");
				paired++; if(query[q]==pivot[r]){matches++;}
				if(q>=query.length-20 && r>=pivot.length-20){terminal++;}
				q++; r++;
			}else if(op=='I'){q++;}
			else if(op=='D'){r++;}
			else{throw new IllegalArgumentException("Unknown tail operation: "+(char)op);}
			require(q<=query.length && r<=pivot.length, "Tail construction trace overran query or pivot");
		}
		require(q==query.length && r==pivot.length, "Tail construction trace did not consume its declared right ends");
		if(alignment.path[alignment.path.length-1]!='m'){return "unpaired_end";}
		if(policy.checkSupport && paired<policy.minPaired(query.length, pivot.length)){return "short_support";}
		if(policy.checkIdentity && 10L*matches<3L*alignment.path.length){return "low_identity";}
		if(policy.checkTerminal && terminal<10){return "short_terminal_support";}
		return null;
	}

	private static void checkCounts(AAGraph graph, long expected){
		long observed=0;
		for(int i=0; i<graph.ref.length; i++){
			observed=Math.addExact(observed, graph.ref[i].countSum-1L);
			for(AAGraphNode node=graph.ref[i].insEdge; node!=null; node=node.insEdge){observed=Math.addExact(observed, node.countSum);}
			for(AAGraphNode node=graph.del[i].insEdge; node!=null; node=node.insEdge){observed=Math.addExact(observed, node.countSum);}
		}
		require(observed==expected, "Tail graph residue conservation failed: observed="+observed+" expected="+expected);
	}

	private static void validate(byte[][] members, String[] ids){
		require(members!=null && ids!=null && members.length==ids.length && members.length>0,
			"Tail family requires matched nonempty member and identifier arrays");
		final HashSet<String> unique=new HashSet<String>();
		for(int i=0; i<members.length; i++){
			require(ids[i]!=null && !ids[i].isEmpty() && unique.add(ids[i]), "Tail member IDs must be unique and nonempty");
			for(int j=0; j<ids[i].length(); j++){require(!Character.isWhitespace(ids[i].charAt(j)), "Tail member IDs must be single FASTA tokens so external clustering cannot collapse distinct copies");}
			require(members[i]!=null && members[i].length>=MIN_LENGTH && members[i].length<=CAP,
				"Tail members must contain100..200 encoded amino acids");
			for(byte residue : members[i]){require(residue>=0 && residue<20, "Tail members cannot contain stops, ambiguity or non-amino-acid codes");}
		}
	}

	private static Result retired(String[] reasons, String reason){
		for(int i=0; i<reasons.length; i++){if(reasons[i]==null){reasons[i]=reason;}}
		return new Result(null, reasons, 0, 0, 0, 0, 0, 0, 0, reason);
	}
	private static Result failed(String[] reasons, String reason){
		for(int i=0; i<reasons.length; i++){if(reasons[i]==null){reasons[i]=reason;}}
		return new Result(null, reasons, 0, 0, 0, 0, 0, 0, 0, reason);
	}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final class Stage {
		Stage(AAGraph graph_){this(graph_, null);}
		Stage(AAGraph graph_, String[] reasons_){graph=graph_; reasons=reasons_;}
		final AAGraph graph;
		final String[] reasons;
		int accepted;
		long residues, placed, clipped;
	}

	/** Null reason means retained. A retired group has no graph and no retained copies. */
	public static final class Result {
		Result(AAGraph graph_, String[] reasons_, int retained_, long residues_, long placed_, long clipped_,
				long seedClipped_, long growthClipped_, int consensusTrim_, String retiredReason_){
			graph=graph_; reasons=reasons_.clone(); retained=retained_;
			retainedInputResidues=residues_; placedResidues=placed_; upstreamClippedResidues=clipped_;
			seedClippedResidues=seedClipped_; growthClippedResidues=growthClipped_;
			consensusPrefixTrim=consensusTrim_; retiredReason=retiredReason_;
			int count=0; for(String reason : reasons){if(reason==null){count++;}}
			require(count==retained && (graph!=null)==(retained>=MIN_COPIES), "Published tail membership must equal the graph's accepted-copy roster");
		}
		public final AAGraph graph;
		public final String[] reasons;
		public final int retained, consensusPrefixTrim;
		public final long retainedInputResidues, placedResidues, upstreamClippedResidues, seedClippedResidues, growthClippedResidues;
		public final String retiredReason;
	}
	private HbmTailProfileRefiner(){}
	public static final int MIN_LENGTH=100, CAP=200, MIN_COPIES=3;

	public static enum TailPolicy {
		LEGACY_CONSTRUCTION(false, true, true, true),
		FINAL_BIDIRECTIONAL(true, true, true, true),
		PLACEMENT_ONLY(false, false, false, false);

		TailPolicy(boolean bidirectional_, boolean checkSupport_, boolean checkIdentity_, boolean checkTerminal_){
			bidirectional=bidirectional_; checkSupport=checkSupport_; checkIdentity=checkIdentity_;
			checkTerminal=checkTerminal_;
		}
		int minPaired(int queryLength, int pivotLength){
			require(queryLength>0 && pivotLength>0, "Tail policy lengths must be positive");
			final int length=bidirectional ? Math.max(queryLength, pivotLength) : Math.min(queryLength, pivotLength);
			return Math.max(70, (7*length+9)/10);
		}
		private final boolean bidirectional, checkSupport, checkIdentity, checkTerminal;
	}

	public static final class MembershipGate {
		public MembershipGate(float identityLower_, float identityUpper_, float blosumFinal_){
			require(Float.isFinite(identityLower_) && identityLower_>=0 && identityLower_<=100 &&
				Float.isFinite(identityUpper_) && identityUpper_>=identityLower_ && identityUpper_<=Math.nextUp(110f) &&
				Float.isFinite(blosumFinal_), "Invalid membership gate");
			identityLower=identityLower_; identityUpper=identityUpper_; blosumFinal=blosumFinal_;
		}
		String rejection(byte[] query, byte[] pivot, HbmPositionModel.TailResult alignment){
			if(alignment==null){return "no_alignment";}
			require(alignment.path!=null && alignment.path.length>0
				&& alignment.queryStart>=0 && alignment.start>=0
				&& alignment.queryEnd==query.length-1 && alignment.end==pivot.length-1,
				"Membership gate requires a recorded positive trace ending at both sequence endpoints");
			final int self=selfScore(query);
			int q=alignment.queryStart, r=alignment.start, matches=0, score=0;
			for(byte op : alignment.path){
				if(op=='m'){
					require(q<query.length && r<pivot.length, "Membership paired column exceeds its source sequences");
					if(query[q]==pivot[r]){matches++;}
					score+=(int)(HbmPositionModel.SCALE*Blosum62.score(query[q++], pivot[r++]));
				}else if(op=='I'){q++; score-=HbmPositionModel.GAP;}
				else if(op=='D'){r++; score-=HbmPositionModel.GAP;}
				else{throw new IllegalArgumentException("Unknown BLOSUM trace operation: "+op);}
				require(q<=query.length && r<=pivot.length, "Membership trace overran query or consensus");
			}
			require(q==query.length && r==pivot.length && score==alignment.score,
				"Membership trace consumption/score differs from bottom-right BLOSUM alignment");
			final float identity=(float)(100.0*matches/alignment.path.length);
			final float blosum=(float)(score/(HbmPositionModel.SCALE*(double)self));
			if(identity>identityUpper){return null;}
			if(identity<identityLower){return "identity_reject";}
			return blosum>blosumFinal ? null : "blosum_reject";
		}
		private static int selfScore(byte[] query){
			int self=0; for(byte aa : query){self+=Blosum62.score(aa, aa);}
			require(self>0, "Membership BLOSUM denominator must be positive");
			return self;
		}
		private final float identityLower, identityUpper, blosumFinal;
	}
}
