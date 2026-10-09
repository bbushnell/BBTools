package prot;

/** Two frozen-profile passes for one family; no assignment, serialization, or corpus scheduling. */
public final class HbmProfileRefiner {
	private HbmProfileRefiner(){}

	/**
	 * Realigns members to the old profile, derives a trimmed consensus/profile pair,
	 * then accumulates a fresh pad0 final graph on that frozen new consensus.
	 * Inputs are borrowed read-only for the call. Background probabilities are fixed
	 * across both passes; scores use beta .01, clipping [-4,11], and linear gap4.
	 * The growth pass uses caller-specified neutral padding and depth trim .1.
	 * Insufficient growth padding fails loudly; the caller can retry with more.
	 * Final pad0 overhang exclusions preserve the existing HBM format's edge policy
	 * and are counted explicitly. The returned graph is owned by the caller.
	 */
	public static Result refine(final AAGraph initial, final byte[][] members,
			final double[] background, final int padding){
		if(initial.pad!=0 || members.length==0 || padding<0){
			throw new IllegalArgumentException("Refinement needs an unpadded initial HBM, members, and nonnegative growth padding");
		}
		final byte[] original=initial.pivot.clone();
		final HbmPositionModel oldProfile=HbmPositionModel.deriveOne(initial,"logodds",0.01,true,background,-4);
		oldProfile.requireConsensus(original);
		final HbmPositionModel growthProfile=oldProfile.padded(padding);
		final AAGraph growth=new AAGraph(original,padding);
		growth.trimDepthFraction=0.1f;
		growthProfile.requireConsensus(growth.pivot);
		long residues=0;
		for(byte[] member : members){
			final HbmPositionModel.Result alignment=growthProfile.align(member,true);
			final int dropped=growth.addTrace(member,alignment.start,alignment.path);
			if(dropped!=0){throw new PaddingException(padding,dropped);}
			residues=Math.addExact(residues,member.length);
		}
		checkCounts(growth,residues);
		final AAGraph.Traversal selected=growth.traverseProfile();
		final byte[] consensus=selected.consensus();
		if(consensus.length==0){throw new IllegalArgumentException("Profile refinement produced an empty consensus");}
		final HbmPositionModel next=HbmPositionModel.deriveTraversal(selected,background);
		next.requireConsensus(consensus);
		final AAGraph result=new AAGraph(consensus,0);
		long excluded=0;
		for(byte[] member : members){
			final HbmPositionModel.Result alignment=next.align(member,true);
			excluded=Math.addExact(excluded,result.addTrace(member,alignment.start,alignment.path));
		}
		checkCounts(result,residues-excluded);
		HbmBundleBuilder.validateGraph(new HbmBundleBuilder.FamilyInput("refinement",consensus,result));
		return new Result(result,selected,members.length,residues,excluded);
	}

	/** Every placed query residue contributes once, in addition to one scaffold observation per column. */
	private static void checkCounts(final AAGraph graph, final long expected){
		long count=0;
		for(int i=0; i<graph.ref.length; i++){
			count=Math.addExact(count,graph.ref[i].countSum-1L);
			for(AAGraphNode node=graph.ref[i].insEdge; node!=null; node=node.insEdge){count=Math.addExact(count,node.countSum);}
			for(AAGraphNode node=graph.del[i].insEdge; node!=null; node=node.insEdge){count=Math.addExact(count,node.countSum);}
		}
		if(count!=expected){throw new AssertionError("Refined graph residue conservation failed: placed="+count+" expected="+expected);}
	}

	/** The final graph is mutable caller-owned output; the selected growth columns are immutable. */
	public static final class Result{
		private Result(AAGraph graph_,AAGraph.Traversal selected_,int members_,long residues_,long excluded_){
			graph=graph_;selected=selected_;members=members_;residues=residues_;excludedTerminalResidues=excluded_;
		}
		public final AAGraph graph;
		public final AAGraph.Traversal selected;
		public final int members;
		public final long residues,excludedTerminalResidues;
	}
	/** A caller may retry with larger neutral padding; other failures must not be mistaken for this condition. */
	public static final class PaddingException extends IllegalArgumentException{
		private static final long serialVersionUID=1L;
		PaddingException(int padding,int excluded_){
			super("Growth padding="+padding+" excluded "+excluded_+" terminal residues; retry this family with more padding");
			excluded=excluded_;
		}
		public final int excluded;
	}
}
