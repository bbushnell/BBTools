package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

/** Exact orientation/endpoint propagation checks for multi-K postpass products. @author Fischl */
public final class MultiKBranchRescueTest {

	/** Exercises every pair orientation and unmerged closed components. */
	public static void main(final String[] args){
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Tests require -ea.");}
		if(args.length>0 && args[0].startsWith("scope-")){
			testCoverageScope(args);
			return;
		}
		for(int mask=0; mask<4; mask++){testJoin((mask&1)!=0, (mask&2)!=0);}
		testCycle();
		testCoverageArguments();
		System.out.println("MULTIK_BRANCH_METADATA_PASS orientations=4 cycle_unitigs=2");
	}

	/** Invalid initial ratios fail during parsing, before any input or writer opens. */
	private static void testCoverageArguments(){
		final String[] valid={"0", "1", "1.75", "3"};
		for(String value : valid){new MultiKBranchRescue(coverageArgs(value));}
		final String[] invalid={"-1", "0.5", "NaN", "Infinity", "-Infinity", "bad"};
		for(String value : invalid){
			boolean failed=false;
			try{new MultiKBranchRescue(coverageArgs(value));}
			catch(IllegalArgumentException expected){failed=true;}
			assert(failed) : "Invalid initial coverage ratio escaped the pre-I/O validation: "+value;
		}
		System.out.println("MULTIK_INITIAL_COVERAGE_ARGUMENTS_PASS valid="+valid.length+" invalid="+invalid.length);
	}

	/** Supplies parse-only paths; this helper never opens files or loads a count table. */
	private static String[] coverageArgs(final String ratio){
		assert(ratio!=null) : "Argument fixtures must supply an explicit ratio spelling.";
		return new String[]{"in=unused-input.fa", "out=unused-output.fa", "k=96,124,64,32",
				"assemblek=96", "initialmergecoverageratio="+ratio};
	}

	/** Runs one supplied fixture in-process so static-policy leakage cannot hide at JVM exit. */
	private static void testCoverageScope(final String[] args){
		final boolean expectFailure=args[0].equals("scope-failure");
		if(!expectFailure && !args[0].equals("scope-success")){
			throw new IllegalArgumentException("Expected scope-success or scope-failure followed by driver arguments.");
		}
		final BubblePopper.CoverageGate original=BubblePopper.directCoverageGate;
		final BubblePopper.CoverageGate sentinel=new BubblePopper.CoverageGate(17);
		final GraphBranchGate originalBranch=BubblePopper.directBranchGate;
		final GraphBranchGate branchSentinel=new GraphBranchGate();
		BubblePopper.directCoverageGate=sentinel;
		BubblePopper.directBranchGate=branchSentinel;
		boolean failed=false;
		try{
			MultiKBranchRescue.main(Arrays.copyOfRange(args, 1, args.length));
		}catch(IllegalArgumentException expected){
			if(!expectFailure || !expected.getMessage().contains("Fusion neural gate requires 109 inputs")){throw expected;}
			failed=true;
		}finally{
			final BubblePopper.CoverageGate restored=BubblePopper.directCoverageGate;
			final GraphBranchGate restoredBranch=BubblePopper.directBranchGate;
			BubblePopper.directCoverageGate=original;
			BubblePopper.directBranchGate=originalBranch;
			assert(restored==sentinel) : "An initial cleanup gate leaked out of the driver on success or failure.";
			assert(restoredBranch==branchSentinel) : "Initial/final branch gates must restore their caller on success or failure.";
		}
		assert(failed==expectFailure) : "The failure fixture must reach fusion-model validation after prepare() installs the initial gate.";
		System.out.println("MULTIK_INITIAL_COVERAGE_SCOPE_PASS failure="+failed);
		System.out.println("MULTIK_INITIAL_BRANCH_SCOPE_PASS failure="+failed);
	}

	/** The start's physical outside end and finish's outside end must survive inversion. */
	private static void testJoin(final boolean reverseA, final boolean reverseB){
		final Contig a=contig("AAAACCCC", 0), b=contig("CCCCGGGG", 1);
		a.leftCode=Tadpole.F_BRANCH; a.rightCode=Tadpole.B_BRANCH;
		b.leftCode=Tadpole.D_BRANCH; b.rightCode=Tadpole.LOOP;
		a.leftRatio=7; a.rightRatio=2; b.leftRatio=3; b.rightRatio=11;
		a.leftBridgeEndpoint=true; b.rightFusionEndpoint=true;
		if(reverseA){a.rcomp();}
		if(reverseB){b.rcomp();}
		a.graphLeftLimited=!reverseA; a.graphRightLimited=reverseA;
		b.graphLeftLimited=reverseB; b.graphRightLimited=!reverseB;
		final byte[] originalA=a.bases.clone(), originalB=b.bases.clone();
		final ArrayList<Contig> input=new ArrayList<Contig>(Arrays.asList(a, b));
		final int[] mates={-1,-1,-1,-1};
		final int ae=reverseA ? 0 : 1, be=reverseB ? 3 : 2;
		mates[ae]=be; mates[be]=ae;
		final SameKJoinMerge.Result result=SameKJoinMerge.merge(input, mates, 5);
		assert(result.joins==1 && result.contigs.size()==1) : "One reciprocal connection must join exactly two unitigs.";
		final Contig c=result.contigs.get(0);
		assert(Arrays.equals(c.bases, "AAAACCCCGGGG".getBytes(StandardCharsets.US_ASCII))) : "K-1 overlap spelling changed with orientation.";
		assert(c.leftCode==Tadpole.F_BRANCH && c.rightCode==Tadpole.LOOP && c.leftRatio==7 && c.rightRatio==11) :
				"Surviving outside branch codes and ratios must not become default dead ends.";
		assert(c.leftBridgeEndpoint && !c.rightBridgeEndpoint && !c.leftFusionEndpoint && c.rightFusionEndpoint) :
				"Bridge/fusion eligibility belongs to physical surviving ends, not source orientation.";
		assert(c.graphLeftLimited && c.graphRightLimited && c.leftEdges==null && c.rightEdges==null) :
				"Endpoint limits survive; graph edges must be rebuilt for the new contig identities.";
		assert(Arrays.equals(a.bases, originalA) && Arrays.equals(b.bases, originalB)) : "Spelling must not mutate the original unitigs.";
	}

	/** Closed mate components remain original unitigs, including their endpoint semantics. */
	private static void testCycle(){
		final Contig a=contig("AAAACCCC", 0), b=contig("CCCCAAAA", 1);
		a.leftCode=Tadpole.LOOP; a.rightCode=Tadpole.LOOP; a.leftRatio=13;
		b.leftCode=Tadpole.D_BRANCH; b.rightCode=Tadpole.B_BRANCH; b.rightRatio=17;
		final SameKJoinMerge.Result result=SameKJoinMerge.merge(new ArrayList<Contig>(Arrays.asList(a, b)), new int[]{3,2,1,0}, 5);
		assert(result.joins==0 && result.cycleUnitigs==2 && result.contigs.size()==2) : "Closed components must not be arbitrarily opened.";
		assert(result.contigs.get(0).leftCode==a.leftCode && result.contigs.get(0).leftRatio==13 &&
				result.contigs.get(1).rightCode==b.rightCode && result.contigs.get(1).rightRatio==17) : "Unmerged cycle metadata must remain original.";
	}

	/** Creates a small exact-sequence unitig without adding graph edges. */
	private static Contig contig(final String bases, final int id){
		assert(bases!=null && bases.length()>5) : "Fixture unitigs must exceed the K-1 overlap.";
		return new Contig(bases.getBytes(StandardCharsets.US_ASCII), id);
	}
}
