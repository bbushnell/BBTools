package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

/** Boundary, orientation, call-site and restoration tests for the final branch veto. @author Fischl */
public final class GraphBranchGateTest {

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs synthetic checks; optional ordinary driver arguments also exercise final-phase restoration. */
	public static void main(final String[] args){
		boolean enabled=false;
		assert(enabled=true);
		if(!enabled){throw new IllegalStateException("Tests require -ea.");}
		final byte[] reference=randomBases(1000);
		final int[] codes={Tadpole.KEEP_GOING, Tadpole.B_BRANCH, Tadpole.F_BRANCH, Tadpole.D_BRANCH};
		int cases=0;
		for(int left : codes){
			for(int right : codes){
				for(int length : new int[]{0, 1, 19, 31, 32, 33, 165}){
					for(int mask=0; mask<4; mask++){
						for(int mode=0; mode<3; mode++){
							testGraph(reference, left, right, length, 0, mask, mode);
							cases++;
						}
					}
				}
			}
		}
		for(int overlap : new int[]{13, 32, 95}){
			for(int mask=0; mask<4; mask++){
				for(int mode=0; mode<3; mode++){
					testGraph(reference, Tadpole.F_BRANCH, Tadpole.D_BRANCH, 0, overlap, mask, mode);
					cases++;
				}
			}
		}
		testConfig();
		// The standalone suite has no input requirement; the delegated driver owns
		// standard config/output parsing only when actual assembly args are supplied.
		if(args.length>0){testRestoration(args);}
		System.out.println("GRAPH_BRANCH_GATE_PASS cases="+cases+" lifecycle="+(args.length>0));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Static Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Exercises real direct merging; mode0 is disabled, mode1 guarded, mode2 guarded cross-K. */
	private static void testGraph(final byte[] reference, final int leftCode, final int rightCode,
			final int pathLength, final int overlap, final int mask, final int mode){
		final int k=32, parentLength=200;
		assert(pathLength>=0 && overlap>=0 && mask>=0 && mask<4 && mode>=0 && mode<3) :
				"Fixture geometry and mode must describe a supported BubblePopper direct join.";
		final int rightStart=parentLength+(overlap>0 ? -overlap : pathLength-k);
		final Contig a=new Contig(Arrays.copyOf(reference, parentLength), 0);
		final Contig b=new Contig(Arrays.copyOfRange(reference, rightStart, rightStart+parentLength), 1);
		a.coverage=b.coverage=10;
		a.leftCode=b.rightCode=Tadpole.DEAD_END;
		a.rightCode=leftCode; b.leftCode=rightCode;
		final byte[] forward=pathLength==0 ? null : Arrays.copyOfRange(reference, parentLength, parentLength+pathLength);
		// Left-exit payloads are stored in contig strand, not outward reverse-complement order.
		final byte[] backward=pathLength==0 ? null : Arrays.copyOfRange(reference, parentLength-k, parentLength-k+pathLength);
		final Edge out=new Edge(0, 1, pathLength, 1, 10, forward, overlap);
		final Edge back=new Edge(1, 0, pathLength, 2, 10, backward, overlap);
		a.addRightEdge(out); b.addLeftEdge(back);
		final HashMap<Integer, ArrayList<Edge>> incoming=new HashMap<Integer, ArrayList<Edge>>();
		incoming.put(0, new ArrayList<Edge>(Arrays.asList(back)));
		incoming.put(1, new ArrayList<Edge>(Arrays.asList(out)));
		if((mask&1)!=0){a.flip(incoming.get(0));}
		if((mask&2)!=0){b.flip(incoming.get(1));}
		final byte[] originalB=b.bases.clone();
		final boolean branch=leftCode==Tadpole.F_BRANCH || leftCode==Tadpole.D_BRANCH
				|| rightCode==Tadpole.F_BRANCH || rightCode==Tadpole.D_BRANCH;
		final boolean expectMerge=mode!=1 || overlap>0 || !branch;
		final GraphBranchGate oldGate=BubblePopper.directBranchGate;
		final GraphOverlapGate oldOverlap=BubblePopper.directOverlapGate;
		final BubblePopper.CoverageGate oldCoverage=BubblePopper.directCoverageGate;
		final FusionKmerSupport oldSupport=BubblePopper.crossKSupport;
		final boolean oldDirect=BubblePopper.popDirect, oldIndirect=BubblePopper.popIndirect;
		final boolean oldUnzip=BubblePopper.unzipBubbles, oldCross=BubblePopper.crossKMerge;
		final boolean oldValidate=BubblePopper.validateGraph;
		final int oldMismatches=BubblePopper.crossKMaxMismatches;
		final GraphBranchGate gate=mode==0 ? null : new GraphBranchGate();
		try{
			BubblePopper.directBranchGate=gate;
			BubblePopper.directOverlapGate=null; BubblePopper.directCoverageGate=null;
			BubblePopper.crossKSupport=null; BubblePopper.crossKMaxMismatches=-1;
			BubblePopper.popDirect=true; BubblePopper.popIndirect=false;
			BubblePopper.unzipBubbles=false; BubblePopper.crossKMerge=mode==2;
			BubblePopper.validateGraph=true;
			final BubblePopper popper=new BubblePopper(new ArrayList<Contig>(Arrays.asList(a, b)), incoming, k);
			final int merges=popper.expand(a);
			assert(merges==(expectMerge ? 1 : 0) && b.used()==expectMerge) :
					"Only guarded ordinary implicit paths at outward branches may be vetoed: codes="+
					leftCode+","+rightCode+", path="+pathLength+", mask="+mask+", mode="+mode;
			if(expectMerge){
				assert(Arrays.equals(a.bases, Arrays.copyOf(reference, rightStart+parentLength))) :
						"Allowed joins must preserve exact native spelling, including zero-gap and overlap boundaries.";
			}else{
				assert(a.length()==parentLength && Arrays.equals(b.bases, originalB) && !b.associate()) :
						"A veto must not consume, trim or flip its destination.";
				assert(popper.validate(a) && popper.validate(b)) : "Vetoed joins retain their reciprocal graph edges.";
			}
			if(gate!=null){
				assert(gate.evaluations==(mode==1 && overlap==0 ? 1 : 0)
						&& gate.rejected==(expectMerge ? 0 : 1)) :
						"Counters exclude explicit overlaps/cross-K paths but include implicit zero/short paths exactly once.";
			}
		}finally{
			BubblePopper.directBranchGate=oldGate;
			BubblePopper.directOverlapGate=oldOverlap; BubblePopper.directCoverageGate=oldCoverage;
			BubblePopper.crossKSupport=oldSupport; BubblePopper.crossKMaxMismatches=oldMismatches;
			BubblePopper.popDirect=oldDirect; BubblePopper.popIndirect=oldIndirect;
			BubblePopper.unzipBubbles=oldUnzip; BubblePopper.crossKMerge=oldCross;
			BubblePopper.validateGraph=oldValidate;
		}
	}

	/** The development flag must neither alter public defaults nor create an unrequested table phase. */
	private static void testConfig(){
		final String[] base={"k=96,124,64", "assemblek=96", "out=unused.fa", "graphmergecoverageratio=1.75"};
		assert(!new TadpoleMulti.Config(base).graphMergeBranchGuard) : "Public/default configuration must leave this experiment off.";
		final ArrayList<String> enabled=new ArrayList<String>(Arrays.asList(base));
		enabled.add("graphbranchguard=t");
		new MultiKBranchRescue(enabled.toArray(new String[0]));
		boolean refused=false;
		try{
			new MultiKBranchRescue(new String[]{"k=96,64", "out=unused.fa", "sweeplen=0", "graphbranchguard=t"});
		}catch(IllegalArgumentException expected){refused=true;}
		assert(refused) : "The branch guard cannot silently add an unrequested final graph table.";
	}

	/** Runs a harmless assembly and verifies restoration of an already installed caller policy. */
	private static void testRestoration(final String[] args){
		assert(args.length>0) : "Lifecycle test requires a complete native fixture invocation.";
		final GraphBranchGate old=BubblePopper.directBranchGate, sentinel=new GraphBranchGate();
		try{
			BubblePopper.directBranchGate=sentinel;
			MultiKBranchRescue.main(args);
			assert(BubblePopper.directBranchGate==sentinel) :
					"TadpoleMulti.extractFinalGraph must restore its caller's branch policy after final simplification.";
		}finally{BubblePopper.directBranchGate=old;}
	}

	/** Creates deterministic sequence for exact spelling tests. */
	private static byte[] randomBases(final int length){
		assert(length>0) : "Join fixtures need nonempty reference context.";
		final byte[] bases=new byte[length], alphabet={'A','C','G','T'};
		final Random random=new Random(85731);
		for(int i=0; i<length; i++){bases[i]=alphabet[random.nextInt(4)];}
		return bases;
	}
}
