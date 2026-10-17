package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import dna.AminoAcid;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;

/** Exact spelling, orientation and no-mutation tests for the optional graph veto. @author Fischl */
public final class GraphOverlapGateTest {

	/** Runs synthetic cases and optionally the preserved three-record native proof fixture. */
	public static void main(final String[] args){
		boolean enabled=false;
		assert(enabled=true);
		if(!enabled){throw new IllegalStateException("Tests require -ea.");}
		final PreParser pp=new PreParser(args, GraphOverlapGateTest.class, false);
		try{
			final byte[] reference=randomBases(1100, 1234);
			final byte[] left=Arrays.copyOfRange(reference, 0, 500);
			final byte[] right=Arrays.copyOfRange(reference, 405, 905);
			final byte[] edge=Arrays.copyOfRange(right, 95, 268);
			testPair(left, right, edge);
			testBoundaries(reference);
			testConfig();
			if(pp.args.length>0){
				if(pp.args.length!=1 || !pp.args[0].startsWith("in=")){
					throw new IllegalArgumentException("Optional in=parents.fa contains left, right and edge records.");
				}
				final ArrayList<Read> reads=ReadInputStream.toReads(pp.args[0].substring(3), FileFormat.FA, -1);
				assert(reads.size()==3 && reads.get(0).id.equals("left") && reads.get(1).id.equals("right")
						&& reads.get(2).id.equals("edge")) : "Native proof fixture must preserve its three named inputs.";
				assert(reads.get(0).length()==500 && reads.get(1).length()==500 && reads.get(2).length()==173) :
						"The proved duplication uses 500bp parent windows and a 173bp edge.";
				testPair(reads.get(0).bases, reads.get(1).bases, reads.get(2).bases);
			}
			System.out.println("GRAPH_OVERLAP_GATE_PASS native="+(pp.args.length>0));
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** The veto must act in every stored orientation; disabled and cross-K paths remain unchanged. */
	private static void testPair(final byte[] left, final byte[] right, final byte[] edge){
		assert(Arrays.equals(Arrays.copyOfRange(left, 405, 500), Arrays.copyOf(right, 95))) :
				"This regression reproduces an existing 95bp parent overlap.";
		for(int mask=0; mask<4; mask++){
			testGraph(left, right, edge, mask, true, false, false);
			testGraph(left, right, edge, mask, false, false, true);
			testGraph(left, right, edge, mask, true, true, true);
		}
	}

	/** Exercises the actual BubblePopper call site with reciprocal edges and optional flips. */
	private static void testGraph(final byte[] left, final byte[] right, final byte[] path,
			final int mask, final boolean enabled, final boolean crossK, final boolean expectMerge){
		final Contig a=new Contig(left.clone(), 0), b=new Contig(right.clone(), 1);
		a.coverage=b.coverage=10;
		a.leftCode=b.rightCode=Tadpole.DEAD_END;
		a.rightCode=b.leftCode=Tadpole.B_BRANCH;
		final Edge forward=new Edge(0, 1, path.length, 1, 10, path.clone());
		final Edge backward=new Edge(1, 0, path.length, 2, 10, reverse(path));
		a.addRightEdge(forward); b.addLeftEdge(backward);
		final HashMap<Integer, ArrayList<Edge>> incoming=new HashMap<Integer, ArrayList<Edge>>();
		incoming.put(0, new ArrayList<Edge>(Arrays.asList(backward)));
		incoming.put(1, new ArrayList<Edge>(Arrays.asList(forward)));
		if((mask&1)!=0){a.flip(incoming.get(0));}
		if((mask&2)!=0){b.flip(incoming.get(1));}
		final byte[] originalB=b.bases.clone();
		final GraphOverlapGate oldGate=BubblePopper.directOverlapGate;
		final BubblePopper.CoverageGate oldCoverage=BubblePopper.directCoverageGate;
		final boolean oldDirect=BubblePopper.popDirect, oldIndirect=BubblePopper.popIndirect;
		final boolean oldUnzip=BubblePopper.unzipBubbles, oldCross=BubblePopper.crossKMerge;
		final boolean oldValidate=BubblePopper.validateGraph;
		final GraphOverlapGate gate=enabled ? new GraphOverlapGate(95) : null;
		try{
			BubblePopper.directOverlapGate=gate; BubblePopper.directCoverageGate=null;
			BubblePopper.popDirect=true; BubblePopper.popIndirect=false;
			BubblePopper.unzipBubbles=false; BubblePopper.crossKMerge=crossK;
			BubblePopper.validateGraph=true;
			final BubblePopper popper=new BubblePopper(new ArrayList<Contig>(Arrays.asList(a, b)), incoming, 64);
			final int merges=popper.expand(a);
			assert(merges==(expectMerge ? 1 : 0) && b.used()==expectMerge) :
					"Only the enabled ordinary-graph policy may veto this walk; mask="+mask+", crossK="+crossK;
			if(expectMerge){
				assert(a.length()==left.length+path.length+right.length-64) : "Disabled/cross-K merge spelling must remain unchanged.";
			}else{
				assert(Arrays.equals(originalB, b.bases) && a.length()==left.length && gate.rejected==1) :
						"Rejected destinations must not be flipped, consumed or truncated.";
				assert(popper.validate(a) && popper.validate(b)) : "Vetoes must retain valid reciprocal graph connectivity.";
			}
		}finally{
			BubblePopper.directOverlapGate=oldGate; BubblePopper.directCoverageGate=oldCoverage;
			BubblePopper.popDirect=oldDirect; BubblePopper.popIndirect=oldIndirect;
			BubblePopper.unzipBubbles=oldUnzip; BubblePopper.crossKMerge=oldCross;
			BubblePopper.validateGraph=oldValidate;
		}
	}

	/** Preserves true gaps, ordinary K-1 overlaps and matching spacing at the lower boundary. */
	private static void testBoundaries(final byte[] reference){
		final byte[] left=Arrays.copyOf(reference, 500);
		final byte[] gapRight=Arrays.copyOfRange(reference, 520, 1020);
		final byte[] gapPath=Arrays.copyOfRange(reference, 500, 584);
		for(int mask=0; mask<4; mask++){testGraph(left, gapRight, gapPath, mask, true, false, true);}
		final GraphOverlapGate gate=new GraphOverlapGate(95);
		for(int overlap : new int[]{63,64,95,96}){
			final Contig a=new Contig(left, 0);
			final Contig b=new Contig(Arrays.copyOfRange(reference, 500-overlap, 1000-overlap), 1);
			final Edge e=new Edge(0, 1, 1, 1, 10, new byte[]{'A'});
			assert(gate.accept(a, b, e, 64)==(overlap<64 || overlap>95)) : "Exact overlap must obey both search bounds: "+overlap;
			if(overlap==64){
				e.bases=null; e.length=0;
				assert(gate.accept(a, b, e, 64)) : "An identical implied displacement is not a conflict.";
			}
			if(overlap==95){
				b.bases[40]='N';
				assert(gate.accept(a, b, e, 64)) : "Undefined bases cannot establish an exact overlap veto.";
			}
		}
		final Contig a=new Contig(left, 0), b=new Contig(Arrays.copyOfRange(reference, 405, 905), 1);
		final Edge exact=new Edge(0, 1, 0, 1, 10, null, 95);
		assert(gate.accept(a, b, exact, 64)) : "Explicit fusion overlaps must retain their own policy.";
		assert(gate.accept(a, b, new Edge(0, 1, 1, 1, 10, new byte[]{'A'}), 96)) :
				"Initial assemble-K cannot match an overlap shorter than itself.";
	}

	/** The experimental flag is consumed by the driver, not leaked into per-K native arguments. */
	private static void testConfig(){
		final String[] args={"k=96,124,64", "assemblek=96", "out=unused.fa", "graphmergecoverageratio=1.75", "graphoverlapguard=t"};
		new MultiKBranchRescue(args);
		boolean refused=false;
		try{
			new MultiKBranchRescue(new String[]{"k=96,64", "out=unused.fa", "sweeplen=0", "graphoverlapguard=t"});
		}catch(IllegalArgumentException expected){refused=true;}
		assert(refused) : "The veto must not silently create an unrequested final table load.";
	}

	/** Generates repeat-free-in-practice fixture context with a deterministic seed. */
	private static byte[] randomBases(final int length, final long seed){
		assert(length>0) : "Fixtures require nonempty parental sequence.";
		final Random random=new Random(seed);
		final byte[] bases=new byte[length], alphabet={'A','C','G','T'};
		for(int i=0; i<length; i++){bases[i]=alphabet[random.nextInt(4)];}
		return bases;
	}

	/** Returns a separate reverse complement for stored-orientation fixtures. */
	private static byte[] reverse(final byte[] bases){
		assert(bases!=null) : "Orientation fixtures require an explicit sequence payload.";
		final byte[] reversed=new byte[bases.length];
		for(int i=0; i<bases.length; i++){reversed[i]=AminoAcid.baseToComplementExtended[bases[bases.length-1-i]];}
		return reversed;
	}
}
