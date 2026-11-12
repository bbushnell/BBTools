package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import idaligner.AlignmentStats;
import idaligner.ContigEdgeAlignment;
import idaligner.EndClippedAligner;
import parse.PreParser;
import prok.NcrnaAlignmentCache;
import prok.NcrnaAlignmentWork;
import prok.NcrnaCpuProfile;
import prok.TrnaBoundaryFeatures;
import prok.TrnaConsensusBuilder;
import shared.Shared;
import stream.Read;

/** Actual kernel census fixtures: cache hits, hidden HBM/boundary traces,
 * no-op edge refinement and independent workers. Counts and operand lengths
 * have literal expectations; elapsed CPU is not a throughput benchmark.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaAlignmentWorkTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaAlignmentWorkTest.class, false);
		try{
			check(pp.args.length==1 && pp.args[0].startsWith("out="), "Supply one fresh out=profile prefix");
			final String prefix=pp.args[0].substring(4);
			System.setProperty("bbtools.diagnostic.ncrnaCpuOut", prefix);
			run(prefix);
		}finally{Shared.closeStream(pp.outstream);}
		System.out.println("NCRNA_ALIGNMENT_WORK_TEST_PASS");
	}

	private static void run(final String prefix){
		final byte[] model=bytes("ACGTGCAATGCTAGCATGACCTGATCGTAGCTAGTCGATGCATGACTGACGTACGATCGTAGCTAGCTGACGATCGTACGATGCATCGATGCTAGCTA");
		final byte[] candidate=new byte[model.length+20];
		Arrays.fill(candidate, (byte)'A');
		System.arraycopy(model, 0, candidate, 10, model.length);
		final ArrayList<Read> members=new ArrayList<Read>();
		members.add(new Read(model.clone(), null, "member", 0));
		final BaseGraph graph=TrnaConsensusBuilder.buildGraphFromReference(members, model, "fixture");
		final NcrnaCpuProfile profile=NcrnaCpuProfile.create(1);
		check(profile!=null, "The property must initialize the real CPU clock before profile creation");
		profile.begin(0, NcrnaCpuProfile.SCAVENGE);
		final NcrnaAlignmentWork work=profile.alignmentWork(0);
		try{
			final NcrnaAlignmentCache cache=new NcrnaAlignmentCache(new byte[][]{model}, 2, 512);
			final int[] positions=new int[4];
			cache.quantum(0, candidate, positions);
			Arrays.fill(positions, 0);
			cache.quantum(0, candidate, positions);
			check(work.logical==2 && work.physical==1 && work.hits==1,
				"Two identical requested kernels execute once and retain one logical cache hit");
			NcrnaAlignmentWork.scrabbleTrace(model, candidate, new AlignmentStats(true), model.length);
			TrnaConsensusBuilder.scoreAgainstModel(model, graph);
			TrnaBoundaryFeatures.alignToModelFrame(model, graph, new AlignmentStats(true));
			check(work.logical==5 && work.physical==4 && work.hbmCalls==1,
				"Trim, HBM traceback and boundary traceback each add one kernel; HBM likelihood is counted separately");

			final byte[] partial=Arrays.copyOfRange(model, 16, model.length);
			final EndClippedAligner.Result edge=new EndClippedAligner.Result();
			edge.qLen=model.length;
			edge.qStart=16;
			edge.qStop=model.length-1;
			edge.rLen=partial.length;
			edge.rStart=0;
			edge.rStop=partial.length-1;
			edge.physicalClips=3;
			ContigEdgeAlignment.refineInterior(model, partial, edge, 0, 0, NcrnaAlignmentWork.refinementObserver());
			check(work.logical==5, "Two fixed physical ends return without issuing a kernel");
			edge.physicalClips=1;
			ContigEdgeAlignment.refineInterior(model, partial, edge, 0, 0, NcrnaAlignmentWork.refinementObserver());
			check(work.logical==6 && work.physical==5 && work.hits==1,
				"One interior-facing end issues exactly one observed-slice alignment");
			check(work.modelBases==5L*model.length+partial.length
				&& work.physicalModelBases==4L*model.length+partial.length,
				"Model bases count each actual model operand, including the shorter observed slice and logical hit");
			check(work.cpuNanos>=0 && work.hbmCpuNanos>=0, "Enabled worker CPU measurements must be nonnegative");
		}finally{profile.end();}
		check(NcrnaAlignmentWork.current()==null, "Completed family scopes release their thread-local counter");
		profile.write(prefix, new String[]{"fixture"});
		workers(model);
	}

	private static void workers(final byte[] model){
		final NcrnaCpuProfile[] profiles=new NcrnaCpuProfile[2];
		final Throwable[] errors=new Throwable[2];
		final Thread[] workers=new Thread[2];
		for(int i=0; i<workers.length; i++){
			final int id=i;
			workers[i]=new Thread(() -> {
				try{
					final NcrnaCpuProfile profile=NcrnaCpuProfile.create(1);
					profiles[id]=profile;
					profile.begin(0, NcrnaCpuProfile.SCAVENGE);
					try{
						for(int j=0; j<=id; j++){
							NcrnaAlignmentWork.scrabble(model, model, null, model.length);
						}
					}finally{profile.end();}
					check(NcrnaAlignmentWork.current()==null, "Worker scope must not leak");
				}catch(Throwable failure){errors[id]=failure;}
			});
			workers[i].start();
		}
		for(final Thread worker:workers){
			try{worker.join();}
			catch(InterruptedException failure){
				Thread.currentThread().interrupt();
				throw new IllegalStateException(failure);
			}
		}
		for(final Throwable error:errors){if(error!=null){throw new AssertionError("Profile worker failed", error);}}
		check(profiles[0].alignmentWork(0).logical==1 && profiles[1].alignmentWork(0).logical==2,
			"Concurrent workers retain separate one-request and two-request counters");
		final NcrnaCpuProfile total=NcrnaCpuProfile.create(1);
		total.merge(profiles[0]);
		total.merge(profiles[1]);
		check(total.alignmentWork(0).logical==3 && total.alignmentWork(0).modelBases==3L*model.length,
			"Joined workers sum exactly once without per-request atomics or shared increments");
	}

	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean condition, final String reason){if(!condition){throw new AssertionError(reason);}}
}
