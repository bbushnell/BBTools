package test.prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import ddl.FamilySketchBuilder;
import prot.ProteinSequence;
import prot.cluster.ClusterSketchCodec;
import prot.cluster.ClusteringMetrics;
import prot.cluster.ConsensusRefiner;
import prot.cluster.Family;
import prot.cluster.GreedyClusterer;
import prot.cluster.TailAlignmentDecision;
import prot.cluster.TailAlignmentGate;

/** Fixed-seed retention and greedy-worker limits, using real gate/index/codec paths.
 * Controlled updates isolate recovery; native refiners cover empty parallel rebuilds.
 * These bounded fixtures make no biological-quality or performance claim.
 * @author Fischl
 */
public final class FinalPruneTest{

	/** Runs inside the canonical Dori driver, with a new output directory. */
	public static void run(final Path root) throws Exception{
		Files.createDirectory(root);
		smallFamilyRecovery(root);
		emptyReferenceRecovery(root);
		nativeRetention(root);
		invalidConfigurations(root);
		greedyWorkers(root);
		System.out.println("FINAL_PRUNE_TEST_PASS checks="+checks);
	}

	/** Four members survive until a valid bridge model recruits the fifth on pass two. */
	private static void smallFamilyRecovery(final Path root) throws Exception{
		final ProteinSequence[] input=sequences("AWCFDY","AWCFDY","AWCFDY","AWCFDY","CWDFAY");
		final ClusterSketchCodec.Loaded seeds=seeds(root.resolve("small.pddl"),new int[]{7},"AWCFDY");
		final byte[] bridge=new ProteinSequence("bridge","AWCFAY").enc;
		check(!gate().decide(input[4].enc,input[0].enc).accepted(),"Fifth input must initially miss the only seed");
		for(ProteinSequence s:input){check(gate().decide(s.enc,bridge).accepted(),"Replacement must admit every recovered member");}
		for(int threads:new int[]{1,4}){
			for(boolean deferred:new boolean[]{false,true}){
				final int[] updates={0};
				final GreedyClusterer.ConsensusUpdater updater=(family,seq,g,m,online)->{
					check(!online && !family.isEmpty(),"Fixed-seed rebuild may only update assigned families");
					if(updates[0]++==0){
						check(family.size()==4,"All four initially matching members must survive the first pass");
						family.setConsensus(bridge);
					}
				};
				final GreedyClusterer.Result result=run(root,"small_"+threads+"_"+deferred,input,seeds,2,threads,deferred,updater);
				check(result.assigned==(deferred ? 5 : 0) && result.families.length==(deferred ? 1 : 0),
					"Only deferred pruning can recover a four-member seed before the final five-member cutoff");
				check(updates[0]==(deferred ? 2 : 0),"Pruned references must not be rebuilt or queried later");
			}
		}
	}

	/** A weaker seed gets zero members, then wins a deterministic tie after a valid model update. */
	private static void emptyReferenceRecovery(final Path root) throws Exception{
		final ProteinSequence[] input=sequences("AWCFAY","AWCFAY","AWCFAY","AWCFAY","AWCFAY");
		final ClusterSketchCodec.Loaded seeds=seeds(root.resolve("empty_recovery.pddl"),new int[]{2,19},"CWDFAY","AWCFDY");
		final TailAlignmentDecision weak=gate().decide(input[0].enc,seeds.families[0].consensus);
		final TailAlignmentDecision strong=gate().decide(input[0].enc,seeds.families[1].consensus);
		check(weak.accepted() && strong.accepted() && strong.scoreRatio>weak.scoreRatio,
			"The later sparse-ID seed must strictly beat the first before the controlled update");
		for(int threads:new int[]{1,4}){
			final int[] updates={0};
			final GreedyClusterer.ConsensusUpdater updater=(family,seq,g,m,online)->{
				check(!family.isEmpty(),"A retained empty seed must keep its old model without entering the updater");
				check(family.familyId==(updates[0]==0 ? 19 : 2),"Pass two must recruit into the previously empty reference");
				if(updates[0]++==0){family.setConsensus(seeds.families[0].consensus);}
			};
			final GreedyClusterer.Result result=run(root,"empty_recovery_"+threads,input,seeds,2,threads,true,updater);
			check(result.assigned==5 && result.families.length==1 && result.families[0].familyId==2 && updates[0]==2,
				"Previously empty seed must recover all five members; the newly empty seed is removed only at the end");
			for(int i=0; i<input.length; i++){check(result.familyFor(i)==2,"Recovered membership retains stable sparse family ID");}
			phase(root.resolve("empty_recovery_"+threads+".checkpoints.tsv"),"consensus",1,2,5);
		}
	}

	/** Native serial/parallel paths retain zero-member models without rebuilding an empty family. */
	private static void nativeRetention(final Path root) throws Exception{
		final String a="AWCFDY",b="RRRRRR",c="WWWWWW";
		final ProteinSequence[] input=sequences(a,a,a,a,a,b,b,b,b);
		final ClusterSketchCodec.Loaded seeds=seeds(root.resolve("native.pddl"),new int[]{2,7,19},a,b,c);
		GreedyClusterer.Result serial=null;
		for(int threads:new int[]{1,4}){
			try(ConsensusRefiner updater=new ConsensusRefiner(1,false,0,null,false)){
				final String name="native_"+threads;
				final GreedyClusterer.Result result=run(root,name,input,seeds,3,threads,true,updater);
				check(result.assigned==5 && result.families.length==1 && result.families[0].familyId==2,
					"Final cutoff keeps exactly five members and removes four-member and empty references");
				for(int pass=1; pass<3; pass++){phase(root.resolve(name+".checkpoints.tsv"),"consensus",pass,3,9);}
				phase(root.resolve(name+".checkpoints.tsv"),"consensus",3,1,5);
				if(serial==null){serial=result;}else{same(serial,result);}
			}
		}
		for(int threads:new int[]{1,4}){
			try(ConsensusRefiner updater=new ConsensusRefiner(1,false,0,null,false)){
				final GreedyClusterer.Result regular=run(root,"one_off_"+threads,input,seeds,1,threads,false,updater);
				final GreedyClusterer.Result deferred=run(root,"one_on_"+threads,input,seeds,1,threads,true,updater);
				same(regular,deferred);
			}
		}
		final ClusterSketchCodec.Loaded noSeeds=seeds(root.resolve("no_seeds.pddl"),new int[0]);
		final ProteinSequence[][] populations={new ProteinSequence[0],sequences("ACD"),sequences(c,c)};
		for(int threads:new int[]{1,4}){
			for(int i=0; i<populations.length; i++){
				try(ConsensusRefiner updater=new ConsensusRefiner(1,false,0,null,false)){
					final GreedyClusterer.Result result=run(root,"no_survivors_"+threads+"_"+i,populations[i],seeds,2,threads,true,updater);
					check(result.assigned==0 && result.families.length==0,"Empty, filtered or subthreshold inputs must finish with no families");
				}
			}
			try(ConsensusRefiner updater=new ConsensusRefiner(1,false,0,null,false)){
				final GreedyClusterer.Result result=run(root,"no_models_"+threads,input,noSeeds,2,threads,true,updater);
				check(result.assigned==0 && result.families.length==0,"An empty seed codec must not create new families");
			}
		}
	}

	/** Configuration is checked at execution, so fluent setter order cannot bypass the scope. */
	private static void invalidConfigurations(final Path root) throws Exception{
		final ProteinSequence[] input=sequences("AWCFDY");
		final ClusterSketchCodec.Loaded seeds=seeds(root.resolve("guards.pddl"),new int[]{7},"AWCFDY");
		for(int scenario=0; scenario<4; scenario++){
			try(ClusteringMetrics metrics=metrics(null); ConsensusRefiner updater=new ConsensusRefiner(1,false,0,null,false)){
				final GreedyClusterer engine=engine(input,scenario==1 ? 0 : 1,4,updater,metrics);
				if(scenario==2){engine.fixedSeeds(seeds).merge(true).pruneFinal(true);}
				else if(scenario==3){engine.pruneFinal(true).fixedSeeds(seeds).merge(true);}
				else{engine.pruneFinal(true);}
				reject(()->engine.run(),"prunefinal requires");
			}
		}
		for(int n:new int[]{0,-1,5}){
			try(ClusteringMetrics metrics=metrics(null); ConsensusRefiner updater=new ConsensusRefiner(1,false,0,null,false)){
				reject(()->engine(input,1,4,updater,metrics).greedyThreads(n),"Greedy workers");
			}
		}
	}

	/** Configured and actual worker counts must both honor the cap or explicit override. */
	private static void greedyWorkers(final Path root) throws Exception{
		final String base="AWCFDY";
		final String[] raw=new String[400]; Arrays.fill(raw,base);
		final ProteinSequence[] input=sequences(raw);
		for(int override:new int[]{0,2,20}){
			try(ClusteringMetrics metrics=metrics(null); ConsensusRefiner updater=new ConsensusRefiner(1,false,0,null,false)){
				final GreedyClusterer engine=engine(input,1,32,updater,metrics).concurrentGreedy(.02);
				if(override>0){engine.greedyThreads(override);}
				final GreedyClusterer.Result result=engine.run();
				check(result.assigned==input.length && result.families.length==1,"Worker configuration must preserve identical-input ownership");
				final int expected=override==0 ? 16 : override;
				final String stats=String.join("\n",engine.concurrentStatistics());
				check(stats.contains("greedy_threads\t"+expected+"\n") && stats.contains("greedy_workers_started\t"+expected+"\n"),
					"Actual asynchronous workers must follow default16 or explicit1..t, independently of deterministic t32");
				reject(()->engine.greedyThreads(1),"Greedy workers");
			}
		}
	}

	/** Uses a separate model file: seed counts are source metadata, never fabricated input assignments. */
	private static ClusterSketchCodec.Loaded seeds(final Path path,final int[] ids,final String...raw){
		check(ids.length==raw.length,"Seed ID/model arrays must have identical cardinality");
		final ProteinSequence[] input=sequences(raw);
		final Family[] families=new Family[input.length];
		for(int i=0; i<input.length; i++){
			families[i]=new Family(ids[i],"seed_"+ids[i],i,input[i].enc);
			families[i].refreshStats(input,sketches());
		}
		ClusterSketchCodec.write(path.toString(),false,families,12,64,7,false,0,0);
		return ClusterSketchCodec.read(path.toString());
	}

	/** Runs with independent metrics and records intermediate family populations. */
	private static GreedyClusterer.Result run(final Path root,final String name,final ProteinSequence[] input,
			final ClusterSketchCodec.Loaded seeds,final int passes,final int threads,final boolean deferred,
			final GreedyClusterer.ConsensusUpdater updater){
		try(ClusteringMetrics metrics=metrics(root.resolve(name+".checkpoints.tsv"))){
			return engine(input,passes,threads,updater,metrics).pruneFinal(deferred).fixedSeeds(seeds).run();
		}
	}

	/** All controlled six-residue candidates have empty K12 sketches, retaining the actual gate. */
	private static GreedyClusterer engine(final ProteinSequence[] input,final int passes,final int threads,
			final GreedyClusterer.ConsensusUpdater updater,final ClusteringMetrics metrics){
		return new GreedyClusterer(input,6,8,5,passes,0,17,sketches(),gate(),updater,metrics,true,threads);
	}
	private static FamilySketchBuilder sketches(){return new FamilySketchBuilder(64,12,7,false,0,0,null);}
	private static TailAlignmentGate gate(){return new TailAlignmentGate(TailAlignmentGate.Mode.PIN_RIGHT_FULL_QUERY,
		TailAlignmentGate.Denominator.FULL_QUERY,6,6,.6,.1);}
	private static ClusteringMetrics metrics(final Path path){return new ClusteringMetrics(path==null ? null : path.toString(),null,false,8,Long.MAX_VALUE);}

	/** Generates unique IDs while retaining caller-specified sequence multiplicities. */
	private static ProteinSequence[] sequences(final String...raw){
		final ProteinSequence[] out=new ProteinSequence[raw.length];
		for(int i=0; i<raw.length; i++){out[i]=new ProteinSequence("q"+i,raw[i]);}
		return out;
	}

	/** Checks the terminal checkpoint for one phase, including retained empty model slots. */
	private static void phase(final Path path,final String stage,final int pass,final int families,final int assigned) throws Exception{
		String[] last=null;
		for(String line:Files.readAllLines(path,StandardCharsets.US_ASCII)){
			final String[] cells=line.split("\t");
			if(cells[0].equals(stage) && cells[1].equals(Integer.toString(pass))){last=cells;}
		}
		check(last!=null && Integer.parseInt(last[6])==families && Integer.parseInt(last[4])==assigned,
			"Unexpected terminal population for "+stage+":"+pass+" in "+path);
	}

	/** Compares final assignments and exported model payloads across worker counts or one-pass modes. */
	private static void same(final GreedyClusterer.Result a,final GreedyClusterer.Result b){
		check(a.assigned==b.assigned && a.inputs()==b.inputs() && a.families.length==b.families.length,"Final populations must agree");
		for(int i=0; i<a.inputs(); i++){check(a.familyFor(i)==b.familyFor(i),"Every final assignment must agree");}
		for(int i=0; i<a.families.length; i++){
			check(a.families[i].familyId==b.families[i].familyId && Arrays.equals(a.families[i].consensus(),b.families[i].consensus()) &&
				Arrays.equals(a.families[i].sketch(),b.families[i].sketch()),"Final model identity, sequence and sketch must agree");
		}
	}

	/** Requires the intended configuration failure rather than accepting any unrelated exception. */
	private static void reject(final Runnable action,final String message){
		try{action.run();}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(message),"Wrong rejection: "+expected); return;
		}
		throw new AssertionError("Expected rejection: "+message);
	}
	private static void check(final boolean ok,final String why){checks++; if(!ok){throw new AssertionError(why);}}
	private static int checks;
}
