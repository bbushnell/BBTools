package prot.cluster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;

import ddl.FamilySketchBuilder;
import ddl.MutableFamilyIndex;
import fileIO.ByteStreamWriter;
import prot.HbmPositionModel;
import structures.ByteBuilder;
import template.Accumulator;
import template.ThreadWaiter;

/** Scores each retained final consensus against the best accepted families in its retained-family shortlist. */
public final class ConsensusNeighborQualityReport{

	private ConsensusNeighborQualityReport(){}

	/** Computes neighbor quality and writes the deterministic paths returned by outputPaths(prefix). */
	public static void write(Family[] retained,ConsensusQualityReport cohesion,
			FamilySketchBuilder sketchBuilder,TailAlignmentGate gate,int threads,String prefix,
			boolean overwrite){
		final String[] paths=outputPaths(prefix);
		for(int i=0; i<paths.length; i++){
			for(int j=i+1; j<paths.length; j++){
				if(paths[i].equals(paths[j])){
					throw new IllegalArgumentException("Neighbor quality outputs alias each other");
				}
			}
		}
		final long start=System.nanoTime();
		final Prepared prepared=new Prepared(retained,cohesion,sketchBuilder,gate,threads);
		prepared.run();
		final long elapsed=System.nanoTime()-start;
		writeFamilies(paths[0],overwrite,prepared.families,prepared.results,cohesion);
		writeCohorts(paths[1],overwrite,prepared.families,prepared.results,prepared.resultByFamilyId);
		writeMethod(paths[2],overwrite,prepared,elapsed);
	}

	/** Output paths are deterministic so the driver can preflight them beside other fixed outputs. */
	public static String[] outputPaths(String prefix){
		if(prefix==null || prefix.length()<1){throw new IllegalArgumentException("Neighbor quality output prefix is required");}
		return new String[]{prefix+".neighbor_quality.tsv",prefix+".neighbor_cohorts.tsv",
			prefix+".neighbor_method.tsv"};
	}

	private static void writeFamilies(String path,boolean overwrite,Family[] families,
			Result[] results,ConsensusQualityReport cohesion){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("family_id\tmembers\tcohesion_available\tcohesion_identity\tcohesion_ratio");
		bb.append("\tmax_identity_available\tmax_identity\tmax_identity_neighbor_id\tmax_identity_neighbor_ratio");
		bb.append("\tmax_ratio_available\tmax_ratio\tmax_ratio_neighbor_id\tmax_ratio_neighbor_identity");
		bb.append("\traw_delta_available\traw_delta\tnormalized_delta_available\tnormalized_delta");
		bb.append("\tnormalized_status\tmissing_reason\tattempted_pairs\tpositive_span_pairs\n");
		for(int i=0; i<families.length; i++){
			final Family family=families[i];
			final ConsensusQualityReport.Snapshot a=cohesion.snapshot(family.familyId);
			final Result r=results[i];
			bb.append(family.familyId).tab().append(family.size()).tab().append(a.available() ? 1 : 0).tab();
			value(bb,a.available(),a.identityMean()).tab();
			value(bb,a.available(),a.scoreRatioMean()).tab();
			bb.append(r.hasIdentity ? 1 : 0).tab();
			value(bb,r.hasIdentity,r.identity).tab();
			value(bb,r.hasIdentity,r.identityNeighbor).tab();
			value(bb,r.hasIdentity,r.identityRatio).tab();
			bb.append(r.hasRatio ? 1 : 0).tab();
			value(bb,r.hasRatio,r.ratio).tab();
			value(bb,r.hasRatio,r.ratioNeighbor).tab();
			value(bb,r.hasRatio,r.ratioIdentity).tab();
			final boolean raw=a.available() && r.hasIdentity;
			final double delta=raw ? a.identityMean()-r.identity : 0;
			bb.append(raw ? 1 : 0).tab();
			value(bb,raw,delta).tab();
			final double denominator=raw ? Math.max(a.identityMean(),r.identity) : 0;
			final boolean normalized=raw && denominator>0;
			bb.append(normalized ? 1 : 0).tab();
			value(bb,normalized,normalized ? delta/denominator : 0).tab();
			bb.append(raw ? (normalized ? "ok" : "zero_denominator") : "no_raw_delta").tab();
			bb.append(a.available() ? (r.hasIdentity ? "ok" : r.missingReason) : a.reason()).tab();
			bb.append(r.attemptedPairs).tab().append(r.positiveSpanPairs).nl();
		}
		write(path,overwrite,bb);
	}

	private static void writeMethod(String path,boolean overwrite,Prepared prepared,long elapsed){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("mode\tdenominator\tmax_length\tmin_span\tneighbor_count\tthreads\tretained");
		bb.append("\tattempted_pairs\tpositive_span_pairs\telapsed_ns\telapsed_seconds\n");
		bb.append(prepared.gate.mode.toString()).tab().append(prepared.gate.denominator.toString()).tab();
		bb.append(prepared.gate.maxLength).tab().append(prepared.gate.minSpan).tab().append(10).tab();
		bb.append(prepared.workerCount).tab().append(prepared.families.length).tab();
		bb.append(prepared.attemptedPairs).tab().append(prepared.positiveSpanPairs).tab();
		bb.append(elapsed).tab().appendSlow(elapsed*1e-9).nl();
		write(path,overwrite,bb);
	}

	private static void writeCohorts(String path,boolean overwrite,Family[] families,
			Result[] results,int[] resultByFamilyId){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("#delta_percentiles\tnearest_rank\tceil(p*n)-1\n");
		bb.append("#normalized_delta\t(a-b)/max(a,b)\tzero denominator excluded from normalized summaries\n");
		bb.append("cohort\tretained\traw_delta_count\tnormalized_delta_count\tcohesion_missing\tneighbor_missing\tzero_denominator");
		bb.append("\traw_p10\traw_p50\traw_p90\traw_fraction_negative");
		bb.append("\tnormalized_p10\tnormalized_p50\tnormalized_p90\tnormalized_fraction_negative\n");
		cohort(bb,"all",families,results,resultByFamilyId);
		final Family[] top=families.clone();
		Arrays.sort(top,SIZE_ORDER);
		cohort(bb,"top_80000_by_size",top,Math.min(80000,top.length),results,resultByFamilyId);
		int large=0;
		for(Family family : top){if(family.size()>=100){large++;}}
		cohort(bb,"size_at_least_100",top,large,results,resultByFamilyId);
		write(path,overwrite,bb);
	}

	private static void cohort(ByteBuilder bb,String name,Family[] families,Result[] results,
			int[] resultByFamilyId){
		cohort(bb,name,families,families.length,results,resultByFamilyId);
	}

	private static void cohort(ByteBuilder bb,String name,Family[] families,int count,
			Result[] results,int[] resultByFamilyId){
		final double[] raw=new double[count], normalized=new double[count];
		int rawCount=0, normalizedCount=0, missingA=0, missingB=0, zero=0;
		for(int i=0; i<count; i++){
			final Result r=results[resultByFamilyId[families[i].familyId]];
			if(!r.hasA){missingA++;}
			if(!r.hasIdentity){missingB++;}
			if(r.hasA && r.hasIdentity){
				raw[rawCount++]=r.rawDelta;
				if(r.normalizedAvailable){normalized[normalizedCount++]=r.normalizedDelta;}
				else{zero++;}
			}
		}
		Arrays.sort(raw,0,rawCount); Arrays.sort(normalized,0,normalizedCount);
		bb.append(name).tab().append(count).tab().append(rawCount).tab().append(normalizedCount).tab()
			.append(missingA).tab().append(missingB).tab().append(zero).tab();
		metrics(bb,raw,rawCount); bb.tab(); metrics(bb,normalized,normalizedCount); bb.nl();
	}

	private static void metrics(ByteBuilder bb,double[] sorted,int count){
		if(count<1){bb.append("NA\tNA\tNA\tNA"); return;}
		value(bb,true,percentile(sorted,count,0.10)).tab();
		value(bb,true,percentile(sorted,count,0.50)).tab();
		value(bb,true,percentile(sorted,count,0.90)).tab();
		int negative=0;
		for(int i=0; i<count; i++){if(sorted[i]<0){negative++;}}
		value(bb,true,negative/(double)count);
	}

	private static double percentile(double[] sorted,int count,double p){
		final int index=Math.max(0,Math.min(count-1,(int)Math.ceil(p*count)-1));
		return sorted[index];
	}

	private static ByteBuilder value(ByteBuilder bb,boolean available,double value){return available ? bb.appendSlow(value) : bb.append("NA");}
	private static ByteBuilder value(ByteBuilder bb,boolean available,int value){return available ? bb.append(value) : bb.append("NA");}

	private static void write(String path,boolean overwrite,ByteBuilder bb){
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,false);
		writer.start(); writer.print(bb);
		if(writer.poisonAndWait()){throw new IllegalStateException("Consensus neighbor quality output failed: "+path);}
	}

	private static final Comparator<Family> SIZE_ORDER=new Comparator<Family>(){
		@Override
		public int compare(Family a,Family b){
			final int size=Integer.compare(b.size(),a.size());
			return size!=0 ? size : Integer.compare(a.familyId,b.familyId);
		}
	};

	private static final class Prepared implements Accumulator<Worker>{
		Prepared(Family[] retained,ConsensusQualityReport cohesion,FamilySketchBuilder builder,
				TailAlignmentGate gate,int threads){
			if(retained==null || cohesion==null || builder==null || gate==null || threads<1){
				throw new IllegalArgumentException("Neighbor quality needs families, cohesion, sketcher, gate and positive threads");
			}
			families=retained.clone(); this.cohesion=cohesion; this.builder=builder; this.gate=gate;
			final int width=builder.width();
			int max=-1;
			for(Family family : families){
				if(family==null || family.size()<1 || family.familyId<0){
					throw new IllegalArgumentException("Every retained family must be nonempty with a nonnegative stable ID");
				}
				max=Math.max(max,family.familyId);
			}
			byId=new Family[max+1]; resultByFamilyId=new int[max+1]; Arrays.fill(resultByFamilyId,-1);
			final MutableFamilyIndex mutable=new MutableFamilyIndex(width,true);
			for(int i=0; i<families.length; i++){
				final Family family=families[i];
				if(byId[family.familyId]!=null){
					throw new IllegalArgumentException("Duplicate retained family ID: "+family.familyId);
				}
				final char[] sketch=family.sketch();
				if(sketch==null || sketch.length!=width){
					throw new IllegalArgumentException("Retained family lacks an exact-width sketch: "+family.familyId);
				}
				if(family.consensusUnsafe()==null || family.consensusUnsafe().length<1){
					throw new IllegalArgumentException("Retained family lacks a final consensus: "+family.familyId);
				}
				final ConsensusQualityReport.Snapshot snapshot=cohesion.snapshot(family.familyId);
				if(snapshot.members()!=family.size()){
					throw new IllegalArgumentException("Consensus quality member count mismatch for family "+
						family.familyId+": "+snapshot.members()+" != "+family.size());
				}
				byId[family.familyId]=family; resultByFamilyId[family.familyId]=i;
				mutable.upsert(family.familyId,sketch);
			}
			view=mutable.freeze(); results=new Result[families.length];
			workerCount=Math.min(threads,Math.max(1,families.length));
		}
		void run(){
			final Worker[] workers=new Worker[workerCount];
			for(int i=0; i<workers.length; i++){workers[i]=new Worker(this);}
			try{
				if(!ThreadWaiter.startAndWait(Arrays.asList(workers),this)){
					throw new IllegalStateException("Consensus neighbor workers failed");
				}
			}catch(Throwable cause){
				if(failure==null){failure=cause;}
				final ArrayList<Worker> started=new ArrayList<Worker>();
				for(Worker worker : workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null){throw new IllegalStateException("Consensus neighbor workers failed",failure);}
			if(done!=families.length){throw new IllegalStateException("Consensus neighbor workers finished "+done+" of "+families.length+" families");}
		}
		int claim(){
			final long index=next.getAndIncrement();
			return index>=families.length ? families.length : (int)index;
		}
		@Override public void accumulate(Worker worker){
			if(worker.failure!=null && failure==null){failure=worker.failure;}
			done+=worker.done;
			attemptedPairs+=worker.attemptedPairs;
			positiveSpanPairs+=worker.positiveSpanPairs;
		}
		@Override public ReadWriteLock rwlock(){return null;}
		@Override public boolean success(){return failure==null && done==families.length;}
		final Family[] families,byId;
		final int[] resultByFamilyId;
		final ConsensusQualityReport cohesion;
		final FamilySketchBuilder builder;
		final TailAlignmentGate gate;
		final MutableFamilyIndex.FrozenView view;
		final Result[] results;
		final AtomicLong next=new AtomicLong();
		final int workerCount;
		long attemptedPairs,positiveSpanPairs;
		int done;
		Throwable failure;
	}

	private static final class Worker extends Thread implements MutableFamilyIndex.Cosine{
		Worker(Prepared parent){
			this.parent=parent; scorer=new HbmPositionModel.ConsensusScorer(parent.gate.matrix);
			builder=parent.builder.fork(); querySketch=new char[builder.width()];
			candidateIds=new int[11]; neighborGate=new TailAlignmentGate(parent.gate.mode,
				parent.gate.denominator,parent.gate.maxLength,parent.gate.minSpan,0,0,null,parent.gate.matrix,parent.gate.minPaired);
		}
		@Override public void run(){
			try{
				for(int i=parent.claim(); i<parent.families.length; i=parent.claim()){
					parent.results[i]=score(parent.families[i]); done++;
				}
			}catch(Throwable problem){failure=problem;}
		}
		private Result score(Family query){
			final byte[] consensus=query.consensusUnsafe();
			builder.sketchEncoded(consensus,querySketch); setMonomers(consensus);
			final int candidates=parent.view.rankIndexed(querySketch,candidateIds.length,candidateIds,null,this,workspace);
			final Result result=new Result(parent.cohesion.snapshot(query.familyId));
			int kept=0;
			for(int i=0; i<candidates && kept<10; i++){
				final int id=candidateIds[i];
				if(id==query.familyId){continue;}
				kept++;
				final byte[] reference=parent.byId[id].consensusUnsafe();
				scorer.resetConsensus(reference);
				result.observe(id,neighborGate.decide(consensus,reference,scorer));
			}
			attemptedPairs+=result.attemptedPairs;
			positiveSpanPairs+=result.positiveSpanPairs;
			result.finish(); return result;
		}
		private void setMonomers(byte[] consensus){
			Arrays.fill(queryMonomers,0); queryNorm=0;
			for(byte residue : consensus){if(residue>=0 && residue<20){queryMonomers[residue]++;}}
			for(double count : queryMonomers){queryNorm+=count*count;}
		}
		@Override public double score(int familyId){
			final float[] target=parent.byId[familyId].monomersUnsafe();
			double dot=0,norm=0;
			for(int i=0; i<20; i++){dot+=queryMonomers[i]*target[i]; norm+=target[i]*(double)target[i];}
			return queryNorm==0 || norm==0 ? 0 : Math.min(1,Math.max(0,dot/Math.sqrt(queryNorm*norm)));
		}
		final Prepared parent;
		final FamilySketchBuilder builder;
		final TailAlignmentGate neighborGate;
		final HbmPositionModel.ConsensusScorer scorer;
		final MutableFamilyIndex.Workspace workspace=new MutableFamilyIndex.Workspace();
		final int[] candidateIds;
		final char[] querySketch;
		final double[] queryMonomers=new double[20];
		double queryNorm;
		long attemptedPairs,positiveSpanPairs;
		int done;
		Throwable failure;
	}

	private static final class Result{
		Result(ConsensusQualityReport.Snapshot cohesion){
			hasA=cohesion.available(); a=cohesion.identityMean();
		}
		void observe(int familyId,TailAlignmentDecision decision){
			attemptedPairs++;
			if(!decision.accepted()){return;}
			positiveSpanPairs++;
			if(!hasIdentity || decision.identity>identity || decision.identity==identity &&
					(decision.scoreRatio>identityRatio || decision.scoreRatio==identityRatio && familyId<identityNeighbor)){
				hasIdentity=true; identity=decision.identity; identityRatio=decision.scoreRatio; identityNeighbor=familyId;
			}
			if(!hasRatio || decision.scoreRatio>ratio || decision.scoreRatio==ratio && familyId<ratioNeighbor){
				hasRatio=true; ratio=decision.scoreRatio; ratioIdentity=decision.identity; ratioNeighbor=familyId;
			}
		}
		void finish(){
			if(!hasIdentity){
				missingReason=attemptedPairs<1 ? "no_shortlisted_neighbor" :
					positiveSpanPairs<1 ? "no_positive_span_valid_neighbor" : "no_accepted_neighbor";
				return;
			}
			missingReason="ok";
			if(hasA){
				rawDelta=a-identity;
				final double denominator=Math.max(a,identity);
				if(denominator>0){normalizedAvailable=true; normalizedDelta=rawDelta/denominator;}
			}
		}
		final boolean hasA;
		final double a;
		boolean hasIdentity,hasRatio,normalizedAvailable;
		int attemptedPairs,positiveSpanPairs;
		double identity,identityRatio,ratio,ratioIdentity,rawDelta,normalizedDelta;
		int identityNeighbor=-1,ratioNeighbor=-1;
		String missingReason;
	}
}
