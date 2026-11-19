package prot.cluster;

import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.Arrays;

import ddl.MutableFamilyIndex;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Phase-local counters with driver-only worker merges; timers never influence assignment.
 * Heap values are samples, not a minimum-memory estimate.
 * @author Keqing
 */
public final class ClusteringMetrics implements AutoCloseable {

	public ClusteringMetrics(String checkpoints,String histograms,boolean overwrite,int shortlist,long interval){
		this(checkpoints, histograms, overwrite, shortlist, interval, null);
	}

	/** The optional progress stream is borrowed; each diagnostic line is explicitly flushed. */
	public ClusteringMetrics(String checkpoints,String histograms,boolean overwrite,int shortlist,long interval,PrintStream progress){
		if(shortlist<1 || interval<1){throw new IllegalArgumentException("Positive shortlist and checkpoint interval required");}
		this.interval=interval;
		this.progress=progress;
		shortlistHistogram=new long[Math.addExact(shortlist,1)];
		rankHistogram=new long[shortlistHistogram.length];
		checkpointWriter=writer(checkpoints,overwrite,CHECKPOINT_HEADER);
		try{histogramWriter=writer(histograms,overwrite,HISTOGRAM_HEADER);}
		catch(RuntimeException | Error failure){
			if(checkpointWriter!=null){checkpointWriter.poisonAndWait();}
			throw failure;
		}
	}

	/** Attaches the engine-owned index for observational counters, never for mutation. */
	void attachIndex(MutableFamilyIndex index){
		if(active || closed || index==null || rankingIndex!=null){throw new IllegalStateException("Attach one ranking index before starting clustering");}
		rankingIndex=index;
	}

	/** Starts one unique stage/pass. The caller closes the preceding stage first. */
	public void begin(String stage,int pass,long eligible){
		if(active || closed || stage==null || eligible<0 || pass<0 ||
			!(stage.equals("greedy") || stage.equals("seeds") || stage.equals("scavenge") || stage.equals("reassign") ||
				stage.equals("consensus") || stage.equals("merge") || stage.equals("final")) ||
			((stage.equals("scavenge") || stage.equals("seeds")) && pass!=0) || (stage.equals("merge") && pass==0)){
			throw new IllegalStateException("Invalid clustering stage transition: "+stage+" pass="+pass);
		}
		this.stage=stage; this.pass=pass; this.eligible=eligible;
		assignmentAlignments=constructionAlignments=shortlistQueries=shortlistSum=0;
		hitFiltered=hydroFiltered=0;
		seedCandidates=seedAlignments=seedHitFiltered=seedHydroFiltered=seedMatches=0;
		sketchUpdates=consensusUpdates=0;
		Arrays.fill(shortlistHistogram,0); Arrays.fill(rankHistogram,0);
		started=System.nanoTime(); cpuStarted=cpuTime();
		nextCheckpoint=interval; lastProcessed=-1; active=true;
		phaseThread=Thread.currentThread(); begunPhases=Math.incrementExact(begunPhases);
	}

	/** Records all candidates, including zero-hit candidates, and rank zero for no admitted match. */
	public void query(int candidates,int winningRank){
		requireActive();
		if(stage.equals("merge")){throw new IllegalStateException("Merge has no original-input shortlist observations");}
		if(candidates<0 || candidates>=shortlistHistogram.length || winningRank<0 || winningRank>candidates){
			throw new IllegalArgumentException("Winning rank must fit the actual shortlist: "+winningRank+" / "+candidates);
		}
		shortlistQueries=Math.addExact(shortlistQueries,1);
		shortlistSum=Math.addExact(shortlistSum,candidates);
		shortlistHistogram[candidates]=Math.addExact(shortlistHistogram[candidates],1);
		rankHistogram[winningRank]=Math.addExact(rankHistogram[winningRank],1);
	}

	public void alignment(boolean construction){
		alignments(construction, 1);
	}

	/** Disjoint screens: hit failures are counted before hydrophobicity is examined. */
	public void filtered(long hits,long hydro){
		requireActive();
		if(hits<0 || hydro<0){throw new IllegalArgumentException("Prefilter counts cannot be negative");}
		hitFiltered=Math.addExact(hitFiltered,hits); hydroFiltered=Math.addExact(hydroFiltered,hydro);
		totalHitFiltered=Math.addExact(totalHitFiltered,hits); totalHydroFiltered=Math.addExact(totalHydroFiltered,hydro);
	}

	/** Records one competitive-seed candidate outside the initial shortlist.
	 * PASS records its actual alignment call; the caller must not count that call again.
	 * Screened candidates enter both the ordinary work totals and the seed-only subsets.
	 */
	public void seedCandidate(final int screened){
		requireActive();
		if(screened!=CandidateFilter.PASS && screened!=CandidateFilter.LOW_HITS && screened!=CandidateFilter.HYDRO){
			throw new IllegalArgumentException("Unknown competitive-seed screening outcome: "+screened);
		}
		seedCandidates=Math.addExact(seedCandidates,1);
		totalSeedCandidates=Math.addExact(totalSeedCandidates,1);
		if(screened==CandidateFilter.PASS){
			alignment(false);
			seedAlignments=Math.addExact(seedAlignments,1); totalSeedAlignments=Math.addExact(totalSeedAlignments,1);
		}else if(screened==CandidateFilter.LOW_HITS){
			filtered(1,0);
			seedHitFiltered=Math.addExact(seedHitFiltered,1); totalSeedHitFiltered=Math.addExact(totalSeedHitFiltered,1);
		}else{
			filtered(0,1);
			seedHydroFiltered=Math.addExact(seedHydroFiltered,1); totalSeedHydroFiltered=Math.addExact(totalSeedHydroFiltered,1);
		}
	}

	/** Records an admitted competitive-seed candidate, never a newly created family. */
	public void seedMatch(){
		requireActive();
		if(seedMatches>=seedAlignments){throw new IllegalStateException("A seed match requires a recorded seed alignment");}
		seedMatches=Math.addExact(seedMatches,1); totalSeedMatches=Math.addExact(totalSeedMatches,1);
	}

	/** Driver-only merge after a successful worker join and close.
	 * Workers begin one matching phase with no output sinks; close freezes their counters
	 * without requiring a known-at-start partition size. Phase clocks, heap samples,
	 * eligible counts and checkpoint positions belong only to the driver and are not summed.
	 */
	public void mergeCompletedWorker(final ClusteringMetrics local){
		requireActive();
		if(Thread.currentThread()!=phaseThread || local==null || local==this || !local.closed || local.merged ||
			local.begunPhases!=1 || local.checkpointWriter!=null || local.histogramWriter!=null || local.progress!=null ||
			!stage.equals(local.stage) || pass!=local.pass || shortlistHistogram.length!=local.shortlistHistogram.length){
			throw new IllegalStateException("Merge requires one closed outputless worker from this driver phase and shortlist geometry");
		}
		local.validateCounts(false);
		alignments(false,local.assignmentAlignments); alignments(true,local.constructionAlignments);
		filtered(local.hitFiltered,local.hydroFiltered);
		consensusUpdated(local.consensusUpdates);
		sketchUpdates=Math.addExact(sketchUpdates,local.sketchUpdates);
		totalSketchUpdates=Math.addExact(totalSketchUpdates,local.sketchUpdates);
		shortlistQueries=Math.addExact(shortlistQueries,local.shortlistQueries);
		shortlistSum=Math.addExact(shortlistSum,local.shortlistSum);
		for(int i=0; i<shortlistHistogram.length; i++){
			shortlistHistogram[i]=Math.addExact(shortlistHistogram[i],local.shortlistHistogram[i]);
			rankHistogram[i]=Math.addExact(rankHistogram[i],local.rankHistogram[i]);
		}
		seedCandidates=Math.addExact(seedCandidates,local.seedCandidates);
		seedAlignments=Math.addExact(seedAlignments,local.seedAlignments);
		seedHitFiltered=Math.addExact(seedHitFiltered,local.seedHitFiltered);
		seedHydroFiltered=Math.addExact(seedHydroFiltered,local.seedHydroFiltered);
		seedMatches=Math.addExact(seedMatches,local.seedMatches);
		totalSeedCandidates=Math.addExact(totalSeedCandidates,local.seedCandidates);
		totalSeedAlignments=Math.addExact(totalSeedAlignments,local.seedAlignments);
		totalSeedHitFiltered=Math.addExact(totalSeedHitFiltered,local.seedHitFiltered);
		totalSeedHydroFiltered=Math.addExact(totalSeedHydroFiltered,local.seedHydroFiltered);
		totalSeedMatches=Math.addExact(totalSeedMatches,local.seedMatches);
		local.merged=true;
	}

	/** Merges completed worker work on the driver thread, in serial visitation order. */
	public void alignments(boolean construction,long count){
		requireActive();
		if(count<0){throw new IllegalArgumentException("Alignment count cannot be negative");}
		if(stage.equals("merge") && !construction && count>0){
			throw new IllegalArgumentException("Nested-model and union alignments belong to merge construction work");
		}
		if(construction){
			constructionAlignments=Math.addExact(constructionAlignments,count);
			totalConstructionAlignments=Math.addExact(totalConstructionAlignments,count);
		}else{
			assignmentAlignments=Math.addExact(assignmentAlignments,count);
			totalAssignmentAlignments=Math.addExact(totalAssignmentAlignments,count);
		}
	}

	public void sketchUpdated(){
		requireActive(); sketchUpdates=Math.addExact(sketchUpdates,1);
		totalSketchUpdates=Math.addExact(totalSketchUpdates,1);
	}

	public void consensusUpdated(){
		consensusUpdated(1);
	}

	/** Accepted proposals only; rejected and unchanged proposals do not increment this counter. */
	public void consensusUpdated(long count){
		requireActive();
		if(count<0){throw new IllegalArgumentException("Consensus update count cannot be negative");}
		consensusUpdates=Math.addExact(consensusUpdates,count);
		totalConsensusUpdates=Math.addExact(totalConsensusUpdates,count);
	}

	/** Assigned/unassigned partition the whole eligible population, including not-yet-visited inputs. */
	public void checkpoint(long processed,long assigned,long unassigned,int live,long created){
		checkpoint(processed,assigned,unassigned,live,created,false);
	}

	/** Emits driver-observed progress while concurrent workers are still running.
	 * These independently read values are not a coherent scientific checkpoint.
	 * No counters, timers, heap samples, checkpoint positions or TSV files are changed.
	 */
	public void concurrentProgress(final long processed, final int publishedFamilies){
		requireActive();
		if(Thread.currentThread()!=phaseThread || !stage.equals("greedy") || processed<0 || processed>eligible ||
			publishedFamilies<0 || publishedFamilies>eligible){
			throw new IllegalArgumentException("Concurrent progress needs driver-owned greedy observations within the eligible population");
		}
		if(progress==null){return;}
		final ByteBuilder line=new ByteBuilder(256).append("clustering_concurrent_progress\tstage=").append(stage)
			.append("\tpass=").append(pass).append("\tcompleted_inputs=").append(processed)
			.append("\tpublished_families=").append(publishedFamilies)
			.append("\telapsed_seconds=").append((System.nanoTime()-started)*1e-9,6)
			.append("\tobservation=live_not_coherent_checkpoint");
		progress.println(line.toString()); progress.flush();
		if(progress.checkError()){throw new IllegalStateException("Concurrent clustering progress failed to write");}
	}

	/** Completes a stage; processed counts accounted input ordinals, not alignment invocations. */
	public void end(long processed,long assigned,long unassigned,int live,long created){
		requireActive();
		if(processed!=eligible){throw new IllegalStateException("Incomplete stage "+stage+": "+processed+" / "+eligible);}
		checkpoint(processed,assigned,unassigned,live,created,true);
		validateCounts(true);
		active=false;
	}

	/** Validate disjoint shortlist and competitive-seed work; only the driver writes histograms. */
	private void validateCounts(final boolean write){
		long queries=0,ranks=0,weighted=0;
		for(int i=0; i<shortlistHistogram.length; i++){
			queries=Math.addExact(queries,shortlistHistogram[i]);
			ranks=Math.addExact(ranks,rankHistogram[i]);
			weighted=Math.addExact(weighted,Math.multiplyExact((long)i,shortlistHistogram[i]));
			if(write){
				writeHistogram("shortlist_size",i,shortlistHistogram[i]);
				writeHistogram("match_rank",i,rankHistogram[i]);
			}
		}
		if(queries!=shortlistQueries || ranks!=queries || weighted!=shortlistSum){
			throw new IllegalStateException("Shortlist/rank histograms must conserve every classification query");
		}
		if(seedAlignments>assignmentAlignments || seedHitFiltered>hitFiltered || seedHydroFiltered>hydroFiltered ||
			seedMatches>seedAlignments || Math.addExact(seedAlignments,Math.addExact(seedHitFiltered,seedHydroFiltered))!=seedCandidates){
			throw new IllegalStateException("Competitive-seed subsets must conserve candidates and actual alignment work");
		}
		if(Math.addExact(assignmentAlignments-seedAlignments,
			Math.addExact(hitFiltered-seedHitFiltered,hydroFiltered-seedHydroFiltered))!=shortlistSum){
			throw new IllegalStateException("Every ranked candidate must be aligned or rejected by exactly one prefilter in "+stage+"/"+pass);
		}
	}

	private void checkpoint(long processed,long assigned,long unassigned,int live,long created,boolean force){
		requireActive();
		if(processed<0 || processed>eligible || processed<lastProcessed || assigned<0 || unassigned<0 ||
			Math.addExact(assigned,unassigned)!=eligible || live<0 || created<live){
			throw new IllegalStateException("Checkpoint must conserve eligible input ordinals in "+stage);
		}
		if(!force && processed<nextCheckpoint){return;}
		final Runtime runtime=Runtime.getRuntime();
		final long committed=runtime.totalMemory(),used=committed-runtime.freeMemory();
		sampledPeak=Math.max(sampledPeak,used);
		final long now=System.nanoTime(),cpu=cpuTime();
		final boolean measuredCpu=cpuStarted>=0 && cpu>=cpuStarted;
		if(checkpointWriter!=null){
			final ByteBuilder row=new ByteBuilder(256);
			row.append(stage).tab().append(pass).tab().append(processed).tab().append(eligible).tab()
				.append(assigned).tab().append(unassigned).tab().append(live).tab().append(created).tab()
				.append(Math.addExact(assignmentAlignments,constructionAlignments)).tab()
				.append(assignmentAlignments).tab().append(constructionAlignments).tab()
				.append(shortlistQueries).tab().append(shortlistSum).tab().append(sketchUpdates).tab()
				.append(consensusUpdates).tab().append((now-started)*1e-9, 6).tab()
				.append(measuredCpu ? (cpu-cpuStarted)*1e-9 : 0, 6).tab().append(measuredCpu ? 1 : 0).tab()
				.append(used).tab().append(committed).tab().append(runtime.maxMemory()).tab().append(sampledPeak);
			checkpointWriter.println(row);
		}
		if(progress!=null){
			final ByteBuilder line=new ByteBuilder(256).append("clustering_progress\tstage=").append(stage)
				.append("\tpass=").append(pass).append("\tprocessed=").append(processed)
				.append("\teligible=").append(eligible).append("\tassigned=").append(assigned)
				.append("\tclusters_live=").append(live).append("\tclusters_created=").append(created)
				.append("\telapsed_seconds=").append((now-started)*1e-9, 6)
				.append("\theap_used_bytes=").append(used).append("\theap_max_bytes=").append(runtime.maxMemory());
			if(rankingIndex!=null){
				line.append("\tindex_live_postings=").append(rankingIndex.livePostings())
					.append("\tindex_delta_records=").append(rankingIndex.deltaRecords())
					.append("\tindex_compactions=").append(rankingIndex.compactions())
					.append("\tindex_primitive_bytes=").append(rankingIndex.primitivePayloadBytes());
			}
			progress.println(line.toString()); progress.flush();
			if(progress.checkError()){throw new IllegalStateException("Clustering progress diagnostics failed to write");}
		}
		lastProcessed=processed;
		nextCheckpoint=processed>Long.MAX_VALUE-interval ? Long.MAX_VALUE : processed+interval;
	}

	private void writeHistogram(String metric,int bin,long count){
		if(histogramWriter!=null && count>0){
			histogramWriter.println(new ByteBuilder().append(stage).tab().append(pass).tab()
				.append(metric).tab().append(bin).tab().append(count));
		}
	}

	private void requireActive(){
		if(!active || closed){throw new IllegalStateException("Counters require an open clustering stage");}
	}

	private static long cpuTime(){
		return CPU.isCurrentThreadCpuTimeSupported() && CPU.isThreadCpuTimeEnabled() ? CPU.getCurrentThreadCpuTime() : -1;
	}

	private static ByteStreamWriter writer(String path,boolean overwrite,String header){
		if(path==null){return null;}
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,true);
		writer.start(); writer.println(header); return writer;
	}

	@Override
	public void close(){
		if(closed){return;}
		closed=true;
		boolean failed=false;
		if(checkpointWriter!=null){failed|=checkpointWriter.poisonAndWait();}
		if(histogramWriter!=null){failed|=histogramWriter.poisonAndWait();}
		if(failed){throw new IllegalStateException("Clustering diagnostic output failed");}
	}

	public long totalAlignments(){return Math.addExact(totalAssignmentAlignments,totalConstructionAlignments);}
	/** All assignment alignment attempts, including the competitive-seed subset. */
	public long assignmentAlignments(){return totalAssignmentAlignments;}
	public long constructionAlignments(){return totalConstructionAlignments;}
	public long sketchUpdates(){return totalSketchUpdates;}
	public long consensusUpdates(){return totalConsensusUpdates;}
	/** All hit-screened candidates, including seedHitFiltered(). */
	public long hitFiltered(){return totalHitFiltered;}
	/** All hydrophobicity-screened candidates, including seedHydroFiltered(). */
	public long hydroFiltered(){return totalHydroFiltered;}
	/** Competitive-seed totals are subsets of ordinary alignment/filter work, except matches/candidates. */
	public long seedCandidates(){return totalSeedCandidates;}
	public long seedAlignments(){return totalSeedAlignments;}
	public long seedHitFiltered(){return totalSeedHitFiltered;}
	public long seedHydroFiltered(){return totalSeedHydroFiltered;}
	public long seedMatches(){return totalSeedMatches;}

	public static final String CHECKPOINT_HEADER="stage\tpass\tprocessed\teligible\tassigned\tunassigned\tclusters_live\tclusters_created\talignments\tassignment_alignments\tconstruction_alignments\tshortlist_queries\tshortlist_sum\tsketch_updates\tconsensus_updates\telapsed_seconds\tdriver_cpu_seconds\tcpu_available\theap_used_bytes\theap_committed_bytes\theap_max_bytes\tsampled_peak_used_heap_bytes";
	public static final String HISTOGRAM_HEADER="stage\tpass\tmetric\tbin\tcount";
	private static final ThreadMXBean CPU=ManagementFactory.getThreadMXBean();
	private final ByteStreamWriter checkpointWriter,histogramWriter;
	private final PrintStream progress;
	private MutableFamilyIndex rankingIndex;
	private final long[] shortlistHistogram,rankHistogram;
	private final long interval;
	private String stage;
	private int pass,begunPhases;
	private boolean active,closed,merged;
	private Thread phaseThread;
	private long eligible,started,cpuStarted,nextCheckpoint,lastProcessed,sampledPeak;
	private long assignmentAlignments,constructionAlignments,shortlistQueries,shortlistSum,sketchUpdates,consensusUpdates;
	private long totalAssignmentAlignments,totalConstructionAlignments,totalSketchUpdates,totalConsensusUpdates;
	private long hitFiltered,hydroFiltered,totalHitFiltered,totalHydroFiltered;
	private long seedCandidates,seedAlignments,seedHitFiltered,seedHydroFiltered,seedMatches;
	private long totalSeedCandidates,totalSeedAlignments,totalSeedHitFiltered,totalSeedHydroFiltered,totalSeedMatches;
}
