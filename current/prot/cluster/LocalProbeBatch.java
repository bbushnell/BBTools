package prot.cluster;

import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;
import template.Accumulator;
import template.ThreadWaiter;

/** Independent local responses from one frozen baseline and persisted stable IDs.
 * Workers construct and evaluate proposals concurrently with private caches and
 * scratch over shared immutable family and pool indexes.
 * Only small result rows survive a probe; no per-probe genome reassignment or
 * replacement codec export occurs. Output is opened only after all workers join.
 * @author Brian Bushnell, Keqing
 */
final class LocalProbeBatch{

	static void run(final StructuralProbeBatch.Plan plan,final LocalProbeBaseline baseline,
			final LayerCakeClustering engine,final LocalProbeEvaluator evaluator,final int threads,
			final double temperature,final String output,final String rebuildOutput,
			final boolean overwrite){
		run(plan,baseline,engine,evaluator,threads,temperature,output,rebuildOutput,null,overwrite);
	}

	static void run(final StructuralProbeBatch.Plan plan,final LocalProbeBaseline baseline,
			final LayerCakeClustering engine,final LocalProbeEvaluator evaluator,final int threads,
			final double temperature,final String output,final String rebuildOutput,final String auditOutput,
			final boolean overwrite){
		run(plan,baseline,engine,evaluator,threads,temperature,output,rebuildOutput,auditOutput,null,overwrite);
	}

	static void run(final StructuralProbeBatch.Plan plan,final LocalProbeBaseline baseline,
			final LayerCakeClustering engine,final LocalProbeEvaluator evaluator,final int threads,
			final double temperature,final String output,final String rebuildOutput,final String auditOutput,
			final String featuresOutput,final boolean overwrite){
		run(plan,baseline,engine,evaluator,threads,temperature,output,rebuildOutput,auditOutput,featuresOutput,null,overwrite);
	}

	static void run(final StructuralProbeBatch.Plan plan,final LocalProbeBaseline baseline,
			final LayerCakeClustering engine,final LocalProbeEvaluator evaluator,final int threads,
			final double temperature,final String output,final String rebuildOutput,final String auditOutput,
			final String featuresOutput,final String refinementOutput,final boolean overwrite){
		if(plan==null || baseline==null || engine==null || evaluator==null || threads<1 ||
			!Double.isFinite(temperature) || temperature<0 || output==null){
			throw new IllegalArgumentException("Local batch requires a frozen baseline, engine, evaluator, output and valid worker/temperature settings");
		}
		if(featuresOutput!=null && auditOutput==null){throw new IllegalArgumentException("Feature export requires explicit audit outcomes");}
		validate(plan,baseline);
		final long setupStart=System.nanoTime();
		final Work work=new Work(plan,baseline,engine,evaluator,threads,temperature,rebuildOutput!=null,auditOutput!=null,featuresOutput!=null);
		final long evaluationStart=System.nanoTime();
		work.run();
		final long evaluationEnd=System.nanoTime();
		System.err.println("Local probe worker/index preparation seconds: "+((evaluationStart-setupStart)*1e-9));
		System.err.println("Local probe construction, scoring and features seconds: "+((evaluationEnd-evaluationStart)*1e-9));
		System.err.println("Local probe concurrent workers: "+work.workers.size()+"; requests: "+plan.size());
		final ByteStreamWriter writer=new ByteStreamWriter(output,overwrite,false,true); writer.start();
		try{
			writer.println(HEADER); final ByteBuilder row=new ByteBuilder();
			for(int i=0;i<plan.rows.length;i++){
				final StructuralProbeBatch.Request request=plan.rows[i]; final Outcome result=work.results[i];
				row.clear().append("local_replacements_v1").tab().append(request.id).tab().append(request.seed).tab()
					.append(request.phase.name().toLowerCase(Locale.ROOT)).tab().append(request.left).tab().append(request.right).tab()
					.append(result.features==null ? 0 : 1).tab().append(result.reason).tab().append(result.affected).tab()
					.append(evaluator.poolSize()).tab().append(result.slotCost).tab().appendSlow(temperature);
				if(result.features==null){
					for(int missing=0;missing<6;missing++){row.tab().append("NA");}
				}else{
					final LayerCakeClustering.ProbeFeatures f=result.features;
					row.tab().appendSlow(f.balance).tab().appendSlow(f.cohesionDelta).tab().appendSlow(f.separationDelta)
						.tab().appendSlow(f.sizePressure).tab().appendSlow(f.slotDelta).tab().append(result.response.netGain);
				}
				if(result.response==null){
					for(int missing=0;missing<10;missing++){row.tab().append("NA");}
				}else{
					final LocalProbeEvaluator.Result r=result.response;
					row.tab().append(r.poolVisited).tab().append(r.poolAccepted).tab().append(r.affectedAccepted).tab()
						.append(r.accepted).tab().append(r.candidates).tab().append(r.alignments).tab()
						.append(r.lowHitRejections).tab().append(r.hydroRejections).tab().append(r.postingsVisited).tab()
						.append(r.affected-r.affectedAccepted);
				}
				row.tab().append(result.elapsedNanos); writer.println(row);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Local response output failed: "+output);}}
		if(rebuildOutput!=null){writeRebuildReport(plan,work.results,rebuildOutput,overwrite);}
		if(auditOutput!=null){writeAuditReport(plan,work.results,auditOutput,overwrite);}
		if(featuresOutput!=null){writeFeatures(plan,work.results,featuresOutput,overwrite);}
		if(refinementOutput!=null){writeRefinement(plan,work.results,refinementOutput,overwrite);}
	}

	/** Forced construction membership is separate from successful final gene-tail assignment. */
	private static void writeRefinement(final StructuralProbeBatch.Plan plan,final Outcome[] results,
			final String output,final boolean overwrite){
		assert(results.length==plan.rows.length) : "Refinement report identities must stay aligned with every requested outcome";
		final ByteStreamWriter writer=new ByteStreamWriter(output,overwrite,false,true);writer.start();
		try{
			writer.println(REFINEMENT_HEADER);final ByteBuilder row=new ByteBuilder();
			for(int i=0;i<plan.rows.length;i++){
				final StructuralProbeBatch.Request request=plan.rows[i];final Outcome result=results[i];
				final CompetitiveRefinement.Summary s=result.refinement;
				final boolean available=result.audit!=null && result.audit.response!=null;
				row.clear().append("local_refinement_v1").tab().append(request.id).tab().append(request.seed).tab()
					.append(request.phase.name().toLowerCase(Locale.ROOT)).tab().append(request.left).tab().append(request.right)
					.tab().append(s==null ? "not_attempted" : available ? "completed" : "unavailable")
					.tab().append(result.reason).tab().append(s==null ? 0 : s.rounds).tab().append(s==null ? 0 : s.moved)
					.tab().append(s==null ? 0 : s.forcedRetention).tab().append(s==null ? 0 : s.assignmentAlignments)
					.tab().append(s==null ? 0 : s.rebuildAlignments).tab().append(result.elapsedNanos);
				writer.println(row);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Refinement diagnostic output failed: "+output);}}
	}

	/** One row per request, including honest missingness for unconstructible candidates. */
	private static void writeFeatures(final StructuralProbeBatch.Plan plan,final Outcome[] results,
			final String output,final boolean overwrite){
		final ByteStreamWriter writer=new ByteStreamWriter(output,overwrite,false,true);writer.start();
		try{
			writer.println(FEATURE_HEADER);final ByteBuilder row=new ByteBuilder();
			for(int i=0;i<plan.rows.length;i++){
				final StructuralProbeBatch.Request request=plan.rows[i];final Outcome result=results[i];
				final LayerCakeClustering.ProbeFeatures f=result.measured;
				assert(result.audit!=null) : "Feature rows require the corresponding actual audit outcome";
				row.clear().append("local_features_v1").tab().append(request.id).tab().append(request.seed).tab()
					.append(request.phase.name().toLowerCase(Locale.ROOT)).tab().append(request.left).tab().append(request.right)
					.tab().append(result.audit.status).tab().append(f==null ? "unavailable" : "available").tab()
					.append(f==null && result.audit.response!=null ? "quality_unavailable" : result.audit.reason);
				if(f==null){for(int col=9;col<20;col++){row.tab().append("NA");}}
				else{
					row.tab().append(f.parent.members).tab().append(f.left.members).tab().append(f.right.members)
						.tab().appendSlow(f.parent.cohesion).tab().appendSlow(f.parent.neighbor)
						.tab().appendSlow(f.left.cohesion).tab().appendSlow(f.left.neighbor)
						.tab().appendSlow(f.right.cohesion).tab().appendSlow(f.right.neighbor)
						.tab().appendSlow(f.cohesionDelta).tab().appendSlow(f.separationDelta);
				}
				row.tab();if(result.audit.response==null){row.append("NA");}else{row.append(result.audit.response.netGain);}
				pair(row,f==null ? null : f.leftPair);pair(row,f==null ? null : f.rightPair);
				row.tab().append(result.elapsedNanos);writer.println(row);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Local feature output failed: "+output);}}
	}

	private static void pair(final ByteBuilder row,final FamilyNeighborEvidence.Result pair){
		assert(row!=null) : "Directed pair fields must append to their existing request row";
		final boolean available=pair!=null && pair.explicitAvailable;
		row.tab().append(available ? 1 : 0);
		if(!available){for(int i=0;i<4;i++){row.tab().append("NA");}}
		else{row.tab().appendSlow(pair.explicitIdentity).tab().appendSlow(pair.explicitScoreRatio)
			.tab().appendSlow(pair.explicitSharedFraction).tab().appendSlow(pair.explicitQuality);}
	}

	/** Writes scalar audit outcomes only; no candidate model survives worker evaluation. */
	private static void writeAuditReport(final StructuralProbeBatch.Plan plan,final Outcome[] results,
			final String output,final boolean overwrite){
		final ByteStreamWriter writer=new ByteStreamWriter(output,overwrite,false,true); writer.start();
		try{
			writer.println(AUDIT_HEADER); final ByteBuilder row=new ByteBuilder();
			for(int i=0;i<plan.rows.length;i++){
				final StructuralProbeBatch.Request request=plan.rows[i]; final Outcome result=results[i];
				final AuditOutcome audit=result.audit;
				if(audit==null){throw new IllegalStateException("Every requested audit needs an explicit scalar outcome");}
				row.clear().append("local_reject_audit_v1").tab().append(request.id).tab().append(request.seed).tab()
					.append(request.phase.name().toLowerCase(Locale.ROOT)).tab().append(request.left).tab().append(request.right)
					.tab().append(result.reason).tab().append(audit.status).tab().append(audit.reason)
					.tab().append(audit.leftState).tab().append(audit.rightState).tab().append(audit.unionState);
				optional(row,audit.leftMembers); optional(row,audit.rightMembers); optional(row,audit.unionMembers);
				optional(row,audit.leftLength); optional(row,audit.rightLength); optional(row,audit.unionLength);
				final LocalProbeEvaluator.Result r=audit.response;
				if(r==null){for(int j=0;j<11;j++){row.tab().append("NA");}}
				else{
					row.tab().append(r.poolVisited).tab().append(r.poolAccepted).tab().append(r.affectedAccepted)
						.tab().append(r.accepted).tab().append(r.netGain).tab().append(r.candidates).tab().append(r.alignments)
						.tab().append(r.lowHitRejections).tab().append(r.hydroRejections).tab().append(r.postingsVisited)
						.tab().append(r.affected-r.affectedAccepted);
				}
				row.tab().append(audit.rebuilds).tab().append(audit.accepted).tab().append(audit.rejected)
					.tab().append(audit.alignments).tab().append(result.elapsedNanos);
				writer.println(row);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Local audit output failed: "+output);}}
	}

	private static void optional(final ByteBuilder row,final int value){
		assert(value>=-1) : "Audit lengths/member counts use -1 only for explicit missingness";
		row.tab(); if(value<0){row.append("NA");}else{row.append(value);}
	}

	private static void writeRebuildReport(final StructuralProbeBatch.Plan plan,final Outcome[] results,
			final String output,final boolean overwrite){
		final ByteStreamWriter writer=new ByteStreamWriter(output,overwrite,false,true); writer.start();
		try{
			writer.println(REBUILD_HEADER); final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<plan.rows.length; i++){
				final StructuralProbeBatch.Request request=plan.rows[i];
				final Outcome result=results[i];
				if(request.phase==LayerCakeSchedule.Phase.FISSION){
					writeRebuildRow(writer,row,request,result,"left",result.leftRebuild);
					writeRebuildRow(writer,row,request,result,"right",result.rightRebuild);
				}else{
					writeRebuildRow(writer,row,request,result,"union",result.unionRebuild);
				}
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Local rebuild report failed: "+output);}}
	}

	private static void writeRebuildRow(final ByteStreamWriter writer,final ByteBuilder row,
			final StructuralProbeBatch.Request request,final Outcome result,
			final String side,final ConsensusRefiner.RebuildSummary summary){
		row.clear().append("local_rebuild_v1").tab().append(request.id).tab().append(request.seed).tab()
			.append(request.phase.name().toLowerCase(Locale.ROOT)).tab().append(request.left).tab().append(request.right)
			.tab().append(result.reason).tab().append(side).tab().append(summary==null ? "skipped" : "attempted");
		if(summary==null){
			row.tab().append("NA").tab().append("NA").tab().append("NA")
				.tab().append("NA").tab().append("NA").tab().append("NA").tab().append("NA")
				.tab().append("NA").tab().append("NA").tab().append("unavailable");
		}else{
			row.tab().append(summary.members).tab().append(summary.oldLength).tab();
			if(summary.proposedLength<0){row.append("NA");}
			else{row.append(summary.proposedLength);}
			row.tab().append(summary.proposals).tab().append(summary.accepted).tab().append(summary.unchanged)
				.tab().append(summary.rejected).tab().append(summary.rejectedMemberChecks)
				.tab().append(summary.incompatibleSnapshots).tab();
			summary.appendReasonCounts(row);
		}
		writer.println(row);
	}

	/** Reject stale or empty endpoints before the driver builds the large pool index. */
	static void validate(final StructuralProbeBatch.Plan plan,final LocalProbeBaseline baseline){
		if(plan==null || baseline==null){throw new IllegalArgumentException("Local requests require a selected baseline");}
		for(StructuralProbeBatch.Request request:plan.rows){
			endpoint(baseline,request.left);
			if(request.right>=0){endpoint(baseline,request.right);}
		}
	}

	private static Family endpoint(final LocalProbeBaseline baseline,final int stableId){
		final int at=baseline.selectedIndex(stableId);
		if(at<0 || baseline.selected[at].size()==0){
			throw new IllegalArgumentException("A probe endpoint needs current members in the selected baseline: "+stableId);
		}
		return baseline.selected[at];
	}

	private static final class Work implements Accumulator<Worker>{
		Work(StructuralProbeBatch.Plan plan_,LocalProbeBaseline baseline_,LayerCakeClustering engine_,
				LocalProbeEvaluator evaluator,int threads,double temperature_,boolean rebuildReport_,boolean auditReport_,boolean featuresReport_){
			plan=plan_;baseline=baseline_;engine=engine_;temperature=temperature_;rebuildReport=rebuildReport_;auditReport=auditReport_;
			featuresReport=featuresReport_;
			results=new Outcome[plan.size()];
			for(int i=0;i<Math.min(threads,plan.size());i++){workers.add(new Worker(this,evaluator.newWorker(),engine.forkLocalProbes()));}
		}
		void run(){
			try{
				if(!ThreadWaiter.startAndWait(workers,this)){fail(new IllegalStateException("Local batch workers did not finish"));}
			}catch(Throwable cause){
				fail(cause); final ArrayList<Worker> started=new ArrayList<Worker>();
				for(Worker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null || completed!=results.length){throw new IllegalStateException("Local batch failed before result publication",failure);}
			for(Outcome result:results){if(result==null){throw new IllegalStateException("Local batch omitted a requested probe");}}
		}
		int claim(){return (int)Math.min(next.getAndIncrement(),results.length);}
		synchronized void fail(Throwable cause){if(failure==null){failure=cause;}}
		@Override public void accumulate(Worker worker){
			if(!worker.success){fail(new IllegalStateException("Local response worker stopped early"));}
			completed=Math.addExact(completed,worker.completed);
		}
		@Override public boolean success(){return failure==null && completed==results.length;}
		@Override public ReadWriteLock rwlock(){return null;}
		final StructuralProbeBatch.Plan plan;
		final LocalProbeBaseline baseline;
		final LayerCakeClustering engine;
		final double temperature;
		final boolean rebuildReport,auditReport,featuresReport;
		final Outcome[] results;
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		final AtomicLong next=new AtomicLong();
		volatile Throwable failure;
		int completed;
	}

	private static final class Worker extends Thread{
		Worker(Work owner_,LocalProbeEvaluator.Worker evaluator_,LayerCakeClustering engine_){owner=owner_;evaluator=evaluator_;engine=engine_;}
		@Override public void run(){
			try{
				for(int index=owner.claim();index<owner.results.length && owner.failure==null;index=owner.claim()){
					final long start=System.nanoTime();
					final StructuralProbeBatch.Request request=owner.plan.rows[index];
					final Family a=endpoint(owner.baseline,request.left),b=request.right<0 ? null : endpoint(owner.baseline,request.right);
					final int affected=Math.addExact(a.size(),b==null ? 0 : b.size());
					final LayerCakeClustering.Probe proposal=engine.probeByFamilyId(request.phase,
						request.left,request.right,owner.temperature,owner.rebuildReport,owner.auditReport);
					LocalProbeEvaluator.Result response=null;
					final LayerCakeClustering.Audit candidate=proposal.audit;
					final boolean auditAvailable=owner.auditReport && candidate!=null && candidate.available();
					int[] members=null;
					if(proposal.feasible() || auditAvailable){
						members=new int[affected];
						for(int i=0;i<a.size();i++){members[i]=a.memberIndex(i);}
						if(b!=null){for(int i=0;i<b.size();i++){members[a.size()+i]=b.memberIndex(i);}}
					}
					if(proposal.feasible()){
						final boolean split=request.phase==LayerCakeSchedule.Phase.FISSION;
						response=evaluator.evaluate(split ? proposal.left : proposal.parent,split ? proposal.right : null,members,null);
						if(response.affected!=affected){throw new IllegalStateException("Local response must charge exactly the removed current members");}
					}
					LocalProbeEvaluator.Result auditResponse=response;
					if(!proposal.feasible() && auditAvailable){
						auditResponse=evaluator.evaluate(candidate.fission ? candidate.left : candidate.union,
							candidate.fission ? candidate.right : null,members,null);
						if(auditResponse.affected!=affected){throw new IllegalStateException("Audit response must charge the same removed current members");}
					}
					final AuditOutcome audit=owner.auditReport ? new AuditOutcome(proposal,auditResponse) : null;
					final LayerCakeClustering.ProbeFeatures measured=owner.featuresReport ? engine.measureLocalFeatures(proposal) : null;
					owner.results[index]=new Outcome(proposal,response,audit,measured,affected,System.nanoTime()-start); completed++;
				}
				success=owner.failure==null;
			}catch(Throwable cause){owner.fail(cause);}
		}
		final Work owner;
		final LocalProbeEvaluator.Worker evaluator;
		final LayerCakeClustering engine;
		int completed;
		boolean success;
	}

	private static final class Outcome{
		Outcome(LayerCakeClustering.Probe proposal,LocalProbeEvaluator.Result response_,AuditOutcome audit_,
				LayerCakeClustering.ProbeFeatures measured_,int affected_,long elapsed_){
			features=proposal.features;reason=proposal.reason;slotCost=proposal.slotCost;
			leftRebuild=proposal.leftRebuild;rightRebuild=proposal.rightRebuild;unionRebuild=proposal.unionRebuild;
			response=response_;audit=audit_;affected=affected_;elapsedNanos=elapsed_;
			measured=measured_;
			refinement=proposal.refinement;
			if((features==null)!=(response==null) || elapsedNanos<0){throw new IllegalStateException("Local feasibility and response availability must agree");}
		}
		final LayerCakeClustering.ProbeFeatures features;
		final LayerCakeClustering.ProbeFeatures measured;
		final CompetitiveRefinement.Summary refinement;
		final LocalProbeEvaluator.Result response;
		final AuditOutcome audit;
		final ConsensusRefiner.RebuildSummary leftRebuild,rightRebuild,unionRebuild;
		final String reason;
		final int affected,slotCost;
		final long elapsedNanos;
	}

	/** Copies metadata and immutable counts; deliberately contains no Family or sequence fields. */
	private static final class AuditOutcome{
		AuditOutcome(final LayerCakeClustering.Probe proposal,final LocalProbeEvaluator.Result response_){
			final LayerCakeClustering.Audit candidate=proposal.audit;
			final boolean split=proposal.phase==LayerCakeSchedule.Phase.FISSION;
			status=proposal.feasible() ? "feasible_control" : candidate!=null && candidate.available() ? "audited" : "unavailable";
			reason=proposal.feasible() ? "ordinary_feasible" : candidate==null ? proposal.reason :
				candidate.available() ? "strict_rebuild_veto" : candidate.unavailableReason;
			response=response_;
			if((response==null)!=status.equals("unavailable") || reason==null || reason.isEmpty()){
				throw new IllegalStateException("Audit response availability must agree with its explicit status and reason");
			}
			final ConsensusRefiner.RebuildSummary left=candidate==null ? proposal.leftRebuild : candidate.leftRebuild;
			final ConsensusRefiner.RebuildSummary right=candidate==null ? proposal.rightRebuild : candidate.rightRebuild;
			final ConsensusRefiner.RebuildSummary union=candidate==null ? proposal.unionRebuild : candidate.unionRebuild;
			leftState=split ? state(left) : "NA";rightState=split ? state(right) : "NA";unionState=split ? "NA" : state(union);
			leftMembers=left==null ? -1 : left.members;rightMembers=right==null ? -1 : right.members;unionMembers=union==null ? -1 : union.members;
			leftLength=left==null ? -1 : left.proposedLength;rightLength=right==null ? -1 : right.proposedLength;
			unionLength=union==null ? -1 : union.proposedLength;
			rebuilds=candidate==null ? 0 : candidate.rebuilds;accepted=candidate==null ? 0 : candidate.accepted;
			rejected=candidate==null ? 0 : candidate.rejected;alignments=candidate==null ? 0 : candidate.alignments;
		}
		private static String state(final ConsensusRefiner.RebuildSummary summary){
			if(summary==null){return "skipped";}
			if(summary.proposals!=1 || summary.accepted+summary.unchanged+summary.rejected!=1){
				throw new IllegalStateException("Audit side state requires one completed native rebuild proposal");
			}
			return summary.accepted==1 ? "accepted" : summary.unchanged==1 ? "unchanged" : "rejected";
		}
		final String status,reason,leftState,rightState,unionState;
		final int leftMembers,rightMembers,unionMembers,leftLength,rightLength,unionLength;
		final LocalProbeEvaluator.Result response;
		final long rebuilds,accepted,rejected,alignments;
	}

	static final String HEADER="response_kind\tcandidate_index\tsample_seed\toperation\tleft_id\tright_id\tfeasible\treason\taffected_members\tpool_size\tslot_cost\ttemperature\tbalance\tweighted_cohesion_delta\tweighted_separation_delta\tsize_pressure\tslot_delta\tnet_gain\tpool_visited\tpool_accepted\taffected_accepted\taccepted\tcandidates\talignments\tlow_hit_rejections\thydro_rejections\tpostings_visited\taffected_lost\telapsed_ns";
	static final String FEATURE_HEADER="report_kind\tcandidate_index\tsample_seed\toperation\tleft_id\tright_id\tstatus\tfeature_status\treason\tparent_members\tleft_members\tright_members\tparent_q\tparent_n\tleft_q\tleft_n\tright_q\tright_n\tweighted_cohesion_delta\tweighted_separation_delta\tnet_gain\tleft_pair_available\tleft_pair_identity\tleft_pair_score_ratio\tleft_pair_shared_fraction\tleft_pair_quality\tright_pair_available\tright_pair_identity\tright_pair_score_ratio\tright_pair_shared_fraction\tright_pair_quality\telapsed_ns";
	static final String REFINEMENT_HEADER="report_kind\tcandidate_index\tsample_seed\toperation\tleft_id\tright_id\tstatus\treason\trounds\tmoved\tforced_retention\tassignment_alignments\trebuild_alignments\telapsed_ns";
	static final String REBUILD_HEADER="report_kind\tcandidate_index\tsample_seed\toperation\tleft_id\tright_id\ttop_reason\tside\tstate\tconstruction_members\told_length\tproposed_length\tproposals\taccepted\tunchanged\trejected\trejected_member_checks\tincompatible_snapshots\treason_counts";
	static final String AUDIT_HEADER="report_kind\tcandidate_index\tsample_seed\toperation\tleft_id\tright_id\ttop_reason\taudit_status\taudit_reason\tleft_state\tright_state\tunion_state\tleft_members\tright_members\tunion_members\tleft_length\tright_length\tunion_length\tpool_visited\tpool_accepted\taffected_accepted\taccepted\tnet_gain\tcandidates\talignments\tlow_hit_rejections\thydro_rejections\tpostings_visited\taffected_lost\taudit_rebuilds\taudit_accepted\taudit_rejected\taudit_alignments\telapsed_ns";
}
