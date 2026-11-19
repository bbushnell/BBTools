package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.locks.ReadWriteLock;

import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import idaligner.QuantumAligner;
import map.LongHashSet;
import stream.Read;
import structures.IntList;
import template.Accumulator;
import template.ThreadWaiter;

/** Serial growing-family admission followed by frozen parallel reassignment.
 * Member ordinals and family IDs stay stable. The three policies share the
 * same input order, Quantum gate, rebuild schedule and final membership census.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlClusterer implements Accumulator<NcrnaDdlClusterer.ProcessThread>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlClusterer(final ArrayList<Read> input_, final NcrnaClusterAssignment.Policy policy_,
			final int patience_, final int bins_, final long seed_, final LongHashSet blacklist_,
			final float initialIdentity_, final float recruitIdentity_, final int rounds_,
			final int passes_, final boolean trim_, final boolean recruit_, final int threads_){
		if(input_==null || input_.isEmpty() || policy_==null || policy_!=NcrnaClusterAssignment.Policy.ALL && blacklist_==null || threads_<1 || patience_<0
				|| rounds_<0 || passes_<0 || !validIdentity(initialIdentity_) || !validIdentity(recruitIdentity_)){
			throw new IllegalArgumentException("Nucleotide clustering needs real input, explicit policies and finite identity/rebuild bounds");
		}
		input=input_;
		policy=policy_;
		patience=patience_;
		bins=policy==NcrnaClusterAssignment.Policy.ALL ? 1 : bins_;
		seed=seed_;
		blacklist=blacklist_;
		initialIdentity=initialIdentity_;
		recruitIdentity=recruitIdentity_;
		rounds=rounds_;
		passes=passes_;
		trim=trim_;
		recruit=recruit_;
		threads=threads_;
		assignments=new int[input.size()];
		Arrays.fill(assignments, -1);
		sketches=policy==NcrnaClusterAssignment.Policy.ALL ? null : new NcrnaDdlSketchBuilder(bins, seed, blacklist);
		voter=sketches==null ? null : sketches.newSketch();
		votes=sketches==null ? null : new NcrnaDdlLayoutVotes(bins);
		prefix=sketches==null ? null : new NcrnaDdlSketchBuilder.Prefix(bins);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Clustering         ----------------*/
	/*--------------------------------------------------------------*/

	/** Input order is the existing builder's stable decreasing-length order.
	 * Fixed seeds never acquire invented members; unmatched records stay explicit. */
	public ArrayList<ArrayList<Read>> cluster(final ArrayList<byte[]> seeds){
		if(started){throw new IllegalStateException("A clusterer consumes its original input once");}
		started=true;
		final boolean fixed=seeds!=null && !seeds.isEmpty();
		membershipIdentity=fixed ? recruitIdentity : initialIdentity;
		if(fixed){for(final byte[] bases:seeds){families.add(new Family(families.size(), bases));}}
		Snapshot snapshot=snapshot();
		final ProcessThread greedy=new ProcessThread(this, snapshot, 0, 0);
		progress("admission", 0, 0, 0);
		for(int ordinal=0; ordinal<input.size(); ordinal++){
			final Read read=input.get(ordinal);
			int id=greedy.best(read, snapshot, fixed ? recruitIdentity : initialIdentity);
			if(id<0 && !fixed){
				id=families.size();
				families.add(new Family(id, read.bases));
			}
			assignments[ordinal]=id;
			if((ordinal+1)%1000==0 && ordinal+1<input.size()){
				progress("admission", 0, ordinal+1, greedy.matcher.alignments);
			}
			if(id<0){continue;}
			final Family family=families.get(id);
			family.members.add(ordinal);
			final int size=family.members.size;
			final boolean consensusDue=size>=family.nextConsensus;
			if(sketches!=null){
				sketches.fill(read.bases, voter);
				if(family.liveVotes.addAndSelect(voter.maxArray(), family.livePrefix, 100, family.profile.maxArray())){
					family.profile=DynamicDemiLog.fromArray(family.profile.maxArray(), family.id, null, 10, -1);
				}
				profileUpdates++;
			}
			if(consensusDue){
				refine(family, membershipIdentity);
				if(sketches!=null){
					sketches.fillPrefix(family.consensus, family.livePrefix);
					family.profile=DynamicDemiLog.fromArray(family.liveVotes.finish(family.livePrefix, 100, true).keys, family.id, null, 10, -1);
					profileRebuilds++;
				}
				while(family.nextConsensus<=size){family.nextConsensus=Math.addExact(Math.multiplyExact(family.nextConsensus-1, 2), 1);}
			}
			snapshot=publish(snapshot, family);
		}
		accumulate(greedy);
		progress("admission", 0, input.size(), 0);
		for(final Family family:families){family.liveVotes=null;family.livePrefix=null;}
		refreshAll(membershipIdentity);
		for(int round=0; round<rounds; round++){
			progress("reassignment", round+1, 0, 0);
			snapshot=snapshot();
			final int[] next=new int[input.size()];
			Arrays.fill(next, -1);
			final int count=Math.min(threads, input.size());
			final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(count);
			for(int worker=0; worker<count; worker++){
				final ProcessThread thread=new ProcessThread(this, snapshot,
					(int)((long)input.size()*worker/count), (int)((long)input.size()*(worker+1)/count));
				thread.output=next;
				workers.add(thread);
			}
			if(!ThreadWaiter.startAndWait(workers, this)){throw new IllegalStateException("A nucleotide reassignment worker failed");}
			progress("reassignment", round+1, input.size(), 0);
			int changed=0;
			for(final Family family:families){family.members.clear();}
			for(int ordinal=0; ordinal<next.length; ordinal++){
				if(next[ordinal]!=assignments[ordinal]){changed++;}
				if(next[ordinal]>=0){families.get(next[ordinal]).members.add(ordinal);}
			}
			assignments=next;
			membershipIdentity=recruitIdentity;
			refreshAll(recruitIdentity);
			if(changed==0){break;}
		}
		if(recruit){recruitOrphans();}
		final ArrayList<ArrayList<Read>> result=new ArrayList<ArrayList<Read>>(families.size());
		int accounted=0;
		for(final Family family:families){
			final ArrayList<Read> members=members(family);
			result.add(members);
			accounted+=members.size();
		}
		for(final int id:assignments){if(id<0){accounted++;}}
		assert(accounted==input.size()) : "Every input must occur once in a returned cluster or remain explicitly unassigned";
		return result;
	}

	public byte[] consensus(final int family){return families.get(family).consensus;}
	public float membershipIdentity(){return membershipIdentity;}

	/** Owner-only coarse progress. Admission's matcher is merged only at the
	 * phase end; include its live count without changing assignment or totals.
	 * Reassignment reports surround the joined workers, so counters are complete
	 * and no per-gene lock, atomic update or worker-output ordering is introduced. */
	private void progress(final String phase, final int round, final int processed, final long liveAlignments){
		assert(round>=0 && processed>=0 && processed<=input.size() && liveAlignments>=0)
			: "Construction progress must retain the full input denominator and report only completed model comparisons";
		System.err.println("Clustering "+phase+" round "+round+": genes processed "+processed+" of "+input.size()
			+"; model alignments "+(alignments+liveAlignments)+"; consensus alignments "+consensusAlignments);
	}

	/** Rechecks genuine orphans against the final rebuilt models. There is no
	 * additional k-mer veto: a new-policy no-pass result already exhausted models. */
	private void recruitOrphans(){
		final Snapshot snapshot=snapshot();
		final ProcessThread worker=new ProcessThread(this, snapshot, 0, 0);
		int recruited=0;
		for(int ordinal=0; ordinal<assignments.length; ordinal++){
			if(assignments[ordinal]>=0){continue;}
			final int family=worker.best(input.get(ordinal), snapshot, recruitIdentity);
			if(family>=0){
				assignments[ordinal]=family;
				families.get(family).members.add(ordinal);
				recruited++;
			}
		}
		accumulate(worker);
		if(recruited>0){
			membershipIdentity=Math.min(membershipIdentity, recruitIdentity);
			refreshAll(membershipIdentity);
		}
	}

	private void refreshAll(final float minimum){
		for(final Family family:families){refresh(family, true, minimum);}
	}

	private void refresh(final Family family, final boolean refine, final float minimum){
		if(family.members.isEmpty()){
			family.consensus=null;
			if(family.profile!=null){family.profile.clear();}
			return;
		}
		if(refine){refine(family, minimum);}
		if(sketches==null){return;}
		votes.reset();
		for(int i=0; i<family.members.size; i++){
			sketches.fill(input.get(family.members.get(i)).bases, voter);
			votes.add(voter.maxArray());
		}
		sketches.fillPrefix(family.consensus, prefix);
		final char[] selected=votes.finish(prefix, 100, true).keys;
		family.profile=DynamicDemiLog.fromArray(selected, family.id, null, 10, -1);
		while(family.nextConsensus<=family.members.size){family.nextConsensus=Math.addExact(Math.multiplyExact(family.nextConsensus-1, 2), 1);}
		profileRebuilds++;
	}

	private void refine(final Family family, final float minimum){
		final NcrnaDdlConsensus.Fit fit=NcrnaDdlConsensus.refine(members(family), family.consensus, minimum, passes, trim);
		family.consensus=fit.bases;
		consensusAlignments+=fit.alignments;
	}

	private ArrayList<Read> members(final Family family){
		final ArrayList<Read> reads=new ArrayList<Read>(family.members.size);
		for(int i=0; i<family.members.size; i++){reads.add(input.get(family.members.get(i)));}
		return reads;
	}

	private Snapshot snapshot(){
		if(sketches!=null && families.size()>65536){throw new IllegalStateException("Small DDL clustering index supports at most65536stable family slots");}
		final byte[][] models=new byte[families.size()][];
		final DynamicDemiLog[] profiles=sketches==null ? null : new DynamicDemiLog[families.size()];
		for(int id=0; id<families.size(); id++){
			final Family family=families.get(id);
			models[id]=family.consensus;
			if(profiles!=null){profiles[id]=family.profile;}
		}
		final DDLIndexSmall index=profiles==null || profiles.length==0 ? null : new DDLIndexSmall(profiles);
		if(index!=null){indexRebuilds++;}
		return new Snapshot(models, profiles, index);
	}

	/** Serial greedy publication overlays exact current rows on a frozen prefix.
	 * New families rebuild that prefix; at64changed IDs it compacts. Old posting
	 * scores are replaced, never added to a changed row's current score. */
	private Snapshot publish(final Snapshot current, final Family family){
		if(current.models.length!=families.size()){return snapshot();}
		current.models[family.id]=family.consensus;
		if(current.profiles!=null){
			current.profiles[family.id]=family.profile;
			if(!current.dirty[family.id]){
				current.dirty[family.id]=true;
				current.changed[current.changedCount++]=family.id;
			}
			if(current.changedCount>=64){return snapshot();}
		}
		return current;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Workers          ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void accumulate(final ProcessThread worker){
		if(worker.output!=null && !worker.success){errorState=true;}
		queries+=worker.matcher.queries;
		alignments+=worker.matcher.alignments;
		passing+=worker.matcher.passing;
		earlyStops+=worker.matcher.earlyStops;
	}
	@Override
	public boolean success(){return !errorState;}
	@Override
	public ReadWriteLock rwlock(){return null;}

	public static final class ProcessThread extends Thread{
		private ProcessThread(final NcrnaDdlClusterer owner_, final Snapshot snapshot_, final int from_, final int to_){
			owner=owner_;
			snapshot=snapshot_;
			from=from_;
			to=to_;
			matcher=new NcrnaClusterAssignment(owner.policy, owner.patience, owner.bins);
			sketches=owner.sketches==null ? null : new NcrnaDdlSketchBuilder(owner.bins, owner.seed, owner.blacklist);
			query=sketches==null ? null : sketches.newSketch();
			setDaemon(true);
		}

		@Override
		public void run(){
			for(int i=from; i<to; i++){output[i]=best(owner.input.get(i), snapshot, owner.recruitIdentity);}
			success=true;
		}

		private int best(final Read read, final Snapshot current, final float minimum){
			if(counts.length!=current.models.length){counts=new int[current.models.length];}
			if(current.index!=null){
				sketches.fill(read.bases, query);
				current.index.queryReplacing(query, counts, current.profiles, current.changed, current.changedCount);
			}
			return matcher.bestFamily(read.bases, current.models, counts, minimum, ALIGNMENT);
		}

		private final NcrnaDdlClusterer owner;
		private final Snapshot snapshot;
		private final int from, to;
		private final NcrnaClusterAssignment matcher;
		private final NcrnaDdlSketchBuilder sketches;
		private final DynamicDemiLog query;
		private int[] counts=new int[0], output;
		private boolean success;
	}

	private final class Family{
		private Family(final int id_, final byte[] bases){
			if(bases==null || bases.length==0){throw new IllegalArgumentException("A new family needs a real nonempty representative");}
			id=id_;
			consensus=bases;
			if(sketches!=null){
				profile=sketches.newSketch();
				sketches.fill(bases, profile);
				liveVotes=new NcrnaDdlLayoutVotes(bins);
				livePrefix=new NcrnaDdlSketchBuilder.Prefix(bins);
				sketches.fillPrefix(bases, livePrefix);
			}
		}
		private final int id;
		private final IntList members=new IntList();
		private DynamicDemiLog profile;
		private NcrnaDdlLayoutVotes liveVotes;
		private NcrnaDdlSketchBuilder.Prefix livePrefix;
		private byte[] consensus;
		private int nextConsensus=5;
	}

	private static final class Snapshot{
		private Snapshot(final byte[][] models_, final DynamicDemiLog[] profiles_, final DDLIndexSmall index_){
			models=models_;
			profiles=profiles_;
			index=index_;
			changed=new int[profiles==null ? 0 : models.length];
			dirty=new boolean[changed.length];
		}
		private final byte[][] models;
		private final DynamicDemiLog[] profiles;
		private final DDLIndexSmall index;
		private final int[] changed;
		private final boolean[] dirty;
		private int changedCount;
	}

	private static boolean validIdentity(final float value){return Float.isFinite(value) && value>=0 && value<=1;}
	private static final NcrnaClusterAssignment.Alignment ALIGNMENT=(query, reference)->QuantumAligner.alignStatic(query, reference, null);
	public long queries, alignments, consensusAlignments, passing, earlyStops, profileRebuilds, profileUpdates, indexRebuilds;
	private final ArrayList<Read> input;
	private final ArrayList<Family> families=new ArrayList<Family>();
	private final NcrnaClusterAssignment.Policy policy;
	private final LongHashSet blacklist;
	private final int patience, bins, rounds, passes, threads;
	private final long seed;
	private final float initialIdentity, recruitIdentity;
	private float membershipIdentity;
	private final boolean trim, recruit;
	private final NcrnaDdlSketchBuilder sketches;
	private final DynamicDemiLog voter;
	private final NcrnaDdlLayoutVotes votes;
	private final NcrnaDdlSketchBuilder.Prefix prefix;
	private int[] assignments;
	private boolean started, errorState;
}
