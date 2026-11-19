package prot.cluster;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteStreamWriter;
import prot.ProteinSequence;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Read-only coverage, size and member-to-consensus quality for consecutive live states.
 * Callers must treat passed Family objects and input sequences as immutable. The
 * live state publisher keeps unchanged families by identity and constructs new
 * objects for every changed model/membership, allowing exact quality reuse.
 * @author Brian Bushnell, Keqing
 */
public final class LiveProbeMetrics{

	public LiveProbeMetrics(final ProteinSequence[] sequences_,final TailAlignmentGate gate_,
			final int threads_,final int budget_){
		if(sequences_==null || gate_==null || threads_<1 || budget_<1){
			throw new IllegalArgumentException("Live metrics need immutable inputs, an alignment gate, positive threads and budget");
		}
		sequences=sequences_;gate=gate_;threads=threads_;budget=budget_;
	}

	/** Measures only new Family objects; same IDs with new objects are always rescored.
	 * Identity and raw BLOSUM score-ratio means use the existing descriptive scorer,
	 * which ignores acceptance thresholds. Missing placements remain explicit.
	 */
	public Snapshot measure(final Family[] models){
		if(models==null){throw new IllegalArgumentException("Live metrics require a complete model roster");}
		final long start=System.nanoTime();
		final Family[] ordered=models.clone();
		for(Family model:ordered){
			if(model==null || model.familyId<0 || model.size()<0){throw new IllegalArgumentException("Invalid live metric model");}
		}
		Arrays.sort(ordered,(a,b)->Integer.compare(a.familyId,b.familyId));
		final int[] sizes=new int[ordered.length];
		final ArrayList<Family> activeList=new ArrayList<Family>();
		long covered=0;
		for(int i=0;i<ordered.length;i++){
			final Family model=ordered[i];
			if(i>0 && ordered[i-1].familyId==model.familyId){throw new IllegalArgumentException("Duplicate metric family ID: "+model.familyId);}
			sizes[i]=model.size();covered=Math.addExact(covered,model.size());
			if(model.size()>0){activeList.add(model);}
		}
		if(covered>sequences.length){throw new IllegalArgumentException("Metric memberships exceed the original input population");}
		final Family[] active=activeList.toArray(new Family[activeList.size()]);
		final ConsensusQualityReport.Snapshot[] quality=new ConsensusQualityReport.Snapshot[active.length];
		final ArrayList<Family> changed=new ArrayList<Family>();
		int old=0,reused=0;
		long changedMembers=0;
		for(int i=0;i<active.length;i++){
			final Family family=active[i];
			while(old<previous.length && previous[old].familyId<family.familyId){old++;}
			if(old<previous.length && previous[old]==family){
				quality[i]=previousQuality[old];reused++;
				if(quality[i].members()!=family.size()){
					throw new IllegalStateException("A previously measured Family changed membership in place: "+family.familyId);
				}
			}else{changed.add(family);changedMembers=Math.addExact(changedMembers,family.size());}
		}
		if(!changed.isEmpty()){
			final ConsensusQualityReport measured=SavedConsensusQuality.measure(sequences,
				changed.toArray(new Family[changed.size()]),gate,threads);
			for(int i=0;i<quality.length;i++){
				if(quality[i]==null){quality[i]=measured.snapshot(active[i].familyId);}
			}
		}
		final IntHashMap index=new IntHashMap(active.length*2+1);
		final QualityTotals all=new QualityTotals(),top=new QualityTotals();
		for(int i=0;i<active.length;i++){index.put(active[i].familyId,i);all.add(quality[i]);}
		final Family[] ranked=active.clone();
		Arrays.sort(ranked,(a,b)->{
			final int cmp=Integer.compare(b.size(),a.size());
			return cmp==0 ? Integer.compare(a.familyId,b.familyId) : cmp;
		});
		long topMembers=0;
		for(int i=0;i<Math.min(budget,ranked.length);i++){
			final Family family=ranked[i];
			topMembers=Math.addExact(topMembers,family.size());top.add(quality[index.get(family.familyId)]);
		}
		final int slot=ranked.length>=budget ? ranked[budget-1].size() : 0;
		if(all.members+all.missing!=covered || top.members+top.missing!=topMembers){
			throw new IllegalStateException("Quality availability must partition the full and top-budget member censuses");
		}
		Arrays.sort(sizes);
		previous=active;previousQuality=quality;
		return new Snapshot(sizes,active.length,covered,topMembers,slot,all,top,
			changed.size(),reused,changedMembers,System.nanoTime()-start);
	}

	private static final class QualityTotals{
		void add(final ConsensusQualityReport.Snapshot snapshot){
			assert(snapshot!=null && snapshot.members()>0) : "Each active family must have one descriptive quality measurement";
			if(snapshot.available()){
				members=Math.addExact(members,snapshot.members());
				identity+=snapshot.members()*snapshot.identityMean();
				ratio+=snapshot.members()*snapshot.scoreRatioMean();
				if(!Double.isFinite(identity) || !Double.isFinite(ratio)){throw new IllegalStateException("Quality sums must remain finite");}
			}else{missing=Math.addExact(missing,snapshot.members());}
		}
		void append(final ByteBuilder line){
			line.append(members).tab().append(missing).tab();
			if(members==0){line.append("NA\tNA");}
			else{line.appendSlow(identity/members).tab().appendSlow(ratio/members);}
		}
		long members,missing;
		double identity,ratio;
	}

	/** Counts describe current memberships; top-budget coverage is a census, not a reassignment experiment. */
	public static final class Snapshot{
		private Snapshot(int[] sizes_,int active_,long covered_,long topMembers_,int slot_,
				QualityTotals all_,QualityTotals top_,int measured_,int reused_,long changedMembers_,long elapsed_){
			sizes=sizes_;active=active_;covered=covered_;topMembers=topMembers_;slot=slot_;
			all=all_;top=top_;measured=measured_;reused=reused_;changedMembers=changedMembers_;elapsed=elapsed_;
			assert(measured+reused==active && topMembers<=covered && elapsed>=0) :
				"Metrics cover every active family once and a subset of current memberships";
		}
		public void writeSummary(final ByteStreamWriter writer,final String stage){
			validate(writer,stage);
			final ByteBuilder line=new ByteBuilder();
			line.append(stage).tab().append(sizes.length).tab().append(active).tab().append(covered).tab()
				.append(topMembers).tab().append(slot).tab();
			all.append(line);line.tab();top.append(line);
			line.tab().append(measured).tab().append(reused).tab().append(changedMembers).tab().append(elapsed);
			writer.println(line);
		}
		/** Exact size-frequency rows include zero-member payloads, with zero contributed tails. */
		public void writeHistogram(final ByteStreamWriter writer,final String stage){
			validate(writer,stage);
			final ByteBuilder line=new ByteBuilder();
			for(int i=0;i<sizes.length;){
				int end=i+1;while(end<sizes.length && sizes[end]==sizes[i]){end++;}
				line.clear().append(stage).tab().append(sizes[i]).tab().append(end-i).tab().append((long)sizes[i]*(end-i));
				writer.println(line);i=end;
			}
		}
		public int measuredFamilies(){return measured;}
		public int reusedFamilies(){return reused;}
		public long covered(){return covered;}
		public long topMembers(){return topMembers;}
		public int slotCost(){return slot;}
		public long qualityMembers(){return all.members;}
		public long missingQualityMembers(){return all.missing;}
		public long elapsedNanos(){return elapsed;}
		private static void validate(final ByteStreamWriter writer,final String stage){
			if(writer==null || stage==null || stage.isEmpty()){throw new IllegalArgumentException("Metrics output needs a writer and stage token");}
			for(int i=0;i<stage.length();i++){
				if(stage.charAt(i)<=32 || stage.charAt(i)>=127){throw new IllegalArgumentException("Metric stage must be one ASCII token");}
			}
		}
		private final int[] sizes;
		private final int active,slot,measured,reused;
		private final long covered,topMembers,changedMembers,elapsed;
		private final QualityTotals all,top;
	}

	public static final String SUMMARY_HEADER="stage\tmodels\tactive_models\tcovered_members\ttop_budget_members\tslot_cost\tquality_members\tquality_missing_members\tidentity_mean\tscore_ratio_mean\ttop_quality_members\ttop_quality_missing_members\ttop_identity_mean\ttop_score_ratio_mean\trescored_families\treused_families\trescored_member_population\telapsed_ns";
	public static final String HISTOGRAM_HEADER="stage\tmembers\tfamilies\ttails";
	private final ProteinSequence[] sequences;
	private final TailAlignmentGate gate;
	private final int threads,budget;
	private Family[] previous=new Family[0];
	private ConsensusQualityReport.Snapshot[] previousQuality=new ConsensusQualityReport.Snapshot[0];
}
