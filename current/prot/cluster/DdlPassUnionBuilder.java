package prot.cluster;

import java.util.Arrays;
import java.util.Comparator;

import ddl.FamilySketchBuilder;
import prot.ProteinSequence;
import structures.IntList;

/** Builds post-pass DDL merge unions in memory without mutating engine state. */
public final class DdlPassUnionBuilder{

	public DdlPassUnionBuilder(ProteinSequence[] sequences_,FamilySketchBuilder sketches_,
			TailAlignmentGate gate_,ConsensusRefiner refiner_,ClusteringMetrics metrics_){
		if(sequences_==null || sketches_==null || gate_==null || refiner_==null || metrics_==null){
			throw new IllegalArgumentException("Pass unions need sequences, sketch scratch, gate, refiner and metrics");
		}
		if(refiner_.rounds()!=1){throw new IllegalArgumentException("Pass unions need a one-round worker refiner");}
		sequences=sequences_; sketches=sketches_; gate=gate_; refiner=refiner_; metrics=metrics_;
	}

	/** Collects live family consensuses in stable family-ID order for the nested model run. */
	public static ModelInput collectModels(Family[] familiesById){
		if(familiesById==null){throw new IllegalArgumentException("Missing family table");}
		final IntList ids=new IntList();
		for(int id=0; id<familiesById.length; id++){
			final Family family=familiesById[id];
			if(family==null){continue;}
			if(family.familyId!=id || family.size()<1){
				throw new IllegalArgumentException("Live family table must be addressed by nonempty stable family ID");
			}
			ids.add(id);
		}
		final int[] familyIds=ids.toArray();
		final ProteinSequence[] models=new ProteinSequence[familyIds.length];
		for(int i=0; i<familyIds.length; i++){
			final Family family=familiesById[familyIds[i]];
			models[i]=new ProteinSequence("family_"+family.familyId, ClusteringInput.decode(family.consensusUnsafe()));
		}
		return new ModelInput(models, familyIds);
	}

	/** Converts a nested consensus-clustering result into multi-family merge groups only. */
	public static Group[] mergeGroups(ModelInput input,GreedyClusterer.Result modelResult){
		if(input==null || modelResult==null){throw new IllegalArgumentException("Missing model input or result");}
		if(modelResult.inputs()!=input.size() || modelResult.eligible!=input.size() ||
				modelResult.assigned!=input.size()){
			throw new IllegalArgumentException("Nested model clustering must retain every collected family");
		}
		final int[] assignments=new int[input.size()];
		for(int i=0; i<assignments.length; i++){assignments[i]=modelResult.familyFor(i);}
		return mergeGroups(input, assignments, modelResult.families);
	}

	/** Converts dense model assignments into deterministic multi-family merge groups. */
	public static Group[] mergeGroups(ModelInput input,int[] modelFamilyByOrdinal,Family[] modelSeeds){
		if(input==null || modelFamilyByOrdinal==null || modelSeeds==null ||
				modelFamilyByOrdinal.length!=input.size()){
			throw new IllegalArgumentException("Model assignments must match the collected families");
		}
		int max=-1;
		for(Family seed : modelSeeds){
			if(seed==null){throw new IllegalArgumentException("Model seed list contains null");}
			if(seed.familyId>=input.size()){
				throw new IllegalArgumentException("Nested model family ID exceeds the dense model ordinal range");
			}
			if(seed.size()<1){throw new IllegalArgumentException("Nested model family contains no consensuses");}
			max=Math.max(max, seed.familyId);
		}
		final Family[] seedById=new Family[max+1];
		final IntList[] sourceIds=new IntList[max+1];
		for(Family seed : modelSeeds){
			if(seedById[seed.familyId]!=null){throw new IllegalArgumentException("Duplicate model family ID "+seed.familyId);}
			seedById[seed.familyId]=seed; sourceIds[seed.familyId]=new IntList();
		}
		for(int ordinal=0; ordinal<modelFamilyByOrdinal.length; ordinal++){
			final int modelFamilyId=modelFamilyByOrdinal[ordinal];
			if(modelFamilyId<0 || modelFamilyId>=seedById.length || seedById[modelFamilyId]==null){
				throw new IllegalArgumentException("Model assignment references missing family "+modelFamilyId);
			}
			sourceIds[modelFamilyId].add(input.familyIds[ordinal]);
		}
		for(Family seed : modelSeeds){
			if(sourceIds[seed.familyId].size!=seed.size()){
				throw new IllegalArgumentException("Nested model membership count differs from the dense assignments");
			}
		}
		int groups=0;
		for(IntList list : sourceIds){if(list!=null && list.size>1){groups++;}}
		final Group[] out=new Group[groups];
		for(int id=0, next=0; id<sourceIds.length; id++){
			final IntList list=sourceIds[id];
			if(list!=null && list.size>1){out[next++]=group(seedById[id], list.toArray());}
		}
		Arrays.sort(out, SURVIVOR_ORDER);
		return out;
	}

	/** Defines one replacement union: one model seed and two or more source families. */
	public static Group group(Family seed,int[] sourceFamilyIds){
		return new Group(seed, sourceFamilyIds);
	}

	/** Builds one replacement without editing the caller's family table, index or assignments. */
	public Replacement build(Group group,Family[] familiesById,boolean captureQuality){
		if(group==null || familiesById==null){throw new IllegalArgumentException("Missing union group or family table");}
		final IntList members=new IntList();
		for(int id : group.sourceFamilyIds){
			if(id<0 || id>=familiesById.length || familiesById[id]==null || familiesById[id].familyId!=id){
				throw new IllegalArgumentException("Union references missing source family "+id);
			}
			final Family source=familiesById[id];
			if(source.size()<1){throw new IllegalArgumentException("Union references empty source family "+id);}
			for(int i=0; i<source.size(); i++){members.add(source.memberIndex(i));}
		}
		members.sort();
		checkMembers(members);
		final int survivorId=group.survivorId();
		final Family survivor=familiesById[survivorId];
		final Family replacement=new Family(survivorId,survivor.id,members.array[0],sequences[members.array[0]].enc);
		replacement.setConsensus(group.seed.consensusUnsafe());
		for(int i=1; i<members.size; i++){replacement.addMember(members.array[i]);}
		final ConsensusRefiner.Update update=refiner.updateAfterConcurrentGreedy(replacement,
			sequences, gate, metrics, captureQuality);
		replacement.refreshStats(sequences, sketches);
		replacement.markSketchRebuilt();
		return new Replacement(replacement, group.sourceFamilyIds, status(update), reason(update), update);
	}

	private void checkMembers(IntList members){
		if(members.size<1){throw new IllegalStateException("Pass union did not retain any source members");}
		int previous=-1;
		for(int i=0; i<members.size; i++){
			final int member=members.array[i];
			if(member<=previous || member>=sequences.length || sequences[member]==null){
				throw new IllegalStateException("Union contains an invalid or duplicate member ordinal: "+member);
			}
			previous=member;
		}
	}

	private static Status status(ConsensusRefiner.Update update){
		return update.accepted()>0 ? Status.UNION_ACCEPTED :
			update.unchanged()>0 ? Status.UNION_UNCHANGED : Status.UNION_REJECTED;
	}

	private static String reason(ConsensusRefiner.Update update){
		return update.accepted()>0 ? "all_members_pass" :
			update.unchanged()>0 ? "same_sequence" :
			update.incompatibleSnapshots()>0 ? "snapshot_no_positive_placement" :
			update.rejected()>0 ? "strict_new_validation_failed" : "not_run";
	}

	/** Dense consensus inputs plus the stable family ID represented by each ordinal. */
	public static final class ModelInput{
		private ModelInput(ProteinSequence[] sequences_,int[] familyIds_){
			sequences=sequences_; familyIds=familyIds_;
		}
		public int size(){return sequences.length;}
		public final ProteinSequence[] sequences;
		public final int[] familyIds;
	}

	public static final class Group{
		private Group(Family seed_,int[] sourceFamilyIds_){
			if(seed_==null || sourceFamilyIds_==null || sourceFamilyIds_.length<2){
				throw new IllegalArgumentException("A pass union needs one model seed and at least two source families");
			}
			sourceFamilyIds=sourceFamilyIds_.clone();
			Arrays.sort(sourceFamilyIds);
			for(int i=0; i<sourceFamilyIds.length; i++){
				if(sourceFamilyIds[i]<0 || (i>0 && sourceFamilyIds[i]==sourceFamilyIds[i-1])){
					throw new IllegalArgumentException("Source family IDs must be unique and nonnegative");
				}
			}
			seed=seed_;
		}
		public int survivorId(){return sourceFamilyIds[0];}
		public final Family seed;
		public final int[] sourceFamilyIds;
	}

	public static final class Replacement{
		private Replacement(Family family_,int[] sourceFamilyIds_,Status status_,String reason_,
				ConsensusRefiner.Update update_){
			family=family_; sourceFamilyIds=sourceFamilyIds_.clone();
			status=status_; reason=reason_; update=update_;
			quality=update.quality();
		}
		public int survivorId(){return family.familyId;}
		public final Family family;
		public final int[] sourceFamilyIds;
		public final Status status;
		public final String reason;
		public final ConsensusRefiner.Update update;
		public final ConsensusQualityReport.Snapshot quality;
	}

	public enum Status {UNION_ACCEPTED, UNION_UNCHANGED, UNION_REJECTED}

	private final ProteinSequence[] sequences;
	private final FamilySketchBuilder sketches;
	private final TailAlignmentGate gate;
	private final ConsensusRefiner refiner;
	private final ClusteringMetrics metrics;

	private static final Comparator<Group> SURVIVOR_ORDER=new Comparator<Group>(){
		@Override
		public int compare(Group a,Group b){return Integer.compare(a.survivorId(), b.survivorId());}
	};
}
