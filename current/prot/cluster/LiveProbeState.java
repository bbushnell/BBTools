package prot.cluster;

import java.util.Arrays;
import java.util.HashMap;

import prot.ProteinSequence;
import structures.IntHashMap;

/** Frozen model roster and exact covered/residual partition between live-probe passes.
 * Arrays and models returned by accessors are BORROWED READ-ONLY. Baseline and
 * unaffected models are shared read-only; published replacements own fresh
 * membership and model arrays. Evaluated model payloads are not recomputed at
 * publication; optional ranking descriptors describe the actual winning members.
 * @author Brian Bushnell, Fischl
 */
public final class LiveProbeState{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Takes a read-only baseline lease and reconstructs one owner per eligible ordinal. */
	public static LiveProbeState fromBaseline(final LocalProbeBaseline baseline, final int sequenceCount){
		require(baseline!=null && baseline.sourceHeader!=null && baseline.selected!=null &&
			baseline.residualOrdinals!=null && sequenceCount>=0 && baseline.firstNewId>=0,
			"Live state needs a complete frozen baseline and original sequence count");
		final Family[] models=baseline.selected.clone();
		final IntHashMap index=index(models, baseline.sourceHeader.width);
		final int[] owners=new int[sequenceCount]; Arrays.fill(owners, FILTERED);
		int assigned=0;
		for(Family model:models){
			require(model.familyId<baseline.firstNewId, "Fresh family IDs must exceed all baseline selected IDs");
			for(int i=0; i<model.size(); i++){
				final int ordinal=model.memberIndex(i);
				require(validOrdinal(ordinal, owners) && owners[ordinal]==FILTERED,
					"Baseline family memberships must be unique original input ordinals");
				owners[ordinal]=model.familyId; assigned++;
			}
		}
		final int[] residual=baseline.residualOrdinals.clone();
		for(int i=0; i<residual.length; i++){
			final int ordinal=residual[i];
			require(validOrdinal(ordinal, owners) && owners[ordinal]==FILTERED && (i==0 || ordinal>residual[i-1]),
				"Baseline residuals must be sorted unique eligible ordinals outside every assigned family");
			owners[ordinal]=RESIDUAL;
		}
		final int filtered=sequenceCount-assigned-residual.length;
		require(assigned==baseline.assigned && filtered==baseline.filtered,
			"Baseline assigned, residual and filtered counts must partition the original input");
		return new LiveProbeState(models, owners, residual, index, baseline.firstNewId,
			assigned, filtered, baseline.sourceHeader, 0, new boolean[models.length],null,0);
	}

	/** Binds publication to the original population after its family descriptors are hydrated.
	 * This returns a new state envelope without modifying borrowed models. Rebinding at a
	 * later layer must retain the same cache; it never attaches a second family accumulator.
	 */
	public LiveProbeState withMemberRanking(final ProteinSequence[] population,
			final MemberRankInputs inputs,final int layer){
		require(population!=null && inputs!=null && inputs.belongsTo(population) &&
			population.length==owners.length && layer>=0 &&
			(memberInputs==null || memberInputs==inputs && layer>=memberLayer),
			"Ranked publication requires the identical original population/cache and a nondecreasing measurement layer");
		final MemberRankStatistics.View view=new MemberRankStatistics.View();
		for(Family model:models){
			require(model.memberRankingUses(inputs),"Every frozen model must use the same original ordinal cache before ranked publication");
			model.readMemberRanking(view);
			require(view.layer<=layer,"Publication layer cannot precede an existing member measurement");
		}
		return new LiveProbeState(models,owners,residual,modelIndex,firstNewId,assigned,filtered,
			sourceHeader,lastMarginalGain,fusionProducts,inputs,layer);
	}

	private LiveProbeState(final Family[] models_, final int[] owners_, final int[] residual_,
			final IntHashMap index_, final int firstNewId_, final int assigned_, final int filtered_,
			final ClusterSketchCodec.Header sourceHeader_, final long delta, final boolean[] fusionProducts_,
			final MemberRankInputs inputs_,final int layer_){
		models=models_; owners=owners_; residual=residual_; modelIndex=index_; fusionProducts=fusionProducts_;
		firstNewId=firstNewId_; assigned=assigned_; filtered=filtered_; sourceHeader=sourceHeader_;
		lastMarginalGain=delta;
		memberInputs=inputs_;memberLayer=layer_;
		int count=0;
		for(Family model:models){if(model.size()>0){count++;}}
		active=new Family[count]; count=0;
		for(Family model:models){if(model.size()>0){active[count++]=model;}}
		require(fusionProducts!=null && fusionProducts.length==models.length,
			"Fusion-product tags must be parallel to the frozen model roster");
		assert((long)assigned+filtered+residual.length==owners.length) : "Every original ordinal has exactly one coverage status";
	}

	/*--------------------------------------------------------------*/
	/*----------------         Publication          ----------------*/
	/*--------------------------------------------------------------*/

	/** Publishes evaluated winners in accepted order, validating their exact marginal gains.
	 * Unaffected models keep roster order; replacements follow in accepted order,
	 * first then second. Membership in each replacement is ascending input ordinal.
	 * nextUnusedFamilyId reserves IDs allocated while making this pass's proposals,
	 * including rejected proposals. A failure cannot mutate this state or proposals.
	 * Reused parent IDs retain their original names. A fresh model's colliding name
	 * receives a deterministic metadata-only suffix; its numeric ID and arrays stay exact.
	 */
	public LiveProbeState commit(final Proposal[] proposals, final LiveProbeSelection.Result selection,
			final int nextUnusedFamilyId){
		require(proposals!=null && selection!=null && nextUnusedFamilyId>=firstNewId,
			"Publication needs evaluated proposals, selection and a nondecreasing fresh-ID cursor");
		for(Family model:models){
			require(memberInputs==null ? !model.hasMemberRanking() : model.memberRankingUses(memberInputs),
				"Publication cannot mix ranked and legacy models or discard an unbound ranking cache");
		}
		final Proposal[] input=proposals.clone();
		final IntHashMap requests=new IntHashMap(mapCapacity(input.length));
		for(int i=0; i<input.length; i++){
			final Proposal proposal=input[i];
			require(proposal!=null && requests.get(proposal.requestId)<0, "Evaluated request identities must be unique and present");
			validate(proposal, nextUnusedFamilyId);
			requests.put(proposal.requestId, i);
		}
		require(selection.size()<=input.length && (long)selection.size()*2<=Integer.MAX_VALUE-8,
			"Selection cannot publish more requests than supplied or overflow the replacement roster");
		final boolean[] consumed=new boolean[models.length], selected=new boolean[input.length];
		final Family[] replacements=new Family[selection.size()*2];
		final boolean[] replacementFusionProducts=new boolean[replacements.length];
		final IntHashMap replacementIndex=new IntHashMap(mapCapacity(replacements.length));
		final HashMap<String,Family> names=names(models);
		final int[] nextOwners=owners.clone();
		int replacementCount=0, removed=0;
		long delta=0;
		for(int step=0; step<selection.size(); step++){
			final int at=requests.get(selection.requestId(step));
			require(at>=0 && !selected[at], "Accepted requests must occur exactly once among the evaluated proposals");
			selected[at]=true;
			final Proposal proposal=input[at];
			for(int parent:proposal.parents){
				final int row=modelIndex.get(parent);
				require(!consumed[row], "Accepted operations cannot share an input parent");
				consumed[row]=true; removed++;
				final Family old=models[row];
				for(int i=0; i<old.size(); i++){nextOwners[old.memberIndex(i)]=RESIDUAL;}
			}
			final int firstAt=replacementCount;
			replacementCount=appendCopy(proposal.first, replacements, replacementCount, replacementIndex, names);
			final int secondAt=replacementCount;
			if(proposal.second!=null){replacementCount=appendCopy(proposal.second, replacements, replacementCount, replacementIndex, names);}
			else{replacementFusionProducts[firstAt]=true;}
			final LiveProbeAssignments.Snapshot assignments=proposal.assignments;
			int duplicateClaims=0;
			for(int r=0, p=0; r<assignments.retainedCount() || p<assignments.rescuedCount();){
				final boolean affected=p>=assignments.rescuedCount() ||
					(r<assignments.retainedCount() && assignments.retainedOrdinal(r)<assignments.rescuedOrdinal(p));
				if(affected){
					final int ordinal=assignments.retainedOrdinal(r), winner=assignments.retainedWinner(r);
					nextOwners[ordinal]=winner;
					addReplacementMember(proposal, replacements, firstAt, secondAt, winner, ordinal);
					r++;
				}else{
					final int ordinal=assignments.rescuedOrdinal(p);
					if(nextOwners[ordinal]==RESIDUAL){
						final int winner=assignments.rescuedWinner(p);
						nextOwners[ordinal]=winner;
						addReplacementMember(proposal, replacements, firstAt, secondAt, winner, ordinal);
					}else{duplicateClaims++;}
					p++;
				}
			}
			final long marginal=Math.subtractExact(proposal.measuredNetGain, duplicateClaims);
			require(marginal==selection.marginalGain(step), "Selected marginal gain differs from actual first-owner rescue accounting");
			delta=Math.addExact(delta, marginal);
		}
		final long rosterSize=(long)models.length-removed+replacementCount;
		require(rosterSize<=Integer.MAX_VALUE-8, "Published roster exceeds primitive array capacity");
		final Family[] nextModels=new Family[(int)rosterSize];
		final boolean[] nextFusionProducts=new boolean[nextModels.length];
		final IntHashMap nextIndex=new IntHashMap(mapCapacity(nextModels.length));
		int out=0;
		for(int i=0; i<models.length; i++){
			if(!consumed[i]){
				final Family model=models[i];
				nextModels[out]=model;
				nextFusionProducts[out]=fusionProducts[i];
				nextIndex.put(model.familyId, out);
				out++;
			}
		}
		for(int i=0; i<replacementCount; i++){
			final Family replacement=replacements[i];
			require(nextIndex.get(replacement.familyId)<0, "Published replacement family IDs must not collide with retained models");
			nextModels[out]=replacement;
			nextFusionProducts[out]=replacementFusionProducts[i];
			nextIndex.put(replacement.familyId, out);
			out++;
		}
		assert(out==nextModels.length) : "Published roster construction must fill the exact new model count";
		final long expectedAssigned=Math.addExact(assigned, delta);
		require(expectedAssigned>=0 && expectedAssigned<=owners.length-filtered,
			"Sum of selected marginal gains must leave a valid assigned count");
		final int[] nextResidual=new int[owners.length-filtered-(int)expectedAssigned];
		int countedAssigned=0, countedResidual=0;
		for(int ordinal=0; ordinal<nextOwners.length; ordinal++){
			final int owner=nextOwners[ordinal];
			require((owners[ordinal]==FILTERED)==(owner==FILTERED), "Publication cannot change the eligible input population");
			if(owner==FILTERED){continue;}
			if(owner==RESIDUAL){
				require(countedResidual<nextResidual.length, "Actual losses exceed the selected marginal-gain ledger");
				nextResidual[countedResidual++]=ordinal;
			}else{
				require(owner>=0, "An assigned ordinal must name a nonnegative model ID");
				countedAssigned++;
			}
		}
		require(countedAssigned==expectedAssigned && countedResidual==nextResidual.length,
			"Actual assigned-count delta must equal the sum of selected marginal gains");
		long memberCount=0;
		for(Family model:nextModels){memberCount+=model.size();}
		require(memberCount==countedAssigned, "Published memberships must cover exactly the assigned ordinal ledger");
		return new LiveProbeState(nextModels, nextOwners, nextResidual, nextIndex,
			nextUnusedFamilyId, countedAssigned, filtered, sourceHeader, delta, nextFusionProducts,memberInputs,memberLayer);
	}

	/** Checks every evaluated record against the original state, not partially committed state. */
	private void validate(final Proposal proposal, final int nextUnused){
		long affected=0;
		for(int parent:proposal.parents){
			final int at=modelIndex.get(parent);
			require(at>=0 && models[at].size()>0, "Probe parents must be active models in this frozen state");
			affected+=models[at].size();
		}
		validateReplacement(proposal.first, proposal.parents, nextUnused);
		if(proposal.second!=null){validateReplacement(proposal.second, proposal.parents, nextUnused);}
		final LiveProbeAssignments.Snapshot rows=proposal.assignments;
		for(int i=0; i<rows.retainedCount(); i++){
			final int ordinal=rows.retainedOrdinal(i);
			require(validOrdinal(ordinal, owners) && Arrays.binarySearch(proposal.parents, owners[ordinal])>=0,
				"Retained ordinals must belong to this proposal's original affected parents");
			require(proposal.hasWinner(rows.retainedWinner(i)), "Retained winner must name an evaluated replacement");
		}
		for(int i=0; i<rows.rescuedCount(); i++){
			final int ordinal=rows.rescuedOrdinal(i);
			require(validOrdinal(ordinal, owners) && owners[ordinal]==RESIDUAL,
				"Rescued ordinals must come from the frozen state's residual pool");
			require(proposal.hasWinner(rows.rescuedWinner(i)), "Rescued winner must name an evaluated replacement");
		}
		final long measured=(long)rows.rescuedCount()+rows.retainedCount()-affected;
		require(measured==proposal.measuredNetGain, "Measured net gain must equal rescued plus retained minus ALL original affected members");
	}

	/** Existing IDs belong only to this operation's parents; fresh IDs lie in the reserved interval. */
	private void validateReplacement(final Family model, final int[] parents, final int nextUnused){
		validateModel(model, sourceHeader.width);
		require(memberInputs==null ? !model.hasMemberRanking() : model.memberRankingUses(memberInputs),
			"Evaluated replacements must retain the frozen state's member-ranking mode and empty-template provenance");
		require(model.representativeIndex>=0 && model.representativeIndex<owners.length,
			"Evaluated replacement needs its real original representative ordinal for an exact detached copy");
		require(Arrays.binarySearch(parents, model.familyId)>=0 ||
			model.familyId>=firstNewId && model.familyId<nextUnused,
			"Replacement IDs must reuse an affected parent or belong to this pass's reserved fresh-ID interval");
	}

	/** Copies payload bytes through defensive setters, clearing ONLY the new object's initial member. */
	private int appendCopy(final Family source, final Family[] copies, final int at,
			final IntHashMap index, final HashMap<String,Family> names){
		require(index.get(source.familyId)<0, "Selected replacement family IDs must be globally unique");
		final String name=publicationName(source, names);
		final Family copy=new Family(source.familyId, name, source.representativeIndex, source.representative());
		copy.clearMembers();
		if(memberInputs!=null){
			final MemberRankStatistics.Snapshot measured=source.memberRankingSnapshot();
			require(measured.layer<=memberLayer,"Publication layer cannot precede the evaluated proposal's member measurement");
			copy.enableMemberRanking(memberInputs,memberLayer,measured);
		}
		copy.setConsensus(source.consensusUnsafe()); copy.setSketch(source.sketchUnsafe());
		copy.setMonomers(source.monomersUnsafe()); copy.setScalars(source.scalarsUnsafe(), source.scalarMask());
		copies[at]=copy; index.put(source.familyId, at);
		names.put(name, copy);
		return at+1;
	}

	/** Appends a retained or first-claim rescue to the copied winning child. */
	private static void addReplacementMember(final Proposal proposal, final Family[] replacements,
			final int firstAt, final int secondAt, final int winner, final int ordinal){
		final Family first=replacements[firstAt];
		if(winner==first.familyId){first.addMember(ordinal); return;}
		final Family second=replacements[secondAt];
		require(second!=null && winner==second.familyId,
			"Published assignment winner must name one of the copied replacements");
		second.addMember(ordinal);
	}

	/** Reserves original names even for consumed parents, then disambiguates fresh models only.
	 * The accepted order and stable numeric ID determine suffix ownership. A colliding
	 * fallback receives an increasing ordinal, so preexisting suffixes are safe too.
	 */
	private String publicationName(final Family source, final HashMap<String,Family> names){
		final int original=modelIndex.get(source.familyId);
		if(original>=0){return models[original].id;}
		if(!names.containsKey(source.id)){return source.id;}
		final String stem=source.id+".live"+source.familyId;
		String name=stem;
		for(long suffix=1; names.containsKey(name); suffix++){name=stem+"."+suffix;}
		return name;
	}

	/** A name maps to its model, not a boxed numeric ID. Codec export requires unique names. */
	private static HashMap<String,Family> names(final Family[] models){
		final HashMap<String,Family> names=new HashMap<String,Family>(mapCapacity(models.length));
		for(Family model:models){
			require(model!=null && names.put(model.id, model)==null, "Frozen roster names must be unique for codec export");
		}
		return names;
	}

	/** Builds a stable-ID lookup while rejecting duplicate IDs/names and incomplete payload geometry. */
	private static IntHashMap index(final Family[] models, final int width){
		names(models);
		final IntHashMap index=new IntHashMap(mapCapacity(models.length));
		for(int i=0; i<models.length; i++){
			validateModel(models[i], width);
			require(index.get(models[i].familyId)<0, "Frozen roster family IDs must be unique");
			index.put(models[i].familyId, i);
		}
		return index;
	}
	private static int mapCapacity(final int entries){
		return Math.max(4, entries>=(Integer.MAX_VALUE-8)/2 ? Integer.MAX_VALUE-8 : entries*2+1);
	}
	private static void validateModel(final Family model, final int width){
		require(model!=null && model.sketchUnsafe()!=null && model.sketchUnsafe().length==width,
			"Frozen model payload must have the source sketch geometry");
	}
	private static boolean validOrdinal(final int ordinal, final int[] owners){return ordinal>=0 && ordinal<owners.length;}

	/*--------------------------------------------------------------*/
	/*----------------       Evaluated proposal     ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable envelope; evaluated Family payloads are borrowed READ-ONLY until commit completes. */
	public static final class Proposal{
		public Proposal(final int requestId_, final int[] parentIds, final Family first_, final Family second_,
				final LiveProbeAssignments.Snapshot assignments_, final long measuredNetGain_){
			require(requestId_>=0 && parentIds!=null && parentIds.length>0 && first_!=null && assignments_!=null &&
				(second_==null || first_.familyId!=second_.familyId), "Proposal needs stable identity, parents and distinct evaluated replacements");
			parents=parentIds.clone();
			for(int i=0; i<parents.length; i++){
				require(parents[i]>=0 && (i==0 || parents[i]>parents[i-1]), "Proposal parent IDs must be sorted unique nonnegative IDs");
			}
			requestId=requestId_; first=first_; second=second_; assignments=assignments_; measuredNetGain=measuredNetGain_;
		}
		/** Exact immutable input to the deterministic overlap selector. */
		public LiveProbeSelection.Operation selectionOperation(){
			return new LiveProbeSelection.Operation(requestId, parents, measuredNetGain, assignments.rescuedOrdinals());
		}
		private boolean hasWinner(final int id){return id==first.familyId || second!=null && id==second.familyId;}
		public int[] parentIds(){return parents.clone();}
		public final int requestId;
		public final Family first, second;
		public final LiveProbeAssignments.Snapshot assignments;
		public final long measuredNetGain;
		private final int[] parents;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Read-only borrows      ----------------*/
	/*--------------------------------------------------------------*/

	/** Marks an imported current-state roster with caller-persisted fusion-product eligibility. */
	public LiveProbeState withFusionProductIds(final int[] familyIds){
		require(familyIds!=null, "Fusion-product family IDs are required");
		final boolean[] next=new boolean[models.length];
		for(int id:familyIds){
			final int at=modelIndex.get(id);
			require(at>=0 && models[at].size()>0, "Fusion-product IDs must name active imported models");
			require(!next[at], "Fusion-product IDs must be unique");
			next[at]=true;
		}
		return new LiveProbeState(models, owners, residual, modelIndex, firstNewId,
			assigned, filtered, sourceHeader, lastMarginalGain, next,memberInputs,memberLayer);
	}

	public int[] fusionProductIds(){
		int count=0;
		for(int i=0; i<models.length; i++){if(fusionProducts[i] && models[i].size()>0){count++;}}
		final int[] ids=new int[count];
		for(int i=0, out=0; i<models.length; i++){
			if(fusionProducts[i] && models[i].size()>0){ids[out++]=models[i].familyId;}
		}
		return ids;
	}

	public Family[] models(){return models;}
	/** The coordinator borrows the original cache; null keeps legacy queries explicit. */
	MemberRankInputs memberRankingInputs(){return memberInputs;}
	/** Measurement layer for this immutable publication envelope. */
	int memberRankingLayer(){return memberLayer;}
	public Family[] active(){return active;}
	public int[] residualOrdinals(){return residual;}
	int rowOf(final int familyId){return modelIndex.get(familyId);}
	public int firstNewId(){return firstNewId;}
	public int assigned(){return assigned;}
	public int filtered(){return filtered;}
	public int sequenceCount(){return owners.length;}
	public ClusterSketchCodec.Header sourceHeader(){return sourceHeader;}
	public long lastMarginalGain(){return lastMarginalGain;}
	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}
	private static final int RESIDUAL=-1, FILTERED=-2;
	private final Family[] models, active;
	private final int[] owners, residual;
	private final boolean[] fusionProducts;
	private final IntHashMap modelIndex;
	private final int firstNewId, assigned, filtered;
	private final ClusterSketchCodec.Header sourceHeader;
	private final long lastMarginalGain;
	private final MemberRankInputs memberInputs;
	private final int memberLayer;
}
