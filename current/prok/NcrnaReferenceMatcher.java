package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import map.IntObjectMap;
import map.ObjectMap;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;
import structures.IntList;
import structures.LongList;

/** Matches calls against scored references while retaining explicit don't-care masks.
 * Input IDs must identify exact assemblies, not taxa or bare contig names. Family
 * names define the scoring view; a subtype union must be explicitly projected by
 * the caller. This class does not infer biological truth or change coordinates.
 *
 * Scored eligibility takes precedence over masks, even for a duplicate or a call
 * left unmatched by one-to-one assignment. Otherwise any positive overlap with a
 * same-input, same-family, same-contig, same-strand mask makes a call unscored.
 * Exact duplicate queries never consume additional scored references. All rows
 * retain their original list indices, including duplicate and masked queries.
 *
 * Matching uses the existing maximum-cardinality matcher within each input/family
 * partition. Duplicate and adjacency construction are quadratic within a partition;
 * augmenting-path assignment can cost more on dense graphs. Partitions are never
 * compared across the genome union. The caller must supply unique biological
 * reference loci; this matcher does not adjudicate duplicate reference annotations.
 * Inputs are borrowed read-only for the duration of match(); results are owned by
 * the caller. Independent invocations share no mutable state. Linear loci only.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReferenceMatcher{

	private NcrnaReferenceMatcher(){}

	/*--------------------------------------------------------------*/
	/*----------------           Matching           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns original-index pairings and unscored-call witnesses.
	 * @param calls Query loci, each with a positive exact-input ID in callInputs
	 * @param callInputs One assembly ID per query, using the reference ID namespace
	 * @param references Scored and masked reference loci in their original order
	 * @param referenceInputs One positive exact-input ID per reference
	 * @param masks True for unscored reference rows; never inferred from length
	 * @param overlap Minimum reciprocal overlap fraction for scored matches, in (0,1]
	 */
	public static Result match(final List<Locus> calls, final int[] callInputs,
		final List<Locus> references, final int[] referenceInputs, final boolean[] masks,
		final float overlap){
		return match(calls, callInputs, references, referenceInputs, masks, overlap, false);
	}

	/** Optionally adds positive-overlap pairs after freezing the strict assignment.
	 * Added pairs consume only unused unique queries and unused scored references;
	 * they are endpoint-inexact, never additional copies of an existing match.
	 * Strict pairs are preserved even when another strict assignment could allow
	 * more residual pairs. This is a strict-first extension, not joint optimization.
	 */
	public static Result match(final List<Locus> calls, final int[] callInputs,
		final List<Locus> references, final int[] referenceInputs, final boolean[] masks,
		final float overlap, final boolean includeInexact){
		require(calls!=null && references!=null && callInputs!=null && referenceInputs!=null && masks!=null,
			"Calls, references and their input/mask arrays are required");
		require(calls.size()==callInputs.length && references.size()==referenceInputs.length
			&& references.size()==masks.length, "Every original row requires exactly one input ID and reference mask flag");
		require(overlap>0 && overlap<=1, "Reciprocal scored overlap must be finite and in (0,1]: "+overlap);
		final IntObjectMap<InputGroup> inputs=new IntObjectMap<InputGroup>(16, InputGroup.class);
		final ArrayList<Partition> partitions=new ArrayList<Partition>();
		final Result result=new Result(calls.size(), references.size());
		for(int i=0; i<references.size(); i++){
			final Locus locus=references.get(i);
			final Partition p=partition(inputs, partitions, referenceInputs[i], locus);
			if(masks[i]){
				p.masks.add(i);
				result.maskedReferences++;
			}else{p.references.add(i);}
		}
		for(int i=0; i<calls.size(); i++){
			partition(inputs, partitions, callInputs[i], calls.get(i)).calls.add(i);
		}
		for(final Partition p:partitions){
			matchPartition(p, calls, references, overlap, result);
			if(includeInexact){matchInexact(p, calls, references, overlap, result);}
		}
		result.scoredReferences=references.size()-result.maskedReferences;
		result.falseNegatives=result.scoredReferences-result.matched;
		result.falsePositives=calls.size()-result.matched-result.maskedCalls;
		assert(result.falseNegatives>=0 && result.falsePositives>=0)
			: "One-to-one matches and masks must conserve original query/reference rows";
		return result;
	}

	/** Preserve strict owners, then maximize residual positive-overlap pairs.
	 * Broad scored eligibility precedes masks even when no unused reference is
	 * available; otherwise a duplicate could disappear into a don't-care region.
	 */
	private static void matchInexact(final Partition p, final List<Locus> calls,
		final List<Locus> references, final float overlap, final Result result){
		final ArrayList<Locus> residualCalls=new ArrayList<Locus>();
		final ArrayList<Locus> residualReferences=new ArrayList<Locus>();
		final IntList callIds=new IntList();
		final IntList referenceIds=new IntList();
		for(int i=0; i<p.references.size; i++){
			final int ri=p.references.get(i);
			if(result.referenceToCall[ri]<0){
				residualReferences.add(references.get(ri));
				referenceIds.add(ri);
			}
		}
		for(int i=0; i<p.calls.size; i++){
			final int ci=p.calls.get(i);
			final Locus call=calls.get(ci);
			boolean eligible=result.scoredEligible[ci];
			for(int j=0; !eligible && j<p.references.size; j++){
				eligible=NcrnaGffFamilyGrader.overlapsReciprocal(call, references.get(p.references.get(j)), 0);
			}
			result.scoredEligible[ci]=eligible;
			if(eligible && result.callMaskReference[ci]>=0){
				result.callMaskReference[ci]=-1;
				result.maskedCalls--;
			}
			if(result.callToReference[ci]<0 && !result.duplicateCalls[ci] && result.callMaskReference[ci]<0){
				residualCalls.add(call);
				callIds.add(ci);
			}
		}
		final MatchResult added=NcrnaGffFamilyGrader.matchOneToOne(residualCalls, residualReferences, 0);
		for(int i=0; i<callIds.size; i++){
			final int local=added.callMatchedTruthIdx[i];
			if(local<0){continue;}
			final int ci=callIds.get(i);
			final int ri=referenceIds.get(local);
			require(result.callToReference[ci]<0 && result.referenceToCall[ri]<0
				&& !NcrnaGffFamilyGrader.overlapsReciprocal(calls.get(ci), references.get(ri), overlap),
				"Residual pairs must be unused endpoint-inexact rows after maximum-cardinality strict matching");
			result.callToReference[ci]=ri;
			result.referenceToCall[ri]=ci;
			result.endpointInexact[ci]=true;
		}
		result.matched+=added.matched;
		result.inexactMatched+=added.matched;
	}

	/** Groups by exact input and family without boxed IDs or concatenated row keys. */
	private static Partition partition(final IntObjectMap<InputGroup> inputs,
		final ArrayList<Partition> order, final int input, final Locus locus){
		require(input>0, "Exact input IDs must be positive: "+input);
		require(locus!=null && locus.seqid!=null && !locus.seqid.isEmpty()
			&& locus.family!=null && !locus.family.isEmpty(), "A locus needs literal contig and family identity");
		require(locus.strand=='+' || locus.strand=='-', "Unknown strand cannot establish scored or masked eligibility");
		InputGroup g=inputs.get(input);
		if(g==null){
			g=new InputGroup();
			inputs.put(input, g);
		}
		Partition p=g.families.get(locus.family);
		if(p==null){
			p=new Partition();
			g.families.put(locus.family, p);
			order.add(p);
		}
		return p;
	}

	/** Scored-edge existence is checked before duplicate suppression or masking.
	 * Duplicate detection uses retained loci within this input/family partition;
	 * there is no temporary composite String or boxed coordinate per comparison.
	 */
	private static void matchPartition(final Partition p, final List<Locus> calls,
		final List<Locus> references, final float overlap, final Result result){
		final ArrayList<Locus> scored=new ArrayList<Locus>(p.references.size);
		for(int i=0; i<p.references.size; i++){scored.add(references.get(p.references.get(i)));}
		final ArrayList<Locus> unique=new ArrayList<Locus>(p.calls.size);
		final IntList originalCalls=new IntList(1);
		for(int i=0; i<p.calls.size; i++){
			final int ci=p.calls.get(i);
			final Locus call=calls.get(ci);
			boolean duplicate=false;
			for(int j=0; j<i; j++){
				final Locus previous=calls.get(p.calls.get(j));
				if(call.start==previous.start && call.stop==previous.stop
					&& sameTarget(call, previous)){
					duplicate=true;
					break;
				}
			}
			result.duplicateCalls[ci]=duplicate;
			if(duplicate){result.duplicates++;}
			boolean eligible=false;
			for(final Locus ref:scored){
				if(NcrnaGffFamilyGrader.overlapsReciprocal(call, ref, overlap)){
					eligible=true;
					break;
				}
			}
			result.scoredEligible[ci]=eligible;
			if(!eligible){
				for(int j=0; j<p.masks.size; j++){
					final int ri=p.masks.get(j);
					final Locus mask=references.get(ri);
					if(sameTarget(call, mask) && call.start<=mask.stop && mask.start<=call.stop){
						result.callMaskReference[ci]=ri;
						result.maskedCalls++;
						break;
					}
				}
			}
			if(!duplicate && result.callMaskReference[ci]<0){
				unique.add(call);
				originalCalls.add(ci);
			}
		}
		final MatchResult match=NcrnaGffFamilyGrader.matchOneToOne(unique, scored, overlap);
		result.matched+=match.matched;
		for(int i=0; i<originalCalls.size; i++){
			final int localReference=match.callMatchedTruthIdx[i];
			if(localReference<0){continue;}
			final int ci=originalCalls.get(i), ri=p.references.get(localReference);
			assert(result.referenceToCall[ri]<0 && result.callToReference[ci]<0
				&& result.callMaskReference[ci]<0) : "A scored pair cannot consume a masked or already assigned original row";
			result.callToReference[ci]=ri;
			result.referenceToCall[ri]=ci;
		}
	}

	/** Family and input equality have already been enforced by the partition. */
	private static boolean sameTarget(final Locus a, final Locus b){
		return a.strand==b.strand && a.seqid.equals(b.seqid);
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}

	/** Fixed input/family partitions, duplicate groups and scored/masked eligibility.
	 * Constructor lists and arrays are borrowed only while preparing the pool.
	 * Repeated matches take distinct pool indices in native selected-call order;
	 * exact-coordinate duplicates are determined within THAT selection, not against
	 * inactive candidates. Scored eligibility always precedes mask handling.
	 * One instance belongs to one thread. Pairing accessors describe only the last
	 * completed trial; copyResult() is the explicit allocating, owned-result path.
	 */
	public static final class Prepared{
		public Prepared(final List<Locus> candidates, final int[] candidateInputs,
				final List<Locus> references, final int[] referenceInputs, final boolean[] masks,
				final float overlap){
			this(candidates, candidateInputs, references, referenceInputs, masks, overlap, false);
		}

		/** Optional strict-first inclusive matching, equivalent to the allocating
		 * match(..., includeInexact) route. Strict owners remain fixed; only unused
		 * unique calls/references enter residual positive-overlap augmentation.
		 */
		public Prepared(final List<Locus> candidates, final int[] candidateInputs,
				final List<Locus> references, final int[] referenceInputs, final boolean[] masks,
				final float overlap, final boolean includeInexact_){
			includeInexact=includeInexact_;
			require(candidates!=null && references!=null && candidateInputs!=null && referenceInputs!=null && masks!=null,
				"Prepared grading requires fixed candidates, references and input/mask arrays");
			require(candidates.size()==candidateInputs.length && references.size()==referenceInputs.length
				&& references.size()==masks.length, "Every prepared row requires exactly one input ID and reference mask flag");
			require(overlap>0 && overlap<=1, "Prepared reciprocal overlap must be finite and in (0,1]");
			final int count=candidates.size();
			result=new Result(count, references.size());
			positions=new int[count];
			localIndices=new int[count];
			duplicateGroups=new int[count];
			groupSeen=new boolean[count];
			eligible=new boolean[count];
			maskReferences=new int[count];
			Arrays.fill(maskReferences, -1);
			owners=new PreparedPartition[count];
			final IntObjectMap<InputGroup> inputs=new IntObjectMap<InputGroup>(16, InputGroup.class);
			final ArrayList<Partition> groups=new ArrayList<Partition>();
			for(int i=0; i<references.size(); i++){
				final Partition p=partition(inputs, groups, referenceInputs[i], references.get(i));
				if(masks[i]){
					p.masks.add(i);
					result.maskedReferences++;
				}else{p.references.add(i);}
			}
			for(int i=0; i<count; i++){
				partition(inputs, groups, candidateInputs[i], candidates.get(i)).calls.add(i);
			}
			partitions=new PreparedPartition[groups.size()];
			for(int pi=0; pi<groups.size(); pi++){
				final Partition p=groups.get(pi);
				final ArrayList<Locus> calls=new ArrayList<Locus>(p.calls.size);
				final ArrayList<Locus> scored=new ArrayList<Locus>(p.references.size);
				for(int i=0; i<p.references.size; i++){scored.add(references.get(p.references.get(i)));}
				for(int i=0; i<p.calls.size; i++){
					final int candidate=p.calls.get(i);
					final Locus call=candidates.get(candidate);
					calls.add(call);
					localIndices[candidate]=i;
					duplicateGroups[candidate]=candidate;
					for(int j=0; j<i; j++){
						final int previous=p.calls.get(j);
						final Locus other=candidates.get(previous);
						if(call.start==other.start && call.stop==other.stop && sameTarget(call, other)){
							duplicateGroups[candidate]=duplicateGroups[previous];
							break;
						}
					}
					for(final Locus ref:scored){
						if(NcrnaGffFamilyGrader.overlapsReciprocal(call, ref, includeInexact ? 0 : overlap)){
							eligible[candidate]=true;
							break;
						}
					}
					if(!eligible[candidate]){
						for(int j=0; j<p.masks.size; j++){
							final int reference=p.masks.get(j);
							final Locus mask=references.get(reference);
							if(sameTarget(call, mask) && call.start<=mask.stop && mask.start<=call.stop){
								maskReferences[candidate]=reference;
								break;
							}
						}
					}
				}
				final PreparedPartition prepared=new PreparedPartition(p, calls, scored, overlap, includeInexact);
				partitions[pi]=prepared;
				for(int i=0; i<p.calls.size; i++){owners[p.calls.get(i)]=prepared;}
			}
			result.scoredReferences=references.size()-result.maskedReferences;
		}

		/** Matches every prepared candidate in original pool order. */
		public int matchAll(){
			begin();
			for(int i=0; i<positions.length; i++){activate(i);}
			return finish();
		}
		/** Matches a selected subset without allocating loci, partitions or results. */
		public int match(final IntList order){
			require(order!=null, "Prepared grading requires an explicit selected-candidate order");
			begin();
			for(int i=0; i<order.size; i++){activate(order.get(i));}
			return finish();
		}
		private void begin(){
			ready=false;
			calls=0;
			Arrays.fill(positions, -1);
			Arrays.fill(groupSeen, false);
			Arrays.fill(result.callToReference, -1);
			Arrays.fill(result.referenceToCall, -1);
			result.matched=0;
			result.inexactMatched=0;
			Arrays.fill(result.endpointInexact, false);
			result.maskedCalls=0;
			result.duplicates=0;
			for(final PreparedPartition p:partitions){p.order.clear();}
		}
		private void activate(final int candidate){
			if(candidate<0 || candidate>=positions.length || positions[candidate]>=0){
				throw new IllegalArgumentException("Selected candidate must be unique and inside the prepared grading pool: "+candidate);
			}
			final int ci=calls++;
			positions[candidate]=ci;
			final int group=duplicateGroups[candidate];
			final boolean duplicate=groupSeen[group];
			groupSeen[group]=true;
			result.duplicateCalls[ci]=duplicate;
			result.scoredEligible[ci]=eligible[candidate];
			result.callMaskReference[ci]=maskReferences[candidate];
			if(duplicate){result.duplicates++;}
			if(maskReferences[candidate]>=0){result.maskedCalls++;}
			if(!duplicate && maskReferences[candidate]<0){owners[candidate].order.add(localIndices[candidate]);}
		}
		private int finish(){
			for(final PreparedPartition p:partitions){
				result.matched+=p.matcher.match(p.order);
				for(int i=0; i<p.order.size; i++){
					final int reference=p.matcher.matchedReference(i);
					if(reference<0){continue;}
					final int ci=positions[p.candidates[p.order.get(i)]], ri=p.references[reference];
					assert(ci>=0 && result.referenceToCall[ri]<0)
						: "Prepared partitions must retain unique original-reference owners and selected query ordinals";
					result.callToReference[ci]=ri;
					result.referenceToCall[ri]=ci;
				}
				if(includeInexact){p.matchInexact(result, positions);}
			}
			result.falseNegatives=result.scoredReferences-result.matched;
			result.falsePositives=calls-result.matched-result.maskedCalls;
			assert(result.falseNegatives>=0 && result.falsePositives>=0)
				: "Prepared one-to-one matches and masks must conserve the selected rows";
			ready=true;
			return result.matched;
		}

		/** Returns an independent exact-sized result, retained across later trials. */
		public Result copyResult(){
			require(ready, "Copying prepared results requires a completed match");
			final Result copy=new Result(calls, result.referenceToCall.length);
			System.arraycopy(result.callToReference, 0, copy.callToReference, 0, calls);
			System.arraycopy(result.callMaskReference, 0, copy.callMaskReference, 0, calls);
			System.arraycopy(result.scoredEligible, 0, copy.scoredEligible, 0, calls);
			System.arraycopy(result.duplicateCalls, 0, copy.duplicateCalls, 0, calls);
			System.arraycopy(result.endpointInexact, 0, copy.endpointInexact, 0, calls);
			System.arraycopy(result.referenceToCall, 0, copy.referenceToCall, 0, copy.referenceToCall.length);
			copy.matched=result.matched;
			copy.inexactMatched=result.inexactMatched;
			copy.falsePositives=result.falsePositives;
			copy.falseNegatives=result.falseNegatives;
			copy.scoredReferences=result.scoredReferences;
			copy.maskedReferences=result.maskedReferences;
			copy.maskedCalls=result.maskedCalls;
			copy.duplicates=result.duplicates;
			return copy;
		}
		/** Original scored reference for a selected query ordinal, or -1. */
		public int matchedReference(final int call){
			checkCall(call);
			return result.callToReference[call];
		}
		/** Original mask witness for a selected query ordinal, or -1. */
		public int maskReference(final int call){
			checkCall(call);
			return result.callMaskReference[call];
		}
		public boolean duplicateCall(final int call){
			checkCall(call);
			return result.duplicateCalls[call];
		}
		public boolean scoredEligible(final int call){
			checkCall(call);
			return result.scoredEligible[call];
		}
		/** Selected query ordinal paired to an original reference, or -1. */
		public int matchedCall(final int reference){
			require(ready && reference>=0 && reference<result.referenceToCall.length,
				"Reference lookup requires a completed prepared match and original reference index");
			return result.referenceToCall[reference];
		}
		public int callCount(){
			checkResult();
			return calls;
		}
		public int matched(){
			checkResult();
			return result.matched;
		}
		/** Added positive-overlap matches, excluded from endpoint accuracy. */
		public int inexactMatched(){
			checkResult();
			return result.inexactMatched;
		}
		public boolean endpointInexact(final int call){
			checkCall(call);
			return result.endpointInexact[call];
		}
		public int falsePositives(){
			checkResult();
			return result.falsePositives;
		}
		public int falseNegatives(){
			checkResult();
			return result.falseNegatives;
		}
		public int maskedCalls(){
			checkResult();
			return result.maskedCalls;
		}
		public int duplicates(){
			checkResult();
			return result.duplicates;
		}
		private void checkResult(){require(ready, "Counters require a completed prepared match");}
		private void checkCall(final int call){
			require(ready && call>=0 && call<calls, "Query lookup requires a completed prepared match and selected ordinal");
		}
		private final PreparedPartition[] partitions, owners;
		private final int[] positions, localIndices, duplicateGroups, maskReferences;
		private final boolean[] groupSeen, eligible;
		private final Result result;
		private final boolean includeInexact;
		private int calls;
		private boolean ready;
	}

	/** One immutable input/family pool and its reusable native-order selection. */
	private static final class PreparedPartition{
		PreparedPartition(final Partition p, final List<Locus> calls, final List<Locus> truth,
				final float overlap, final boolean includeInexact){
			candidates=p.calls.toArray();
			references=p.references.toArray();
			order=new IntList(candidates.length);
			matcher=new NcrnaGffFamilyGrader.PreparedMatcher(calls, truth, overlap, 0);
			if(!includeInexact){
				offsets=edges=residualOwners=visited=null;
				return;
			}
			offsets=new int[Math.addExact(calls.size(), 1)];
			final IntList adjacency=new IntList();
			final LongList row=new LongList();
			for(int ci=0; ci<calls.size(); ci++){
				offsets[ci]=adjacency.size;
				row.clear();
				final Locus call=calls.get(ci);
				for(int ri=0; ri<truth.size(); ri++){
					final Locus reference=truth.get(ri);
					if(!sameTarget(call, reference)){continue;}
					final int shared=Math.min(call.stop, reference.stop)-Math.max(call.start, reference.start)+1;
					if(shared>0){row.add(((long)(Integer.MAX_VALUE-shared)<<32)|ri);}
				}
				// Same overlap-descending/reference-index tie order as PreparedMatcher.
				Arrays.sort(row.array, 0, row.size);
				for(int i=0; i<row.size; i++){adjacency.add((int)row.get(i));}
			}
			offsets[calls.size()]=adjacency.size;
			edges=adjacency.toArray();
			residualOwners=new int[truth.size()];
			visited=new int[truth.size()];
		}
		/** Reuses immutable edges and primitive owners; never reroutes a strict pair. */
		void matchInexact(final Result result, final int[] positions){
			Arrays.fill(residualOwners, -1);
			for(int i=0; i<order.size; i++){
				final int candidate=order.get(i);
				if(result.callToReference[positions[candidates[candidate]]]>=0){continue;}
				if(stamp==Integer.MAX_VALUE){Arrays.fill(visited, 0); stamp=0;}
				stamp++;
				if(augment(candidate, result)){
					result.matched++;
					result.inexactMatched++;
				}
			}
			for(int reference=0; reference<residualOwners.length; reference++){
				final int candidate=residualOwners[reference];
				if(candidate<0){continue;}
				final int ci=positions[candidates[candidate]], ri=references[reference];
				assert(ci>=0 && result.callToReference[ci]<0 && result.referenceToCall[ri]<0)
					: "Residual augmentation may consume only unused selected queries and references; strict owners stay fixed";
				result.callToReference[ci]=ri;
				result.referenceToCall[ri]=ci;
				result.endpointInexact[ci]=true;
			}
		}
		private boolean augment(final int candidate, final Result result){
			assert(candidate>=0 && candidate<candidates.length)
				: "Residual owners use this prepared partition's candidate indices";
			for(int i=offsets[candidate]; i<offsets[candidate+1]; i++){
				final int reference=edges[i];
				if(visited[reference]==stamp || result.referenceToCall[references[reference]]>=0){continue;}
				visited[reference]=stamp;
				if(residualOwners[reference]<0 || augment(residualOwners[reference], result)){
					residualOwners[reference]=candidate;
					return true;
				}
			}
			return false;
		}
		final int[] candidates, references;
		final IntList order;
		final NcrnaGffFamilyGrader.PreparedMatcher matcher;
		final int[] offsets, edges, residualOwners, visited;
		int stamp;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Types             ----------------*/
	/*--------------------------------------------------------------*/

	/** All index arrays refer to original input-list positions. -1 means no pairing.
	 * callMaskReference records the first qualifying mask, a witness rather than a
	 * consuming assignment. Masks never contribute matched pairs or false negatives.
	 * Duplicate flags are independent of scored/masked outcomes.
	 */
	public static final class Result{
		private Result(final int calls, final int references){
			callToReference=new int[calls];
			referenceToCall=new int[references];
			callMaskReference=new int[calls];
			scoredEligible=new boolean[calls];
			duplicateCalls=new boolean[calls];
			endpointInexact=new boolean[calls];
			Arrays.fill(callToReference, -1);
			Arrays.fill(referenceToCall, -1);
			Arrays.fill(callMaskReference, -1);
		}
		public final int[] callToReference, referenceToCall, callMaskReference;
		public final boolean[] scoredEligible, duplicateCalls, endpointInexact;
		/** Additional positive-overlap pairs failing the strict reciprocal criterion. */
		public int inexactMatched;
		public int matched, falsePositives, falseNegatives;
		public int scoredReferences, maskedReferences, maskedCalls, duplicates;
	}

	private static final class InputGroup{
		final ObjectMap<String, Partition> families=new ObjectMap<String, Partition>(4, String.class, Partition.class);
	}
	private static final class Partition{
		final IntList calls=new IntList(1), references=new IntList(1), masks=new IntList(1);
	}
}
