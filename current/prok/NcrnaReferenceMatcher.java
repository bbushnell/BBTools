package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import map.IntObjectMap;
import map.ObjectMap;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;
import structures.IntList;

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
		for(final Partition p:partitions){matchPartition(p, calls, references, overlap, result);}
		result.scoredReferences=references.size()-result.maskedReferences;
		result.falseNegatives=result.scoredReferences-result.matched;
		result.falsePositives=calls.size()-result.matched-result.maskedCalls;
		assert(result.falseNegatives>=0 && result.falsePositives>=0)
			: "One-to-one matches and masks must conserve original query/reference rows";
		return result;
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
			Arrays.fill(callToReference, -1);
			Arrays.fill(referenceToCall, -1);
			Arrays.fill(callMaskReference, -1);
		}
		public final int[] callToReference, referenceToCall, callMaskReference;
		public final boolean[] scoredEligible, duplicateCalls;
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
