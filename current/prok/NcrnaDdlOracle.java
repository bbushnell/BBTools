package prok;

import java.util.Arrays;

import prok.NcrnaGffFamilyGrader.Locus;

/** One worker's exhaustive RF00174 verification result for one captured window.
 * The driver supplies a separate, normally configured observer-free scavenger
 * and the same family's original model-name order. Every model is verified once;
 * primitive accepted spans survive all subsequent sketch geometries for that
 * window. No shortlist, initial identity or sketch score substitutes for the full
 * ordinary verifier. Input bases and seed arrays are borrowed read-only.
 *
 * Results use one-based inclusive ORIGINAL-FORWARD coordinates, matching the
 * reference GFF. They are replaced by verify(), so consume or copy them before
 * the next window. Compatibility links are not one-to-one recoveries: a span
 * compatible with multiple tandem parents must retain that ambiguity in the
 * driver's assignment/reporting, rather than multiplying detections.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlOracle{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlOracle(final String input_, final NcrnaScavenger caller_, final String[] names){
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Exhaustive oracle requires -ea for requested-model identity invariants");}
		if(input_==null || input_.isEmpty() || caller_==null || names==null || names.length<1){
			throw new IllegalArgumentException("Oracle requires an exact input ID, separate configured worker and original model order");
		}
		input=input_;
		caller=caller_;
		modelNames=names.clone();
		caller.validateOracleModels(modelNames);
		starts=new int[names.length];
		stops=new int[names.length];
		partials=new int[names.length];
		scores=new float[names.length];
		parentCounts=new int[names.length];
		uniqueParents=new int[names.length];
	}

	/*--------------------------------------------------------------*/
	/*----------------          Verification        ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces the result with one full verification per original library model.
	 * Bounds and bases are caller-strand oriented, zero-based inclusive. The seed
	 * stream comes from that original oriented contig, not a truth-bounded fragment.
	 * A failed replacement invalidates every accessor, including preceding results.
	 * @return Number of models accepted by the full verifier, not recovered loci */
	public int verify(final String name, final byte[] bases, final int strand_, final int pass,
			final int start, final int stop, final int source, final int[] hits, final long[] keys){
		ready=false;
		attributed=false;
		if(name==null || name.isEmpty() || bases==null || (strand_!=0 && strand_!=1)
			|| (pass!=1 && pass!=2) || start<0 || stop<start || stop>=bases.length || hits==null){
			throw new IllegalArgumentException("Oracle needs a captured native window and its original oriented contig/seed stream");
		}
		if(name!=header){
			contig=seqid(name);
			header=name;
		}
		strand=strand_==0 ? '+' : '-';
		Arrays.fill(starts, -1);
		Arrays.fill(stops, -1);
		Arrays.fill(partials, 0);
		Arrays.fill(scores, Float.NaN);
		int accepted=0;
		for(int model=0; model<modelNames.length; model++){
			singleModel[0]=model;
			final Orf result=caller.verifyOracleModel(name, bases, strand_, pass, start, stop, source, hits, keys, singleModel);
			if(result==null){continue;}
			assert(modelNames[model].equals(result.trnaModel) && name.equals(result.scafName))
				: "Oracle results must preserve the driver's original model order and exact contig: model="+model;
			starts[model]=strand_==0 ? result.start+1 : bases.length-result.stop;
			stops[model]=strand_==0 ? result.stop+1 : bases.length-result.start;
			partials[model]=result.ncrnaPartial;
			scores[model]=result.orfScore;
			assert(starts[model]>0 && stops[model]>=starts[model] && stops[model]<=bases.length)
				: "Forward one-based accepted spans must remain on their original contig";
			accepted++;
		}
		ready=true;
		return accepted;
	}

	/** Writes every model's compatibility with one already-qualified reference parent.
	 * The driver owns subtype/mask qualification; this method cannot admit truth.
	 * Same-input identity is mandatory even when two genomes reuse a contig name.
	 * Output is overwritten and borrowed, with one flag per original model.
	 * @return Number of compatible model spans, not one-to-one detections */
	public int compatible(final String parentInput, final Locus parent, final boolean[] output){
		checkReady();
		if(!input.equals(parentInput) || parent==null || output==null || output.length!=modelNames.length){
			throw new IllegalArgumentException("Oracle compatibility requires this exact input and one output flag per original model");
		}
		Arrays.fill(output, false);
		if(!"cobalamin".equals(parent.family) || !contig.equals(parent.seqid) || strand!=parent.strand){return 0;}
		int count=0;
		for(int model=0; model<modelNames.length; model++){
			if(starts[model]<0){continue;}
			output[model]=overlaps(starts[model], stops[model], parent);
			if(output[model]){count++;}
		}
		assert(count<=modelNames.length) : "Each original model contributes at most one compatibility flag per parent";
		return count;
	}

	/** Attributes spans to qualified SCORED parents, without selecting truth.
	 * Masks/context do not enter this array; native scored eligibility precedes
	 * masks. Multiple compatible parents stay ambiguous, not absent or recovered.
	 * @return Number of accepted model spans compatible with multiple parents */
	public int attribute(final String parentInput, final Locus[] parents){
		checkReady();
		attributed=false;
		if(!input.equals(parentInput)){throw new IllegalArgumentException("Attribution cannot mix exact input namespaces");}
		final int ambiguous=attribute(contig, strand, starts, stops, parents, parentCounts, uniqueParents);
		attributed=true;
		return ambiguous;
	}

	/** Literal primitive-span seam for independent tandem fixtures. Parent indices
	 * refer to the supplied scored list. This computes links, not detection grading.
	 * Inputs are borrowed read-only; the two distinct output arrays are overwritten.
	 */
	public static int attribute(final String contig, final char strand, final int[] starts, final int[] stops,
			final Locus[] parents, final int[] counts, final int[] unique){
		if(contig==null || contig.isEmpty() || (strand!='+' && strand!='-') || starts==null || stops==null
			|| parents==null || counts==null || unique==null || starts.length!=stops.length
			|| starts.length!=counts.length || starts.length!=unique.length || counts==unique
			|| counts==starts || counts==stops || unique==starts || unique==stops){
			throw new IllegalArgumentException("Attribution requires bound primitive spans, scored parents and independent count/index outputs");
		}
		for(final Locus parent:parents){
			if(parent==null){throw new IllegalArgumentException("A scored parent list cannot contain null records");}
		}
		Arrays.fill(counts, 0);
		Arrays.fill(unique, -1);
		int ambiguous=0;
		for(int model=0; model<starts.length; model++){
			if(starts[model]==-1 && stops[model]==-1){continue;}
			if(starts[model]<1 || stops[model]<starts[model]){throw new IllegalArgumentException("Attribution needs complete one-based accepted spans or paired -1 sentinels");}
			for(int parent=0; parent<parents.length; parent++){
				final Locus locus=parents[parent];
				if(!"cobalamin".equals(locus.family) || !contig.equals(locus.seqid) || strand!=locus.strand){continue;}
				if(overlaps(starts[model], stops[model], locus)){
					counts[model]++;
					unique[model]=counts[model]==1 ? parent : -1;
				}
			}
			if(counts[model]>1){ambiguous++;}
		}
		return ambiguous;
	}

	/** Number of scored-parent links, preserving zero versus ambiguous support. */
	public int parentCount(final int model){checkAttribution(model);return parentCounts[model];}
	/** Original supplied parent index only when attribution is unique; otherwise -1. */
	public int uniqueParent(final int model){checkAttribution(model);return uniqueParents[model];}

	/** NcrnaGffFamilyGrader.overlapsReciprocal, without a temporary Locus per pair.
	 * Exact input/contig/strand/family checks precede this primitive predicate. */
	private static boolean overlaps(final int start, final int stop, final Locus parent){
		assert(start>0 && stop>=start) : "Reciprocal overlap consumes already-validated one-based accepted spans";
		final int overlap=Math.min(stop, parent.stop)-Math.max(start, parent.start)+1;
		return overlap>0 && overlap>=.5f*(stop-start+1) && overlap>=.5f*(parent.stop-parent.start+1);
	}

	public int models(){return modelNames.length;}
	public int start(final int model){checkModel(model);return starts[model];}
	public int stop(final int model){checkModel(model);return stops[model];}
	public int partial(final int model){checkModel(model);return partials[model];}
	public float score(final int model){checkModel(model);return scores[model];}

	private void checkModel(final int model){
		checkReady();
		if(model<0 || model>=modelNames.length){throw new IllegalArgumentException("Unknown original oracle model: "+model);}
	}
	private void checkReady(){
		if(!ready){throw new IllegalStateException("Verify one complete window before reading its oracle result");}
	}
	private void checkAttribution(final int model){
		checkModel(model);
		if(!attributed){throw new IllegalStateException("Attribute this verified window before reading scored-parent links");}
	}

	/** Native GFF uses the first space/tab token; retain the full title separately. */
	public static String seqid(final String name){
		if(name==null || name.isEmpty()){throw new IllegalArgumentException("A reference join needs an observed contig title");}
		int end=0;
		while(end<name.length() && name.charAt(end)!=' ' && name.charAt(end)!='\t'){end++;}
		if(end<1){throw new IllegalArgumentException("The observed contig title has an empty native GFF token");}
		return end==name.length() ? name : name.substring(0, end);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String input;
	private final NcrnaScavenger caller;
	private final String[] modelNames;
	private final int[] starts, stops, partials, parentCounts, uniqueParents;
	private final float[] scores;
	private final int[] singleModel=new int[1];
	private String header, contig;
	private char strand;
	private boolean ready, attributed;
}
