package prok;

import java.util.Arrays;

import structures.IntList;

/** Greedy model-subset proposal from gene-to-model support edges.
 * Removes the model with the fewest uniquely supported genes until another
 * removal would breach the requested fraction of ALL input genes. This is a
 * proposal for subsequent whole-genome calling, not a caller-preservation test
 * or a proof of minimum cardinality. Inputs are borrowed and never mutated;
 * all mutable state belongs to one invocation.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaModelSubset{

	/*--------------------------------------------------------------*/
	/*----------------           Reduction          ----------------*/
	/*--------------------------------------------------------------*/

	/** Reduces a sparse support graph without repeating alignments.
	 * @param models Number of models, whose IDs are in [0, models)
	 * @param starts Gene offsets into supporters; length is gene count plus one
	 * @param supporters Distinct model IDs within each gene's slice
	 * @param target Required fraction of all genes with at least one kept supporter
	 * @return Kept-model flags, removal order and original/final coverage counts
	 * @throws IllegalArgumentException If the graph is malformed or target unattainable
	 */
	public static Result reduce(final int models, final int[] starts,
			final int[] supporters, final double target){
		return reduce(models, starts, supporters, target, Integer.MAX_VALUE);
	}

	/** Adds a per-deletion loss cap to the existing total-coverage target.
	 * Loss means genes that become unsupported now, not all genes a model supports.
	 * Incremental counts are updated after every drop, including zero-loss drops.
	 * @param maxModelLoss Maximum newly uncovered genes permitted for one deletion
	 */
	public static Result reduce(final int models, final int[] starts,
			final int[] supporters, final double target, final int maxModelLoss){
		return reduce(models, starts, supporters, target, maxModelLoss, 0);
	}

	/** Stops the same greedy trajectory before fewer than minModels remain.
	 * Coverage and per-deletion loss limits may stop it earlier with more models.
	 * A zero floor preserves the original reduction behavior.
	 * @param minModels Minimum retained model count, in [0, models]
	 */
	public static Result reduce(final int models, final int[] starts,
			final int[] supporters, final double target, final int maxModelLoss,
			final int minModels){
		if(maxModelLoss<0){throw new IllegalArgumentException("maxmodelloss must be nonnegative");}
		if(minModels<0 || minModels>models){
			throw new IllegalArgumentException("minmodels must be between zero and the input model count: "+minModels+", "+models);
		}
		validate(models, starts, supporters, target);
		final int genes=starts.length-1;
		final int required=minimumCovered(genes, target);
		final int[] remaining=new int[genes];
		final int[] unique=new int[models];
		final int[] modelStarts=new int[models+1];
		int covered=0;
		for(int gene=0; gene<genes; gene++){
			final int count=starts[gene+1]-starts[gene];
			remaining[gene]=count;
			if(count>0){covered++;}
			if(count==1){unique[supporters[starts[gene]]]++;}
			for(int edge=starts[gene]; edge<starts[gene+1]; edge++){
				modelStarts[supporters[edge]+1]++;
			}
		}
		if(covered<required){
			throw new IllegalArgumentException("Coverage target is unattainable: "+
				covered+" of "+genes+" genes have support; "+required+" are required.");
		}

		// Reverse adjacency visits only genes affected by a removed model.
		for(int model=0; model<models; model++){
			modelStarts[model+1]+=modelStarts[model];
		}
		final int[] cursor=Arrays.copyOf(modelStarts, models);
		final int[] modelGenes=new int[supporters.length];
		for(int gene=0; gene<genes; gene++){
			for(int edge=starts[gene]; edge<starts[gene+1]; edge++){
				modelGenes[cursor[supporters[edge]]++]=gene;
			}
		}
		final boolean[] kept=new boolean[models];
		Arrays.fill(kept, true);
		final IntList removed=new IntList(models);
		final int originalCovered=covered;
		int retained=models;
		while(retained>minModels){
			final int victim=leastUnique(kept, unique);
			if(victim<0 || unique[victim]>maxModelLoss || covered-unique[victim]<required){break;}
			final int expectedLoss=unique[victim];
			int lost=0;
			kept[victim]=false;
			removed.add(victim);
			retained--;
			for(int edge=modelStarts[victim]; edge<modelStarts[victim+1]; edge++){
				final int gene=modelGenes[edge];
				assert(remaining[gene]>0) : "A kept edge must contribute to its gene's supporter count; gene="+gene;
				remaining[gene]--;
				if(remaining[gene]==0){
					lost++;
				}else if(remaining[gene]==1){
					// A gene reaches this transition at most once during the reduction.
					final int survivor=soleSupporter(gene, starts, supporters, kept);
					unique[survivor]++;
				}
			}
			assert(lost==expectedLoss) : "Unique-support counts predict the exact coverage loss on removal; model="+
				victim+" expected="+expectedLoss+" observed="+lost;
			covered-=lost;
		}
		assert(covered>=required) : "A proposal must preserve the requested corpus coverage before caller testing";
		assert(retained>=minModels && retained+removed.size==models)
			: "Each greedy deletion removes exactly one model and must preserve the requested model-count floor";
		return new Result(kept, removed.toArray(), genes, originalCovered, covered);
	}

	/** Uses the reported coverage comparison, avoiding ceil(product) rounding up
	 * an otherwise attainable decimal fraction, such as seven out of twenty-five. */
	private static int minimumCovered(final int genes, final double target){
		assert(genes>0 && target>0 && target<=1) : "Coverage requires a nonempty population and a fractional target";
		int required=(int)(target*genes);
		while(required/(double)genes<target){required++;}
		while(required>0 && (required-1)/(double)genes>=target){required--;}
		return required;
	}

	/** Finds the least uniquely supporting kept model; library order resolves ties. */
	private static int leastUnique(final boolean[] kept, final int[] unique){
		assert(kept.length==unique.length) : "Model flags and unique-support counters must share their index space";
		int best=-1;
		for(int model=0; model<kept.length; model++){
			if(kept[model] && (best<0 || unique[model]<unique[best])){best=model;}
		}
		return best;
	}

	/** Finds the remaining supporter after a gene's count changes from two to one. */
	private static int soleSupporter(final int gene, final int[] starts,
			final int[] supporters, final boolean[] kept){
		for(int edge=starts[gene]; edge<starts[gene+1]; edge++){
			final int model=supporters[edge];
			if(kept[model]){return model;}
		}
		throw new IllegalStateException("Support count says one model remains, but no kept edge exists for gene "+gene);
	}

	/** Validates framing and distinct edges before primitive-array indexing. */
	private static void validate(final int models, final int[] starts,
			final int[] supporters, final double target){
		if(models<1 || starts==null || starts.length<2 || supporters==null){
			throw new IllegalArgumentException("Reduction requires models and at least one gene in a framed support graph");
		}
		if(!(target>0 && target<=1)){
			throw new IllegalArgumentException("target must be in (0,1]");
		}
		if(starts[0]!=0 || starts[starts.length-1]!=supporters.length){
			throw new IllegalArgumentException("Gene offsets must cover the complete support-edge array");
		}
		final int[] seen=new int[models];
		Arrays.fill(seen, -1);
		for(int gene=0; gene<starts.length-1; gene++){
			final int from=starts[gene], to=starts[gene+1];
			if(from<0 || to<from || to>supporters.length){
				throw new IllegalArgumentException("Invalid support offsets for gene "+gene);
			}
			for(int edge=from; edge<to; edge++){
				final int model=supporters[edge];
				if(model<0 || model>=models || seen[model]==gene){
					throw new IllegalArgumentException("Support IDs must be distinct and within the library; gene="+gene);
				}
				seen[model]=gene;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Result            ----------------*/
	/*--------------------------------------------------------------*/

	/** One proposal, owned by its recipient. Arrays never alias input graph arrays. */
	public static final class Result{
		private Result(final boolean[] kept_, final int[] removed_, final int genes_,
				final int originalCovered_, final int covered_){
			kept=kept_;
			removed=removed_;
			genes=genes_;
			originalCovered=originalCovered_;
			covered=covered_;
		}
		/** Flags in the original library's model order. */
		public final boolean[] kept;
		/** Model IDs in their actual greedy removal order. */
		public final int[] removed;
		/** All input genes, including genes unsupported by the original library. */
		public final int genes;
		/** Number supported before reduction. */
		public final int originalCovered;
		/** Number supported by the proposed subset. */
		public final int covered;
	}

	private NcrnaModelSubset(){}
}
