package test.prok;

import java.util.Arrays;

import prok.NcrnaModelSubset;
import prok.NcrnaModelSubset.Result;
import structures.IntList;

/** Literal and exhaustive small-graph checks for the greedy proposal contract.
 * The oracle counts support directly from the original edges; it does not use
 * the production reverse adjacency or incremental counters.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaModelSubsetTest{

	public static void main(final String[] args){
		literalCases();
		final int checked=smallGraphs();
		System.err.println("Model subset checks passed: literal cases and "+checked+" feasible small graphs.");
	}

	/** Exercises a new unique-support transition and the all-genes denominator. */
	private static void literalCases(){
		final int[] starts={0, 2, 4, 6, 7};
		final int[] edges={0, 1, 0, 1, 1, 2, 2};
		final Result result=NcrnaModelSubset.reduce(3, starts, edges, 0.75);
		assert(Arrays.equals(result.kept, new boolean[]{false, true, false})) :
			"After removing model 0, model 1 uniquely supports two genes; model 2 is the next legal deletion";
		assert(result.covered==3 && result.originalCovered==4) : "Only the fourth gene loses its last supporter";
		final Result lossless=NcrnaModelSubset.reduce(3, starts, edges, 0.25, 0);
		assert(Arrays.equals(lossless.kept, new boolean[]{false, true, true}) && lossless.covered==4)
			: "Explicit zero permits redundant-model removal but must stop before the next newly uncovered gene";
		final Result oneLoss=NcrnaModelSubset.reduce(3, starts, edges, 0.25, 1);
		assert(Arrays.equals(oneLoss.kept, new boolean[]{false, true, false}) && oneLoss.covered==3)
			: "Model 2 supports two genes but loses only one; later model 1 uniquely supports three and exceeds the cap";
		final Result totalFloor=NcrnaModelSubset.reduce(3, starts, edges, 1, 1);
		assert(totalFloor.covered==4 && Arrays.equals(totalFloor.kept, lossless.kept))
			: "A permissive per-drop cap must not override the stricter total coverage target";
		boolean badCap=false;
		try{
			NcrnaModelSubset.reduce(3, starts, edges, 0.25, -1);
		}catch(IllegalArgumentException expected){badCap=expected.getMessage().contains("maxmodelloss");}
		assert(badCap) : "Only omission selects the legacy policy; the explicit API cap cannot be negative";
		final int[] redundantStarts={0, 3, 6};
		final int[] redundantEdges={0, 1, 2, 0, 1, 2};
		final Result bounded=NcrnaModelSubset.reduce(3, redundantStarts, redundantEdges, 1, 0, 2);
		assert(Arrays.equals(bounded.kept, new boolean[]{false, true, true})
			&& Arrays.equals(bounded.removed, new int[]{0}) && bounded.covered==2)
			: "The count floor stops zero-loss deletions without changing the original library-order tie break";
		final Result untouched=NcrnaModelSubset.reduce(3, redundantStarts, redundantEdges, 1, 0, 3);
		assert(untouched.removed.length==0 && untouched.covered==2)
			: "A floor equal to the original library count permits no deletion";
		for(final int floor : new int[]{-1, 4}){
			boolean badFloor=false;
			try{
				NcrnaModelSubset.reduce(3, redundantStarts, redundantEdges, 1, 0, floor);
			}catch(IllegalArgumentException expected){
				badFloor=expected.getMessage().contains("minmodels");
			}
			assert(badFloor) : "A model-count floor outside the supplied library axis must fail before reduction";
		}
		final int[] partialStarts={0, 1, 2, 2, 2};
		final int[] partialEdges={0, 1};
		final Result partial=NcrnaModelSubset.reduce(2, partialStarts, partialEdges, 0.5);
		assert(partial.kept[0] && partial.kept[1] && partial.genes==4) :
			"Originally unsupported genes remain in the target denominator";
		boolean rejected=false;
		try{
			NcrnaModelSubset.reduce(2, partialStarts, partialEdges, 0.75);
		}catch(IllegalArgumentException expected){
			rejected=expected.getMessage().contains("2 of 4");
		}
		assert(rejected) : "Unattainable coverage must report the available count rather than return an invalid proposal";
		final int[] decimalStarts=new int[26];
		final int[] decimalEdges=new int[25];
		for(int gene=0; gene<25; gene++){
			decimalStarts[gene+1]=gene+1;
			decimalEdges[gene]=gene;
		}
		for(final double target : new double[]{Math.nextDown(0.28), 0.28, Math.nextUp(0.28)}){
			final Result decimal=NcrnaModelSubset.reduce(25, decimalStarts, decimalEdges, target);
			final int expected=target>0.28 ? 8 : 7;
			assert(decimal.covered==expected) : "The coverage ratio itself governs the threshold, without an arbitrary epsilon";
		}
	}

	/** Enumerates every support relation between four genes and three models. */
	private static int smallGraphs(){
		int checked=0;
		for(int graph=0; graph<4096; graph++){
			final IntList edges=new IntList();
			final int[] starts=new int[5];
			for(int gene=0; gene<4; gene++){
				for(int model=0; model<3; model++){
					if((graph&(1<<(3*gene+model)))!=0){edges.add(model);}
				}
				starts[gene+1]=edges.size;
			}
			final int[] support=edges.toArray();
			final int[] savedStarts=starts.clone();
			final int[] savedSupport=support.clone();
			final int attainable=count(starts, support, new boolean[]{true, true, true});
			for(int required=1; required<=attainable; required++){
				final Result result=NcrnaModelSubset.reduce(3, starts, support, required/4.0);
				checkTrace(result, starts, support, required, Integer.MAX_VALUE);
				final Result uncapped=NcrnaModelSubset.reduce(3, starts, support, required/4.0, Integer.MAX_VALUE);
				assert(Arrays.equals(result.kept, uncapped.kept) && Arrays.equals(result.removed, uncapped.removed)
					&& result.covered==uncapped.covered) : "The legacy API must preserve its exact deletion order and result";
				for(int cap=0; cap<=2; cap++){
					final Result capped=NcrnaModelSubset.reduce(3, starts, support, required/4.0, cap);
					checkTrace(capped, starts, support, required, cap);
					for(int floor=0; floor<=3; floor++){
						final Result bounded=NcrnaModelSubset.reduce(3, starts, support, required/4.0, cap, floor);
						checkTrace(bounded, starts, support, required, cap, floor);
						assert(bounded.removed.length<=capped.removed.length)
							: "A count floor may shorten, but never extend, the same capped greedy trajectory";
						for(int i=0; i<bounded.removed.length; i++){
							assert(bounded.removed[i]==capped.removed[i])
								: "Stopping early must preserve every preceding greedy deletion and its tie break";
						}
					}
				}
				checked++;
			}
			assert(Arrays.equals(starts, savedStarts) && Arrays.equals(support, savedSupport)) :
				"Reduction borrows its graph, which must remain usable by other consumers";
		}
		assert(checked>0) : "The finite graph exercise must include feasible reductions";
		return checked;
	}

	/** Checks each deletion against directly counted alternatives and final coverage. */
	private static void checkTrace(final Result result, final int[] starts,
			final int[] edges, final int required, final int maxModelLoss){
		checkTrace(result, starts, edges, required, maxModelLoss, 0);
	}

	/** Recounts support independently while also checking the retained-model floor. */
	private static void checkTrace(final Result result, final int[] starts,
			final int[] edges, final int required, final int maxModelLoss, final int minModels){
		final boolean[] kept={true, true, true};
		int retained=kept.length;
		for(final int victim : result.removed){
			assert(retained>minModels) : "A recorded deletion must leave at least the requested number of models";
			assert(kept[victim]) : "A model may only be removed once";
			final int before=count(starts, edges, kept);
			kept[victim]=false;
			final int after=count(starts, edges, kept);
			kept[victim]=true;
			assert(after>=required) : "Every recorded proposal must retain the target coverage";
			assert(before-after<=maxModelLoss) : "Each deletion must satisfy its cap after preceding drops change unique support";
			for(int other=0; other<kept.length; other++){
				if(!kept[other]){continue;}
				kept[other]=false;
				final int alternative=count(starts, edges, kept);
				kept[other]=true;
				assert(before-after<=before-alternative) : "Chosen deletion must have minimum unique-support cost";
			}
			kept[victim]=false;
			retained--;
		}
		assert(Arrays.equals(kept, result.kept)) : "Final flags must agree with the actual deletion history";
		assert(count(starts, edges, kept)==result.covered) : "Reported coverage must equal independently counted support";
		assert(retained>=minModels) : "Final model count must respect the independent stopping floor";
		if(retained==minModels){return;}
		final int before=count(starts, edges, kept);
		for(int model=0; model<kept.length; model++){
			if(!kept[model]){continue;}
			kept[model]=false;
			final int after=count(starts, edges, kept);
			assert(after<required || before-after>maxModelLoss)
				: "A completed proposal has no next deletion satisfying both total coverage and per-drop loss";
			kept[model]=true;
		}
	}

	/** Counts genes with any kept supporter directly from their input slices. */
	private static int count(final int[] starts, final int[] edges, final boolean[] kept){
		assert(starts[starts.length-1]==edges.length) : "The test oracle requires complete framed input";
		int covered=0;
		for(int gene=0; gene<starts.length-1; gene++){
			for(int edge=starts[gene]; edge<starts[gene+1]; edge++){
				if(kept[edges[edge]]){
					covered++;
					break;
				}
			}
		}
		return covered;
	}
}
