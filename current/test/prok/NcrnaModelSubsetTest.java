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
				checkTrace(result, starts, support, required);
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
			final int[] edges, final int required){
		final boolean[] kept={true, true, true};
		for(final int victim : result.removed){
			assert(kept[victim]) : "A model may only be removed once";
			final int before=count(starts, edges, kept);
			kept[victim]=false;
			final int after=count(starts, edges, kept);
			kept[victim]=true;
			assert(after>=required) : "Every recorded proposal must retain the target coverage";
			for(int other=0; other<kept.length; other++){
				if(!kept[other]){continue;}
				kept[other]=false;
				final int alternative=count(starts, edges, kept);
				kept[other]=true;
				assert(before-after<=before-alternative) : "Chosen deletion must have minimum unique-support cost";
			}
			kept[victim]=false;
		}
		assert(Arrays.equals(kept, result.kept)) : "Final flags must agree with the actual deletion history";
		assert(count(starts, edges, kept)==result.covered) : "Reported coverage must equal independently counted support";
		for(int model=0; model<kept.length; model++){
			if(!kept[model]){continue;}
			kept[model]=false;
			assert(count(starts, edges, kept)<required) : "A completed greedy proposal must have no affordable next deletion";
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
