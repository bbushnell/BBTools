package test.prok;

import java.util.Arrays;
import prok.NcrnaLearningSplit;
import prok.NcrnaLearningSplit.Result;

/** Finite population/quota cases, independent of a particular shuffled membership. */
public final class NcrnaLearningSplitTest {
	public static void main(String[] args){
		assert(args.length==0) : "Split fixture has no external population inputs";
		final int[] taxa=new int[61], phyla=new int[61], models=new int[61];
		for(int i=0; i<60; i++){taxa[i]=i/3+1;}
		// An additional copy puts TaxID1 in a rare bin without changing the group target.
		taxa[60]=1;
		models[60]=64;
		final int[] original=taxa.clone();
		final Result result=NcrnaLearningSplit.split(taxa, phyla, models, 1, 65, 42);
		assert(result.target==2 && result.selected==2 && !result.validationByGroup[0])
			: "Twenty groups permit two holdouts, but the unique model64 group must stay in training";
		for(int i=0; i<60; i+=3){
			assert(result.validationByGene[i]==result.validationByGene[i+1]
				&& result.validationByGene[i]==result.validationByGene[i+2]) : "Copies of one TaxID must never straddle the split";
		}
		assert(Arrays.equals(taxa, original)) : "Partitioning must not mutate source identities";
		final int[] reversed=reverse(taxa), reversedModels=reverse(models);
		final Result reordered=NcrnaLearningSplit.split(reversed, phyla, reversedModels, 1, 65, 42);
		assert(Arrays.equals(result.validationByGroup, reordered.validationByGroup))
			: "A fixed seed is applied to sorted TaxIDs, independent of input record order";

		final int[] t30=new int[30], p30=new int[30], m30=new int[30];
		for(int i=0; i<30; i++){
			t30[i]=i+1;
			m30[i]=i/15;
		}
		final Result underfilled=NcrnaLearningSplit.split(t30, p30, m30, 1, 2, 11);
		assert(underfilled.target==3 && underfilled.selected==2)
			: "Two bins with fifteen groups each permit only one holdout each; do not force the global target";
		final int[] pRare=p30.clone();
		pRare[29]=1;
		final Result rare=NcrnaLearningSplit.split(t30, pRare, new int[30], 2, 1, 7);
		assert(!rare.validationByGroup[29] && rare.selected==2)
			: "The singleton phylum stays in training; the other twenty-nine groups permit two holdouts";
		final int[] conflicting=t30.clone(), conflictPhyla=new int[30];
		conflicting[29]=1;
		conflictPhyla[29]=1;
		fails(()->NcrnaLearningSplit.split(conflicting, conflictPhyla, new int[30], 2, 1, 0));
		final int[] missingTaxid=t30.clone();
		missingTaxid[0]=0;
		fails(()->NcrnaLearningSplit.split(missingTaxid, p30, new int[30], 1, 1, 0));
		for(int i=0; i<30; i++){m30[i]=i;}
		fails(()->NcrnaLearningSplit.split(t30, p30, m30, 1, 30, 0));
		System.out.println("Learning split passed: whole TaxIDs, rare phyla/bins, word boundary, underfill and invalid populations");
	}

	private static int[] reverse(int[] input){
		assert(input!=null) : "Order-control fixture requires an existing population";
		final int[] output=new int[input.length];
		for(int i=0; i<input.length; i++){output[i]=input[input.length-1-i];}
		return output;
	}

	private static void fails(Runnable action){
		try{
			action.run();
		}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid population must fail instead of silently relaxing protection limits");
	}
}
