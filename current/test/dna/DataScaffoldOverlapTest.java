package test.dna;

import dna.Data;

/** Literal overlap oracle with two scaffolds separated by twenty padding bases.
 * Exercises the Data API used by RandomReads3's minimum-overlap filter.
 * @author Sucrose
 */
public final class DataScaffoldOverlapTest{

	/** Runs without reference files, restoring the borrowed static scaffold metadata. */
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No arguments expected");}
		final int[][] oldLocs=Data.scaffoldLocs, oldLengths=Data.scaffoldLengths;
		final int oldPadding=Data.interScaffoldPadding;
		try{
			//Inclusive scaffold spans are [10,12] and [33,37]; [13,32] is padding.
			Data.scaffoldLocs=new int[][]{null, {10, 33}};
			Data.scaffoldLengths=new int[][]{null, {3, 5}};
			Data.interScaffoldPadding=20;
			final int[][] queries={{10, 12}, {12, 13}, {13, 14}, {33, 37}, {37, 38},
					{11, 11}, {34, 35}, {9, 10}, {0, 1}, {20, 21}, {38, 39}, {10, 37}};
			//The final query spans both scaffolds; the selected second scaffold contributes five bases.
			final int[] expected={3, 1, 0, 5, 1, 1, 2, 1, 0, 0, 0, 5}, actual=new int[queries.length];
			for(int i=0; i<queries.length; i++){
				actual[i]=Data.scaffoldOverlapLength(1, queries[i][0], queries[i][1]);
				System.out.println("query=["+queries[i][0]+","+queries[i][1]+
						"] expected="+expected[i]+" actual="+actual[i]);
			}
			for(int i=0; i<queries.length; i++){
				if(actual[i]!=expected[i]){
					throw new AssertionError("Overlap includes padding at query "+i);
				}
			}
			System.out.println("PASS scaffold-overlap");
		}finally{
			Data.scaffoldLocs=oldLocs;
			Data.scaffoldLengths=oldLengths;
			Data.interScaffoldPadding=oldPadding;
		}
	}
}
