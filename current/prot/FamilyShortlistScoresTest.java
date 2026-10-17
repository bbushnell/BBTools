package prot;

import java.util.Arrays;

/** Exact small arithmetic oracles for cross-pack normalization and tie ordering.
 * Does not load real models or establish biological annotation accuracy.
 * @author Keqing
 */
public final class FamilyShortlistScoresTest {
	public static void main(String[] args){
		final int[][] k={{0, 2}, {0, 2}};
		final float[][] c={{0, 0}, {1, 1}}, out={new float[2], new float[2]};
		final double[][] length={{100, 100}, {100, 100}}, sd={{0, 0}, {0, 0}};
		FamilyShortlistScores.score(k, c, length, sd, 100, out);
		//Global k mean/sd=1/1; composition mean/sd=.5/.5. No length penalty.
		check(Arrays.equals(out[0], new float[]{-3, -1})
			&& Arrays.equals(out[1], new float[]{1, 3}), "Closed-form global F4 oracle differs");
		check(Arrays.equals(c[0], new float[]{0, 0}) && Arrays.equals(c[1], new float[]{1, 1}), "Raw composition inputs changed");
		final int[] packs=new int[3], ranks=new int[3];
		FamilyShortlistScores.top(out, packs, ranks);
		check(Arrays.equals(packs, new int[]{1, 1, 0}) && Arrays.equals(ranks, new int[]{1, 0, 1}), "Global candidate ordering differs");
		//Local normalization erases the between-pack composition difference.
		final float[][] local={new float[2], new float[2]};
		for(int p=0; p<2; p++){
			FamilyShortlistScores.score(new int[][]{k[p]}, new float[][]{c[p]},
				new double[][]{length[p]}, new double[][]{sd[p]}, 100, new float[][]{local[p]});
		}
		final int[] onePack=new int[1], oneRank=new int[1];
		FamilyShortlistScores.top(local, onePack, oneRank);
		check(onePack[0]==0 && oneRank[0]==1, "Counterexample must expose the wrong pack-local winner");
		final float[][] flat={new float[4]}, split={new float[3], new float[1]};
		FamilyShortlistScores.score(new int[][]{{0, 2, 0, 2}}, new float[][]{{0, 0, 1, 1}},
			new double[][]{{100, 100, 100, 100}}, new double[][]{{0, 0, 0, 0}}, 100, flat);
		FamilyShortlistScores.score(new int[][]{{0, 2, 0}, {2}}, new float[][]{{0, 0, 1}, {1}},
			new double[][]{{100, 100, 100}, {100}}, new double[][]{{0, 0, 0}, {0}}, 100, split);
		check(Arrays.equals(flat[0], new float[]{-3, -1, 1, 3})
			&& Arrays.equals(split[0], new float[]{-3, -1, 1}) && split[1][0]==3,
			"Changing pack boundaries must preserve canonical-order F4 values");
		FamilyShortlistScores.top(new float[][]{{7, 7}, {7}}, packs, ranks);
		check(Arrays.equals(packs, new int[]{0, 0, 1}) && Arrays.equals(ranks, new int[]{0, 1, 0}), "F4 ties must follow canonical pack/local order");
		final float[][] constant={new float[1]};
		FamilyShortlistScores.score(new int[][]{{0}}, new float[][]{{0}},
			new double[][]{{100}}, new double[][]{{0}}, 100, constant);
		check(constant[0][0]==0, "Zero feature variances and equal length must score zero");
		final float[][] allNan={new float[2], new float[1]};
		final float[][] allFinite={new float[2], new float[1]};
		FamilyShortlistScores.score(new int[][]{{0, 2}, {4}},
			new float[][]{{Float.NaN, Float.NaN}, {Float.NaN}},
			new double[][]{{100, 100}, {100}}, new double[][]{{0, 0}, {0}}, 100, allNan);
		FamilyShortlistScores.score(new int[][]{{0, 2}, {4}},
			new float[][]{{7, 7}, {7}},
			new double[][]{{100, 100}, {100}}, new double[][]{{0, 0}, {0}}, 100, allFinite);
		final float expected=(float)Math.sqrt(1.5);
		check(Arrays.equals(allNan[0], new float[]{-expected, 0}) && allNan[1][0]==expected
			&& Arrays.equals(allNan[0], allFinite[0]) && Arrays.equals(allNan[1], allFinite[1]),
			"All-NaN composition must act like a constant composition term");
		reject(()->FamilyShortlistScores.score(k, c, length, sd, 0, out));
		reject(()->FamilyShortlistScores.score(k, c, length, sd, 100, c));
		reject(()->FamilyShortlistScores.top(new float[][]{{Float.NaN}}, onePack, oneRank));
		reject(()->FamilyShortlistScores.top(new float[][]{{1}}, packs, ranks));
		reject(()->FamilyShortlistScores.score(new int[][]{{-1}}, new float[][]{{0}},
			new double[][]{{100}}, new double[][]{{0}}, 100, constant));
		reject(()->FamilyShortlistScores.score(new int[][]{{0}}, new float[][]{{Float.POSITIVE_INFINITY}},
			new double[][]{{100}}, new double[][]{{0}}, 100, constant));
		reject(()->FamilyShortlistScores.score(new int[][]{{0, 1}}, new float[][]{{Float.NaN, 0}},
			new double[][]{{100, 100}}, new double[][]{{0, 0}}, 100, new float[][]{new float[2]}));
		reject(()->FamilyShortlistScores.score(new int[][]{{0}, {1}}, new float[][]{{Float.NaN}, {0}},
			new double[][]{{100}, {100}}, new double[][]{{0}, {0}}, 100, new float[][]{new float[1], new float[1]}));
		System.err.println("FAMILY_SHORTLIST_SCORES_PASS closed_form=true local_normalization_counterexample=true partition_invariant=true ties=true all_nan_composition=true rejection_cases=8");
	}
	private static void reject(Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Expected invalid global shortlist input to fail loudly");
	}
	private static void check(boolean ok, String reason){if(!ok){throw new AssertionError(reason);}}
}
