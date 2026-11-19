package test.prot;

import java.lang.reflect.Method;
import java.util.Arrays;

/** Independent finite full-sort oracle for the public covering selector.
 * Uses primitive score fixtures without files, digests or biological corpora.
 * Does not load the private historical throughput benchmark.
 * @author Sayu, Fischl
 */
public final class HbmTailCoveringSelectionTest{

	/*--------------------------------------------------------------*/
	/*----------------       Finite oracle checks   ----------------*/
	/*--------------------------------------------------------------*/

	/** Enumerates lengths zero through six over scores -1, 0, 1 and 2.
	 * Tests requested widths through length+2, then signed integer boundaries.
	 * No arguments, input files or external resources are required.
	 */
	public static void main(String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("No arguments");}
		final Method product=selector("prot.HbmTailPairedShortlist");
		for(int length=0; length<=6; length++){
			final int patterns=1<<(2*length);
			for(int pattern=0; pattern<patterns; pattern++){
				final int[] scores=new int[length];
				for(int i=0; i<length; i++){
					scores[i]=((pattern>>>(2*i))&3)-1;
				}
				for(int width=1; width<=length+2; width++){
					verify(scores, width, product);
				}
			}
		}
		for(int width=1; width<=8; width++){
			verify(new int[]{Integer.MIN_VALUE, 0, Integer.MAX_VALUE, 1, Integer.MAX_VALUE, -1},
					width, product);
		}
		System.out.println("HBM_TAIL_COVERING_SELECTION_TEST_PASS cases="+cases+
				" files=0 hashes=0 product_oracle=true private_dependency=false");
	}

	/** Resolves the product's private static selector with its two-array signature.
	 * Reflection permits direct comparison without adding a public test-only API.
	 */
	private static Method selector(String name) throws Exception{
		final Method method=Class.forName(name).getDeclaredMethod("topCovering", int[].class, int[].class);
		method.setAccessible(true);
		return method;
	}

	/** Compares with a positive-only full sort: score descending, then index ascending.
	 * The product receives a clone, leaving the supplied scores as the immutable
	 * oracle input. Checks returned count, selected prefix and unchanged scores;
	 * the unused output suffix is deliberately outside this test's contract.
	 */
	private static void verify(final int[] scores, int width, Method product) throws Exception{
		int positives=0;
		for(int score:scores){if(score>0){positives++;}}
		final Integer[] ordered=new Integer[positives];
		for(int i=0, at=0; i<scores.length; i++){
			if(scores[i]>0){ordered[at++]=i;}
		}
		Arrays.sort(ordered, (a, b)->{
			final int scoreOrder=Integer.compare(scores[b], scores[a]);
			return scoreOrder!=0 ? scoreOrder : Integer.compare(a, b);
		});
		final int expected=Math.min(width, positives);
		final int[] productScores=scores.clone();
		final int[] actual=new int[width];
		Arrays.fill(actual, -37);
		int got=-1;
		try{
			got=(Integer)product.invoke(null, productScores, actual);
		}catch(ReflectiveOperationException failure){
			throw new AssertionError("Selector invocation failed"+
					details(scores, width, expected, got, ordered, actual), failure);
		}
		if(!Arrays.equals(scores, productScores)){
			throw new AssertionError("Selector changed score input"+
					details(scores, width, expected, got, ordered, actual)+
					" productScores="+Arrays.toString(productScores));
		}
		if(got!=expected){
			throw new AssertionError("Independent selector count differs"+
					details(scores, width, expected, got, ordered, actual));
		}
		for(int i=0; i<expected; i++){
			if(actual[i]!=ordered[i]){
				throw new AssertionError("Independent score/roster order failed"+
						details(scores, width, expected, got, ordered, actual));
			}
		}
		cases++;
	}

	/** Retains the complete finite witness when a selector assertion fails. */
	private static String details(int[] scores, int width, int expected, int got, Integer[] ordered, int[] actual){
		return " scores="+Arrays.toString(scores)+" width="+width+
				" expectedCount="+expected+" productCount="+got+
				" expectedOrder="+Arrays.toString(ordered)+
				" productRanks="+Arrays.toString(actual);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Completed comparisons in this JVM; the standalone launcher starts a fresh process. */
	private static long cases;
}
