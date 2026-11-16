package consensus;

/** Bridges nucleotide model codecs to graph counts and probability initialization.
 * These operations are unsynchronized. Load counts into a fresh graph, initialize
 * probabilities once, then publish it to scoring readers. Counts do not include
 * insertion chains or graph-level read/quality totals.
 * @author Brian Bushnell
 */
public class BaseGraphHelper{

	/** Initializes reference/insertion probability caches after all counts are loaded.
	 * @param bg Fully populated graph whose probability caches are still absent
	 */
	public static void initForScoring(BaseGraph bg){
		bg.calcProbs();
	}

	/** Replaces one reference position's four base counts and deletion support.
	 * Sets allele weights equal to counts for imported unweighted models. Does not
	 * update insertion chains, graph totals or existing probability caches.
	 * @param bg Fresh graph being populated, before probability initialization
	 * @param pos Zero-based reference position, including any graph padding
	 * @param a Nonnegative A observations
	 * @param c Nonnegative C observations
	 * @param g Nonnegative G observations
	 * @param t Nonnegative T observations; the four-base total must fit an int
	 * @param del Nonnegative deletion observations
	 */
	public static void setCounts(BaseGraph bg, int pos, int a, int c, int g, int t, int del){
		BaseNode rnode=bg.ref[pos];
		rnode.acgtCount[0]=a;
		rnode.acgtCount[1]=c;
		rnode.acgtCount[2]=g;
		rnode.acgtCount[3]=t;
		rnode.countSum=a+c+g+t;
		rnode.weightSum=a+c+g+t;
		rnode.acgtWeight[0]=a;
		rnode.acgtWeight[1]=c;
		rnode.acgtWeight[2]=g;
		rnode.acgtWeight[3]=t;
		bg.del[pos].countSum=del;
		bg.del[pos].weightSum=del;
	}

	/** Copies the A/C/G/T/deletion counts at one reference coordinate.
	 * @param bg Graph whose accumulation or loading has finished
	 * @param pos Zero-based reference position, including padding
	 * @return A new five-element array in A,C,G,T,deletion order
	 */
	public static int[] getCounts(BaseGraph bg, int pos){
		BaseNode rnode=bg.ref[pos];
		int del=bg.del[pos].countSum;
		return new int[]{rnode.acgtCount[0], rnode.acgtCount[1],
			rnode.acgtCount[2], rnode.acgtCount[3], del};
	}

	/** Returns the reference-backbone length, including padding and excluding insertions. */
	public static int length(BaseGraph bg){
		return bg.ref.length;
	}
}
