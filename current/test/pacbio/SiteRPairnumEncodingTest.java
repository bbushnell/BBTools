package test.pacbio;

import java.util.ArrayList;

import pacbio.SiteR;
import shared.Shared;

/** Literal key and accessor oracles for SiteR's packed numericID/pairnum key.
 * @author Brian Bushnell, Jean
 */
public final class SiteRPairnumEncodingTest{

	/** Runs each finite SiteR key row and reports all observed mismatches.
	 * @param args Unused
	 */
	public static void main(String[] args){
		final ArrayList<String> failures=new ArrayList<String>();
		checkPair(0, 0L, 0L, failures);
		checkPair(1, 0L, Long.MIN_VALUE, failures);
		checkPair(0, 1L, 1L, failures);
		checkPair(1, 1L, -1L, failures);
		checkPair(0, Long.MAX_VALUE, Long.MAX_VALUE, failures);
		checkPair(1, Long.MAX_VALUE, Long.MIN_VALUE+1, failures);
		checkMinusStrand(failures);
		if(!failures.isEmpty()){
			throw new AssertionError(failures.size()+" SiteR pairnum failure(s):\n"+
					join(failures));
		}
	}

	/** Checks key, accessor, and text expansion for one strand-plus SiteR row. */
	private static void checkPair(int pairnum, long numericID, long expectedKey,
			ArrayList<String> failures){
		final SiteR site;
		try{
			site=new SiteR(17, 19, 3, Shared.PLUS, numericID, pairnum);
		}catch(Throwable t){
			failures.add("constructor pair="+pairnum+" id="+numericID+
					" threw "+t.getClass().getSimpleName()+": "+t.getMessage());
			return;
		}
		check(site.idPairnum==expectedKey, "key", pairnum, numericID, expectedKey,
				site.idPairnum, failures);
		check(site.numericID()==numericID, "numericID", pairnum, numericID, numericID,
				site.numericID(), failures);
		check(site.pairNum()==pairnum, "pairnum", pairnum, numericID, pairnum,
				site.pairNum(), failures);
		check(site.chrom()==3, "chrom", pairnum, numericID, 3, site.chrom(), failures);
		check(site.strand()==Shared.PLUS, "strand", pairnum, numericID, Shared.PLUS,
				site.strand(), failures);
		final String expected="17,19,3,"+Shared.PLUS+","+numericID+","+pairnum;
		if(!expected.equals(site.toString())){
			failures.add("text pair="+pairnum+" id="+numericID+
					" expected="+expected+" observed="+site.toString());
		}
	}

	/** Verifies the strand bit still roundtrips independently of mate encoding. */
	private static void checkMinusStrand(ArrayList<String> failures){
		final SiteR site=new SiteR(17, 19, 3, Shared.MINUS, 0, 0);
		if(site.chrom()!=3){failures.add("minus chrom expected=3 observed="+site.chrom());}
		if(site.strand()!=Shared.MINUS){
			failures.add("minus strand expected="+Shared.MINUS+" observed="+site.strand());
		}
	}

	/** Appends one labeled scalar mismatch if a literal equality check fails. */
	private static void check(boolean condition, String label, int pairnum, long numericID,
			long expected, long observed, ArrayList<String> failures){
		if(!condition){
			failures.add(label+" pair="+pairnum+" id="+numericID+
					" expected="+expected+" observed="+observed);
		}
	}

	/** Joins failure messages for one aggregated AssertionError. */
	private static String join(ArrayList<String> failures){
		StringBuilder sb=new StringBuilder();
		for(String s : failures){sb.append(s).append('\n');}
		return sb.toString();
	}

	/** Static-only fixture. */
	private SiteRPairnumEncodingTest(){}
}
