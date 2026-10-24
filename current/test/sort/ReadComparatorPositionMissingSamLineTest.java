package test.sort;

import java.util.Comparator;

import sort.ReadComparatorPosition;
import stream.Read;
import stream.SamLine;

/** Baseline and fixed oracle for position sorting roots with missing SamLine metadata.
 * @author Sucrose
 */
public final class ReadComparatorPositionMissingSamLineTest{

	/** Runs expect=baseline or expect=fixed in a fresh test JVM. */
	public static void main(String[] args){
		check(args.length==1, "Expected expect=baseline|fixed");
		check(args[0].startsWith("expect="), "Expected expect=baseline|fixed, got "+args[0]);
		final String expect=args[0].substring("expect=".length());
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline|fixed, got "+args[0]);

		checkMixedMappedUnmappedMissing(expect);
		checkNameTies(expect);
		checkHomogeneousControls(expect);
	}

	/** Checks the approved policy for roots with mapped, unmapped, and missing metadata. */
	private static void checkMixedMappedUnmappedMissing(String expect){
		final Read mapped=read("z-mapped", 1, samLine(0, 100));
		final Read unmapped=read("m-unmapped", 2, new SamLine());
		final Read missing=read("a-missing", 3, null);

		if(fixed(expect)){
			checkOrder(expect, "mixed-ascending", ReadComparatorPosition.ascending,
					mapped, unmapped, missing);
			checkOrder(expect, "mixed-descending", ReadComparatorPosition.descending,
					missing, unmapped, mapped);
		}else{
			checkOrder(expect, "mixed-ascending", ReadComparatorPosition.ascending,
					missing, mapped, unmapped);
			checkOrder(expect, "mixed-descending", ReadComparatorPosition.descending,
					unmapped, mapped, missing);
		}
	}

	/** Checks that position ties and missing-metadata ties still fall back to name. */
	private static void checkNameTies(String expect){
		final Read mappedAlpha=read("alpha", 1, samLine(0, 100));
		final Read mappedBeta=read("beta", 2, samLine(0, 100));
		final Read missingAlpha=read("alpha", 3, null);
		final Read missingBeta=read("beta", 4, null);

		checkOrder(expect, "mapped-name-tie-ascending", ReadComparatorPosition.ascending,
				mappedAlpha, mappedBeta);
		checkOrder(expect, "mapped-name-tie-descending", ReadComparatorPosition.descending,
				mappedBeta, mappedAlpha);
		checkOrder(expect, "missing-name-tie-ascending", ReadComparatorPosition.ascending,
				missingAlpha, missingBeta);
		checkOrder(expect, "missing-name-tie-descending", ReadComparatorPosition.descending,
				missingBeta, missingAlpha);
	}

	/** Checks homogeneous mapped, unmapped, and missing-metadata inputs remain total orders. */
	private static void checkHomogeneousControls(String expect){
		checkOrder(expect, "mapped-homogeneous-ascending", ReadComparatorPosition.ascending,
				read("z-pos300", 1, samLine(0, 300)),
				read("y-pos400", 2, samLine(0, 400)),
				read("x-scaf1", 3, samLine(1, 50)));
		checkOrder(expect, "mapped-homogeneous-descending", ReadComparatorPosition.descending,
				read("x-scaf1", 3, samLine(1, 50)),
				read("y-pos400", 2, samLine(0, 400)),
				read("z-pos300", 1, samLine(0, 300)));

		checkOrder(expect, "unmapped-homogeneous-ascending", ReadComparatorPosition.ascending,
				read("alpha", 1, new SamLine()),
				read("beta", 2, new SamLine()),
				read("gamma", 3, new SamLine()));
		checkOrder(expect, "unmapped-homogeneous-descending", ReadComparatorPosition.descending,
				read("gamma", 3, new SamLine()),
				read("beta", 2, new SamLine()),
				read("alpha", 1, new SamLine()));

		checkOrder(expect, "missing-homogeneous-ascending", ReadComparatorPosition.ascending,
				read("alpha", 1, null),
				read("beta", 2, null),
				read("gamma", 3, null));
		checkOrder(expect, "missing-homogeneous-descending", ReadComparatorPosition.descending,
				read("gamma", 3, null),
				read("beta", 2, null),
				read("alpha", 1, null));
	}

	/** Builds a read with optional SamLine metadata. */
	private static Read read(String id, long numericID, SamLine samLine){
		final Read r=new Read(new byte[] {'A'}, null, id, numericID);
		r.samline=samLine;
		return r;
	}

	/** Builds a mapped SamLine whose public position fields drive the comparator. */
	private static SamLine samLine(int scafnum, int pos){
		final SamLine sl=new SamLine();
		sl.scafnum=scafnum;
		sl.pos=pos;
		return sl;
	}

	/** Requires the listed reads to be in strict comparator order. */
	private static void checkOrder(String expect, String label, Comparator<Read> comparator,
			Read... reads){
		for(int i=0; i<reads.length; i++){
			for(int j=0; j<reads.length; j++){
				checkCompare(expect, label+"-"+i+"-"+j, Integer.compare(i, j),
						comparator, reads[i], reads[j]);
			}
		}
		System.out.println("PASS\t"+expect+"\t"+label);
	}

	/** Checks one comparison sign. */
	private static void checkCompare(String expect, String label, int expected,
			Comparator<Read> comparator, Read a, Read b){
		final int observed=sign(comparator.compare(a, b));
		check(observed==expected, expect+" "+label+" expected "+expected+" observed "+observed);
	}

	/** Returns the sign of a comparator result. */
	private static int sign(int value){return value<0 ? -1 : value>0 ? 1 : 0;}

	/** Returns true when this run should assert the repaired policy. */
	private static boolean fixed(String expect){return "fixed".equals(expect);}

	/** Throws an AssertionError when a test expectation is violated. */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
