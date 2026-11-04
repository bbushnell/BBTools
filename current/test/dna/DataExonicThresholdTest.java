package test.dna;

import dna.Data;
import driver.Search;
import structures.Range;

/** Literal threshold-overlap checks for Search and Data's non-coding exon wrapper.
 * Uses at most three ranges and restores the public exon-cache slot in finally.
 * No annotation resources, reflection, coding-cache access or large arrays are used.
 * Expected booleans come from explicit inclusive coordinates, not production search.
 * @author Sucrose
 */
public final class DataExonicThresholdTest{

	/** Runs case=controls|next-only|bounds in a fresh JVM under either assertion mode. */
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("case="), "Require one case= argument");
		final String mode=args[0].substring(5);
		require("controls".equals(mode) || "next-only".equals(mode) || "bounds".equals(mode),
				"Unknown threshold case: "+mode);
		final Range[] prior=Data.exonRangeMatrix[CHROM];
		try{
			final Range[] ranges={new Range(10, 20), new Range(40, 50), new Range(70, 80)};
			if("controls".equals(mode)){controls(ranges);}
			else if("next-only".equals(mode)){
				require(check("reaches-next-only", ranges, 38, 2, true),
						"[36,40] reaches the next range [40,50], not the preceding range [10,20]");
			}else{bounds();}
		}finally{
			Data.exonRangeMatrix[CHROM]=prior;
		}
		System.out.println("PASS "+mode);
	}

	/** Preserves empty, exact, gap, preceding-range and outer-boundary behavior. */
	private static void controls(Range[] ranges){
		boolean ok=true;
		ok&=check("empty", new Range[0], 30, 3, false);
		ok&=check("exact-hit", ranges, 15, 0, true);
		ok&=check("exact-right-start", ranges, 40, 0, true);
		ok&=check("gap-miss", ranges, 30, 5, false);
		ok&=check("reaches-preceding-only", ranges, 22, 2, true);
		ok&=check("reaches-both", ranges, 30, 10, true);
		ok&=check("before-first-reaches", ranges, 9, 1, true);
		ok&=check("before-first-misses", ranges, 8, 1, false);
		ok&=check("after-last-reaches", ranges, 81, 1, true);
		ok&=check("after-last-misses", ranges, 82, 1, false);
		ok&=check("maximum-exact", new Range[]{new Range(Integer.MAX_VALUE, Integer.MAX_VALUE)},
				Integer.MAX_VALUE, 0, true);
		ok&=check("minimum-exact", new Range[]{new Range(Integer.MIN_VALUE, Integer.MIN_VALUE)},
				Integer.MIN_VALUE, 0, true);
		require(ok, "Finite threshold controls changed");
	}

	/** Mathematical intervals extend past int bounds but still include the supplied extreme endpoint. */
	private static void bounds(){
		boolean ok=true;
		ok&=check("upper-overflow-reaches", new Range[]{new Range(Integer.MAX_VALUE, Integer.MAX_VALUE)},
				Integer.MAX_VALUE-1, 2, true);
		ok&=check("lower-overflow-reaches", new Range[]{new Range(Integer.MIN_VALUE, Integer.MIN_VALUE)},
				Integer.MIN_VALUE+1, 2, true);
		require(ok, "Threshold arithmetic must not wrap away from a covered integer endpoint");
	}

	/** Runs both public paths against the same independent expected answer and reports both observations. */
	private static boolean check(String label, Range[] ranges, int point, int threshold, boolean expected){
		final boolean direct=Search.containsPointBinary(point, ranges, threshold);
		Data.exonRangeMatrix[CHROM]=ranges;
		final boolean wrapper=Data.isExonic(CHROM, point, threshold, false);
		final boolean ok=direct==expected && wrapper==expected;
		System.out.println(label+" point="+point+" threshold="+threshold+" expected="+expected+
				" direct="+direct+" wrapper="+wrapper+" pass="+ok);
		return ok;
	}

	/** Keeps the oracle active when Java assertions are disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Private test use of the public non-coding exon slot, restored even after an assertion failure. */
	private static final byte CHROM=1;
}
