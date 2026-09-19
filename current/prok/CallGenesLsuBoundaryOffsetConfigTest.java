package prok;

import java.util.Arrays;

/** Focused configuration check for the dormant LSU boundary candidate range. */
public class CallGenesLsuBoundaryOffsetConfigTest {
	public static void main(String[] args){
		check("LSU start", CallGenes.boundaryStartOffsets("lsu"), new int[]{-3,-2,-1,0,1,2});
		check("LSU stop", CallGenes.boundaryStopOffsets("lsu"), new int[]{-4,-3,-2,-1,0,1,2,3,4});
		check("R58 stop unchanged", CallGenes.boundaryStopOffsets("r58"), new int[]{-3,-2,-1,0,1,2});
		System.out.println("PASS CallGenesLsuBoundaryOffsetConfigTest");
	}
	private static void check(String label, int[] actual, int[] expected){
		if(!Arrays.equals(actual, expected)){
			throw new AssertionError(label+" expected "+Arrays.toString(expected)+" found "+Arrays.toString(actual));
		}
	}
}
