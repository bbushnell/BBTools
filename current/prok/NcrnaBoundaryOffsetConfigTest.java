package prok;

import java.util.Arrays;

/** Verifies the family-specific boundary candidate configuration loaded by CallGenes. */
public class NcrnaBoundaryOffsetConfigTest {

	public static void main(String[] args){
		check("rnasep", new int[]{-3,-2,-1,0,1,2}, new int[]{-1,0,1,2,3,4});
		check("srp_small", new int[]{-3,-2,-1,0,1,2}, new int[]{-2,-1,0,1,2,3});
		check("srp_large", new int[]{-3,-2,-1,0,1,2}, new int[]{-2,-1,0,1,2,3});
		testFactoryReturnsFreshArrays();
		testFamilyClonesArrays();
		testUnknownFamilyFailsLoud();
		System.out.println("PASS NcrnaBoundaryOffsetConfigTest");
	}

	private static void check(String name, int[] expectedStart, int[] expectedStop){
		final int[] start=CallGenes.boundaryStartOffsets(name);
		final int[] stop=CallGenes.boundaryStopOffsets(name);
		if(!Arrays.equals(expectedStart, start)){
			throw new AssertionError(name+" start offsets: "+Arrays.toString(start));
		}
		if(!Arrays.equals(expectedStop, stop)){
			throw new AssertionError(name+" stop offsets: "+Arrays.toString(stop));
		}
	}

	private static void testFactoryReturnsFreshArrays(){
		final int[] first=CallGenes.boundaryStopOffsets("rnasep");
		first[0]=99;
		if(CallGenes.boundaryStopOffsets("rnasep")[0]!=-1){
			throw new AssertionError("Offset factory returned shared mutable state");
		}
	}

	private static void testFamilyClonesArrays(){
		final int[] start={-3,-2,-1,0,1,2}, stop={-1,0,1,2,3,4};
		final NcrnaFamily family=new NcrnaFamily("test", null, null, null, null, 17, 1, 1,
			7, 60, true, 1f, 1f, 1f, 1, 0f, 1f, 0.7f, 0.6f, 0.7f, 0.85f,
			null, null, null, null, -1, -1, -1, -1, 0f, start, stop);
		start[0]=99; stop[0]=99;
		if(family.boundaryStartOffsets[0]!=-3 || family.boundaryStopOffsets[0]!=-1){
			throw new AssertionError("NcrnaFamily retained caller-owned mutable arrays");
		}
	}

	private static void testUnknownFamilyFailsLoud(){
		try{
			CallGenes.boundaryStartOffsets("future_family");
			throw new AssertionError("Unknown family silently received offsets");
		}catch(IllegalArgumentException expected){/* pass */}
	}
}
