package prok;

import java.util.Arrays;

/** Focused configuration fixture for the default-off LSU 3'-only experiment. */
public class CallGenesLsuThreePrimeOnlyConfigTest {

	private static final int[] START={-3,-2,-1,0,1,2};
	private static final int[] STOP={-4,-3,-2,-1,0,1,2,3,4};

	public static void main(String[] args){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED;
		final boolean oldPair=CallGenes.R58LSU_ENABLED;
		final boolean oldBoundary=CallGenes.NCRNA_BOUNDARY_NN_ENABLED;
		final boolean oldLsuBoundary=CallGenes.LSU_BOUNDARY_NN_ENABLED;
		final boolean oldThreePrime=CallGenes.LSU_BOUNDARY_3PRIME_ONLY;
		final float oldMarginStart=CallGenes.LSU_BOUNDARY_MARGIN_START;
		final float oldMarginStop=CallGenes.LSU_BOUNDARY_MARGIN_STOP;
		final String oldNet=CallGenes.LSU_BOUNDARY_NET_OVERRIDE;
		final String oldStart=CallGenes.LSU_BOUNDARY_START_TABLE_OVERRIDE;
		final String oldStop=CallGenes.LSU_BOUNDARY_STOP_TABLE_OVERRIDE;
		try{
			if(CallGenes.LSU_BOUNDARY_3PRIME_ONLY){throw new AssertionError("3' control is not default-off");}
			check("default start", CallGenes.effectiveBoundaryStartOffsets("lsu"), START);
			check("default stop", CallGenes.boundaryStopOffsets("lsu"), STOP);

			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.R58LSU_ENABLED=true;
			CallGenes.NCRNA_BOUNDARY_NN_ENABLED=false;
			CallGenes.LSU_BOUNDARY_NN_ENABLED=true;
			CallGenes.LSU_BOUNDARY_3PRIME_ONLY=true;
			CallGenes.validateNcrnaGateCombo();
			check("3' start", CallGenes.effectiveBoundaryStartOffsets("lsu"), new int[]{0});
			check("3' stop", CallGenes.boundaryStopOffsets("lsu"), STOP);
			check("R58 start unchanged", CallGenes.effectiveBoundaryStartOffsets("r58"), START);

			CallGenes.LSU_BOUNDARY_3PRIME_ONLY=false;
			check("OFF start parity", CallGenes.effectiveBoundaryStartOffsets("lsu"), START);
			CallGenes.LSU_BOUNDARY_NN_ENABLED=false;
			CallGenes.LSU_BOUNDARY_3PRIME_ONLY=true;
			expectInvalid("3-prime without gate");
			CallGenes.LSU_BOUNDARY_3PRIME_ONLY=false;
			CallGenes.LSU_BOUNDARY_NET_OVERRIDE="unused.bbnet";
			expectInvalid("resource override without gate");
			CallGenes.LSU_BOUNDARY_NET_OVERRIDE=null;
			CallGenes.LSU_BOUNDARY_MARGIN_STOP=0.01f;
			expectInvalid("margin without gate");
			CallGenes.LSU_BOUNDARY_MARGIN_STOP=0f;
			try{CallGenes.parseLsuBoundaryMargin("-0.01"); throw new AssertionError("negative margin accepted");}
			catch(IllegalArgumentException expected){/* pass */}
			System.out.println("PASS CallGenesLsuThreePrimeOnlyConfigTest");
		}finally{
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna;
			CallGenes.R58LSU_ENABLED=oldPair;
			CallGenes.NCRNA_BOUNDARY_NN_ENABLED=oldBoundary;
			CallGenes.LSU_BOUNDARY_NN_ENABLED=oldLsuBoundary;
			CallGenes.LSU_BOUNDARY_3PRIME_ONLY=oldThreePrime;
			CallGenes.LSU_BOUNDARY_MARGIN_START=oldMarginStart;
			CallGenes.LSU_BOUNDARY_MARGIN_STOP=oldMarginStop;
			CallGenes.LSU_BOUNDARY_NET_OVERRIDE=oldNet;
			CallGenes.LSU_BOUNDARY_START_TABLE_OVERRIDE=oldStart;
			CallGenes.LSU_BOUNDARY_STOP_TABLE_OVERRIDE=oldStop;
		}
	}

	private static void expectInvalid(String label){
		try{CallGenes.validateNcrnaGateCombo(); throw new AssertionError(label+" was accepted");}
		catch(IllegalArgumentException expected){/* pass */}
	}

	private static void check(String label, int[] actual, int[] expected){
		if(!Arrays.equals(actual, expected)){
			throw new AssertionError(label+" expected "+Arrays.toString(expected)+" found "+Arrays.toString(actual));
		}
	}
}
