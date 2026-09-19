package prok;

import java.util.Arrays;

/** Focused configuration fixture for the default-off R58 3'-only experiment. */
public class CallGenesR58ThreePrimeOnlyConfigTest {

	private static final int[] SIX_POINT={-3,-2,-1,0,1,2};

	public static void main(String[] args){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED;
		final boolean oldPair=CallGenes.R58LSU_ENABLED;
		final boolean oldBoundary=CallGenes.NCRNA_BOUNDARY_NN_ENABLED;
		final boolean oldR58Boundary=CallGenes.R58_BOUNDARY_NN_ENABLED;
		final boolean oldThreePrime=CallGenes.R58_BOUNDARY_3PRIME_ONLY;
		final float oldMarginStart=CallGenes.R58_BOUNDARY_MARGIN_START;
		final float oldMarginStop=CallGenes.R58_BOUNDARY_MARGIN_STOP;
		final String oldNet=CallGenes.R58_BOUNDARY_NET_OVERRIDE;
		final String oldStart=CallGenes.R58_BOUNDARY_START_TABLE_OVERRIDE;
		final String oldStop=CallGenes.R58_BOUNDARY_STOP_TABLE_OVERRIDE;
		try{
			if(CallGenes.R58_BOUNDARY_3PRIME_ONLY){throw new AssertionError("3' control is not default-off");}
			check("default start", CallGenes.effectiveBoundaryStartOffsets("r58"), SIX_POINT);
			check("default stop", CallGenes.boundaryStopOffsets("r58"), SIX_POINT);

			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.R58LSU_ENABLED=true;
			CallGenes.NCRNA_BOUNDARY_NN_ENABLED=false;
			CallGenes.R58_BOUNDARY_NN_ENABLED=true;
			CallGenes.R58_BOUNDARY_3PRIME_ONLY=true;
			CallGenes.validateNcrnaGateCombo();
			check("3' start", CallGenes.effectiveBoundaryStartOffsets("r58"), new int[]{0});
			check("3' stop unchanged", CallGenes.boundaryStopOffsets("r58"), SIX_POINT);
			check("non-R58 start unchanged", CallGenes.effectiveBoundaryStartOffsets("lsu"), SIX_POINT);

			CallGenes.R58_BOUNDARY_3PRIME_ONLY=false;
			check("OFF start parity", CallGenes.effectiveBoundaryStartOffsets("r58"), SIX_POINT);
			check("OFF stop parity", CallGenes.boundaryStopOffsets("r58"), SIX_POINT);

			CallGenes.R58_BOUNDARY_NN_ENABLED=false;
			CallGenes.R58_BOUNDARY_3PRIME_ONLY=true;
			expectInvalid("3-prime without gate");
			CallGenes.R58_BOUNDARY_3PRIME_ONLY=false;
			CallGenes.R58_BOUNDARY_NET_OVERRIDE="unused.bbnet";
			expectInvalid("resource override without gate");
			CallGenes.R58_BOUNDARY_NET_OVERRIDE=null;
			CallGenes.R58_BOUNDARY_MARGIN_STOP=0.01f;
			expectInvalid("margin without gate");
			CallGenes.R58_BOUNDARY_MARGIN_STOP=0f;
			try{CallGenes.parseR58BoundaryMargin("-0.01"); throw new AssertionError("negative margin accepted");}
			catch(IllegalArgumentException expected){/* pass */}
			System.out.println("PASS CallGenesR58ThreePrimeOnlyConfigTest");
		}finally{
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna;
			CallGenes.R58LSU_ENABLED=oldPair;
			CallGenes.NCRNA_BOUNDARY_NN_ENABLED=oldBoundary;
			CallGenes.R58_BOUNDARY_NN_ENABLED=oldR58Boundary;
			CallGenes.R58_BOUNDARY_3PRIME_ONLY=oldThreePrime;
			CallGenes.R58_BOUNDARY_MARGIN_START=oldMarginStart;
			CallGenes.R58_BOUNDARY_MARGIN_STOP=oldMarginStop;
			CallGenes.R58_BOUNDARY_NET_OVERRIDE=oldNet;
			CallGenes.R58_BOUNDARY_START_TABLE_OVERRIDE=oldStart;
			CallGenes.R58_BOUNDARY_STOP_TABLE_OVERRIDE=oldStop;
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
