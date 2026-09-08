package prok;

/** Focused guards for family-specific ncRNA candidate-padding overrides. */
public class CallGenesNcrnaFamilyPadConfigTest {

	public static void main(String[] args){
		testPadOverridesRequireNcrnaGate();
		System.out.println("PASS CallGenesNcrnaFamilyPadConfigTest");
	}

	private static void testPadOverridesRequireNcrnaGate(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED;
		final int oldRnasep=CallGenes.RNASEP_PAD_OVERRIDE;
		final int oldSrpSmall=CallGenes.SRPSMALL_PAD_OVERRIDE;
		final int oldSrpLarge=CallGenes.SRPLARGE_PAD_OVERRIDE;
		try{
			CallGenes.NCRNA_FAMILIES_ENABLED=false;
			assertRejected("rnasep", new Runnable(){@Override public void run(){CallGenes.RNASEP_PAD_OVERRIDE=1;}});
			CallGenes.RNASEP_PAD_OVERRIDE=-1;
			assertRejected("srp_small", new Runnable(){@Override public void run(){CallGenes.SRPSMALL_PAD_OVERRIDE=1;}});
			CallGenes.SRPSMALL_PAD_OVERRIDE=-1;
			assertRejected("srp_large", new Runnable(){@Override public void run(){CallGenes.SRPLARGE_PAD_OVERRIDE=1;}});
			CallGenes.SRPLARGE_PAD_OVERRIDE=-1;
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.RNASEP_PAD_OVERRIDE=1;
			CallGenes.SRPSMALL_PAD_OVERRIDE=1;
			CallGenes.SRPLARGE_PAD_OVERRIDE=1;
			CallGenes.validateNcrnaSweepOverrides();
		}finally{
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna;
			CallGenes.RNASEP_PAD_OVERRIDE=oldRnasep;
			CallGenes.SRPSMALL_PAD_OVERRIDE=oldSrpSmall;
			CallGenes.SRPLARGE_PAD_OVERRIDE=oldSrpLarge;
		}
	}

	private static void assertRejected(String name, Runnable setter){
		setter.run();
		try{
			CallGenes.validateNcrnaSweepOverrides();
			throw new AssertionError(name+" pad override without ncrna=t was accepted");
		}catch(IllegalArgumentException expected){/* pass */}
	}
}
