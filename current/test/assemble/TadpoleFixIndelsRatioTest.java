package test.assemble;

import java.util.ArrayList;
import java.util.Arrays;

import assemble.Tadpole;
import shared.Shared;
import ukmer.Kmer;

/**
 * Current local-edit parser regression, retaining the historical test name.
 * The ratio API was removed in25457191; verify its rejection rather than
 * referencing the deleted field. No correction algorithm is exercised.
 * Requires enabled assertions for the positive routing cases. Shared settings are
 * changed during construction, so run in a separate JVM with the expected input fixture.
 * @author Fischl
 */
public final class TadpoleFixIndelsRatioTest{

	/** Static test entry point only. */
	private TadpoleFixIndelsRatioTest(){}

	/** Checks three accepted routes under assertions and five rejected option sets.
	 * Constructs small tables without reading or correcting input sequences.
	 * @param args One BBTools root containing the required testdata input
	 */
	public static void main(final String[] args){
		if(args.length!=1){throw new IllegalArgumentException("Expected BBTools root containing testdata.");}
		final String root=args[0];
		assert(!make(root).localEdit) : "Local edits must remain default-off.";
		assert(make(root, "fixindels=t").localEdit) : "Explicit fixindels=t must enable the local edit path.";
		assert(make(root, "fixindels=t", "fixindelspairs=t").localEditPairs) : "Pair lookahead requires enabled local edits.";
		expectFailure(root, "fixindelspairs=t");
		expectFailure(root, "fixindelsmax=0");
		expectFailure(root, "fixindelsstride=0");
		expectFailure(root, "fixindels=t", "fixindelsbatch=t", "fixindelspatchmax=1");
		expectFailure(root, "fixindelsratio=17");
		System.out.println("TADPOLE_FIXINDELS_RATIO_TEST_OK current_routing removed_ratio_rejected");
	}

	/** Resets representation and invocation metadata, then constructs with one thread.
	 * Requires the named fixture path; this helper does not verify its content hash.
	 */
	private static Tadpole make(final String root, final String... extra){
		Kmer.PACKED=false;
		Tadpole.FORCE_TADPOLE2=false;
		Shared.COMMAND_LINE=null;
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(
				"in="+root+"/testdata/crossk_left_bridge_reads.fa", "out=null", "k=31", "t=1",
				"prealloc=f", "prefilter=f", "pop=f", "mode=correct", "ecc=f", "minprob=0"));
		args.addAll(Arrays.asList(extra));
		return Tadpole.makeTadpole(args.toArray(new String[args.size()]), true);
	}

	/** Requires any RuntimeException for an invalid option set; does not inspect its diagnostic.
	 * The positive cases separately establish valid construction when assertions are enabled.
	 */
	private static void expectFailure(final String root, final String... extra){
		try{make(root, extra);}
		catch(final RuntimeException expected){return;}
		throw new AssertionError("Invalid or retired local-edit arguments accepted: "+Arrays.toString(extra));
	}
}
