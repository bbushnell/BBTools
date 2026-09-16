package assemble;

import java.util.Arrays;
import shared.Shared;
import ukmer.Kmer;

/** Config dispatch regressions; constructs tables but never loads reads or assembles.
 * @author Fischl
 */
public class TadpoleConfigUnitTest {

	public static void main(String[] args){
		assert(args.length==1) : "Expected BBTools root for checked-in config and input fixtures.";
		final String root=args[0], dir=root+"/testdata/tadpole_config/";
		final String[] original={"config="+dir+"k63.config", "k=31"};
		Shared.COMMAND_LINE=null;
		final String[] expanded=Tadpole.expandConfigArgs(original);
		assert(Arrays.equals(original, Shared.COMMAND_LINE)) : "Config expansion must preserve original invocation metadata.";
		assert(Tadpole.preparseK(expanded)==31) : "Later CLI k must override config at its position.";
		assert(Tadpole.expandConfigArgs(expanded)==expanded) : "Expanded arguments without config should not be copied/reparsed.";
		checkFactory(root, 63, true, "config="+dir+"k63.config");
		checkFactory(root, 31, false, "config="+dir+"k31.config");
		checkFactory(root, 63, true, "config="+dir+"k31.config", "k=63");
		checkFactory(root, 31, false, "k=63", "config="+dir+"k31.config");
		checkFactory(root, 64, true, "config="+dir+"packed.config");
		checkFactory(root, 31, true, "config="+dir+"forced.config");
		Kmer.PACKED=false;
		Tadpole.FORCE_TADPOLE2=false;
		final String[] multi=Tadpole.expandConfigArgs(new String[]{"config="+dir+"multi.config", "out=unused-config-test.fa"});
		assert(TadpoleMulti.hasMultipleK(multi)) : "Config k-list must select multi-k before factory integer parsing.";
		final TadpoleMulti.Config config=new TadpoleMulti.Config(multi);
		assert(config.assembleK==63 && Arrays.equals(config.fuseKs, new int[]{31})) : "Expanded list must retain multi-k phase semantics.";
		//Missing files can terminate the VM in ReadWrite; test that in the CLI subprocess suite.
		System.out.println("PASS Tadpole config dispatch: six factories, precedence, packed/forced, multi-k, metadata");
	}

	private static void checkFactory(final String root, final int k, final boolean longK, final String... options){
		Kmer.PACKED=false;
		Tadpole.FORCE_TADPOLE2=false;
		final String[] args=Arrays.copyOf(options, options.length+6);
		int i=options.length;
		args[i++]="in="+root+"/testdata/crossk_left_bridge_reads.fa";
		args[i++]="out=null"; args[i++]="t=1"; args[i++]="prealloc=f";
		args[i++]="prefilter=f"; args[i++]="pop=f";
		final Tadpole tad=Tadpole.makeTadpole(args, true);
		assert(tad.kbig==k) : "Factory and constructor must agree on k; expected="+k+", actual="+tad.kbig;
		assert((tad instanceof Tadpole2)==longK) : "Wrong representation for config k="+k+", longK="+longK;
	}
}
