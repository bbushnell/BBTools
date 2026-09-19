package assemble;

import java.util.ArrayList;
import java.util.Arrays;

import shared.Shared;
import ukmer.Kmer;

/** Focused argument-routing checks for the public Tadpole magnitude guard.
 * Constructs tables but does not load, correct, or write reads.
 * @author Fischl */
public class TadpoleFixIndelsRatioTest {

	public static void main(final String[] args){
		assert(args.length==1) : "Expected BBTools root containing testdata.";
		final String root=args[0];
		final Tadpole off=makeRaw(root);
		assert(!off.localEdit && off.localEditMagnitudeFactor==17) : "fixindels must remain default-off while retaining the qualified ratio default.";
		final Tadpole defaults=make(root);
		assert(defaults.localEdit && defaults.localEditMagnitudeFactor==17) : "fixindels must default to the qualified factor17 guard.";
		final Tadpole disabled=make(root,"fixindelsratio=0");
		assert(disabled.localEditMagnitudeFactor==0) : "Explicit factor0 must preserve the diagnostic unguarded path.";
		final Tadpole boundary=make(root,"fixindelsratio=2");
		assert(boundary.localEditMagnitudeFactor==2) : "Factor2 is the minimum enabled ratio.";
		final Tadpole explicit=make(root,"fixindelsratio=23");
		assert(explicit.localEditMagnitudeFactor==23) : "Explicit integer ratios >=2 must route unchanged.";
		final Tadpole pairs=make(root,"fixindelspairs=t","fixindelsratio=0");
		assert(pairs.localEditPairs && pairs.localEditMagnitudeFactor==0) : "Pair lookahead remains available only with the guard explicitly disabled.";
		expectFailure(root,"fixindelsratio=1");
		expectFailure(root,"fixindelsratio=-1");
		expectFailure(root,"fixindelspairs=t");
		System.out.println("TADPOLE_FIXINDELS_RATIO_TEST_OK default_off default17 explicit0 boundary2 explicit23 pair_requires0 invalid1 invalid_negative");
	}

	private static Tadpole make(final String root,final String... extra){
		final String[] enabled=new String[extra.length+1];enabled[0]="fixindels=t";System.arraycopy(extra,0,enabled,1,extra.length);
		return makeRaw(root,enabled);
	}

	private static Tadpole makeRaw(final String root,final String... extra){
		Kmer.PACKED=false;Tadpole.FORCE_TADPOLE2=false;Shared.COMMAND_LINE=null;
		final ArrayList<String> list=new ArrayList<String>(Arrays.asList(
			"in="+root+"/testdata/crossk_left_bridge_reads.fa","out=null","k=31","t=1",
			"prealloc=f","prefilter=f","pop=f","mode=correct","ecc=f","minprob=0"));
		list.addAll(Arrays.asList(extra));
		return Tadpole.makeTadpole(list.toArray(new String[list.size()]),true);
	}

	private static void expectFailure(final String root,final String... extra){
		boolean failed=false;
		try{make(root,extra);}catch(IllegalArgumentException expected){failed=true;}
		assert(failed) : "Invalid fixindelsratio combination was accepted: "+Arrays.toString(extra);
	}
}
