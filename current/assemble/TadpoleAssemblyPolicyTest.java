package assemble;

import java.util.Arrays;

/** Tests public-default selection separately from expensive biological validation.
 * @author Fischl
 */
public final class TadpoleAssemblyPolicyTest {

	/** Checks mode routing, explicit overrides, resource resolution and guard parsing. */
	public static void main(final String[] args){
		assert(args!=null) : "Run with assertions enabled through the test launcher.";
		check(TadpoleAssemblyPolicy.isAssembly(new String[]{"k=96"}), "ordinary assembly");
		for(String option : new String[]{"mode=correct", "mode=extend", "mode=insert", "mode=graph",
				"mode=5", "ecc", "ecct=t", "aecc", "conservative", "markerrors", "hpindel",
				"fixindels", "localedit", "tossjunk", "tu", "tossdepth=1", "el=10", "er=10",
				"maxextension=10", "contigs=assembly.fa"}){
			final TadpoleAssemblyPolicy policy=new TadpoleAssemblyPolicy(new String[]{option});
			check(!policy.assembly && !policy.rescueRequested, option);
			check(Arrays.equals(policy.args, new String[]{option}), "unchanged nonassembly arguments: "+option);
		}
		for(String[] options : new String[][]{{"ecc=t", "ecc=f"}, {"el=10", "el=0"},
				{"mode=correct", "mode=contig"}, {"mode=0"}, {"mode=00"}, {"mode=0k"}, {"tossdepth=0"}}){
			check(TadpoleAssemblyPolicy.isAssembly(options), "last value/explicit mode: "+Arrays.toString(options));
		}
		final String[] original={"k=96,124,64,32", "in=reads.fq", "out=contigs.fa"};
		final TadpoleAssemblyPolicy defaults=new TadpoleAssemblyPolicy(original);
		check(defaults.assembly && defaults.multi && defaults.rescueRequested, "default neural assembly");
		final TadpoleMulti.Config config=defaults.multiConfig();
		check(config.assembleK==96 && config.graphK==32, "user-controlled K/order");
		check(Arrays.equals(config.phaseKs, new int[]{96,124,64,32}), "four phases, no appended K");
		check(config.fuseNet.endsWith("tadpole_fusion.bbnet") && config.fuseCutoff==0.667098f, "bundled fusion model");
		check(config.fuseCoverageRatio==1.75f && config.graphMergeCoverageRatio==1.75f, "depth safeguards");
		check(config.fuseShortEvidence && config.graphMergeBranchGuard, "ambiguity safeguards");
		for(int k : new int[]{24,32,64,96,124,127}){check(config.neuralKAllowed(k), "bundled fusion domain "+k);}
		for(int k : new int[]{1,23,128,300}){check(!config.neuralKAllowed(k), "bundled fusion abstention "+k);}
		final TadpoleAssemblyPolicy disabled=new TadpoleAssemblyPolicy(new String[]{
				"assemblydefaults=f", "k=96,124,64,32", "in=reads.fq", "out=contigs.fa"});
		check(Arrays.equals(original, disabled.args) && !disabled.rescueRequested, "legacy defaults opt-out");
		check(disabled.multiConfig().neuralKAllowed(300), "explicit caller model domain remains caller-controlled");
		for(boolean first : new boolean[]{true, false}){
			final String[] overrides={"k=96,124,64,32", "out=a.fa", "fusenet=null", "branchrescue=f",
					"fusecoverageratio=0", "graphmergecoverageratio=0", "graphbranchguard=f",
					"fuseshortevidence=f", "minprob=0.8", "retainshortcontigs=f", "trimcircular=t"};
			final String[] input=new String[overrides.length+1];
			System.arraycopy(overrides, 0, input, first ? 1 : 0, overrides.length);
			input[first ? 0 : overrides.length]="assemblydefaults=t";
			final TadpoleAssemblyPolicy policy=new TadpoleAssemblyPolicy(input);
			final TadpoleMulti.Config overridden=new TadpoleMulti.Config(policy.args);
			check(!policy.rescueRequested && overridden.fuseNet==null, "neural opt-out in either order");
			check(overridden.fuseCoverageRatio==0 && overridden.graphMergeCoverageRatio==0
					&& !overridden.fuseShortEvidence && !overridden.graphMergeBranchGuard, "individual safeguards disabled");
			check(!overridden.retainShortContigs && overridden.retainShortContigsSet, "explicit short-contig override");
		}
		final TadpoleAssemblyPolicy noPost=new TadpoleAssemblyPolicy(new String[]{"branchrescuenet=null"});
		check(!noPost.rescueRequested, "null postpass model disables");
		for(String[] invalid : new String[][]{{"branchrescuecutoff=NaN"}, {"branchrescuecutoff=1.1"},
				{"initialmergecoverageratio=0.5"}, {"initialmergecoverageratio=NaN"}, {"mode=correct", "branchrescue=t"}}){
			boolean rejected=false;
			try{new TadpoleAssemblyPolicy(invalid);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "reject invalid policy: "+Arrays.toString(invalid));
		}
		System.err.println("TADPOLE_ASSEMBLY_POLICY_TEST_PASS checks="+checks);
	}

	/** Throws under either assertion setting so test failures never become silent. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
}
