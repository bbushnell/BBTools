package test.prot.cluster;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import dna.Data;
import prot.cluster.ClusteringLauncher;

/**
 * Checks six explicit launcher targets and resource-fallback restoration on failure.
 * The unknown tool argument and its stack trace establish dispatch; this does not
 * test successful resource loading, which belongs to separate fresh-JVM fixtures.
 * Creates one fresh empty directory and leaves it in place. Use an isolated JVM:
 * target parsing can touch process-global settings beyond the fallback checked here.
 * Documentation and formatting by Fischl; original authorship is unresolved.
 */
public final class ClusteringLauncherTest{

	/** Checks dispatch, duplicate/invalid bindings and the rejected historical replay target.
	 * @param args One fresh fixture directory whose parent exists; normalized to an absolute path */
	public static void main(String[] args)throws Exception{
		require(args.length==1, "Supply a fresh fixture directory");
		final Path root=Paths.get(args[0]).toAbsolutePath().normalize();
		require(!Files.exists(root), "Launcher fixture directory must be fresh");
		Files.createDirectory(root);
		final String previous=Data.DORI_RESOURCES;
		final String[] targets={"cluster", "audit", "calibrate", "eval", "call", "lean"};
		final String[] classes={"prot.cluster.DdlClust", "prot.cluster.ClusterOutputAudit", "prot.FamilyCalibrate",
			"prot.FamilyEval", "prok.CallGenes", "benchmark.CallGenesLeanReplayBenchmark"};
		for(int i=0; i<targets.length; i++){
			Throwable failure=null;
			try{
				ClusteringLauncher.main(new String[]{"target="+targets[i], "resourceroot="+root,
					"launcher_fixture_unknown=t"});
			}catch(IllegalArgumentException | AssertionError expected){failure=expected;}
			require(failure!=null && failure.getMessage()!=null &&
				failure.getMessage().contains("launcher_fixture_unknown"),
				"Target must receive the unchanged tool argument: "+targets[i]);
			boolean dispatched=false;
			for(StackTraceElement frame:failure.getStackTrace()){dispatched|=frame.getClassName().equals(classes[i]);}
			require(dispatched, "Failure must originate in the selected concrete tool: "+classes[i]);
			require(previous.equals(Data.DORI_RESOURCES), "Resource fallback must be restored after target failure");
		}
		reject(new String[]{"target=java.lang.Runtime", "resourceroot="+root, "x=t"}, "target must be");
		reject(new String[]{"target=replay", "resourceroot="+root, "x=t"}, "target must be");
		reject(new String[]{"target=replay", "resourceroot=relative", "x=t"}, "target must be");
		reject(new String[]{"target=cluster", "target=audit", "resourceroot="+root, "x=t"},
			"Duplicate launcher target");
		reject(new String[]{"target=audit", "resourceroot="+root, "resourceroot="+root, "x=t"},
			"Duplicate launcher resourceroot");
		reject(new String[]{"target=eval", "resourceroot=relative", "x=t"}, "absolute directory");
		reject(new String[]{"target=calibrate", "resourceroot="+root.resolve("absent"), "x=t"},
			"absolute directory");
		reject(new String[]{"target=lean", "resourceroot="+root}, "Supply arguments");
		require(previous.equals(Data.DORI_RESOURCES), "Invalid launcher bindings must preserve fallback");
		System.out.println("CLUSTERING_LAUNCHER_TEST_PASS targets=6 replay_rejected=true fallback_restored=true");
	}

	/** Reject invalid launcher bindings before invoking a tool or changing global lookup. */
	private static void reject(String[] args, String message)throws Exception{
		try{ClusteringLauncher.main(args);}catch(IllegalArgumentException expected){
			require(expected.getMessage().contains(message), "Wrong launcher rejection: "+expected);
			return;
		}
		throw new AssertionError("Invalid launcher binding was accepted");
	}

	/** Throws AssertionError on a failed expectation, independent of -ea. */
	private static void require(boolean valid, String message){
		if(!valid){throw new AssertionError(message);}
	}
}
