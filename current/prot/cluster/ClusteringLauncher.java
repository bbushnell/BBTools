package prot.cluster;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;

import dna.Data;

/** Runs a clustering, family-scoring, caller or lean benchmark entry point with a bound resource fallback.
 * This changes only the current JVM's resource lookup, never installed files.
 * Checkout resources retain precedence under Data.findPath. The selected policy
 * still validates its official resource and any explicitly supplied blacklist.
 * Intended for one invocation per JVM, not concurrent calls with different roots.
 * @author Keqing
 */
public final class ClusteringLauncher{

	private ClusteringLauncher(){}

	/** Launcher options are direct arguments; all other arguments, including config,
	 * reach the selected tool unchanged and in their original order.
	 */
	public static void main(final String[] args) throws Exception{
		require(args!=null, "ClusteringLauncher requires an argument array");
		if(args.length==0 || (args.length==1 && help(args[0]))){
			System.err.println("Usage: ClusteringLauncher target=cluster|audit|calibrate|eval|call|lean resourceroot=/absolute/resources [tool arguments]");
			System.err.println("Runs the named tool with a process-local resource fallback.");
			return;
		}
		String target=null, resourceRoot=null;
		final String[] forwarded=new String[args.length];
		int count=0;
		for(final String arg:args){
			require(arg!=null, "Null argument cannot be forwarded to the clustering tool");
			final int split=arg.indexOf('=');
			final String key=(split<0 ? arg : arg.substring(0, split)).toLowerCase(Locale.ROOT);
			final String value=(split<0 ? null : arg.substring(split+1));
			if(key.equals("target")){
				require(target==null, "Duplicate launcher target");
				require(value!=null && !value.isEmpty(), "target requires an explicit supported tool name");
				target=value.toLowerCase(Locale.ROOT);
			}else if(key.equals("resourceroot")){
				require(resourceRoot==null, "Duplicate launcher resourceroot");
				require(value!=null && !value.isEmpty(), "resourceroot requires an absolute directory");
				resourceRoot=value;
			}else{forwarded[count++]=arg;}
		}
		require("cluster".equals(target) || "audit".equals(target) || "calibrate".equals(target) ||
			"eval".equals(target) || "call".equals(target) || "lean".equals(target),
			"target must be cluster, audit, calibrate, eval, call or lean; received "+target);
		require(resourceRoot!=null, "An explicit resourceroot is required before official blacklist lookup");
		final Path root=Paths.get(resourceRoot).normalize();
		require(root.isAbsolute() && Files.isDirectory(root) && Files.isReadable(root),
			"resourceroot must be an existing readable absolute directory: "+resourceRoot);
		require(count>0, "Supply arguments for the selected "+target+" tool");
		final String[] toolArgs=Arrays.copyOf(forwarded, count);
		final String previous=Data.DORI_RESOURCES;
		Data.DORI_RESOURCES=root.toString()+File.separator;
		try{
			if(target.equals("cluster")){DdlClust.main(toolArgs);}
			else if(target.equals("audit")){ClusterOutputAudit.main(toolArgs);}
			else if(target.equals("calibrate")){prot.FamilyCalibrate.main(toolArgs);}
			else if(target.equals("eval")){prot.FamilyEval.main(toolArgs);}
			else if(target.equals("call")){prok.CallGenes.main(toolArgs);}
			else{benchmark.CallGenesLeanReplayBenchmark.main(toolArgs);}
		}finally{Data.DORI_RESOURCES=previous;}
	}

	/** Help here describes the launcher; tool-specific arguments are delegated. */
	private static boolean help(final String arg){
		return arg!=null && (arg.equalsIgnoreCase("help") || arg.equalsIgnoreCase("-h") || arg.equalsIgnoreCase("--help"));
	}

	/** Invalid launch bindings must fail before starting work or changing lookup. */
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}
}
