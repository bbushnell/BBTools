package prok;

import java.io.File;
import java.util.ArrayList;

import dna.Data;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Runs the ordinary caller through its public launcher in an owned process group.
 * Cancellation terminates that group, including a waiting shell's Java child,
 * and preserves interrupt status. It never signals an account or unrelated job.
 * Caller configuration is an input, not a receipt. Native output is forwarded.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCallerProcess{

	/** Writes and closes a fresh native caller configuration in the existing argument order.
	 * Drivers validate their scientific policy, paths and resource choices before calling.
	 * Repeated in= records preserve the full genome list without an OS argument-size limit.
	 * @param config Fresh configuration file in the driver's output directory
	 * @param arguments Validated scientific policy records
	 * @param genomes Validated genome paths in evaluation order
	 * @param gff Output GFF path, retaining the driver's chosen spelling
	 * @param family Selected registered family
	 * @param threads Positive caller worker budget
	 * @param resources Complete resource-override records owned by this driver
	 */
	public static void writeConfiguration(final File config, final String[] arguments,
			final Iterable<String> genomes, final File gff, final String family, final int threads,
			final String... resources){
		assert(config!=null && arguments!=null && genomes!=null && gff!=null && family!=null && threads>0 && resources!=null)
			: "Caller drivers must establish their input, output, policy and budget before writing a configuration";
		try(Output output=new Output(config.getPath())){
			final ByteBuilder row=new ByteBuilder();
			for(String arg:arguments){
				row.clear();
				output.writer.print(row.append(arg).nl());
			}
			for(String genome:genomes){
				row.clear();
				output.writer.print(row.append("in=").append(genome).nl());
			}
			row.clear();
			row.append("out=").append(gff.getPath()).nl();
			row.append("ncrnafamily=").append(family).nl();
			for(String resource:resources){row.append(resource).nl();}
			row.append("t=").append(threads).append("\now=f\nappend=f\n");
			output.writer.print(row);
		}
	}

	/** Uses the script beside this BBTools installation, with an explicit heap in MiB. */
	public static void runInstalled(final File config, final long heapMiB){
		final String script=Data.findPath("?callgenes.sh", false);
		if(script==null){throw new IllegalStateException("Cannot locate the public callgenes.sh launcher");}
		run(new File(script), config, heapMiB);
	}

	/** Runs one closed configuration file; the script may wait for its native caller.
	 * @param script Public caller script; explicit injection also supports finite process fixtures
	 * @param config Native flag=value configuration, with all input/output routes already assigned
	 * @param heapMiB Maximum child heap; the owner divides its budget across concurrent callers
	 */
	public static void run(final File script, final File config, final long heapMiB){
		if(script==null || !script.isFile() || config==null || !config.isFile() || heapMiB<1){
			throw new IllegalArgumentException("Caller execution needs readable script/configuration files and a positive MiB heap");
		}
		final String path=configurationPath(config);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(script.getAbsolutePath());
		// Match the existing reducer's module selection; an explicit config simd=
		// still follows this default in native argument order and retains precedence.
		command.add("simd="+(Shared.SIMD ? "t" : "f"));
		command.add("config="+path);
		NcrnaProcess.addJvmOptions(command, heapMiB);
		NcrnaProcess.run(command, null, "CallGenes for "+path);
	}

	/** Native config= currently splits its filename on commas and equals signs. */
	public static String configurationPath(final File config){
		if(config==null || !config.isFile() || !config.canRead()){
			throw new IllegalArgumentException("Cannot read native caller configuration: "+config);
		}
		final String path=config.getAbsolutePath();
		if(path.indexOf(',')>=0 || path.indexOf('=')>=0 || path.indexOf('\n')>=0 || path.indexOf('\r')>=0){
			throw new IllegalArgumentException("Native config= parsing cannot represent commas, equals signs or line breaks in its filename: "+path);
		}
		return path;
	}

	private NcrnaCallerProcess(){}
}
