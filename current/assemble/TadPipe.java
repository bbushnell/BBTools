package assemble;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.Parse;
import parse.PreParser;
import shared.Shared;

/**
 * Assembles pretrimmed, filtered, quality-calibrated paired Illumina reads.
 * Every processing stage runs through its own launcher in a fresh JVM. The
 * coordinator owns filenames, resource limits, logs, and intermediate cleanup.
 *
 * @author Brian Bushnell, Fischl
 */
public class TadPipe {

	public static void main(String[] args) throws IOException, InterruptedException{
		final PreParser pp=new PreParser(args, TadPipe.class, false);
		try{
			new TadPipe(new Config(pp.args), pp.outstream).process();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Keeps pipeline state private to this invocation, never in tool statics. */
	TadPipe(final Config config_, final PrintStream outstream_){
		assert(config_!=null && outstream_!=null) : "A pipeline needs configuration and a progress stream.";
		config=config_;
		outstream=outstream_;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Execution           ----------------*/
	/*--------------------------------------------------------------*/

	/** Publishes an assembly only after every child succeeds; failure retains evidence. */
	void process() throws IOException, InterruptedException{
		config.validate();
		Files.createDirectories(config.temp);
		work=Files.createTempDirectory(config.temp, "tadpipe-").toRealPath();
		final ArrayList<Stage> stages=plan();
		validateTools(stages);
		final ArrayList<String> commands=new ArrayList<String>();
		for(Stage stage : stages){commands.add(display(stage.command));}
		Files.write(work.resolve("commands.sh"), commands, StandardCharsets.UTF_8);
		outstream.println("TadPipe work directory: "+work);
		if(config.dryrun){
			for(String command : commands){outstream.println(command);}
			outstream.println("Dry run only; no processing stages launched.");
			return;
		}
		try{
			for(Stage stage : stages){run(stage);}
			if(!Files.isRegularFile(finalAssembly) || !hasSequence(finalAssembly)){
				throw new IOException("Assembly is empty; no final output published. Inspect "+work);
			}
			if(config.overwrite){Files.move(finalAssembly, config.out, StandardCopyOption.REPLACE_EXISTING);}
			else{Files.move(finalAssembly, config.out);}
			if(config.delete){
				for(Path path : intermediates){Files.deleteIfExists(path);}
			}
			outstream.println("TadPipe complete: "+config.out+"; logs retained in "+work);
		}catch(IOException | InterruptedException | RuntimeException e){
			outstream.println("TadPipe failed; intermediates and logs retained in "+work);
			throw e;
		}
	}

	/** Compressed headers alone are nonzero bytes, but are not an assembly. */
	static boolean hasSequence(final Path path) throws IOException{
		assert(path!=null) : "Publication validation needs the exact final stage output.";
		final ByteFile input=ByteFile.makeByteFile(FileFormat.testInput(path.toString(), FileFormat.FASTA, null, true, true));
		boolean header=false, sequence=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='>'){header=true;}
				else if(header){sequence=true;}
				else{throw new IOException("Final output is not FASTA: "+path);}
			}
		}finally{
			if(input.close()){throw new IOException("Error reading completed assembly: "+path);}
		}
		return sequence;
	}

	/** Runs one child without a shell command string or inherited phase state. */
	private void run(final Stage stage) throws IOException, InterruptedException{
		assert(stage!=null && stage.command.size()>2) : "A stage must name a launcher and its explicit arguments.";
		final Path log=work.resolve(stage.name+".log");
		outstream.println("Starting "+stage.name+"; log="+log);
		final long start=System.nanoTime();
		final Process child=new ProcessBuilder(stage.command).directory(work.toFile())
				.redirectErrorStream(true).redirectOutput(log.toFile()).start();
		child.getOutputStream().close();
		final int status;
		try{
			status=child.waitFor();
		}catch(InterruptedException e){
			child.destroy();
			Thread.currentThread().interrupt();
			throw e;
		}
		if(status!=0){throw new IOException("Stage "+stage.name+" exited "+status+"; see "+log);}
		for(Path output : stage.outputs){
			if(!Files.isRegularFile(output)){
				throw new IOException("Stage "+stage.name+" did not create "+output+"; see "+log);
			}
		}
		outstream.printf(Locale.ROOT, "Finished %s in %.3f seconds.%n", stage.name, (System.nanoTime()-start)*1e-9);
	}

	/** Constructs the complete dataflow before launching any expensive processing. */
	ArrayList<Stage> plan(){
		assert(work!=null) : "The unique work directory must exist before stage paths are assigned.";
		final ArrayList<Stage> stages=new ArrayList<Stage>();
		Path paired=reads("paired");
		Stage first=stage(config.dedupe ? "01_dedupe" : "01_pair", config.dedupe ? "clumpify.sh" : "reformat.sh",
				config.dedupe ? "dedupe" : "pair");
		if(config.dedupe){first.defaults("dedupe=t", "optical=f", "ecc=f", "passes=1", "unpair=f", "repair=t");}
		first.defaults("interleaved="+(config.in2==null));
		first.overrides(config);
		first.io("in", config.in);
		if(config.in2!=null){first.io("in2", config.in2);}
		first.output("out", paired);
		stages.add(first);
		if(config.ecco){
			paired=pairedStage(stages, "02_ecco", "bbmerge.sh", "ecco", paired,
					"ecco=t", "mix=t", "kfilter=1", "k=31", "ordered=t");
		}
		if(config.clump){
			paired=pairedStage(stages, "03_clump", "clumpify.sh", "clump", paired,
					"ecc=t", "dedupe=f", "passes=4", "unpair=t", "repair=t");
		}
		if(config.ecc){
			paired=pairedStage(stages, "04_ecc", "tadpole.sh", "ecc", paired,
					"ecc=t", "k=62", "wash=t", "tu=t", "tossdepth=1", "ldf=0.15", "ordered=t");
		}
		Path merged=null, unmerged=paired;
		if(config.merge){
			final ArrayList<Path> parts=new ArrayList<Path>();
			final int[] ks=(config.rem ? new int[] {0,124,145,93} : new int[] {0});
			for(int k : ks){
				final String group="merge"+k;
				final Stage merge=stage("05_"+group, "bbmerge.sh", "merge", group);
				merge.defaults("interleaved=t", "ordered=t");
				if(k==0){merge.defaults("kfilter=1");}
				else{
					merge.defaults("rem=t", "k="+k, "extend2="+(k==124 ? 120 : k==145 ? 140 : 100));
					if(k==93){merge.defaults("strict=t");}
				}
				merge.overrides(config);
				merge.io("in", unmerged);
				if(!parts.isEmpty()){merge.args.put("extra", join(parts));}
				final Path part=reads(group);
				unmerged=reads(group+"_unmerged");
				merge.output("out", part);
				merge.output("outu", unmerged);
				parts.add(part);
				stages.add(merge);
			}
			merged=parts.get(0);
			if(parts.size()>1){
				merged=reads("merged");
				final Stage combine=stage("06_combine", "cat.sh", "combine");
				combine.args.put("in", join(parts));
				combine.output("out", merged);
				stages.add(combine);
			}
		}
		if(config.qtrim){
			unmerged=pairedStage(stages, "07_qtrim", "bbduk.sh", "qtrim", unmerged,
					"qtrim=r", "trimq=15", "maq=14", "minlen=90", "maxns=1", "ordered=t");
		}
		if(config.extend){
			final Path originalMerged=merged, originalUnmerged=unmerged;
			if(merged!=null){
				merged=extension(stages, "08_extend_merged", "extendm", originalMerged, originalUnmerged, false, 145);
			}
			unmerged=extension(stages, "09_extend_unmerged", "extendu", originalUnmerged, originalMerged, true, 124);
		}
		final Stage assembly=stage("10_assemble", "tadpole.sh", "assemble");
		assembly.defaults("assemblek=124", "k=124,300,96,64,32", "minprob=0", "minprobmain=f",
				"fusecoverageratio=1.75", "graphmergecoverageratio=1.75");
		if(config.nn){assembly.defaults("fusenet="+config.bbtools.resolve("networks/tadpole_fusion.bbnet"), "fusencutoff=0.667098");}
		assembly.defaults("interleaved="+(merged==null));
		assembly.overrides(config);
		if("null".equalsIgnoreCase(assembly.args.get("fusenet"))){assembly.args.remove("fusencutoff");}
		assembly.io("in", merged==null ? unmerged : merged);
		if(merged!=null){assembly.io("extra", unmerged);}
		finalAssembly=work.resolve(config.out.toString().endsWith(".gz") ? "assembly.fa.gz" : "assembly.fa");
		assembly.output("out", finalAssembly);
		stages.add(assembly);
		for(Stage stage : stages){stage.finish(config);}
		return stages;
	}

	/** Adds a paired-input transformation and returns its interleaved output. */
	private Path pairedStage(final ArrayList<Stage> stages, final String name, final String tool,
			final String group, final Path input, final String... defaults){
		assert(input!=null) : "Each paired transformation consumes its predecessor's output.";
		final Path output=reads(name);
		final Stage stage=stage(name, tool, group);
		stage.defaults(defaults);
		stage.defaults("interleaved=t");
		stage.overrides(config);
		stage.io("in", input);
		stage.output("out", output);
		stages.add(stage);
		return output;
	}

	/** Uses both unextended read sets as evidence, independent of extension order. */
	private Path extension(final ArrayList<Stage> stages, final String name, final String group,
			final Path input, final Path extra, final boolean paired, final int k){
		assert(input!=null) : "Extension requires its own read stream, even if it is empty.";
		final Path output=reads(name);
		final Stage stage=stage(name, "tadpole.sh", "extend", group);
		stage.defaults("mode=extend", "mce=4", "el=10", "er="+(paired ? 0 : 10), "k="+k,
				"interleaved="+paired, "ordered=t");
		stage.overrides(config);
		stage.io("in", input);
		if(extra!=null){stage.io("extra", extra);}
		stage.output("out", output);
		stages.add(stage);
		return output;
	}

	/** Applies resource defaults without forwarding them through another phase's statics. */
	private Stage stage(final String name, final String tool, final String... groups){
		assert(groups.length>0) : "Every stage has a documented override group.";
		final Stage stage=new Stage(name, tool, groups);
		stage.defaults("t="+config.threads, "zl="+config.zipLevel, "overwrite=f");
		if(tool.equals("tadpole.sh")){stage.defaults("prealloc=f", "prefilter=f", "hashkmers=explicit");}
		return stage;
	}

	/** Tracks only pipeline-created sequence files for nonrecursive success cleanup. */
	private Path reads(final String name){
		assert(work!=null && name.indexOf('/')<0) : "Intermediate names stay inside the unique work directory.";
		final Path path=work.resolve(name+(config.gz ? ".fq.gz" : ".fq"));
		intermediates.add(path);
		return path;
	}

	/** Checks launchers and the selected neural model before any stage runs. */
	private void validateTools(final ArrayList<Stage> stages) throws IOException{
		assert(!stages.isEmpty()) : "A pipeline must contain an assembly stage.";
		for(Stage stage : stages){
			if(!Files.isRegularFile(config.bbtools.resolve(stage.tool))){
				throw new IOException("Missing BBTools launcher: "+config.bbtools.resolve(stage.tool));
			}
			if(Files.exists(config.out) && Files.isSameFile(config.out, config.bbtools.resolve(stage.tool))){
				throw new IOException("Output aliases a required launcher: "+config.out);
			}
			final String net=stage.args.get("fusenet");
			if(net!=null && !net.equalsIgnoreCase("null") && !Files.isReadable(Paths.get(net))){
				throw new IOException("Unreadable fusion model: "+net);
			}
			if(net!=null && !net.equalsIgnoreCase("null") && Files.exists(config.out) && Files.isSameFile(config.out, Paths.get(net))){
				throw new IOException("Output aliases the fusion model: "+config.out);
			}
		}
	}

	/** Joins internal evidence paths for the BBTools comma-list interface. */
	private static String join(final ArrayList<Path> paths){
		assert(!paths.isEmpty()) : "A concatenation or extra-evidence list cannot be empty.";
		final StringBuilder sb=new StringBuilder();
		for(Path path : paths){if(sb.length()>0){sb.append(',');} sb.append(path);}
		return sb.toString();
	}

	/** Shell-quotes the reproduction log only; ProcessBuilder executes the original list. */
	private static String display(final ArrayList<String> command){
		assert(!command.isEmpty()) : "Cannot record an empty subprocess command.";
		final StringBuilder sb=new StringBuilder();
		for(String arg : command){
			if(sb.length()>0){sb.append(' ');}
			sb.append('\'').append(arg.replace("'", "'\\''")).append('\'');
		}
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Configuration          ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-invocation options; prefixes configure tools but cannot redirect pipeline files. */
	static final class Config {
		Config(final String[] args){
			assert(args!=null) : "PreParser supplies a nonnull argument array.";
			for(String arg : args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String key=split[0].toLowerCase(Locale.ROOT), value=(split.length>1 ? split[1] : null);
				if(key.equals("in") || key.equals("in1")){in=path(value);}
				else if(key.equals("in2")){in2=path(value);}
				else if(key.equals("out")){out=path(value);}
				else if(key.equals("bbtools")){bbtools=path(value);}
				else if(key.equals("temp") || key.equals("tempdir") || key.equals("tmpdir")){temp=path(value);}
				else if(key.equals("childheap")){heap=value;}
				else if(key.equals("t") || key.equals("threads")){threads=Integer.parseInt(value);}
				else if(key.equals("zl") || key.equals("ziplevel")){zipLevel=Integer.parseInt(value);}
				else if(key.equals("delete")){delete=Parse.parseBoolean(value);}
				else if(key.equals("gz")){gz=Parse.parseBoolean(value);}
				else if(key.equals("overwrite") || key.equals("ow")){overwrite=Parse.parseBoolean(value);}
				else if(key.equals("dryrun")){dryrun=Parse.parseBoolean(value);}
				else if(key.equals("dedupe")){dedupe=Parse.parseBoolean(value);}
				else if(key.equals("ecco")){ecco=Parse.parseBoolean(value);}
				else if(key.equals("clump") || key.equals("clumpecc")){clump=Parse.parseBoolean(value);}
				else if(key.equals("ecc") || key.equals("correct")){ecc=Parse.parseBoolean(value);}
				else if(key.equals("merge")){merge=Parse.parseBoolean(value);}
				else if(key.equals("rem")){rem=Parse.parseBoolean(value);}
				else if(key.equals("extend")){extend=Parse.parseBoolean(value);}
				else if(key.equals("qtrim")){qtrim=Parse.parseBoolean(value);}
				else if(key.equals("nn")){nn=Parse.parseBoolean(value);}
				else if(key.equals("interleaved") || key.equals("int")){
					interleaved=Parse.parseBoolean(value);
				}else if(key.indexOf('_')>0){addOverride(key, value);}
				else if(key.equals("k") || key.equals("assemblek") || key.equals("graphk")
						|| key.equals("bridgek") || key.equals("fusek")){
					//These convenience aliases configure only the final Tadpole process.
					addOverride("assemble_"+key, value);
				}
				else{throw new IllegalArgumentException("Unknown TadPipe option: "+arg+"; use a documented phase_ prefix.");}
			}
		}

		/** Fails before processing on invalid resources, missing inputs, or output aliases. */
		void validate() throws IOException{
			if(in==null || out==null || bbtools==null || temp==null){throw new IllegalArgumentException("Specify in, out, a temp directory and the BBTools installation.");}
			if(in2==null && Boolean.FALSE.equals(interleaved)){
				throw new IllegalArgumentException("Single-file input must be interleaved pairs; otherwise specify in2.");
			}
			if(threads<1 || zipLevel<1 || zipLevel>9 || heap==null || !heap.matches("[1-9][0-9]*[kKmMgG]?")){
				throw new IllegalArgumentException("Invalid t, zl or childheap; use positive threads and a Java heap such as 16g.");
			}
			if(!Files.isRegularFile(in) || !Files.isReadable(in) || (in2!=null && (!Files.isRegularFile(in2) || !Files.isReadable(in2)))){
				throw new IOException("Inputs must be readable named files; stdin is not supported.");
			}
			in=in.toRealPath();
			if(in2!=null){in2=in2.toRealPath();}
			if(in2!=null && Files.isSameFile(in, in2)){throw new IOException("in and in2 refer to the same file.");}
			if(Files.exists(out)){
				if(Files.isSameFile(in, out) || (in2!=null && Files.isSameFile(in2, out))){throw new IOException("Output aliases an input: "+out);}
				if(!overwrite){throw new IOException("Output already exists; use overwrite=t: "+out);}
			}
			if(!Files.isDirectory(out.getParent())){throw new IOException("Output parent directory does not exist: "+out.getParent());}
			for(Path path : Arrays.asList(in, in2, out, bbtools, temp)){
				if(path!=null && (path.toString().contains(",") || path.toString().contains("\n"))){
					throw new IOException("Pipeline paths cannot contain commas or newlines: "+path);
				}
			}
			final Map<String,String> assembly=overrides.get("assemble");
			if(assembly!=null && assembly.containsKey("fusenet") && !"null".equalsIgnoreCase(assembly.get("fusenet"))){
				if(!assembly.containsKey("fusencutoff")){throw new IllegalArgumentException("A custom assemble_fusenet requires assemble_fusencutoff.");}
				assembly.put("fusenet", path(assembly.get("fusenet")).toString());
			}
		}

		/** Stores overrides in CLI order; stage-specific options follow group-wide ones. */
		private void addOverride(final String key, final String value){
			final int split=key.indexOf('_');
			String group=key.substring(0, split);
			final String option=key.substring(split+1);
			if(group.equals("correct")){group="ecc";}
			if(group.equals("clumpify")){group="clump";}
			if(group.equals("extend1")){group="extend";}
			if(!Arrays.asList("dedupe", "ecco", "clump", "ecc", "merge", "merge0", "merge124", "merge145",
					"merge93", "qtrim", "extend", "extendm", "extendu", "assemble").contains(group)){
				throw new IllegalArgumentException("Unknown phase: "+group+"; inputs must already be adapter-trimmed, filtered and calibrated.");
			}
			if(option.isEmpty() || option.startsWith("in") || option.startsWith("out") || option.equals("extra")
					|| option.equals("config") || option.equals("append") || option.equals("ow") || option.equals("overwrite")
					|| option.equals("interleaved") || option.equals("int") || option.equals("gfa") || option.equals("dot")){
				throw new IllegalArgumentException("TadPipe owns stage inputs, outputs and pairing: "+key);
			}
			Map<String,String> map=overrides.get(group);
			if(map==null){map=new LinkedHashMap<String,String>(); overrides.put(group, map);}
			map.put(option, value==null ? "t" : value);
		}

		/** Resolves user paths before any child changes its working directory. */
		private static Path path(final String value){
			if(value==null || value.isEmpty() || value.equalsIgnoreCase("null")){return null;}
			return Paths.get(value).toAbsolutePath().normalize();
		}

		Path in, in2;
		Path out=path("contigs.fa"), bbtools=path("."), temp=path(System.getProperty("java.io.tmpdir"));
		String heap="14g";
		int threads=Shared.threads(), zipLevel=4;
		boolean delete=true, gz=true, overwrite=false, dryrun=false;
		Boolean interleaved=null;
		boolean dedupe=true, ecco=true, clump=true, ecc=true, merge=true, rem=true, extend=true, qtrim=true, nn=true;
		final Map<String,Map<String,String>> overrides=new LinkedHashMap<String,Map<String,String>>();
	}

	/** One fresh process, with explicit outputs checked before consumers are launched. */
	static final class Stage {
		Stage(final String name_, final String tool_, final String[] groups_){
			assert(name_!=null && tool_.endsWith(".sh")) : "Stages execute BBTools launchers, never raw tool main methods.";
			name=name_;
			tool=tool_;
			groups=groups_;
		}

		/** Later defaults or overrides replace a value without changing argument boundaries. */
		void defaults(final String... values){
			for(String value : values){
				final int split=value.indexOf('=');
				assert(split>0) : "Internal stage defaults must be key=value: "+value;
				args.put(value.substring(0, split), value.substring(split+1));
			}
		}

		/** Applies the generic group, then its specialized stage overrides. */
		void overrides(final Config config){
			assert(config!=null) : "Overrides belong to this pipeline invocation.";
			for(String group : groups){
				final Map<String,String> map=config.overrides.get(group);
				if(map!=null){args.putAll(map);}
			}
		}

		/** Adds an absolute managed input path. */
		void io(final String key, final Path value){
			assert(value!=null) : "Managed stage paths must be resolved before launching: "+key;
			args.put(key, value.toString());
		}

		/** Records an output which must exist before any consumer starts. */
		void output(final String key, final Path value){
			io(key, value);
			outputs.add(value);
		}

		/** Uses a low-starting heap; sequential children share the same maximum budget. */
		void finish(final Config config){
			assert(command.isEmpty()) : "A stage command is finalized exactly once.";
			command.add("bash");
			command.add(config.bbtools.resolve(tool).toString());
			command.add("-Xmx"+config.heap);
			command.add("-Xms32m");
			command.add("-eoom");
			for(Map.Entry<String,String> entry : args.entrySet()){command.add(entry.getKey()+"="+entry.getValue());}
		}

		final String name, tool;
		final String[] groups;
		final LinkedHashMap<String,String> args=new LinkedHashMap<String,String>();
		final ArrayList<Path> outputs=new ArrayList<Path>();
		final ArrayList<String> command=new ArrayList<String>();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Config config;
	private final PrintStream outstream;
	private final ArrayList<Path> intermediates=new ArrayList<Path>();
	Path work;
	private Path finalAssembly;
}
