package test.prok;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import prok.CallGenes;
import prok.FamilyResourceResolver;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaFamilyRegistry;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PackagedBoundaryConfig;
import prok.TrnaConsensusBuilder;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Finite packaged-default precedence and native loader checks using installed SAM-I resources.
 * Temporary links point to read-only resources; no installed file is modified.
 * @author Brian Bushnell, Raiden
 */
public final class PackagedBoundaryConfigTest{

	public static void main(String[] args) throws Exception{
		assert(args.length==1 && args[0].startsWith("installation=")) : "Fixture requires its launcher's installation path";
		final String installation=args[0].substring("installation=".length());
		Shared.setThreads(1);
		ReadWrite.setZipThreads(1);
		final Path root=Files.createTempDirectory("packaged-boundary-");
		final String[] roles={"consensus", "models", "seeds"};
		final Path[] links=new Path[roles.length];
		final Path manifest=root.resolve("manifest.tsv");
		final Path config=root.resolve("boundary.config");
		try(AutoCloseable cleanup=()->{
			Files.deleteIfExists(config);
			Files.deleteIfExists(manifest);
			for(Path link:links){if(link!=null){Files.deleteIfExists(link);}}
			Files.delete(root);
		}){
			final FamilyResourceResolver original=new FamilyResourceResolver();
			final ByteBuilder bb=new ByteBuilder("#ncrna_family\t1\tsam1\nrole\tmodel\tend\tpath\tformat\tbytes\tsha80\n");
			for(int i=0; i<roles.length; i++){
				final File source=new File(original.resolve("sam1", roles[i], null).path).getAbsoluteFile();
				links[i]=root.resolve(roles[i]+"_"+source.getName());
				Files.createSymbolicLink(links[i], source.toPath());
				bb.append(roles[i]).append("\t.\t.\t").append(links[i].getFileName().toString()).tab();
				bb.append(roles[i].equals("models") ? "hbm" : "fasta").append("\t-\t-\n");
			}
			final String core=bb.toString();
			bb.append("boundary_config\t.\t.\tboundary.config\tconfig\t-\t-\n");
			write(manifest, bb);
			final byte[][] library=TrnaConsensusBuilder.loadLibrary(links[0].toString());
			final int mean=PackagedBoundaryConfig.upperMedianLength(library);
			checkRecords(config, mean, library);
			write(config, recipe(false, mean, null, null, null));
			final NcrnaFamilyLoader.Options off=options(root);
			off.boundary=true;
			final NcrnaFamily disabled=load(off);
			assert(disabled.boundary5NetTemplate==null && disabled.boundary3NetTemplate==null)
				: "A packaged OFF selection must replace the historical default when no explicit flag exists";
			assert(disabled.boundaryStartOffsets.length==5 && disabled.boundaryStopOffsets.length==7)
				: "Packaged radii2/3 must reach the native family";
			assert(disabled.boundaryMarginStart==.12f && disabled.boundaryMarginStop==.15f)
				: "Unspecified margins must use packaged defaults";

			write(config, recipe(true, mean, "missing.net", "missing.start", "missing.stop"));
			final NcrnaFamilyLoader.Options explicitOff=options(root);
			explicitOff.boundaryOverride=false;
			explicitOff.boundaryStartOffsets=new int[]{0};
			explicitOff.boundaryMarginStart=0f;
			final NcrnaFamily overridden=load(explicitOff);
			assert(overridden.boundary5NetTemplate==null && overridden.boundaryStartOffsets.length==1
				&& overridden.boundaryMarginStart==0f && overridden.boundaryMarginStop==.15f)
				: "Explicit OFF/radius/zero-margin must win without opening overridden missing payloads";
			write(config, recipe(true, mean+1, "missing.net", "missing.start", "missing.stop"));
			final NcrnaFamilyLoader.Options differentLibraryOff=options(root);
			differentLibraryOff.boundaryOverride=false;
			assert(load(differentLibraryOff).boundary5NetTemplate==null)
				: "Explicit OFF must ignore unused packaged normalization from a different model subset";
			boolean wrongMedian=false;
			try{load(options(root));}catch(IllegalArgumentException e){wrongMedian=e.getMessage().contains("normalization");}
			assert(wrongMedian) : "Enabled refinement must reject normalization from a different accepted library";
			write(config, recipe(true, mean, "missing.net", "missing.start", "missing.stop"));
			checkCli(installation, root);

			final String net=original.resolve("sam1", "boundary", null).path;
			final String start=original.resolve("sam1", "boundary_start", null).path;
			final String stop=original.resolve("sam1", "boundary_stop", null).path;
			write(config, recipe(true, mean, net, start, stop));
			final NcrnaFamily enabled=load(options(root));
			assert(enabled.boundary5NetTemplate!=null && enabled.boundary3NetTemplate!=null)
				: "A packaged ON selection must load its declared native network and tables";
			final NcrnaFamilyLoader.Options badPath=options(root);
			badPath.resources.boundary=root.resolve("explicit-missing.net").toString();
			boolean failed=false;
			try{load(badPath);}catch(IllegalArgumentException e){failed=e.getMessage().contains("explicit-missing.net");}
			assert(failed) : "An explicit path must not silently fall back to the package's working network";

			write(manifest, new ByteBuilder(core));
			final NcrnaFamilyLoader.Options absent=options(root);
			absent.boundary=false;
			absent.boundaryOverride=true;
			final NcrnaFamily unchanged=load(absent);
			assert(unchanged.boundary5NetTemplate==null && unchanged.boundaryMarginStart==0f
				&& unchanged.boundaryMarginStop==0f)
				: "Without a packaged role, historical enable selection and zero margins must be unchanged";
			checkFlagPrecedence();
			System.out.println("Packaged boundary checks passed: ON/OFF, explicit choices, missing payloads, absent package and flag precedence");
		}
	}

	private static void checkRecords(Path path, int mean, byte[][] library){
		assert(mean>0 && library.length>0) : "Record checks use the actual fixture's valid normalization";
		final String base=recipe(false, mean, null, null, null).toString();
		for(final String option:new String[]{"meanlen", "startradius", "stopradius", "marginstart", "marginstop", "cutoff"}){
			for(final String value:new String[]{"", "not_a_number"}){
				write(path, new ByteBuilder(base.replaceFirst(option+"=[^\\n]*", option+"="+value)));
				boolean rejected=false;
				try{PackagedBoundaryConfig.read(path.toString());}catch(IllegalArgumentException failure){
					rejected=failure.getMessage().contains(option);
				}
				assert(rejected) : "Missing/malformed packaged number must identify its option: "+option+"="+value;
			}
		}
		for(String family:new String[]{"rnasep", "tmrna", "fmn", "sam1", "sixs_rf00013", "sixs_rf01685",
			"srp_small", "srp_large", "cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"}){
			assert(NcrnaFamilyRegistry.sharedBoundaryMember(family).id.equals(family))
				: "Explicit shared-boundary members must remain supported: "+family;
			assert(NcrnaFamilyRegistry.modelOverrideMember(family).id.equals(family))
				: "Every current shared-boundary member also permits a paired replacement library: "+family;
			write(path, new ByteBuilder(base.replace("family=sam1", "family="+family)));
			PackagedBoundaryConfig.read(path.toString());
		}
		for(String family:new String[]{"tRNA", "r58", "lsu", "prok5S", "16S", "23S", "euk5S", "euk18S", "cobalamin"}){
			boolean rejected=false;
			try{
				NcrnaFamilyRegistry.sharedBoundaryMember(family);
			}catch(IllegalArgumentException expected){
				rejected=true;
			}
			assert(rejected) : "Specialized adapters and aggregate selectors cannot become shared boundary models: "+family;
		}
		for(String family:new String[]{"r58", "lsu"}){
			assert(NcrnaFamilyRegistry.modelOverrideMember(family).id.equals(family))
				: "Paired model replacement must retain 5.8S/LSU support independently of boundary training: "+family;
		}
		for(String family:new String[]{"tRNA", "prok5S", "16S", "23S", "euk5S", "euk18S", "cobalamin"}){
			boolean rejected=false;
			try{NcrnaFamilyRegistry.modelOverrideMember(family);}catch(IllegalArgumentException expected){rejected=true;}
			assert(rejected) : "Independent verifiers, unpaired libraries and aggregate selectors cannot accept one generic model pair: "+family;
		}
		final String[] bad={base+"unknown=1\n", base+"enabled=f\n",
			base.replace("family=sam1", "family=cobalamin"),
			base.replace("enabled=f", "enabled=maybe"),
			base.replace("features=v1", "features=v2"),
			base.replace("meanlen="+mean, "meanlen=4294967297")};
		for(String text:bad){
			write(path, new ByteBuilder(text));
			boolean rejected=false;
			try{
				PackagedBoundaryConfig.read(path.toString());
			}catch(IllegalArgumentException e){
				rejected=true;
			}
			assert(rejected) : "Unknown/duplicate keys, feature mismatch and overflowing normalization must fail loudly";
		}
		write(path, recipe(true, mean, "missing=network", "missing.start", "missing.stop"));
		final NcrnaFamilyLoader.Options off=new NcrnaFamilyLoader.Options();
		off.boundaryOverride=false;
		PackagedBoundaryConfig.read(path.toString()).applyDefaults("sam1", off, library);
		assert(off.resources.boundary.endsWith("missing=network"))
			: "Only the first equals sign separates the key; the resource name must remain intact";
	}

	/** Missing package payloads distinguish actual ON/OFF parsing before any genomic call. */
	private static void checkCli(String installation, Path root) throws Exception{
		assert(new File(installation, "callgenes.sh").isFile()) : "CLI fixture must use the tested installation's ordinary caller";
		final Path input=root.resolve("input.fa");
		final ByteBuilder sequence=new ByteBuilder(">fixture\n");
		for(int i=0; i<500; i++){sequence.append('A');}
		write(input, sequence.nl());
		try(AutoCloseable inputCleanup=()->Files.deleteIfExists(input)){
			for(int scenario=0; scenario<3; scenario++){
				final Path output=root.resolve("cli"+scenario+".gff");
				final Path log=root.resolve("cli"+scenario+".log");
				try(AutoCloseable trialCleanup=()->{
					Files.deleteIfExists(output);
					Files.deleteIfExists(log);
				}){
					final ArrayList<String> command=new ArrayList<String>();
					command.add("bash");
					command.add(new File(installation, "callgenes.sh").getPath());
					for(String arg:new String[]{"in="+input, "out="+output, "sam1root="+root, "sam1=t",
						"ncrna=t", "fmn=f", "cobalamin=f", "taxonomy=f", "callcds=f", "euk=f", "ow=f", "t=1", "-Xmx1g"}){
						command.add(arg);
					}
					command.add("ncrnaboundarynet="+(scenario==2 ? "t" : "f"));
					if(scenario>0){command.add("sam1boundarynn="+(scenario==1 ? "t" : "f"));}
					final int exit=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log.toFile()).start().waitFor();
					final boolean passed=scenario==1 ? exit!=0 && readLog(log, "missing.net", false) : exit==0;
					if(!passed){readLog(log, null, true);}
					assert(passed) : "CLI per-family/global OFF precedence failed; scenario="+scenario+", exit="+exit;
				}
			}
		}
	}

	private static boolean readLog(Path path, String text, boolean echo){
		assert(path!=null) : "Expected failure must be identified from its own diagnostic log";
		boolean found=false;
		try(Input owner=new Input(path.toString())){
			final ByteFile input=owner.reader;
			for(byte[] row=input.nextLine(); row!=null; row=input.nextLine()){
				if(text!=null){found|=Tools.indexOf(row, text, 0)>=0;}
				if(echo){
					System.err.write(row, 0, row.length);
					System.err.println();
				}
			}
		}
		return found;
	}

	private static NcrnaFamilyLoader.Options options(Path root){
		assert(root!=null) : "Each loader case must bind the isolated manifest";
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.root=root.toString();
		return options;
	}

	private static NcrnaFamily load(NcrnaFamilyLoader.Options options){
		assert(options!=null) : "Fixture must state its explicit overrides";
		return NcrnaFamilyLoader.load("sam1", options, new FamilyResourceResolver());
	}

	private static ByteBuilder recipe(boolean enabled, int mean, String net, String start, String stop){
		assert(mean>0) : "Boundary length normalization uses the installed library's upper median";
		final ByteBuilder bb=new ByteBuilder("family=sam1\nfeatures=v1\nfeatureorder=");
		bb.append(PackagedBoundaryConfig.FEATURE_ORDER).nl();
		bb.append("enabled=").append(enabled ? 't' : 'f').nl();
		bb.append("meanlen=").append(mean).nl();
		bb.append("startradius=2\nstopradius=3\nmarginstart=0.12\nmarginstop=0.15\ncutoff=none\n");
		if(enabled){
			bb.append("net=").append(net).nl();
			bb.append("tablestart=").append(start).nl();
			bb.append("tablestop=").append(stop).nl();
		}
		return bb;
	}

	private static void write(Path path, ByteBuilder text){
		assert(path!=null && text!=null) : "Fixture configuration requires a complete explicit record";
		try(Output owner=new Output(new ByteStreamWriter(path.toString(), true, false, false))){
			owner.writer.print(text);
		}
	}

	private static void checkFlagPrecedence(){
		final Boolean oldGlobal=CallGenes.NCRNA_BOUNDARY_NN_OVERRIDE;
		final Boolean oldFamily=CallGenes.SAM1_BOUNDARY_NN_OVERRIDE;
		try{
			CallGenes.NCRNA_BOUNDARY_NN_OVERRIDE=null;
			CallGenes.SAM1_BOUNDARY_NN_OVERRIDE=null;
			assert(CallGenes.ncrnaBoundaryOverride("sam1")==null) : "Absent flags must leave the package in control";
			CallGenes.NCRNA_BOUNDARY_NN_OVERRIDE=false;
			assert(Boolean.FALSE.equals(CallGenes.ncrnaBoundaryOverride("sam1"))) : "Explicit global OFF beats package ON";
			CallGenes.SAM1_BOUNDARY_NN_OVERRIDE=true;
			assert(Boolean.TRUE.equals(CallGenes.ncrnaBoundaryOverride("sam1"))) : "Per-family ON beats global OFF";
			CallGenes.SAM1_BOUNDARY_NN_OVERRIDE=false;
			CallGenes.NCRNA_BOUNDARY_NN_OVERRIDE=true;
			assert(Boolean.FALSE.equals(CallGenes.ncrnaBoundaryOverride("sam1"))) : "Per-family OFF beats global ON";
		}finally{
			CallGenes.NCRNA_BOUNDARY_NN_OVERRIDE=oldGlobal;
			CallGenes.SAM1_BOUNDARY_NN_OVERRIDE=oldFamily;
		}
	}
}
