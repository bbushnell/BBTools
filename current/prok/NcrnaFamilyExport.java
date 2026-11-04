package prok;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Path;

import map.LongHashSet;
import prok.NcrnaTextIO.ConfigInput;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Writes a fresh installable directory for one shared-boundary family.
 * The builder supplies already validated, immutable accepted-library inputs and
 * learned products. Files are copied once; generated output is never reread here.
 * Group exports preserve manifest-listed siblings and the installed shared seeds.
 * The manifest is written last. Failure leaves an incomplete directory for diagnosis;
 * existing destinations are never overwritten. No global caller state is changed.
 * Text output follows the ByteFile template's native byte reader/writer idioms.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaFamilyExport{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies inputs corresponding to the borrowed accepted family and writes its boundary policy.
	 * @param accepted Loaded family whose library and seeds correspond to the supplied input files
	 * @param consensus Accepted FASTA, plain or gzip
	 * @param models Accepted native .hbm, .hbmc or .hbmc.gz file
	 * @param seeds Accepted forward 17-mer FASTA, plain or gzip
	 * @param boundary Learned policy and products; OFF omits the network and tables
	 * @param output Fresh family directory, or full group directory for a grouped member
	 */
	public static void write(final NcrnaFamily accepted, final String consensus,
			final String models, final String seeds, final Boundary boundary, final File output){
		final NcrnaFamilyExport export=new NcrnaFamilyExport(accepted, consensus, models, seeds, boundary, output);
		export.write();
	}

	private NcrnaFamilyExport(final NcrnaFamily accepted_, final String consensus_, final String models_,
			final String seeds_, final Boundary boundary_, final File output){
		if(accepted_==null || boundary_==null || output==null || output.exists()){
			throw new IllegalArgumentException("Family export needs a loaded family, boundary policy and fresh output directory");
		}
		family=NcrnaFamilyRegistry.sharedBoundaryMember(accepted_.name);
		consensus=readable(consensus_);
		models=readable(models_);
		seeds=readable(seeds_);
		boundary=boundary_;
		destination=output.toPath();
		upperMedianLength=PackagedBoundaryConfig.upperMedianLength(accepted_.library);
		if(boundary.startRadius>upperMedianLength || boundary.stopRadius>upperMedianLength){
			throw new IllegalArgumentException("Boundary radii cannot exceed the accepted library's upper-median length");
		}
		if(accepted_.kLong!=17 || accepted_.kmerSet==null || accepted_.kmerSet.size()==0){
			throw new IllegalArgumentException("Family export requires the accepted forward 17-mer seed set");
		}
		modelName="models"+modelSuffix(models);
		consensusName="consensus.fa"+compressionSuffix(consensus);
		seedName="seeds.fa"+compressionSuffix(seeds);
		netName=boundary.enabled ? "boundary.bbnet"+compressionSuffix(boundary.net) : null;
		if(boundary.enabled){
			readable(boundary.net);
			readable(boundary.startTable);
			readable(boundary.stopTable);
		}
		installed=resolver.installedManifest(family.id);
		sharedSeeds=family.group!=null && family.group.shared("seeds");
		for(FamilyResourceManifest.Entry entry:installed.bindings()){
			if(entry==null || !entry.member.equals(family.member)){continue;}
			if(entry.role.equals("length_bounds")){
				throw new IllegalArgumentException("Selected per-model length bounds need reconciliation with the accepted library before export");
			}
			if(!replaced(entry.role) && !entry.role.equals("cutoffs")){
				throw new IllegalArgumentException("Unsupported selected-family export role: "+entry.role);
			}
		}
		if(sharedSeeds){
			final FamilyResourceManifest.Entry binding=installed.get(".", "seeds", ".", ".");
			final String path=resolver.manifestPayload(family.id, installed, binding).path;
			final LongHashSet original=ProkObject.loadLongKmers(path, 17);
			if(original.size()!=accepted_.kmerSet.size()){
				throw new IllegalArgumentException("Accepted seeds differ from installed shared seeds; group-wide evaluation is required");
			}
			for(long key:original.toArray()){
				if(!accepted_.kmerSet.contains(key)){
					throw new IllegalArgumentException("Accepted seeds differ from installed shared seeds; group-wide evaluation is required");
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Export              ----------------*/
	/*--------------------------------------------------------------*/

	private void write(){
		try{
			Files.createDirectory(destination);
			final ByteBuilder manifest=new ByteBuilder("#ncrna_family\t");
			manifest.append(family.group==null ? 1 : 2).tab().append(family.group==null ? family.id : family.group.id).nl();
			if(family.group!=null){manifest.append("member\t");}
			manifest.append("role\tmodel\tend\tpath\tformat\tbytes\tsha80\n");
			for(FamilyResourceManifest.Entry entry:installed.bindings()){
				if(entry==null || (entry.member.equals(family.member) && replaced(entry.role))){continue;}
				final String source=resolver.manifestPayload(family.id, installed, entry).path;
				if(entry.role.equals("boundary_config")){portableSiblingConfig(source, entry.member);}
				copy(source, entry.member, entry.path);
				row(manifest, entry.member, entry.role, entry.model, entry.end, entry.path, entry.format);
			}
			copy(consensus, family.member, consensusName);
			row(manifest, family.member, "consensus", ".", ".", consensusName, "fasta");
			copy(models, family.member, modelName);
			row(manifest, family.member, "models", ".", ".", modelName, "hbm");
			if(!sharedSeeds){
				copy(seeds, family.member, seedName);
				row(manifest, family.member, "seeds", ".", ".", seedName, "fasta");
			}
			if(boundary.enabled){
				copy(boundary.net, family.member, netName);
				copy(boundary.startTable, family.member, "tables/start.tsv"+compressionSuffix(boundary.startTable));
				copy(boundary.stopTable, family.member, "tables/stop.tsv"+compressionSuffix(boundary.stopTable));
				row(manifest, family.member, "boundary", ".", ".", netName, "bbnet");
				row(manifest, family.member, "boundary_start", ".", ".", "tables/start.tsv"+compressionSuffix(boundary.startTable), "tsv");
				row(manifest, family.member, "boundary_stop", ".", ".", "tables/stop.tsv"+compressionSuffix(boundary.stopTable), "tsv");
			}
			writeText(target(family.member, "boundary.config"), configuration());
			row(manifest, family.member, "boundary_config", ".", ".", "boundary.config", "config");
			writeText(destination.resolve("manifest.tsv"), manifest);
		}catch(IOException failure){
			throw new IllegalStateException("Could not export family directory "+destination, failure);
		}
	}

	private ByteBuilder configuration(){
		final ByteBuilder text=new ByteBuilder("family=").append(family.id).nl();
		text.append("features=v1\nfeatureorder=").append(PackagedBoundaryConfig.FEATURE_ORDER).nl();
		text.append("enabled=").append(boundary.enabled ? 't' : 'f').nl();
		text.append("meanlen=").append(upperMedianLength).nl();
		text.append("startradius=").append(boundary.startRadius).nl();
		text.append("stopradius=").append(boundary.stopRadius).nl();
		text.append("marginstart=").append(decimal(boundary.startMargin)).nl();
		text.append("marginstop=").append(decimal(boundary.stopMargin)).nl();
		text.append("cutoff=");
		if(Float.isNaN(boundary.cutoff)){text.append("none");}else{text.append(decimal(boundary.cutoff));}
		text.nl();
		if(boundary.enabled){
			text.append("net=").append(netName).nl();
			text.append("tablestart=tables/start.tsv").append(compressionSuffix(boundary.startTable)).nl();
			text.append("tablestop=tables/stop.tsv").append(compressionSuffix(boundary.stopTable)).nl();
		}
		return text;
	}

	/** Input packages with external sibling dependencies cannot become portable by copying one directory. */
	private void portableSiblingConfig(final String path, final String member){
		final FamilyResourceManifest.Entry[] entries=installed.bindings();
		final String[] bindings=new String[entries.length];
		for(int i=0; i<entries.length; i++){
			final FamilyResourceManifest.Entry entry=entries[i];
			if(entry!=null && entry.member.equals(member)){
				bindings[i]=resolver.manifestPayload(family.id, installed, entry).path;
			}
		}
		requirePortableConfiguration(path, bindings);
	}

	/** Requires every relative network/table reference to resolve to a same-member manifest payload.
	 * The caller supplies resolver-produced paths from a validated manifest; null cells are ignored.
	 * Relative paths are interpreted from the configuration's parent, exactly as the consumer does.
	 * Absolute references are rejected because copying the directory cannot preserve their portability.
	 */
	public static void requirePortableConfiguration(final String path, final String[] bindings){
		if(path==null || bindings==null){throw new IllegalArgumentException("Portable configuration checking needs its path and manifest bindings");}
		final Path parent=new File(path).getAbsoluteFile().toPath().getParent();
		try(ConfigInput input=new ConfigInput(path)){
			while(input.next()){
				final String a=input.key;
				if(!a.equalsIgnoreCase("net") && !a.equalsIgnoreCase("tablestart") && !a.equalsIgnoreCase("tablestop")){continue;}
				final String relative=input.value;
				if(relative.isEmpty() || new File(relative).isAbsolute()){
					throw new IllegalArgumentException("Sibling configuration dependency must be relative: "+relative);
				}
				final Path resolved=parent.resolve(relative).normalize();
				boolean bound=false;
				for(String binding:bindings){
					bound|=binding!=null && new File(binding).getAbsoluteFile().toPath().normalize().equals(resolved);
				}
				if(!bound){throw new IllegalArgumentException("Sibling configuration dependency is not manifest-bound: "+relative);}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private void row(final ByteBuilder text, final String member, final String role, final String model, final String end, final String path, final String format){
		assert(text!=null && path!=null) : "Every copied or generated payload needs one complete manifest binding";
		if(family.group!=null){text.append(member).tab();}
		text.append(role).tab().append(model).tab().append(end).tab().append(path).tab().append(format).append("\t-\t-\n");
	}

	private Path target(final String member, final String path){
		FamilyResourceManifest.validateRelative(path);
		return (member.isEmpty() || member.equals(".") ? destination : destination.resolve(member)).resolve(path);
	}

	private void copy(final String source, final String member, final String relative) throws IOException{
		final Path target=target(member, relative);
		Files.createDirectories(target.getParent());
		Files.copy(new File(source).toPath(), target);
	}

	private static void writeText(final Path path, final ByteBuilder text){
		try(Output output=new Output(path.toString())){output.writer.print(text);}
	}

	private static String readable(final String path){
		if(path==null || !new File(path).isFile() || !new File(path).canRead()){
			throw new IllegalArgumentException("Missing readable export input: "+path);
		}
		return new File(path).getAbsolutePath();
	}

	private static String compressionSuffix(final String path){
		if(path.endsWith(".gz")){return ".gz";}
		if(path.endsWith(".bz2") || path.endsWith(".xz") || path.endsWith(".zip")){
			throw new IllegalArgumentException("Export inputs support plain or gzip files: "+path);
		}
		return "";
	}

	private static String modelSuffix(final String path){
		if(path.endsWith(".hbmc.gz")){return ".hbmc.gz";}
		if(path.endsWith(".hbmc")){return ".hbmc";}
		if(path.endsWith(".hbm")){return ".hbm";}
		throw new IllegalArgumentException("Native model suffix must be .hbm, .hbmc or .hbmc.gz: "+path);
	}

	/** Three configuration values, not a record loop: preserve every finite float without exponent notation. */
	private static String decimal(final float value){
		assert(Float.isFinite(value)) : "Only validated finite configuration values may become numeric records";
		return new BigDecimal(Float.toString(value)).toPlainString();
	}

	private static boolean replaced(final String role){
		return role.equals("consensus") || role.equals("models") || role.equals("seeds")
			|| role.equals("boundary") || role.equals("boundary_start") || role.equals("boundary_stop") || role.equals("boundary_config");
	}

	/** Immutable selection. A null network selects OFF and requires both tables to be null. */
	public static final class Boundary{
		/** Retains a validated OFF or complete shared-network selection and its endpoint geometry. */
		public Boundary(final String net_, final String startTable_, final String stopTable_, final int startRadius_,
				final int stopRadius_, final float startMargin_, final float stopMargin_, final float cutoff_){
			enabled=net_!=null;
			if(enabled ? startTable_==null || stopTable_==null : startTable_!=null || stopTable_!=null){
				throw new IllegalArgumentException("Enabled export requires one network and both tables; OFF omits all three");
			}
			if(startRadius_<0 || stopRadius_<0 || !Float.isFinite(startMargin_) || startMargin_<0
					|| !Float.isFinite(stopMargin_) || stopMargin_<0 || (!Float.isNaN(cutoff_) && !Float.isFinite(cutoff_))){
				throw new IllegalArgumentException("Boundary export needs nonnegative radii/margins and a finite cutoff or NaN for none");
			}
			net=net_;
			startTable=startTable_;
			stopTable=stopTable_;
			startRadius=startRadius_;
			stopRadius=stopRadius_;
			startMargin=startMargin_;
			stopMargin=stopMargin_;
			cutoff=cutoff_;
		}
		/** True for a complete network/table selection; false represents OFF. */
		public final boolean enabled;
		/** Source network and endpoint-table paths, all null for OFF. */
		public final String net, startTable, stopTable;
		/** Maximum correction distance in bases from each predicted endpoint. */
		public final int startRadius, stopRadius;
		/** Minimum score gain required to move the respective endpoint. */
		public final float startMargin, stopMargin;
		/** Absolute boundary score floor, or NaN when no floor is selected. */
		public final float cutoff;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final NcrnaFamilyRegistry.Entry family;
	private final String consensus, models, seeds, consensusName, modelName, seedName, netName;
	private final Boundary boundary;
	private final Path destination;
	private final int upperMedianLength;
	private final boolean sharedSeeds;
	private final FamilyResourceResolver resolver=new FamilyResourceResolver();
	private final FamilyResourceManifest installed;
}
