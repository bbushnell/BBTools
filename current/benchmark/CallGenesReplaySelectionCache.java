package benchmark;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import ddl.DDLFamilySearch.StructuralInputs;
import fileIO.ByteFile;
import parse.LineParser1;
import prok.CallGenes;
import prot.DigestSuffix;
import structures.ByteBuilder;

/** Small, explicitly pinned per-genome taxonomy/model manifest for lean replay.
 * FASTA/model identities are inherited from packet-level guards, not rehashed
 * per genome. Only these small manifest files are hashed here. Ordinary model
 * resolution remains authoritative and must reproduce each saved selection.
 * Construct inventories/caches serially before publishing them to replay workers;
 * loading a relocated cache extends its inventory's structural input tracker.
 * @author Sayu, Fischl
 */
public final class CallGenesReplaySelectionCache{

	/** Own the newly parsed ordered entries without relocation state. */
	private CallGenesReplaySelectionCache(Entry[] entries_){this(entries_, null, null);}
	/** Retain parsed entries and optional relocation provenance without copying.
	 * Only the loader supplies these internally owned objects. */
	private CallGenesReplaySelectionCache(Entry[] entries_, Inventory relocated_, Path original_){
		assert(entries_!=null) : "A cache must retain its complete ordered cohort";
		entries=entries_;
		relocated=relocated_;
		original=original_;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Manifest Loading       ----------------*/
	/*--------------------------------------------------------------*/

	/** Read the pinned declared model inventory, including unselected models.
	 * Checks manifest bytes before/after parsing, positive sizes and unique absolute
	 * paths. Does not enumerate a directory or read/hash model payloads. */
	public static Inventory readInventory(Path path, String pin)throws IOException{
		checkPin(path, pin);
		final HashMap<String, String> models=new HashMap<String, String>();
		final ByteFile in=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(Arrays.equals(in.nextLine(), "path\tsha80\tsize_bytes".getBytes(StandardCharsets.US_ASCII)),
				"Model inventory header");
			for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
				lp.set(line);
				require(lp.terms()==3, "Model inventory row width");
				final String name=lp.parseString(0), hash=lp.parseString(1);
				absolute(name);
				pin(hash);
				require(lp.parseLong(2)>0 && models.put(name, hash)==null, "Duplicate or empty model inventory entry");
			}
		}finally{
			require(!in.close(), "Model inventory read failed");
		}
		require(!models.isEmpty(), "Empty model inventory");
		checkPin(path, pin);
		return new Inventory(pin, models);
	}

	/** Read and validate a pinned cache into an owned ordered entry array.
	 * The caller must match every row with requireGenome before selecting models. */
	public static CallGenesReplaySelectionCache read(Path path, String pin, Inventory inventory)throws IOException{
		checkPin(path, pin);
		final ArrayList<Entry> rows=new ArrayList<Entry>();
		final HashSet<String> ids=new HashSet<String>(), fastas=new HashSet<String>();
		final ByteFile in=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(Arrays.equals(in.nextLine(), binding(inventory).getBytes(StandardCharsets.US_ASCII)),
				"Cache model-set binding changed");
			require(Arrays.equals(in.nextLine(), HEADER.getBytes(StandardCharsets.US_ASCII)), "Selection cache header");
			for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
				lp.set(line);
				require(lp.terms()==11, "Selection cache row width");
				final String[] fields=new String[11];
				for(int i=0; i<fields.length; i++){
					fields[i]=lp.parseString(i);
					token(fields[i]);
				}
				final Entry row=new Entry(fields, inventory);
				require(ids.add(row.id) && fastas.add(row.fasta), "Duplicate cached genome or FASTA");
				rows.add(row);
			}
		}finally{
			require(!in.close(), "Selection cache read failed");
		}
		require(!rows.isEmpty(), "Empty selection cache");
		checkPin(path, pin);
		return new CallGenesReplaySelectionCache(rows.toArray(new Entry[rows.size()]));
	}

	/** Explicit private relocation: identities describe unchanged original manifests.
	 * The execution wrapper must compare originals and mapped payload bytes before
	 * and after use. Structural tracking detects metadata changes; it does not certify bytes. */
	public static Inventory readRelocatedInventory(Path original, String inheritedPin, Path map, Path root)throws IOException{
		pin(inheritedPin);
		require(root!=null && root.isAbsolute() && root.normalize().equals(root)
			&& Files.isDirectory(root) && root.toRealPath().equals(root),
			"Relocation root must be a canonical absolute directory");
		final StructuralInputs inputs=new StructuralInputs();
		inputs.track(original);
		inputs.track(map);
		final HashMap<String, String> models=new HashMap<String, String>();
		final ArrayList<ModelFile> order=new ArrayList<ModelFile>();
		final ByteBuilder inverse=new ByteBuilder("path\tsha80\tsize_bytes\n");
		final ByteFile in=ByteFile.makeByteFile1(original.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(Arrays.equals(in.nextLine(), "path\tsha80\tsize_bytes".getBytes(StandardCharsets.US_ASCII)),
				"Model inventory header");
			for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
				lp.set(line);
				require(lp.terms()==3, "Model inventory row width");
				final String name=lp.parseString(0), identity=lp.parseString(1);
				final long size=lp.parseLong(2);
				token(name);
				absolute(name);
				pin(identity);
				require(size>0 && models.put(name, identity)==null, "Duplicate or empty model inventory entry");
				order.add(new ModelFile(name, identity, size));
			}
		}finally{require(!in.close(), "Model inventory read failed");}
		require(!order.isEmpty(), "Empty model inventory");
		final HashMap<String, String> forward=new HashMap<String, String>(), reverse=new HashMap<String, String>();
		final ByteFile mapping=ByteFile.makeByteFile1(map.toString(), false);
		int row=0;
		try{
			require(Arrays.equals(mapping.nextLine(), "#CallGenesReplayModelRelocation\t1".getBytes(StandardCharsets.US_ASCII)),
				"Model relocation version header");
			require(Arrays.equals(mapping.nextLine(), "original_path\trelative_path".getBytes(StandardCharsets.US_ASCII)),
				"Model relocation columns");
			for(byte[] line=mapping.nextLine(); line!=null; line=mapping.nextLine()){
				lp.set(line);
				require(lp.terms()==2 && row<order.size(), "Relocation map row width or extra row");
				final String from=lp.parseString(0), relative=lp.parseString(1);
				token(from);
				token(relative);
				final ModelFile model=order.get(row++);
				require(from.equals(model.name), "Relocation map must cover inventory exactly once in inventory order: "+from);
				final Path name=Paths.get(relative);
				require(!name.isAbsolute() && name.normalize().equals(name) && name.toString().equals(relative),
					"Relocation destination must be a normalized relative path");
				final Path target=root.resolve(name).normalize();
				require(!target.equals(root) && target.startsWith(root), "Relocation destination escapes resource root");
				Path part=root;
				for(Path segment:name){
					part=part.resolve(segment);
					require(!Files.isSymbolicLink(part), "Relocation destination contains a symlink: "+part);
				}
				require(Files.isRegularFile(target) && target.toRealPath().equals(target),
					"Missing or aliased relocation target: "+target);
				require(Files.size(target)==model.size, "Relocation target size differs: "+target);
				inputs.track(target);
				require(forward.put(from, target.toString())==null && reverse.put(target.toString(), from)==null,
					"Duplicate relocation source or destination");
				inverse.append(reverse.get(target.toString())).tab().append(model.pin).tab().append(model.size).nl();
			}
		}finally{require(!mapping.close(), "Model relocation map read failed");}
		require(row==order.size(), "Relocation map is missing inventory models");
		require(Arrays.equals(inverse.toBytes(), Files.readAllBytes(original)),
			"Relocation inventory inverse differs from original bytes or LF format");
		requireLf(map);
		inputs.verify();
		return new Inventory(inheritedPin, models, forward, reverse, inputs);
	}

	/** Read unchanged cache metadata without measuring new digests in relocation mode. */
	public static CallGenesReplaySelectionCache readRelocatedCache(Path original, String inheritedPin,
		Inventory inventory)throws IOException{
		require(inventory!=null && inventory.relocated(), "Relocated cache requires an explicit mapped inventory");
		pin(inheritedPin);
		inventory.inputs.track(original);
		final ArrayList<Entry> rows=new ArrayList<Entry>();
		final HashSet<String> ids=new HashSet<String>(), fastas=new HashSet<String>();
		final ByteFile in=ByteFile.makeByteFile1(original.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(Arrays.equals(in.nextLine(), binding(inventory).getBytes(StandardCharsets.US_ASCII)),
				"Cache model-set binding changed");
			require(Arrays.equals(in.nextLine(), HEADER.getBytes(StandardCharsets.US_ASCII)), "Selection cache header");
			for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
				lp.set(line);
				require(lp.terms()==11, "Selection cache row width");
				final String[] fields=new String[11];
				for(int i=0; i<fields.length; i++){
					fields[i]=lp.parseString(i);
					token(fields[i]);
				}
				final Entry entry=new Entry(fields, inventory);
				require(ids.add(entry.id) && fastas.add(entry.fasta), "Duplicate cached genome or FASTA");
				inventory.resolvedPath(entry.modelPath);
				rows.add(entry);
			}
		}finally{require(!in.close(), "Selection cache read failed");}
		require(!rows.isEmpty(), "Empty selection cache");
		requireLf(original);
		inventory.verify();
		return new CallGenesReplaySelectionCache(rows.toArray(new Entry[rows.size()]), inventory, original);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Bound Selections       ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of ordered genome bindings, checked against the complete input roster. */
	public int size(){return entries.length;}

	/** Check one ordinal against its exact genome, FASTA path and inherited identity.
	 * The caller must check the whole ordered cohort before any selections. */
	public void requireGenome(int ordinal, String id, String fasta, String fastaPin){
		require(ordinal>=0 && ordinal<entries.length, "Selection cache ordinal");
		final Entry row=entries[ordinal];
		require(row.id.equals(id) && row.fasta.equals(fasta) && row.fastaPin.equals(fastaPin),
			"Selection cache roster/FASTA identity changed at "+id);
	}

	/** Resolve a previously matched ordinal through the caller-owned model selector.
	 * Verification deliberately calls the fresh selector; ordinary reuse never does.
	 * Relocation requires cached taxonomy and rejects fresh verification. The caller
	 * owns the mutable CallGenes instance and must serialize its selections. */
	public CallGenes.FileModelSelection select(CallGenes owner, int ordinal, boolean verify){
		final Entry row=entries[ordinal];
		if(relocated!=null){
			require(!verify, "Relocation must preserve cached taxonomy, not perform fresh classification");
			final CallGenes.FileModelSelection selected=owner.selectModelForReplayCached(row.fasta, row.domain, row.phylum,
				row.source, row.fallback, relocated.resolvedPath(row.modelPath));
			final String[] paths=selected.modelPaths();
			require(selected.taxonomy && paths.length==1, "Relocated cache requires one taxonomy-selected model");
			validateResolved(ordinal, paths[0], selected.domain, selected.phylum, selected.classificationSource,
				selected.fallback, selected.cacheHit);
			return selected;
		}
		final CallGenes.FileModelSelection selected=verify ? owner.selectModelForReplayFresh(row.fasta) :
			owner.selectModelForReplayCached(row.fasta, row.domain, row.phylum, row.source, row.fallback, row.modelPath);
		final String[] paths=selected.modelPaths();
		require(selected.taxonomy && paths.length==1 && paths[0].equals(row.modelPath)
			&& same(selected.domain, row.domain) && same(selected.phylum, row.phylum)
			&& same(selected.classificationSource, row.source) && selected.fallback.equals(row.fallback)
			&& selected.cacheHit==row.cacheHit, "Cached and selected model metadata differ for "+row.id);
		return selected;
	}

	/** Emit the fresh or verified choices, with no FASTA/model-content hash walk. */
	public static ByteBuilder begin(Inventory inventory){
		return new ByteBuilder(binding(inventory)).nl().append(HEADER).nl();
	}

	/** Append one actual selection after validating its policy, inherited FASTA
	 * identity and membership in the pinned complete model inventory. */
	public static void append(ByteBuilder out, Inventory inventory, String id, String fasta, String fastaPin,
		CallGenes.FileModelSelection selection){
		final String[] paths=selection.modelPaths();
		require(selection.taxonomy && paths.length==1, "Cache requires per-phylum single-model selections");
		final String modelPin=inventory.models.get(paths[0]);
		require(modelPin!=null, "Resolved model is not in the bound inventory: "+paths[0]);
		final String[] f={id, fasta, fastaPin, "PER_PHYLUM", text(selection.domain), text(selection.phylum),
			text(selection.classificationSource), selection.fallback, selection.cacheHit ? "1" : "0", paths[0], modelPin};
		new Entry(f, inventory);
		for(int i=0; i<f.length; i++){
			token(f[i]);
			if(i>0){out.tab();}
			out.append(f[i]);
		}
		out.nl();
	}

	/** Shared selection boundary for native use and metadata-only rejection fixtures. */
	public void validateResolved(int ordinal, String path, String domain, String phylum, String source, String fallback,
		boolean cacheHit){
		require(relocated!=null && ordinal>=0 && ordinal<entries.length, "Relocated selection ordinal");
		final Entry row=entries[ordinal];
		require(relocated.resolvedPath(row.modelPath).equals(path) && same(domain, row.domain) && same(phylum, row.phylum)
			&& same(source, row.source) && row.fallback.equals(fallback) && cacheHit==row.cacheHit,
			"Cached and relocated selected model metadata differ for "+row.id);
	}

	/** The resolution table is not a cache bearing an inherited digest. */
	public static ByteBuilder beginResolution(){
		return new ByteBuilder("#CallGenesReplayModelResolution\t1\n")
			.append(HEADER.replace("\tmodel_path\t", "\toriginal_model_path\t")).append("\tresolved_model_path\n");
	}

	/** Build inverse evidence from actual selected fields, not the expected cache. */
	public void appendResolved(ByteBuilder inverse, ByteBuilder resolution, int ordinal, String path, String domain,
		String phylum, String source, String fallback, boolean cacheHit){
		validateResolved(ordinal, path, domain, phylum, source, fallback, cacheHit);
		final Entry row=entries[ordinal];
		final String originalPath=relocated.reverse.get(path);
		require(originalPath!=null, "Selected model has no inverse relocation");
		final String[] fields={row.id, row.fasta, row.fastaPin, "PER_PHYLUM", text(domain), text(phylum), text(source), fallback,
			cacheHit ? "1" : "0", originalPath, relocated.models.get(originalPath)};
		new Entry(fields, relocated);
		for(int i=0; i<fields.length; i++){
			token(fields[i]);
			if(i>0){
				inverse.tab();
				resolution.tab();
			}
			inverse.append(fields[i]);
			resolution.append(fields[i]);
		}
		inverse.nl();
		resolution.tab().append(path).nl();
	}

	/** Exact original bytes are recovered from all actual ordered selections before workers run. */
	public void verifyInverse(ByteBuilder inverse)throws IOException{
		require(relocated!=null, "Inverse verification requires relocation mode");
		relocated.verify();
		require(Arrays.equals(inverse.toBytes(), Files.readAllBytes(original)),
			"Actual selection inverse differs from original cache bytes");
		relocated.verify();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Bound Records         ----------------*/
	/*--------------------------------------------------------------*/

	/** Manifest identities bind payloads protected by the external packet guard.
	 * Maps are owned here and populated only during construction. Relocated cache
	 * loading extends inputs, so finish all loading before concurrent read-only use. */
	public static final class Inventory{
		/** Retain internally built ordinary model identities. */
		private Inventory(String pin_, HashMap<String, String> models_){
			this(pin_, models_, null, null, null);
		}
		/** Retain the loader-owned maps and optional structural tracker without copying. */
		private Inventory(String pin_, HashMap<String, String> models_, HashMap<String, String> forward_,
			HashMap<String, String> reverse_, StructuralInputs inputs_){
			assert(models_!=null) : "Model identities must remain separate from relocation paths";
			pin=pin_;
			models=models_;
			forward=forward_;
			reverse=reverse_;
			inputs=inputs_;
		}
		/** The configured ordinary fallback must be a member of this model set. */
		public void requireFallback(String path){
			require(relocated() ? reverse.containsKey(path) : models.containsKey(path),
				"Fallback is absent from model inventory");
		}
		/** Whether explicit forward and inverse relocation mappings are present. */
		public boolean relocated(){return forward!=null;}
		/** Number of declared models, including those unused by this cohort. */
		public int size(){return models.size();}
		/** Ordinary caller lookup must select exactly this mapped absolute path. */
		public String resolvedPath(String original){
			require(relocated() && forward.containsKey(original), "Model is absent from explicit relocation map: "+original);
			return forward.get(original);
		}
		/** Verify tracked structural metadata in relocation mode; no-op otherwise.
		 * Does not hash payloads or establish their byte identity. */
		public void verify()throws IOException{
			if(inputs!=null){inputs.verify();}
		}
		/** Inventory manifest identity, measured by ordinary loading or inherited in relocation mode. */
		public final String pin;
		/** Original absolute model paths mapped to declared payload identities. */
		private final HashMap<String, String> models;
		/** Original-to-resolved and resolved-to-original path maps; null outside relocation mode. */
		private final HashMap<String, String> forward, reverse;
		/** Relocation metadata/payload stability tracker; null for an ordinary inventory. */
		private final StructuralInputs inputs;
	}
	/** Immutable original inventory row used to require ordered, complete relocation. */
	private static final class ModelFile{
		/** Retain the validated path, declared identity and positive stored byte size. */
		ModelFile(String name_, String pin_, long size_){
			assert(size_>0) : "Mapped model size comes from the positive inventory field";
			name=name_;
			pin=pin_;
			size=size_;
		}
		/** Original absolute path and inherited payload identity. */
		final String name, pin;
		/** Declared positive stored byte size, checked against the relocated target. */
		final long size;
	}

	/** One immutable roster/classification/model binding. Nullable taxonomy values
	 * retain unavailable classifications; they are never inferred or repaired. */
	private static final class Entry{
		/** Validate the eleven parsed fields against the bound inventory.
		 * The loader separately enforces field token syntax and cohort uniqueness. */
		Entry(String[] f, Inventory inventory){
			id=f[0];
			fasta=f[1];
			fastaPin=f[2];
			domain=value(f[4]);
			phylum=value(f[5]);
			source=value(f[6]);
			fallback=f[7];
			modelPath=f[9];
			cacheHit=f[8].equals("1");
			require(id.matches("g[0-9]{4}") && f[3].equals("PER_PHYLUM"), "Cached genome/policy");
			absolute(fasta);
			absolute(modelPath);
			pin(fastaPin);
			pin(f[10]);
			require(f[8].equals("0") || f[8].equals("1"), "Cached model cache-hit field");
			require(f[10].equals(inventory.models.get(modelPath)), "Cached model identity differs from inventory");
			require(phylum==null ? fallback.equals("classification_unavailable") :
				(fallback.equals("none") || fallback.equals("phylum_model_missing")), "Cached fallback/taxonomy mismatch");
		}
		/** Immutable genome/FASTA binding and nullable taxonomy, fallback and original model path. */
		final String id, fasta, fastaPin, domain, phylum, source, fallback, modelPath;
		/** Expected model-cache state, distinct from whether this selection cache was used. */
		final boolean cacheHit;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Manifest Checks        ----------------*/
	/*--------------------------------------------------------------*/

	/** Bind even unselected models so additions/removals invalidate an old cache. */
	private static String binding(Inventory inventory){
		return "#CallGenesReplaySelectionCache\t1\tmodel_inventory_sha80\t"+inventory.pin;
	}

	/** Check only small manifest bytes before and after parsing. Payload integrity
	 * remains the responsibility of the packet-level input/model guards. */
	private static void checkPin(Path file, String expected)throws IOException{
		pin(expected);
		require(Files.isRegularFile(file) && DigestSuffix.file(file.toString()).equals(expected),
			"Manifest sha80 mismatch: "+file);
	}
	/** Validate the spelling of a 20-lowercase-hex identity; do not measure content. */
	private static void pin(String value){
		require(value!=null && value.matches("[0-9a-f]{20}"), "Measured sha80 required for selection cache");
	}
	/** A successful inverse comparison must preserve exact text representation. */
	private static void requireLf(Path file)throws IOException{
		final byte[] bytes=Files.readAllBytes(file);
		require(bytes.length>0 && bytes[bytes.length-1]=='\n', "Relocation metadata requires terminal LF");
		for(byte b:bytes){require(b!='\r', "Relocation metadata requires LF, not CRLF");}
	}
	/** Require a normalized absolute path without opening or canonicalizing its target. */
	private static void absolute(String value){
		final Path p=Paths.get(value);
		require(p.isAbsolute() && p.normalize().equals(p), "Normalized absolute cache path required");
	}
	/** TSV fields may contain ordinary spaces, including Candidatus phylum names. */
	private static void token(String value){
		require(value!=null && !value.isEmpty(), "Empty cache value");
		for(int i=0; i<value.length(); i++){
			require(value.charAt(i)>=32 && value.charAt(i)!=127, "Control character in cache value");
		}
	}
	/** Decode the reserved none token as unavailable taxonomy. */
	private static String value(String value){return value.equals("none") ? null : value;}
	/** Encode unavailable taxonomy as the reserved none token. */
	private static String text(String value){return value==null ? "none" : value;}
	/** Compare nullable taxonomy values without inferring missing classifications. */
	private static boolean same(String a, String b){return a==null ? b==null : a.equals(b);}
	/** Reject malformed or mismatched manifests independently of assertions. */
	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	/** Owned ordered cohort; immutable entries are never exposed or replaced. */
	private final Entry[] entries;
	/** Explicit relocated inventory, or null for ordinary model resolution. */
	private final Inventory relocated;
	/** Original cache file used by the inverse byte comparison; null outside relocation mode. */
	private final Path original;
	/** Exact version-one cache column order used by loading and serialization. */
	private static final String HEADER="genome\tfasta\tfasta_sha80\tpolicy\tdomain\tphylum\tsource\tfallback\tcache_hit\tmodel_path\tmodel_sha80";
}
