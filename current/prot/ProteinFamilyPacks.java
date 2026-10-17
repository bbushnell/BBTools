package prot;

import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

/** Discovers immutable BBGAP reference packs without changing their installation.
 * Legacy models/labels stay in place; explicit packs live under packs/ID.
 * Discovery reads only small manifests. Loading verifies model and label pins;
 * packaging validation owns the bulk public-evidence payloads.
 * @author Keqing
 */
public final class ProteinFamilyPacks {

	/** Identifies the complete competition space, independent of installation path.
	 * Call after loading/validating every pack; callers keep package files immutable.
	 * Artifact manifests bind all scientific inputs, including shortlist sidecars.
	 */
	public static String referencePin(final Pack[] packs){
		require(packs!=null && packs.length>0, "Reference provenance requires at least one pack");
		final StringBuilder text=new StringBuilder("bbgap_competition_v1\n");
		int previous=0;
		for(Pack pack : packs){
			require(pack!=null && pack.loadOrder>previous, "Reference provenance requires canonical pack order");
			text.append(pack.id).append('\t').append(pack.loadOrder).append('\t')
				.append(DigestSuffix.file(pack.models.resolve("artifacts.tsv").toString())).append('\n');
			previous=pack.loadOrder;
		}
		return DigestSuffix.bytes(text.toString().getBytes(java.nio.charset.StandardCharsets.UTF_8));
	}

	/** Returns canonical load order, independent of directory enumeration order.
	 * Test packs require explicit inclusion. An installation may contain only later
	 * packs; the legacy first pack is not a prerequisite for their discovery.
	 */
	public static Pack[] discover(final Path root, final boolean includeTest) throws IOException{
		require(root!=null && Files.isDirectory(root), "BBGAP reference root is not a directory: "+root);
		final Path directory=root.toRealPath();
		final ArrayList<Pack> found=new ArrayList<Pack>();
		if(Files.exists(directory.resolve("models"))){
			require(Files.isDirectory(directory.resolve("models")), "Legacy models path is not a directory");
			found.add(new Pack(directory, "part1", 1, null));
		}
		final Path packs=directory.resolve("packs");
		if(Files.exists(packs)){
			require(Files.isDirectory(packs), "BBGAP packs path is not a directory");
			try(DirectoryStream<Path> entries=Files.newDirectoryStream(packs)){
				for(Path entry : entries){
					if(!Files.isDirectory(entry)){continue;}
					final Pack pack=read(entry);
					if(pack.test && !includeTest){continue;}
					require(Files.isDirectory(pack.models), "Missing models directory in pack "+pack.id);
					require(Files.isRegularFile(pack.labels), "Missing compressed label file in pack "+pack.id);
					found.add(pack);
				}
			}
		}
		require(!found.isEmpty(), "No enabled BBGAP reference packs in "+root);
		Collections.sort(found, Comparator.comparingInt(p->p.loadOrder));
		for(int i=1; i<found.size(); i++){
			require(found.get(i-1).loadOrder!=found.get(i).loadOrder,
				"Duplicate reference pack load_order: "+found.get(i).loadOrder);
		}
		return found.toArray(new Pack[found.size()]);
	}

	private static Pack read(final Path root){
		final Path manifest=member(root, "pack_manifest.tsv");
		require(Files.isRegularFile(manifest), "Missing pack_manifest.tsv in "+root);
		final List<String[]> rows=HmmComparisonData.rows(manifest.toString());
		final HashMap<String,String> fields=new HashMap<String,String>();
		for(String[] row : rows){
			require(row.length==2 && !row[1].isEmpty(), "Pack manifest requires nonempty key/value rows: "+root);
			require(Arrays.asList(KEYS).contains(row[0]), "Unknown pack manifest field: "+row[0]);
			require(fields.put(row[0], row[1])==null, "Duplicate pack manifest field: "+row[0]);
		}
		for(String key : KEYS){require(fields.containsKey(key), "Missing pack manifest field: "+key);}
		require("1".equals(fields.get("schema_version")) && "bbgap_reference_pack".equals(fields.get("kind")),
			"Unsupported BBGAP pack manifest schema or kind: "+root);
		final String id=fields.get("pack_id"), role=fields.get("pack_role"), state=fields.get("release_state");
		require(id.matches("[A-Za-z0-9][A-Za-z0-9._-]*") && !id.equals("part1")
			&& id.equals(root.getFileName().toString()), "Pack ID must match its directory and cannot replace part1: "+id);
		require((role.equals("test") && state.equals("development"))
			|| (role.equals("production") && (state.equals("released") || state.equals("immutable"))),
			"Invalid pack role/release state: "+role+"/"+state);
		require(!role.equals("test") || id.toLowerCase(java.util.Locale.ROOT).contains("test"),
			"Development pack ID must visibly identify its test role: "+id);
		final int order=Integer.parseInt(fields.get("load_order"));
		require(order>1, "Explicit pack load_order must follow reserved legacy order 1");
		require(Integer.parseInt(fields.get("family_count"))>=50, "Each pack must support the full top50 rule");
		require(ProteinFamilyAssigner.RULE.equals(fields.get("rule")), "Unsupported competitive rule in pack "+id);
		for(String key : PINS){DigestSuffix.requireSuffix(fields.get(key), key);}
		return new Pack(root, id, order, fields);
	}

	/** Metadata and paths for one pack. Local family indexes never change here.
	 * Legacy labels may be absent: the consumer decides whether names are required
	 * or explicitly disabled. Explicit packs always bind a compressed label file.
	 */
	public static final class Pack {
		private Pack(final Path root_, final String id_, final int order_, final HashMap<String,String> fields_){
			root=root_; id=id_; loadOrder=order_; fields=fields_; legacy=fields==null;
			test=!legacy && fields.get("pack_role").equals("test");
			models=root.resolve("models"); labels=legacy ? member(root.resolve("labels"), "family_labels.tsv")
				: root.resolve("labels/family_labels.tsv.gz");
			labelVersion=legacy ? null : fields.get("label_version");
		}

		/** Verifies stored package pins and decoded model bindings before returning
		 * a searchable library. Does not read or decompress bulk annotation evidence.
		 */
		public ProteinFamilyAssigner load() throws IOException{
			if(!legacy){
				pin(models.resolve("artifacts.tsv"), fields.get("model_artifacts_sha80"));
				pin(labels, fields.get("label_sha80"));
				pin(member(root, "manifest.tsv"), fields.get("payload_manifest_sha80"));
			}
			final ProteinFamilyAssigner library=ProteinFamilyAssigner.load(models);
			if(!legacy){
				require(library.familyCount()==Integer.parseInt(fields.get("family_count")), "Pack family count differs: "+id);
				require(library.identityPin().equals(fields.get("family_identity_sha80")), "Pack decoded identity pin differs: "+id);
				require(library.consensusPin().equals(fields.get("consensus_sha80")), "Pack decoded consensus pin differs: "+id);
			}
			return library;
		}

		/** Exact explicit-pack inventory for output-alias protection, without walking
		 * directories or reading bulk evidence. Stored manifest pins are verified by
		 * load(); this preflight only validates the declared paths and table shape.
		 * Legacy installs keep their existing model/selected-label protection rules.
		 */
		public Path[] payloadFiles(){
			if(legacy){return new Path[0];}
			final Path manifest=member(root, "manifest.tsv");
			require(Files.isRegularFile(manifest), "Missing pack payload manifest: "+manifest);
			final List<String[]> rows=HmmComparisonData.rows(manifest.toString());
			require(!rows.isEmpty() && Arrays.equals(rows.get(0), PAYLOAD_HEADER), "Unsupported pack payload manifest header: "+manifest);
			final Path directory=root.toAbsolutePath().normalize();
			final ArrayList<Path> paths=new ArrayList<Path>(); final HashSet<Path> seen=new HashSet<Path>();
			for(int i=1; i<rows.size(); i++){
				final String[] row=rows.get(i);
				require(row.length==PAYLOAD_HEADER.length && !row[1].isEmpty(), "Malformed pack payload row: "+i);
				final Path relative=Paths.get(row[1]), file=directory.resolve(relative).normalize();
				// Removing '..' lexically can change the target when a preceding component is a symlink.
				require(relative.equals(relative.normalize()), "Pack payload paths must use canonical relative spelling: "+row[1]);
				require(!relative.isAbsolute() && !file.equals(directory) && file.startsWith(directory),
					"Pack payload paths must stay relative to the pack root: "+row[1]);
				require(seen.add(file), "Duplicate pack payload path: "+row[1]); paths.add(file);
			}
			return paths.toArray(new Path[paths.size()]);
		}
		public final Path root, models, labels;
		public final String id, labelVersion;
		public final int loadOrder;
		public final boolean legacy, test;
		private final HashMap<String,String> fields;
	}

	private static void pin(final Path file, final String expected){
		require(Files.isRegularFile(file), "Missing pack artifact: "+file);
		require(expected.equals(DigestSuffix.file(file.toString())), "Pack stored-byte pin differs: "+file);
	}
	private static Path member(final Path root, final String name){
		final Path gzip=root.resolve(name+".gz");
		return Files.isRegularFile(gzip) ? gzip : root.resolve(name);
	}
	private static void require(final boolean ok, final String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	private ProteinFamilyPacks(){}
	private static final String[] KEYS={"schema_version", "kind", "pack_id", "pack_role", "release_state", "load_order", "family_count", "rule",
		"common_background_id", "profile_recipe_id", "shortlist_recipe_id", "model_artifacts_sha80", "family_identity_sha80", "consensus_sha80",
		"label_version", "label_sha80", "payload_manifest_sha80"};
	private static final String[] PINS={"model_artifacts_sha80", "family_identity_sha80", "consensus_sha80", "label_sha80", "payload_manifest_sha80"};
	private static final String[] PAYLOAD_HEADER={"package_role", "path", "bytes", "sha80", "source"};
}
