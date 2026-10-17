package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/** Metadata-only installation fixtures; model/label content pins are checked by
 * the separate real-package integration run, not by these dummy payloads.
 * @author Keqing
 */
public final class ProteinFamilyPacksTest {

	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Specify a fresh fixture directory");
		final Path base=Files.createDirectory(Paths.get(args[0]));
		final Path root=Files.createDirectory(base.resolve("installed"));
		Files.createDirectory(root.resolve("models"));
		pack(root, "part3", 3, "production", "released", true, "");
		pack(root, "part2", 2, "production", "immutable", true, "");
		//Disabled test packs may be under construction: do not open their model payloads.
		final Path test=pack(root, "testonly100", 10, "test", "development", false, "");
		final ProteinFamilyPacks.Pack[] production=ProteinFamilyPacks.discover(root, false);
		check(production.length==3 && production[0].id.equals("part1") && production[0].legacy
			&& production[1].id.equals("part2") && production[2].id.equals("part3"),
			"Discovery must retain legacy part1 and sort explicit production packs by load_order");
		reject(()->ProteinFamilyPacks.discover(root, true), "Missing models directory");
		payload(test);
		final ProteinFamilyPacks.Pack[] enabled=ProteinFamilyPacks.discover(root, true);
		check(enabled.length==4 && enabled[3].test && enabled[3].id.equals("testonly100"), "Explicit test opt-in was ignored");
		final Path later=Files.createDirectory(base.resolve("later_only"));
		pack(later, "part2", 2, "production", "released", true, "");
		check(ProteinFamilyPacks.discover(later, false).length==1, "A later pack must not require installation of part1");
		final Path duplicate=Files.createDirectory(base.resolve("duplicate_order"));
		pack(duplicate, "part2", 2, "production", "released", true, "");
		pack(duplicate, "part3", 2, "production", "released", true, "");
		reject(()->ProteinFamilyPacks.discover(duplicate, false), "Duplicate reference pack load_order");
		final Path role=Files.createDirectory(base.resolve("bad_role"));
		pack(role, "testonly100", 10, "test", "released", true, "");
		reject(()->ProteinFamilyPacks.discover(role, false), "Invalid pack role/release state");
		final Path unknown=Files.createDirectory(base.resolve("unknown_key"));
		pack(unknown, "part2", 2, "production", "released", true, "typo\tvalue\n");
		reject(()->ProteinFamilyPacks.discover(unknown, false), "Unknown pack manifest field");
		final Path repeated=Files.createDirectory(base.resolve("duplicate_key"));
		pack(repeated, "part2", 2, "production", "released", true, "load_order\t2\n");
		reject(()->ProteinFamilyPacks.discover(repeated, false), "Duplicate pack manifest field");
		final Path reserved=Files.createDirectory(base.resolve("reserved_id"));
		pack(reserved, "part1", 2, "production", "released", true, "");
		reject(()->ProteinFamilyPacks.discover(reserved, false), "cannot replace part1");
		final Path hidden=Files.createDirectory(base.resolve("hidden_test"));
		pack(hidden, "part2", 2, "test", "development", true, "");
		reject(()->ProteinFamilyPacks.discover(hidden, false), "visibly identify its test role");
		final Path empty=Files.createDirectory(base.resolve("empty"));
		reject(()->ProteinFamilyPacks.discover(empty, false), "No enabled BBGAP reference packs");
		System.err.println("FAMILY_PACK_DISCOVERY_PASS canonical_order=true legacy=true later_only=true test_optin=true rejection_cases=8");
	}

	private static Path pack(final Path root, final String id, final int order, final String role,
			final String state, final boolean data, final String extra) throws Exception{
		final Path pack=Files.createDirectories(root.resolve("packs").resolve(id));
		final String text="schema_version\t1\nkind\tbbgap_reference_pack\npack_id\t"+id
			+"\npack_role\t"+role+"\nrelease_state\t"+state+"\nload_order\t"+order
			+"\nfamily_count\t100\nrule\t"+ProteinFamilyAssigner.RULE
			+"\ncommon_background_id\tshared_v1\nprofile_recipe_id\tprofile_v1\nshortlist_recipe_id\tf4_none_v1\n"
			+"model_artifacts_sha80\t"+PIN+"\nfamily_identity_sha80\t"+PIN+"\nconsensus_sha80\t"+PIN
			+"\nlabel_version\ttest_v1\nlabel_sha80\t"+PIN+"\npayload_manifest_sha80\t"+PIN+"\n"+extra;
		Files.write(pack.resolve("pack_manifest.tsv"), text.getBytes(StandardCharsets.US_ASCII));
		if(data){payload(pack);}
		return pack;
	}
	private static void payload(final Path pack) throws Exception{
		Files.createDirectory(pack.resolve("models"));
		Files.createDirectory(pack.resolve("labels"));
		Files.write(pack.resolve("labels/family_labels.tsv.gz"), new byte[]{0});
	}
	private interface CheckedAction {void run() throws Exception;}
	private static void reject(final CheckedAction action, final String text) throws Exception{
		try{action.run();}catch(IllegalArgumentException e){check(e.getMessage().contains(text), "Unexpected rejection: "+e); return;}
		throw new AssertionError("Expected rejection: "+text);
	}
	private static void check(final boolean ok, final String reason){if(!ok){throw new AssertionError(reason);}}
	private static final String PIN="00000000000000000000";
}
