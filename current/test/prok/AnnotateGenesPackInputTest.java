package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import prot.ProteinFamilyAssigner;
import prok.AnnotateGenes;
import shared.Resources;

/** CLI/discovery/alias fixtures without loading HBM data or claiming assignment
 * parity. Metadata uses dummy pins; real package loading is an integration test.
 * @author Keqing
 */
public final class AnnotateGenesPackInputTest {

	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Specify a fresh fixture directory");
		final Path base=Files.createDirectory(Paths.get(args[0]));
		final Path fasta=write(base.resolve("genome.fna"), ">ctg\nATGAAATAA\n");
		final Path gff=write(base.resolve("genes.gff"), "ctg\tx\tCDS\t1\t9\t.\t+\t0\tID=a;partial=00\n");
		final Path root=Files.createDirectory(base.resolve("reference"));
		final Path second=pack(root, "part2", 2, false, true);
		final Path out=base.resolve("output.gff");
		construct(fasta, gff, out, "reference="+root);
		construct(fasta, gff, out, "reference="+root, "labels=null");
		construct(fasta, gff, out, "models="+second.resolve("models"));
		reject(()->construct(fasta, gff, out, "reference="+root, "models="+second.resolve("models")), "either models=DIR");
		reject(()->construct(fasta, gff, out, "models="+second.resolve("models"), "includetest=t"), "includetest applies");
		final Path third=pack(root, "part3", 3, false, true);
		construct(fasta, gff, out, "reference="+root);
		reject(()->construct(fasta, gff, third.resolve("evidence.tsv"), "reference="+root), "immutable model package");
		final Path model=write(third.resolve("models/model.hbmc"), "model fixture\n");
		final Path alias=base.resolve("model_alias.gff"); Files.createLink(alias, model);
		reject(()->construct(fasta, gff, alias, "reference="+root), "aliases a model package member");
		final Path manifestAlias=base.resolve("manifest_alias.gff");
		Files.createLink(manifestAlias, third.resolve("pack_manifest.tsv"));
		reject(()->construct(fasta, gff, manifestAlias, "reference="+root), "aliases a reference pack member");
		reject(()->construct(fasta, gff, out, "reference="+root, "labels="+second.resolve("labels/family_labels.tsv.gz")),
			"labels=FILE override requires a single pack");
		construct(fasta, gff, out, "reference="+root, "labels=null");
		final Path test=pack(root, "testonly100", 10, true, false);
		construct(fasta, gff, out, "reference="+root);
		missing(()->construct(fasta, gff, out, "reference="+root, "includetest=t"));
		payload(test);
		construct(fasta, gff, out, "reference="+root, "includetest=t");
		final Path legacy=Files.createDirectory(base.resolve("legacy")); Files.createDirectory(legacy.resolve("models"));
		missing(()->construct(fasta, gff, out, "reference="+legacy, "labels=null"));
		missing(()->construct(fasta, gff, out, "models="+legacy.resolve("models"), "labels=null"));
		write(legacy.resolve("models/artifacts.tsv"), "metadata fixture; not a loadable HBM manifest\n");
		missing(()->construct(fasta, gff, out, "reference="+legacy));
		construct(fasta, gff, out, "reference="+legacy, "labels=null");
		construct(fasta, gff, out, "models="+legacy.resolve("models"));
		final Path emptyPack=Files.createDirectories(base.resolve("empty_modular/packs/part2"));
		missing(()->construct(fasta, gff, out, "reference="+emptyPack.getParent().getParent()));
		final Path badArtifacts=Files.createDirectories(base.resolve("bad_artifacts/models/artifacts.tsv"));
		reject(()->construct(fasta, gff, out, "models="+badArtifacts.getParent(), "labels=null"),
			"Legacy model artifacts path is not a regular file");
		final Path danglingModels=base.resolve("dangling_models");
		Files.createSymbolicLink(danglingModels, Paths.get("absent_model_target"));
		reject(()->construct(fasta, gff, out, "models="+danglingModels, "labels=null"), "Model package directory not found");
		final Path labels=Files.createDirectory(legacy.resolve("labels"));
		Files.createSymbolicLink(labels.resolve("family_labels.tsv.gz"), Paths.get("absent_labels_target"));
		reject(()->construct(fasta, gff, out, "reference="+legacy), "BBGAP label path is not a regular file");
		final Path danglingPacks=Files.createDirectories(base.resolve("dangling_pack/packs"));
		Files.createSymbolicLink(danglingPacks.resolve("part2"), Paths.get("absent_pack_target"));
		reject(()->construct(fasta, gff, out, "reference="+danglingPacks.getParent()),
			"BBGAP pack entry is a dangling or non-directory symbolic link");
		final Path partialRoot=Files.createDirectory(base.resolve("missing_payload_manifest"));
		final Path partial=pack(partialRoot, "part2", 2, false, true);
		Files.delete(partial.resolve("manifest.tsv"));
		missing(()->construct(fasta, gff, out, "reference="+partialRoot));
		Files.createSymbolicLink(partial.resolve("manifest.tsv.gz"), Paths.get("absent_manifest_target"));
		reject(()->construct(fasta, gff, out, "reference="+partialRoot), "Pack payload manifest is not a regular file");
		Files.delete(partial.resolve("manifest.tsv.gz"));
		Files.createDirectory(partial.resolve("manifest.tsv"));
		reject(()->construct(fasta, gff, out, "reference="+partialRoot), "Pack payload manifest is not a regular file");
		reject(()->construct(fasta, gff, out, "reference=null"), "reference directory not found");
		check(!Files.exists(out) && new String(Files.readAllBytes(model), StandardCharsets.US_ASCII).equals("model fixture\n"),
			"Constructor validation must not write outputs or aliased inputs");
		System.err.println("ANNOTATE_GENES_PACK_INPUT_PASS later_only=true test_optin=true legacy_optional_labels=true aliases_preserved=true partial_pack_absence=true malformed_paths_strict=true rejection_cases=19");
	}

	private static void construct(final Path fasta, final Path gff, final Path out, final String... options) throws Exception{
		final String[] args=new String[options.length+4];
		args[0]=fasta.toString(); args[1]=gff.toString(); args[2]="out="+out; args[3]="t=1";
		System.arraycopy(options, 0, args, 4, options.length); new AnnotateGenes(args);
	}
	private static Path pack(final Path root, final String id, final int order, final boolean test, final boolean data) throws Exception{
		final Path pack=Files.createDirectories(root.resolve("packs").resolve(id));
		write(pack.resolve("pack_manifest.tsv"), "schema_version\t1\nkind\tbbgap_reference_pack\npack_id\t"+id
			+"\npack_role\t"+(test ? "test" : "production")+"\nrelease_state\t"+(test ? "development" : "released")
			+"\nload_order\t"+order+"\nfamily_count\t100\nrule\t"+ProteinFamilyAssigner.RULE
			+"\ncommon_background_id\tshared_v1\nprofile_recipe_id\tprofile_v1\nshortlist_recipe_id\tf4_none_v1\n"
			+"model_artifacts_sha80\t"+PIN+"\nfamily_identity_sha80\t"+PIN+"\nconsensus_sha80\t"+PIN
			+"\nlabel_version\ttest_v1\nlabel_sha80\t"+PIN+"\npayload_manifest_sha80\t"+PIN+"\n");
		if(data){payload(pack);} return pack;
	}
	private static void payload(final Path pack) throws Exception{
		Files.createDirectory(pack.resolve("models")); Files.createDirectory(pack.resolve("labels"));
		write(pack.resolve("models/artifacts.tsv"), "metadata fixture; not a loadable HBM manifest\n");
		try(java.io.OutputStream out=new java.util.zip.GZIPOutputStream(Files.newOutputStream(pack.resolve("labels/family_labels.tsv.gz")))){
			out.write("metadata fixture\n".getBytes(StandardCharsets.US_ASCII));
		}
		Files.createDirectory(pack.resolve("evidence"));
		write(pack.resolve("evidence/decisions.tsv"), "retained evidence fixture\n");
		write(pack.resolve("manifest.tsv"), "package_role\tpath\tbytes\tsha80\tsource\n"
			+"selected_labels\tlabels/family_labels.tsv.gz\t1\t"+PIN+"\tfixture\n"
			+"decisions_raw_provisional\tevidence/decisions.tsv\t1\t"+PIN+"\tfixture\n");
	}
	private static Path write(final Path path, final String text) throws Exception{return Files.write(path, text.getBytes(StandardCharsets.US_ASCII));}
	private interface CheckedAction {void run() throws Exception;}
	/** Missing installed components must use the typed acquisition diagnostic, not a generic input error. */
	private static void missing(final CheckedAction action) throws Exception{
		try{action.run();}
		catch(Resources.MissingResourceException e){
			check(e.getMessage().contains("Required AnnotateGenes HBM reference pack not found:") &&
				e.getMessage().contains("The AnnotateGenes HBM reference pack is not distributed yet."),
				"Missing installation must identify the unavailable HBM pack: "+e);
			return;
		}
		throw new AssertionError("Expected missing HBM installation guidance");
	}
	private static void reject(final CheckedAction action, final String text) throws Exception{
		try{action.run();}catch(IllegalArgumentException e){check(e.getMessage().contains(text), "Unexpected rejection: "+e); return;}
		throw new AssertionError("Expected rejection: "+text);
	}
	private static void check(final boolean ok, final String reason){if(!ok){throw new AssertionError(reason);}}
	private static final String PIN="00000000000000000000";
}
