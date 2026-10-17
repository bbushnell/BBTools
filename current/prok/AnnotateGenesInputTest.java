package prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;

import structures.ByteBuilder;

/** Independent file/metadata fixtures; deliberately does not need an HBM library.
 * Full-model assignment and threaded output parity require the integration pilot.
 * @author Keqing
 */
public final class AnnotateGenesInputTest {

	public static void main(String[] args) throws Exception{
		check(args.length==1, "Specify one fresh fixture directory");
		final Path dir=Files.createDirectory(Paths.get(args[0]));
		final Path fasta=write(dir, "input.fna", ">ctg|one description\nATGNNNTAA\n>table4\nATGTGATAA\n");
		final HashMap<String,byte[]> sequences=AnnotateGenes.readScaffolds(fasta.toString());
		check(sequences.size()==2 && sequences.containsKey("ctg|one"), "FASTA first-token identifiers must match GFF identifiers");
		final String feature="ctg%7Cone\tcaller\tCDS\t1\t9\t7.25\t+\t0\tID=a;partial=00;product=original;Note=bbgap_status=inside_note";
		final String rna="ctg%7Cone\tcaller\ttRNA\t1\t8\t.\t-\t.\tID=rna;anticodon=TCC;Note=keep%3Bthis";
		final String[] original={"##gff-version 3", "# arbitrary comment", feature, rna, "##FASTA", ">ctg|one", "ATGNNNTAA"};
		final Path gff=write(dir, "mixed.gff", String.join("\r\n", original)+"\r\n");
		final AnnotationGff doc=AnnotationGff.load(gff.toString(), sequences, 11);
		check(doc.lines.size()==original.length && doc.cds.size()==1, "Mixed GFF/embedded FASTA rows must stay in order");
		for(int i=0; i<original.length; i++){
			check(Arrays.equals(doc.lines.get(i), original[i].getBytes(StandardCharsets.UTF_8)), "Changed original GFF bytes on line "+(i+1));
		}
		check("ambiguous_codon".equals(doc.cds.get(0).reason), "Ambiguous CDS must be unsupported, not searched as no_hit");
		final ByteBuilder row=new ByteBuilder();
		AnnotateGenes.annotatedRow(row, doc.lines.get(2), doc.cds.get(0), "01234567890123456789", null);
		check(row.toString().startsWith(feature+";bbgap_status=unsupported;"), "Annotation changed existing CDS fields or attributes");
		check(row.toString().contains(";bbgap_reason=ambiguous_codon"), "Unsupported annotation needs its specific reason");
		final byte[] empty="ctg\tx\tCDS\t1\t9\t.\t+\t0\t.".getBytes(StandardCharsets.US_ASCII);
		AnnotationGff.beginAnnotation(row, empty).append("bbgap_status=no_hit");
		check(row.toString().endsWith("\tbbgap_status=no_hit"), "Empty attribute placeholder was not replaced");
		row.clear(); AnnotationGff.escaped(row, "iron; x=y%,&β");
		check(row.toString().equals("iron%3B%20x%3Dy%25%2C%26%CE%B2"), "New attribute values need reversible UTF-8 escaping");
		check("a+b|c".equals(AnnotationGff.decode("a+b%7Cc")), "GFF '+' is not URL-encoded whitespace");
		final Path header=write(dir, "table.gff", "##transl_table 4\ntable4\tx\tCDS\t1\t9\t.\t+\t0\tID=w\n");
		check("MW".equals(new String(AnnotationGff.load(header.toString(), sequences, 11).cds.get(0).translation.peptide, StandardCharsets.US_ASCII)), "File translation-table directive ignored");
		final String segmented="table4\tx\tCDS\t1\t9\t.\t+\t0\tID=x;Parent=gene\n"
			+"table4\tx\tCDS\t1\t9\t.\t+\t0\tID=y;Parent=gene\n";
		final AnnotationGff compound=AnnotationGff.load(write(dir, "compound.gff", segmented).toString(), sequences, 4);
		check(compound.cds.size()==2 && "compound_cds".equals(compound.cds.get(0).reason)
			&& "compound_cds".equals(compound.cds.get(1).reason), "Both segments of compound CDS must remain explicitly unsupported");
		reject(()->AnnotationGff.load(writeUnchecked(dir, "owned.gff", feature+";bbgap_status=assigned\n"), sequences, 11), "already has BBGAP");
		reject(()->AnnotationGff.load(writeUnchecked(dir, "short.gff", "ctg%7Cone\tx\tCDS\t1\t9\n"), sequences, 11), "nine GFF");
		reject(()->AnnotationGff.load(writeUnchecked(dir, "missing.gff", feature.replace("ctg%7Cone", "absent")+"\n"), sequences, 11), "Missing FASTA");
		reject(()->AnnotationGff.load(writeUnchecked(dir, "bounds.gff", feature.replace("\t1\t9\t", "\t1\t10\t")+"\n"), sequences, 11), "outside contig");

		final String meta="#schema_version\t1\n#kind\tbbgap_family_labels\n#label_version\tpilot_v1\n"
			+"#family_identity_sha80\t01234567890123456789\n#consensus_sha80\t98765432109876543210\n";
		final String columns="family_id\tactive_index\trep_id\tlabel_version\tdecision_class\tproduct\n";
		final String label="42\t0\trepA\tpilot_v1\thomology_supported_broad_name\tIron-interacting β protein\n"
			+"99\t1\trepB\tpilot_v1\tdomain_supported_broad_name\tBeta-interacting iron domain\n";
		final int[] ids={42, 99, 123}; final String[] reps={"repA", "repB", "repC"};
		final FamilyLabels labels=FamilyLabels.load(write(dir, "labels.tsv", meta+columns+label).toString(), ids, reps,
			"01234567890123456789", "98765432109876543210");
		check(labels.get(0).product.equals("Iron-interacting β protein")
			&& labels.get(1).product.equals("Beta-interacting iron domain")
			&& labels.version.equals("pilot_v1"),
			"Optional labels lost UTF-8, broad classes or release semantics");
		check(labels.get(2)==null, "An omitted family must stay unnamed");
		row.clear().append(feature); AnnotateGenes.appendFamilyLabel(row, labels, 0);
		check(row.toString().equals(feature+";bbgap_label_class=homology_supported_broad_name;bbgap_product=Iron-interacting%20%CE%B2%20protein"),
			"Homology-supported broad labels must retain original attributes and escape the selected product");
		row.clear().append(feature); AnnotateGenes.appendFamilyLabel(row, labels, 1);
		check(row.toString().equals(feature+";bbgap_label_class=domain_supported_broad_name;bbgap_product=Beta-interacting%20iron%20domain"),
			"Domain-supported broad labels must retain their evidence class and selected product");
		row.clear().append(feature); AnnotateGenes.appendFamilyLabel(row, labels, 2);
		check(row.toString().equals(feature+";bbgap_label_class=unnamed"), "An omitted label must not emit a product");
		row.clear().append(feature); AnnotateGenes.appendFamilyLabel(row, null, 0);
		check(row.toString().equals(feature+";bbgap_label_class=unnamed"), "An absent sidecar must not emit a product");
		final String meta2="#schema_version\t2\n#kind\tbbgap_family_labels\n"
			+"#label_version\tpilot_v2\n#family_identity_sha80\t01234567890123456789\n"
			+"#consensus_sha80\t98765432109876543210\n";
		final String columns2="family_id\tactive_index\trep_id\tdecision_class\tproduct\t"
			+"label_version\tdisplay_label\tsource_methods\tsupport_tag\tcoverage_tag\t"
			+"score_tag\tproduct_text_tag\treview_tag\tselected_hit_id\n";
		final String label2="42\t0\trepA\tprovisional_full\tIron; β=2 protein\t"
			+"pilot_v2\tprovisional_full: Iron; β=2 protein "
			+"[NCBIfam/HMMER,Swiss-Prot/DIAMOND]\t"
			+"NCBIfam/HMMER,Swiss-Prot/DIAMOND\tnative_full\tcoverage_full\t"
			+"score_passed\tsame_product_text\tnone\tNF;42\n"
			+"99\t1\trepB\tweak_provisional\ttentative homology to Beta=iron\t"
			+"pilot_v2\tweak_provisional: tentative homology to Beta=iron "
			+"[Swiss-Prot/DIAMOND]\tSwiss-Prot/DIAMOND\tswiss_partial\tquery_low\t"
			+"retained_search_hit\tswiss_product_text_disagreement\tnone\tP99\n"
			+"123\t2\trepC\tconflicting_evidence\t\tpilot_v2\t\tunknown\tunknown\t"
			+"unknown\tunknown\tnot_applicable\tnone\t\n";
		final FamilyLabels labels2=FamilyLabels.load(write(dir, "labels_v2.tsv",
			meta2+columns2+label2).toString(), ids, reps, "01234567890123456789",
			"98765432109876543210");
		check(labels2.get(0).product.equals("Iron; β=2 protein")
			&& labels2.get(1).score.equals("retained_search_hit")
			&& labels2.get(2).product==null, "Schema-2 labels lost provisional fields");
		row.clear().append(feature); AnnotateGenes.appendFamilyLabel(row, labels2, 0);
		check(row.toString().equals(feature
			+";bbgap_label_class=provisional_full;bbgap_product=Iron%3B%20%CE%B2%3D2%20protein"
			+";bbgap_label_display=provisional_full:%20Iron%3B%20%CE%B2%3D2%20protein%20[NCBIfam/HMMER%2CSwiss-Prot/DIAMOND]"
			+";bbgap_label_source=NCBIfam/HMMER%2CSwiss-Prot/DIAMOND"
			+";bbgap_label_support=native_full;bbgap_label_coverage=coverage_full"
			+";bbgap_label_score=score_passed;bbgap_label_product_text=same_product_text"
			+";bbgap_label_review=none;bbgap_label_hit=NF%3B42"),
			"Schema-2 labels should retain provenance fields and escape them");
		row.clear().append(feature); AnnotateGenes.appendFamilyLabel(row, labels2, 1);
		check(row.toString().contains(";bbgap_label_class=weak_provisional;"
			+"bbgap_product=tentative%20homology%20to%20Beta%3Diron"),
			"Weak provisional labels should emit their visible weak product");
		row.clear().append(feature); AnnotateGenes.appendFamilyLabel(row, labels2, 2);
		check(row.toString().equals(feature+";bbgap_label_class=conflicting_evidence"
			+";bbgap_label_display=;bbgap_label_source=unknown;bbgap_label_support=unknown"
			+";bbgap_label_coverage=unknown;bbgap_label_score=unknown"
			+";bbgap_label_product_text=not_applicable;bbgap_label_review=none"
			+";bbgap_label_hit="),
			"Schema-2 reviewed abstentions should stay productless but load");
		reject(()->FamilyLabels.load(writeUnchecked(dir, "labels_v2_missing.tsv",
			meta2+columns2.replace("\treview_tag", "")+label2), ids, reps,
			"01234567890123456789", "98765432109876543210"), "Missing label column");
		reject(()->FamilyLabels.load(writeUnchecked(dir, "labels_v1_provisional.tsv",
			meta+columns+"42\t0\trepA\tpilot_v1\tprovisional_full\tfoo\n"), ids, reps,
			"01234567890123456789", "98765432109876543210"), "schema 2");
		reject(()->FamilyLabels.load(writeUnchecked(dir, "labels_v2_bad_provisional.tsv",
			meta2+columns2+"42\t0\trepA\tprovisional_full\tfoo\tpilot_v2\t\t"
				+"NCBIfam/HMMER\tnative_full\tcoverage_full\tscore_passed\t"
				+"same_product_text\tnone\tncbi:0\n"), ids, reps,
			"01234567890123456789", "98765432109876543210"), "complete evidence");
		reject(()->FamilyLabels.load(writeUnchecked(dir, "stale.tsv", meta+columns+label), ids, reps, "bad", "98765432109876543210"), "identity pin");
		reject(()->FamilyLabels.load(writeUnchecked(dir, "stale_consensus.tsv", meta+columns+label), ids, reps, "01234567890123456789", "bad"), "consensus pin");
		reject(()->FamilyLabels.load(writeUnchecked(dir, "duplicate.tsv", meta+columns+label+label), ids, reps, "01234567890123456789", "98765432109876543210"), "Duplicate");
		reject(()->FamilyLabels.load(writeUnchecked(dir, "wrong_family.tsv", meta+columns+label.replace("42\t", "43\t")), ids, reps, "01234567890123456789", "98765432109876543210"), "family ID");
		reject(()->FamilyLabels.load(writeUnchecked(dir, "wrong_version.tsv", meta+columns+label.replace("pilot_v1", "pilot_v2")), ids, reps, "01234567890123456789", "98765432109876543210"), "Mixed label");
		final Path models=Files.createDirectory(dir.resolve("models"));
		final Path alias=dir.resolve("alias.gff"); Files.createLink(alias, gff);
		reject(()->{try{new AnnotateGenes(new String[]{"in="+fasta, "gff="+gff, "out="+alias, "models="+models, "t=1"});}
			catch(java.io.IOException e){throw new RuntimeException(e);}}, "aliases an input");
		check(Arrays.equals(Files.readAllBytes(gff), (String.join("\r\n", original)+"\r\n").getBytes(StandardCharsets.UTF_8)), "Rejected output alias modified the input");
		final Path model=write(models, "model.hbmc", "immutable model fixture\n"), modelAlias=dir.resolve("model_alias.gff");
		Files.createLink(modelAlias, model);
		reject(()->{try{new AnnotateGenes(new String[]{"in="+fasta, "gff="+gff, "out="+modelAlias, "models="+models, "t=1"});}
			catch(java.io.IOException e){throw new RuntimeException(e);}}, "aliases a model");
		reject(()->{try{new AnnotateGenes(new String[]{"in="+fasta, "gff="+gff, "out=stdout.gff", "evidence=stdout.tsv", "models="+models, "t=1"});}
			catch(java.io.IOException e){throw new RuntimeException(e);}}, "standard stream");
		check(new String(Files.readAllBytes(model), StandardCharsets.UTF_8).equals("immutable model fixture\n"), "Rejected alias modified the model");
		new AnnotateGenes(new String[]{fasta.toString(), gff.toString(), "out="+dir.resolve("positional.gff"), "models="+models, "labels=null", "t=1"});
		new AnnotateGenes(new String[]{gff.toString(), "in="+fasta, "out="+dir.resolve("mixed_args.gff"), "models="+models, "labels=null", "t=1"});
		final Path zippedFasta=writeGzip(dir, "genome.FNA.gz", ">table4\nATGTGATAA\n");
		final Path zippedGff=writeGzip(dir, "genes.gff3.gz", "table4\tx\tCDS\t1\t9\t.\t+\t0\tID=w\n");
		new AnnotateGenes(new String[]{zippedGff.toString(), zippedFasta.toString(), "out="+dir.resolve("compressed_args.gff"), "models="+models, "labels=null", "t=1"});
		reject(()->construct(fasta.toString(), "input="+fasta), "Multiple FASTA inputs");
		reject(()->construct("in1="+fasta, fasta.toString()), "Multiple FASTA inputs");
		reject(()->construct(gff.toString(), "gff="+gff), "Multiple GFF inputs");
		reject(()->construct("gff="+gff, gff.toString()), "Multiple GFF inputs");
		reject(()->construct("in=null", fasta.toString()), "requires a filename");
		reject(()->construct("gff=", gff.toString()), "requires a filename");
		reject(()->construct("unknown_argument.txt"), "Unknown AnnotateGenes parameter");
		final Path sibling=Files.createDirectory(dir.resolve("labels"));
		final Path siblingLabels=write(sibling, "family_labels.tsv", meta+columns+label);
		reject(()->construct(fasta.toString(), gff.toString(), "models="+models,
			"out="+siblingLabels, "t=1"), "aliases an input");
		new AnnotateGenes(new String[]{fasta.toString(), gff.toString(), "models="+models,
			"labels=null", "out="+siblingLabels, "t=1"});
		check(new String(Files.readAllBytes(siblingLabels), StandardCharsets.UTF_8).equals(meta+columns+label),
			"Constructor checks must not write the automatically discovered label file");
		new AnnotateGenes(new String[]{fasta.toString(), gff.toString(), "models="+models,
			"labels="+dir.resolve("labels.tsv"), "out="+siblingLabels, "t=1"});
		final Path siblingGzip=writeGzip(sibling, "family_labels.tsv.gz", meta+columns+label);
		reject(()->construct(fasta.toString(), gff.toString(), "models="+models,
			"out="+siblingGzip, "t=1"), "aliases an input");
		new AnnotateGenes(new String[]{fasta.toString(), gff.toString(), "models="+models,
			"out="+siblingLabels, "t=1"});
		final String oldDori=dna.Data.DORI_RESOURCES, oldNersc=dna.Data.NERSC_RESOURCES;
		final String resourceKey="annotategenes_fixture_"+dir.getFileName();
		final Path resourceModels=Files.createDirectories(dir.resolve(resourceKey).resolve("models"));
		try{
			dna.Data.DORI_RESOURCES=dir.toString()+"/";
			dna.Data.NERSC_RESOURCES=dir.resolve("absent").toString()+"/";
			reject(()->construct(fasta.toString(), gff.toString(), "models=?"+resourceKey+"/models",
				"labels=null", "out="+resourceModels.resolve("output.gff"), "t=1"), "immutable model package");
			final Path defaultModels=Files.createDirectories(dir.resolve("bbgap/models"));
			check(Paths.get(dna.Data.findPath("?bbgap/models", false)).toRealPath().equals(defaultModels.toRealPath()),
				"Default-resource fixtures require an isolated runtime without an installed BBGAP package");
			reject(()->construct(fasta.toString(), gff.toString(), "out="+dir.resolve("missing_labels.gff"), "t=1"),
				"installed BBGAP package is missing labels");
			final Path defaultLabels=write(Files.createDirectories(dir.resolve("bbgap/labels")), "family_labels.tsv", meta+columns+label);
			reject(()->construct(fasta.toString(), gff.toString(), "out="+defaultLabels, "t=1"), "aliases an input");
		}finally{dna.Data.DORI_RESOURCES=oldDori; dna.Data.NERSC_RESOURCES=oldNersc;}
		System.err.println("ANNOTATE_GENES_CLI_PASS positional=true compressed=true duplicates=true sibling_labels=true resource_lookup=true");
		System.err.println("ANNOTATE_GENES_INPUT_PASS mixed_gff=true escaping=true table_header=true compound=true labels=true rejection_cases=12 alias_preserved=true");
	}

	private static void construct(String... args){
		try{new AnnotateGenes(args);}catch(java.io.IOException e){throw new RuntimeException(e);}
	}
	private static Path writeGzip(Path dir, String name, String text) throws java.io.IOException{
		final Path path=dir.resolve(name);
		try(java.io.OutputStream out=new java.util.zip.GZIPOutputStream(Files.newOutputStream(path))){
			out.write(text.getBytes(StandardCharsets.UTF_8));
		}
		return path;
	}

	private static Path write(Path dir, String name, String text) throws java.io.IOException{
		return Files.write(dir.resolve(name), text.getBytes(StandardCharsets.UTF_8));
	}

	public static void checkProducedBbgapLabel(final String labelPath, final String familyPath,
			final int rank, final String expected) throws java.io.IOException{
		final String[] lines=new String(Files.readAllBytes(Paths.get(familyPath)),
			StandardCharsets.UTF_8).split("\n");
		int max=rank;
		for(int i=1; i<lines.length; i++){
			if(lines[i].length()<1){continue;}
			final String[] row=lines[i].split("\t", -1);
			max=Math.max(max, Integer.parseInt(row[0]));
		}
		final int[] ids=new int[max+1];
		final String[] reps=new String[max+1];
		Arrays.fill(reps, "");
		for(int i=1; i<lines.length; i++){
			if(lines[i].length()<1){continue;}
			final String[] row=lines[i].split("\t", -1);
			final int activeIndex=Integer.parseInt(row[0]);
			ids[activeIndex]=Integer.parseInt(row[1]);
			reps[activeIndex]=row[3];
		}
		final FamilyLabels labels=FamilyLabels.load(labelPath, ids, reps, "idpin", "conpin");
		final ByteBuilder row=new ByteBuilder();
		row.append("ctg\tx\tCDS\t1\t9\t.\t+\t0\tID=x;partial=10;anticodon=TCC");
		AnnotateGenes.appendFamilyLabel(row, labels, rank);
		check(row.toString().equals(expected), "Produced BBGAP label did not round-trip into GFF");
	}

	private static String writeUnchecked(Path dir, String name, String text){
		try{return write(dir, name, text).toString();}catch(java.io.IOException e){throw new RuntimeException(e);}
	}
	private static void check(boolean ok, String reason){if(!ok){throw new AssertionError(reason);}}
	private static void reject(Runnable action, String text){
		try{action.run();}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(text), "Wrong failure: "+expected); return;
		}
		throw new AssertionError("Expected rejection: "+text);
	}
}
