package test.prok;

import prok.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.charset.StandardCharsets;

/** Unknown names must fail even when a row carries a valid family label.
 * @author Brian Bushnell, Raiden */
public final class Euk18sCallAttributionTest {
	public static void main(String[] args)throws Exception{
		final Path dir=Files.createTempDirectory("euk18s_attribution_");
		final String models=dir.resolve("models.fa").toString();
		Files.write(dir.resolve("models.fa"), ">euk18S_mitochondria\nACGTACGTACGTACGTACGT\n".getBytes(StandardCharsets.US_ASCII));
		final String valid=row("rRNA", "18S,model:euk18S_mitochondria,family:euk18S")
			+row("RNA", "RNA,model:euk5S_test,family:euk5S")+row("rRNA", "16S,len:100")+row("CDS", "CDS,len:100");
		final Path calls=dir.resolve("calls.gff"), out=dir.resolve("audit.tsv");
		Files.write(calls, valid.getBytes(StandardCharsets.US_ASCII));
		new Euk18sCallAttribution().run(models, calls.toString(), out.toString());
		final String observed=new String(Files.readAllBytes(out), StandardCharsets.US_ASCII);
		check(observed.contains("TOTAL\tRNA_AND_RRNA\t3\n") && observed.contains("TOTAL\tGENERIC_RNA\t1\n")
			&& observed.contains("FAMILY\tLEGACY_RRNA\t1\n") && observed.contains("MODEL\teuk18S_mitochondria\t1\n"), "Mixed known families and legacy rows must conserve the literal three-RNA fixture");
		reject(dir, models, row("rRNA", "18S,model:euk18S_unknown,family:euk18S"), "absent from its actual library");
		reject(dir, models, row("RNA", "RNA,model:euk18S_mitochondria"), "lacks explicit family attribution");
		reject(dir, models, row("RNA", "RNA,model:euk18S_mitochondria,family:unregistered"), "Unknown emitted RNA family");
		reject(dir, models, row("RNA", "RNA,model:euk18S_mitochondria,family:euk5S"), "cannot be attributed to a different family");
		reject(dir, models, row("RNA", "RNA,family:euk18S,family:euk5S"), "Missing/duplicate family value");
		final String[] names={"euk18S_universal", "euk18S_fungi", "euk18S_plant", "euk18S_invertebrate", "euk18S_vertebrate_other", "euk18S_vertebrate_mammalian", "euk18S_protist", "euk18S_protozoa_component_001", "euk18S_protozoa_component_003"};
		Euk18sRuntimeConfig.validateDefaultNames(names);
		names[5]="euk18S_animal";
		Euk18sRuntimeConfigTest.reject(()->Euk18sRuntimeConfig.validateDefaultNames(names), "roster mismatch");
		Euk18sRuntimeConfigTest.reject(()->Euk18sRuntimeConfig.validateDefaultNames(new String[10]), "complete declared roster");
		System.out.println("EUK18S_ATTRIBUTION_TEST_PASS mixed_families legacy unknown_model unknown_family contradictory_family duplicate_family default_roster");
	}
	private static void reject(Path dir, String models, String row, String expected)throws Exception{
		final Path file=Files.createTempFile(dir, "invalid", ".gff");Files.write(file, row.getBytes(StandardCharsets.US_ASCII));
		Euk18sRuntimeConfigTest.reject(()->new Euk18sCallAttribution().run(models, file.toString(), file+".audit"), expected);
	}
	private static String row(String type, String attrs){return "seq\tfixture\t"+type+"\t1\t100\t.\t+\t.\t"+attrs+"\n";}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
