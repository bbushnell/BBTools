package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/** Small fixtures for Java-owned BBGAP label-pilot compaction.
 * @author Yelan
 */
public final class BbgapLabelEvidenceTest {

	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Specify one fresh fixture directory");
		final Path dir=Files.createDirectory(Paths.get(args[0]));
		write(dir, "pilot_families.tsv",
			"active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues\tconsensus_aa\tpilot_stratum\n"
			+"7\t42\t7\trepA\t100\t1000\t100\thigh_count\n");
		write(dir, "hmm_PGAP.tsv",
			"#ncbi_accession\tsource_identifier\tlabel\tsequence_cutoff\tdomain_cutoff\thmm_length\t"
			+"family_type\tfor_structural_annotation\tfor_naming\tfor_AMRFinder\tproduct_name\t"
			+"gene_symbol\tgene_synonyms\tec_numbers\tgo_terms\tpmids\ttaxonomic_range\t"
			+"taxonomic_range_name\ttaxonomic_rank_name\tn_refseq_protein_hits\tsource\tname_orig\thmm_name\tcomment\n"
			+"NF000035.1\t\tDnaB\t500\t500\t274\tequivalog\tY\tY\tY\tDNA helicase DnaB\tdnaB\t\t3.6.4.12\t\t\t\t\t\t0\tNCBIFAM\tDnaB\tDnaB\t\n");
		write(dir, "pilot.ncbi.tblout",
			"repA - HMM_A NF000035.1 1e-04 45.2 0.0 2e-04 44.0 0.0 3.0 3 0 0 3 3 2 1 family_id=42 active_index=7 has spaces\n");
		write(dir, "pilot.ncbi.domtblout",
			"repA - 100 HMM_A NF000035.1 274 1e-04 45.2 0.0 1 2 5e-05 7e-05 44.0 0.0 10 100 5 95 4 96 0.96 family_id=42 active_index=7 first\n"
			+"repA - 100 HMM_A NF000035.1 274 1e-04 45.2 0.0 2 2 6e-05 8e-05 25.0 0.0 150 200 40 90 38 92 0.91 family_id=42 active_index=7 second\n");
		write(dir, "pilot.swissprot.tsv",
			"repA\tsp|P12345|DNAB_ECOLI\t99.0\t80\t100\t200\t1\t80\t11\t90\t1e-40\t200\t"
			+"Replicative DNA helicase DnaB OS=Escherichia coli\n");
		final Path out=Files.createDirectory(dir.resolve("out"));
		BbgapLabelEvidence.main(new String[]{"in="+dir, "out="+out, "mode=both",
			"expectedtargets=1", "expecteddomains=2", "expectedswissprot=1"});
		final String targets=read(out, "pilot.ncbi.targets.compact.tsv");
		final String domains=read(out, "pilot.ncbi.domains.compact.tsv");
		final String swiss=read(out, "pilot.swissprot.compact.tsv");
		final String validation=read(out, "compact_validation.tsv");
		check(!targets.contains("sequence_cutoff_bits") && !domains.contains("domain_cutoff_bits"),
			"NCBI native cutoffs must not be divided into synthetic bit columns");
		check(targets.contains("\t45.2\t2e-04\t44.0\t2\t500\t500\t274\t"),
			"Target row did not retain native NCBI cutoffs");
		check(validation.contains("ncbi_target_sequence_cutoff_pass_rows\t0\trecorded\tOK"),
			"NF000035-like 45.2 score must not pass a native 500 cutoff");
		check(domains.contains("\t0.332116788321\t0.910000000000\t500\t500\t"),
			"Domain coverage math changed");
		check(swiss.contains("42\t7\trepA\tP12345\tDNAB_ECOLI\t99.0\t80\t100\t200\t0.800000000000\t0.400000000000\t"),
			"Swiss-Prot join or coverage math changed");
		write(dir, "bad.tblout", "repA - HMM_A NF404 1e-04 45.2 0.0 2e-04 44.0 0.0 1.0 1 0 0 1 1 1 1 family_id=42 active_index=7\n");
		reject(()->BbgapLabelEvidence.main(new String[]{"in="+dir, "tbl="+dir.resolve("bad.tblout"),
			"out="+dir.resolve("bad"), "mode=compact"}), "Missing NCBI metadata");
		write(dir, "bad.domtblout", "repA - 100 HMM_A NF000035.1 274 1e-04 45.2 0.0 1 1 5e-05 7e-05 44.0 0.0 10 100 5 95 4 96 0.96 active_index=7\n");
		reject(()->BbgapLabelEvidence.main(new String[]{"in="+dir, "dom="+dir.resolve("bad.domtblout"),
			"out="+dir.resolve("bad"), "mode=compact"}), "Missing family_id");
		System.err.println("BBGAP_LABEL_EVIDENCE_PASS native_cutoffs=true compact_tables=true validation=true");
	}

	private static void write(final Path dir, final String name, final String text) throws Exception{
		Files.write(dir.resolve(name), text.getBytes(StandardCharsets.UTF_8));
	}

	private static String read(final Path dir, final String name) throws Exception{
		return new String(Files.readAllBytes(dir.resolve(name)), StandardCharsets.UTF_8);
	}

	private static void check(final boolean ok, final String reason){
		if(!ok){throw new AssertionError(reason);}
	}

	private static void reject(final ThrowingRunnable action, final String text) throws Exception{
		try{
			action.run();
		}catch(Throwable expected){
			check(expected.toString().contains(text), "Wrong failure: "+expected);
			return;
		}
		throw new AssertionError("Expected rejection: "+text);
	}

	private interface ThrowingRunnable {
		void run() throws Exception;
	}
}
