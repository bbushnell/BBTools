package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/** Small fixtures for BBGAP full-evidence triage.
 * @author Yelan
 */
public final class BbgapFullEvidenceTriageTest {

	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Specify one fresh fixture directory");
		final Path dir=Files.createDirectory(Paths.get(args[0]));
		write(dir, "family_identity.tsv",
			"active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues\n"
			+"7\t42\t7\trepA\t100\t1000\n"
			+"8\t43\t8\trepB\t50\t500\n"
			+"9\t44\t9\trepC\t10\t100\n"
			+"10\t45\t10\trepD\t5\t50\n"
			+"11\t46\t11\trepE\t4\t40\n");
		write(dir, "consensus.faa",
			">repA active_index=7 family_id=42\nMKTIIALSY\n"
			+">repB active_index=8 family_id=43\nMKTIIALSYF\n"
			+">repC active_index=9 family_id=44\nMKTIIALSYFG\n"
			+">repD active_index=10 family_id=45\nMKTIIALSYFGH\n"
			+">repE active_index=11 family_id=46\nMKTIIALSYFGHI\n");
		write(dir, "targets.tsv",
			targetHeader()
			+"42\t7\trepA\tNF_A\tHMM_A\t1e-5\t200\t1e-5\t190\t1\t100\t90\t100\tequivalog\tY\tNCBIFAM\tgeneA\t\tproduct A\n"
			+"43\t8\trepB\tNF_B\tHMM_B\t1e-5\t200\t1e-5\t190\t1\t100\t90\t100\tequivalog\tY\tNCBIFAM\tgeneB\t\tproduct B\n"
			+"45\t10\trepD\tNF_D1\tHMM_D1\t1e-5\t200\t1e-5\t190\t1\t100\t90\t100\tequivalog\tY\tNCBIFAM\tgeneD1\t\tproduct D1\n"
			+"45\t10\trepD\tNF_D2\tHMM_D2\t1e-5\t150\t1e-5\t140\t1\t100\t90\t100\tequivalog\tY\tNCBIFAM\tgeneD2\t\tproduct D2\n"
			+"46\t11\trepE\tNF_E1\tHMM_E1\t1e-5\t200\t1e-5\t190\t1\t250\t90\t100\tequivalog\tY\tNCBIFAM\tgeneE1\t\tproduct E1\n"
			+"46\t11\trepE\tNF_E2\tHMM_E2\t1e-5\t150\t1e-5\t140\t1\t100\t90\t100\tequivalog\tY\tNCBIFAM\tgeneE2\t\tproduct E2\n");
		write(dir, "domains.tsv",
			domainHeader()
			+"42\t7\trepA\t9\tNF_A\tHMM_A\t100\t1\t1\t1e-5\t200\t1e-5\t1e-5\t195\t1\t40\t1\t5\t1\t5\t0.400000000000\t0.555555555556\t100\t90\tequivalog\tNCBIFAM\tgeneA\tshort A\n"
			+"42\t7\trepA\t9\tNF_A\tHMM_A\t100\t2\t2\t1e-5\t200\t1e-5\t1e-5\t190\t1\t90\t1\t9\t1\t9\t0.900000000000\t1.000000000000\t100\t90\tequivalog\tNCBIFAM\tgeneA\tproduct A\n"
			+"43\t8\trepB\t10\tNF_B\tHMM_B\t100\t1\t1\t1e-5\t200\t1e-5\t1e-5\t190\t1\t40\t1\t5\t1\t5\t0.400000000000\t0.500000000000\t100\t90\tequivalog\tNCBIFAM\tgeneB\tproduct B\n"
			+"45\t10\trepD\t12\tNF_D1\tHMM_D1\t100\t1\t1\t1e-5\t200\t1e-5\t1e-5\t190\t1\t90\t1\t9\t1\t9\t0.900000000000\t0.818181818182\t100\t90\tequivalog\tNCBIFAM\tgeneD1\tproduct D1\n"
			+"45\t10\trepD\t12\tNF_D2\tHMM_D2\t100\t1\t1\t1e-5\t150\t1e-5\t1e-5\t140\t1\t90\t1\t9\t1\t9\t0.900000000000\t0.818181818182\t100\t90\tequivalog\tNCBIFAM\tgeneD2\tproduct D2\n"
			+"46\t11\trepE\t13\tNF_E1\tHMM_E1\t100\t1\t1\t1e-5\t200\t1e-5\t1e-5\t190\t1\t90\t1\t9\t1\t9\t0.900000000000\t0.818181818182\t250\t90\tequivalog\tNCBIFAM\tgeneE1\tproduct E1\n"
			+"46\t11\trepE\t13\tNF_E2\tHMM_E2\t100\t1\t1\t1e-5\t150\t1e-5\t1e-5\t140\t1\t90\t1\t9\t1\t9\t0.900000000000\t0.818181818182\t100\t90\tequivalog\tNCBIFAM\tgeneE2\tproduct E2\n");
		write(dir, "swiss.tsv",
			swissHeader()
			+"43\t8\trepB\tP1\tFULL_B\t70\t10\t10\t10\t1.000000000000\t1.000000000000\t1\t10\t1\t10\t1e-9\t99\tsp|P1|FULL_B reviewed B\n");
		BbgapFullEvidenceTriage.main(new String[]{"families="+dir.resolve("family_identity.tsv"),
			"consensus="+dir.resolve("consensus.faa"), "targets="+dir.resolve("targets.tsv"),
			"domains="+dir.resolve("domains.tsv"), "swiss="+dir.resolve("swiss.tsv"),
			"out="+dir.resolve("full"), "audit="+dir.resolve("audit.tsv"), "expectedfamilies=5"});
		final String review=read(dir, "full.review_candidates.tsv");
		final String models=read(dir, "full.ncbi_model_candidates.tsv");
		final String summary=read(dir, "full.review_summary.tsv");
		final String audit=read(dir, "audit.tsv");
		checkWidths(review, 23, "review");
		checkWidths(models, 20, "model");
		checkWidths(audit, 10, "audit");
		checkReviewPositions(review);
		checkModelPositions(models);
		check(review.contains("42\t7\trepA\t100\t1000\t9\tgated_full_length\t1\t1\t0\t0\tNF_A"),
			"family 42 should be a native gated full-length candidate");
		check(review.contains("43\t8\trepB\t50\t500\t10\tswissprot_full_length_review\t1\t0\t1\t1"),
			"family 43 should retain full-length Swiss-Prot support without native eligibility");
		check(review.contains("44\t9\trepC\t10\t100\t11\tno_evidence\t0\t0\t0\t0"),
			"family 44 should be explicit no-evidence");
		check(review.contains("45\t10\trepD\t5\t50\t12\tnative_multimodel_review\t2\t2\t0\t0"),
			"family 45 should be a multimodel review candidate");
		check(review.contains("46\t11\trepE\t4\t40\t13\tgated_full_length\t2\t1\t0\t0\tNF_E2"),
			"family 46 should find the lower-scoring eligible model");
		check(models.contains("43\t8\trepB\tNF_B\t200\t100\tY\t190\t90\tY\t0.400000000000\t0.500000000000\tN\tN"),
			"same-hit native cutoff pass must still fail the mutual coverage gate");
		check(summary.contains("total_families\t5\n")
			&& summary.contains("gated_full_length\t2\n")
			&& summary.contains("native_multimodel_review\t1\n")
			&& summary.contains("swissprot_full_length_review\t1\n")
			&& summary.contains("no_evidence\t1\n"), "summary counts changed");
		check(audit.contains("\n46\t11\tNF_E1\t200\t250\tN\tN\t1\tNF_E2\tproduct E2\n"),
			"top-raw masking audit should report the hidden lower-scoring eligible model");
		checkBadDomainLength(dir);
		checkMissingTarget(dir);
		checkDecimalBoundaries(dir);
		System.err.println("BBGAP_FULL_EVIDENCE_TRIAGE_PASS rows=true same_hit=true no_hit=true decimals=true");
	}

	private static void checkDecimalBoundaries(final Path dir) throws Exception{
		final String domains=read(dir, "domains.tsv").replace("0.900000000000", "0.800000000000");
		final String swiss=read(dir, "swiss.tsv").replace("1.000000000000", "0.800000000000");
		write(dir, "boundary.domains.tsv", domains);
		write(dir, "boundary.swiss.tsv", swiss);
		runBoundary(dir, "boundary");
		check(read(dir, "boundary.review_summary.tsv").equals(read(dir, "full.review_summary.tsv")),
			"Exact 0.800000000000 must pass both native-model and Swiss coverage gates");
		check(rowByModel(read(dir, "boundary.ncbi_model_candidates.tsv"), "42", "NF_A")[13].equals("Y"),
			"Exact coverage lost the lower-scoring eligible domain");
		write(dir, "below.domains.tsv", domains.replace("0.800000000000", "0.799999999999"));
		write(dir, "below.swiss.tsv", swiss.replace("0.800000000000", "0.799999999999"));
		runBoundary(dir, "below");
		final String below=read(dir, "below.review_summary.tsv");
		check(below.contains("gated_full_length\t0\n") && below.contains("native_multimodel_review\t0\n")
			&& below.contains("swissprot_full_length_review\t0\n") && below.contains("partial_or_below_cutoff\t4\n"),
			"Values strictly below the coverage cutoff must fail; do not use an epsilon");
	}

	private static void runBoundary(final Path dir, final String prefix){
		BbgapFullEvidenceTriage.main(new String[]{"families="+dir.resolve("family_identity.tsv"),
			"consensus="+dir.resolve("consensus.faa"), "targets="+dir.resolve("targets.tsv"),
			"domains="+dir.resolve(prefix+".domains.tsv"), "swiss="+dir.resolve(prefix+".swiss.tsv"),
			"out="+dir.resolve(prefix), "expectedfamilies=5"});
	}

	private static void checkReviewPositions(final String review){
		final String[] row42=rowByFamily(review, "42"), row43=rowByFamily(review, "43");
		final String[] row44=rowByFamily(review, "44");
		check(row42[11].equals("NF_A") && row42[13].equals("190")
			&& row42[16].equals("product A") && row42[17].equals("NF_A"),
			"family 42 eligible and top-model fields shifted");
		check(row43[11].equals("") && row43[16].equals("") && row43[17].equals("NF_B")
			&& row43[19].equals("P1") && row43[22].equals("sp|P1|FULL_B reviewed B"),
			"family 43 null eligible-model fields shifted");
		for(int i=11; i<row44.length; i++){
			check(row44[i].equals(""), "family 44 no-evidence column "+i+" was not empty");
		}
	}

	private static void checkModelPositions(final String models){
		final String[] rowA=rowByModel(models, "42", "NF_A");
		check(rowA[7].equals("190") && rowA[10].equals("0.900000000000")
			&& rowA[11].equals("1.000000000000") && rowA[13].equals("Y")
			&& rowA[19].equals("product A"),
			"NF_A should use the lower-scoring domain that passes coverage");
		final String[] rowB=rowByModel(models, "43", "NF_B");
		check(rowB[7].equals("190") && rowB[10].equals("0.400000000000")
			&& rowB[12].equals("N") && rowB[13].equals("N"),
			"NF_B should fail coverage and same-hit native eligibility");
	}

	private static void checkBadDomainLength(final Path dir) throws Exception{
		write(dir, "domains_bad_len.tsv",
			domainHeader()
			+"45\t10\trepD\t11\tNF_D1\tHMM_D1\t100\t1\t1\t1e-5\t200\t1e-5\t1e-5\t190\t1\t90\t1\t9\t1\t9\t0.900000000000\t0.818181818182\t100\t90\tequivalog\tNCBIFAM\tgeneD1\tproduct D1\n");
		boolean ok=false;
		try{
			BbgapFullEvidenceTriage.main(new String[]{"families="+dir.resolve("family_identity.tsv"),
				"consensus="+dir.resolve("consensus.faa"), "targets="+dir.resolve("targets.tsv"),
				"domains="+dir.resolve("domains_bad_len.tsv"), "swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("bad_len"), "expectedfamilies=5"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Consensus length mismatch");
		}
		check(ok, "Bad domain consensus_len should crash loudly");
	}

	private static void checkMissingTarget(final Path dir) throws Exception{
		write(dir, "domains_missing_target.tsv",
			domainHeader()
			+"45\t10\trepD\t12\tNF_MISSING\tHMM_MISSING\t100\t1\t1\t1e-5\t200\t1e-5\t1e-5\t190\t1\t90\t1\t9\t1\t9\t0.900000000000\t0.818181818182\t100\t90\tequivalog\tNCBIFAM\tgeneM\tproduct M\n");
		boolean ok=false;
		try{
			BbgapFullEvidenceTriage.main(new String[]{"families="+dir.resolve("family_identity.tsv"),
				"consensus="+dir.resolve("consensus.faa"), "targets="+dir.resolve("targets.tsv"),
				"domains="+dir.resolve("domains_missing_target.tsv"), "swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("missing_target"), "expectedfamilies=5"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Domain row without target row");
		}
		check(ok, "Domain without a compact target row should crash loudly");
	}

	private static String targetHeader(){
		return "family_id\tactive_index\tconsensus_id\tncbi_accession\thmm_name\t"
			+"full_evalue\tfull_score\tbest_domain_evalue\tbest_domain_score\treported_domains\t"
			+"sequence_cutoff\tdomain_cutoff\thmm_length\tfamily_type\tfor_naming\tsource\t"
			+"gene_symbol\tec_numbers\tproduct_name\n";
	}

	private static String domainHeader(){
		return "family_id\tactive_index\tconsensus_id\tconsensus_len\tncbi_accession\t"
			+"hmm_name\thmm_len\tdomain_index\tdomain_count\tfull_evalue\tfull_score\t"
			+"conditional_evalue\tindependent_evalue\tdomain_score\thmm_from\thmm_to\t"
			+"ali_from\tali_to\tenv_from\tenv_to\thmm_cov\tconsensus_cov\tsequence_cutoff\t"
			+"domain_cutoff\tfamily_type\tsource\tgene_symbol\tproduct_name\n";
	}

	private static String swissHeader(){
		return "family_id\tactive_index\tconsensus_id\tuniprot_accession\tuniprot_entry\t"
			+"identity\taligned_len\tquery_len\tsubject_len\tquery_cov\tsubject_cov\tqstart\t"
			+"qend\tsstart\tsend\tevalue\tbitscore\ttitle\n";
	}

	private static void write(final Path dir, final String name, final String text) throws Exception{
		Files.write(dir.resolve(name), text.getBytes(StandardCharsets.UTF_8));
	}

	private static String read(final Path dir, final String name) throws Exception{
		return new String(Files.readAllBytes(dir.resolve(name)), StandardCharsets.UTF_8);
	}

	private static void checkWidths(final String text, final int width, final String name){
		final String[] lines=text.split("\n");
		check(lines.length>1, "No data rows in "+name);
		for(int i=0; i<lines.length; i++){
			final int found=lines[i].split("\t", -1).length;
			check(found==width, name+" line "+(i+1)+" width "+found+" != "+width);
		}
	}

	private static String[] rowByFamily(final String text, final String familyId){
		final String[] lines=text.split("\n");
		for(int i=1; i<lines.length; i++){
			final String[] split=lines[i].split("\t", -1);
			if(split[0].equals(familyId)){return split;}
		}
		throw new AssertionError("Missing family "+familyId);
	}

	private static String[] rowByModel(final String text, final String familyId, final String accession){
		final String[] lines=text.split("\n");
		for(int i=1; i<lines.length; i++){
			final String[] split=lines[i].split("\t", -1);
			if(split[0].equals(familyId) && split[3].equals(accession)){return split;}
		}
		throw new AssertionError("Missing family/model "+familyId+"/"+accession);
	}

	private static void check(final boolean ok, final String reason){
		if(!ok){throw new AssertionError(reason);}
	}
}
