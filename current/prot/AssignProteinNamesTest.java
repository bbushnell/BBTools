package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import prok.AnnotateGenesInputTest;

/** Fixtures for BBGAP automatic protein naming.
 * @author Yelan
 */
public final class AssignProteinNamesTest {

	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Specify one fresh fixture directory");
		final Path dir=Files.createDirectory(Paths.get(args[0]));
		writeInputs(dir);
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss.tsv"),
			"reviewed="+dir.resolve("reviewed.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("selected.tsv"),
			"evidence="+dir.resolve("evidence.tsv"),
			"unresolved="+dir.resolve("unresolved.tsv"),
			"summary="+dir.resolve("summary.tsv"),
			"audit="+dir.resolve("audit.tsv"),
			"labelversion=test_v1",
			"identitypin=idpin",
			"consensuspin=conpin",
			"maxaudit=10",
			"seed=7"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss.tsv"),
			"reviewed="+dir.resolve("reviewed.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("case_selected.tsv"),
			"evidence="+dir.resolve("case_evidence.tsv"),
			"unresolved="+dir.resolve("case_unresolved.tsv"),
			"summary="+dir.resolve("case_summary.tsv"),
			"audit="+dir.resolve("case_audit.tsv"),
			"labelversion=test_case",
			"identitypin=idpin",
			"consensuspin=conpin",
			"maxaudit=10",
			"seed=7",
			"ignorecase=t"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss.tsv"),
			"reviewed="+dir.resolve("reviewed.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("hierarchy_selected.tsv"),
			"evidence="+dir.resolve("hierarchy_evidence.tsv"),
			"unresolved="+dir.resolve("hierarchy_unresolved.tsv"),
			"summary="+dir.resolve("hierarchy_summary.tsv"),
			"membercoverage="+dir.resolve("hierarchy_member_coverage.tsv"),
			"memberunresolved="+dir.resolve("hierarchy_member_unresolved.tsv"),
			"labelversion=test_hierarchy",
			"identitypin=idpin",
			"consensuspin=conpin",
			"maxaudit=10",
			"conflictgroups="+dir.resolve("hierarchy_conflict_groups.tsv"),
			"mode=hierarchy",
			"seed=7"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss.tsv"),
			"reviewed="+dir.resolve("reviewed.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("nearquery_selected.tsv"),
			"evidence="+dir.resolve("nearquery_evidence.tsv"),
			"unresolved="+dir.resolve("nearquery_unresolved.tsv"),
			"summary="+dir.resolve("nearquery_summary.tsv"),
			"labelversion=test_nearquery",
			"identitypin=idpin",
			"consensuspin=conpin",
			"maxaudit=10",
			"mode=hierarchy",
			"seed=7",
			"ncbinearquerydomains="+dir.resolve("coverage_domains.tsv"),
			"ncbinearquery=t"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss.tsv"),
			"reviewed="+dir.resolve("reviewed.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("nearquery_provisional_selected.tsv"),
			"evidence="+dir.resolve("nearquery_provisional_evidence.tsv"),
			"unresolved="+dir.resolve("nearquery_provisional_unresolved.tsv"),
			"summary="+dir.resolve("nearquery_provisional_summary.tsv"),
			"labelversion=test_nearquery_provisional",
			"identitypin=idpin",
			"consensuspin=conpin",
			"maxaudit=10",
			"mode=hierarchy",
			"seed=7",
			"ncbinearquerydomains="+dir.resolve("coverage_domains.tsv"),
			"ncbinearquery=t",
			"provisional=t"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss.tsv"),
			"reviewed="+dir.resolve("reviewed.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("provisional_selected.tsv"),
			"evidence="+dir.resolve("provisional_evidence.tsv"),
			"unresolved="+dir.resolve("provisional_unresolved.tsv"),
			"summary="+dir.resolve("provisional_summary.tsv"),
			"labelversion=test_provisional",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy",
			"seed=7",
			"provisional=t"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss.tsv"),
			"reviewed="+dir.resolve("reviewed.tsv"),
			"reviewedevidence="+dir.resolve("reviewed_evidence.tsv"),
			"reviewedevidenceout="+dir.resolve("reviewed_evidence_validated.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("reviewed_evidence_selected.tsv"),
			"evidence="+dir.resolve("reviewed_evidence_evidence.tsv"),
			"unresolved="+dir.resolve("reviewed_evidence_unresolved.tsv"),
			"summary="+dir.resolve("reviewed_evidence_summary.tsv"),
			"labelversion=test_reviewed_evidence",
			"identitypin=idpin",
			"consensuspin=conpin",
			"seed=7"});
		write(dir, "swiss_reviewed.tsv",
			read(dir, "reviewed.tsv")
			+"7\t7\tgated_full_length\tY\tdomain_supported_broad_name\t"
			+"LysR family transcriptional regulator\tprimary\tcompeting\t"
			+"manual exception\tpreserve Swiss review\n");
		write(dir, "swiss_reviewed_hits.tsv", read(dir, "swiss.tsv")
			+"7\t7\trep7\tP77333\tP77333_LYSR\t49.7\t100\t100\t100\t0.951298701299\t0.983277591973\t1\t100\t1\t100\t1e-50\t276\tsp|P77333|P77333_LYSR Putative transcriptional regulator OS=Species OX=1 GN=lysR PE=1 SV=1\n");
		write(dir, "swiss_reviewed_evidence.tsv",
			reviewedEvidenceSidecar(reviewedEvidenceRow(6, 6, "rep6", "85.5")
				+swissReviewedEvidenceRow(7, 7, "rep7")));
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss_reviewed_hits.tsv"),
			"reviewed="+dir.resolve("swiss_reviewed.tsv"),
			"reviewedevidence="+dir.resolve("swiss_reviewed_evidence.tsv"),
			"reviewedevidenceout="+dir.resolve("swiss_reviewed_evidence_validated.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("swiss_reviewed_evidence_selected.tsv"),
			"evidence="+dir.resolve("swiss_reviewed_evidence_evidence.tsv"),
			"unresolved="+dir.resolve("swiss_reviewed_evidence_unresolved.tsv"),
			"summary="+dir.resolve("swiss_reviewed_evidence_summary.tsv"),
			"labelversion=test_swiss_reviewed_evidence",
			"identitypin=idpin",
			"consensuspin=conpin",
			"seed=7"});
		write(dir, "native_reviewed.tsv",
			read(dir, "reviewed.tsv")
			+"7\t7\tgated_full_length\tY\tdomain_supported_broad_name\t"
			+"native curated protein\tprimary\tcompeting\t"
			+"manual exception\tpreserve native review\n");
		write(dir, "native_reviewed_evidence.tsv",
			reviewedEvidenceSidecar(reviewedEvidenceRow(6, 6, "rep6", "85.5")
				+nativeReviewedEvidenceRow(7, 7, "rep7")));
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss.tsv"),
			"reviewed="+dir.resolve("native_reviewed.tsv"),
			"reviewedevidence="+dir.resolve("native_reviewed_evidence.tsv"),
			"reviewedevidenceout="+dir.resolve("native_reviewed_evidence_validated.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("native_reviewed_evidence_selected.tsv"),
			"evidence="+dir.resolve("native_reviewed_evidence_evidence.tsv"),
			"unresolved="+dir.resolve("native_reviewed_evidence_unresolved.tsv"),
			"summary="+dir.resolve("native_reviewed_evidence_summary.tsv"),
			"labelversion=test_native_reviewed_evidence",
			"identitypin=idpin",
			"consensuspin=conpin",
			"seed=7"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("family_identity.tsv"),
			"ncbi="+dir.resolve("ncbi.tsv"),
			"swiss="+dir.resolve("swiss.tsv"),
			"reviewed="+dir.resolve("native_reviewed.tsv"),
			"reviewedevidence="+dir.resolve("native_reviewed_evidence.tsv"),
			"reviewedevidenceout="+dir.resolve("native_reviewed_evidence_provisional_validated.tsv"),
			"synonyms="+dir.resolve("synonyms.tsv"),
			"out="+dir.resolve("native_reviewed_evidence_provisional_selected.tsv"),
			"evidence="+dir.resolve("native_reviewed_evidence_provisional_evidence.tsv"),
			"unresolved="+dir.resolve("native_reviewed_evidence_provisional_unresolved.tsv"),
			"summary="+dir.resolve("native_reviewed_evidence_provisional_summary.tsv"),
			"labelversion=test_native_reviewed_evidence_provisional",
			"identitypin=idpin",
			"consensuspin=conpin",
			"seed=7",
			"mode=hierarchy",
			"provisional=t"});
		writeEggnogInputs(dir);
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("eggnog_family_identity.tsv"),
			"ncbi="+dir.resolve("eggnog_ncbi.tsv"),
			"swiss="+dir.resolve("eggnog_swiss.tsv"),
			"reviewed="+dir.resolve("eggnog_reviewed.tsv"),
			"out="+dir.resolve("eggnog_selected.tsv"),
			"evidence="+dir.resolve("eggnog_evidence.tsv"),
			"unresolved="+dir.resolve("eggnog_unresolved.tsv"),
			"summary="+dir.resolve("eggnog_summary.tsv"),
			"eggnoghomology="+dir.resolve("eggnog_homology.tsv"),
			"eggnoghomologyout="+dir.resolve("eggnog_homology_validated.tsv"),
			"labelversion=test_eggnog",
			"identitypin=idpin",
			"consensuspin=conpin",
			"seed=7",
			"mode=hierarchy",
			"provisional=t"});
		write(dir, "eggnog_homology_relaxed.tsv", eggnogHomologySidecar(
			eggnogHomologyRow(4, "DUF1815", "DUF1815_2_101", "DUF1815@2|A-1",
				"-", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80", "29609749")
			+eggnogHomologyRow(6, "DUF1,DUF2", "DUF1_1_10,DUF2_1_10",
				"DUF2_1_10,DUF1_1_10", "DUF1_1_10,DUF2_1_10",
				"DUF1@2|A-1,DUF2@2|A-1", "-", "-", "-", "-", "f",
				"seed_pfam_width_fraction_0.50_0.80", "6")
			+eggnogHomologyRow(7, "DUF3", "DUF3_1_10", "DUF3@2|A-1",
				"yghQ", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80", "7")
			+eggnogHomologyRow(8, "ActD", "ActD_1_10", "ActD@2|A-1",
				"actD", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80", "8")
			+eggnogHomologyRow(9, "DUF6", "DUF6_1_10", "DUF6@2|A-1",
				"enzyme", "K00001", "1.1.1.1", "K00001=kinase", "f",
				"seed_pfam_width_fraction_lt_0.50", "9")
			+eggnogHomologyRow(10, "DUF1,DUF2", "DUF1_1_10,DUF2_1_10",
				"DUF1_1_11,DUF2_1_10", "DUF1_1_10,DUF2_1_10",
				"DUF1@2|A-1,DUF2@2|A-1", "-", "-", "-", "-", "f",
				"seed_pfam_width_fraction_0.50_0.80", "10")
			+eggnogHomologyRow(11, "DUF1,DUF2", "DUF1_1_10,DUF2_1_10",
				"DUF1_1_10", "DUF1_1_10,DUF2_1_10",
				"DUF1@2|A-1,DUF2@2|A-1", "-", "-", "-", "-", "f",
				"seed_pfam_width_fraction_0.50_0.80", "11")
			+eggnogHomologyRow(12, "DUF1,DUF2", "DUF1_1_10,DUF2_1_10",
				"DUF1_1_10,DUF1_1_10,DUF2_1_10", "DUF1_1_10,DUF2_1_10",
				"DUF1@2|A-1,DUF2@2|A-1", "-", "-", "-", "-", "f",
				"seed_pfam_width_fraction_0.50_0.80", "12")));
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("eggnog_family_identity.tsv"),
			"ncbi="+dir.resolve("eggnog_ncbi.tsv"),
			"swiss="+dir.resolve("eggnog_swiss.tsv"),
			"reviewed="+dir.resolve("eggnog_reviewed.tsv"),
			"out="+dir.resolve("eggnog_relaxed_selected.tsv"),
			"evidence="+dir.resolve("eggnog_relaxed_evidence.tsv"),
			"unresolved="+dir.resolve("eggnog_relaxed_unresolved.tsv"),
			"summary="+dir.resolve("eggnog_relaxed_summary.tsv"),
			"eggnoghomology="+dir.resolve("eggnog_homology_relaxed.tsv"),
			"eggnoghomologyout="+dir.resolve("eggnog_homology_relaxed_validated.tsv"),
			"labelversion=test_eggnog_relaxed",
			"identitypin=idpin",
			"consensuspin=conpin",
			"seed=7",
			"mode=hierarchy",
			"provisional=t",
			"eggnogdomainrelaxed=t"});
		final String selected=read(dir, "selected.tsv");
		final String evidence=read(dir, "evidence.tsv");
		final String unresolved=read(dir, "unresolved.tsv");
		final String summary=read(dir, "summary.tsv");
		final String caseSelected=read(dir, "case_selected.tsv");
		final String caseEvidence=read(dir, "case_evidence.tsv");
		final String caseUnresolved=read(dir, "case_unresolved.tsv");
		final String caseSummary=read(dir, "case_summary.tsv");
		final String hierarchySelected=read(dir, "hierarchy_selected.tsv");
		final String hierarchyEvidence=read(dir, "hierarchy_evidence.tsv");
		final String hierarchyUnresolved=read(dir, "hierarchy_unresolved.tsv");
		final String hierarchySummary=read(dir, "hierarchy_summary.tsv");
		final String hierarchyConflictGroups=read(dir, "hierarchy_conflict_groups.tsv");
		final String hierarchyMemberCoverage=read(dir, "hierarchy_member_coverage.tsv");
		final String hierarchyMemberUnresolved=read(dir, "hierarchy_member_unresolved.tsv");
		final String nearquerySelected=read(dir, "nearquery_selected.tsv");
		final String nearqueryEvidence=read(dir, "nearquery_evidence.tsv");
		final String nearqueryUnresolved=read(dir, "nearquery_unresolved.tsv");
		final String nearquerySummary=read(dir, "nearquery_summary.tsv");
		final String nearqueryProvisionalSelected=read(dir,
			"nearquery_provisional_selected.tsv");
		final String provisionalSelected=read(dir, "provisional_selected.tsv");
		final String provisionalEvidence=read(dir, "provisional_evidence.tsv");
		final String provisionalSummary=read(dir, "provisional_summary.tsv");
		final String reviewedEvidenceSelected=read(dir, "reviewed_evidence_selected.tsv");
		final String swissReviewedEvidenceSelected=read(dir,
			"swiss_reviewed_evidence_selected.tsv");
		final String nativeReviewedEvidenceSelected=read(dir,
			"native_reviewed_evidence_selected.tsv");
		final String eggnogSelected=read(dir, "eggnog_selected.tsv");
		final String eggnogEvidence=read(dir, "eggnog_evidence.tsv");
		final String eggnogHomology=read(dir, "eggnog_homology_validated.tsv");
		final String eggnogSummary=read(dir, "eggnog_summary.tsv");
		final String eggnogRelaxedSelected=read(dir, "eggnog_relaxed_selected.tsv");
		final String eggnogRelaxedHomology=read(dir,
			"eggnog_homology_relaxed_validated.tsv");
		final String eggnogRelaxedSummary=read(dir, "eggnog_relaxed_summary.tsv");
		check(selected.contains("\n1\t1\trep1\tspecific_full_length\tproduct A\ttest_v1\n"),
			"exact NCBI/Swiss-Prot agreement should select product A");
		check(selected.contains("\n2\t2\trep2\tspecific_full_length\tDNA gyrase subunit A\ttest_v1\n"),
			"explicit synonym map should select the canonical product");
		check(selected.contains("\n5\t5\trep5\tconflicting_evidence\t\ttest_v1\n"),
			"reviewed conflicts should be preserved with an empty product");
		check(selected.contains("\n6\t6\trep6\tdomain_supported_broad_name\tpeptide ABC transporter permease\ttest_v1\n"),
			"reviewed broad names should be preserved");
		check(selected.contains("\n18\t18\trep18\tspecific_full_length\tproduct R\ttest_v1\n"),
			"exact coverage-boundary Swiss-Prot rows should be retained");
		check(unresolved.contains("\n19\t19\trep19\tunnamed\t\t"
			+"insufficient_full_length_agreement\t1\t0\tNF_S\t\tproduct S\t\tnot_selected\t"),
			"below-boundary Swiss-Prot rows should be filtered");
		check(evidence.contains("\n3\t3\trep3\tconflicting_evidence\t\tcompeting_ncbi_products\t2\t1\tNF_C1,NF_C2\tP3\tproduct C1|product C2\tproduct C1\tnot_selected\tP3\tproduct C1"),
			"NCBI product conflicts should retain accessions and products");
		check(unresolved.contains("\n4\t4\trep4\tconflicting_evidence\t\tcompeting_swissprot_products\t1\t2\tNF_D\tP4A,P4B\tproduct D\tproduct D1|product D2\tnot_selected\tP4A\tproduct D1"),
			"Swiss-Prot product conflicts should be unresolved");
		check(unresolved.contains("\n7\t7\trep7\tunnamed\t\tinsufficient_full_length_agreement\t0\t1\t\tP7\t\tproduct G\tnot_selected\tP7\tproduct G"),
			"for_naming=N NCBI products should not count");
		check(summary.contains("total_families\t34\n")
			&& summary.contains("selected_auto\t4\n")
			&& summary.contains("selected_auto_exact\t4\n")
			&& summary.contains("selected_auto_specific\t4\n")
			&& summary.contains("selected_auto_with_swissprot\t4\n")
			&& summary.contains("reviewed_overrides\t2\n")
			&& summary.contains("reviewed_named\t1\n")
			&& summary.contains("reviewed_unresolved\t1\n")
			&& summary.contains("unresolved_conflict\t14\n")
			&& summary.contains("unresolved_insufficient\t14\n"),
			"summary counts changed");
		check(caseSelected.contains("\n8\t8\trep8\tspecific_full_length\tproduct H\ttest_case\n"),
			"ignorecase=t should select case-only NCBI/Swiss-Prot agreement");
		check(caseEvidence.contains("\n8\t8\trep8\tspecific_full_length\tproduct H\t"
			+"case_insensitive_ncbi_swissprot_full_length_agreement\t1\t1\tNF_H\tP8\t"
			+"product H\tProduct H\tagreement\tP8\tProduct H"),
			"case-folded exact calls should preserve original product text");
		check(!caseUnresolved.contains("\n8\t8\trep8\t"),
			"case-folded agreement should not remain unresolved");
		check(caseSummary.contains("total_families\t34\n")
			&& caseSummary.contains("selected_auto\t6\n")
			&& caseSummary.contains("selected_auto_exact\t6\n")
			&& caseSummary.contains("selected_auto_specific\t6\n")
			&& caseSummary.contains("selected_auto_with_swissprot\t6\n")
			&& caseSummary.contains("unresolved_conflict\t12\n")
			&& caseSummary.contains("unresolved_insufficient\t14\n")
			&& caseSummary.contains("ignore_case\ttrue\n"),
			"ignorecase summary counts changed");
		check(caseSelected.contains("\n9\t9\trep9\tspecific_full_length\tDNA gyrase subunit A\ttest_case\n"),
			"ignorecase=t should fold explicit synonym sources consistently");
		check(caseEvidence.contains("\n9\t9\trep9\tspecific_full_length\tDNA gyrase subunit A\t"
			+"case_insensitive_ncbi_swissprot_full_length_agreement\t1\t1\tNF_I\tP9\t"
			+"DNA TOPOISOMERASE (ATP-HYDROLYZING) SUBUNIT A\tDNA gyrase subunit A\t"
			+"agreement\tP9\tDNA gyrase subunit A"),
			"case-folded synonyms should preserve original mapped product text");
		check(evidence.contains("\n2\t2\trep2\tspecific_full_length\tDNA gyrase subunit A\t"
			+"exact_ncbi_swissprot_full_length_agreement\t1\t1\tNF_B\tP2\t"
			+"DNA topoisomerase (ATP-hydrolyzing) subunit A\tDNA gyrase subunit A\t"
			+"agreement\tP2\tDNA gyrase subunit A"),
			"synonym-collapsed exact calls should preserve original product text");
		check(hierarchySelected.contains("\n10\t10\trep10\tspecific_full_length\t"
			+"product J\ttest_hierarchy\n"),
			"hierarchy mode should select single specific NCBI products without Swiss-Prot");
		check(hierarchySelected.contains("\n11\t11\trep11\tdomain_supported_broad_name\t"
			+"subfamily protein\ttest_hierarchy\n"),
			"hierarchy mode should select broad names from subfamily_domain rows");
		check(hierarchyEvidence.contains("\n11\t11\trep11\tdomain_supported_broad_name\t"
			+"subfamily protein\tncbi_broad_rank55_native_name\t1\t1\tNF_K\tP11\t"
			+"subfamily protein\tsubfamily offshoot\tmismatch\tP11\tsubfamily offshoot"),
			"hierarchy mode should report same-query Swiss-Prot mismatches");
		check(hierarchyEvidence.contains("\n17\t17\trep17\tspecific_full_length\t"
			+"product Q\tncbi_specific_rank70_native_name\t1\t3\tNF_Q\tP17B,P17C,P17A\t"
			+"product Q\tproduct R|product S|product Q\tmixed\tP17B\tproduct R\t98\t"
			+"1.000000000000\t1.000000000000\t200\tP17C\tproduct S\t97\t"
			+"1.000000000000\t1.000000000000\t190\t10.000\t1.052632\tP17A\t"
			+"product Q\t96\t1.000000000000\t1.000000000000\t100\tP17B\tproduct R\t"
			+"98\t1.000000000000\t1.000000000000\t200\t-100.000\t0.500000\n"),
			"hierarchy mode should report mixed Swiss-Prot matches and stronger mismatches");
		check(hierarchyEvidence.contains("\n12\t12\trep12\tspecific_full_length\t"
			+"exception product\tncbi_specific_rank77_native_name\t2\t0\tNF_L1,NF_L2\t\t"
			+"exception product|domain product\t\tabsent\t"),
			"hierarchy mode should prefer higher-precedence NCBI products");
		check(hierarchyUnresolved.contains("\n13\t13\trep13\tconflicting_evidence\t\t"
			+"competing_top_ncbi_products\t2\t0\tNF_M1,NF_M2\t\tproduct M1|product M2\t\t"
			+"not_selected\t"),
			"hierarchy mode should leave same-precedence product splits unresolved");
		check(hierarchySelected.contains("\n14\t14\trep14\tdomain_supported_broad_name\t"
			+"DUF1234 domain-containing protein\ttest_hierarchy\n")
			&& hierarchySelected.contains("\n16\t16\trep16\tdomain_supported_broad_name\t"
			+"DUF4321 domain-containing protein\ttest_hierarchy\n"),
			"hierarchy mode should deterministically broaden hypothetical top rows");
		check(hierarchyUnresolved.contains("\n15\t15\trep15\tunnamed\t\t"
			+"uninformative_native_ncbi_product\t1\t0\tNF_O1,NF_O2\t\t"
			+"hypothetical protein\t\tnot_selected\t"),
			"hierarchy mode should not count generic hypothetical protein as named");
		check(hierarchySelected.contains("\n29\t29\trep29\tspecific_full_length\t"
			+"product AC\ttest_hierarchy\n"),
			"hierarchy mode should keep exact-boundary full-coverage NCBI rows eligible");
		check(hierarchySummary.contains("selected_auto\t14\n")
			&& hierarchySummary.contains("selected_auto_exact\t0\n")
			&& hierarchySummary.contains("selected_auto_specific\t11\n")
			&& hierarchySummary.contains("selected_auto_broad\t3\n")
			&& hierarchySummary.contains("selected_auto_swissonly\t0\n")
			&& hierarchySummary.contains("selected_auto_with_swissprot\t9\n")
			&& hierarchySummary.contains("selected_auto_ncbi_only\t5\n")
			&& hierarchySummary.contains("selected_auto_subfamily_domain\t1\n")
			&& hierarchySummary.contains("unresolved_conflict\t7\n")
			&& hierarchySummary.contains("unresolved_insufficient\t11\n")
			&& hierarchySummary.contains("mode\thierarchy\n"),
			"hierarchy summary counts changed");
		check(hierarchyConflictGroups.contains("rank\tconflict_source\treason\t"
			+"family_count\tproduct_count\tfamily_ids\tnormalized_products\t"
			+"display_products\n")
			&& hierarchyConflictGroups.contains("\tncbi_top_tier\t"
				+"competing_top_ncbi_products\t2\t2\t13,30\tproduct M1|product M2\t"
				+"product M1|product M2\n")
			&& hierarchyConflictGroups.contains("\tswissprot\t"
				+"competing_swissprot_products\t2\t2\t28,31\tproduct AB|product AC\t"
				+"product AB|product AC\n")
			&& hierarchyConflictGroups.contains("\tncbi_top_tier\t"
				+"competing_top_ncbi_products\t1\t2\t32\ta|b|c\ta|b|c\n")
			&& hierarchyConflictGroups.contains("\tncbi_top_tier\t"
				+"competing_top_ncbi_products\t1\t2\t33\ta|b|c\ta|b|c\n"),
			"hierarchy conflict grouping should count recurring effective keysets");
		check(hierarchyMemberCoverage.contains("bucket\tdecision_class\trule_or_reason\t"
			+"family_count\tmembers\traw_residues\ttotal_families\ttotal_members\t"
			+"total_raw_residues\tmember_fraction\traw_residue_fraction\n")
			&& hierarchyMemberCoverage.contains("named_reviewed\tdomain_supported_broad_name"
				+"\t\t1\t10\t100\t34\t550\t5500\t0.018181818182\t0.018181818182\n")
			&& hierarchyMemberCoverage.contains("unresolved_conflict\tconflicting_evidence"
				+"\tcompeting_top_ncbi_products\t5\t130\t1300\t34\t550\t5500\t"
				+"0.236363636364\t0.236363636364\n")
			&& hierarchyMemberCoverage.contains("unresolved_reviewed\tconflicting_evidence"
				+"\treviewed_override\t1\t70\t700\t34\t550\t5500\t0.127272727273"
				+"\t0.127272727273\n"),
			"reference-package member coverage should retain weighted denominators");
		checkMemberCoverageTotals(hierarchyMemberCoverage, 34, 550, 5500);
		check(hierarchyMemberUnresolved.contains("rank\tfamily_id\tactive_index\trep_id\t"
			+"decision_class\trule_or_reason\tmembers\traw_residues\tncbi_products\t"
			+"swissprot_products\n1\t13\t13\trep13\tconflicting_evidence\t"
			+"competing_top_ncbi_products\t90\t900")
			&& hierarchyMemberUnresolved.contains("\n2\t15\t15\trep15\tunnamed\t"
				+"uninformative_native_ncbi_product\t80\t800")
			&& hierarchyMemberUnresolved.contains("\n3\t5\t5\trep5\tconflicting_evidence\t"
				+"reviewed_override\t70\t700")
			&& hierarchyMemberUnresolved.contains("\n4\t3\t3\trep3\tconflicting_evidence\t"
				+"competing_top_ncbi_products\t10\t100")
			&& hierarchyMemberUnresolved.contains("\n5\t7\t7\trep7\tunnamed\t"
				+"insufficient_native_ncbi_evidence\t10\t100"),
			"member unresolved report should sort by reference-package members");
		check(nearquerySelected.contains("\n20\t20\trep20\tspecific_full_length\t"
			+"product T\ttest_nearquery\n"),
			"near-query fallback should rescue one Swiss-supported equivalog/exception key");
		check(nearqueryEvidence.contains("\n20\t20\trep20\tspecific_full_length\t"
			+"product T\tncbi_near_query_swiss_supported_native_name\t1\t1\t"
			+"NF_T,NF_T2\tP20\tproduct T\tproduct T\tagreement\t"),
			"near-query fallback should report all same-product NCBI alternatives");
		check(nearqueryUnresolved.contains("\n21\t21\trep21\tunnamed\t\t"
			+"insufficient_native_ncbi_evidence\t0\t1\t\tP21")
			&& nearqueryUnresolved.contains("\n22\t22\trep22\tunnamed\t\t"
			+"insufficient_native_ncbi_evidence\t0\t1\t\tP22")
			&& nearqueryUnresolved.contains("\n23\t23\trep23\tunnamed\t\t"
			+"insufficient_native_ncbi_evidence\t0\t1\t\tP23")
			&& nearqueryUnresolved.contains("\n24\t24\trep24\tunnamed\t\t"
			+"insufficient_native_ncbi_evidence\t0\t1\t\tP24")
			&& nearqueryUnresolved.contains("\n25\t25\trep25\tunnamed\t\t"
			+"insufficient_native_ncbi_evidence\t0\t1\t\tP25")
			&& nearqueryUnresolved.contains("\n26\t26\trep26\tunnamed\t\t"
			+"insufficient_native_ncbi_evidence\t0\t1\t\tP26")
			&& nearqueryUnresolved.contains("\n27\t27\trep27\tunnamed\t\t"
			+"insufficient_native_ncbi_evidence\t0\t1\t\tP27")
			&& nearqueryUnresolved.contains("\n28\t28\trep28\tconflicting_evidence\t\t"
			+"competing_swissprot_products\t0\t2\t\tP28A,P28B"),
			"near-query fallback should preserve coverage, family-type, naming and Swiss gates");
		check(nearquerySummary.contains("selected_auto\t15\n")
			&& nearquerySummary.contains("selected_auto_specific\t12\n")
			&& nearquerySummary.contains("selected_auto_ncbinearquery\t1\n")
			&& nearquerySummary.contains("selected_auto_with_swissprot\t10\n")
			&& nearquerySummary.contains("selected_auto_ncbi_only\t5\n")
			&& nearquerySummary.contains("unresolved_conflict\t7\n")
			&& nearquerySummary.contains("unresolved_insufficient\t10\n"),
			"near-query summary counts changed");
		check(nearqueryProvisionalSelected.contains("\n20\t20\trep20\t"
			+"specific_full_length\tproduct T\ttest_nearquery_provisional\t"
			+"product T\tNCBIfam/HMMER,Swiss-Prot/DIAMOND\tnative_near_query\tquery_low\t"
			+"score_passed\tnot_applicable\tnone\t\n"),
			"Accepted near-query rows should remain visibly low-query in schema 2");
		check(!selected.contains("display_label") &&
			provisionalSelected.contains("#schema_version\t2\n")
			&& provisionalSelected.contains("display_label\tsource_methods\t"
				+"support_tag\tcoverage_tag\tscore_tag\tproduct_text_tag\treview_tag\t"
				+"selected_hit_id\n"),
			"provisional labels must be schema-2 opt-in only");
		check(provisionalSelected.contains("\n7\t7\trep7\tprovisional_full\t"
			+"product G\ttest_provisional\tprovisional_full: product G "
			+"[NCBIfam/HMMER,Swiss-Prot/DIAMOND]\t"
			+"NCBIfam/HMMER,Swiss-Prot/DIAMOND\tswiss_full\tcoverage_full\t"
			+"retained_search_hit\t"),
			"Swiss-only unresolved rows should become explicit provisional labels");
		check(provisionalSelected.contains("\n15\t15\trep15\t"
			+"weak_provisional_placeholder\thypothetical protein\ttest_provisional"),
			"uninformative placeholders should stay visibly weak");
		check(provisionalSelected.contains("\n15\t15\trep15\t"
			+"weak_provisional_placeholder\thypothetical protein\ttest_provisional\t"
			+"weak_provisional_placeholder: hypothetical protein [NCBIfam/HMMER]\t"
			+"NCBIfam/HMMER\tnative_placeholder\t"),
			"Native placeholders should not be promoted to partial support");
		check(provisionalSelected.contains("\n20\t20\trep20\tweak_provisional\t"
			+"product T\ttest_provisional\tweak_provisional: product T "
			+"[NCBIfam/HMMER,Swiss-Prot/DIAMOND]\tNCBIfam/HMMER,Swiss-Prot/DIAMOND\t"
			+"native_partial\tquery_low\tscore_passed\tsame_product_text\t"
			+"none\tncbi:1"),
			"partial NCBI rows should outrank agreeing full Swiss hits without a "
			+"source-priority tag");
		check(provisionalSelected.contains("\n34\t34\trep34\t"
			+"weak_provisional_placeholder\thypothetical protein\ttest_provisional\t"
			+"weak_provisional_placeholder: hypothetical protein [NCBIfam/HMMER]\t"
			+"NCBIfam/HMMER\tnative_placeholder\tcoverage_full\tscore_passed\t"
			+"ncbi_and_swiss_product_text_disagreement\tnone\tncbi:1"),
			"Generic fallback should not prefer a 7.9-bit named hit over a "
			+"175.2-bit placeholder hit");
		check(provisionalEvidence.contains("\n34\t34\trep34\tncbi:0\tN\t"
			+"NCBIfam/HMMER\tNF608_BAD\tstress response YncL\t")
			&& provisionalEvidence.contains("\n34\t34\trep34\tncbi:1\tY\t"
				+"NCBIfam/HMMER\tNF608_HIGH\thypothetical protein\t")
			&& provisionalEvidence.contains("\n34\t34\trep34\tswiss:0\tN\t"
				+"Swiss-Prot/DIAMOND\tP608\tUPF0126 family protein\t"),
			"Generic fallback should keep selected and alternate placeholder evidence");
		check(provisionalSelected.contains("\n5\t5\trep5\tconflicting_evidence\t\t"
			+"test_provisional\t"),
			"reviewed abstentions should remain abstentions in schema 2");
		check(provisionalSummary.contains("provisional_labels\t18\n")
			&& provisionalSummary.contains("provisional_informative\t16\n")
			&& provisionalSummary.contains("provisional_placeholder\t2\n")
			&& provisionalSummary.contains("provisional_source_priority_partial_ncbi_over_full_swiss\t3\n"),
			"provisional summary counts changed");
		check(reviewedEvidenceSelected.contains("#schema_version\t2\n")
			&& !reviewedEvidenceSelected.contains("provisional_full")
			&& !reviewedEvidenceSelected.contains("weak_provisional"),
			"reviewed evidence should opt selected labels into schema 2 without provisional rows");
		check(reviewedEvidenceSelected.contains("\n6\t6\trep6\t"
			+"domain_supported_broad_name\tpeptide ABC transporter permease\t"
			+"test_reviewed_evidence\tpeptide ABC transporter permease\teggNOG/DIAMOND\t"
			+"reviewed_eggnog_homology_domain\tcoverage_full\tselected_seed_passed\t"
			+"eggnog_domain_family\tmanual_noevidence_reviewβ\teggnog:38557254;DUF3097\n"),
			"reviewed evidence should project raw UTF-8 evidence into selected labels");
		check(reviewedEvidenceSelected.contains("\n5\t5\trep5\tconflicting_evidence\t\t"
			+"test_reviewed_evidence\t\tunknown\tunknown\tunknown\tunknown\t"
			+"not_applicable\tnone\t\n"),
			"reviewed abstentions without a sidecar row should remain schema-2 empty rows");
		check(read(dir, "reviewed_evidence_validated.tsv").equals(
			read(dir, "reviewed_evidence.tsv")),
			"validated reviewed evidence must round-trip the sidecar bytes");
		check(swissReviewedEvidenceSelected.contains("\n7\t7\trep7\t"
			+"domain_supported_broad_name\tLysR family transcriptional regulator\t"
			+"test_swiss_reviewed_evidence\tLysR family transcriptional regulator\t"
			+"Swiss-Prot/DIAMOND\tretained_full_swissprot\tcoverage_full\t"
			+"retained_bitscore\tfull_coverage_family\tscientific_review\t"
			+"swiss:P77333\n"),
			"Swiss-Prot reviewed evidence should project into schema-2 labels");
		check(read(dir, "swiss_reviewed_evidence_validated.tsv").equals(
			read(dir, "swiss_reviewed_evidence.tsv")),
			"Swiss-Prot reviewed evidence must round-trip the sidecar bytes");
		check(nativeReviewedEvidenceSelected.contains("\n7\t7\trep7\t"
			+"domain_supported_broad_name\tnative curated protein\t"
			+"test_native_reviewed_evidence\tnative curated protein\t"
			+"NCBIfam/HMMER\tretained_full_native\tcoverage_full\t"
			+"native_bitscore\tnative_family\tscientific_review\tncbi:0\n"),
			"NCBIfam/HMMER reviewed evidence should project into schema-2 labels");
		check(read(dir, "native_reviewed_evidence_validated.tsv").equals(
			read(dir, "native_reviewed_evidence.tsv")),
			"NCBIfam/HMMER reviewed evidence must round-trip the sidecar bytes");
		AnnotateGenesInputTest.checkProducedBbgapLabel(
			dir.resolve("reviewed_evidence_selected.tsv").toString(),
			dir.resolve("family_identity.tsv").toString(),
			6,
			"ctg\tx\tCDS\t1\t9\t.\t+\t0\tID=x;partial=10;anticodon=TCC"
			+";bbgap_label_class=domain_supported_broad_name"
			+";bbgap_product=peptide%20ABC%20transporter%20permease"
			+";bbgap_label_display=peptide%20ABC%20transporter%20permease"
			+";bbgap_label_source=eggNOG/DIAMOND"
			+";bbgap_label_support=reviewed_eggnog_homology_domain"
			+";bbgap_label_coverage=coverage_full"
			+";bbgap_label_score=selected_seed_passed"
			+";bbgap_label_product_text=eggnog_domain_family"
			+";bbgap_label_review=manual_noevidence_review%CE%B2"
			+";bbgap_label_hit=eggnog:38557254%3BDUF3097");
		AnnotateGenesInputTest.checkProducedBbgapLabel(
			dir.resolve("swiss_reviewed_evidence_selected.tsv").toString(),
			dir.resolve("family_identity.tsv").toString(),
			7,
			"ctg\tx\tCDS\t1\t9\t.\t+\t0\tID=x;partial=10;anticodon=TCC"
			+";bbgap_label_class=domain_supported_broad_name"
			+";bbgap_product=LysR%20family%20transcriptional%20regulator"
			+";bbgap_label_display=LysR%20family%20transcriptional%20regulator"
			+";bbgap_label_source=Swiss-Prot/DIAMOND"
			+";bbgap_label_support=retained_full_swissprot"
			+";bbgap_label_coverage=coverage_full"
			+";bbgap_label_score=retained_bitscore"
			+";bbgap_label_product_text=full_coverage_family"
			+";bbgap_label_review=scientific_review"
			+";bbgap_label_hit=swiss:P77333");
		AnnotateGenesInputTest.checkProducedBbgapLabel(
			dir.resolve("native_reviewed_evidence_selected.tsv").toString(),
			dir.resolve("family_identity.tsv").toString(),
			7,
			"ctg\tx\tCDS\t1\t9\t.\t+\t0\tID=x;partial=10;anticodon=TCC"
			+";bbgap_label_class=domain_supported_broad_name"
			+";bbgap_product=native%20curated%20protein"
			+";bbgap_label_display=native%20curated%20protein"
			+";bbgap_label_source=NCBIfam/HMMER"
			+";bbgap_label_support=retained_full_native"
			+";bbgap_label_coverage=coverage_full"
			+";bbgap_label_score=native_bitscore"
			+";bbgap_label_product_text=native_family"
			+";bbgap_label_review=scientific_review"
			+";bbgap_label_hit=ncbi:0");
		check(eggnogSelected.contains("\n1\t1\terep1\tspecific_full_length\tproduct A\ttest_eggnog\t")
			&& eggnogSelected.contains("\n2\t2\terep2\tdomain_supported_broad_name\treviewed protein\ttest_eggnog\t")
			&& eggnogSelected.contains("\n3\t3\terep3\tweak_provisional\tproduct C\ttest_eggnog\t"),
			"eggNOG homology must preserve existing accepted, reviewed and provisional rows");
		check(eggnogSelected.contains("\n4\t4\terep4\tweak_provisional\t"
				+"putative DUF1815 domain-containing protein\ttest_eggnog\t"
				+"weak_provisional: same-seed homology to putative DUF1815 "
				+"domain-containing protein [eggNOG/DIAMOND]\teggNOG/DIAMOND\t"
				+"same_seed_domain_family_homology\tseed_pfam_width_fraction_ge_0.80\t"
				+"selected_seed_60_80_80_pass\tpfam_og_root_match_nonfunctional_lead\t"
				+"none\teggnog:29609749;DUF1815\n")
			&& eggnogSelected.contains("\n5\t5\terep5\tweak_provisional\t"
				+"putative 12TM_1 domain-containing protein\ttest_eggnog\t")
			&& eggnogSelected.contains("\n12\t12\terep12\tweak_provisional\t"
				+"putative DUF9 domain-containing protein\ttest_eggnog\t")
			&& eggnogSelected.contains("\n13\t13\terep13\tweak_provisional\t"
				+"putative DUF10 domain-containing protein\ttest_eggnog\t")
			&& eggnogSelected.contains("\n14\t14\terep14\tweak_provisional\t"
				+"putative DUF11 domain-containing protein\ttest_eggnog\t"),
			"eligible same-seed Pfam/OG roots should allow locus-like Preferred_name placeholders");
		check(!eggnogSelected.contains("\n6\t6\terep6\t")
			&& !eggnogSelected.contains("\n7\t7\terep7\t")
			&& !eggnogSelected.contains("\n8\t8\terep8\t")
			&& !eggnogSelected.contains("\n9\t9\terep9\t")
			&& !eggnogSelected.contains("\n10\t10\terep10\t")
			&& !eggnogSelected.contains("\n11\t11\terep11\t"),
			"valid ineligible eggNOG rows should not emit labels");
		check(!eggnogEvidence.contains("raw_eggnog"),
			"raw eggNOG rows belong in the separate evidence output");
		check(eggnogSummary.contains("provisional_labels\t6\n")
			&& eggnogSummary.contains("provisional_eggnog_homology\t5\n"),
			"eggNOG homology summary counts changed");
		check(eggnogRelaxedSelected.contains("\n6\t6\terep6\tweak_provisional\t"
				+"putative DUF1/DUF2 domain-containing protein\ttest_eggnog_relaxed\t")
			&& eggnogRelaxedSelected.contains(
				"\tpfam_og_multi_root_match_nonfunctional_lead\tnone\teggnog:6;DUF1,DUF2\n")
			&& eggnogRelaxedSelected.contains("\n7\t7\terep7\tweak_provisional\t"
				+"putative DUF3 domain-containing protein\ttest_eggnog_relaxed\t")
			&& eggnogRelaxedSelected.contains(
				"\tpfam_og_root_match_symbol_ambiguity\tunresolved_symbol_ambiguity\t"
				+"eggnog:7;DUF3\n")
			&& eggnogRelaxedSelected.contains("\n8\t8\terep8\tweak_provisional\t"
				+"putative ActD domain-containing protein\ttest_eggnog_relaxed\t")
				&& !eggnogRelaxedSelected.contains("\n9\t9\terep9\t")
				&& !eggnogRelaxedSelected.contains("\n10\t10\terep10\t")
				&& !eggnogRelaxedSelected.contains("\n11\t11\terep11\t")
				&& !eggnogRelaxedSelected.contains("\n12\t12\terep12\t"),
				"relaxed eggNOG domain rules should require exact Pfam token multisets");
		check(eggnogRelaxedHomology.contains("\n6\t6\terep6\t")
			&& eggnogRelaxedHomology.contains("\tY\tunnamed\t"
				+"insufficient_native_ncbi_evidence\t.\n")
			&& eggnogRelaxedHomology.contains("\n7\t7\terep7\t")
			&& eggnogRelaxedHomology.contains("\n9\t9\terep9\t")
			&& eggnogRelaxedHomology.contains("\tN\tunnamed\t"
				+"insufficient_native_ncbi_evidence\tpreferred_name_functional_lead\n"),
			"relaxed eggNOG raw evidence should retain selected and rejected rows");
		check(eggnogRelaxedSummary.contains("provisional_labels\t5\n")
			&& eggnogRelaxedSummary.contains("provisional_eggnog_homology\t4\n"),
			"relaxed eggNOG summary counts changed");
		check(eggnogHomology.contains("\n4\t4\terep4\t")
			&& eggnogHomology.contains("\tY\tunnamed\tinsufficient_native_ncbi_evidence\t.\n")
			&& eggnogHomology.contains("\n1\t1\terep1\t")
			&& eggnogHomology.contains("\tN\tspecific_full_length\t\texisting_product\n")
			&& eggnogHomology.contains("enzyme β-chain note\tK00001\t1.1.1.1")
			&& eggnogHomology.contains("\tpreferred_name_functional_lead\n")
			&& eggnogHomology.contains("\tK00001,K00002\t-\t"
				+"K00001=uncharacterized protein,K00002=kinase\tf\t")
			&& eggnogHomology.contains("\tK00003\t-\t"
				+"K00003=uncharacterized methyltransferase\tf\t")
			&& eggnogHomology.contains("\tko_functional_lead\n")
			&& eggnogHomology.contains("\tN\tunnamed\tinsufficient_native_ncbi_evidence\t"
				+"absent_pfam\n"),
			"separate eggNOG homology evidence should retain raw UTF-8 rows and reasons");
		AnnotateGenesInputTest.checkProducedBbgapLabel(
			dir.resolve("eggnog_selected.tsv").toString(),
			dir.resolve("eggnog_family_identity.tsv").toString(),
			4,
			"ctg\tx\tCDS\t1\t9\t.\t+\t0\tID=x;partial=10;anticodon=TCC"
			+";bbgap_label_class=weak_provisional"
			+";bbgap_product=putative%20DUF1815%20domain-containing%20protein"
			+";bbgap_label_display=weak_provisional:%20same-seed%20homology%20to%20"
			+"putative%20DUF1815%20domain-containing%20protein%20[eggNOG/DIAMOND]"
			+";bbgap_label_source=eggNOG/DIAMOND"
			+";bbgap_label_support=same_seed_domain_family_homology"
			+";bbgap_label_coverage=seed_pfam_width_fraction_ge_0.80"
			+";bbgap_label_score=selected_seed_60_80_80_pass"
			+";bbgap_label_product_text=pfam_og_root_match_nonfunctional_lead"
			+";bbgap_label_review=none"
			+";bbgap_label_hit=eggnog:29609749%3BDUF1815");
		checkWidths(read(dir, "audit.tsv"), 42, "audit");
		checkReviewedEvidenceRejections(dir);
		checkEggnogHomologyRejections(dir);
		checkDuplicateOutput(dir);
		checkDuplicateSynonym(dir);
		checkDuplicateFoldedSynonym(dir);
		checkBadSwissTitle(dir);
		checkBadSwissMetric(dir);
		checkBadNcbiNearHeader(dir);
		checkBadNcbiNearMetric(dir);
		checkBadNcbiNearFamilyIdentity(dir);
		checkZeroSwissBits(dir);
		checkSwissOnlyFallback(dir);
		checkBroadConflictRules(dir);
		checkNcbiNearRequiresHierarchy(dir);
		checkUnknownFamilyType(dir);
		checkBadMemberCounts(dir);
		checkZeroRawMemberCoverage(dir);
		System.err.println("ASSIGN_PROTEIN_NAMES_PASS exact=true synonyms=true reviewed=true unresolved=true");
	}

	private static void writeInputs(final Path dir) throws Exception{
		write(dir, "family_identity.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\n"
			+"1\t1\t0\trep1\t10\t100\n"
			+"2\t2\t0\trep2\t10\t100\n"
			+"3\t3\t0\trep3\t10\t100\n"
			+"4\t4\t0\trep4\t10\t100\n"
			+"5\t5\t0\trep5\t70\t700\n"
			+"6\t6\t0\trep6\t10\t100\n"
			+"7\t7\t0\trep7\t10\t100\n"
			+"8\t8\t0\trep8\t10\t100\n"
			+"9\t9\t0\trep9\t10\t100\n"
			+"10\t10\t0\trep10\t10\t100\n"
			+"11\t11\t0\trep11\t10\t100\n"
			+"12\t12\t0\trep12\t10\t100\n"
			+"13\t13\t0\trep13\t90\t900\n"
			+"14\t14\t0\trep14\t10\t100\n"
			+"15\t15\t0\trep15\t80\t800\n"
			+"16\t16\t0\trep16\t10\t100\n"
			+"17\t17\t0\trep17\t10\t100\n"
			+"18\t18\t0\trep18\t10\t100\n"
			+"19\t19\t0\trep19\t10\t100\n"
			+"20\t20\t0\trep20\t10\t100\n"
			+"21\t21\t0\trep21\t10\t100\n"
			+"22\t22\t0\trep22\t10\t100\n"
			+"23\t23\t0\trep23\t10\t100\n"
			+"24\t24\t0\trep24\t10\t100\n"
			+"25\t25\t0\trep25\t10\t100\n"
			+"26\t26\t0\trep26\t10\t100\n"
			+"27\t27\t0\trep27\t10\t100\n"
			+"28\t28\t0\trep28\t10\t100\n"
			+"29\t29\t0\trep29\t10\t100\n"
			+"30\t30\t0\trep30\t10\t100\n"
			+"31\t31\t0\trep31\t10\t100\n"
			+"32\t32\t0\trep32\t10\t100\n"
			+"33\t33\t0\trep33\t10\t100\n"
			+"34\t34\t0\trep34\t10\t100\n");
		write(dir, "ncbi.tsv", modelHeader()
			+"1\t1\trep1\tNF_A\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneA\tproduct A\n"
			+"2\t2\trep2\tNF_B\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneB\tDNA topoisomerase (ATP-hydrolyzing) subunit A\n"
			+"3\t3\trep3\tNF_C1\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneC1\tproduct C1\n"
			+"3\t3\trep3\tNF_C2\t190\t100\tY\t180\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneC2\tproduct C2\n"
			+"4\t4\trep4\tNF_D\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneD\tproduct D\n"
			+"5\t5\trep5\tNF_E\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneE\tproduct E\n"
			+"6\t6\trep6\tNF_F\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneF\toligopeptide ABC transporter permease\n"
			+"7\t7\trep7\tNF_G\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tN\tNCBIFAM\tgeneG\tproduct G\n"
			+"8\t8\trep8\tNF_H\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneH\tproduct H\n"
			+"9\t9\trep9\tNF_I\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneI\tDNA TOPOISOMERASE (ATP-HYDROLYZING) SUBUNIT A\n"
			+"10\t10\trep10\tNF_J\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog_domain\tY\tNCBIFAM\tgeneJ\tproduct J\n"
			+"11\t11\trep11\tNF_K\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tsubfamily_domain\tY\tNCBIFAM\tgeneK\tsubfamily protein\n"
			+"12\t12\trep12\tNF_L1\t190\t100\tY\t180\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\texception\tY\tNCBIFAM\tgeneL1\texception product\n"
			+"12\t12\trep12\tNF_L2\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tdomain\tY\tNCBIFAM\tgeneL2\tdomain product\n"
			+"13\t13\trep13\tNF_M1\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneM1\tproduct M1\n"
			+"13\t13\trep13\tNF_M2\t190\t100\tY\t180\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\thypoth_equivalog\tY\tNCBIFAM\tgeneM2\tproduct M2\n"
			+"14\t14\trep14\tNF_N1\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\thypoth_equivalog\tY\tNCBIFAM\tgeneN1\tDUF1234 domain-containing protein\n"
			+"14\t14\trep14\tNF_N2\t190\t100\tY\t180\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneN2\tDUF1234 domain-containing protein\n"
			+"15\t15\trep15\tNF_O1\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\thypoth_equivalog\tY\tNCBIFAM\tgeneO1\thypothetical protein\n"
			+"15\t15\trep15\tNF_O2\t190\t100\tY\t180\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneO2\thypothetical protein\n"
			+"16\t16\trep16\tNF_P1\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneP1\tDUF4321 domain-containing protein\n"
			+"16\t16\trep16\tNF_P2\t190\t100\tY\t180\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\thypoth_equivalog\tY\tNCBIFAM\tgeneP2\tDUF4321 domain-containing protein\n"
			+"17\t17\trep17\tNF_Q\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneQ\tproduct Q\n"
			+"18\t18\trep18\tNF_R\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneR\tproduct R\n"
			+"19\t19\trep19\tNF_S\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneS\tproduct S\n"
			+"20\t20\trep20\tNF_T1\t200\t100\tY\t190\t100\tY\t0.800000000000\t0.750000000000\tN\tN\thmm\tequivalog\tY\tNCBIFAM\tgeneT1\tproduct T\n"
			+"20\t20\trep20\tNF_T2\t210\t100\tY\t195\t100\tY\t1.000000000000\t0.799999999999\tN\tN\thmm\texception\tY\tNCBIFAM\tgeneT2\tproduct T\n"
			+"21\t21\trep21\tNF_U\t200\t100\tY\t190\t100\tY\t0.900000000000\t0.800000000000\tN\tN\thmm\tequivalog\tY\tNCBIFAM\tgeneU\tproduct U\n"
			+"22\t22\trep22\tNF_V\t200\t100\tY\t190\t100\tY\t0.900000000000\t0.749999999999\tN\tN\thmm\tequivalog\tY\tNCBIFAM\tgeneV\tproduct V\n"
			+"23\t23\trep23\tNF_W\t200\t100\tY\t190\t100\tY\t0.799999999999\t0.750000000000\tN\tN\thmm\tequivalog\tY\tNCBIFAM\tgeneW\tproduct W\n"
			+"24\t24\trep24\tNF_X\t200\t100\tY\t190\t100\tY\t0.900000000000\t0.750000000000\tN\tN\thmm\tsubfamily_domain\tY\tNCBIFAM\tgeneX\tproduct X\n"
			+"25\t25\trep25\tNF_Y\t200\t100\tY\t190\t100\tY\t0.900000000000\t0.750000000000\tN\tN\thmm\tequivalog\tN\tNCBIFAM\tgeneY\tproduct Y\n"
			+"26\t26\trep26\tNF_Z1\t200\t100\tY\t190\t100\tY\t0.900000000000\t0.750000000000\tN\tN\thmm\tequivalog\tY\tNCBIFAM\tgeneZ1\tproduct Z1\n"
			+"26\t26\trep26\tNF_Z2\t210\t100\tY\t195\t100\tY\t0.900000000000\t0.750000000000\tN\tN\thmm\texception\tY\tNCBIFAM\tgeneZ2\tproduct Z2\n"
			+"27\t27\trep27\tNF_AA\t200\t100\tY\t190\t100\tY\t0.900000000000\t0.750000000000\tN\tN\thmm\tequivalog\tY\tNCBIFAM\tgeneAA\tproduct AA\n"
			+"28\t28\trep28\tNF_AB\t200\t100\tY\t190\t100\tY\t0.900000000000\t0.750000000000\tN\tN\thmm\tequivalog\tY\tNCBIFAM\tgeneAB\tproduct AB\n"
			+"29\t29\trep29\tNF_AC\t200\t100\tY\t190\t100\tY\t0.800000000000\t0.800000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneAC\tproduct AC\n"
			+"30\t30\trep30\tNF_AD2\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneAD2\tproduct M2\n"
			+"30\t30\trep30\tNF_AD1\t190\t100\tY\t180\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneAD1\tproduct M1\n"
			+"30\t30\trep30\tNF_LOW\t500\t100\tY\t500\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tdomain\tY\tNCBIFAM\tgeneLOW\tproduct lower\n"
			+"32\t32\trep32\tNF_A_B\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneAB\ta|b\n"
			+"32\t32\trep32\tNF_C\t190\t100\tY\t180\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneC\tc\n"
			+"33\t33\trep33\tNF_A\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneA\ta\n"
			+"33\t33\trep33\tNF_B_C\t190\t100\tY\t180\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneBC\tb|c\n"
			+"34\t34\trep34\tNF608_BAD\t7.9\t100\tN\t7.9\t100\tN\t0.049300000000\t0.366700000000\tN\tN\thmm\tequivalog\tY\tNCBIFAM\tyncL\tstress response YncL\n"
			+"34\t34\trep34\tNF608_HIGH\t175.2\t100\tY\t175.2\t100\tY\t0.887900000000\t0.977000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgene608\thypothetical protein\n");
		write(dir, "coverage_domains.tsv", domainHeader()
			+"20\t20\trep20\t100\tNF_T\thmmT\t100\t1\t2\t1e-40\t200\t1e-40\t1e-40\t195\t1\t50\t1\t50\t1\t50\t0.500000000000\t0.500000000000\t100\t100\tequivalog\tNCBIFAM\tgeneT\tproduct T\tY\tboth\n"
			+"20\t20\trep20\t100\tNF_T\thmmT\t100\t2\t2\t1e-40\t200\t1e-35\t1e-35\t190\t1\t90\t1\t76\t1\t76\t0.900000000000\t0.760000000000\t100\t100\tequivalog\tNCBIFAM\tgeneT\tproduct T\tY\tquery_only\n"
			+"20\t20\trep20\t100\tNF_T2\thmmT2\t100\t1\t1\t1e-42\t210\t1e-37\t1e-37\t195\t1\t100\t1\t79\t1\t79\t1.000000000000\t0.790000000000\t100\t100\texception\tNCBIFAM\tgeneT2\tproduct T\tY\tquery_only\n"
			+"21\t21\trep21\t100\tNF_U\thmmU\t100\t1\t1\t1e-40\t200\t1e-40\t1e-40\t190\t1\t90\t1\t80\t1\t80\t0.900000000000\t0.800000000000\t100\t100\tequivalog\tNCBIFAM\tgeneU\tproduct U\tY\tcoverage_pass\n"
			+"22\t22\trep22\t100\tNF_V\thmmV\t100\t1\t1\t1e-40\t200\t1e-40\t1e-40\t190\t1\t90\t1\t74\t1\t74\t0.900000000000\t0.749999999999\t100\t100\tequivalog\tNCBIFAM\tgeneV\tproduct V\tY\tquery_only\n"
			+"23\t23\trep23\t100\tNF_W\thmmW\t100\t1\t1\t1e-40\t200\t1e-40\t1e-40\t190\t1\t79\t1\t75\t1\t75\t0.799999999999\t0.750000000000\t100\t100\tequivalog\tNCBIFAM\tgeneW\tproduct W\tY\tboth\n"
			+"24\t24\trep24\t100\tNF_X\thmmX\t100\t1\t1\t1e-40\t200\t1e-40\t1e-40\t190\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t100\t100\tsubfamily_domain\tNCBIFAM\tgeneX\tproduct X\tY\tquery_only\n"
			+"25\t25\trep25\t100\tNF_Y\thmmY\t100\t1\t1\t1e-40\t200\t1e-40\t1e-40\t190\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t100\t100\tequivalog\tNCBIFAM\tgeneY\tproduct Y\tN\tquery_only\n"
			+"26\t26\trep26\t100\tNF_Z1\thmmZ1\t100\t1\t1\t1e-40\t200\t1e-40\t1e-40\t190\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t100\t100\tequivalog\tNCBIFAM\tgeneZ1\tproduct Z1\tY\tquery_only\n"
			+"26\t26\trep26\t100\tNF_Z2\thmmZ2\t100\t1\t1\t1e-42\t210\t1e-42\t1e-42\t195\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t100\t100\texception\tNCBIFAM\tgeneZ2\tproduct Z2\tY\tquery_only\n"
			+"27\t27\trep27\t100\tNF_AA\thmmAA\t100\t1\t1\t1e-40\t200\t1e-40\t1e-40\t190\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t100\t100\tequivalog\tNCBIFAM\tgeneAA\tproduct AA\tY\tquery_only\n"
			+"28\t28\trep28\t100\tNF_AB\thmmAB\t100\t1\t1\t1e-40\t200\t1e-40\t1e-40\t190\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t100\t100\tequivalog\tNCBIFAM\tgeneAB\tproduct AB\tY\tquery_only\n");
		write(dir, "swiss.tsv", swissHeader()
			+"1\t1\trep1\tP1\tP1_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P1|P1_ENTRY product A OS=Species OX=1 GN=a PE=1 SV=1\n"
			+"2\t2\trep2\tP2\tP2_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P2|P2_ENTRY DNA gyrase subunit A OS=Species OX=1 GN=b PE=1 SV=1\n"
			+"3\t3\trep3\tP3\tP3_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P3|P3_ENTRY product C1 OS=Species OX=1 GN=c PE=1 SV=1\n"
			+"4\t4\trep4\tP4A\tP4A_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P4A|P4A_ENTRY product D1 OS=Species OX=1 GN=d PE=1 SV=1\n"
			+"4\t4\trep4\tP4B\tP4B_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t190\tsp|P4B|P4B_ENTRY product D2 OS=Species OX=1 GN=d PE=1 SV=1\n"
			+"5\t5\trep5\tP5\tP5_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P5|P5_ENTRY product E OS=Species OX=1 GN=e PE=1 SV=1\n"
			+"6\t6\trep6\tP6\tP6_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P6|P6_ENTRY dipeptide ABC transporter permease OS=Species OX=1 GN=f PE=1 SV=1\n"
			+"7\t7\trep7\tP7\tP7_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P7|P7_ENTRY product G OS=Species OX=1 GN=g PE=1 SV=1\n"
			+"8\t8\trep8\tP8\tP8_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P8|P8_ENTRY Product H OS=Species OX=1 GN=h PE=1 SV=1\n"
			+"9\t9\trep9\tP9\tP9_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P9|P9_ENTRY DNA gyrase subunit A OS=Species OX=1 GN=i PE=1 SV=1\n"
			+"10\t10\trep10\tP10\tP10_ENTRY\t92\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-40\t180\tsp|P10|P10_ENTRY product J OS=Species OX=1 GN=j PE=1 SV=1\n"
			+"11\t11\trep11\tP11\tP11_ENTRY\t88\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-30\t150\tsp|P11|P11_ENTRY subfamily offshoot OS=Species OX=1 GN=k PE=1 SV=1\n"
			+"17\t17\trep17\tP17B\tP17B_ENTRY\t98\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-60\t200\tsp|P17B|P17B_ENTRY product R OS=Species OX=1 GN=q PE=1 SV=1\n"
			+"17\t17\trep17\tP17C\tP17C_ENTRY\t97\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-55\t190\tsp|P17C|P17C_ENTRY product S OS=Species OX=1 GN=q PE=1 SV=1\n"
			+"17\t17\trep17\tP17A\tP17A_ENTRY\t96\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-20\t100\tsp|P17A|P17A_ENTRY product Q OS=Species OX=1 GN=q PE=1 SV=1\n"
			+"18\t18\trep18\tP18\tP18_ENTRY\t99\t100\t100\t100\t0.800000000000\t0.800000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P18|P18_ENTRY product R OS=Species OX=1 GN=r PE=1 SV=1\n"
			+"19\t19\trep19\tP19\tP19_ENTRY\t99\t100\t100\t100\t0.799999999999\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P19|P19_ENTRY product S OS=Species OX=1 GN=s PE=1 SV=1\n"
			+"20\t20\trep20\tP20\tP20_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P20|P20_ENTRY product T OS=Species OX=1 GN=t PE=1 SV=1\n"
			+"21\t21\trep21\tP21\tP21_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P21|P21_ENTRY product U OS=Species OX=1 GN=u PE=1 SV=1\n"
			+"22\t22\trep22\tP22\tP22_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P22|P22_ENTRY product V OS=Species OX=1 GN=v PE=1 SV=1\n"
			+"23\t23\trep23\tP23\tP23_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P23|P23_ENTRY product W OS=Species OX=1 GN=w PE=1 SV=1\n"
			+"24\t24\trep24\tP24\tP24_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P24|P24_ENTRY product X OS=Species OX=1 GN=x PE=1 SV=1\n"
			+"25\t25\trep25\tP25\tP25_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P25|P25_ENTRY product Y OS=Species OX=1 GN=y PE=1 SV=1\n"
			+"26\t26\trep26\tP26\tP26_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P26|P26_ENTRY product Z1 OS=Species OX=1 GN=z PE=1 SV=1\n"
			+"27\t27\trep27\tP27\tP27_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P27|P27_ENTRY product BB OS=Species OX=1 GN=aa PE=1 SV=1\n"
			+"28\t28\trep28\tP28A\tP28A_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P28A|P28A_ENTRY product AB OS=Species OX=1 GN=ab PE=1 SV=1\n"
			+"28\t28\trep28\tP28B\tP28B_ENTRY\t98\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-40\t180\tsp|P28B|P28B_ENTRY product AC OS=Species OX=1 GN=ab PE=1 SV=1\n"
			+"31\t31\trep31\tP31B\tP31B_ENTRY\t98\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-40\t180\tsp|P31B|P31B_ENTRY product AC OS=Species OX=1 GN=ab PE=1 SV=1\n"
			+"31\t31\trep31\tP31A\tP31A_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|P31A|P31A_ENTRY product AB OS=Species OX=1 GN=ab PE=1 SV=1\n"
			+"34\t34\trep34\tP608\tP608_UPF\t51.9\t100\t100\t100\t0.955000000000\t0.977000000000\t1\t100\t1\t100\t1e-45\t185\tsp|P608|P608_UPF UPF0126 family protein OS=Species OX=1 GN=upf PE=1 SV=1\n");
		write(dir, "reviewed.tsv",
			"family_id\tactive_index\tcandidate_class\tpilot_anchor\t"
			+"provisional_decision_class\tprovisional_product\tprimary_support\t"
			+"competing_support\ttransfer_rule_refinement\treviewer_notes\n"
			+"5\t5\tgated_full_length\tY\tconflicting_evidence\t.\t"
			+"primary\tcompeting\tmanual exception\tpreserve abstention\n"
			+"6\t6\tgated_full_length\tY\tdomain_supported_broad_name\t"
			+"peptide ABC transporter permease\tprimary\tcompeting\tmanual exception\tpreserve broad name\n");
		write(dir, "reviewed_evidence.tsv", reviewedEvidenceSidecar(reviewedEvidenceRow(6, 6,
			"rep6", "85.5")));
		write(dir, "synonyms.tsv",
			"from_product\tcanonical_product\tnote\n"
			+"DNA topoisomerase (ATP-hydrolyzing) subunit A\tDNA gyrase subunit A\treviewed synonym\n");
	}

	private static void checkReviewedEvidenceRejections(final Path dir) throws Exception{
		final String evidence=reviewedEvidenceSidecar(reviewedEvidenceRow(6, 6, "rep6",
			"85.5"));
		write(dir, "bad_reviewed_evidence_header.tsv",
			evidence.replace("source_methods\tsupport_tag",
				"support_tag\tsource_methods"));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_header.tsv",
			"bad_reviewed_evidence_header", "Unexpected reviewed-evidence column");
		write(dir, "bad_reviewed_evidence_meta.tsv",
			evidence.replace("#schema_version\tbbgap_reviewed_evidence_v1\n",
				"#schema_version\tbbgap_reviewed_evidence_v1\n"
				+"#schema_version\tbbgap_reviewed_evidence_v1\n"));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_meta.tsv",
			"bad_reviewed_evidence_meta", "Duplicate schema");
		write(dir, "bad_reviewed_evidence_pin.tsv",
			evidence.replace("#consensus_sha80\tconpin", "#consensus_sha80\tbad"));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_pin.tsv",
			"bad_reviewed_evidence_pin", "consensus pin");
		write(dir, "bad_reviewed_evidence_qcov.tsv",
			reviewedEvidenceSidecar(reviewedEvidenceRow(6, 6, "rep6", "101")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_qcov.tsv",
			"bad_reviewed_evidence_qcov", "Bad seed_query_cov");
		write(dir, "bad_reviewed_evidence_unreviewed.tsv",
			reviewedEvidenceSidecar(reviewedEvidenceRow(7, 7, "rep7", "85.5")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_unreviewed.tsv",
			"bad_reviewed_evidence_unreviewed", "lacks reviewed=");
		write(dir, "bad_reviewed_evidence_rep.tsv",
			reviewedEvidenceSidecar(reviewedEvidenceRow(6, 6, "rep_bad", "85.5")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_rep.tsv",
			"bad_reviewed_evidence_rep", "Representative mismatch");
		write(dir, "bad_reviewed_evidence_swiss_seed.tsv",
			reviewedEvidenceSidecar(swissReviewedEvidenceRow(6, 6, "rep6")
				.replace("scientific_review\t.\t.",
					"scientific_review\tfake_ortholog\t.")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_swiss_seed.tsv",
			"bad_reviewed_evidence_swiss_seed", "carries eggNOG field seed_ortholog");
		write(dir, "bad_reviewed_evidence_swiss_hit.tsv",
			reviewedEvidenceSidecar(swissReviewedEvidenceRow(6, 6, "rep6")
				.replace("swiss:P77333", "swiss:BAD")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_swiss_hit.tsv",
			"bad_reviewed_evidence_swiss_hit", "lacks retained hit");
		write(dir, "bad_reviewed_evidence_swiss_wrong_family.tsv",
			reviewedEvidenceSidecar(swissReviewedEvidenceRow(6, 6, "rep6")
				.replace("swiss:P77333", "swiss:P7")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_swiss_wrong_family.tsv",
			"bad_reviewed_evidence_swiss_wrong_family", "lacks retained hit");
		write(dir, "bad_reviewed_evidence_swiss_primary.tsv",
			reviewedEvidenceSidecar(validSwissReviewedEvidenceRow(6, 6, "rep6")
				.replace("Swiss-Prot P77333 identity=49.7 qcov=0.951298701299 "
					+"scov=0.983277591973 bitscore=276", ".")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_swiss_primary.tsv",
			"bad_reviewed_evidence_swiss_primary", "lacks primary_support");
		write(dir, "bad_reviewed_evidence_swiss_competing.tsv",
			reviewedEvidenceSidecar(validSwissReviewedEvidenceRow(6, 6, "rep6")
				.replace("Q01610 and P72131 remain LysR-family alternatives", ". . .")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_swiss_competing.tsv",
			"bad_reviewed_evidence_swiss_competing", "lacks competing_support");
		write(dir, "bad_reviewed_evidence_swiss_notes.tsv",
			reviewedEvidenceSidecar(validSwissReviewedEvidenceRow(6, 6, "rep6")
				.replace("PF00126/PF03466/PS50931 support LysR-family homology only",
					"   ")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_swiss_notes.tsv",
			"bad_reviewed_evidence_swiss_notes", "lacks reviewer_notes");
		write(dir, "bad_reviewed_evidence_ncbi_seed.tsv",
			reviewedEvidenceSidecar(nativeReviewedEvidenceRow(6, 6, "rep6")
				.replace("scientific_review\t.\t.",
					"scientific_review\tfake_ortholog\t.")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_ncbi_seed.tsv",
			"bad_reviewed_evidence_ncbi_seed", "carries eggNOG field seed_ortholog");
		write(dir, "bad_reviewed_evidence_ncbi_hit.tsv",
			reviewedEvidenceSidecar(nativeReviewedEvidenceRow(6, 6, "rep6")
				.replace("ncbi:0", "ncbi:1")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_ncbi_hit.tsv",
			"bad_reviewed_evidence_ncbi_hit", "lacks retained hit");
		write(dir, "bad_reviewed_evidence_dotted_eggnog.tsv",
			reviewedEvidenceSidecar(reviewedEvidenceRow(6, 6, "rep6", "85.5")
				.replace("\t38557254\t73.2\t85.5\t99.0\t1e-42",
					"\t.\t73.2\t85.5\t99.0\t1e-42")));
		checkReviewedEvidenceRejection(dir, "bad_reviewed_evidence_dotted_eggnog.tsv",
			"bad_reviewed_evidence_dotted_eggnog", "Bad selected_seed_numeric");
	}

	private static void writeEggnogInputs(final Path dir) throws Exception{
		write(dir, "eggnog_family_identity.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\n"
			+"1\t1\t0\terep1\t10\t100\n"
			+"2\t2\t0\terep2\t10\t100\n"
			+"3\t3\t0\terep3\t10\t100\n"
			+"4\t4\t0\terep4\t10\t100\n"
			+"5\t5\t0\terep5\t10\t100\n"
			+"6\t6\t0\terep6\t10\t100\n"
			+"7\t7\t0\terep7\t10\t100\n"
			+"8\t8\t0\terep8\t10\t100\n"
			+"9\t9\t0\terep9\t10\t100\n"
			+"10\t10\t0\terep10\t10\t100\n"
			+"11\t11\t0\terep11\t10\t100\n"
			+"12\t12\t0\terep12\t10\t100\n"
			+"13\t13\t0\terep13\t10\t100\n"
			+"14\t14\t0\terep14\t10\t100\n");
		write(dir, "eggnog_ncbi.tsv", modelHeader()
			+"1\t1\terep1\tNF_EA\t200\t100\tY\t190\t100\tY\t1.000000000000\t1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneA\tproduct A\n"
			+"3\t3\terep3\tNF_EC\t200\t100\tN\t190\t100\tN\t1.000000000000\t1.000000000000\tY\tN\thmm\tequivalog\tN\tNCBIFAM\tgeneC\tproduct C\n");
		write(dir, "eggnog_swiss.tsv", swissHeader()
			+"1\t1\terep1\tEP1\tEP1_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tsp|EP1|EP1_ENTRY product A OS=Species OX=1 GN=a PE=1 SV=1\n");
		write(dir, "eggnog_reviewed.tsv",
			"family_id\tactive_index\tcandidate_class\tpilot_anchor\t"
			+"provisional_decision_class\tprovisional_product\tprimary_support\t"
			+"competing_support\ttransfer_rule_refinement\treviewer_notes\n"
			+"2\t2\tgated_full_length\tY\tdomain_supported_broad_name\treviewed protein\t"
			+"primary\tcompeting\tmanual exception\tpreserve reviewed\n");
		write(dir, "eggnog_homology.tsv", eggnogHomologySidecar(
			eggnogHomologyRow(1, "DUF1815", "DUF1815_2_101", "DUF1815@2|A-1",
				"-", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80", "29609749")
			+eggnogHomologyRow(2, "DUF1815", "DUF1815_2_101", "DUF1815@2|A-1",
				"-", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80", "29609749")
			+eggnogHomologyRow(3, "DUF1815", "DUF1815_2_101", "DUF1815@2|A-1",
				"-", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80", "29609749")
			+eggnogHomologyRow(4, "DUF1815", "DUF1815_2_101", "DUF1815@2|A-1",
				"-", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80", "29609749")
			+eggnogHomologyRow(5, "12TM_1", "12TM_1_40_476",
				"12TM_1@2|A-1,12TM_1@131567|A-2", "-", "-", "-", "-", "f",
				"seed_pfam_width_fraction_ge_0.80", "24219188")
			+eggnogHomologyRow(6, "DUF1,DUF2", "DUF1_1_10,DUF2_1_10",
				"DUF1@2|A-1,DUF2@2|A-1", "-", "-", "-", "-", "f",
				"seed_pfam_width_fraction_0.50_0.80", "6")
			+eggnogHomologyRow(7, "DUF3", "DUF3_1_10", "DUF4@2|A-1",
				"-", "-", "-", "-", "f", "seed_pfam_width_lt_80", "7")
			+eggnogHomologyRow(8, "DUF5", ".", "UNK",
				"-", "-", "-", "-", "f", "not_applicable", "8")
			+eggnogHomologyRow(9, "DUF6", "DUF6_1_10", "DUF6@2|A-1",
				"enzyme β-chain note", "K00001", "1.1.1.1", "alcohol dehydrogenase",
				"f", "seed_pfam_width_fraction_lt_0.50", "9")
			+eggnogHomologyRow(10, "DUF7", "DUF7_1_10", "DUF7@2|A-1",
				"-", "K00001,K00002", "-", "K00001=uncharacterized protein,K00002=kinase",
				"f", "seed_pfam_width_fraction_ge_0.80", "10")
			+eggnogHomologyRow(11, "DUF8", "DUF8_1_10", "DUF8@2|A-1",
				"-", "K00003", "-", "K00003=uncharacterized methyltransferase",
				"f", "seed_pfam_width_fraction_ge_0.80", "11")
			+eggnogHomologyRow(12, "DUF9", "DUF9_1_10", "DUF9@2|A-1",
				"orf", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80", "12")
			+eggnogHomologyRow(13, "DUF10", "DUF10_1_10", "DUF10@2|A-1",
				"LOC123456", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80",
				"13")
			+eggnogHomologyRow(14, "DUF11", "DUF11_1_10", "DUF11@2|A-1",
				"1d0326", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80",
				"14")));
	}

	private static void checkEggnogHomologyRejections(final Path dir) throws Exception{
		checkEggnogHomologyRejection(dir, "eggnog_requires_provisional",
			eggnogHomologySidecar(eggnogHomologyRow(4, "DUF1815", "DUF1815_2_101",
				"DUF1815@2|A-1", "-", "-", "-", "-", "f",
				"seed_pfam_width_fraction_ge_0.80", "29609749")),
			"eggnoghomology requires provisional", false);
		checkEggnogHomologyRejection(dir, "eggnog_duplicate",
			eggnogHomologySidecar(
				eggnogHomologyRow(4, "DUF1815", "DUF1815_2_101", "DUF1815@2|A-1",
					"-", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80",
					"29609749")
				+eggnogHomologyRow(4, "DUF1815", "DUF1815_2_101", "DUF1815@2|A-1",
					"-", "-", "-", "-", "f", "seed_pfam_width_fraction_ge_0.80",
					"29609749")),
			"Duplicate eggNOG-homology", true);
		checkEggnogHomologyRejection(dir, "eggnog_bad_suffix",
			eggnogHomologySidecar(eggnogHomologyRow(4, "DUF1815", "DUF1815_x_101",
				"DUF1815@2|A-1", "-", "-", "-", "-", "f",
				"seed_pfam_width_fraction_ge_0.80", "29609749")),
			"Bad pfam_raw_start", true);
	}

	private static void checkEggnogHomologyRejection(final Path dir, final String tag,
			final String evidence, final String message, final boolean provisional)
			throws Exception{
		write(dir, tag+".tsv", evidence);
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("eggnog_family_identity.tsv"),
				"ncbi="+dir.resolve("eggnog_ncbi.tsv"),
				"swiss="+dir.resolve("eggnog_swiss.tsv"),
				"reviewed="+dir.resolve("eggnog_reviewed.tsv"),
				"out="+dir.resolve(tag+"_selected.tsv"),
				"evidence="+dir.resolve(tag+"_evidence.tsv"),
				"unresolved="+dir.resolve(tag+"_unresolved.tsv"),
				"summary="+dir.resolve(tag+"_summary.tsv"),
				"eggnoghomology="+dir.resolve(tag+".tsv"),
				"eggnoghomologyout="+dir.resolve(tag+"_validated.tsv"),
				"labelversion="+tag,
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy",
				"provisional="+provisional});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains(message);
		}
		check(ok, "eggNOG-homology rejection failed for "+tag);
	}

	private static String eggnogHomologySidecar(final String rows){
		return "#schema_version\tbbgap_eggnog_homology_v2\n"
			+"#family_identity_sha80\tidpin\n"
			+"#consensus_sha80\tconpin\n"
			+"family_id\tactive_index\trep_id\tseed_ortholog\tselected_seed_numeric\t"
			+"seed_identity\tseed_query_cov\tseed_subject_cov\tseed_evalue\tseed_score\t"
			+"pfam_root\tpfam_raw\tseed_db_pfam\tannotation_pfam\teggNOG_OGs\t"
			+"Preferred_name\tKOs\tEC\tKO_definition\tfunctional_lead\ttax_ceiling\t"
			+"farthest_donor_lineage\tannotation_confidence\tspan_bin\n"
			+rows;
	}

	private static String eggnogHomologyRow(final int familyId, final String pfamRoot,
			final String pfamRaw, final String ogs, final String preferredName,
			final String kos, final String ec, final String koDefinition,
			final String functionalLead, final String spanBin, final String seedNumeric){
		return eggnogHomologyRow(familyId, pfamRoot, pfamRaw, pfamRaw, pfamRaw, ogs,
			preferredName, kos, ec, koDefinition, functionalLead, spanBin, seedNumeric);
	}

	private static String eggnogHomologyRow(final int familyId, final String pfamRoot,
			final String pfamRaw, final String seedDbPfam, final String annotationPfam,
			final String ogs, final String preferredName, final String kos, final String ec,
			final String koDefinition, final String functionalLead, final String spanBin,
			final String seedNumeric){
		return familyId+"\t"+familyId+"\terep"+familyId+"\tseed_"+familyId+"\t"
			+seedNumeric+"\t85.5\t91.0\t92.0\t1e-20\t99.5\t"+pfamRoot+"\t"
			+pfamRaw+"\t"+seedDbPfam+"\t"+annotationPfam+"\t"+ogs+"\t"
			+preferredName+"\t"+kos+"\t"+ec+"\t"+koDefinition+"\t"+functionalLead
			+"\t2\tBacteria\tlow\t"+spanBin+"\n";
	}

	private static void checkReviewedEvidenceRejection(final Path dir, final String sidecar,
			final String tag, final String message){
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"reviewed="+dir.resolve("reviewed.tsv"),
				"reviewedevidence="+dir.resolve(sidecar),
				"reviewedevidenceout="+dir.resolve(tag+"_validated.tsv"),
				"out="+dir.resolve(tag+"_selected.tsv"),
				"evidence="+dir.resolve(tag+"_evidence.tsv"),
				"unresolved="+dir.resolve(tag+"_unresolved.tsv"),
				"summary="+dir.resolve(tag+"_summary.tsv"),
				"labelversion="+tag,
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains(message);
		}
		check(ok, "reviewed-evidence rejection failed for "+tag);
	}

	private static String reviewedEvidenceSidecar(final String rows){
		return "#schema_version\tbbgap_reviewed_evidence_v1\n"
			+"#family_identity_sha80\tidpin\n"
			+"#consensus_sha80\tconpin\n"
			+"family_id\tactive_index\trep_id\thit_id\tsource_methods\tsupport_tag\t"
			+"coverage_tag\tscore_tag\tproduct_text_tag\treview_tag\tseed_ortholog\t"
			+"selected_seed_numeric\tseed_identity\tseed_query_cov\tseed_subject_cov\t"
			+"seed_evalue\tseed_score\teggNOG_OGs\ttax_ceiling\t"
			+"farthest_donor_lineage\traw_Preferred_name\traw_KEGG_ko\traw_EC\t"
			+"raw_PFAMs\tannotation_confidence\tprimary_support\tcompeting_support\t"
			+"warnings\treviewer_notes\n"
			+rows;
	}

	private static String reviewedEvidenceRow(final int familyId, final int activeIndex,
			final String repId, final String qcov){
		return familyId+"\t"+activeIndex+"\t"+repId+"\teggnog:38557254;DUF3097\t"
			+"eggNOG/DIAMOND\treviewed_eggnog_homology_domain\tcoverage_full\t"
			+"selected_seed_passed\teggnog_domain_family\tmanual_noevidence_reviewβ\t"
			+"38557254.LLCGMKFL_04186\t38557254\t73.2\t"+qcov+"\t99.0\t"
			+"1e-42\t254.5\tCOG5432@1|root\t2759\tBacteria\tDUF3097\tK02004\t\t"
			+"PF03733\tmedium\tOG/Pfam support\tno material competitor\t"
			+"non-ASCII warning β\tpunctuation; hit id retained\n";
	}

	private static String swissReviewedEvidenceRow(final int familyId, final int activeIndex,
			final String repId){
		return familyId+"\t"+activeIndex+"\t"+repId+"\tswiss:P77333\t"
			+"Swiss-Prot/DIAMOND\tretained_full_swissprot\tcoverage_full\t"
			+"retained_bitscore\tfull_coverage_family\tscientific_review\t.\t.\t.\t.\t"
			+".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t"
			+"Swiss-Prot P77333 identity=49.7 qcov=0.951298701299 "
			+"scov=0.983277591973 bitscore=276\t"
			+"Q01610 and P72131 remain LysR-family alternatives\t.\t"
			+"PF00126/PF03466/PS50931 support LysR-family homology only\n";
	}

	private static String validSwissReviewedEvidenceRow(final int familyId,
			final int activeIndex, final String repId){
		return swissReviewedEvidenceRow(familyId, activeIndex, repId)
			.replace("swiss:P77333", "swiss:P6");
	}

	private static String nativeReviewedEvidenceRow(final int familyId, final int activeIndex,
			final String repId){
		return familyId+"\t"+activeIndex+"\t"+repId+"\tncbi:0\t"
			+"NCBIfam/HMMER\tretained_full_native\tcoverage_full\tnative_bitscore\t"
			+"native_family\tscientific_review\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t."
			+"\t.\t.\t.\tNCBIfam NF_G retained native hit\t"
			+"No material competitor\t.\tNative homology only\n";
	}

	private static void checkBadSwissTitle(final Path dir) throws Exception{
		write(dir, "bad_swiss.tsv", swissHeader()
			+"1\t1\trep1\tP1\tP1_ENTRY\t99\t100\t100\t100\t1.000000000000\t1.000000000000\t1\t100\t1\t100\t1e-50\t200\tproduct A\n");
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("bad_swiss.tsv"),
				"out="+dir.resolve("bad_selected.tsv"),
				"evidence="+dir.resolve("bad_evidence.tsv"),
				"unresolved="+dir.resolve("bad_unresolved.tsv"),
				"summary="+dir.resolve("bad_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Bad Swiss-Prot title");
		}
		check(ok, "malformed Swiss-Prot title should crash loudly");
	}

	private static void checkBadSwissMetric(final Path dir) throws Exception{
		write(dir, "bad_swiss_metric.tsv", swissHeader()
			+"1\t1\trep1\tP1\tP1_ENTRY\tNaN\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P1|P1_ENTRY product A OS=Species OX=1 GN=a PE=1 SV=1\n");
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("bad_swiss_metric.tsv"),
				"out="+dir.resolve("bad_swiss_metric_selected.tsv"),
				"evidence="+dir.resolve("bad_swiss_metric_evidence.tsv"),
				"unresolved="+dir.resolve("bad_swiss_metric_unresolved.tsv"),
				"summary="+dir.resolve("bad_swiss_metric_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Bad Swiss-Prot identity");
		}
		check(ok, "non-finite Swiss-Prot identity should crash loudly");
	}

	private static void checkBadNcbiNearMetric(final Path dir) throws Exception{
		write(dir, "bad_ncbi_near.tsv", domainHeader()
			+"1\t1\trep1\t100\tNF_BAD\thmmBad\t100\t1\t1\t1e-40\tNaN\t1e-40\t"
			+"1e-40\t190\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t"
			+"100\t100\tequivalog\tNCBIFAM\tgeneBad\tproduct Bad\tY\tquery_only\n");
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"ncbinearquerydomains="+dir.resolve("bad_ncbi_near.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("bad_ncbi_near_selected.tsv"),
				"evidence="+dir.resolve("bad_ncbi_near_evidence.tsv"),
				"unresolved="+dir.resolve("bad_ncbi_near_unresolved.tsv"),
				"summary="+dir.resolve("bad_ncbi_near_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy",
				"ncbinearquery=t"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Bad NCBI full_score");
		}
		check(ok, "non-finite NCBI near-query metrics should crash loudly");
	}

	private static void checkBadNcbiNearHeader(final Path dir) throws Exception{
		write(dir, "bad_ncbi_near_header.tsv", swappedDomainHeader()
			+"1\t1\trep1\t100\tNF_BAD\thmmBad\t100\t1\t1\t1e-40\t200\t1e-40\t"
			+"1e-40\t190\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t"
			+"100\t100\tequivalog\tNCBIFAM\tgeneBad\tproduct Bad\tY\tquery_only\n");
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"ncbinearquerydomains="+dir.resolve("bad_ncbi_near_header.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("bad_ncbi_near_header_selected.tsv"),
				"evidence="+dir.resolve("bad_ncbi_near_header_evidence.tsv"),
				"unresolved="+dir.resolve("bad_ncbi_near_header_unresolved.tsv"),
				"summary="+dir.resolve("bad_ncbi_near_header_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy",
				"ncbinearquery=t"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Unexpected NCBI near-query domain header");
		}
		check(ok, "swapped NCBI near-query domain headers should crash loudly");
	}

	private static void checkBadNcbiNearFamilyIdentity(final Path dir) throws Exception{
		checkBadNcbiNearFamilyIdentity(dir, "active",
			"1\t2\trep1\t100\tNF_BAD\thmmBad\t100\t1\t1\t1e-40\t200\t1e-40\t"
			+"1e-40\t190\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t"
			+"100\t100\tequivalog\tNCBIFAM\tgeneBad\tproduct Bad\tY\tquery_only\n",
			"Active-index mismatch");
		checkBadNcbiNearFamilyIdentity(dir, "rep",
			"1\t1\trep2\t100\tNF_BAD\thmmBad\t100\t1\t1\t1e-40\t200\t1e-40\t"
			+"1e-40\t190\t1\t90\t1\t75\t1\t75\t0.900000000000\t0.750000000000\t"
			+"100\t100\tequivalog\tNCBIFAM\tgeneBad\tproduct Bad\tY\tquery_only\n",
			"Representative mismatch");
	}

	private static void checkBadNcbiNearFamilyIdentity(final Path dir, final String tag,
			final String row, final String message) throws Exception{
		write(dir, "bad_ncbi_near_"+tag+".tsv", domainHeader()+row);
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"ncbinearquerydomains="+dir.resolve("bad_ncbi_near_"+tag+".tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("bad_ncbi_near_"+tag+"_selected.tsv"),
				"evidence="+dir.resolve("bad_ncbi_near_"+tag+"_evidence.tsv"),
				"unresolved="+dir.resolve("bad_ncbi_near_"+tag+"_unresolved.tsv"),
				"summary="+dir.resolve("bad_ncbi_near_"+tag+"_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy",
				"ncbinearquery=t"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains(message);
		}
		check(ok, "bad NCBI near-query family identity should crash loudly for "+tag);
	}

	private static void checkNcbiNearRequiresHierarchy(final Path dir) throws Exception{
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("near_exact_selected.tsv"),
				"evidence="+dir.resolve("near_exact_evidence.tsv"),
				"unresolved="+dir.resolve("near_exact_unresolved.tsv"),
				"summary="+dir.resolve("near_exact_summary.tsv"),
				"labelversion=test_near_exact",
				"identitypin=idpin",
				"consensuspin=conpin",
				"ncbinearquery=t"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("ncbinearquery requires mode=hierarchy");
		}
		check(ok, "near-query fallback should require hierarchy mode");
	}

	private static void checkZeroSwissBits(final Path dir) throws Exception{
		write(dir, "zero_family.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\n"
			+"1\t1\t0\trep1\t10\t100\n");
		write(dir, "zero_ncbi.tsv", modelHeader()
			+"1\t1\trep1\tNF_X\t200\t100\tY\t190\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneX\tproduct X\n");
		write(dir, "zero_swiss.tsv", swissHeader()
			+"1\t1\trep1\tP1\tP1_ENTRY\t99\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t100\t"
			+"sp|P1|P1_ENTRY product X OS=Species OX=1 GN=x PE=1 SV=1\n"
			+"1\t1\trep1\tP2\tP2_ENTRY\t98\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1\t0\t"
			+"sp|P2|P2_ENTRY product Y OS=Species OX=1 GN=y PE=1 SV=1\n");
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("zero_family.tsv"),
			"ncbi="+dir.resolve("zero_ncbi.tsv"),
			"swiss="+dir.resolve("zero_swiss.tsv"),
			"out="+dir.resolve("zero_selected.tsv"),
			"evidence="+dir.resolve("zero_evidence.tsv"),
			"unresolved="+dir.resolve("zero_unresolved.tsv"),
			"summary="+dir.resolve("zero_summary.tsv"),
			"labelversion=test_v1",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy"});
		check(read(dir, "zero_evidence.tsv").contains("\n1\t1\trep1\tspecific_full_length\t"
			+"product X\tncbi_specific_rank70_native_name\t1\t2\tNF_X\tP1,P2\t"
			+"product X\tproduct X|product Y\tmixed\tP1\tproduct X\t99\t"
			+"1.000000000000\t1.000000000000\t100\tP2\tproduct Y\t98\t"
			+"1.000000000000\t1.000000000000\t0\t100.000\t\tP1\tproduct X\t"
			+"99\t1.000000000000\t1.000000000000\t100\tP2\tproduct Y\t98\t"
			+"1.000000000000\t1.000000000000\t0\t100.000\t\n"),
			"zero Swiss-Prot bits should retain deltas and emit undefined ratios");
	}

	private static void checkSwissOnlyFallback(final Path dir) throws Exception{
		write(dir, "swiss_only_family.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\n"
			+"1\t1\t0\trep1\t10\t100\n"
			+"2\t2\t0\trep2\t10\t100\n"
			+"3\t3\t0\trep3\t10\t100\n"
			+"4\t4\t0\trep4\t10\t100\n"
			+"5\t5\t0\trep5\t10\t100\n"
			+"6\t6\t0\trep6\t10\t100\n");
		write(dir, "swiss_only_ncbi.tsv", modelHeader());
		write(dir, "swiss_only_swiss.tsv", swissHeader()
			+"1\t1\trep1\tP1\tP1_ENTRY\t95\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P1|P1_ENTRY product A OS=Species OX=1 GN=a PE=1 SV=1\n"
			+"2\t2\trep2\tP2\tP2_ENTRY\t89\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-40\t180\t"
			+"sp|P2|P2_ENTRY product B OS=Species OX=1 GN=b PE=1 SV=1\n"
			+"3\t3\trep3\tP3A\tP3A_ENTRY\t96\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-60\t220\t"
			+"sp|P3A|P3A_ENTRY product C OS=Species OX=1 GN=c PE=1 SV=1\n"
			+"3\t3\trep3\tP3B\tP3B_ENTRY\t97\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-70\t240\t"
			+"sp|P3B|P3B_ENTRY product D OS=Species OX=1 GN=d PE=1 SV=1\n"
			+"4\t4\trep4\tP4\tP4_ENTRY\t99\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-80\t260\t"
			+"sp|P4|P4_ENTRY Uncharacterized glycosyltransferase MK0977 "
			+"OS=Species OX=1 GN=e PE=1 SV=1\n"
			+"5\t5\trep5\tP5A\tP5A_ENTRY\t89\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-90\t300\t"
			+"sp|P5A|P5A_ENTRY product E OS=Species OX=1 GN=f PE=1 SV=1\n"
			+"5\t5\trep5\tP5B\tP5B_ENTRY\t95\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-30\t100\t"
			+"sp|P5B|P5B_ENTRY product E OS=Species OX=1 GN=f PE=1 SV=1\n"
			+"6\t6\trep6\tP6A\tP6A_ENTRY\t95\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P6A|P6A_ENTRY product F OS=Species OX=1 GN=g PE=1 SV=1\n"
			+"6\t6\trep6\tP6B\tP6B_ENTRY\t99\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-80\t260\t"
			+"sp|P6B|P6B_ENTRY Uncharacterized protein YaaB OS=Species OX=1 GN=h PE=1 SV=1\n");
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("swiss_only_family.tsv"),
			"ncbi="+dir.resolve("swiss_only_ncbi.tsv"),
			"swiss="+dir.resolve("swiss_only_swiss.tsv"),
			"out="+dir.resolve("swiss_only_selected.tsv"),
			"evidence="+dir.resolve("swiss_only_evidence.tsv"),
			"unresolved="+dir.resolve("swiss_only_unresolved.tsv"),
			"summary="+dir.resolve("swiss_only_summary.tsv"),
			"labelversion=test_swiss_only",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy",
			"swissonly=t",
			"minswissidentity=90"});
		final String selected=read(dir, "swiss_only_selected.tsv");
		final String evidence=read(dir, "swiss_only_evidence.tsv");
		final String summary=read(dir, "swiss_only_summary.tsv");
		check(selected.contains("\n1\t1\trep1\tswissprot_full_length\tproduct A\t"
			+"test_swiss_only\n"),
			"Swiss-only singleton above the identity threshold should be selected");
		check(!selected.contains("\n2\t") && !selected.contains("\n3\t") &&
			!selected.contains("\n4\t") && !selected.contains("\n5\t") &&
			!selected.contains("\n6\t"),
			"Swiss-only fallback should reject low-identity, conflicting and placeholder rows");
		check(evidence.contains("\n1\t1\trep1\tswissprot_full_length\tproduct A\t"
			+"swissprot_singleton_identity_90.0\t0\t1"),
			"Swiss-only fallback should record its rule");
		check(evidence.contains("\n5\t5\trep5\tunnamed\t\t"
			+"insufficient_native_ncbi_evidence\t0\t1\t\tP5A,P5B")
			&& evidence.contains("\n6\t6\trep6\tconflicting_evidence\t\t"
				+"competing_swissprot_products\t0\t2\t\tP6A,P6B"),
			"Swiss-only fallback should threshold the best same-name hit and retain placeholders");
		check(summary.contains("selected_auto_swissonly\t1\n")
			&& summary.contains("selected_auto\t1\n")
			&& summary.contains("unresolved_conflict\t2\n")
			&& summary.contains("unresolved_insufficient\t3\n")
			&& summary.contains("swiss_only_fallback\ttrue\n")
			&& summary.contains("min_swiss_identity\t90.0\n"),
			"Swiss-only fallback summary counts changed");
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("swiss_only_family.tsv"),
			"ncbi="+dir.resolve("swiss_only_ncbi.tsv"),
			"swiss="+dir.resolve("swiss_only_swiss.tsv"),
			"out="+dir.resolve("swiss_only_default_selected.tsv"),
			"evidence="+dir.resolve("swiss_only_default_evidence.tsv"),
			"unresolved="+dir.resolve("swiss_only_default_unresolved.tsv"),
			"summary="+dir.resolve("swiss_only_default_summary.tsv"),
			"labelversion=test_swiss_only",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy"});
		check(!read(dir, "swiss_only_default_selected.tsv").contains("\n1\t"),
			"Swiss-only fallback should be off by default");
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("swiss_only_family.tsv"),
				"ncbi="+dir.resolve("swiss_only_ncbi.tsv"),
				"swiss="+dir.resolve("swiss_only_swiss.tsv"),
				"out="+dir.resolve("swiss_only_exact_selected.tsv"),
				"evidence="+dir.resolve("swiss_only_exact_evidence.tsv"),
				"unresolved="+dir.resolve("swiss_only_exact_unresolved.tsv"),
				"summary="+dir.resolve("swiss_only_exact_summary.tsv"),
				"labelversion=test_swiss_only",
				"identitypin=idpin",
				"consensuspin=conpin",
				"swissonly=t"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("swissonly requires mode=hierarchy");
		}
		check(ok, "Swiss-only fallback should require hierarchy mode");
	}

	private static void checkDuplicateOutput(final Path dir) throws Exception{
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("dup.tsv"),
				"evidence="+dir.resolve("dup.tsv"),
				"unresolved="+dir.resolve("dup_unresolved.tsv"),
				"summary="+dir.resolve("dup_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(Throwable t){
			ok=true;
		}
		check(ok, "duplicate output paths should be rejected before writing");

		final Path a=dir.resolve("hardlink_member_coverage.tsv");
		final Path b=dir.resolve("hardlink_member_unresolved.tsv");
		Files.deleteIfExists(a);
		Files.deleteIfExists(b);
		Files.write(a, "preimage\n".getBytes(StandardCharsets.UTF_8));
		Files.createLink(b, a);
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("hardlink_selected.tsv"),
				"evidence="+dir.resolve("hardlink_evidence.tsv"),
				"unresolved="+dir.resolve("hardlink_unresolved.tsv"),
				"summary="+dir.resolve("hardlink_summary.tsv"),
				"membercoverage="+a,
				"memberunresolved="+b,
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin",
				"ow=t"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Output aliases output");
		}
		check(ok && read(dir, "hardlink_member_coverage.tsv").equals("preimage\n"),
			"hard-linked outputs should be rejected before overwriting their preimage");

		final Path c=dir.resolve("hardlink_reviewed_evidence.tsv");
		final Path d=dir.resolve("hardlink_reviewed_evidence_validated.tsv");
		Files.deleteIfExists(c);
		Files.deleteIfExists(d);
		Files.write(c, reviewedEvidenceSidecar(reviewedEvidenceRow(6, 6, "rep6",
			"85.5")).getBytes(StandardCharsets.UTF_8));
		Files.createLink(d, c);
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"reviewed="+dir.resolve("reviewed.tsv"),
				"reviewedevidence="+c,
				"reviewedevidenceout="+d,
				"out="+dir.resolve("hardlink_reviewed_evidence_selected.tsv"),
				"evidence="+dir.resolve("hardlink_reviewed_evidence_evidence.tsv"),
				"unresolved="+dir.resolve("hardlink_reviewed_evidence_unresolved.tsv"),
				"summary="+dir.resolve("hardlink_reviewed_evidence_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin",
				"ow=t"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Output aliases input");
		}
		check(ok && read(dir, "hardlink_reviewed_evidence.tsv").contains(
			"manual_noevidence_reviewβ"),
			"reviewed-evidence input/output aliases should be rejected");
	}

	private static void checkBroadConflictRules(final Path dir) throws Exception{
		write(dir, "broad_family.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\n"
			+"1\t1\t0\trep1\t10\t100\n"
			+"2\t2\t0\trep2\t10\t100\n"
			+"3\t3\t0\trep3\t10\t100\n"
			+"4\t4\t0\trep4\t10\t100\n"
			+"5\t5\t0\trep5\t10\t100\n"
			+"6\t6\t0\trep6\t10\t100\n"
			+"7\t7\t0\trep7\t10\t100\n"
			+"8\t8\t0\trep8\t10\t100\n"
			+"9\t9\t0\trep9\t10\t100\n"
			+"10\t10\t0\trep10\t10\t100\n"
			+"11\t11\t0\trep11\t10\t100\n"
			+"12\t12\t0\trep12\t10\t100\n");
		write(dir, "broad_ncbi.tsv", modelHeader()
			+"1\t1\trep1\tNF_RND1\t200\t100\tY\t190\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tsubfamily\tY\tNCBIFAM\t\t"
			+"efflux RND transporter permease subunit\n"
			+"1\t1\trep1\tNF_RND2\t190\t100\tY\t180\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tsubfamily\tY\tNCBIFAM\t\t"
			+"multidrug efflux RND transporter permease subunit\n"
			+"3\t3\trep3\tNF_RND3A\t200\t100\tY\t190\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tsubfamily\tY\tNCBIFAM\t\t"
			+"efflux RND transporter permease subunit\n"
			+"3\t3\trep3\tNF_RND3B\t190\t100\tY\t180\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tsubfamily\tY\tNCBIFAM\t\t"
			+"multidrug efflux RND transporter permease subunit\n"
			+"3\t3\trep3\tNF_RND3C\t180\t100\tY\t170\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tsubfamily\tY\tNCBIFAM\t\t"
			+"heavy metal efflux RND transporter permease subunit\n"
			+"6\t6\trep6\tNF_KEEP\t200\t100\tY\t190\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneKeep\t"
			+"product Keep\n"
			+"8\t8\trep8\tNF_RND8\t200\t100\tY\t190\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tsubfamily\tY\tNCBIFAM\t\t"
			+"efflux RND transporter permease subunit\n"
			+"8\t8\trep8\tNF_ABC8\t190\t100\tY\t180\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tsubfamily\tY\tNCBIFAM\t\t"
			+"peptide ABC transporter permease\n"
			+"9\t9\trep9\tNF_RND9\t200\t100\tY\t190\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tsubfamily\tY\tNCBIFAM\t\t"
			+"efflux RND transporter permease subunit\n"
			+"9\t9\trep9\tNF_MISSING9\t190\t100\tY\t180\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tsubfamily\tY\tNCBIFAM\t\t"
			+"copper efflux RND transporter permease subunit\n");
		write(dir, "broad_swiss.tsv", swissHeader()
			+"1\t1\trep1\tP1\tP1_ENTRY\t95\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P1|P1_ENTRY Multidrug efflux pump subunit AcrB OS=Species OX=1 GN=a PE=1 SV=1\n"
			+"2\t2\trep2\tP2A\tP2A_ENTRY\t45\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P2A|P2A_ENTRY HTH-type transcriptional regulator CrgA OS=Species OX=1 GN=b PE=1 SV=1\n"
			+"2\t2\trep2\tP2B\tP2B_ENTRY\t44\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-40\t180\t"
			+"sp|P2B|P2B_ENTRY HTH-type transcriptional activator AaeR OS=Species OX=1 GN=b PE=1 SV=1\n"
			+"4\t4\trep4\tP4A\tP4A_ENTRY\t45\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P4A|P4A_ENTRY HTH-type transcriptional regulator CrgA OS=Species OX=1 GN=d PE=1 SV=1\n"
			+"4\t4\trep4\tP4B\tP4B_ENTRY\t44\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-40\t180\t"
			+"sp|P4B|P4B_ENTRY Other HTH-type transcriptional regulator OS=Species OX=1 GN=d PE=1 SV=1\n"
			+"5\t5\trep5\tP5A\tP5A_ENTRY\t45\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P5A|P5A_ENTRY HTH-type transcriptional regulator CrgA OS=Species OX=1 GN=e PE=1 SV=1\n"
			+"5\t5\trep5\tP5B\tP5B_ENTRY\t44\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-40\t180\t"
			+"sp|P5B|P5B_ENTRY HTH-type transcriptional activator AaeR OS=Species OX=1 GN=e PE=1 SV=1\n"
			+"10\t10\trep10\tP10A\tP10A_ENTRY\t45\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P10A|P10A_ENTRY HTH-type transcriptional regulator CrgA OS=Species OX=1 GN=j PE=1 SV=1\n"
			+"10\t10\trep10\tP10B\tP10B_ENTRY\t44\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-40\t180\t"
			+"sp|P10B|P10B_ENTRY HTH-type transcriptional activator AaeR OS=Species OX=1 GN=j PE=1 SV=1\n"
			+"11\t11\trep11\tP11A\tP11A_ENTRY\t45\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P11A|P11A_ENTRY HTH-type transcriptional regulator CrgA OS=Species OX=1 GN=k PE=1 SV=1\n"
			+"11\t11\trep11\tP11B\tP11B_ENTRY\t44\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-40\t180\t"
			+"sp|P11B|P11B_ENTRY Response regulator protein PmrA OS=Species OX=1 GN=k PE=1 SV=1\n"
			+"12\t12\trep12\tP12A\tP12A_ENTRY\t45\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P12A|P12A_ENTRY HTH-type transcriptional regulator CrgA OS=Species OX=1 GN=l PE=1 SV=1\n"
			+"12\t12\trep12\tP12B\tP12B_ENTRY\t44\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-40\t180\t"
			+"sp|P12B|P12B_ENTRY Unlisted HTH-type transcriptional regulator OS=Species OX=1 GN=l PE=1 SV=1\n");
		write(dir, "broad_reviewed.tsv",
			"family_id\tactive_index\tcandidate_class\tpilot_anchor\t"
			+"provisional_decision_class\tprovisional_product\tprimary_support\t"
			+"competing_support\ttransfer_rule_refinement\treviewer_notes\n"
			+"5\t5\tswissprot_full_length_review\tY\tconflicting_evidence\t.\t"
			+"primary\tcompeting\tmanual exception\tpreserve reviewed conflict\n"
			+"10\t10\tswissprot_full_length_review\tY\tconflicting_evidence\t.\t"
			+"primary\tcompeting\tmanual exception\tpreserve reviewed parent conflict\n");
		write(dir, "broad_synonyms.tsv",
			"from_product\tcanonical_product\tnote\n"
			+"RND efflux membrane protein\tefflux RND transporter permease subunit\t"
			+"exercise curated broad normalization\n");
		write(dir, "broad_rules.tsv",
			"rule_id\tconflict_source\tdecision_class\tproduct\tconflict_product\tevidence\n"
			+"rnd_permease\tncbi_top_tier\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\tRND efflux membrane protein\t"
			+"synonym-normalized conflict key\n"
			+"rnd_permease\tncbi_top_tier\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\t"
			+"multidrug efflux RND transporter permease subunit\texact conflict key\n"
			+"hth_regulator\tswissprot\thomology_supported_broad_name\t"
			+"HTH-type transcriptional regulator\t"
			+"HTH-type transcriptional regulator CrgA\texact conflict key\n"
			+"hth_regulator\tswissprot\thomology_supported_broad_name\t"
			+"HTH-type transcriptional regulator\t"
			+"HTH-type transcriptional activator AaeR\texact conflict key\n");
		write(dir, "broad_parent_rules.tsv",
			"parent_id\tconflict_source\tdecision_class\tparent_product\t"
			+"child_product\tevidence\n"
			+"rnd_parent\tncbi_top_tier\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\t"
			+"RND efflux membrane protein\tsynonym-normalized parent key\n"
			+"rnd_parent\tncbi_top_tier\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\t"
			+"multidrug efflux RND transporter permease subunit\texplicit parent key\n"
			+"rnd_parent\tncbi_top_tier\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\t"
			+"heavy metal efflux RND transporter permease subunit\texplicit parent key\n"
			+"abc_permease\tncbi_top_tier\tdomain_supported_broad_name\t"
			+"ABC transporter permease\tpeptide ABC transporter permease\t"
			+"mixed-parent control\n"
			+"abc_permease\tncbi_top_tier\tdomain_supported_broad_name\t"
			+"ABC transporter permease\tamino acid ABC transporter permease\t"
			+"unused second child\n"
			+"hth_parent\tswissprot\thomology_supported_broad_name\t"
			+"HTH-type transcriptional regulator\t"
			+"HTH-type transcriptional regulator CrgA\texplicit parent key\n"
			+"hth_parent\tswissprot\thomology_supported_broad_name\t"
			+"HTH-type transcriptional regulator\t"
			+"HTH-type transcriptional activator AaeR\texplicit parent key\n"
			+"hth_parent\tswissprot\thomology_supported_broad_name\t"
			+"HTH-type transcriptional regulator\t"
			+"Other HTH-type transcriptional regulator\texplicit parent key\n"
			+"response_parent\tswissprot\thomology_supported_broad_name\t"
			+"response regulator\tResponse regulator protein PmrA\t"
			+"mixed-parent control\n"
			+"response_parent\tswissprot\thomology_supported_broad_name\t"
			+"response regulator\tResponse regulator protein QseB\t"
			+"unused second child\n");
		write(dir, "bad_broad_rules.tsv",
			"rule_id\tconflict_source\tdecision_class\tproduct\tconflict_product\tnote\n"
			+"bad\tswissprot\thomology_supported_broad_name\tp1\tp1\tbad sixth field\n"
			+"bad\tswissprot\thomology_supported_broad_name\tp1\tp2\tbad sixth field\n");
		write(dir, "broad_rules_5col.tsv",
			"rule_id\tconflict_source\tdecision_class\tproduct\tconflict_product\n"
			+"rnd5\tncbi_top_tier\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\t"
			+"RND efflux membrane protein\n"
			+"rnd5\tncbi_top_tier\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\t"
			+"multidrug efflux RND transporter permease subunit\n");
		write(dir, "bad_broad_parent_header.tsv",
			"parent_id\tconflict_source\tdecision_class\tparent_product\t"
			+"child_product\tevidence\textra\n"
			+"bad\tswissprot\thomology_supported_broad_name\tp1\tp1\tbad\tbad\n");
		write(dir, "bad_broad_parent_missing_evidence.tsv",
			"parent_id\tconflict_source\tdecision_class\tparent_product\tchild_product\n"
			+"bad\tswissprot\thomology_supported_broad_name\tp1\tp1\n"
			+"bad\tswissprot\thomology_supported_broad_name\tp1\tp2\n");
		write(dir, "bad_broad_parent_duplicate.tsv",
			"parent_id\tconflict_source\tdecision_class\tparent_product\t"
			+"child_product\tevidence\n"
			+"parent1\tswissprot\thomology_supported_broad_name\tp1\tp1\tfirst\n"
			+"parent2\tswissprot\thomology_supported_broad_name\tp2\tp1\tduplicate\n");
		write(dir, "bad_broad_parent_decision.tsv",
			"parent_id\tconflict_source\tdecision_class\tparent_product\t"
			+"child_product\tevidence\n"
			+"bad\tswissprot\tdomain_supported_broad_name\tp1\tp1\twrong decision\n"
			+"bad\tswissprot\tdomain_supported_broad_name\tp1\tp2\twrong decision\n");
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("broad_family.tsv"),
			"ncbi="+dir.resolve("broad_ncbi.tsv"),
			"swiss="+dir.resolve("broad_swiss.tsv"),
			"reviewed="+dir.resolve("broad_reviewed.tsv"),
			"out="+dir.resolve("broad_off_selected.tsv"),
			"evidence="+dir.resolve("broad_off_evidence.tsv"),
			"unresolved="+dir.resolve("broad_off_unresolved.tsv"),
			"summary="+dir.resolve("broad_off_summary.tsv"),
			"labelversion=broad_off",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("broad_family.tsv"),
			"ncbi="+dir.resolve("broad_ncbi.tsv"),
			"swiss="+dir.resolve("broad_swiss.tsv"),
			"reviewed="+dir.resolve("broad_reviewed.tsv"),
			"synonyms="+dir.resolve("broad_synonyms.tsv"),
			"broadconflicts="+dir.resolve("broad_rules.tsv"),
			"out="+dir.resolve("broad_on_provisional_selected.tsv"),
			"evidence="+dir.resolve("broad_on_provisional_evidence.tsv"),
			"unresolved="+dir.resolve("broad_on_provisional_unresolved.tsv"),
			"summary="+dir.resolve("broad_on_provisional_summary.tsv"),
			"labelversion=broad_on_provisional",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy",
			"provisional=t"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("broad_family.tsv"),
			"ncbi="+dir.resolve("broad_ncbi.tsv"),
			"swiss="+dir.resolve("broad_swiss.tsv"),
			"reviewed="+dir.resolve("broad_reviewed.tsv"),
			"synonyms="+dir.resolve("broad_synonyms.tsv"),
			"broadconflicts="+dir.resolve("broad_rules_5col.tsv"),
			"out="+dir.resolve("broad_5col_selected.tsv"),
			"evidence="+dir.resolve("broad_5col_evidence.tsv"),
			"unresolved="+dir.resolve("broad_5col_unresolved.tsv"),
			"summary="+dir.resolve("broad_5col_summary.tsv"),
			"labelversion=broad_5col",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("broad_family.tsv"),
			"ncbi="+dir.resolve("broad_ncbi.tsv"),
			"swiss="+dir.resolve("broad_swiss.tsv"),
			"reviewed="+dir.resolve("broad_reviewed.tsv"),
			"synonyms="+dir.resolve("broad_synonyms.tsv"),
			"broadparents="+dir.resolve("broad_parent_rules.tsv"),
			"out="+dir.resolve("broad_parent_selected.tsv"),
			"evidence="+dir.resolve("broad_parent_evidence.tsv"),
			"unresolved="+dir.resolve("broad_parent_unresolved.tsv"),
			"summary="+dir.resolve("broad_parent_summary.tsv"),
			"labelversion=broad_parent",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("broad_family.tsv"),
			"ncbi="+dir.resolve("broad_ncbi.tsv"),
			"swiss="+dir.resolve("broad_swiss.tsv"),
			"reviewed="+dir.resolve("broad_reviewed.tsv"),
			"synonyms="+dir.resolve("broad_synonyms.tsv"),
			"broadconflicts="+dir.resolve("broad_rules.tsv"),
			"broadparents="+dir.resolve("broad_parent_rules.tsv"),
			"out="+dir.resolve("broad_both_selected.tsv"),
			"evidence="+dir.resolve("broad_both_evidence.tsv"),
			"unresolved="+dir.resolve("broad_both_unresolved.tsv"),
			"summary="+dir.resolve("broad_both_summary.tsv"),
			"labelversion=broad_both",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy"});
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("broad_family.tsv"),
			"ncbi="+dir.resolve("broad_ncbi.tsv"),
			"swiss="+dir.resolve("broad_swiss.tsv"),
			"reviewed="+dir.resolve("broad_reviewed.tsv"),
			"synonyms="+dir.resolve("broad_synonyms.tsv"),
			"broadconflicts="+dir.resolve("broad_rules.tsv"),
			"out="+dir.resolve("broad_on_selected.tsv"),
			"evidence="+dir.resolve("broad_on_evidence.tsv"),
			"unresolved="+dir.resolve("broad_on_unresolved.tsv"),
			"summary="+dir.resolve("broad_on_summary.tsv"),
			"labelversion=broad_on",
			"identitypin=idpin",
			"consensuspin=conpin",
			"mode=hierarchy"});
		final String offSelected=read(dir, "broad_off_selected.tsv");
		final String offSummary=read(dir, "broad_off_summary.tsv");
		check(!offSelected.contains("\n1\t") && !offSelected.contains("\n2\t"),
			"Curated broad conflict fallback should be off by default");
		check(!offSummary.contains("selected_auto_broad_conflict"),
			"Curated broad counters should be off by default for summary parity");
		check(!offSummary.contains("selected_auto_broad_parent"),
			"Curated broad parent counters should be off by default for summary parity");
		final String selected=read(dir, "broad_on_selected.tsv");
		final String evidence=read(dir, "broad_on_evidence.tsv");
		final String unresolved=read(dir, "broad_on_unresolved.tsv");
		final String summary=read(dir, "broad_on_summary.tsv");
		final String provisionalSelected=read(dir, "broad_on_provisional_selected.tsv");
		final String parentSelected=read(dir, "broad_parent_selected.tsv");
		final String parentEvidence=read(dir, "broad_parent_evidence.tsv");
		final String parentUnresolved=read(dir, "broad_parent_unresolved.tsv");
		final String parentSummary=read(dir, "broad_parent_summary.tsv");
		final String bothEvidence=read(dir, "broad_both_evidence.tsv");
		final String bothSummary=read(dir, "broad_both_summary.tsv");
		final String fiveColumnSelected=read(dir, "broad_5col_selected.tsv");
		check(selected.contains("\n1\t1\trep1\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\tbroad_on\n"),
			"NCBI top-tier exact conflict set should select the curated RND broad product");
		check(selected.contains("\n2\t2\trep2\thomology_supported_broad_name\t"
			+"HTH-type transcriptional regulator\tbroad_on\n"),
			"Swiss-Prot exact conflict set should select a homology-supported class");
		check(selected.contains("\n5\t5\trep5\tconflicting_evidence\t\tbroad_on\n")
			&& selected.contains("\n6\t6\trep6\tspecific_full_length\tproduct Keep\tbroad_on\n")
			&& selected.contains("\n10\t10\trep10\tconflicting_evidence\t\tbroad_on\n"),
			"Curated broad conflicts should preserve reviews and existing names");
		check(provisionalSelected.contains("\n2\t2\trep2\t"
			+"homology_supported_broad_name\tHTH-type transcriptional regulator\t"
			+"broad_on_provisional\tHTH-type transcriptional regulator\t"
			+"Swiss-Prot/DIAMOND\tswiss_full\tcoverage_full\tretained_search_hit\t"
			+"not_applicable\tnone\t\n")
			&& provisionalSelected.contains("\n5\t5\trep5\tconflicting_evidence\t\t"
				+"broad_on_provisional\t\tunknown\tunknown\tunknown\tunknown\t"
				+"not_applicable\tnone\t\n"),
			"Schema-2 accepted broad labels and reviewed abstentions should expose "
			+"their actual support");
		check(unresolved.contains("\n3\t3\trep3\tconflicting_evidence\t\t"
			+"competing_top_ncbi_products\t3\t0\tNF_RND3A,NF_RND3B,NF_RND3C")
			&& unresolved.contains("\n4\t4\trep4\tconflicting_evidence\t\t"
				+"competing_swissprot_products\t0\t2"),
			"Exact conflict sets should not match missing or extra products");
		check(evidence.contains("\n1\t1\trep1\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\t"
			+"rnd_permease_ncbi_top_tier_exact_conflict_set\t2\t1")
			&& evidence.contains("\n2\t2\trep2\thomology_supported_broad_name\t"
				+"HTH-type transcriptional regulator\t"
				+"hth_regulator_swissprot_exact_conflict_set\t0\t2"),
			"Curated broad decisions should retain source and rule provenance");
		check(summary.contains("selected_auto\t3\n")
			&& summary.contains("selected_auto_specific\t1\n")
			&& summary.contains("selected_auto_broad\t2\n")
			&& summary.contains("selected_auto_broad_conflict\t2\n")
			&& summary.contains("unresolved_conflict\t6\n")
			&& summary.contains("unresolved_insufficient\t1\n"),
			"Curated broad summary counts changed");
		check(parentSelected.contains("\n1\t1\trep1\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\tbroad_parent\n")
			&& parentSelected.contains("\n2\t2\trep2\thomology_supported_broad_name\t"
				+"HTH-type transcriptional regulator\tbroad_parent\n")
			&& parentSelected.contains("\n3\t3\trep3\tdomain_supported_broad_name\t"
				+"efflux RND transporter permease subunit\tbroad_parent\n")
			&& parentSelected.contains("\n4\t4\trep4\thomology_supported_broad_name\t"
				+"HTH-type transcriptional regulator\tbroad_parent\n"),
			"Broad parents alone should map every listed product to one parent");
		check(parentSelected.contains("\n5\t5\trep5\tconflicting_evidence\t\tbroad_parent\n")
			&& parentSelected.contains("\n6\t6\trep6\tspecific_full_length\tproduct Keep\tbroad_parent\n")
			&& parentSelected.contains("\n10\t10\trep10\tconflicting_evidence\t\tbroad_parent\n"),
			"Broad parents should preserve reviews and existing names");
		check(parentUnresolved.contains("\n8\t8\trep8\tconflicting_evidence\t\t"
			+"competing_top_ncbi_products\t2\t0\tNF_RND8,NF_ABC8")
			&& parentUnresolved.contains("\n9\t9\trep9\tconflicting_evidence\t\t"
				+"competing_top_ncbi_products\t2\t0\tNF_RND9,NF_MISSING9")
			&& parentUnresolved.contains("\n11\t11\trep11\tconflicting_evidence\t\t"
				+"competing_swissprot_products\t0\t2")
			&& parentUnresolved.contains("\n12\t12\trep12\tconflicting_evidence\t\t"
				+"competing_swissprot_products\t0\t2"),
			"Broad parents should abstain on mixed parents and unlisted children");
		check(parentEvidence.contains("\n1\t1\trep1\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\trnd_parent_ncbi_top_tier_parent_vocabulary\t")
			&& parentEvidence.contains("\n2\t2\trep2\thomology_supported_broad_name\t"
				+"HTH-type transcriptional regulator\thth_parent_swissprot_parent_vocabulary\t"),
			"Broad parents should retain parent source and rule provenance");
		check(parentSummary.contains("selected_auto\t5\n")
			&& parentSummary.contains("selected_auto_specific\t1\n")
			&& parentSummary.contains("selected_auto_broad\t4\n")
			&& parentSummary.contains("selected_auto_broad_parent\t4\n")
			&& !parentSummary.contains("selected_auto_broad_conflict")
			&& parentSummary.contains("unresolved_conflict\t4\n")
			&& parentSummary.contains("unresolved_insufficient\t1\n"),
			"Broad-parent summary counts changed");
		check(bothEvidence.contains("rnd_permease_ncbi_top_tier_exact_conflict_set")
			&& bothEvidence.contains("hth_regulator_swissprot_exact_conflict_set")
			&& bothEvidence.contains("rnd_parent_ncbi_top_tier_parent_vocabulary")
			&& bothEvidence.contains("hth_parent_swissprot_parent_vocabulary"),
			"Exact broad conflicts should have precedence while broad parents catch extras");
		check(bothSummary.contains("selected_auto\t5\n")
			&& bothSummary.contains("selected_auto_broad_conflict\t2\n")
			&& bothSummary.contains("selected_auto_broad_parent\t2\n")
			&& bothSummary.contains("unresolved_conflict\t4\n"),
			"Combined exact-set and parent-vocabulary summary counts changed");
		check(fiveColumnSelected.contains("\n1\t1\trep1\tdomain_supported_broad_name\t"
			+"efflux RND transporter permease subunit\tbroad_5col\n"),
			"Five-column broad-conflict rules should remain supported");
		final String originalBroadRules=read(dir, "broad_rules.tsv");
		final Path alias=dir.resolve("broad_rules_alias.tsv");
		Files.createLink(alias, dir.resolve("broad_rules.tsv"));
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("broad_family.tsv"),
				"ncbi="+dir.resolve("broad_ncbi.tsv"),
				"swiss="+dir.resolve("broad_swiss.tsv"),
				"reviewed="+dir.resolve("broad_reviewed.tsv"),
				"broadconflicts="+dir.resolve("broad_rules.tsv"),
				"out="+alias,
				"evidence="+dir.resolve("broad_alias_evidence.tsv"),
				"unresolved="+dir.resolve("broad_alias_unresolved.tsv"),
				"summary="+dir.resolve("broad_alias_summary.tsv"),
				"labelversion=broad_alias",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Output aliases input");
		}
		check(ok, "Broad-conflict input hardlink aliases should be rejected");
		check(read(dir, "broad_rules.tsv").equals(originalBroadRules),
			"Broad-conflict input hardlink must be rejected before writing");
		final String originalBroadParentRules=read(dir, "broad_parent_rules.tsv");
		final Path parentAlias=dir.resolve("broad_parent_rules_alias.tsv");
		Files.createLink(parentAlias, dir.resolve("broad_parent_rules.tsv"));
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("broad_family.tsv"),
				"ncbi="+dir.resolve("broad_ncbi.tsv"),
				"swiss="+dir.resolve("broad_swiss.tsv"),
				"reviewed="+dir.resolve("broad_reviewed.tsv"),
				"broadparents="+dir.resolve("broad_parent_rules.tsv"),
				"out="+parentAlias,
				"evidence="+dir.resolve("broad_parent_alias_evidence.tsv"),
				"unresolved="+dir.resolve("broad_parent_alias_unresolved.tsv"),
				"summary="+dir.resolve("broad_parent_alias_summary.tsv"),
				"labelversion=broad_parent_alias",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Output aliases input");
		}
		check(ok, "Broad-parent input hardlink aliases should be rejected");
		check(read(dir, "broad_parent_rules.tsv").equals(originalBroadParentRules),
			"Broad-parent input hardlink must be rejected before writing");
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("broad_family.tsv"),
				"ncbi="+dir.resolve("broad_ncbi.tsv"),
				"swiss="+dir.resolve("broad_swiss.tsv"),
				"reviewed="+dir.resolve("broad_reviewed.tsv"),
				"broadconflicts="+dir.resolve("bad_broad_rules.tsv"),
				"out="+dir.resolve("bad_broad_selected.tsv"),
				"evidence="+dir.resolve("bad_broad_evidence.tsv"),
				"unresolved="+dir.resolve("bad_broad_unresolved.tsv"),
				"summary="+dir.resolve("bad_broad_summary.tsv"),
				"labelversion=bad_broad",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Unexpected broad-conflict header");
		}
		check(ok, "Broad-conflict rules should reject unexpected optional columns");
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("broad_family.tsv"),
				"ncbi="+dir.resolve("broad_ncbi.tsv"),
				"swiss="+dir.resolve("broad_swiss.tsv"),
				"broadparents="+dir.resolve("bad_broad_parent_header.tsv"),
				"out="+dir.resolve("bad_broad_parent_header_selected.tsv"),
				"evidence="+dir.resolve("bad_broad_parent_header_evidence.tsv"),
				"unresolved="+dir.resolve("bad_broad_parent_header_unresolved.tsv"),
				"summary="+dir.resolve("bad_broad_parent_header_summary.tsv"),
				"labelversion=bad_broad_parent_header",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Unexpected broad-parent header");
		}
		check(ok, "Broad-parent rules should reject unexpected optional columns");
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("broad_family.tsv"),
				"ncbi="+dir.resolve("broad_ncbi.tsv"),
				"swiss="+dir.resolve("broad_swiss.tsv"),
				"broadparents="+dir.resolve("bad_broad_parent_missing_evidence.tsv"),
				"out="+dir.resolve("bad_broad_parent_missing_selected.tsv"),
				"evidence="+dir.resolve("bad_broad_parent_missing_evidence.tsv.out"),
				"unresolved="+dir.resolve("bad_broad_parent_missing_unresolved.tsv"),
				"summary="+dir.resolve("bad_broad_parent_missing_summary.tsv"),
				"labelversion=bad_broad_parent_missing",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Unexpected broad-parent header");
		}
		check(ok, "Broad-parent rules should reject missing evidence columns");
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("broad_family.tsv"),
				"ncbi="+dir.resolve("broad_ncbi.tsv"),
				"swiss="+dir.resolve("broad_swiss.tsv"),
				"broadparents="+dir.resolve("bad_broad_parent_duplicate.tsv"),
				"out="+dir.resolve("bad_broad_parent_duplicate_selected.tsv"),
				"evidence="+dir.resolve("bad_broad_parent_duplicate_evidence.tsv"),
				"unresolved="+dir.resolve("bad_broad_parent_duplicate_unresolved.tsv"),
				"summary="+dir.resolve("bad_broad_parent_duplicate_summary.tsv"),
				"labelversion=bad_broad_parent_duplicate",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Duplicate broad-parent product");
		}
		check(ok, "Broad-parent rules should reject duplicate source/product children");
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("broad_family.tsv"),
				"ncbi="+dir.resolve("broad_ncbi.tsv"),
				"swiss="+dir.resolve("broad_swiss.tsv"),
				"broadparents="+dir.resolve("bad_broad_parent_decision.tsv"),
				"out="+dir.resolve("bad_broad_parent_decision_selected.tsv"),
				"evidence="+dir.resolve("bad_broad_parent_decision_evidence.tsv"),
				"unresolved="+dir.resolve("bad_broad_parent_decision_unresolved.tsv"),
				"summary="+dir.resolve("bad_broad_parent_decision_summary.tsv"),
				"labelversion=bad_broad_parent_decision",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Bad broad-parent decision");
		}
		check(ok, "Broad-parent rules should reject source/decision mismatches");
	}

	private static void checkDuplicateSynonym(final Path dir) throws Exception{
		write(dir, "bad_synonyms.tsv",
			"from_product\tcanonical_product\tnote\n"
			+"DNA topoisomerase (ATP-hydrolyzing) subunit A\tDNA gyrase subunit A\ta\n"
			+"DNA topoisomerase (ATP-hydrolyzing) subunit A\tother\tb\n");
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"synonyms="+dir.resolve("bad_synonyms.tsv"),
				"out="+dir.resolve("bad_syn_selected.tsv"),
				"evidence="+dir.resolve("bad_syn_evidence.tsv"),
				"unresolved="+dir.resolve("bad_syn_unresolved.tsv"),
				"summary="+dir.resolve("bad_syn_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Duplicate synonym source");
		}
		check(ok, "duplicate synonym sources should crash loudly");
	}

	private static void checkDuplicateFoldedSynonym(final Path dir) throws Exception{
		write(dir, "bad_folded_synonyms.tsv",
			"from_product\tcanonical_product\tnote\n"
			+"product X\tproduct Y\ta\n"
			+"PRODUCT X\tproduct Z\tb\n");
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("family_identity.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"synonyms="+dir.resolve("bad_folded_synonyms.tsv"),
				"out="+dir.resolve("bad_fold_selected.tsv"),
				"evidence="+dir.resolve("bad_fold_evidence.tsv"),
				"unresolved="+dir.resolve("bad_fold_unresolved.tsv"),
				"summary="+dir.resolve("bad_fold_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin",
				"ignorecase=t"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Duplicate synonym source");
		}
		check(ok, "duplicate folded synonym sources should crash loudly");
	}

	private static void checkUnknownFamilyType(final Path dir) throws Exception{
		write(dir, "unknown_family.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\n"
			+"1\t1\t0\trep1\t10\t100\n");
		write(dir, "unknown_ncbi.tsv", modelHeader()
			+"1\t1\trep1\tNF_X\t200\t100\tY\t190\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tmystery\tY\tNCBIFAM\tgeneX\tproduct X\n");
		write(dir, "unknown_swiss.tsv", swissHeader()
			+"1\t1\trep1\tP1\tP1_ENTRY\t99\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P1|P1_ENTRY product X OS=Species OX=1 GN=x PE=1 SV=1\n");
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("unknown_family.tsv"),
			"ncbi="+dir.resolve("unknown_ncbi.tsv"),
			"swiss="+dir.resolve("unknown_swiss.tsv"),
			"out="+dir.resolve("unknown_selected.tsv"),
			"evidence="+dir.resolve("unknown_evidence.tsv"),
			"unresolved="+dir.resolve("unknown_unresolved.tsv"),
			"summary="+dir.resolve("unknown_summary.tsv"),
			"labelversion=test_v1",
			"identitypin=idpin",
			"consensuspin=conpin"});
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("unknown_family.tsv"),
				"ncbi="+dir.resolve("unknown_ncbi.tsv"),
				"swiss="+dir.resolve("unknown_swiss.tsv"),
				"out="+dir.resolve("bad_unknown_selected.tsv"),
				"evidence="+dir.resolve("bad_unknown_evidence.tsv"),
				"unresolved="+dir.resolve("bad_unknown_unresolved.tsv"),
				"summary="+dir.resolve("bad_unknown_summary.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin",
				"mode=hierarchy"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Unsupported NCBI naming family_type mystery");
		}
		check(ok, "hierarchy mode should reject unknown NCBI family types");
	}

	private static void checkBadMemberCounts(final Path dir) throws Exception{
		write(dir, "missing_members.tsv",
			"active_index\tfamily_id\tignored\trep_id\traw_residues\n"
			+"1\t1\t0\trep1\t100\n");
		boolean ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("missing_members.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("missing_members_selected.tsv"),
				"evidence="+dir.resolve("missing_members_evidence.tsv"),
				"unresolved="+dir.resolve("missing_members_unresolved.tsv"),
				"summary="+dir.resolve("missing_members_summary.tsv"),
				"membercoverage="+dir.resolve("missing_members_coverage.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("family_identity requires members and raw_residues");
		}
		check(ok, "member reports should require members in family_identity");

		write(dir, "missing_member_raw.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\n"
			+"1\t1\t0\trep1\t10\n");
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("missing_member_raw.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("missing_member_selected.tsv"),
				"evidence="+dir.resolve("missing_member_evidence.tsv"),
				"unresolved="+dir.resolve("missing_member_unresolved.tsv"),
				"summary="+dir.resolve("missing_member_summary.tsv"),
				"membercoverage="+dir.resolve("missing_member_coverage.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("family_identity requires members and raw_residues");
		}
		check(ok, "member reports should require raw_residues in family_identity");

		write(dir, "bad_members.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\n"
			+"1\t1\t0\trep1\t0\t100\n");
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("bad_members.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("bad_members_selected.tsv"),
				"evidence="+dir.resolve("bad_members_evidence.tsv"),
				"unresolved="+dir.resolve("bad_members_unresolved.tsv"),
				"summary="+dir.resolve("bad_members_summary.tsv"),
				"memberunresolved="+dir.resolve("bad_members_unresolved_members.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Bad members");
		}
		check(ok, "member reports should reject nonpositive members");

		write(dir, "bad_raw_residues.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\n"
			+"1\t1\t0\trep1\t10\t-1\n");
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("bad_raw_residues.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("bad_raw_selected.tsv"),
				"evidence="+dir.resolve("bad_raw_evidence.tsv"),
				"unresolved="+dir.resolve("bad_raw_unresolved.tsv"),
				"summary="+dir.resolve("bad_raw_summary.tsv"),
				"membercoverage="+dir.resolve("bad_raw_coverage.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Bad raw_residues");
		}
		check(ok, "member reports should reject negative raw_residues");

		write(dir, "duplicate_member_column.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\tmembers\n"
			+"1\t1\t0\trep1\t10\t100\t10\n");
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("duplicate_member_column.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("dup_member_selected.tsv"),
				"evidence="+dir.resolve("dup_member_evidence.tsv"),
				"unresolved="+dir.resolve("dup_member_unresolved.tsv"),
				"summary="+dir.resolve("dup_member_summary.tsv"),
				"membercoverage="+dir.resolve("dup_member_coverage.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(IllegalArgumentException iae){
			ok=iae.getMessage().contains("Duplicate family_identity column members");
		}
		check(ok, "member reports should reject duplicate member columns");

		final String overflow=read(dir, "family_identity.tsv")
			.replace("1\t1\t0\trep1\t10\t100\n",
				"1\t1\t0\trep1\t9223372036854775807\t100\n")
			.replace("2\t2\t0\trep2\t10\t100\n",
				"2\t2\t0\trep2\t9223372036854775807\t100\n");
		write(dir, "overflow_members.tsv", overflow);
		ok=false;
		try{
			AssignProteinNames.main(new String[]{
				"families="+dir.resolve("overflow_members.tsv"),
				"ncbi="+dir.resolve("ncbi.tsv"),
				"swiss="+dir.resolve("swiss.tsv"),
				"out="+dir.resolve("overflow_selected.tsv"),
				"evidence="+dir.resolve("overflow_evidence.tsv"),
				"unresolved="+dir.resolve("overflow_unresolved.tsv"),
				"summary="+dir.resolve("overflow_summary.tsv"),
				"membercoverage="+dir.resolve("overflow_member_coverage.tsv"),
				"labelversion=test_v1",
				"identitypin=idpin",
				"consensuspin=conpin"});
		}catch(ArithmeticException ae){
			ok=true;
		}
		check(ok, "member reports should crash loudly on long overflow");
	}

	private static void checkZeroRawMemberCoverage(final Path dir) throws Exception{
		write(dir, "zero_raw_family.tsv",
			"active_index\tfamily_id\tignored\trep_id\tmembers\traw_residues\n"
			+"1\t1\t0\trep1\t3\t0\n");
		write(dir, "zero_raw_ncbi.tsv", modelHeader()
			+"1\t1\trep1\tNF_A\t200\t100\tY\t190\t100\tY\t1.000000000000\t"
			+"1.000000000000\tY\tY\thmm\tequivalog\tY\tNCBIFAM\tgeneA\t"
			+"product A\n");
		write(dir, "zero_raw_swiss.tsv", swissHeader()
			+"1\t1\trep1\tP1\tP1_ENTRY\t99\t100\t100\t100\t1.000000000000\t"
			+"1.000000000000\t1\t100\t1\t100\t1e-50\t200\t"
			+"sp|P1|P1_ENTRY product A OS=Species OX=1 GN=a PE=1 SV=1\n");
		AssignProteinNames.main(new String[]{
			"families="+dir.resolve("zero_raw_family.tsv"),
			"ncbi="+dir.resolve("zero_raw_ncbi.tsv"),
			"swiss="+dir.resolve("zero_raw_swiss.tsv"),
			"out="+dir.resolve("zero_raw_selected.tsv"),
			"evidence="+dir.resolve("zero_raw_evidence.tsv"),
			"unresolved="+dir.resolve("zero_raw_unresolved.tsv"),
			"summary="+dir.resolve("zero_raw_summary.tsv"),
			"membercoverage="+dir.resolve("zero_raw_member_coverage.tsv"),
			"labelversion=test_v1",
			"identitypin=idpin",
			"consensuspin=conpin"});
		check(read(dir, "zero_raw_member_coverage.tsv").contains(
			"named_auto\tspecific_full_length\t\t1\t3\t0\t1\t3\t0\t1.000000000000\t.\n"),
			"zero total raw residues should emit a missing raw-residue fraction");
	}

	private static String modelHeader(){
		return "family_id\tactive_index\trep_id\tncbi_accession\tfull_score\t"
			+"sequence_cutoff\tsequence_pass\tdomain_score\tdomain_cutoff\tdomain_pass\t"
			+"hmm_cov\tconsensus_cov\tmutual_coverage_pass\tsame_hit_native_eligible\t"
			+"hmm_name\tfamily_type\tfor_naming\tsource\tgene_symbol\tproduct_name\n";
	}

	private static String swissHeader(){
		return "family_id\tactive_index\tconsensus_id\tuniprot_accession\tuniprot_entry\t"
			+"identity\taligned_len\tquery_len\tsubject_len\tquery_cov\tsubject_cov\t"
			+"qstart\tqend\tsstart\tsend\tevalue\tbitscore\ttitle\n";
	}

	private static String domainHeader(){
		return BbgapCoverageAudit.DOMAIN_HEADER+"\tfor_naming\tcoverage_failure\n";
	}

	private static String swappedDomainHeader(){
		final String[] terms=BbgapCoverageAudit.DOMAIN_HEADER.split("\t");
		final String temp=terms[10];
		terms[10]=terms[13];
		terms[13]=temp;
		return String.join("\t", terms)+"\tfor_naming\tcoverage_failure\n";
	}

	private static void write(final Path dir, final String name, final String text) throws Exception{
		Files.write(dir.resolve(name), text.getBytes(StandardCharsets.UTF_8));
	}

	private static String read(final Path dir, final String name) throws Exception{
		return new String(Files.readAllBytes(dir.resolve(name)), StandardCharsets.UTF_8);
	}

	private static void checkWidths(final String text, final int width, final String name){
		final String[] lines=text.split("\n");
		for(int i=1; i<lines.length; i++){
			final int found=lines[i].split("\t", -1).length;
			check(found==width, name+" line "+(i+1)+" width "+found+" != "+width);
		}
	}

	private static void checkMemberCoverageTotals(final String text, final int families,
			final long members, final long rawResidues){
		final String[] lines=text.split("\n");
		long familySum=0, memberSum=0, rawSum=0;
		for(int i=1; i<lines.length; i++){
			final String[] row=lines[i].split("\t", -1);
			check(row.length==11, "member coverage line "+(i+1)+" width "+row.length+" != 11");
			familySum+=Integer.parseInt(row[3]);
			memberSum+=Long.parseLong(row[4]);
			rawSum+=Long.parseLong(row[5]);
			check(Integer.parseInt(row[6])==families, "bad member coverage family denominator");
			check(Long.parseLong(row[7])==members, "bad member coverage member denominator");
			check(Long.parseLong(row[8])==rawResidues,
				"bad member coverage raw-residue denominator");
		}
		check(familySum==families, "member coverage families do not sum");
		check(memberSum==members, "member coverage members do not sum");
		check(rawSum==rawResidues, "member coverage raw residues do not sum");
	}

	private static void check(final boolean ok, final String reason){
		if(!ok){throw new AssertionError(reason);}
	}
}
