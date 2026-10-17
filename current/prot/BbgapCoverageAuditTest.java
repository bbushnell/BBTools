package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/** Literal evidence fixtures for lossless domain inventory and strict joins.
 * @author Brian Bushnell, Keqing
 */
public final class BbgapCoverageAuditTest {

	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Required: new fixture directory");
		final Path dir=Paths.get(args[0]);
		Files.createDirectory(dir);
		final String decisions="family_id\tactive_index\trep_id\tdecision_class\tproduct\trule_or_reason\n"
			+"10\t1\trep10\tunnamed\t\tinsufficient_native_ncbi_evidence\n"
			+"20\t2\trep20\tunnamed\t\tinsufficient_native_ncbi_evidence\n"
			+"30\t3\trep30\tspecific_full_length\talready named\texisting_rule\n"
			+"#raw_provisional_evidence\n"
			+"family_id\tactive_index\trep_id\thit_id\tselected\tsource_methods\t"
			+"accession\tproduct\tsupport_tag\tcoverage_tag\tscore_tag\tsequence_score\t"
			+"domain_score\tquery_cov\ttarget_cov\tsequence_pass\tdomain_pass\t"
			+"coverage_pass\tnative_eligible\tfor_naming\toriginal_decision\t"
			+"original_reason\n"
			+"10\t1\trep10\tncbi:0\tN\tNCBIfam/HMMER\tNF_IGNORED\tignored\t"
			+"native_partial\tcoverage_full\tsequence_and_domain_score_low\t1\t1\t1\t"
			+"1\tN\tN\tY\tN\tY\tunnamed\tinsufficient_native_ncbi_evidence\n";
		final String targets=BbgapCoverageAudit.TARGET_HEADER+"\n"
			+target(10, 1, "NF_A", "equivalog", "Y", "alpha", 2)
			+target(10, 1, "NF_B", "equivalog", "Y", "beta", 1)
			+target(10, 1, "NF_C", "exception", "N", "non-naming", 1)
			+target(20, 2, "NF_D", "equivalog", "Y", "below-cutoff", 1)
			+target(20, 2, "NF_ZERO", "equivalog", "Y", "no reported domains", 0)
			+target(30, 3, "NF_E", "equivalog", "Y", "already named", 1);
		final String a1=domain(10, 1, "NF_A", 1, 2, 90, 0.2, 0.4, "equivalog", "alpha");
		final String a2=domain(10, 1, "NF_A", 2, 2, 80, 0.9, 0.76, "equivalog", "alpha");
		final String b=domain(10, 1, "NF_B", 1, 1, 85, 0.76, 0.9, "equivalog", "beta");
		final String c=domain(10, 1, "NF_C", 1, 1, 95, 1, 0.8, "exception", "non-naming");
		final String domains=BbgapCoverageAudit.DOMAIN_HEADER+"\n"+a1+a2+b+c
			+domain(20, 2, "NF_D", 1, 1, 40, 1, 1, "equivalog", "below-cutoff")
			+domain(30, 3, "NF_E", 1, 1, 95, 1, 1, "equivalog", "already named");
		write(dir, "decisions.tsv", decisions);
		write(dir, "targets.tsv", targets);
		write(dir, "domains.tsv", domains);
		run(dir, "domains.tsv", "audit");
		final String rows=read(dir, "audit.domains.tsv");
		check(rows.contains(withoutNewline(a1)+"\tY\tboth\n"), "Highest-score partial domain lost");
		check(rows.contains(withoutNewline(a2)+"\tY\tquery_only\n"), "Lower-score alternative domain lost");
		check(rows.contains(withoutNewline(b)+"\tY\tmodel_only\n"), "Equal-type alternative model lost");
		check(rows.contains(withoutNewline(c)+"\tN\tnone\n"), "Non-naming evidence was discarded");
		check(!rows.contains("NF_D") && !rows.contains("NF_E"), "Failed score or already-named family leaked into inventory");
		checkWidths(rows, 30);
		final String familyRows=read(dir, "audit.families.tsv");
		check(familyRows.contains("\n10\t1\trep10\tunnamed\t\tinsufficient_native_ncbi_evidence"
			+"\t3\t4\t1\t1\t1\t1\t0\t1\t1\t1\n"), "Family flags do not conserve all domain alternatives");
		check(familyRows.contains("\n20\t2\trep20\tunnamed\t\tinsufficient_native_ncbi_evidence"
			+"\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"), "Zero-support candidate family disappeared");
		checkWidths(familyRows, 16);
		check(read(dir, "audit.summary.tsv").equals("metric\tvalue\ninput_families\t3\n"
			+"candidate_families\t2\nfamilies_with_scorepass_domains\t1\n"
			+"families_with_coverage_fail_domains\t1\nscorepass_domain_rows\t4\nmin_coverage\t0.8\n"),
			"Inventory summary differs from the literal oracle");
		// Formatting the same decimal with more zeroes must not change its gate.
		write(dir, "boundary.tsv", domains.replace("\t0.8\t", "\t0.800000000000\t"));
		run(dir, "boundary.tsv", "boundary");
		check(read(dir, "boundary.domains.tsv").contains(
			withoutNewline(c).replace("\t0.8\t", "\t0.800000000000\t")+"\tN\tnone\n"),
			"Exact 0.8 coverage was rejected after adding decimal zeroes");
		check(read(dir, "boundary.summary.tsv").equals(read(dir, "audit.summary.tsv")),
			"Equivalent decimal formatting changed summary counts");
		write(dir, "below.tsv", domains.replace("\t0.8\t", "\t0.799999999999\t"));
		run(dir, "below.tsv", "below");
		check(read(dir, "below.domains.tsv").contains(
			withoutNewline(c).replace("\t0.8\t", "\t0.799999999999\t")+"\tN\tquery_only\n"),
			"A value strictly below 0.8 must still fail; do not add an epsilon");
		write(dir, "missing.tsv", domains.replace("NF_A", "NF_MISSING"));
		fails(dir, "missing.tsv", "missing", "Domain lacks target");
		write(dir, "duplicate.tsv", domains+a1);
		fails(dir, "duplicate.tsv", "duplicate", "Duplicate/invalid domain index");
		write(dir, "truncated.tsv", domains.replace(a2, ""));
		fails(dir, "truncated.tsv", "truncated", "Incomplete reported domains");
		write(dir, "missingmodel.tsv", domains.replace(a1, "").replace(a2, ""));
		fails(dir, "missingmodel.tsv", "missingmodel", "Incomplete reported domains");
		write(dir, "count.tsv", domains.replace(a2,
			domain(10, 1, "NF_A", 2, 3, 80, 0.9, 0.76, "equivalog", "alpha")));
		fails(dir, "count.tsv", "count", "Reported-domain count mismatch");
		write(dir, "missingfailed.tsv", domains.replace(
			domain(20, 2, "NF_D", 1, 1, 40, 1, 1, "equivalog", "below-cutoff"), ""));
		fails(dir, "missingfailed.tsv", "missingfailed", "Incomplete reported domains");
		write(dir, "identity.tsv", domains.replace("rep10", "wrong_rep"));
		fails(dir, "identity.tsv", "identity", "Family identity mismatch");
		write(dir, "nan.tsv", domains.replace("\t0.76\t", "\tNaN\t"));
		fails(dir, "nan.tsv", "nan", "Non-finite");
		// An existing output alias must be rejected even with overwrite enabled.
		write(dir, "alias.domains.tsv", domains);
		boolean rejected=false;
		try{
			BbgapCoverageAudit.main(new String[]{"decisions="+dir.resolve("decisions.tsv"),
				"targets="+dir.resolve("targets.tsv"), "domains="+dir.resolve("alias.domains.tsv"),
				"out="+dir.resolve("alias"), "ow=t"});
		}catch(IllegalArgumentException e){rejected=e.getMessage().contains("aliases");}
		check(rejected && read(dir, "alias.domains.tsv").equals(domains), "Input alias was not preserved");
		System.err.println("BBGAP_COVERAGE_AUDIT_PASS alternatives=true joins=true counts=true aliases=true completeness=true decimals=true");
	}

	private static String target(final int family, final int active, final String accession,
			final String type, final String naming, final String product, final int reported){
		return family+"\t"+active+"\trep"+family+"\t"+accession+"\tmodel\t1e-20\t100\t1e-20"
			+"\t95\t"+reported+"\t50\t50\t100\t"+type+"\t"+naming+"\ttest\tgene\t\t"+product+"\n";
	}

	private static String domain(final int family, final int active, final String accession,
			final int index, final int count, final int score, final double hmmCov,
			final double queryCov, final String type, final String product){
		return family+"\t"+active+"\trep"+family+"\t100\t"+accession+"\tmodel\t100\t"
			+index+"\t"+count+"\t1e-20\t100\t1e-20\t1e-20\t"+score
			+"\t1\t"+Math.round(100*hmmCov)+"\t1\t"+Math.round(100*queryCov)
			+"\t1\t"+Math.round(100*queryCov)+"\t"+hmmCov+"\t"+queryCov
			+"\t50\t50\t"+type+"\ttest\tgene\t"+product+"\n";
	}
	private static void run(final Path dir, final String domainFile, final String prefix) throws Exception{
		BbgapCoverageAudit.main(new String[]{"decisions="+dir.resolve("decisions.tsv"),
			"targets="+dir.resolve("targets.tsv"), "domains="+dir.resolve(domainFile),
			"out="+dir.resolve(prefix)});
	}
	private static void fails(final Path dir, final String input, final String prefix,
			final String diagnostic) throws Exception{
		boolean rejected=false;
		try{run(dir, input, prefix);}
		catch(IllegalArgumentException e){rejected=e.getMessage().contains(diagnostic);}
		check(rejected && !Files.exists(dir.resolve(prefix+".domains.tsv")),
			"Expected diagnostic without output: "+diagnostic);
	}
	private static void checkWidths(final String text, final int width){
		for(String line : text.split("\n")){
			check(line.split("\t", -1).length==width, "Wrong output width: "+line);
		}
	}
	private static String withoutNewline(final String text){return text.substring(0, text.length()-1);}
	private static void write(final Path dir, final String file, final String text) throws Exception{
		Files.write(dir.resolve(file), text.getBytes(StandardCharsets.UTF_8));
	}
	private static String read(final Path dir, final String file) throws Exception{
		return new String(Files.readAllBytes(dir.resolve(file)), StandardCharsets.UTF_8);
	}
	private static void check(final boolean valid, final String diagnostic){
		if(!valid){throw new AssertionError(diagnostic);}
	}
}
