package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.LinkedHashMap;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.IntHashMap2;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Exports every native-score-passing domain for unnamed BBGAP families.
 * This is an evidence inventory, not a naming rule. No domain or model is
 * chosen by score or precedence; the original metadata remains available.
 * @author Brian Bushnell, Keqing
 */
public final class BbgapCoverageAudit {

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, BbgapCoverageAudit.class, false);
		try{new BbgapCoverageAudit(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	BbgapCoverageAudit(final String[] args) throws Exception{
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			final String value=pair.length==2 ? pair[1] : null;
			if(key.equals("decisions")){decisions=value;}
			else if(key.equals("targets")){targets=value;}
			else if(key.equals("domains")){domains=value;}
			else if(key.equals("out")){out=value;}
			else if(key.equals("mincov")){minCoverage=Double.parseDouble(value);}
			else if(key.equals("ow") || key.equals("overwrite")){overwrite=Parse.parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		require(decisions!=null && targets!=null && domains!=null && out!=null,
			"Required: decisions= targets= domains= out=PREFIX");
		require(minCoverage>0 && minCoverage<=1, "mincov must be finite in (0,1]: "+minCoverage);
		outputPaths=new String[]{out+".domains.tsv", out+".families.tsv", out+".summary.tsv"};
		final ArrayList<Path> paths=new ArrayList<Path>();
		for(String input : new String[]{decisions, targets, domains}){
			final Path path=Paths.get(input).toRealPath();
			require(Files.isRegularFile(path) && Files.isReadable(path), "Unreadable input: "+input);
			paths.add(path);
		}
		for(String output : outputPaths){
			final Path path=Paths.get(output).toAbsolutePath().normalize();
			require(overwrite || !Files.exists(path), "Output exists; use a new prefix: "+output);
			for(Path other : paths){
				require(!path.equals(other) && !(Files.exists(path) && Files.isSameFile(path, other)),
					"Output aliases an input or another output: "+output+" and "+other);
			}
			paths.add(path);
		}
	}

	void process(){
		loadDecisions();
		loadTargets();
		loadDomains();
		validateDomainCounts();
		// All joins and numeric inputs are validated before output writers start.
		writeDomains();
		writeFamilies();
		writeSummary();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Evidence           ----------------*/
	/*--------------------------------------------------------------*/

	private void loadDecisions(){
		final ByteFile bf=ByteFile.makeByteFile(decisions, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNo=0;
		int width=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line); lineNo++;
				if(lineNo==1){
					require(lp.terms()>=6 && lp.termEquals("family_id", 0) &&
						lp.termEquals("active_index", 1) && lp.termEquals("rep_id", 2) &&
						lp.termEquals("decision_class", 3) && lp.termEquals("product", 4) &&
						lp.termEquals("rule_or_reason", 5), "Unexpected decisions header: "+decisions);
					width=lp.terms(); decisionHeader=line.clone(); continue;
				}
				if(line.length>0 && line[0]=='#'){
					require(lp.terms()==1 && lp.termEquals("#raw_provisional_evidence", 0),
						"Unexpected decisions section marker at line "+lineNo);
					break;
				}
				require(lp.terms()==width, "Wrong decisions width at line "+lineNo);
				final int id=lp.parseInt(0), active=lp.parseInt(1);
				require(id>=0 && active>=0 && familyIndex.get(id)<0 && activeIndex.get(active)<0,
					"Duplicate/negative family or active ID at decisions line "+lineNo);
				final boolean candidate=lp.termEquals("unnamed", 3) &&
					lp.termEquals("insufficient_native_ncbi_evidence", 5);
				require(!candidate || lp.termEquals("", 4), "Unnamed candidate has a product: "+id);
				final Family family=new Family(id, active, lp.parseString(2), candidate,
					candidate ? line.clone() : null);
				require(!family.rep.isEmpty(), "Missing representative for family "+id);
				familyIndex.put(id, families.size()); activeIndex.put(active, families.size());
				families.add(family);
				if(candidate){candidateFamilies++;}
			}
		}finally{require(!bf.close(), "I/O error reading "+decisions);}
		require(lineNo>1, "Decisions file has no family records: "+decisions);
	}

	private void loadTargets(){
		final ByteFile bf=ByteFile.makeByteFile(targets, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNo=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line); lineNo++;
				if(lineNo==1){checkHeader(lp, TARGET_HEADER, targets); continue;}
				require(lp.terms()==19, "Wrong target width at line "+lineNo);
				final Family family=family(lp, targets, lineNo);
				if(!family.candidate){continue;}
				final String accession=lp.parseString(3);
				require(!accession.isEmpty(), "Missing target accession at line "+lineNo);
				final Model model=new Model(accession, finite(lp, 6), finite(lp, 10),
					finite(lp, 11), lp.parseInt(12), lp.parseString(13), lp.parseString(14), lp.parseInt(9));
				require(model.reportedDomains>=0, "Negative reported-domain count for family "+family.id+" model "+accession);
				require(model.hmmLength>0 && (model.naming.equals("Y") || model.naming.equals("N")),
					"Invalid target length/naming flag for family "+family.id+" model "+accession);
				require(family.models.put(accession, model)==null,
					"Duplicate target for family "+family.id+" model "+accession);
			}
		}finally{require(!bf.close(), "I/O error reading "+targets);}
		require(lineNo>0, "Empty targets file: "+targets);
	}

	private void loadDomains(){
		final ByteFile bf=ByteFile.makeByteFile(domains, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNo=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line); lineNo++;
				if(lineNo==1){checkHeader(lp, DOMAIN_HEADER, domains); continue;}
				require(lp.terms()==28, "Wrong domain width at line "+lineNo);
				final Family family=family(lp, domains, lineNo);
				if(!family.candidate){continue;}
				final String accession=lp.parseString(4);
				final Model model=family.models.get(accession);
				require(model!=null, "Domain lacks target: family "+family.id+" model "+accession);
				final int length=lp.parseInt(3), index=lp.parseInt(7), count=lp.parseInt(8);
				require(length>0 && (family.length==0 || family.length==length),
					"Inconsistent consensus length for family "+family.id+" at domain line "+lineNo);
				family.length=length;
				require(count==model.reportedDomains,
					"Reported-domain count mismatch for family "+family.id+" model "+accession+
					": target="+model.reportedDomains+" domain="+count+"; use corrected compact targets");
				require(index>0 && index<=count && !model.seen.get(index),
					"Duplicate/invalid domain index "+index+" for family "+family.id+" model "+accession);
				model.seen.set(index);
				require(model.hmmLength==lp.parseInt(6) && lp.termEquals(model.type, 24) &&
					model.score==finite(lp, 10) && model.cutoff==finite(lp, 22) &&
					model.domainCutoff==finite(lp, 23),
					"Target/domain metadata disagree for family "+family.id+" model "+accession);
				final double score=finite(lp, 13), hmmCov=finite(lp, 20), queryCov=finite(lp, 21);
				require(hmmCov>=0 && hmmCov<=1 && queryCov>=0 && queryCov<=1,
					"Invalid coverage for family "+family.id+" model "+accession);
				if(model.score<model.cutoff || score<model.domainCutoff){continue;}
				final int flags=(hmmCov<minCoverage ? 1 : 0)|(queryCov<minCoverage ? 2 : 0);
				family.counts[flags]++; model.passing++;
				if(model.naming.equals("Y")){family.namingCounts[flags]++;}
				retained.add(new Domain(family, model, flags, line.clone()));
			}
		}finally{require(!bf.close(), "I/O error reading "+domains);}
		require(lineNo>0, "Empty domains file: "+domains);
	}

	/** HMMER 3.4 p7_tophits_TabularDomains emits sequential indices for reported
	 * domains and nreported as their count (p7_tophits.c:1619-1645). The target
	 * compactor copies tblout's rep field, so every promised row must be present,
	 * including rows below the native score gate. Zero reported domains is valid.
	 */
	private void validateDomainCounts(){
		for(Family family : families){
			for(Model model : family.models.values()){
				final int observed=model.seen.cardinality();
				require(observed==model.reportedDomains,
					"Incomplete reported domains for family "+family.id+" model "+model.accession+
					": observed="+observed+" expected="+model.reportedDomains+
					"; HMMER reports one row per declared domain");
			}
		}
	}

	private Family family(final LineParser1 lp, final String file, final long lineNo){
		final int index=familyIndex.get(lp.parseInt(0));
		require(index>=0, "Unknown family in "+file+" line "+lineNo);
		final Family family=families.get(index);
		require(family.active==lp.parseInt(1) && lp.termEquals(family.rep, 2),
			"Family identity mismatch in "+file+" line "+lineNo);
		return family;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeDomains(){
		final ByteStreamWriter writer=writer(outputPaths[0]);
		final ByteBuilder bb=new ByteBuilder(4096);
		try{
			writer.println(DOMAIN_HEADER+"\tfor_naming\tcoverage_failure");
			for(Domain domain : retained){
				bb.clear().append(domain.raw).tab().append(domain.model.naming).tab()
					.append(FAILURES[domain.flags]).nl();
				writer.print(bb);
			}
		}finally{close(writer, outputPaths[0]);}
	}

	private void writeFamilies(){
		final ByteStreamWriter writer=writer(outputPaths[1]);
		final ByteBuilder bb=new ByteBuilder(8192);
		try{
			bb.append(decisionHeader).append("\tscorepass_models\tscorepass_domains\tcoverage_pass_domains")
				.append("\tmodel_only_fail_domains\tquery_only_fail_domains\tboth_fail_domains")
				.append("\tnaming_coverage_pass_domains\tnaming_model_only_fail_domains")
				.append("\tnaming_query_only_fail_domains\tnaming_both_fail_domains\n");
			writer.print(bb);
			for(Family family : families){
				if(!family.candidate){continue;}
				int models=0, total=0;
				for(Model model : family.models.values()){if(model.passing>0){models++;}}
				for(int count : family.counts){total+=count;}
				bb.clear().append(family.raw).tab().append(models).tab().append(total);
				for(int count : family.counts){bb.tab().append(count);}
				for(int count : family.namingCounts){bb.tab().append(count);}
				writer.print(bb.nl());
			}
		}finally{close(writer, outputPaths[1]);}
	}

	private void writeSummary(){
		final ByteStreamWriter writer=writer(outputPaths[2]);
		long withDomains=0, withFailure=0;
		for(Family family : families){
			final long fail=(long)family.counts[1]+family.counts[2]+family.counts[3];
			if(fail+family.counts[0]>0){withDomains++;}
			if(fail>0){withFailure++;}
		}
		try{
			writer.println("metric\tvalue");
			writer.println("input_families\t"+families.size());
			writer.println("candidate_families\t"+candidateFamilies);
			writer.println("families_with_scorepass_domains\t"+withDomains);
			writer.println("families_with_coverage_fail_domains\t"+withFailure);
			writer.println("scorepass_domain_rows\t"+retained.size());
			writer.println("min_coverage\t"+minCoverage);
		}finally{close(writer, outputPaths[2]);}
	}

	private ByteStreamWriter writer(final String file){
		final ByteStreamWriter writer=new ByteStreamWriter(file, overwrite, false, true);
		writer.start(); return writer;
	}
	private static void close(final ByteStreamWriter writer, final String file){
		require(!writer.poisonAndWait(), "I/O error writing "+file);
	}
	private static double finite(final LineParser1 lp, final int field){
		// Parse.parseDouble's reciprocal multiplication rounds 0.800000000000
		// below 0.8. Gate comparisons need the same correctly rounded parser
		// as CLI thresholds, regardless of the input's trailing decimal zeroes.
		final double value=Double.parseDouble(lp.parseString(field));
		require(Double.isFinite(value), "Non-finite numeric field "+field+": "+lp.parseString(field));
		return value;
	}
	private static void checkHeader(final LineParser1 lp, final String expected, final String file){
		final String[] terms=expected.split("\t", -1);
		require(lp.terms()==terms.length, "Wrong header width: "+file);
		for(int i=0; i<terms.length; i++){
			require(lp.termEquals(terms[i], i), "Unexpected column "+i+" in "+file);
		}
	}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	private static final class Family {
		Family(final int id_, final int active_, final String rep_, final boolean candidate_, final byte[] raw_){
			id=id_; active=active_; rep=rep_; candidate=candidate_; raw=raw_;
		}
		final int id, active;
		final String rep;
		final boolean candidate;
		final byte[] raw;
		final LinkedHashMap<String, Model> models=new LinkedHashMap<String, Model>();
		final int[] counts=new int[4], namingCounts=new int[4];
		int length;
	}
	private static final class Model {
		Model(final String accession_, final double score_, final double cutoff_, final double domainCutoff_,
				final int hmmLength_, final String type_, final String naming_, final int reportedDomains_){
			accession=accession_; score=score_; cutoff=cutoff_; domainCutoff=domainCutoff_;
			hmmLength=hmmLength_; type=type_; naming=naming_;
			reportedDomains=reportedDomains_;
		}
		final String accession, type, naming;
		final double score, cutoff, domainCutoff;
		final int hmmLength, reportedDomains;
		final BitSet seen=new BitSet();
		int passing;
	}
	private static final class Domain {
		Domain(final Family family_, final Model model_, final int flags_, final byte[] raw_){
			family=family_; model=model_; flags=flags_; raw=raw_;
		}
		final Family family;
		final Model model;
		final int flags;
		final byte[] raw;
	}

	private String decisions, targets, domains, out;
	private String[] outputPaths;
	private byte[] decisionHeader;
	private boolean overwrite=false;
	private double minCoverage=0.8;
	private long candidateFamilies;
	private final IntHashMap2 familyIndex=new IntHashMap2(), activeIndex=new IntHashMap2();
	private final ArrayList<Family> families=new ArrayList<Family>();
	private final ArrayList<Domain> retained=new ArrayList<Domain>();
	private static final String[] FAILURES={"none", "model_only", "query_only", "both"};
	static final String TARGET_HEADER="family_id\tactive_index\tconsensus_id\tncbi_accession\thmm_name\tfull_evalue\tfull_score\tbest_domain_evalue\tbest_domain_score\treported_domains\tsequence_cutoff\tdomain_cutoff\thmm_length\tfamily_type\tfor_naming\tsource\tgene_symbol\tec_numbers\tproduct_name";
	static final String DOMAIN_HEADER="family_id\tactive_index\tconsensus_id\tconsensus_len\tncbi_accession\thmm_name\thmm_len\tdomain_index\tdomain_count\tfull_evalue\tfull_score\tconditional_evalue\tindependent_evalue\tdomain_score\thmm_from\thmm_to\tali_from\tali_to\tenv_from\tenv_to\thmm_cov\tconsensus_cov\tsequence_cutoff\tdomain_cutoff\tfamily_type\tsource\tgene_symbol\tproduct_name";
}
