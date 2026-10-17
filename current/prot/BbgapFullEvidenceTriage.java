package prot;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Rolls compact full-cohort BBGAP public evidence into one row per family.
 * @author Yelan
 */
public final class BbgapFullEvidenceTriage {

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BbgapFullEvidenceTriage.class, false);
		try{
			new BbgapFullEvidenceTriage(pp.args).process();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	BbgapFullEvidenceTriage(final String[] args){
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("families")){familiesFile=b;}
			else if(a.equals("consensus")){consensusFile=b;}
			else if(a.equals("targets")){targetsFile=b;}
			else if(a.equals("domains")){domainsFile=b;}
			else if(a.equals("swiss") || a.equals("swissprot")){swissFile=b;}
			else if(a.equals("out")){outPrefix=b;}
			else if(a.equals("toprawaudit") || a.equals("audit")){topRawAuditFile=b;}
			else if(a.equals("mincov")){minCoverage=Double.parseDouble(b);}
			else if(a.equals("expectedfamilies")){expectedFamilies=Parse.parseInt(b, 0);}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		require(familiesFile!=null && consensusFile!=null && targetsFile!=null &&
			domainsFile!=null && swissFile!=null && outPrefix!=null,
			"Required: families= consensus= targets= domains= swiss= out=");
		require(minCoverage>0 && minCoverage<=1, "mincov must be in (0, 1]: "+minCoverage);
	}

	void process(){
		loadFamilies();
		loadConsensus();
		loadTargets();
		loadDomains();
		loadSwiss();
		choose();
		writeModelCandidates();
		writeReview();
		writeSummary();
		if(topRawAuditFile!=null){writeTopRawAudit();}
	}

	private void loadFamilies(){
		final ByteFile bf=ByteFile.makeByteFile(familiesFile, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0){continue;}
				lp.set(line);
				if(lineNumber==1){
					require(lp.parseString(0).equals("active_index") && lp.parseString(3).equals("rep_id"),
						"Unexpected family_identity header in "+familiesFile);
					continue;
				}
				require(lp.terms()>=6, "Too few family_identity columns in "+familiesFile+" line "+lineNumber);
				final int active=lp.parseInt(0);
				final int familyId=lp.parseInt(1);
				final Family family=new Family(familyId, active, lp.parseString(3),
					lp.parseLong(4), lp.parseLong(5));
				require(families.put(familyId, family)==null, "Duplicate family_id "+familyId);
				ordered.add(family);
			}
		}finally{
			require(!bf.close(), "I/O error reading "+familiesFile);
		}
		require(ordered.size()==expectedFamilies, "Expected "+expectedFamilies+" families; found "+ordered.size());
	}

	private void loadConsensus(){
		final ByteFile bf=ByteFile.makeByteFile(consensusFile, true);
		Family current=null;
		int len=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='>'){
					if(current!=null){finishConsensus(current, len);}
					final String header=new String(line, 1, line.length-1);
					final String rep=firstToken(header);
					final int familyId=parseDescInt(header, "family_id");
					final int active=parseDescInt(header, "active_index");
					current=families.get(familyId);
					require(current!=null, "Consensus for unknown family_id "+familyId+": "+rep);
					require(current.activeIndex==active, "Active-index mismatch for family "+familyId);
					require(current.repId.equals(rep), "Representative mismatch for family "+familyId);
					len=0;
				}else{
					require(current!=null, "Sequence before first consensus header");
					len+=line.length;
				}
			}
		}finally{
			require(!bf.close(), "I/O error reading "+consensusFile);
		}
		if(current!=null){finishConsensus(current, len);}
		for(Family family : ordered){
			require(family.consensusAa>0, "Missing consensus for family_id "+family.familyId);
		}
	}

	private void finishConsensus(final Family family, final int len){
		require(family.consensusAa==0, "Duplicate consensus for family_id "+family.familyId);
		require(len>0, "Empty consensus for family_id "+family.familyId);
		family.consensusAa=len;
	}

	private void loadTargets(){
		final ByteFile bf=ByteFile.makeByteFile(targetsFile, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				lp.set(line);
				if(lineNumber==1){
					require(lp.terms()==19, "Unexpected target header width in "+targetsFile);
					continue;
				}
				require(lp.terms()==19, "Unexpected target width in "+targetsFile+" line "+lineNumber);
				final Family family=requiredFamily(lp, targetsFile, lineNumber);
				require(family.repId.equals(lp.parseString(2)),
					"Consensus ID mismatch in "+targetsFile+" line "+lineNumber);
				final String accession=lp.parseString(3);
				final ModelEvidence model=model(family, accession);
				model.hmmName=lp.parseString(4);
				model.fullScore=decimal(lp, 6);
				model.bestReportedDomainScore=decimal(lp, 8);
				model.reportedDomains=lp.parseInt(9);
				model.sequenceCutoff=decimal(lp, 10);
				model.hmmLength=lp.parseInt(12);
				model.familyType=lp.parseString(13);
				model.forNaming=lp.parseString(14);
				model.source=lp.parseString(15);
				model.gene=lp.parseString(16);
				model.product=lp.parseString(18);
				model.hasTarget=true;
			}
		}finally{
			require(!bf.close(), "I/O error reading "+targetsFile);
		}
	}

	private void loadDomains(){
		final ByteFile bf=ByteFile.makeByteFile(domainsFile, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				lp.set(line);
				if(lineNumber==1){
					require(lp.terms()==28, "Unexpected domain header width in "+domainsFile);
					continue;
				}
				require(lp.terms()==28, "Unexpected domain width in "+domainsFile+" line "+lineNumber);
				final Family family=requiredFamily(lp, domainsFile, lineNumber);
				require(family.repId.equals(lp.parseString(2)),
					"Consensus ID mismatch in "+domainsFile+" line "+lineNumber);
				require(family.consensusAa==lp.parseInt(3),
					"Consensus length mismatch in "+domainsFile+" line "+lineNumber);
				final ModelEvidence model=model(family, lp.parseString(4));
				require(model.hasTarget, "Domain row without target row for family_id "
					+family.familyId+" accession "+model.accession+" in "+domainsFile
					+" line "+lineNumber);
				final DomainEvidence domain=new DomainEvidence();
				domain.hmmName=lp.parseString(5);
				domain.score=decimal(lp, 13);
				domain.hmmCov=decimal(lp, 20);
				domain.consensusCov=decimal(lp, 21);
				domain.cutoff=decimal(lp, 23);
				domain.source=lp.parseString(25);
				domain.gene=lp.parseString(26);
				domain.product=lp.parseString(27);
				model.bestDomain=betterDomain(domain, model.bestDomain);
				if(domain.eligible(minCoverage)){
					model.bestEligibleDomain=betterDomain(domain, model.bestEligibleDomain);
				}
				if(model.bestDomain==domain){
					model.hmmName=domain.hmmName;
					model.source=domain.source;
					model.gene=domain.gene;
					model.product=domain.product;
					model.hasDomain=true;
				}
			}
		}finally{
			require(!bf.close(), "I/O error reading "+domainsFile);
		}
	}

	private void loadSwiss(){
		final ByteFile bf=ByteFile.makeByteFile(swissFile, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				lp.set(line);
				if(lineNumber==1){
					require(lp.terms()==18, "Unexpected Swiss-Prot header width in "+swissFile);
					continue;
				}
				require(lp.terms()==18, "Unexpected Swiss-Prot width in "+swissFile+" line "+lineNumber);
				final Family family=requiredFamily(lp, swissFile, lineNumber);
				require(family.repId.equals(lp.parseString(2)),
					"Consensus ID mismatch in "+swissFile+" line "+lineNumber);
				require(family.consensusAa==lp.parseInt(7),
					"Consensus length mismatch in "+swissFile+" line "+lineNumber);
				final SwissEvidence swiss=new SwissEvidence();
				swiss.accession=lp.parseString(3);
				swiss.entry=lp.parseString(4);
				swiss.identity=decimal(lp, 5);
				swiss.queryCov=decimal(lp, 9);
				swiss.subjectCov=decimal(lp, 10);
				swiss.bitscore=decimal(lp, 16);
				swiss.title=lp.parseString(17);
				family.swiss.add(swiss);
			}
		}finally{
			require(!bf.close(), "I/O error reading "+swissFile);
		}
	}

	private Family requiredFamily(final LineParser1 lp, final String path, final long lineNumber){
		final int familyId=lp.parseInt(0);
		final int active=lp.parseInt(1);
		final Family family=families.get(familyId);
		require(family!=null, "Unknown family_id "+familyId+" in "+path+" line "+lineNumber);
		require(family.activeIndex==active, "Active-index mismatch in "+path+" line "+lineNumber);
		return family;
	}

	private ModelEvidence model(final Family family, final String accession){
		ModelEvidence model=family.models.get(accession);
		if(model==null){
			model=new ModelEvidence();
			model.family=family;
			model.accession=accession;
			family.models.put(accession, model);
		}
		return model;
	}

	private void choose(){
		for(Family family : ordered){
			for(ModelEvidence model : family.models.values()){
				if(model.eligible(minCoverage)){
					family.eligibleModels++;
					family.bestEligible=betterNative(model, family.bestEligible);
				}
				family.topModel=betterRaw(model, family.topModel);
			}
			for(SwissEvidence swiss : family.swiss){
				if(swiss.fullLength(minCoverage)){family.fullSwiss++;}
				family.topSwiss=betterSwiss(swiss, family.topSwiss);
			}
			if(family.models.isEmpty() && family.swiss.isEmpty()){family.candidateClass="no_evidence";}
			else if(family.eligibleModels>1){family.candidateClass="native_multimodel_review";}
			else if(family.eligibleModels==1){family.candidateClass="gated_full_length";}
			else if(family.fullSwiss>0){family.candidateClass="swissprot_full_length_review";}
			else{family.candidateClass="partial_or_below_cutoff";}
		}
	}

	private void writeReview(){
		final ByteStreamWriter bsw=writer(outPrefix+".review_candidates.tsv");
		final ByteBuilder bb=new ByteBuilder(1<<16);
		try{
			bsw.println("family_id\tactive_index\trep_id\tmembers\traw_residues\tconsensus_aa\t"
				+"candidate_class\tncbi_model_count\tncbi_eligible_count\tswissprot_hit_count\t"
				+"swissprot_full_length_count\teligible_ncbi_accession\teligible_ncbi_score\t"
				+"eligible_domain_score\teligible_hmm_cov\teligible_consensus_cov\teligible_product\t"
				+"top_ncbi_accession\ttop_ncbi_score\ttop_swissprot_accession\t"
				+"top_swissprot_query_cov\ttop_swissprot_subject_cov\ttop_swissprot_title");
			for(Family family : ordered){
				bb.clear();
				family.append(bb).tab().append(family.candidateClass).tab()
					.append(family.models.size()).tab().append(family.eligibleModels).tab()
					.append(family.swiss.size()).tab().append(family.fullSwiss).tab();
				appendModel(bb, family.bestEligible);
				bb.tab();
				appendModelSummary(bb, family.topModel);
				bb.tab();
				appendSwiss(bb, family.topSwiss);
				bb.nl();
				bsw.print(bb);
			}
		}finally{
			close(bsw, outPrefix+".review_candidates.tsv");
		}
	}

	private void writeModelCandidates(){
		final ArrayList<ModelEvidence> models=new ArrayList<ModelEvidence>();
		for(Family family : ordered){models.addAll(family.models.values());}
		Collections.sort(models, MODEL_ORDER);
		final ByteStreamWriter bsw=writer(outPrefix+".ncbi_model_candidates.tsv");
		final ByteBuilder bb=new ByteBuilder(1<<16);
		try{
			bsw.println("family_id\tactive_index\trep_id\tncbi_accession\tfull_score\t"
				+"sequence_cutoff\tsequence_pass\tdomain_score\tdomain_cutoff\tdomain_pass\t"
				+"hmm_cov\tconsensus_cov\tmutual_coverage_pass\tsame_hit_native_eligible\t"
				+"hmm_name\tfamily_type\tfor_naming\tsource\tgene_symbol\tproduct_name");
			for(ModelEvidence model : models){
				bb.clear();
				final DomainEvidence domain=model.displayDomain();
				model.family.appendIds(bb).tab().append(model.accession).tab()
					.append(model.hasTarget ? model.fullScore : Double.NaN, 3, true).tab()
					.append(model.hasTarget ? model.sequenceCutoff : Double.NaN, 3, true).tab()
					.append(model.sequencePass() ? 'Y' : 'N').tab()
					.append(domain==null ? Double.NaN : domain.score, 3, true).tab()
					.append(domain==null ? Double.NaN : domain.cutoff, 3, true).tab()
					.append(domain!=null && domain.domainPass() ? 'Y' : 'N').tab();
				appendDecimal(bb, domain==null ? Double.NaN : domain.hmmCov).tab();
				appendDecimal(bb, domain==null ? Double.NaN : domain.consensusCov).tab()
					.append(domain!=null && domain.coveragePass(minCoverage) ? 'Y' : 'N').tab()
					.append(model.eligible(minCoverage) ? 'Y' : 'N').tab()
					.append(model.hmmName).tab().append(model.familyType).tab()
					.append(model.forNaming).tab().append(domain==null ? model.source : domain.source).tab()
					.append(domain==null ? model.gene : domain.gene).tab()
					.append(domain==null ? model.product : domain.product).nl();
				bsw.print(bb);
			}
		}finally{
			close(bsw, outPrefix+".ncbi_model_candidates.tsv");
		}
	}

	private void writeSummary(){
		long noEvidence=0, gated=0, multi=0, swissOnly=0, partial=0;
		for(Family family : ordered){
			if(family.candidateClass.equals("no_evidence")){noEvidence++;}
			else if(family.candidateClass.equals("gated_full_length")){gated++;}
			else if(family.candidateClass.equals("native_multimodel_review")){multi++;}
			else if(family.candidateClass.equals("swissprot_full_length_review")){swissOnly++;}
			else if(family.candidateClass.equals("partial_or_below_cutoff")){partial++;}
		}
		final ByteStreamWriter bsw=writer(outPrefix+".review_summary.tsv");
		try{
			bsw.println("metric\tvalue");
			bsw.println("total_families\t"+ordered.size());
			bsw.println("gated_full_length\t"+gated);
			bsw.println("native_multimodel_review\t"+multi);
			bsw.println("swissprot_full_length_review\t"+swissOnly);
			bsw.println("partial_or_below_cutoff\t"+partial);
			bsw.println("no_evidence\t"+noEvidence);
		}finally{
			close(bsw, outPrefix+".review_summary.tsv");
		}
		require(noEvidence+gated+multi+swissOnly+partial==ordered.size(),
			"Summary categories do not conserve families");
	}

	private void writeTopRawAudit(){
		final ByteStreamWriter bsw=writer(topRawAuditFile);
		final ByteBuilder bb=new ByteBuilder(1<<16);
		try{
			bsw.println("family_id\tactive_index\ttop_ncbi_accession\ttop_full_score\t"
				+"top_sequence_cutoff\ttop_sequence_pass\ttop_same_hit_native_eligible\t"
				+"eligible_ncbi_count\teligible_ncbi_accessions\teligible_products");
			for(Family family : ordered){
				if(family.eligibleModels<1 || family.topModel==null ||
						family.topModel.eligible(minCoverage)){continue;}
				final ArrayList<ModelEvidence> models=eligibleModels(family);
				bb.clear();
				bb.append(family.familyId).tab().append(family.activeIndex).tab()
					.append(family.topModel.accession).tab()
					.append(family.topModel.fullScore, 3, true).tab()
					.append(family.topModel.sequenceCutoff, 3, true).tab()
					.append(family.topModel.sequencePass() ? 'Y' : 'N').tab()
					.append(family.topModel.eligible(minCoverage) ? 'Y' : 'N').tab()
					.append(models.size()).tab();
				appendAccessions(bb, models).tab();
				appendProducts(bb, models).nl();
				bsw.print(bb);
			}
		}finally{
			close(bsw, topRawAuditFile);
		}
	}

	private static ArrayList<ModelEvidence> eligibleModels(final Family family){
		final ArrayList<ModelEvidence> list=new ArrayList<ModelEvidence>();
		for(ModelEvidence model : family.models.values()){
			if(model.eligible(minCoverage)){list.add(model);}
		}
		Collections.sort(list, MODEL_ORDER);
		require(list.size()==family.eligibleModels, "Eligible-count mismatch for family_id "
			+family.familyId);
		return list;
	}

	private static ByteBuilder appendAccessions(final ByteBuilder bb,
			final ArrayList<ModelEvidence> models){
		for(int i=0; i<models.size(); i++){
			if(i>0){bb.append(',');}
			bb.append(models.get(i).accession);
		}
		return bb;
	}

	private static ByteBuilder appendProducts(final ByteBuilder bb,
			final ArrayList<ModelEvidence> models){
		for(int i=0; i<models.size(); i++){
			if(i>0){bb.append('|');}
			bb.append(models.get(i).product);
		}
		return bb;
	}

	private static ModelEvidence betterNative(final ModelEvidence a, final ModelEvidence b){
		if(b==null){return a;}
		if(a.fullScore!=b.fullScore){return a.fullScore>b.fullScore ? a : b;}
		if(a.eligibleDomainScore()!=b.eligibleDomainScore()){
			return a.eligibleDomainScore()>b.eligibleDomainScore() ? a : b;
		}
		return a.accession.compareTo(b.accession)<0 ? a : b;
	}

	private static ModelEvidence betterRaw(final ModelEvidence a, final ModelEvidence b){
		if(b==null){return a;}
		final double as=a.hasTarget ? a.fullScore : a.domainScore();
		final double bs=b.hasTarget ? b.fullScore : b.domainScore();
		if(as!=bs){return as>bs ? a : b;}
		return a.accession.compareTo(b.accession)<0 ? a : b;
	}

	private static SwissEvidence betterSwiss(final SwissEvidence a, final SwissEvidence b){
		if(b==null){return a;}
		if(a.bitscore!=b.bitscore){return a.bitscore>b.bitscore ? a : b;}
		return a.accession.compareTo(b.accession)<0 ? a : b;
	}

	private static DomainEvidence betterDomain(final DomainEvidence a, final DomainEvidence b){
		if(b==null){return a;}
		return a.score>b.score ? a : b;
	}

	private static void appendModel(final ByteBuilder bb, final ModelEvidence model){
		if(model==null){bb.append("\t\t\t\t\t"); return;}
		final DomainEvidence domain=model.displayDomain();
		bb.append(model.accession).tab().append(model.fullScore, 3, true).tab()
			.append(domain.score, 3, true).tab();
		appendDecimal(bb, domain.hmmCov).tab();
		appendDecimal(bb, domain.consensusCov).tab().append(domain.product);
	}

	private static void appendModelSummary(final ByteBuilder bb, final ModelEvidence model){
		if(model==null){bb.append('\t'); return;}
		bb.append(model.accession).tab().append(model.hasTarget ? model.fullScore : model.domainScore(), 3, true);
	}

	private static void appendSwiss(final ByteBuilder bb, final SwissEvidence swiss){
		if(swiss==null){bb.append("\t\t\t"); return;}
		bb.append(swiss.accession).tab();
		appendDecimal(bb, swiss.queryCov).tab();
		appendDecimal(bb, swiss.subjectCov).tab().append(swiss.title);
	}

	private static ByteBuilder appendDecimal(final ByteBuilder bb, final double value){
		return Double.isNaN(value) ? bb : bb.append(String.format(java.util.Locale.ROOT, "%.12f", value));
	}

	/** Uses the same correctly rounded decimal parsing as CLI thresholds.
	 * The fast reciprocal-based parser places 0.800000000000 below 0.8.
	 */
	private static double decimal(final LineParser1 lp, final int field){
		final double value=Double.parseDouble(lp.parseString(field));
		require(Double.isFinite(value), "Non-finite evidence field "+field+": "+lp.parseString(field));
		return value;
	}

	private static String firstToken(final String s){
		int end=0;
		while(end<s.length() && s.charAt(end)>32){end++;}
		return s.substring(0, end);
	}

	private static int parseDescInt(final String desc, final String key){
		final String prefix=key+"=";
		int start=0;
		while(start<desc.length()){
			while(start<desc.length() && desc.charAt(start)<=' '){start++;}
			int stop=start;
			while(stop<desc.length() && desc.charAt(stop)>' '){stop++;}
			if(desc.startsWith(prefix, start)){return Parse.parseInt(desc.substring(start+prefix.length(), stop), 0);}
			start=stop+1;
		}
		throw new IllegalArgumentException("Missing "+key+" in "+desc);
	}

	private static ByteStreamWriter writer(final String path){
		final ByteStreamWriter bsw=new ByteStreamWriter(path, true, false, true);
		bsw.start();
		return bsw;
	}

	private static void close(final ByteStreamWriter bsw, final String out){
		if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+out);}
	}

	private static void require(final boolean ok, final String message){
		if(!ok){throw new IllegalArgumentException(message);}
	}

	private static final Comparator<ModelEvidence> MODEL_ORDER=new Comparator<ModelEvidence>() {
		@Override
		public int compare(final ModelEvidence a, final ModelEvidence b){
			if(a.family.activeIndex!=b.family.activeIndex){return a.family.activeIndex-b.family.activeIndex;}
			if(a.eligible(minCoverage)!=b.eligible(minCoverage)){return a.eligible(minCoverage) ? -1 : 1;}
			final int score=Double.compare(b.hasTarget ? b.fullScore : b.domainScore(),
				a.hasTarget ? a.fullScore : a.domainScore());
			return score!=0 ? score : a.accession.compareTo(b.accession);
		}
	};

	private static final class Family {
		Family(final int familyId_, final int activeIndex_, final String repId_,
				final long members_, final long rawResidues_){
			familyId=familyId_;
			activeIndex=activeIndex_;
			repId=repId_;
			members=members_;
			rawResidues=rawResidues_;
		}
		ByteBuilder appendIds(final ByteBuilder bb){
			return bb.append(familyId).tab().append(activeIndex).tab().append(repId);
		}
		ByteBuilder append(final ByteBuilder bb){
			return appendIds(bb).tab().append(members).tab().append(rawResidues).tab().append(consensusAa);
		}
		final int familyId, activeIndex;
		final long members, rawResidues;
		final String repId;
		int consensusAa, eligibleModels, fullSwiss;
		String candidateClass;
		ModelEvidence bestEligible, topModel;
		SwissEvidence topSwiss;
		final HashMap<String, ModelEvidence> models=new HashMap<String, ModelEvidence>();
		final ArrayList<SwissEvidence> swiss=new ArrayList<SwissEvidence>();
	}

	private static final class ModelEvidence {
		boolean sequencePass(){return hasTarget && fullScore>=sequenceCutoff;}
		boolean domainPass(){return bestEligibleDomain!=null;}
		boolean coveragePass(final double minCoverage){
			return bestEligibleDomain!=null;
		}
		boolean eligible(final double minCoverage){
			return sequencePass() && domainPass() && coveragePass(minCoverage);
		}
		DomainEvidence displayDomain(){
			return bestEligibleDomain==null ? bestDomain : bestEligibleDomain;
		}
		double eligibleDomainScore(){
			return bestEligibleDomain==null ? Double.NaN : bestEligibleDomain.score;
		}
		double domainScore(){
			return bestDomain==null ? Double.NaN : bestDomain.score;
		}
		Family family;
		String accession="", hmmName="", familyType="", forNaming="", source="", gene="", product="";
		double fullScore, bestReportedDomainScore, sequenceCutoff;
		DomainEvidence bestDomain, bestEligibleDomain;
		int reportedDomains, hmmLength;
		boolean hasTarget, hasDomain;
	}

	private static final class DomainEvidence {
		boolean domainPass(){return score>=cutoff;}
		boolean coveragePass(final double minCoverage){
			return hmmCov>=minCoverage && consensusCov>=minCoverage;
		}
		boolean eligible(final double minCoverage){return domainPass() && coveragePass(minCoverage);}
		String hmmName="", source="", gene="", product="";
		double score, cutoff, hmmCov, consensusCov;
	}

	private static final class SwissEvidence {
		boolean fullLength(final double minCoverage){return queryCov>=minCoverage && subjectCov>=minCoverage;}
		String accession, entry, title;
		double identity, queryCov, subjectCov, bitscore;
	}

	private final LinkedHashMap<Integer, Family> families=new LinkedHashMap<Integer, Family>();
	private final ArrayList<Family> ordered=new ArrayList<Family>(8537);
	private String familiesFile, consensusFile, targetsFile, domainsFile, swissFile, outPrefix;
	private String topRawAuditFile;
	private int expectedFamilies=8537;
	private static double minCoverage=0.80;
}
