package prot;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/** Assigns BBGAP protein-family names from NCBI/Swiss-Prot agreement and reviewed overrides.
 * @author Yelan
 */
public final class AssignProteinNames {

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, AssignProteinNames.class, false);
		try{
			new AssignProteinNames(pp.args).process();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	AssignProteinNames(final String[] args){
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("families")){familiesFile=b;}
			else if(a.equals("ncbi")){ncbiFile=b;}
			else if(a.equals("ncbinearquerydomains")){ncbiNearQueryDomainsFile=b;}
			else if(a.equals("swiss") || a.equals("swissprot")){swissFile=b;}
			else if(a.equals("reviewed")){reviewedFile=b;}
			else if(a.equals("reviewedevidence")){reviewedEvidenceFile=b;}
			else if(a.equals("reviewedevidenceout")){reviewedEvidenceOutFile=b;}
			else if(a.equals("eggnoghomology")){eggnogHomologyFile=b;}
			else if(a.equals("eggnoghomologyout")){eggnogHomologyOutFile=b;}
			else if(a.equals("eggnogdomainrelaxed")){
				eggnogDomainRelaxed=Parse.parseBoolean(b);
			}
			else if(a.equals("synonyms")){synonymsFile=b;}
			else if(a.equals("broadconflicts")){broadConflictRulesFile=b;}
			else if(a.equals("broadparents")){broadParentRulesFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("evidence")){evidenceFile=b;}
			else if(a.equals("unresolved")){unresolvedFile=b;}
			else if(a.equals("summary")){summaryFile=b;}
			else if(a.equals("audit")){auditFile=b;}
			else if(a.equals("conflictgroups")){conflictGroupsFile=b;}
			else if(a.equals("membercoverage")){memberCoverageFile=b;}
			else if(a.equals("memberunresolved")){memberUnresolvedFile=b;}
			else if(a.equals("labelversion")){labelVersion=b;}
			else if(a.equals("identitypin")){identityPin=b;}
			else if(a.equals("consensuspin")){consensusPin=b;}
			else if(a.equals("mincov")){minCoverage=Double.parseDouble(b);}
			else if(a.equals("minswissidentity")){minSwissIdentity=Double.parseDouble(b);}
			else if(a.equals("swissonly") || a.equals("swissonlyfallback")){
				swissOnlyFallback=Parse.parseBoolean(b);
			}
			else if(a.equals("ncbinearquery") || a.equals("ncbinearqueryfallback")){
				ncbiNearQueryFallback=Parse.parseBoolean(b);
			}
			else if(a.equals("provisional") || a.equals("provisionallabels")){
				provisionalLabels=Parse.parseBoolean(b);
			}
			else if(a.equals("maxaudit")){maxAudit=Parse.parseInt(b, 0);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("ignorecase")){ignoreCase=Parse.parseBoolean(b);}
			else if(a.equals("mode") || a.equals("rule")){mode=b.toLowerCase();}
			else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		require(familiesFile!=null && ncbiFile!=null && swissFile!=null && outFile!=null
			&& evidenceFile!=null && unresolvedFile!=null && summaryFile!=null,
			"Required: families= ncbi= swiss= out= evidence= unresolved= summary=");
		require(labelVersion!=null && labelVersion.matches("[A-Za-z0-9._-]+"),
			"Required: labelversion=ASCII_TOKEN");
		require(identityPin!=null && identityPin.length()>0, "Required: identitypin=");
		require(consensusPin!=null && consensusPin.length()>0, "Required: consensuspin=");
		require(mode.equals("exact") || mode.equals("hierarchy"),
			"mode must be exact or hierarchy: "+mode);
		require(conflictGroupsFile==null || mode.equals("hierarchy"),
			"conflictgroups requires mode=hierarchy");
		require(!swissOnlyFallback || mode.equals("hierarchy"),
			"swissonly requires mode=hierarchy");
		require(!ncbiNearQueryFallback || mode.equals("hierarchy"),
			"ncbinearquery requires mode=hierarchy");
		require(!ncbiNearQueryFallback || ncbiNearQueryDomainsFile!=null,
			"ncbinearquery requires ncbinearquerydomains=");
		require(broadConflictRulesFile==null || mode.equals("hierarchy"),
			"broadconflicts requires mode=hierarchy");
		require(broadParentRulesFile==null || mode.equals("hierarchy"),
			"broadparents requires mode=hierarchy");
		require(minCoverage>0 && minCoverage<=1, "mincov must be in (0, 1]: "+minCoverage);
		require(minSwissIdentity>=0 && minSwissIdentity<=100,
			"minswissidentity must be in [0, 100]: "+minSwissIdentity);
		require(maxAudit>=0, "maxaudit must be nonnegative: "+maxAudit);
		require(reviewedEvidenceFile==null || reviewedFile!=null,
			"reviewedevidence requires reviewed=");
		require((reviewedEvidenceFile==null)==(reviewedEvidenceOutFile==null),
			"reviewedevidence and reviewedevidenceout must be supplied together");
		require((eggnogHomologyFile==null)==(eggnogHomologyOutFile==null),
			"eggnoghomology and eggnoghomologyout must be supplied together");
		require(eggnogHomologyFile==null || provisionalLabels,
			"eggnoghomology requires provisional=t");
		checkFiles();
	}

	void process(){
		loadFamilies();
		if(synonymsFile!=null){loadSynonyms();}
		if(broadConflictRulesFile!=null){loadBroadConflictRules();}
		if(broadParentRulesFile!=null){loadBroadParentRules();}
		loadNcbi();
		if(ncbiNearQueryFallback){loadNcbiNearQueryDomains();}
		loadSwiss();
		if(reviewedFile!=null){loadReviewed();}
		if(reviewedEvidenceFile!=null){loadReviewedEvidence();}
		if(eggnogHomologyFile!=null){loadEggnogHomology();}
		choose();
		if(provisionalLabels){chooseProvisionalLabels();}
		if(eggnogHomologyFile!=null){chooseEggnogHomologyLabels();}
		writeLabels();
		writeEvidence();
		writeUnresolved();
		writeSummary();
		if(auditFile!=null){writeAudit();}
		if(conflictGroupsFile!=null){writeConflictGroups();}
		if(memberCoverageFile!=null){writeMemberCoverage();}
		if(memberUnresolvedFile!=null){writeMemberUnresolved();}
		if(reviewedEvidenceOutFile!=null){writeReviewedEvidence();}
		if(eggnogHomologyOutFile!=null){writeEggnogHomologyEvidence();}
	}

	private void loadFamilies(){
		final ByteFile bf=ByteFile.makeByteFile(familiesFile, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		int membersColumn=-1, rawResiduesColumn=-1;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0){continue;}
				lp.set(line);
				if(lineNumber==1){
					require(lp.parseString(0).equals("active_index") && lp.parseString(3).equals("rep_id"),
						"Unexpected family_identity header in "+familiesFile);
					if(memberReportsRequested()){
						membersColumn=headerIndex(lp, "members", familiesFile);
						rawResiduesColumn=headerIndex(lp, "raw_residues", familiesFile);
						require(membersColumn>=0 && rawResiduesColumn>=0,
							"family_identity requires members and raw_residues for member reports");
					}
					continue;
				}
				final int minColumns=memberReportsRequested() ?
					Math.max(membersColumn, rawResiduesColumn)+1 : 4;
				require(lp.terms()>=minColumns, "Too few family_identity columns in "+familiesFile
					+" line "+lineNumber);
				final Family family=new Family(lp.parseInt(1), lp.parseInt(0), lp.parseString(3));
				if(memberReportsRequested()){
					family.members=parsePositiveLong(lp, membersColumn, familiesFile, lineNumber,
						"members");
					family.rawResidues=parseNonnegativeLong(lp, rawResiduesColumn, familiesFile,
						lineNumber, "raw_residues");
				}
				require(families.put(family.familyId, family)==null,
					"Duplicate family_id "+family.familyId+" in "+familiesFile);
				ordered.add(family);
			}
		}finally{
			require(!bf.close(), "I/O error reading "+familiesFile);
		}
		require(!ordered.isEmpty(), "No families loaded from "+familiesFile);
	}

	private void loadSynonyms(){
		final ByteFile bf=ByteFile.makeByteFile(synonymsFile, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lineNumber==1 && lp.parseString(0).equals("from_product")){continue;}
				require(lp.terms()>=2, "Too few synonym columns in "+synonymsFile
					+" line "+lineNumber);
				final String from=comparisonKey(lp.parseString(0));
				final String to=lp.parseString(1);
				require(from.length()>0 && to.length()>0, "Empty synonym in "
					+synonymsFile+" line "+lineNumber);
				require(synonyms.put(from, to)==null, "Duplicate synonym source in "
					+synonymsFile+" line "+lineNumber+": "+from);
			}
		}finally{
			require(!bf.close(), "I/O error reading "+synonymsFile);
		}
	}

	private void loadBroadConflictRules(){
		final ByteFile bf=ByteFile.makeByteFile(broadConflictRulesFile, true);
		final LineParser1 lp=new LineParser1('\t');
		final LinkedHashMap<String, CuratedBroadConflict> byRule=
			new LinkedHashMap<String, CuratedBroadConflict>();
		long lineNumber=0;
		boolean sawHeader=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(!sawHeader){
					sawHeader=true;
					require((lp.terms()==5 || lp.terms()==6)
						&& lp.parseString(0).equals("rule_id")
						&& lp.parseString(1).equals("conflict_source")
						&& lp.parseString(2).equals("decision_class")
						&& lp.parseString(3).equals("product")
						&& lp.parseString(4).equals("conflict_product")
						&& (lp.terms()==5 || lp.parseString(5).equals("evidence")),
						"Unexpected broad-conflict header in "+broadConflictRulesFile);
					continue;
				}
				require(lp.terms()==5 || lp.terms()==6, "Unexpected broad-conflict width in "
					+broadConflictRulesFile+" line "+lineNumber);
				final String ruleId=lp.parseString(0);
				final String source=lp.parseString(1);
				final String decision=lp.parseString(2);
				final String product=lp.parseString(3);
				final String conflictProduct=lp.parseString(4);
				require(ruleId.matches("[A-Za-z0-9._-]+"), "Bad broad-conflict rule_id "
					+ruleId+" line "+lineNumber);
				require(source.equals("ncbi_top_tier") || source.equals("swissprot"),
					"Bad broad-conflict source "+source+" line "+lineNumber);
				require(decision.equals("domain_supported_broad_name") ||
					decision.equals("homology_supported_broad_name"),
					"Bad broad-conflict decision "+decision+" line "+lineNumber);
				require(product.length()>0 && conflictProduct.length()>0,
					"Empty broad-conflict product line "+lineNumber);
				CuratedBroadConflict rule=byRule.get(ruleId);
				if(rule==null){
					rule=new CuratedBroadConflict(ruleId, source, decision, product);
					byRule.put(ruleId, rule);
				}
				rule.add(source, decision, product, conflictProduct, lineNumber);
			}
		}finally{
			require(!bf.close(), "I/O error reading "+broadConflictRulesFile);
		}
		require(sawHeader, "Missing broad-conflict header in "+broadConflictRulesFile);
		for(CuratedBroadConflict rule : byRule.values()){
			require(rule.keys.size()>1, "Broad-conflict rule has fewer than two products: "
				+rule.ruleId);
			final CuratedBroadConflict old=curatedBroadConflicts.put(rule.mapKey(), rule);
			if(old!=null){
				require(false, "Duplicate broad-conflict product set: "+rule.ruleId
					+" and "+old.ruleId);
			}
		}
	}

	private void loadBroadParentRules(){
		final ByteFile bf=ByteFile.makeByteFile(broadParentRulesFile, true);
		final LineParser1 lp=new LineParser1('\t');
		final LinkedHashMap<String, CuratedBroadParent> byId=
			new LinkedHashMap<String, CuratedBroadParent>();
		long lineNumber=0;
		boolean sawHeader=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(!sawHeader){
					sawHeader=true;
					require(lp.terms()==6
						&& lp.parseString(0).equals("parent_id")
						&& lp.parseString(1).equals("conflict_source")
						&& lp.parseString(2).equals("decision_class")
						&& lp.parseString(3).equals("parent_product")
						&& lp.parseString(4).equals("child_product")
						&& lp.parseString(5).equals("evidence"),
						"Unexpected broad-parent header in "+broadParentRulesFile);
					continue;
				}
				require(lp.terms()==6, "Unexpected broad-parent width in "
					+broadParentRulesFile+" line "+lineNumber);
				final String parentId=lp.parseString(0);
				final String source=lp.parseString(1);
				final String decision=lp.parseString(2);
				final String product=lp.parseString(3);
				final String childProduct=lp.parseString(4);
				require(parentId.matches("[A-Za-z0-9._-]+"), "Bad broad-parent parent_id "
					+parentId+" line "+lineNumber);
				require(source.equals("ncbi_top_tier") || source.equals("swissprot"),
					"Bad broad-parent source "+source+" line "+lineNumber);
				require(decision.equals(source.equals("ncbi_top_tier") ?
					"domain_supported_broad_name" : "homology_supported_broad_name"),
					"Bad broad-parent decision "+decision+" line "+lineNumber);
				require(product.length()>0 && childProduct.length()>0,
					"Empty broad-parent product line "+lineNumber);
				CuratedBroadParent parent=byId.get(parentId);
				if(parent==null){
					parent=new CuratedBroadParent(parentId, source, decision, product);
					byId.put(parentId, parent);
				}
				parent.add(source, decision, product, childProduct, lineNumber);
			}
		}finally{
			require(!bf.close(), "I/O error reading "+broadParentRulesFile);
		}
		require(sawHeader, "Missing broad-parent header in "+broadParentRulesFile);
		for(CuratedBroadParent parent : byId.values()){
			require(parent.childCount>=2, "Broad-parent rule has fewer than two products: "
				+parent.parentId);
		}
	}

	private void loadNcbi(){
		final ByteFile bf=ByteFile.makeByteFile(ncbiFile, true);
		final LineParser1 lp=new LineParser1('\t');
		final boolean retainRawNcbi=provisionalLabels || reviewedEvidenceFile!=null;
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				lp.set(line);
				if(lineNumber==1){
					require(lp.terms()==20 && lp.parseString(0).equals("family_id")
						&& lp.parseString(19).equals("product_name"),
						"Unexpected NCBI candidate header in "+ncbiFile);
					continue;
				}
				require(lp.terms()==20, "Unexpected NCBI candidate width in "
					+ncbiFile+" line "+lineNumber);
				final int familyId=lp.parseInt(0);
				final Family family=requiredFamily(familyId, lp.parseInt(1), lp.parseString(2),
					ncbiFile, lineNumber);
				if(retainRawNcbi){
					family.addRawNcbi(lp.parseString(3), lp.parseString(19),
						lp.parseString(15), lp.parseString(4), lp.parseString(7),
						lp.parseString(11), lp.parseString(10), lp.termEquals("Y", 6),
						lp.termEquals("Y", 9), lp.termEquals("Y", 12),
						lp.termEquals("Y", 13), lp.termEquals("Y", 16));
				}
				if(!lp.termEquals("Y", 13)){continue;}
				if(lp.termEquals("Y", 16)){
					family.addNcbi(lp.parseString(3), lp.parseString(19),
						mode.equals("hierarchy") ? lp.parseString(15) : null);
				}
			}
		}finally{
			require(!bf.close(), "I/O error reading "+ncbiFile);
		}
	}

	private void loadNcbiNearQueryDomains(){
		final ByteFile bf=ByteFile.makeByteFile(ncbiNearQueryDomainsFile, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				lp.set(line);
				if(lineNumber==1){
					require(Tools.equals(NCBI_NEAR_QUERY_DOMAIN_HEADER, line),
						"Unexpected NCBI near-query domain header in "
						+ncbiNearQueryDomainsFile);
					continue;
				}
				require(lp.terms()==30, "Unexpected NCBI near-query domain width in "
					+ncbiNearQueryDomainsFile+" line "+lineNumber);
				if(!nearQueryEligible(lp, ncbiNearQueryDomainsFile, lineNumber)){continue;}
				final int familyId=lp.parseInt(0);
				final Family family=requiredFamily(familyId, lp.parseInt(1), lp.parseString(2),
					ncbiNearQueryDomainsFile, lineNumber);
				family.addNcbi(lp.parseString(4), lp.parseString(27), lp.parseString(24), true);
			}
		}finally{
			require(!bf.close(), "I/O error reading "+ncbiNearQueryDomainsFile);
		}
	}

	private boolean nearQueryEligible(final LineParser1 lp, final String path,
			final long lineNumber){
		if(!lp.termEquals("Y", 28)){return false;}
		final String familyType=lp.parseString(24);
		if(!familyType.equals("equivalog") && !familyType.equals("exception")){return false;}
		final double fullScore=parseFinite(lp.parseString(10), path, lineNumber, "full_score");
		final double sequenceCutoff=parseFinite(lp.parseString(22), path, lineNumber,
			"sequence_cutoff");
		final double domainScore=parseFinite(lp.parseString(13), path, lineNumber,
			"domain_score");
		final double domainCutoff=parseFinite(lp.parseString(23), path, lineNumber,
			"domain_cutoff");
		final double hmmCov=parseCoverage(lp.parseString(20), path, lineNumber, "hmm_cov");
		final double consensusCov=parseCoverage(lp.parseString(21), path, lineNumber,
			"consensus_cov");
		require(fullScore>=sequenceCutoff, "NCBI near-query sequence_pass below cutoff in "
			+path+" line "+lineNumber);
		require(domainScore>=domainCutoff, "NCBI near-query domain_pass below cutoff in "
			+path+" line "+lineNumber);
		return hmmCov>=minCoverage && consensusCov>=NEAR_QUERY_MIN_COVERAGE
			&& consensusCov<minCoverage && lp.termEquals("query_only", 29);
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
					require(lp.terms()==18 && lp.parseString(0).equals("family_id")
						&& lp.parseString(17).equals("title"),
						"Unexpected Swiss-Prot compact header in "+swissFile);
					continue;
				}
				require(lp.terms()==18, "Unexpected Swiss-Prot compact width in "
					+swissFile+" line "+lineNumber);
				final int familyId=lp.parseInt(0);
				final Family family=requiredFamily(familyId, lp.parseInt(1), lp.parseString(2),
					swissFile, lineNumber);
				final String accession=lp.parseString(3);
				final String identity=lp.parseString(5);
				final String queryCov=lp.parseString(9);
				final String subjectCov=lp.parseString(10);
				final String bitscore=lp.parseString(16);
				final String product=swissProduct(lp.parseString(17), swissFile, lineNumber);
				if(provisionalLabels){
					family.addRawSwiss(accession, product, identity, queryCov, subjectCov, bitscore);
				}
				if(Double.parseDouble(queryCov)<minCoverage ||
						Double.parseDouble(subjectCov)<minCoverage){continue;}
				family.addSwiss(accession, product, identity, queryCov, subjectCov, bitscore);
			}
		}finally{
			require(!bf.close(), "I/O error reading "+swissFile);
		}
	}

	private void loadReviewed(){
		final ByteFile bf=ByteFile.makeByteFile(reviewedFile, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0){continue;}
				lp.set(line);
				if(lineNumber==1){
					require(lp.terms()==10 && lp.parseString(0).equals("family_id")
						&& lp.parseString(5).equals("provisional_product"),
						"Unexpected reviewed header in "+reviewedFile);
					continue;
				}
				require(lp.terms()==10, "Unexpected reviewed width in "+reviewedFile
					+" line "+lineNumber);
				final Family family=requiredFamily(lp.parseInt(0), lp.parseInt(1), null,
					reviewedFile, lineNumber);
				require(!family.reviewed, "Duplicate reviewed family_id "+family.familyId);
				family.reviewed=true;
				family.decision=lp.parseString(4);
				family.product=blankDot(lp.parseString(5));
				family.rule="reviewed_override";
			}
		}finally{
			require(!bf.close(), "I/O error reading "+reviewedFile);
		}
	}

	private void loadEggnogHomology(){
		final ByteFile bf=ByteFile.makeByteFile(eggnogHomologyFile, true);
		final LineParser1 lp=new LineParser1('\t');
		boolean sawSchema=false, sawIdentity=false, sawConsensus=false, sawHeader=false;
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0){continue;}
				lp.set(line);
				if(line[0]=='#'){
					require(lp.terms()==2, "Bad eggNOG-homology metadata in "
						+eggnogHomologyFile+" line "+lineNumber);
					final String key=lp.parseString(0), value=lp.parseString(1);
					if(key.equals("#schema_version")){
						require(!sawSchema, "Duplicate eggNOG-homology schema in "
							+eggnogHomologyFile);
						require(value.equals("bbgap_eggnog_homology_v2"),
							"Unsupported eggNOG-homology schema in "+eggnogHomologyFile
							+": "+value);
						sawSchema=true;
					}else if(key.equals("#family_identity_sha80")){
						require(!sawIdentity, "Duplicate family_identity metadata in "
							+eggnogHomologyFile);
						require(value.equals(identityPin),
							"eggNOG-homology identity pin differs in "+eggnogHomologyFile);
						sawIdentity=true;
					}else if(key.equals("#consensus_sha80")){
						require(!sawConsensus, "Duplicate consensus metadata in "
							+eggnogHomologyFile);
						require(value.equals(consensusPin),
							"eggNOG-homology consensus pin differs in "+eggnogHomologyFile);
						sawConsensus=true;
					}else{throw new IllegalArgumentException("Unknown eggNOG-homology metadata "
						+key+" in "+eggnogHomologyFile);}
					continue;
				}
				if(!sawHeader){
					require(sawSchema && sawIdentity && sawConsensus,
						"eggNOG-homology rows require schema and pin metadata in "
						+eggnogHomologyFile);
					require(lp.terms()==EGGNOG_HOMOLOGY_COLUMNS.length,
						"Unexpected eggNOG-homology width in "+eggnogHomologyFile);
					for(int i=0; i<EGGNOG_HOMOLOGY_COLUMNS.length; i++){
						require(lp.parseString(i).equals(EGGNOG_HOMOLOGY_COLUMNS[i]),
							"Unexpected eggNOG-homology column "+i+" in "
							+eggnogHomologyFile+": "+lp.parseString(i));
					}
					sawHeader=true;
					continue;
				}
				require(lp.terms()==EGGNOG_HOMOLOGY_COLUMNS.length,
					"Unexpected eggNOG-homology width in "
					+eggnogHomologyFile+" line "+lineNumber);
				final Family family=requiredFamily(lp.parseInt(0), lp.parseInt(1),
					lp.parseString(2), eggnogHomologyFile, lineNumber);
				require(family.eggnog.isEmpty(),
					"Duplicate eggNOG-homology family_id "+family.familyId);
				family.eggnog.add(new EggnogHomology(lp, eggnogHomologyFile, lineNumber));
			}
		}finally{
			require(!bf.close(), "I/O error reading "+eggnogHomologyFile);
		}
		require(sawHeader, "Missing eggNOG-homology header in "+eggnogHomologyFile);
	}

	private void loadReviewedEvidence(){
		final ByteFile bf=ByteFile.makeByteFile(reviewedEvidenceFile, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		boolean sawHeader=false, sawSchema=false, sawIdentity=false, sawConsensus=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0){continue;}
				lp.set(line);
				if(line[0]=='#'){
					require(!sawHeader && lp.terms()==2,
						"Bad reviewed-evidence metadata in "+reviewedEvidenceFile
						+" line "+lineNumber);
					final String key=lp.parseString(0);
					final String value=lp.parseString(1);
					if(key.equals("#schema_version")){
						require(!sawSchema, "Duplicate schema metadata in "+reviewedEvidenceFile);
						require(value.equals("bbgap_reviewed_evidence_v1"),
							"Bad reviewed-evidence schema in "+reviewedEvidenceFile);
						sawSchema=true;
					}else if(key.equals("#family_identity_sha80")){
						require(!sawIdentity, "Duplicate identity metadata in "
							+reviewedEvidenceFile);
						require(value.equals(identityPin),
							"Reviewed-evidence identity pin differs in "+reviewedEvidenceFile);
						sawIdentity=true;
					}else if(key.equals("#consensus_sha80")){
						require(!sawConsensus, "Duplicate consensus metadata in "
							+reviewedEvidenceFile);
						require(value.equals(consensusPin),
							"Reviewed-evidence consensus pin differs in "+reviewedEvidenceFile);
						sawConsensus=true;
					}else{throw new IllegalArgumentException("Unknown reviewed-evidence metadata "
						+key+" in "+reviewedEvidenceFile);}
					continue;
				}
				if(!sawHeader){
					require(sawSchema && sawIdentity && sawConsensus,
						"Reviewed-evidence rows require schema and pin metadata in "
						+reviewedEvidenceFile);
					require(lp.terms()==REVIEWED_EVIDENCE_COLUMNS.length,
						"Unexpected reviewed-evidence width in "+reviewedEvidenceFile);
					for(int i=0; i<REVIEWED_EVIDENCE_COLUMNS.length; i++){
						require(lp.parseString(i).equals(REVIEWED_EVIDENCE_COLUMNS[i]),
							"Unexpected reviewed-evidence column "+i+" in "
							+reviewedEvidenceFile+": "+lp.parseString(i));
					}
					sawHeader=true;
					continue;
				}
				require(lp.terms()==REVIEWED_EVIDENCE_COLUMNS.length,
					"Unexpected reviewed-evidence width in "
					+reviewedEvidenceFile+" line "+lineNumber);
				final Family family=requiredFamily(lp.parseInt(0), lp.parseInt(1),
					lp.parseString(2), reviewedEvidenceFile, lineNumber);
				require(family.reviewed, "Reviewed-evidence row lacks reviewed= row for family_id "
					+family.familyId);
				require(family.reviewedEvidence==null,
					"Duplicate reviewed-evidence family_id "+family.familyId);
				final ReviewedEvidence evidence=new ReviewedEvidence(lp);
				evidence.validate(family, reviewedEvidenceFile, lineNumber);
				require(evidence.sourceMethods.length>0 && evidence.support.length>0
					&& evidence.coverage.length>0 && evidence.score.length>0
					&& evidence.productText.length>0 && evidence.review.length>0
					&& evidence.hitId.length>0,
					"Empty reviewed-evidence schema2 field in "+reviewedEvidenceFile
					+" line "+lineNumber);
				family.reviewedEvidence=evidence;
			}
		}finally{
			require(!bf.close(), "I/O error reading "+reviewedEvidenceFile);
		}
		require(sawHeader, "Missing reviewed-evidence header in "+reviewedEvidenceFile);
	}

	private void choose(){
		for(Family family : ordered){
			if(family.reviewed){
				reviewed++;
				if(family.product.length()>0){reviewedNamed++;}
				else{reviewedEmpty++;}
			}else if(mode.equals("hierarchy")){
				chooseHierarchy(family);
			}else if(family.ncbi.productCount()==1 && family.swiss.productCount()==1){
				final String ncbi=family.ncbi.firstKey();
				final String swiss=family.swiss.firstKey();
				if(ncbi.equals(swiss)){
					family.decision="specific_full_length";
					family.product=family.ncbi.firstProduct();
					family.rule=ignoreCase ? "case_insensitive_ncbi_swissprot_full_length_agreement" :
						"exact_ncbi_swissprot_full_length_agreement";
					autoSelected++;
					autoExact++;
					autoSpecific++;
					autoWithSwiss++;
				}else{
					family.decision="conflicting_evidence";
					family.reason="ncbi_swissprot_product_disagreement";
					unresolvedConflict++;
				}
			}else if(family.ncbi.productCount()>1 || family.swiss.productCount()>1){
				family.decision="conflicting_evidence";
				family.reason=family.ncbi.productCount()>1 ?
					"competing_ncbi_products" : "competing_swissprot_products";
				unresolvedConflict++;
			}else{
				family.decision="unnamed";
				family.reason="insufficient_full_length_agreement";
				unresolvedInsufficient++;
			}
		}
	}

	private void chooseHierarchy(final Family family){
		final ArrayList<NcbiModel> top=topNcbiTier(family);
		final Evidence products=topNcbiProducts(top);
		if(products.productCount()==1){
			final String product=products.firstProduct();
			if(isGenericUnknown(product)){
				family.decision="unnamed";
				family.reason="uninformative_native_ncbi_product";
				unresolvedInsufficient++;
				return;
			}
			final int specificity=supportedSpecificity(top);
			family.decision=specificity==SPECIFIC ? "specific_full_length" :
				"domain_supported_broad_name";
			family.product=product;
			family.rule=(specificity==SPECIFIC ? "ncbi_specific_" : "ncbi_broad_")
				+"rank"+top.get(0).precedence+"_native_name";
			autoSelected++;
			if(specificity==SPECIFIC){autoSpecific++;}
			else{autoBroad++;}
			if(family.swiss.productCount()>0){autoWithSwiss++;}
			else{autoNcbiOnly++;}
			if(hasFamilyType(top, "subfamily_domain")){autoSubfamilyDomain++;}
		}else if(products.productCount()>1){
			if(chooseCuratedBroadConflict(family, "ncbi_top_tier", products)){return;}
			if(chooseCuratedBroadParent(family, "ncbi_top_tier", products)){return;}
			family.decision="conflicting_evidence";
			family.reason="competing_top_ncbi_products";
			unresolvedConflict++;
		}else if(chooseSwissOnly(family)){
			return;
		}else if(family.swiss.productCount()>1){
			if(chooseCuratedBroadConflict(family, "swissprot", family.swiss)){return;}
			if(chooseCuratedBroadParent(family, "swissprot", family.swiss)){return;}
			family.decision="conflicting_evidence";
			family.reason="competing_swissprot_products";
			unresolvedConflict++;
		}else if(chooseNcbiNearQuery(family)){
			return;
		}else{
			family.decision="unnamed";
			family.reason="insufficient_native_ncbi_evidence";
			unresolvedInsufficient++;
		}
	}

	private boolean chooseSwissOnly(final Family family){
		if(!swissOnlyFallback){return false;}
		final SwissHit best=family.singletonInformativeSwiss();
		if(best==null || best.identity<minSwissIdentity){return false;}
		family.decision="swissprot_full_length";
		family.product=best.display;
		family.rule="swissprot_singleton_identity_"+minSwissIdentity;
		autoSelected++;
		autoSpecific++;
		autoWithSwiss++;
		autoSwissOnly++;
		return true;
	}

	private boolean chooseNcbiNearQuery(final Family family){
		if(!ncbiNearQueryFallback || family.swiss.productCount()!=1 ||
				family.ncbiNearModels.isEmpty()){
			return false;
		}
		final Evidence products=new Evidence();
		for(NcbiModel model : family.ncbiNearModels){
			products.add(model.accession, model.key, model.display, model.original);
		}
		if(products.productCount()!=1 || !products.firstKey().equals(family.swiss.firstKey())){
			return false;
		}
		final String product=products.firstProduct();
		if(isGenericUnknown(product)){return false;}
		for(NcbiModel model : family.ncbiNearModels){
			family.ncbi.add(model.accession, model.key, model.display, model.original);
		}
		family.decision="specific_full_length";
		family.product=product;
		family.rule="ncbi_near_query_swiss_supported_native_name";
		autoSelected++;
		autoSpecific++;
		autoWithSwiss++;
		autoNcbiNearQuery++;
		return true;
	}

	private boolean chooseCuratedBroadConflict(final Family family, final String source,
			final Evidence products){
		if(curatedBroadConflicts.isEmpty()){return false;}
		final ArrayList<String> keys=products.sortedKeys();
		final CuratedBroadConflict rule=curatedBroadConflicts.get(source+"\t"
			+Evidence.productSetKey(keys));
		if(rule==null){return false;}
		family.decision=rule.decision;
		family.product=rule.product;
		family.rule=rule.ruleId+"_"+source+"_exact_conflict_set";
		autoSelected++;
		autoBroad++;
		autoBroadConflict++;
		if(family.swiss.productCount()>0){autoWithSwiss++;}
		else{autoNcbiOnly++;}
		return true;
	}

	private boolean chooseCuratedBroadParent(final Family family, final String source,
			final Evidence products){
		if(curatedBroadParents.isEmpty()){return false;}
		CuratedBroadParent expected=null;
		final ArrayList<String> keys=products.sortedKeys();
		for(String key : keys){
			final CuratedBroadParent parent=curatedBroadParents.get(source+"\t"+key);
			if(parent==null){return false;}
			if(expected==null){expected=parent;}
			else if(expected!=parent){return false;}
		}
		if(expected==null){return false;}
		family.decision=expected.decision;
		family.product=expected.product;
		family.rule=expected.parentId+"_"+source+"_parent_vocabulary";
		autoSelected++;
		autoBroad++;
		autoBroadParent++;
		if(family.swiss.productCount()>0){autoWithSwiss++;}
		else{autoNcbiOnly++;}
		return true;
	}

	private void chooseProvisionalLabels(){
		for(Family family : ordered){
			if(family.product!=null || family.reviewed){continue;}
			final RawHit selected=chooseProvisional(family);
			if(selected==null){continue;}
			family.provisional=new Provisional(family, selected);
			provisionalLabelCount++;
			if(selected.informative){provisionalInformative++;}
			else{provisionalPlaceholder++;}
			if(family.provisional.sourcePriority){provisionalSourcePriority++;}
		}
	}

	private void chooseEggnogHomologyLabels(){
		for(Family family : ordered){
			if(family.product!=null || family.reviewed || family.provisional!=null ||
					family.eggnog.isEmpty()){
				continue;
			}
			assert(family.eggnog.size()==1) : "loadEggnogHomology rejects duplicates";
			final EggnogHomology selected=family.eggnog.get(0);
			if(!selected.eligible(eggnogDomainRelaxed)){continue;}
			family.provisional=new Provisional(selected, eggnogDomainRelaxed);
			provisionalLabelCount++;
			provisionalInformative++;
			provisionalEggnogHomology++;
		}
	}

	private RawHit chooseProvisional(final Family family){
		RawHit best=bestRawNcbi(family, true, true, true);
		if(best!=null){return best;}
		best=bestRawNcbi(family, false, true, true);
		if(best!=null){return best;}
		best=bestRawSwiss(family, true);
		if(best!=null){return best;}
		best=bestRawNcbi(family, false, false, false);
		return best==null ? bestRawSwiss(family, false) : best;
	}

	private RawHit bestRawNcbi(final Family family, final boolean nativeOnly,
			final boolean requireScores, final boolean requireInformative){
		RawHit best=null;
		for(RawHit hit : family.rawNcbi){
			if(nativeOnly && !hit.nativeEligible){continue;}
			if(requireScores && (!hit.sequencePass || !hit.domainPass || !hit.forNaming)){
				continue;
			}
			if(requireInformative && !hit.informative){continue;}
			if(hit.betterThan(best)){best=hit;}
		}
		return best;
	}

	private RawHit bestRawSwiss(final Family family, final boolean requireInformative){
		RawHit best=null;
		for(RawHit hit : family.rawSwiss){
			if(requireInformative && !hit.informative){continue;}
			if(hit.betterThan(best)){best=hit;}
		}
		return best;
	}

	private RawHit bestDifferentFullRawSwiss(final Family family, final RawHit selected){
		RawHit best=null;
		for(RawHit hit : family.rawSwiss){
			if(!hit.coveragePass || !hit.informative || hit.key.equals(selected.key)){
				continue;
			}
			if(hit.betterThan(best)){best=hit;}
		}
		return best;
	}

	private static Evidence topNcbiProducts(final ArrayList<NcbiModel> top){
		final Evidence products=new Evidence();
		for(NcbiModel model : top){
			products.add(model.accession, model.key, model.display, model.original);
		}
		return products;
	}

	private static ArrayList<NcbiModel> topNcbiTier(final Family family){
		int topPrecedence=Integer.MIN_VALUE;
		for(NcbiModel model : family.ncbiModels){
			if(model.precedence>topPrecedence){topPrecedence=model.precedence;}
		}
		final ArrayList<NcbiModel> top=new ArrayList<NcbiModel>();
		for(NcbiModel model : family.ncbiModels){
			if(model.precedence==topPrecedence){top.add(model);}
		}
		return top;
	}

	private static int supportedSpecificity(final ArrayList<NcbiModel> top){
		int specificity=SPECIFIC;
		for(NcbiModel model : top){
			if(model.specificity==BROAD){specificity=BROAD;}
		}
		return specificity;
	}

	private static boolean hasFamilyType(final ArrayList<NcbiModel> top, final String type){
		for(NcbiModel model : top){
			if(model.familyType.equals(type)){return true;}
		}
		return false;
	}

	private static boolean isGenericUnknown(final String product){
		final String p=toLowerCaseAscii(product);
		return p.equals("hypothetical protein") || p.equals("conserved hypothetical protein") ||
			p.equals("conserverd hypothetical protein") || p.equals("uncharacterized protein");
	}

	private static boolean isSwissPlaceholder(final String product){
		return !isInformative(product);
	}

	private static boolean isInformative(final String product){
		final String p=toLowerCaseAscii(product);
		return p.length()>0 &&
			!p.equals("hypothetical protein") &&
			!p.equals("conserved hypothetical protein") &&
			!p.startsWith("uncharacterized") &&
			!p.startsWith("upf") &&
			!p.startsWith("duf") &&
			!p.startsWith("putative uncharacterized ") &&
			!p.startsWith("probable uncharacterized ");
	}

	private void writeLabels(){
		final ByteStreamWriter bsw=writer(outFile);
		final ByteBuilder bb=new ByteBuilder(1<<14);
		try{
			writeLabelMetadata(bsw);
			for(Family family : ordered){
				if(family.product!=null || family.reviewed || family.provisional!=null){
					bb.clear();
					family.appendLabel(bb, labelVersion).nl();
					bsw.print(bb);
				}
			}
		}finally{
			close(bsw, outFile);
		}
	}

	private void writeEvidence(){
		final ByteStreamWriter bsw=writer(evidenceFile);
		final ByteBuilder bb=new ByteBuilder(1<<16);
		try{
			writeEvidenceHeader(bsw);
			for(Family family : ordered){
				bb.clear();
				family.appendEvidence(bb).nl();
				bsw.print(bb);
			}
			if(provisionalLabels){
				writeRawEvidenceHeader(bsw);
				for(Family family : ordered){
					if(family.provisional==null){continue;}
					for(RawHit hit : family.rawNcbi){
						bb.clear();
						family.appendRawEvidence(bb, hit).nl();
						bsw.print(bb);
					}
					for(RawHit hit : family.rawSwiss){
						bb.clear();
						family.appendRawEvidence(bb, hit).nl();
						bsw.print(bb);
					}
				}
			}
		}finally{
			close(bsw, evidenceFile);
		}
	}

	private void writeEggnogHomologyEvidence(){
		final ByteStreamWriter bsw=writer(eggnogHomologyOutFile);
		final ByteBuilder bb=new ByteBuilder(1<<14);
		try{
			bsw.println("#schema_version\tbbgap_eggnog_homology_v2");
			bsw.println("#family_identity_sha80\t"+identityPin);
			bsw.println("#consensus_sha80\t"+consensusPin);
			bsw.println(EGGNOG_HOMOLOGY_EVIDENCE_HEADER);
			for(Family family : ordered){
				for(EggnogHomology hit : family.eggnog){
					bb.clear();
					family.appendEggnogHomologyEvidence(bb, hit).nl();
					bsw.print(bb);
				}
			}
		}finally{
			close(bsw, eggnogHomologyOutFile);
		}
	}

	private void writeReviewedEvidence(){
		final ByteStreamWriter bsw=writer(reviewedEvidenceOutFile);
		final ByteBuilder bb=new ByteBuilder(1<<14);
		try{
			bsw.println("#schema_version\tbbgap_reviewed_evidence_v1");
			bsw.println("#family_identity_sha80\t"+identityPin);
			bsw.println("#consensus_sha80\t"+consensusPin);
			bsw.println(REVIEWED_EVIDENCE_HEADER);
			for(Family family : ordered){
				if(family.reviewedEvidence==null){continue;}
				bb.clear();
				family.reviewedEvidence.append(bb, family).nl();
				bsw.print(bb);
			}
		}finally{
			close(bsw, reviewedEvidenceOutFile);
		}
	}

	private void writeUnresolved(){
		final ByteStreamWriter bsw=writer(unresolvedFile);
		final ByteBuilder bb=new ByteBuilder(1<<16);
		try{
			writeEvidenceHeader(bsw);
			for(Family family : ordered){
				if(family.product==null || family.product.length()<1){
					bb.clear();
					family.appendEvidence(bb).nl();
					bsw.print(bb);
				}
			}
		}finally{
			close(bsw, unresolvedFile);
		}
	}

	private void writeSummary(){
		final ByteStreamWriter bsw=writer(summaryFile);
		try{
			bsw.println("metric\tvalue");
			bsw.println("total_families\t"+ordered.size());
			bsw.println("selected_auto\t"+autoSelected);
			bsw.println("selected_auto_exact\t"+autoExact);
			bsw.println("reviewed_overrides\t"+reviewed);
			bsw.println("reviewed_named\t"+reviewedNamed);
			bsw.println("reviewed_unresolved\t"+reviewedEmpty);
			bsw.println("selected_auto_specific\t"+autoSpecific);
			bsw.println("selected_auto_broad\t"+autoBroad);
			bsw.println("selected_auto_swissonly\t"+autoSwissOnly);
			if(broadConflictRulesFile!=null){
				bsw.println("selected_auto_broad_conflict\t"+autoBroadConflict);
			}
			if(broadParentRulesFile!=null){
				bsw.println("selected_auto_broad_parent\t"+autoBroadParent);
			}
			if(ncbiNearQueryFallback){
				bsw.println("selected_auto_ncbinearquery\t"+autoNcbiNearQuery);
			}
			bsw.println("selected_auto_with_swissprot\t"+autoWithSwiss);
			bsw.println("selected_auto_ncbi_only\t"+autoNcbiOnly);
			bsw.println("selected_auto_subfamily_domain\t"+autoSubfamilyDomain);
			bsw.println("unresolved_conflict\t"+unresolvedConflict);
			bsw.println("unresolved_insufficient\t"+unresolvedInsufficient);
			bsw.println("mode\t"+mode);
			bsw.println("ignore_case\t"+(ignoreCase ? "true" : "false"));
			bsw.println("swiss_only_fallback\t"+(swissOnlyFallback ? "true" : "false"));
			bsw.println("min_swiss_identity\t"+minSwissIdentity);
			bsw.println("label_rows\t"+(autoSelected+reviewed));
			bsw.println("named_products\t"+(autoSelected+reviewedNamed));
			bsw.println("unresolved\t"+(unresolvedConflict+unresolvedInsufficient+reviewedEmpty));
			if(provisionalLabels){
				bsw.println("label_rows_with_provisional\t"+
					(autoSelected+reviewed+provisionalLabelCount));
				bsw.println("provisional_labels\t"+provisionalLabelCount);
				bsw.println("provisional_informative\t"+provisionalInformative);
				bsw.println("provisional_placeholder\t"+provisionalPlaceholder);
				bsw.println("provisional_source_priority_partial_ncbi_over_full_swiss\t"+
					provisionalSourcePriority);
				if(eggnogHomologyFile!=null){
					bsw.println("provisional_eggnog_homology\t"+
						provisionalEggnogHomology);
				}
			}
		}finally{
			close(bsw, summaryFile);
		}
		require(autoSelected+reviewed+unresolvedConflict+unresolvedInsufficient==ordered.size(),
			"Naming decisions do not conserve families");
	}

	private void checkFiles(){
		require(Tools.testInputFiles(false, true, familiesFile, ncbiFile,
			ncbiNearQueryDomainsFile, swissFile, reviewedFile, synonymsFile,
				broadConflictRulesFile, broadParentRulesFile, reviewedEvidenceFile,
				eggnogHomologyFile),
			"Can't read input files");
		require(Tools.testOutputFiles(overwrite, false, false, outFile, evidenceFile,
				unresolvedFile, summaryFile, auditFile, conflictGroupsFile, memberCoverageFile,
				memberUnresolvedFile, reviewedEvidenceOutFile, eggnogHomologyOutFile),
			"Can't write output files");
		require(Tools.testForDuplicateFiles(true, familiesFile, ncbiFile,
			ncbiNearQueryDomainsFile, swissFile, reviewedFile, synonymsFile,
			broadConflictRulesFile, broadParentRulesFile, outFile, evidenceFile,
			unresolvedFile, summaryFile, auditFile, conflictGroupsFile, memberCoverageFile,
				memberUnresolvedFile, reviewedEvidenceFile, reviewedEvidenceOutFile,
				eggnogHomologyFile, eggnogHomologyOutFile),
			"Duplicate input or output path");
		checkExistingAliases();
	}

	private void checkExistingAliases(){
		final String[] inputs={familiesFile, ncbiFile, ncbiNearQueryDomainsFile, swissFile,
			reviewedFile, synonymsFile, broadConflictRulesFile, broadParentRulesFile,
				reviewedEvidenceFile, eggnogHomologyFile};
		final String[] outputs={outFile, evidenceFile, unresolvedFile, summaryFile, auditFile,
			conflictGroupsFile, memberCoverageFile, memberUnresolvedFile,
			reviewedEvidenceOutFile, eggnogHomologyOutFile};
		for(String in : inputs){
			if(in==null){continue;}
			final Path ip=Paths.get(in);
			for(String out : outputs){
				if(out==null){continue;}
				final Path op=Paths.get(out);
				try{
					require(!Files.exists(op) || !Files.isSameFile(ip, op),
						"Output aliases input: "+out+" -> "+in);
				}catch(IOException e){
					throw new IllegalArgumentException("Could not check file alias: "
						+in+" vs "+out, e);
				}
			}
		}
		for(int i=0; i<outputs.length; i++){
			if(outputs[i]==null){continue;}
			final Path a=Paths.get(outputs[i]);
			if(!Files.exists(a)){continue;}
			for(int j=i+1; j<outputs.length; j++){
				if(outputs[j]==null){continue;}
				final Path b=Paths.get(outputs[j]);
				try{
					require(!Files.exists(b) || !Files.isSameFile(a, b),
						"Output aliases output: "+outputs[i]+" -> "+outputs[j]);
				}catch(IOException e){
					throw new IllegalArgumentException("Could not check file alias: "
						+outputs[i]+" vs "+outputs[j], e);
				}
			}
		}
	}

	private void writeAudit(){
		final ArrayList<Family> auto=new ArrayList<Family>(), conflicts=new ArrayList<Family>();
		final ArrayList<Family> insufficient=new ArrayList<Family>();
		for(Family family : ordered){
			if(!family.reviewed && family.product!=null){auto.add(family);}
			else if(family.reason.startsWith("competing") ||
					family.reason.equals("ncbi_swissprot_product_disagreement")){
				conflicts.add(family);
			}else if(!family.reviewed){insufficient.add(family);}
		}
		final ByteStreamWriter bsw=writer(auditFile);
		final ByteBuilder bb=new ByteBuilder(1<<16);
		try{
			bsw.print("audit_stratum\t");
			writeEvidenceHeader(bsw);
			writeAuditStratum(bsw, bb, auto, "selected_auto");
			writeAuditStratum(bsw, bb, conflicts, "conflict");
			writeAuditStratum(bsw, bb, insufficient, "insufficient");
		}finally{
			close(bsw, auditFile);
		}
	}

	private void writeAuditStratum(final ByteStreamWriter bsw, final ByteBuilder bb,
			final ArrayList<Family> list, final String stratum){
		Collections.sort(list, auditOrder);
		final int limit=Math.min(maxAudit, list.size());
		for(int i=0; i<limit; i++){
			final Family family=list.get(i);
			bb.clear().append(stratum).tab();
			family.appendEvidence(bb).nl();
			bsw.print(bb);
		}
	}

	private void writeConflictGroups(){
		final HashMap<String, ConflictGroup> map=new HashMap<String, ConflictGroup>();
		for(Family family : ordered){
			final boolean ncbi=family.reason.equals("competing_top_ncbi_products");
			final boolean swiss=family.reason.equals("competing_swissprot_products");
			if(!ncbi && !swiss){continue;}
			final String source=ncbi ? "ncbi_top_tier" : "swissprot";
			final Evidence evidence=ncbi ? topNcbiProducts(topNcbiTier(family)) : family.swiss;
			require(evidence.productCount()>1, "Conflict without competing products for family "
				+family.familyId);
			final ArrayList<String> keys=evidence.sortedKeys();
			final String normalized=Evidence.join(keys, '|');
			final String products=evidence.productsForKeys(keys);
			final String key=source+"\t"+Evidence.productSetKey(keys);
			ConflictGroup group=map.get(key);
			if(group==null){
				group=new ConflictGroup(source, family.reason, keys.size(), normalized, products);
				map.put(key, group);
			}
			group.add(family.familyId);
		}
		final ArrayList<ConflictGroup> groups=new ArrayList<ConflictGroup>(map.values());
		Collections.sort(groups, new Comparator<ConflictGroup>() {
			@Override
			public int compare(final ConflictGroup a, final ConflictGroup b){
				if(a.familyIds.size()!=b.familyIds.size()){
					return b.familyIds.size()-a.familyIds.size();
				}
				final int x=a.source.compareTo(b.source);
				return x==0 ? a.normalizedProducts.compareTo(b.normalizedProducts) : x;
			}
		});
		final ByteStreamWriter bsw=writer(conflictGroupsFile);
		final ByteBuilder bb=new ByteBuilder(1<<14);
		try{
			bsw.println("rank\tconflict_source\treason\tfamily_count\tproduct_count\t"
				+"family_ids\tnormalized_products\tdisplay_products");
			for(int i=0; i<groups.size(); i++){
				final ConflictGroup group=groups.get(i);
				bb.clear().append(i+1).tab().append(group.source).tab().append(group.reason)
					.tab().append(group.familyIds.size()).tab().append(group.productCount)
					.tab().append(group.familyIds()).tab().append(group.normalizedProducts)
					.tab().append(group.displayProducts).nl();
				bsw.print(bb);
			}
		}finally{
			close(bsw, conflictGroupsFile);
		}
	}

	private void writeMemberCoverage(){
		final LinkedHashMap<String, MemberBucket> buckets=
			new LinkedHashMap<String, MemberBucket>();
		long totalMembers=0, totalRawResidues=0;
		for(Family family : ordered){
			totalMembers=Math.addExact(totalMembers, family.members);
			totalRawResidues=Math.addExact(totalRawResidues, family.rawResidues);
			final String bucket=memberBucket(family);
			final String reason=memberReason(family);
			final String key=bucket+"\t"+family.decision+"\t"+reason;
			MemberBucket mb=buckets.get(key);
			if(mb==null){
				mb=new MemberBucket(bucket, family.decision, reason);
				buckets.put(key, mb);
			}
			mb.add(family);
		}
		final ArrayList<MemberBucket> list=new ArrayList<MemberBucket>(buckets.values());
		Collections.sort(list, memberBucketOrder);
		final ByteStreamWriter bsw=writer(memberCoverageFile);
		final ByteBuilder bb=new ByteBuilder(1<<14);
		try{
			bsw.println("bucket\tdecision_class\trule_or_reason\tfamily_count\tmembers\t"
				+"raw_residues\ttotal_families\ttotal_members\ttotal_raw_residues\t"
				+"member_fraction\traw_residue_fraction");
			for(MemberBucket bucket : list){
				bb.clear().append(bucket.bucket).tab().append(bucket.decision).tab()
					.append(bucket.reason).tab().append(bucket.families).tab()
					.append(bucket.members).tab().append(bucket.rawResidues).tab()
					.append(ordered.size()).tab().append(totalMembers).tab()
					.append(totalRawResidues).tab();
				bb.appendSlow(bucket.members/(double)totalMembers, 12).tab();
				if(totalRawResidues>0){bb.appendSlow(bucket.rawResidues/(double)totalRawResidues, 12);}
				else{bb.append('.');}
				bb.nl();
				bsw.print(bb);
			}
		}finally{
			close(bsw, memberCoverageFile);
		}
	}

	private void writeMemberUnresolved(){
		final ArrayList<Family> list=new ArrayList<Family>();
		for(Family family : ordered){
			if(family.product==null || family.product.length()<1){list.add(family);}
		}
		Collections.sort(list, memberUnresolvedOrder);
		final ByteStreamWriter bsw=writer(memberUnresolvedFile);
		final ByteBuilder bb=new ByteBuilder(1<<14);
		try{
			bsw.println("rank\tfamily_id\tactive_index\trep_id\tdecision_class\t"
				+"rule_or_reason\tmembers\traw_residues\tncbi_products\t"
				+"swissprot_products");
			for(int i=0; i<list.size(); i++){
				final Family family=list.get(i);
				bb.clear().append(i+1).tab().append(family.familyId).tab()
					.append(family.activeIndex).tab().append(family.repId).tab()
					.append(family.decision).tab().append(memberReason(family)).tab()
					.append(family.members).tab().append(family.rawResidues).tab()
					.append(memberProducts(family.ncbi)).tab()
					.append(memberProducts(family.swiss)).nl();
				bsw.print(bb);
			}
		}finally{
			close(bsw, memberUnresolvedFile);
		}
	}

	private Family requiredFamily(final int familyId, final int active, final String rep,
			final String path, final long lineNumber){
		final Family family=families.get(familyId);
		require(family!=null, "Unknown family_id "+familyId+" in "+path+" line "+lineNumber);
		require(family.activeIndex==active, "Active-index mismatch in "+path+" line "+lineNumber);
		require(rep==null || family.repId.equals(rep),
			"Representative mismatch in "+path+" line "+lineNumber);
		return family;
	}

	private static String swissProduct(final String title, final String path, final long lineNumber){
		final int start=title.startsWith("sp|") ? title.indexOf(' ')+1 : 0;
		require(start>0, "Bad Swiss-Prot title in "+path+" line "+lineNumber+": "+title);
		int stop=title.indexOf(" OS=", start);
		if(stop<0){stop=title.length();}
		final String product=title.substring(start, stop);
		require(product.length()>0, "Empty Swiss-Prot product in "+path+" line "+lineNumber);
		return product;
	}

	private String displayProduct(final String product){
		final String x=synonyms.get(comparisonKey(product));
		return x==null ? product : x;
	}

	private String comparisonKey(final String display){
		return ignoreCase ? toLowerCaseAscii(display) : display;
	}

	private static String toLowerCaseAscii(final String s){
		char[] chars=null;
		for(int i=0; i<s.length(); i++){
			final char c=s.charAt(i);
			if(c>='A' && c<='Z'){
				if(chars==null){chars=s.toCharArray();}
				chars[i]=(char)(c+32);
			}
		}
		return chars==null ? s : new String(chars);
	}

	private static String blankDot(final String s){return s.equals(".") ? "" : s;}

	private static int headerIndex(final LineParser1 lp, final String name, final String path){
		int found=-1;
		for(int i=0; i<lp.terms(); i++){
			if(lp.parseString(i).equals(name)){
				require(found<0, "Duplicate family_identity column "+name+" in "+path);
				found=i;
			}
		}
		return found;
	}

	private static long parsePositiveLong(final LineParser1 lp, final int column,
			final String path, final long lineNumber, final String field){
		final long x=parseLong(lp, column, path, lineNumber, field);
		require(x>0, "Bad "+field+" in "+path+" line "+lineNumber+": "+x);
		return x;
	}

	private static long parseNonnegativeLong(final LineParser1 lp, final int column,
			final String path, final long lineNumber, final String field){
		final long x=parseLong(lp, column, path, lineNumber, field);
		require(x>=0, "Bad "+field+" in "+path+" line "+lineNumber+": "+x);
		return x;
	}

	private static long parseLong(final LineParser1 lp, final int column,
			final String path, final long lineNumber, final String field){
		final String s=lp.parseString(column);
		try{
			return Long.parseLong(s);
		}catch(NumberFormatException nfe){
			throw new IllegalArgumentException("Bad "+field+" in "+path+" line "
				+lineNumber+": "+s, nfe);
		}
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

	private static double parseFinite(final String text, final String path,
			final long lineNumber, final String field){
		final double value=Double.parseDouble(text);
		require(Double.isFinite(value), "Bad NCBI "+field+" in "+path+" line "+lineNumber
			+": "+text);
		return value;
	}

	private static double parseCoverage(final String text, final String path,
			final long lineNumber, final String field){
		final double value=parseFinite(text, path, lineNumber, field);
		require(value>=0 && value<=1, "Bad NCBI "+field+" in "+path+" line "+lineNumber
			+": "+text);
		return value;
	}

	private void writeLabelMetadata(final ByteStreamWriter bsw){
		bsw.println("#schema_version\t"+(schema2Labels() ? "2" : "1"));
		bsw.println("#kind\tbbgap_family_labels");
		bsw.println("#label_version\t"+labelVersion);
		bsw.println("#family_identity_sha80\t"+identityPin);
		bsw.println("#consensus_sha80\t"+consensusPin);
		bsw.print("family_id\tactive_index\trep_id\tdecision_class\tproduct\tlabel_version");
		if(schema2Labels()){
			bsw.print("\tdisplay_label\tsource_methods\tsupport_tag\tcoverage_tag\t"+
				"score_tag\tproduct_text_tag\treview_tag\tselected_hit_id");
		}
		bsw.println();
	}

	private boolean schema2Labels(){return provisionalLabels || reviewedEvidenceFile!=null;}

	private static void writeEvidenceHeader(final ByteStreamWriter bsw){
		bsw.println("family_id\tactive_index\trep_id\tdecision_class\tproduct\trule_or_reason\t"
			+"ncbi_product_count\tswissprot_product_count\tncbi_accessions\t"
			+"swissprot_accessions\tncbi_products\tswissprot_products\tswissprot_overlay\t"
			+"swissprot_best_accession\tswissprot_best_product\tswissprot_best_identity\t"
			+"swissprot_best_query_cov\tswissprot_best_subject_cov\t"
			+"swissprot_best_bitscore\tswissprot_second_accession\t"
			+"swissprot_second_product\tswissprot_second_identity\t"
			+"swissprot_second_query_cov\tswissprot_second_subject_cov\t"
			+"swissprot_second_bitscore\tswissprot_best_second_bits_delta\t"
			+"swissprot_best_second_bits_ratio\tswissprot_match_accession\t"
			+"swissprot_match_product\tswissprot_match_identity\t"
			+"swissprot_match_query_cov\tswissprot_match_subject_cov\t"
			+"swissprot_match_bitscore\tswissprot_diff_accession\t"
			+"swissprot_diff_product\tswissprot_diff_identity\t"
			+"swissprot_diff_query_cov\tswissprot_diff_subject_cov\t"
			+"swissprot_diff_bitscore\tswissprot_match_diff_bits_delta\t"
			+"swissprot_match_diff_bits_ratio");
	}

	private static void writeRawEvidenceHeader(final ByteStreamWriter bsw){
		bsw.println("#raw_provisional_evidence");
		bsw.println("family_id\tactive_index\trep_id\thit_id\tselected\tsource_methods\t"
			+"accession\tproduct\tsupport_tag\tcoverage_tag\tscore_tag\tsequence_score\t"
			+"domain_score\tquery_cov\ttarget_cov\tsequence_pass\tdomain_pass\t"
			+"coverage_pass\tnative_eligible\tfor_naming\toriginal_decision\t"
			+"original_reason");
	}

	private final Comparator<Family> auditOrder=new Comparator<Family>() {
		@Override
		public int compare(final Family a, final Family b){
			final int x=Long.compare(auditKey(a.familyId), auditKey(b.familyId));
			return x!=0 ? x : a.activeIndex-b.activeIndex;
		}
	};

	private long auditKey(final int familyId){
		long x=familyId+seed+0x9E3779B97F4A7C15L;
		x=(x^(x>>>30))*0xBF58476D1CE4E5B9L;
		x=(x^(x>>>27))*0x94D049BB133111EBL;
		return (x^(x>>>31))&Long.MAX_VALUE;
	}

	private boolean memberReportsRequested(){
		return memberCoverageFile!=null || memberUnresolvedFile!=null;
	}

	private static String memberBucket(final Family family){
		if(family.product!=null && family.product.length()>0){
			return family.reviewed ? "named_reviewed" : "named_auto";
		}
		if(family.reviewed){return "unresolved_reviewed";}
		return family.decision.equals("conflicting_evidence") ?
			"unresolved_conflict" : "unresolved_insufficient";
	}

	private static String memberReason(final Family family){
		if(family.product!=null && family.product.length()>0){return "";}
		return family.reviewed ? "reviewed_override" : family.reason;
	}

	private static String memberProducts(final Evidence evidence){
		final String products=evidence.products();
		return products.length()<1 ? "." : products;
	}

	private final Comparator<MemberBucket> memberBucketOrder=new Comparator<MemberBucket>() {
		@Override
		public int compare(final MemberBucket a, final MemberBucket b){
			final int x=bucketRank(a.bucket)-bucketRank(b.bucket);
			if(x!=0){return x;}
			final int y=a.decision.compareTo(b.decision);
			return y==0 ? a.reason.compareTo(b.reason) : y;
		}
	};

	private final Comparator<Family> memberUnresolvedOrder=new Comparator<Family>() {
		@Override
		public int compare(final Family a, final Family b){
			final int x=Long.compare(b.members, a.members);
			return x==0 ? Integer.compare(a.familyId, b.familyId) : x;
		}
	};

	private static int bucketRank(final String bucket){
		if(bucket.equals("named_auto")){return 0;}
		if(bucket.equals("named_reviewed")){return 1;}
		if(bucket.equals("unresolved_conflict")){return 2;}
		if(bucket.equals("unresolved_insufficient")){return 3;}
		if(bucket.equals("unresolved_reviewed")){return 4;}
		return 5;
	}

	private final class Family {
		Family(final int familyId_, final int activeIndex_, final String repId_){
			familyId=familyId_;
			activeIndex=activeIndex_;
			repId=repId_;
		}

		void addNcbi(final String accession, final String product, final String familyType){
			addNcbi(accession, product, familyType, false);
		}

		void addNcbi(final String accession, final String product, final String familyType,
				final boolean nearQuery){
			require(product.length()>0, "Empty NCBI product for "+accession);
			final String display=displayProduct(product);
			final String key=comparisonKey(display);
			if(nearQuery){ncbiNearModels.add(new NcbiModel(accession, key, display,
				product, familyType));}
			else{
				ncbi.add(accession, key, display, product);
				if(familyType!=null){
					ncbiModels.add(new NcbiModel(accession, key, display, product, familyType));
				}
			}
		}

		void addSwiss(final String accession, final String product, final String identity,
				final String queryCov, final String subjectCov, final String bitscore){
			final String display=displayProduct(product);
			final String key=comparisonKey(display);
			swiss.add(accession, key, display, product);
			swissHits.add(new SwissHit(accession, key, display, identity, queryCov,
				subjectCov, bitscore));
		}

		void addRawNcbi(final String accession, final String product, final String familyType,
				final String fullScore, final String domainScore, final String queryCov,
				final String targetCov, final boolean sequencePass, final boolean domainPass,
				final boolean coveragePass, final boolean nativeEligible,
				final boolean forNaming){
			final String display=displayProduct(product);
			rawNcbi.add(new RawHit("ncbi:"+rawNcbi.size(), "NCBIfam/HMMER", accession,
				comparisonKey(display), display, familyType, fullScore, domainScore,
				queryCov, targetCov,
				sequencePass, domainPass, coveragePass, nativeEligible, forNaming,
				fullScore, minCoverage, isInformative(display)));
		}

		void addRawSwiss(final String accession, final String product, final String identity,
				final String queryCov, final String subjectCov, final String bitscore){
			final String display=displayProduct(product);
			final double q=Double.parseDouble(queryCov);
			final double s=Double.parseDouble(subjectCov);
			rawSwiss.add(new RawHit("swiss:"+rawSwiss.size(), "Swiss-Prot/DIAMOND",
				accession, comparisonKey(display), display, ".", bitscore, bitscore,
				queryCov, subjectCov,
				true, true, q>=minCoverage && s>=minCoverage, false, true,
				bitscore, minCoverage, !isSwissPlaceholder(display)));
		}

		boolean hasSwissHit(final String accession){
			for(SwissHit hit : swissHits){
				if(hit.accession.equals(accession)){return true;}
			}
			return false;
		}

		boolean hasRawNcbiHit(final String id){
			for(RawHit hit : rawNcbi){
				if(hit.id.equals(id)){return true;}
			}
			return false;
		}

		ByteBuilder appendLabel(final ByteBuilder bb, final String labelVersion){
			final Provisional p=product==null ? provisional : null;
			bb.append(familyId).tab().append(activeIndex).tab().append(repId).tab()
				.append(p==null ? decision : p.decision).tab()
				.append(p==null ? product==null ? "" : product : p.cleanProduct).tab()
				.append(labelVersion);
			if(schema2Labels()){
				if(p==null && reviewedEvidence!=null){
					reviewedEvidence.appendLabel(bb, product==null ? "" : product);
				}else{
					bb.tab().append(p==null ? displayLabel() : p.display).tab()
						.append(sourceMethods(p)).tab()
						.append(p==null ? acceptedSupportTag() : p.support).tab()
						.append(p==null ? acceptedCoverageTag() : p.coverage).tab()
						.append(p==null ? acceptedScoreTag() : p.score).tab()
						.append(p==null ? productTextTag() : p.productText).tab()
						.append(p==null ? reviewTag() : p.review).tab()
						.append(p==null ? selectedHitId() : p.hitId);
				}
			}
			return bb;
		}

		String displayLabel(){
			return product==null ? "" : product;
		}

		String sourceMethods(final Provisional p){
			if(p!=null){return p.sourceMethod;}
			if(product==null || product.length()<1){return "unknown";}
			final String key=comparisonKey(product);
			final boolean n=ncbi.productMap.containsKey(key);
			final boolean s=swiss.productMap.containsKey(key);
			if(n && s){return "NCBIfam/HMMER,Swiss-Prot/DIAMOND";}
			if(rule.startsWith("reviewed_")){return "unknown";}
			if(decision.equals("swissprot_full_length") || rule.indexOf("swissprot")>=0){
				return "Swiss-Prot/DIAMOND";
			}
			if(rule.indexOf("ncbi")>=0){return "NCBIfam/HMMER";}
			return "unknown";
		}

		String acceptedSupportTag(){
			if(product==null || product.length()<1){return "unknown";}
			if(rule.startsWith("reviewed_")){return "reviewed";}
			if(rule.indexOf("swissprot")>=0){return "swiss_full";}
			if(rule.startsWith("ncbi_near_query_")){return "native_near_query";}
			if(sourceMethods(null).indexOf(',')>=0){return "native_and_swiss_full";}
			if(rule.indexOf("ncbi")>=0){return "native_full";}
			return "unknown";
		}

		String acceptedCoverageTag(){
			if(product==null || product.length()<1 || rule.startsWith("reviewed_")){
				return "unknown";
			}
			if(rule.startsWith("ncbi_near_query_")){return "query_low";}
			if(sourceMethods(null).equals("unknown")){return "unknown";}
			return "coverage_full";
		}

		String acceptedScoreTag(){
			if(product==null || product.length()<1 || rule.startsWith("reviewed_")){
				return "unknown";
			}
			return sourceMethods(null).equals("Swiss-Prot/DIAMOND") ?
				"retained_search_hit" : "score_passed";
		}

		String productTextTag(){
			return "not_applicable";
		}

		String reviewTag(){
			return "none";
		}

		String selectedHitId(){
			return "";
		}

		ByteBuilder appendEvidence(final ByteBuilder bb){
			return appendCoreEvidence(bb.append(familyId).tab());
		}

		ByteBuilder appendRawEvidence(final ByteBuilder bb, final RawHit hit){
			final Provisional p=provisional;
			return bb.append(familyId).tab().append(activeIndex).tab().append(repId).tab()
				.append(hit.id).tab().append(p!=null && p.selected==hit ? 'Y' : 'N')
				.tab().append(hit.method).tab().append(hit.accession).tab()
				.append(hit.display).tab().append(hit.supportTag()).tab()
				.append(hit.coverageTag()).tab().append(hit.scoreTag()).tab()
				.append(hit.score).tab().append(hit.domainScore).tab()
				.appendSlow(hit.queryCov, 12).tab().appendSlow(hit.targetCov, 12).tab()
				.append(hit.sequencePass ? 'Y' : 'N').tab()
				.append(hit.domainPass ? 'Y' : 'N').tab()
				.append(hit.coveragePass ? 'Y' : 'N').tab()
				.append(hit.nativeEligible ? 'Y' : 'N').tab()
				.append(hit.forNaming ? 'Y' : 'N').tab()
				.append(decision).tab().append(reason);
		}

		ByteBuilder appendEggnogHomologyEvidence(final ByteBuilder bb,
				final EggnogHomology hit){
			final Provisional p=provisional;
			return hit.append(bb.append(familyId).tab().append(activeIndex).tab()
				.append(repId).tab()).tab()
				.append(p!=null && p.eggnog==hit ? 'Y' : 'N').tab()
				.append(decision).tab().append(reason).tab()
				.append(eggnogHomologyReason(hit));
		}

		String eggnogHomologyReason(final EggnogHomology hit){
			if(provisional!=null && provisional.eggnog==hit){return ".";}
			if(reviewed){return "reviewed_override";}
			if(product!=null){return "existing_product";}
			if(provisional!=null){return "existing_provisional";}
			return hit.ineligibleReason;
		}

		ByteBuilder appendCoreEvidence(final ByteBuilder bb){
			final SwissHit best=bestSwiss(null), second=bestSwiss(best==null ? null : best.key);
			final String productKey=(product==null || product.length()<1) ? null :
				comparisonKey(product);
			final SwissHit match=productKey==null ? null : bestSwiss(productKey, true);
			final SwissHit diff=productKey==null ? null : bestSwiss(productKey, false);
			bb.append(activeIndex).tab().append(repId).tab().append(decision).tab()
				.append(product==null ? "" : product).tab()
				.append(product==null ? reason : rule).tab()
				.append(ncbi.productCount()).tab().append(swiss.productCount()).tab()
				.append(ncbi.accessions()).tab().append(swiss.accessions()).tab()
				.append(ncbi.products()).tab().append(swiss.products()).tab()
				.append(swissOverlay()).tab();
			appendSwissHit(bb, best);
			bb.tab();
			appendSwissHit(bb, second);
			bb.tab();
			if(best!=null && second!=null){bb.appendSlow(best.bits-second.bits, 3);}
			bb.tab();
			if(best!=null && second!=null && second.bits!=0){
				bb.appendSlow(best.bits/second.bits, 6);
			}
			bb.tab();
			appendSwissHit(bb, match);
			bb.tab();
			appendSwissHit(bb, diff);
			bb.tab();
			if(match!=null && diff!=null){bb.appendSlow(match.bits-diff.bits, 3);}
			bb.tab();
			if(match!=null && diff!=null && diff.bits!=0){
				bb.appendSlow(match.bits/diff.bits, 6);
			}
			return bb;
		}

		String swissOverlay(){
			if(product==null || product.length()<1){return "not_selected";}
			if(swissHits.isEmpty()){return "absent";}
			final String key=comparisonKey(product);
			boolean matches=false, differs=false;
			for(SwissHit hit : swissHits){
				if(hit.key.equals(key)){matches=true;}
				else{differs=true;}
			}
			return matches ? (differs ? "mixed" : "agreement") : "mismatch";
		}

		SwissHit bestSwiss(final String excludedKey){
			SwissHit best=null;
			for(SwissHit hit : swissHits){
				if(hit.key.equals(excludedKey)){continue;}
				if(better(hit, best)){best=hit;}
			}
			return best;
		}

		SwissHit bestSwiss(final String key, final boolean include){
			SwissHit best=null;
			for(SwissHit hit : swissHits){
				if(hit.key.equals(key)!=include){continue;}
				if(better(hit, best)){best=hit;}
			}
			return best;
		}

		SwissHit singletonInformativeSwiss(){
			String key=null;
			SwissHit best=null;
			boolean placeholder=false;
			for(SwissHit hit : swissHits){
				if(key==null){key=hit.key;}
				else if(!key.equals(hit.key)){return null;}
				if(isSwissPlaceholder(hit.display)){placeholder=true;}
				if(better(hit, best)){best=hit;}
			}
			return placeholder ? null : best;
		}

		final int familyId, activeIndex;
		final String repId;
		final Evidence ncbi=new Evidence(), swiss=new Evidence();
		final ArrayList<NcbiModel> ncbiModels=new ArrayList<NcbiModel>();
		final ArrayList<NcbiModel> ncbiNearModels=new ArrayList<NcbiModel>();
		final ArrayList<SwissHit> swissHits=new ArrayList<SwissHit>();
		final ArrayList<RawHit> rawNcbi=new ArrayList<RawHit>();
		final ArrayList<RawHit> rawSwiss=new ArrayList<RawHit>();
		final ArrayList<EggnogHomology> eggnog=new ArrayList<EggnogHomology>();
		Provisional provisional;
		ReviewedEvidence reviewedEvidence;
		String decision, product, rule="", reason="";
		long members=-1, rawResidues=-1;
		boolean reviewed;
	}

	private static final class MemberBucket {
		MemberBucket(final String bucket_, final String decision_, final String reason_){
			bucket=bucket_;
			decision=decision_;
			reason=reason_;
		}

		void add(final Family family){
			families++;
			members=Math.addExact(members, family.members);
			rawResidues=Math.addExact(rawResidues, family.rawResidues);
		}

		final String bucket, decision, reason;
		long families, members, rawResidues;
	}

	private static void appendSwissHit(final ByteBuilder bb, final SwissHit hit){
		if(hit==null){bb.tab().tab().tab().tab().tab();}
		else{
			bb.append(hit.accession).tab().append(hit.display).tab()
				.append(hit.identityText).tab().append(hit.queryCovText).tab()
				.append(hit.subjectCovText).tab().append(hit.bitsText);
		}
	}

	private static boolean better(final SwissHit a, final SwissHit b){
		if(b==null){return true;}
		if(a.bits!=b.bits){return a.bits>b.bits;}
		if(a.identity!=b.identity){return a.identity>b.identity;}
		return a.accession.compareTo(b.accession)<0;
	}

	private static final class SwissHit {
		SwissHit(final String accession_, final String key_, final String display_,
				final String identity_, final String queryCov_, final String subjectCov_,
				final String bits_){
			accession=accession_;
			key=key_;
			display=display_;
			identityText=identity_;
			queryCovText=queryCov_;
			subjectCovText=subjectCov_;
			bitsText=bits_;
			identity=Double.parseDouble(identity_);
			queryCov=Double.parseDouble(queryCov_);
			subjectCov=Double.parseDouble(subjectCov_);
			bits=Double.parseDouble(bits_);
			require(Double.isFinite(identity) && identity>=0 && identity<=100,
				"Bad Swiss-Prot identity for "+accession+": "+identityText);
			require(Double.isFinite(queryCov) && queryCov>=0 && queryCov<=1,
				"Bad Swiss-Prot query coverage for "+accession+": "+queryCovText);
			require(Double.isFinite(subjectCov) && subjectCov>=0 && subjectCov<=1,
				"Bad Swiss-Prot subject coverage for "+accession+": "+subjectCovText);
			require(Double.isFinite(bits), "Bad Swiss-Prot bitscore for "+accession+": "
				+bitsText);
		}
		final String accession, key, display, identityText, queryCovText, subjectCovText;
		final String bitsText;
		final double identity, queryCov, subjectCov, bits;
	}

	private final class Provisional {
		Provisional(final Family family, final RawHit selected_){
			selected=selected_;
			eggnog=null;
			decision=selected.informative && selected.coveragePass && selected.sequencePass
				&& selected.domainPass ? "provisional_full" :
				selected.informative ? "weak_provisional" : "weak_provisional_placeholder";
			cleanProduct=selected.source.equals("swiss") && decision.equals("weak_provisional")
				? "tentative homology to "+selected.display : selected.display;
			sourceMethod=provisionalSourceMethods(family, selected);
			support=selected.supportTag();
			coverage=selected.coverageTag();
			score=selected.scoreTag();
			productText=productTextTag(family, selected);
			sourcePriority=selected.source.equals("ncbi") && !selected.nativeEligible &&
				bestDifferentFullRawSwiss(family, selected)!=null;
			review=sourcePriority ? "source_priority_partial_ncbi_over_full_swiss" : "none";
			hitId=selected.id;
			display=decision+": "+cleanProduct+" ["+sourceMethod+"]";
		}

		Provisional(final EggnogHomology selected_, final boolean relaxed){
			selected=null;
			eggnog=selected_;
			decision="weak_provisional";
			cleanProduct=selected_.product();
			sourceMethod="eggNOG/DIAMOND";
			support="same_seed_domain_family_homology";
			coverage=selected_.spanBin;
			score="selected_seed_60_80_80_pass";
			productText=selected_.productTextTag(relaxed);
			review=selected_.reviewTag(relaxed);
			hitId=selected_.hitId();
			display=decision+": same-seed homology to "+cleanProduct+" ["+sourceMethod+"]";
			sourcePriority=false;
		}

		private String provisionalSourceMethods(final Family family, final RawHit selected){
			boolean ncbi=selected.source.equals("ncbi"), swiss=selected.source.equals("swiss");
			for(RawHit hit : family.rawNcbi){
				if(hit.informative && hit.key.equals(selected.key)){ncbi=true;}
			}
			for(RawHit hit : family.rawSwiss){
				if(hit.informative && hit.key.equals(selected.key)){swiss=true;}
			}
			return ncbi && swiss ? "NCBIfam/HMMER,Swiss-Prot/DIAMOND" :
				ncbi ? "NCBIfam/HMMER" : "Swiss-Prot/DIAMOND";
		}

		private String productTextTag(final Family family, final RawHit selected){
			boolean ncbi=false, swiss=false;
			for(RawHit hit : family.rawNcbi){
				if(!hit.key.equals(selected.key)){ncbi=true;}
			}
			for(RawHit hit : family.rawSwiss){
				if(!hit.key.equals(selected.key)){swiss=true;}
			}
			return ncbi && swiss ? "ncbi_and_swiss_product_text_disagreement" :
				ncbi ? "ncbi_product_text_disagreement" :
				swiss ? "swiss_product_text_disagreement" : "same_product_text";
		}

		final RawHit selected;
		final EggnogHomology eggnog;
		final String decision, cleanProduct, display, sourceMethod, support, coverage;
		final String score, productText, review, hitId;
		final boolean sourcePriority;
	}

	private static final class RawHit {
		RawHit(final String id_, final String method_, final String accession_,
				final String key_, final String display_, final String familyType_,
				final String score_, final String domainScore_, final String queryCov_,
				final String targetCov_, final boolean sequencePass_,
				final boolean domainPass_, final boolean coveragePass_,
				final boolean nativeEligible_, final boolean forNaming_,
				final String sortScore_, final double minCoverage_,
				final boolean informative_){
			id=id_;
			method=method_;
			accession=accession_;
			key=key_;
			display=display_;
			familyType=familyType_;
			score=score_;
			domainScore=domainScore_;
			queryCov=Double.parseDouble(queryCov_);
			targetCov=Double.parseDouble(targetCov_);
			queryLow=queryCov<minCoverage_;
			targetLow=targetCov<minCoverage_;
			sortScore=Double.parseDouble(sortScore_);
			sequencePass=sequencePass_;
			domainPass=domainPass_;
			coveragePass=coveragePass_;
			nativeEligible=nativeEligible_;
			forNaming=forNaming_;
			informative=informative_;
			source=id.startsWith("ncbi:") ? "ncbi" : "swiss";
		}

		boolean betterThan(final RawHit b){
			if(b==null){return true;}
			final double aScore=scoreValue(sortScore), bScore=scoreValue(b.sortScore);
			if(aScore!=bScore){return aScore>bScore;}
			if(coveragePass!=b.coveragePass){return coveragePass;}
			if(sequencePass!=b.sequencePass){return sequencePass;}
			if(domainPass!=b.domainPass){return domainPass;}
			return accession.compareTo(b.accession)<0;
		}

		String supportTag(){
			if(!informative){return source.equals("ncbi") ? "native_placeholder" :
				"swiss_placeholder";}
			return source.equals("ncbi") ?
				(nativeEligible ? "native_full" :
					"native_partial") :
				(coveragePass ? "swiss_full" :
					"swiss_partial");
		}

		String coverageTag(){
			if(coveragePass){return "coverage_full";}
			return queryLow && targetLow ? "query_and_target_low" :
				queryLow ? "query_low" : targetLow ? "target_low" : "coverage_low";
		}

		String scoreTag(){
			if(source.equals("swiss")){return "retained_search_hit";}
			return !sequencePass && !domainPass ? "sequence_and_domain_score_low" :
				!sequencePass ? "sequence_score_low" :
				!domainPass ? "domain_score_low" : "score_passed";
		}

		private static double scoreValue(final double x){
			return Double.isNaN(x) ? Double.NEGATIVE_INFINITY : x;
		}

		final String id, source, method, accession, key, display, familyType;
		final String score, domainScore;
		final double queryCov, targetCov, sortScore;
		final boolean sequencePass, domainPass, coveragePass, nativeEligible, forNaming;
		final boolean queryLow, targetLow;
		final boolean informative;
	}

	private static final class EggnogHomology {
		EggnogHomology(final LineParser1 lp, final String path, final long lineNumber){
			seedOrtholog=lp.parseByteArray(3);
			selectedSeedNumeric=lp.parseByteArray(4);
			seedIdentity=lp.parseByteArray(5);
			seedQueryCov=lp.parseByteArray(6);
			seedSubjectCov=lp.parseByteArray(7);
			seedEvalue=lp.parseByteArray(8);
			seedScore=lp.parseByteArray(9);
			pfamRoot=lp.parseByteArray(10);
			pfamRaw=lp.parseByteArray(11);
			seedDbPfam=lp.parseByteArray(12);
			annotationPfam=lp.parseByteArray(13);
			eggnogOgs=lp.parseByteArray(14);
			preferredName=lp.parseByteArray(15);
			kos=lp.parseByteArray(16);
			ec=lp.parseByteArray(17);
			koDefinition=lp.parseByteArray(18);
			functionalLead=lp.parseByteArray(19);
			taxCeiling=lp.parseByteArray(20);
			farthestDonorLineage=lp.parseByteArray(21);
			annotationConfidence=lp.parseByteArray(22);
			spanBin=lp.parseString(23);
			validate(path, lineNumber);
		}

		private void validate(final String path, final long lineNumber){
			try{
				Long.parseLong(ascii(selectedSeedNumeric));
			}catch(NumberFormatException nfe){
				throw new IllegalArgumentException("Bad selected_seed_numeric in "
					+path+" line "+lineNumber+": "+ascii(selectedSeedNumeric), nfe);
			}
			final double identity=parseEggnogDouble(seedIdentity, path, lineNumber,
				"seed_identity");
			final double queryCov=parseEggnogDouble(seedQueryCov, path, lineNumber,
				"seed_query_cov");
			final double subjectCov=parseEggnogDouble(seedSubjectCov, path, lineNumber,
				"seed_subject_cov");
			require(identity>=0 && identity<=100,
				"Bad seed_identity in "+path+" line "+lineNumber+": "+ascii(seedIdentity));
			require(queryCov>=0 && queryCov<=100,
				"Bad seed_query_cov in "+path+" line "+lineNumber+": "+ascii(seedQueryCov));
			require(subjectCov>=0 && subjectCov<=100,
				"Bad seed_subject_cov in "+path+" line "+lineNumber+": "
				+ascii(seedSubjectCov));
			require(parseEggnogDouble(seedEvalue, path, lineNumber, "seed_evalue")>=0,
				"Bad seed_evalue in "+path+" line "+lineNumber+": "+ascii(seedEvalue));
			parseEggnogDouble(seedScore, path, lineNumber, "seed_score");
			require(knownSpanBin(spanBin), "Bad eggNOG span_bin in "+path+" line "
				+lineNumber+": "+spanBin);
			final String pRoot=ascii(pfamRoot);
			final String rawText=ascii(pfamRaw);
			final String seedText=ascii(seedDbPfam);
			final String annotationText=ascii(annotationPfam);
			final LinkedHashSet<String> expected=roots(pRoot);
			final LinkedHashSet<String> raw=coordinateRoots(pfamRaw, path, lineNumber,
				"pfam_raw");
			final LinkedHashSet<String> seedDb=coordinateRoots(seedDbPfam, path, lineNumber,
				"seed_db_pfam");
			final LinkedHashSet<String> annotation=coordinateRoots(annotationPfam, path,
				lineNumber, "annotation_pfam");
			final ArrayList<String> rawTerms=coordinateTerms(pfamRaw, path, lineNumber,
				"pfam_raw");
			final ArrayList<String> seedDbTerms=coordinateTerms(seedDbPfam, path,
				lineNumber, "seed_db_pfam");
			final ArrayList<String> annotationTerms=coordinateTerms(annotationPfam, path,
				lineNumber, "annotation_pfam");
			final LinkedHashSet<String> og=ogRoots(ascii(eggnogOgs));
			ineligibleReason=ineligibleReason(identity, queryCov, subjectCov, expected,
				raw, seedDb, annotation, og, rawText, seedText, annotationText,
				parseBoolean(functionalLead, path, lineNumber), functionalRawReason());
			eligible=ineligibleReason.equals(".");
			relaxedType=relaxedType(identity, queryCov, subjectCov, expected, raw, seedDb,
				annotation, rawTerms, seedDbTerms, annotationTerms, og,
				parseBoolean(functionalLead, path, lineNumber), functionalRawReason(),
				ascii(preferredName), ascii(kos), ascii(ec));
		}

		private static String ineligibleReason(final double identity, final double queryCov,
				final double subjectCov, final LinkedHashSet<String> expected,
				final LinkedHashSet<String> raw, final LinkedHashSet<String> seedDb,
				final LinkedHashSet<String> annotation, final LinkedHashSet<String> og,
				final String rawText, final String seedText, final String annotationText,
				final boolean functionLead, final String functionalRawReason){
			if(expected.isEmpty() || raw.isEmpty() || seedDb.isEmpty() || annotation.isEmpty()){
				return "absent_pfam";
			}
			if(og.isEmpty()){return "absent_og";}
			if(!rawText.equals(seedText) || !rawText.equals(annotationText)){
				return "pfam_raw_mismatch";
			}
			if(expected.size()!=1 || raw.size()!=1 || seedDb.size()!=1 ||
					annotation.size()!=1 || og.size()!=1){
				return "multi_domain";
			}
			final String root=expected.iterator().next();
			if(!raw.contains(root) || !seedDb.contains(root) || !annotation.contains(root) ||
					!og.contains(root)){
				return "pfam_og_mismatch";
			}
			if(functionLead){return "functional_lead";}
			if(functionalRawReason!=null){return functionalRawReason;}
			if(identity<60){return "identity_lt_60";}
			if(queryCov<80){return "query_cov_lt_80";}
			if(subjectCov<80){return "subject_cov_lt_80";}
			return ".";
		}

		private static String relaxedType(final double identity, final double queryCov,
				final double subjectCov, final LinkedHashSet<String> expected,
				final LinkedHashSet<String> raw, final LinkedHashSet<String> seedDb,
				final LinkedHashSet<String> annotation, final ArrayList<String> rawTerms,
				final ArrayList<String> seedDbTerms, final ArrayList<String> annotationTerms,
				final LinkedHashSet<String> og, final boolean functionLead,
				final String functionalRawReason, final String preferred,
				final String koText, final String ecText){
			if(expected.isEmpty() || raw.isEmpty() || seedDb.isEmpty() || annotation.isEmpty()){
				return null;
			}
			if(og.isEmpty()){return null;}
			if(!sameTerms(rawTerms, seedDbTerms) || !sameTerms(rawTerms, annotationTerms)){
				return null;
			}
			if(!sameRoots(expected, raw) || !sameRoots(expected, seedDb) ||
					!sameRoots(expected, annotation) || !sameRoots(expected, og)){
				return null;
			}
			if(functionLead || identity<60 || queryCov<80 || subjectCov<80){return null;}
			if(!blankFunctional(koText) || !blankFunctional(ecText)){return null;}
			if(functionalRawReason==null){return expected.size()>1 ? "multi" : null;}
			return expected.size()==1 && functionalRawReason.equals("preferred_name_functional_lead")
				&& unresolvedEggnogPreferredSymbol(preferred) ? "symbol" : null;
		}

		private String functionalRawReason(){
			final String root=ascii(pfamRoot);
			final String preferred=ascii(preferredName);
			if(informativeEggnogPreferred(preferred) &&
					!toLowerCaseAscii(preferred).equals(toLowerCaseAscii(root))){
				return "preferred_name_functional_lead";
			}
			if(!blankFunctional(ascii(ec))){return "ec_functional_lead";}
			final String koText=ascii(kos);
			if(blankFunctional(koText)){return null;}
			final String koDef=ascii(koDefinition);
			if(blankFunctional(koDef) || containsAscii(koDef, "unresolved") ||
					containsAscii(koDef, "missing")){
				return "ko_definition_unresolved";
			}
			final LinkedHashSet<String> resolved=new LinkedHashSet<String>();
			for(String rawMapping : koDef.split(",")){
				final int equals=rawMapping.indexOf('=');
				if(equals<=0){return "ko_definition_unresolved";}
				final String ko=rawMapping.substring(0, equals);
				final String def=rawMapping.substring(equals+1);
				if(blankFunctional(def) || containsAscii(def, "unresolved") ||
						containsAscii(def, "missing")){
					return "ko_definition_unresolved";
				}
				if(!explicitNonfunctionalKoDefinition(def)){
					return "ko_functional_lead";
				}
				resolved.add(ko);
			}
			for(String ko : koText.split(",")){
				if(!blankFunctional(ko) && !resolved.contains(ko)){
					return "ko_definition_unresolved";
				}
			}
			return null;
		}

		private static boolean blankFunctional(final String s){
			return s.length()<1 || s.equals(".") || s.equals("-");
		}

		private static boolean containsAscii(final String text, final String term){
			return text.toLowerCase().indexOf(term)>=0;
		}

		private static boolean explicitNonfunctionalKoDefinition(final String def){
			final String lower=toLowerCaseAscii(def);
			return lower.equals("uncharacterized protein") ||
				lower.equals("hypothetical protein") ||
				lower.equals("unknown protein") ||
				lower.equals("unknown");
		}

		private static boolean informativeEggnogPreferred(final String s){
			if(blankFunctional(s)){return false;}
			if(toLowerCaseAscii(s).equals("orf")){return false;}
			if(isLocPlaceholder(s)){return false;}
			if(isNumericLocusPlaceholder(s)){return false;}
			return isInformative(s);
		}

		private static boolean unresolvedEggnogPreferredSymbol(final String s){
			if(!informativeEggnogPreferred(s)){return false;}
			for(int i=0; i<s.length(); i++){
				final char c=s.charAt(i);
				if(!((c>='A' && c<='Z') || (c>='a' && c<='z') ||
						(c>='0' && c<='9') || c=='_')){
					return false;
				}
			}
			return true;
		}

		private static boolean isLocPlaceholder(final String s){
			if(!s.startsWith("LOC") || s.length()<4){return false;}
			for(int i=3; i<s.length(); i++){
				if(!Character.isDigit(s.charAt(i))){return false;}
			}
			return true;
		}

		private static boolean isNumericLocusPlaceholder(final String s){
			if(s.length()<1 || !Character.isDigit(s.charAt(0))){return false;}
			for(int i=1; i<s.length(); i++){
				final char c=s.charAt(i);
				if(!(Character.isDigit(c) || (c>='A' && c<='Z') ||
						(c>='a' && c<='z') || c=='_')){
					return false;
				}
			}
			return true;
		}

		private static LinkedHashSet<String> coordinateRoots(final byte[] raw,
				final String path, final long lineNumber, final String field){
			final LinkedHashSet<String> roots=new LinkedHashSet<String>();
			final String text=ascii(raw);
			if(text.length()<1 || text.equals(".") || text.equals("-")){return roots;}
			for(String term : text.split(",")){
				if(term.length()<1 || term.equals(".")){continue;}
				roots.add(rootOf(term, path, lineNumber, field));
			}
			return roots;
		}

		private static ArrayList<String> coordinateTerms(final byte[] raw,
				final String path, final long lineNumber, final String field){
			final ArrayList<String> terms=new ArrayList<String>();
			final String text=ascii(raw);
			if(text.length()<1 || text.equals(".") || text.equals("-")){return terms;}
			for(String term : text.split(",")){
				if(term.length()<1 || term.equals(".")){continue;}
				rootOf(term, path, lineNumber, field);
				terms.add(term);
			}
			Collections.sort(terms);
			return terms;
		}

		private static LinkedHashSet<String> ogRoots(final String raw){
			final LinkedHashSet<String> roots=new LinkedHashSet<String>();
			if(raw.length()<1 || raw.equals(".")){return roots;}
			for(String term : raw.split(",")){
				final int at=term.indexOf('@');
				final String root=at<0 ? term : term.substring(0, at);
				if(root.length()>0 && !root.equals("UNK")){roots.add(root);}
			}
			return roots;
		}

		private static String rootOf(final String raw, final String path,
				final long lineNumber, final String field){
			final int last=raw.lastIndexOf('_');
			final int prev=last<0 ? -1 : raw.lastIndexOf('_', last-1);
			require(prev>0 && last>prev+1 && last<raw.length()-1,
				"Bad coordinate-suffixed Pfam "+field+" in "+path+" line "
				+lineNumber+": "+raw);
			final int start=parsePositiveInt(raw.substring(prev+1, last), path,
				lineNumber, field+"_start");
			final int end=parsePositiveInt(raw.substring(last+1), path, lineNumber,
				field+"_end");
			require(start<=end, "Bad coordinate-suffixed Pfam "+field+" in "+path
				+" line "+lineNumber+": "+raw);
			return raw.substring(0, prev);
		}

		private static LinkedHashSet<String> roots(final String raw){
			final LinkedHashSet<String> roots=new LinkedHashSet<String>();
			if(raw.length()<1 || raw.equals(".")){return roots;}
			for(String term : raw.split(",")){
				if(term.length()>0 && !term.equals(".")){roots.add(term);}
			}
			return roots;
		}

		private static boolean sameRoots(final LinkedHashSet<String> a,
				final LinkedHashSet<String> b){
			return a.size()==b.size() && a.containsAll(b);
		}

		private static boolean sameTerms(final ArrayList<String> a,
				final ArrayList<String> b){
			if(a.size()!=b.size()){return false;}
			for(int i=0; i<a.size(); i++){
				if(!a.get(i).equals(b.get(i))){return false;}
			}
			return true;
		}

		private static boolean parseBoolean(final byte[] field, final String path,
				final long lineNumber){
			final String s=ascii(field);
			if(s.equals("t") || s.equals("true") || s.equals("1") || s.equals("Y")){
				return true;
			}
			if(s.equals("f") || s.equals("false") || s.equals("0") || s.equals("N")){
				return false;
			}
			throw new IllegalArgumentException("Bad functional_lead in "+path
				+" line "+lineNumber+": "+s);
		}

		private static int parsePositiveInt(final String s, final String path,
				final long lineNumber, final String field){
			final int x;
			try{
				x=Integer.parseInt(s);
			}catch(NumberFormatException nfe){
				throw new IllegalArgumentException("Bad "+field+" in "+path+" line "
					+lineNumber+": "+s, nfe);
			}
			require(x>0, "Bad "+field+" in "+path+" line "+lineNumber+": "+s);
			return x;
		}

		private static boolean knownSpanBin(final String s){
			return s.equals("seed_pfam_width_fraction_ge_0.80") ||
				s.equals("seed_pfam_width_fraction_0.50_0.80") ||
				s.equals("seed_pfam_width_fraction_lt_0.50") ||
				s.equals("seed_pfam_width_lt_80") ||
				s.equals("not_assessed") ||
				s.equals("not_applicable");
		}

		private static double parseEggnogDouble(final byte[] value, final String path,
				final long lineNumber, final String field){
			final String s=ascii(value);
			try{
				final double x=Double.parseDouble(s);
				require(Double.isFinite(x), "Bad "+field+" in "+path+" line "
					+lineNumber+": "+s);
				return x;
			}catch(NumberFormatException nfe){
				throw new IllegalArgumentException("Bad "+field+" in "+path
					+" line "+lineNumber+": "+s, nfe);
			}
		}

		boolean eligible(final boolean relaxed){
			return eligible || (relaxed && relaxedType!=null);
		}

		String productTextTag(final boolean relaxed){
			return !eligible && relaxed && relaxedType!=null ?
				relaxedType.equals("symbol") ? "pfam_og_root_match_symbol_ambiguity" :
				"pfam_og_multi_root_match_nonfunctional_lead" :
				"pfam_og_root_match_nonfunctional_lead";
		}

		String reviewTag(final boolean relaxed){
			return !eligible && relaxed && "symbol".equals(relaxedType) ?
				"unresolved_symbol_ambiguity" : "none";
		}

		String product(){
			return "putative "+ascii(pfamRoot).replace(',', '/')+" domain-containing protein";
		}
		String hitId(){return "eggnog:"+ascii(selectedSeedNumeric)+";"+ascii(pfamRoot);}

		ByteBuilder append(final ByteBuilder bb){
			return bb.append(seedOrtholog).tab().append(selectedSeedNumeric).tab()
				.append(seedIdentity).tab().append(seedQueryCov).tab()
				.append(seedSubjectCov).tab().append(seedEvalue).tab()
				.append(seedScore).tab().append(pfamRoot).tab().append(pfamRaw).tab()
				.append(seedDbPfam).tab().append(annotationPfam).tab()
				.append(eggnogOgs).tab().append(preferredName).tab()
				.append(kos).tab().append(ec).tab().append(koDefinition).tab()
				.append(functionalLead).tab().append(taxCeiling).tab()
				.append(farthestDonorLineage).tab().append(annotationConfidence).tab()
				.append(spanBin);
		}

		private static String ascii(final byte[] b){
			return new String(b, StandardCharsets.US_ASCII);
		}

		final byte[] seedOrtholog, selectedSeedNumeric, seedIdentity, seedQueryCov;
		final byte[] seedSubjectCov, seedEvalue, seedScore, pfamRoot, pfamRaw;
		final byte[] seedDbPfam, annotationPfam, eggnogOgs, preferredName, kos, ec;
		final byte[] koDefinition, functionalLead, taxCeiling, farthestDonorLineage;
		final byte[] annotationConfidence;
		final String spanBin;
		String ineligibleReason;
		String relaxedType;
		boolean eligible;
	}

	private final class CuratedBroadConflict {
		CuratedBroadConflict(final String ruleId_, final String source_,
				final String decision_, final String product_){
			ruleId=ruleId_;
			source=source_;
			decision=decision_;
			product=product_;
		}
		void add(final String source_, final String decision_, final String product_,
				final String conflictProduct, final long lineNumber){
			require(source.equals(source_) && decision.equals(decision_) &&
				product.equals(product_), "Inconsistent broad-conflict rule "+ruleId
				+" line "+lineNumber);
			final String key=comparisonKey(displayProduct(conflictProduct));
			require(key.length()>0, "Empty broad-conflict key line "+lineNumber);
			require(keys.putIfAbsent(key, conflictProduct)==null,
				"Duplicate broad-conflict product for "+ruleId+" line "+lineNumber
				+": "+conflictProduct);
		}
		String mapKey(){
			final ArrayList<String> sorted=new ArrayList<String>(keys.keySet());
			Collections.sort(sorted);
			return source+"\t"+Evidence.productSetKey(sorted);
		}
		final String ruleId, source, decision, product;
		final LinkedHashMap<String, String> keys=new LinkedHashMap<String, String>();
	}

	private final class CuratedBroadParent {
		CuratedBroadParent(final String parentId_, final String source_,
				final String decision_, final String product_){
			parentId=parentId_;
			source=source_;
			decision=decision_;
			product=product_;
		}
		void add(final String source_, final String decision_, final String product_,
				final String childProduct, final long lineNumber){
			require(source.equals(source_) && decision.equals(decision_) &&
				product.equals(product_), "Inconsistent broad-parent rule "+parentId
				+" line "+lineNumber);
			final String key=comparisonKey(displayProduct(childProduct));
			require(key.length()>0, "Empty broad-parent key line "+lineNumber);
			final String mapKey=source+"\t"+key;
			final CuratedBroadParent old=curatedBroadParents.putIfAbsent(mapKey, this);
			require(old==null, "Duplicate broad-parent product for source "+source
				+" line "+lineNumber+": "+childProduct);
			childCount++;
		}
		final String parentId, source, decision, product;
		int childCount;
	}

	private static final class NcbiModel {
		NcbiModel(final String accession_, final String key_, final String display_,
				final String original_, final String familyType_){
			accession=accession_;
			key=key_;
			display=display_;
			original=original_;
			familyType=familyType_;
			precedence=precedence(familyType);
			specificity=specificity(familyType);
			require(precedence>=0, "Unsupported NCBI naming family_type "+familyType
				+" for "+accession);
		}
		final String accession, key, display, original, familyType;
		final int precedence, specificity;
	}

	private static final class ReviewedEvidence {
		ReviewedEvidence(final LineParser1 lp){
			hitId=lp.parseByteArray(3);
			sourceMethods=lp.parseByteArray(4);
			support=lp.parseByteArray(5);
			coverage=lp.parseByteArray(6);
			score=lp.parseByteArray(7);
			productText=lp.parseByteArray(8);
			review=lp.parseByteArray(9);
			seedOrtholog=lp.parseByteArray(10);
			selectedSeedNumeric=lp.parseByteArray(11);
			seedIdentity=lp.parseByteArray(12);
			seedQueryCov=lp.parseByteArray(13);
			seedSubjectCov=lp.parseByteArray(14);
			seedEvalue=lp.parseByteArray(15);
			seedScore=lp.parseByteArray(16);
			eggnogOgs=lp.parseByteArray(17);
			taxCeiling=lp.parseByteArray(18);
			farthestDonorLineage=lp.parseByteArray(19);
			rawPreferredName=lp.parseByteArray(20);
			rawKeggKo=lp.parseByteArray(21);
			rawEc=lp.parseByteArray(22);
			rawPfams=lp.parseByteArray(23);
			annotationConfidence=lp.parseByteArray(24);
			primarySupport=lp.parseByteArray(25);
			competingSupport=lp.parseByteArray(26);
			warnings=lp.parseByteArray(27);
			reviewerNotes=lp.parseByteArray(28);
		}

		ByteBuilder appendLabel(final ByteBuilder bb, final String product){
			return bb.tab().append(product).tab().append(sourceMethods).tab()
				.append(support).tab().append(coverage).tab().append(score).tab()
				.append(productText).tab().append(review).tab().append(hitId);
		}

		void validate(final Family family, final String path, final long lineNumber){
			final String method=ascii(sourceMethods);
			if(method.equals("Swiss-Prot/DIAMOND")){
				validateSwissProt(family, path, lineNumber);
			}else if(method.equals("NCBIfam/HMMER")){
				validateNcbi(family, path, lineNumber);
			}else{
				validateEggnog(path, lineNumber);
			}
		}

		void validateSwissProt(final Family family, final String path, final long lineNumber){
			requireDot(seedOrtholog, "seed_ortholog", path, lineNumber);
			requireDot(selectedSeedNumeric, "selected_seed_numeric", path, lineNumber);
			requireDot(seedIdentity, "seed_identity", path, lineNumber);
			requireDot(seedQueryCov, "seed_query_cov", path, lineNumber);
			requireDot(seedSubjectCov, "seed_subject_cov", path, lineNumber);
			requireDot(seedEvalue, "seed_evalue", path, lineNumber);
			requireDot(seedScore, "seed_score", path, lineNumber);
			requireDot(eggnogOgs, "eggNOG_OGs", path, lineNumber);
			requireDot(taxCeiling, "tax_ceiling", path, lineNumber);
			requireDot(farthestDonorLineage, "farthest_donor_lineage", path, lineNumber);
			requireDot(rawPreferredName, "raw_Preferred_name", path, lineNumber);
			requireDot(rawKeggKo, "raw_KEGG_ko", path, lineNumber);
			requireDot(rawEc, "raw_EC", path, lineNumber);
			requireDot(rawPfams, "raw_PFAMs", path, lineNumber);
			requireDot(annotationConfidence, "annotation_confidence", path, lineNumber);
			final String id=ascii(hitId);
			require(id.startsWith("swiss:") && id.length()>6,
				"Bad Swiss-Prot reviewed-evidence hit_id in "+path+" line "+lineNumber
				+": "+id);
			final String accession=id.substring(6);
			require(family.hasSwissHit(accession),
				"Swiss-Prot reviewed-evidence hit_id lacks retained hit in "+path
				+" line "+lineNumber+": "+id);
			requireMeaningful(primarySupport, "primary_support", path, lineNumber);
			requireMeaningful(competingSupport, "competing_support", path, lineNumber);
			requireMeaningful(reviewerNotes, "reviewer_notes", path, lineNumber);
		}

		void validateNcbi(final Family family, final String path, final long lineNumber){
			requireDot(seedOrtholog, "seed_ortholog", path, lineNumber);
			requireDot(selectedSeedNumeric, "selected_seed_numeric", path, lineNumber);
			requireDot(seedIdentity, "seed_identity", path, lineNumber);
			requireDot(seedQueryCov, "seed_query_cov", path, lineNumber);
			requireDot(seedSubjectCov, "seed_subject_cov", path, lineNumber);
			requireDot(seedEvalue, "seed_evalue", path, lineNumber);
			requireDot(seedScore, "seed_score", path, lineNumber);
			requireDot(eggnogOgs, "eggNOG_OGs", path, lineNumber);
			requireDot(taxCeiling, "tax_ceiling", path, lineNumber);
			requireDot(farthestDonorLineage, "farthest_donor_lineage", path, lineNumber);
			requireDot(rawPreferredName, "raw_Preferred_name", path, lineNumber);
			requireDot(rawKeggKo, "raw_KEGG_ko", path, lineNumber);
			requireDot(rawEc, "raw_EC", path, lineNumber);
			requireDot(rawPfams, "raw_PFAMs", path, lineNumber);
			requireDot(annotationConfidence, "annotation_confidence", path, lineNumber);
			final String id=ascii(hitId);
			require(id.startsWith("ncbi:") && family.hasRawNcbiHit(id),
				"NCBIfam/HMMER reviewed-evidence hit_id lacks retained hit in "
				+path+" line "+lineNumber+": "+id);
			requireMeaningful(primarySupport, "primary_support", path, lineNumber);
			requireMeaningful(competingSupport, "competing_support", path, lineNumber);
			requireMeaningful(reviewerNotes, "reviewer_notes", path, lineNumber);
		}

		void validateEggnog(final String path, final long lineNumber){
			final double identity=parseReviewedDouble(seedIdentity, path, lineNumber,
				"seed_identity");
			final double queryCov=parseReviewedDouble(seedQueryCov, path, lineNumber,
				"seed_query_cov");
			final double subjectCov=parseReviewedDouble(seedSubjectCov, path, lineNumber,
				"seed_subject_cov");
			require(identity>=0 && identity<=100,
				"Bad seed_identity in "+path+" line "+lineNumber);
			require(queryCov>=0 && queryCov<=100,
				"Bad seed_query_cov in "+path+" line "+lineNumber);
			require(subjectCov>=0 && subjectCov<=100,
				"Bad seed_subject_cov in "+path+" line "+lineNumber);
			require(parseReviewedDouble(seedEvalue, path, lineNumber, "seed_evalue")>=0,
				"Bad seed_evalue in "+path+" line "+lineNumber);
			parseReviewedDouble(seedScore, path, lineNumber, "seed_score");
			try{
				Long.parseLong(new String(selectedSeedNumeric, StandardCharsets.US_ASCII));
			}catch(NumberFormatException nfe){
				throw new IllegalArgumentException("Bad selected_seed_numeric in "+path
					+" line "+lineNumber, nfe);
			}
		}

		private void requireDot(final byte[] field, final String name, final String path,
				final long lineNumber){
			require(field.length==1 && field[0]=='.',
				"Reviewed-evidence row carries eggNOG field "+name+" in "
				+path+" line "+lineNumber);
		}

		private void requireMeaningful(final byte[] field, final String name,
				final String path, final long lineNumber){
			boolean meaningful=false;
			for(byte b : field){
				if(!Character.isWhitespace((char)b) && b!='.'){meaningful=true;}
			}
			require(meaningful, "Reviewed-evidence row lacks "+name+" in "
				+path+" line "+lineNumber);
		}

		private String ascii(final byte[] field){
			return new String(field, StandardCharsets.US_ASCII);
		}

		ByteBuilder append(final ByteBuilder bb, final Family family){
			return bb.append(family.familyId).tab().append(family.activeIndex).tab()
				.append(family.repId).tab().append(hitId).tab().append(sourceMethods).tab()
				.append(support).tab().append(coverage).tab().append(score).tab()
				.append(productText).tab().append(review).tab().append(seedOrtholog).tab()
				.append(selectedSeedNumeric).tab().append(seedIdentity).tab()
				.append(seedQueryCov).tab().append(seedSubjectCov).tab()
				.append(seedEvalue).tab().append(seedScore).tab().append(eggnogOgs).tab()
				.append(taxCeiling).tab().append(farthestDonorLineage).tab()
				.append(rawPreferredName).tab().append(rawKeggKo).tab().append(rawEc).tab()
				.append(rawPfams).tab().append(annotationConfidence).tab()
				.append(primarySupport).tab().append(competingSupport).tab().append(warnings).tab()
				.append(reviewerNotes);
		}

		final byte[] hitId, sourceMethods, support, coverage, score, productText, review;
		final byte[] seedOrtholog, selectedSeedNumeric, seedIdentity, seedQueryCov;
		final byte[] seedSubjectCov, seedEvalue, seedScore, eggnogOgs, taxCeiling;
		final byte[] farthestDonorLineage, rawPreferredName, rawKeggKo, rawEc, rawPfams;
		final byte[] annotationConfidence, primarySupport, competingSupport, warnings;
		final byte[] reviewerNotes;
	}

	private static double parseReviewedDouble(final byte[] field, final String path,
			final long lineNumber, final String name){
		final String s=new String(field, StandardCharsets.US_ASCII);
		try{
			final double x=Double.parseDouble(s);
			require(Double.isFinite(x), "Bad "+name+" in "+path+" line "+lineNumber+": "+s);
			return x;
		}catch(NumberFormatException nfe){
			throw new IllegalArgumentException("Bad "+name+" in "+path
				+" line "+lineNumber+": "+s, nfe);
		}
	}

	private static int precedence(final String type){
		// *_domain subfamily/paralog ranks are provisional BBGAP policy for NCBIFAM rows.
		if(type.equals("exception")){return 77;}
		if(type.equals("equivalog") || type.equals("hypoth_equivalog") ||
				type.equals("equivalog_domain") || type.equals("hypoth_equivalog_domain")){
			return 70;
		}
		if(type.equals("PfamEq")){return 57;}
		if(type.equals("subfamily") || type.equals("subfamily_domain")){return 55;}
		if(type.equals("PfamAutoEq")){return 37;}
		if(type.equals("paralog") || type.equals("paralog_domain")){return 35;}
		if(type.equals("superfamily")){return 33;}
		if(type.equals("domain")){return 30;}
		return -1;
	}

	private static int specificity(final String type){
		if(type.equals("exception") || type.equals("equivalog") ||
				type.equals("equivalog_domain")){
			return SPECIFIC;
		}
		return BROAD;
	}

	private static final class Evidence {
		void add(final String accession, final String key, final String display,
				final String original){
			accessionList.add(accession);
			productMap.putIfAbsent(key, display);
			productOriginals.putIfAbsent(original, original);
		}
		int productCount(){return productMap.size();}
		String firstKey(){
			for(String product : productMap.keySet()){return product;}
			throw new AssertionError("No products");
		}
		String firstProduct(){
			for(String product : productMap.values()){return product;}
			throw new AssertionError("No products");
		}
		String accessions(){return join(accessionList, ',');}
		String products(){
			final ArrayList<String> products=new ArrayList<String>(productOriginals.keySet());
			return join(products, '|');
		}
		private static String join(final ArrayList<String> list, final char delimiter){
			final ByteBuilder bb=new ByteBuilder(256);
			for(int i=0; i<list.size(); i++){
				if(i>0){bb.append(delimiter);}
				bb.append(list.get(i));
			}
			return bb.toString();
		}
		private static String productSetKey(final ArrayList<String> keys){
			final ByteBuilder bb=new ByteBuilder(256);
			for(String key : keys){
				bb.append(key.length()).append(':').append(key);
			}
			return bb.toString();
		}
		ArrayList<String> sortedKeys(){
			final ArrayList<String> keys=new ArrayList<String>(productMap.keySet());
			Collections.sort(keys);
			return keys;
		}
		String productsForKeys(final ArrayList<String> keys){
			final ByteBuilder bb=new ByteBuilder(256);
			for(int i=0; i<keys.size(); i++){
				if(i>0){bb.append('|');}
				bb.append(productMap.get(keys.get(i)));
			}
			return bb.toString();
		}
		final ArrayList<String> accessionList=new ArrayList<String>();
		final LinkedHashMap<String, String> productMap=new LinkedHashMap<String, String>();
		final LinkedHashMap<String, String> productOriginals=new LinkedHashMap<String, String>();
	}

	private static final class ConflictGroup {
		ConflictGroup(final String source_, final String reason_, final int productCount_,
				final String normalizedProducts_, final String displayProducts_){
			source=source_;
			reason=reason_;
			productCount=productCount_;
			normalizedProducts=normalizedProducts_;
			displayProducts=displayProducts_;
		}
		void add(final int familyId){familyIds.add(familyId);}
		String familyIds(){
			final ByteBuilder bb=new ByteBuilder(128);
			for(int i=0; i<familyIds.size(); i++){
				if(i>0){bb.append(',');}
				bb.append(familyIds.get(i));
			}
			return bb.toString();
		}
		final String source, reason, normalizedProducts, displayProducts;
		final int productCount;
		final IntList familyIds=new IntList();
	}

	private final LinkedHashMap<Integer, Family> families=new LinkedHashMap<Integer, Family>();
	private final ArrayList<Family> ordered=new ArrayList<Family>(8537);
	private final HashMap<String, String> synonyms=new HashMap<String, String>();
	private final HashMap<String, CuratedBroadConflict> curatedBroadConflicts=
		new HashMap<String, CuratedBroadConflict>();
	private final HashMap<String, CuratedBroadParent> curatedBroadParents=
		new HashMap<String, CuratedBroadParent>();
	private String familiesFile, ncbiFile, ncbiNearQueryDomainsFile, swissFile, reviewedFile;
	private String reviewedEvidenceFile, reviewedEvidenceOutFile;
	private String eggnogHomologyFile, eggnogHomologyOutFile;
	private String synonymsFile, broadConflictRulesFile, broadParentRulesFile;
	private String outFile, evidenceFile, unresolvedFile, summaryFile, auditFile;
	private String conflictGroupsFile, memberCoverageFile, memberUnresolvedFile;
	private String labelVersion, identityPin, consensusPin;
	private String mode="exact";
	private long seed=1;
	private int maxAudit=20;
	private double minCoverage=0.80, minSwissIdentity=90;
	private long autoSelected, autoExact, reviewed, reviewedNamed, reviewedEmpty;
	private long autoSpecific, autoBroad, autoSwissOnly, autoWithSwiss, autoNcbiOnly;
	private long autoSubfamilyDomain, autoNcbiNearQuery, autoBroadConflict, autoBroadParent;
	private long unresolvedConflict, unresolvedInsufficient;
	private long provisionalLabelCount, provisionalInformative, provisionalPlaceholder;
	private long provisionalSourcePriority;
	private long provisionalEggnogHomology;
	private boolean overwrite=true, ignoreCase=false, swissOnlyFallback=false;
	private boolean ncbiNearQueryFallback=false, provisionalLabels=false;
	private boolean eggnogDomainRelaxed=false;

	private static final double NEAR_QUERY_MIN_COVERAGE=0.75;
	private static final String REVIEWED_EVIDENCE_HEADER=
		"family_id\tactive_index\trep_id\thit_id\tsource_methods\tsupport_tag\t"
		+"coverage_tag\tscore_tag\tproduct_text_tag\treview_tag\tseed_ortholog\t"
		+"selected_seed_numeric\tseed_identity\tseed_query_cov\tseed_subject_cov\t"
		+"seed_evalue\tseed_score\teggNOG_OGs\ttax_ceiling\t"
		+"farthest_donor_lineage\traw_Preferred_name\traw_KEGG_ko\traw_EC\t"
		+"raw_PFAMs\tannotation_confidence\tprimary_support\tcompeting_support\t"
		+"warnings\treviewer_notes";
	private static final String[] REVIEWED_EVIDENCE_COLUMNS=REVIEWED_EVIDENCE_HEADER.split("\t");
	private static final String EGGNOG_HOMOLOGY_HEADER=
		"family_id\tactive_index\trep_id\tseed_ortholog\tselected_seed_numeric\t"
		+"seed_identity\tseed_query_cov\tseed_subject_cov\tseed_evalue\tseed_score\t"
		+"pfam_root\tpfam_raw\tseed_db_pfam\tannotation_pfam\teggNOG_OGs\t"
		+"Preferred_name\tKOs\tEC\tKO_definition\tfunctional_lead\ttax_ceiling\t"
		+"farthest_donor_lineage\tannotation_confidence\tspan_bin";
	private static final String[] EGGNOG_HOMOLOGY_COLUMNS=
		EGGNOG_HOMOLOGY_HEADER.split("\t");
	private static final String EGGNOG_HOMOLOGY_EVIDENCE_HEADER=
		EGGNOG_HOMOLOGY_HEADER+"\tselected\toriginal_decision\toriginal_reason\t"
		+"ineligible_reason";
	private static final String NCBI_NEAR_QUERY_DOMAIN_HEADER=
		BbgapCoverageAudit.DOMAIN_HEADER+"\tfor_naming\tcoverage_failure";
	private static final int SPECIFIC=1, BROAD=2;
}
