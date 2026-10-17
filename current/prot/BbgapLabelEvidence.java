package prot;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Joins raw BBGAP label-pilot searches with retained family/model metadata.
 * @author Yelan
 */
public final class BbgapLabelEvidence {

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BbgapLabelEvidence.class, false);
		try{
			new BbgapLabelEvidence(pp.args).process();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	BbgapLabelEvidence(final String[] args){
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("in") || a.equals("root") || a.equals("pilot")){inDir=b;}
			else if(a.equals("out")){outDir=b;}
			else if(a.equals("mode")){mode=b.toLowerCase();}
			else if(a.equals("meta")){metaFile=b;}
			else if(a.equals("families")){familiesFile=b;}
			else if(a.equals("compact") || a.equals("compactdir")){compactDir=b;}
			else if(a.equals("tbl") || a.equals("targets")){tbloutFile=b;}
			else if(a.equals("dom") || a.equals("domains")){domtbloutFile=b;}
			else if(a.equals("swiss") || a.equals("swissprot")){swissFile=b;}
			else if(a.equals("expectedtargets")){expectedTargets=Parse.parseInt(b, 0);}
			else if(a.equals("expecteddomains")){expectedDomains=Parse.parseInt(b, 0);}
			else if(a.equals("expectedswissprot")){expectedSwiss=Parse.parseInt(b, 0);}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(outDir==null){outDir=inDir;}
	}

	void process(){
		resolveDefaults();
		if(mode.equals("compact") || mode.equals("both")){compact();}
		if(mode.equals("validate") || mode.equals("both")){validate();}
		if(!mode.equals("compact") && !mode.equals("validate") && !mode.equals("both")){
			throw new IllegalArgumentException("Unknown mode: "+mode);
		}
	}

	private void compact(){
		final HashMap<String, NcbifamModel> models=loadNcbifamModels(metaFile);
		final HashMap<String, PilotFamily> families=loadPilotFamilies(familiesFile);
		writeNcbiTargets(models, tbloutFile, outDir+"/pilot.ncbi.targets.compact.tsv");
		writeNcbiDomains(models, domtbloutFile, outDir+"/pilot.ncbi.domains.compact.tsv");
		writeSwiss(families, swissFile, outDir+"/pilot.swissprot.compact.tsv");
	}

	private void writeNcbiTargets(final HashMap<String, NcbifamModel> models, final String in,
			final String out){
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		final ByteStreamWriter bsw=writer(out);
		final ByteBuilder bb=new ByteBuilder(1<<14);
		final WhiteLine wl=new WhiteLine();
		long rows=0;
		bsw.println("family_id\tactive_index\tconsensus_id\tncbi_accession\thmm_name\t"
			+"full_evalue\tfull_score\tbest_domain_evalue\tbest_domain_score\treported_domains\t"
			+"sequence_cutoff\tdomain_cutoff\thmm_length\tfamily_type\tfor_naming\tsource\t"
			+"gene_symbol\tec_numbers\tproduct_name");
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				wl.set(line, 18, in, rows+1);
				final String accession=wl.token(3);
				final NcbifamModel model=requiredModel(models, accession, in, rows+1);
				final String desc=wl.description();
				bb.clear();
				appendDescValue(bb, desc, "family_id", in, rows+1).tab();
				appendDescValue(bb, desc, "active_index", in, rows+1).tab();
				wl.appendToken(bb, 0).tab();
				bb.append(accession).tab();
				wl.appendToken(bb, 2).tab();
				wl.appendToken(bb, 4).tab();
				wl.appendToken(bb, 5).tab();
				wl.appendToken(bb, 7).tab();
				wl.appendToken(bb, 8).tab();
				wl.appendToken(bb, 16).tab();// HMMER 3.4 --tblout rep/nreported; token 15 is detected domains.
				model.appendTarget(bb).nl();
				bsw.print(bb);
				rows++;
			}
		}finally{
			close(bf, bsw, out);
		}
	}

	private void writeNcbiDomains(final HashMap<String, NcbifamModel> models, final String in,
			final String out){
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		final ByteStreamWriter bsw=writer(out);
		final ByteBuilder bb=new ByteBuilder(1<<14);
		final WhiteLine wl=new WhiteLine();
		long rows=0;
		bsw.println("family_id\tactive_index\tconsensus_id\tconsensus_len\tncbi_accession\t"
			+"hmm_name\thmm_len\tdomain_index\tdomain_count\tfull_evalue\tfull_score\t"
			+"conditional_evalue\tindependent_evalue\tdomain_score\thmm_from\thmm_to\t"
			+"ali_from\tali_to\tenv_from\tenv_to\thmm_cov\tconsensus_cov\tsequence_cutoff\t"
			+"domain_cutoff\tfamily_type\tsource\tgene_symbol\tproduct_name");
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				wl.set(line, 22, in, rows+1);
				final String accession=wl.token(4);
				final NcbifamModel model=requiredModel(models, accession, in, rows+1);
				final String desc=wl.description();
				final int consensusLen=wl.parseInt(2);
				final int hmmLen=wl.parseInt(5);
				final int hmmFrom=wl.parseInt(15);
				final int hmmTo=wl.parseInt(16);
				final int aliFrom=wl.parseInt(17);
				final int aliTo=wl.parseInt(18);
				final double hmmCov=(hmmTo-hmmFrom+1)/(double)hmmLen;
				final double consensusCov=(aliTo-aliFrom+1)/(double)consensusLen;
				require(coverageOk(hmmCov) && coverageOk(consensusCov),
					"Bad coverage derived from "+in+" line "+(rows+1));
				bb.clear();
				appendDescValue(bb, desc, "family_id", in, rows+1).tab();
				appendDescValue(bb, desc, "active_index", in, rows+1).tab();
				wl.appendToken(bb, 0).tab();
				wl.appendToken(bb, 2).tab();
				bb.append(accession).tab();
				wl.appendToken(bb, 3).tab();
				wl.appendToken(bb, 5).tab();
				wl.appendToken(bb, 9).tab();
				wl.appendToken(bb, 10).tab();
				wl.appendToken(bb, 6).tab();
				wl.appendToken(bb, 7).tab();
				wl.appendToken(bb, 11).tab();
				wl.appendToken(bb, 12).tab();
				wl.appendToken(bb, 13).tab();
				wl.appendToken(bb, 15).tab();
				wl.appendToken(bb, 16).tab();
				wl.appendToken(bb, 17).tab();
				wl.appendToken(bb, 18).tab();
				wl.appendToken(bb, 19).tab();
				wl.appendToken(bb, 20).tab();
				appendDecimal(bb, hmmCov).tab();
				appendDecimal(bb, consensusCov).tab();
				model.appendDomain(bb).nl();
				bsw.print(bb);
				rows++;
			}
		}finally{
			close(bf, bsw, out);
		}
	}

	private void writeSwiss(final HashMap<String, PilotFamily> families, final String in,
			final String out){
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		final ByteStreamWriter bsw=writer(out);
		final LineParser1 lp=new LineParser1('\t');
		final ByteBuilder bb=new ByteBuilder(1<<14);
		long lineNumber=0;
		bsw.println("family_id\tactive_index\tconsensus_id\tuniprot_accession\tuniprot_entry\t"
			+"identity\taligned_len\tquery_len\tsubject_len\tquery_cov\tsubject_cov\tqstart\t"
			+"qend\tsstart\tsend\tevalue\tbitscore\ttitle");
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0){continue;}
				lp.set(line);
				require(lp.terms()==13, "Expected 13 DIAMOND columns in "+in+" line "+lineNumber);
				final String consensus=lp.parseString(0);
				final PilotFamily family=families.get(consensus);
				require(family!=null, "Swiss-Prot row for unknown consensus "+consensus+" in "+in+" line "+lineNumber);
				final String[] uniprot=lp.parseString(1).split("\\|", -1);
				require(uniprot.length==3 && uniprot[1].length()>0 && uniprot[2].length()>0,
					"Bad UniProt sseqid in "+in+" line "+lineNumber+": "+lp.parseString(1));
				final int queryLen=lp.parseInt(4);
				final int subjectLen=lp.parseInt(5);
				final int qStart=lp.parseInt(6);
				final int qEnd=lp.parseInt(7);
				final int sStart=lp.parseInt(8);
				final int sEnd=lp.parseInt(9);
				final double qCov=(qEnd-qStart+1)/(double)queryLen;
				final double sCov=(sEnd-sStart+1)/(double)subjectLen;
				require(coverageOk(qCov) && coverageOk(sCov), "Bad DIAMOND coverage in "+in+" line "+lineNumber);
				bb.clear();
				bb.append(family.familyId).tab();
				bb.append(family.activeIndex).tab();
				bb.append(consensus).tab();
				bb.append(uniprot[1]).tab();
				bb.append(uniprot[2]).tab();
				lp.appendTerm(bb, 2).tab();
				lp.appendTerm(bb, 3).tab();
				lp.appendTerm(bb, 4).tab();
				lp.appendTerm(bb, 5).tab();
				appendDecimal(bb, qCov).tab();
				appendDecimal(bb, sCov).tab();
				lp.appendTerm(bb, 6).tab();
				lp.appendTerm(bb, 7).tab();
				lp.appendTerm(bb, 8).tab();
				lp.appendTerm(bb, 9).tab();
				lp.appendTerm(bb, 10).tab();
				lp.appendTerm(bb, 11).tab();
				lp.appendTerm(bb, 12).nl();
				bsw.print(bb);
			}
		}finally{
			close(bf, bsw, out);
		}
	}

	private void validate(){
		final Validation v=new Validation();
		v.rawTargets=countHmmerRows(tbloutFile);
		v.rawDomains=countHmmerRows(domtbloutFile);
		v.rawSwiss=countRows(swissFile);
		validateTargets(v, compactDir+"/pilot.ncbi.targets.compact.tsv");
		validateDomains(v, compactDir+"/pilot.ncbi.domains.compact.tsv");
		validateSwiss(v, compactDir+"/pilot.swissprot.compact.tsv");
		writeValidation(v, outDir+"/compact_validation.tsv");
	}

	private void validateTargets(final Validation v, final String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				lp.set(line);
				if(lineNumber==1){
					require(lp.terms()==19, "Unexpected target compact header width in "+path);
					continue;
				}
				require(lp.terms()==19, "Unexpected target compact width in "+path+" line "+lineNumber);
				v.compactTargets++;
				addFamily(v.targetFamilies, lp, path, lineNumber);
				if(lp.parseString(10).length()<1 || lp.parseString(11).length()<1){v.missingCutoff++;}
				if(lp.parseDouble(6)>=lp.parseDouble(10)){v.targetSeqPass++;}
				if(lp.parseDouble(8)>=lp.parseDouble(11)){v.targetDomPass++;}
			}
		}finally{
			bf.close();
		}
	}

	private void validateDomains(final Validation v, final String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				lp.set(line);
				if(lineNumber==1){
					require(lp.terms()==28, "Unexpected domain compact header width in "+path);
					continue;
				}
				require(lp.terms()==28, "Unexpected domain compact width in "+path+" line "+lineNumber);
				v.compactDomains++;
				addFamily(v.domainFamilies, lp, path, lineNumber);
				final double hmmCov=lp.parseDouble(20);
				final double consensusCov=lp.parseDouble(21);
				if(!coverageOk(hmmCov) || !coverageOk(consensusCov)){v.badCoverage++;}
				if(lp.parseString(22).length()<1 || lp.parseString(23).length()<1){v.missingCutoff++;}
				if(lp.parseDouble(10)>=lp.parseDouble(22)){v.domainSeqPass++;}
				if(lp.parseDouble(13)>=lp.parseDouble(23)){v.domainDomPass++;}
			}
		}finally{
			bf.close();
		}
	}

	private void validateSwiss(final Validation v, final String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				lp.set(line);
				if(lineNumber==1){
					require(lp.terms()==18, "Unexpected Swiss-Prot compact header width in "+path);
					continue;
				}
				require(lp.terms()==18, "Unexpected Swiss-Prot compact width in "+path+" line "+lineNumber);
				v.compactSwiss++;
				addFamily(v.swissFamilies, lp, path, lineNumber);
				final double queryCov=lp.parseDouble(9);
				final double subjectCov=lp.parseDouble(10);
				if(!coverageOk(queryCov) || !coverageOk(subjectCov)){v.badCoverage++;}
			}
		}finally{
			bf.close();
		}
	}

	private void writeValidation(final Validation v, final String out){
		final ByteStreamWriter bsw=writer(out);
		int failures=0;
		try{
			final ByteBuilder bb=new ByteBuilder(1<<12);
			bsw.println("metric\tvalue\texpected\tstatus");
			failures+=validationRow(bsw, bb, "raw_ncbi_target_rows", v.rawTargets, expectedTargets);
			failures+=validationRow(bsw, bb, "compact_ncbi_target_rows", v.compactTargets, v.rawTargets);
			failures+=validationRow(bsw, bb, "raw_ncbi_domain_rows", v.rawDomains, expectedDomains);
			failures+=validationRow(bsw, bb, "compact_ncbi_domain_rows", v.compactDomains, v.rawDomains);
			failures+=validationRow(bsw, bb, "raw_swissprot_rows", v.rawSwiss, expectedSwiss);
			failures+=validationRow(bsw, bb, "compact_swissprot_rows", v.compactSwiss, v.rawSwiss);
			failures+=recordedRow(bsw, bb, "unique_ncbi_target_families", v.targetFamilies.size());
			failures+=recordedRow(bsw, bb, "unique_ncbi_domain_families", v.domainFamilies.size());
			failures+=recordedRow(bsw, bb, "unique_swissprot_families", v.swissFamilies.size());
			failures+=validationRow(bsw, bb, "missing_bbgap_ids", v.missingIds, 0);
			failures+=validationRow(bsw, bb, "missing_ncbi_cutoffs", v.missingCutoff, 0);
			failures+=validationRow(bsw, bb, "bad_coverage_values", v.badCoverage, 0);
			failures+=recordedRow(bsw, bb, "ncbi_target_sequence_cutoff_pass_rows", v.targetSeqPass);
			failures+=recordedRow(bsw, bb, "ncbi_target_domain_cutoff_pass_rows", v.targetDomPass);
			failures+=recordedRow(bsw, bb, "ncbi_domain_sequence_cutoff_pass_rows", v.domainSeqPass);
			failures+=recordedRow(bsw, bb, "ncbi_domain_domain_cutoff_pass_rows", v.domainDomPass);
		}finally{
			if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+out);}
		}
		if(failures>0){throw new IllegalStateException("Validation failed with "+failures+" failing metric(s): "+out);}
	}

	private static int recordedRow(final ByteStreamWriter bsw, final ByteBuilder bb, final String metric, final long value){
		bb.clear().append(metric).tab().append(value).append("\trecorded\tOK");
		bsw.println(bb);
		return 0;
	}

	private static int validationRow(final ByteStreamWriter bsw, final ByteBuilder bb, final String metric,
			final long observed, final long expected){
		final boolean ok=expected<0 || observed==expected;
		bb.clear().append(metric).tab().append(observed).tab();
		if(expected<0){bb.append("recorded");}
		else{bb.append(expected);}
		bb.tab().append(ok ? "OK" : "FAIL");
		bsw.println(bb);
		return ok ? 0 : 1;
	}

	private static void addFamily(final HashSet<String> set, final LineParser1 lp, final String path,
			final long lineNumber){
		final String family=lp.parseString(0);
		final String active=lp.parseString(1);
		if(family.length()<1 || active.length()<1){throw new IllegalStateException("Missing BBGAP ids in "+path+" line "+lineNumber);}
		set.add(family);
	}

	private static HashMap<String, NcbifamModel> loadNcbifamModels(final String path){
		final HashMap<String, NcbifamModel> map=new HashMap<String, NcbifamModel>(1<<15);
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				require(lp.terms()>=21, "Too few NCBI metadata columns in "+path+" line "+lineNumber);
				final NcbifamModel model=new NcbifamModel(
					required(lp, 0, path, lineNumber, "ncbi_accession"),
					requiredDouble(lp, 3, path, lineNumber, "sequence_cutoff"),
					requiredDouble(lp, 4, path, lineNumber, "domain_cutoff"),
					requiredInt(lp, 5, path, lineNumber, "hmm_length"),
					required(lp, 6, path, lineNumber, "family_type"),
					lp.parseString(8), lp.parseString(10), lp.parseString(11),
					lp.parseString(13), lp.parseString(20));
				require(map.put(model.accession, model)==null, "Duplicate NCBI accession "+model.accession);
			}
		}finally{
			bf.close();
		}
		require(!map.isEmpty(), "No NCBI metadata rows loaded from "+path);
		return map;
	}

	private static HashMap<String, PilotFamily> loadPilotFamilies(final String path){
		final HashMap<String, PilotFamily> map=new HashMap<String, PilotFamily>(128);
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lineNumber==1){
					require(lp.parseString(0).equals("active_index") && lp.parseString(3).equals("rep_id"),
						"Unexpected pilot family header in "+path);
					continue;
				}
				require(lp.terms()>=4, "Too few pilot family columns in "+path+" line "+lineNumber);
				final PilotFamily family=new PilotFamily(lp.parseInt(0), lp.parseInt(1), lp.parseString(3));
				require(map.put(family.repId, family)==null, "Duplicate pilot rep_id "+family.repId);
			}
		}finally{
			bf.close();
		}
		require(!map.isEmpty(), "No pilot families loaded from "+path);
		return map;
	}

	private static NcbifamModel requiredModel(final HashMap<String, NcbifamModel> models,
			final String accession, final String file, final long row){
		final NcbifamModel model=models.get(accession);
		require(model!=null, "Missing NCBI metadata for "+accession+" in "+file+" line "+row);
		return model;
	}

	private static String required(final LineParser1 lp, final int col, final String path,
			final long lineNumber, final String name){
		final String value=lp.parseString(col);
		require(value.length()>0, "Missing "+name+" in "+path+" line "+lineNumber);
		return value;
	}

	private static int requiredInt(final LineParser1 lp, final int col, final String path,
			final long lineNumber, final String name){
		final String value=required(lp, col, path, lineNumber, name);
		try{return Parse.parseInt(value, 0);}
		catch(Throwable t){throw new IllegalArgumentException("Bad integer "+name+" in "+path+" line "+lineNumber+": "+value, t);}
	}

	private static double requiredDouble(final LineParser1 lp, final int col, final String path,
			final long lineNumber, final String name){
		final String value=required(lp, col, path, lineNumber, name);
		try{return Double.parseDouble(value);}
		catch(Throwable t){throw new IllegalArgumentException("Bad score "+name+" in "+path+" line "+lineNumber+": "+value, t);}
	}

	private static long countHmmerRows(final String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		long rows=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length>0 && line[0]!='#'){rows++;}
			}
		}finally{
			bf.close();
		}
		return rows;
	}

	private static long countRows(final String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		long rows=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length>0){rows++;}
			}
		}finally{
			bf.close();
		}
		return rows;
	}

	private static ByteStreamWriter writer(final String path){
		final File parent=new File(path).getParentFile();
		if(parent!=null){parent.mkdirs();}
		final ByteStreamWriter bsw=new ByteStreamWriter(path, true, false, true);
		bsw.start();
		return bsw;
	}

	private static void close(final ByteFile bf, final ByteStreamWriter bsw, final String out){
		try{
			bf.close();
		}finally{
			if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+out);}
		}
	}

	private void resolveDefaults(){
		if(inDir==null){
			inDir="/clusterfs/jgi/groups/gentech/homes/bbushnell/Keqing/bbgap/data/labels/pilot_80";
		}
		if(outDir==null){outDir=inDir;}
		if(compactDir==null){compactDir=outDir;}
		if(metaFile==null){metaFile=inDir+"/hmm_PGAP.tsv";}
		if(familiesFile==null){familiesFile=inDir+"/pilot_families.tsv";}
		if(tbloutFile==null){tbloutFile=inDir+"/pilot.ncbi.tblout";}
		if(domtbloutFile==null){domtbloutFile=inDir+"/pilot.ncbi.domtblout";}
		if(swissFile==null){swissFile=inDir+"/pilot.swissprot.tsv";}
	}

	private static boolean coverageOk(final double x){return x>=0 && x<=1 && Double.isFinite(x);}

	private static ByteBuilder appendDescValue(final ByteBuilder bb, final String desc, final String key,
			final String file, final long row){
		final String value=descValue(desc, key);
		require(value.length()>0, "Missing "+key+" in HMMER description from "+file+" line "+row);
		return bb.append(value);
	}

	private static String descValue(final String desc, final String key){
		final String prefix=key+"=";
		int start=0;
		while(start<desc.length()){
			while(start<desc.length() && desc.charAt(start)<=' '){start++;}
			int stop=start;
			while(stop<desc.length() && desc.charAt(stop)>' '){stop++;}
			if(desc.startsWith(prefix, start)){return desc.substring(start+prefix.length(), stop);}
			start=stop+1;
		}
		return "";
	}

	private static ByteBuilder appendDecimal(final ByteBuilder bb, final double value){
		return bb.append(String.format(java.util.Locale.ROOT, "%.12f", value));
	}

	private static void require(final boolean ok, final String message){
		if(!ok){throw new IllegalArgumentException(message);}
	}

	private static final class NcbifamModel {
		NcbifamModel(final String accession_, final double sequenceCutoff_, final double domainCutoff_,
				final int hmmLength_, final String familyType_, final String forNaming_,
				final String product_, final String gene_, final String ec_, final String source_){
			accession=accession_;
			sequenceCutoff=sequenceCutoff_;
			domainCutoff=domainCutoff_;
			hmmLength=hmmLength_;
			familyType=familyType_;
			forNaming=forNaming_;
			product=product_;
			gene=gene_;
			ec=ec_;
			source=source_;
		}
		ByteBuilder appendTarget(final ByteBuilder bb){
			return appendShared(bb).tab().append(forNaming).tab().append(source).tab()
				.append(gene).tab().append(ec).tab().append(product);
		}
		ByteBuilder appendDomain(final ByteBuilder bb){
			return bb.append(sequenceCutoff, 3, true).tab().append(domainCutoff, 3, true)
				.tab().append(familyType).tab().append(source).tab().append(gene).tab().append(product);
		}
		private ByteBuilder appendShared(final ByteBuilder bb){
			return bb.append(sequenceCutoff, 3, true).tab().append(domainCutoff, 3, true).tab()
				.append(hmmLength).tab().append(familyType);
		}
		final String accession, familyType, forNaming, product, gene, ec, source;
		final double sequenceCutoff, domainCutoff;
		final int hmmLength;
	}

	private static final class PilotFamily {
		PilotFamily(final int activeIndex_, final int familyId_, final String repId_){
			activeIndex=activeIndex_;
			familyId=familyId_;
			repId=repId_;
		}
		final int activeIndex, familyId;
		final String repId;
	}

	private static final class Validation {
		long rawTargets, rawDomains, rawSwiss;
		long compactTargets, compactDomains, compactSwiss;
		long missingIds, missingCutoff, badCoverage;
		long targetSeqPass, targetDomPass, domainSeqPass, domainDomPass;
		final HashSet<String> targetFamilies=new HashSet<String>();
		final HashSet<String> domainFamilies=new HashSet<String>();
		final HashSet<String> swissFamilies=new HashSet<String>();
	}

	private static final class WhiteLine {
		WhiteLine set(final byte[] line_, final int tokens, final String file, final long row){
			line=line_;
			ensure(tokens);
			count=0;
			int pos=0;
			while(count<tokens && pos<line.length){
				while(pos<line.length && line[pos]<=' '){pos++;}
				if(pos>=line.length){break;}
				starts[count]=pos;
				while(pos<line.length && line[pos]>' '){pos++;}
				stops[count]=pos;
				count++;
			}
			while(pos<line.length && line[pos]<=' '){pos++;}
			descStart=pos;
			require(count==tokens, "Expected at least "+tokens+" HMMER columns in "+file+" line "+row);
			return this;
		}
		String token(final int term){return new String(line, starts[term], stops[term]-starts[term]);}
		int parseInt(final int term){return Parse.parseInt(line, starts[term], stops[term]);}
		ByteBuilder appendToken(final ByteBuilder bb, final int term){
			return bb.append(line, starts[term], stops[term]-starts[term]);
		}
		String description(){
			return descStart>=line.length ? "" : new String(line, descStart, line.length-descStart);
		}
		private void ensure(final int len){
			if(starts.length<len){
				starts=new int[len];
				stops=new int[len];
			}
		}
		private byte[] line;
		private int[] starts=new int[24], stops=new int[24];
		private int count, descStart;
	}

	private String inDir, outDir, compactDir, metaFile, familiesFile, tbloutFile, domtbloutFile, swissFile;
	private String mode="both";
	private long expectedTargets=-1, expectedDomains=-1, expectedSwiss=-1;
}
