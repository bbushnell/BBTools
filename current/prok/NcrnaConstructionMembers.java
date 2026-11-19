package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import map.IntHashMap;
import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

/** Resolves FMN member admission from saved original, native and assembly evidence.
 * Keeps original members and native events in separate, complete provenance tables.
 * One unambiguous, physically corresponding clean FMN core may replace old bounds;
 * weak/multiple-core associations preserve the old prior for review. Unsupported old members
 * remain labelled priors. Assembly exclusion is independent of admission.
 * No sequence payload, gene-identity filter or automatic rare-lineage exception.
 * Outputs require review and extraction validation before resource construction.
 * @author Ganyu
 */
public final class NcrnaConstructionMembers{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionMembers.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String genomes=null;
			String bindings=null;
			String evidence=null;
			String groups=null;
			String raw=null;
			String imported=null;
			String models=null;
			String taxonomy=null;
			int crossCount=20;
			int batchCount=-1;
			long rawCount=-1;
			String mode="members";
			int coreCount=-1;
			int oldCount=-1;
			int genomeCount=-1;
			int evalCount=-1;
			int newCount=-1;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("genomes")){genomes=b;}
				else if(a.equals("mode")){mode=b;}
				else if(a.equals("groups")){groups=b;}
				else if(a.equals("raw")){raw=b;}
				else if(a.equals("imported")){imported=b;}
				else if(a.equals("models")){models=b;}
				else if(a.equals("taxonomy")){taxonomy=b;}
				else if(a.equals("crosscount")){crossCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedbatches")){batchCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedraw")){rawCount=Parse.parseKMG(b);}
				else if(a.equals("expectedcores")){coreCount=Parse.parseIntKMG(b);}
				else if(a.equals("bindings")){bindings=b;}
				else if(a.equals("evidence")){evidence=b;}
				else if(a.equals("expectedmembers")){oldCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedgenomes")){genomeCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedeval")){evalCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectednew")){newCount=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown construction-member parameter: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			require(assertions && parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use -ea in=groups.inputs.tsv.gz genomes= bindings=prefix evidence=prefix out=fresh_prefix; no sampling or overwrite");
			if("procedure-report".equals(mode)){
				require(bindings==null && evidence==null && genomes==null && groups==null && taxonomy==null && imported==null
					&& models!=null && batchCount==-1 && rawCount==-1, "Procedure comparison consumes selected inputs and both completed native procedures");
				procedureReport(parser.in1, raw, models, parser.out1, crossCount);
				pp.outstream.println("NCRNA_FMN_PROCEDURE_COMPARISON_COMPLETE_NOT_SPLIT");
			}else if("cross-report".equals(mode)){
				require(bindings==null && evidence==null && genomes==null && groups==null && taxonomy==null && imported==null
					&& models!=null && batchCount==-1 && rawCount==-1, "Cross report reads selected inputs, expected cores, native collection and the original model roster");
				crossReport(parser.in1, raw, models, parser.out1, crossCount);
				pp.outstream.println("NCRNA_FMN_CROSS_REPORT_COMPLETE_NO_SPLIT_CHANGE");
			}else if("lineages".equals(mode)){
				require(bindings==null && evidence==null && raw==null && imported==null && models==null
					&& batchCount==-1 && rawCount==-1, "Lineage comparison reads saved metadata only");
				lineages(parser.in1, groups, genomes, taxonomy, parser.out1,
					genomeCount, oldCount, evalCount, newCount, coreCount);
				pp.outstream.println("NCRNA_FMN_LINEAGES_COMPLETE_NOT_SPLIT_FROZEN");
			}else if("cross-select".equals(mode)){
				require(bindings==null && evidence==null && raw==null && imported==null && models==null
					&& batchCount==-1 && rawCount==-1, "Cross selection consumes only saved member/group/genome metadata and accession-bound taxonomy");
				crossSelect(parser.in1, groups, genomes, taxonomy, parser.out1,
					genomeCount, oldCount, evalCount, newCount, coreCount, crossCount);
				pp.outstream.println("NCRNA_FMN_CROSS_SELECTION_COMPLETE_NOT_SEARCHED");
			}else if("stages".equals(mode)){
				require(taxonomy==null, "Full taxonomy is consumed only by cross-select");
				require(bindings==null, "Stage audit uses saved admission, not original bindings input");
				stageAudit(parser.in1, groups, genomes, evidence, imported, raw, models, parser.out1,
					genomeCount, oldCount, evalCount, newCount, coreCount, batchCount, rawCount);
				pp.outstream.println("NCRNA_CONSTRUCTION_STAGE_AUDIT_COMPLETE_NO_SPLIT_CHANGE");
			}else if("coverage".equals(mode)){
				require(taxonomy==null, "Full taxonomy is consumed only by cross-select");
				require(raw==null && imported==null && models==null && batchCount==-1 && rawCount==-1, "Raw-stage parameters require mode=stages");
				require(bindings==null && evidence==null, "Coverage mode consumes saved members; bindings/evidence would be ignored");
				coverageAudit(parser.in1, groups, genomes, parser.out1, genomeCount, oldCount, evalCount, newCount, coreCount);
				pp.outstream.println("NCRNA_CONSTRUCTION_COVERAGE_AUDIT_COMPLETE_NO_SPLIT_CHANGE");
			}else{
				require(taxonomy==null, "Full taxonomy is consumed only by cross-select");
				require(raw==null && imported==null && models==null && batchCount==-1 && rawCount==-1, "Raw-stage parameters require mode=stages");
				require("members".equals(mode) && groups==null && coreCount==-1, "Use mode=members or coverage with mode-specific inputs");
				process(parser.in1, genomes, bindings, evidence, parser.out1, genomeCount, oldCount, evalCount, newCount);
				pp.outstream.println("NCRNA_CONSTRUCTION_MEMBER_LEDGER_COMPLETE_REVIEW_AND_EXTRACTION_REQUIRED");
			}
		}finally{Shared.closeStream(pp.outstream);}
	}

	public static void process(final String groups, final String genomes, final String bindings, final String evidence,
			final String output, final int genomeCount, final int oldCount, final int evalCount, final int newCount){
		require(groups!=null && genomes!=null && bindings!=null && evidence!=null && output!=null
			&& genomeCount>0 && oldCount>0 && evalCount>0 && newCount>0,
			"Complete metadata inputs and positive whole-cohort counts are required");
		RrnaResourceIO.fresh(output, ".old.tsv.gz", ".events.tsv.gz", ".extract.tsv.gz", ".coverage.tsv", ".summary.tsv");
		final NcrnaConstructionMembers ledger=new NcrnaConstructionMembers();
		ledger.genomes(genomes, genomeCount);
		ledger.groups(groups, evalCount, newCount);
		ledger.members(bindings+".members.tsv.gz", oldCount);
		ledger.bindings(bindings+".bindings.tsv", genomeCount);
		ledger.expected(evidence+".summary.tsv");
		ledger.events(evidence+".fmn_events.tsv.gz");
		ledger.links(evidence+".links.tsv.gz");
		ledger.states(evidence+".old_evidence.tsv.gz");
		ledger.resolve();
		ledger.write(output);
	}

	private NcrnaConstructionMembers(){}

	/** Audit the saved admission result without redoing admission or changing labels. */
	public static void coverageAudit(final String members, final String groups, final String genomes, final String output,
			final int genomeCount, final int oldCount, final int evalCount, final int newCount, final int coreCount){
		require(members!=null && groups!=null && genomes!=null && output!=null && genomeCount>0
			&& oldCount>=0 && evalCount>0 && newCount>0 && coreCount>=0, "Coverage audit needs complete inputs and explicit conserved counts");
		RrnaResourceIO.fresh(output, ".coverage.tsv", ".unlabelled.tsv.gz", ".summary.tsv");
		final NcrnaConstructionMembers ledger=new NcrnaConstructionMembers();
		ledger.genomes(genomes, genomeCount);
		ledger.groups(groups, evalCount, newCount);
		ledger.auditBindings();
		ledger.auditCores(members+".extract.tsv.gz", coreCount);
		ledger.auditOld(members+".old.tsv.gz", oldCount);
		ledger.writeAudit(output, oldCount, coreCount);
	}

	private void auditBindings(){
		for(final Node node:nodes){
			final Meta meta=node.meta;
			if(meta==null){continue;}
			require(meta.auditNode==null || meta.auditNode.group==node.group && meta.auditNode.split.equals(node.split),
				"One literal original FASTA cannot have conflicting resolved group/split bindings");
			if(meta.auditNode==null || node.cohort<meta.auditNode.cohort
				|| node.cohort==meta.auditNode.cohort && node.id<meta.auditNode.id){meta.auditNode=node;}
			meta.auditDiscovery|=node.cohort==1;
		}
		for(final Meta meta:metadata){require(meta.auditNode!=null, "Coverage denominator lost an original genome, including zero-locus genomes: "+meta.genome);}
	}

	private void auditCores(final String path, final int expected){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, EXTRACT_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 21);
				final Node node=node(cohort(1), p.parseInt(2));
				final Meta meta=node.meta;
				require(ids.add(text(0)) && !text(0).isEmpty(), "Saved extraction identities must remain unique and nonempty");
				require(text(3).equals(node.fasta) && p.parseInt(4)==node.group && text(13).equals(node.split),
					"Saved extraction source, assembly group and split must agree with resolved inputs");
				require(text(5).equals(meta==null ? "NA" : meta.taxid) && text(6).equals(meta==null ? "UNKNOWN" : meta.domain)
					&& text(7).equals(meta==null ? "UNKNOWN" : meta.phylum), "Saved taxonomy must come from the same original literal-path roster; do not silently inherit labels");
				final boolean intact=flag(14);
				final boolean seed=flag(15);
				final boolean train=node.split.equals(TRAIN);
				require(seed && flag(16)==(train && intact) && flag(17)==train
					&& p.termEquals("PENDING_SEQUENCE_AND_FLANK_VALIDATION", 20), "Saved extraction roles or validation state changed before the audit");
				final boolean nativeCore=token(text(12), "CM_FULL_CLEAN");
				require(!nativeCore || intact, "A clean native full core must retain its intact metadata role");
				node.auditCores++;
				node.auditIntact|=intact;
				node.auditNative|=nativeCore;
				if(nativeCore){node.auditNativeCores++;}
				if(intact && !nativeCore){node.auditPriorCores++;}
			}
		}
		require(ids.size()==expected, "Coverage audit must conserve every saved candidate core");
	}

	private void auditOld(final String path, final int expected){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, OLD_OUTPUT_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 39);
				final Meta meta=byGenome.get(text(6));
				require(meta!=null && text(8).equals(meta.taxid) && text(9).equals(meta.domain) && text(10).equals(meta.phylum),
					"Old member taxonomy must agree with the same original genome roster");
				require(flag(0) && text(5).startsWith("RF00050_") && ids.add(text(5)), "Every original FMN member must occur exactly once");
				require(p.parseInt(34)==meta.auditNode.group && text(33).equals(meta.auditNode.split),
					"Old-member assembly split must agree with its resolved original-genome binding");
				meta.auditOld++;
			}
		}
		require(ids.size()==expected, "Coverage audit must conserve every original member, including review-only priors");
	}

	/** Select a small crossed-procedure diagnostic, never a new training split. */
	public static void crossSelect(final String members, final String groups, final String genomes, final String taxonomy,
			final String output, final int genomeCount, final int oldCount, final int evalCount, final int newCount,
			final int coreCount, final int selectedCount){
		require(members!=null && groups!=null && genomes!=null && taxonomy!=null && output!=null && selectedCount>0,
			"Cross selection requires saved populations, original accession taxonomy and a positive per-arm size");
		RrnaResourceIO.fresh(output, ".selected.tsv", ".composition.tsv", ".expected.tsv", ".summary.tsv",
			".evaluation.paths.tsv", ".discovery.paths.tsv");
		final NcrnaConstructionMembers ledger=new NcrnaConstructionMembers();
		ledger.genomes(genomes, genomeCount);
		ledger.groups(groups, evalCount, newCount);
		ledger.auditBindings();
		ledger.auditCores(members+".extract.tsv.gz", coreCount);
		ledger.auditOld(members+".old.tsv.gz", oldCount);
		ledger.crossTaxonomy(taxonomy);
		ledger.crossWrite(members, output, selectedCount);
	}

	/** Read the original GTDB accession table without allocating keys for unrelated rows. */
	private void crossTaxonomy(final String path){
		final ObjectMap<NcrnaReferenceGeometry.Target, CrossTax> wanted=
			new ObjectMap<NcrnaReferenceGeometry.Target, CrossTax>(NcrnaReferenceGeometry.Target.class, CrossTax.class);
		final NcrnaReferenceGeometry.Target probe=new NcrnaReferenceGeometry.Target();
		for(final Meta meta:metadata){
			if(!meta.phylum.equals("Bacteroidota")){continue;}
			require(meta.domain.equals("bacteria") && !meta.assembly.isEmpty(), "Bacteroidota taxonomy needs the original bacterial assembly accession");
			final byte[] accession=meta.assembly.getBytes(StandardCharsets.UTF_8);
			probe.set(accession, 0, accession.length);
			CrossTax tax=wanted.get(probe);
			if(tax==null){
				tax=new CrossTax();
				wanted.put(probe.copy(), tax);
			}
			meta.crossTax=tax;
		}
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 2);
				p.setBounds(0);
				int start=p.a();
				if(p.b()-start>3 && ((row[start]=='R' && row[start+1]=='S') || (row[start]=='G' && row[start+1]=='B')) && row[start+2]=='_'){start+=3;}
				probe.set(row, start, p.b());
				final CrossTax tax=wanted.get(probe);
				if(tax==null){continue;}
				require(!tax.seen, "Duplicate original GTDB accession would make class/order selection ambiguous");
				final String lineage=text(1);
				require(rank(lineage, "d__").equals("Bacteria") && rank(lineage, "p__").equals("Bacteroidota"),
					"Accession-bound original GTDB taxonomy must agree with the frozen construction domain/phylum");
				tax.clazz=rank(lineage, "c__");
				tax.order=rank(lineage, "o__");
				tax.family=rank(lineage, "f__");
				tax.genus=rank(lineage, "g__");
				tax.seen=true;
			}
		}
	}

	/** Compare the complete saved Bacteroidota populations, without selecting or changing a split. */
	public static void lineages(final String members, final String groups, final String genomes, final String taxonomy,
			final String output, final int genomeCount, final int oldCount, final int evalCount,
			final int newCount, final int coreCount){
		require(members!=null && groups!=null && genomes!=null && taxonomy!=null && output!=null
			&& genomeCount>0 && evalCount>0 && newCount>0 && oldCount>=0 && coreCount>=0,
			"Lineage comparison needs complete saved metadata, positive genome/cohort counts and nonnegative old/core counts");
		RrnaResourceIO.fresh(output, ".lineages.tsv", ".lineage_inputs.tsv.gz", ".summary.tsv");
		final NcrnaConstructionMembers ledger=new NcrnaConstructionMembers();
		ledger.genomes(genomes, genomeCount);
		ledger.groups(groups, evalCount, newCount);
		ledger.auditBindings();
		ledger.auditCores(members+".extract.tsv.gz", coreCount);
		ledger.auditOld(members+".old.tsv.gz", oldCount);
		ledger.crossTaxonomy(taxonomy);
		ledger.writeLineages(output);
	}

	private void writeLineages(final String prefix){
		final ObjectMap<String, LineageCategory> byRank=new ObjectMap<String, LineageCategory>(String.class, LineageCategory.class);
		final ObjectMap<String, String> groupRanks=new ObjectMap<String, String>(String.class, String.class);
		final ArrayList<LineageCategory> categories=new ArrayList<LineageCategory>();
		final IntHashMap positiveGroups=new IntHashMap();
		final int[] inputs=new int[2];
		final int[] missing=new int[2];
		final int[] unbound=new int[2];
		try(Output detail=new Output(prefix+".lineage_inputs.tsv.gz")){
			detail.writer.println("cohort\tinput_id\tfasta\tassembly_group\tassembly\tclass\torder\tfamily\tgenus\ttaxonomy_status\tCM_admitted_cores");
			for(final Node node:nodes){
				if(node.meta==null){
					unbound[node.cohort]++;
					continue;
				}
				final CrossTax tax=node.meta.crossTax;
				if(tax==null){continue;}
				inputs[node.cohort]++;
				final String status=tax.lineageStatus(true);
				if(!status.equals("GTDB_ACCESSION_MATCH")){missing[node.cohort]++;}
				if(node.cohort==0 && node.auditNative){positiveGroups.put(node.group, 1);}
				detail.writer.print(b.clear().append(NAMES[node.cohort]).tab().append(node.id).tab().append(node.fastaBytes)
					.tab().append(node.group).tab().append(node.meta.assembly.getBytes(StandardCharsets.UTF_8)).tab()
					.append(tax.clazz.getBytes(StandardCharsets.UTF_8)).tab().append(tax.order.getBytes(StandardCharsets.UTF_8)).tab()
					.append(tax.family.getBytes(StandardCharsets.UTF_8)).tab().append(tax.genus.getBytes(StandardCharsets.UTF_8)).tab()
					.append(status).tab().append(node.auditNativeCores).nl());
				for(int rank=0; rank<2; rank++){
					final String key=(rank==0 ? "FAMILY" : "GENUS")+'\t'+tax.family+'\t'
						+(rank==0 ? "NA" : tax.genus)+'\t'+tax.lineageStatus(rank==1);
					final String groupKey=rank+"\t"+node.group;
					final String previous=groupRanks.get(groupKey);
					require(previous==null || previous.equals(key), "One resolved assembly group cannot have conflicting GTDB lineage bindings: "+groupKey);
					groupRanks.put(groupKey, key);
					LineageCategory category=byRank.get(key);
					if(category==null){
						category=new LineageCategory(key, rank, tax.lineageStatus(rank==1).equals("GTDB_ACCESSION_MATCH"));
						byRank.put(key, category);
						categories.add(category);
					}
					if(node.cohort==0){
						category.evaluationInputs++;
						category.evaluationGroups.put(node.group, 1);
						category.cores+=node.auditNativeCores;
						if(node.auditNative){category.positive.put(node.group, 1);}
					}else{
						category.discoveryInputs++;
						category.discoveryGroups.put(node.group, 1);
					}
				}
			}
		}
		java.util.Collections.sort(categories, (a, z)->a.key.compareTo(z.key));
		final int[] absent=new int[2];
		try(Output table=new Output(prefix+".lineages.tsv"); Output summary=new Output(prefix+".summary.tsv")){
			table.writer.println("rank\tfamily\tgenus\ttaxonomy_status\tevaluation_inputs\tevaluation_groups\tevaluation_CM_positive_groups\tevaluation_CM_core_records\tdiscovery_inputs\tdiscovery_groups\tdiscovery_coverage");
			for(final LineageCategory category:categories){
				if(category.known && category.discoveryInputs==0){absent[category.rank]+=category.positive.size();}
				table.writer.print(b.clear().append(category.key.getBytes(StandardCharsets.UTF_8)).tab()
					.append(category.evaluationInputs).tab().append(category.evaluationGroups.size()).tab()
					.append(category.positive.size()).tab().append(category.cores).tab()
					.append(category.discoveryInputs).tab().append(category.discoveryGroups.size()).tab()
					.append(!category.known ? "UNKNOWN_TAXONOMY" : category.discoveryInputs==0 ? "NO_MATCHED_DISCOVERY_INPUT" : "DISCOVERY_PRESENT").nl());
			}
			summary.writer.println("evaluation_Bacteroidota_inputs\tdiscovery_Bacteroidota_inputs\tevaluation_CM_positive_groups\tevaluation_missing_family_or_genus\tdiscovery_missing_family_or_genus\tpositive_groups_no_matched_discovery_family\tpositive_groups_no_matched_discovery_genus\tall_evaluation_unbound_inputs\tall_discovery_unbound_inputs\tstatus");
			summary.writer.print(b.clear().append(inputs[0]).tab().append(inputs[1]).tab().append(positiveGroups.size()).tab()
				.append(missing[0]).tab().append(missing[1]).tab().append(absent[0]).tab().append(absent[1]).tab()
				.append(unbound[0]).tab().append(unbound[1]).append("\tSAVED_METADATA_NOT_BIOLOGICAL_ABSENCE_OR_FROZEN_SPLIT\n"));
		}
	}

	private static final class LineageCategory{
		LineageCategory(final String key_, final int rank_, final boolean known_){
			assert(key_!=null && rank_>=0 && rank_<2) : "Each lineage category belongs to exactly one reported taxonomic rank";
			key=key_;
			rank=rank_;
			known=known_;
		}
		final String key;
		final int rank;
		final boolean known;
		final IntHashMap evaluationGroups=new IntHashMap();
		final IntHashMap discoveryGroups=new IntHashMap();
		final IntHashMap positive=new IntHashMap();
		int evaluationInputs;
		int discoveryInputs;
		long cores;
	}

	private static String rank(final String lineage, final String prefix){
		require(lineage!=null && prefix.length()==3, "Taxonomy ranks need an explicit GTDB rank prefix");
		int from=0;
		String result=null;
		while(from<lineage.length()){
			int to=lineage.indexOf(';', from);
			if(to<0){to=lineage.length();}
			if(lineage.startsWith(prefix, from)){
				require(result==null, "One accession lineage cannot repeat a taxonomy rank");
				result=lineage.substring(from+3, to);
			}
			from=to+1;
		}
		return result==null || result.isEmpty() ? "UNKNOWN" : result;
	}

	private void crossWrite(final String members, final String prefix, final int count){
		final ObjectMap<String, CrossCategory> byRank=new ObjectMap<String, CrossCategory>(String.class, CrossCategory.class);
		final ArrayList<CrossCategory> categories=new ArrayList<CrossCategory>();
		final long[] inputs=new long[2];
		final long[] missing=new long[2];
		for(final Node node:nodes){
			if(node.meta==null || node.meta.crossTax==null){continue;}
			final CrossTax tax=node.meta.crossTax;
			final String key=NAMES[node.cohort]+'\t'+tax.clazz+'\t'+tax.order+'\t'+tax.status();
			CrossCategory category=byRank.get(key);
			if(category==null){
				category=new CrossCategory(key, node.cohort);
				byRank.put(key, category);
				categories.add(category);
			}
			category.inputs++;
			category.cores+=node.auditNativeCores;
			category.groups.put(node.group, 1);
			if(node.auditNative){category.positive.put(node.group, 1);}
			if(node.cohort==1 || node.auditNative){category.candidates.add(node);}
			inputs[node.cohort]++;
			if(!tax.status().equals("GTDB_ACCESSION_MATCH")){missing[node.cohort]++;}
		}
		java.util.Collections.sort(categories, (a, z)->a.key.compareTo(z.key));
		for(final CrossCategory category:categories){java.util.Collections.sort(category.candidates, (a, z)->Integer.compare(a.id, z.id));}
		final ArrayList<Node> chosen=new ArrayList<Node>();
		final ObjectSet<String> selected=new ObjectSet<String>(String.class);
		for(int cohort=0; cohort<2; cohort++){
			final IntHashMap groups=new IntHashMap();
			while(groups.size()<count){
				final int before=groups.size();
				for(final CrossCategory category:categories){
					if(category.cohort!=cohort || groups.size()==count){continue;}
					while(category.next<category.candidates.size() && groups.contains(category.candidates.get(category.next).group)){category.next++;}
					if(category.next==category.candidates.size()){continue;}
					final Node node=category.candidates.get(category.next++);
					groups.put(node.group, 1);
					require(selected.add(node.key), "A crossed diagnostic cannot select one cohort/input twice");
					chosen.add(node);
				}
				require(groups.size()>before, "Insufficient distinct eligible assemblies for the requested crossed arm");
			}
		}
		try(Output out=new Output(prefix+".selected.tsv"); Output composition=new Output(prefix+".composition.tsv");
			Output eval=new Output(prefix+".evaluation.paths.tsv"); Output discovery=new Output(prefix+".discovery.paths.tsv");
			Output expected=new Output(prefix+".expected.tsv"); Output summary=new Output(prefix+".summary.tsv")){
			out.writer.println(CROSS_HEADER);
			for(final Node node:chosen){
				final Meta meta=node.meta;
				final CrossTax tax=meta.crossTax;
				out.writer.print(b.clear().append(NAMES[node.cohort]).tab().append(node.id).tab().append(node.fastaBytes)
					.tab().append(node.group).tab().append(meta.assembly.getBytes(StandardCharsets.UTF_8)).tab().append(meta.taxid)
					.append("\tbacteria\tBacteroidota\t").append(tax.clazz.getBytes(StandardCharsets.UTF_8)).tab()
					.append(tax.order.getBytes(StandardCharsets.UTF_8)).tab().append(tax.status()).tab().append(node.auditNativeCores).nl());
				(node.cohort==0 ? eval : discovery).writer.print(b.clear().append(node.fastaBytes).tab().append(node.id).nl());
			}
			composition.writer.println("cohort\tclass\torder\ttaxonomy_status\tinputs\tassembly_groups\tCM_positive_groups\tCM_admitted_cores");
			for(final CrossCategory category:categories){
				composition.writer.print(b.clear().append(category.key.getBytes(StandardCharsets.UTF_8)).tab().append(category.inputs)
					.tab().append(category.groups.size()).tab().append(category.positive.size()).tab().append(category.cores).nl());
			}
			expected.writer.println(EXTRACT_HEADER);
			int expectedCores=0;
			try(Input in=new Input(members+".extract.tsv.gz")){
				RrnaResourceIO.header(in.reader, EXTRACT_HEADER, members);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					parse(row, 21);
					if(selected.contains(key(cohort(1), p.parseInt(2))) && token(text(12), "CM_FULL_CLEAN")){
						expected.writer.print(b.clear().append(row).nl());
						expectedCores++;
					}
				}
			}
			require(expectedCores>=count, "Every selected evaluation assembly needs its actual admitted native FMN locus preserved");
			summary.writer.println("evaluation_inputs\tdiscovery_inputs\tevaluation_missing_ranks\tdiscovery_missing_ranks\tselected_per_arm\tselected_native_cores\tstatus");
			summary.writer.print(b.clear().append(inputs[0]).tab().append(inputs[1]).tab().append(missing[0]).tab().append(missing[1])
				.tab().append(count).tab().append(expectedCores).append("\tCROSSED_PROCEDURE_SELECTION_NOT_SEARCHED_OR_SPLIT_CHANGE\n"));
		}
	}

	private static final class CrossTax{
		String status(){return !seen ? "NO_ACCESSION_MATCH" : clazz.equals("UNKNOWN") || order.equals("UNKNOWN") ? "MISSING_RANK" : "GTDB_ACCESSION_MATCH";}
		String lineageStatus(final boolean withGenus){return !seen ? "NO_ACCESSION_MATCH" : family.equals("UNKNOWN") || withGenus && genus.equals("UNKNOWN") ? "MISSING_RANK" : "GTDB_ACCESSION_MATCH";}
		String clazz="UNKNOWN";
		String order="UNKNOWN";
		String family="UNKNOWN";
		String genus="UNKNOWN";
		boolean seen;
	}
	private static final class CrossCategory{
		CrossCategory(final String key_, final int cohort_){
			assert(key_!=null && cohort_>=0 && cohort_<2) : "Cross categories retain one declared cohort and rank key";
			key=key_;
			cohort=cohort_;
		}
		final String key;
		final int cohort;
		final ArrayList<Node> candidates=new ArrayList<Node>();
		final IntHashMap groups=new IntHashMap(), positive=new IntHashMap();
		int inputs, next;
		long cores;
	}

	/** Counts native GA hits and reports known-locus presence, not one-to-one grading. */
	public static void crossReport(final String selected, final String collection, final String models, final String output, final int perArm){
		require(selected!=null && collection!=null && models!=null && output!=null && perArm>0, "Cross report requires selected metadata, original models and a complete native collection");
		RrnaResourceIO.fresh(output, ".inputs.tsv", ".hits.tsv", ".expected.tsv", ".summary.tsv");
		final NcrnaConstructionMembers reader=new NcrnaConstructionMembers();
		reader.crossResults(selected, collection, models, output, perArm);
	}

	private void crossResults(final String selected, final String collection, final String models, final String prefix, final int perArm){
		final StageAudit roster=new StageAudit();
		roster.readModels(models);
		final ObjectMap<String, CrossPoint> points=new ObjectMap<String, CrossPoint>(String.class, CrossPoint.class);
		final ArrayList<CrossPoint> ordered=new ArrayList<CrossPoint>();
		final int[] sizes=new int[2];
		try(Input in=new Input(selected+".selected.tsv")){
			RrnaResourceIO.header(in.reader, CROSS_HEADER, selected);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 12);
				final int cohort=cohort(0);
				final CrossPoint point=new CrossPoint(cohort, p.parseInt(1), text(2), p.parseInt(11));
				require(points.get(point.key)==null && point.id>0 && point.fasta.startsWith("/")
					&& p.termEquals("Bacteroidota", 7) && (cohort==1 || point.expected>0), "Cross results must bind unique selected Bacteroidota inputs and positive evaluation cores");
				points.put(point.key, point);
				ordered.add(point);
				sizes[cohort]++;
			}
		}
		require(sizes[0]==perArm && sizes[1]==perArm, "The completed cross report requires the full declared selected population in both arms");
		try(Output hits=new Output(prefix+".hits.tsv"); Input in=new Input(collection)){
			hits.writer.println("cohort\tinput_id\tseqid\tstart1\tstop1\tstrand\tnative_GA\tfull_model_untruncated\tnative_table\tnative_line\traw_tblout_row");
			RrnaResourceIO.header(in.reader, "cohort\tinput_id\tprepared\tnative", collection);
			final ObjectSet<String> directories=new ObjectSet<String>(String.class);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 4);
				final CrossPoint point=points.get(key(cohort(0), p.parseInt(1)));
				final String prepared=text(2);
				final String nativePath=text(3);
				require(point!=null && !point.complete && prepared.startsWith("/") && nativePath.startsWith("/")
					&& directories.add(nativePath), "Each selected input must own one distinct complete native result");
				crossPoint(point, prepared, nativePath, hits, roster);
				point.complete=true;
			}
		}
		try(Output expected=new Output(prefix+".expected.tsv"); Input in=new Input(selected+".expected.tsv")){
			expected.writer.println("cohort\tinput_id\tmember_id\tseqid\tstrand\tstart1\tstop1\tGA_full_overlap_present\texact_GA_full_present\tinterpretation");
			RrnaResourceIO.header(in.reader, EXTRACT_HEADER, selected);
			final ObjectSet<String> ids=new ObjectSet<String>(String.class);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 21);
				final CrossPoint point=points.get(key(cohort(1), p.parseInt(2)));
				require(point!=null && ids.add(text(0)) && point.fasta.equals(text(3)) && token(text(12), "CM_FULL_CLEAN"),
					"Expected physical cores must uniquely bind selected input metadata");
				final String contig=text(8);
				final char strand=strand(9);
				final int[] bounds=intervals(10);
				require(bounds.length==2, "A clean native core must have one contiguous physical interval");
				boolean overlap=false;
				boolean exact=false;
				for(final CrossHit hit:point.hits){
					if(!hit.included || !hit.full || hit.strand!=strand || !hit.contig.equals(contig)){continue;}
					final long shared=Math.min(bounds[1], hit.stop)-Math.max(bounds[0], hit.start)+1L;
					overlap|=2*shared>bounds[1]-bounds[0]+1L && 2*shared>hit.stop-hit.start+1L;
					exact|=bounds[0]==hit.start && bounds[1]==hit.stop;
				}
				point.expectedRead++;
				if(overlap){point.recovered++;}
				if(exact){point.exact++;}
				expected.writer.print(b.clear().append(NAMES[point.cohort]).tab().append(point.id).tab().append(text(0).getBytes(StandardCharsets.UTF_8))
					.tab().append(contig.getBytes(StandardCharsets.UTF_8)).tab().append(strand).tab().append(bounds[0]).tab().append(bounds[1])
					.tab().append(overlap).tab().append(exact).append("\tPRESENCE_ONLY_NOT_ONE_TO_ONE_GRADING\n"));
			}
		}
		final long[][] totals=new long[2][7];
		try(Output out=new Output(prefix+".inputs.tsv"); Output summary=new Output(prefix+".summary.tsv")){
			out.writer.println("cohort\tinput_id\tFMN_reported\tFMN_native_GA\tFMN_GA_full\texpected_cores\texpected_overlap_present\texpected_exact_present\tGA_positive_input");
			for(final CrossPoint point:ordered){
				require(point.complete && point.expectedRead==point.expected, "No failed, missing or partially bound search may become a zero-hit result");
				out.writer.print(b.clear().append(NAMES[point.cohort]).tab().append(point.id).tab().append(point.hits.size())
					.tab().append(point.ga).tab().append(point.full).tab().append(point.expected).tab().append(point.recovered)
					.tab().append(point.exact).tab().append(point.ga>0).nl());
				final long[] total=totals[point.cohort];
				total[0]+=point.hits.size();
				total[1]+=point.ga;
				total[2]+=point.full;
				total[3]+=point.expected;
				total[4]+=point.recovered;
				total[5]+=point.exact;
				if(point.ga>0){total[6]++;}
			}
			summary.writer.println("cohort\tinputs\tFMN_reported\tFMN_native_GA\tFMN_GA_full\texpected_cores\texpected_overlap_present\texpected_exact_present\tGA_positive_inputs\tstatus");
			for(int cohort=0; cohort<2; cohort++){
				b.clear().append(NAMES[cohort]).tab().append(perArm);
				for(final long value:totals[cohort]){b.tab().append(value);}
				summary.writer.print(b.append("\tTARGETED_CROSSED_PROCEDURE_NOT_GENERALISATION_OR_SPLIT_AUTHORITY\n"));
			}
		}
	}

	/** Compare two procedures on the same discovery inputs without relabelling their cohort. */
	public static void procedureReport(final String selected, final String collection, final String models, final String output, final int expected){
		require(selected!=null && collection!=null && models!=null && output!=null && expected>0, "Declare the complete selected population and paired procedures");
		RrnaResourceIO.fresh(output, ".hits.tsv", ".inputs.tsv", ".summary.tsv");
		new NcrnaConstructionMembers().procedureResults(selected, collection, models, output, expected);
	}
	private void procedureResults(final String selected, final String collection, final String models, final String output, final int expected){
		final StageAudit roster=new StageAudit();
		roster.readModels(models);
		final ObjectMap<String, CrossPoint> points=new ObjectMap<String, CrossPoint>(String.class, CrossPoint.class);
		final ArrayList<CrossPoint> ordered=new ArrayList<CrossPoint>();
		final IntHashMap groups=new IntHashMap();
		try(Input in=new Input(selected)){
			RrnaResourceIO.header(in.reader, "cohort\tinput_id\tfasta\tassembly_group\tassembly\tgenus\tCM_admitted_cores", selected);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 7);
				if(p.termEquals("evaluation", 0)){continue;}
				require(p.termEquals("discovery", 0) && p.parseInt(1)>0 && p.parseInt(3)>=0 && !groups.contains(p.parseInt(3)),
					"Paired procedures need unique discovery assembly groups");
				groups.put(p.parseInt(3), 1);
				for(int procedure=0; procedure<2; procedure++){
					final CrossPoint point=new CrossPoint(1, p.parseInt(1), text(2), 0, procedure, PROCEDURES[procedure]);
					point.genus=text(5);
					final String key=point.label+'\t'+point.id;
					require(points.get(key)==null && point.fasta.startsWith("/"), "Duplicate or nonliteral discovery source");
					points.put(key, point);
					ordered.add(point);
				}
			}
		}
		require(groups.size()==expected && ordered.size()==2*expected, "Every selected discovery input must have two procedure slots");
		try(Input in=new Input(collection); Output hits=new Output(output+".hits.tsv")){
			RrnaResourceIO.header(in.reader, "procedure\tinput_id\tprepared\tnative", collection);
			hits.writer.println("procedure\tinput_id\tseqid\tstart1\tstop1\tstrand\tnative_GA\tfull_model_untruncated\tnative_table\tnative_line\traw_tblout_row");
			final ObjectSet<String> directories=new ObjectSet<String>(String.class);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 4);
				final CrossPoint point=points.get(text(0)+'\t'+p.parseInt(1));
				final String prepared=text(2);
				final String nativePath=text(3);
				require(point!=null && !point.complete && prepared.startsWith("/") && nativePath.startsWith("/") && directories.add(nativePath),
					"Every paired procedure requires its own bound successful native search");
				crossPoint(point, prepared, nativePath, hits, roster);
				point.complete=true;
			}
		}
		final long[][] sums=new long[2][4];
		try(Output inputs=new Output(output+".inputs.tsv"); Output summary=new Output(output+".summary.tsv")){
			inputs.writer.println("procedure\tinput_id\tgenus\tFMN_reported\tFMN_native_GA\tFMN_GA_full\tGA_positive_input");
			for(final CrossPoint point:ordered){
				require(point.complete, "A missing or failed procedure cannot be reported as zero hits");
				inputs.writer.print(b.clear().append(point.label).tab().append(point.id).tab().append(point.genus.getBytes(StandardCharsets.UTF_8))
					.tab().append(point.hits.size()).tab().append(point.ga).tab().append(point.full).tab().append(point.ga>0).nl());
				sums[point.procedure][0]+=point.hits.size();
				sums[point.procedure][1]+=point.ga;
				sums[point.procedure][2]+=point.full;
				if(point.ga>0){sums[point.procedure][3]++;}
			}
			summary.writer.println("procedure\tdiscovery_inputs\tFMN_reported\tFMN_native_GA\tFMN_GA_full\tGA_positive_inputs\tstatus");
			for(int procedure=0; procedure<2; procedure++){
				b.clear().append(PROCEDURES[procedure]).tab().append(expected);
				for(final long value:sums[procedure]){b.tab().append(value);}
				summary.writer.print(b.append("\tPAIRED_DISCOVERY_INPUT_PROCEDURES_NOT_SPLIT_AUTHORITY\n"));
			}
		}
	}
	private static final String[] PROCEDURES={"DISCOVERY_BUNDLE", "EVALUATION_FMN"};

	private void crossPoint(final CrossPoint point, final String prepared, final String nativePath, final Output out, final StageAudit roster){
		try(Input exit=new Input(nativePath+"/exit.txt")){
			final byte[] line=exit.reader.nextLine();
			require(line!=null && line.length==1 && line[0]=='0' && exit.reader.nextLine()==null, "Only a completed successful native search may enter cross results");
		}
		try(Input members=new Input(prepared+"/members.tsv")){
			RrnaResourceIO.header(members.reader, NcrnaEvidenceCatalog.CENSUS_HEADER, prepared);
			final byte[] row=members.reader.nextLine();
			require(row!=null, "Crossed native batch must retain its original input binding");
			parse(row, 5);
			require(p.parseInt(0)==point.id && text(1).equals(point.fasta) && members.reader.nextLine()==null,
				"Each cross point is exactly the selected complete genome, not another input or a pooled batch");
			point.bases=p.parseLong(3);
			point.records=p.parseLong(4);
		}
		final ObjectMap<NcrnaReferenceGeometry.Target, CrossContig> targets=
			new ObjectMap<NcrnaReferenceGeometry.Target, CrossContig>(NcrnaReferenceGeometry.Target.class, CrossContig.class);
		final NcrnaReferenceGeometry.Target probe=new NcrnaReferenceGeometry.Target();
		long bases=0;
		try(Input map=new Input(prepared+"/targets.tsv")){
			RrnaResourceIO.header(map.reader, PrepareCmSearchBatches.TARGET_HEADER, prepared);
			for(byte[] row=map.reader.nextLine(); row!=null; row=map.reader.nextLine()){
				parse(row, 6);
				probe.set(p, 1);
				require(p.parseInt(2)==point.id && p.parseInt(5)>0 && !targets.contains(probe), "Native target map must uniquely bind a complete selected input");
				targets.put(probe.copy(), new CrossContig(text(4), p.parseInt(5)));
				bases+=p.parseInt(5);
			}
		}
		require(targets.size()==point.records && bases==point.bases && bases>0, "Crossed map must conserve the complete original genome dimensions");
		long all=0;
		long number=0;
		boolean done=false;
		boolean version=false;
		final RrnaCmTruthBuilder.Fields fields=new RrnaCmTruthBuilder.Fields();
		final String table=nativePath+"/hits.tbl";
		int query=-1;
		final long[] observed=new long[roster.accessions.size()];
		try(Input in=new Input(table)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				number++;
				if(row.length==0){continue;}
				require(!done, "Native cross tblout cannot continue after its completion footer");
				if(row[0]=='#'){
					if(shared.Tools.startsWith(row, "# Version:")){version=new String(row, StandardCharsets.US_ASCII).contains("1.1.5 ");}
					if(java.util.Arrays.equals(row, "# [ok]".getBytes(StandardCharsets.US_ASCII))){done=true;}
					continue;
				}
				fields.set(row);
				if(query<0 || !fields.eq(3, roster.accessions.get(query))){query=roster.modelIndex.get(fields.text(3));}
				require(query>=0 && fields.eq(2, roster.queryNames.get(query)), "Every crossed raw model must retain its original accession/name binding");
				all++;
				observed[query]++;
				probe.set(row, fields.start[0], fields.end[0]);
				final CrossContig contig=targets.get(probe);
				require(contig!=null && fields.eq(4, "cm") && (fields.eq(16, "!") || fields.eq(16, "?")), "Crossed hits require bound targets, CM scoring and original native inclusion");
				if(point.procedure==0){require(fields.eq(16, "!"), "Discovery cut_ga procedure cannot report native-excluded hits");}
				if(!fields.eq(3, "RF00050")){
					require(point.procedure==0, "Evaluation FMN procedure searches only the original FMN model");
					continue;
				}
				require(fields.eq(2, "FMN"), "Crossed FMN accession and native query name must agree");
				final int from=fields.number(7);
				final int to=fields.number(8);
				final int start=Math.min(from, to);
				final int stop=Math.max(from, to);
				final char strand=fields.eq(9, "+") ? '+' : fields.eq(9, "-") ? '-' : '?';
				require(start>0 && stop<=contig.length && (strand=='+' && from<=to || strand=='-' && from>=to), "Crossed native coordinates must be within the captured contig and ordered for their strand");
				final boolean included=fields.eq(16, "!");
				require(fields.number(5)>0 && fields.number(6)>=fields.number(5) && fields.number(6)<=roster.modelLengths.get(query),
					"Native crossed model intervals must stay within the unchanged FMN model");
				final boolean full=fields.number(5)==1 && fields.number(6)==roster.modelLengths.get(query) && fields.eq(10, "no");
				point.hits.add(new CrossHit(contig.name, start, stop, strand, included, full));
				if(included){point.ga++;}
				if(included && full){point.full++;}
				out.writer.print(b.clear().append(point.label).tab().append(point.id).tab().append(contig.name.getBytes(StandardCharsets.UTF_8))
					.tab().append(start).tab().append(stop).tab().append(strand).tab().append(included).tab().append(full)
					.tab().append(table.getBytes(StandardCharsets.UTF_8)).tab().append(number).tab().append(row).nl());
			}
		}
		require(done && version, "Crossed native tblout must have Infernal1.1.5 and completion evidence");
		crossSearchReport(point, nativePath+"/search.out", all, observed, roster);
	}

	/** Cross-check every native model report, including zero-hit searches. */
	private static void crossSearchReport(final CrossPoint point, final String path, final long raw, final long[] observed, final StageAudit roster){
		int queries=0;
		int accessions=0;
		int dimensions=0;
		int totals=0;
		long hits=0;
		boolean fmn=false;
		boolean done=false;
		String name=null;
		int length=-1;
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		final boolean[] counted=new boolean[observed.length];
		int activeQuery=-1;
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				if(shared.Tools.startsWith(row, "Query:")){
					final String[] fields=new String(row, StandardCharsets.US_ASCII).trim().split("\\s+");
					require(fields.length==3 && fields[2].startsWith("[CLEN=") && fields[2].endsWith("]"), "Native cross query must identify its name and CM length");
					name=fields[1];
					length=Integer.parseInt(fields[2].substring(6, fields[2].length()-1));
					queries++;
				}
				else if(shared.Tools.startsWith(row, "Accession:")){
					final String accession=new String(row, StandardCharsets.US_ASCII).substring(10).trim();
					final int query=roster.modelIndex.get(accession);
					require(query>=0 && seen.add(accession) && roster.queryNames.get(query).equals(name)
						&& roster.modelLengths.get(query)==length && (point.procedure==0 || accession.equals("RF00050")),
						"Crossed query name/accession/length must match each original model exactly once");
					accessions++;
					activeQuery=query;
					fmn|=accession.equals("RF00050");
				}else if(shared.Tools.startsWith(row, "Target sequences:")){
					final String[] tokens=new String(row, StandardCharsets.US_ASCII).trim().split("\\s+");
					require(Long.parseLong(tokens[2])==point.records && Long.parseLong(tokens[3].substring(1))==2*point.bases,
						"Every native model must have searched the complete selected genome on both strands");
					dimensions++;
				}else if(shared.Tools.startsWith(row, "Total CM hits reported:")){
					final String value=new String(row, StandardCharsets.US_ASCII).substring(23).trim();
					final long reported=Long.parseLong(value.split("\\s+")[0]);
					require(activeQuery>=0 && !counted[activeQuery] && reported==observed[activeQuery],
						"Each native model report must independently conserve its raw rows, including zero-hit models");
					counted[activeQuery]=true;
					hits+=reported;
					totals++;
				}else if(java.util.Arrays.equals(row, "[ok]".getBytes(StandardCharsets.US_ASCII))){done=true;}
			}
		}
		final int expected=point.procedure==0 ? roster.accessions.size() : 1;
		require(queries==expected && accessions==expected && dimensions==expected && totals==expected && fmn && done && hits==raw,
			"Crossed native report must conserve every expected model, complete target dimensions and raw hit count");
	}

	private static final class CrossPoint{
		CrossPoint(final int cohort_, final int id_, final String fasta_, final int expected_){
			this(cohort_, id_, fasta_, expected_, cohort_, NAMES[cohort_]);
		}
		CrossPoint(final int cohort_, final int id_, final String fasta_, final int expected_, final int procedure_, final String label_){
			assert(cohort_>=0 && cohort_<2 && expected_>=0) : "Cross points retain a declared cohort and nonnegative saved core count";
			assert(procedure_>=0 && procedure_<2 && label_!=null) : "Search procedure is independent of the original input cohort";
			cohort=cohort_;
			procedure=procedure_;
			label=label_;
			id=id_;
			fasta=fasta_;
			expected=expected_;
			key=key(cohort, id);
		}
		final int cohort, id, expected, procedure;
		final String key, fasta, label;
		String genus="NA";
		final ArrayList<CrossHit> hits=new ArrayList<CrossHit>();
		long bases, records;
		int ga, full, expectedRead, recovered, exact;
		boolean complete;
	}
	private static final class CrossContig{
		CrossContig(final String name_, final int length_){
			assert(name_!=null && length_>0) : "A captured contig retains its name and positive length";
			name=name_;
			length=length_;
		}
		final String name;
		final int length;
	}
	private static final class CrossHit{
		CrossHit(final String contig_, final int start_, final int stop_, final char strand_, final boolean included_, final boolean full_){
			assert(start_>0 && stop_>=start_) : "Retained native hits use normalized inclusive coordinates";
			contig=contig_;
			start=start_;
			stop=stop_;
			strand=strand_;
			included=included_;
			full=full_;
		}
		final String contig;
		final int start, stop;
		final char strand;
		final boolean included, full;
	}
	public static final String CROSS_HEADER="cohort\tinput_id\tfasta\tassembly_group\tassembly\ttaxid\tdomain\tphylum\tclass\torder\ttaxonomy_status\tCM_admitted_cores";

	/** Compare native discovery rows with saved import, classification and admission. */
	public static void stageAudit(final String members, final String groups, final String genomes, final String evidence,
			final String imported, final String collection, final String models, final String output,
			final int genomeCount, final int oldCount, final int evalCount, final int newCount,
			final int coreCount, final int batchCount, final long rawCount){
		require(members!=null && groups!=null && genomes!=null && evidence!=null && imported!=null && collection!=null
			&& models!=null && output!=null && genomeCount>0 && oldCount>=0 && evalCount>0 && newCount>0
			&& coreCount>=0 && batchCount>0 && rawCount>=0, "Stage audit needs all saved populations and raw batch bindings");
		RrnaResourceIO.fresh(output, ".stages.tsv", ".examples.tsv.gz", ".summary.tsv");
		final NcrnaConstructionMembers ledger=new NcrnaConstructionMembers();
		ledger.genomes(genomes, genomeCount);
		ledger.groups(groups, evalCount, newCount);
		ledger.auditBindings();
		ledger.auditCores(members+".extract.tsv.gz", coreCount);
		ledger.auditOld(members+".old.tsv.gz", oldCount);
		for(final Node node:ledger.nodes){node.stages=new long[10];}
		final StageAudit audit=ledger.new StageAudit();
		audit.readModels(models);
		audit.readImport(imported+".hits.tsv.gz");
		ledger.expected(evidence+".summary.tsv");
		audit.readEvents(evidence+".fmn_events.tsv.gz");
		try(Output examples=new Output(output+".examples.tsv.gz")){
			examples.writer.println("domain\tphylum\tinput_id\tdropping_stage\tnative_table\tnative_line\traw_tblout_row");
			audit.readRaw(collection, batchCount, newCount, rawCount, examples);
		}
		audit.write(output, rawCount, batchCount);
	}

	private final class StageAudit{
		StageAudit(){
			for(int i=0; i<nodes.size(); i++){
				final Node node=nodes.get(i);
				if(node.cohort==1){discoveryIndex.put(node.id, i);}
			}
		}
		Node discovery(final int id){
			final int index=discoveryIndex.get(id);
			require(index>=0, "Stage evidence refers to an absent discovery input");
			return nodes.get(index);
		}
		void readModels(final String path){
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, NcrnaReferenceEvidence.QUERY_HEADER, path);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					parse(row, 3);
					final String accession=text(0);
					require(!modelIndex.contains(accession) && p.parseInt(2)>0, "Raw model roster needs unique accessions and positive CLEN");
					modelIndex.put(accession, modelIndex.size());
					accessions.add(accession);
					queryNames.add(text(1));
					modelLengths.add(p.parseInt(2));
				}
			}
			require(modelIndex.contains("RF00050"), "Zero FMN rows are interpretable only if RF00050 was searched");
		}

		void readImport(final String path){
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, BuildNcrnaReference.BATCH_HIT_HEADER, path);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					parse(row, 24);
					final Node node=discovery(p.parseInt(0));
					node.stages[4]++;
					importRows++;
					if(!p.termEquals("RF00050", 7)){continue;}
					require(flag(17)==flag(22) && p.termEquals(BuildNcrnaReference.BOUND_GA_BASIS, 23),
						"Imported FMN GA classification must retain the original native inclusion decision");
					require(p.parseInt(13)==modelLengths.get(modelIndex.get("RF00050")), "Imported FMN CLEN must match the searched model roster");
					node.stages[5]++;
					final String key=witnessKey(text(1), p.parseLong(2));
					require(witnesses.get(key)==null, "One physical native FMN row cannot be imported twice");
					final StageWitness witness=new StageWitness(node, text(3), text(21), p.parseInt(4), p.parseInt(5), strand(6),
						p.parseInt(11), p.parseInt(12), p.parseInt(14), flag(22));
					witnesses.put(key, witness);
					witnessList.add(witness);
				}
			}
		}

		void readEvents(final String path){
			final ObjectSet<String> eventIds=new ObjectSet<String>(String.class);
			final long[] counts=new long[2];
			final long[] fulls=new long[2];
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, FMN_HEADER, path);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					parse(row, 32);
					final int cohort=cohort(0);
					final Node node=node(cohort, p.parseInt(3));
					require(p.termEquals("RF00050", 8) && p.termEquals("FMN", 9)
						&& p.parseInt(2)>0 && eventIds.add(key(cohort, p.parseInt(2))), "Saved FMN events need unique cohort-qualified IDs and consistent family");
					final boolean included=flag(19);
					require(p.parseInt(12)==modelLengths.get(modelIndex.get("RF00050")), "Both cohort classifications must refer to the same FMN model length");
					final boolean full=full(10, 11, 12, 13, 19);
					require(full==flag(1), "Saved native full-model flag must match its original measurements");
					final boolean clean=full && p.termEquals("CONFIDENT", 24);
					require(!p.termEquals("CONFIDENT", 24) || full, "A confident saved FMN event requires full native inclusion support");
					counts[cohort]++;
					if(full){fulls[cohort]++;}
					node.stages[6]++;
					if(included){node.stages[7]++;}
					if(full){node.stages[8]++;}
					if(clean){node.stages[9]++;}
					if(cohort==1){
						final StageWitness witness=witnesses.get(witnessKey(text(28), p.parseLong(29)));
						require(witness!=null && witness.node==node && !witness.classified,
							"A discovery classification must bind exactly one imported native FMN row and input");
						require(witness.from==p.parseInt(10) && witness.to==p.parseInt(11) && witness.truncation==p.parseInt(13)
							&& witness.included==included, "Classification must preserve imported native model span, truncation and inclusion");
						require(witness.target.equals(text(4)) && witness.start==p.parseInt(5) && witness.stop==p.parseInt(6)
							&& witness.strand==strand(7), "Classification must preserve the imported contig, coordinates and strand of the same physical locus");
						witness.classified=true;
						witness.reason=clean ? null : !included ? "NATIVE_NOT_INCLUDED" : !full
							? "NOT_FULL_MODEL_OR_NATIVE_TRUNCATION" : "TIER_"+text(24)+":"+text(25);
					}
				}
			}
			for(int c=0; c<2; c++){require(counts[c]==expectedEvents[c] && fulls[c]==expectedFull[c], "The stage audit must read the complete accepted FMN event populations");}
		}

		void readRaw(final String path, final int expectedBatches, final int expectedInputs, final long expectedRows, final Output examples){
			final ObjectSet<String> batchPaths=new ObjectSet<String>(String.class);
			final IntHashMap seenInputs=new IntHashMap();
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, "scan_group\tbatch_id\tprepared_metadata\tnative_result\tgenomes\tbases\trecords\tstate", path);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					parse(row, 8);
					final String batch=text(1);
					final String prepared=text(2);
					final String nativeRoot=text(3);
					final int inputCount=p.parseInt(4);
					final long bases=p.parseLong(5);
					final long records=p.parseLong(6);
					require(batchPaths.add(nativeRoot) && inputCount>0 && bases>0 && records>0
						&& p.termEquals("RAW_NATIVE_GA_COMPLETE", 7), "Read each complete declared raw batch exactly once");
					final ObjectMap<NcrnaReferenceGeometry.Target, StageTarget> targets=new ObjectMap<NcrnaReferenceGeometry.Target, StageTarget>(NcrnaReferenceGeometry.Target.class, StageTarget.class);
					final IntHashMap batchInputs=new IntHashMap();
					long mapped=0;
					try(Input map=new Input(prepared+"/targets.tsv")){
						RrnaResourceIO.header(map.reader, PrepareCmSearchBatches.TARGET_HEADER, prepared);
						for(byte[] target=map.reader.nextLine(); target!=null; target=map.reader.nextLine()){
							parse(target, 6);
							final Node node=discovery(p.parseInt(2));
							require(p.termEquals(batch, 0) && p.parseLong(3)>0 && p.parseInt(5)>0, "Captured target map must bind this batch and a positive contig");
							probe.set(p, 1);
							require(!targets.contains(probe), "Synthetic native targets must be unique inside their batch");
							targets.put(probe.copy(), new StageTarget(node, text(4), p.parseInt(5)));
							if(!batchInputs.contains(node.id)){
								require(!seenInputs.contains(node.id), "A discovery input must not be counted in two raw batches");
								batchInputs.put(node.id, 1);
								seenInputs.put(node.id, 1);
							}
							mapped++;
						}
					}
					require(mapped==records && batchInputs.size()==inputCount, "Raw batch map must conserve every declared contig and input");
					final long[] expected=completion(nativeRoot, batch, records, bases);
					final long[] observed=new long[expected.length];
					final String table=nativeRoot+"/hits.tbl";
					long line=0;
					int model=-1;
					boolean program=false, version=false, cutGa=false, finished=false;
					try(Input hits=new Input(table)){
						for(byte[] hit=hits.reader.nextLine(); hit!=null; hit=hits.reader.nextLine()){
							line++;
							if(hit.length==0){continue;}
							require(!finished, "Raw tblout content follows its completion footer");
							if(hit[0]=='#'){
								final String footer=new String(hit, StandardCharsets.UTF_8);
								if(footer.startsWith("# Program:")){program=footer.substring(10).trim().equals("cmsearch");}
								if(footer.startsWith("# Version:")){version=footer.substring(10).trim().startsWith("1.1.5 ");}
								if(footer.startsWith("# Option settings:")){cutGa=footer.contains(" --cut_ga ");}
								if(footer.equals("# [ok]")){finished=true;}
								continue;
							}
							fields.set(hit);
							if(model<0 || !fields.eq(3, accessions.get(model))){model=modelIndex.get(fields.text(3));}
							require(model>=0 && fields.eq(2, queryNames.get(model)) && fields.eq(4, "cm") && fields.eq(16, "!"),
								"Native GA bundle rows must use the declared CM accession/name and native inclusion");
							probe.set(hit, fields.start[0], fields.end[0]);
							final StageTarget target=targets.get(probe);
							require(target!=null, "Raw native target lacks its captured input binding");
							final Node node=target.node;
							observed[model]++;
							rawRows++;
							node.stages[0]++;
							if(fields.eq(3, "RF00010")){node.stages[2]++;}
							if(fields.eq(3, "RF00174")){node.stages[3]++;}
							if(!fields.eq(3, "RF00050")){continue;}
							node.stages[1]++;
							final int seqFrom=fields.number(7);
							final int seqTo=fields.number(8);
							final int start=Math.min(seqFrom, seqTo);
							final int stop=Math.max(seqFrom, seqTo);
							final char strand=fields.eq(9, "+") ? '+' : fields.eq(9, "-") ? '-' : '?';
							require(start>0 && stop<=target.length && (strand=='+' && seqFrom<=seqTo || strand=='-' && seqFrom>=seqTo),
								"Raw FMN coordinates must be ordered for their strand and contained in the captured contig");
							final StageWitness witness=witnesses.get(witnessKey(table, line));
							if(witness!=null){
								require(witness.node==node && !witness.rawSeen, "Imported FMN witness must bind the same raw input and unique physical line");
								require(fields.eq(0, witness.nativeTarget) && witness.target.equals(target.contig)
									&& witness.start==start && witness.stop==stop && witness.strand==strand,
									"Raw/import native target, mapped contig, coordinates and strand must identify the same physical locus");
								// Same native flag convention as NcrnaReferenceEvidence.Hit.
								final int truncation=fields.eq(10, "no") ? 0 : fields.eq(10, "5'") ? 1 : fields.eq(10, "3'") ? 2 : fields.eq(10, "5'&3'") ? 3 : -1;
								require(witness.from==fields.number(5) && witness.to==fields.number(6) && witness.truncation==truncation
									&& witness.included, "Raw/import model spans and native flags must agree before assigning a dropping stage");
								witness.rawSeen=true;
							}
							final String reason=witness==null ? "ABSENT_FROM_IMPORT" : !witness.classified ? "ABSENT_FROM_TIER_EVENTS" : witness.reason;
							if(reason!=null){example(examples, node, table, line, reason, hit);}
						}
					}
					require(program && version && cutGa && finished, "Raw tblout needs complete cmsearch 1.1.5 GA execution evidence");
					for(int i=0; i<expected.length; i++){require(expected[i]==observed[i], "Every raw model count, including RF00050 zeroes, must agree with native completion");}
				}
			}
			require(batchPaths.size()==expectedBatches && seenInputs.size()==expectedInputs && rawRows==expectedRows,
				"Raw audit must cover every declared completed batch, discovery input and native row");
			for(final StageWitness witness:witnessList){require(witness.rawSeen, "An imported FMN witness is absent from the complete native scan");}
		}

		long[] completion(final String directory, final String batch, final long records, final long bases){
			try(Input exit=new Input(directory+"/exit.txt")){
				final byte[] row=exit.reader.nextLine();
				require(row!=null && row.length==1 && row[0]=='0' && exit.reader.nextLine()==null, "Raw batch search must have exited successfully");
			}
			final long[] result=new long[modelIndex.size()];
			java.util.Arrays.fill(result, -1);
			try(Input in=new Input(directory+"/completion.tsv")){
				RrnaResourceIO.header(in.reader, NcrnaReferenceBatch.BUNDLED_COMPLETION_HEADER, directory);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					parse(row, 7);
					final int model=modelIndex.get(text(1));
					require(model>=0 && result[model]<0 && text(0).equals(batch) && p.parseLong(2)==records
						&& p.parseLong(3)==bases && p.parseLong(4)==Math.multiplyExact(2L, bases)
						&& p.parseLong(5)>=0 && p.termEquals("COMPLETE", 6), "Native completion must cover this exact batch and each searched model once");
					result[model]=p.parseLong(5);
				}
			}
			for(final long count:result){require(count>=0, "Missing model completion cannot be interpreted as a zero-hit search");}
			return result;
		}

		void example(final Output out, final Node node, final String table, final long line, final String reason, final byte[] raw){
			require(node.meta!=null, "Discovery examples require the same original taxonomy binding as their denominator");
			final String key=node.meta.domain+'\t'+node.meta.phylum+'\t'+reason;
			final int previous=exampleCounts.get(key);
			final int count=previous<0 ? 0 : previous;
			if(count>=5){return;}
			exampleCounts.put(key, count+1);
			out.writer.print(b.clear().append(node.meta.domain.getBytes(StandardCharsets.UTF_8)).tab()
				.append(node.meta.phylum.getBytes(StandardCharsets.UTF_8)).tab().append(node.id).tab().append(reason.getBytes(StandardCharsets.UTF_8))
				.tab().append(table.getBytes(StandardCharsets.UTF_8)).tab().append(line).tab().append(raw).nl());
		}

		void write(final String prefix, final long expectedRaw, final int batches){
			final ObjectMap<String, StageCoverage> labels=new ObjectMap<String, StageCoverage>(String.class, StageCoverage.class);
			final ArrayList<StageCoverage> ordered=new ArrayList<StageCoverage>();
			for(final Node node:nodes){
				require(node.auditNativeCores<=node.stages[9], "Unique CM-admitted cores cannot exceed their complete clean native event witnesses");
				final String key=NAMES[node.cohort]+'\t'+(node.meta==null ? "UNKNOWN\tUNKNOWN" : node.meta.domain+'\t'+node.meta.phylum);
				StageCoverage row=labels.get(key);
				if(row==null){
					row=new StageCoverage(key, node.cohort);
					labels.put(key, row);
					ordered.add(row);
				}
				row.inputs++;
				if(node.split.equals(TRAIN)){row.training++;}
				for(int i=0; i<node.stages.length; i++){row.counts[i]+=node.stages[i];}
				row.nativeCores+=node.auditNativeCores;
				row.priorCores+=node.auditPriorCores;
				if(node.split.equals(TRAIN)){
					row.trainingNativeCores+=node.auditNativeCores;
					row.trainingPriorCores+=node.auditPriorCores;
				}
				if(node.auditIntact){row.intactGroups.put(node.group, 1);}
				if(node.stages[7]>0){row.gaGroups.put(node.group, 1);}
				if(node.auditNative){row.nativeGroups.put(node.group, 1);}
				if(node.meta!=null && node.meta.auditNode==node){row.oldMembers+=node.meta.auditOld;}
			}
			try(Output out=new Output(prefix+".stages.tsv"); Output summary=new Output(prefix+".summary.tsv")){
				out.writer.println("cohort\tdomain\tphylum\tinputs\ttraining_inputs\traw_all_events\traw_FMN_events\traw_RF00010_events\traw_RF00174_events\timport_all_events\timport_FMN_events\ttier_FMN_events\tnative_GA_FMN_events\tnative_GA_full_untruncated_events\tfull_CONFIDENT_events\tCM_admitted_cores\told_only_intact_cores\told_members\tintact_candidate_groups\tnative_GA_positive_groups\tCM_admitted_groups\ttraining_CM_admitted_cores\ttraining_old_only_intact_cores\tintact_groups_with_native_GA\tintact_groups_without_native_GA");
				for(final StageCoverage row:ordered){
					int supported=0;
					for(final int group:row.intactGroups.keys()){
						if(group>=0 && row.gaGroups.contains(group)){supported++;}
					}
					b.clear().append(row.key.getBytes(StandardCharsets.UTF_8)).tab().append(row.inputs).tab().append(row.training);
					for(int i=0; i<row.counts.length; i++){
						b.tab();
						if(row.cohort==0 && i<6){b.append("NA");}else{b.append(row.counts[i]);}
					}
					out.writer.print(b.tab().append(row.nativeCores).tab().append(row.priorCores).tab().append(row.oldMembers)
						.tab().append(row.intactGroups.size()).tab().append(row.gaGroups.size()).tab().append(row.nativeGroups.size())
						.tab().append(row.trainingNativeCores).tab().append(row.trainingPriorCores)
						.tab().append(supported).tab().append(row.intactGroups.size()-supported).nl());
				}
				require(rawRows==expectedRaw, "Reported raw totals must retain the declared native population");
				summary.writer.println("raw_batches\traw_rows\timport_rows\timport_FMN_witnesses\tstatus");
				summary.writer.print(b.clear().append(batches).tab().append(rawRows).tab().append(importRows).tab().append(witnessList.size())
					.append("\tSAVED_STAGE_COUNTS_NOT_BIOLOGICAL_ABSENCE_OR_SPLIT_AUTHORITY\n"));
			}
		}
		private String witnessKey(final String table, final long line){
			require(table.startsWith("/") && line>0, "A stage witness needs its literal absolute native path and positive physical line");
			return table+'\t'+line;
		}
		final ObjectIntMap<String> modelIndex=new ObjectIntMap<String>(String.class);
		final ArrayList<String> accessions=new ArrayList<String>();
		final ArrayList<String> queryNames=new ArrayList<String>();
		final IntList modelLengths=new IntList();
		final IntHashMap discoveryIndex=new IntHashMap();
		final ObjectIntMap<String> exampleCounts=new ObjectIntMap<String>(String.class);
		final ObjectMap<String, StageWitness> witnesses=new ObjectMap<String, StageWitness>(String.class, StageWitness.class);
		final ArrayList<StageWitness> witnessList=new ArrayList<StageWitness>();
		final NcrnaReferenceGeometry.Target probe=new NcrnaReferenceGeometry.Target();
		final RrnaCmTruthBuilder.Fields fields=new RrnaCmTruthBuilder.Fields();
		long rawRows, importRows;
	}

	private static final class StageWitness{
		StageWitness(final Node node_, final String target_, final String nativeTarget_, final int start_, final int stop_, final char strand_,
				final int from_, final int to_, final int truncation_, final boolean included_){
			require(!target_.isEmpty() && !nativeTarget_.isEmpty() && start_>0 && stop_>=start_,
				"An imported witness must retain complete contig and normalized locus identity");
			node=node_;
			target=target_;
			nativeTarget=nativeTarget_;
			start=start_;
			stop=stop_;
			strand=strand_;
			from=from_;
			to=to_;
			truncation=truncation_;
			included=included_;
		}
		final Node node;
		final String target, nativeTarget;
		final int start, stop;
		final char strand;
		final int from, to, truncation;
		final boolean included;
		boolean classified, rawSeen;
		String reason;
	}
	private static final class StageTarget{
		StageTarget(final Node node_, final String contig_, final int length_){
			require(!contig_.isEmpty() && length_>0, "Captured native targets require their literal contig and length");
			node=node_;
			contig=contig_;
			length=length_;
		}
		final Node node;
		final String contig;
		final int length;
	}
	private static final class StageCoverage{
		StageCoverage(final String key_, final int cohort_){key=key_; cohort=cohort_;}
		final String key;
		final int cohort;
		final long[] counts=new long[10];
		final IntHashMap intactGroups=new IntHashMap(), gaGroups=new IntHashMap(), nativeGroups=new IntHashMap();
		long inputs, training, nativeCores, priorCores, oldMembers, trainingNativeCores, trainingPriorCores;
	}

	private void writeAudit(final String prefix, final int oldCount, final int coreCount){
		final ObjectMap<String, AuditCoverage> byLabel=new ObjectMap<String, AuditCoverage>(String.class, AuditCoverage.class);
		final ArrayList<AuditCoverage> labels=new ArrayList<AuditCoverage>();
		for(final Meta meta:metadata){
			final AuditCoverage coverage=auditLabel(meta.domain+'\t'+meta.phylum, byLabel, labels);
			coverage.counts[7]+=meta.auditOld;
			if(meta.auditDiscovery){coverage.counts[8]+=meta.auditOld;}
			if(meta.auditNode.split.equals(TRAIN)){coverage.counts[9]+=meta.auditOld;}
			if(meta.auditOld>0){
				coverage.counts[10]++;
				if(meta.auditNode.split.equals(TRAIN)){coverage.counts[11]++;}
			}
		}
		int missing=0;
		int training=0;
		int cores=0;
		try(Output unknown=new Output(prefix+".unlabelled.tsv.gz"); Output out=new Output(prefix+".coverage.tsv");
				Output summary=new Output(prefix+".summary.tsv")){
			unknown.writer.println("cohort\tinput_id\tfasta\tassembly_group\tassembly_split\tcandidate_cores\tintact_metadata_present\tnative_clean_core_present\treason");
			for(final Node node:nodes){
				final Meta meta=node.meta;
				final AuditCoverage coverage=auditLabel(meta==null ? "UNKNOWN\tUNKNOWN" : meta.domain+'\t'+meta.phylum, byLabel, labels);
				final boolean train=node.split.equals(TRAIN);
				coverage.counts[node.cohort]++;
				if(train){
					training++;
					coverage.counts[2]++;
				}
				if(node.cohort==1 && node.auditIntact){coverage.counts[3]++;}
				if(train && node.auditIntact){coverage.counts[4]++;}
				if(node.cohort==1 && node.auditNative){coverage.counts[5]++;}
				if(train && node.auditNative){coverage.counts[6]++;}
				if(node.auditIntact){coverage.intact.put(node.group, 1);}
				if(train && node.auditIntact){coverage.trainIntact.put(node.group, 1);}
				cores+=node.auditCores;
				if(meta==null){
					missing++;
					unknown.writer.print(b.clear().append(NAMES[node.cohort]).tab().append(node.id).tab().append(node.fastaBytes)
						.tab().append(node.group).tab().append(node.split).tab().append(node.auditCores).tab().append(node.auditIntact)
						.tab().append(node.auditNative).append("\tNO_LITERAL_FASTA_BINDING_IN_ORIGINAL_ROSTER_NOT_PROOF_OF_UNSEEN_ASSEMBLY\n"));
				}
			}
			out.writer.println(AUDIT_HEADER);
			long evalSum=0;
			long discoverySum=0;
			long oldSum=0;
			for(final AuditCoverage coverage:labels){
				b.clear().append(coverage.key.getBytes(StandardCharsets.UTF_8));
				for(final long value:coverage.counts){b.tab().append(value);}
				out.writer.print(b.tab().append(coverage.intact.size()).tab().append(coverage.trainIntact.size()).nl());
				evalSum+=coverage.counts[0];
				discoverySum+=coverage.counts[1];
				oldSum+=coverage.counts[7];
			}
			require(evalSum+discoverySum==nodes.size() && oldSum==oldCount && cores==coreCount,
				"Per-label denominators must conserve all inputs, original members and saved candidate cores");
			summary.writer.println("original_genomes\tevaluation_inputs\tdiscovery_inputs\ttraining_inputs\told_members\tcandidate_cores\tunlabelled_inputs\tstatus");
			summary.writer.print(b.clear().append(metadata.size()).tab().append(evalSum).tab().append(discoverySum).tab().append(training)
				.tab().append(oldCount).tab().append(cores).tab().append(missing).append("\tORIGINAL_LABELS_ONLY_NO_TAXONOMY_OR_SPLIT_CHANGE\n"));
		}
	}

	private static AuditCoverage auditLabel(final String key, final ObjectMap<String, AuditCoverage> map, final ArrayList<AuditCoverage> list){
		require(key!=null && !key.isEmpty(), "Coverage labels must retain domain and phylum, including explicit unknowns");
		AuditCoverage result=map.get(key);
		if(result==null){
			result=new AuditCoverage(key);
			map.put(key, result);
			list.add(result);
		}
		return result;
	}
	private static boolean token(final String values, final String token){
		require(values!=null && token!=null && !token.isEmpty(), "Admission token matching requires explicit values");
		int from=0;
		while(from<values.length()){
			int end=values.indexOf(';', from);
			if(end<0){end=values.length();}
			if(end-from==token.length() && values.regionMatches(from, token, 0, token.length())){return true;}
			from=end+1;
		}
		return false;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Original Identity      ----------------*/
	/*--------------------------------------------------------------*/

	private void genomes(final String path, final int expected){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, GENOME_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 11);
				require(p.parseInt(0)==metadata.size(), "Original genome ordinals must remain consecutive");
				final Meta meta=new Meta(text(1), text(7), text(6), text(3), text(4), text(5), text(2));
				require(meta.fasta.startsWith("/") && meta.gff.startsWith("/") && !meta.genome.isEmpty()
					&& !meta.taxid.isEmpty() && !meta.domain.isEmpty() && !meta.phylum.isEmpty(), "Original genome provenance is incomplete");
				require(byGenome.get(meta.genome)==null && byFasta.get(meta.fasta)==null,
					"Original genome names and literal FASTA paths must remain unique");
				byGenome.put(meta.genome, meta);
				byFasta.put(meta.fasta, meta);
				metadata.add(meta);
			}
		}
		require(metadata.size()==expected, "Original genome roster was truncated or expanded");
	}

	private void groups(final String path, final int expectedEval, final int expectedNew){
		final IntHashMap trainGroups=new IntHashMap();
		final IntHashMap evalGroups=new IntHashMap();
		final int[] counts=new int[2];
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, GROUP_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 14);
				final int cohort=cohort(0);
				final int id=p.parseInt(1);
				final String fasta=text(2);
				final int group=p.parseInt(10);
				final boolean evaluation=flag(11);
				final String disposition=text(13);
				require(id>0 && fasta.startsWith("/") && group>=0, "Integrated input identity is invalid");
				final boolean train=disposition.equals("TRAIN_CANDIDATE_ASSEMBLY_IDENTITY");
				final boolean uncertain=disposition.equals("ASSEMBLY_SEQUENCE_IDENTITY_UNRESOLVED")
					|| disposition.equals("ASSEMBLY_IDENTITY_UNRESOLVED");
				require(cohort==0 ? evaluation && disposition.equals("EVALUATION") : evaluation
					? disposition.equals("EXCLUDE_KNOWN_EVALUATION_ASSEMBLY") || disposition.equals("EXCLUDE_VERIFIED_EVALUATION_ASSEMBLY_ALIAS")
					: train || uncertain, "Resolved disposition conflicts with cohort/evaluation-group evidence");
				final Node node=new Node(cohort, id, fasta, group, evaluation ? HELD : train ? TRAIN : UNRESOLVED, byFasta.get(fasta));
				require(nodesById.get(node.key)==null, "Repeated cohort-qualified assembly input: "+node.key);
				if(cohort==1){require(node.meta!=null, "Every discovery input must bind the complete original PGAP taxonomy roster: "+fasta);}
				if(train){
					require(!trainGroups.contains(group), "Multiple training paths share one assembly group; resolve representative locus provenance before weighting them");
					trainGroups.put(group, id);
				}
				if(cohort==0){evalGroups.put(group, id);}
				nodesById.put(node.key, node);
				nodes.add(node);
				counts[cohort]++;
			}
		}
		require(counts[0]==expectedEval && counts[1]==expectedNew, "Complete integrated cohort sizes differ from the declared counts");
		for(final Node node:nodes){require(node.split.equals(HELD)==evalGroups.contains(node.group),
			"An evaluation-labelled group needs an actual evaluation input, and may never supply training aliases");}
	}

	private void members(final String path, final int expected){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionEvidence.MEMBER_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 32);
				require(flag(0) && p.termEquals("UNVERIFIED", 3) && p.termEquals("NOT_EVALUATED", 4), "Consume original, unmodified construction bindings");
				final Meta meta=byGenome.get(text(6));
				require(meta!=null && text(8).equals(meta.taxid) && text(9).equals(meta.domain) && text(10).equals(meta.phylum),
					"Original locus taxonomy must agree with its original genome roster");
				final Old old=new Old(row, text(5), meta, text(11), strand(12), intervals(13), p.parseInt(14),
					flag(17), flag(18), p.parseInt(22), p.parseInt(20), p.termEquals("ANNOTATION_FULL_COPY", 19));
				require(old.id.startsWith("RF00050_") && oldById.get(old.id)==null && old.undefined>=0,
					"Original FMN members must be unique, with explicit ambiguity counts");
				oldById.put(old.id, old);
				olds.add(old);
				meta.members.add(old);
			}
		}
		require(olds.size()==expected, "Every original consensus/seed input must survive admission accounting");
	}

	private void bindings(final String path, final int expectedGenomes){
		final ObjectSet<String> observed=new ObjectSet<String>(String.class);
		final ObjectSet<String> seenBindings=new ObjectSet<String>(String.class);
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionEvidence.BINDING_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 6);
				final Meta meta=byGenome.get(text(0));
				final Node node=node(cohort(3), p.parseInt(4));
				require(meta!=null && text(1).equals(meta.fasta) && text(2).equals(meta.gff)
					&& node.fasta.equals(meta.fasta) && p.termEquals("LITERAL_CATALOG_PATH_NOT_SEQUENCE_PROOF", 5),
					"Original binding must preserve its actual genome, FASTA, GFF and cohort-qualified input");
				require(seenBindings.add(meta.genome+'\t'+node.key), "Repeated original genome/input binding would duplicate provenance");
				observed.add(meta.genome);
				for(final Old old:meta.members){
					final String key=old.id+'\t'+node.key;
					require(stateById.get(key)==null, "Repeated original member/input binding");
					final State state=new State(old, node);
					stateById.put(key, state);
					allStates.add(state);
					old.states.add(state);
					if(old.node==null || node.cohort<old.node.cohort || node.cohort==old.node.cohort && node.id<old.node.id){old.node=node;}
				}
			}
		}
		require(observed.size()==expectedGenomes, "Every original genome, including zero-locus genomes, needs a binding");
		for(final Old old:olds){
			require(old.node!=null, "Original member lacks a catalog binding: "+old.id);
			for(final State state:old.states){require(state.node.group==old.node.group && state.node.split.equals(old.node.split),
				"Literal aliases of one original member cannot disagree on assembly grouping or exclusion");}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Native Evidence      ----------------*/
	/*--------------------------------------------------------------*/

	private void expected(final String path){
		final boolean[] seen=new boolean[2];
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, SUMMARY_INPUT_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 8);
				final int cohort=cohort(0);
				require(!seen[cohort] && p.termEquals("EVIDENCE_JOIN_NOT_ADMISSION_OR_ALIAS_PROOF", 7), "Need one complete diagnostic summary per cohort");
				seen[cohort]=true;
				expectedStates[cohort]=p.parseLong(1);
				expectedLinks[cohort]=p.parseLong(3);
				expectedEvents[cohort]=p.parseLong(5);
				expectedFull[cohort]=p.parseLong(6);
				require(expectedStates[cohort]>=0 && expectedLinks[cohort]>=0 && expectedEvents[cohort]>=0
					&& expectedFull[cohort]>=0 && expectedFull[cohort]<=expectedEvents[cohort], "Saved cohort counts cannot be negative or exceed their denominators");
			}
		}
		require(seen[0] && seen[1], "Both cohort summaries are required, including zero-event cohorts");
	}

	private void events(final String path){
		final long[] counts=new long[2];
		final long[] fulls=new long[2];
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, FMN_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 32);
				final int cohort=cohort(0);
				final int id=p.parseInt(2);
				final Node node=node(cohort, p.parseInt(3));
				final boolean full=full(10, 11, 12, 13, 19);
				require(id>0 && flag(1)==full && p.termEquals("RF00050", 8) && p.termEquals("FMN", 9),
					"Retained native FMN identity and complete-model flag must agree with the native fields");
				final int start=p.parseInt(5);
				final int stop=p.parseInt(6);
				require(start>0 && stop>=start, "Native FMN cores require positive inclusive coordinates");
				final Event event=new Event(row, node, id, text(4), strand(7), start, stop,
					full && p.termEquals("CONFIDENT", 24));
				require(!p.termEquals("CONFIDENT", 24) || full, "A confident FMN event must retain its full native model/inclusion support");
				require(eventById.get(event.key)==null, "Repeated cohort-qualified native FMN event");
				eventById.put(event.key, event);
				events.add(event);
				counts[cohort]++;
				if(full){fulls[cohort]++;}
			}
		}
		for(int c=0; c<2; c++){require(counts[c]==expectedEvents[c] && fulls[c]==expectedFull[c], "FMN event output must conserve its accepted cohort summary");}
	}

	private void links(final String path){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		final long[] counts=new long[2];
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, MEMBER_LINK_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 35);
				final int cohort=cohort(0);
				final String oldId=text(1);
				final Node node=node(cohort, p.parseInt(6));
				final State state=stateById.get(oldId+'\t'+node.key);
				final int eventId=p.parseInt(5);
				require(state!=null && eventId>0 && seen.add(state.old.id+'\t'+node.key+'\t'+eventId), "Every native link must name one unique original member/input/event association");
				final Old old=state.old;
				final int start=p.parseInt(8);
				final int stop=p.parseInt(9);
				final int overlap=old.overlap(start, stop);
				require(p.parseInt(2)==old.length && p.parseInt(3)==overlap && text(7).equals(old.contig) && strand(10)==old.strand,
					"Native association must bind the literal contig, strand and disjoint old-core overlap");
				require(p.termEquals("SOURCE_LINE", 4) || p.termEquals("BOUND_CORE_COORDINATE", 4)
					|| p.termEquals("SOURCE_LINE_AND_BOUND_CORE", 4), "Saved link basis must remain explicit");
				final boolean full=full(13, 14, 15, 16, 22);
				final boolean own=p.termEquals("RF00050", 11);
				require(own==p.termEquals("FMN", 12), "Native model/family identity must agree on every member link");
				if(own){
					final Event event=eventById.get(key(cohort, eventId));
					require(event!=null, "FMN member association has no retained native event");
					p.setBounds(5);
					final int begin=p.a();
					eventParser.set(event.raw);
					eventParser.setBounds(2);
					require(equalTail(row, begin, event.raw, eventParser.a()), "FMN association fields differ from the retained native event");
					if(overlap>0 && event.eligible){
						old.nativeAssociations.add(event);
						event.oldIds.add(old.id);
					}
				}
				if(overlap>0 && full){
					if(own){state.own++;}
					else{
						state.other++;
						state.maxOther=Math.max(state.maxOther, overlap);
						old.competing=true;
					}
				}
				state.links++;
				counts[cohort]++;
			}
		}
		for(int c=0; c<2; c++){require(counts[c]==expectedLinks[c], "Member-event associations must conserve their accepted cohort summary");}
	}

	private void states(final String path){
		final long[] counts=new long[2];
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, OLD_EVIDENCE_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 12);
				final int cohort=cohort(0);
				final State state=stateById.get(text(1)+'\t'+key(cohort, p.parseInt(2)));
				require(state!=null && !state.seen && p.termEquals("NOT_EVALUATED", 11), "Every saved old-member state must be unique and bound");
				require(p.parseInt(5)==state.links && p.parseInt(6)==state.own && p.parseInt(7)==state.other
					&& p.parseInt(8)==state.maxOther && p.parseInt(9)==state.old.length,
					"Reconstructed native link counts/overlap must equal the independent old-member diagnostic");
				state.seen=true;
				counts[cohort]++;
			}
		}
		for(int c=0; c<2; c++){require(counts[c]==expectedStates[c], "Saved old-member state table must conserve its complete cohort summary");}
		for(final State state:allStates){require(state.seen, "Original member/input binding is missing its saved diagnostic");}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Admission and Reports   ----------------*/
	/*--------------------------------------------------------------*/

	private void resolve(){
		for(final Event event:events){
			if(!event.eligible){continue;}
			final String interval=event.start+"-"+event.stop;
			final Candidate candidate=add(new Candidate("CM_"+NAMES[event.node.cohort]+'_'+event.id, event.node,
				event.contig, event.strand, interval, event.stop-event.start+1, true, true, "CM_FULL_CLEAN", NAMES[event.node.cohort]+':'+event.id));
			candidate.oldIds.addAll(event.oldIds);
			event.candidate=candidate;
		}
		for(final Old old:olds){
			if(old.competing){competingOld++;}
			final String issue=replacementIssue(old);
			if(old.competing){old.admission="COMPETING_FULL_OTHER_FAMILY_REVIEW";}
			else if(!old.nativeAssociations.isEmpty() && issue==null){old.admission="REPLACED_BY_FULL_CLEAN_NATIVE_CORE";}
			else{
				final boolean intact=old.intactPrior();
				old.admission=intact ? "RETAINED_OLD_PRIOR_NOT_NEWLY_CM_CHECKED" : "RETAINED_OLD_FRAGMENT_OR_AMBIGUITY_SEED_ONLY";
				if(!old.nativeAssociations.isEmpty()){
					ambiguousOld++;
					old.admission="RETAINED_OLD_PRIOR_"+issue;
				}
				p.set(old.raw);
				old.candidate=add(new Candidate("OLD_"+old.id, old.node, old.contig, old.strand, text(13), old.length,
					intact, true, old.admission, old.id));
				old.candidate.oldIds.add(old.id);
			}
		}
	}

	/** Reuse the declared prior/own-hit correspondence rule from
	 * NcrnaReferenceAnchors.observe (single own hit, reciprocal-half physical overlap),
	 * not an invented similarity threshold. Exact duplicate native witnesses name
	 * one physical core; distinct cores or multipart priors remain review cases.
	 * Same-family proximity alone is never positive contradiction of an old prior.
	 */
	private String replacementIssue(final Old old){
		if(old.nativeAssociations.isEmpty()){return null;}
		if(old.intervals.length!=2){return "MULTIPART_NATIVE_ASSOCIATION";}
		final Event single=old.nativeAssociations.get(0);
		for(final Event event:old.nativeAssociations){
			if(event.start!=single.start || event.stop!=single.stop || event.strand!=single.strand
				|| !event.contig.equals(single.contig) || !event.node.fasta.equals(single.node.fasta)){
				return "MULTIPLE_NATIVE_CORES";
			}
		}
		final long overlap=old.overlap(single.start, single.stop);
		final long nativeLength=(long)single.stop-single.start+1;
		if(overlap<NcrnaGffFamilyGrader.MIN_OVERLAP_FRAC*old.length
			|| overlap<NcrnaGffFamilyGrader.MIN_OVERLAP_FRAC*nativeLength){return "WEAK_NATIVE_ASSOCIATION";}
		return null;
	}

	private Candidate add(final Candidate candidate){
		final String core=candidate.node.key+'\t'+candidate.contig+'\t'+candidate.strand+'\t'+candidate.intervals;
		final Candidate existing=byCore.get(core);
		if(existing!=null){
			// Extraction roles belong to the shared physical sequence. Independent
			// complete CM support can make that sequence an intact candidate without
			// changing a partial old source's own metadata or its review disposition.
			existing.intact|=candidate.intact;
			existing.seed|=candidate.seed;
			existing.sources.add(candidate.source);
			existing.bases.add(candidate.basis);
			return existing;
		}
		candidate.sources.add(candidate.source);
		candidate.bases.add(candidate.basis);
		byCore.put(core, candidate);
		candidates.add(candidate);
		return candidate;
	}

	private void write(final String prefix){
		long training=0;
		long trainingSeeds=0;
		long held=0;
		long uncertain=0;
		try(Output oldOut=new Output(prefix+".old.tsv.gz"); Output eventOut=new Output(prefix+".events.tsv.gz");
				Output extracts=new Output(prefix+".extract.tsv.gz"); Output coverageOut=new Output(prefix+".coverage.tsv");
				Output summary=new Output(prefix+".summary.tsv")){
			oldOut.writer.println(OLD_OUTPUT_HEADER);
			for(final Old old:olds){
				b.clear().append(old.raw).tab().append(old.admission).tab().append(old.node.split).tab().append(old.node.group).tab();
				final ObjectSet<String> ids=new ObjectSet<String>(String.class);
				if(old.candidate!=null){ids.add(old.candidate.id);}
				for(final Event event:old.nativeAssociations){ids.add(event.candidate.id);}
				appendSorted(ids);
				oldOut.writer.print(b.tab().append(old.competing).tab().append(!old.competing && old.intactPrior()).tab().append(!old.competing).nl());
			}
			eventOut.writer.println(FMN_HEADER+"\tadmission\tassembly_split\tassembly_group\textraction_member");
			for(final Event event:events){
				eventOut.writer.print(b.clear().append(event.raw).tab().append(event.eligible ? "FULL_CLEAN_NATIVE_CORE" : "NATIVE_NOT_INTACT_ADMISSION")
					.tab().append(event.node.split).tab().append(event.node.group).tab().append(event.candidate==null ? "NA" : event.candidate.id).nl());
			}
			extracts.writer.println(EXTRACT_HEADER);
			for(final Candidate candidate:candidates){
				final Node node=candidate.node;
				final Meta meta=node.meta;
				final String taxid=meta==null ? "NA" : meta.taxid;
				final String domain=meta==null ? "UNKNOWN" : meta.domain;
				final String phylum=meta==null ? "UNKNOWN" : meta.phylum;
				final boolean train=node.split.equals(TRAIN);
				if(train){require(meta!=null && positiveTaxid(taxid), "Training candidates need explicit positive taxonomy IDs before FASTA/header construction");}
				if(train && candidate.intact){training++;}
				if(train && candidate.seed){trainingSeeds++;}
				if(node.split.equals(HELD)){held++;}
				if(node.split.equals(UNRESOLVED)){uncertain++;}
				b.clear().append(candidate.id.getBytes(StandardCharsets.UTF_8)).tab().append(NAMES[node.cohort]).tab().append(node.id).tab().append(node.fastaBytes)
					.tab().append(node.group).tab().append(taxid.getBytes(StandardCharsets.UTF_8))
					.tab().append(domain.getBytes(StandardCharsets.UTF_8)).tab().append(phylum.getBytes(StandardCharsets.UTF_8))
					.tab().append(candidate.contigBytes).tab().append(candidate.strand).tab().append(candidate.intervals).tab().append(candidate.length).tab();
				appendSorted(candidate.bases);
				b.tab().append(node.split).tab().append(candidate.intact).tab().append(candidate.seed)
					.tab().append(train && candidate.intact).tab().append(train && candidate.seed).tab();
				appendSorted(candidate.oldIds);
				b.tab();
				appendSorted(candidate.sources);
				extracts.writer.print(b.append("\tPENDING_SEQUENCE_AND_FLANK_VALIDATION\n"));
				if(candidate.intact){
					final String key=domain+'\t'+phylum;
					Coverage coverage=byLineage.get(key);
					if(coverage==null){
						coverage=new Coverage(key);
						byLineage.put(key, coverage);
					lineages.add(coverage);
					}
					coverage.all.put(node.group, 1);
					if(train){coverage.train.put(node.group, 1);}
				}
			}
			coverageOut.writer.println("domain\tphylum\tintact_assembly_groups\ttraining_intact_assembly_groups\tstatus");
			int deficits=0;
			for(final Coverage coverage:lineages){
				final boolean deficit=coverage.train.size()==0;
				if(deficit){deficits++;}
				coverageOut.writer.print(b.clear().append(coverage.key.getBytes(StandardCharsets.UTF_8)).tab().append(coverage.all.size()).tab().append(coverage.train.size())
					.tab().append(deficit ? "REVIEW_NO_TRAINING_REPRESENTATIVE_NOT_AUTOMATIC_EXCEPTION" : "REPRESENTED").nl());
			}
			require(training<=trainingSeeds && trainingSeeds+held+uncertain==candidates.size(), "Disjoint assembly splits must conserve all selected cores, including seed-only priors");
			summary.writer.println("old_members\tnative_events\tunique_candidate_cores\ttraining_intact_cores\ttraining_seed_cores\theldout_cores\tunresolved_assembly_cores\told_competing_evidence\told_ambiguous_native_associations\tlineages_without_training\tstatus");
			summary.writer.print(b.clear().append(olds.size()).tab().append(events.size()).tab().append(candidates.size()).tab()
				.append(training).tab().append(trainingSeeds).tab().append(held).tab().append(uncertain).tab().append(competingOld).tab()
				.append(ambiguousOld).tab().append(deficits).append("\tREVIEW_AND_EXTRACTION_REQUIRED_NO_AUTOMATIC_COVERAGE_EXCEPTION\n"));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Parsing             ----------------*/
	/*--------------------------------------------------------------*/

	private void parse(final byte[] row, final int fields){
		p.set(row);
		require(p.terms()==fields, "Construction table field count changed: expected "+fields);
	}
	private String text(final int field){
		p.setBounds(field);
		final String result=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(result.indexOf('\ufffd')<0, "Invalid UTF-8 would change member provenance");
		return result;
	}
	private int cohort(final int field){
		require(p.termEquals("evaluation", field) || p.termEquals("discovery", field), "Unknown source cohort");
		return p.termEquals("evaluation", field) ? 0 : 1;
	}
	private Node node(final int cohort, final int id){
		final Node result=nodesById.get(key(cohort, id));
		require(result!=null, "Evidence names an absent cohort-qualified assembly input: "+key(cohort, id));
		return result;
	}
	private boolean flag(final int field){
		require(p.termEquals("true", field) || p.termEquals("false", field), "Saved boolean is not literal true/false");
		return p.termEquals("true", field);
	}
	private char strand(final int field){
		require(p.termEquals("+", field) || p.termEquals("-", field), "Construction strand must be explicit");
		return p.termEquals("+", field) ? '+' : '-';
	}
	private boolean full(final int fromField, final int toField, final int lengthField, final int truncationField, final int includedField){
		final int from=p.parseInt(fromField);
		final int to=p.parseInt(toField);
		final int length=p.parseInt(lengthField);
		final int truncation=p.parseInt(truncationField);
		require(from>0 && to>=from && to<=length && truncation>=0 && truncation<=3, "Native model span or truncation is invalid");
		return flag(includedField) && from==1 && to==length && truncation==0;
	}
	private int[] intervals(final int field){
		p.setBounds(field);
		final byte[] row=p.line();
		final IntList result=new IntList();
		long value=0;
		boolean digits=false;
		for(int i=p.a(); i<=p.b(); i++){
			final int c=i==p.b() ? ',' : row[i];
			if(c=='-' || c==','){
				require(digits && value>0 && value<=Integer.MAX_VALUE && (c=='-' ? result.size()%2==0 : result.size()%2==1), "Expected positive start-stop interval pairs");
				result.add((int)value);
				value=0;
				digits=false;
			}else{
				require(c>='0' && c<='9' && value<=Integer.MAX_VALUE/10L, "Invalid or overflowing interval coordinate");
				value=value*10+c-'0';
				digits=true;
			}
		}
		return result.toArray();
	}
	private void appendSorted(final ObjectSet<String> values){
		final ArrayList<String> sorted=new ArrayList<String>();
		for(final String value:values.keys()){if(value!=null){sorted.add(value);}}
		java.util.Collections.sort(sorted);
		if(sorted.isEmpty()){b.append("NA");}
		for(int i=0; i<sorted.size(); i++){
			if(i>0){b.append(';');}
			b.append(sorted.get(i).getBytes(StandardCharsets.UTF_8));
		}
	}
	private static boolean equalTail(final byte[] left, final int a, final byte[] right, final int b){
		if(left.length-a!=right.length-b){return false;}
		for(int i=0; i<left.length-a; i++){if(left[a+i]!=right[b+i]){return false;}}
		return true;
	}
	private static boolean positiveTaxid(final String value){
		if(value.isEmpty() || value.charAt(0)<'1' || value.charAt(0)>'9'){return false;}
		for(int i=1; i<value.length(); i++){if(value.charAt(i)<'0' || value.charAt(i)>'9'){return false;}}
		return true;
	}
	private static String key(final int cohort, final int id){return NAMES[cohort]+'\t'+id;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------        Retained State        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Meta{
		Meta(final String genome_, final String fasta_, final String gff_, final String taxid_, final String domain_, final String phylum_, final String assembly_){
			genome=genome_;
			fasta=fasta_;
			gff=gff_;
			taxid=taxid_;
			domain=domain_;
			phylum=phylum_;
			assembly=assembly_;
		}
		final String genome, fasta, gff, taxid, domain, phylum, assembly;
		CrossTax crossTax;
		final ArrayList<Old> members=new ArrayList<Old>();
		Node auditNode;
		boolean auditDiscovery;
		int auditOld;
	}
	private static final class Node{
		Node(final int cohort_, final int id_, final String fasta_, final int group_, final String split_, final Meta meta_){
			cohort=cohort_;
			id=id_;
			fasta=fasta_;
			group=group_;
			fastaBytes=fasta.getBytes(StandardCharsets.UTF_8);
			split=split_;
			meta=meta_;
			key=key(cohort, id);
		}
		final int cohort, id, group;
		final String fasta, split, key;
		final byte[] fastaBytes;
		final Meta meta;
		int auditCores, auditNativeCores, auditPriorCores;
		boolean auditIntact, auditNative;
		long[] stages;
	}
	private static final class Old{
		Old(final byte[] raw_, final String id_, final Meta meta_, final String contig_, final char strand_, final int[] intervals_,
				final int length_, final boolean partial_, final boolean pseudo_, final int undefined_, final int limit, final boolean fullAnnotation_){
			raw=raw_;
			id=id_;
			meta=meta_;
			contig=contig_;
			strand=strand_;
			intervals=intervals_;
			length=length_;
			partial=partial_;
			pseudo=pseudo_;
			undefined=undefined_;
			fullAnnotation=fullAnnotation_;
			int total=0;
			int previous=0;
			for(int i=0; i<intervals.length; i+=2){
				require(intervals[i]>previous && intervals[i+1]>=intervals[i] && intervals[i+1]<=limit, "Original intervals must be disjoint, ordered and inside their contig");
				total=Math.addExact(total, intervals[i+1]-intervals[i]+1);
				previous=intervals[i+1];
			}
			require(total==length && total>0 && !contig.isEmpty(), "Old core length must equal its physical disjoint intervals");
		}
		int overlap(final int start, final int stop){
			require(start>0 && stop>=start, "Linked native interval must be valid before measuring physical overlap");
			int result=0;
			for(int i=0; i<intervals.length; i+=2){result+=Math.max(0, Math.min(stop, intervals[i+1])-Math.max(start, intervals[i])+1);}
			assert(result<=length) : "Disjoint old-core intersections cannot exceed the old core";
			return result;
		}
		boolean intactPrior(){
			return fullAnnotation && intervals.length==2 && !partial && !pseudo && undefined==0;
		}
		final byte[] raw;
		final String id, contig;
		final Meta meta;
		final char strand;
		final int[] intervals;
		final int length, undefined;
		final boolean partial, pseudo, fullAnnotation;
		final ArrayList<State> states=new ArrayList<State>();
		final ArrayList<Event> nativeAssociations=new ArrayList<Event>();
		Node node;
		Candidate candidate;
		boolean competing;
		String admission;
	}
	private static final class State{
		State(final Old old_, final Node node_){
			old=old_;
			node=node_;
		}
		final Old old;
		final Node node;
		int links, own, other, maxOther;
		boolean seen;
	}
	private static final class Event{
		Event(final byte[] raw_, final Node node_, final int id_, final String contig_, final char strand_, final int start_, final int stop_, final boolean eligible_){
			raw=raw_;
			node=node_;
			id=id_;
			contig=contig_;
			strand=strand_;
			start=start_;
			stop=stop_;
			eligible=eligible_;
			key=key(node.cohort, id);
		}
		final byte[] raw;
		final Node node;
		final int id, start, stop;
		final String contig, key;
		final char strand;
		final boolean eligible;
		final ObjectSet<String> oldIds=new ObjectSet<String>(String.class);
		Candidate candidate;
	}
	private static final class Candidate{
		Candidate(final String id_, final Node node_, final String contig_, final char strand_, final String intervals_,
				final int length_, final boolean intact_, final boolean seed_, final String basis_, final String source_){
			id=id_;
			node=node_;
			contig=contig_;
			strand=strand_;
			intervals=intervals_;
			length=length_;
			intact=intact_;
			contigBytes=contig.getBytes(StandardCharsets.UTF_8);
			seed=seed_;
			basis=basis_;
			source=source_;
		}
		final String id, contig, intervals, basis, source;
		final byte[] contigBytes;
		final Node node;
		final char strand;
		final int length;
		boolean intact, seed;
		final ObjectSet<String> oldIds=new ObjectSet<String>(String.class);
		final ObjectSet<String> sources=new ObjectSet<String>(String.class);
		final ObjectSet<String> bases=new ObjectSet<String>(String.class);
	}
	private static final class Coverage{
		Coverage(final String key_){key=key_;}
		final String key;
		final IntHashMap all=new IntHashMap();
		final IntHashMap train=new IntHashMap();
	}
	private static final class AuditCoverage{
		AuditCoverage(final String key_){key=key_;}
		final String key;
		final long[] counts=new long[12];
		final IntHashMap intact=new IntHashMap();
		final IntHashMap trainIntact=new IntHashMap();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fields and Schemas     ----------------*/
	/*--------------------------------------------------------------*/

	private final LineParser1 p=new LineParser1('\t');
	private final LineParser1 eventParser=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private final ArrayList<Meta> metadata=new ArrayList<Meta>();
	private final ArrayList<Node> nodes=new ArrayList<Node>();
	private final ArrayList<Old> olds=new ArrayList<Old>();
	private final ArrayList<State> allStates=new ArrayList<State>();
	private final ArrayList<Event> events=new ArrayList<Event>();
	private final ArrayList<Candidate> candidates=new ArrayList<Candidate>();
	private final ArrayList<Coverage> lineages=new ArrayList<Coverage>();
	private final ObjectMap<String, Meta> byGenome=new ObjectMap<String, Meta>(String.class, Meta.class);
	private final ObjectMap<String, Meta> byFasta=new ObjectMap<String, Meta>(String.class, Meta.class);
	private final ObjectMap<String, Node> nodesById=new ObjectMap<String, Node>(String.class, Node.class);
	private final ObjectMap<String, Old> oldById=new ObjectMap<String, Old>(String.class, Old.class);
	private final ObjectMap<String, State> stateById=new ObjectMap<String, State>(String.class, State.class);
	private final ObjectMap<String, Event> eventById=new ObjectMap<String, Event>(String.class, Event.class);
	private final ObjectMap<String, Candidate> byCore=new ObjectMap<String, Candidate>(String.class, Candidate.class);
	private final ObjectMap<String, Coverage> byLineage=new ObjectMap<String, Coverage>(String.class, Coverage.class);
	private final long[] expectedStates=new long[2];
	private final long[] expectedLinks=new long[2];
	private final long[] expectedEvents=new long[2];
	private final long[] expectedFull=new long[2];
	private int competingOld;
	private int ambiguousOld;
	private static final String[] NAMES={"evaluation", "discovery"};
	private static final String TRAIN="TRAIN";
	private static final String HELD="EVALUATION_ASSEMBLY_HELD_OUT";
	private static final String UNRESOLVED="UNRESOLVED_ASSEMBLY_NOT_TRAINING";
	public static final String GENOME_HEADER="index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs";
	public static final String GROUP_HEADER=NcrnaConstructionAssemblyGroups.INPUT_HEADER+"\tresolved_group_id\tresolved_evaluation_group\tgroup_has_unresolved_pair\tresolved_assembly_disposition";
	public static final String FMN_HEADER="cohort\tnative_included_full_model\t"+NcrnaConstructionEvidence.EVENT_HEADER;
	public static final String MEMBER_LINK_HEADER="cohort\told_member_id\told_core_length\tcore_overlap_bp\tlink_basis\t"+NcrnaConstructionEvidence.EVENT_HEADER;
	public static final String SUMMARY_INPUT_HEADER="cohort\told_member_bindings\tmissing_source_bindings\tdistinct_member_event_links\tall_native_events\tFMN_native_events\tnative_full_GA_FMN_events\tstatus";
	public static final String OLD_EVIDENCE_HEADER="cohort\told_member_id\tinput_id\tsource_rows_expected\tsource_rows_observed\tdistinct_linked_events\tfull_GA_FMN\tfull_GA_other_family\tmax_other_core_overlap_bp\told_core_length\tclassification\teligibility_status";
	public static final String EXTRACT_HEADER="member_id\tcohort\tinput_id\tfasta\tassembly_group\ttaxid\tdomain\tphylum\tseqid\tstrand\tintervals1\tlength\tadmission\tassembly_split\tintact_metadata\tseed_metadata\ttrain_intact\ttrain_seed\told_members\tevidence_sources\tsequence_validation";
	public static final String OLD_OUTPUT_HEADER=NcrnaConstructionEvidence.MEMBER_HEADER+"\tadmission\tassembly_split\tassembly_group\textraction_members\tcompeting_full_other_family\told_prior_intact_metadata\told_prior_seed_metadata";
	public static final String AUDIT_HEADER="domain\tphylum\tevaluation_inputs\tdiscovery_scanned_inputs\ttraining_inputs\tdiscovery_inputs_with_intact_candidate\ttraining_inputs_with_intact_candidate\tdiscovery_inputs_with_native_clean_core\ttraining_inputs_with_native_clean_core\told_members\told_members_on_discovery_inputs\told_members_on_training_inputs\told_member_genomes\told_member_training_genomes\tintact_assembly_groups\ttraining_intact_assembly_groups";
}
