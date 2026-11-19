package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaConstructionEvidence;
import prok.NcrnaConstructionMembers;
import prok.BuildNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal original/native conservation, corrected cores, aliases and prior retention.
 * Paths deliberately do not exist: admission must not read genome payloads.
 * @author Ganyu
 */
public final class NcrnaConstructionMembersTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionMembersTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown member fixture argument: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_directory without sampling or overwrite");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Member fixture requires a fresh output directory");
			new NcrnaConstructionMembersTest(root.getPath()).run();
			pp.outstream.println("NCRNA_CONSTRUCTION_MEMBERS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaConstructionMembersTest(final String root_){root=root_;}

	private void run(){
		String genomes=NcrnaConstructionMembers.GENOME_HEADER+'\n';
		String bindings=NcrnaConstructionEvidence.BINDING_HEADER+'\n';
		for(int i=0; i<7; i++){
			genomes+=i+"\tg"+i+"\tGCF_"+i+".1\t"+(100+i)+"\tbacteria\t"+phylum(i)
				+'\t'+fasta(i)+".gff\t"+fasta(i)+"\t0\t0\t0\n";
			bindings+="g"+i+'\t'+fasta(i)+'\t'+fasta(i)+".gff\t"+cohort(i)+'\t'+input(i)+"\tLITERAL_CATALOG_PATH_NOT_SEQUENCE_PROOF\n";
		}
		final String groups=NcrnaConstructionMembers.GROUP_HEADER+'\n'
			+group(0, 0, true, "EVALUATION")+group(5, 1, true, "EVALUATION")
			+group(1, 2, false, "TRAIN_CANDIDATE_ASSEMBLY_IDENTITY")
			+group(2, 0, true, "EXCLUDE_VERIFIED_EVALUATION_ASSEMBLY_ALIAS")
			+group(3, 4, false, "ASSEMBLY_SEQUENCE_IDENTITY_UNRESOLVED")
			+group(4, 5, false, "TRAIN_CANDIDATE_ASSEMBLY_IDENTITY")
			+group(6, 7, false, "TRAIN_CANDIDATE_ASSEMBLY_IDENTITY")
			+"evaluation\t3\t/never-opened/unlabelled.fa\t/never-opened/unlabelled.fa\tNA\t6\ttrue\tfalse\tfalse\tORIGINAL_DIAGNOSTIC\t6\ttrue\tfalse\tEVALUATION\n";
		final String old=NcrnaConstructionEvidence.MEMBER_HEADER+'\n'
			+old("e", 0, 10, 19, false)+old("a", 1, 10, 19, false)
			+old("b", 2, 10, 19, false)+old("u", 3, 30, 39, false)
				.replace("\t30-39\t10\t", "\t30-34,40-44\t10\t").replace("ANNOTATION_FULL_COPY", "MULTIPART_UNRESOLVED")
			+old("p", 1, 40, 44, true)+old("c", 4, 50, 59, false)+old("r", 5, 60, 69, false)
			+old("w", 1, 200, 209, false)+old("m", 4, 300, 319, false)+old("q", 4, 400, 419, true);
		final String e=event(1, 0, 9, 20, false, true);
		final String a=event(1, 1, 9, 20, false, true);
		final String duplicate=event(2, 1, 9, 20, false, true);
		final String partial=event(3, 1, 80, 84, true, true);
		final String addition=event(4, 4, 90, 99, false, true);
		final String alias=event(5, 2, 100, 109, false, true);
		final String unresolved=event(6, 3, 100, 109, false, true);
		final String gap=event(7, 3, 35, 39, false, true);
		final String weak=event(8, 1, 209, 220, false, true);
		final String multiOne=event(9, 4, 299, 320, false, true);
		final String multiTwo=event(10, 4, 300, 321, false, true);
		final String partialExact=event(11, 4, 400, 419, false, true);
		final String partialOther=event(12, 4, 399, 420, false, true);
		final String competitor=event(100, 4, 50, 59, false, false);
		final String events=NcrnaConstructionMembers.FMN_HEADER+'\n'
			+"evaluation\ttrue\t"+e+"discovery\ttrue\t"+a+"discovery\ttrue\t"+duplicate
			+"discovery\tfalse\t"+partial+"discovery\ttrue\t"+addition
			+"discovery\ttrue\t"+alias+"discovery\ttrue\t"+unresolved+"discovery\ttrue\t"+gap
			+"discovery\ttrue\t"+weak+"discovery\ttrue\t"+multiOne+"discovery\ttrue\t"+multiTwo
			+"discovery\ttrue\t"+partialExact+"discovery\ttrue\t"+partialOther;
		final String links=NcrnaConstructionMembers.MEMBER_LINK_HEADER+'\n'
			+link("e", 0, 10, 10, e)+link("a", 1, 10, 10, a)+link("a", 1, 10, 10, duplicate)
			+link("c", 4, 10, 10, competitor)+link("u", 3, 10, 0, gap)
			+link("w", 1, 10, 1, weak)+link("m", 4, 20, 20, multiOne)+link("m", 4, 20, 20, multiTwo)
			+link("q", 4, 20, 20, partialExact)+link("q", 4, 20, 20, partialOther);
		final String states=NcrnaConstructionMembers.OLD_EVIDENCE_HEADER+'\n'
			+state("e", 0, 10, 1, 1, 0)+state("a", 1, 10, 2, 2, 0)
			+state("b", 2, 10, 0, 0, 0)+state("u", 3, 10, 1, 0, 0)
			+state("p", 1, 5, 0, 0, 0)+state("c", 4, 10, 1, 0, 1)+state("r", 5, 10, 0, 0, 0)
			+state("w", 1, 10, 1, 1, 0)+state("m", 4, 20, 2, 2, 0)+state("q", 4, 20, 2, 2, 0);
		final String summary=NcrnaConstructionMembers.SUMMARY_INPUT_HEADER+'\n'
			+"evaluation\t2\t0\t1\t1\t1\t1\tEVIDENCE_JOIN_NOT_ADMISSION_OR_ALIAS_PROOF\n"
			+"discovery\t8\t0\t9\t13\t12\t11\tEVIDENCE_JOIN_NOT_ADMISSION_OR_ALIAS_PROOF\n";
		write("genomes.tsv", genomes);
		write("groups.tsv", groups);
		write("original.members.tsv.gz", old);
		write("original.bindings.tsv", bindings);
		write("evidence.fmn_events.tsv.gz", events);
		write("evidence.links.tsv.gz", links);
		write("evidence.old_evidence.tsv.gz", states);
		write("evidence.summary.tsv", summary);
		process("result");
		checkResult("result");
		NcrnaConstructionMembers.main(new String[]{"in="+path("groups.tsv"), "genomes="+path("genomes.tsv"),
			"bindings="+path("original"), "evidence="+path("evidence"), "out="+path("cli"),
			"expectedgenomes=7", "expectedmembers=10", "expectedeval=3", "expectednew=5", "t=1"});
		checkResult("cli");
		check(rows(path("result.summary.tsv")).equals(rows(path("cli.summary.tsv"))), "CLI and library must produce the same complete ledger accounting");
		audit("audit");
		final ArrayList<String> coverage=rows(path("audit.coverage.tsv"));
		check(coverage.size()==5, "Coverage audit includes both zero-positive and unlabelled categories");
		check(coverage.contains("bacteria\tCommon\t1\t4\t2\t4\t2\t4\t2\t9\t8\t6\t5\t2\t4\t2"),
			"Input denominators, old members, old-positive genomes and resolved intact groups have independently specified units");
		check(coverage.contains("bacteria\tRare\t1\t0\t0\t0\t0\t0\t0\t1\t0\t0\t1\t0\t1\t0"),
			"The held-out rare category retains its old member without creating training evidence");
		check(coverage.contains("bacteria\tQuiet\t0\t1\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"),
			"A scanned training genome with no FMN signal remains in its category's denominator");
		check(coverage.contains("UNKNOWN\tUNKNOWN\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"),
			"Missing taxonomy is explicit even when the input has no candidate core");
		check(rows(path("audit.summary.tsv")).get(1).equals("7\t3\t5\t3\t10\t17\t1\tORIGINAL_LABELS_ONLY_NO_TAXONOMY_OR_SPLIT_CHANGE"),
			"Audit conservation covers every input, original member and saved candidate");
		check(rows(path("audit.unlabelled.tsv.gz")).get(1).equals("evaluation\t3\t/never-opened/unlabelled.fa\t6\tEVALUATION_ASSEMBLY_HELD_OUT\t0\tfalse\tfalse\tNO_LITERAL_FASTA_BINDING_IN_ORIGINAL_ROSTER_NOT_PROOF_OF_UNSEEN_ASSEMBLY"),
			"Unlabelled input inventory preserves its exact source and group without inferring novelty");
		NcrnaConstructionMembers.main(new String[]{"mode=coverage", "in="+path("result"), "groups="+path("groups.tsv"),
			"genomes="+path("genomes.tsv"), "out="+path("audit_cli"), "expectedgenomes=7", "expectedmembers=10",
			"expectedeval=3", "expectednew=5", "expectedcores=17", "t=1"});
		check(coverage.equals(rows(path("audit_cli.coverage.tsv")))
			&& rows(path("audit.summary.tsv")).equals(rows(path("audit_cli.summary.tsv"))), "Coverage CLI and library reports agree");
		final String savedCores=String.join("\n", rows(path("result.extract.tsv.gz")))+'\n';
		final String savedOld=String.join("\n", rows(path("result.old.tsv.gz")))+'\n';
		auditFault("audit_taxonomy", "result.extract.tsv.gz", savedCores, savedCores.replace("\tbacteria\tCommon\t", "\tbacteria\tOther\t"));
		auditFault("audit_source", "result.extract.tsv.gz", savedCores, savedCores.replace(fasta(1), "/changed/α.fa"));
		auditFault("audit_truncation", "result.extract.tsv.gz", savedCores, savedCores.substring(0, savedCores.lastIndexOf('\n', savedCores.length()-2)+1));
		auditFault("audit_old_duplicate", "result.old.tsv.gz", savedOld, savedOld+rows(path("result.old.tsv.gz")).get(1)+'\n');
		auditFault("audit_old_split", "result.old.tsv.gz", savedOld, savedOld.replace("\tTRAIN\t", "\tEVALUATION_ASSEMBLY_HELD_OUT\t"));

		fault("missing_event", "evidence.fmn_events.tsv.gz", events, events.replace("discovery\ttrue\t"+unresolved, ""));
		fault("wrong_full_flag", "evidence.fmn_events.tsv.gz", events, events.replace("evaluation\ttrue\t", "evaluation\tfalse\t"));
		fault("duplicate_event", "evidence.fmn_events.tsv.gz", events, events+"evaluation\ttrue\t"+e);
		fault("wrong_binding_path", "original.bindings.tsv", bindings, bindings.replace(fasta(1), "/changed/α.fa"));
		fault("missing_link", "evidence.links.tsv.gz", links, links.replace(link("a", 1, 10, 10, duplicate), ""));
		fault("changed_native_fields", "evidence.links.tsv.gz", links, links.replace("\t100.0\t66.0\t", "\t101.0\t66.0\t"));
		fault("wrong_diagnostic", "evidence.old_evidence.tsv.gz", states,
			states.replace(state("a", 1, 10, 2, 2, 0), state("a", 1, 10, 1, 1, 0)));
		fault("missing_diagnostic", "evidence.old_evidence.tsv.gz", states, states.replace(state("u", 3, 10, 1, 0, 0), ""));
		fault("gap_as_core", "evidence.links.tsv.gz", links, links.replace(link("u", 3, 10, 0, gap), link("u", 3, 10, 5, gap)));
		fault("missing_old", "original.members.tsv.gz", old, old.replace(old("p", 1, 40, 44, true), ""));
		fault("invalid_interval", "original.members.tsv.gz", old, old.replace("\t10-19\t", "\t19-10\t"));
		fault("false_training_alias", "groups.tsv", groups,
			groups.replace("EXCLUDE_VERIFIED_EVALUATION_ASSEMBLY_ALIAS", "TRAIN_CANDIDATE_ASSEMBLY_IDENTITY"));
		fault("missing_cohort_summary", "evidence.summary.tsv", summary, summary.substring(0, summary.indexOf("discovery\t")));
		stages(events, summary);
		crossSelection(genomes);
	}

	private void crossSelection(final String genomes){
		write("cross.genomes.tsv", bacteroidota(genomes));
		write("cross.old.tsv.gz", bacteroidota(String.join("\n", rows(path("result.old.tsv.gz")))+'\n'));
		write("cross.extract.tsv.gz", bacteroidota(String.join("\n", rows(path("result.extract.tsv.gz")))+'\n'));
		String taxonomy="";
		for(int i=0; i<6; i++){
			final String clazz=i==0 || i==2 ? "A" : i==4 ? "C" : i==5 ? "" : "Béta";
			taxonomy+="RS_GCF_"+i+".1\td__Bacteria;p__Bacteroidota;c__"+clazz+";o__Order"+clazz+";f__F;g__G;s__S\n";
		}
		write("cross.taxonomy.tsv", taxonomy);
		cross("cross_selection", 1);
		final ArrayList<String> summary=rows(path("cross_selection.summary.tsv"));
		check(summary.get(1).equals("2\t5\t1\t1\t1\t2\tCROSSED_PROCEDURE_SELECTION_NOT_SEARCHED_OR_SPLIT_CHANGE"),
			"Full Bacteroidota denominators include missing rank/accession inputs; each arm selects one assembly with two native cores total");
		final ArrayList<String> selected=rows(path("cross_selection.selected.tsv"));
		check(selected.size()==3 && selected.get(1).startsWith("evaluation\t1\t"+fasta(0)+'\t')
			&& selected.get(2).startsWith("discovery\t2\t"+fasta(2)+'\t'),
			"Class/order rounds choose the stable lowest input ID and never select an old-only evaluation prior");
		check(rows(path("cross_selection.evaluation.paths.tsv")).equals(java.util.Arrays.asList(fasta(0)+"\t1"))
			&& rows(path("cross_selection.discovery.paths.tsv")).equals(java.util.Arrays.asList(fasta(2)+"\t2")),
			"Native batch preparation consumes exact selected paths with original input IDs, not guessed names");
		check(rows(path("cross_selection.expected.tsv")).size()==3, "Both selected native cores preserve their complete physical evidence");
		boolean unicode=false;
		for(final String row:rows(path("cross_selection.composition.tsv"))){unicode|=row.contains("Béta");}
		check(unicode, "GTDB class/order text must retain UTF-8 through the native output writer");
		NcrnaConstructionMembers.main(new String[]{"mode=cross-select", "in="+path("cross"), "groups="+path("groups.tsv"),
			"genomes="+path("cross.genomes.tsv"), "taxonomy="+path("cross.taxonomy.tsv"), "out="+path("cross_cli"),
			"expectedgenomes=7", "expectedmembers=10", "expectedeval=3", "expectednew=5", "expectedcores=17", "crosscount=1", "t=1"});
		for(final String suffix:new String[]{".selected.tsv", ".composition.tsv", ".expected.tsv", ".summary.tsv", ".evaluation.paths.tsv", ".discovery.paths.tsv"}){
			check(rows(path("cross_selection"+suffix)).equals(rows(path("cross_cli"+suffix))), "Crossed selection CLI parity: "+suffix);
		}
		crossFault("cross_duplicate_taxonomy", taxonomy+taxonomy.substring(0, taxonomy.indexOf('\n')+1));
		crossFault("cross_wrong_phylum", taxonomy.replace("p__Bacteroidota", "p__Foreign"));
		write("cross.taxonomy.tsv", taxonomy);
		boolean rejected=false;
		try{cross("cross_insufficient", 2);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Insufficient native-positive evaluation groups must not be filled with old-only priors");
		lineageComparison(taxonomy);
		crossReport();
	}
	private void lineageComparison(final String originalTaxonomy){
		final String taxonomy=originalTaxonomy.replace("g__G", "g__Géna")
			.replace("c__;o__Order;f__F;g__Géna", "c__;o__Order;f__OldOnly;g__Absent");
		write("cross.taxonomy.tsv", taxonomy);
		lineage("lineages");
		check(rows(path("lineages.summary.tsv")).get(1).equals("2\t5\t1\t0\t1\t0\t0\t1\t0\tSAVED_METADATA_NOT_BIOLOGICAL_ABSENCE_OR_FROZEN_SPLIT"),
			"Lineage denominators include five discovery inputs, one native-positive evaluation group, one unmatched accession and one unbound evaluation input");
		final ArrayList<String> table=rows(path("lineages.lineages.tsv"));
		check(table.contains("GENUS\tF\tGéna\tGTDB_ACCESSION_MATCH\t1\t1\t1\t1\t4\t4\tDISCOVERY_PRESENT"),
			"Shared genus preserves UTF8, one evaluation native core/group (the other selected core is discovery), and all four matched discovery groups");
		check(table.contains("GENUS\tOldOnly\tAbsent\tGTDB_ACCESSION_MATCH\t1\t1\t0\t0\t0\t0\tNO_MATCHED_DISCOVERY_INPUT"),
			"A discovery-absent old-only lineage must not become a native-positive group or justify a training exception");
		check(table.contains("GENUS\tUNKNOWN\tUNKNOWN\tNO_ACCESSION_MATCH\t0\t0\t0\t0\t1\t1\tUNKNOWN_TAXONOMY"),
			"Missing accession taxonomy remains unknown rather than biological or corpus absence");
		check(rows(path("lineages.lineage_inputs.tsv.gz")).size()==8, "All seven labelled Bacteroidota input records must remain auditable");
		NcrnaConstructionMembers.main(new String[]{"mode=lineages", "in="+path("cross"), "groups="+path("groups.tsv"),
			"genomes="+path("cross.genomes.tsv"), "taxonomy="+path("cross.taxonomy.tsv"), "out="+path("lineages_cli"),
			"expectedgenomes=7", "expectedmembers=10", "expectedeval=3", "expectednew=5", "expectedcores=17", "t=1"});
		for(final String suffix:new String[]{".lineages.tsv", ".lineage_inputs.tsv.gz", ".summary.tsv"}){
			check(rows(path("lineages"+suffix)).equals(rows(path("lineages_cli"+suffix))), "Lineage CLI must preserve the same complete metadata census: "+suffix);
		}
		write("cross.taxonomy.tsv", taxonomy.replace("RS_GCF_2.1\td__Bacteria;p__Bacteroidota;c__A;o__OrderA;f__F;g__Géna", "RS_GCF_2.1\td__Bacteria;p__Bacteroidota;c__A;o__OrderA;f__F;g__Conflicting"));
		boolean rejected=false;
		try{lineage("lineage_conflict");}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("conflicting GTDB lineage");}
		check(rejected, "Aliases of one assembly cannot silently enter conflicting genus bins");
		write("cross.taxonomy.tsv", originalTaxonomy);
	}
	private void lineage(final String output){
		NcrnaConstructionMembers.lineages(path("cross"), path("groups.tsv"), path("cross.genomes.tsv"),
			path("cross.taxonomy.tsv"), path(output), 7, 10, 3, 5, 17);
	}
	private void crossReport(){
		write("cross.models.tsv", "accession\tname\tCLEN\nRF00050\tFMN\t138\nRF00010\tRNaseP\t10\n");
		String collection="cohort\tinput_id\tprepared\tnative\n";
		for(final int genome:new int[]{0,2}){
			final String prepared="cross_prepared_"+genome;
			final String nativeDir="cross_native_"+genome;
			check(new File(path(prepared)).mkdir() && new File(path(nativeDir)).mkdir(), "Crossed native fixture directories must be fresh");
			write(prepared+"/members.tsv", "input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records\n"
				+input(genome)+'\t'+fasta(genome)+"\t100\t1000\t1\n");
			write(prepared+"/targets.tsv", "batch_id\ttarget_id\tinput_id\tcontig_ordinal\toriginal_target\tlength\n"
				+"0\ti"+genome+"c1\t"+input(genome)+"\t1\tcontig"+genome+"\t1000\n");
			write(nativeDir+"/exit.txt", "0\n");
			final int from=genome==0 ? 9 : 100;
			final int to=genome==0 ? 20 : 109;
			String hits=rawHit(genome, from, to, false, "RF00050", "FMN").replace(" cm 1 10 ", " cm 1 138 ");
			if(genome==2){hits+=rawHit(genome, 400, 509, false, "RF00050", "FMN").replace(" cm 1 10 ", " cm 1 138 ").replace(" ! ", " ? ");}
			write(nativeDir+"/hits.tbl", "# Version: 1.1.5 fixture\n"+hits+"# [ok]\n");
			write(nativeDir+"/search.out", crossNativeReport(genome==0 ? 1 : 2, genome==0));
			collection+=cohort(genome)+'\t'+input(genome)+'\t'+path(prepared)+'\t'+path(nativeDir)+'\n';
		}
		write("cross.collection.tsv", collection);
		report("cross_report");
		final ArrayList<String> summary=rows(path("cross_report.summary.tsv"));
		check(summary.size()==3 && summary.get(1).startsWith("evaluation\t1\t1\t1\t1\t1\t1\t1\t1\t")
			&& summary.get(2).startsWith("discovery\t1\t2\t1\t1\t1\t1\t1\t1\t"),
			"Native sub-GA reports stay separate from original GA inclusion; expected physical cores remain bound by cohort and coordinates");
		NcrnaConstructionMembers.main(new String[]{"mode=cross-report", "in="+path("cross_selection"),
			"raw="+path("cross.collection.tsv"), "models="+path("cross.models.tsv"), "out="+path("cross_report_cli"), "crosscount=1", "t=1"});
		for(final String suffix:new String[]{".inputs.tsv", ".hits.tsv", ".expected.tsv", ".summary.tsv"}){
			check(rows(path("cross_report"+suffix)).equals(rows(path("cross_report_cli"+suffix))), "Cross report native CLI parity: "+suffix);
		}
		final String report=crossNativeReport(2, false);
		crossReportFault("cross_missing_model", "cross_native_0/search.out", crossNativeReport(1, true), crossNativeReport(1, false));
		crossReportFault("cross_wrong_size", "cross_native_2/search.out", report, report.replace("2000 residues", "1999 residues"));
		crossReportFault("cross_wrong_hit_count", "cross_native_2/search.out", report, crossNativeReport(3, false));
		crossReportFault("cross_missing_complete", "cross_native_2/search.out", report, report.replace("[ok]\n", ""));
		crossReportFault("cross_failed_search", "cross_native_2/exit.txt", "0\n", "1\n");
		crossReportFault("cross_missing_point", "cross.collection.tsv", collection, collection.substring(0, collection.lastIndexOf("discovery\t")));
		final String target="batch_id\ttarget_id\tinput_id\tcontig_ordinal\toriginal_target\tlength\n0\ti2c1\t2\t1\tcontig2\t1000\n";
		crossReportFault("cross_wrong_input", "cross_prepared_2/targets.tsv", target, target.replace("\t2\t1\t", "\t99\t1\t"));
		final String hits=String.join("\n", rows(path("cross_native_2/hits.tbl")))+'\n';
		write("cross_native_2/hits.tbl", "# Version: 1.1.5 fixture\n# [ok]\n");
		write("cross_native_2/search.out", crossNativeReport(0, false));
		report("cross_zero");
		check(rows(path("cross_zero.summary.tsv")).get(2).startsWith("discovery\t1\t0\t0\t0\t1\t0\t0\t0\t"),
			"A completed zero-hit search is retained as zero with its unrecovered expected core, never dropped");
		write("cross_native_2/hits.tbl", hits);
		write("cross_native_2/search.out", report);
		procedureReport();
	}
	private void procedureReport(){
		check(new File(path("procedure_bundle")).mkdir(), "Paired procedure fixture needs a distinct native directory");
		write("procedure_bundle/exit.txt", "0\n");
		write("procedure_bundle/hits.tbl", "# Version: 1.1.5 fixture\n"
			+rawHit(2, 100, 109, false, "RF00050", "FMN").replace(" cm 1 10 ", " cm 1 138 ")+"# [ok]\n");
		write("procedure_bundle/search.out", crossNativeReport(1, true));
		write("procedure_selected.tsv", "cohort\tinput_id\tfasta\tassembly_group\tassembly\tgenus\tCM_admitted_cores\n"
			+"discovery\t2\t"+fasta(2)+"\t0\tGCF_2.1\tChryseobacterium\t1\n");
		final String collection="procedure\tinput_id\tprepared\tnative\nDISCOVERY_BUNDLE\t2\t"+path("cross_prepared_2")+'\t'+path("procedure_bundle")
			+"\nEVALUATION_FMN\t2\t"+path("cross_prepared_2")+'\t'+path("cross_native_2")+'\n';
		write("procedure_collection.tsv", collection);
		NcrnaConstructionMembers.procedureReport(path("procedure_selected.tsv"), path("procedure_collection.tsv"), path("cross.models.tsv"), path("procedures"), 1);
		final ArrayList<String> summary=rows(path("procedures.summary.tsv"));
		check(summary.get(1).equals("DISCOVERY_BUNDLE\t1\t1\t1\t1\t1\tPAIRED_DISCOVERY_INPUT_PROCEDURES_NOT_SPLIT_AUTHORITY")
			&& summary.get(2).equals("EVALUATION_FMN\t1\t2\t1\t1\t1\tPAIRED_DISCOVERY_INPUT_PROCEDURES_NOT_SPLIT_AUTHORITY"),
			"Both native procedures must report the same discovery input without relabelling its cohort; sub-GA stays separate");
		NcrnaConstructionMembers.main(new String[]{"mode=procedure-report", "in="+path("procedure_selected.tsv"), "raw="+path("procedure_collection.tsv"),
			"models="+path("cross.models.tsv"), "out="+path("procedures_cli"), "crosscount=1", "t=1"});
		check(rows(path("procedures.inputs.tsv")).equals(rows(path("procedures_cli.inputs.tsv"))), "Paired-procedure CLI must retain actual source identity and counts");
		write("procedure_missing.tsv", collection.substring(0, collection.indexOf("EVALUATION_FMN")));
		boolean rejected=false;
		try{NcrnaConstructionMembers.procedureReport(path("procedure_selected.tsv"), path("procedure_missing.tsv"), path("cross.models.tsv"), path("procedures_missing"), 1);}
		catch(IllegalArgumentException failure){rejected=failure.getMessage().contains("missing or failed procedure");}
		check(rejected, "An omitted arm must not become a zero-hit result");
	}
	private static String crossNativeReport(final int hits, final boolean bundle){
		check(hits>=0, "Native fixture reports conserve nonnegative hit counts");
		return "Query: FMN [CLEN=138]\nAccession: RF00050\nTarget sequences: 1 (2000 residues searched)\nTotal CM hits reported: "+hits+" (fixture)\n"
			+(bundle ? "Query: RNaseP [CLEN=10]\nAccession: RF00010\nTarget sequences: 1 (2000 residues searched)\nTotal CM hits reported: 0 (fixture)\n" : "")+"[ok]\n";
	}
	private void report(final String output){
		NcrnaConstructionMembers.crossReport(path("cross_selection"), path("cross.collection.tsv"), path("cross.models.tsv"), path(output), 1);
	}
	private void crossReportFault(final String name, final String file, final String valid, final String invalid){
		check(!valid.equals(invalid), "Cross report fault must change the native evidence");
		write(file, invalid);
		boolean rejected=false;
		try{report(name);}
		catch(IllegalArgumentException expected){rejected=true;}
		finally{write(file, valid);}
		check(rejected, "Incomplete native cross evidence must not become a false zero: "+name);
	}
	private static String bacteroidota(final String text){
		check(text!=null, "Fixture relabelling requires actual saved metadata");
		return text.replace("\tCommon\t", "\tBacteroidota\t").replace("\tRare\t", "\tBacteroidota\t").replace("\tQuiet\t", "\tBacteroidota\t");
	}
	private void cross(final String output, final int count){
		NcrnaConstructionMembers.crossSelect(path("cross"), path("groups.tsv"), path("cross.genomes.tsv"),
			path("cross.taxonomy.tsv"), path(output), 7, 10, 3, 5, 17, count);
	}
	private void crossFault(final String name, final String taxonomy){
		write("cross.taxonomy.tsv", taxonomy);
		boolean rejected=false;
		try{cross(name, 1);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Original accession taxonomy inconsistency must fail: "+name);
	}

	private void stages(final String events, final String summary){
		check(new File(path("raw")).mkdir() && new File(path("prepared")).mkdir(), "Stage fixture needs fresh native/map directories");
		final String table=path("raw/hits.tbl");
		final int[] genomes={1,1,1,4,2,3,3,1,4,4,4,4};
		final int[] starts={9,9,80,90,100,100,35,209,299,300,400,399};
		final int[] stops={20,20,84,99,109,109,39,220,320,321,419,420};
		String raw="# native fixture\n";
		String imported=BuildNcrnaReference.BATCH_HIT_HEADER+'\n';
		for(int i=0; i<genomes.length; i++){
			raw+=rawHit(genomes[i], starts[i], stops[i], i==2, "RF00050", "FMN");
			imported+=importHit(table, i+2, genomes[i], starts[i], stops[i], i==2, "RF00050", "FMN");
		}
		// The pre-existing competitor witness has native physical line101.
		for(int line=14; line<=100; line++){raw+='\n';}
		raw+=rawHit(4, 50, 59, false, "RF00010", "RNaseP");
		imported+=importHit(table, 101, 4, 50, 59, false, "RF00010", "RNaseP");
		raw+="# Program: cmsearch\n# Version: 1.1.5 fixture\n# Option settings: cmsearch --cut_ga --noali models.cm batch.fa\n# [ok]\n";
		String targets="batch_id\ttarget_id\tinput_id\tcontig_ordinal\toriginal_target\tlength\n";
		for(final int genome:new int[]{1,2,3,4,6}){
			targets+="0\ti"+genome+"c1\t"+genome+"\t1\tcontig"+genome+"\t1000\n";
		}
		final String complete="batch_id\tmodel\tFASTA_records\tinput_bases\tsearched_residues\treported_hits\tstatus\n"
			+"0\tRF00050\t5\t5000\t10000\t12\tCOMPLETE\n"
			+"0\tRF00010\t5\t5000\t10000\t1\tCOMPLETE\n"
			+"0\tRF00174\t5\t5000\t10000\t0\tCOMPLETE\n";
		write("models.tsv", "accession\tname\tCLEN\nRF00050\tFMN\t10\nRF00010\tRNaseP\t10\nRF00174\tcobalamin\t10\n");
		write("prepared/targets.tsv", targets);
		write("raw/hits.tbl", raw);
		write("raw/completion.tsv", complete);
		write("raw/exit.txt", "0\n");
		write("collection.tsv", "scan_group\tbatch_id\tprepared_metadata\tnative_result\tgenomes\tbases\trecords\tstate\n"
			+"fixture\t0\t"+path("prepared")+'\t'+path("raw")+"\t5\t5000\t5\tRAW_NATIVE_GA_COMPLETE\n");
		write("imported.hits.tsv.gz", imported);
		final String joinedEvents=events.replace("/saved/native.tbl", table);
		write("stage_evidence.fmn_events.tsv.gz", joinedEvents);
		write("stage_evidence.summary.tsv", summary);
		stageProcess("stages", 13);
		final ArrayList<String> counts=rows(path("stages.stages.tsv"));
		check(counts.contains("discovery\tbacteria\tCommon\t4\t2\t13\t12\t1\t0\t13\t12\t12\t12\t11\t11\t10\t3\t8\t4\t4\t4\t7\t2\t4\t0"),
			"Raw, imported, classified events, collapsed native cores and prior-only cores retain distinct literal denominators");
		check(counts.contains("discovery\tbacteria\tQuiet\t1\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"),
			"A mapped zero-hit input remains present despite no raw native rows");
		check(counts.contains("evaluation\tbacteria\tRare\t1\t0\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\t0\t0\t0\t1\t1\t1\t0\t0\t0\t0\t0\t1"),
			"Old-only evaluation priors are not relabelled as native-supported, and uninspected evaluation raw counts are NA");
		check(counts.contains("evaluation\tbacteria\tCommon\t1\t0\tNA\tNA\tNA\tNA\tNA\tNA\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t0\t1\t0"),
			"Actual evaluation native support is reported independently of retained old priors");
		final ArrayList<String> examples=rows(path("stages.examples.tsv.gz"));
		check(examples.size()==2 && examples.get(1).contains("\t1\tNOT_FULL_MODEL_OR_NATIVE_TRUNCATION\t"+table+"\t4\t"),
			"The native included partial hit drops at full-span/truncation and retains its physical raw-row witness");
		check(rows(path("stages.summary.tsv")).get(1).equals("1\t13\t13\t12\tSAVED_STAGE_COUNTS_NOT_BIOLOGICAL_ABSENCE_OR_SPLIT_AUTHORITY"),
			"Raw batch and imported populations conserve their independently specified rows");
		NcrnaConstructionMembers.main(new String[]{"mode=stages", "in="+path("result"), "groups="+path("groups.tsv"),
			"genomes="+path("genomes.tsv"), "evidence="+path("stage_evidence"), "imported="+path("imported"),
			"raw="+path("collection.tsv"), "models="+path("models.tsv"), "out="+path("stages_cli"),
			"expectedgenomes=7", "expectedmembers=10", "expectedeval=3", "expectednew=5", "expectedcores=17",
			"expectedbatches=1", "expectedraw=13", "t=1"});
		check(counts.equals(rows(path("stages_cli.stages.tsv"))), "Stage CLI and library must expose the same counts");
		stageFault("stage_model_missing", "raw/completion.tsv", complete, complete.replace("0\tRF00050\t5\t5000\t10000\t12\tCOMPLETE\n", ""));
		stageFault("stage_raw_count", "raw/completion.tsv", complete, complete.replace("\t12\tCOMPLETE", "\t11\tCOMPLETE"));
		stageFault("stage_raw_target", "raw/hits.tbl", raw, raw.replace("i1c1", "absent"));
		stageFault("stage_missing_footer", "raw/hits.tbl", raw, raw.replace("# [ok]\n", ""));
		stageFault("stage_failed_native", "raw/exit.txt", "0\n", "1\n");
		stageFault("stage_changed_import_flags", "imported.hits.tsv.gz", imported, imported.replace("\t1\t5\t10\t2\t", "\t1\t5\t10\t0\t"));
		stageFault("stage_changed_native_target", "imported.hits.tsv.gz", imported, imported.replace("\ti4c1\ttrue\t", "\ti4c2\ttrue\t"));
		stageFault("stage_changed_event_contig", "stage_evidence.fmn_events.tsv.gz", joinedEvents, joinedEvents.replace("\tcontig4\t", "\twrong_contig\t"));
		stageFault("stage_changed_event_strand", "stage_evidence.fmn_events.tsv.gz", joinedEvents, joinedEvents.replace("\t+\tRF00050\t", "\t-\tRF00050\t"));
		stageFault("stage_changed_map_contig", "prepared/targets.tsv", targets, targets.replace("\tcontig4\t", "\twrong_contig\t"));
		final String importA=importHit(table, 10, 4, 299, 320, false, "RF00050", "FMN");
		final String importB=importHit(table, 11, 4, 300, 321, false, "RF00050", "FMN");
		final String swappedImport=imported.replace(importA, importHit(table, 11, 4, 299, 320, false, "RF00050", "FMN"))
			.replace(importB, importHit(table, 10, 4, 300, 321, false, "RF00050", "FMN"));
		stageLocusFault("stage_import_line_swap", imported, swappedImport, joinedEvents, joinedEvents,
			"Classification must preserve the imported contig, coordinates and strand");
		final String eventA=event(9, 4, 299, 320, false, true).replace("/saved/native.tbl", table);
		final String eventB=event(10, 4, 300, 321, false, true).replace("/saved/native.tbl", table);
		final String swappedEvents=joinedEvents.replace(eventA, eventA.replace(table+"\t10\t", table+"\t11\t"))
			.replace(eventB, eventB.replace(table+"\t11\t", table+"\t10\t"));
		stageLocusFault("stage_raw_line_swap", imported, swappedImport, joinedEvents, swappedEvents,
			"Raw/import native target, mapped contig, coordinates and strand");
		// A raw-only FMN row must be localised before import, not silently disappear.
		final String reverseRaw=rawHit(6, 10, 19, false, "RF00050", "FMN").replace(" 10 19 + ", " 19 10 - ");
		check(reverseRaw.contains(" 19 10 - "), "Reverse-strand native example must use descending raw coordinates");
		write("raw/hits.tbl", raw.replace("# Program:", reverseRaw+"# Program:"));
		write("raw/completion.tsv", complete.replace("\t12\tCOMPLETE", "\t13\tCOMPLETE"));
		stageProcess("stage_raw_only", 14);
		check(rows(path("stage_raw_only.examples.tsv.gz")).get(2).contains("\tQuiet\t6\tABSENT_FROM_IMPORT\t"),
			"A real raw row missing from the imported ledger is retained as an explicit dropping stage example");
		write("raw/hits.tbl", raw);
		write("raw/completion.tsv", complete);
	}

	private static String rawHit(final int genome, final int start, final int stop, final boolean partial, final String model, final String family){
		return "i"+genome+"c1 - "+family+' '+model+" cm 1 "+(partial ? 5 : 10)+' '+start+' '+stop
			+" + "+(partial ? "3'" : "no")+" 1 0.5 0.0 100.0 1e-9 ! -\n";
	}
	private static String importHit(final String table, final int line, final int genome, final int start, final int stop,
			final boolean partial, final String model, final String family){
		return genome+"\t"+table+'\t'+line+"\tcontig"+genome+'\t'+start+'\t'+stop+"\t+\t"+model+"\t15.1\t"+family
			+"\tTARGET\t1\t"+(partial ? 5 : 10)+"\t10\t"+(partial ? 2 : 0)
			+"\t100.0\t1e-9\ttrue\tVALID_INTERVAL\tfalse\tfalse\ti"+genome+"c1\ttrue\t"+BuildNcrnaReference.BOUND_GA_BASIS+'\n';
	}
	private void stageProcess(final String output, final int rawCount){
		NcrnaConstructionMembers.stageAudit(path("result"), path("groups.tsv"), path("genomes.tsv"), path("stage_evidence"),
			path("imported"), path("collection.tsv"), path("models.tsv"), path(output), 7, 10, 3, 5, 17, 1, rawCount);
	}
	private void stageFault(final String name, final String file, final String valid, final String invalid){
		check(!valid.equals(invalid), "Stage fault must actually change the fixture: "+name);
		write(file, invalid);
		boolean rejected=false;
		try{stageProcess(name, 13);}
		catch(IllegalArgumentException expected){rejected=true;}
		finally{write(file, valid);}
		check(rejected, "Malformed stage evidence must fail loudly: "+name);
	}
	private void stageLocusFault(final String name, final String validImport, final String changedImport,
			final String validEvents, final String changedEvents, final String diagnostic){
		check(!validImport.equals(changedImport), "Swapped witnesses must really move between two distinct full-span loci");
		if(name.equals("stage_raw_line_swap")){check(!validEvents.equals(changedEvents), "Raw-only mismatch test must move both import and event witnesses together");}
		write("imported.hits.tsv.gz", changedImport);
		write("stage_evidence.fmn_events.tsv.gz", changedEvents);
		boolean rejected=false;
		try{stageProcess(name, 13);}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(diagnostic), "Wrong rejection would not prove locus binding: "+expected.getMessage());
			rejected=true;
		}finally{
			write("imported.hits.tsv.gz", validImport);
			write("stage_evidence.fmn_events.tsv.gz", validEvents);
		}
		check(rejected, "Same-genome, same-model-span native line swaps must fail at the physical locus identity check");
	}

	private void checkResult(final String prefix){
		final ArrayList<String> old=rows(path(prefix+".old.tsv.gz"));
		final ArrayList<String> events=rows(path(prefix+".events.tsv.gz"));
		final ArrayList<String> extracts=rows(path(prefix+".extract.tsv.gz"));
		check(old.size()==11 && events.size()==14 && extracts.size()==18, "Ten original members and thirteen events survive; shared-core provenance yields seventeen extraction cores");
		check(rows(path(prefix+".summary.tsv")).get(1).equals("10\t13\t17\t9\t10\t4\t3\t1\t3\t1\tREVIEW_AND_EXTRACTION_REQUIRED_NO_AUTOMATIC_COVERAGE_EXCEPTION"),
			"Admission, seed-only retention, heldout aliases, unresolved assembly and rare-lineage deficits have distinct conserved counts");
		final String[] corrected=find(extracts, "CM_discovery_1");
		check(corrected[3].equals(fasta(1)) && corrected[10].equals("9-20") && corrected[11].equals("12")
			&& corrected[16].equals("true") && corrected[18].equals("RF00050_a")
			&& corrected[19].equals("discovery:1;discovery:2"), "Corrected native bounds and duplicate witnesses preserve UTF-8 paths and original-member provenance");
		final String[] fragment=find(extracts, "OLD_RF00050_p");
		check(fragment[14].equals("false") && fragment[15].equals("true") && fragment[16].equals("false")
			&& fragment[17].equals("true"), "An old fragment survives as a seed prior without being declared intact training material");
		check(find(extracts, "CM_discovery_5")[16].equals("false") && find(extracts, "CM_discovery_6")[16].equals("false"),
			"Neither a positive evaluation-assembly alias nor an unresolved nonmatch may train");
		check(find(old, "true", "RF00050_c", 5)[32].equals("COMPETING_FULL_OTHER_FAMILY_REVIEW"),
			"Competing evidence stays explicitly reviewable rather than silently deleting the old record");
		check(find(old, "true", "RF00050_a", 5)[35].equals("CM_discovery_1"), "Two native witnesses replace one old member only once");
		final String[] multipart=find(extracts, "OLD_RF00050_u");
		check(multipart[10].equals("30-34,40-44") && multipart[14].equals("false")
			&& find(old, "true", "RF00050_u", 5)[35].equals("OLD_RF00050_u"),
			"A native event entirely inside an old multipart gap cannot replace its core or make the multipart prior intact");
		check(find(extracts, "OLD_RF00050_w")[10].equals("200-209")
			&& find(old, "true", "RF00050_w", 5)[32].equals("RETAINED_OLD_PRIOR_WEAK_NATIVE_ASSOCIATION"),
			"A one-base neighboring overlap is not a replacement of the old prior");
		check(find(extracts, "OLD_RF00050_m")[10].equals("300-319")
			&& find(old, "true", "RF00050_m", 5)[32].equals("RETAINED_OLD_PRIOR_MULTIPLE_NATIVE_CORES"),
			"Two distinct clean cores retain the original prior as an explicit review case even when each overlap is strong");
		final String[] partialPrior=find(old, "true", "RF00050_q", 5);
		final String[] sharedCore=find(extracts, "CM_discovery_11");
		check(partialPrior[32].equals("RETAINED_OLD_PRIOR_MULTIPLE_NATIVE_CORES")
			&& partialPrior[37].equals("false") && partialPrior[38].equals("true"),
			"A partial old prior keeps its own seed-only eligibility and ambiguous association state");
		check(sharedCore[10].equals("400-419") && sharedCore[14].equals("true") && sharedCore[16].equals("true")
			&& sharedCore[12].equals("CM_FULL_CLEAN;RETAINED_OLD_PRIOR_MULTIPLE_NATIVE_CORES")
			&& sharedCore[18].equals("RF00050_q") && sharedCore[19].contains("RF00050_q")
			&& sharedCore[19].contains("discovery:11"),
			"The one shared physical core retains both support bases; native complete support does not relabel its partial original source");
		check(rows(path(prefix+".coverage.tsv")).contains("bacteria\tRare\t1\t0\tREVIEW_NO_TRAINING_REPRESENTATIVE_NOT_AUTOMATIC_EXCEPTION"),
			"Rare-lineage loss is reported without silently admitting an evaluation assembly");
		for(int i=1; i<extracts.size(); i++){
			check(extracts.get(i).split("\t", -1).length==21, "Provenance values must not inject tabs into extraction columns");
			check(!extracts.get(i).startsWith("OLD_RF00050_c\t"), "Contradicted original prior is held for review, not extracted as a positive");
		}
	}

	private void fault(final String name, final String file, final String valid, final String invalid){
		check(!valid.equals(invalid), "Fault mutation must actually change the fixture: "+name);
		write(file, invalid);
		boolean rejected=false;
		try{process(name);}
		catch(IllegalArgumentException expected){rejected=true;}
		finally{write(file, valid);}
		check(rejected, "Malformed evidence must fail loudly: "+name);
	}
	private void process(final String output){
		NcrnaConstructionMembers.process(path("groups.tsv"), path("genomes.tsv"), path("original"), path("evidence"), path(output), 7, 10, 3, 5);
	}
	private void audit(final String output){
		NcrnaConstructionMembers.coverageAudit(path("result"), path("groups.tsv"), path("genomes.tsv"), path(output), 7, 10, 3, 5, 17);
	}
	private void auditFault(final String name, final String file, final String valid, final String invalid){
		check(!valid.equals(invalid), "Audit fault must actually change the fixture: "+name);
		write(file, invalid);
		boolean rejected=false;
		try{audit(name);}
		catch(IllegalArgumentException expected){rejected=true;}
		finally{write(file, valid);}
		check(rejected, "Malformed saved admission evidence must fail loudly: "+name);
	}
	private String path(final String file){return new File(root, file).getPath();}
	private void write(final String file, final String text){
		final File previous=new File(path(file));
		check(!previous.exists() || previous.delete(), "Synthetic input replacement must succeed");
		try(Output out=new Output(path(file))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private static ArrayList<String> rows(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){result.add(new String(row, StandardCharsets.UTF_8));}
		}
		return result;
	}
	private static String[] find(final ArrayList<String> rows, final String id){return find(rows, id, id, 0);}
	private static String[] find(final ArrayList<String> rows, final String first, final String id, final int field){
		for(final String row:rows){
			final String[] split=row.split("\t", -1);
			if(split[0].equals(first) && split[field].equals(id)){return split;}
		}
		throw new AssertionError("Missing literal fixture identity: "+id);
	}
	private static String fasta(final int index){return index==1 ? "/never-opened/α.fa" : "/never-opened/g"+index+".fa";}
	private static String phylum(final int index){return index==5 ? "Rare" : index==6 ? "Quiet" : "Common";}
	private static String cohort(final int index){return index==0 || index==5 ? "evaluation" : "discovery";}
	private static int input(final int index){return index==0 ? 1 : index==5 ? 2 : index;}
	private static String group(final int index, final int group, final boolean eval, final String disposition){
		return cohort(index)+'\t'+input(index)+'\t'+fasta(index)+'\t'+fasta(index)+"\tGCF_"+index+'\t'+group+'\t'+eval
			+"\tfalse\tfalse\tORIGINAL_DIAGNOSTIC\t"+group+'\t'+eval+'\t'+(index==3)+'\t'+disposition+'\n';
	}
	private static String old(final String id, final int genome, final int start, final int stop, final boolean partial){
		return "true\tfalse\tORIGINAL_PATH\tUNVERIFIED\tNOT_EVALUATED\tRF00050_"+id+"\tg"+genome+"\tGCF_"+genome+".1\t"+(100+genome)
			+"\tbacteria\t"+phylum(genome)+"\tcontig"+genome+"\t+\t"+start+'-'+stop+'\t'+(stop-start+1)
			+"\t1\tfeature_"+id+'\t'+partial+"\tfalse\t"+(partial ? "SOURCE_PARTIAL" : "ANNOTATION_FULL_COPY")
			+"\t1000\tfalse\t0\t0\t0\t1\told_record\t"+(100+genome)+"\tunknown\torphan\t\torphan\n";
	}
	private static String event(final int id, final int genome, final int start, final int stop, final boolean partial, final boolean fmn){
		return id+"\t"+input(genome)+"\tcontig"+genome+'\t'+start+'\t'+stop+"\t+\t"+(fmn ? "RF00050\tFMN" : "RF00010\tRNaseP")
			+"\t1\t"+(partial ? 5 : 10)+"\t10\t"+(partial ? 2 : 0)+"\tfalse\tfalse\t100.0\t66.0\t151.5\ttrue\tfalse\t90\ttrue\tfalse\t"
			+(partial ? "MARGINAL\tpartial_model_support" : "CONFIDENT\tnative_GA_full_model_clean_family")
			+"\tNATIVE_FLAGS\t0\t/saved/native.tbl\t"+(id+1)+"\tPRINTED_0.1_BIT_NOT_INTERNAL_PRECISION\t6_DECIMAL_PLACES_SOURCE_PAIR_RETAINED\n";
	}
	private static String link(final String old, final int genome, final int length, final int overlap, final String event){
		return cohort(genome)+"\tRF00050_"+old+'\t'+length+'\t'+overlap+"\tSOURCE_LINE\t"+event;
	}
	private static String state(final String old, final int genome, final int length, final int links, final int own, final int other){
		return cohort(genome)+"\tRF00050_"+old+'\t'+input(genome)+"\t1\t1\t"+links+'\t'+own+'\t'+other+'\t'+(other>0 ? length : 0)
			+'\t'+length+"\tDIAGNOSTIC\tNOT_EVALUATED\n";
	}
	private static void check(final boolean value, final String why){if(!value){throw new AssertionError(why);}}
	private final String root;
}
