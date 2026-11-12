package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import prok.BuildNcrnaReference;
import prok.DraftNcrnaScenarios;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceEvidence;
import prok.WriteNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Independent finite geometry/precedence expectations over literal saved tables.
 * Genomes, native reports and source GFF payloads intentionally do not exist.
 * @author Brian Bushnell, Raiden
 */
final class DraftNcrnaScenariosTest{
	DraftNcrnaScenariosTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Scenario-construction fixture must be new");
		prepare();
		invoke("result");
		final String q='\t'+DraftNcrnaScenarios.QUALIFICATION;
		final String[] summary=lines("result/summary.tsv");
		check(summary.length==3 && summary[1].equals("2\t2\t16\t15\t1\t3\t4\t2\t6\t3\tGA\t6\t4\t7\t1\tEVIDENCE_COMPLETE_NOT_ADJUDICATED"+q)
			&& summary[2].equals("2\t2\t16\t15\t1\t3\t4\t2\t6\t3\tGA_PLUS_MARGIN\t4\t6\t7\t1\tEVIDENCE_COMPLETE_NOT_ADJUDICATED"+q),
			"Sixteen events yield fifteen model intervals and fourteen model decisions; four prior rows stay fixed");
		final String[] low=lines("result/ga.tsv.gz"), high=lines("result/margin.tsv.gz");
		check(low.length==19 && high.length==19, "Both scenarios retain the same eighteen declared decisions");
		check(find(low, "event_1\t1\t").contains("\tRESOLVED\tscored\t")
			&& find(high, "event_1\t1\t").contains("\tRESOLVED\tdont_care\t"), "Only the numeric margin changes this clear full candidate");
		check(find(low, "event_1\t2\t").contains("\tRESOLVED\tscored\t"), "Equal native line and coordinates in another input stay independent");
		check(find(low, "event_2\t1\t").contains("\tRESOLVED\tscored\t"), "A false partial claim and a context pseudo claim cannot block the RNA candidate");
		check(find(low, "event_3\t1\t").contains("\tscored\t") && absent(low, "event_4\t1\t"), "Identical same-model events collapse without creating a competing locus");
		check(find(low, "event_5\t1\t").contains("\tUNRESOLVED\tunresolved\t")
			&& find(low, "event_6\t1\t").contains("\tUNRESOLVED\tunresolved\t"), "Different-family overlap creates no masking exemption on either side");
		check(find(low, "event_7\t1\t").contains("\t130\t139\t+\tFMN\tRESOLVED\tdont_care\t")
			&& find(low, "event_8\t1\t").contains("\t135\t144\t+\tFMN\tRESOLVED\tdont_care\t"), "Same-family competition retains two original intervals, never their union");
		check(find(low, "event_9\t1\t").contains("overlaps_preserved_prior"), "A new candidate cannot widen or duplicate a protected prior");
		check(find(low, "event_10\t1\t").contains("unresolved_cross_family"), "Direct other-family RNA source claim blocks automatic score and mask");
		check(find(low, "event_11\t1\t").contains("\tRNaseP\tRESOLVED\tdont_care\tdiagnostic_declared_policy\t.\t")
			&& absent(low, "event_12\t1\t"), "Identical physical masks collapse across models while model counts stay separate");
		check(find(low, "event_13\t1\t").contains("no_automatic_addition")
			&& find(low, "event_14\t1\t").contains("model_context_not_mature_tRNA_truth"), "No-add and model-only context stay fixed outside scored additions");
		check(find(low, "event_15\t1\t").contains("fixed_pending_structure") && find(high, "event_15\t1\t").contains("fixed_pending_structure"),
			"An isolated fixed-mask model remains masked even with full high-scoring support");
		check(find(low, "event_15\t1\t").equals(find(high, "event_15\t1\t")), "Fixed dispositions retain the same policy and witness text on both sides");
		check(find(low, "event_1\t1\t").endsWith("#1/event_1") && find(low, "event_1\t2\t").endsWith("#2/event_1"),
			"Witness navigation must include the input, since native event numbers repeat between genomes");
		for(final String id:new String[]{"prior_2\t", "prior_3\t", "prior_4\t", "prior_5\t"}){
			check(find(low, id).equals(find(high, id)), "Prior import is identical between scenarios");
		}
		check(find(low, "prior_3\t").contains("\t320\t329\t+\tFMN\tRESOLVED\tscored\t"), "Accepted prior without any current hit remains scored");
		check(find(low, "prior_4\t").contains("\tUNRESOLVED\tunresolved\t"), "Unreviewed source row gets no inherited mask and does not freeze a new hit at its position");
		final String[] witnesses=lines("result/witnesses.tsv.gz");
		check(witnesses.length==17 && find(witnesses, "event_3\t1\t/missing/hits1.tbl\t4\t").contains("\t3\t"), "Every duplicate original event keeps its own native-line witness");
		check(find(witnesses, "event_11\t1\t/missing/hits1.tbl\t12\t").contains("RF00011"), "Collapsed masks preserve the second model's original evidence");
		check(find(lines("result/counts.tsv"), "RF00050\tFMN\tGA\t").equals("RF00050\tFMN\tGA\t9\t4\t2\t3\t0"+q), "Model counts use unique model intervals, not collapsed physical masks or raw duplicates");
		check(find(lines("result/family_counts.tsv"), "FMN\tGA\t").equals("FMN\tGA\t6\t2\t4\t0"+q)
			&& find(lines("result/family_counts.tsv"), "RNaseP\tGA\t").equals("RNaseP\tGA\t0\t2\t0\t0"+q), "Family rows include preserved priors and count coincident masks once");
		check(lines("result/prior_witnesses.tsv.gz").length==5 && find(lines("result/prior_witnesses.tsv.gz"), "prior_5\t").contains("\tB12clade\t190\t199\t192\t197\t"),
			"Clan-only original core and classification survive without becoming scored or masked");
		check(!new File(path("result/reference.tsv")).exists() && !new File(path("result/input_1.gff.gz")).exists(), "Draft tables do not publish biological reference GFFs/manifests");
		final String warning=write("warning.links.tsv", links.replace("SOURCE_PARTIAL_FALSE", "SOURCE_PARTIAL_TRUE"));
		invoke("warning", "links="+warning);
		check(find(lines("warning/ga.tsv.gz"), "event_2\t1\t").contains("direct_family_claim_or_RNA_attribute_warning"), "A true direct-RNA partial claim blocks automatic scoring");
		final String regulatory=write("regulatory.links.tsv", links.replace(contextLink,
			link(3, 2, "regulatory_region", "OTHER_ROW", "SOURCE_OTHER_KNOWN_FAMILY", "SOURCE_PSEUDO_TRUE")));
		invoke("regulatory", "links="+regulatory);
		check(find(lines("regulatory/ga.tsv.gz"), "event_2\t1\t").contains("unresolved_cross_family"),
			"Direct regulatory-row RFAM conflict stays visible without inheriting that row's pseudo attribute");
		final String duplicateLinks=links.replace(contextLink, crossLink);
		bad("duplicate_link", "Repeated input/GFF", "links="+write("duplicate.links.tsv", duplicateLinks));
		bad("unknown_flag", "saved source flag", "links="+write("unknown.links.tsv", links.replace("SOURCE_PARTIAL_FALSE", "FUTURE_WARNING")));
		bad("repeated_flag", "saved source flag", "links="+write("repeated.links.tsv", links.replace("SOURCE_PARTIAL_FALSE", "SOURCE_PARTIAL_FALSE;SOURCE_PARTIAL_FALSE")));
		bad("policy_identity", "captured metadata", "policy="+write("bad.policy.tsv", policies.replace("00000000000000000001", "00000000000000000002")));
		bad("prior_bounds", "contig observation", "priors="+write("bad.prior.tsv", priorText.replace("\tseq\t400\t", "\tseq\t401\t")));
		variant("conflict", hits.replace(hit(1, 4, "RF00050", "FMN", 70, 79, 60), hit(1, 4, "RF00050", "FMN", 70, 79, 61)), null);
		invoke("conflict_result", "in="+path("conflict"));
		check(find(lines("conflict_result/ga.tsv.gz"), "event_3\t1\t").contains("conflicting_duplicate_observations"), "Conflicting duplicate scores must not choose a convenient winner");
		final String first=hit(1, 1, "RF00050", "FMN", 10, 19, 40);
		variant("edge", hits.replace(first, hit(1, 1, "RF00050", "FMN", 1, 10, 60).replace("\tfalse\tfalse\n", "\ttrue\tfalse\n")), null);
		invoke("edge_result", "in="+path("edge"));
		check(find(lines("edge_result/ga.tsv.gz"), "event_1\t1\t").contains("physical_contig_edge"), "Physical clipping stays masked despite a full high score");
		final String partial=first.replace("\tTARGET_CLASS\t1\t10\t10\t0\t", "\tTARGET_CLASS\t2\t9\t10\t0\t");
		variant("partial", hits.replace(first, partial), "RF00050\tfixture\t10\t10\t0\t9\t10\t2\nRF00050\tfixture\t10\t8\t0\t1\t10\t2\n");
		invoke("partial_result", "in="+path("partial"));
		check(find(lines("partial_result/ga.tsv.gz"), "event_1\t1\t").contains("partial_model_support"), "Partial model span cannot acquire full-gene scoring status");
		variant("same_below_ga", hits.replace(hit(1, 8, "RF00050", "FMN", 135, 144, 60), hit(1, 8, "RF00050", "FMN", 135, 144, 19).replace("\ttrue\tVALID_INTERVAL", "\tfalse\tVALID_INTERVAL")), null);
		bad("same_below_result", "requires GA-passing native events", "in="+path("same_below_ga"));
		variant("other_below_ga", hits.replace(hit(1, 6, "RF00162", "SAM_I", 105, 114, 70), hit(1, 6, "RF00162", "SAM_I", 105, 114, 19).replace("\ttrue\tVALID_INTERVAL", "\tfalse\tVALID_INTERVAL")), null);
		bad("other_below_result", "requires GA-passing native events", "in="+path("other_below_ga"));
		variant("bad_ga_flag", hits.replace(first, first.replace("\ttrue\tVALID_INTERVAL", "\tfalse\tVALID_INTERVAL")), null);
		bad("bad_ga_flag_result", "Saved GA flag disagrees", "in="+path("bad_ga_flag"));
		variant("exact_ga", hits.replace(first, hit(1, 1, "RF00050", "FMN", 10, 19, 20)), null);
		invoke("exact_ga_result", "in="+path("exact_ga"));
		check(find(lines("exact_ga_result/ga.tsv.gz"), "event_1\t1\t").contains("\tRESOLVED\tscored\t"), "The GA boundary is inclusive, not strictly greater-than");
		final String unicodePolicy=write("policy%Ω.tsv", policies);
		final String unicodeLinks=write("links%Ω.tsv", links);
		invoke("unicode%Ω", "policy="+unicodePolicy, "links="+unicodeLinks);
		final String encodedPolicy=utf8(BuildNcrnaReference.field(new ByteBuilder(), unicodePolicy));
		final String encodedLinks=utf8(BuildNcrnaReference.field(new ByteBuilder(), unicodeLinks));
		final String encodedWitness=utf8(BuildNcrnaReference.field(new ByteBuilder(), path("unicode%Ω/witnesses.tsv.gz")));
		check(find(lines("unicode%Ω/ga.tsv.gz"), "event_1\t1\t").endsWith(encodedPolicy+":2\t"+encodedWitness+"#1/event_1"),
			"Cached provenance bytes preserve percent escaping and UTF-8 without double encoding");
		check(find(lines("unicode%Ω/witnesses.tsv.gz"), "event_1\t1\t").endsWith(encodedLinks+q),
			"Every native witness preserves the original Unicode companion path");
		System.err.println("DRAFT_NCRNA_SCENARIOS_TEST_PASS");
	}

	private void prepare(){
		write("census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t/missing/one.fa\t410\t400\t1\n2\t/missing/two.fa\t410\t400\t1\n");
		write("members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+"\nfixture\t1\t1\t/missing/one.fa\tNA\tNA\tNA\tUNKNOWN\t/missing/source.gff\tNA\nfixture\t2\t2\t/missing/two.fa\tNA\tNA\tNA\tUNKNOWN\tNA\tNA\n");
		final String[] accessions={"RF00050", "RF00162", "RF00010", "RF00011", "RF01851", "RF00005", "RF00174"};
		final String[] families={"FMN", "SAM_I", "RNaseP", "RNaseP", "tmRNA", "tRNA", "cobalamin"};
		final ByteBuilder scope=new ByteBuilder(NcrnaReferenceEvidence.SCOPE_HEADER+"\n");
		final ByteBuilder assets=new ByteBuilder(NcrnaReferenceEvidence.ASSET_HEADER+"\n");
		final ByteBuilder queries=new ByteBuilder(NcrnaReferenceEvidence.QUERY_HEADER+"\n");
		final ByteBuilder policy=new ByteBuilder(DraftNcrnaScenarios.POLICY_HEADER+"\n");
		for(int i=0; i<accessions.length; i++){
			scope.append(accessions[i]).tab().append(families[i]).append("\tprimary\tfixture\t").append(families[i]).append("\tTARGET_CLASS\n");
			assets.append(accessions[i]).append("\tfixture\t10\t20\t00000000000000000001\t/missing/model.cm\n");
			if(i==6){continue;}
			queries.append(accessions[i]).append("\tmodel_").append(i).append("\t10\n");
			policy.append(accessions[i]).tab().append(families[i]).append("\tfixture\t00000000000000000001\t20\t10\t0\t0\tfixture\t1\t1\t1\t30\t30\t30\t1\t1\t")
				.append((i<3) ? "30\tBOUNDARY_UNSELECTED" : i==3 ? "NA\tMASK_PENDING_STRUCTURE_REVIEW" : i==4 ? "NA\tNO_AUTOMATIC_ADDITION" : "NA\tCONTEXT_PENDING_BRIAN")
				.append("\tfixture\tDRAFT_NOT_APPLIED\tfixture\t1\n");
		}
		write("scope.tsv", utf8(scope));
		write("assets.tsv", utf8(assets));
		write("queries.tsv", utf8(queries));
		policies=utf8(policy);
		write("policy.tsv", policies);
		write("evidence.summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+"\n2\t2\t800\t2\t16\t2\t1\t3\t3\t0\t3\t4\t6\t12\t5\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n");
		final ByteBuilder coverage=new ByteBuilder(BuildNcrnaReference.COVERAGE_HEADER+"\n");
		for(int i=0; i<6; i++){coverage.append(accessions[i]).append("\tfixture\t10\t10\t0\t").append(i==0 ? 10 : i==3 ? 2 : 1).tab().append(i==0 ? 10 : i==3 ? 2 : 1).append("\t2\n");}
		write("evidence.coverage.tsv", utf8(coverage));
		write("evidence.sources.tsv.gz", BuildNcrnaReference.SOURCE_HEADER+"\n1\t/missing/source.gff\tPARSED_GFF\t1\t3\t3\t0\t0\t0\t0\t0\t0\n2\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n");
		hits=hit(1, 1, "RF00050", "FMN", 10, 19, 40)+hit(1, 2, "RF00050", "FMN", 40, 49, 60)
			+hit(1, 3, "RF00050", "FMN", 70, 79, 60)+hit(1, 4, "RF00050", "FMN", 70, 79, 60)
			+hit(1, 5, "RF00050", "FMN", 100, 109, 60)+hit(1, 6, "RF00162", "SAM_I", 105, 114, 70)
			+hit(1, 7, "RF00050", "FMN", 130, 139, 60)+hit(1, 8, "RF00050", "FMN", 135, 144, 60)
			+hit(1, 9, "RF00050", "FMN", 170, 179, 60)+hit(1, 10, "RF00050", "FMN", 200, 209, 60)
			+hit(1, 11, "RF00010", "RNaseP", 230, 239, 60)+hit(1, 12, "RF00011", "RNaseP", 230, 239, 60)
			+hit(1, 13, "RF01851", "tmRNA", 260, 269, 60)+hit(1, 14, "RF00005", "tRNA", 290, 299, 60)
			+hit(1, 15, "RF00011", "RNaseP", 350, 359, 60)+hit(2, 1, "RF00050", "FMN", 10, 19, 40);
		write("evidence.hits.tsv.gz", BuildNcrnaReference.HIT_HEADER+'\n'+hits);
		write("evidence.overlaps.tsv.gz", BuildNcrnaReference.OVERLAP_HEADER+'\n'+pair(3, 4, 10)+pair(5, 6, 5)+pair(7, 8, 5)+pair(11, 12, 10)
			+sourcePair(1, 10)+sourcePair(2, 2)+sourcePair(3, 2));
		write("evidence.contigs.tsv.gz", WriteNcrnaReference.CONTIG_HEADER+"\n1\t/missing/one.fa\tseq\t400\n2\t/missing/two.fa\tseq\t400\n");
		crossLink=link(1, 10, "ncRNA", "RNA_ROW", "SOURCE_OTHER_KNOWN_FAMILY", "NONE");
		contextLink=link(3, 2, "gene", "GENE_OR_TRANSCRIPT_CONTEXT", "SOURCE_SAME_ACCESSION", "SOURCE_PSEUDO_TRUE");
		links=DraftNcrnaScenarios.LINKS_HEADER+'\n'+crossLink+link(2, 2, "ncRNA", "RNA_ROW", "SOURCE_SAME_ACCESSION", "SOURCE_PARTIAL_FALSE")+contextLink;
		write("links.tsv", links);
		priorText=PRIOR_HEADER+'\n'+prior(1, "RF00050", "ACCEPTED_NOVEL_FMN", 170, 179, 0, 0)
			+prior(2, "RF00050", "ACCEPTED_NOVEL_FMN", 320, 329, 0, 0)
			+prior(3, "RF00050", "UNRESOLVED", 10, 19, 0, 0)+prior(4, "RF00174", "B12clade", 190, 199, 192, 197);
		write("priors.tsv", priorText);
	}
	private static String hit(final int input, final int line, final String accession, final String family, final int first, final int last, final int score){
		return input+"\t/missing/hits"+input+".tbl\t"+line+"\tseq\t"+first+'\t'+last+"\t+\t"+accession+"\tfixture\t"+family+"\tTARGET_CLASS\t1\t10\t10\t0\t"+score+"\t1e-8\ttrue\tVALID_INTERVAL\tfalse\tfalse\n";
	}
	private static String pair(final int a, final int b, final int overlap){return "1\tMODEL_MODEL\t/missing/hits1.tbl\t"+a+"\t/missing/hits1.tbl\t"+b+'\t'+overlap+"\tSAME_STRAND_GEOMETRY\tfalse\n";}
	private static String sourcePair(final int line, final int hit){return "1\tHISTORICAL_MODEL\t/missing/source.gff\t"+line+"\t/missing/hits1.tbl\t"+hit+"\t10\tSAME_STRAND_CANDIDATE\tfalse\n";}
	private static String link(final int line, final int hit, final String type, final String kind, final String comparison, final String flags){return "1\t/missing/source.gff\t"+line+"\t/missing/hits1.tbl\t"+hit+'\t'+type+'\t'+kind+"\tSAME_STRAND_OVERLAP\t"+comparison+'\t'+flags+'\n';}
	private static String prior(final int line, final String accession, final String assignment, final int first, final int last, final int coreFirst, final int coreLast){return "fixture\t/missing/prior.tsv\t"+line+"\tlocus\tgenome\t/missing/one.fa\t1\t/missing/one.fa\tseq\t400\t"+accession+'\t'+assignment+'\t'+first+'\t'+last+'\t'+coreFirst+'\t'+coreLast+"\t+\t"+(accession.equals("RF00050") ? "PRIOR_REVIEWED_WHOLE_INTERVAL_NOT_NEW_APTAMER_BOUNDS" : "SOURCE_AND_REVIEWED_CORE")+"\tfixture_authority\tfixture_row\toriginal_row\n";}
	/** Tiny saved-table variants are data fixtures; no program/runtime copy exists. */
	private void variant(final String name, final String changedHits, final String fmnCoverage){
		for(final String suffix:new String[]{".summary.tsv", ".sources.tsv.gz", ".overlaps.tsv.gz", ".contigs.tsv.gz"}){
			try(Input in=new Input(path("evidence"+suffix)); Output out=new Output(path(name+suffix))){
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){out.writer.println(row);}
			}
		}
		final ByteBuilder coverage=new ByteBuilder();
		for(final String row:lines("evidence.coverage.tsv")){
			if(fmnCoverage!=null && row.startsWith("RF00050\t")){coverage.append(fmnCoverage);}
			else{coverage.append(row).nl();}
		}
		write(name+".coverage.tsv", utf8(coverage));
		write(name+".hits.tsv.gz", BuildNcrnaReference.HIT_HEADER+'\n'+changedHits);
	}
	private void invoke(final String name, final String... extra){
		final String[] base={"in="+path("evidence"), "census="+path("census.tsv"), "memberships="+path("members.tsv"), "scope="+path("scope.tsv"),
			"assets="+path("assets.tsv"), "queries="+path("queries.tsv"), "policy="+path("policy.tsv"), "priors="+path("priors.tsv"), "links="+path("links.tsv"),
			"out="+path(name), "t=1", "pigz=f", "bgzip=f", "zipthreads=1"};
		final String[] args=Arrays.copyOf(base, base.length+extra.length);
		System.arraycopy(extra, 0, args, base.length, extra.length);
		DraftNcrnaScenarios.main(args);
	}
	private void bad(final String name, final String expected, final String... args){
		try{invoke(name, args);}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(expected), "Unexpected construction failure: "+failure);
			check(!new File(path(name+"/summary.tsv")).exists(), "Invalid construction cannot publish completion");
			return;
		}
		throw new AssertionError("Expected construction failure: "+expected);
	}
	private String path(final String file){return new File(root, file).getAbsolutePath();}
	private String write(final String name, final String text){
		final String file=path(name);
		try(Output out=new Output(file)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		return file;
	}
	private String[] lines(final String name){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				result.add(new String(row, StandardCharsets.UTF_8));
			}
		}
		return result.toArray(new String[result.size()]);
	}
	private static String utf8(final ByteBuilder b){return new String(b.array, 0, b.length(), StandardCharsets.UTF_8);}
	private static String find(final String[] rows, final String prefix){
		for(final String row:rows){if(row.startsWith(prefix)){return row;}}
		throw new AssertionError("Missing expected row: "+prefix);
	}
	private static boolean absent(final String[] rows, final String prefix){
		for(final String row:rows){if(row.startsWith(prefix)){return false;}}
		return true;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private String policies, links, priorText, contextLink, crossLink, hits;
	private static final String PRIOR_HEADER="source_id\tdecision_file\tdecision_line\toriginal_locus\toriginal_genome\toriginal_fasta\tinput_id\tresolved_fasta\tliteral_target\tcontig_length\tsource_family\tassignment\tsource_start1\tsource_stop1\tcore_start1\tcore_stop1\tstrand\tinterval_semantics\tsource_authority\trow_authority\toriginal_row";
}
