package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDataReconciliation;
import prok.NcrnaDataReconciliation.Assembly;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Metadata fixture: subtype witnesses, protected priors, Unicode and unmapped
 * renamed/reversed aliases. Empty source files prove no genome reader is needed.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDataReconciliationTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaDataReconciliationTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, split[0], split.length>1 ? split[1] : null), "Unknown fixture argument: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Supply out=fresh_native_fixture_directory");
			fixture(new File(parser.out1).toPath());
			System.err.println("NCRNA_DATA_RECONCILIATION_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void fixture(final Path root) throws Exception{
		Files.createDirectory(root);
		final Path sources=Files.createDirectory(root.resolve("söurces"));
		final String left=Files.createFile(sources.resolve("left.fna.gz")).toString();
		final String right=Files.createFile(sources.resolve("right.fna.gz")).toString();
		final String third=Files.createFile(sources.resolve("newly_scanned.fna.gz")).toString();
		Files.createFile(sources.resolve("unscanned.fna.gz"));
		write(root, "groups.tsv", GROUP_HEADER
			+"evaluation\t1\t"+left+"\t"+left+"\tGCF_1\t0\ttrue\tfalse\tfalse\tEVALUATION\t0\ttrue\tfalse\tEVALUATION\n"
			+"discovery\t1\t"+right+"\t"+right+"\tGCF_1\t1\tfalse\tfalse\tfalse\tALIAS\t0\ttrue\tfalse\tALIAS\n");
		write(root, "mapping.tsv", "family\tmember\treference_model\ncobalamin\tRF00174\tRF00174\ncobalamin\tRF01482\tRF01482\ncobalamin\tRF01689\tRF01689\n");
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ArrayList<Row> oldRows=new ArrayList<Row>();
		oldRows.add(reference(codec, root, "old.witness", 100, 199, false, "RF00174"));
		codec.write(root.resolve("old.gff").toString(), "1", oldRows);
		final String exact=oldRows.get(0).id;
		final ArrayList<Row> newRows=new ArrayList<Row>();
		newRows.add(reference(codec, root, "new.witness", 100, 199, false, "RF01482"));
		newRows.add(reference(codec, root, "new.witness", 400, 499, true, "RF00174"));
		newRows.add(reference(codec, root, "new.witness", 700, 799, false, "RF00174"));
		newRows.add(reference(codec, root, "new.witness", 900, 999, false, "RF00174"));
		codec.write(root.resolve("new.gff").toString(), "1", newRows);
		codec.write(root.resolve("third.gff").toString(), "1", new ArrayList<Row>());
		final String weak=newRows.get(1).id, prior=newRows.get(2).id, priorOnly=newRows.get(3).id;
		write(root, "old.census", CENSUS_HEADER+"1\t"+left+"\t0\t1000\t1\n");
		write(root, "new.census", CENSUS_HEADER+"1\t"+right+"\t0\t1000\t1\n");
		write(root, "third.census", CENSUS_HEADER+"1\t"+third+"\t0\t1000\t1\n");
		write(root, "old.reference", "input_id\tfasta\tgff\n1\t"+left+"\t"+root.resolve("old.gff")+"\n");
		write(root, "new.reference", "input_id\tfasta\tgff\n1\t"+right+"\t"+root.resolve("new.gff")+"\n");
		write(root, "third.reference", "input_id\tfasta\tgff\n1\t"+third+"\t"+root.resolve("third.gff")+"\n");
		write(root, "old.events", EVENT_HEADER+event(1, 100, 199, "RF00174", "CONFIDENT", true, 0)
			+event(2, 250, 299, "RF00174", "MARGINAL", false, 0));
		write(root, "new.events", EVENT_HEADER+event(1, 100, 199, "RF00174", "CONFIDENT", true, 0)
			+event(2, 400, 499, "RF00174", "BELOW_MARGINAL", false, 0)
			+event(3, 400, 499, "RF00174", "CONFIDENT", true, 16)
			+event(4, 100, 199, "RF01482", "CONFIDENT", true, 0)
			+event(5, 400, 499, "RF00174", "MARGINAL", true, 0)
			+event(6, 700, 799, "RF00174", "COVERED_BY_PROTECTED_PRIOR", true, 0));
		write(root, "third.events", EVENT_HEADER+event(1, 100, 199, "RF00174", "CONFIDENT", true, 0)
			+event(2, 190, 250, "RF00174", "CONFIDENT", true, 0)
			+event(3, 240, 299, "RF00174", "CONFIDENT", true, 0)
			+event(4, 300, 320, "RF00174", "CONFIDENT", true, 0)
			+event(5, 100, 199, "RF00174", "CONFIDENT", true, 0).replace("\t+\tRF00174", "\t-\tRF00174"));
		write(root, "third.witness", WITNESS_HEADER);
		final String validOld=WITNESS_HEADER+witness(exact, 1, 100, 199, "RF00174", -1, "QUALIFYING_SAME_FAMILY_GROUP_MEMBER");
		write(root, "old.witness", validOld);
		write(root, "new.witness", WITNESS_HEADER+witness(exact, 1, 100, 199, "RF00174", -1, "QUALIFYING_SAME_FAMILY_GROUP_MEMBER")
			+witness(exact, 4, 100, 199, "RF01482", -1, "QUALIFYING_SAME_FAMILY_GROUP_MEMBER")
			+witness(weak, 2, 400, 499, "RF00174", -1, "WEAK_OVERLAP_WITNESS_NOT_GROUPING_SUPPORT")
			+witness(weak, 3, 400, 499, "RF00174", -1, "TANDEM_BRIDGE_WITNESS_NOT_REPRESENTATIVE")
			+witness(weak, 5, 400, 499, "RF00174", -1, "QUALIFYING_SAME_FAMILY_GROUP_MEMBER")
			+prior(prior, 700, 799, 1)+witness(prior, 6, 700, 799, "RF00174", 1, "SAME_FAMILY_OVERLAP_NOT_RESCORING_PRIOR_BOUNDS")
			+prior(priorOnly, 900, 999, 2));
		write(root, "cohorts.tsv", "cohort\tcensus\treference\tevents\twitnesses\texpected_inputs\texpected_scored_RF00174\n"
			+"original\t"+root.resolve("old.census")+"\t"+root.resolve("old.reference")+"\t"+root.resolve("old.events")+"\t"+root.resolve("old.witness")+"\t1\t1\n"
			+"import\t"+root.resolve("new.census")+"\t"+root.resolve("new.reference")+"\t"+root.resolve("new.events")+"\t"+root.resolve("new.witness")+"\t1\t3\n"
			+"completion\t"+root.resolve("third.census")+"\t"+root.resolve("third.reference")+"\t"+root.resolve("third.events")+"\t"+root.resolve("third.witness")+"\t1\t0\n");
		write(root, "genomes.tsv", GENOME_HEADER
			+"0\tdonör\tGCF_1.1\t1\tbacteria\tp\tfixture\t"+left+"\t1\t1000\t1\n"
			+"1\tdonor_alias\tGCF_1.1\t1\tbacteria\tp\tfixture\t"+right+"\t1\t1000\t1\n");
		final String header="RF00174_1_0 source=donör", renamed="RF00174_2_0 source=donor_alias";
		write(root, "assignments.tsv", "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n"
			+header+" member_family=RF00174 build_group=RF00174\t1\tRF00174\tRF00174_0\tcobalamin_RF00174_0\tkept\n"
			+renamed+" member_family=RF00174 build_group=RF00174\t1\tRF00174\tRF00174_0\tcobalamin_RF00174_0\tkept\n"
			+"orphan\t1\tRF00174\tRF00174_1\t\torphan\n");
		write(root, "ledger.tsv", "locus_id\tsource_family\tmember_family\tbuild_group\trouting_rule\trouting_bits\tlength\tinput_header\n"
			+"RF00174_1_0\tRF00174\tRF00174\tRF00174\tfixture\t0\t100\t"+header+"\n"
			+"RF00174_2_0\tRF00174\tRF00174\tRF00174\tfixture\t0\t100\t"+renamed+"\n");
		write(root, "proxy.tsv", prok.NcrnaDdlLocusOverlap.RAW_HEADER+"\n"
			+"RF00174_1_0\tdonör\tGCF_1.1\t1\tbacteria\tp\tNC_000001.1\t+\t100-199\t100\t1\tid\tfalse\tfalse\tfixture\t1000\tfalse\t0\t0\t0\t1\n"
			+"RF00174_2_0\tdonor_alias\tGCF_1.1\t1\tbacteria\tp\trenamed\t-\t800-899\t100\t1\tid\tfalse\tfalse\tfixture\t1000\tfalse\t0\t0\t0\t1\n");
		write(root, "recovered.tsv", prok.NcrnaDdlLocusOverlap.RECOVERED_HEADER+"\n");
		write(root, "overlay.tsv", prok.NcrnaDdlLocusOverlap.OVERLAY_HEADER+"\n");
		final String[] parameters={"in="+root.resolve("cohorts.tsv"), "out="+root.resolve("result"),
			"mapping="+root.resolve("mapping.tsv"), "groups="+root.resolve("groups.tsv"), "roots="+sources,
			"genomes="+root.resolve("genomes.tsv"), "assignments="+root.resolve("assignments.tsv"),
			"ledger="+root.resolve("ledger.tsv"), "proxy="+root.resolve("proxy.tsv"), "recovered="+root.resolve("recovered.tsv"),
			"overlays="+root.resolve("overlay.tsv"), "expectedgroups=1", "expectedreference=1", "expectedconfident=2",
			"expectedkept=2", "expectedassigned=3", "expectedgenomes=2", "expectedmodels=1", "expectedcohorts=3", "t=1", "pigz=f", "unpigz=f"};
		NcrnaDataReconciliation.main(parameters);
		check(value(root, "qualified_exact_coordinate_groups")==9, "All three cohort censuses contribute qualified coordinates; unproved source frames stay separate");
		check(value(root, "qualified_same_frame_overlap_loci")==7 && value(root, "coordinate_groups_collapsed_by_overlap")==2,
			"The three-interval chain merges transitively; adjacency, opposite strand and other coordinate frames must remain separate");
		check(value(root, "accepted_alias_audit_assembly_groups")==1 && value(root, "scanned_assembly_groups")==2,
			"The old alias-audit check is separate from newly scanned census coverage");
		check(value(root, "distinct_scanned_source_frames")==3 && value(root, "coordinate_comparison_source_frames")==2
			&& value(root, "coordinate_comparison_frame_pairs")==1 && value(root, "coordinate_comparison_nucleotide_bases")==2000,
			"Duplicate-resolution cost inputs must count distinct saved source frames and within-group pairs, not aliases or guessed byte totals");
		check(value(root, "saved_scored_RF00174_annotation_rows")==4, "RF00174 support must survive another representative subtype and preserve prior authority");
		check(value(root, "RF00174_prior_annotations_without_native_GA_hit")==1, "Prior-only annotations remain outside the CM roster");
		check(value(root, "native_GA_events_with_prior_annotation")==1, "Native GA hits covered by a prior remain included");
		check(value(root, "marginal_RF00174_events")==2 && value(root, "marginal_native_GA_events")==1, "Saved quality tiers and native GA are overlapping classifications");
		check(value(root, "exact_member_matches")==1 && value(root, "unmapped_members")==1, "Renamed/reversed cross-file aliases must not become absent or exact by guessed names");
		check(value(root, "current_live_assembly_groups")==3 && value(root, "unscanned_live_assembly_groups")==1
			&& value(root, "unscanned_live_paths")==1, "The newly scanned third source must be covered before inventory; retain the genuinely unscanned source");
		boolean chain=false, adjacent=false, reverse=false;
		int mergedRows=0;
		final LineParser1 components=new LineParser1('\t');
		try(Input input=new Input(root.resolve("result/merged_loci.tsv.gz").toString())){
			check(input.reader.nextLine()!=null, "The merged roster must retain an explicit header");
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				components.set(row);
				check(components.terms()==10, "Merged roster must retain bounds and member-group provenance");
				mergedRows++;
				if(!new String(row, StandardCharsets.UTF_8).contains("\t"+third+"\t")){continue;}
				final int start=components.parseInt(4), stop=components.parseInt(5), groups=components.parseInt(7);
				chain|=start==100 && stop==299 && groups==3 && components.termEquals("+", 6);
				adjacent|=start==300 && stop==320 && groups==1 && components.termEquals("+", 6);
				reverse|=start==100 && stop==199 && groups==1 && components.termEquals("-", 6);
			}
		}
		check(mergedRows==7 && chain && adjacent && reverse,
			"Actual roster rows must retain the transitive union bounds and separate adjacent/opposite-strand observations");
		check(NcrnaDataReconciliation.reciprocal80(80, 100, 100) && !NcrnaDataReconciliation.reciprocal80(79, 100, 100), "Reciprocal eighty-percent boundary must be exact");
		check(NcrnaDataReconciliation.coordinateFrame(new Assembly(0, left, left, false), "CONTIG1.1").equals(left), "Accession-looking strings do not prove a coordinate map");
		boolean unicode=false;
		try(Input input=new Input(root.resolve("result/assemblies.tsv.gz").toString())){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				unicode|=new String(line, StandardCharsets.UTF_8).contains("söurces");
			}
		}
		check(unicode, "UTF-8 source paths must survive decoding and output");
		Files.move(root.resolve("old.witness"), root.resolve("old.witness.valid"));
		write(root, "old.witness", WITNESS_HEADER+witness(exact, 1, 101, 199, "RF00174", -1, "QUALIFYING_SAME_FAMILY_GROUP_MEMBER"));
		parameters[1]="out="+root.resolve("malformed_result");
		boolean rejected=false;
		try{NcrnaDataReconciliation.main(parameters);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("Witness event/input/model/interval/native binding mismatch");}
		check(rejected, "A malformed witness interval must fail against the actual native event");
	}

	private static Row reference(final NcrnaReferenceGff codec, final Path root, final String evidence,
			final int start, final int stop, final boolean mask, final String model){
		return codec.create("1", new Locus("NC_000001.1", "cobalamin", start, stop, '+'), mask,
			mask ? "fixture_marginal" : null, model, root.resolve(evidence).toString());
	}
	private static String event(final int id, final int start, final int stop, final String model,
			final String tier, final boolean included, final int flags){
		return id+"\t1\tNC_000001.1\t"+start+"\t"+stop+"\t+\t"+model+"\tcobalamin\t1\t100\t100\t0\tfalse\tfalse\t60.0\t40.0\t150.0\t"+included+"\tfalse\t100\ttrue\tfalse\t"+tier+"\tfixture\t0\t"+flags+"\tnative.tbl\t"+id+"\tfixture\t6\n";
	}
	private static String witness(final String id, final int event, final int start, final int stop,
			final String model, final int prior, final String association){
		return id+"\t1\tCM_HIT\t"+event+"\tnative.tbl\t"+event+"\t"+model+"\t"+start+"\t"+stop+"\t"+prior+"\t"+association+"\n";
	}
	private static String prior(final String id, final int start, final int stop, final int line){
		return id+"\t1\tPROTECTED_PRIOR\t.\t.\t.\tRF00174\t"+start+"\t"+stop+"\t"+line+"\tORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY\n";
	}
	private static void write(final Path root, final String name, final String text){
		try(Output output=new Output(root.resolve(name).toString())){output.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private static long value(final Path root, final String key){
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(root.resolve("result/summary.tsv").toString())){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.termEquals(key, 0)){return parser.parseLong(1);}
			}
		}
		throw new AssertionError("Missing fixture summary metric: "+key);
	}
	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
	private static final String GROUP_HEADER="cohort\tinput_id\tfasta\tcanonical_fasta\tassembly_record\tgroup_id\tevaluation_group\tsequence_comparison_pending\tunresolved_claim\tdisposition\tresolved_group_id\tresolved_evaluation_group\tgroup_has_unresolved_pair\tresolved_assembly_disposition\n";
	private static final String CENSUS_HEADER="input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records\n";
	private static final String GENOME_HEADER="index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs\n";
	private static final String WITNESS_HEADER="annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation\n";
	private static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision\n";
}
