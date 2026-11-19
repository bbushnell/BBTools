package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import fileIO.FileFormat;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaRosterMembers;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import prok.RrnaResourceIO;
import shared.Shared;
import stream.Read;
import stream.StreamerFactory;
import structures.ListNum;

/** Full native extraction/ledger fixture, including construction equality scope,
 * reverse inclusive bounds, ambiguity and conservation of original called genes.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaRosterMembersTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaRosterMembersTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, split[0], split.length>1 ? split[1] : null), "Unknown fixture argument");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Supply out=fresh_fixture_directory");
			fixture(Files.createDirectory(new File(parser.out1).toPath()));
			System.err.println("NCRNA_ROSTER_MEMBERS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void fixture(final Path root) throws Exception{
		final Path sources=Files.createDirectory(root.resolve("söurces"));
		final Path a=sources.resolve("a.fa"), b=sources.resolve("b.fa"), c=sources.resolve("c.fa");
		write(a, ">alpha description\naAcCgGtT\n>dup\nAAAA\n>dup repeated\nCCCC\n");
		write(b, ">beta\nAACCGGTT\n>variant\nAANAATCGATT\n>ambiguous\nNCGT\n>defined\nAAAAANACGTAAAAAA\n");
		write(c, ">gamma\nAACCGGTT\n>trunc\nAACC\n>span\nGGTT\n");
		final Path alias=sources.resolve("a_link.fa"), other=sources.resolve("other.fa");
		Files.createSymbolicLink(alias, a.getFileName());
		final Path input=Files.createDirectory(root.resolve("metadata"));
		write(input.resolve("summary.tsv"), "metric\tcount\tqualification\n"
			+"raw_RF00174_events\t13\tfixture\nnative_GA_RF00174_events\t12\tfixture\n"
			+"qualified_exact_coordinate_groups\t12\tfixture\nqualified_same_frame_overlap_loci\t12\tfixture\n");
		write(input.resolve("source_frames.tsv.gz"), FRAME_HEADER+frame(0, a, alias, 16, 3)+frame(0, b, 39, 4)+frame(1, c, 16, 3));
		final String[] contigs={"alpha","beta","variant","gamma","missing","dup","alpha","beta","ambiguous","defined","trunc","span"};
		final Path[] paths={a,b,b,c,a,a,a,b,b,b,c,c};
		final int[] starts={2,4,7,4,1,1,1,3,1,7,1,1}, stops={5,7,10,7,3,4,99,7,4,10,4,4};
		final StringBuilder merged=new StringBuilder(MERGED_HEADER), links=new StringBuilder("exact_coordinate_group_id\tmerged_locus_id\n");
		final StringBuilder originals=new StringBuilder(ORIGINAL_HEADER);
		final StringBuilder events=new StringBuilder(EVENT_HEADER);
		for(int id=0; id<12; id++){
			final int group=id==3 || id>=10 ? 1 : 0;
			final char strand=id==0 || id==9 ? '-' : '+';
			merged.append(id).append('\t').append(group).append('\t').append(paths[id]).append('\t').append(contigs[id]).append('\t')
				.append(starts[id]).append('\t').append(stops[id]).append('\t').append(strand).append("\t1\t1\tfalse\n");
			links.append(id).append('\t').append(id).append('\n');
			originals.append(original(id, group, paths[id], contigs[id], starts[id], stops[id], strand));
			final int source=paths[id].equals(a) ? 1 : paths[id].equals(b) ? 2 : 3;
			events.append(event(id, source, group, contigs[id], starts[id], stops[id], strand));
		}
		write(input.resolve("merged_loci.tsv.gz"), merged.toString());
		write(input.resolve("locus_components.tsv.gz"), links.toString());
		write(input.resolve("loci.tsv.gz"), originals.toString());
		write(input.resolve("events.tsv.gz"), events.toString());
		write(input.resolve("census.tsv"), "input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records\n1\t"+alias+"\t0\t16\t3\n2\t"+b+"\t0\t39\t4\n3\t"+c+"\t0\t16\t3\n");
		final Path cohorts=input.resolve("cohorts.tsv");
		write(cohorts, "cohort\tcensus\treference\tevents\twitnesses\texpected_inputs\texpected_scored_RF00174\nfixture\t"
			+input.resolve("census.tsv")+"\tDO_NOT_READ_ANNOTATIONS\tDO_NOT_READ_OTHER_EVENTS\tDO_NOT_READ_WITNESSES\t3\t-1\n");
		write(other, ">alpha description\nTTTTTTTT\n>dup\nGGGG\n>dup repeated\nCCCC\n");
		Files.move(alias, sources.resolve("a_link.before"));
		Files.createSymbolicLink(alias, other.getFileName());
		check(Files.isSameFile(alias, other) && !Files.isSameFile(alias, a),
			"The fixture must actually retarget the supplied alias to different bases with the same record/base dimensions");
		final Path metadata=root.resolve("eligibility");
		NcrnaRosterMembers.main(new String[]{"in="+input, "cohorts="+cohorts, "out="+metadata, "expectedloci=12", "mode=metadata"});
		check(!Files.exists(metadata.resolve("members.fa.gz")) && value(metadata, "eligibility.tsv", "original_called_coordinate_groups")==12,
			"Metadata mode must conserve the whole roster without creating sequence outputs");
		check(value(metadata, "eligibility.tsv", "infernal_truncated_only")==0
			&& value(metadata, "eligibility.tsv", "incomplete_model_coverage_only")==3
			&& value(metadata, "eligibility.tsv", "infernal_truncated_and_incomplete_model_coverage")==1
			&& value(metadata, "eligibility.tsv", "neither_partial_diagnostic")==8
			&& value(metadata, "eligibility.tsv", "construction_candidates")==11,
			"Literal model start>5/end<length-5 cutoffs must allow start5/end15 of20, including contig-edge contact");
		final Path strict=root.resolve("strict");
		NcrnaRosterMembers.main(new String[]{"in="+input, "cohorts="+cohorts, "out="+strict, "expectedloci=12", "mode=metadata", "modeltolerance=1", "coveragefilter=t"});
		check(value(strict, "eligibility.tsv", "infernal_truncated_and_incomplete_model_coverage")==1
			&& value(strict, "eligibility.tsv", "incomplete_model_coverage_total")==5,
			"A changed tolerance must reach actual eligibility, including overlap between truncation and coverage");
		final Path stopped=root.resolve("stopped_before_sequences");
		boolean rejected=false;
		try{
			NcrnaRosterMembers.main(new String[]{"in="+input, "cohorts="+cohorts, "out="+stopped, "expectedloci=12", "maxcoveragefraction=0", "coveragefilter=t"});
		}catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("no sequence was read");}
		check(rejected && Files.exists(stopped.resolve("eligibility.tsv")) && !Files.exists(stopped.resolve("members.fa.gz")),
			"Coverage-fraction refusal must save its metadata evidence before any FASTA extraction");
		final Path full=root.resolve("full");
		NcrnaRosterMembers.main(new String[]{"in="+input, "cohorts="+cohorts, "out="+full, "expectedloci=12", "t=2", "pigz=f", "unpigz=f", "coveragefilter=f", "maxcoveragefraction=0"});
		check(value(full, "original_called_coordinate_groups")==12 && value(full, "merged_loci")==12,
			"Every original native called gene must retain its overlap-component identity");
		check(value(full, "extracted_loci")==9 && value(full, "native_eligible_extracted_loci")==8 && value(full, "construction_representatives")==7
			&& value(full, "removed_construction_sequence_copies")==1 && value(full, "unresolved_merged_loci")==3,
			"Only same-group exact full sequences collapse; missing, ambiguous and out-of-bounds sources remain unresolved");
		check(value(full, "unresolved_original_calls")==3, "Unresolved original calls must remain explicitly counted");
		check(value(full, "partial_or_truncated_original_calls")==1 && value(full, "construction_ineligible_components")==1,
			"With coverage filtering off, only native truncation excludes construction; incomplete model coverage stays diagnostic");
		check(value(full, "eligible_identical_sequence_classes_across_groups")==1
			&& value(full, "eligible_loci_in_cross_group_identical_classes")==3
			&& value(full, "retained_representatives_in_cross_group_identical_classes")==2,
			"Cross-group identical sequence is counted but both accepted groups keep their representative");
		final ObjectMap<String, byte[]> sequences=sequences(full.resolve("members.fa.gz"));
		check(sequences.size()==7 && string(sequences.get("RF00174_X_0")).equals("CGGT") && !sequences.contains("RF00174_X_1"),
			"Reverse-strand inclusive bounds and same-group cross-frame equality must choose the stable original representative");
		check(string(sequences.get("RF00174_X_2")).equals("CGAT") && string(sequences.get("RF00174_X_3")).equals("CGGT")
			&& string(sequences.get("RF00174_X_7")).equals("CCGGT"), "Different sequences/lengths and equal genes in other assembly groups remain separate");
		check(string(sequences.get("RF00174_X_8")).equals("NCGT") && string(sequences.get("RF00174_X_9")).equals("ACGT"),
			"Ambiguous symbols are exact bytes, not wildcard matches");
		check(string(sequences.get("RF00174_X_11")).equals("GGTT") && !sequences.contains("RF00174_X_10"),
			"Coverage-only native genes remain in construction while Infernal-truncated genes do not; an off coverage filter also disables its fraction stop");
		final LineParser1 parser=new LineParser1('\t');
		int rows=0, linkedOriginals=0;
		try(Input ledger=new Input(full.resolve("original_calls.tsv.gz").toString())){
			check(ledger.reader.nextLine()!=null, "Original-call ledger requires its source schema");
			for(byte[] row=ledger.reader.nextLine(); row!=null; row=ledger.reader.nextLine()){
				parser.set(row);
				check(parser.terms()==29, "Original calls must retain native eligibility, strand-aware flanks, adjacent N runs and bit scores");
				if(parser.parseInt(0)==2){check(parser.termEquals("INTERNAL_BOTH_FLANKS", 27), "Internal incomplete sequence must retain both actual genomic flanks");}
				if(parser.parseInt(0)==11){check(parser.termEquals("CONTIG_END_REACH", 27), "Missing native model positions exceeding available sequence must retain the contig-end reach class");}
				if(parser.parseInt(0)==9){
					check(parser.termEquals("ADJACENT_N_RUN", 27) && parser.parseInt(25)==0 && parser.parseInt(26)==1,
						"Minus-strand model-right context is the lower-coordinate adjacent genomic N run");
				}
				rows++;
				if(parser.parseInt(0)==0 || parser.parseInt(0)==1){
					check(parser.parseInt(13)==parser.parseInt(0) && parser.parseInt(14)==0,
						"Both original calls must retain their own native component while sharing only the construction representative");
					linkedOriginals++;
				}
			}
		}
		check(rows==12 && linkedOriginals==2, "No original call may vanish during construction sequence deduplication");
		check(value(full, "assembly_loss.tsv", "groups_left_zero_removing_all_incomplete")==1
			&& value(full, "assembly_loss.tsv", "groups_left_zero_removing_internal_incomplete")==0,
			"Assembly loss uses every original call, including retained coverage variants and truncated evidence");
		try(Input scores=new Input(full.resolve("bit_scores.tsv").toString())){
			boolean found=false;
			for(byte[] row=scores.reader.nextLine(); row!=null; row=scores.reader.nextLine()){
				parser.set(row);
				if(parser.termEquals("INTERNAL_INCOMPLETE", 0)){
					check(parser.parseInt(1)==1 && parser.parseDouble(4)==30, "Internal incomplete bit-score median must use the retained native score");
					found=true;
				}
			}
			check(found, "The internal incomplete bit-score distribution must be reported explicitly");
		}
		int unresolvedReasons=0;
		try(Input ledger=new Input(full.resolve("merged_members.tsv.gz").toString())){
			check(ledger.reader.nextLine()!=null, "Merged extraction ledger needs its complete schema");
			for(byte[] row=ledger.reader.nextLine(); row!=null; row=ledger.reader.nextLine()){
				parser.set(row);
				final int id=parser.parseInt(0);
				if(id==4 || id==5 || id==6){
					check(parser.parseInt(8)==-1 && parser.termEquals(id==4 ? "CONTIG_NOT_FOUND" : id==5 ? "AMBIGUOUS_CONTIG_ID" : "OUT_OF_RANGE_COORDINATES", 9),
						"Different unresolved mapping causes must retain explicit reasons and no invented representative");
					unresolvedReasons++;
				}
			}
		}
		check(unresolvedReasons==3, "All missing, ambiguous and out-of-range calls need their own ledger rows");
		final Path canary=root.resolve("canary");
		NcrnaRosterMembers.main(new String[]{"in="+input, "cohorts="+cohorts, "out="+canary, "expectedloci=12", "mode=canary", "maxsources=1", "t=1", "pigz=f", "unpigz=f"});
		check(value(canary, "skipped_operational_canary_loci")==8 && value(canary, "original_called_coordinate_groups")==12,
			"An operational canary must label its skipped sources while preserving the complete input population");
		Files.move(c, sources.resolve("c.fa.saved"));
		final Path missing=root.resolve("missing_source");
		NcrnaRosterMembers.main(new String[]{"in="+input, "cohorts="+cohorts, "out="+missing, "expectedloci=12", "t=2", "pigz=f", "unpigz=f"});
		check(value(missing, "unresolved_merged_loci")==6 && value(missing, "unresolved_original_calls")==6
			&& value(missing, "original_called_coordinate_groups")==12, "Unavailable source files must retain their original called genes as unresolved");
	}

	private static String frame(final int group, final Path path, final long bases, final int records){
		return frame(group, path, path, bases, records);
	}
	private static String frame(final int group, final Path canonical, final Path supplied, final long bases, final int records){
		return group+"\t"+canonical+"\t"+supplied+"\t0\t"+bases+"\t"+records+"\t1\t"+(group==0 ? 2 : 1)+"\t1\ttrue\t"+(group==0)+"\ttrue\tfixture\n";
	}
	private static String original(final int id, final int group, final Path path, final String contig, final int start, final int stop, final char strand){
		return id+"\t"+group+"\t"+path+"\t"+contig+"\t"+start+"\t"+stop+"\t"+strand+"\ttrue\t"+(id==0)+"\t1\t0\t0\tfixture\n";
	}
	private static String event(final int id, final int source, final int group, final String contig, final int start, final int stop, final char strand){
		final String[] fields=new String[35];
		Arrays.fill(fields, "0");
		fields[0]="fixture";
		fields[1]=""+group;
		fields[2]="true";
		fields[3]=""+(id==0);
		fields[5]=""+(id+1);
		fields[6]=""+source;
		fields[7]=contig;
		fields[8]=""+start;
		fields[9]=""+stop;
		fields[10]=""+strand;
		fields[11]="RF00174";
		fields[12]="cobalamin";
		fields[13]=""+(id==2 || id==3 || id==10 || id==11 ? 6 : id==9 ? 5 : 1);
		fields[14]=""+(id==9 || id==10 ? 15 : 20);
		fields[15]="20";
		fields[16]=""+(id==10 ? 1 : 0);
		fields[17]=fields[18]=""+(id==9);
		fields[19]=""+(10+id*10);
		fields[22]="true";
		return String.join("\t", fields)+"\n";
	}
	private static ObjectMap<String, byte[]> sequences(final Path path){
		final ObjectMap<String, byte[]> values=new ObjectMap<String, byte[]>(String.class, byte[].class);
		final FileFormat format=FileFormat.testInput(path.toString(), FileFormat.FASTA, null, true, true);
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 0))){
			input.stream.start();
			for(ListNum<Read> batch=input.stream.nextList(); batch!=null; batch=input.stream.nextList()){
				for(final Read read:batch.list){
					final String id=RrnaResourceIO.first(read.id);
					check(!values.contains(id), "Construction FASTA IDs must remain unique");
					values.put(id, read.bases);
				}
				input.stream.returnList(batch);
			}
		}
		return values;
	}
	private static long value(final Path root, final String key){
		return value(root, "summary.tsv", key);
	}
	private static long value(final Path root, final String table, final String key){
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(root.resolve(table).toString())){
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				parser.set(row);
				if(parser.termEquals(key, 0)){return parser.parseLong(1);}
			}
		}
		throw new AssertionError("Missing extraction metric: "+key);
	}
	private static String string(final byte[] bases){
		check(bases!=null, "Expected construction representative missing");
		return new String(bases, StandardCharsets.US_ASCII);
	}
	private static void write(final Path path, final String text){try(Output output=new Output(path.toString())){output.writer.print(text.getBytes(StandardCharsets.UTF_8));}}
	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
	private static final String FRAME_HEADER="assembly_group\tcoordinate_frame\tsupplied_fasta\tfile_bytes\tnucleotide_bases\tfasta_records\tcohort_input_rows\tframes_in_group\tnative_GA_exact_groups\tassembly_has_native_GA\tcoordinate_comparison_needed\tsource_available\treason\n";
	private static final String MERGED_HEADER="locus_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\texact_coordinate_groups\tnative_GA_events\tprotected_prior_witness\n";
	private static final String ORIGINAL_HEADER="locus_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\tqualified\tprior_annotation_overlap\traw_events\texact_members\treciprocal80_members\tdisposition\n";
	private static final String EVENT_HEADER="cohort\tassembly_group\tnative_GA_qualified\tprior_annotation_overlap\tqualifying_annotation_links\tevent_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision\n";
}
