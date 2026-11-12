package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import prok.BuildNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Literal interval/edge oracles for saved geometry, with no readable genome payload.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceGeometryTest{

	NcrnaReferenceGeometryTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Geometry fixture root must be new");
		prepare();
		final String hits=hit(1, 1, TARGET, 10, 19, '+', A, 1, 10, "VALID_INTERVAL")
			+hit(1, 2, TARGET, 10, 19, '+', B, 1, 10, "VALID_INTERVAL")
			+hit(1, 3, TARGET, 30, 39, '+', A, 1, 10, "VALID_INTERVAL")
			+hit(1, 4, TARGET, 18, 31, '+', A, 2, 8, "VALID_INTERVAL")
			+hit(1, 5, TARGET, 5, 45, '+', B, 1, 10, "VALID_INTERVAL")
			+hit(1, 6, TARGET, 12, 22, '+', A, 1, 10, "VALID_INTERVAL")
			+hit(1, 7, TARGET, 20, 25, '-', C, 1, 6, "VALID_INTERVAL")
			+hit(1, 8, TARGET, 90, 110, '+', A, 1, 10, "OUT_OF_RANGE")
			+hit(1, 9, TARGET, 95, 105, '+', A, 1, 10, "OUT_OF_RANGE")
			+hit(1, 10, "Aa", 1, 10, '+', A, 1, 10, "VALID_INTERVAL")
			+hit(1, 11, "BB", 1, 10, '+', A, 1, 10, "VALID_INTERVAL")
			+hit(2, 1, TARGET, 10, 19, '+', A, 1, 10, "VALID_INTERVAL");
		final String pairs=pair(1, 2, 10, false)+pair(1, 4, 2, false)+pair(1, 5, 10, false)+pair(1, 6, 8, false)
			+pair(2, 4, 2, false)+pair(2, 5, 10, false)+pair(2, 6, 8, false)+pair(3, 4, 2, false)
			+pair(3, 5, 10, false)+pair(4, 5, 14, false)+pair(4, 6, 5, false)+pair(5, 6, 11, false)
			+pair(4, 7, 6, true)+pair(5, 7, 6, true)+pair(6, 7, 3, true)+pair(8, 9, 11, false);
		writeEvidence("good", SUMMARY, COVERAGE, SOURCES, hits, pairs);
		invoke("good", "result");
		check(lines("result.geometry_summary.tsv")[1].equals("3\t3\t3\t12\t10\t2\t5\t16\t12\t3\t1\t2\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tDESCRIPTIVE_GEOMETRY_NOT_GENES"),
			"Twelve events form five review groups, with two invalid events and two bridge witnesses");
		final String[] members=lines("result.group_members.tsv.gz");
		check(members.length==13, "Every original event must survive, including invalid intervals");
		final String main=find(members, "\t1\t/missing/search1/hits.tbl\t1\t");
		final String same=find(members, "\t1\t/missing/search1/hits.tbl\t3\t");
		final String otherInput=find(members, "\t2\t/missing/search2/hits.tbl\t1\t");
		check(group(main).equals(group(same)) && !group(main).equals(group(otherInput)),
			"Chaining joins disjoint events for review; equal native lines/targets in another input remain separate");
		check(main.contains(TARGET) && !main.contains("%2525"), "Saved UTF-8/percent framing must survive exactly once");
		check(!group(find(members, "\tAa\t")).equals(group(find(members, "\tBB\t"))),
			"Equal hash codes for Aa/BB must not merge distinct target bytes");
		check(find(members, "\t1\t/missing/search1/hits.tbl\t8\t").startsWith("NA\t"), "Invalid intervals have no group assignment");
		final String[] geometry=lines("result.model_pairs.tsv.gz");
		check(geometry.length==17, "Every original model edge must survive, including opposite-strand and invalid geometry");
		// ByteBuilder.append(double,int) emits integral values without fractional zeros.
		check(find(geometry, "\t1\t2\t"+A+"\t"+B+"\t").endsWith("\tEQUAL\t10\t10\t10\t1\t1\tSAME_STRAND\ttrue"),
			"Co-located full-model alternatives are equal geometry, not tandem evidence");
		check(find(geometry, "\t1\t6\t"+A+"\t"+A+"\t").endsWith("\tPARTIAL_OVERLAP\t8\t10\t11\t0.800000\t0.727273\tSAME_STRAND\ttrue"),
			"Partially overlapping full hits retain exact overlap numerators/denominators and descriptive fractions");
		check(find(geometry, "\t8\t9\t"+A+"\t"+A+"\t").contains("\tUNASSESSABLE\t"), "Out-of-range pairs are not used for grouping");
		final String[] bridges=lines("result.bridges.tsv.gz");
		check(bridges.length==3, "Only the shorter touching hit and outer covering hit bridge disjoint full events");
		check(find(bridges, "\t4\t1\t3\t").endsWith("\tFMN\t2\t2\t10\tfalse\tCANDIDATE_BRIDGE_NOT_FUSION"),
			"A bridge touching two portions does not cover their outer span");
		check(find(bridges, "\t5\t1\t3\t").endsWith("\tFMN\t10\t10\t10\ttrue\tCANDIDATE_BRIDGE_NOT_FUSION"),
			"Outer-span coverage is separate from the decision to keep both original hits");
		strataAlternatives();
		failure("duplicate", hits, pairs+pair(1, 2, 10, false), "Repeated saved model pair");
		failure("missing_edge", hits, pairs.replace(pair(1, 2, 10, false), ""), "overlap counts disagree");
		failure("unknown_hit", hits, pairs.replace(pair(1, 2, 10, false), pair(1, 99, 10, false)), "two distinct retained events");
		failure("wrong_target", hits, pairs.replace(pair(1, 2, 10, false), pair(10, 11, 10, false)), "target/overlap disagrees");
		failure("wrong_overlap", hits, pairs.replace(pair(1, 2, 10, false), pair(1, 2, 9, false)), "target/overlap disagrees");
		failure("duplicate_hit", hits.replace("1\t/missing/search1/hits.tbl\t2\t", "1\t/missing/search1/hits.tbl\t1\t"), pairs, "Repeated input/table/native-line");
		writeEvidence("sample", SUMMARY.replace("3\t3\t300\t3", "2\t3\t200\t2")
			.replace("\t12\t3\t3\t", "\t12\t2\t2\t").replace("\t9\t5\tEVIDENCE_COMPLETE", "\t6\t2\tSAMPLE_EVIDENCE"),
			COVERAGE.replace("\t3\n", "\t2\n"), SOURCES.substring(0, SOURCES.indexOf("3\tNA")), hits, pairs);
		invoke("sample", "sample_result");
		check(lines("sample_result.geometry_summary.tsv")[1].equals("2\t3\t3\t12\t10\t2\t5\t16\t12\t3\t1\t2\tSAMPLE_EVIDENCE_NOT_ADJUDICATED\tDESCRIPTIVE_GEOMETRY_NOT_GENES"),
			"A sampled geometry report must keep the original catalog denominator and scope");
		writeEvidence("zero", "3\t3\t300\t3\t0\t3\t3\t0\t0\t0\t0\t0\t3\t9\t9\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n",
			A+"\tfixture\t10\tNA\tNA\t0\t0\t3\n"+B+"\tfixture\t10\tNA\tNA\t0\t0\t3\n"+C+"\tfixture\t6\tNA\tNA\t0\t0\t3\n", SOURCES, "", "");
		invoke("zero", "zero_result");
		check(lines("zero_result.geometry_summary.tsv")[1].equals("3\t3\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tDESCRIPTIVE_GEOMETRY_NOT_GENES"),
			"Completed zero-hit inputs/models must survive without invented groups");
	}

	private void strataAlternatives(){
		write("good.historical.tsv.gz", BuildNcrnaReference.HISTORICAL_HEADER+'\n');
		write("good.relations.tsv.gz", BuildNcrnaReference.RELATION_HEADER+'\n');
		write("capture.attributes.tsv.gz", "input_id\tgff\tline\traw_attributes\n");
		write("capture.attribute_summary.tsv", "inputs_observed\tcatalog_inputs\tretained_source_rows\tmatched_attribute_rows\tsource_paths_read\tunique_selected_lines\tsource_lines_read\tevidence_status\tstatus\n"
			+"3\t3\t0\t0\t0\t0\t0\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED\n");
		BuildNcrnaReference.main(new String[]{"mode=strata", "in="+path("good"), "attributes="+path("capture"), "out="+path("strata"),
			"census="+path("census.tsv"), "memberships="+path("members.tsv"), "scope="+path("scope.tsv"),
			"assets="+path("assets.tsv"), "queries="+path("queries.tsv"), "t=1", "bgzip=f", "zipthreads=1"});
		final String event=find(lines("strata.event_strata.tsv.gz"), "1\t/missing/search1/hits.tbl\t1\t");
		check(event.contains("OVERLAPPING_SAME_FAMILY_OTHER_MODEL") && event.contains("OVERLAPPING_SAME_MODEL"),
			"Same-model overlaps and same-family alternative CMs are independent observations");
		check(find(lines("strata.joint_strata.tsv"), "RF00050\t100_percent\t0\t20_to_below50\tNO_AVAILABLE_GFF\tNONE\tNONE\ttrue\ttrue\tfalse\tfalse\tfalse\t")
			.endsWith("\t2\t1\tDESCRIPTIVE_BINS_NONADDITIVE_GROUPS"), "Two same-bin events in the complex review group must count as one group");
		check(find(lines("strata.joint_strata.tsv"), "RF00050\t100_percent\t0\t20_to_below50\tUNASSESSABLE_EVENT\t")
			.endsWith("\t2\t0\tDESCRIPTIVE_BINS_NONADDITIVE_GROUPS"), "Invalid events have no invented review-group membership");
	}

	private void prepare(){
		final StringBuilder census=new StringBuilder(NcrnaEvidenceCatalog.CENSUS_HEADER).append('\n');
		final StringBuilder members=new StringBuilder(NcrnaEvidenceCatalog.MEMBERSHIP_HEADER).append('\n');
		for(int i=1; i<=3; i++){
			final String fasta=path("absent"+i+".fa");
			check(!new File(fasta).exists(), "Geometry must work without genomic payloads");
			census.append(i).append('\t').append(fasta).append("\t110\t100\t1\n");
			members.append("fixture\t").append(i).append('\t').append(i).append('\t').append(fasta).append("\tNA\tNA\tNA\tUNKNOWN\tNA\tNA\n");
		}
		write("census.tsv", census.toString());
		write("members.tsv", members.toString());
		write("scope.tsv", NcrnaReferenceEvidence.SCOPE_HEADER+'\n'+A+"\tvariant_A\tprimary\tfixture\tFMN\tTARGET_CLASS\n"
			+B+"\tvariant_B\tprimary\tfixture\tFMN\tTARGET_CLASS\n"+C+"\tother\tprimary\tfixture\t5S\tRELATED_CLASS\n");
		write("assets.tsv", NcrnaReferenceEvidence.ASSET_HEADER+'\n'+A+"\tfixture\t10\t20\t00000000000000000001\t/missing/a.cm\n"
			+B+"\tfixture\t10\t20\t00000000000000000002\t/missing/b.cm\n"+C+"\tfixture\t6\t20\t00000000000000000003\t/missing/c.cm\n");
		write("queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+'\n'+A+"\ta\t10\n"+B+"\tb\t10\n"+C+"\tc\t6\n");
	}
	private static String hit(final int input, final int line, final String target, final int start, final int stop,
		final char strand, final String model, final int from, final int to, final String range){
		return input+"\t/missing/search"+input+"/hits.tbl\t"+line+'\t'+target+'\t'+start+'\t'+stop+'\t'+strand+'\t'+model
			+"\tfixture\t"+(model.equals(C) ? "5S\tRELATED_CLASS" : "FMN\tTARGET_CLASS")+'\t'+from+'\t'+to+'\t'
			+(model.equals(C) ? 6 : 10)+"\t0\t50.000000\t1e-8\ttrue\t"+range+'\t'+(start==1)+"\tfalse\n";
	}
	private static String pair(final int a, final int b, final int overlap, final boolean opposite){
		return "1\tMODEL_MODEL\t/missing/search1/hits.tbl\t"+a+"\t/missing/search1/hits.tbl\t"+b+'\t'+overlap+'\t'
			+(opposite ? "OPPOSITE_STRAND_GEOMETRY" : "SAME_STRAND_GEOMETRY")+"\tfalse\n";
	}
	private void failure(final String name, final String hits, final String pairs, final String reason){
		writeEvidence(name, SUMMARY, COVERAGE, SOURCES, hits, pairs);
		boolean failed=false;
		try{invoke(name, name+"_result");}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(reason), "Unexpected fixture failure: "+expected);
			failed=true;
		}
		check(failed && !new File(path(name+"_result.groups.tsv.gz")).exists(), "Invalid event joins must fail before writing groups: "+name);
	}
	private void writeEvidence(final String name, final String summary, final String coverage, final String sources, final String hits, final String pairs){
		write(name+".summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+'\n'+summary);
		write(name+".coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+'\n'+coverage);
		write(name+".sources.tsv.gz", BuildNcrnaReference.SOURCE_HEADER+'\n'+sources);
		write(name+".hits.tsv.gz", BuildNcrnaReference.HIT_HEADER+'\n'+hits);
		write(name+".overlaps.tsv.gz", BuildNcrnaReference.OVERLAP_HEADER+'\n'+pairs);
	}
	private void invoke(final String input, final String output){
		BuildNcrnaReference.main(new String[]{"mode=geometry", "in="+path(input), "out="+path(output), "census="+path("census.tsv"),
			"memberships="+path("members.tsv"), "scope="+path("scope.tsv"), "assets="+path("assets.tsv"), "queries="+path("queries.tsv"),
			"t=1", "bgzip=f", "zipthreads=1"});
	}
	private void write(final String name, final String text){try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}}
	private String[] lines(final String name){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){rows.add(new String(row, StandardCharsets.UTF_8));}
		}
		return rows.toArray(new String[rows.size()]);
	}
	private static String find(final String[] rows, final String fragment){
		String found=null;
		for(final String row:rows){
			if(row.contains(fragment)){
				check(found==null, "Fixture key must identify exactly one row: "+fragment);
				found=row;
			}
		}
		check(found!=null, "Missing fixture row: "+fragment);
		return found;
	}
	private static String group(final String row){return row.substring(0, row.indexOf('\t'));}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
	private static final String A="RF00050", B="RF00051", C="RF00001", TARGET="Ω%25seq";
	private static final String SUMMARY="3\t3\t300\t3\t12\t3\t3\t0\t0\t0\t0\t16\t3\t9\t5\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n";
	private static final String COVERAGE=A+"\tfixture\t10\t10\t0\t8\t9\t3\n"+A+"\tfixture\t10\t7\t0\t1\t9\t3\n"
		+B+"\tfixture\t10\t10\t0\t2\t2\t3\n"+C+"\tfixture\t6\t6\t0\t1\t1\t3\n";
	private static final String SOURCES="1\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n"
		+"2\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n"+"3\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n";
}
