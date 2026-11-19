package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaReference;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Finite hand-specified correspondence graph; nonexistent FASTAs prove the
 * comparison reads metadata only. Tests include a long source/short core pair,
 * duplicate aliases, both strands, one-to-many/many-to-one/many-to-many paths,
 * masks, rejected/unresolved sources, unassigned families and hitless priors.
 * @author Ganyu
 */
public final class NcrnaSourceComparisonTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSourceComparisonTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Supply out=fresh_fixture_directory, without input or overwrite");
			new NcrnaSourceComparisonTest(new File(parser.out1).getAbsoluteFile()).run();
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaSourceComparisonTest(final File root_){
		root=root_;
		check(root.mkdir(), "Fixture directory must be fresh");
	}

	private void run(){
		final Row same=parent(100, 199, '+', false);
		final int sameEvent=event(same, "CONFIDENT");
		final int duplicateModel=event(same, "CONFIDENT");
		final int a=source(100, 199, '+', "FMN", 3);
		final int alias=source(100, 199, '+', "FMN", 3);
		link(a, sameEvent, "SAME_FAMILY", "CONFIDENT");
		link(a, duplicateModel, "SAME_FAMILY", "CONFIDENT");
		link(alias, sameEvent, "SAME_FAMILY", "CONFIDENT");
		final Row minus=parent(300, 399, '-', false);
		link(source(280, 419, '-', "FMN", 3), event(minus, "CONFIDENT"), "SAME_FAMILY", "CONFIDENT");
		final Row core=parent(500, 509, '+', false);
		link(source(450, 649, '+', "FMN", 3), event(core, "CONFIDENT"), "SAME_FAMILY", "CONFIDENT");
		final Row left=parent(700, 749, '+', false);
		final Row right=parent(800, 849, '+', false);
		final int split=source(690, 859, '+', "FMN", 3);
		link(split, event(left, "CONFIDENT"), "SAME_FAMILY", "CONFIDENT");
		link(split, event(right, "CONFIDENT"), "SAME_FAMILY", "CONFIDENT");
		final int merged=event(parent(900, 999, '+', false), "CONFIDENT");
		link(source(900, 949, '+', "FMN", 3), merged, "SAME_FAMILY", "CONFIDENT");
		link(source(950, 999, '+', "FMN", 3), merged, "SAME_FAMILY", "CONFIDENT");
		link(source(1090, 1209, '+', "FMN", 4), event(parent(1100, 1199, '+', true), "MARGINAL"), "SAME_FAMILY", "MARGINAL");
		event(parent(1300, 1399, '+', false), "CONFIDENT");
		link(source(1500, 1599, '+', "UNASSIGNED", 3), event(parent(1500, 1599, '+', false), "CONFIDENT"), "UNASSIGNED_SOURCE", "CONFIDENT");
		prior(parent(1700, 1799, '+', false));
		source(1700, 1799, '+', "FMN", 2);
		final Row protectedHit=parent(1900, 1999, '+', false);
		prior(protectedHit);
		event(protectedHit, "COVERED_BY_PROTECTED_PRIOR");
		source(2100, 2199, '+', "FMN", 9);
		source(2300, 2399, '+', "FMN", 3);
		final int e1=event(parent(2500, 2599, '+', false), "CONFIDENT");
		final int e2=event(parent(2600, 2699, '+', false), "CONFIDENT");
		final int s1=source(2500, 2699, '+', "FMN", 3);
		final int s2=source(2510, 2709, '+', "FMN", 3);
		link(s1, e1, "SAME_FAMILY", "CONFIDENT");
		link(s1, e2, "SAME_FAMILY", "CONFIDENT");
		link(s2, e1, "SAME_FAMILY", "CONFIDENT");
		link(s2, e2, "SAME_FAMILY", "CONFIDENT");
		final Row tandemLeft=parent(3000, 3099, '+', false);
		final Row tandemRight=parent(3200, 3299, '+', false);
		final int tandemSourceLeft=source(3000, 3099, '+', "FMN", 3);
		final int tandemSourceRight=source(3200, 3299, '+', "FMN", 3);
		link(tandemSourceLeft, event(tandemLeft, "CONFIDENT"), "SAME_FAMILY", "CONFIDENT");
		link(tandemSourceRight, event(tandemRight, "CONFIDENT"), "SAME_FAMILY", "CONFIDENT");
		final int bridge=event(tandemRight, "CONFIDENT", 3000, 3299, "TANDEM_BRIDGE_WITNESS_NOT_REPRESENTATIVE");
		link(tandemSourceLeft, bridge, "SAME_FAMILY", "CONFIDENT");
		final Row weakParent=parent(3500, 3599, '+', false);
		event(weakParent, "CONFIDENT");
		final int weak=event(weakParent, "CONFIDENT", 3490, 3509, "WEAK_OVERLAP_WITNESS_NOT_GROUPING_SUPPORT");
		link(source(3480, 3509, '+', "FMN", 3), weak, "SAME_FAMILY", "CONFIDENT");
		link(source(3485, 3509, '+', "UNASSIGNED", 3), weak, "UNASSIGNED_SOURCE", "CONFIDENT");
		final int maskAndWeak=source(3540, 3659, '+', "FMN", 4);
		link(maskAndWeak, event(parent(3550, 3649, '+', true), "MARGINAL"), "SAME_FAMILY", "MARGINAL");
		link(maskAndWeak, event(weakParent, "CONFIDENT", 3540, 3599, "WEAK_OVERLAP_WITNESS_NOT_GROUPING_SUPPORT"), "SAME_FAMILY", "CONFIDENT");
		final Row priorSame=parent(3800, 3899, '+', false);
		prior(priorSame);
		link(source(3800, 3899, '+', "FMN", 2), event(priorSame, "COVERED_BY_PROTECTED_PRIOR"), "SAME_FAMILY", "COVERED_BY_PROTECTED_PRIOR");
		final Row priorChanged=parent(4100, 4199, '-', false);
		prior(priorChanged);
		link(source(4000, 4299, '-', "FMN", 2), event(priorChanged, "COVERED_BY_PROTECTED_PRIOR"), "SAME_FAMILY", "COVERED_BY_PROTECTED_PRIOR");
		link(source(3000, 3099, '+', "FMN", 0), bridge, "SAME_FAMILY", "CONFIDENT");
		codec.write(path("reference.gff"), "1", parents);
		codec.write(path("empty.gff"), "2", new ArrayList<Row>());
		write("reference.tsv", "input_id\tfasta\tgff\n1\t"+path("unopened1.fa")+"\t"+path("reference.gff")+"\n2\t"+path("unopened2.fa")+"\t"+path("empty.gff")+"\n");
		write("source_provenance.tsv.gz", sources);
		write("source_event_links.tsv.gz", links);
		write("event_provenance.tsv.gz", events);
		write("annotation_witnesses.tsv.gz", witnesses);
		write("summary.tsv", "input_id\tscored_annotations\tmasked_annotations\tprotected_priors\tstatus\n1\t17\t2\t4\tDECLARED_POLICY_APPLIED_NOT_INDEPENDENT_BIOLOGICAL_VALIDATION\n2\t0\t0\t0\tDECLARED_POLICY_APPLIED_NOT_INDEPENDENT_BIOLOGICAL_VALIDATION\n");
		write("family_counts.tsv", "family\traw_events\tCONFIDENT\tMARGINAL\tBELOW_MARGINAL\tUNCERTAIN_FLOOR\tREJECTED_WRONG_FAMILY\tCOVERED_BY_PROTECTED_PRIOR\tINVALID_INTERVAL\tscored_annotations\tmasked_annotations\tprotected_priors\nFMN\t22\t17\t2\t0\t0\t0\t3\t0\t17\t2\t4\n");
		final ByteBuilder counts=new ByteBuilder("declared_family\tsource_rows");
		for(final String state:STATES){counts.tab().append(state);}
		counts.nl();
		for(int family=0; family<2; family++){
			counts.append(family==0 ? "FMN" : "UNASSIGNED").tab().append(family==0 ? 20 : 2);
			for(final long count:sourceCounts[family]){counts.tab().append(count);}
			counts.nl();
		}
		write("source_family_counts.tsv", counts);
		compare("result");
		final ArrayList<String> sourceResults=lines("result.sources.tsv.gz");
		check(sourceResults.size()==21, "Twenty-one selected claims collapse to twenty physical sources; context stays separate");
		check(count(sourceResults, "\tSAME_BOUNDS_QUALIFYING_GROUP\t")==3, "One duplicate-collapsed locus plus two independent tandem loci");
		check(count(sourceResults, "\tSAME_BOUNDS_PROTECTED_PRIOR\t")==1
			&& count(sourceResults, "\tCHANGED_BOUNDS_PROTECTED_PRIOR\t")==1, "Prior-bound same/changed pairs stay separate from CM qualifying-group pairs");
		final String leftResult=find(sourceResults, "\t3000\t3099\t+");
		final String rightResult=find(sourceResults, "\t3200\t3299\t+");
		check(!leftResult.split("\t")[8].equals(rightResult.split("\t")[8]), "A tandem bridge must not join the two direct correspondence components");
		check(count(lines("result.paths.tsv.gz"), "\tTANDEM_BRIDGE_WITNESS_NOT_REPRESENTATIVE")==1,
			"The tandem bridge survives as exactly one diagnostic path");
		check(count(sourceResults, "\tONLY_WEAK_OR_BRIDGE_SCORED_PATH\t")==1,
			"A source with only weak scored-parent correspondence is neither kept nor removed");
		check(find(sourceResults, "\t3485\t3509\t+").contains("\tUNASSIGNED_SOURCE\t"),
			"Unknown-family source remains unassigned even with only weak scored correspondence");
		check(find(sourceResults, "\t3540\t3659\t+").contains("\tMASKED_NO_SCORED_PARENT\t"),
			"A direct mask takes precedence over weak scored correspondence");
		check(find(lines("result.summary.tsv"), "FMN\t").startsWith("FMN\t2\t19\t18\t1\t1\t1\t"),
			"Selected claims, physical loci, context rows and context links conserve their separate populations");
		check(count(sourceResults, "\tCHANGED_BOUNDS_QUALIFYING_GROUP\t")==2, "Two unique changed CM-group pairs, excluding changed prior bounds");
		check(count(sourceResults, "\tSPLIT\t")==1 && count(sourceResults, "\tMERGE\t")==2
			&& count(sourceResults, "\tCOMPLEX\t")==2, "Graph components, not arbitrary pairing, define split/merge/complex");
		check(find(sourceResults, "\t280\t419\t-").endsWith("\t20\t-20\t40"), "Minus-strand end orientation");
		check(find(sourceResults, "\t450\t649\t+").endsWith("\t50\t-140\t190"), "Short CM core remains linked and measures full source shrinkage");
		check(count(sourceResults, "\tMASKED_NO_SCORED_PARENT\t")==2 && count(sourceResults, "\tREJECTED_NO_SCORED_PARENT\t")==1
			&& count(sourceResults, "\tUNRESOLVED_NO_SCORED_PARENT\t")==1, "Mask, rejection and missing correspondence are distinct");
		final ArrayList<String> parentResults=lines("result.parents.tsv.gz");
		check(count(parentResults, "\tCM_PARENT_WITHOUT_SOURCE_LINK\t")==1, "Only the unlinked CM parent is added");
		check(count(parentResults, "\tWEAK_OR_BRIDGE_SOURCE_LINK_ONLY\t")==1,
			"Diagnostic-only source correspondence is not a CM-added claim");
		check(count(parentResults, "\tHITLESS_PROTECTED_PRIOR\t")==1 && count(parentResults, "\tUNLINKED_PROTECTED_PRIOR\t")==1
			&& count(parentResults, "\tUNASSIGNED_SOURCE_LINK_ONLY\t")==1, "Priors and unassigned-source links are never CM-added");
		checkShifts("5prime", new double[]{20, 35, 50, 50, 35, 50, 50});
		checkShifts("3prime", new double[]{-140, -80, -20, -20, 80, 140, 140});
		expectFailure("result", null); // Fresh-output guard.
		expectFailure("sample", "maxinputs=1");
		// A missing witness fails before any successful summary can be published.
		write("broken_reference.tsv", "input_id\tfasta\tgff\n2\t"+path("unopened2.fa")+"\t"+path("empty.gff")+"\n");
		boolean failed=false;
		try{BuildNcrnaReference.main(new String[]{"mode=sourcecomparison", "in="+root, "reference="+path("broken_reference.tsv"), "out="+path("broken")});}
		catch(IllegalArgumentException expected){failed=true;}
		check(failed && !new File(path("broken.summary.tsv")).exists(), "Unknown saved input must fail without publishing completion");
		System.err.println("NcrnaSourceComparisonTest PASS");
	}

	private Row parent(final int start, final int stop, final char strand, final boolean mask){
		final Row row=codec.create("1", new Locus("contig_α%", "FMN", start, stop, strand), mask, mask ? "fixture" : null, "RF00050", path("annotation_witnesses.tsv.gz"));
		parents.add(row);
		return row;
	}
	private int event(final Row row, final String tier){
		return event(row, tier, row.locus.start, row.locus.stop,
			tier.equals("COVERED_BY_PROTECTED_PRIOR") ? "SAME_FAMILY_OVERLAP_NOT_RESCORING_PRIOR_BOUNDS" : "QUALIFYING_SAME_FAMILY_GROUP_MEMBER");
	}
	private int event(final Row row, final String tier, final int start, final int stop, final String association){
		final int id=++eventId;
		events.append(id).append("\t1\tcontig_α%25\t".getBytes(StandardCharsets.UTF_8)).append(start).tab().append(stop).tab().append(row.locus.strand)
			.append("\tRF00050\tFMN\t1\t100\t100\t0\tfalse\tfalse\t50\t40\t125\ttrue\tfalse\t80\ttrue\tfalse\t")
			.append(tier).append("\tfixture\tNONE\t0\t/native.tbl\t").append(id+10).append("\tPRINTED_0.1_BIT_NOT_INTERNAL_PRECISION\t6_DECIMAL_PLACES_SOURCE_PAIR_RETAINED\n");
		BuildNcrnaReference.field(witnesses, row.id).append("\t1\tCM_HIT\t").append(id).append("\t/native.tbl\t").append(id+10)
			.append("\tRF00050\t").append(start).tab().append(stop).tab()
			.append(tier.equals("COVERED_BY_PROTECTED_PRIOR") ? 9 : -1).tab().append(association).nl();
		return id;
	}
	private void prior(final Row row){
		BuildNcrnaReference.field(witnesses, row.id).append("\t1\tPROTECTED_PRIOR\t.\t.\t.\tRF00050\t").append(row.locus.start).tab()
			.append(row.locus.stop).append("\t9\tORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY\n");
	}
	private int source(final int start, final int stop, final char strand, final String family, final int disposition){
		final int line=++sourceLine;
		sourceCounts[family.equals("FMN") ? 0 : 1][disposition]++;
		sources.append("1\t/source.gff\t").append(line).append("\tcontig_α%25\tncRNA\t".getBytes(StandardCharsets.UTF_8)).append(start).tab().append(stop).tab().append(strand).tab()
			.append(family).tab().append(disposition==0 ? 0 : 1).tab().append(STATES[disposition]).append("\tNONE\tSOURCE_CLAIM_NOT_INDEPENDENT_BIOLOGICAL_TRUTH\n");
		return line;
	}
	private void link(final int source, final int event, final String agreement, final String tier){
		links.append("1\t/source.gff\t").append(source).tab().append(event).append("\ttrue\t").append(agreement).tab().append(tier).nl();
	}
	private void compare(final String prefix){
		BuildNcrnaReference.main(new String[]{"mode=sourcecomparison", "in="+root, "reference="+path("reference.tsv"), "out="+path(prefix)});
	}
	private void expectFailure(final String prefix, final String extra){
		boolean failed=false;
		try{
			if(extra==null){compare(prefix);}
			else{BuildNcrnaReference.main(new String[]{"mode=sourcecomparison", "in="+root, "reference="+path("reference.tsv"), "out="+path(prefix), extra});}
		}catch(IllegalArgumentException expected){failed=true;}
		catch(AssertionError expected){
			// Tools.testOutputFiles uses this assertion under the required -ea.
			check(extra==null && ("File "+path(prefix+".source_rows.tsv.gz")+" exists and overwrite=false").equals(expected.getMessage()),
				"Unexpected assertion while checking rejected output/sampling: "+expected.getMessage());
			failed=true;
		}
		check(failed, "Expected rejection of overwrite or sampling");
	}
	/** Native numeric formatting may omit .0; assert values, not decimal padding. */
	private void checkShifts(final String end, final double[] expected){
		final LineParser1 parser=new LineParser1('\t');
		parser.set(find(lines("result.shifts.tsv"), "FMN\t"+end+"\t").getBytes(StandardCharsets.UTF_8));
		check(parser.terms()==11 && parser.parseInt(2)==2 && expected.length==7,
			"Shift statistics require two changed qualifying-group pairs and seven numeric summaries");
		for(int i=0; i<expected.length; i++){
			check(Double.compare(parser.parseDouble(i+3), expected[i])==0,
				"Incorrect "+end+" shift statistic at field "+(i+3));
		}
	}
	private String path(final String name){return new File(root, name).getPath();}
	private void write(final String name, final String text){
		try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private void write(final String name, final ByteBuilder bytes){
		try(Output out=new Output(path(name))){out.writer.print(bytes);}
	}
	private ArrayList<String> lines(final String name){
		final ArrayList<String> lines=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){lines.add(new String(line, StandardCharsets.UTF_8));}
		}
		return lines;
	}
	private static String find(final ArrayList<String> lines, final String term){
		for(final String line:lines){if(line.contains(term)){return line;}}
		throw new AssertionError("Missing expected row: "+term);
	}
	private static int count(final ArrayList<String> lines, final String term){
		int n=0;
		for(final String line:lines){if(line.contains(term)){n++;}}
		return n;
	}
	private static void check(final boolean value, final String message){if(!value){throw new AssertionError(message);}}

	private final File root;
	private final NcrnaReferenceGff codec=new NcrnaReferenceGff();
	private final ArrayList<Row> parents=new ArrayList<Row>();
	private int eventId, sourceLine;
	private final long[][] sourceCounts=new long[2][11];
	private static final String[] STATES={"CONTEXT_PROVENANCE", "UNASSESSABLE_SOURCE", "PROTECTED_PRIOR_SOURCE", "SUPPORTED_CONFIDENT_EVENT", "SUPPORTED_MARGINAL_EVENT", "UNCERTAINTY_UNMASKED", "REMOVED_REVIEWED_FAMILY_CONFLICT", "PROVENANCE_NEAR_PROTECTED_PRIOR", "REMOVED_BELOW_MARGINAL", "REMOVED_NO_CM_SUPPORT", "REMOVED_OTHER_FAMILY_ONLY"};
	private final ByteBuilder sources=new ByteBuilder("input_id\tgff\tsource_line\ttarget\tsource_type\tstart\tstop\tstrand\tdeclared_family\tselection_bits\tdisposition\tsource_flags\tqualification\n");
	private final ByteBuilder links=new ByteBuilder("input_id\tgff\tsource_line\tevent_id\tsame_strand_valid\tfamily_agreement\tevent_tier\n");
	private final ByteBuilder witnesses=new ByteBuilder("annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation\n");
	private final ByteBuilder events=new ByteBuilder("event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision\n");
}
