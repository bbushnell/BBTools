package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaConstructionEvidence;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal saved-evidence joins: namespace separation, disjoint cores and loud
 * corruption failures. No training/admission conclusion follows these fixtures.
 * @author Ganyu
 */
public final class NcrnaConstructionEvidenceTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionEvidenceTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown evidence fixture parameter: "+arg);
			}
			boolean assertions=false;
			assert(assertions=true);
			check(assertions && parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use -ea and out=fresh_fixture_directory without sampling/overwrite overrides");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture directory must be fresh");
			run(root.getPath());
			pp.outstream.println("NCRNA_CONSTRUCTION_EVIDENCE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final String root){
		final Fixture good=new Fixture(root, "good");
		good.run();
		checkResult(good.dir+"/result");
		NcrnaConstructionEvidence.main(new String[]{"bindings="+good.dir+"/original", "evaluation="+good.dir+"/evaluation",
			"discovery="+good.dir+"/discovery", "expectedmembers=3", "out="+good.dir+"/cli", "t=1"});
		checkResult(good.dir+"/cli");
		final Fixture mixed=new Fixture(root, "mixed");
		mixed.run();
		check(rows(mixed.dir+"/result.old_evidence.tsv.gz").get(1).equals(
			"evaluation\tRF00050_0_0\t1\t2\t2\t2\t1\t1\t10\t20\tMIXED_FULL_MODEL_EVIDENCE\tNOT_EVALUATED"),
			"Conflicting full-model evidence must remain visible without a deletion/admission decision");
		final Fixture coordinate=new Fixture(root, "coordinate_only");
		coordinate.run();
		check(rows(coordinate.dir+"/result.old_evidence.tsv.gz").get(3).equals(
			"discovery\tRF00050_2_0\t2\t1\t0\t2\t1\t1\t4\t20\tMISSING_SOURCE_WITNESS\tNOT_EVALUATED"),
			"Missing saved GFF rows cannot hide native same-family or competing-family coordinate support");
		final ArrayList<String> coordinateLinks=rows(coordinate.dir+"/result.links.tsv.gz");
		check(coordinateLinks.size()==7
			&& coordinateLinks.get(5).startsWith("discovery\tRF00050_2_0\t20\t20\tBOUND_CORE_COORDINATE\t5\t2\tcontig\t51\t80\t+")
			&& coordinateLinks.get(6).startsWith("discovery\tRF00050_2_0\t20\t4\tBOUND_CORE_COORDINATE\t6\t2\tcontig\t55\t58\t+"),
			"Coordinate-only links retain exact input identity and an explicit basis, never invented source witnesses");
		check(rows(coordinate.dir+"/result.fmn_events.tsv.gz").get(6).startsWith("discovery\ttrue\t7\t2\tcontig\t61\t70\t+"),
			"An FMN event wholly in the multipart gap must remain in native events without becoming a member link");
		final Fixture partial=new Fixture(root, "partial_witness");
		partial.run();
		check(rows(partial.dir+"/result.old_evidence.tsv.gz").get(1).equals(
			"evaluation\tRF00050_0_0\t1\t2\t1\t1\t1\t0\t0\t20\tMISSING_SOURCE_WITNESS\tNOT_EVALUATED")
			&& rows(partial.dir+"/result.links.tsv.gz").get(1).startsWith(
				"evaluation\tRF00050_0_0\t20\t20\tSOURCE_LINE_AND_BOUND_CORE\t1\t"),
			"Partial witness loss must union source and coordinate support without double-counting the same member/event");
		final String[] faults={"missing_event", "truncated_events", "bad_source", "wrong_input", "wrong_contig", "wrong_strand",
			"duplicate_source", "dangling_source", "bad_boolean", "invalid_interval", "bad_member_count", "bad_intervals", "bad_line"};
		final String[] reasons={"absent native event", "conserve", "no longer overlaps", "retain their input", "retain their input", "retain their input",
			"witness is duplicated", "missing original source", "literal true/false", "invalid native interval", "declared complete population",
			"disjoint", "positive"};
		check(faults.length==reasons.length, "Every malformed fixture requires a named failure");
		for(int i=0; i<faults.length; i++){
			final Fixture bad=new Fixture(root, faults[i]);
			try{bad.run();}
			catch(IllegalArgumentException error){
				check(error.getMessage().contains(reasons[i]), "Wrong rejection for "+faults[i]+": "+error.getMessage());
				continue;
			}
			throw new AssertionError("Malformed evidence was accepted: "+faults[i]);
		}
	}

	private static void checkResult(final String prefix){
		final ArrayList<String> old=rows(prefix+".old_evidence.tsv.gz");
		check(old.size()==4, "Every bound original member, including missing sources, must be retained");
		check(old.get(1).equals("evaluation\tRF00050_0_0\t1\t2\t2\t1\t1\t0\t0\t20\tOWN_FULL_FMN\tNOT_EVALUATED"),
			"Two source witnesses and duplicate links must yield one member/event association");
		check(old.get(2).equals("discovery\tRF00050_1_0\t1\t1\t1\t3\t0\t1\t4\t20\tOTHER_FULL_MODEL_EVIDENCE\tNOT_EVALUATED"),
			"Same numeric input/event/source IDs in another cohort must remain separate; a core gap contributes zero overlap");
		check(old.get(3).equals("discovery\tRF00050_2_0\t2\t1\t0\t0\t0\t0\t0\t20\tMISSING_SOURCE_WITNESS\tNOT_EVALUATED"),
			"Missing source evidence must remain unresolved, never silently dropped or labelled a negative");
		final ArrayList<String> summary=rows(prefix+".summary.tsv");
		check(summary.size()==3 && summary.get(1).equals("evaluation\t1\t0\t1\t2\t1\t1\tEVIDENCE_JOIN_NOT_ADMISSION_OR_ALIAS_PROOF")
			&& summary.get(2).equals("discovery\t2\t1\t3\t4\t3\t2\tEVIDENCE_JOIN_NOT_ADMISSION_OR_ALIAS_PROOF"),
			"Both complete event populations and binding totals must be conserved");
		final ArrayList<String> links=rows(prefix+".links.tsv.gz");
		check(links.size()==5 && links.get(2).startsWith("discovery\tRF00050_1_0\t20\t4\tSOURCE_LINE\t1\t1\tcontig\t15\t18\t+\tRF00005\ttRNA\t")
			&& links.get(4).startsWith("discovery\tRF00050_1_0\t20\t0\tSOURCE_LINE\t3\t"),
			"Links retain literal native identities and measured disjoint-core overlaps");
		final ArrayList<String> fmn=rows(prefix+".fmn_events.tsv.gz");
		check(fmn.size()==5 && fmn.get(2).startsWith("discovery\tfalse\t2\t") && fmn.get(4).contains("\t-\tRF00050\tFMN\t"),
			"All FMN events survive, including partial hits and events on the opposite strand without a valid source link");
	}

	private static final class Fixture{
		Fixture(final String root, final String fault_){
			fault=fault_;
			dir=root+'/'+fault;
			check(new File(dir).mkdir() && new File(dir+"/evaluation").mkdir() && new File(dir+"/discovery").mkdir(),
				"Each case requires fresh cohort directories");
			write("original.members.tsv.gz", NcrnaConstructionEvidence.MEMBER_HEADER+'\n'
				+member(0, "eval", "11-30", fault.equals("bad_line") ? "0" : "7,8")
				+member(1, "train", fault.equals("bad_intervals") ? "11-20,20-29" : "11-20,31-40", "7")
				+member(2, "missing", fault.equals("coordinate_only") ? "51-60,71-80" : "51-70", "9"));
			write("original.bindings.tsv", NcrnaConstructionEvidence.BINDING_HEADER+'\n'
				+binding("eval", "evaluation", 1)+binding("train", "discovery", 1)+binding("missing", "discovery", 2));
			final String source7=source(7, fault.equals("bad_source") ? 71 : 11, fault.equals("bad_source") ? 90 : 30);
			write("evaluation/source_provenance.tsv.gz", NcrnaConstructionEvidence.SOURCE_HEADER+'\n'+source7
				+(fault.equals("partial_witness") ? "" : source(8, 11, 30))
				+(fault.equals("duplicate_source") ? source7 : ""));
			write("discovery/source_provenance.tsv.gz", NcrnaConstructionEvidence.SOURCE_HEADER+'\n'+source(7, 11, 40));
			write("evaluation/source_event_links.tsv.gz", NcrnaConstructionEvidence.LINK_HEADER+'\n'
				+link(7, 1, "true")+link(7, 1, "true")+(fault.equals("partial_witness") ? "" : link(8, 1, "true"))
				+link(7, 2, fault.equals("mixed") ? "true" : "false")
				+(fault.equals("missing_event") ? link(7, 99, "true") : "")
				+(fault.equals("bad_boolean") ? link(7, 2, "yes") : ""));
			write("discovery/source_event_links.tsv.gz", NcrnaConstructionEvidence.LINK_HEADER+'\n'
				+link(7, 1, "true")+link(7, 2, "true")+link(7, 3, "true")+link(7, 4, "false")
				+(fault.equals("dangling_source") ? "2\t/original.gff\t9\t4\ttrue\tSAME_FAMILY\tCONFIDENT\n" : ""));
			String first=event(1, "RF00050", "FMN", 11, 30, '+', 20, 0, true);
			if(fault.equals("wrong_input")){first=first.replace("1\t1\tcontig", "1\t2\tcontig");}
			if(fault.equals("wrong_contig")){first=first.replace("contig", "another");}
			if(fault.equals("wrong_strand")){first=first.replace("\t+\t", "\t-\t");}
			if(fault.equals("invalid_interval")){first=first.replace("CONFIDENT", "INVALID_INTERVAL");}
			write("evaluation/event_provenance.tsv.gz", NcrnaConstructionEvidence.EVENT_HEADER+'\n'+first
				+(fault.equals("truncated_events") ? "" : event(2, "RF00005", "tRNA", 16, 25, fault.equals("mixed") ? '+' : '-', 20, 0, true)));
			write("discovery/event_provenance.tsv.gz", NcrnaConstructionEvidence.EVENT_HEADER+'\n'
				+event(1, "RF00005", "tRNA", 15, 18, '+', 20, 0, true)
				+event(2, "RF00050", "FMN", 11, 20, '+', 10, 1, false)
				+event(3, "RF00050", "FMN", 21, 30, '+', 20, 0, true)
				+event(4, "RF00050", "FMN", 11, 30, '-', 20, 0, true)
				+(fault.equals("coordinate_only") ? event(5, "RF00050", "FMN", 51, 80, '+', 20, 0, true).replace("5\t1\tcontig", "5\t2\tcontig")
					+event(6, "RF00005", "tRNA", 55, 58, '+', 20, 0, true).replace("6\t1\tcontig", "6\t2\tcontig")
					+event(7, "RF00050", "FMN", 61, 70, '+', 20, 0, true).replace("7\t1\tcontig", "7\t2\tcontig") : ""));
			write("evaluation/family_counts.tsv", NcrnaConstructionEvidence.FAMILY_HEADER+'\n'+family("FMN", 1)+family("tRNA", 1));
			write("discovery/family_counts.tsv", NcrnaConstructionEvidence.FAMILY_HEADER+'\n'
				+family("FMN", fault.equals("coordinate_only") ? 5 : 3)+family("tRNA", fault.equals("coordinate_only") ? 2 : 1));
		}
		void run(){
			NcrnaConstructionEvidence.process(dir+"/original", dir+"/evaluation", dir+"/discovery", dir+"/result",
				fault.equals("bad_member_count") ? 4 : 3);
		}
		void write(final String path, final String text){
			try(Output out=new Output(dir+'/'+path)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		}
		final String dir, fault;
	}

	private static String member(final int index, final String genome, final String intervals, final String lines){
		return "true\tfalse\tNO_RECORDED_EVALUATION_IDENTITY\tUNVERIFIED\tNOT_EVALUATED\tRF00050_"+index+"_0\t"+genome
			+"\tGCF_TEST.1\t1\tbacteria\tP\tcontig\t+\t"+intervals+"\t20\t7\tfeature\tfalse\tfalse\tANNOTATION_FULL_COPY\t100\tfalse\t0\t0\t0\t"
			+lines+"\tfull header\t1\tunknown\tcluster\t\torphan\n";
	}
	private static String binding(final String genome, final String cohort, final int input){
		return genome+'\t'+"/"+genome+".fa\t/original.gff\t"+cohort+'\t'+input+"\tLITERAL_CATALOG_PATH_NOT_SEQUENCE_PROOF\n";
	}
	private static String source(final int line, final int start, final int stop){
		return "1\t/original.gff\t"+line+"\tcontig\tncRNA\t"+start+'\t'+stop
			+"\t+\tFMN\t1\tCONFIDENT\tNONE\tSOURCE_CLAIM_NOT_INDEPENDENT_BIOLOGICAL_TRUTH\n";
	}
	private static String link(final int line, final int event, final String valid){
		return "1\t/original.gff\t"+line+'\t'+event+'\t'+valid+"\tSAME_FAMILY\tCONFIDENT\n";
	}
	private static String event(final int id, final String model, final String family, final int start, final int stop,
			final char strand, final int modelTo, final int truncation, final boolean included){
		return id+"\t1\tcontig\t"+start+'\t'+stop+'\t'+strand+'\t'+model+'\t'+family+"\t1\t"+modelTo
			+"\t20\t"+truncation+"\tfalse\tfalse\t50.0\t40\t125.000000\t"+included
			+"\tfalse\t100\ttrue\tfalse\tCONFIDENT\tFIXTURE\tNONE\t0\t/native.tbl\t"+id
			+"\tPRINTED_0.1_BIT_NOT_INTERNAL_PRECISION\t6_DECIMAL_PLACES_SOURCE_PAIR_RETAINED\n";
	}
	private static String family(final String name, final int count){return name+'\t'+count+"\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";}
	private static ArrayList<String> rows(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){result.add(new String(row, StandardCharsets.UTF_8));}
		}
		check(!result.isEmpty(), "Fixture result must include its header: "+path);
		return result;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
