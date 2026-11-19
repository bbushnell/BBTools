package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaConstructionMembers;
import prok.NcrnaDdlLocusOverlap;
import prok.NcrnaDdlLocusOverlap.Config;
import prok.NcrnaDdlMembers;
import prok.NcrnaMemberLengthBounds;
import prok.NcrnaSavedRankSummary;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal saved-window, corrected-source and actual-kept-sequence joins.
 * Every genome path is nonexistent; all sequence evidence is saved synthetic FASTA.
 * @author Raiden
 */
public final class NcrnaDdlLocusOverlapTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlLocusOverlapTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown locus-overlap fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Supply out=fresh_fixture_directory without sampling or overwrite");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Locus-overlap fixture directory must be fresh");
			Shared.TRIM_READ_DESCRIPTION=false;
			NcrnaExactSequenceSetTest.run();
			geometry();
			new NcrnaDdlLocusOverlapTest(root.getPath()).run();
			pp.outstream.println("NCRNA_DDL_LOCUS_OVERLAP_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaDdlLocusOverlapTest(final String root_){root=root_;}

	private static void geometry(){
		final byte[] forward=dna("AGTCGA"), reverse=dna("TCGACT");
		check(Arrays.equals(NcrnaDdlLocusOverlap.extract(forward, 10, 15, 0, 12, 14, 0), dna("GTC")), "Inclusive forward extraction uses one-based query bounds");
		check(Arrays.equals(NcrnaDdlLocusOverlap.extract(forward, 10, 15, 0, 12, 14, 1), dna("GAC")), "Minus query is reverse complemented from a forward capture");
		check(Arrays.equals(NcrnaDdlLocusOverlap.extract(reverse, 10, 15, 1, 12, 14, 0), dna("GTC")), "Forward query is reconstructed from a reverse capture");
		check(Arrays.equals(NcrnaDdlLocusOverlap.extract(reverse, 10, 15, 1, 12, 14, 1), dna("GAC")), "Minus query preserves a reverse capture's orientation");
		check(Arrays.equals(NcrnaDdlLocusOverlap.extract(forward, 10, 15, 0, 11, 16, 0), forward), "Both inclusive window endpoints belong to the captured sequence");
		check(NcrnaDdlLocusOverlap.extract(forward, 10, 15, 0, 10, 16, 0)==null
			&& NcrnaDdlLocusOverlap.extract(forward, 10, 15, 0, 11, 17, 0)==null, "One missing base at either end is incomplete query sequence");
		check(NcrnaDdlLocusOverlap.extract(dna("AGNCGA"), 10, 15, 0, 12, 14, 0)==null, "Undefined query bases retain unknown exact-sequence status");
		check(Arrays.equals(NcrnaDdlLocusOverlap.extract(dna("NGTCGA"), 10, 15, 0, 12, 14, 0), dna("GTC")), "Undefined flanking bases do not invalidate a complete defined gene");
		check(!NcrnaDdlLocusOverlap.sameLocus(10, 19, 0, 15, 24, 0)
			&& NcrnaDdlLocusOverlap.sameLocus(10, 19, 0, 14, 23, 0), "Same-locus evidence needs strictly greater than half overlap of both spans");
		check(!NcrnaDdlLocusOverlap.sameLocus(10, 19, 0, 10, 19, 1), "Opposite-strand intervals are not the same source locus");
		final int[][] expected={{0,1,2}, {1,1,1}, {2,1,2}};
		for(int a=0; a<3; a++){
			for(int b=0; b<3; b++){
				check(NcrnaDdlLocusOverlap.classify(a, b)==expected[a][b], "Positive evidence wins; negative status needs both complete negative components");
			}
		}
	}

	private void run(){
		final String rawA="RF00174_a source=gA status=ANNOTATION_FULL_COPY";
		final String rawB="RF00174_b source=gA contig=ctgA strand=+ start=500 stop=519 historical_start=100 historical_stop=119 recovery=CM_FULL_SPAN_SUPPORTED";
		final String rawC="RF01482_c source=gB source_family=RF01482 member_family=B12clade contig=ctgB strand=- start=700 stop=719 coordinate_model=RF00174";
		final String rawO="RF00174_o source=gC status=ANNOTATION_FULL_COPY";
		final String rawS="RF00174_s source=gC status=ANNOTATION_FULL_COPY";
		final String a=rawA+" member_family=RF00174 build_group=RF00174";
		final String b=rawB+" member_family=RF00174 build_group=RF00174";
		final String c=rawC+" build_group=RF00174";
		final String orphan=rawO+" member_family=RF00174 build_group=RF00174";
		final String small=rawS+" member_family=RF00174 build_group=RF00174";
		write("members.fa", fasta(a, A)+fasta(b, C)+fasta(c, G)+fasta(orphan, T)+fasta(small, "ACGTACGTACGTACGTACGT"));
		write("references.fa", fasta(MODEL+" n=3", A));
		write("ledger.tsv", "locus_id\tsource_family\tmember_family\tbuild_group\trouting_rule\trouting_bits\tlength\tinput_header\n"
			+ledger("RF00174_a", "RF00174", "RF00174", rawA)+ledger("RF00174_b", "RF00174", "RF00174", rawB)
			+ledger("RF01482_c", "RF01482", "B12clade", rawC)+ledger("RF00174_o", "RF00174", "RF00174", rawO)
			+ledger("RF00174_s", "RF00174", "RF00174", rawS));
		write("assignments.tsv", "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n"
			+assignment(a, MODEL, "kept")+assignment(b, MODEL, "kept")+assignment(c, MODEL, "kept")
			+assignment(orphan, "", "orphan")+assignment(small, "", "small"));
		write("census.tsv", NcrnaMemberLengthBounds.HEADER+'\n'+MODEL+"\tRF00174\t20\t3\t20\t20\t20\t20\t20\t0.6666667\tMIXED\t2\t0.6666667\t0\t0\t0\t0\t1\t0.3333333\n");
		write("summary.tsv", "emitted_family\tinput_members\toriginal_kept\toriginal_small\toriginal_orphan\tadded_memberships\tconstruction_memberships\tmodels\tpure\tmixed\tempty\tpurity_cut\tgroup_min\tgroup_p10\tgroup_median\tgroup_p90\tgroup_max\n"
			+"RF00174\t5\t3\t1\t1\t0\t3\t1\t0\t1\t0\t0.8\t20\t20\t20\t20\t20\n");
		String genomes=NcrnaConstructionMembers.GENOME_HEADER+'\n';
		String groups=NcrnaConstructionMembers.GROUP_HEADER+'\n';
		String assemblies=NcrnaDdlLocusOverlap.ASSEMBLY_HEADER+'\n';
		for(int i=0; i<3; i++){
			genomes+=i+"\tg"+(char)('A'+i)+"\tasm"+i+"\t7\tbacteria\tP\t/never/g.gff\t"+PATHS[i]+"\t1\t1000\t1\n";
			groups+="evaluation\t"+(i+1)+'\t'+PATHS[i]+'\t'+PATHS[i]+"\trecord\t"+i+"\ttrue\tfalse\tfalse\tACCEPTED\t"+i+"\ttrue\tfalse\tACCEPTED\n";
			assemblies+=i+"\t"+(i+1)+'\t'+PATHS[i]+'\t'+i+'\t'+(i==0 ? 2 : i==1 ? 1 : 0)+"\tSTATUS\tASSEMBLY_ONLY\n";
		}
		write("genomes.tsv", genomes);
		write("groups.tsv", groups);
		write("assemblies.tsv", assemblies);
		final String rawBRow=raw("RF00174_b", "gA", "ctgA", "+", "100-119", 20);
		write("proxy.tsv.gz", NcrnaDdlLocusOverlap.RAW_HEADER+'\n'+raw("RF00174_a", "gA", "ctgA", "+", "11-30", 20)
			+rawBRow+raw("RF00174_o", "gC", "ctgC", "+", "101-120", 20));
		final String recovered=NcrnaDdlLocusOverlap.RECOVERED_HEADER+'\n'+rawBRow.trim()
			+"\tCM_FULL_SPAN_SUPPORTED\t1\t1\t1\t20\t1\t20\tno\t50\t20\t0\t20\tHIT_WITHIN_ANNOTATION\tRECOVERED_CM_CORE\t500\t519\t20\n";
		write("recovered.tsv", recovered);
		final String overlay=NcrnaDdlLocusOverlap.OVERLAY_HEADER+'\n'+raw("RF01482_c", "gB", "ctgB", "-", "100-119", 20).trim()
			+"\tRF01482\tRF00174\ttrue\tRF00174\tRF00174\tCURRENT\t700\t719\t20\tWITHIN_BAND\tB12clade\tCORRECTED\n";
		write("overlay.tsv.gz", overlay);
		final int[] slots={0,0,0,1,2,2,2};
		final int[] starts={11,500,100,700,11,101,201};
		final int[] stops={30,519,119,719,30,120,220};
		final String[] sequences={A,T,T,G,A,T,T.substring(1),A,A};
		String parents=NcrnaSavedRankSummary.PARENT_HEADER+'\n';
		String loci=NcrnaSavedRankSummary.LOCUS_HEADER+'\n';
		String windows=NcrnaDdlLocusOverlap.WINDOW_HEADER+'\n';
		String captured="", forwardCaptured="";
		for(int i=0; i<slots.length; i++){
			final int slot=slots[i];
			parents+=slot+"\t"+(slot+1)+'\t'+PATHS[slot]+"\tp"+i+"\tctg"+(char)('A'+slot)+'\t'+starts[i]+'\t'+stops[i]+'\t'+(i==3 ? "-" : "+")+"\tfalse\ttrue\tfalse\tfalse\tRF00174\t.\n";
			for(int method=0; method<2; method++){
				final int count=i==0 ? 3 : 1;
				loci+=slot+"\t"+(slot+1)+"\tp"+i+'\t'+method+'\t'+count+'\t'+count+"\t0\t"+count+"\t1\tRANKED\t1\t1\n";
			}
		}
		for(int event=0; event<9; event++){
			final int parent=event<7 ? event : 0;
			final int slot=slots[parent], start=starts[parent]-1, stop=stops[parent]-1-(parent==6 ? 1 : 0), strand=parent==3 ? 1 : 0;
			windows+=slot+"\t"+(slot+1)+'\t'+event+"\t1\tctg"+(char)('A'+slot)+"\theader\t"+strand+"\t0\t"
				+(strand==0 ? start : 999-stop)+'\t'+(strand==0 ? stop : 999-start)+'\t'+start+'\t'+stop+"\t1000\tSAVED\t"+(stop-start+1)+'\n';
			final String record=fasta(slot+"_"+event, sequences[event]);
			captured=record+captured;
			forwardCaptured+=record;
		}
		writeRank("rank", parents, loci, windows, captured);
		final Config config=config("result", "rank");
		NcrnaDdlLocusOverlap.process(config);
		final ArrayList<String> labels=rows("result.labels.tsv");
		check(labels.size()==8 && labels.get(1).equals("0\t1\tp0\tSEEN") && labels.get(2).equals("0\t1\tp1\tSEEN")
			&& labels.get(3).equals("0\t1\tp2\tNOT_SEEN") && labels.get(4).equals("1\t2\tp3\tSEEN")
			&& labels.get(5).equals("2\t3\tp4\tSEEN") && labels.get(6).equals("2\t3\tp5\tNOT_SEEN")
			&& labels.get(7).equals("2\t3\tp6\tUNKNOWN"), "Literal source/sequence membership includes cross-assembly sequence presence and retains partial capture as unknown");
		final ArrayList<String> result=rows("result.summary.tsv");
		check(result.get(1).equals("NOT_SEEN\t2\t2\t0\t0\t0\t0\t0\t2")
			&& result.get(2).equals("SEEN\t4\t4\t3\t3\t2\t0\t0\t4")
			&& result.get(3).equals("UNKNOWN\t1\t1\t0\t0\t0\t1\t0\t0"), "Exact-sequence, corrected-source and final common-locus denominators are independently specified");
		final ArrayList<String> detail=rows("result.loci.tsv");
		check(detail.get(2).endsWith("\tNOT_SEEN\tSEEN\tSEEN\tRF00174_b"), "Recovered final bounds, not historical annotation bounds, identify the source locus");
		check(detail.get(4).endsWith("\tSEEN\tSEEN\tSEEN\tRF01482_c"), "Current admitted clan bounds and minus-strand capture preserve the actual cross-labelled source");
		NcrnaDdlLocusOverlap.main(new String[]{"in="+path("rank"), "out="+path("cli"), "members="+path("members.fa"),
			"ledger="+path("ledger.tsv"), "assignments="+path("assignments.tsv"), "census="+path("census.tsv"),
			"summary="+path("summary.tsv"), "ref="+path("references.fa"), "genomes="+path("genomes.tsv"),
			"groups="+path("groups.tsv"), "assemblies="+path("assemblies.tsv"), "proxy="+path("proxy.tsv.gz"),
			"recovered="+path("recovered.tsv"), "overlays="+path("overlay.tsv.gz"), "expectedinputs=3", "expectedgenomes=3",
			"expectedtargets=7", "expectedcommon=7", "expectedmodels=1", "expectedmethods=2", "t=1"});
		check(labels.equals(rows("cli.labels.tsv")) && result.equals(rows("cli.summary.tsv"))
			&& detail.equals(rows("cli.loci.tsv")), "Every locus-audit CLI argument reaches the same complete saved-data computation");
		NcrnaSavedRankSummary.process(path("rank"), path("rank_strata"), path("result.labels.tsv"), new int[]{0,1}, 2, 1, 7, 7, 0, 0);
		final ArrayList<String> ranks=rows("rank_strata.summary.tsv");
		check(ranks.size()==9 && ranks.get(3).startsWith("0\t4\t4\t1\t") && ranks.get(5).startsWith("0\t2\t2\t1\t")
			&& ranks.get(7).startsWith("0\t1\t1\t1\t"), "The actual saved-rank reporter consumes generated SEEN/NOT_SEEN/UNKNOWN labels and conserves all common loci");
		observer();
		write("recovered_missing.tsv", NcrnaDdlLocusOverlap.RECOVERED_HEADER+'\n');
		final Config missingCorrection=config("missing_correction", "rank");
		missingCorrection.recovered=path("recovered_missing.tsv");
		NcrnaDdlLocusOverlap.process(missingCorrection);
		check(rows("missing_correction.labels.tsv").get(2).equals("0\t1\tp1\tUNKNOWN")
			&& rows("missing_correction.labels.tsv").get(3).equals("0\t1\tp2\tUNKNOWN"),
			"A missing same-length shifted correction cannot fall back to historical source coordinates");
		write("recovered_shifted.tsv", recovered.replace("\t500\t519\t20\n", "\t600\t619\t20\n"));
		final Config shifted=config("shifted_correction", "rank");
		shifted.recovered=path("recovered_shifted.tsv");
		reject(shifted, "retained FASTA header coordinates");
		final String undefined=fasta("0_8", 'N'+A.substring(1));
		writeRank("undefined_first", parents, loci, windows, captured.replace(fasta("0_8", A), undefined));
		NcrnaDdlLocusOverlap.process(config("undefined_first_result", "undefined_first"));
		writeRank("undefined_last", parents, loci, windows, forwardCaptured.replace(fasta("0_8", A), undefined));
		NcrnaDdlLocusOverlap.process(config("undefined_last_result", "undefined_last"));
		check(result.equals(rows("undefined_first_result.summary.tsv")) && result.equals(rows("undefined_last_result.summary.tsv"))
			&& labels.equals(rows("undefined_first_result.labels.tsv")) && labels.equals(rows("undefined_last_result.labels.tsv")),
			"Defined agreeing captures provide the same sequence and labels before or after an undefined containing capture");
		writeRank("conflicting", parents, loci, windows, captured.replace(fasta("0_8", A), fasta("0_8", C)));
		reject(config("conflicting_result", "conflicting"), "Conflicting complete captured sequences");
		write("recovered_bad.tsv", recovered.replace("\t500\t519\t20\n", "\t500\t520\t20\n"));
		final Config bad=config("bad_coordinates", "rank");
		bad.recovered=path("recovered_bad.tsv");
		reject(bad, "Final corrected source coordinates");
		writeRank("duplicate_window", parents, loci, windows, captured+fasta("0_0", A));
		reject(config("bad_duplicate", "duplicate_window"), "Duplicate or truncated captured window");
		writeRank("missing_window", parents, loci, windows, captured.replace(fasta("0_0", A), ""));
		reject(config("bad_missing", "missing_window"), "without its FASTA record");
	}

	/** Check both the unchanged null observer API and loud callback failure paths. */
	private void observer(){
		final NcrnaDdlMembers plain=load(null);
		check(plain.keptMembers()==3 && plain.members(0).length==3, "The original loader API retains its ordinary member arrays");
		final ArrayList<String> seen=new ArrayList<String>();
		load((header, model, bases)->{seen.add(header); check(model==0 && bases.length==20, "Observer borrows complete checked members");});
		check(seen.size()==3 && seen.get(0).startsWith("RF00174_a ") && seen.get(2).startsWith("RF01482_c "), "Orphan and small records never reach the kept-member observer");
		final RuntimeException runtime=new IllegalStateException("literal observer runtime failure");
		try{load((header, model, bases)->{throw runtime;}); throw new AssertionError("Observer runtime failure was swallowed");}
		catch(RuntimeException actual){check(actual==runtime, "Draining preserves the original observer runtime failure");}
		final AssertionError error=new AssertionError("literal observer assertion failure");
		try{load((header, model, bases)->{throw error;}); throw new AssertionError("Observer assertion was swallowed");}
		catch(AssertionError actual){check(actual==error, "Draining preserves the original observer assertion failure");}
	}

	private NcrnaDdlMembers load(final NcrnaDdlMembers.KeptMemberSink sink){
		if(sink==null){return NcrnaDdlMembers.load(path("members.fa"), path("ledger.tsv"), path("assignments.tsv"), path("census.tsv"), path("summary.tsv"), new String[]{MODEL}, new byte[][]{dna(A)});}
		return NcrnaDdlMembers.load(path("members.fa"), path("ledger.tsv"), path("assignments.tsv"), path("census.tsv"), path("summary.tsv"), new String[]{MODEL}, new byte[][]{dna(A)}, sink);
	}

	private Config config(final String output, final String prefix){
		final Config c=new Config();
		c.input=path(prefix);
		c.output=path(output);
		c.members=path("members.fa");
		c.ledger=path("ledger.tsv");
		c.assignments=path("assignments.tsv");
		c.census=path("census.tsv");
		c.summary=path("summary.tsv");
		c.references=path("references.fa");
		c.genomes=path("genomes.tsv");
		c.groups=path("groups.tsv");
		c.assemblies=path("assemblies.tsv");
		c.proxy=path("proxy.tsv.gz");
		c.recovered=path("recovered.tsv");
		c.overlays=new String[]{path("overlay.tsv.gz")};
		c.inputs=3;
		c.genomeCount=3;
		c.targets=7;
		c.common=7;
		c.models=1;
		c.methods=2;
		return c;
	}

	private void reject(final Config config, final String reason){
		try{NcrnaDdlLocusOverlap.process(config);}
		catch(IllegalArgumentException expected){check(expected.getMessage().contains(reason), "Wrong rejection: "+expected.getMessage()); return;}
		throw new AssertionError("Malformed source/window input was accepted");
	}
	private void writeRank(final String prefix, final String parents, final String loci, final String windows, final String fasta){
		write(prefix+".parents.tsv.gz", parents);
		write(prefix+".loci.tsv.gz", loci);
		write(prefix+".windows.tsv.gz", windows);
		write(prefix+".windows.fa.gz", fasta);
	}
	private static String raw(final String id, final String source, final String seqid, final String strand, final String span, final int length){
		return id+'\t'+source+"\tasm\t7\tbacteria\tP\t"+seqid+'\t'+strand+'\t'+span+'\t'+length+"\t1\tids\tfalse\tfalse\tANNOTATION_FULL_COPY\t1000\tfalse\t0\t0\t0\t1\n";
	}
	private static String ledger(final String id, final String family, final String biological, final String raw){return id+'\t'+family+'\t'+biological+"\tRF00174\tRULE\tNA\t20\t"+raw+'\n';}
	private static String assignment(final String header, final String model, final String status){return header+"\t7\tunknown\tRF00174_fixture\t"+model+'\t'+status+'\n';}
	private static String fasta(final String header, final String sequence){return '>'+header+'\n'+sequence+'\n';}
	private static byte[] dna(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private String path(final String name){return root+'/'+name;}
	private void write(final String name, final String text){try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}}
	private ArrayList<String> rows(final String name){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path(name))){for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){result.add(new String(row, StandardCharsets.UTF_8));}}
		return result;
	}
	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
	private final String root;
	private static final String MODEL="cobalamin_RF00174_fixture";
	private static final String A="AAAAAAAAAAAAAAAAAAAA", C="CCCCCCCCCCCCCCCCCCCC", G="GGGGGGGGGGGGGGGGGGGG", T="TTTTTTTTTTTTTTTTTTTT";
	private static final String[] PATHS={"/never-opened/α.fa", "/never-opened/B.fa", "/never-opened/C.fa"};
}
