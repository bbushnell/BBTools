package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import parse.LineParser1;
import prok.NcrnaSavedRankSummary;
import prok.NcrnaDdlLocusOverlap;
import prok.NcrnaDdlTrainingOverlap;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal rank sums and joined named cases, plus incomplete/misbound-row rejection.
 * No sequence inputs or caller initialization.
 * @author Raiden
 */
public final class NcrnaSavedRankSummaryTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSavedRankSummaryTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown saved-rank fixture argument: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"The saved-rank fixture requires out=fresh_directory without sampling or overwrite");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture directory must be fresh");
			new NcrnaSavedRankSummaryTest(root.getPath()).run();
			pp.outstream.println("NCRNA_SAVED_RANK_SUMMARY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaSavedRankSummaryTest(final String root_){root=root_;}

	private void run(){
		String parents=NcrnaSavedRankSummary.PARENT_HEADER+'\n';
		for(int i=0; i<5; i++){
			if(i==1){parents+=parent("other", false);}
			parents+=parent("p"+i, true);
		}
		final int[][] ranks={{1,2,1,3}, {3,1,2,2}, {-1,4,2,1}, {-1,-1,-1,-1}, {-1,-1,-1,-1}};
		String loci=NcrnaSavedRankSummary.LOCUS_HEADER+'\n';
		for(int p=0; p<ranks.length; p++){
			for(int m=0; m<4; m++){
				final String reason=p==3 ? "NO_COMPATIBLE_MODEL" : p==4 ? "NO_CAPTURED_WINDOW"
					: p==2 && m==0 ? "ZERO_NATIVE_INDEX_SUPPORT" : "RANKED";
				loci+=row("p"+p, m, ranks[p][m], reason);
			}
		}
		write("input.parents.tsv.gz", parents);
		write("input.loci.tsv.gz", loci);
		final String strata=NcrnaSavedRankSummary.STRATA_HEADER+"\n0\tinput\tp0\tSEEN\n0\tinput\tp1\tNOT_SEEN\n"
			+"0\tinput\tp2\tUNKNOWN\n0\tinput\tp3\tUNKNOWN\n0\tinput\tp4\tUNKNOWN\n";
		write("strata.tsv", strata);
		process("result", 2);
		final ArrayList<String> summary=rows("result.summary.tsv");
		check(summary.size()==17, "One native and three selected methods retain ALL, SEEN, NOT_SEEN and UNKNOWN strata");
		checkNumbers(summary.get(1), 0, 2, 4, 2.0);
		check(summary.get(1).endsWith("ALL_COMMON\t2\t2\t2\t3\t3\t3"), "Coverage and nearest-rank quantiles use the same literal common set");
		checkNumbers(summary.get(2), 1, 2, 3, 1.5);
		checkNumbers(summary.get(3), 2, 2, 3, 1.5);
		checkNumbers(summary.get(4), 3, 2, 5, 2.5);
		checkNumbers(summary.get(5), 0, 1, 1, 1.0);
		checkNumbers(summary.get(9), 0, 1, 3, 3.0);
		check(summary.get(13).startsWith("0\t0\t0\tNA\t"), "An empty UNKNOWN common stratum has no invented mean");
		final ArrayList<String> zero=rows("result.zero_support.tsv");
		final ArrayList<String> gaps=rows("result.no_compatible.tsv");
		check(zero.size()==2 && zero.get(1).contains("\tp2\t") && zero.get(1).contains("/never-opened/α.fa")
			&& zero.get(1).contains("\t-1\t1,2,3\t4,2,1\t"), "Zero-support list joins full parent provenance and all selected ranks");
		check(gaps.size()==2 && gaps.get(1).contains("\tp3\t"), "Only the named no-compatible-model parent enters the rebuild list");
		NcrnaSavedRankSummary.main(new String[]{"in="+path("input"), "out="+path("cli"), "methods=0,1,2,3",
			"strata="+path("strata.tsv"), "expectedmethods=4", "expectedmodels=8", "expectedloci=5", "expectedcommon=2", "expectedzero=1", "expectedgaps=1", "t=1"});
		check(summary.equals(rows("cli.summary.tsv")), "CLI flags must reach the same saved-rank computation");
		NcrnaSavedRankSummary.process(path("input"), path("unassessed"), null, new int[]{0,1,2,3}, 4, 8, 5, 2, 1, 1);
		check(rows("unassessed.summary.tsv").get(9).startsWith("0\t0\t0\tNA\t"),
			"Omitting labels must not silently declare any locus unseen");
		checkNumbers(rows("unassessed.summary.tsv").get(13), 0, 2, 4, 2.0);
		fault("truncated", loci, loci.substring(0, loci.lastIndexOf("0\tinput\tp4\t3\t")));
		fault("duplicate_method", loci, loci.replace(row("p0", 1, 2, "RANKED"), row("p0", 0, 2, "RANKED")));
		fault("wrong_parent", loci, loci.replace("0\tinput\tp0\t1\t", "0\tinput\twrong\t1\t"));
		fault("bad_rank", loci, loci.replace(row("p0", 1, 2, "RANKED"), row("p0", 1, 9, "RANKED")));
		fault("missing_sketch", loci, loci.replace(row("p0", 1, 2, "RANKED"), row("p0", 1, -1, "RANKED")));
		fault("contradictory_sketch_reason", loci, loci.replace(row("p0", 1, 2, "RANKED"),
			row("p0", 1, 2, "RANKED").replace("\tRANKED\t", "\tNO_COMPATIBLE_MODEL\t")));
		fault("contradictory_gap_rank", loci, loci.replace(row("p3", 1, -1, "NO_COMPATIBLE_MODEL"), row("p3", 1, 2, "RANKED")));
		faultLabels("truncated_labels", strata, strata.substring(0, strata.lastIndexOf("0\tinput\tp4\t")));
		faultLabels("reordered_labels", strata, strata.replace("0\tinput\tp0\tSEEN\n0\tinput\tp1\tNOT_SEEN\n",
			"0\tinput\tp1\tNOT_SEEN\n0\tinput\tp0\tSEEN\n"));
		faultLabels("invalid_stratum", strata, strata.replace("\tp0\tSEEN\n", "\tp0\tMAYBE\n"));
		boolean rejected=false;
		try{process("wrong_count", 3);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A wrong common-set denominator must prevent a complete report");
		combined(parents, loci, strata);
	}

	/** Excluded loci remain validated; selected inputs without targets remain visible. */
	private void combined(final String parents, final String loci, final String strata){
		final String excluded="1\texcluded\t/never-opened/β.fa";
		write("subset.parents.tsv.gz", parents.replace(parent("p1", true),
			parent("p1", true).replace("0\tinput\t/never-opened/α.fa", excluded)));
		write("subset.loci.tsv.gz", loci.replace("0\tinput\tp1\t", "1\texcluded\tp1\t"));
		write("subset.labels.tsv", strata.replace("0\tinput\tp1\t", "1\texcluded\tp1\t"));
		final String assembly=NcrnaDdlLocusOverlap.ASSEMBLY_HEADER+'\n'
			+assembly(0, "input", "/never-opened/α.fa", 1)
			+assembly(1, "excluded", "/never-opened/β.fa", 2)
			+assembly(2, "empty", "/never-opened/empty.fa", 3);
		write("subset.assemblies.tsv", assembly);
		final String captured=NcrnaDdlTrainingOverlap.SAMPLE_HEADER+'\n'
			+captured(0, "input", "/never-opened/α.fa")
			+captured(1, "excluded", "/never-opened/β.fa")
			+captured(2, "empty", "/never-opened/empty.fa");
		write("subset.inputs.tsv", captured);
		final String reuse=NcrnaSavedRankSummary.REUSE_HEADER+"\n1\t0\tinput\t/never-opened/α.fa\n3\t2\tempty\t/never-opened/empty.fa\n";
		write("reuse.tsv", reuse);
		write("extra.parents.tsv.gz", NcrnaSavedRankSummary.PARENT_HEADER+'\n'
			+parent("p5", true).replace("0\tinput\t/never-opened/α.fa", "0\tnew\t/never-opened/γ.fa"));
		final int[] extraRanks={3, 2, 3, 1};
		final StringBuilder extra=new StringBuilder(NcrnaSavedRankSummary.LOCUS_HEADER).append('\n');
		for(int m=0; m<4; m++){
			extra.append(row("p5", m, extraRanks[m], "RANKED").replace("0\tinput\t", "0\tnew\t"));
		}
		write("extra.loci.tsv.gz", extra.toString());
		write("extra.labels.tsv", NcrnaSavedRankSummary.STRATA_HEADER+"\n0\tnew\tp5\tNOT_SEEN\n");
		write("extra.assemblies.tsv", NcrnaDdlLocusOverlap.ASSEMBLY_HEADER+'\n'+assembly(0, "new", "/never-opened/γ.fa", 4));
		write("extra.inputs.tsv", NcrnaDdlTrainingOverlap.SAMPLE_HEADER+'\n'+captured(0, "new", "/never-opened/γ.fa"));
		final String panels=NcrnaSavedRankSummary.PANEL_HEADER+"\nreused\t"+path("subset")+'\t'+path("subset.labels.tsv")
			+'\t'+path("subset.assemblies.tsv")+'\t'+path("reuse.tsv")+"\t3\t2\t5\t2\t1\t1\nnew\t"
			+path("extra")+'\t'+path("extra.labels.tsv")+'\t'+path("extra.assemblies.tsv")+"\t.\t1\t1\t1\t1\t0\t0\n";
		write("panels.tsv", panels);
		NcrnaSavedRankSummary.combine(path("panels.tsv"), path("combined"), new int[]{0,1,2,3}, 4, 8, 3);
		final ArrayList<String> summary=rows("combined.summary.tsv");
		check(summary.size()==13 && summary.get(1).contains("\tNOT_SEEN\t") && summary.get(5).contains("\tSEEN\t"),
			"Combined output reports not-seen first, seen second and unknown separately, with no pooled rows");
		checkNumbers(summary.get(1), 0, 1, 3, 3);
		checkNumbers(summary.get(5), 0, 1, 1, 1);
		check(summary.get(9).startsWith("0\t0\t0\tNA\t"), "An excluded not-seen locus cannot leak into another panel's sum");
		final ArrayList<String> inventory=rows("combined.inputs.tsv");
		check(inventory.size()==4 && inventory.get(1).endsWith("\t1\t4\t1\t1\t1\t1\t0")
			&& inventory.get(2).contains("\tempty\t/never-opened/empty.fa\t3\t0\t0\t0\t0\t0\t0"),
			"Selected group inventory preserves UTF-8 paths, all dispositions and the selected zero-target input");
		check(rows("combined.zero_support.tsv").size()==2 && rows("combined.zero_support.tsv").get(1).startsWith("reused\t"),
			"Named zero-support cases retain their source panel and bound parent metadata");
		NcrnaSavedRankSummary.main(new String[]{"panels="+path("panels.tsv"), "out="+path("combined_cli"),
			"methods=0,1,2,3", "expectedmethods=4", "expectedmodels=8", "expectedgroups=3", "t=1"});
		check(summary.equals(rows("combined_cli.summary.tsv")), "Combined CLI must use the same subset and histogram computation");
		write("wrong_reuse.tsv", reuse.replace("1\t0\tinput", "1\t1\tinput"));
		rejectCombined("wrong_reuse", panels.replace(path("reuse.tsv"), path("wrong_reuse.tsv")));
		write("wrong_path.tsv", reuse.replace("/never-opened/α.fa", "/never-opened/other.fa"));
		rejectCombined("wrong_path", panels.replace(path("reuse.tsv"), path("wrong_path.tsv")));
		write("duplicate_groups.tsv", NcrnaDdlLocusOverlap.ASSEMBLY_HEADER+'\n'+assembly(0, "new", "/never-opened/γ.fa", 1));
		rejectCombined("duplicate_group", panels.replace(path("extra.assemblies.tsv"), path("duplicate_groups.tsv")));
		write("excluded_bad_labels.tsv", strata.replace("0\tinput\tp1\tNOT_SEEN\n", ""));
		rejectCombined("excluded_bad_labels", panels.replace(path("subset.labels.tsv"), path("excluded_bad_labels.tsv")));
		write("subset.inputs.tsv", captured.replace(captured(2, "empty", "/never-opened/empty.fa"), ""));
		rejectCombined("missing_zero_target_input", panels);
		write("subset.inputs.tsv", captured);
		rejectCombined("wrong_source_count", panels.replace("\t3\t2\t5\t2\t1\t1", "\t3\t2\t5\t3\t1\t1"));
		fullCapture(panels, assembly);
	}

	/** Restore contributors; an alias-only panel contributes no duplicate locus. */
	private void fullCapture(final String panels, final String assembly){
		final String contributor=assembly.replace("\t2\t0\tNO_KNOWN_MEMBER_CONTRIBUTOR\t",
			"\t2\t1\tKNOWN_MEMBER_CONTRIBUTOR\t");
		check(!contributor.equals(assembly), "The full-capture fixture must actually restore a contributing assembly");
		write("contributors.tsv", contributor);
		final String full=panels.replace(path("reuse.tsv")+"\t3\t2", ".\t3\t3")
			.replace(path("subset.assemblies.tsv"), path("contributors.tsv"))
			.replace(path("extra.assemblies.tsv"), path("duplicate_groups.tsv"));
		write("full_sources.panels.tsv", full);
		NcrnaSavedRankSummary.main(new String[]{"panels="+path("full_sources.panels.tsv"), "out="+path("full"),
			"fullcapture=t", "methods=0,1,2,3", "expectedmethods=4", "expectedmodels=8", "t=1"});
		check(rows("full.population.tsv").get(1).equals("4\t3\t1"),
			"Every captured input is either one unique assembly group or a declared discarded alias");
		check(rows("full.panels.tsv").size()==2 && rows("full.panels.tsv").get(1).endsWith("\t3\t3\t5\t2\t1\t1"),
			"Effective panels retain actual selected counts and omit a panel containing aliases only");
		final ArrayList<String> inventory=rows("full.inputs.tsv");
		check(inventory.size()==4 && inventory.get(2).contains("\texcluded\t/never-opened/β.fa\t2\t1\t1\t0\t0\t0\t0"),
			"Full capture restores a contributor and preserves its original target disposition");
		checkNumbers(rows("full.summary.tsv").get(1), 0, 1, 3, 3);
		check(rows("full.aliases.tsv").size()==2 && rows("full.aliases.tsv").get(1).startsWith("new\t0\tnew\t"),
			"The duplicated assembly's second capture is named, never silently counted twice");
		rejectCombined("legacy_contributor", full);
		write("subset_request.tsv", panels);
		boolean rejected=false;
		try{NcrnaSavedRankSummary.combine(path("subset_request.tsv"), path("full_subset"), new int[]{0,1,2,3}, 4, 8, 0, true);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A subset manifest cannot masquerade as a full-capture census");
	}

	private void rejectCombined(final String name, final String panels){
		write(name+".panels.tsv", panels);
		boolean rejected=false;
		try{NcrnaSavedRankSummary.combine(path(name+".panels.tsv"), path(name), new int[]{0,1,2,3}, 4, 8, 3);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Combined source/subset disagreement must fail loudly: "+name);
	}

	private static String assembly(final int slot, final String id, final String path, final int group){
		return slot+"\t"+id+'\t'+path+'\t'+group+"\t0\tNO_KNOWN_MEMBER_CONTRIBUTOR\tFIXTURE\n";
	}
	private static String captured(final int slot, final String id, final String path){
		return slot+"\t"+id+'\t'+path+"\t1\t10\t0\t0\t0\t0\t0\t0\t0\t0\n";
	}

	private void fault(final String name, final String valid, final String invalid){
		check(!valid.equals(invalid), "A fixture fault must actually change saved rows: "+name);
		write("input.loci.tsv.gz", invalid);
		boolean rejected=false;
		try{process(name, 2);}catch(IllegalArgumentException expected){rejected=true;}
		finally{write("input.loci.tsv.gz", valid);}
		check(rejected, "Malformed saved report must fail loudly: "+name);
	}
	/** ByteBuilder intentionally prints integral doubles without a decimal suffix. */
	private static void checkNumbers(final String row, final int method, final int count, final long sum, final double mean){
		final LineParser1 p=new LineParser1('\t');
		p.set(row.getBytes(StandardCharsets.UTF_8));
		check(p.parseInt(0)==method && p.parseInt(1)==count && p.parseLong(2)==sum && p.parseDouble(3)==mean,
			"Literal method/common-count/sum/mean must agree numerically, including compact integral means: "+row);
	}
	private void process(final String output, final int common){
		NcrnaSavedRankSummary.process(path("input"), path(output), path("strata.tsv"), new int[]{0,1,2,3}, 4, 8, 5, common, 1, 1);
	}
	private void faultLabels(final String name, final String valid, final String invalid){
		check(!valid.equals(invalid), "A label fault must actually change the supplied metadata: "+name);
		write("strata.tsv", invalid);
		boolean rejected=false;
		try{process(name, 2);}catch(IllegalArgumentException expected){rejected=true;}
		finally{write("strata.tsv", valid);}
		check(rejected, "Incomplete or invalid seen-status labels must fail loudly: "+name);
	}
	private void write(final String name, final String text){
		final File previous=new File(path(name));
		check(!previous.exists() || previous.delete(), "Only this fixture's synthetic input may be replaced");
		try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private ArrayList<String> rows(final String name){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				result.add(new String(row, StandardCharsets.UTF_8));
			}
		}
		return result;
	}
	private static String parent(final String id, final boolean target){
		return "0\tinput\t/never-opened/α.fa\t"+id+"\tcontig\t1\t10\t+\tfalse\t"+target+"\tfalse\tfalse\t1\t0\n";
	}
	private static String row(final String id, final int method, final int rank, final String reason){
		final int windows=reason.equals("NO_CAPTURED_WINDOW") ? 0 : 1;
		final int compatible=windows==0 || reason.equals("NO_COMPATIBLE_MODEL") ? 0 : 1;
		return "0\tinput\t"+id+'\t'+method+'\t'+windows+'\t'+compatible+"\t0\t"+compatible+'\t'+rank+'\t'+reason+"\tfalse\tfalse\n";
	}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
	private final String root;
}
