package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaConstructionMembers;
import prok.NcrnaDdlTrainingOverlap;
import prok.NcrnaMemberLengthBounds;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal assembly aliases, actual kept votes and explicit unresolved bindings.
 * Metadata paths intentionally do not exist; no sequence input is opened.
 * @author Raiden
 */
public final class NcrnaDdlTrainingOverlapTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlTrainingOverlapTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown assembly-overlap fixture argument: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_directory without overwrite, append or sampling");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Assembly-overlap fixture directory must be fresh");
			new NcrnaDdlTrainingOverlapTest(root.getPath()).run();
			pp.outstream.println("NCRNA_DDL_TRAINING_OVERLAP_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaDdlTrainingOverlapTest(final String root_){root=root_;}

	private void run(){
		final String[] paths={"/never-opened/α.fa", "/never-opened/alias.fa", "/never-opened/orphan.fa",
			"/never-opened/small.fa", "/never-opened/no_member.fa", "/never-opened/missing.fa"};
		String sample=NcrnaDdlTrainingOverlap.SAMPLE_HEADER+'\n';
		for(int i=0; i<paths.length; i++){
			sample+=i+"\t"+(10+i)+"\t"+paths[i]+"\t1\t100\t0\t0\t0\t0\t0\t0\t0\t0\n";
		}
		String groups=NcrnaConstructionMembers.GROUP_HEADER+'\n';
		groups+=group(10, paths[0], "/never-opened/sourceA.fa", 0);
		groups+=group(11, paths[1], "/never-opened/sourceA.fa", 0);
		for(int i=2; i<5; i++){groups+=group(10+i, paths[i], paths[i], i-1);}
		String genomes=NcrnaConstructionMembers.GENOME_HEADER+'\n';
		final String[] donorPaths={"/never-opened/sourceA.fa", paths[2], paths[3], paths[4], "/never-opened/unmapped_donor.fa"};
		for(int i=0; i<donorPaths.length; i++){
			final String genome=i==0 ? "gα" : "g"+i;
			genomes+=i+"\t"+genome+"\tassembly"+i+"\t7\tbacteria\tP\t/never-opened/"+i+".gff\t"+donorPaths[i]+"\t0\t100\t1\n";
		}
		final String census=NcrnaMemberLengthBounds.HEADER+'\n'+model("m0", 2)+model("mβ", 1);
		final String assignments="#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n"
			+assignment("a", "gα", "m0", "kept")+assignment("b", "gα", "m0", "kept")
			+assignment("c", "gα", "mβ", "kept")+assignment("d", "g1", "", "orphan")
			+assignment("e", "g2", "", "small");
		write("sample.tsv", sample);
		write("groups.tsv", groups);
		write("genomes.tsv", genomes);
		write("census.tsv", census);
		write("assignments.tsv", assignments);
		process("result");
		final ArrayList<String> summary=rows("result.summary.tsv");
		check(summary.size()==2 && summary.get(1).equals("6\t2\t3\t1\t1\t3\t0\t1\t3\t1\t0\t0\tKEPT_ASSIGNMENT_ASSEMBLY_CENSUS_NOT_EXACT_SEQUENCE_OR_HOLDOUT_PROOF"),
			"Two sample aliases share one contributing assembly; three actual kept votes come from one donor source");
		final ArrayList<String> inputs=rows("result.inputs.tsv");
		check(inputs.size()==7 && inputs.get(1).startsWith("0\t10\t"+paths[0]+"\t0\t3\tKNOWN_MEMBER_CONTRIBUTOR\t"),
			"Canonical path aliases retain the original UTF-8 sample path and actual vote count");
		check(inputs.get(2).contains("\t0\t3\tKNOWN_MEMBER_CONTRIBUTOR\t"), "Known aliases must have the same contribution status");
		for(int i=3; i<=5; i++){
			check(inputs.get(i).contains("\t0\tNO_KNOWN_MEMBER_CONTRIBUTOR\t"), "Orphans, small clusters and a roster entry alone do not contribute votes");
		}
		check(inputs.get(6).contains("\t-1\t0\tUNRESOLVED_ASSEMBLY_CONTRIBUTION\t"), "A missing binding is unresolved rather than an unseen claim");
		NcrnaDdlTrainingOverlap.main(new String[]{"in="+path("assignments.tsv"), "genomes="+path("genomes.tsv"),
			"groups="+path("groups.tsv"), "sample="+path("sample.tsv"), "census="+path("census.tsv"), "out="+path("cli"),
			"expectedinputs=6", "expectedgenomes=5", "expectedmodels=2", "expectedmembers=5", "expectedorphans=1", "expectedsmall=1", "t=1"});
		check(summary.equals(rows("cli.summary.tsv")) && inputs.equals(rows("cli.inputs.tsv")), "CLI and direct metadata audit preserve the same literal result");
		write("assignments_unknown.tsv", assignments.replace("c source=gα", "c source=g4"));
		NcrnaDdlTrainingOverlap.process(path("assignments_unknown.tsv"), path("genomes.tsv"), path("groups.tsv"),
			path("sample.tsv"), path("census.tsv"), path("unknown"), 6, 5, 2, 5, 1, 1);
		check(rows("unknown.summary.tsv").get(1).equals("6\t2\t0\t4\t1\t0\t3\t1\t3\t2\t1\t0\tKEPT_ASSIGNMENT_ASSEMBLY_CENSUS_NOT_EXACT_SEQUENCE_OR_HOLDOUT_PROOF"),
			"An unmapped actual donor prevents negative assembly conclusions while preserving positive evidence");
		write("groups_uncertain.tsv", groups.replace(group(12, paths[2], paths[2], 1),
			group(12, paths[2], paths[2], 1).replace("\tfalse\tACCEPTED\n", "\ttrue\tACCEPTED\n")));
		NcrnaDdlTrainingOverlap.process(path("assignments.tsv"), path("genomes.tsv"), path("groups_uncertain.tsv"),
			path("sample.tsv"), path("census.tsv"), path("uncertain"), 6, 5, 2, 5, 1, 1);
		check(rows("uncertain.summary.tsv").get(1).startsWith("6\t2\t2\t2\t1\t2\t1\t1\t"),
			"An unresolved assembly-pair nomination remains visible for that noncontributing group");
		fault("duplicate_assignment", "assignments.tsv", assignments, assignments+assignment("a", "gα", "m0", "kept"));
		fault("wrong_model", "assignments.tsv", assignments, assignments.replace("\tm0\tkept\n", "\tforeign\tkept\n"));
		fault("nonvoting_model", "assignments.tsv", assignments, assignments.replace("\t\torphan\n", "\tm0\torphan\n"));
		fault("wrong_census", "census.tsv", census, census.replace(model("mβ", 1), model("mβ", 2)));
		fault("invalid_identity", "sample.tsv", sample, sample.replace(paths[0], "/never-opened/\ufffd.fa"));
		fault("path_conflict", "groups.tsv", groups, groups.replace(group(12, paths[2], paths[2], 1), group(12, paths[2], donorPaths[0], 1)));
		fault("evaluation_misbind", "groups.tsv", groups, groups.replace(paths[1], "/never-opened/changed.fa"));
		fault("missing_sample", "sample.tsv", sample, sample.substring(0, sample.lastIndexOf("5\t15\t")));
		selection(paths, groups);
	}

	/** Captured aliases win over earlier uncaptured paths in the same assembly. */
	private void selection(final String[] paths, final String groups){
		final String alias="/never-opened/alias_orphan.fa";
		String selected=NcrnaDdlTrainingOverlap.SELECTION_HEADER+'\n';
		for(int i=0; i<paths.length; i++){selected+=selectionRow(10+i, paths[i]);}
		selected+=selectionRow(16, alias);
		write("selection.tsv", selected);
		write("selection_groups.tsv", groups+group(16, alias, paths[2], 1));
		final String captured=NcrnaDdlTrainingOverlap.SAMPLE_HEADER+'\n'
			+captureRow(0, 11, paths[1])+captureRow(1, 16, alias)+captureRow(2, 14, paths[4]);
		write("captured.tsv", captured);
		select("expansion", "captured.tsv");
		final ArrayList<String> summary=rows("expansion.expansion.tsv");
		check(summary.size()==2 && summary.get(1).equals("3\t2\t1\t100\t1\tASSEMBLY_SELECTION_NOT_EXACT_SEQUENCE_HOLDOUT_PROOF"),
			"Exactly three negative assembly groups split into two reused captures and one new genome; unknown and contributor are excluded");
		final ArrayList<String> reuse=rows("expansion.reuse.tsv");
		check(reuse.size()==3 && reuse.contains("1\t1\t16\t"+alias) && reuse.contains("3\t2\t14\t"+paths[4]),
			"Later captured alias replaces an earlier uncaptured path and retains its original capture slot");
		check(rows("expansion.uncaptured.selection.tsv").equals(java.util.Arrays.asList(
			NcrnaDdlTrainingOverlap.SELECTION_HEADER, selectionRow(13, paths[3]).trim())),
			"The uncaptured representative retains its literal frozen selection row without invented caller accounting");
		check(rows("expansion.uncaptured.config").equals(java.util.Arrays.asList("in="+paths[3])),
			"Only one uncaptured representative is sent to the future caller, never the reused or unknown inputs");
		NcrnaDdlTrainingOverlap.main(new String[]{"in="+path("assignments.tsv"), "genomes="+path("genomes.tsv"),
			"groups="+path("selection_groups.tsv"), "sample="+path("selection.tsv"), "census="+path("census.tsv"),
			"captured="+path("captured.tsv"), "out="+path("expansion_cli"), "expectedcaptured=3", "expectedinputs=7",
			"expectedgenomes=5", "expectedmodels=2", "expectedmembers=5", "expectedorphans=1", "expectedsmall=1", "t=1"});
		for(final String suffix:new String[]{".inputs.tsv", ".summary.tsv", ".reuse.tsv", ".uncaptured.selection.tsv", ".uncaptured.config", ".expansion.tsv"}){
			check(rows("expansion"+suffix).equals(rows("expansion_cli"+suffix)), "Expansion CLI must preserve each selected/reused population: "+suffix);
		}
		selectionFault("changed_capture_path", captured.replace(alias, "/never-opened/foreign.fa"));
		selectionFault("changed_capture_size", captured.replace("\t1\t100\t", "\t1\t101\t"));
		selectionFault("duplicate_capture", captured.replace(captureRow(2, 14, paths[4]), captureRow(2, 16, alias)));
		selectionFault("missing_capture", captured.replace(captureRow(2, 14, paths[4]), ""));
	}

	private static String selectionRow(final int id, final String fasta){
		return id+"\t"+fasta+"\t100\t1\tbacteria\tP\t6\t0\n";
	}

	private static String captureRow(final int slot, final int id, final String fasta){
		return slot+"\t"+id+'\t'+fasta+"\t1\t100\t0\t0\t0\t0\t0\t0\t0\t0\n";
	}

	private void select(final String output, final String captured){
		NcrnaDdlTrainingOverlap.process(path("assignments.tsv"), path("genomes.tsv"), path("selection_groups.tsv"),
			path("selection.tsv"), path("census.tsv"), path(output), 7, 5, 2, 5, 1, 1, path(captured), 3);
	}

	private void selectionFault(final String name, final String captured){
		write(name+".captured.tsv", captured);
		boolean rejected=false;
		try{select(name, name+".captured.tsv");}
		catch(final IllegalArgumentException expected){rejected=true;}
		check(rejected, "A changed captured population must not cause recapture or suppress an uncaptured group: "+name);
	}

	private static String group(final int id, final String fasta, final String canonical, final int group){
		return "evaluation\t"+id+'\t'+fasta+'\t'+canonical+"\trecord\t"+group+"\ttrue\tfalse\tfalse\tACCEPTED\t"+group+"\ttrue\tfalse\tACCEPTED\n";
	}

	private static String model(final String id, final int members){
		return id+"\tRF00174\t100\t"+members+"\t100\t100\t100\t100\t100\t1\tPURE\t"+members+"\t1\t0\t0\t0\t0\t0\t0\n";
	}

	private static String assignment(final String id, final String source, final String model, final String status){
		return id+" source="+source+"\t7\tgroup\tcluster\t"+model+'\t'+status+'\n';
	}

	private void process(final String output){
		NcrnaDdlTrainingOverlap.process(path("assignments.tsv"), path("genomes.tsv"), path("groups.tsv"),
			path("sample.tsv"), path("census.tsv"), path(output), 6, 5, 2, 5, 1, 1);
	}

	private void fault(final String name, final String input, final String original, final String changed){
		check(!original.equals(changed), "A malformed fixture must actually change its input: "+name);
		write("fault_"+name+".tsv", changed);
		boolean rejected=false;
		try{
			NcrnaDdlTrainingOverlap.process(path(input.equals("assignments.tsv") ? "fault_"+name+".tsv" : "assignments.tsv"),
				path(input.equals("genomes.tsv") ? "fault_"+name+".tsv" : "genomes.tsv"),
				path(input.equals("groups.tsv") ? "fault_"+name+".tsv" : "groups.tsv"),
				path(input.equals("sample.tsv") ? "fault_"+name+".tsv" : "sample.tsv"),
				path(input.equals("census.tsv") ? "fault_"+name+".tsv" : "census.tsv"), path(name), 6, 5, 2, 5, 1, 1);
		}catch(final IllegalArgumentException expected){rejected=true;}
		check(rejected, "Malformed training-overlap input was accepted: "+name);
	}

	private void write(final String name, final String text){
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

	private String path(final String name){return root+'/'+name;}
	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
	private final String root;
}
