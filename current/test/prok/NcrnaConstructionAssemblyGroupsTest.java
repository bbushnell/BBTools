package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaConstructionAssemblyCompare;
import prok.NcrnaConstructionAssemblyGroups;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal metadata-only alias propagation, conservative nonmatches and conservation.
 * FASTA paths deliberately need not exist: this stage must never reopen payloads.
 * @author Ganyu
 */
public final class NcrnaConstructionAssemblyGroupsTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionAssemblyGroupsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown assembly-group fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Supply out=fresh_directory without inputs, overwrite or sampling");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Group fixtures need a fresh directory");
			new NcrnaConstructionAssemblyGroupsTest(root.getPath()).run();
			pp.outstream.println("NCRNA_CONSTRUCTION_ASSEMBLY_GROUPS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaConstructionAssemblyGroupsTest(final String root_){root=root_;}

	private void run(){
		final String inputs=NcrnaConstructionAssemblyGroups.INPUT_HEADER+'\n'
			+node(true, "1", 0, true, false, "EVALUATION")
			+node(true, "2", 1, true, false, "EVALUATION")
			+node(true, "3", 2, true, false, "EVALUATION")
			+node(false, "1", 3, false, true, "ASSEMBLY_SEQUENCE_CHECK_PENDING")
			+node(false, "2", 3, false, true, "ASSEMBLY_SEQUENCE_CHECK_PENDING")
			+node(false, "3", 5, false, true, "ASSEMBLY_SEQUENCE_CHECK_PENDING")
			+node(false, "4", 6, false, false, "TRAIN_CANDIDATE_ASSEMBLY_IDENTITY")
			+node(false, "5", 0, true, false, "EXCLUDE_KNOWN_EVALUATION_ASSEMBLY");
		final String one=pair("1", "1");
		final String two=pair("2", "1");
		final String three=pair("3", "1");
		final String four=pair("3", "3");
		final String nominations=NcrnaConstructionAssemblyCompare.HEADER+'\n'+one+two+three+four;
		final String comparisons=NcrnaConstructionAssemblyCompare.HEADER+"\tsequence_identity\tcriterion\n"
			+result(one, true)+result(two, true)+result(three, false)+result(four, false);
		write("inputs.tsv", inputs);
		write("nominations.tsv", nominations);
		write("comparisons.tsv", comparisons);
		process("inputs.tsv", "nominations.tsv", "comparisons.tsv", "result", 4);
		checkResult("result");
		NcrnaConstructionAssemblyGroups.main(new String[]{"in="+path("inputs.tsv"),
			"nominations="+path("nominations.tsv"), "comparisons="+path("comparisons.tsv"),
			"out="+path("cli"), "expectedeval=3", "expectednew=5", "expectedpairs=4", "t=1"});
		checkResult("cli");

		fault("missing_result", inputs, nominations, comparisons.replace(result(four, false), ""), 4);
		fault("extra_result", inputs, nominations, comparisons+result(one, true), 4);
		fault("duplicate_pair", inputs, nominations+one, comparisons+result(one, true), 5);
		fault("wrong_count", inputs, nominations, comparisons, 5);
		fault("wrong_basis", inputs, nominations.replace(NcrnaConstructionAssemblyCompare.BASIS, "GUESS"), comparisons, 4);
		fault("wrong_criterion", inputs, nominations, comparisons.replace(NcrnaConstructionAssemblyCompare.CRITERION, "SHORT_GENE_MATCH"), 4);
		fault("unknown_outcome", inputs, nominations, comparisons.replace(NcrnaConstructionAssemblyCompare.MATCH, "INDEPENDENT"), 4);
		fault("wrong_path", inputs.replace(fasta(false, "1"), "/changed/input.fa"), nominations, comparisons, 4);
		fault("orphan_pending", inputs.replace(node(false, "4", 6, false, false, "TRAIN_CANDIDATE_ASSEMBLY_IDENTITY"),
			node(false, "4", 6, false, true, "ASSEMBLY_SEQUENCE_CHECK_PENDING")), nominations, comparisons, 4);
		fault("group_flags", inputs.replace(node(false, "2", 3, false, true, "ASSEMBLY_SEQUENCE_CHECK_PENDING"),
			node(false, "2", 3, false, false, "TRAIN_CANDIDATE_ASSEMBLY_IDENTITY")), nominations, comparisons, 4);
		fault("duplicate_input", inputs.replace(node(false, "2", 3, false, true, "ASSEMBLY_SEQUENCE_CHECK_PENDING"),
			node(false, "1", 3, false, true, "ASSEMBLY_SEQUENCE_CHECK_PENDING")), nominations, comparisons, 4);
		fault("mislabelled_disposition", inputs.replace("TRAIN_CANDIDATE_ASSEMBLY_IDENTITY", "EXCLUDE_KNOWN_EVALUATION_ASSEMBLY"), nominations, comparisons, 4);
	}

	private void checkResult(final String name){
		final ArrayList<String> inputs=rows(path(name+".inputs.tsv.gz"));
		final ArrayList<String> pairs=rows(path(name+".pairs.tsv.gz"));
		check(inputs.size()==9 && pairs.size()==5, "All original inputs and pair outcomes must survive integration");
		check(inputs.get(4).endsWith("\t0\ttrue\ttrue\tEXCLUDE_VERIFIED_EVALUATION_ASSEMBLY_ALIAS")
			&& inputs.get(5).endsWith("\t0\ttrue\ttrue\tEXCLUDE_VERIFIED_EVALUATION_ASSEMBLY_ALIAS"),
			"A positive match excludes its original group alias too; a separate nonmatch stays visible");
		check(inputs.get(6).endsWith("\t5\tfalse\ttrue\tASSEMBLY_SEQUENCE_IDENTITY_UNRESOLVED"),
			"A discovery input with only a nonmatch must never become an independent training candidate");
		check(inputs.get(2).endsWith("\t0\ttrue\ttrue\tEVALUATION")
			&& inputs.get(3).endsWith("\t2\ttrue\ttrue\tEVALUATION"),
			"Two positive evaluation links join transitively; the nonmatching evaluation group stays separate");
		check(inputs.get(7).endsWith("\t6\tfalse\tfalse\tTRAIN_CANDIDATE_ASSEMBLY_IDENTITY")
			&& inputs.get(8).endsWith("\t0\ttrue\ttrue\tEXCLUDE_KNOWN_EVALUATION_ASSEMBLY"),
			"Unchanged candidates and preexisting exclusions retain distinct dispositions");
		check(pairs.get(3).contains(NcrnaConstructionAssemblyCompare.UNRESOLVED)
			&& pairs.get(3).endsWith("\t2\t0") && pairs.get(4).endsWith("\t2\t5"),
			"Neither unresolved pair becomes a positive edge, even beside independent positive evidence");
		check(inputs.get(4).contains(fasta(false, "1")), "UTF-8 input provenance must survive verbatim");
		check(rows(path(name+".summary.tsv")).get(1).equals("3\t1\t2\t1\t0\t1\t2\t2\t4\tASSEMBLY_GROUPS_NOT_FROZEN_MEMBER_SPLIT"),
			"Unique input counts differ from pair counts and conserve three evaluation plus five discovery inputs");
	}

	private void fault(final String name, final String inputs, final String nominations, final String comparisons, final int count){
		write(name+".inputs.tsv", inputs);
		write(name+".nominations.tsv", nominations);
		write(name+".comparisons.tsv", comparisons);
		boolean rejected=false;
		try{process(name+".inputs.tsv", name+".nominations.tsv", name+".comparisons.tsv", name, count);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Malformed saved evidence must fail before group output: "+name);
		check(!new File(path(name+".summary.tsv")).exists(), "Rejected evidence cannot publish a complete summary");
	}
	private void process(final String inputs, final String nominations, final String comparisons, final String out, final int count){
		NcrnaConstructionAssemblyGroups.process(path(inputs), path(nominations), path(comparisons), path(out), 3, 5, count);
	}
	private static String node(final boolean eval, final String id, final int group, final boolean prior, final boolean pending, final String disposition){
		final String path=fasta(eval, id);
		return (eval ? "evaluation" : "discovery")+'\t'+id+'\t'+path+'\t'+path+"\tGCF_"+group+'\t'
			+group+'\t'+prior+'\t'+pending+"\tfalse\t"+disposition+'\n';
	}
	private static String pair(final String eval, final String discovery){
		return eval+'\t'+discovery+'\t'+fasta(true, eval)+'\t'+fasta(false, discovery)
			+"\t20\t4\t"+NcrnaConstructionAssemblyCompare.BASIS+'\n';
	}
	private static String result(final String row, final boolean match){
		return row.substring(0, row.length()-1)+'\t'+(match ? NcrnaConstructionAssemblyCompare.MATCH
			: NcrnaConstructionAssemblyCompare.UNRESOLVED)+'\t'+NcrnaConstructionAssemblyCompare.CRITERION+'\n';
	}
	private static String fasta(final boolean evaluation, final String id){return "/metadata_only_\u00e9/"+(evaluation ? "E" : "D")+id+".fa";}
	private String path(final String name){return root+'/'+name;}
	private void write(final String name, final String text){try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}}
	private static ArrayList<String> rows(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				result.add(new String(row, StandardCharsets.UTF_8));
			}
		}
		return result;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final String root;
}
