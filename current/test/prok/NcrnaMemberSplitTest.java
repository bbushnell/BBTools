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
import prok.NcrnaMemberSplit;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Conserved historical rows with group-based aliases, heldout and outside roles.
 * @author Ganyu
 */
public final class NcrnaMemberSplitTest{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaMemberSplitTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Unknown ledger fixture option");
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append, "Use out=fresh_fixture_directory");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Ledger fixture directory must be fresh");
			fixture(root.getPath());
			pp.outstream.println("NCRNA_MEMBER_SPLIT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private static void fixture(final String root){
		final int[] group={0, 0, 1, 1, 2, 3};
		final String[] cohort={"evaluation", "discovery", "discovery", "evaluation", "evaluation", "evaluation"};
		final ArrayList<String> nodes=new ArrayList<String>();
		for(int i=0; i<6; i++){
			nodes.add(cohort[i]+'\t'+(i+1)+"\t/source"+i+".fa\t/source"+i+".fa\tassembly"+i+'\t'+group[i]
				+"\tfalse\tfalse\tfalse\tHISTORICAL\t"+group[i]+"\tfalse\tfalse\tHISTORICAL");
		}
		write(root+"/nodes.tsv", NcrnaConstructionMembers.GROUP_HEADER, nodes);
		write(root+"/split.groups.tsv", "assembly_group\tstratum\tpool_inputs\tassignment",
			Arrays.asList("0\tG0\t1\tTRAIN", "1\tG1\t1\tHELD_OUT", "2\tG2\t1\tCOVERAGE_CRITICAL_TRAIN"));
		write(root+"/split.inputs.tsv.gz", "assembly\tfasta\tassembly_group\tstratum\tassignment",
			Arrays.asList("A0\t/source0.fa\t0\tG0\tTRAIN", "A2\t/source2.fa\t1\tG1\tHELD_OUT", "A4\t/source4.fa\t2\tG2\tCOVERAGE_CRITICAL_TRAIN"));
		write(root+"/split.coassigned.tsv", "pool_group\tpartner_group\tinput_group\tfasta\tin_PGAP_pool\tassignment\treason", new ArrayList<String>());
		final ArrayList<String> olds=new ArrayList<String>();
		for(int i=0; i<5; i++){
			final String[] f=fields(39);
			f[5]="OLD_"+i;
			f[33]="historical_side";
			f[34]=Integer.toString(Math.min(i, 3));
			f[35]=i==0 ? "CORE_0;CORE_1" : i==1 ? "CORE_2;CORE_3" : i==2 ? "CORE_4" : i==3 ? "CORE_5" : "NA";
			f[36]=Boolean.toString(i==4);
			f[37]=Boolean.toString(i!=4);
			f[38]=Boolean.toString(i!=4);
			olds.add(String.join("\t", f));
		}
		write(root+"/members.old.tsv.gz", NcrnaConstructionMembers.OLD_OUTPUT_HEADER, olds);
		final ArrayList<String> events=new ArrayList<String>();
		final ArrayList<String> cores=new ArrayList<String>();
		for(int i=0; i<7; i++){
			final int node=Math.min(i, 5);
			final String[] e=fields(36);
			e[0]=cohort[node];
			e[2]=Integer.toString(i+1);
			e[3]=Integer.toString(node+1);
			e[8]="RF00050";
			e[33]="historical_side";
			e[34]=Integer.toString(group[node]);
			e[35]="CORE_"+i;
			events.add(String.join("\t", e));
			final String[] f=fields(21);
			f[0]="CORE_"+i;
			f[1]=cohort[node];
			f[2]=Integer.toString(node+1);
			f[3]="/source"+node+".fa";
			f[4]=Integer.toString(group[node]);
			f[11]="138";
			f[13]="historical_side";
			f[14]="true";
			f[15]="true";
			f[16]="false";
			f[17]="false";
			f[18]=i==6 ? "NA" : "OLD_"+group[node];
			f[20]="PENDING_SEQUENCE_AND_FLANK_VALIDATION";
			cores.add(String.join("\t", f));
		}
		write(root+"/members.events.tsv.gz", NcrnaConstructionMembers.FMN_HEADER+"\tadmission\tassembly_split\tassembly_group\textraction_member", events);
		write(root+"/members.extract.tsv.gz", NcrnaConstructionMembers.EXTRACT_HEADER, cores);
		final NcrnaMemberSplit.Config c=new NcrnaMemberSplit.Config();
		c.members=root+"/members";
		c.groups=root+"/nodes.tsv";
		c.split=root+"/split";
		c.output=root+"/result";
		c.inputs=6;
		c.pool=3;
		c.old=5;
		c.events=7;
		c.cores=7;
		NcrnaMemberSplit.run(c);
		preserved(olds, rows(c.output+".old.tsv.gz"));
		preserved(events, rows(c.output+".events.tsv.gz"));
		final ArrayList<String> out=rows(c.output+".extract.tsv.gz");
		preserved(cores, out);
		check(out.get(2).endsWith("\tTRAIN\ttrue\ttrue\tNA\tGROUP_SEQUENCE_DEDUP_PENDING"), "An outside path in a pool training group follows that group");
		check(out.get(4).endsWith("\tHELD_OUT\tfalse\tfalse\tNA\tGROUP_SEQUENCE_DEDUP_PENDING"), "An outside path in a heldout group cannot leak into training");
		check(out.get(6).endsWith("\tOUTSIDE_POOL\ttrue\ttrue\tNA\tGROUP_SEQUENCE_DEDUP_PENDING"), "Uncontradicted outside-pool old members remain eligible");
		check(out.get(7).endsWith("\tOUTSIDE_POOL\tfalse\tfalse\tNA\tGROUP_SEQUENCE_DEDUP_PENDING"), "Outside-pool native-only cores remain recorded without expanding the training pool");
		check(rows(c.output+".old.tsv.gz").get(5).endsWith("\tOUTSIDE_POOL\tfalse\tfalse\tNA"), "Positive contradiction remains excluded under the outside-pool rule");
		check(rows(c.output+".labels.tsv").contains("OUTSIDE_POOL\t2\t1\t1\t1\t1\t2\t2\t2\t2\t1\t1"), "Outside-pool candidate counts and methodology eligibility must stay distinct");
		check(rows(c.output+".summary.tsv").get(1).equals("5\t7\t7\t3\t6\tMETADATA_REBOUND_GROUP_SEQUENCE_DEDUP_REQUIRED_BEFORE_LEARNING"), "All existing row populations must be conserved independently");
		NcrnaMemberSplit.main(new String[]{"in="+c.members, "groups="+c.groups, "split="+c.split, "out="+root+"/cli",
			"expectedinputs=6", "expectedpool=3", "expectedold=5", "expectedevents=7", "expectedcores=7", "t=1"});
		check(rows(root+"/cli.labels.tsv").equals(rows(c.output+".labels.tsv")), "CLI and direct label counts must match");
		final ArrayList<String> pairedNodes=new ArrayList<String>();
		for(final String row:nodes){
			final String[] f=row.split("\t", -1);
			if(f[10].equals("0") || f[10].equals("2")){f[12]="true";}
			pairedNodes.add(String.join("\t", f));
		}
		write(root+"/paired_nodes.tsv", NcrnaConstructionMembers.GROUP_HEADER, pairedNodes);
		write(root+"/paired.groups.tsv", "assembly_group\tstratum\tpool_inputs\tassignment",
			Arrays.asList("0\tG0\t1\tCOVERAGE_CRITICAL_TRAIN", "1\tG1\t1\tHELD_OUT", "2\tG0\t1\tCOVERAGE_CRITICAL_TRAIN"));
		write(root+"/paired.inputs.tsv.gz", "assembly\tfasta\tassembly_group\tstratum\tassignment",
			Arrays.asList("A0\t/source0.fa\t0\tG0\tCOVERAGE_CRITICAL_TRAIN", "A2\t/source2.fa\t1\tG1\tHELD_OUT", "A4\t/source4.fa\t2\tG0\tCOVERAGE_CRITICAL_TRAIN"));
		final ArrayList<String> pairRows=new ArrayList<String>(Arrays.asList(
			"0\t2\t0\t/source0.fa\ttrue\tCOVERAGE_CRITICAL_TRAIN\tunresolved, co-assigned to prevent leakage",
			"0\t2\t0\t/source1.fa\tfalse\tCOVERAGE_CRITICAL_TRAIN\tunresolved, co-assigned to prevent leakage",
			"0\t2\t2\t/source4.fa\ttrue\tCOVERAGE_CRITICAL_TRAIN\tunresolved, co-assigned to prevent leakage"));
		final String pairHeader="pool_group\tpartner_group\tinput_group\tfasta\tin_PGAP_pool\tassignment\treason";
		write(root+"/paired.coassigned.tsv", pairHeader, pairRows);
		c.groups=root+"/paired_nodes.tsv";
		c.split=root+"/paired";
		c.output=root+"/paired_result";
		NcrnaMemberSplit.run(c);
		check(rows(c.output+".old.tsv.gz").get(1).endsWith("\tCOVERAGE_CRITICAL_TRAIN\ttrue\ttrue\tunresolved, co-assigned to prevent leakage"),
			"Unresolved same-side provenance must accompany eligible old-member methodology labels");
		pairRows.remove(1);
		write(root+"/missing_pair.groups.tsv", "assembly_group\tstratum\tpool_inputs\tassignment", rows(root+"/paired.groups.tsv").subList(1, 4));
		write(root+"/missing_pair.inputs.tsv.gz", "assembly\tfasta\tassembly_group\tstratum\tassignment", rows(root+"/paired.inputs.tsv.gz").subList(1, 4));
		write(root+"/missing_pair.coassigned.tsv", pairHeader, pairRows);
		c.split=root+"/missing_pair";
		c.output=root+"/missing_alias";
		boolean aliasRejected=false;
		try{NcrnaMemberSplit.run(c);}
		catch(IllegalArgumentException failure){aliasRejected=failure.getMessage().contains("every original alias path");}
		check(aliasRejected && !new File(c.output+".old.tsv.gz").exists(), "A truncated unresolved-pair table must fail before publishing member assignments");
		c.groups=root+"/nodes.tsv";
		c.split=root+"/split";
		final String[] wrong=cores.get(3).split("\t", -1);
		wrong[4]="0";
		cores.set(3, String.join("\t", wrong));
		write(root+"/bad.old.tsv.gz", NcrnaConstructionMembers.OLD_OUTPUT_HEADER, olds);
		write(root+"/bad.events.tsv.gz", NcrnaConstructionMembers.FMN_HEADER+"\tadmission\tassembly_split\tassembly_group\textraction_member", events);
		write(root+"/bad.extract.tsv.gz", NcrnaConstructionMembers.EXTRACT_HEADER, cores);
		c.members=root+"/bad";
		c.output=root+"/rejected";
		boolean rejected=false;
		try{NcrnaMemberSplit.run(c);}
		catch(IllegalArgumentException failure){rejected=failure.getMessage().contains("original source path/group");}
		check(rejected, "A heldout alias cannot gain training eligibility by changing only its saved core group");
	}
	private static String[] fields(final int n){
		final String[] f=new String[n];
		Arrays.fill(f, "preserved_α");
		return f;
	}
	private static void preserved(final ArrayList<String> source, final ArrayList<String> result){
		check(result.size()==source.size()+1, "Appending methodology labels must conserve every original row");
		for(int i=0; i<source.size(); i++){check(result.get(i+1).startsWith(source.get(i)+'\t'), "All original evidence and historical eligibility fields must remain byte-equivalent, including UTF8");}
	}
	private static void write(final String path, final String header, final java.util.List<String> rows){
		try(Output out=new Output(path)){
			out.writer.println(header);
			for(final String row:rows){out.writer.print((row+'\n').getBytes(StandardCharsets.UTF_8));}
		}
	}
	private static ArrayList<String> rows(final String path){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path)){for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){rows.add(new String(row, StandardCharsets.UTF_8));}}
		return rows;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
}
