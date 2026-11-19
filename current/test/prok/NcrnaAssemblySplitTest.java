package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;

import map.ObjectMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaAssemblySplit;
import prok.NcrnaConstructionMembers;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal rare/ordinary/family-fallback populations and order/annotation invariance.
 * @author Ganyu
 */
public final class NcrnaAssemblySplitTest{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaAssemblySplitTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Unknown split test option");
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append, "Use out=fresh_fixture_directory");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Split fixtures require a fresh output directory");
			fixture(root.getPath());
			pp.outstream.println("NCRNA_ASSEMBLY_SPLIT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private static void fixture(final String root){
		final ArrayList<String> roster=new ArrayList<String>();
		final ArrayList<String> groups=new ArrayList<String>();
		final ArrayList<String> taxonomy=new ArrayList<String>();
		for(int i=0; i<20; i++){
			final int group=i==19 ? 0 : i;
			final String accession="GCF_"+i+".1";
			final String fasta="/original/"+accession+".fa.gz";
			final String family=group<4 ? "A" : group<9 ? "B" : "C";
			final String genus=group<4 ? "A" : group<9 ? "B" : "";
			roster.add(i+"\tgenome"+i+'\t'+accession+"\t1\tbacteria\tP\t/original/gff"+i+'\t'+fasta+"\t0\t0\t0");
			groups.add("discovery\t"+(i+1)+'\t'+fasta+'\t'+fasta+'\t'+accession+'\t'+group
				+"\tfalse\tfalse\tfalse\tTRAIN_CANDIDATE_ASSEMBLY_IDENTITY\t"+group+"\tfalse\tfalse\tTRAIN_CANDIDATE_ASSEMBLY_IDENTITY");
			taxonomy.add("RS_"+accession+"\td__Bacteria;p__P;c__C;o__O;f__"+family+";g__"+genus+";s__species"+i);
		}
		write(root+"/roster.tsv", NcrnaConstructionMembers.GENOME_HEADER, roster);
		write(root+"/groups.tsv", NcrnaConstructionMembers.GROUP_HEADER, groups);
		write(root+"/taxonomy.tsv", null, taxonomy);
		NcrnaAssemblySplit.run(root+"/roster.tsv", root+"/groups.tsv", root+"/taxonomy.tsv", root+"/first", 20261113, 20);
		final ArrayList<String> output=rows(root+"/first.groups.tsv");
		check(output.size()==20, "Twenty inputs with one alias must make nineteen assembly groups plus header");
		final ObjectMap<String, String> assignment=new ObjectMap<String, String>(String.class, String.class);
		int rare=0, heldB=0, heldC=0;
		for(int i=1; i<output.size(); i++){
			final String[] f=output.get(i).split("\t", -1);
			check(f.length==4 && assignment.put(f[0], f[3])==null, "No group may have two assignments");
			final int group=Integer.parseInt(f[0]);
			if(group<4){check(f[3].equals("COVERAGE_CRITICAL_TRAIN"), "All four rare groups train, even with five alias paths"); rare++;}
			else if(f[3].equals("HELD_OUT")){
				if(group<9){heldB++;}else{heldC++;}
			}
		}
		check(rare==4 && heldB==1 && heldC==2, "Literal four/five/ten-group strata must conserve zero/one/two heldout groups");
		for(final String row:rows(root+"/first.inputs.tsv.gz").subList(1, 21)){
			final String[] f=row.split("\t", -1);
			check(f[4].equals(assignment.get(f[2])), "Every alias path must inherit its group's one side");
		}
		check(rows(root+"/first.summary.tsv").get(1).startsWith("20\t19\t3\t3\t12\t4\t20261113\t"),
			"Summary must conserve input, group and side units independently");
		Collections.reverse(roster);
		Collections.reverse(groups);
		Collections.reverse(taxonomy);
		for(int i=0; i<roster.size(); i++){
			final String old=roster.get(i);
			roster.set(i, old.substring(0, old.length()-6)+"\t999\t0\t0");
		}
		write(root+"/shuffled.tsv", NcrnaConstructionMembers.GENOME_HEADER, roster);
		write(root+"/shuffled_groups.tsv", NcrnaConstructionMembers.GROUP_HEADER, groups);
		write(root+"/shuffled_taxonomy.tsv", null, taxonomy);
		NcrnaAssemblySplit.main(new String[]{"in="+root+"/shuffled.tsv", "groups="+root+"/shuffled_groups.tsv",
			"taxonomy="+root+"/shuffled_taxonomy.tsv", "out="+root+"/second", "seed=20261113", "expectedgenomes=20", "t=1"});
		for(final String suffix:new String[]{".inputs.tsv.gz", ".groups.tsv", ".strata.tsv", ".summary.tsv"}){
			check(rows(root+"/first"+suffix).equals(rows(root+"/second"+suffix)),
				"Reversed input/group/taxonomy order and different locus counts must preserve every split output via CLI");
		}
		groups.set(0, groups.get(0).replace("\tfalse\tTRAIN_CANDIDATE_ASSEMBLY_IDENTITY", "\ttrue\tASSEMBLY_SEQUENCE_IDENTITY_UNRESOLVED"));
		write(root+"/unresolved.tsv", NcrnaConstructionMembers.GROUP_HEADER, groups);
		boolean rejected=false;
		try{NcrnaAssemblySplit.run(root+"/roster.tsv", root+"/unresolved.tsv", root+"/taxonomy.tsv", root+"/bad", 20261113, 20);}
		catch(IllegalArgumentException failure){rejected=failure.getMessage().contains("Resolve every pool assembly");}
		check(rejected, "An unresolved assembly must stop the split rather than silently leak across sides");
		// The reversed first row is alias input19 in group0. Its external partner
		// is deliberately outside the20-input pool and must not become a20th group.
		groups.add("evaluation\t100\t/external.fa\t/external.fa\tNA\t90\ttrue\tfalse\tfalse\tEVALUATION\t90\ttrue\ttrue\tEVALUATION");
		write(root+"/coassigned_groups.tsv", NcrnaConstructionMembers.GROUP_HEADER, groups);
		final ArrayList<String> exception=new ArrayList<String>();
		exception.add("0\t90\tunresolved, co-assigned to prevent leakage");
		write(root+"/coassign.tsv", "pool_group\tpartner_group\treason", exception);
		NcrnaAssemblySplit.run(root+"/roster.tsv", root+"/coassigned_groups.tsv", root+"/taxonomy.tsv", root+"/allowed", 20261113, 20, root+"/coassign.tsv");
		check(rows(root+"/allowed.summary.tsv").equals(rows(root+"/first.summary.tsv")), "An external co-assignment changes no PGAP population or quota");
		check(rows(root+"/allowed.coassigned.tsv").contains("0\t90\t90\t/external.fa\tfalse\tCOVERAGE_CRITICAL_TRAIN\tunresolved, co-assigned to prevent leakage"),
			"The reviewed outside-pool partner inherits the pool group's side without an identity claim");
		// Both comparison groups may also contain pool genomes. Retain both quota
		// units and record every external alias without counting it as a third unit.
		groups.add("evaluation\t101\t/external_pair.fa\t/external_pair.fa\tNA\t1\ttrue\tfalse\tfalse\tEVALUATION\t1\ttrue\ttrue\tEVALUATION");
		write(root+"/both_pool_groups.tsv", NcrnaConstructionMembers.GROUP_HEADER, groups);
		exception.set(0, "0\t1\tunresolved, co-assigned to prevent leakage");
		write(root+"/both_pool_pair.tsv", "pool_group\tpartner_group\treason", exception);
		NcrnaAssemblySplit.run(root+"/roster.tsv", root+"/both_pool_groups.tsv", root+"/taxonomy.tsv", root+"/both", 20261113, 20, root+"/both_pool_pair.tsv");
		check(rows(root+"/both.summary.tsv").equals(rows(root+"/first.summary.tsv")), "Co-assignment preserves both distinct PGAP groups and all stratum quotas");
		check(rows(root+"/both.coassigned.tsv").contains("0\t1\t1\t/original/GCF_1.1.fa.gz\ttrue\tCOVERAGE_CRITICAL_TRAIN\tunresolved, co-assigned to prevent leakage")
			&& rows(root+"/both.coassigned.tsv").contains("0\t1\t1\t/external_pair.fa\tfalse\tCOVERAGE_CRITICAL_TRAIN\tunresolved, co-assigned to prevent leakage"),
			"A PGAP comparison partner and its external alias share the assignment without erasing the partner's pool membership");
		check(assignment.get("4").equals("TRAIN") && assignment.get("6").equals("HELD_OUT"),
			"The frozen seed's literal five-group stratum puts group4 in training and group6 in holdout");
		final ArrayList<String> conflict=rows(root+"/shuffled_groups.tsv");
		conflict.remove(0);
		for(int i=0; i<conflict.size(); i++){
			final String[] fields=conflict.get(i).split("\t", -1);
			if(fields[10].equals("4")){
				fields[12]="true";
				fields[13]="ASSEMBLY_SEQUENCE_IDENTITY_UNRESOLVED";
				conflict.set(i, String.join("\t", fields));
			}
		}
		write(root+"/conflict_groups.tsv", NcrnaConstructionMembers.GROUP_HEADER, conflict);
		exception.set(0, "4\t6\tunresolved, co-assigned to prevent leakage");
		write(root+"/conflict_pair.tsv", "pool_group\tpartner_group\treason", exception);
		boolean sideRejected=false;
		try{NcrnaAssemblySplit.run(root+"/roster.tsv", root+"/conflict_groups.tsv", root+"/taxonomy.tsv", root+"/conflict", 20261113, 20, root+"/conflict_pair.tsv");}
		catch(IllegalArgumentException failure){sideRejected=failure.getMessage().contains("share a stratum and side without altering group quotas");}
		check(sideRejected && !new File(root+"/conflict.groups.tsv").exists() && !new File(root+"/conflict.inputs.tsv.gz").exists(),
			"Opposite-side co-assignment must reach its specific guard and fail before writing any assignment, never move a group to repair quota");
		check(NcrnaAssemblySplit.stratum("d__B;p__P;c__C;o__O;f__F;g__;s__").startsWith("FAMILY:"), "Missing genus uses its explicit family");
	}
	private static void write(final String path, final String header, final ArrayList<String> rows){
		try(Output out=new Output(path)){
			if(header!=null){out.writer.println(header);}
			for(final String row:rows){out.writer.print((row+'\n').getBytes(StandardCharsets.UTF_8));}
		}
	}
	private static ArrayList<String> rows(final String path){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){rows.add(new String(row, StandardCharsets.UTF_8));}
		}
		return rows;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
}
