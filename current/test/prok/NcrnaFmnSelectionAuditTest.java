package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaReference;
import prok.NcrnaConstructionMembers;
import prok.NcrnaFmnSelectionAudit;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Original-path intersection, species overlap and independent native inclusion.
 * @author Ganyu
 */
public final class NcrnaFmnSelectionAuditTest{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaFmnSelectionAuditTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Unknown selection test option");
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append, "Use out=fresh_fixture_directory");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Selection fixture directory must be fresh");
			fixture(root.getPath());
			pp.outstream.println("NCRNA_FMN_SELECTION_AUDIT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private static void fixture(final String root){
		String roster=NcrnaConstructionMembers.GENOME_HEADER+'\n';
		for(int i=0; i<5; i++){roster+=i+"\tg"+i+"\tGCF_"+i+".1\t1\tbacteria\tBacteroidota\t/gff"+i+"\t/pool/"+i+".fa\t0\t0\t0\n";}
		write(root+"/roster.tsv", roster);
		write(root+"/catalog.tsv", "input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records\n"
			+"1\t/pool/0.fa\t1\t1\t1\n2\t/pool/1.fa\t1\t1\t1\n3\t/external.fa\t1\t1\t1\n");
		String lineages="cohort\tinput_id\tfasta\tassembly_group\tassembly\tclass\torder\tfamily\tgenus\ttaxonomy_status\tCM_admitted_cores\n";
		final int[] accessions={0, 1, 5, 2, 3};
		final String[] species={"A", "B", "External", "B", "A"};
		String taxonomy="";
		for(int i=0; i<5; i++){
			final String cohort=i<3 ? "evaluation" : "discovery";
			final String path=i==2 ? "/external.fa" : "/pool/"+accessions[i]+".fa";
			final String genus=i==2 ? "Other" : "Chryseobacterium";
			lineages+=cohort+'\t'+(i+1)+'\t'+path+'\t'+i+"\tGCF_"+accessions[i]+".1\tC\tO\tF\t"+genus+"\tGTDB_ACCESSION_MATCH\t"+(i==0 || i==2 ? 1 : 0)+'\n';
			taxonomy+="RS_GCF_"+accessions[i]+".1\td__Bacteria;p__Bacteroidota;c__C;o__O;f__F;g__"+genus+";s__"+species[i]+'\n';
		}
		write(root+"/lineages.tsv", lineages);
		write(root+"/taxonomy.tsv", taxonomy);
		final String raw="# fixture\nchrA - FMN RF00050 cm 1 138 10 147 + no 1 0.4 0 70 1e-8 ! fixture\n"
			+"# Option settings: cmsearch --cut_ga model target\n# [ok]\n";
		write(root+"/native.tbl", raw);
		final String hits=BuildNcrnaReference.HIT_HEADER+'\n'+"1\t"+root+"/native.tbl\t2\tchrA\t10\t147\t+\tRF00050\t15.1\tFMN\tTARGET_CLASS\t1\t138\t138\t0\t70\t1e-8\tfalse\tVALID_INTERVAL\tfalse\tfalse\n";
		write(root+"/hits.tsv", hits);
		final NcrnaFmnSelectionAudit.Config c=new NcrnaFmnSelectionAudit.Config();
		c.roster=root+"/roster.tsv";
		c.catalog=root+"/catalog.tsv";
		c.lineages=root+"/lineages.tsv";
		c.taxonomy=root+"/taxonomy.tsv";
		c.hits=root+"/hits.tsv";
		c.output=root+"/audit";
		c.pool=5;
		c.evaluation=3;
		c.overlap=2;
		c.labels=5;
		c.hitRows=1;
		NcrnaFmnSelectionAudit.run(c);
		check(rows(c.output+".summary.tsv").get(1).equals("5\t3\t2\t3\t2\t1\t2\t1\t1\t1\tSAVED_SELECTION_DIAGNOSTIC_NOT_SPLIT"),
			"Original PGAP intersection excludes the external evaluation path; native ! overrides the misleading legacy printed-score boolean");
		check(rows(c.output+".species.tsv").contains("s__B\t0\t1\t1\ttrue"), "The negative evaluation group shares its species with one discovery group");
		check(rows(c.output+".witnesses.tsv").size()==2, "The counted positive must retain its exact original native witness");
		NcrnaFmnSelectionAudit.main(new String[]{"in="+c.roster, "catalog="+c.catalog, "lineages="+c.lineages,
			"taxonomy="+c.taxonomy, "hits="+c.hits, "out="+root+"/cli", "expectedgenomes=5", "expectedevaluation=3",
			"expectedoverlap=2", "expectedlineages=5", "expectedhits=1", "t=1"});
		check(rows(root+"/cli.summary.tsv").equals(rows(c.output+".summary.tsv")), "Actual CLI must preserve literal diagnostic counts");
		write(root+"/excluded.tbl", raw.replace(" ! fixture", " ? fixture"));
		write(root+"/bad_hits.tsv", hits.replace("/native.tbl", "/excluded.tbl"));
		c.hits=root+"/bad_hits.tsv";
		c.output=root+"/bad";
		boolean rejected=false;
		try{NcrnaFmnSelectionAudit.run(c);}
		catch(IllegalArgumentException failure){rejected=failure.getMessage().contains("physical row");}
		check(rejected, "Native exclusion cannot be converted to GA inclusion by a saved imported row");
	}
	private static void write(final String path, final String value){try(Output out=new Output(path)){out.writer.print(value.getBytes(StandardCharsets.UTF_8));}}
	private static ArrayList<String> rows(final String path){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path)){for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){rows.add(new String(row, StandardCharsets.UTF_8));}}
		return rows;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
}
