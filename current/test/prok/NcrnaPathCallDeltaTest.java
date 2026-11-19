package test.prok;

import java.io.File;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaPathCallDelta;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal source, multiplicity, split and coordinate-change accounting.
 * Saved miniature GFFs only; no genome calling or accuracy inference.
 * @author Ganyu
 */
public final class NcrnaPathCallDeltaTest{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPathCallDeltaTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Use native fixture flag=value options");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append,
				"Supply out=fresh_fixture_directory");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture directory must be fresh");
			run(root);
			pp.outstream.println("NCRNA_PATH_CALL_DELTA_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final File root){
		final String a="##gff-version 3\n##input_ordinal0\t0\n"
			+cds(10, 90)+rna(100, 150, "RF00174")+cds(160, 220)+rna(300, 400, "RF01482")
			+rna(500, 550, "RF00174")+rna(700, 750, "RF00174")+rna(900, 950, "RF00174")+rna(900, 950, "RF00174")
			+"##input_ordinal0\t1\n"+trna()+"##input_ordinal0\t2\n";
		final String b="##gff-version 3\n##input_ordinal0\t0\n"
			+cds(10, 90)+rna(102, 152, "RF00174")+cds(165, 220)+rna(300, 360, "RF01482")+rna(340, 400, "RF01689")
			+rna(500, 550, "RF00174")+rna(510, 540, "RF01689")+rna(800, 850, "RF00174")+rna(900, 950, "RF00174")
			+"##input_ordinal0\t1\n"+trna()+rna(100, 150, "RF00174")+"##input_ordinal0\t2\n";
		write(root, "before.gff", a);
		write(root, "after.gff", b);
		write(root, "before.config", config(root, "before.gff"));
		write(root, "after.config", config(root, "after.gff"));
		final ByteBuilder sources=new ByteBuilder("input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff\n");
		for(int i=0; i<3; i++){
			sources.append(i).append("\tassembly").append(i).append("\tid").append(i).tab().append(path(root, "id"+i+".fna.gz"))
				.tab().append(path(root, "truth"+i+".gff")).nl();
		}
		write(root, "sources.tsv", sources.toString());
		write(root, "genomes.tsv", "input_id\tcalls\tmatched\tmasked\tunmatched\nid0\t6\t3\t2\t1\nid1\t0\t0\t0\t0\nid2\t0\t0\t0\t0\n");
		write(root, "references.tsv", "input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored\n"
			+"id0\tr1\tctg\t100\t150\t+\ttrue\n"
			+"id0\tr2\tctg\t300\t400\t+\ttrue\n"
			+"id0\tlost\tctg\t700\t750\t+\ttrue\n"
			+"id0\tgained\tctg\t800\t850\t+\ttrue\n"
			+"id0\tduplicate\tctg\t900\t950\t+\tfalse\n"
			+"id1\tsame_contig_other_source\tctg\t100\t150\t+\ttrue\n");
		new NcrnaPathCallDelta(args(root, "result")).process();
		final LineParser1 clan=row(root, "result/summary.tsv", "CLAN");
		check(clan.terms()==20 && clan.parseLong(1)==6 && clan.parseLong(2)==8 && clan.parseLong(3)==2
			&& clan.parseLong(4)==4 && clan.parseLong(5)==6, "Exact differences retain all six/eight clan calls and duplicate multiplicity");
		check(clan.parseLong(6)==1 && clan.parseLong(7)==1 && clan.parseLong(10)==1 && clan.parseLong(11)==2,
			"One unique coordinate pair is separate from one lost and two gained non-overlapping calls");
		check(clan.parseLong(14)==1 && clan.parseLong(16)==1 && clan.parseLong(18)==2 && clan.parseLong(19)==1,
			"Exact/spatial duplicates and split/merge candidates remain visible, including an unchanged split anchor");
		final LineParser1 beforeOverlap=row(root, "result/clan_overlap_summary.tsv", "BEFORE");
		final LineParser1 afterOverlap=row(root, "result/clan_overlap_summary.tsv", "AFTER");
		check(beforeOverlap.parseLong(1)==6 && beforeOverlap.parseLong(2)==2 && beforeOverlap.parseLong(3)==1
			&& beforeOverlap.parseLong(4)==51 && beforeOverlap.parseLong(5)==51,
			"Two duplicate physical calls form one unordered overlap pair and affect both calls");
		check(afterOverlap.parseLong(1)==8 && afterOverlap.parseLong(2)==4 && afterOverlap.parseLong(3)==2
			&& afterOverlap.parseLong(4)==52 && afterOverlap.parseLong(5)==31 && afterOverlap.parseDouble(6)==26,
			"All-path overlap includes the unchanged split anchor, with mean distinct from maximum");
		final LineParser1 cds=row(root, "result/summary.tsv", "CDS");
		check(cds.parseLong(1)==2 && cds.parseLong(2)==2 && cds.parseLong(3)==1 && cds.parseLong(4)==1 && cds.parseLong(5)==1
			&& cds.parseLong(6)==1 && cds.parseLong(7)==1, "A changed CDS start is one unique coordinate pair, with both raw rows retained");
		final LineParser1 refs=row(root, "result/reference_summary.tsv", "ALL");
		check(refs.parseLong(1)==6 && refs.parseLong(2)==4 && refs.parseLong(3)==5 && refs.parseLong(4)==2 && refs.parseLong(5)==1,
			"Frozen reference IDs retain both coverage directions without crossing source ownership");
		try(Input input=new Input(path(root, "result/reference_coverage.tsv"))){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			boolean found=false;
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				if(p.termEquals("same_contig_other_source", 2)){
					check(p.termEquals("id1", 1) && p.parseInt(8)==0 && p.parseInt(9)==1, "Repeated contig names cannot borrow another genome's calls");
					found=true;
				}
			}
			check(found, "The repeated-name reference must actually be checked");
		}
		check(read(root, "result/neighbours.tsv.gz").contains("\tLEFT\t0\t")
			&& read(root, "result/neighbours.tsv.gz").contains("\tOVERLAP\t"), "Changed loci retain flanking and overlapping before/after feature evidence");
		check(read(root, "result/neighbours.tsv.gz").contains("\tREMOVED\t4\tAFTER\t4\tCOUNTERPART\t49\t"),
			"The unique opposite-path coordinate partner is identified explicitly, not mislabelled as a replacement neighbour");
		final String[] legacy=args(root, "legacy");
		legacy[4]="legacyfilenames=t";
		new NcrnaPathCallDelta(legacy).process();
		for(final String file:new String[]{"summary.tsv", "reference_summary.tsv", "reference_coverage.tsv", "bindings.tsv", "clan_overlap_summary.tsv"}){
			check(read(root, "result/"+file).equals(read(root, "legacy/"+file)), "Explicit historical filename binding must match the fully supplied source roster: "+file);
		}
		bad(root, "overflow", a.replace("\t100\t150\t", "\t4294967396\t150\t"));
		bad(root, "missing_input", a.replace("##input_ordinal0\t2\n", ""));
		bad(root, "bad_strand", a.replace("\t+\t0\t", "\t?\t0\t"));
		bad(root, "nonfinite_score", a.replace("\t0.0\t+", "\tNaN\t+"));
		bad(root, "overflow_score", a.replace("\t0.0\t+", "\t1e100\t+"));
		bad(root, "bad_score", a.replace("\t0.0\t+", "\t1.2.3\t+"));
		bad(root, "empty_score", a.replace("\t0.0\t+", "\t\t+"));
		bad(root, "unknown_member", a.replace("cobalamin_RF00174_", "cobalamin_UNKNOWN_"));
		bad(root, "ambiguous_member", a.replace("model:cobalamin_RF00174_fixture", "model:cobalamin_RF00174_fixture,model:cobalamin_RF01482_fixture"));
		bad(root, "wrong_member_type", a.replace("\tRNA\t", "\tCDS\t"));
		bad(root, "sequence_payload", a+"##FASTA\n>ctg\nAAAA\n");
		final String exponent="+1.25e-4";
		pair(root, "exponent", a.replace("0.0", exponent), b.replace("0.0", exponent));
		check(read(root, "result/summary.tsv").equals(read(root, "exponent/summary.tsv")),
			"Valid exponent scores retain the same feature relationships without entering numeric matching");
		nested(root);
	}
	/** A very long earlier interval must survive the index's lower-bound search;
	 * the reverse strand cannot borrow an overlapping forward-strand counterpart. */
	private static void nested(final File root){
		final String head="##gff-version 3\n##input_ordinal0\t0\n";
		final String tail="##input_ordinal0\t1\n##input_ordinal0\t2\n";
		final String common=(rna(1, 2000, "RF01482")+rna(300, 320, "RF00174")).replace("\t+\t", "\t-\t");
		final String a=head+common+rna(700, 800, "RF00174").replace("\t+\t", "\t-\t")+tail;
		final String b=head+common+rna(710, 810, "RF00174").replace("\t+\t", "\t-\t")+rna(1999, 2001, "RF00174")+tail;
		pair(root, "nested", a, b);
		final LineParser1 clan=row(root, "nested/summary.tsv", "CLAN");
		check(clan.parseLong(6)==1 && clan.parseLong(7)==1 && clan.parseLong(11)==1,
			"Nested reverse-strand coordinate pair stays unique; opposite-strand gain stays unlinked");
		check(read(root, "nested/overlap_links.tsv.gz").contains("\t3\t5\t101\tfalse\ttrue\tfalse\tfalse\n"),
			"Prefix maximum stops must retain the unchanged enclosing interval's link to the changed feature");
		final LineParser1 after=row(root, "nested/clan_overlap_summary.tsv", "AFTER");
		check(after.parseLong(2)==4 && after.parseLong(3)==3 && after.parseLong(4)==124,
			"Whole-path overlaps retain enclosing intervals and the two-base opposite-strand overlap");
	}
	private static void pair(final File root, final String name, final String before, final String after){
		write(root, name+".before.gff", before);
		write(root, name+".after.gff", after);
		write(root, name+".before.config", config(root, name+".before.gff"));
		write(root, name+".after.config", config(root, name+".after.gff"));
		final String[] args=args(root, name);
		args[0]="in1="+path(root, name+".before.gff");
		args[1]="in2="+path(root, name+".after.gff");
		args[2]="beforeconfig="+path(root, name+".before.config");
		args[3]="afterconfig="+path(root, name+".after.config");
		new NcrnaPathCallDelta(args).process();
	}
	private static void bad(final File root, final String name, final String data){
		write(root, name+".gff", data);
		write(root, name+".config", config(root, name+".gff"));
		final String[] args=args(root, "bad_"+name);
		args[0]="in1="+path(root, name+".gff");
		args[2]="beforeconfig="+path(root, name+".config");
		boolean failed=false;
		try{new NcrnaPathCallDelta(args).process();}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !new File(root, "bad_"+name+"/complete.tsv").exists(), "Broken saved data cannot publish completion: "+name);
	}
	private static String[] args(final File root, final String out){
		return new String[]{"in1="+path(root, "before.gff"), "in2="+path(root, "after.gff"),
			"beforeconfig="+path(root, "before.config"), "afterconfig="+path(root, "after.config"),
			"sources="+path(root, "sources.tsv"), "genomes="+path(root, "genomes.tsv"),
			"reference="+path(root, "references.tsv"), "expectedinputs=3", "out="+path(root, out)};
	}
	private static String config(final File root, final String output){
		return "ordered=t\ngffinputordinal=t\nin="+path(root, "id0.fna.gz")+"\nin="+path(root, "id1.fna.gz")
			+"\nin="+path(root, "id2.fna.gz")+"\nout="+path(root, output)+"\n";
	}
	private static String rna(final int start, final int stop, final String member){
		return "ctg\tBBTools\tRNA\t"+start+'\t'+stop+"\t0.0\t+\t0\tmodel:cobalamin_"+member+"_fixture\n";
	}
	private static String cds(final int start, final int stop){return "ctg\tBBTools\tCDS\t"+start+'\t'+stop+"\t0.0\t+\t0\tCDS\n";}
	private static String trna(){return "ctg\tBBTools\ttRNA\t20\t30\t0.0\t-\t0\ttRNA\n";}
	private static LineParser1 row(final File root, final String file, final String key){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path(root, file))){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				if(p.termEquals(key, 0)){return p;}
			}
		}
		throw new AssertionError("Missing literal result row: "+key);
	}
	private static String read(final File root, final String name){
		final ByteBuilder text=new ByteBuilder();
		try(Input input=new Input(path(root, name))){
			for(byte[] row; (row=input.reader.nextLine())!=null;){text.append(row).nl();}
		}
		return text.toString();
	}
	private static void write(final File root, final String name, final String text){try(Output output=new Output(path(root, name))){output.writer.print(text);}}
	private static String path(final File root, final String name){return new File(root, name).getAbsolutePath();}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
}
