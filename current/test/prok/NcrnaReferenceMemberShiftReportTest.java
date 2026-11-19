package test.prok;

import java.io.File;
import java.util.ArrayList;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaReferenceMemberShiftReport;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal source/strand/CM and missing-key oracles; no biological result implied.
 * @author Brian, Ganyu
 */
public final class NcrnaReferenceMemberShiftReportTest{

	public static void main(final String[] args){
		final PreParser pre=new PreParser(args, NcrnaReferenceMemberShiftReportTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pre.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(java.util.Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				require(parser.parse(arg, a, b), "Unknown test argument");
			}
			require(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Require out=fresh_directory without filtering or overwrite");
			final File root=new File(parser.out1).getAbsoluteFile();
			require(root.mkdir(), "Test directory must be fresh");
			run(root);
			pre.outstream.println("NCRNA_REFERENCE_MEMBER_SHIFT_REPORT_TEST_PASS");
		}finally{Shared.closeStream(pre.outstream);}
	}

	private static void run(final File root){
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ArrayList<Row> left=new ArrayList<Row>(), right=new ArrayList<Row>();
		left.add(row(codec, "left", "same", '+', 100, 199, "RF00174"));
		left.add(codec.create("left", new Locus("second", "cobalamin", 300, 399, '+'), false, null,
			"RF00174", "literal_fixture", codec.scores("RF01482,RF00174", new long[]{1000,500}, new String[]{"10.00","10.00"})));
		left.add(row(codec, "left", "near", '+', 701, 800, "RF01482"));//Overlap is deliberately not an exact key.
		left.add(row(codec, "left", "third", '-', 500, 599, "RF00174,RF01689"));
		right.add(row(codec, "right", "same", '-', 100, 199, "RF00174"));
		right.add(row(codec, "right", "second", '+', 300, 399, "RF01482"));//Same contig/bounds as left; different source.
		right.add(row(codec, "right", "third", '+', 900, 999, null));
		final String leftPath=new File(root, "left.gff").getPath(), rightPath=new File(root, "right.gff").getPath();
		codec.write(leftPath, "left", left);
		codec.write(rightPath, "right", right);
		final String bindings=write(root, "sources.tsv", "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff\n"
			+"0\tassembly_left\tleft\t/unused.fa\t"+leftPath+"\n"
			+"1\tassembly_right\tright\t/unused.fa\t"+rightPath+"\n"
			+"2\tassembly_omit\tomit\t/unused.fa\t/must_not_read_excluded.gff\n");
		final String auditText="# HELD OUT (panel genomes excluded from tables and training): exclusion bindings; not a fitted result\n"
			+"panel\trow1\tfasta\tassembly\tstatus\n/panel.txt\t1\t/omitted.fa\tassembly_omit\tEXCLUDED_ASSEMBLY\n";
		final String audit=write(root, "exclude.tsv", auditText);
		final String matchText="#family\tseqid\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3\n"
			+"cobalamin\t4:left/same\t+\t100\t199\t100\t199\t0\t0\n" //Clan alias is not an extra vote.
			+"RF00174\t4:left/same\t+\t100\t199\t100\t199\t0\t0\n"
			+"RF00174\t5:right/same\t-\t100\t199\t90\t189\t10\t10\n"
			+"RF01482\t4:left/second\t+\t300\t399\t319\t428\t19\t29\n"
			+"RF01482\t5:right/second\t+\t300\t399\t300\t399\t0\t0\n"
			+"RF01482\t4:left/near\t+\t700\t799\t698\t795\t-2\t-4\n"
			+"RF01482\t4:omit/second\t+\t300\t399\t319\t428\t19\t29\n"
			+"RF01689\t4:left/third\t-\t500\t599\t493\t612\t-13\t7\n"
			+"RF01689\t5:right/third\t+\t900\t999\t904\t1005\t4\t6\n";
		final String matches=write(root, "matches.tsv", matchText);
		final String out=new File(root, "result").getPath();
		execute(matches, bindings, audit, out, "2,4,2", "2,3,2");
		final LineParser1 cross=find(out+".summary.tsv", "RF01482", "RF00174", "5prime");
		require(cross.terms()==19 && cross.parseLong(3)==1 && cross.parseInt(9)==19 && cross.parseLong(4)==0,
			"Recorded RF00174 attribution retains +19 even though RF01482 is first and higher-scoring in cm_models");
		require(find(out+".summary.tsv", "RF01482", "RF00174", "3prime").parseInt(9)==29,
			"Same cross-model match retains its literal +29 stop");
		final LineParser1 own=find(out+".summary.tsv", "RF01482", "RF01482", "5prime");
		require(own.parseLong(3)==1 && own.parseInt(9)==0 && own.parseDouble(5)==100,
			"Identical contig/bounds in the other source join that source's own model, never the first source");
		require(find(out+".summary.tsv", "RF01482", "UNRESOLVED_NO_EXACT_REFERENCE", "5prime").parseLong(3)==1,
			"One-base shifted overlapping reference must remain unresolved");
		require(find(out+".summary.tsv", "RF01689", "AMBIGUOUS_MODEL", "5prime").parseInt(9)==-13,
			"Comma-valued model provenance remains ambiguous, with the reverse-strand error unchanged");
		require(find(out+".summary.tsv", "RF01689", "UNRESOLVED_MODEL", "3prime").parseInt(9)==6,
			"Missing reference model remains an explicit unresolved cohort");
		final LineParser1 main=find(out+".summary.tsv", "RF00174", "RF00174", "5prime");
		require(main.parseLong(3)==2 && main.parseLong(4)==2 && main.parseInt(9)==0 && main.parseInt(13)==10
			&& main.parseInt(14)==10 && main.parseInt(15)==0 && main.parseLong(16)==1 && main.parseInt(18)==2,
			"Literal values0,10 give nearest-rank median0, central80/IQR10, tied modes and both within ten");
		try(Input in=new Input(out+".summary.tsv")){
			final LineParser1 p=new LineParser1('\t');
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				if(p.termEquals("query_member", 0)){continue;}
				require(p.terms()==19 && p.parseLong(3)>0, "Only existing nonempty distributions are reported");
				require(!(p.termEquals("RF00174", 0) && p.termEquals("RF01482", 1)),
					"An empty cohort has no summary row; its nonexistent statistics are not encoded as placeholders");
				for(int field=3; field<19; field++){require(Double.isFinite(p.parseDouble(field)), "Every reported numeric statistic must be defined and finite");}
			}
		}
		try(Input in=new Input(out+".population.tsv")){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(int member=0; member<3; member++){
				p.set(in.reader.nextLine());
				require(p.terms()==10 && p.parseInt(4)==3 && p.parseInt(5)==1 && p.parseInt(6)==1
					&& p.parseInt(7)==1 && p.parseInt(8)==0 && p.parseInt(9)==2,
					"One excluded assembly is never read; only two reference files serve the retained population");
				if(member==1){require(p.parseInt(1)==4 && p.parseInt(2)==3 && p.parseInt(3)==1, "Second-member selected=retained+excluded");}
			}
			require(in.reader.nextLine()==null, "Population report has exactly three declared members");
		}
		int loci=0, evidence=0;
		try(Input in=new Input(out+".loci.tsv")){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==18, "The locus ledger must include recorded attribution and reference_evidence");
				if(p.termEquals("literal_fixture", 15)){evidence++;}
				else{require(p.termEquals(".", 15) && p.termEquals("UNRESOLVED_NO_EXACT_REFERENCE", 11), "Only the missing exact reference lacks evidence in this fixture");}
				loci++;
			}
		}
		require(loci==7 && evidence==6, "All retained loci and every available evidence pointer must survive, including ambiguous/missing models");
		final String badOffset=write(root, "bad_offset.tsv", matchText.replace("319\t428\t19\t29", "319\t428\t18\t29"));
		fails(()->execute(badOffset, bindings, audit, new File(root, "bad_offset").getPath(), "2,4,2", "2,3,2"));
		final String duplicate=write(root, "duplicate.tsv", matchText+"RF01482\t4:left/same\t+\t100\t199\t100\t199\t0\t0\n");
		fails(()->execute(duplicate, bindings, audit, new File(root, "duplicate").getPath(), "2,5,2", "2,4,2"));
		final String wrongIdentity=write(root, "wrong_source.tsv", read(bindings).replace(leftPath, rightPath));
		fails(()->execute(matches, wrongIdentity, audit, new File(root, "wrong_source").getPath(), "2,4,2", "2,3,2"));
		final String badAudit=write(root, "bad_audit.tsv", auditText.replace("/panel.txt\t1", "/panel.txt\t2"));
		fails(()->execute(matches, bindings, badAudit, new File(root, "bad_audit").getPath(), "2,4,2", "2,3,2"));
		fails(()->execute(matches, bindings, audit, new File(root, "bad_counts").getPath(), "2,4,2", "2,4,2"));
		require(!new File(root, "bad_counts.summary.tsv").exists() && !new File(root, "wrong_source.summary.tsv").exists(),
			"Malformed population or source identity fails before any output is published");
	}

	private static Row row(final NcrnaReferenceGff codec, final String source, final String contig, final char strand,
			final int start, final int stop, final String model){
		return codec.create(source, new Locus(contig, "cobalamin", start, stop, strand), false, null, model, "literal_fixture");
	}
	private static void execute(final String matches, final String bindings, final String audit, final String out,
			final String before, final String after){
		new NcrnaReferenceMemberShiftReport(new String[]{"in="+matches, "bindings="+bindings, "exclude="+audit,
			"out="+out, "expectedall="+before, "expected="+after}).process();
	}
	private static LineParser1 find(final String path, final String member, final String group, final String end){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				if(p.termEquals(member, 0) && p.termEquals(group, 1) && p.termEquals(end, 2)){
					return new LineParser1('\t').set(row.clone());
				}
			}
		}
		throw new AssertionError("Missing literal summary row: "+member+"/"+group+"/"+end);
	}
	private static String write(final File root, final String name, final String text){
		final String path=new File(root, name).getPath();
		try(Output out=new Output(path)){out.writer.print(text.getBytes(java.nio.charset.StandardCharsets.UTF_8));}
		return path;
	}
	private static String read(final String path){
		final structures.ByteBuilder out=new structures.ByteBuilder();
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){out.append(row).nl();}
		}
		return out.toString();
	}
	private static void fails(final Runnable action){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		require(failed, "Invalid saved evidence must fail loudly");
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
