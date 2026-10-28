package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaModelMembershipReport;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Checks inclusive purity threshold, separate added membership and population conservation.
 * @author Brian, Raiden
 */
public final class NcrnaModelMembershipReportTest {

	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final ByteBuilder fasta=new ByteBuilder(), assignments=new ByteBuilder();
		assignments.append("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n");
		for(int i=0; i<5; i++){
			final String header="q"+i+" tid_1 member_family="+(i==4?"B12clade":"RF00174")+" build_group=RF00174";
			fasta.append('>').append(header).nl().append(bases(i==4?'G':'A', 20+i)).nl();
			assignments.append(header).append("\t1\tunknown\tRF00174_kept_base\tcobalamin_RF00174_base n=5\tkept\n");
		}
		write(root.resolve("input.fa"), fasta.toString());write(root.resolve("assignments.tsv"), assignments.toString());
		final String ref=">cobalamin_RF00174_base n=5\n"+bases('A', 20)+"\n>cobalamin_RF00174_empty\n"+bases('T', 25)+
			"\n>cobalamin_RF00174_extra_q4 n=1\n"+bases('G', 24)+'\n';
		write(root.resolve("ref.fa"), ref);
		final String head="locus_id\tmodel_id\temitted_family\tmember_family\tlength\tconstruction\n";
		final String added="q4\tcobalamin_RF00174_extra_q4\tRF00174\tB12clade\t24\tEXACT_SINGLE_MEMBER_NO_REASSIGNMENT\n";
		write(root.resolve("additions.tsv"), head+added);
		run(root, "input.fa", "ref.fa", "additions.tsv", "good");
		final ByteFile rows=RrnaResourceIO.open(root.resolve("good.models.tsv").toString());final LineParser1 p=new LineParser1('\t');
		try{
			rows.nextLine();p.set(rows.nextLine());
			check(p.terms()==19 && p.termEquals("cobalamin_RF00174_base", 0) && p.parseInt(3)==5 && p.parseInt(4)==20 && p.parseInt(8)==24
				&& p.termEquals("0.800000", 9) && p.termEquals("PURE", 10) && p.parseInt(11)==4 && p.parseInt(17)==1,
				"Four own plus one clan-only member must report0.800000/PURE and assigned bounds20..24; fraction="+p.parseString(9)+" purity="+p.parseString(10)+" fields="+p.terms());
			p.set(rows.nextLine());check(p.parseInt(3)==0 && p.termEquals("NA", 4) && p.termEquals("NA", 9) && p.termEquals("EMPTY", 10),
				"A model with no documented members must expose missing composition and bounds, not pretend purity");
			p.set(rows.nextLine());check(p.parseInt(3)==1 && p.parseInt(4)==24 && p.parseInt(8)==24 && p.parseDouble(9)==0 && p.termEquals("MIXED", 10)
				&& p.parseInt(17)==1 && p.parseDouble(18)==1 && rows.nextLine()==null,
				"Added clan-only singleton contributes to its new model without changing its biological label");
		}finally{check(!rows.close(), "Model report read failed");}
		final ByteFile summary=RrnaResourceIO.open(root.resolve("good.summary.tsv").toString());
		try{
			summary.nextLine();p.set(summary.nextLine());check(p.terms()==17 && p.parseInt(1)==5 && p.parseInt(2)==5 && p.parseInt(5)==1
				&& p.parseInt(6)==6 && p.parseInt(7)==3 && p.parseInt(8)==1 && p.parseInt(9)==1 && p.parseInt(10)==1
				&& p.parseInt(12)==20 && p.parseInt(14)==22 && p.parseInt(16)==24,
				"Five unique loci create six model-membership events; the summary must keep those units separate");
		}finally{check(!summary.close(), "Summary read failed");}
		write(root.resolve("duplicate.tsv"), head+added+added);reject(root, "input.fa", "ref.fa", "duplicate.tsv", "duplicate");
		write(root.resolve("relabel.tsv"), head+added.replace("\tB12clade\t", "\tRF00174\t"));reject(root, "input.fa", "ref.fa", "relabel.tsv", "relabel");
		write(root.resolve("changedref.fa"), ref.replace(bases('G', 24), bases('C', 24)));reject(root, "input.fa", "changedref.fa", "additions.tsv", "changed");
		final String complete=fasta.toString();write(root.resolve("missing.fa"), complete.substring(0, complete.indexOf(">q4 ")));
		reject(root, "missing.fa", "ref.fa", "additions.tsv", "missing");
		System.out.println("NCRNA_MODEL_MEMBERSHIP_REPORT_TEST_PASS");
	}
	private static void run(Path root, String input, String ref, String additions, String out){
		NcrnaModelMembershipReport.main(new String[]{"in="+root.resolve(input), "ref="+root.resolve(ref), "assignments="+root.resolve("assignments.tsv"),
			"additions="+root.resolve(additions), "outprefix="+root.resolve(out), "family=RF00174", "labels=RF00174,RF01482,RF01689,B12clade",
			"modelprefix=cobalamin_RF00174_", "expected=5", "purity=0.8"});
	}
	private static void reject(Path root, String input, String ref, String additions, String out){
		boolean failed=false;try{run(root, input, ref, additions, out);}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve(out+".models.tsv")), "Invalid contributor provenance must fail before publishing a report: "+out);
	}
	private static String bases(char base, int length){
		check(length>0, "Fixture sequences need positive lengths");final ByteBuilder b=new ByteBuilder();for(int i=0; i<length; i++){b.append(base);}return b.toString();
	}
	private static void write(Path path, String text){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
