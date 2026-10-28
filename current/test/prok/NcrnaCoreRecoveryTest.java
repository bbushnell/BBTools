package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaCoreRecovery;
import prok.RrnaResourceIO;
import stream.Read;

/** Tests sense-coordinate conversion, ambiguity location and unresolved-row conservation.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCoreRecoveryTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final StringBuilder sequence=new StringBuilder(50);for(int i=0; i<50; i++){sequence.append("ACGT".charAt(i&3));}
		final String seq=sequence.toString();
		check(seq.length()==50, "Fixture is40source bases plus5flanks on each side");
		final ByteStreamWriter audit=RrnaResourceIO.writer(root.resolve("audit.tsv").toString()), contexts=RrnaResourceIO.writer(root.resolve("contexts.fa").toString());
		try{
			audit.println(NcrnaCoreRecovery.AUDIT_HEADER);
			for(int i=0; i<6; i++){
				final String id="RF00174_"+i, strand=i==1 || i==5?"-":"+", state=i==4?"CM_SUPPORTED_INCOMPLETE_OR_UNRESOLVED":"CM_FULL_SPAN_SUPPORTED";
				audit.println(id+"\tg"+i+"\tGCF_1.1\t1\tbacteria\tP\tcontig\t"+strand+"\t101-140\t40\t1\tf\tfalse\tfalse\tANNOTATION_FULL_COPY\t1000\tfalse\t"+(i==2 || i==3?1:0)+"\t5\t5\t1\t"+state+
					'\t'+(i==5?2:1)+"\t1\t"+(i==4?2:1)+(i==5?"\t189\t5\t45\tno\t50\t40\t-1\t0\tANNOTATION_WITHIN_HIT":"\t189\t9\t28\tno\t50\t20\t3\t17\tHIT_WITHIN_ANNOTATION"));
				final String bases=i==2?seq.substring(0, 6)+'N'+seq.substring(7):i==3?seq.substring(0, 10)+'N'+seq.substring(11):seq;
				contexts.println('>'+id+" source=g"+i+" contig=contig strand="+strand+" start=101 stop=140 lflank=5 rflank=5\n"+bases);
			}
		}finally{RrnaResourceIO.finish(audit, contexts);}
		NcrnaCoreRecovery.main(new String[]{"audit="+root.resolve("audit.tsv"), "contexts="+root.resolve("contexts.fa"), "outprefix="+root.resolve("recovery"),
			"expected=6", "minsource=40", "clen=189", "gabits=31.5"});
		final ArrayList<Read> cores=new ArrayList<Read>();RrnaResourceIO.stream(root.resolve("recovery.cores.fa").toString(), cores::add);
		check(cores.size()==4, "Plus, minus, outside-core ambiguity and supported one-base extension are recoverable");
		for(int i=0; i<cores.size(); i++){check(new String(cores.get(i).bases, java.nio.charset.StandardCharsets.US_ASCII).equals(i==3?seq.substring(4, 45):seq.substring(8, 28)), "Sense contexts must be sliced exactly once, never reverse-complemented again");}
		check(RrnaResourceIO.field(cores.get(3).id, "start").equals("101") && RrnaResourceIO.field(cores.get(3).id, "stop").equals("141") && cores.get(3).length()==41,
			"A minus-strand full-model hit with one overlapping and one neighboring hit must extend genomic stop, preserving the41bp sense slice");
		check(RrnaResourceIO.field(cores.get(0).id, "start").equals("104") && RrnaResourceIO.field(cores.get(0).id, "stop").equals("123"), "Plus-strand genomic core must reflect one-based inclusive context coordinates");
		check(RrnaResourceIO.field(cores.get(1).id, "start").equals("118") && RrnaResourceIO.field(cores.get(1).id, "stop").equals("137"), "Minus-strand conversion must reverse coordinate offsets, not sequence orientation");
		final ByteFile in=RrnaResourceIO.open(root.resolve("recovery.ledger.tsv").toString());final LineParser1 p=new LineParser1('\t');int rows=0;
		try{
			in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);check(p.terms()==38, "Recovery must preserve all original audit fields plus disposition");
				if(rows==3){check(p.termEquals("AMBIGUOUS_CORE", 34), "Ambiguity inside core must exclude the core");}
				if(rows==4){check(p.termEquals("CM_SUPPORTED_INCOMPLETE_OR_UNRESOLVED", 34), "Incomplete-model support must remain unresolved");}rows++;
			}
		}finally{check(!in.close(), "Recovery fixture ledger read failed");}
		check(rows==6, "Rejected and unresolved source rows must remain in the ledger");
		final ByteStreamWriter changed=RrnaResourceIO.writer(root.resolve("changed.fa").toString());
		try{changed.println(">RF00174_0 source=g0 contig=contig strand=+ start=102 stop=140 lflank=5 rflank=5\n"+seq);}
		finally{RrnaResourceIO.finish(changed);}
		boolean failed=false;
		try{NcrnaCoreRecovery.main(new String[]{"audit="+root.resolve("audit.tsv"), "contexts="+root.resolve("changed.fa"), "outprefix="+root.resolve("changed"),
			"expected=6", "minsource=40", "clen=189", "gabits=31.5"});}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Changed source coordinates must fail before a plausible core can be emitted");
		System.out.println("NCRNA_CORE_RECOVERY_TEST_PASS");
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
