package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaClanEvidence;
import prok.NcrnaCoreRecovery;
import prok.RrnaResourceIO;

/** Checks complete source joins, conflicting model support and minus-strand projection.
 * @author Brian, Raiden
 */
public final class NcrnaClanEvidenceTest {
	public static void main(String[] args){
		require(args.length==1, "Supply a fresh fixture directory");final File dir=new File(args[0]);
		require(!dir.exists() && dir.mkdirs(), "Fixture must not replace prior evidence");final String b=dir.getAbsolutePath();
		final String[] source={source(0, 20, '+'), source(1, 30, '+'), source(2, 18, '-'), source(3, 18, '+')};
		write(b+"/rf1.tsv", NcrnaCoreRecovery.AUDIT_HEADER, source[0]+none(), source[1]+full(30, 30, 30), source[2]+none(), source[3]+none());
		write(b+"/rf14.tsv", NcrnaCoreRecovery.AUDIT_HEADER, source[0]+full(20, 20, 20), source[1]+full(20, 30, 30), source[2]+none(),
			source[3]+"\tCM_SUPPORTED_INCOMPLETE_OR_UNRESOLVED\t1\t1\t5\t20\t6\t20\tno\t45\t15\t0\t3\tHIT_WITHIN_ANNOTATION");
		write(b+"/rf16.tsv", NcrnaCoreRecovery.AUDIT_HEADER, source[0]+none(), source[1]+none(), source[2]+full(10, 10, 18), source[3]+none());
		write(b+"/models.tsv", "model_id\tclen\tga_bits\taudit", "RF00174\t30\t25\t"+b+"/rf1.tsv", "RF01482\t20\t35\t"+b+"/rf14.tsv", "RF01689\t10\t20\t"+b+"/rf16.tsv");
		run(b, "models.tsv", "good");final ArrayList<String> summary=lines(b+"/good.summary.tsv"), rows=lines(b+"/good.loci.tsv.gz"), bounds=lines(b+"/good.bounds.tsv");
		require(summary.size()==6 && rows.size()==5 && bounds.size()==2, "All source loci and all five support classes remain explicit");
		for(int i=1; i<=4; i++){require(summary.get(i).contains("\t1\t"), "Fixture must retain no-full, own-only, own-plus-other and other-only cases separately");}
		require(summary.get(5).contains("\t0\t"), "Unobserved multiple-other class must be zero, not absent");
		require(rows.get(2).contains("\tRF01482\t2\tNA\ttrue\tOWN_AND_OTHER_FULL\t"), "Two full CMs must not silently become one biological subtype");
		require(rows.get(3).contains("\tRF01482\t1\tRF01689\tfalse\tOTHER_FULL_ONLY\t") && rows.get(3).contains("\tHIT_WITHIN_ANNOTATION\t109\t118\t"), "Sole other-model support remains a proposal; reverse-strand CM coordinates project correctly");
		require(bounds.get(1).equals("RF01482\t20\t2\t20\t20\t20\t20\t30\tRF01482_0_0\tRF01482_1_0\tPROPOSED_EXTREMA_REVIEW_NOT_A_CALLER_CUTOFF"), "Only exact own-model/source-full evidence determines proposed extrema; lower-rank quantiles remain explicit");
		write(b+"/rf14none.tsv", NcrnaCoreRecovery.AUDIT_HEADER, source[0]+none(), source[1]+none(), source[2]+none(), source[3]+none());
		write(b+"/empty_models.tsv", "model_id\tclen\tga_bits\taudit", "RF00174\t30\t25\t"+b+"/rf1.tsv", "RF01482\t20\t35\t"+b+"/rf14none.tsv", "RF01689\t10\t20\t"+b+"/rf16.tsv");
		run(b, "empty_models.tsv", "empty");
		require(lines(b+"/empty.bounds.tsv").get(1).equals("RF01482\t20\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNO_EXACT_FULL_SOURCE_EVIDENCE"), "No own-model exact support must not invent a length bound or discard the original population");
		write(b+"/bad.tsv", NcrnaCoreRecovery.AUDIT_HEADER, source[0]+none(), source[1].replace("\tg1\t", "\twrong\t")+full(30, 30, 30), source[2]+none(), source[3]+none());
		write(b+"/bad_models.tsv", "model_id\tclen\tga_bits\taudit", "RF00174\t30\t25\t"+b+"/bad.tsv", "RF01482\t20\t35\t"+b+"/rf14.tsv", "RF01689\t10\t20\t"+b+"/rf16.tsv");
		boolean failed=false;try{run(b, "bad_models.tsv", "bad");}catch(IllegalArgumentException e){failed=e.getMessage().contains("changed source metadata");}
		require(failed && !new File(b+"/bad.summary.tsv").exists(), "Mismatched provenance must fail before publishing a complete clan report");
		System.out.println("NCRNA_CLAN_EVIDENCE_TEST_PASS");
	}
	private static void run(String b, String manifest, String out){
		assert(!manifest.isEmpty()):"Fixture explicitly binds all three model audits";
		NcrnaClanEvidence.main(new String[]{"models="+b+'/'+manifest, "sourcefamily=RF01482", "expected=4", "modelcount=3", "outprefix="+b+'/'+out});
	}
	private static String source(int i, int length, char strand){
		assert(length>0):"Fixture source coordinates are inclusive and nonempty";
		return "RF01482_"+i+"_0\tg"+i+"\tGCF_1.1\t1\tbacteria\tP\tchr\t"+strand+"\t101-"+(100+length)+'\t'+length+"\t1\tfeature\tfalse\tfalse\tANNOTATION_FULL_COPY\t1000\tfalse\t0\t5\t5\t1";
	}
	private static String none(){return "\tNO_REPORTED_HIT\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\t0\tNA\tNA\tNONE";}
	private static String full(int modelLength, int hitLength, int sourceLength){
		assert(hitLength<=sourceLength):"Fixture full-model hits are contained within their source annotation";
		return "\tCM_FULL_SPAN_SUPPORTED\t1\t1\t1\t"+modelLength+"\t6\t"+(5+hitLength)+"\tno\t60\t"+hitLength+"\t0\t"+(sourceLength-hitLength)+'\t'+(hitLength==sourceLength?"EXACT":"HIT_WITHIN_ANNOTATION");
	}
	private static void write(String path, String... rows){
		assert(rows.length>0):"Every fixture file includes its schema";
		final ByteStreamWriter w=RrnaResourceIO.writer(path);try{for(String row:rows){w.println(row);}}finally{RrnaResourceIO.finish(w);}
	}
	private static ArrayList<String> lines(String path){
		final ByteFile in=RrnaResourceIO.open(path);final ArrayList<String> rows=new ArrayList<String>();
		try{for(byte[] row; (row=in.nextLine())!=null; ){rows.add(new String(row, StandardCharsets.UTF_8));}}
		finally{require(!in.close(), "Fixture output read failed");}return rows;
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
