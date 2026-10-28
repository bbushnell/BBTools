package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaClanEvidence;
import prok.NcrnaClanOverlay;
import prok.NcrnaClanRawEvidence;
import prok.NcrnaCmLocusAudit;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Exercises raw admission through the native join and actual core extraction.
 * @author Brian, Raiden
 */
public final class NcrnaClanRawEvidenceTest {

	public static void main(String[] args){
		require(args.length==1, "Supply a fresh fixture directory");final File dir=new File(args[0]);
		require(!dir.exists() && dir.mkdirs(), "Fixture must preserve earlier evidence");final String b=dir.getAbsolutePath();
		final String[] sources={source(0, '+'), source(1, '-'), source(2, '+'), source(3, '+')};
		write(b+"/loci.tsv", NcrnaCmLocusAudit.HEADER, sources[0], sources[1], sources[2], sources[3]);
		table(b+"/rf1.default.tbl", hit(1, "RF00174", 1, 189, 201, 320, 35), hit(3, "RF00174", 1, 189, 1, 150, 100), hit(3, "RF00174", 1, 189, 201, 350, 30));
		table(b+"/rf14.default.tbl", hit(0, "RF01482", 1, 149, 198, 275, 60), hit(0, "RF01482", 50, 60, 210, 220, 2), hit(2, "RF01482", 1, 149, 201, 350, 45));
		table(b+"/rf16.default.tbl", hit(1, "RF01689", 1, 81, 201, 281, 60));
		table(b+"/rf1.sensitive.tbl");table(b+"/rf16.sensitive.tbl");
		table(b+"/rf14.sensitive.tbl", hit(0, "RF01482", 1, 149, 201, 350, 80), hit(2, "RF01482", 1, 149, 201, 350, 60));
		write(b+"/rawmodels.tsv", "model_id\tclen\tga_bits\tstage\thits",
			"RF00174\t189\t31.5\tDEFAULT\t"+b+"/rf1.default.tbl", "RF00174\t189\t31.5\tSENSITIVE\t"+b+"/rf1.sensitive.tbl",
			"RF01482\t149\t50\tDEFAULT\t"+b+"/rf14.default.tbl", "RF01482\t149\t50\tSENSITIVE\t"+b+"/rf14.sensitive.tbl",
			"RF01689\t81\t36.5\tDEFAULT\t"+b+"/rf16.default.tbl", "RF01689\t81\t36.5\tSENSITIVE\t"+b+"/rf16.sensitive.tbl");
		raw(b, "raw");
		NcrnaClanEvidence.main(new String[]{"models="+b+"/raw.models.tsv", "sourcefamily=RF01482", "expected=4", "modelcount=3", "outprefix="+b+"/join"});
		final ByteBuilder bases=new ByteBuilder();for(int i=0; i<550; i++){bases.append("ACGTTGCA".charAt(i%8));}
		final ByteStreamWriter contexts=RrnaResourceIO.writer(b+"/contexts.fa");
		try{for(int i=0; i<4; i++){contexts.println(">RF01482_"+i+"_0 source=g"+i+" contig=chr strand="+(i==1?'-':'+')+" start=1001 stop=1150 lflank=200 rflank=200");contexts.println(bases);}}
		finally{RrnaResourceIO.finish(contexts);}
		NcrnaClanOverlay.main(new String[]{"joined="+b+"/join.loci.tsv.gz", "contexts="+b+"/contexts.fa", "outprefix="+b+"/overlay",
			"sourcefamily=RF01482", "expected=4", "minlen=109", "maxlen=187", "evidencestage=RAW_ADJUDICATED"});
		require(lines(b+"/overlay.summary.tsv").get(1).equals("RF01482\t4\t3\t2\t1\t1\t1\t109\t187\tREVIEW_REQUIRED_BEFORE_FITTING"),
			"Geometry correction, sensitive rescue and cross-label admission must conserve original loci; a strong neighbor and sub-GA full hit must not admit the fourth source");
		final ArrayList<String> cores=lines(b+"/overlay.admitted.fa"), rows=lines(b+"/overlay.overlay.tsv.gz"), selected=lines(b+"/raw.selection.tsv");
		require(cores.size()==6 && cores.get(1).equals(bases.toString().substring(197, 275))
			&& cores.get(3).equals(bases.toString().substring(200, 320)) && cores.get(5).equals(bases.toString().substring(200, 350)),
			"Actual corrected, minus-source and rescued core bases must match the chosen raw hits exactly");
		require(cores.get(0).contains("start=998 stop=1075") && cores.get(2).contains("start=1031 stop=1150"), "Both coordinate directions must survive extraction");
		require(rows.get(2).contains("\tCROSS_LABELED\ttrue\tRF00174,RF01689\tRF00174\tRAW_ADJUDICATED\t1031\t1150\t120\tWITHIN_BAND\tB12clade\tGEOMETRY-CORRECTED"),
			"Two other-only models retain both support IDs and clan membership; greater source overlap supplies coordinates despite lower bit score");
		require(rows.get(1).endsWith("\tGEOMETRY-CORRECTED") && rows.get(3).endsWith("\tSOURCE_EXACT") && rows.get(4).endsWith("\tUNRESOLVED"),
			"Truth views must expose corrected geometry separately from original exact and unresolved loci");
		boolean preserved=false;for(String row:selected){if(row.startsWith("RF01482_0_0\tRF01482\t")){
			preserved=row.contains("\t2\t1\t2\t1\t1\t1\tDEFAULT\t3\t198\t275\tPARTIAL_OVERLAP\tFULL_SUPPORT");}}
		require(preserved, "A weak additional hit must not hide full support; default/sensitive event counts and selected default coordinates stay explicit");
		table(b+"/rf16.cross_tied.tbl", hit(1, "RF01689", 1, 81, 230, 349, 60));
		final ArrayList<String> crossModels=lines(b+"/rawmodels.tsv");
		for(int i=0; i<crossModels.size(); i++){crossModels.set(i, crossModels.get(i).replace("/rf16.default.tbl", "/rf16.cross_tied.tbl"));}
		write(b+"/rawmodels_cross_tied.tsv", crossModels.toArray(new String[crossModels.size()]));raw(b, "cross_tied", "rawmodels_cross_tied.tsv");
		NcrnaClanEvidence.main(new String[]{"models="+b+"/cross_tied.models.tsv", "sourcefamily=RF01482", "expected=4", "modelcount=3", "outprefix="+b+"/cross_tied_join"});
		boolean crossFailed=false;
		try{NcrnaClanOverlay.main(new String[]{"joined="+b+"/cross_tied_join.loci.tsv.gz", "contexts="+b+"/contexts.fa", "outprefix="+b+"/cross_tied_overlay",
			"sourcefamily=RF01482", "expected=4", "minlen=109", "maxlen=187", "evidencestage=RAW_ADJUDICATED"});}
		catch(IllegalArgumentException e){crossFailed=e.getMessage().contains("Equal-overlap other-only full cores");}
		require(crossFailed && !new File(b+"/cross_tied_overlay.summary.tsv").exists(), "Actual join-to-overlay must reject equally overlapping different cores from two other-only models");
		// Equal-overlap, different qualified coordinates require a reviewed choice,
		// rather than implicitly selecting a second biological copy by bit score.
		table(b+"/rf14.tied.tbl", hit(0, "RF01482", 1, 149, 201, 275, 60), hit(0, "RF01482", 1, 149, 276, 350, 90));
		final ArrayList<String> tiedModels=lines(b+"/rawmodels.tsv");
		for(int i=0; i<tiedModels.size(); i++){tiedModels.set(i, tiedModels.get(i).replace("/rf14.default.tbl", "/rf14.tied.tbl"));}
		write(b+"/rawmodels_tied.tsv", tiedModels.toArray(new String[tiedModels.size()]));
		boolean failed=false;try{raw(b, "tied", "rawmodels_tied.tsv");}catch(IllegalArgumentException e){failed=e.getMessage().contains("Equal-overlap qualified cores");}
		require(failed && !new File(b+"/tied.summary.tsv").exists(), "Coordinate ambiguity must stop publication before any fitted population is asserted");
		System.out.println("NCRNA_CLAN_RAW_EVIDENCE_TEST_PASS");
	}
	private static void raw(String b, String out){raw(b, out, "rawmodels.tsv");}
	private static void raw(String b, String out, String manifest){
		NcrnaClanRawEvidence.main(new String[]{"loci="+b+"/loci.tsv", "rawmodels="+b+'/'+manifest, "outprefix="+b+'/'+out, "sourcefamily=RF01482", "expected=4"});
	}
	private static String source(int i, char strand){
		return "RF01482_"+i+"_0\tg"+i+"\tGCF_1.1\t1\tbacteria\tP\tchr\t"+strand+"\t1001-1150\t150\t1\tfeature\tfalse\tfalse\tANNOTATION_FULL_COPY\t2000\tfalse\t0\t200\t200\t1";
	}
	private static String hit(int i, String model, int from, int to, int start, int stop, int score){
		return "RF01482_"+i+"_0 - model "+model+" cm "+from+' '+to+' '+start+' '+stop+" + no 1 0.5 0.0 "+score+" 1e-10 ! fixture";
	}
	private static void table(String path, String... rows){
		final ByteStreamWriter w=RrnaResourceIO.writer(path);
		try{w.println("# Program: cmsearch");w.println("# Version: 1.1.5 (fixture)");for(String row:rows){w.println(row);}w.println("# [ok]");}
		finally{RrnaResourceIO.finish(w);}
	}
	private static void write(String path, String... rows){
		assert(rows.length>0):"Fixture ledgers always include the complete schema";
		final ByteStreamWriter w=RrnaResourceIO.writer(path);try{for(String row:rows){w.println(row);}}finally{RrnaResourceIO.finish(w);}
	}
	private static ArrayList<String> lines(String path){
		final ByteFile in=RrnaResourceIO.open(path);final ArrayList<String> rows=new ArrayList<String>();
		try{for(byte[] row; (row=in.nextLine())!=null; ){rows.add(new String(row, StandardCharsets.UTF_8));}}
		finally{require(!in.close(), "Fixture output read failed");}return rows;
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
