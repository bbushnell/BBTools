package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaClanEvidence;
import prok.NcrnaClanOverlay;
import prok.NcrnaCoreRecovery;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Checks clan admission, source preservation, strand projection and targeted rescue.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaClanOverlayTest {

	public static void main(String[] args){
		require(args.length==1, "Supply a fresh fixture directory");final File dir=new File(args[0]);
		require(!dir.exists() && dir.mkdirs(), "Fixture must preserve earlier evidence");final String b=dir.getAbsolutePath();
		final String[] sources={source(0, '+'), source(1, '+'), source(2, '-'), source(3, '+')};
		write(b+"/rf1.tsv", NcrnaCoreRecovery.AUDIT_HEADER, sources[0]+none(), sources[1]+full(189, 150), sources[2]+none(), sources[3]+none());
		write(b+"/rf14.tsv", NcrnaCoreRecovery.AUDIT_HEADER, sources[0]+full(149, 150), sources[1]+full(149, 150), sources[2]+none(), sources[3]+none());
		write(b+"/rf16.tsv", NcrnaCoreRecovery.AUDIT_HEADER, sources[0]+none(), sources[1]+none(), sources[2]+full(81, 81), sources[3]+none());
		write(b+"/models.tsv", "model_id\tclen\tga_bits\taudit", "RF00174\t189\t31.5\t"+b+"/rf1.tsv", "RF01482\t149\t50\t"+b+"/rf14.tsv", "RF01689\t81\t36.5\t"+b+"/rf16.tsv");
		NcrnaClanEvidence.main(new String[]{"models="+b+"/models.tsv", "sourcefamily=RF01482", "expected=4", "modelcount=3", "outprefix="+b+"/join"});
		final ByteStreamWriter contexts=RrnaResourceIO.writer(b+"/contexts.fa");final ByteBuilder sequence=new ByteBuilder();
		for(int i=0; i<170; i++){sequence.append("ACGTTGCA".charAt(i%8));}
		try{for(int i=0; i<4; i++){contexts.println(">RF01482_"+i+"_0 source=g"+i+" contig=chr strand="+(i==2?'-':'+')+" start=1001 stop=1150 lflank=10 rflank=10");contexts.print(sequence);contexts.println();}}
		finally{RrnaResourceIO.finish(contexts);}
		run(b, "first", null);
		require(lines(b+"/first.summary.tsv").get(1).equals("RF01482\t4\t3\t2\t1\t1\t1\t109\t187\tREVIEW_REQUIRED_BEFORE_FITTING"),
			"Full-supported cross-label core remains admitted outside the diagnostic band; unresolved stays distinct");
		final ArrayList<String> overlay=lines(b+"/first.overlay.tsv.gz"), joined=lines(b+"/join.loci.tsv.gz");
		require(overlay.size()==5, "Overlay must conserve all original source rows");
		boolean minus=false;for(String row:overlay){if(row.startsWith("RF01482_2_0\t")){minus=row.contains("\tRF01482\tCROSS_LABELED\ttrue\tRF01689\tRF01689\tFIRST_PASS\t1070\t1150\t81\tOUTSIDE_BAND_SUPPORTED\tB12clade");}}
		require(minus, "Cross-labeled minus-strand core must retain original source and genomic coordinates");
		final ArrayList<String> cores=lines(b+"/first.admitted.fa"), followup=lines(b+"/first.followup.contexts.fa");
		final String bases=sequence.toString();require(cores.size()==6, "Exactly three admitted core records must be emitted");
		for(int i=0; i<3; i++){
			final boolean cross=i==2;
			final String expectedHeader=">RF01482_"+i+"_0 source_family=RF01482 member_family="+(cross?"B12clade":"RF01482")
				+" source=g"+i+" contig=chr strand="+(cross?'-':'+')+" start="+(cross?1070:1001)+" stop=1150 coordinate_model="
				+(cross?"RF01689":"RF01482")+" lflank=0 rflank=0";
			require(cores.get(2*i).equals(expectedHeader) && cores.get(2*i+1).equals(bases.substring(10, cross?91:160)),
				"Supported core header and exact sense-oriented slice must agree with selected CM coordinates, including the minus-source81bp core");
		}
		require(lines(b+"/first.followup.loci.tsv").size()==2 && followup.size()==2
			&& followup.get(0).equals(">RF01482_3_0 source=g3 contig=chr strand=+ start=1001 stop=1150 lflank=10 rflank=10") && followup.get(1).equals(bases),
			"Followup contains exactly the unresolved original context, including flanks");
		// Reuse an actual native join row with identical source geometry, changing only
		// the complete source identity to the unresolved locus for a rescue fixture.
		String rescued=joined.get(1).replace("RF01482_0_0\tg0\t", "RF01482_3_0\tg3\t");
		write(b+"/rescue.tsv", joined.get(0), rescued);run(b, "rescued", b+"/rescue.tsv");
		require(lines(b+"/rescued.summary.tsv").get(1).equals("RF01482\t4\t4\t3\t1\t0\t1\t109\t187\tREVIEW_REQUIRED_BEFORE_FITTING"),
			"Rescue may add supported loci without changing pre-existing admissions or dropping source rows");
		require(new File(b+"/rescued.followup.contexts.fa").length()==0 && lines(b+"/rescued.followup.loci.tsv").size()==1,
			"Empty unresolved cohort keeps the ledger header and has no invented contexts");
		write(b+"/badrescue.tsv", joined.get(0), joined.get(1));boolean failed=false;
		try{run(b, "bad", b+"/badrescue.tsv");}catch(IllegalArgumentException e){failed=e.getMessage().contains("originally unresolved");}
		require(failed && !new File(b+"/bad.summary.tsv").exists(), "A rescue cannot replace an already admitted source");
		final ArrayList<String> originalContexts=lines(b+"/contexts.fa");
		originalContexts.set(1, bases.substring(0, 15)+'N'+bases.substring(16));
		write(b+"/ambiguous.fa", originalContexts.toArray(new String[originalContexts.size()]));failed=false;
		try{run(b, "ambiguous", null, "ambiguous.fa");}catch(IllegalArgumentException e){failed=e.getMessage().contains("Ambiguous supported core");}
		require(failed && !new File(b+"/ambiguous.summary.tsv").exists(), "An N inside a supported core must fail instead of publishing a successful fitting population");
		System.out.println("NCRNA_CLAN_OVERLAY_TEST_PASS");
	}
	private static void run(String b, String out, String rescue){run(b, out, rescue, "contexts.fa");}
	private static void run(String b, String out, String rescue, String contexts){
		final ArrayList<String> args=new ArrayList<String>();
		args.add("joined="+b+"/join.loci.tsv.gz");args.add("contexts="+b+'/'+contexts);args.add("outprefix="+b+'/'+out);
		args.add("sourcefamily=RF01482");args.add("expected=4");args.add("minlen=109");args.add("maxlen=187");if(rescue!=null){args.add("rescue="+rescue);}
		NcrnaClanOverlay.main(args.toArray(new String[args.size()]));
	}
	private static String source(int i, char strand){
		return "RF01482_"+i+"_0\tg"+i+"\tGCF_1.1\t1\tbacteria\tP\tchr\t"+strand+"\t1001-1150\t150\t1\tfeature\tfalse\tfalse\tANNOTATION_FULL_COPY\t2000\tfalse\t0\t10\t10\t1";
	}
	private static String none(){return "\tNO_REPORTED_HIT\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\t0\tNA\tNA\tNONE";}
	private static String full(int clen, int length){
		assert(length<=150):"Fixture CM hit is within the source interval";
		return "\tCM_FULL_SPAN_SUPPORTED\t1\t1\t1\t"+clen+"\t11\t"+(10+length)+"\tno\t70\t"+length+"\t0\t"+(150-length)+'\t'+(length==150?"EXACT":"HIT_WITHIN_ANNOTATION");
	}
	private static void write(String path, String... rows){
		assert(rows.length>0):"Every fixture ledger has a schema";
		final ByteStreamWriter w=RrnaResourceIO.writer(path);try{for(String row:rows){w.println(row);}}finally{RrnaResourceIO.finish(w);}
	}
	private static ArrayList<String> lines(String path){
		final ByteFile in=RrnaResourceIO.open(path);final ArrayList<String> rows=new ArrayList<String>();
		try{for(byte[] row; (row=in.nextLine())!=null; ){rows.add(new String(row, StandardCharsets.UTF_8));}}
		finally{require(!in.close(), "Fixture output read failed");}return rows;
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
