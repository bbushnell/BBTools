package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.CobalaminProductionPrecision;
import prok.CobalaminGenomicContext;
import prok.CobalaminPrecisionReview;
import prok.NcrnaReviewJoin;
import prok.NcrnaReviewUnion;
import prok.RrnaResourceIO;

/** Shifting sorted IDs cannot transfer dispositions; reverse-context CM hits
 * must map back to the correct genome and genomic interval.
 * @author Brian Bushnell, Ganyu
 */
public final class CobalaminProductionPrecisionTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply one fresh fixture directory");
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		write(root, "old.tsv", NcrnaReviewUnion.LOCUS_HEADER+'\n'+locus("review_0", 10, 20, "+", 5)+locus("review_1", 100, 110, "+", 5));
		write(root, "new.tsv", NcrnaReviewUnion.LOCUS_HEADER+'\n'+locus("review_0", 10, 20, "+", 5)+locus("review_1", 200, 220, "-", 5)+locus("review_2", 100, 110, "+", 5));
		write(root, "review.tsv", CobalaminProductionPrecision.REVIEW_HEADER+"\nreview_0\tPARTIAL/BELOW-GA\toriginal partial\t1\t1\t0\tfalse\ttrue\tfalse\nreview_1\tREAL-B12\toriginal inspected call\t1\t0\t1\tfalse\tfalse\tfalse\n");
		write(root, "source.tsv", NcrnaReviewUnion.LOCUS_HEADER+'\n'+locus("source_minus", 100, 300, "-", 10));
		write(root, "hits.tbl", "source_minus - Cobalamin RF00174 cm 1 189 91 111 + no 1 0.5 0 40 1e-10 ! fixture\n# Program: cmsearch\n# Version: 1.1.5 (fixture)\n# [ok]\n");
		write(root, "sources.tsv", CobalaminProductionPrecision.SOURCE_HEADER+"\nreview_1\tsource_minus\tREVIEW\t"+root.resolve("source.tsv")+'\t'+root.resolve("hits.tbl")+"\tRF00174\tDEFAULT\tRAW\n");
		write(root, "decisions.tsv", CobalaminPrecisionReview.DECISION_HEADER+"\nreview_1\tANNOTATED-GEOMETRY\tInspected same-strand annotated core\n");
		write(root, "context.tsv", NcrnaReviewJoin.CONTEXT_HEADER+'\n'+context("review_0", 10, 20, "+", 0)+context("review_1", 200, 220, "-", 1)+context("review_2", 100, 110, "+", 0));
		write(root, "features.tsv", NcrnaReviewJoin.OVERLAP_HEADER+'\n'+feature());
		write(root, "neighborhood.tsv", CobalaminGenomicContext.HEADER+'\n'+neighborhood("review_0")+neighborhood("review_1")+neighborhood("review_2"));
		final String union=call("union", "a", 10, 20, "+", "APPARENT_FP")+call("union", "b", 200, 220, "-", "APPARENT_FP")+call("union", "c", 100, 110, "+", "TP");
		write(root, "calls.tsv", NcrnaReviewUnion.PANEL_CALL_HEADER+'\n'+union+union.replace("\tunion\t", "\tRF00174\t"));
		write(root, "summary.tsv", CobalaminProductionPrecision.SUMMARY_HEADER+'\n'+summary("RF00174", 3, 1, 2)+summary("RF01482", 0, 0, 0)+summary("RF01689", 0, 0, 0)+summary("union", 3, 1, 2));
		run(root, "valid");
		check(contains(root.resolve("valid.mapping.tsv"), "review_1\treview_2\tg\tseq\t+\t100\t110"), "Original ID must move by exact coordinates when a new locus shifts sorting");
		check(contains(root.resolve("valid.loci.tsv"), "review_1\tNA\tANNOTATED-GEOMETRY"), "Added locus must use independent evidence, never the old disposition at its new row number");
		check(contains(root.resolve("valid.precision.tsv"), "union\t1\t2\t0\t0\t1\t1\t0.33333333\t0.50000000"), "Only the supported FP is excluded optimistically; TP stays unchanged");
		check(contains(root.resolve("valid.hits.tsv"), "review_1\tsource_minus\t"+root.resolve("hits.tbl")+"\t1\tRF00174\tDEFAULT\tseq\t200\t220\t-"), "Reverse-context positions91..111 must project to genomic200..220, not forward coordinates");
		write(root, "projected.tsv", NcrnaReviewJoin.HIT_HEADER+"\nsource_minus\t1\tseq\t200\t220\t-\t1\t189\tno\t40\t1e-10\ttrue\t21\tfalse\n");
		write(root, "sources.tsv", CobalaminProductionPrecision.SOURCE_HEADER+"\nreview_1\tsource_minus\tREVIEW\t"+root.resolve("source.tsv")+'\t'+root.resolve("projected.tsv")+"\tRF00174\tDEFAULT\tPROJECTED\n");
		run(root, "projected");
		check(contains(root.resolve("projected.precision.tsv"), "union\t1\t2\t0\t0\t1\t1"), "Existing genomic hits must be evaluated at the new exact call span, not reuse the old reciprocal-overlap decision");
		write(root, "pilot.tsv", prok.NcrnaCmLocusAudit.HEADER+"\nsource_minus\tg\tassembly\t1\tbacteria\tfixture\tseq\t-\t100-300\t201\t1\tID=fixture\tfalse\tfalse\tANNOTATION_FULL_COPY\t1000\tfalse\t0\t10\t10\t1\n");
		write(root, "sources.tsv", CobalaminProductionPrecision.SOURCE_HEADER+"\nreview_1\tsource_minus\tPILOT\t"+root.resolve("pilot.tsv")+'\t'+root.resolve("hits.tbl")+"\tRF00174\tHISTORICAL_GA\tRAW\n");
		run(root, "historical");
		check(contains(root.resolve("historical.sources.tsv"), "review_1\tsource_minus\tPILOT"), "Historical full-GA support stays labeled separately from DEFAULT/MAX scans that were never performed");
		write(root, "sources.tsv", CobalaminProductionPrecision.SOURCE_HEADER+"\nreview_1\tsource_minus\tREVIEW\t"+root.resolve("source.tsv")+'\t'+root.resolve("hits.tbl")+"\tRF00174\tDEFAULT\tRAW\n");
		write(root, "source.tsv", NcrnaReviewUnion.LOCUS_HEADER+'\n'+locus("source_minus", 400, 500, "-", 10));
		fails(root, "uncovered", "A distant source window cannot lend its evidence to the call");
		write(root, "source.tsv", NcrnaReviewUnion.LOCUS_HEADER+'\n'+locus("source_minus", 100, 300, "-", 10).replace("\tg\t", "\tother_genome\t"));
		fails(root, "wrong_genome", "Shared contig text is not permission to reuse evidence from another genome");
		write(root, "source.tsv", NcrnaReviewUnion.LOCUS_HEADER+'\n'+locus("source_minus", 100, 300, "-", 10));
		write(root, "features.tsv", NcrnaReviewJoin.OVERLAP_HEADER+'\n'+feature().replace("Dbxref=Rfam:RF00174", "Dbxref=Rfam:RF001740"));
		fails(root, "wrong_annotation", "Only the exact Rfam token can support annotated geometry");
		write(root, "features.tsv", NcrnaReviewJoin.OVERLAP_HEADER+'\n'+feature());
		write(root, "decisions.tsv", CobalaminPrecisionReview.DECISION_HEADER+'\n');
		final String[] models={"RF00174", "RF01482", "RF01689"};final int[] lengths={189, 149, 81};
		String six=CobalaminProductionPrecision.SOURCE_HEADER+'\n';
		for(int i=0; i<models.length; i++){
			write(root, models[i]+".tbl", cm(models[i], lengths[i], 10));
			for(String stage:new String[]{"DEFAULT", "MAX"}){six+="review_1\tsource_minus\tREVIEW\t"+root.resolve("source.tsv")+'\t'+root.resolve(models[i]+".tbl")+'\t'+models[i]+'\t'+stage+"\tRAW\n";}
		}
		write(root, "sources.tsv", six);run(root, "partial");
		check(contains(root.resolve("partial.loci.tsv"), "review_1\tNA\tPARTIAL/BELOW-GA")
			&& contains(root.resolve("partial.precision.tsv"), "union\t1\t2\t0\t0\t0\t2\t0.33333333\t0.33333333"),
			"Six below-GA scans close evidence coverage without an optimistic exclusion or a new truth call");
		write(root, "sources.tsv", six.substring(0, six.lastIndexOf("review_1\t")));
		fails(root, "missing_scan", "Five scans cannot silently masquerade as all-model DEFAULT/MAX coverage");
		write(root, "sources.tsv", six);write(root, "RF00174.tbl", cm("RF00174", 189, 40));
		fails(root, "undecided_full", "New full-GA support requires an explicit reviewed disposition, even with all six scans");
		System.out.println("COBALAMIN_PRODUCTION_PRECISION_TEST_PASS");
	}
	private static String cm(String model, int length, int bits){return "source_minus - fixture "+model+" cm 1 "+length+" 91 111 + no 1 0.5 0 "+bits+" 1e-10 ! fixture\n# Program: cmsearch\n# Version: 1.1.5 (fixture)\n# [ok]\n";}
	private static String locus(String id, int start, int stop, String strand, int flank){return id+"\tg\tassembly\tgenome.fa\tseq\tcobalamin\t"+strand+'\t'+start+'\t'+stop+"\t1000\t"+flank+'\t'+flank+'\t'+(stop-start+1+2*flank)+"\tPREDICTION_UNRESOLVED\n";}
	private static String context(String id, int start, int stop, String strand, int features){
		final String[] f=new String[26];java.util.Arrays.fill(f, "0");
		f[0]=id;f[1]="g";f[2]="assembly";f[3]="seq";f[4]=""+start;f[5]=""+stop;f[6]=strand;f[7]=""+(stop-start+1);f[8]="189";f[11]=""+features;f[22]="source.gff";f[23]="00000000000000000000";f[24]="genome.fa";f[25]="11111111111111111111";
		return String.join("\t", f)+'\n';
	}
	private static String feature(){return "review_1\tg\tassembly\tgenome.fa\t11111111111111111111\tsource.gff\t00000000000000000000\t1\tTARGET_RNA_ANNOTATION\ttrue\t21\tSAME\tfalse\tseq\tNCBI\tncRNA\t200\t220\t.\t-\t.\tID=test;Dbxref=Rfam:RF00174\n";}
	private static String neighborhood(String id){return id+"\tNO_CDS_OVERLAP\t0\t0\t0\t0\tfalse\tNA\tfalse\tNA\tNA\t0\tNA\tfalse\tNA\tNA\t0\tNA\n";}
	private static String call(String scope, String id, int start, int stop, String strand, String outcome){return "FULL_CM\t"+scope+'\t'+id+"\tg\tseq\tRF00174\t"+start+'\t'+stop+'\t'+strand+"\tcobalamin_RF00174_fixture\t"+outcome+"\tNA\tNA\tNA\n";}
	private static String summary(String scope, int count, int tp, int fp){return "FULL_CM\t"+scope+"\t1\t"+tp+'\t'+count+'\t'+tp+'\t'+fp+"\t0\t0\t0\tNA\tNA\tNA\t0\t0\t0\n";}
	private static void run(Path root, String output){
		CobalaminProductionPrecision.main(new String[]{"old="+root.resolve("old.tsv"), "new="+root.resolve("new.tsv"), "review="+root.resolve("review.tsv"), "sources="+root.resolve("sources.tsv"), "decisions="+root.resolve("decisions.tsv"), "calls="+root.resolve("calls.tsv"), "summary="+root.resolve("summary.tsv"), "context="+root.resolve("context.tsv"), "features="+root.resolve("features.tsv"), "neighborhood="+root.resolve("neighborhood.tsv"), "outprefix="+root.resolve(output), "oldcount=2", "newcount=3", "fpcount=2"});
	}
	private static void fails(Path root, String output, String why){boolean failed=false;try{run(root, output);}catch(IllegalArgumentException expected){failed=true;}check(failed && !Files.exists(root.resolve(output+".precision.tsv")), why);}
	private static boolean contains(Path path, String prefix){final ByteFile in=RrnaResourceIO.open(path.toString());boolean found=false;try{for(byte[] row;(row=in.nextLine())!=null;){if(new String(row, java.nio.charset.StandardCharsets.UTF_8).startsWith(prefix)){found=true;}}}finally{check(!in.close(), "Fixture read failed");}return found;}
	/** Rewrites are deliberate fixture mutations inside the newly created test directory. */
	private static void write(Path root, String name, String value){
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(root.resolve(name).toString(), FileFormat.TEXT, null, true, true, false, false));
		out.start();try{out.print(value);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
