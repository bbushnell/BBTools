package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.CobalaminPrecisionReview;
import prok.CobalaminGenomicContext;
import prok.NcrnaReviewJoin;
import prok.NcrnaReviewUnion;
import prok.NcrnaStageDiagLogSink;
import prok.RrnaResourceIO;

/** No added TPs, explicit geometry/conflict decisions, complete sensitive evidence,
 * and seed-window event accounting are independent fixture expectations.
 * @author Brian Bushnell, Ganyu
 */
public final class CobalaminPrecisionReviewTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");
		final Path dir=Paths.get(args[0]);Files.createDirectory(dir);
		final String[] ids={"real", "geometry", "conflict", "partial", "none", "newfull"};
		final String hash="00000000000000000000";
		final StringBuilder loci=new StringBuilder(NcrnaReviewUnion.LOCUS_HEADER).append('\n');
		final StringBuilder context=new StringBuilder(NcrnaReviewJoin.CONTEXT_HEADER).append('\n');
		final StringBuilder features=new StringBuilder(NcrnaReviewJoin.OVERLAP_HEADER).append('\n');
		final StringBuilder members=new StringBuilder(NcrnaReviewJoin.MEMBER_HEADER).append('\n');
		final StringBuilder calls=new StringBuilder(NcrnaReviewUnion.PANEL_CALL_HEADER).append('\n');
		final StringBuilder seeds=new StringBuilder("call_index\temitted_family\toutcome\tproducing_windows\tminimum_distinct_keys\tmaximum_distinct_keys\n");
		final StringBuilder windows=new StringBuilder("call_index\tgenome\tseqid\temitted_family\toutcome\tstrand\tstage_line\tpass\twindow_start0_oriented\twindow_stop0_oriented\tseed_occurrences\tdistinct_seed_keys\n");
		final StringBuilder stages=new StringBuilder(NcrnaStageDiagLogSink.HEADER).append('\n');int stageLine=1;
		final StringBuilder neighborhood=new StringBuilder(CobalaminGenomicContext.HEADER).append('\n');
		for(String id:ids){
			loci.append(id).append("\tg\tA\tgenome.fa\t").append(id).append("\tcobalamin\t+\t100\t199\t1000\t50\t250\t400\tPREDICTION_UNRESOLVED\n");
			final boolean feature=id.equals("real") || id.equals("geometry") || id.equals("conflict");
			context.append(id).append("\tg\tA\t").append(id).append("\t100\t199\t+\t100\t189\t-89\t0.529101\t")
				.append(feature?1:0).append("\t1\t0\t").append(id.equals("geometry")?1:0).append("\t0\t")
				.append(feature && !id.equals("geometry")?1:0).append("\t0\t").append(feature?1:0).append("\t0\t0\t")
				.append(feature?"FEATURE_OVERLAP":"INTERGENIC").append("\tgenome.gff\t").append(hash).append("\tgenome.fa\t").append(hash).append('\n');
			features.append(feature(id, hash, 1, "REGION_METADATA", "region", 1, 1000, 100, false));
			if(feature){
				final boolean annotated=id.equals("geometry");final int start=id.equals("real")?193:id.equals("conflict")?140:95;
				features.append(feature(id, hash, 2, annotated?"TARGET_RNA_ANNOTATION":"CDS", annotated?"riboswitch":"CDS",
					start, 300, id.equals("real")?7:id.equals("conflict")?60:100, annotated));
			}
			members.append("point\t").append(id).append('\t').append(id).append("\tAPPARENT_FP\tAPPARENT_FP\n");
			calls.append(call(id, "APPARENT_FP"));
			stageLine=seedRows(seeds, windows, stages, stageLine, id, "APPARENT_FP", id.equals("none"));
			final int cds=id.equals("real")?7:id.equals("conflict")?60:0;
			neighborhood.append(id).append('\t').append(cds>0?"CDS_SAME":id.equals("geometry")?"OTHER_NCRNA":"INTERGENIC")
				.append('\t').append(cds).append('\t').append(cds/100.0).append("\t0\t0\tfalse")
				.append(id.equals("none")?"\t10\ttrue\tmetH\tCobalamin dependent enzyme\t2":"\tNA\tfalse\tNA\tNA\t0")
				.append(id.equals("none")?"\t10\ttrue\tmetH\tCobalamin dependent enzyme\t2\t+\n":"\tNA\tfalse\tNA\tNA\t0\tNA\n");
		}
		for(String id:new String[]{"tp1", "tp2", "unscored"}){
			final String outcome=id.equals("unscored")?"UNSCORED":"TP";calls.append(call(id, outcome));stageLine=seedRows(seeds, windows, stages, stageLine, id, outcome, false);
		}
		write(dir, "loci.tsv", loci.toString());write(dir, "context.tsv", context.toString());write(dir, "features.tsv", features.toString());
		write(dir, "members.tsv", members.toString());write(dir, "calls.tsv", calls.toString());write(dir, "seeds.tsv", seeds.toString());write(dir, "windows.tsv", windows.toString());
		write(dir, "stages.tsv", stages.toString());
		write(dir, "neighborhood.tsv", neighborhood.toString());
		final StringBuilder modelRows=new StringBuilder(CobalaminPrecisionReview.MODEL_HEADER).append('\n');
		for(String stage:new String[]{"DEFAULT", "MAX"}){for(String model:new String[]{"RF00174", "RF01482", "RF01689"}){
			final String name=stage+"_"+model;final boolean first=model.equals("RF00174");
			final StringBuilder summaries=new StringBuilder(NcrnaReviewJoin.CM_HEADER).append('\n');
			final StringBuilder hits=new StringBuilder(NcrnaReviewJoin.HIT_HEADER).append('\n');
			int line=0;
			for(String id:ids){
				if(stage.equals("MAX") && (id.equals("real") || id.equals("conflict"))){continue;}
				final boolean has=first && (stage.equals("DEFAULT")?(id.equals("real") || id.equals("conflict") || id.equals("partial")):(id.equals("geometry") || id.equals("newfull")));
				final boolean strict=has && !id.equals("partial") && !id.equals("newfull");
				summaries.append(id).append(has?"\t1\t1\t":"\t0\t0\t").append(strict?1:0)
					.append(strict?"\tONE_FULL_GA_OVERLAPPING_HIT\n":"\tNO_FULL_GA_OVERLAPPING_HIT\n");
				if(has){
					final int start=id.equals("real")?70:id.equals("geometry")?90:id.equals("newfull")?190:100;
					final int stop=id.equals("real")?170:id.equals("geometry")?210:id.equals("newfull")?379:199;
					final int overlap=id.equals("real")?71:id.equals("newfull")?10:100;
					hits.append(id).append('\t').append(++line).append('\t').append(id).append('\t').append(start).append('\t').append(stop).append("\t+\t1\t")
						.append(id.equals("partial")?100:189).append("\tno\t").append(id.equals("partial")?15:60)
						.append("\t1e-10\ttrue\t").append(overlap).append('\t').append(strict).append('\n');
				}
			}
			write(dir, name+".loci.tsv", summaries.toString());write(dir, name+".hits.tsv", hits.toString());
			modelRows.append(model).append('\t').append(stage).append('\t').append(dir.resolve(name+".loci.tsv")).append('\t').append(dir.resolve(name+".hits.tsv")).append('\n');
		}}
		write(dir, "models.tsv", modelRows.toString());
		final String decisions=CobalaminPrecisionReview.DECISION_HEADER+"\nreal\tREAL-B12\tCore is upstream of overlapping CDS.\n"+
			"geometry\tANNOTATED-GEOMETRY\tShort prediction overlaps an annotated B12 locus.\nconflict\tFEATURE-CONFLICT\tCore overlaps a source feature.\n";
		write(dir, "decisions.tsv", decisions);
		final String runPrefix=CobalaminPrecisionReview.RUN_HEADER+"\npoint\t"+dir.resolve("calls.tsv")+'\t'+dir.resolve("seeds.tsv")+'\t'+dir.resolve("windows.tsv")+"\t9\t";
		write(dir, "runs.tsv", runPrefix+dir.resolve("stages.tsv")+'\n');
		run(dir, "pending", "decisions.tsv");
		LineParser1 p=second(dir.resolve("pending.precision.tsv"));
		check(p.parseInt(1)==2 && p.parseInt(2)==6 && p.parseInt(9)==1 && p.termEquals("NA", 12),
			"Unknown full-GA geometry blocks final optimistic precision without altering original TP/FP counts");
		check(contains(dir.resolve("pending.loci.tsv"), "newfull\tFULL-GA-REVIEW-REQUIRED"), "A non-reciprocal full-GA hit must be routed to inspection");
		write(dir, "accepted.tsv", decisions+"newfull\tREAL-B12\tExplicitly inspected non-reciprocal full-GA geometry.\n");
		run(dir, "accepted", "accepted.tsv");p=second(dir.resolve("accepted.precision.tsv"));
		check(p.parseInt(1)==2 && p.parseInt(2)==6 && p.parseInt(3)==1 && p.parseInt(4)==2 && p.parseInt(5)==1
			&& p.parseInt(6)==1 && p.parseInt(7)==1 && p.parseInt(8)==1 && p.parseInt(10)==3
			&& Math.abs(p.parseDouble(11)-.25)<1e-9 && Math.abs(p.parseDouble(12)-.4)<1e-9,
			"Only REAL and GEOMETRY leave the optimistic FP denominator; feature conflict stays FP and no extra TP is invented");
		check(contains(dir.resolve("accepted.no_support.tsv"), "point\tRF01689\t1\t1\t1\t0"),
			"One NO-SUPPORT call with two producing windows remains one call");
		check(contains(dir.resolve("accepted.no_support_histogram.tsv"), "point\tRF01689\t1\t1\tCOMPLETE")
			&& contains(dir.resolve("accepted.no_support_histogram.tsv"), "point\tRF01689\t2\t1\tCOMPLETE"),
			"Distinct-key histogram counts both windows separately");
		check(contains(dir.resolve("accepted.hits.tsv"), "RF00174\tDEFAULT\treal\t1\treal\t70\t170\t+\t1\t189\tno\t60\t1e-10\ttrue\t71\ttrue\ttrue\t0\t0"),
			"Seven-base caller/CDS overlap is not a CM-core conflict when the core ends23 bases upstream");
		check(contains(dir.resolve("accepted.loci.tsv"), "partial\tPARTIAL/BELOW-GA\tReported overlapping partial-span or below-GA evidence from either scan\t1\t0\t0\tfalse\ttrue\ttrue"),
			"DEFAULT weak overlap remains PARTIAL even when MAX has none; both flags stay visible");
		check(contains(dir.resolve("accepted.calls.tsv"), "point\tnone\tRF01689\tAPPARENT_FP\tnone\tNO-SUPPORT\tfalse\t2\t1\t2\tcobalamin_RF01689_fixture\t0.67\t0.71\tMULTIPLE_MATCHING_VALUES"),
			"Different matching-window identities must remain an explicit range, never a guessed selected value");
		check(contains(dir.resolve("accepted.calls.tsv"), "point\ttp2\tRF00174\tTP\tNA\tTP\tfalse\t1\t1\t1\tcobalamin_RF00174_fixture\t0.67\t0.67\tONE_MATCHING_VALUE"),
			"Minus-strand identity joins must convert oriented coordinates back to the exact emitted span");
		check(contains(dir.resolve("accepted.context_summary.tsv"), "point\tRF01689\tNO-SUPPORT\tINTERGENIC\ttrue\ttrue\t1"),
			"A nearby B12-name hint groups the unsupported call but cannot admit it or change precision");
		write(dir, "wrong_stages.tsv", stages.toString().replace("cobalamin_RF01689_fixture", "cobalamin_RF01689_wrong"));
		write(dir, "wrong_runs.tsv", runPrefix+dir.resolve("wrong_stages.tsv")+'\n');boolean wrongModel=false;
		try{run(dir, "wrong_identity", "accepted.tsv", "models.tsv", "wrong_runs.tsv");}catch(IllegalArgumentException e){wrongModel=true;}
		check(wrongModel && !Files.exists(dir.resolve("wrong_identity.precision.tsv")), "Matching coordinates with a different model cannot establish identity");
		// Remove a whole model/stage row, never interpret its absence as no support.
		final String missing=modelRows.substring(0, modelRows.lastIndexOf("RF01689\tMAX"));
		write(dir, "missing_models.tsv", missing);boolean failed=false;
		try{run(dir, "missing", "accepted.tsv", "missing_models.tsv");}catch(IllegalArgumentException e){failed=true;}
		check(failed && !Files.exists(dir.resolve("missing.precision.tsv")), "Missing sensitive evidence must fail before publishing precision");
		System.out.println("COBALAMIN_PRECISION_REVIEW_TEST_PASS");
	}
	private static String family(String id){return id.equals("none")?"RF01689":id.equals("partial")?"RF01482":"RF00174";}
	private static String call(String id, String outcome){
		return "FULL_CM\tunion\t"+id+"\tg\t"+id+'\t'+family(id)+"\t100\t199\t"+(id.equals("tp2")?"-":"+")+"\tcobalamin_"+family(id)+"_fixture\t"+outcome+"\tNA\tNA\tNA\n";
	}
	private static int seedRows(StringBuilder calls, StringBuilder windows, StringBuilder stages, int line, String id, String outcome, boolean two){
		final boolean minus=id.equals("tp2");final int a=minus?750:50, z=minus?999:448;
		calls.append(id).append('\t').append(family(id)).append('\t').append(outcome).append(two?"\t2\t1\t2\n":"\t1\t1\t1\n");
		for(int n=1; n<=(two?2:1); n++){windows.append(id).append("\tg\t").append(id).append('\t').append(family(id)).append('\t').append(outcome)
			.append(minus?"\t-\t":"\t+\t").append(++line).append("\t1\t").append(a).append('\t').append(z).append("\t2\t").append(n).append('\n');
			stages.append("WINDOW\tcobalamin_").append(family(id).toLowerCase(java.util.Locale.ROOT)).append('\t').append(id).append(" description\t")
				.append(minus?1:0).append("\t1000\t1\t").append(a).append('\t').append(z).append("\tACCEPT_IDPASS\t2\t1\t0\tcobalamin_")
				.append(family(id)).append("_fixture\t").append(n==1?"0.67":"0.71").append("\t0.0\t0.0\t")
				.append(minus?801:99).append('\t').append(minus?900:198).append("\tfalse\t.\t.\n");
		}return line;
	}
	private static String feature(String id, String hash, int line, String kind, String type, int start, int stop, int overlap, boolean annotated){
		return id+"\tg\tA\tgenome.fa\t"+hash+"\tgenome.gff\t"+hash+'\t'+line+'\t'+kind+'\t'+annotated+'\t'+overlap+
			"\tSAME\tfalse\t"+id+"\tsource\t"+type+'\t'+start+'\t'+stop+"\t.\t+\t.\t"+(annotated?"Dbxref=Rfam:RF00174;":"ID=feature;")+'\n';
	}
	private static void run(Path dir, String out, String decisions){run(dir, out, decisions, "models.tsv");}
	private static void run(Path dir, String out, String decisions, String models){
		run(dir, out, decisions, models, "runs.tsv");
	}
	private static void run(Path dir, String out, String decisions, String models, String runs){
		CobalaminPrecisionReview.main(new String[]{"loci="+dir.resolve("loci.tsv"), "context="+dir.resolve("context.tsv"),
			"features="+dir.resolve("features.tsv"), "members="+dir.resolve("members.tsv"), "models="+dir.resolve(models),
			"decisions="+dir.resolve(decisions), "runs="+dir.resolve(runs), "outprefix="+dir.resolve(out),
			"neighborhood="+dir.resolve("neighborhood.tsv"),
			"expected=6", "events=6", "sensitive=4", "runcount=1"});
	}
	private static LineParser1 second(Path path){
		final ByteFile in=RrnaResourceIO.open(path.toString());final LineParser1 p=new LineParser1('\t');
		try{check(in.nextLine()!=null, "Expected output header");final byte[] row=in.nextLine();check(row!=null, "Expected output data");p.set(row.clone());}
		finally{check(!in.close(), "Fixture read failed");}return p;
	}
	private static boolean contains(Path path, String prefix){
		final ByteFile in=RrnaResourceIO.open(path.toString());boolean found=false;
		try{for(byte[] row; (row=in.nextLine())!=null; ){if(new String(row, java.nio.charset.StandardCharsets.UTF_8).startsWith(prefix)){found=true;}}}
		finally{check(!in.close(), "Fixture read failed");}return found;
	}
	private static void write(Path dir, String name, String text){final ByteStreamWriter w=RrnaResourceIO.writer(dir.resolve(name).toString());try{w.print(text);}finally{RrnaResourceIO.finish(w);}}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
