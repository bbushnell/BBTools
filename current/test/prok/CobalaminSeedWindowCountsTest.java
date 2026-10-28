package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import prok.CobalaminPanelPrep;
import prok.CobalaminSeedWindowCounts;
import prok.NcrnaReviewUnion;
import prok.NcrnaStageDiagLogSink;
import prok.RrnaResourceIO;

/** Exercises saved-call/stage joins, both strands and event-versus-locus units.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminSeedWindowCountsTest {

	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String a="AAAAAAAAAAAAAAAAA", g="GGGGGGGGGGGGGGGGG", t="TTTTTTTTTTTTTTTTTT";
		write(root.resolve("seeds.fa"), ">a\n"+a+"\n>g\n"+g+"\n");
		write(root.resolve("genome.fa"), ">one description\n"+a+"A\n>two\n"+a+"N"+g+"\n>minus\n"+t+"\n");
		write(root.resolve("inputs.txt"), root.resolve("genome.fa").toString()+"\n");
		write(root.resolve("contigs.tsv"), CobalaminPanelPrep.CONTIG_HEADER+"\ngenome\tGCF_1.1\tone\t18\ngenome\tGCF_1.1\ttwo\t35\ngenome\tGCF_1.1\tminus\t18\n");
		final String calls=NcrnaReviewUnion.PANEL_CALL_HEADER+'\n'+
			call("1", "one", 18, "+", "TP")+call("2", "two", 35, "+", "APPARENT_FP")+call("3", "minus", 18, "-", "APPARENT_FP");
		write(root.resolve("calls.tsv"), calls);
		final String stages=NcrnaStageDiagLogSink.HEADER+'\n'+window("one description", 18, 0, 1, 2)+window("one description", 18, 0, 2, 2)
			+window("two", 35, 0, 1, 2)+window("minus", 18, 1, 1, 2);
		write(root.resolve("stages.tsv"), stages);run(root, "good", "stages.tsv");
		final ByteFile in=RrnaResourceIO.open(root.resolve("good.calls.tsv").toString());final LineParser1 p=new LineParser1('\t');int seen=0;
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);seen++;final int id=p.parseInt(0), distinct=id==2?2:1, events=id==1?2:1;
			check(id>=1 && id<=3 && p.parseInt(3)==events && p.parseInt(4)==distinct && p.parseInt(5)==distinct,
				"Repeats count once per window, different keys count twice and multiple producing windows remain events");
		}}finally{check(!in.close(), "Fixture call table read failed");}
		check(seen==3, "Three emitted calls must survive the complete join");
		final ByteFile hist=RrnaResourceIO.open(root.resolve("good.histogram.tsv").toString());long total=0;
		try{hist.nextLine();for(byte[] row; (row=hist.nextLine())!=null; ){
			p.set(row);if(!p.termEquals("union", 0)){continue;}final long n=p.parseLong(3);total+=n;
			if(p.termEquals("TP", 1) && p.termEquals("1", 2)){check(n==2, "Two TP-producing windows are events for one emitted TP");}
			if(p.termEquals("APPARENT_FP", 1) && (p.termEquals("1", 2) || p.termEquals("2", 2))){check(n==1, "Minus-strand repeats and two-key forward windows must occupy separate FP buckets");}
		}}finally{check(!hist.close(), "Fixture histogram read failed");}
		check(total==4, "Histogram denominator is four windows, not three calls");
		write(root.resolve("wrong.tsv"), stages.replace("ACCEPT_IDPASS\t2", "ACCEPT_IDPASS\t3"));
		boolean failed=false;try{run(root, "wrong", "wrong.tsv");}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve("wrong.calls.tsv")), "Incorrect saved occurrence counts must fail before publishing a histogram");
		write(root.resolve("missing.tsv"), NcrnaStageDiagLogSink.HEADER+'\n'+window("one description", 18, 0, 1, 2));
		failed=false;try{run(root, "missing", "missing.tsv");}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Missing exact accepted windows must fail instead of inventing a proximity match");
		identityComparison(root.resolve("identity"));
		System.out.println("COBALAMIN_SEED_WINDOW_COUNTS_TEST_PASS");
	}
	private static void identityComparison(Path root)throws Exception{
		Files.createDirectory(root);
		final String model="cobalamin_RF00174_fixture n=3", a="AAAAAAAAAAAAAAAAAA", t="TTTTTTTTTTTTTTTTTT";
		final char[] padding=new char[132];java.util.Arrays.fill(padding, 'N');final String ns=new String(padding);
		write(root.resolve("models.fa"), ">"+model+"\n"+a.substring(1)+"\n");
		write(root.resolve("seeds.fa"), ">seed\n"+a.substring(1)+"\n");
		write(root.resolve("genome.fa"), ">short\n"+a+"\n>wide\n"+a+ns+"\n>reverse\n"+ns+t+"\n");
		write(root.resolve("inputs.txt"), root.resolve("genome.fa")+"\n");
		final StringBuilder contigs=new StringBuilder(CobalaminPanelPrep.CONTIG_HEADER).append('\n');
		final StringBuilder calls=new StringBuilder(NcrnaReviewUnion.PANEL_CALL_HEADER).append('\n');
		final StringBuilder panel=new StringBuilder(CobalaminPanelPrep.HEADER).append('\n');
		final StringBuilder stages=new StringBuilder(NcrnaStageDiagLogSink.HEADER).append('\n');
		final String[] names={"short", "wide", "reverse"};
		for(int i=0; i<3; i++){
			final int length=i==0 ? 18 : 150, strand=i==2 ? 1 : 0, start=i==2 ? 133 : 1, stop=i==2 ? 150 : 18;
			final String seq=names[i], sign=strand==0 ? "+" : "-";
			contigs.append("genome\tGCF_1.1\t").append(seq).append('\t').append(length).append('\n');
			calls.append("FULL_CM\tunion\t").append(i).append("\tgenome\t").append(seq).append("\tRF00174\t")
				.append(start).append('\t').append(stop).append('\t').append(sign).append('\t').append(model).append("\tTP\tt").append(i).append("\t0\t0\n");
			panel.append('t').append(i).append("\tgenome\t").append(seq).append("\tRF00174\tRF00174\t")
				.append(start).append('\t').append(stop).append('\t').append(start).append('\t').append(stop).append('\t').append(sign).append('\n');
			final String prefix="cobalamin_rf00174\t"+seq+'\t'+strand+'\t'+length+"\t1\t0\t"+(length-1);
			stages.append("MODEL\t").append(prefix).append("\tALIGNMENT\t.\t.\t0\t").append(model).append("\t1.0\t.\t.\t.\t.\t.\t.\t.\n");
			stages.append("WINDOW\t").append(prefix).append("\tACCEPT_IDPASS\t2\t1\t0\t").append(model).append("\t1.0\t0\t0\t0\t17\tfalse\t.\t.\n");
		}
		write(root.resolve("contigs.tsv"), contigs.toString());write(root.resolve("calls.tsv"), calls.toString());
		write(root.resolve("panel.tsv"), panel.toString());write(root.resolve("stages.tsv"), stages.toString());
		paired(root, "paired");
		final ByteFile in=RrnaResourceIO.open(root.resolve("paired.identities.tsv").toString());final LineParser1 p=new LineParser1('\t');int count=0;
		try{in.nextLine();for(byte[] row;(row=in.nextLine())!=null;){p.set(row);count++;
			check(p.terms()==12 && p.parseInt(6)==18 && p.parseFloat(7)==1 && p.parseFloat(8)==1 && p.parseFloat(9)==1
				&& p.parseFloat(10)==0 && p.parseFloat(11)==0, "Shorter-query and wide Quantum pairs must preserve exact truth geometry on both strands");
		}}finally{check(!in.close(), "Paired fixture output read failed");}
		check(count==3, "Every exact emitted-model/window pair requires a diagnostic row");
		write(root.resolve("stages.tsv"), stages.toString().replace("\t1.0\t", "\t0.7\t"));
		boolean failed=false;try{paired(root, "bad");}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve("bad.identities.tsv")), "A saved identity that cannot replay must never become a measured cutoff offset");
	}
	private static void paired(Path root, String out){
		CobalaminSeedWindowCounts.main(new String[]{"calls="+root.resolve("calls.tsv"), "stages="+root.resolve("stages.tsv"), "contigs="+root.resolve("contigs.tsv"),
			"inputs="+root.resolve("inputs.txt"), "seeds="+root.resolve("seeds.fa"), "outprefix="+root.resolve(out), "expected=3", "refs="+root.resolve("models.fa"), "panel="+root.resolve("panel.tsv"), "quantum=120"});
	}
	private static String call(String id, String contig, int length, String strand, String outcome){
		return "FULL_CM\tunion\t"+id+"\tgenome\t"+contig+"\tRF00174\t1\t"+length+'\t'+strand+"\tcobalamin_RF00174_fixture\t"+outcome+"\tNA\tNA\tNA\n";
	}
	private static String window(String contig, int length, int strand, int pass, int raw){
		return "WINDOW\tcobalamin_rf00174\t"+contig+'\t'+strand+'\t'+length+'\t'+pass+"\t0\t"+(length-1)+"\tACCEPT_IDPASS\t"+raw
			+"\t1\t0\tcobalamin_RF00174_fixture\t0.8\t0.0\t0.0\t0\t"+(length-1)+"\tfalse\t.\t.\n";
	}
	private static void run(Path root, String out, String stages){
		CobalaminSeedWindowCounts.main(new String[]{"calls="+root.resolve("calls.tsv"), "stages="+root.resolve(stages), "contigs="+root.resolve("contigs.tsv"),
			"inputs="+root.resolve("inputs.txt"), "seeds="+root.resolve("seeds.fa"), "outprefix="+root.resolve(out), "expected=3"});
	}
	private static void write(Path path, String text){
		// Only fixture inputs in this freshly created directory are replaced between cases.
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path.toString(), FileFormat.TEXT, null, true, true, false, false));
		out.start();try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
