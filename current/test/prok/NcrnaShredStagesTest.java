package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaShredStages;
import prok.NcrnaStageDiagLogSink;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Distinguishes seed occurrences, attempted models, accepted candidates and final calls.
 * @author Raiden
 */
public final class NcrnaShredStagesTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture root");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final ByteBuilder fasta=new ByteBuilder();for(int i=0; i<5; i++){fasta.append(">q").append(i).nl();for(int j=0; j<(i==2?100:200); j++){fasta.append('A');}fasta.nl();}
		write(root.resolve("input.fa"), fasta.toString());stages(root.resolve("stages.tsv"), false);
		write(root.resolve("calls.gff"), "##gff-version 3\nq0\tCallGenes\tRNA\t51\t150\t.\t+\t.\tRNA,model:cobalamin_consensus_0\n"+
			"q4\tCallGenes\tRNA\t51\t150\t.\t+\t.\tRNA,model:cobalaminX_consensus_0\n");
		run(root, "good", "stages.tsv", ".65");
		final ByteFile in=RrnaResourceIO.open(root.resolve("good.summary.tsv").toString());final LineParser1 p=new LineParser1('\t');int count=0, hbm=0;
		final long[] events={8, 4, 4, 5, 5, 3, 3, 3, 1, 1}, reached={4, 3, 3, 3, 3, 2, 2, 2, 1, 1};
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);check(p.terms()==11 && p.parseInt(4)==5 && p.parseInt(5)==900, "All five records, including unobserved input, must remain in both denominators");
			if(count<10){check(p.parseLong(6)==reached[count] && p.parseLong(7)==events[count], "Actual events and unique reached inputs differ at stage="+p.parseString(3));}
			else{check(p.termEquals("NA", 6) && p.termEquals("NOT_OBSERVED_IN_STAGE_DIAGNOSTICS", 10), "alignWindowOnce can score HBMs at equal thresholds; terminal scores do not supply evaluation counts");hbm++;}count++;
		}}finally{check(!in.close(), "Stage summary read failed");}
		check(count==12 && hbm==2, "Every requested stage must remain explicit");
		final ByteFile rows=RrnaResourceIO.open(root.resolve("good.loci.tsv").toString());int loci=0;
		try{rows.nextLine();for(byte[] row; (row=rows.nextLine())!=null; ){p.set(row);
			check(p.terms()==16 && p.termEquals("NA_NOT_OBSERVED", 14) && p.termEquals("NA_NOT_OBSERVED", 15), "Per-record HBM evaluation/pass counts require dedicated observations, even when a terminal best score exists");
			if(p.termEquals("q2", 0)){check(p.termEquals("NOT_OBSERVED", 3), "No SEED callback cannot be called an observed zero-hit sequence");}loci++;}}
		finally{check(!rows.close(), "Per-record stage read failed");}check(loci==5, "Zero-hit/unobserved records must not disappear");
		final ByteFile called=RrnaResourceIO.open(root.resolve("good.called.ids.txt").toString());
		try{check(java.util.Arrays.equals(called.nextLine(), new byte[]{'q','0'}) && called.nextLine()==null, "Final GFF, not verifier acceptance or a false model prefix, drives the CM review union");}
		finally{check(!called.close(), "Called-ID read failed");}
		stages(root.resolve("rescue.tsv"), true);reject(root, "rescue", "rescue.tsv", ".65");reject(root, "wrong_band", "stages.tsv", ".60");
		System.out.println("NCRNA_SHRED_STAGES_TEST_PASS");
	}
	private static void stages(Path path, boolean rescue){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());final NcrnaStageDiagLogSink log=new NcrnaStageDiagLogSink(out);
		try{
			out.println(NcrnaStageDiagLogSink.HEADER);
			for(int i:new int[]{0, 1, 3, 4}){log.seed("cobalamin", "q"+i, 0, 200, i==0?new int[]{20, 30}:i==1?new int[0]:new int[]{20});log.seed("cobalamin", "q"+i, 1, 200, new int[0]);}
			log.model("cobalamin", "q0", 0, 200, 1, 0, 199, 0, "cobalamin_0", .60f);
			log.model("cobalamin", "q0", 0, 200, 1, 0, 199, 1, "cobalamin_1", .65f);
			// Edge verification may compute HBM even though equal thresholds select IDPASS.
			log.window("cobalamin", "q0", 0, 200, 1, 0, 199, "ACCEPT_IDPASS", 2, 2, 1, "cobalamin_1", .65f, 0, .80f, 50, 149, false);
			log.model("cobalamin", "q0", 0, 200, 2, 0, 199, 0, "cobalamin_0", .70f);
			log.window("cobalamin", "q0", 0, 200, 2, 0, 199, "ACCEPT_IDPASS", 2, 1, 0, "cobalamin_0", .70f, 0, 0, 50, 149, false);
			log.path("cobalamin", "q0", 0, 200, "cobalamin_1", 50, 149, false);log.path("cobalamin", "q0", 0, 200, "cobalamin_0", 50, 149, false);log.path("cobalamin", "q0", 0, 200, "cobalamin_0", 50, 149, true);
			// Native appendSlow emits this finite value with an exponent.
			log.model("cobalamin", "q3", 0, 200, 1, 0, 199, 0, "cobalamin_0", .0001f);
			log.window("cobalamin", "q3", 0, 200, 1, 0, 199, "REJECT_ID", 1, 1, 0, "cobalamin_0", .0001f, 0, 0, -1, -1, false);
			log.model("cobalamin", "q4", 0, 200, 1, 0, 199, 0, "cobalamin_0", .66f);
			log.window("cobalamin", "q4", 0, 200, 1, 0, 199, rescue?"ACCEPT_RESCUE":"ACCEPT_IDPASS", 1, 1, 0, "cobalamin_0", .66f, 0, 0, 50, 149, false);
			log.path("cobalamin", "q4", 0, 200, "cobalamin_0", 50, 149, false);
		}finally{RrnaResourceIO.finish(out);}
	}
	private static void run(Path root, String out, String stages, String border){
		NcrnaShredStages.main(new String[]{"in="+root.resolve("input.fa"), "stages="+root.resolve(stages), "gff="+root.resolve("calls.gff"),
			"outprefix="+root.resolve(out), "expected=5", "family=cobalamin", "prefix=cobalamin", "label=PANEL", "idpass=.65", "idborderline="+border});
	}
	private static void reject(Path root, String out, String stages, String border){
		boolean failed=false;try{run(root, out, stages, border);}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve(out+".summary.tsv")), "Rescue evidence or a rescue band must not publish as raw-only accounting");
	}
	private static void write(Path path, String text){final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
