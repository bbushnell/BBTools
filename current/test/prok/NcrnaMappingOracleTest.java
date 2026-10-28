package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaCalibrationPrep;
import prok.NcrnaMappingOracle;
import prok.RrnaResourceIO;
import prok.TrnaKmerIndex;
import structures.ByteBuilder;

/** Checks native fallback, strand-aware locus unions and deterministic threaded output.
 * @author Brian, Raiden
 */
public final class NcrnaMappingOracleTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String model="ACGTACGATCGGTACCTAGCTACGATCGTACGTTACCGATCGTACGTACGA";
		final byte[] bases=model.getBytes(StandardCharsets.US_ASCII), query=new byte[200];Arrays.fill(query, (byte)'A');
		System.arraycopy(bases, 0, query, 75, bases.length);
		write(root.resolve("ref.fa"), ">cobalamin_0\n"+model+"\n>cobalamin_1\n"+model+'\n');
		write(root.resolve("truth.tsv"), NcrnaCalibrationPrep.TRUTH_HEADER+"\nctg\tcobalamin\t176\t"+(175+bases.length)+"\t+\nctg\tcobalamin\t"+(826-bases.length)+"\t825\t-\n");
		final ByteBuilder fasta=new ByteBuilder(), windows=new ByteBuilder();windows.append(NcrnaCalibrationPrep.WINDOW_HEADER).nl();
		for(int i=0; i<80; i++){
			final boolean edge=i>=78;final int strand=i<39?0:1;
			fasta.append(">w").append(i).nl().append(query).nl();
			windows.append('w').append(i).append("\tg\tctg\t").append(strand).append("\t1000\t1\t")
				.append(edge?0:100).tab().append(edge?199:299).append("\tACCEPT_IDPASS\t2\t0\t1\t").append(edge).nl();
		}
		write(root.resolve("windows.fa"), fasta.toString());write(root.resolve("windows.tsv"), windows.toString());
		for(int threads:new int[]{1,16}){
			NcrnaMappingOracle.main(new String[]{"in="+root.resolve("windows.fa"), "windows="+root.resolve("windows.tsv"), "truth="+root.resolve("truth.tsv"),
				"ref="+root.resolve("ref.fa"), "outprefix="+root.resolve("t"+threads), "expected=80", "minlen=20", "quantum=120", "idpass=.65", "t="+threads});
		}
		for(String suffix:new String[]{".curves.tsv", ".windows.tsv", ".models.tsv.gz", ".losses.tsv", ".loss_summary.tsv"}){
			final ByteFile a=RrnaResourceIO.open(root.resolve("t1"+suffix).toString()), b=RrnaResourceIO.open(root.resolve("t16"+suffix).toString());
			try{
				for(byte[] line; (line=a.nextLine())!=null; ){check(Arrays.equals(line, b.nextLine()), "Oracle thread count changed ordered output: "+suffix);}
				check(b.nextLine()==null, "Parallel oracle emitted extra rows");
			}finally{final boolean error=a.close()|b.close();check(!error, "Parity input read failed");}
		}
		final ByteFile curves=RrnaResourceIO.open(root.resolve("t16.curves.tsv").toString());final LineParser1 p=new LineParser1('\t');int rows=0;
		try{
			curves.nextLine();
			for(byte[] row; (row=curves.nextLine())!=null; ){
				p.set(row);check(p.terms()==14 && p.parseInt(3)==80 && p.parseInt(4)==78 && p.parseInt(5)==78 && p.parseInt(6)==78,
					"Physical-edge events must stay separate; all78 interior exact-model events pass");
				check(p.parseInt(7)==2 && p.parseInt(8)==2 && p.parseInt(9)==2 && p.parseInt(12)==156,
					"Repeated windows must union to two strand-specific loci, never78 alleged genes");
				if(p.parseInt(1)==1 || p.parseInt(2)==80){check(p.parseInt(10)==78, "Top1 and keep-best fallback each schedule one model per event");}
				rows++;
			}
		}finally{check(!curves.close(), "Curve fixture read failed");}
		check(rows==196, "Every declared k/topN/minHits point, including larger and all-model caps, needs one result row");
		final TrnaKmerIndex nativeIndex=new TrnaKmerIndex(new byte[][]{bases,bases}, 7, false, 0, 0, 0, 80);
		check(nativeIndex.shortlist(query, 3).length==1, "High-cutoff oracle behavior must match native multi-model keep-best fallback");
		NcrnaMappingOracle.main(new String[]{"in="+root.resolve("windows.fa"), "windows="+root.resolve("windows.tsv"), "truth="+root.resolve("truth.tsv"),
			"ref="+root.resolve("ref.fa"), "scores="+root.resolve("t1.models.tsv.gz"), "outprefix="+root.resolve("reuse"), "expected=80", "minlen=20", "quantum=120", "idpass=.65", "t=16"});
		final ByteFile original=RrnaResourceIO.open(root.resolve("t1.curves.tsv").toString()), replay=RrnaResourceIO.open(root.resolve("reuse.curves.tsv").toString());
		final LineParser1 a=new LineParser1('\t'), b=new LineParser1('\t');
		try{
			check(Arrays.equals(original.nextLine(), replay.nextLine()), "Replay curve schema must match fresh verification");
			for(byte[] row; (row=original.nextLine())!=null; ){
				a.set(row);final byte[] next=replay.nextLine();check(next!=null, "Replay must conserve every configuration");b.set(next);
				check(a.terms()==14 && b.terms()==14 && b.parseInt(12)==0, "Reused scores must not be reported as new alignments");
				for(int field=0; field<14; field++){if(field!=12){check(a.termEquals(b.parseString(field), field), "Cached scores changed coverage or shortlist counts");}}
			}
			check(replay.nextLine()==null, "Replay must not append extra configurations");
		}finally{check(!(original.close()|replay.close()), "Replay comparison read failed");}
		capFixture(root, model, bases, query);
		System.out.println("NCRNA_MAPPING_ORACLE_TEST_PASS");
	}
	/** Synthetic cached scores isolate cap eviction; these are not observed aligner outputs. */
	private static void capFixture(Path root, String model, byte[] bases, byte[] query)throws Exception{
		final ByteBuilder refs=new ByteBuilder(), scores=new ByteBuilder();final byte[][] library=new byte[32][];
		for(int i=0; i<32; i++){library[i]=bases;refs.append(">cobalamin_").append(i).nl().append(model).nl();}
		write(root.resolve("capref.fa"), refs.toString());
		write(root.resolve("cap.fa"), ">w0\n"+new String(query, StandardCharsets.US_ASCII)+'\n');
		write(root.resolve("cap.tsv"), NcrnaCalibrationPrep.WINDOW_HEADER+"\nw0\tg\tctg\t0\t1000\t1\t100\t299\tREJECT_ID\t30\t0\t0.6\tfalse\n");
		final int[][] shared=new int[4][32];
		for(int k=6; k<=9; k++){
			final TrnaKmerIndex index=new TrnaKmerIndex(library, k, false, 0, 0, 0, 1);final int[] order=index.shortlist(query, 32);
			check(order.length==32 && order[31]==31, "Native count ties must place the last model below the top30 cap");
			for(int m=0; m<32; m++){shared[k-6][m]=index.lastSharedCount(m);}
		}
		scores.append("window\tmodel\tidentity\taligned_start0\taligned_stop0\tshared6\tshared7\tshared8\tshared9\n");
		for(int m=0; m<32; m++){
			scores.append("w0\t").append(m).tab().append(m==31?"1.0":"0.6").append("\t75\t").append(74+bases.length);
			for(int k=0; k<4; k++){scores.tab().append(shared[k][m]);}scores.nl();
		}
		write(root.resolve("capscores.tsv"), scores.toString());
		NcrnaMappingOracle.main(new String[]{"in="+root.resolve("cap.fa"), "windows="+root.resolve("cap.tsv"), "truth="+root.resolve("truth.tsv"),
			"ref="+root.resolve("capref.fa"), "scores="+root.resolve("capscores.tsv"), "outprefix="+root.resolve("capout"), "expected=1", "minlen=20", "idpass=.65", "t=1"});
		final ByteFile losses=RrnaResourceIO.open(root.resolve("capout.losses.tsv").toString());
		try{
			losses.nextLine();check(new String(losses.nextLine(), StandardCharsets.UTF_8).endsWith("\tCAP_EVICTED"), "An indexed acceptable model below rank30 is cap-evicted, not unseeded");
			check(new String(losses.nextLine(), StandardCharsets.UTF_8).endsWith("\tNO_SCORED_OVERLAPPING_WINDOW"), "A target on the other strand has no scored window, not a mapping failure");
		}finally{check(!losses.close(), "Cap classification read failed");}
		final ByteFile curves=RrnaResourceIO.open(root.resolve("capout.curves.tsv").toString());final LineParser1 p=new LineParser1('\t');int checked=0;
		try{
			curves.nextLine();for(byte[] row; (row=curves.nextLine())!=null; ){
				p.set(row);if(p.parseInt(0)!=7 || p.parseInt(2)!=1){continue;}
				final int top=p.parseInt(1);check(p.parseInt(9)==(top<32?0:1) && p.parseInt(12)==0, "Larger caps recover the accepted model without any new alignment");checked++;
			}
		}finally{check(!curves.close(), "Cap recovery curve read failed");}
		check(checked==7, "Every cap must participate in the controlled K7 comparison");
	}
	private static void write(final Path path, final String text){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
