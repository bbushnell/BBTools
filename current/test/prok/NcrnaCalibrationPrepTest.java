package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaCalibrationPrep;
import prok.NcrnaCoreRecovery;
import prok.NcrnaStageDiagLogSink;
import prok.RrnaResourceIO;
import stream.Read;
import structures.ByteBuilder;

/** Verifies sense-to-genome CM coordinates, real window slices and ambiguity guards.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCalibrationPrepTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String alphabet="ACGTTGCACCTAGGATCTAA";final ByteBuilder seq=new ByteBuilder();
		for(int i=0; i<10; i++){seq.append(alphabet);}
		write(root.resolve("genome.fa"), ">ctg\n"+seq+'\n');
		write(root.resolve("roster.tsv"), NcrnaCalibrationPrep.ROSTER_HEADER+"\n0\tg0\tasm0\t1\tbacteria\tphy\tunused.gff\t"+root.resolve("genome.fa")+"\t3\t200\t1\n");
		write(root.resolve("audit.tsv"), NcrnaCoreRecovery.AUDIT_HEADER+'\n'+row("plus", '+', true)+'\n'+row("minus", '-', true)+'\n'+row("unknown", '+', false)+'\n');
		write(root.resolve("stages.tsv"), NcrnaStageDiagLogSink.HEADER+'\n'+
			"WINDOW\tcobalamin\tctg\t0\t200\t1\t89\t129\tACCEPT_IDPASS\t1\t2\t0\tcobalamin_0\t0.8\t0\t0\t95\t121\ttrue\t.\t.\n"+
			"WINDOW\tcobalamin\tctg\t1\t200\t1\t70\t110\tACCEPT_IDPASS\t1\t2\t0\tcobalamin_0\t0.8\t0\t0\t76\t102\ttrue\t.\t.\n");
		run(root, "good", 1);
		final ArrayList<String> full=lines(root.resolve("good.full.tsv"));
		check(full.size()==3 && full.get(1).equals("ctg\tcobalamin\t96\t122\t+") && full.get(2).equals("ctg\tcobalamin\t97\t123\t-"),
			"CM projections must use the original sense flank on each strand, preserving27bases");
		check(lines(root.resolve("good.source.tsv")).size()==4 && lines(root.resolve("good.unresolved.tsv")).size()==2,
			"Historical3 and supported2/unresolved1 views must conserve every locus");
		final ArrayList<Read> windows=new ArrayList<Read>();RrnaResourceIO.stream(root.resolve("good.windows.fa.gz").toString(), windows::add);
		final byte[] original=seq.toString().getBytes(StandardCharsets.US_ASCII), plus=Arrays.copyOfRange(original, 89, 130), minus=plus.clone();
		AminoAcid.reverseComplementBasesInPlace(minus);
		check(windows.size()==2 && windows.get(0).id.equals("w0") && windows.get(1).id.equals("w1")
			&& Arrays.equals(windows.get(0).bases, plus) && Arrays.equals(windows.get(1).bases, minus),
			"Window extraction must recover the exact forward and reverse-complemented genomic slices");
		write(root.resolve("duplicate.tsv"), NcrnaCalibrationPrep.ROSTER_HEADER+"\n0\tg0\tasm0\t1\tbacteria\tphy\tunused.gff\t"+root.resolve("genome.fa")+
			"\t3\t200\t1\n1\tg1\tasm1\t2\tbacteria\tphy\tunused.gff\t"+root.resolve("genome.fa")+"\t0\t200\t1\n");
		boolean failed=false;
		try{run(root, "duplicate", 2);}catch(IllegalArgumentException expected){failed=expected.getMessage().contains("globally unique");}
		check(failed && !Files.exists(root.resolve("duplicate.summary.tsv")), "Ambiguous contig ownership must fail before a complete calibration view is published");
		System.out.println("NCRNA_CALIBRATION_PREP_TEST_PASS");
	}
	private static void run(final Path root, final String prefix, final int genomes){
		NcrnaCalibrationPrep.main(new String[]{"roster="+root.resolve(genomes==1?"roster.tsv":"duplicate.tsv"), "audit="+root.resolve("audit.tsv"),
			"stages="+root.resolve("stages.tsv"), "outprefix="+root.resolve(prefix), "family=cobalamin", "genomes="+genomes,
			"loci=3", "windows=2", "contigs=1", "clen=18", "gabits=1"});
	}
	private static String row(final String id, final char strand, final boolean supported){
		final String[] row=new String[34];Arrays.fill(row, "NA");
		row[0]=id;row[1]="g0";row[2]="asm0";row[3]="1";row[4]="bacteria";row[5]="phy";row[6]="ctg";row[7]=String.valueOf(strand);
		row[8]=supported?"100-119":"20-39";row[9]="20";row[10]="1";row[11]=id;row[12]="false";row[13]="false";
		row[14]="ANNOTATION_FULL_COPY";row[15]="200";row[16]="false";row[17]="0";row[18]="5";row[19]="5";row[20]="1";
		row[21]=supported?"CM_FULL_SPAN_SUPPORTED":"NO_REPORTED_HIT";
		if(supported){
			row[22]="1";row[23]="1";row[24]="1";row[25]="18";row[26]="2";row[27]="28";row[28]="no";row[29]="10";
			row[30]="20";row[31]="-4";row[32]="-3";row[33]="ANNOTATION_WITHIN_HIT";
		}
		return String.join("\t", row);
	}
	private static ArrayList<String> lines(final Path path){
		final ArrayList<String> result=new ArrayList<String>();final ByteFile in=RrnaResourceIO.open(path.toString());
		try{for(byte[] row; (row=in.nextLine())!=null; ){result.add(new String(row, StandardCharsets.US_ASCII));}}
		finally{check(!in.close(), "Fixture result read failed");}return result;
	}
	private static void write(final Path path, final String text){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
