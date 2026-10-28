package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.CobalaminShredSampler;
import prok.CobalaminShredScreen;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Core-overlap exclusion is strand independent; context-only/partial/low/truncated hits stay explicit.
 * @author Raiden
 */
public final class CobalaminShredScreenTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture root");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final ByteBuilder manifest=new ByteBuilder(CobalaminShredSampler.HEADER+'\n'), cores=new ByteBuilder(), scanned=new ByteBuilder();
		for(int i=0; i<8; i++){
			manifest.append('q').append(i).append("\tpanel\tg\tGCF_1.1\t1\tc\t151\t170\t101\t200\t51\t70\t20\n");
			cores.append(">q").append(i).append("\nACGTACGTACGTACGTACGT\n");if(i<7){scanned.append('q').append(i).nl();}
		}
		write(root.resolve("manifest.tsv"), manifest.toString());write(root.resolve("cores.fa"), cores.toString());write(root.resolve("scanned.txt"), scanned.toString());
		final String footer="# Program: cmsearch\n# Version: 1.1.5 (Nov 2020)\n# [ok]\n";
		write(root.resolve("rf1.tbl"), hit("q0", "RF00174", 1, 10, 45, 54, '+', "no", 50)+hit("q2", "RF00174", 1, 10, 10, 19, '+', "no", 50)
			+hit("q3", "RF00174", 2, 10, 51, 59, '+', "no", 50)+hit("q4", "RF00174", 1, 10, 51, 60, '+', "no", 10)
			+hit("q5", "RF00174", 1, 10, 51, 60, '+', "5'", 50)+footer);
		write(root.resolve("rf14.tbl"), hit("q1", "RF01482", 1, 8, 75, 68, '-', "no", 50)+footer);
		write(root.resolve("rf16.tbl"), hit("q0", "RF01689", 1, 6, 69, 74, '+', "no", 50)+footer);
		final String models=CobalaminShredScreen.MODEL_HEADER+'\n'+"RF00174\t10\t20\t"+root.resolve("rf1.tbl")+"\t00000000000000000001\n"+
			"RF01482\t8\t20\t"+root.resolve("rf14.tbl")+"\t00000000000000000002\n"+"RF01689\t6\t20\t"+root.resolve("rf16.tbl")+"\t00000000000000000003\n";
		write(root.resolve("models.tsv"), models);run(root, "good", "models.tsv", "scanned.txt");
		final ByteFile in=RrnaResourceIO.open(root.resolve("good.disposition.tsv").toString());final LineParser1 p=new LineParser1('\t');int rows=0, excluded=0;
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);check(p.terms()==6, "Disposition must retain every original core and explicit evidence status");
			if(rows<2){check(p.termStartsWith("EXCLUDED", 2), "Full-span hits crossing either core edge must exclude on either strand");excluded++;}
			else if(rows==7){check(p.termEquals("SOURCE_MASKED_BACKGROUND_NOT_CM_ADJUDICATED", 2) && p.parseInt(4)==0, "Unscanned must never become screened-negative by absence of a hit");}
			else{check(p.termEquals("RETAINED_CM_SCREENED_BACKGROUND", 2), "Context-only, partial, belowGA, truncated and zero-hit cores remain distinct retained evidence");}
			if(rows==0){check(p.termEquals("RF00174,RF01689", 3), "All qualifying clan models must remain visible without double-excluding a core");}rows++;
		}}finally{check(!in.close(), "Disposition read failed");}
		check(rows==8 && excluded==2, "Screen must conserve all eight cores and exactly two exclusions");
		final int[] kept={0};RrnaResourceIO.stream(root.resolve("good.cores.fa.gz").toString(), r->{check(!r.id.equals("q0") && !r.id.equals("q1") && r.length()==20, "Only excluded cores leave the retained FASTA");kept[0]++;});
		check(kept[0]==6, "Retained FASTA denominator must include the unadjudicated core");
		final ByteFile summary=RrnaResourceIO.open(root.resolve("good.summary.tsv").toString());
		try{summary.nextLine();p.set(summary.nextLine());check(p.parseInt(0)==8 && p.parseInt(1)==160 && p.parseInt(2)==6 && p.parseInt(3)==120
			&& p.parseInt(4)==2 && p.parseInt(5)==7 && p.parseInt(6)==7 && p.parseInt(7)==1, "Summary must separately conserve cores, bases, hits, scanned and unadjudicated populations");}
		finally{check(!summary.close(), "Summary read failed");}
		write(root.resolve("unscanned.tbl"), hit("q7", "RF00174", 1, 10, 51, 60, '+', "no", 50)+footer);
		write(root.resolve("badmodels.tsv"), models.replace(root.resolve("rf1.tbl").toString(), root.resolve("unscanned.tbl").toString()));reject(root, "outside_subset", "badmodels.tsv");
		write(root.resolve("incomplete.tbl"), "# Program: cmsearch\n# Version: 1.1.5 (Nov 2020)\n");
		write(root.resolve("incomplete_models.tsv"), models.replace(root.resolve("rf1.tbl").toString(), root.resolve("incomplete.tbl").toString()));reject(root, "unfinished_scan", "incomplete_models.tsv");
		write(root.resolve("missing_models.tsv"), models.substring(0, models.indexOf("RF01689\t")));reject(root, "missing_model", "missing_models.tsv");
		write(root.resolve("empty.ids"), "");write(root.resolve("not_scanned.tsv"), CobalaminShredScreen.MODEL_HEADER+"\nRF00174\t10\t20\t.\tNA\nRF01482\t8\t20\t.\tNA\nRF01689\t6\t20\t.\tNA\n");
		CobalaminShredScreen.main(new String[]{"manifest="+root.resolve("manifest.tsv"), "in="+root.resolve("cores.fa"), "models="+root.resolve("not_scanned.tsv"),
			"scanned="+root.resolve("empty.ids"), "expectedscanned=0", "outprefix="+root.resolve("none"), "expected=8"});
		final ByteFile none=RrnaResourceIO.open(root.resolve("none.summary.tsv").toString());
		try{none.nextLine();p.set(none.nextLine());check(p.parseInt(2)==8 && p.parseInt(4)==0 && p.parseInt(5)==0 && p.parseInt(6)==0 && p.parseInt(7)==8
			&& p.termStartsWith("CM_NOT_RUN_EMPTY_SELECTION", 8), "Empty selection must retain all inputs as unadjudicated and explicitly report that no CM ran");}
		finally{check(!none.close(), "Empty-selection summary read failed");}
		System.out.println("COBALAMIN_SHRED_SCREEN_TEST_PASS");
	}
	private static String hit(String id, String model, int a, int z, int from, int to, char strand, String trunc, int bits){
		return id+" - model "+model+" cm "+a+' '+z+' '+from+' '+to+' '+strand+' '+trunc+" 1 0.50 0.0 "+bits+" 1e-8 ! fixture\n";
	}
	private static void run(Path root, String out, String models, String scanned){
		CobalaminShredScreen.main(new String[]{"manifest="+root.resolve("manifest.tsv"), "in="+root.resolve("cores.fa"), "models="+root.resolve(models),
			"scanned="+root.resolve(scanned), "expectedscanned=7", "outprefix="+root.resolve(out), "expected=8"});
	}
	private static void reject(Path root, String out, String models){
		boolean failed=false;try{run(root, out, models, "scanned.txt");}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve(out+".summary.tsv")), "Incomplete or misattributed CM evidence must fail before publishing dispositions");
	}
	private static void write(Path path, String text){final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
