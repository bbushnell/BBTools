package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaNullIdentity;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Exact forward/reverse controls on both sides of caller dispatch, and deterministic parallel output.
 * @author Raiden
 */
public final class NcrnaNullIdentityTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String shortA=bases('A', 120), longC=bases('C', 121), shortT=bases('T', 120), longG=bases('G', 121);
		write(root.resolve("models.fa"), ">a\n"+shortA+"\n>c\n"+longC+'\n');
		final ByteBuilder input=new ByteBuilder();for(int i=0; i<20; i++){
			input.append(">q").append(i).nl();switch(i%5){case 0:input.append(shortA);break;case 1:input.append(shortT);break;case 2:input.append(longC);break;case 3:input.append(longG);break;default:for(int j=0; j<121; j++){input.append("AGAT".charAt(j&3));}}
			input.nl();
		}
		write(root.resolve("input.fa"), input.toString());run(root, "serial", "input.fa", 1, 20, 2);run(root, "parallel", "input.fa", 3, 20, 2);
		for(String suffix:new String[]{".loci.tsv", ".quantiles.tsv", ".thresholds.tsv", ".cm_review.ids.txt", ".summary.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("serial"+suffix)), Files.readAllBytes(root.resolve("parallel"+suffix))), "Worker count must not change scores, order, quantiles or denominator: "+suffix);
		}
		final ByteFile in=RrnaResourceIO.open(root.resolve("serial.loci.tsv").toString());final LineParser1 p=new LineParser1('\t');int count=0;
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);check(p.terms()==11 && p.termEquals("q"+count, 0) && p.termEquals(NcrnaNullIdentity.INTERPRETATION, 10), "Each input needs one explicitly seed-bypassed result in original order");
			final float forward=Float.parseFloat(p.parseString(3)), reverse=Float.parseFloat(p.parseString(4)), best=Float.parseFloat(p.parseString(5));
			check(best==Math.max(forward, reverse), "Best-of-both orientation score must agree with the reported individual maxima");
			if(count%5<4){check(best==1f, "An exact full model in either orientation must score one across short/long dispatch");}
			else{check(Float.isFinite(best) && best<1f, "A mixed sequence cannot become an exact homopolymer model");}
			if(count%5==1 || count%5==3){check(p.termEquals("-", 6) && reverse==1f, "Reverse-oriented exact controls must retain their winning orientation");}count++;
		}}finally{check(!in.close(), "Identity fixture read failed");}
		check(count==20, "Every sample must be scored exactly once");
		final ByteFile summary=RrnaResourceIO.open(root.resolve("serial.summary.tsv").toString());
		try{summary.nextLine();p.set(summary.nextLine());check(p.parseInt(1)==20 && p.parseInt(3)==2 && p.parseInt(4)==2 && p.parseInt(5)==80 && p.parseInt(6)==120,
			"No seed or shortlist may suppress any of the80 required initial alignments");}
		finally{check(!summary.close(), "Summary fixture read failed");}
		write(root.resolve("duplicate.fa"), ">same\n"+shortA+"\n>same\n"+shortT+'\n');reject(root, "duplicate", "duplicate.fa", 2, 2);
		reject(root, "missing", "input.fa", 21, 2);reject(root, "models", "input.fa", 20, 3);
		System.out.println("NCRNA_NULL_IDENTITY_TEST_PASS");
	}
	private static void run(Path root, String out, String input, int threads, int count, int models){
		NcrnaNullIdentity.main(new String[]{"in="+root.resolve(input), "ref="+root.resolve("models.fa"), "outprefix="+root.resolve(out),
			"label=NULL_PANEL", "expected="+count, "expectedmodels="+models, "t="+threads, "quantum=120"});
	}
	private static void reject(Path root, String out, String input, int count, int models){
		boolean failed=false;try{run(root, out, input, 1, count, models);}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve(out+".summary.tsv")), "Wrong denominators or duplicate inputs must fail before publishing a distribution");
	}
	private static String bases(char c, int n){final ByteBuilder b=new ByteBuilder();for(int i=0; i<n; i++){b.append(c);}return b.toString();}
	private static void write(Path path, String text){final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
