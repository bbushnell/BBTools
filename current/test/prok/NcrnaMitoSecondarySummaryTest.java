package test.prok;

import java.io.File;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMitoSecondarySummary;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Hand means distinguish copy weighting, unique weighting and overlapping flags.
 * @author Brian, Raiden */
public final class NcrnaMitoSecondarySummaryTest {
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="),"Use one fresh fixture directory");final String root=args[0].substring(7);check(new File(root).mkdir(),"Fixture directory must not exist");
		final String rows="a\tx\t100\tfalse\tfalse\t0.2\t100\t0\t99\n"+"b\tx\t100\ttrue\tfalse\t0.2\t100\t0\t99\n"+"c\ty\t100\tfalse\tfalse\t0.8\t100\t0\t99\n"+"d\tz\t100\ttrue\ttrue\t0.9\t100\t0\t99\n";
		write(root+"/members.tsv",rows);NcrnaMitoSecondarySummary.summarize(root+"/members.tsv",root+"/summary.tsv","fixture",4);
		final String[] result=read(root+"/summary.tsv").split("\n");check(result.length==5,"Four strata plus header required");
		row(result[1],"ALL",4,3,.525,1.9/3);row(result[2],"SOURCE_NONPARTIAL_NONPSEUDO",2,2,.5,.5);row(result[3],"PARTIAL",2,2,.55,.55);row(result[4],"PSEUDO",1,1,.9,.9);
		int bad=0;bad+=bad(root,rows+rows.split("\n")[0]+'\n',"duplicate",5);
		bad+=bad(root,rows.replace("b\tx\t100\ttrue\tfalse\t0.2","b\tx\t100\ttrue\tfalse\t0.3"),"identity_conflict",4);
		bad+=bad(root,rows.replace("\t0.8\t","\tNaN\t"),"nonfinite",4);bad+=bad(root,rows.replace("\ttrue\ttrue\t","\tmaybe\ttrue\t"),"flag",4);bad+=bad(root,rows,"count",5);
		check(bad==5,"Invalid evidence must fail before emitting aggregate means");
		System.out.println("MITO_SECONDARY_SUMMARY_TEST_PASS copy_unique=true overlapping_flags=true exact_denominators=true negatives=5");
	}
	private static void row(String text,String name,int copies,int unique,double copyMean,double uniqueMean){final String[] t=text.split("\t");check(t.length==7 && t[1].equals(name) && Integer.parseInt(t[2])==copies && Integer.parseInt(t[3])==unique,"Hand counts/stratum differ");check(Math.abs(Double.parseDouble(t[4])-copyMean)<1e-9 && Math.abs(Double.parseDouble(t[5])-uniqueMean)<1e-9,"Copy/unique mean differs from hand arithmetic");}
	private static int bad(String root,String rows,String name,int count){write(root+"/members.tsv",rows);boolean failed=false;try{NcrnaMitoSecondarySummary.summarize(root+"/members.tsv",root+'/'+name+".tsv","fixture",count);}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+'/'+name+".tsv").exists(),"Bad evidence must leave no aggregate: "+name);return 1;}
	private static void write(String path,String rows){final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,true,false,false));out.start();out.println(NcrnaMitoSecondarySummary.HEADER);out.print(rows);RrnaResourceIO.finish(out);}
	private static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{check(!in.close(),"Fixture read failed");}return b.toString();}
	private static void check(boolean ok,String message){if(!ok){throw new AssertionError(message);}}
}
