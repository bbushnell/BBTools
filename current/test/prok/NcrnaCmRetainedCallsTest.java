package test.prok;

import java.io.File;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaCmRetainedCalls;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** A changed E-value may match; endpoints, score, model span and duplicates may not.
 * @author Ganyu */
public final class NcrnaCmRetainedCallsTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="),"Require outdir=NEW_DIRECTORY");final String root=args[0].substring(7);require(new File(root).mkdir(),"Fixture output must be new");
		final String a=row("a",1,20,10,29,"50.0","1e-8"),z=row("a",1,20,10,29,"50.0","0.1");
		write(root+"/reference.tbl",a+a+row("b",1,20,40,59,"51.0","1e-9")+row("c",1,20,70,89,"52.0","1e-9")+row("d",1,20,100,119,"53.0","1e-9"));
		write(root+"/candidate.tbl",row("d",2,20,100,119,"53.0","0.01")+row("c",1,20,70,89,"52.1","0.01")+z+row("b",1,20,41,60,"51.0","0.01"));
		final long[] n=NcrnaCmRetainedCalls.compare(root+"/reference.tbl",root+"/candidate.tbl","RF00001",100,40,root+"/changed");
		require(n[0]==5 && n[1]==4 && n[2]==1 && n[3]==4 && n[4]==3,"One E-value-only row matches; a duplicate, endpoint, score and model-span changes must remain visible");
		write(root+"/same.tbl",z+z+row("d",1,20,100,119,"53.0","0.01")+row("b",1,20,40,59,"51.0","0.01")+row("c",1,20,70,89,"52.0","0.01"));
		NcrnaCmRetainedCalls.main(new String[]{"REFERENCE="+root+"/reference.tbl","candidate="+root+"/same.tbl","model=RF00001","clen=100","ga=40","outprefix="+root+"/exact"});
		require(read(root+"/exact.summary.tsv").contains("\n5\t5\t5\t0\t0\ttrue\n"),"Order and E-values alone must not fail exact fixed-field parity");
		write(root+"/empty.tbl","");final long[] empty=NcrnaCmRetainedCalls.compare(root+"/empty.tbl",root+"/empty.tbl","RF00001",100,40,root+"/zero");require(empty[0]==0 && empty[1]==0 && empty[2]==0,"Complete zero-hit arms remain valid comparisons");
		boolean rejected=false;try{NcrnaCmRetainedCalls.compare(root+"/reference.tbl",root+"/same.tbl","RF00002",100,40,root+"/wrong");}catch(IllegalArgumentException e){rejected=true;}require(rejected,"Wrong model table must fail before a success summary");
		System.out.println("CM_RETAINED_CALLS_TEST_PASS evalue_only=true permutation=true duplicates=true endpoint_score_model_span=true zero_hit=true wrong_model_negative=true");
	}
	private static String row(String name,int mf,int mt,int from,int to,String score,String e){return name+" - fixture RF00001 cm "+mf+" "+mt+" "+from+" "+to+" + no 1 0.50 0.0 "+score+" "+e+" ! description\n";}
	private static void write(String path,String rows){final ByteStreamWriter w=RrnaResourceIO.writer(path);try{w.print("# Program: cmsearch\n# Version: 1.1.5 fixture\n"+rows+"# [ok]\n");}finally{RrnaResourceIO.finish(w);}}
	private static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{require(!in.close(),"Fixture read failed");}return b.toString();}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
