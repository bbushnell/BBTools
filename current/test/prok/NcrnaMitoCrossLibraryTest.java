package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Random;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMitoCrossLibrary;
import prok.RrnaResourceIO;
import structures.ByteBuilder;
import tax.TaxTree;

/** Exhaustive-search fixture: the exact winner is second, duplicated copies span
 * phyla, and a known internal substitution separates copy/unique weighting.
 * @author Brian Bushnell, Raiden */
public final class NcrnaMitoCrossLibraryTest {

	public static void main(String[] args) throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical hash executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir() && new File(sha).canExecute(),"Fixture directory must be new");
		final String a=sequence(400,17),b=sequence(400,29);final char[] changed=a.toCharArray();changed[200]=changed[200]=='A'?'C':'A';final String mutant=new String(changed);
		final String fasta=record(0,4,"u0",false,false,a)+record(1,999,"u0",true,false,a)+record(2,4,"u1",false,true,mutant)+record(3,4,"u2",false,false,b);
		final String models=">B n=1\n"+b+"\n>A n=3\n"+a+'\n';
		write(root+"/input.fa",fasta);write(root+"/models.fa",models);tree(root+"/tree.tsv");
		final String pins=bind(root,sha,"input.fa","models.fa","tree.tsv");run(root,pins,"serial",1,4,2);run(root,pins,"parallel",3,4,2);
		for(String suffix:new String[]{".members.tsv.gz",".scores.tsv.gz",".summary.tsv",".histogram.tsv"}){check(read(root+"/serial"+suffix).equals(read(root+"/parallel"+suffix)),"Serial/parallel complete outputs must match: "+suffix);}
		final String[] rows=read(root+"/serial.members.tsv.gz").split("\n");check(rows.length==5,"Every original copy must survive exact-sequence caching");
		for(int i=1;i<=3;i++){check(rows[i].split("\t",-1)[8].equals("A n=3"),"The second library model must beat the first for all A copies");}
		check(rows[4].split("\t",-1)[8].equals("B n=1"),"Different queries can choose different models");
		check(Math.abs(Double.parseDouble(rows[3].split("\t",-1)[9])-.9975)<1e-6,"One internal substitution in 400 bases has the hand-known identity");
		check(read(root+"/serial.scores.tsv.gz").split("\n").length==7,"Three unique sequences times two models: no early acceptance or duplicate alignment");
		final String[] all=find(read(root+"/serial.summary.tsv"),"fixture\tALL\tALL\t");
		check(all[3].equals("4") && all[4].equals("3") && all[5].equals("1") && all[6].equals("1"),"Copy, unique, partial and pseudo denominators are independent");
		check(Math.abs(Double.parseDouble(all[7])-.999375)<1e-6 && Math.abs(Double.parseDouble(all[8])-(2+.9975)/3)<1e-6,"Copy and unique means differ by the repeated exact copy");
		for(int i=9;i<12;i++){check(all[i].equals("4") && all[i+3].equals("3"),"Inclusive thresholds retain each expected copy and unique sequence");}
		final String[] phylum=find(read(root+"/serial.summary.tsv"),"fixture\tPHYLUM\t2:Fixture phylum\t"),missing=find(read(root+"/serial.summary.tsv"),"fixture\tPHYLUM\tMISSING_TAXID\t");
		check(phylum[3].equals("3") && phylum[4].equals("3") && missing[3].equals("1") && missing[4].equals("1"),"A unique sequence shared between phyla belongs once in each stratum, once globally");
		final String[] h99=find(read(root+"/serial.histogram.tsv"),"fixture\tALL\tALL\t0.99\t"),h100=find(read(root+"/serial.histogram.tsv"),"fixture\tALL\tALL\t1\t");
		check(h99[5].equals("1") && h99[6].equals("1") && h100[5].equals("3") && h100[6].equals("2"),"Histogram distinguishes the known substitution from exact duplicates");
		check(NcrnaMitoCrossLibrary.thresholdBits(Math.nextDown(.70f))==0 && NcrnaMitoCrossLibrary.thresholdBits(.70f)==1 && NcrnaMitoCrossLibrary.thresholdBits(Math.nextDown(.75f))==1 && NcrnaMitoCrossLibrary.thresholdBits(.75f)==3 && NcrnaMitoCrossLibrary.thresholdBits(Math.nextDown(.80f))==3 && NcrnaMitoCrossLibrary.thresholdBits(.80f)==7,"Each native-float cut is inclusive, not rounded decimal comparison");
		boolean failed=false;try{NcrnaMitoCrossLibrary.thresholdBits(Float.NaN);}catch(IllegalArgumentException e){failed=true;}check(failed,"Nonfinite identities cannot be counted");
		bad(root,pins,"wrong_count",5,2);bad(root,pins,"wrong_models",4,3);
		write(root+"/input.fa",fasta+record(0,4,"u0",false,false,a));bad(root,bind(root,sha,"input.fa","models.fa","tree.tsv"),"duplicate_locus",5,2);
		write(root+"/input.fa",fasta.replace(record(1,999,"u0",true,false,a),record(1,999,"u0",true,false,b)));bad(root,bind(root,sha,"input.fa","models.fa","tree.tsv"),"conflicting_uid",4,2);
		write(root+"/input.fa",fasta);write(root+"/models.fa",models.replace(">A n=3",">B n=3"));bad(root,bind(root,sha,"input.fa","models.fa","tree.tsv"),"duplicate_model_token",4,2);
		write(root+"/models.fa",models);bad(root,bind(root,sha,"input.fa","models.fa"),"missing_tree_pin",4,2);
		final ByteBuilder many=new ByteBuilder();for(int i=0;i<12;i++){final char[] bases=a.toCharArray();bases[100+i]=bases[100+i]=='A'?'C':'A';many.append(record(i,4,"v"+i,false,false,new String(bases)));}write(root+"/input.fa",many.toString());
		final String manyPins=bind(root,sha,"input.fa","models.fa","tree.tsv");run(root,manyPins,"many_serial",1,12,2);run(root,manyPins,"many_parallel",3,12,2);
		for(String suffix:new String[]{".members.tsv.gz",".scores.tsv.gz",".summary.tsv",".histogram.tsv"}){check(read(root+"/many_serial"+suffix).equals(read(root+"/many_parallel"+suffix)),"Multiple independent worker chunks preserve exact output: "+suffix);}
		System.out.println("MITO_CROSS_LIBRARY_TEST_PASS exhaustive=true copy_unique=true taxonomy_gaps=true native_thresholds=true histogram=true parallel=true negatives=7");
	}
	private static void run(String root,String pins,String name,int threads,int expected,int models){NcrnaMitoCrossLibrary.main(new String[]{"in="+root+"/input.fa","ref="+root+"/models.fa","tree="+root+"/tree.tsv","verified="+pins,"outprefix="+root+'/'+name,"arm=fixture","expected="+expected,"expectedmodels="+models,"t="+threads});}
	private static void bad(String root,String pins,String name,int copies,int models){boolean failed=false;try{run(root,pins,name,1,copies,models);}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+'/'+name+".members.tsv.gz").exists(),"Malformed evidence must fail before output: "+name);}
	public static String sequence(int n,long seed){check(n>0,"Fixture sequence must be nonempty");final Random r=new Random(seed);final char[] b=new char[n];for(int i=0;i<n;i++){b[i]="ACGT".charAt(r.nextInt(4));}return new String(b);}
	public static String record(int locus,int tid,String uid,boolean partial,boolean pseudo,String bases){return ">tid|"+tid+"|locus"+locus+" unique="+uid+" partial="+partial+" pseudo="+pseudo+"\n"+bases+'\n';}
	public static void tree(String path){final ByteBuilder b=new ByteBuilder().append("##TaxTree2\t1\n##Encoding\tA48\n##Tree\t5\t4\n##Config\t0\t0\t0\t0\t0\n##MergedMap\t0\t0\n##LevelExtended\tTID\tPID\tNumChildren\tMinParentLevelExtended\tMaxChildLevelExtended\tFlag\tName\n");
		for(int i=1;i<=4;i++){final int level=i==1?TaxTree.LIFE_E:i==2?TaxTree.PHYLUM_E:i==3?TaxTree.CLASS_E:TaxTree.SPECIES_E;b.appendA48(level).tab().appendA48(i).tab().appendA48(i==1?1:i-1).tab().appendA48(i==4?0:1).tab().appendA48(level).tab().appendA48(0).tab().appendA48(0).tab().append(i==1?"Life":i==2?"Fixture phylum":i==3?"Fixture class":"Fixture species").nl();}write(path,b.toString());}
	public static String bind(String root,String sha,String... files) throws Exception{check(files.length>0,"Fixture pins must bind actual input files");final ByteBuilder b=new ByteBuilder();for(String file:files){final String path=root+'/'+file;final Process p=new ProcessBuilder(sha,path).redirectError(ProcessBuilder.Redirect.INHERIT).start();final java.io.ByteArrayOutputStream bytes=new java.io.ByteArrayOutputStream();final byte[] buffer=new byte[128];int n;while((n=p.getInputStream().read(buffer))>=0){bytes.write(buffer,0,n);}check(p.waitFor()==0,"Canonical sha80 tool failed");final String digest=new String(bytes.toByteArray(),StandardCharsets.US_ASCII).trim();check(digest.matches("[0-9a-f]{20}"),"Expected canonical sha80");b.append(digest).tab().append(path).nl();}final String path=root+"/pins_"+(pinCount++)+".tsv";write(path,b.toString());return path;}
	public static String[] find(String text,String prefix){String[] found=null;for(String row:text.split("\n")){if(row.startsWith(prefix)){check(found==null,"Duplicate expected row");found=row.split("\t",-1);}}check(found!=null,"Missing row: "+prefix);return found;}
	public static void write(String path,String text){final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,true,false,false));out.start();out.print(text);RrnaResourceIO.finish(out);}
	public static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{check(!in.close(),"Fixture read failed");}return b.toString();}
	public static void check(boolean b,String message){if(!b){throw new AssertionError(message);}}
	private static int pinCount;
}
