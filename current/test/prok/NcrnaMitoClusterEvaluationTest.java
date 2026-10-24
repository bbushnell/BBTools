package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMitoCharacterization;
import prok.NcrnaMitoClusterEvaluation;
import prok.RrnaResourceIO;
import structures.ByteBuilder;
import tax.TaxTree;

/** Complete tiny arm: duplicate copies, overlapping source flags, weak assignment,
 * orphan, explicit taxonomy gaps, distinct files sharing a reference ID, and
 * hand-calculated summary denominators. No reference ancestry is inferred.
 * @author Ganyu */
public final class NcrnaMitoClusterEvaluationTest {
	public static void main(String[] args) throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical sha80 executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir() && new File(sha).canExecute(),"Fixture directory must be new and hashing executable");
		final String a=dna(80,17),b=dna(80,29);final String c=(a.charAt(0)=='A'?'C':'A')+a.substring(1);final String weak=repeat('C',80);
		final String[] ids={"a","b","c","d","e","f"},unique={"uA","uA","uC","uD","uE","uF"},bases={a,a,c,b,a.substring(0,30),weak};
		final int[] tids={3,3,3,999,1,3};final boolean[] partial={false,true,false,false,true,false},pseudo={false,false,false,true,true,false};
		final String[] headers=new String[6];final ByteBuilder fa=new ByteBuilder(),audit=new ByteBuilder().append(NcrnaMitoCharacterization.HEADER).nl(),assign=new ByteBuilder().append("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n");
		for(int i=0;i<6;i++){headers[i]=ids[i]+" unique="+unique[i]+" tid|"+tids[i]+" partial="+partial[i]+" pseudo="+pseudo[i];fa.append('>').append(headers[i]).nl().append(bases[i]).nl();
			final String model=i==4?"":i==3?"clusterB n=1":"clusterA n=4",status=i==4?"orphan":"kept";
			assign.append(headers[i]).tab().append(tids[i]).append("\tunknown\tkey_").append(i==4?"orphan":model).tab().append(model).tab().append(status).nl();
			final String[] row=new String[30];Arrays.fill(row,"NA");row[0]="ANNOTATION";row[1]=ids[i];row[3]=""+tids[i];row[6]=""+bases[i].length();row[7]=""+partial[i];row[8]=""+pseudo[i];row[9]="source_product";row[12]=unique[i];row[15]="16S";row[18]="RF00177";row[23]="CONCORDANT_NOMINAL_MODEL_CLASS";audit.append(String.join("\t",row)).nl();
		}
		write(root+"/arm.fa",fa.toString());write(root+"/audit.tsv",audit.toString());write(root+"/assignments.tsv",assign.toString());
		write(root+"/clusters.fa",">clusterA n=4\n"+a+"\n>clusterB n=1\n"+b+'\n');
		final String hbm="#HBM\t1\n"+model("clusterA n=4",a)+model("clusterB n=1",b);write(root+"/models.hbm",hbm);
		write(root+"/r16.fa",">universal\n"+a+'\n');write(root+"/r18.fa",">euk18S_universal\n"+b+'\n');write(root+"/r23.fa",">universal\n"+b+'\n');write(root+"/rlsu.fa",">lsu_universal_final\n"+a+'\n');writeTree(root+"/taxonomy.tsv");
		final String pins=bind(root,sha);run(root,pins,"serial",1);run(root,pins,"parallel",3);
		for(String suffix:new String[]{".members.tsv.gz",".pairs.tsv.gz",".summary.tsv",".composition.tsv"}){check(read(root+"/serial"+suffix).equals(read(root+"/parallel"+suffix)),"Thread count must preserve every ordered output: "+suffix);}
		final String[] members=read(root+"/serial.members.tsv.gz").split("\n"),pairs=read(root+"/serial.pairs.tsv.gz").split("\n"),summaries=read(root+"/serial.summary.tsv").split("\n");
		check(members.length==7 && pairs.length==12,"Six original copies and 2*4 universals+3 unordered/self pairs must survive");
		final String[] all=find(summaries,"fixture\tARM\tALL\tALL\t");check(all[5].equals("6") && all[6].equals("5") && all[7].equals("1") && all[8].equals("5") && all[9].equals("4") && all[10].equals("1"),"Copy/unique assigned and unassigned denominators are independently known");
		final String[] mc=find(members,"fixture\tc\t"),mf=find(members,"fixture\tf\t");final double ci=Double.parseDouble(mc[17]),wi=Double.parseDouble(mf[17]);
		near(all[11],(3+ci+wi)/5);near(all[12],(2+ci+wi)/4);
		final String[] lineage=find(summaries,"fixture\tLINEAGE\t2:Fixture phylum\tALL\t");
		check(lineage[5].equals("4") && lineage[6].equals("4") && lineage[7].equals("0") && lineage[8].equals("3") && lineage[9].equals("3") && lineage[14].equals("phylum"),"Lineage quality rows retain copy/unique denominators and explicit rank");
		near(lineage[11],(2+ci+wi)/4);near(lineage[12],(1+ci+wi)/3);
		final String[] missing=find(summaries,"fixture\tLINEAGE\tMISSING_TAXID\tALL\t"),unranked=find(summaries,"fixture\tLINEAGE\tRANK_ABSENT\tALL\t");
		check(missing[5].equals("1") && missing[6].equals("1") && missing[8].equals("1"),"Unknown source taxids remain in lineage assignment and identity denominators");near(missing[11],1);near(missing[12],1);
		check(unranked[5].equals("1") && unranked[6].equals("0") && unranked[7].equals("1") && unranked[8].equals("1") && unranked[10].equals("1") && unranked[11].equals("NA") && unranked[12].equals("NA"),"Unassigned lineage rows preserve the denominator without invented identity");
		check(find(summaries,"fixture\tARM\tALL\tSOURCE_NONPARTIAL_NONPSEUDO\t")[5].equals("3"),"Source nonpartial/nonpseudo is a flag stratum, not inferred intactness");
		check(find(summaries,"fixture\tARM\tALL\tPARTIAL\t")[5].equals("2") && find(summaries,"fixture\tARM\tALL\tPSEUDO\t")[5].equals("2"),"Partial and pseudo overlap explicitly; neither is discarded");
		final String[] orphan=find(members,"fixture\te\t");check(orphan[14].equals("orphan") && orphan[16].equals("UNASSIGNED") && orphan[17].equals("NA") && orphan[10].equals("RANK_ABSENT"),"Orphan and absent rank remain explicit without fabricated identity");
		check(find(members,"fixture\td\t")[10].equals("MISSING_TAXID") && find(members,"fixture\ta\t")[10].equals("2:Fixture phylum"),"Pinned taxonomy resolves a real parent chain and reports unknown taxids");
		check(wi<.8 && mf[14].equals("kept"),"A weak member remains in its original assignment without a new threshold");
		near(find(pairs,"fixture\tUNIVERSAL\tclusterA n=4\t16S:universal\t")[7],1);
		check(Double.parseDouble(find(pairs,"fixture\tUNIVERSAL\tclusterA n=4\t23S:universal\t")[7])<1,"Two literal universal IDs in different files must not collapse");
		int bad=0;
		write(root+"/assignments.tsv",assign.toString()+assign.toString().split("\n")[1]+'\n');bad+=bad(root,bind(root,sha),"duplicate_assignment");
		write(root+"/assignments.tsv",assign.toString().replace(assign.toString().split("\n")[1]+'\n',""));bad+=bad(root,bind(root,sha),"missing_assignment");
		write(root+"/assignments.tsv",assign.toString().replace("clusterA n=4","missing_model"));bad+=bad(root,bind(root,sha),"unknown_model");
		final String cAssignment=assign.toString().split("\n")[3];
		write(root+"/assignments.tsv",assign.toString().replace(cAssignment,cAssignment.replace("clusterA n=4","clusterB n=1")));bad+=bad(root,bind(root,sha),"valid_model_wrong_count");
		write(root+"/assignments.tsv",assign.toString());write(root+"/clusters.fa",">clusterA n=4\n"+a+"\n>clusterA n=1\n"+b+'\n');bad+=bad(root,bind(root,sha),"duplicate_model_token");
		write(root+"/clusters.fa",">clusterA n=4\n"+a+"\n>clusterB n=1\n"+b+'\n');
		write(root+"/assignments.tsv",assign.toString());write(root+"/models.hbm",hbm.replaceFirst(a,c));bad+=bad(root,bind(root,sha),"hbm_bases");
		write(root+"/models.hbm",hbm);write(root+"/arm.fa",fa.toString().replace(">"+headers[1]+"\n"+a,">"+headers[1]+"\n"+c));bad+=bad(root,bind(root,sha),"unique_collision");
		check(bad==7,"Malformed memberships/model counts/coordinates must fail before reports are created");
		System.out.println("MITO_CLUSTER_EVALUATION_TEST_PASS conservation=true model_counts=true model_tokens=true copy_unique=true lineage_means=true overlapping_flags=true weak_retained=true unassigned_NA=true taxonomy_gaps=true universal_ids=true all_pairs=true parallel_exact=true negatives=7");
	}
	private static void run(String root,String pins,String name,int threads){NcrnaMitoClusterEvaluation.main(new String[]{"in="+root+"/arm.fa","assignments="+root+"/assignments.tsv","consensus="+root+"/clusters.fa","hbm="+root+"/models.hbm","audit="+root+"/audit.tsv","tree="+root+"/taxonomy.tsv","lineagerank=phylum","arm=fixture","expected=6","verified="+pins,"outprefix="+root+'/'+name,"ref16="+root+"/r16.fa","ref16id=universal","ref18="+root+"/r18.fa","ref18id=euk18S_universal","ref23="+root+"/r23.fa","ref23id=universal","reflsu="+root+"/rlsu.fa","reflsuid=lsu_universal_final","lengthbins=50,100","t="+threads});}
	private static int bad(String root,String pins,String name){boolean failed=false;try{run(root,pins,name,1);}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+'/'+name+".members.tsv.gz").exists(),"Invalid input cannot leave a valid-looking report: "+name);return 1;}
	private static String model(String id,String bases){final ByteBuilder b=new ByteBuilder().append('>').append(id).append("\tdepth=1\n").append(bases).nl();for(int i=0;i<bases.length();i++){b.append("1\t1\t1\t1\t0\n");}return b.toString();}
	private static void writeTree(String path){final ByteBuilder b=new ByteBuilder().append("##TaxTree2\t1\n##Encoding\tA48\n##Tree\t4\t3\n##Config\t0\t0\t0\t0\t0\n##MergedMap\t0\t0\n##LevelExtended\tTID\tPID\tNumChildren\tMinParentLevelExtended\tMaxChildLevelExtended\tFlag\tName\n");
		for(int i=1;i<=3;i++){final int level=i==1?TaxTree.LIFE_E:i==2?TaxTree.PHYLUM_E:TaxTree.SPECIES_E;b.appendA48(level).tab().appendA48(i).tab().appendA48(i==1?1:i-1).tab().appendA48(i==3?0:1).tab().appendA48(level).tab().appendA48(0).tab().appendA48(0).tab().append(i==1?"Life":i==2?"Fixture phylum":"Fixture species").nl();}write(path,b.toString());}
	private static String bind(String root,String sha) throws Exception{final ByteBuilder b=new ByteBuilder();for(String file:new String[]{"arm.fa","audit.tsv","assignments.tsv","clusters.fa","models.hbm","r16.fa","r18.fa","r23.fa","rlsu.fa","taxonomy.tsv"}){
		final String path=root+'/'+file;final Process p=new ProcessBuilder(sha,path).redirectError(ProcessBuilder.Redirect.INHERIT).start();final java.io.ByteArrayOutputStream bytes=new java.io.ByteArrayOutputStream();final byte[] buffer=new byte[128];int n;while((n=p.getInputStream().read(buffer))>=0){bytes.write(buffer,0,n);}check(p.waitFor()==0,"Canonical fixture hash failed");final String digest=new String(bytes.toByteArray(),StandardCharsets.US_ASCII).trim();check(digest.matches("[0-9a-f]{20}"),"Expected sha80 digest");b.append(digest).tab().append(path).nl();}
		final String path=root+"/pins_"+(pinsCreated++)+".tsv";write(path,b.toString());return path;}
	private static String[] find(String[] rows,String prefix){String[] found=null;for(String row:rows){if(row.startsWith(prefix)){check(found==null,"Duplicate expected row");found=row.split("\t",-1);}}check(found!=null,"Missing expected row: "+prefix);return found;}
	private static void near(String value,double expected){check(Math.abs(Double.parseDouble(value)-expected)<1e-8,"Incorrect hand-calculated mean: "+value+" vs "+expected);}
	private static String dna(int n,long seed){final Random r=new Random(seed);final StringBuilder b=new StringBuilder();for(int i=0;i<n;i++){b.append("ACGT".charAt(r.nextInt(4)));}return b.toString();}
	private static String repeat(char c,int n){final char[] b=new char[n];Arrays.fill(b,c);return new String(b);}
	private static void write(String path,String text){final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,true,false,false));out.start();out.print(text);RrnaResourceIO.finish(out);}
	private static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final ByteBuilder b=new ByteBuilder();try{for(byte[] row;(row=in.nextLine())!=null;){b.append(row).nl();}}finally{check(!in.close(),"Fixture output read failed");}return b.toString();}
	private static void check(boolean ok,String message){if(!ok){throw new AssertionError(message);}}
	private static int pinsCreated;
}
