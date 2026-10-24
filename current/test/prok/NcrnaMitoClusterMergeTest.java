package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMitoClusterMerge;
import prok.RrnaResourceIO;
import structures.ByteBuilder;
import tax.TaxTree;

/** Hand-constructed single-linkage bridge and composition table, including a
 * sequence duplicated across two original clusters, ties and an unassigned row.
 * @author Ganyu */
public final class NcrnaMitoClusterMergeTest {
	public static void main(String[] args) throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical hash executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir() && new File(sha).canExecute(),"Fixture directory must be new");
		final double[][] matrix={{1,.85,.61,.4},{.85,1,.80,.4},{.61,.80,1,.4},{.4,.4,.4,1}};
		final int[] bridge=NcrnaMitoClusterMerge.components(matrix,.80);check(bridge[0]==bridge[1] && bridge[1]==bridge[2] && bridge[2]!=bridge[3],"A-B-C bridge merges transitively even when A-C is below the cut");
		final int[] edge=NcrnaMitoClusterMerge.components(matrix,.85);check(edge[0]==edge[1] && edge[1]!=edge[2],"Identity equal to the cut joins; a lower edge does not");
		final String[] labels={"A n=2","B n=1","C n=1","D n=2"};
		final int[][] forest=NcrnaMitoClusterMerge.maximumSpanningForest(matrix,labels);check(Arrays.deepEquals(forest,new int[][]{{0,1},{1,2},{0,3}}),"Maximum spanning forest uses AB, BC, then lexical AD among equal .4 edges");
		final int[] permutation={3,1,0,2};final double[][] reordered=new double[4][4];final String[] renamed=new String[4];for(int i=0;i<4;i++){renamed[i]=labels[permutation[i]];for(int j=0;j<4;j++){reordered[i][j]=matrix[permutation[i]][permutation[j]];}}
		final int[][] permuted=NcrnaMitoClusterMerge.maximumSpanningForest(reordered,renamed);for(int i=0;i<3;i++){check(renamed[permuted[i][0]].equals(labels[forest[i][0]]) && renamed[permuted[i][1]].equals(labels[forest[i][1]]),"Named union-edge order must not depend on input model order");}
		for(double cut:new double[]{.6,.65,.7,.75,.8,.85,.9}){final double[][] sparse=new double[4][4];for(int[] e:forest){sparse[e[0]][e[1]]=sparse[e[1]][e[0]]=matrix[e[0]][e[1]];}check(Arrays.equals(NcrnaMitoClusterMerge.components(matrix,cut),NcrnaMitoClusterMerge.components(sparse,cut)),"Saved union edges must reproduce every reported cut component");}
		final ByteBuilder pairs=new ByteBuilder().append(NcrnaMitoClusterMerge.PAIR_HEADER).nl();
		for(int i=0;i<4;i++){for(int j=i;j<4;j++){pairs.append("fixture\tCLUSTER_PAIR\t").append(labels[i]).tab().append(labels[j]).append("\t80\t80\tleft\t").append(matrix[i][j],3).append("\t1\t0\t79\t80\t1\t1\n");}}
		final ByteBuilder members=new ByteBuilder().append(NcrnaMitoClusterMerge.MEMBER_HEADER).nl();
		final String[] unique={"u0","u1","u0","u2","u3","u4","u5"},products={"pA","pA","pB","pB","pD","pD","unassigned"};
		final int[] model={0,0,1,2,3,3,-1},tids={4,4,4,999,4,999,1};
		for(int i=0;i<7;i++){final String[] row=new String[22];Arrays.fill(row,"NA");row[0]="fixture";row[1]="locus"+i;row[2]=unique[i];row[3]=""+tids[i];row[4]="80";row[5]=row[6]="false";row[7]=products[i];row[12]=i<2?"RF00177":i<4?"RF02542":"NA";row[14]=model[i]<0?"orphan":"kept";row[16]=model[i]<0?"UNASSIGNED":labels[model[i]];row[17]=model[i]<0?"NA":"0.9";members.append(String.join("\t",row)).nl();}
		write(root+"/pairs.tsv",pairs.toString());write(root+"/members.tsv",members.toString());tree(root+"/taxonomy.tsv");run(root,bind(root,sha),"ok");
		final String[] rows=read(root+"/ok.tsv").split("\n");check(rows.length==17,"Seven cuts produce five two-group cuts and two top-three cuts, plus header");
		final String[] treeRows=read(root+"/ok.tree.tsv").split("\n");check(treeRows.length==4 && treeRows[1].contains("\tA n=2\tB n=1\t0.850000000\t") && treeRows[3].contains("\tA n=2\tD n=2\t0.400000000\t"),"Persisted merge tree contains actual named union edges and measured identities");
		for(String row:rows){check(row.split("\t",-1).length==26,"Every output row has a stable complete schema");}
		final String[] merged=find(rows,"fixture\t0.80\t2\t1\t");check(merged[4].equals("A n=2") && merged[5].equals("3") && merged[6].equals("4") && merged[7].equals("3"),"Merged bridge has three models, four copies and only three unique sequences across original clusters");
		check(merged[8].equals("7") && merged[9].equals("6") && merged[10].equals("1") && merged[11].equals("1"),"Unassigned source copy/unique denominator remains explicit");
		check(merged[12].equals("pA | pB") && merged[13].equals("2") && merged[14].equals("2"),"Product mode preserves both equal copy-weighted labels and their support");
		check(merged[15].equals("2:Fixture phylum") && merged[16].equals("3") && merged[18].equals("3:Fixture class") && merged[19].equals("3"),"Actual pinned taxonomy gives independently counted phylum and class composition");
		check(merged[21].equals("RF00177 | RF02542") && merged[22].equals("2") && merged[23].equals("2"),"Best-CM family modes retain ties rather than guessing ancestry");
		final String[] high=find(rows,"fixture\t0.90\t4\t1\t");check(high[4].equals("A n=2") && high[6].equals("2") && high[7].equals("2"),"Equal member/unique group sizes use the lexical model key deterministically");
		final String[] intermediate=find(rows,"fixture\t0.85\t3\t1\t");check(intermediate[6].equals("3") && intermediate[7].equals("2"),"Inclusive AB edge merges duplicate unique IDs across models exactly once");
		final String[] gap=find(rows,"fixture\t0.80\t2\t2\t");check(gap[15].equals("2:Fixture phylum | MISSING_TAXID") && gap[17].equals("2"),"Missing taxonomy participates explicitly in dominant-label ties");
		int bad=0;final String firstPair=pairs.toString().split("\n")[1];
		write(root+"/pairs.tsv",pairs.toString()+firstPair+'\n');bad+=bad(root,bind(root,sha),"duplicate_pair");
		write(root+"/pairs.tsv",pairs.toString().replace(firstPair+'\n',""));bad+=bad(root,bind(root,sha),"missing_pair");
		write(root+"/pairs.tsv",pairs.toString());final String row=members.toString().split("\n")[1];
		write(root+"/members.tsv",members.toString()+row+'\n');bad+=bad(root,bind(root,sha),"duplicate_member");
		write(root+"/members.tsv",members.toString().replace(row,row.replace("A n=2","B n=1")));bad+=bad(root,bind(root,sha),"wrong_model_count");
		write(root+"/members.tsv",members.toString());write(root+"/pairs.tsv",pairs.toString().replace("\t0.850\t","\tNaN\t"));bad+=bad(root,bind(root,sha),"nonfinite_identity");
		check(bad==5,"Malformed pair/member evidence fails before a summary is emitted");
		System.out.println("MITO_CLUSTER_MERGE_TEST_PASS bridge=true inclusive_cuts=true copies_unique=true unassigned=true label_ties=true phylum_class=true deterministic_top3=true all_seven_cuts=true named_forest=true forest_cut_parity=true negatives=5");
	}
	private static void run(String root,String pins,String name){NcrnaMitoClusterMerge.main(new String[]{"in="+root+"/members.tsv","pairs="+root+"/pairs.tsv","tree="+root+"/taxonomy.tsv","arm=fixture","expected=7","expectedmodels=4","verified="+pins,"out="+root+'/'+name+".tsv","outtree="+root+'/'+name+".tree.tsv","cuts=0.60,0.65,0.70,0.75,0.80,0.85,0.90","t=1"});}
	private static int bad(String root,String pins,String name){boolean failed=false;try{run(root,pins,name);}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+'/'+name+".tsv").exists() && !new File(root+'/'+name+".tree.tsv").exists(),"Malformed input must not leave a result: "+name);return 1;}
	private static void tree(String path){final ByteBuilder b=new ByteBuilder().append("##TaxTree2\t1\n##Encoding\tA48\n##Tree\t5\t4\n##Config\t0\t0\t0\t0\t0\n##MergedMap\t0\t0\n##LevelExtended\tTID\tPID\tNumChildren\tMinParentLevelExtended\tMaxChildLevelExtended\tFlag\tName\n");
		for(int i=1;i<=4;i++){final int level=i==1?TaxTree.LIFE_E:i==2?TaxTree.PHYLUM_E:i==3?TaxTree.CLASS_E:TaxTree.SPECIES_E;b.appendA48(level).tab().appendA48(i).tab().appendA48(i==1?1:i-1).tab().appendA48(i==4?0:1).tab().appendA48(level).tab().appendA48(0).tab().appendA48(0).tab().append(i==1?"Life":i==2?"Fixture phylum":i==3?"Fixture class":"Fixture species").nl();}write(path,b.toString());}
	private static String bind(String root,String sha) throws Exception{final ByteBuilder b=new ByteBuilder();for(String file:new String[]{"pairs.tsv","members.tsv","taxonomy.tsv"}){final String path=root+'/'+file;final Process p=new ProcessBuilder(sha,path).redirectError(ProcessBuilder.Redirect.INHERIT).start();final java.io.ByteArrayOutputStream bytes=new java.io.ByteArrayOutputStream();final byte[] buffer=new byte[128];int n;while((n=p.getInputStream().read(buffer))>=0){bytes.write(buffer,0,n);}check(p.waitFor()==0,"Fixture canonical hash failed");final String digest=new String(bytes.toByteArray(),StandardCharsets.US_ASCII).trim();check(digest.matches("[0-9a-f]{20}"),"Expected canonical sha80");b.append(digest).tab().append(path).nl();}final String path=root+"/pins_"+(pinCount++)+".tsv";write(path,b.toString());return path;}
	private static String[] find(String[] rows,String prefix){String[] found=null;for(String row:rows){if(row.startsWith(prefix)){check(found==null,"Duplicate expected summary row");found=row.split("\t",-1);}}check(found!=null,"Missing summary row: "+prefix);return found;}
	private static void write(String path,String text){final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,true,false,false));out.start();out.print(text);RrnaResourceIO.finish(out);}
	private static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{check(!in.close(),"Fixture output read failed");}return b.toString();}
	private static void check(boolean b,String message){if(!b){throw new AssertionError(message);}}
	private static int pinCount;
}
