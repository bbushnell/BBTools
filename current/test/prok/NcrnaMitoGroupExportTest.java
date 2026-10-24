package test.prok;

import static test.prok.NcrnaMitoCrossLibraryTest.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import prok.NcrnaMitoGroupExport;
import prok.NcrnaMitoClusterMerge;
import prok.RrnaResourceIO;
import stream.Read;
import structures.ByteBuilder;

/** Exact connected-component export: transitive bridge, minority taxa, duplicate
 * sequence copies, source flags and original headers all survive.
 * @author Brian Bushnell, Raiden */
public final class NcrnaMitoGroupExportTest {

	public static void main(String[] args) throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical hash executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir() && new File(sha).canExecute(),"Fixture directory must be new");
		final String[] names={"A n=3","B n=1","C n=1","D n=1"},uids={"u0","u1","u2","u0","u3","u4"};final int[] model={0,0,0,1,2,3},tid={4,4,4,999,999,4};
		final double[][] matrix={{1,.8,.6,.3},{.8,1,.7,.3},{.6,.7,1,.3},{.3,.3,.3,1}};
		final ByteBuilder fasta=new ByteBuilder(),models=new ByteBuilder(),members=new ByteBuilder().append(NcrnaMitoClusterMerge.MEMBER_HEADER).nl(),pairs=new ByteBuilder().append(NcrnaMitoClusterMerge.PAIR_HEADER).nl();
		for(int i=0;i<6;i++){final String seq=sequence(80,i==3?17:17+i);fasta.append(record(i,tid[i],uids[i],i==3,i==4,seq));final String[] row=new String[22];Arrays.fill(row,"NA");row[0]="fixture";row[1]="tid|"+tid[i]+"|locus"+i;row[2]=uids[i];row[3]=""+tid[i];row[4]="80";row[5]=""+(i==3);row[6]=""+(i==4);row[14]="kept";row[16]=names[model[i]];row[17]="0.9";members.append(String.join("\t",row)).nl();}
		// Deliberately nonlexical FASTA order; canonical component key is still A.
		for(int i:new int[]{3,1,0,2}){models.append('>').append(names[i]).nl().append(sequence(80,17+i)).nl();}
		for(int i=0;i<4;i++){for(int j=i;j<4;j++){pairs.append("fixture\tCLUSTER_PAIR\t").append(names[i]).tab().append(names[j]).append("\t80\t80\tleft\t").append(matrix[i][j],3).append("\t1\t0\t79\t80\t1\t1\n");}}
		write(root+"/input.fa",fasta.toString());write(root+"/models.fa",models.toString());write(root+"/members.tsv",members.toString());write(root+"/pairs.tsv",pairs.toString());
		String pins=bind(root,sha,FILES);run(root,pins,"ok","A",5);
		final ArrayList<Read> source=new ArrayList<Read>(),exported=new ArrayList<Read>(),seeds=new ArrayList<Read>();RrnaResourceIO.stream(root+"/input.fa",source::add);RrnaResourceIO.stream(root+"/ok.fa.gz",exported::add);RrnaResourceIO.stream(root+"/ok.seeds.fa.gz",seeds::add);
		check(exported.size()==5 && seeds.size()==3,"ABC bridge includes every source copy and all three native seeds");for(int i=0;i<5;i++){check(exported.get(i).id.equals(source.get(i).id) && Arrays.equals(exported.get(i).bases,source.get(i).bases),"Original full headers, order and sequences remain exact");}
		final String[] row=find(read(root+"/ok.summary.tsv"),"group\tfixture\t");check(row[2].equals("A n=3") && row[5].equals("5") && row[6].equals("4") && row[7].equals("3") && row[8].equals("1") && row[9].equals("1") && row[13].equals("A n=3") && row[14].equals("false"),"Selected counts, duplicated UID, flags and seed choice are known independently");
		check(read(root+"/ok.members.tsv").contains("tid|999|locus4"),"Minority taxonomy is included; group label is not a taxon filter");
		bad(root,pins,"noncanonical","B",5);bad(root,pins,"wrong_selected","A",4);
		final String first=pairs.toString().split("\n")[1];write(root+"/pairs.tsv",pairs.toString().replace(first+'\n',""));bad(root,bind(root,sha,FILES),"missing_pair","A",5);
		write(root+"/pairs.tsv",pairs.toString()+first+'\n');bad(root,bind(root,sha,FILES),"duplicate_pair","A",5);write(root+"/pairs.tsv",pairs.toString());
		final String member=members.toString().split("\n")[1];write(root+"/members.tsv",members.toString()+member+'\n');bad(root,bind(root,sha,FILES),"duplicate_member","A",5);
		write(root+"/members.tsv",members.toString().replace(member,member.replace("A n=3","B n=1")));bad(root,bind(root,sha,FILES),"wrong_assignment_count","A",5);write(root+"/members.tsv",members.toString());
		write(root+"/input.fa",fasta.toString().replace("partial=true","partial=false"));bad(root,bind(root,sha,FILES),"source_flag_mismatch","A",5);
		System.out.println("MITO_GROUP_EXPORT_TEST_PASS exact_component=true transitive=true copies_unique=true minority_taxa=true flags=true full_headers=true native_seeds=true negatives=7");
	}
	private static void run(String root,String pins,String name,String target,int selected){NcrnaMitoGroupExport.main(new String[]{"in="+root+"/input.fa","models="+root+"/models.fa","members="+root+"/members.tsv","pairs="+root+"/pairs.tsv","arm=fixture","group=group","target="+target,"cut=.70","expected=6","expectedmodels=4","selected="+selected,"verified="+pins,"outprefix="+root+'/'+name,"t=1"});}
	private static void bad(String root,String pins,String name,String target,int selected){boolean failed=false;try{run(root,pins,name,target,selected);}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+'/'+name+".fa.gz").exists(),"Invalid source evidence must fail before FASTA emission: "+name);}
	private static final String[] FILES={"input.fa","models.fa","members.tsv","pairs.tsv"};
}
