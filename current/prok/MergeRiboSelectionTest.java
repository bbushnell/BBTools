package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import aligner.SingleStateAlignerFlat2;
import consensus.BaseGraph;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import tax.TaxTree;
import tax.GiToTaxid;

/** Black-box fixtures for ranked member selection and duplicate consensus votes.
 * @author Raiden
 */
public final class MergeRiboSelectionTest {
	public static void main(String[] args){
		if(args.length==3 && args[0].equals("compare")){compare(args[1],args[2]);return;}
		check(args.length==2,"Expected create|check and fixture/output directory, or compare baseline output");
		if(args[0].equals("create")){create(args[1]);}else{check(args[0].equals("check"),"Unknown test action");verify(args[1]);}
	}
	static byte[][] sequences(){
		final byte[] b=new byte[120],alphabet={'A','C','G','T'};final Random random=new Random(9272026);
		for(int i=0;i<b.length;i++){b[i]=alphabet[random.nextInt(4)];}
		final byte[] a=b.clone(),c=b.clone(),d=b.clone();
		for(int i=10;i<=100;i+=10){a[i]=different(a[i]);}c[110]=different(c[110]);d[111]=different(d[111]);
		return new byte[][]{a,b,c,d};
	}
	static byte different(byte b){return b=='A'?(byte)'C':(byte)'A';}
	static void create(String dir){
		check(new File(dir).mkdirs(),"Fixture directory must be fresh");final byte[][] seq=sequences();
		final ByteStreamWriter input=writer(dir+"/input.fa"),unique=writer(dir+"/unique.fa"),ref=writer(dir+"/ref.fa");
		try{
			for(int i=0;i<7;i++){fasta(input,"tid|11|A"+i,seq[0]);}
			for(int i=0;i<2;i++){fasta(input,"tid|11|B"+i,seq[1]);}
			fasta(input,"tid|11|C",seq[2]);fasta(input,"tid|11|D",seq[3]);
			for(int i=0;i<4;i++){fasta(unique,"tid|11|unique"+i,seq[i]);}
			fasta(input,"tid|12|A0",seq[0]);fasta(input,"tid|12|A1",seq[0]);
			fasta(input,"tid|13|B",seq[1]);
			fasta(input,"tid|14|C",seq[2]);fasta(input,"tid|14|B",seq[1]);
			fasta(input,"tid|15|D0",seq[3]);fasta(input,"tid|15|D1",seq[3]);
			fasta(ref,"global_seed",seq[1]);
		}finally{finish(input,unique,ref);}
		final ByteStreamWriter names=writer(dir+"/names.dmp"),nodes=writer(dir+"/nodes.dmp");
		try{
			names.println("1\t|\troot\t|\t\t|\tscientific name\t|");nodes.println("1\t|\t1\t|\tno rank\t|\t0");
			names.println("2\t|\tBacteria\t|\t\t|\tscientific name\t|");nodes.println("2\t|\t1\t|\tsuperkingdom\t|\t0");
			names.println("2157\t|\tArchaea\t|\t\t|\tscientific name\t|");nodes.println("2157\t|\t1\t|\tsuperkingdom\t|\t0");
			for(int i=11;i<=14;i++){names.println(i+"\t|\tFixture species "+i+"\t|\t\t|\tscientific name\t|");nodes.println(i+"\t|\t1\t|\tspecies\t|\t0");}
		}finally{finish(names,nodes);}
		final TaxTree tree=TaxTree.loadTaxTree(null,dir+"/names.dmp",dir+"/nodes.dmp",null,System.err,false,false);
		TaxTree.writeTaxTree(tree,dir+"/fixture.taxtree.gz",false);
		System.out.println("MERGERIBO_FIXTURES_CREATED members=18 taxids=5 species11_duplicates_vote=7A_2B_1C_1D");
	}
	static void verify(String dir){
		final byte[][] s=sequences();
		final ArrayList<Read> distinct=load(dir+"/distinct.fa"),repeated=load(dir+"/repeated.fa"),all=load(dir+"/all_distinct.fa"),fast=load(dir+"/fast.fa");
		check(distinct.size()==7&&repeated.size()==9&&all.size()==9,"N=2 distinct/multiplicity and N=10 distinct cardinalities across all five taxa");
		group(distinct,11,s[0],s[1]);group(repeated,11,s[0],s[0]);rankedDistinct(all,11,s);
		group(distinct,12,s[0]);group(repeated,12,s[0],s[0]);group(distinct,13,s[1]);group(distinct,14,s[1],s[2]);group(distinct,15,s[3]);
		group(fast,11,s[1],s[3]);group(fast,12,s[0]);group(fast,13,s[1]);group(fast,14,s[1],s[2]);
		group(load(dir+"/single.fa"),11,s[0]);group(load(dir+"/single_dedupe.fa"),11,s[0]);
		group(load(dir+"/weighted_consensus.fa"),11,s[0]);group(load(dir+"/unique_consensus.fa"),11,s[1]);
		check(!Arrays.equals(s[0],s[1]),"The duplicates-vote fixture must actually change its consensus");
		final ArrayList<Read> dada=load(dir+"/dada2.fa");check(dada.size()==6,"DADA2 outputs every selected member and omits the undefined TaxID15 lineage");
		int named11=0;for(Read r:dada){check(r.id.startsWith("k__NA;")&&r.id.contains("s__Fixture species ")&&!r.id.contains("tid|"),"Every selected DADA2 header is converted");if(r.id.endsWith("species 11")){named11++;}}
		check(named11==2,"DADA2 conversion covers both representatives of one species");
		seedAccounting(s[1]);
		System.out.println("MERGERIBO_SELECTION_TEST_PASS distinct7 repeated9 allDistinct9 fast_order small_single_two duplicate_votes_consensus_A_vs_B dada2_6");
	}
	/** Test the graph API contract: a seed prior is not an already-added member. */
	static void seedAccounting(byte[] seed){
		final BaseGraph graph=new BaseGraph("seed",seed,null,0,10);
		check(graph.readTotal==0 && graph.ref[10].weightSum==1 && graph.ref[10].countSum==1,"Constructor installs one weight-1 prior, zero real reads");
		final SingleStateAlignerFlat2 aligner=new SingleStateAlignerFlat2();
		for(int i=0;i<3;i++){
			final Read r=new Read(seed.clone(),null,"copy"+i,i);
			graph.alignAndGenerateMatch(r,aligner);graph.add(r);
		}
		check(graph.readTotal==3 && graph.ref[10].countSum==4 && graph.ref[10].weightSum==1+3*(BaseGraph.fakeQuality+1),"Three real copies must contribute three full weights, in addition to the reference prior");
		System.out.println("SEED_ACCOUNTING_PASS real_reads=3 node_prior_count=1 node_count=4 node_weight="+graph.ref[10].weightSum);
	}
	/** Assert the requested ranking without imposing an unspecified order on ties. */
	static void rankedDistinct(ArrayList<Read> reads,int taxid,byte[][] expected){
		final boolean[] seen=new boolean[expected.length];final SingleStateAlignerFlat2 aligner=new SingleStateAlignerFlat2();
		float previous=Float.POSITIVE_INFINITY;int n=0;
		for(Read r:reads){if(r.id.startsWith("tid|"+taxid+"|")){
			int which=-1;for(int i=0;i<expected.length;i++){if(Arrays.equals(r.bases,expected[i])){which=i;break;}}
			check(which>=0 && !seen[which],"All distinct output must contain each expected sequence once");seen[which]=true;
			final float score=aligner.align(r.bases,expected[0]);
			check(score<=previous,"Equal-length members must rank by identity to the known majority consensus");previous=score;n++;
			System.out.println("RANK_CHECK taxid="+taxid+" sequence="+which+" identity="+score);
		}}
		check(n==expected.length,"All distinct members must survive N above the distinct count");
	}
	static void compare(String baseline,String output){
		for(String arm:new String[]{"best","fast","consensus"}){
			final ArrayList<Read> before=load(baseline+"/"+arm+".fa"),after=load(output+"/baseline_"+arm+".fa");
			check(before.size()==after.size(),"Voting fix must retain one output per accepted taxon");int basesChanged=0,headersChanged=0;
			for(Read a:before){
				final Integer tid=GiToTaxid.parseTaxidNumber(a.id,'|');Read match=null;
				for(Read b:after){if(tid.equals(GiToTaxid.parseTaxidNumber(b.id,'|'))){check(match==null,"Duplicate output taxon");match=b;}}
				check(match!=null,"Missing output taxon after voting fix: "+tid);
				if(!Arrays.equals(a.bases,match.bases)){basesChanged++;}if(!a.id.equals(match.id)){headersChanged++;}
			}
			System.out.println("REGRESSION_COMPARE mode="+arm+" taxa="+before.size()+" changed_sequences="+basesChanged+" changed_headers="+headersChanged);
		}
	}
	static void group(ArrayList<Read> reads,int taxid,byte[]... expected){
		int n=0;final String prefix="tid|"+taxid+"|";
		for(Read r:reads){if(r.id.startsWith(prefix)){check(n<expected.length&&Arrays.equals(r.bases,expected[n]),"Wrong sequence or ranked order for taxid="+taxid+" rank="+(n+1)+" header="+r.id);n++;}}
		check(n==expected.length,"Wrong output count for taxid="+taxid+": "+n+" != "+expected.length);
	}
	static ArrayList<Read> load(String file){
		final ArrayList<Read> result=new ArrayList<Read>();final Streamer in=StreamerFactory.makeStreamer(FileFormat.testInput(file,FileFormat.FASTA,null,true,true),null,true,-1,false,true,1);in.start();boolean done=false;
		try{for(ListNum<Read> ln;(ln=in.nextList())!=null;){result.addAll(ln.list);in.returnList(ln);}check(!in.errorState(),"Output read failed");done=true;}finally{if(!done){in.close();}}
		return result;
	}
	static ByteStreamWriter writer(String path){final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));out.start();return out;}
	static void fasta(ByteStreamWriter out,String id,byte[] seq){out.print(new ByteBuilder().append('>').append(id).nl().append(seq).nl());}
	static void finish(ByteStreamWriter... writers){for(ByteStreamWriter w:writers){check(!w.poisonAndWait(),"Fixture write failed");}}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
