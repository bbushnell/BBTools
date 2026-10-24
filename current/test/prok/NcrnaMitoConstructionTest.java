package test.prok;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMitoConstruction;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Hand-specified evidence separates locked membership from descriptive ranking.
 * @author Brian, Raiden */
public final class NcrnaMitoConstructionTest {
	public static void main(String[] args) throws Exception {
		require(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical sha80 executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);require(new File(root).mkdir() && new File(sha).canExecute(),"Fixture must be isolated and use actual provenance hashes");
		require(NcrnaMitoConstruction.label5s("\"5S ribosomal RNA\"") && NcrnaMitoConstruction.label5s("5s_rRNA"),"5S source tokens tolerate case and delimiters");
		for(String s:new String[]{"5.8S ribosomal RNA","15S","25S","ribosomal RNA"}){require(!NcrnaMitoConstruction.label5s(s),"Do not mistake a different ribosomal product for 5S: "+s);}
		require(NcrnaMitoConstruction.unsupportedCategory("5S",10,18,true)==0,"5S label takes precedence over length/opposite-strand evidence");
		require(NcrnaMitoConstruction.unsupportedCategory("rRNA",10,18,true)==1,"Short descriptive category precedes the no-hit category");
		final String[] ids={"a","b","c","d","e","f","g","h","i","j","k"};
		final int[] lengths={30,30,20,12,10,25,25,22,18,25,22},starts={1,1,51,81,101,121,151,181,211,241,281};
		final String[] status={SSU,SSU,SSU,NO,NO,NO,NO,SSU,CONFLICT,LSU,SSU};
		final char[] bases={'A','A','C','G','T','A','C','G','T','G','T'};
		final double[] first={.8,.8,.6,.4,.5,.5,.6,.6,.8,.5,.6},second={.7,.7,.8,.6,.7,.7,.7,.6,.7,.6,.8};
		final ByteBuilder loci=new ByteBuilder(header(28,"locus_id")),pairs=new ByteBuilder(header(18,"locus_id")),fasta=new ByteBuilder();
		for(int i=0;i<ids.length;i++){
			final String unique=i==1?"u0":"u"+i,geometry=starts[i]+".."+(starts[i]+lengths[i]-1)+":+";
			final String[] row=blank(28);row[0]=ids[i];row[1]="A";row[2]=unique;row[8]=status[i];row[9]=ids[i];row[10]="A";row[11]="1";row[12]="mitochondrion";
			row[15]=i==3?"5S ribosomal RNA":"small subunit ribosomal RNA";row[16]="EXTRACTED";row[17]=row[18]=i==1?"true":"false";row[23]=geometry;row[24]="false";row[25]=""+lengths[i];row[27]="SEARCHED";
			line(loci,row);
			for(int m=0;m<2;m++){final String[] pair=blank(18);pair[0]=ids[i];pair[1]=unique;pair[2]=m==0?"universal":"euk18S_universal";pair[3]=""+lengths[i];pair[6]=""+(m==0?first[i]:second[i]);line(pairs,pair);}
			fasta.append('>').append(ids[i]).nl();for(int j=0;j<lengths[i];j++){fasta.append(bases[i]);}fasta.nl();
		}
		final ByteBuilder genomes=new ByteBuilder(header(9,"accession"));for(int i=0;i<3;i++){final String[] row=blank(9);row[0]=""+(char)('A'+i);row[1]=""+(i+1);line(genomes,row);}
		final ByteBuilder hits=new ByteBuilder(header(28,"run_id")),links=new ByteBuilder(header(9,"run_id"));
		final String[] models={"RF00177","RF00177","RF00177","RF01960","RF00177","RF02543","RF02543","RF00177","RF01960","RF00177","RF01960","RF01960"};
		final int[] source={0,2,6,7,8,8,9,10,10,-1,-1,-1};
		final int[] scores={80,90,60,70,75,85,100,90,90,100,80,85};
		for(int g=0;g<models.length;g++){
			final int s=source[g];final boolean genomic=s<0;final String accession=genomic?(g==11?"C":"B"):"A";
			final int start=genomic?(g==11?51:11):starts[s],length=genomic?20:lengths[s];final String geometry=start+".."+(start+length-1)+(g==2?":-":":+");
			final String[] row=hit(g,models[g],accession,geometry,genomic,scores[g],true);line(hits,row);
			if(s>=0){association(links,g,ids[s],length,g==2);if(g==0){association(links,g,"b",length,false);}}
			if(g==9){final String[] duplicate=row.clone();duplicate[2]="99";duplicate[21]="false";line(hits,duplicate);}
		}
		write(root+"/loci.tsv",loci.toString());write(root+"/pairs.tsv",pairs.toString());write(root+"/genomes.tsv",genomes.toString());
		write(root+"/hits.tsv",hits.toString());write(root+"/associations.tsv",links.toString());write(root+"/input.fa",fasta.toString());
		final String pins=bind(root,sha);run(root,pins,"ok");
		final String summary=read(root+"/ok.summary.tsv");require(summary.contains("m16S\t2\t1\t1\t1\t1\t30\t30\t18\n"),"Both exact copies survive; partial/pseudo flags do not remove the second vote; conflict-supported intact length18 defines N");
		require(summary.contains("m18S\t2\t2\t2\t0\t0\t20\t22\t18\n"),"Discordance and tied models cannot change the selected m18S members");
		final String audit=read(root+"/ok.audit.tsv.gz");
		for(String row:audit.split("\n")){require(row.split("\t",-1).length==30,"Audit records retain one schema across annotations and CM-only rows");}
		require(audit.contains("DISCORDANT_NOMINAL_MODEL_CLASS") && audit.contains("RF00177,RF01960\tTIED_BEST_MODELS"),"Audit must expose nominal discordance and full tied-model membership");
		require(audit.contains("CM_ONLY_EXACT_GEOMETRY\tCM_ONLY_0\tB\t2") && audit.contains("\t11..30:+\t2\tNA"),"Identical CM-only genomic geometry is one audit row with two model evidence members");
		require(audit.split("\n").length==14,"Eleven annotations plus two exact-geometry CM-only rows and a header are conserved");
		final String reasons=read(root+"/ok.unsupported.tsv");for(String category:new String[]{"SOURCE_5S_LABEL","SHORTER_THAN_OBSERVED_INTACT_SSU_MINIMUM","NO_GA_HIT_EITHER_STRAND","OPPOSITE_STRAND_ONLY_EVIDENCE"}){require(reasons.contains(category+"\t1\t18\tNONE"),"Four descriptive categories are exclusive with fixed precedence: "+category);}
		final ArrayList<String> names=new ArrayList<String>();RrnaResourceIO.stream(root+"/ok.m16S.fa.gz",r->{require(tax.TaxTree.parseTaxID(r.id)==1,"TrnaConsensusBuilder.writeAssignments must parse the actual exported taxid header");names.add(RrnaResourceIO.first(r.id));});require(names.toString().equals("[a, b]"),"m16S FASTA must preserve both original copies in census order");
		names.clear();RrnaResourceIO.stream(root+"/ok.m18S.fa.gz",r->{require(tax.TaxTree.parseTaxID(r.id)==1,"Exported taxid must survive the actual downstream parser");names.add(RrnaResourceIO.first(r.id));});require(names.toString().equals("[c, k]"),"m18S FASTA keeps the discordant and tied-model members");
		write(root+"/input.fa",fasta.toString()+">a\n"+repeat('A',30)+"\n");bad(root,bind(root,sha),"duplicate");
		write(root+"/input.fa",fasta.toString().replace(">b\n"+repeat('A',30),">b\n"+repeat('C',30)));bad(root,bind(root,sha),"wrong_unique");
		write(root+"/input.fa",fasta.toString());write(root+"/loci.tsv",loci.toString().replace("\ta\tA\t1\t","\ta\tA\t99\t"));bad(root,bind(root,sha),"wrong_taxid");
		write(root+"/loci.tsv",loci.toString().replace("\t"+CONFLICT+"\t","\t"+SSU+"\t"));bad(root,bind(root,sha),"wrong_status");
		System.out.println("MITO_CONSTRUCTION_TEST_PASS copies=true partial_pseudo_retained=true discordance_retained=true model_ties=true observed_N=true descriptive_precedence=true cm_only_geometry=true taxid_headers=true negatives=4");
	}
	private static String[] hit(int group,String model,String accession,String geometry,boolean genomic,int score,boolean rep){
		final String[] r=blank(28);r[0]="run"+group;r[2]="1";r[3]=model;r[4]=model.equals("RF02543")?"LSU":"SSU";r[5]=genomic?"genomic":"candidate";r[7]=accession;r[12]="1.0";r[17]=geometry;r[18]="false";r[20]=""+group;r[21]=""+rep;r[22]="no";r[24]=""+score;return r;
	}
	private static void association(ByteBuilder b,int group,String locus,int length,boolean opposite){line(b,new String[]{"run"+group,"1",""+group,"true",locus,opposite?"0":""+length,opposite?""+length:"0",""+length,""+length});}
	private static void run(String root,String pins,String output){NcrnaMitoConstruction.main(new String[]{"in="+root+"/input.fa","loci="+root+"/loci.tsv","pairs="+root+"/pairs.tsv","genomes="+root+"/genomes.tsv","hits="+root+"/hits.tsv","associations="+root+"/associations.tsv","verified="+pins,"outprefix="+root+"/"+output,"expectedloci=11","expectedgenomes=3","expectedhits=13","expectedgroups=12","expected16=2","expected18=2","expectedunsupported=4","t=1"});}
	private static void bad(String root,String pins,String output){boolean failed=false;try{run(root,pins,output);}catch(IllegalArgumentException e){failed=true;}require(failed,"Malformed provenance/content must fail: "+output);require(!new File(root+"/"+output+".summary.tsv").exists(),"Failed inputs cannot emit a successful-looking construction summary");}
	private static String bind(String root,String sha) throws IOException,InterruptedException {
		final ByteBuilder b=new ByteBuilder();for(String name:new String[]{"loci.tsv","pairs.tsv","genomes.tsv","hits.tsv","associations.tsv","input.fa"}){final String file=root+"/"+name;final Process p=new ProcessBuilder(sha,file).redirectError(ProcessBuilder.Redirect.INHERIT).start();final java.io.ByteArrayOutputStream bytes=new java.io.ByteArrayOutputStream();byte[] buffer=new byte[256];int n;while((n=p.getInputStream().read(buffer))>=0){bytes.write(buffer,0,n);}require(p.waitFor()==0,"Canonical fixture hashing failed");final String digest=new String(bytes.toByteArray(),StandardCharsets.US_ASCII).trim();require(digest.matches("[0-9a-f]{20}"),"Canonical hash must be sha80");b.append(digest).tab().append(file).nl();}
		final String path=root+"/pins_"+(pinNumber++)+".tsv";write(path,b.toString());return path;
	}
	private static String header(int n,String first){final String[] cells=blank(n);cells[0]=first;return String.join("\t",cells)+"\n";}
	private static String[] blank(int n){final String[] a=new String[n];java.util.Arrays.fill(a,"NA");return a;}
	private static void line(ByteBuilder b,String[] cells){b.append(String.join("\t",cells)).nl();}
	private static String repeat(char c,int count){final char[] a=new char[count];java.util.Arrays.fill(a,c);return new String(a);}
	private static void write(String path,String text){final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,true,false,false));out.start();out.print(text);RrnaResourceIO.finish(out);}
	private static String read(String path){final ByteBuilder b=new ByteBuilder();final ByteFile in=RrnaResourceIO.open(path);try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{require(!in.close(),"Fixture output read failed");}return b.toString();}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static int pinNumber;
	private static final String SSU="SSU_SUPPORT_REVIEW",LSU="LSU_SUPPORT_REVIEW",CONFLICT="SSU_LSU_CONFLICT_REVIEW",NO="NO_SAME_STRAND_GA_HIT_UNRESOLVED";
}
