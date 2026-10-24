package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMitoCharacterization;
import prok.RrnaResourceIO;
import structures.ByteBuilder;
import structures.DoubleList;

/** Hand-calculated margins and boundary cases; no searches or generated biology.
 * @author Ganyu */
public final class NcrnaMitoCharacterizationTest {
	public static void main(String[] args) throws Exception {
		require(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require a fresh fixture directory and canonical sha80 executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);require(new File(root).mkdir() && new File(sha).canExecute(),"Isolated fixture and canonical hashing are required");
		final int[] bins=NcrnaMitoCharacterization.parseBins("200,500");
		require(NcrnaMitoCharacterization.bin(199,bins)==0 && NcrnaMitoCharacterization.bin(200,bins)==1 && NcrnaMitoCharacterization.bin(499,bins)==1 && NcrnaMitoCharacterization.bin(500,bins)==2,"Length bins are inclusive without gaps or overlap at their boundaries");
		final DoubleList values=new DoubleList(4);for(double v:new double[]{1,2,3,4}){values.add(v);}
		require(NcrnaMitoCharacterization.quantile(values,.5)==2 && NcrnaMitoCharacterization.quantile(values,.95)==4,"Nearest-rank median is lower middle for even N, not interpolated");
		boolean failed=false;try{NcrnaMitoCharacterization.parseBins("500,200");}catch(IllegalArgumentException e){failed=true;}require(failed,"Unordered bins must fail");
		final String a=row("a","m16S",1000,.70,.60,"CONCORDANT_NOMINAL_MODEL_CLASS","NA","rRNA");
		final String b=row("b","m16S",1000,.70,.69,"DISCORDANT_NOMINAL_MODEL_CLASS","NA","rRNA");
		final String c=row("c","m18S",1000,.50,.55,"DISCORDANT_NOMINAL_MODEL_CLASS","NA","rRNA");
		final String d=row("d","m18S",1000,.60,.61,"CONCORDANT_NOMINAL_MODEL_CLASS","NA","rRNA");
		final String e=row("e","UNRESOLVED_UNSUPPORTED",199,.50,.60,"NA","NO_GA_HIT_EITHER_STRAND","rRNA");
		final String f=row("f","UNRESOLVED_UNSUPPORTED",200,.40,.50,"NA","NO_GA_HIT_EITHER_STRAND","rRNA");
		final String g=row("g","UNRESOLVED_UNSUPPORTED",200,.40,.50,"NA","NO_GA_HIT_EITHER_STRAND","large subunit");
		final String h=row("h","UNRESOLVED_UNSUPPORTED",150,.40,.50,"NA","SOURCE_5S_LABEL","5S");
		final String i=row("i","LSU_ONLY",2000,.50,.60,"NA","NA","LSU");
		final String j=row("j","UNRESOLVED_TIE",950,.60,.60,"UNIVERSAL_TIE","NA","SSU");
		final String[] cm=blank();cm[0]="CM_ONLY_EXACT_GEOMETRY";cm[1]="CM0";cm[4]="AUDIT_ONLY";
		final String input=NcrnaMitoCharacterization.HEADER+'\n'+a+b+c+d+e+f+g+h+i+j+String.join("\t",cm)+'\n';
		final String[] ids={"a","b","c","d","e","f","g","h","i","j"};final int[] lengths={1000,1000,1000,1000,199,200,200,150,2000,950};
		final ByteBuilder fasta=new ByteBuilder();for(int idx=0;idx<ids.length;idx++){fasta.append('>').append(ids[idx]).nl();for(int k=0;k<lengths[idx];k++){fasta.append('A');}fasta.nl();}write(root+"/source.fa",fasta.toString());
		write(root+"/audit.tsv",input);run(root,bind(root,sha),"ok",10);
		final String[] margins=read(root+"/ok.margins.tsv").split("\n");require(margins.length==11,"Both cohorts have ALL plus four explicit concordance strata, including empty rows");
		for(String line:margins){require(line.split("\t",-1).length==19,"Every margin row, including empty strata, has the exact schema");}
		final String[] m16=find(margins,"m16S\tALL\t"),m18=find(margins,"m18S\tALL\t");
		near(m16[9],.01);near(m16[13],.055);near(m18[9],-.05);near(m18[13],-.03);
		require(m16[2].equals("2") && m16[14].equals("1") && m16[15].equals("1") && m16[16].equals("1"),"Hand-counted positive margins conserve copies and threshold inclusions");
		require(m18[14].equals("1") && m18[15].equals("1") && m18[16].equals("2"),"Negative margins use absolute thresholds without reversing the signed distribution");
		final String[] discord=find(margins,"m18S\tDISCORDANT_NOMINAL_MODEL_CLASS\t");require(discord[2].equals("1") && discord[3].equals("2"),"CM discordance stratum keeps its own numerator and full cohort denominator");
		final String[] empty=find(margins,"m18S\tSPECIALIZED_MITO_MODEL\t");require(empty[2].equals("0") && empty[4].equals("NA") && empty[14].equals("0"),"Empty strata have NA distributions, not fabricated zeros");
		final String[] nohit=read(root+"/ok.nohit.tsv").split("\n");require(nohit.length==4,"Three occupied product/bin cells retain exactly the three no-GA annotations");
		for(int n=1;n<nohit.length;n++){final String[] cells=nohit[n].split("\t",-1);require(cells.length==21 && cells[2].equals("1") && cells[3].equals("3"),"Product/length table uses the no-hit denominator, not all unsupported annotations");require(cells[6].equals("NA") && cells[7].equals("NA") && cells[8].equals("NOT_RETAINED_BY_CUT_GA"),"Missing below-GA scores must remain explicitly unavailable");}
		find(nohit,"rRNA\t1-199\t");find(nohit,"rRNA\t200-499\t");find(nohit,"large subunit\t200-499\t");
		final java.util.ArrayList<String> exported=new java.util.ArrayList<String>();RrnaResourceIO.stream(root+"/ok.nohit.fa.gz",r->{
			final String id=RrnaResourceIO.first(r.id);exported.add(id);require(tax.TaxTree.parseTaxID(r.id)==1,"Actual downstream taxonomy parser must read every exported header");
			require(r.length()==(id.equals("e")?199:200),"Third-cohort export must not trim source lengths");for(byte base:r.bases){require(base=='A',"Export must preserve exact original bases");}
			if(id.equals("e")){require(r.id.contains("partial=true pseudo=true"),"Flagged third-cohort members remain included");}
			else{require(r.id.contains("unique=u_fg"),"Two identical source copies retain their shared unique ID without collapsing");}
		});require(exported.toString().equals("[e, f, g]"),"Third-cohort membership is exactly the three no-GA IDs in audit order");
		exported.clear();RrnaResourceIO.stream(root+"/ok.ssu.fa.gz",r->{exported.add(RrnaResourceIO.first(r.id));require(tax.TaxTree.parseTaxID(r.id)==1,"Unsplit union retains parseable host taxids");for(byte base:r.bases){require(base=='A',"Unsplit union retains exact source bases");}});
		require(exported.toString().equals("[a, b, c, d, j]"),"Unsplit union contains both locked cohorts plus the preference tie, no unsupported or LSU rows");
		write(root+"/audit.tsv",input+a);bad(root,bind(root,sha),"duplicate",11);
		write(root+"/audit.tsv",input.replace(a,row("a","m16S",1000,.50,.60,"CONCORDANT_NOMINAL_MODEL_CLASS","NA","rRNA")));bad(root,bind(root,sha),"sign",10);
		final String[] wrong=e.trim().split("\t",-1);wrong[18]="RF00177";wrong[19]="49";
		write(root+"/audit.tsv",input.replace(e,String.join("\t",wrong)+'\n'));bad(root,bind(root,sha),"nohit_score",10);
		write(root+"/audit.tsv",input);bad(root,bind(root,sha),"count",11);
		write(root+"/source.fa",fasta.toString()+">e\nA\n");bad(root,bind(root,sha),"duplicate_fasta",10);
		System.out.println("MITO_CHARACTERIZATION_TEST_PASS bins=true quantiles=true signed_margins=true concordance_strata=true explicit_NA=true exact_subset=true source_copies=true taxid_headers=true flagged_retained=true unsplit_tie=true negatives=5");
	}
	private static String row(String id,String cohort,int length,double a,double b,String concordance,String reason,String product){
		final String[] r=blank();r[0]="ANNOTATION";r[1]=id;r[4]=cohort;r[5]=cohort.startsWith("m") || cohort.equals("UNRESOLVED_TIE")?"SSU_SUPPORT_REVIEW":cohort.equals("LSU_ONLY")?"LSU_SUPPORT_REVIEW":"NO_SAME_STRAND_GA_HIT_UNRESOLVED";
		r[3]="1";r[6]=""+length;r[7]=r[8]=id.equals("e")?"true":"false";r[9]=product;r[12]=id.equals("f") || id.equals("g")?"u_fg":"u_"+id;r[16]=""+a;r[17]=""+b;r[23]=concordance;r[29]=reason;return String.join("\t",r)+'\n';
	}
	private static String[] blank(){final String[] r=new String[30];Arrays.fill(r,"NA");return r;}
	private static String[] find(String[] lines,String prefix){String[] found=null;for(String line:lines){if(line.startsWith(prefix)){require(found==null,"Duplicate output group: "+prefix);found=line.split("\t",-1);}}require(found!=null,"Missing output group: "+prefix);return found;}
	private static void near(String value,double expected){require(Math.abs(Double.parseDouble(value)-expected)<1e-9,"Wrong independently calculated summary: "+value+" versus "+expected);}
	private static void run(String root,String pins,String name,int count){NcrnaMitoCharacterization.main(new String[]{"in="+root+"/audit.tsv","source="+root+"/source.fa","verified="+pins,"outprefix="+root+"/"+name,"expected="+count,"expectedcm=1","expectednohit=3","expected16=2","expected18=2","expectedssu=5","lengthbins=200,500","t=1"});}
	private static void bad(String root,String pins,String name,int count){boolean failed=false;try{run(root,pins,name,count);}catch(IllegalArgumentException e){failed=true;}require(failed && !new File(root+"/"+name+".margins.tsv").exists(),"Malformed audit cannot produce an accepted characterization: "+name);}
	private static String bind(String root,String sha) throws Exception {
		final ByteBuilder pins=new ByteBuilder();for(String file:new String[]{"audit.tsv","source.fa"}){
			final String input=root+'/'+file;final Process p=new ProcessBuilder(sha,input).redirectError(ProcessBuilder.Redirect.INHERIT).start();
			final java.io.ByteArrayOutputStream bytes=new java.io.ByteArrayOutputStream();final byte[] buffer=new byte[256];int n;
			while((n=p.getInputStream().read(buffer))>=0){bytes.write(buffer,0,n);}require(p.waitFor()==0,"Canonical fixture hashing failed");
			final String digest=new String(bytes.toByteArray(),StandardCharsets.US_ASCII).trim();require(digest.matches("[0-9a-f]{20}"),"Expected canonical sha80");pins.append(digest).tab().append(input).nl();
		}
		final String path=root+"/pins_"+(pinNumber++)+".tsv";write(path,pins.toString());return path;
	}
	private static void write(String path,String text){final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,true,false,false));out.start();out.print(text);RrnaResourceIO.finish(out);}
	private static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{require(!in.close(),"Fixture output read failed");}return b.toString();}
	private static void require(boolean value,String message){if(!value){throw new IllegalArgumentException(message);}}
	private static int pinNumber;
}
