package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import structures.ByteBuilder;

/** Converts frozen genome-wide Infernal1.1.5 tables to a model-specific truth view.
 * The launcher must verify each TABLES.tsv content hash before execution.
 * Keeps every raw hit in the ledger; complete labels require !, no truncation,
 * and model1..120, matching the Fungi112 predicate. No source GFF is modified.
 * RF00001 defaults preserve the legacy inclusion/full-model predicate. Other
 * models require an explicit GA; their family label and span remain explicit.
 * @author Ganyu, Raiden
 */
public final class RrnaCmTruthBuilder {
	public static void main(String[] args){
		String tables=null,prefix=null,model="RF00001",family="5S";int expected=-1,clen=120,minSpan=120,specified=0;double ga=Double.NEGATIVE_INFINITY;
		for(String arg:args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),value=arg.substring(eq+1);
			if(key.equals("tables")){tables=value;}
			else if(key.equals("outprefix")){prefix=value;}
			else if(key.equals("expectedgenomes")){expected=Integer.parseInt(value);}
			else if(key.equals("modelid")){model=value;}
			else if(key.equals("family")){family=value;specified|=1;}
			else if(key.equals("clen")){clen=Integer.parseInt(value);specified|=2;}
			else if(key.equals("minmodelspan")){minSpan=Integer.parseInt(value);specified|=4;}
			else if(key.equals("gabits")){ga=Double.parseDouble(value);specified|=8;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(tables==null||prefix==null||expected<1){throw new IllegalArgumentException("Required: tables= outprefix= expectedgenomes=; optional modelid= family= clen= minmodelspan= gabits=");}
		if(!model.equals("RF00001") && specified!=15){throw new IllegalArgumentException("A nonlegacy model requires explicit family, clen, minmodelspan and gabits; never inherit 5S constants");}
		run(tables,prefix,expected,new Criteria(model,family,clen,minSpan,ga));
	}

	public static void run(String tables,String prefix,int expected){
		run(tables,prefix,expected,LEGACY);
	}
	public static void run(String tables,String prefix,int expected,Criteria criteria){
		assert(tables!=null&&prefix!=null&&expected>0):"Truth construction requires a named, nonempty expected genome corpus";
		for(String suffix:new String[]{".complete.tsv",".hits.tsv.gz",".genomes.tsv"}){
			if(new File(prefix+suffix).exists()){throw new IllegalArgumentException("Output exists: "+prefix+suffix);}
		}
		final ArrayList<Genome> genomes=loadTables(tables,expected);
		final ArrayList<Hit> all=new ArrayList<Hit>(),qualified=new ArrayList<Hit>();
		final ObjectSet<String> truthKeys=new ObjectSet<String>(String.class);
		for(Genome genome:genomes){
			final ArrayList<Hit> hits=readHits(genome,criteria);
			for(Hit hit:hits){
				genome.raw++;
				if(hit.included){genome.included++;}
				if(hit.status.equals("QUALIFIED")){
					if(!truthKeys.add(hit.truthLine())){throw new IllegalArgumentException("Duplicate truth coordinates, including across genomes: "+hit.truthLine());}
					genome.qualified++;qualified.add(hit);
				}
				all.add(hit);
			}
		}
		Collections.sort(qualified);
		// All input tables and terminal markers pass before opening output writers.
		final ByteStreamWriter truth=writer(prefix+".complete.tsv"),ledger=writer(prefix+".hits.tsv.gz"),counts=writer(prefix+".genomes.tsv");
		final ByteBuilder bb=new ByteBuilder();
		try{
			truth.println("seqid\tfamily\tstart\tstop\tstrand");
			for(Hit hit:qualified){truth.println(hit.truthLine());}
			ledger.println("domain\tgenome\tseqid\ttableLine\tmodelFrom\tmodelTo\tseqFrom\tseqTo\tstrand\ttrunc\tinclusion\tscore\tevalue\tstatus\ttblout\ttblout_sha80");
			for(Hit hit:all){
				bb.clear();bb.append(hit.genome.domain).tab().append(hit.genome.name).tab().append(hit.seqid).tab().append(hit.line).tab().append(hit.mdlFrom).tab().append(hit.mdlTo);
				bb.tab().append(hit.seqFrom).tab().append(hit.seqTo).tab().append(hit.strand).tab().append(hit.trunc).tab().append(hit.included?'!':'?');
				bb.tab().append(hit.score).tab().append(hit.evalue).tab().append(hit.status).tab().append(hit.genome.table).tab().append(hit.genome.tableHash).nl();ledger.print(bb);
			}
			counts.println("domain\tgenome\trawHits\tinclusionHits\tqualifiedHits\ttblout\ttblout_sha80");
			for(Genome genome:genomes){
				bb.clear();bb.append(genome.domain).tab().append(genome.name).tab().append(genome.raw).tab().append(genome.included).tab().append(genome.qualified).tab().append(genome.table).tab().append(genome.tableHash).nl();counts.print(bb);
			}
		}finally{finish(truth,ledger,counts);}
		System.err.println("CM_TRUTH_DONE genomes="+genomes.size()+" rawHits="+all.size()+" qualified="+qualified.size()+"; model="+criteria.model+" inclusion !, trunc=no, modelSpan>="+criteria.minSpan+"/"+criteria.clen+" ga="+criteria.ga+"; no source annotation edit");
	}

	static ArrayList<Genome> loadTables(String path,int expected){
		assert(expected>0):"A missing cohort must not pass as an empty reference";
		final ArrayList<Genome> result=new ArrayList<Genome>();
		final ObjectSet<String> names=new ObjectSet<String>(String.class),paths=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		final ByteFile input=open(path);
		try{
			final byte[] header=input.nextLine();
			if(header==null||!new String(header,StandardCharsets.UTF_8).equals("domain\tgenome\tfasta\tfasta_sha80\ttblout\ttblout_sha80\torigin")){throw new IllegalArgumentException("Wrong TABLES.tsv schema");}
			for(byte[] line;(line=input.nextLine())!=null;){
				p.set(line);if(p.terms()!=7){throw new IllegalArgumentException("Expected seven manifest fields");}
				final Genome genome=new Genome(p.parseString(0),p.parseString(1),p.parseString(4),p.parseString(5));
				if(genome.domain.isEmpty()||genome.name.isEmpty()||!names.add(genome.name)||!paths.add(genome.table)){throw new IllegalArgumentException("Empty or duplicate genome/table identity: "+genome.name);}
				if(!p.parseString(3).matches("[0-9a-f]{20}")||!genome.tableHash.matches("[0-9a-f]{20}")){throw new IllegalArgumentException("Invalid declared sha80: "+genome.name);}
				if(!p.termEquals("SCANNED",6)&&!p.termEquals("REUSED_CANARY",6)&&!p.termEquals("CUT_GA_NOHMMONLY",6)){throw new IllegalArgumentException("Unknown table origin: "+genome.name);}
				result.add(genome);
			}
		}finally{if(input.close()){throw new IllegalStateException("Manifest read failed: "+path);}}
		if(result.size()!=expected){throw new IllegalArgumentException("Genome table count "+result.size()+" != expected "+expected);}
		return result;
	}

	public static ArrayList<Hit> readHits(Genome genome){
		return readHits(genome,LEGACY);
	}
	public static ArrayList<Hit> readHits(Genome genome,Criteria criteria){
		assert(genome!=null&&genome.table!=null):"Every hit must retain its genome and source-table identity";
		final ArrayList<Hit> result=new ArrayList<Hit>();
		final Fields f=new Fields();final ByteFile input=open(genome.table);
		boolean program=false,version=false,done=false;int lineNumber=0;
		try{
			for(byte[] line;(line=input.nextLine())!=null;){
				lineNumber++;
				if(line.length==0){continue;}
				if(done){throw new IllegalArgumentException("Content after CM completion marker: "+genome.table);}
				if(line[0]=='#'){
					final String text=new String(line,StandardCharsets.US_ASCII);
					if(text.startsWith("# Program:")){if(!text.substring(10).trim().equals("cmsearch")){throw new IllegalArgumentException("Not cmsearch output");}program=true;}
					if(text.startsWith("# Version:")){if(!text.substring(10).trim().startsWith("1.1.5 ")){throw new IllegalArgumentException("Not Infernal1.1.5 output");}version=true;}
					if(text.equals("# [ok]")){done=true;}
					continue;
				}
				f.set(line);
				if(!f.eq(3,criteria.model)||!f.eq(4,"cm")){throw new IllegalArgumentException("Expected "+criteria.model+" CM hit at "+genome.table+":"+lineNumber);}
				final Hit hit=new Hit(genome,lineNumber,f,criteria);
				if(hit.mdlFrom<1||hit.mdlTo<hit.mdlFrom||hit.mdlTo>criteria.clen||hit.seqFrom<1||hit.seqTo<1){throw new IllegalArgumentException("Invalid model/sequence bounds at "+genome.table+":"+lineNumber);}
				if((hit.strand!='+'&&hit.strand!='-')||(hit.strand=='+'&&hit.seqFrom>hit.seqTo)||(hit.strand=='-'&&hit.seqFrom<hit.seqTo)){throw new IllegalArgumentException("Invalid strand/coordinate order");}
				if(!f.eq(16,"!")&&!f.eq(16,"?")){throw new IllegalArgumentException("Unknown inclusion marker");}
				if(!hit.trunc.equals("no")&&!hit.trunc.equals("5'")&&!hit.trunc.equals("3'")&&!hit.trunc.equals("5'&3'")){throw new IllegalArgumentException("Unknown Infernal truncation flag: "+hit.trunc);}
				final double score=Double.parseDouble(hit.score),evalue=Double.parseDouble(hit.evalue);
				if(!Double.isFinite(score)||!Double.isFinite(evalue)||evalue<0){throw new IllegalArgumentException("Invalid CM score/E-value");}
				result.add(hit);
			}
		}finally{if(input.close()){throw new IllegalStateException("CM table read failed: "+genome.table);}}
		if(!program||!version||!done){throw new IllegalArgumentException("Missing cmsearch/version/completion footer: "+genome.table);}
		return result;
	}

	public static final class Genome {
		public Genome(String domain_,String name_,String table_,String hash_){domain=domain_;name=name_;table=table_;tableHash=hash_;}
		final String domain,name,table,tableHash;long raw,included,qualified;
	}
	public static final class Hit implements Comparable<Hit> {
		Hit(Genome genome_,int line_,Fields f){
			this(genome_,line_,f,LEGACY);
		}
		Hit(Genome genome_,int line_,Fields f,Criteria criteria){
			genome=genome_;line=line_;seqid=f.text(0);mdlFrom=f.number(5);mdlTo=f.number(6);seqFrom=f.number(7);seqTo=f.number(8);
			final String strandText=f.text(9);if(strandText.length()!=1){throw new IllegalArgumentException("Invalid strand token");}strand=strandText.charAt(0);
			trunc=f.text(10);included=f.eq(16,"!");score=f.text(14);evalue=f.text(15);
			family=criteria.family;
			status=!included?"NOT_INCLUDED":Double.parseDouble(score)<criteria.ga?"BELOW_GA":!trunc.equals("no")?"TRUNCATED":mdlTo-mdlFrom+1<criteria.minSpan?"PARTIAL_MODEL":"QUALIFIED";
		}
		public String truthLine(){return seqid+'\t'+family+'\t'+Math.min(seqFrom,seqTo)+"\t"+Math.max(seqFrom,seqTo)+"\t"+strand;}
		@Override public int compareTo(Hit other){int c=seqid.compareTo(other.seqid);if(c==0){c=Integer.compare(Math.min(seqFrom,seqTo),Math.min(other.seqFrom,other.seqTo));}if(c==0){c=Integer.compare(Math.max(seqFrom,seqTo),Math.max(other.seqFrom,other.seqTo));}return c==0?Character.compare(strand,other.strand):c;}
		final Genome genome;public final int line,mdlFrom,mdlTo,seqFrom,seqTo;public final String seqid,trunc,score,evalue,status,family;public final char strand;final boolean included;
	}
	/** Explicit model contract. Legacy truth used Infernal inclusion alone, not a second GA cutoff. */
	public static final class Criteria {
		public Criteria(String model_,String family_,int clen_,int minSpan_,double ga_){
			if(model_==null || !model_.matches("RF[0-9]{5}") || family_==null || !family_.matches("[A-Za-z0-9_]+")
				|| clen_<1 || minSpan_<1 || minSpan_>clen_ || (!Double.isFinite(ga_) && ga_!=Double.NEGATIVE_INFINITY)
				|| (!"RF00001".equals(model_) && !Double.isFinite(ga_))){
				throw new IllegalArgumentException("CM admission requires a model ID, family, valid CLEN/span and explicit finite GA for nonlegacy models");
			}
			model=model_;family=family_;clen=clen_;minSpan=minSpan_;ga=ga_;
		}
		final String model,family;final int clen,minSpan;final double ga;
	}
	static final Criteria LEGACY=new Criteria("RF00001","5S",120,120,Double.NEGATIVE_INFINITY);
	/** Infernal tblout has17 whitespace-separated fixed fields plus optional description. */
	static final class Fields {
		void set(byte[] line){
			assert(line!=null):"Tblout token bounds require an actual data row";
			bytes=line;int pos=0;
			for(int i=0;i<17;i++){
				while(pos<line.length&&line[pos]<=32){pos++;}start[i]=pos;
				while(pos<line.length&&line[pos]>32){pos++;}end[i]=pos;
				if(start[i]==end[i]){throw new IllegalArgumentException("CM row has fewer than17 fields");}
			}
		}
		int number(int field){assert(field>=0&&field<17):"Fixed tblout field range";return Parse.parseInt(bytes,start[field],end[field]);}
		String text(int field){assert(field>=0&&field<17):"Fixed tblout field range";return new String(bytes,start[field],end[field]-start[field],StandardCharsets.UTF_8);}
		boolean eq(int field,String value){assert(field>=0&&field<17):"Fixed tblout field range";if(end[field]-start[field]!=value.length()){return false;}for(int i=0;i<value.length();i++){if(bytes[start[field]+i]!=value.charAt(i)){return false;}}return true;}
		byte[] bytes;final int[] start=new int[17],end=new int[17];
	}
	private static ByteFile open(String path){return ByteFile.makeByteFile(FileFormat.testInput(path,FileFormat.TEXT,null,true,true));}
	private static ByteStreamWriter writer(String path){ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));out.start();return out;}
	private static void finish(ByteStreamWriter... writers){boolean error=false;for(ByteStreamWriter writer:writers){error|=writer.poisonAndWait();}if(error){throw new IllegalStateException("CM truth output failed");}}
}
