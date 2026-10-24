package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Exact nuclear 5S/SSU membership join to attributed source GFFs.
 * Assembly identity comes from the declared source GFF, never from TaxID.
 * @author Raiden
 */
public final class RrnaSourceTargetPrep {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,RrnaSourceTargetPrep.class,false);args=pp.args;
		String clade=null,five=null,ssu=null,gff=null,sources=null,genomes=null,bulk=null,prefix=null;
		for(String arg:args){
			final int eq=arg.indexOf('=');if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("clade")){clade=v;}else if(k.equals("five")){five=v;}else if(k.equals("ssu")){ssu=v;}
			else if(k.equals("gff")){gff=v;}else if(k.equals("sources")){sources=v;}else if(k.equals("genomes")){genomes=v;}
			else if(k.equals("bulk")){bulk=v;}else if(k.equals("outprefix")){prefix=v;}else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		if(clade==null||five==null||ssu==null||gff==null||sources==null||prefix==null||(genomes==null)==(bulk==null)){
			throw new IllegalArgumentException("Required: clade= five= ssu= gff= sources= outprefix= and exactly one of genomes=<paired path/taxid TSV> or bulk=<genome stream path>");
		}
		new RrnaSourceTargetPrep(clade,sources,genomes,bulk).run(five,ssu,gff,prefix);
	}
	public RrnaSourceTargetPrep(String clade_,String sourceList,String genomeList,String bulk){
		if(!clade_.matches("[a-z_]+")){throw new IllegalArgumentException("Clade must be a stable lowercase identifier");}clade=clade_;
		final ObjectMap<String,Source> genomeByStem=new ObjectMap<String,Source>(String.class,Source.class);
		if(genomeList!=null){
			final ByteFile f=open(genomeList);final LineParser1 p=new LineParser1('\t');
			try{for(byte[] line;(line=f.nextLine())!=null;){
				p.set(line);require(p.terms()==2,"Genome inventory needs path and taxid");
				final String path=p.parseString(0);require(path.startsWith("/")&&path.endsWith(".fna.gz"),"Genome inventory requires absolute .fna.gz paths");
				final Source s=new Source();s.fasta=path;s.taxid=p.parseInt(1);require(s.taxid>0,"Genome taxonomy must be explicit");
				require(genomeByStem.put(stem(path,".fna.gz"),s)==null,"Duplicate paired genome stem");
			}}finally{close(f);}
		}
		final ByteFile f=open(sourceList);
		try{for(byte[] line;(line=f.nextLine())!=null;){
			final String path=new String(line,StandardCharsets.UTF_8);require(path.startsWith("/")&&path.endsWith(".gff.gz"),"Source list requires absolute .gff.gz paths");
			final String id=stem(path,".gff.gz");final Source s=genomeList==null?new Source():genomeByStem.get(id);
			require(s!=null,"No exactly paired genome file for source GFF: "+path);
			s.gff=path;s.id=id;if(genomeList==null){s.fasta=bulk;}
			require(sources.put(new File(path).getName(),s)==null,"Source basename is ambiguous: "+path);
		}}finally{close(f);}
		require(sources.size()>0,"Empty source list");
	}
	public void run(String five,String ssu,String gff,String prefix){
		for(String suffix:new String[]{".targets.gff.gz",".targets.tsv.gz",".counts.tsv",".assemblies.tsv"}){require(!new File(prefix+suffix).exists(),"Output exists: "+prefix+suffix);}
		Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		loadMembers(five,0);loadMembers(ssu,1);
		require(!ordered.isEmpty(),"No input members");
		final ByteStreamWriter gout=writer(prefix+".targets.gff.gz"),ledger=writer(prefix+".targets.tsv.gz");
		final ByteFile input=open(gff);final LineParser1 p=new LineParser1('\t');final Attrs attrs=new Attrs();final ByteBuilder b=new ByteBuilder(),key=new ByteBuilder();
		int matched=0;
		try{
			ledger.println("target_id\tclade\tfamily\tassembly_id\tsource_gff\tgenome_input\ttaxid\tfeature_id\tseqid\tstart\tstop\tstrand\tlength\tsource_product\tmember_input");
			for(byte[] line;(line=input.nextLine())!=null;){
				if(line.length==0||line[0]=='#'){continue;}p.set(line);require(p.terms()==9,"GFF must have nine fields");if(!p.termEquals("rRNA",2)){continue;}
				p.setBounds(8);attrs.set(line,p.a(),p.b());
				// Other rRNA families remain outside this exact-member join.
				final Member m=members.get(attrs.key(key));if(m==null){continue;}
				require(!m.found,"Duplicate source GFF row for selected member: "+m.id);
				require(p.termEquals(m.contig,0)&&p.parseInt(3)==m.start&&p.parseInt(4)==m.stop&&p.termEquals(m.strand,6),"Source GFF coordinates disagree with attributed member: "+m.id);
				require(attrs.eq(6,"nuclear")&&attrs.eq(7,clade)&&attrs.eq(8,m.product),"Source compartment/clade/product disagreement: "+m.id);
				require(!attrs.hasReserved(),"Source attributes already contain task27_* output fields");
				b.clear().append(line).append(";task27_target=").append(m.id).append(";task27_family=").append(FAMILIES[m.family]).append(";task27_assembly=").append(m.source.id).append(";task27_taxid=").append(m.taxid).nl();gout.print(b);
				b.clear().append(m.id).tab().append(clade).tab().append(FAMILIES[m.family]).tab().append(m.source.id).tab().append(m.source.gff).tab().append(m.source.fasta).tab().append(m.taxid).tab().append(m.feature).tab().append(m.contig).tab().append(m.start).tab().append(m.stop).tab().append(m.strand).tab().append(m.length).tab().append(m.product).tab().append(m.input).nl();ledger.print(b);
				m.found=true;matched++;
			}
		}finally{
			final boolean inputError=input.close();boolean outputError=gout.poisonAndWait();outputError|=ledger.poisonAndWait();require(!inputError&&!outputError,"Target I/O failed");
		}
		for(Member m:ordered){require(m.found,"Selected member absent from source GFF: "+m.id+" "+m.source.id+" "+m.feature);}
		require(matched==ordered.size(),"Exact membership join did not conserve inputs");
		final ByteStreamWriter counts=writer(prefix+".counts.tsv"),assemblies=writer(prefix+".assemblies.tsv");
		try{
			counts.println("clade\tfamily\tloci\tassemblies");assemblies.println("clade\tassembly_id\tsource_gff\tgenome_input\t5S_loci\tSSU_loci");
			for(int family=0;family<2;family++){b.clear().append(clade).tab().append(FAMILIES[family]).tab().append(locusCounts[family]).tab().append(assemblySets.get(family).size()).nl();counts.print(b);}
			for(Source s:usedSources){b.clear().append(clade).tab().append(s.id).tab().append(s.gff).tab().append(s.fasta).tab().append(s.counts[0]).tab().append(s.counts[1]).nl();assemblies.print(b);}
		}finally{boolean error=counts.poisonAndWait();error|=assemblies.poisonAndWait();require(!error,"Summary output failed");}
		System.err.println("RRNA_SOURCE_TARGET_PREP_PASS clade="+clade+" members="+matched+" five="+locusCounts[0]+" ssu="+locusCounts[1]+"; no genome scan or biological label change");
	}
	void loadMembers(String path,int family){
		assert(family>=0&&family<2):"Only the two declared membership sets are accepted";
		final Streamer st=StreamerFactory.makeStreamer(FileFormat.testInput(path,FileFormat.FASTA,null,true,true),null,true,-1,false,true,1);
		final Attrs a=new Attrs();final ByteBuilder key=new ByteBuilder();st.start();boolean completed=false;
		try{
			for(ListNum<Read> ln;(ln=st.nextList())!=null;){for(Read read:ln.list){
				final String header=read.id;require(header!=null&&header.startsWith("tid|"),"Member lacks explicit taxid prefix");final int pipe=header.indexOf('|',4);require(pipe>4,"Malformed taxid prefix");
				final int taxid=Integer.parseInt(header.substring(4,pipe));require(taxid>0,"Invalid member taxid");
				final byte[] bytes=header.substring(pipe+1).getBytes(StandardCharsets.UTF_8);a.set(bytes,0,bytes.length);
				final String sourceToken=a.text(1);final Source source=sources.get(new File(sourceToken).getName());
				require(source!=null&&(sourceToken.equals(new File(source.gff).getName())||sourceToken.equals(source.gff)),"Member source GFF absent or ambiguous: "+sourceToken);
				require(source.taxid<1||source.taxid==taxid,"Paired genome taxid conflicts with member metadata");
				require(a.eq(6,"nuclear")&&a.eq(7,clade),"Member has wrong compartment or clade");
				final Member m=new Member();m.source=source;m.taxid=taxid;m.family=family;m.input=path;m.feature=a.text(0);m.contig=a.text(2);m.start=a.number(3);m.stop=a.number(4);m.strand=a.text(5);m.product=a.text(8);
				m.length=read.length();require(m.start>0&&m.stop>=m.start&&m.length==m.stop-m.start+1&&(m.strand.equals("+")||m.strand.equals("-")),"Member sequence length/coordinates/strand invalid: "+m.feature);
				final String prior=contigAssembly.put(m.contig,source.id);require(prior==null||prior.equals(source.id),"Contig accession occurs in multiple assemblies; streaming source identity is unresolved: "+m.contig);
				m.id="task27_"+clade+"_"+ordered.size();require(members.put(a.key(key),m)==null,"Duplicate or cross-family member identity: "+m.feature);ordered.add(m);
				if(source.counts[0]+source.counts[1]==0){usedSources.add(source);}source.counts[family]++;locusCounts[family]++;assemblySets.get(family).add(source.id);
			}st.returnList(ln);}
			completed=true;
		}finally{if(!completed){st.close();}require(!st.errorState(),"Member FASTA read failed: "+path);}
	}
	static final class Attrs{
		void set(byte[] line,int a,int b){
			assert(a>=0&&b>=a&&b<=line.length):"Attribute slice must stay inside its input record";bytes=line;from=a;to=b;Arrays.fill(starts,-1);Arrays.fill(ends,-1);
			for(int pos=a;pos<b;){int end=pos;while(end<b&&line[end]!=';'){end++;}int eq=pos;while(eq<end&&line[eq]!='='){eq++;}
				if(eq<end){for(int k=0;k<KEYS.length;k++){if(eq-pos==KEYS[k].length()&&region(line,pos,eq,KEYS[k])){require(starts[k]<0,"Duplicate source attribute: "+KEYS[k]);starts[k]=eq+1;ends[k]=end;break;}}}pos=end+1;
			}
		}
		String key(ByteBuilder b){b.clear();for(int k:IDENTITY){require(starts[k]>=0&&ends[k]>starts[k],"Missing identity attribute: "+KEYS[k]);if(b.length()>0){b.tab();}b.append(bytes,starts[k],ends[k]-starts[k]);}return b.toString();}
		String text(int k){require(starts[k]>=0&&ends[k]>starts[k],"Missing source attribute: "+KEYS[k]);return new String(bytes,starts[k],ends[k]-starts[k],StandardCharsets.UTF_8);}
		int number(int k){require(starts[k]>=0,"Missing numeric source attribute");return Parse.parseInt(bytes,starts[k],ends[k]);}
		boolean eq(int k,String value){return starts[k]>=0&&region(bytes,starts[k],ends[k],value);}
		boolean hasReserved(){for(int p=from;p<to;p++){if((p==from||bytes[p-1]==';')&&p+7<=to&&region(bytes,p,p+7,"task27_")){return true;}}return false;}
		final int[] starts=new int[KEYS.length],ends=new int[KEYS.length];byte[] bytes;int from,to;
	}
	static boolean region(byte[] b,int a,int end,String s){if(end-a!=s.length()){return false;}for(int i=0;i<s.length();i++){if(b[a+i]!=(byte)s.charAt(i)){return false;}}return true;}
	static String stem(String path,String suffix){final String name=new File(path).getName();require(name.endsWith(suffix)&&name.length()>suffix.length(),"Invalid source filename");return name.substring(0,name.length()-suffix.length());}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static ByteFile open(String path){return ByteFile.makeByteFile(FileFormat.testInput(path,FileFormat.TEXT,null,true,true));}
	static void close(ByteFile f){require(!f.close(),"Source input read failed");}
	public static ByteStreamWriter writer(String path){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));w.start();return w;}
	static final class Source{String id,gff,fasta;int taxid;final long[] counts=new long[2];}
	static final class Member{String id,feature,contig,strand,product,input;Source source;int family,taxid,start,stop,length;boolean found;}
	final String clade;final ObjectMap<String,Source> sources=new ObjectMap<String,Source>(String.class,Source.class);
	final ObjectMap<String,Member> members=new ObjectMap<String,Member>(String.class,Member.class);
	final ObjectMap<String,String> contigAssembly=new ObjectMap<String,String>(String.class,String.class);
	final ArrayList<Member> ordered=new ArrayList<Member>();final ArrayList<Source> usedSources=new ArrayList<Source>();
	final ArrayList<ObjectSet<String>> assemblySets=new ArrayList<ObjectSet<String>>(Arrays.asList(new ObjectSet<String>(String.class),new ObjectSet<String>(String.class)));
	final long[] locusCounts=new long[2];
	static final String[] FAMILIES={"5S","SSU"};
	static final String[] KEYS={"ID","source_gff","source_contig","source_start","source_end","source_strand","source_compartment","source_group","product"};
	static final int[] IDENTITY={1,2,3,4,5,0};
}
