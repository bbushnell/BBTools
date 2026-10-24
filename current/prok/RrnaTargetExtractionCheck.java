package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
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
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Verify an exact-member extraction and emit unique source-attributed FASTA headers.
 * No endpoint correction, locus filtering, or sequence deduplication occurs.
 * @author Raiden
 */
public final class RrnaTargetExtractionCheck {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,RrnaTargetExtractionCheck.class,false);args=pp.args;
		String targets=null,five=null,ssu=null,in=null,out=null,summary=null;int flank=-1;
		for(String arg:args){int eq=arg.indexOf('=');require(eq>0,"Expected key=value");String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("targets")){targets=v;}else if(k.equals("five")){five=v;}else if(k.equals("ssu")){ssu=v;}else if(k.equals("in")){in=v;}else if(k.equals("out")){out=v;}else if(k.equals("summary")){summary=v;}else if(k.equals("flank")){flank=Parse.parseIntKMG(v);}else{throw new IllegalArgumentException("Unknown option: "+k);}}
		require(targets!=null&&five!=null&&ssu!=null&&in!=null&&out!=null&&summary!=null&&flank>=0,"Required targets= five= ssu= in= out= summary= flank=");
		new RrnaTargetExtractionCheck().run(targets,five,ssu,in,out,summary,flank);
	}
	public void run(String targets,String five,String ssu,String in,String out,String summary,int flank){
		require(!new File(out).exists()&&!new File(summary).exists(),"Refusing existing extraction-check output");Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		loadTargets(targets);loadExact(five,"5S");loadExact(ssu,"SSU");
		for(Target t:byId.values()){if(t!=null){require(t.exact!=null,"Target has no exact member: "+t.id);}}
		final Streamer st=stream(in);final Writer writer=WriterFactory.makeWriter(FileFormat.testOutput(out,FileFormat.FASTA,null,true,false,false,false),null,1,null,false);
		final ByteBuilder header=new ByteBuilder();int records=0;boolean completed=false;writer.start();st.start();
		try{
			for(ListNum<Read> ln;(ln=st.nextList())!=null;){for(Read r:ln.list){
				final String id=field(r.id,"task27_target");final Target t=byId.get(id);require(t!=null&&!t.seen,"Foreign or duplicate extracted target: "+id);
				require(field(r.id,"task27_family").equals(t.family)&&field(r.id,"task27_assembly").equals(t.assembly)&&Integer.parseInt(field(r.id,"task27_taxid"))==t.taxid,"Extraction lost family/assembly/taxid attribution: "+id);
				require(new File(field(r.id,"source_gff")).getName().equals(new File(t.source).getName())&&field(r.id,"source_contig").equals(t.seq)&&Integer.parseInt(field(r.id,"source_start"))==t.start&&Integer.parseInt(field(r.id,"source_end"))==t.stop&&field(r.id,"source_strand").equals(t.strand),"Extraction source coordinates changed: "+id);
				final int left=Integer.parseInt(field(r.id,"lflank")),right=Integer.parseInt(field(r.id,"rflank"));final String gcText=field(r.id,"contig_gc");final double gc=Double.parseDouble(gcText);
				require(left>=0&&right>=0&&left<=flank&&right<=flank&&Double.isFinite(gc)&&gc>=0&&gc<=1,"Invalid actual flanks or GC: "+id);
				require((t.strand.equals("+")?left:right)==Math.min(flank,t.start-1),"Genomic-left clipping disagrees with source coordinates: "+id);
				require(r.bases!=null&&r.length()==t.exact.length+left+right,"Extraction length does not match exact member plus actual flanks: "+id);
				for(int i=0;i<r.bases.length;i++){r.bases[i]=upper(r.bases[i]);}
				for(int i=0;i<t.exact.length;i++){if(r.bases[left+i]!=t.exact[i]){throw new IllegalArgumentException("Exact/de-flanked sequence disagreement: "+id+" base "+i);}}
				header.clear().append("tid|").append(t.taxid).append('|').append(t.id).append(" family=").append(t.family).append(" assembly_id=").append(t.assembly).append(" source_gff=").append(t.source).append(" source_contig=").append(t.seq).append(" source_start=").append(t.start).append(" source_end=").append(t.stop).append(" source_strand=").append(t.strand).append(" lflank=").append(left).append(" rflank=").append(right).append(" contig_gc=").append(gcText);
				r.id=header.toString();t.seen=true;records++;
			}writer.addReads(ln);}
			for(Target t:byId.values()){if(t!=null){require(t.seen,"Target missing from extraction: "+t.id);}}
			require(records==byId.size()&&!st.errorState(),"Extracted membership/input stream did not conserve targets");completed=true;
		}finally{
			if(!completed){st.close();}
			// One producer submits contiguous complete batches; there is no missing writer job.
			// Drain on validation failure too, so compression streams receive their terminal.
			//TODO: Probable bug - FastqWriter.finishError leaves BGZF threads alive on external
			// failure (Task27 fixture thread dump, 2026-09-25); avoid that path here. Core fix separate.
			require(!writer.poisonAndWait(),"Validated FASTA output failed");
		}
		final ObjectSet<String> fiveAssemblies=new ObjectSet<String>(String.class),ssuAssemblies=new ObjectSet<String>(String.class);int n5=0,ns=0;
		for(Target t:byId.values()){if(t==null){continue;}if(t.family.equals("5S")){n5++;fiveAssemblies.add(t.assembly);}else{ns++;ssuAssemblies.add(t.assembly);}}
		final ByteStreamWriter count=new ByteStreamWriter(FileFormat.testOutput(summary,FileFormat.TEXT,null,true,false,false,false));count.start();
		try{count.println("family\tloci\tassemblies\tdeFlankParity");header.clear().append("5S\t").append(n5).tab().append(fiveAssemblies.size()).append("\tPASS\nSSU\t").append(ns).tab().append(ssuAssemblies.size()).append("\tPASS\n");count.print(header);}finally{require(!count.poisonAndWait(),"Parity summary output failed");}
		System.err.println("RRNA_TARGET_EXTRACTION_CHECK_PASS targets="+records+" five="+n5+" ssu="+ns+"; exact membership and de-flanked sequence parity; no label change");
	}
	public void loadTargets(String path){
		final ByteFile f=ByteFile.makeByteFile(FileFormat.testInput(path,FileFormat.TEXT,null,true,true));final LineParser1 p=new LineParser1('\t');
		try{final byte[] h=f.nextLine();require(h!=null&&new String(h,StandardCharsets.UTF_8).equals("target_id\tclade\tfamily\tassembly_id\tsource_gff\tgenome_input\ttaxid\tfeature_id\tseqid\tstart\tstop\tstrand\tlength\tsource_product\tmember_input"),"Unexpected target-ledger header");
			for(byte[] line;(line=f.nextLine())!=null;){p.set(line);require(p.terms()==15,"Target schema changed");Target t=new Target();t.id=p.parseString(0);t.family=p.parseString(2);t.assembly=p.parseString(3);t.source=p.parseString(4);t.taxid=p.parseInt(6);t.feature=p.parseString(7);t.seq=p.parseString(8);t.start=p.parseInt(9);t.stop=p.parseInt(10);t.strand=p.parseString(11);
				require(!t.id.isEmpty()&&!t.assembly.isEmpty()&&t.source.startsWith("/")&&t.taxid>0&&(t.family.equals("5S")||t.family.equals("SSU"))&&t.start>0&&t.stop>=t.start&&p.parseInt(12)==t.stop-t.start+1&&(t.strand.equals("+")||t.strand.equals("-")),"Malformed attributed target");
				require(byId.put(t.id,t)==null&&byKey.put(key(t.source,t.seq,t.start,t.stop,t.strand,t.feature),t)==null,"Duplicate target identity");
			}
		}finally{require(!f.close(),"Target ledger read failed");}require(byId.size()>0,"Empty targets");
	}
	public void loadExact(String path,String family){
		final Streamer st=stream(path);final RrnaSourceTargetPrep.Attrs a=new RrnaSourceTargetPrep.Attrs();st.start();boolean done=false;
		try{for(ListNum<Read> ln;(ln=st.nextList())!=null;){for(Read r:ln.list){
			final int pipe=r.id.indexOf('|',4);require(r.id.startsWith("tid|")&&pipe>4,"Exact member needs explicit taxonomy");final int taxid=Integer.parseInt(r.id.substring(4,pipe));byte[] attr=r.id.substring(pipe+1).getBytes(StandardCharsets.UTF_8);a.set(attr,0,attr.length);
			final Target t=byKey.get(key(a.text(1),a.text(2),a.number(3),a.number(4),a.text(5),a.text(0)));
			require(t!=null&&t.exact==null&&t.family.equals(family)&&t.taxid==taxid&&r.bases!=null&&r.length()==t.stop-t.start+1,"Exact member is foreign/duplicate or has changed metadata/length");t.exact=r.bases;for(int i=0;i<t.exact.length;i++){t.exact[i]=upper(t.exact[i]);}
		}st.returnList(ln);}require(!st.errorState(),"Exact member input failed");done=true;}finally{if(!done){st.close();}}
	}
	static String field(String header,String key){
		assert(header!=null&&key!=null):"Scalar metadata requires an intact FASTA header";
		final String prefix=key+"=";int found=-1,end=-1;
		for(int pos=header.indexOf(prefix);pos>=0;pos=header.indexOf(prefix,pos+prefix.length())){if(pos>0&&header.charAt(pos-1)!=';'&&header.charAt(pos-1)!=' '){continue;}require(found<0,"Duplicate metadata: "+key);found=pos+prefix.length();end=found;while(end<header.length()&&header.charAt(end)!=';'&&!Character.isWhitespace(header.charAt(end))){end++;}}
		require(found>=0&&end>found,"Missing metadata: "+key);return header.substring(found,end);
	}
	static byte upper(byte b){return b>='a'&&b<='z'?(byte)(b-32):b;}
	static String key(String source,String seq,int start,int stop,String strand,String id){return new File(source).getName()+"\t"+seq+"\t"+start+"\t"+stop+"\t"+strand+"\t"+id;}
	static Streamer stream(String path){return StreamerFactory.makeStreamer(FileFormat.testInput(path,FileFormat.FASTA,null,true,true),null,true,-1,false,true,1);}
	static void require(boolean condition,String why){if(!condition){throw new IllegalArgumentException(why);}}
	public static class Target{String id,family,assembly,source,feature,seq,strand;int taxid,start,stop;public byte[] exact;boolean seen;}
	public final ObjectMap<String,Target> byId=new ObjectMap<String,Target>(String.class,Target.class),byKey=new ObjectMap<String,Target>(String.class,Target.class);
}
