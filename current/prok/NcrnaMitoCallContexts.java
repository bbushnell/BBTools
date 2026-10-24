package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;

/** Exports merged call neighborhoods from the previously sealed whole-contig FASTA.
 * Contexts retain forward source orientation, including calls on the minus strand.
 * This is extraction only, with no origin classification or new biological cutoff.
 * @author Brian, Ganyu */
public final class NcrnaMitoCallContexts {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoCallContexts.class,false);
		try{final NcrnaMitoCallContexts x=new NcrnaMitoCallContexts();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		for(String arg:args){final String[] p=Parse.splitOnFirst(arg,'=');require(p.length==2,"Use flag=value");final String key=p[0].toLowerCase(Locale.ROOT);
			require(KNOWN.contains("|"+key+"|") && flags.put(key,p[1])==null,"Unknown/repeated flag: "+key);}
		for(String key:new String[]{"in","contigs","cases","verified","outprefix"}){require(flags.get(key)!=null,"Missing input binding: "+key);}
		flank=integer("flank",20000);wholeBelow=integer("wholebelow",200000);expected=integer("expected",55);expectedContigs=integer("expectedcontigs",44);
		require(flank>=0 && wholeBelow>0 && expected>0 && expectedContigs>0,"Require nonnegative flank and positive source/case denominators");
		pins=EukPanelGffValidation.hashes(flags.get("verified"));for(String key:new String[]{"in","contigs","cases"}){require(pins.get(flags.get(key))!=null,"Producer must verify exact input: "+key);}
		Shared.setThreads(1);ReadWrite.setZipThreads(1);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=false;
	}
	private void run(){
		final String prefix=flags.get("outprefix");RrnaResourceIO.fresh(prefix,".fa.gz",".contexts.tsv",".calls.tsv",".summary.tsv");
		loadContigs();loadCalls();plan();
		final Writer out=WriterFactory.makeWriter(FileFormat.testOutput(prefix+".fa.gz",FileFormat.FASTA,null,true,false,false,true),null,1);out.start();final long[] batchId={0};
		try{RrnaResourceIO.stream(flags.get("in"),r->{final String id=RrnaResourceIO.first(r.id);final Source s=sources.get(id);
			require(s!=null && !s.seen && r.bases!=null && r.length()==s.length,"Frozen exported FASTA must contain each bound source exactly once: "+id);
			require(RrnaResourceIO.field(r.id,"assembly").equals(s.assembly),"Exported sequence assembly must match contig metadata: "+id);
			s.seen=true;long acgt=0;for(byte base:r.bases){if(RrnaResourceIO.base(base)>=0){acgt++;}}require(acgt==s.acgt,"Whole-source ACGT denominator changed: "+id);
			final ArrayList<Read> batch=new ArrayList<Read>(s.contexts.size());
			for(Context c:s.contexts){final byte[] bases=Arrays.copyOfRange(r.bases,c.start,c.stop+1);for(byte base:bases){if(RrnaResourceIO.base(base)>=0){c.acgt++;}}
				batch.add(new Read(bases,null,c.id+" source="+s.id+" start0="+c.start+" stop0="+c.stop+" orientation=FORWARD_SOURCE",batchId[0],false));}
			out.add(batch,batchId[0]++);
		});}finally{require(!out.poisonAndWait(),"Complete context FASTA output is required before reporting success");}
		for(Source s:sources.values()){if(s!=null){require(s.seen,"Missing frozen whole-contig source: "+s.id);}}
		write(prefix);
	}
	private void loadContigs(){
		final ByteFile in=RrnaResourceIO.open(flags.get("contigs"));final LineParser1 p=new LineParser1('\t');
		try{p.set(in.nextLine());require(p.terms()==12 && p.termEquals("exported_id",6) && p.termEquals("source_sha80",8),"Require sealed whole-contig evidence schema");
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==12,"Truncated contig evidence");if(p.termEquals("NA",6)){continue;}
				final Source s=new Source(p);require(sources.put(s.id,s)==null && byOriginal.put(s.assembly+'\t'+s.original,s)==null,"Export IDs and assembly/contig keys must be unique");}
		}finally{require(!in.close(),"Contig metadata read failed");}
		require(sources.size()==expectedContigs,"Retain the complete exported contig population");
	}
	private void loadCalls(){
		final ByteFile in=RrnaResourceIO.open(flags.get("cases"));final LineParser1 p=new LineParser1('\t');final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try{p.set(in.nextLine());require(p.terms()==42 && p.termEquals("case_id",0) && p.termEquals("audit_id",29),"Require sealed source-bound case schema");
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==42,"Truncated case evidence");if(!p.termEquals("NUCLEAR_NEW_UNMATCHED",1)){continue;}
				final Source s=byOriginal.get(p.parseString(2)+'\t'+p.parseString(3));require(s!=null,"Every selected call must bind an exported source");
				final Call c=new Call(p,s);require(ids.add(c.id) && p.parseInt(12)==s.length && p.parseLong(13)==s.acgt && p.termEquals(s.fasta,26) && p.termEquals(s.hash,27),"Case source provenance must exactly match whole-contig metadata: "+c.id);
				s.calls.add(c);calls.add(c);}
		}finally{require(!in.close(),"Case evidence read failed");}
		require(calls.size()==expected,"Every selected call must survive extraction exactly once");
	}
	private void plan(){
		final ArrayList<String> ids=new ArrayList<String>();for(String id:sources.keys()){if(id!=null){ids.add(id);}}Collections.sort(ids);
		for(String id:ids){final Source s=sources.get(id);require(!s.calls.isEmpty(),"Every exported candidate contig must contain at least one selected call");
			s.calls.sort((a,b)->{int n=Integer.compare(a.start,b.start);if(n==0){n=Integer.compare(a.stop,b.stop);}return n==0?a.id.compareTo(b.id):n;});Context previous=null;
			for(Call call:s.calls){final boolean whole=s.length<wholeBelow;
				final int start=whole?0:(int)Math.max(0L,(long)call.start-flank),stop=whole?s.length-1:(int)Math.min((long)s.length-1,(long)call.stop+flank);
				if(previous==null || (long)start>1L+previous.stop){previous=new Context("MCTX"+(contexts.size()+1),s,start,stop,whole);contexts.add(previous);s.contexts.add(previous);}
				else{previous.stop=Math.max(previous.stop,stop);}call.context=previous;
			}
		}
		for(Call c:calls){require(c.context.start<=c.start && c.context.stop>=c.stop,"Merged context must contain the complete call");}
	}
	private void write(String prefix){
		final ByteStreamWriter ct=RrnaResourceIO.writer(prefix+".contexts.tsv"),members=RrnaResourceIO.writer(prefix+".calls.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv");final ByteBuilder b=new ByteBuilder();long bases=0,originalBases=0;int whole=0;
		try{ct.println(CONTEXT_HEADER);members.println(CALL_HEADER);
			for(Context c:contexts){final Source s=c.source;bases+=c.stop-c.start+1;whole+=c.whole?1:0;
				ct.print(b.clear().append(c.id).tab().append(s.id).tab().append(s.assembly).tab().append(s.original).tab().append(c.start).tab().append(c.stop)
					.tab().append(c.stop-c.start+1).tab().append(c.acgt).tab().append(s.length).tab().append(s.acgt).tab().append(c.whole?"WHOLE_CONTIG":"CALL_NEIGHBORHOOD")
					.tab().append(s.fasta).tab().append(s.hash).tab().append(pins.get(flags.get("in"))).tab().append(pins.get(flags.get("contigs"))).tab().append(pins.get(flags.get("cases"))).nl());}
			calls.sort((a,z)->a.id.compareTo(z.id));for(Call c:calls){final Source s=c.source;
				members.print(b.clear().append(c.id).tab().append(c.context.id).tab().append(s.id).tab().append(s.assembly).tab().append(s.original).tab().append(c.start+1).tab().append(c.stop+1)
					.tab().append(c.strand).tab().append(c.start-c.context.start).tab().append(c.stop-c.context.start).tab().append(c.audit).tab().append(c.model).nl());}
			for(Source s:sources.values()){if(s!=null){originalBases+=s.length;}}
			summary.println("calls\tsource_contigs\tsource_bases\tcontexts\tcontext_bases\twhole_contig_contexts\tflank\twhole_below\torientation\tnew_alignments");
			summary.print(b.clear().append(calls.size()).tab().append(sources.size()).tab().append(originalBases).tab().append(contexts.size()).tab().append(bases).tab().append(whole).tab().append(flank).tab().append(wholeBelow).append("\tFORWARD_SOURCE\t0\n"));
		}finally{RrnaResourceIO.finish(ct,members,summary);}
		System.err.println("MITO_CALL_CONTEXTS_PASS calls="+calls.size()+" sources="+sources.size()+" contexts="+contexts.size()+" bases="+bases+" no_alignment=true");
	}
	private int integer(String key,int fallback){final String s=flags.get(key);return s==null?fallback:Parse.parseIntKMG(s);}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Source{
		Source(LineParser1 p){assembly=p.parseString(0);original=p.parseString(1);length=p.parseInt(2);acgt=p.parseLong(3);id=p.parseString(6);fasta=p.parseString(7);hash=p.parseString(8);require(length>0 && acgt>=0 && acgt<=length && hash.matches("[0-9a-f]{20}"),"Valid source length, ACGT and hash required");}
		final String assembly,original,id,fasta,hash;final int length;final long acgt;boolean seen;final ArrayList<Call> calls=new ArrayList<Call>();final ArrayList<Context> contexts=new ArrayList<Context>();
	}
	private static final class Call{
		Call(LineParser1 p,Source s){source=s;id=p.parseString(0);start=p.parseInt(4)-1;stop=p.parseInt(5)-1;strand=p.parseString(6);model=p.parseString(8);audit=p.parseString(29);require(start>=0 && stop>=start && stop<s.length && p.parseInt(7)==stop-start+1 && (strand.equals("+") || strand.equals("-")),"Call must fit its bound forward-source contig");}
		final Source source;final String id,strand,model,audit;final int start,stop;Context context;
	}
	private static final class Context{Context(String i,Source s,int a,int z,boolean w){id=i;source=s;start=a;stop=z;whole=w;}final String id;final Source source;final int start;int stop;final boolean whole;long acgt;}
	public static final String CONTEXT_HEADER="context_id\texported_contig\tassembly\toriginal_contig\tsource_start0\tsource_stop0\tcontext_length\tcontext_ACGT\tsource_length\tsource_ACGT\tscope\tsource_fasta\tsource_sha80\texport_sha80\tcontigs_sha80\tcases_sha80";
	public static final String CALL_HEADER="case_id\tcontext_id\texported_contig\tassembly\toriginal_contig\tcall_start1\tcall_stop1\tstrand\tcontext_call_start0\tcontext_call_stop0\taudit_id\twinning_model";
	private static final String KNOWN="|in|contigs|cases|verified|outprefix|flank|wholebelow|expected|expectedcontigs|";
	private final ObjectMap<String,String> flags=new ObjectMap<String,String>(String.class,String.class);private ObjectMap<String,String> pins;
	private final ObjectMap<String,Source> sources=new ObjectMap<String,Source>(String.class,Source.class),byOriginal=new ObjectMap<String,Source>(String.class,Source.class);
	private final ArrayList<Call> calls=new ArrayList<Call>();private final ArrayList<Context> contexts=new ArrayList<Context>();private int flank,wholeBelow,expected,expectedContigs;
}
