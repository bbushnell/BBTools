package prok;

import java.io.File;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Locale;
import java.util.TreeMap;
import java.util.TreeSet;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectIntMap;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Prepare assembly-specific 5S windows and compare genome/window RF00001 labels.
 * Raw CM parsing and eligibility are shared with RrnaCmTruthBuilder. No annotation
 * is edited and no mismatch is discarded to manufacture parity.
 * @author Ganyu
 */
public final class RrnaCmWindowParity {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,RrnaCmWindowParity.class,false);
		String mode=null,assembly=null,windows=null,genome=null,whole=null,flank=null,prefix=null;int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0,"Expected key=value: "+arg);
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(key.equals("mode")){mode=v;}else if(key.equals("assembly")){assembly=v;}
			else if(key.equals("windows")){windows=v;}else if(key.equals("genome")){genome=v;}
			else if(key.equals("whole")){whole=v;}else if(key.equals("flank")){flank=v;}
			else if(key.equals("outprefix")){prefix=v;}else if(key.equals("expected")){expected=Parse.parseIntKMG(v);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(assembly!=null&&windows!=null&&prefix!=null&&expected>0,"Required assembly= windows= outprefix= expected=positive-window-count");
		Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		if("prepare".equals(mode)){require(genome!=null,"prepare requires genome=");prepare(assembly,windows,genome,prefix,expected);}
		else if("compare".equals(mode)){require(whole!=null&&flank!=null,"compare requires whole= flank=");compare(assembly,windows,whole,flank,prefix,expected);}
		else{throw new IllegalArgumentException("mode must be prepare or compare");}
	}

	public static void prepare(String assembly,String windows,String genome,String prefix,int expected){
		assert(expected>0):"The reviewed canary requires known, nonempty 5S membership";
		fresh(prefix,".windows.fa.gz",".windows.tsv",".space.tsv");
		final TreeMap<String,Window> selected=loadWindows(windows,assembly,expected);
		final ObjectIntMap<String> lengths=new ObjectIntMap<String>(String.class);
		final Streamer st=RrnaTargetExtractionCheck.stream(genome);long residues=0;boolean done=false;st.start();
		try{
			for(ListNum<Read> ln;(ln=st.nextList())!=null;){
				for(Read r:ln.list){
					final String id=token(r.id);require(r.bases!=null&&r.length()>0,"Empty genome record: "+id);
					require(!lengths.contains(id),"Duplicate genome seqid: "+id);lengths.put(id,r.length());
					residues=Math.addExact(residues,r.length());
				}st.returnList(ln);
			}
			require(!st.errorState()&&residues>0,"Genome read failed or empty");done=true;
		}finally{if(!done){st.close();}}
		for(Window w:selected.values()){
			require(lengths.contains(w.seqid)&&w.low>=1&&w.high<=lengths.get(w.seqid),"Window outside exact genome contig: "+w.id);
		}
		final String z=zArgument(residues);
		final Writer writer=WriterFactory.makeWriter(FileFormat.testOutput(prefix+".windows.fa.gz",FileFormat.FASTA,null,true,false,false,false),null,1,null,false);
		writer.start();final ArrayList<Read> reads=new ArrayList<Read>();
		for(Window w:selected.values()){reads.add(new Read(w.bases,null,w.header,reads.size()));}
		writer.addReads(new ListNum<Read>(reads,0));require(!writer.poisonAndWait(),"Window FASTA output failed");
		final ByteStreamWriter ledger=writer(prefix+".windows.tsv"),space=writer(prefix+".space.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			ledger.println("assembly\twindow\tseqid\textract_start\textract_stop\tsource_strand\tlength");
			for(Window w:selected.values()){b.clear().append(assembly).tab().append(w.id).tab().append(w.seqid).tab().append(w.low).tab().append(w.high).tab().append(w.strand).tab().append(w.length).nl();ledger.print(b);}
			space.println("assembly\tgenome\tcontigs\tresidues\teffective_residues\tz_megabases\tz_roundtrip_residues\twindows");
			b.clear().append(assembly).tab().append(genome).tab().append(lengths.size()).tab().append(residues).tab().append(Math.multiplyExact(residues,2));
			b.tab().append(z).tab().append((long)(Double.parseDouble(z)*1000000.0)).tab().append(selected.size()).nl();space.print(b);
		}finally{finish(ledger,space);}
		System.err.println("CM_WINDOW_PREPARE_PASS assembly="+assembly+" windows="+selected.size()+" residues="+residues+" Z="+z);
	}

	/** Infernal1.1.5 cmsearch.c:541 casts parsed -Z*1e6 to int64; explicit Z is not doubled. */
	public static String zArgument(long residues){
		require(residues>0&&residues<=(1L<<51),"Residues outside exact double/integer search-space range");
		final long target=Math.multiplyExact(residues,2);
		String text=BigDecimal.valueOf(target).movePointLeft(6).toPlainString();
		double value=Double.parseDouble(text);
		for(int i=0;(long)(value*1000000.0)!=target&&i<8;i++){
			value=(long)(value*1000000.0)<target?Math.nextUp(value):Math.nextDown(value);text=Double.toString(value);
		}
		require((long)(Double.parseDouble(text)*1000000.0)==target,"Cannot express exact Infernal search space for "+residues);
		return text;
	}

	static TreeMap<String,Window> loadWindows(String path,String assembly,int expected){
		assert(path!=null&&assembly!=null):"Window origin requires a concrete input and exact assembly identity";
		final TreeMap<String,Window> result=new TreeMap<String,Window>();
		final Streamer st=RrnaTargetExtractionCheck.stream(path);boolean done=false;st.start();
		try{
			for(ListNum<Read> ln;(ln=st.nextList())!=null;){
				for(Read r:ln.list){
					final String family=field(r.id,"family");require(family.equals("5S")||family.equals("SSU"),"Unknown checked family");
					if(!family.equals("5S")||!field(r.id,"assembly_id").equals(assembly)){continue;}
					final Window w=new Window(r);require(result.put(w.id,w)==null,"Duplicate window identity: "+w.id);
				}st.returnList(ln);
			}
			require(!st.errorState(),"Window input failed");done=true;
		}finally{if(!done){st.close();}}
		require(result.size()==expected,"Selected5S windows "+result.size()+" != "+expected+" for "+assembly);return result;
	}

	public static void compare(String assembly,String windows,String whole,String flank,String prefix,int expected){
		fresh(prefix,".hits.tsv.gz",".loci.tsv",".summary.tsv");
		final TreeMap<String,Window> ws=loadWindows(windows,assembly,expected);
		final ArrayList<RrnaCmTruthBuilder.Hit> wh=RrnaCmTruthBuilder.readHits(new RrnaCmTruthBuilder.Genome("whole",assembly,whole,"launcher_verified"));
		final ArrayList<RrnaCmTruthBuilder.Hit> fh=RrnaCmTruthBuilder.readHits(new RrnaCmTruthBuilder.Genome("window",assembly,flank,"launcher_verified"));
		final ArrayList<Mapped> all=new ArrayList<Mapped>();
		final TreeMap<String,ArrayList<Mapped>> g=new TreeMap<String,ArrayList<Mapped>>(),f=new TreeMap<String,ArrayList<Mapped>>();
		for(RrnaCmTruthBuilder.Hit h:wh){final Mapped m=new Mapped(h,null);all.add(m);if(h.status.equals("QUALIFIED")){add(g,m);}}
		for(RrnaCmTruthBuilder.Hit h:fh){final Window w=ws.get(h.seqid);require(w!=null,"CM hit names an unknown window: "+h.seqid);final Mapped m=new Mapped(h,w);all.add(m);if(h.status.equals("QUALIFIED")){add(f,m);}}
		final TreeSet<String> keys=new TreeSet<String>(g.keySet());keys.addAll(f.keySet());
		long matched=0,genomeOnly=0,windowOnly=0,duplicateOrigins=0,scoreDifferences=0,evalueDifferences=0;
		final ByteStreamWriter hits=writer(prefix+".hits.tsv.gz"),loci=writer(prefix+".loci.tsv"),summary=writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			hits.println("assembly\tarm\torigin\ttableLine\tseqid\tstart\tstop\tstrand\tmdl_from\tmdl_to\ttrunc\tinclusion\tscore\tevalue\tstatus");
			for(Mapped m:all){
				final RrnaCmTruthBuilder.Hit h=m.hit;b.clear().append(assembly).tab().append(m.window==null?"whole":"window").tab().append(h.seqid).tab().append(h.line);
				b.tab().append(m.seqid).tab().append(m.start).tab().append(m.stop).tab().append(m.strand).tab().append(h.mdlFrom).tab().append(h.mdlTo).tab().append(h.trunc).tab().append(h.included?'!':'?').tab().append(h.score).tab().append(h.evalue).tab().append(h.status).nl();hits.print(b);
			}
			loci.println("assembly\tseqid\tstart\tstop\tstrand\tmdl_from\tmdl_to\ttrunc\tinclusion\tcomparison\twindowGeometry\tgenomeOrigins\twindowOrigins\tgenomeScoreEvalue\twindowScoreEvalue");
			for(String key:keys){
				final ArrayList<Mapped> a=g.get(key),c=f.get(key);final Mapped m=a!=null?a.get(0):c.get(0);
				final String status=a!=null&&c!=null?"MATCH":a!=null?"GENOME_ONLY":"WINDOW_ONLY";
				if(a!=null&&c!=null){matched++;for(Mapped w:c){if(!a.get(0).hit.score.equals(w.hit.score)){scoreDifferences++;}if(!a.get(0).hit.evalue.equals(w.hit.evalue)){evalueDifferences++;}}}
				else if(a!=null){genomeOnly++;}else{windowOnly++;}
				if(c!=null){duplicateOrigins+=c.size()-1;}
				b.clear().append(assembly).tab().append(key).tab().append(status).tab().append(geometry(m,ws));
				b.tab().append(a==null?0:a.size()).tab().append(c==null?0:c.size()).tab().append(scores(a)).tab().append(scores(c)).nl();loci.print(b);
			}
			summary.println("assembly\tgenome_raw\twindow_raw\tgenome_qualified_loci\twindow_qualified_loci\tmatched\tgenome_only\twindow_only\tduplicate_window_origins\tmatched_origin_score_differences\tmatched_origin_evalue_differences\tgate");
			b.clear().append(assembly).tab().append(wh.size()).tab().append(fh.size()).tab().append(g.size()).tab().append(f.size()).tab().append(matched).tab().append(genomeOnly).tab().append(windowOnly).tab().append(duplicateOrigins).tab().append(scoreDifferences).tab().append(evalueDifferences).tab().append(genomeOnly+windowOnly==0?"PASS":"FAIL").nl();summary.print(b);
		}finally{finish(hits,loci,summary);}
		require(genomeOnly+windowOnly==0,"RF00001 label mismatch; inspect retained locus ledger: genomeOnly="+genomeOnly+" windowOnly="+windowOnly);
		System.err.println("CM_WINDOW_PARITY_PASS assembly="+assembly+" loci="+matched+" duplicateWindowOrigins="+duplicateOrigins);
	}

	static String geometry(Mapped m,TreeMap<String,Window> windows){
		assert(m!=null):"Geometry diagnoses a concrete mapped locus";boolean overlap=false;
		for(Window w:windows.values()){
			if(!w.seqid.equals(m.seqid)){continue;}
			if(m.start>=w.low&&m.stop<=w.high){return "ENCLOSED";}
			if(m.start<=w.high&&m.stop>=w.low){overlap=true;}
		}return overlap?"BOUNDARY_CENSORED":"OUTSIDE_WINDOWS";
	}
	static void add(TreeMap<String,ArrayList<Mapped>> map,Mapped m){
		assert(m.hit.status.equals("QUALIFIED")):"Only identical eligibility predicates enter the compared sets";
		final String key=m.key();ArrayList<Mapped> list=map.get(key);
		if(list==null){list=new ArrayList<Mapped>();map.put(key,list);}list.add(m);
	}
	static String scores(ArrayList<Mapped> list){
		if(list==null){return "NA";}final StringBuilder s=new StringBuilder();
		for(Mapped m:list){if(s.length()>0){s.append(';');}s.append(m.hit.seqid).append(':').append(m.hit.line).append(':').append(m.hit.score).append(':').append(m.hit.evalue);}return s.toString();
	}
	static final class Window {
		Window(Read r){
			header=r.id;id=token(r.id);bases=r.bases;length=r.length();seqid=field(r.id,"source_contig");
			final String strandText=field(r.id,"source_strand");require(strandText.equals("+")||strandText.equals("-"),"Bad source strand");strand=strandText.charAt(0);
			final int start=Integer.parseInt(field(r.id,"source_start")),stop=Integer.parseInt(field(r.id,"source_end"));
			final int left=Integer.parseInt(field(r.id,"lflank")),right=Integer.parseInt(field(r.id,"rflank"));
			require(start>0&&stop>=start&&left>=0&&right>=0&&bases!=null,"Invalid checked window metadata");
			final long lo=(long)start-(strand=='+'?left:right),hi=(long)stop+(strand=='+'?right:left);
			require(lo>=1&&hi<=Integer.MAX_VALUE&&hi-lo+1==length,"Window length/bounds disagree: "+id);low=(int)lo;high=(int)hi;
		}
		final String header,id,seqid;final byte[] bases;final int low,high,length;final char strand;
	}
	static final class Mapped {
		Mapped(RrnaCmTruthBuilder.Hit h,Window w){
			hit=h;window=w;final int lo=Math.min(h.seqFrom,h.seqTo),hi=Math.max(h.seqFrom,h.seqTo);
			if(w==null){seqid=h.seqid;start=lo;stop=hi;strand=h.strand;}
			else{
				require(lo>=1&&hi<=w.length,"CM hit exceeds window: "+h.seqid);
				seqid=w.seqid;start=w.strand=='+'?w.low+lo-1:w.high-hi+1;stop=w.strand=='+'?w.low+hi-1:w.high-lo+1;
				strand=w.strand=='+'?h.strand:h.strand=='+'?'-':'+';
			}
			assert(start>=1&&stop>=start):"Infernal and extraction coordinates are1-based inclusive";
		}
		String key(){return seqid+"\t"+start+"\t"+stop+"\t"+strand+"\t"+hit.mdlFrom+"\t"+hit.mdlTo+"\t"+hit.trunc+"\t"+(hit.included?'!':'?');}
		final RrnaCmTruthBuilder.Hit hit;final Window window;final String seqid;final int start,stop;final char strand;
	}
	static String token(String header){int n=0;while(n<header.length()&&!Character.isWhitespace(header.charAt(n))){n++;}require(n>0,"Empty FASTA ID");return header.substring(0,n);}
	static String field(String header,String key){return RrnaTargetExtractionCheck.field(header,key);}
	static void fresh(String prefix,String... suffixes){for(String s:suffixes){require(!new File(prefix+s).exists(),"Output exists: "+prefix+s);}}
	static ByteStreamWriter writer(String path){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));w.start();return w;}
	public static void finish(ByteStreamWriter... writers){boolean error=false;for(ByteStreamWriter w:writers){error|=w.poisonAndWait();}require(!error,"Parity output failed");}
	static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}
}
