package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Review-only addendum to frozen truth. Sensitive CM support is not automatic admission.
 * @author Brian, Raiden */
public final class EukPanelHistoricalEvidence {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,EukPanelHistoricalEvidence.class,false);
		String windows=null,models=null,scans=null,truth=null,hash=null,out=null;int expected=724;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("windows")){windows=v;}else if(k.equals("models")){models=v;}else if(k.equals("scans")){scans=v;}else if(k.equals("truth")){truth=v;}
			else if(k.equals("truthsha80")){hash=v;}else if(k.equals("outprefix")){out=v;}else if(k.equals("expected")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException(k);}}
		require(windows!=null && models!=null && scans!=null && truth!=null && SeedOffsetTable.validSha80(hash) && out!=null && expected>0,"Bound windows/models/complete scans/frozen truth and fresh output required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		final EukPanelHistoricalEvidence x=new EukPanelHistoricalEvidence();x.loadModels(models);x.loadTruth(truth);x.loadWindows(windows,expected);
		x.write(scans,out,hash);Shared.closeStream(pp.outstream);
	}
	void loadModels(String path){
		assert(models.isEmpty()):"One independent model roster defines all folds and GA thresholds";
		read(path,EukPanelHistoricalReview.MODEL_HEADER,p->{require(p.terms()==8,"Model schema changed");final Model m=new Model(p);
			require(models.put(m.acc,m)==null && m.clen>0 && Double.isFinite(m.ga),"Model IDs/lengths/thresholds must be unambiguous");});
	}
	void loadTruth(String path){
		assert(truth.isEmpty()):"Candidate addition checks use one unchanged frozen reference";
		read(path,EukPanelTruthFreeze.HEADER,p->{require(p.terms()==35 && p.termEquals("TRUTH_V1",32),"Frozen baseline required");
			if(!p.termEquals("ADMITTED",26)){return;}
			final String key=key(p.parseString(1),p.parseString(2),p.parseString(3),p.parseString(6));
			ArrayList<Locus> list=truth.get(key);if(list==null){list=new ArrayList<Locus>();truth.put(key,list);}
			list.add(new Locus(p.parseString(0),p.parseLong(4),p.parseLong(5)));});
	}
	void loadWindows(String path,int expected){
		assert(rows.isEmpty()):"Each source annotation appears exactly once in the review addendum";
		final ObjectSet<String> ids=new ObjectSet<String>(String.class),sources=new ObjectSet<String>(String.class);
		read(path,EukPanelHistoricalReview.WINDOW_HEADER,p->{require(p.terms()==20,"Window schema changed");final Row r=new Row(p);
			require(ids.add(r.id) && sources.add(r.assembly+'\t'+p.parseString(9)+'\t'+p.parseLong(11)),"Duplicate mask/source identity");
			require(r.start>=1 && r.end>=r.start && r.left>=1 && r.left<=r.start && r.right>=r.end && r.right<=p.parseLong(15),"Window must contain the exact original interval");
			require(r.strand.equals("+") || r.strand.equals("-"),"Historical strand unresolved");rows.add(r);});
		require(rows.size()==expected,"Mask conservation failed");
	}
	void hits(Row r,String path){
		assert(r.hits.isEmpty()):"A sensitive scan may be consumed only once";
		final ByteFile in=RrnaResourceIO.open(path);final RrnaCmTruthBuilder.Fields f=new RrnaCmTruthBuilder.Fields();long lineNo=0;boolean program=false,version=false,done=false;
		try{for(byte[] line;(line=in.nextLine())!=null;){lineNo++;if(line.length==0){continue;}require(!done,"Content after completed CM table");
			if(line[0]=='#'){final String s=Euk18sShredCmSummary.text(line);if(s.startsWith("# Program:")){require(s.substring(10).trim().equals("cmsearch"),"Expected cmsearch");program=true;}
				if(s.startsWith("# Version:")){require(s.substring(10).trim().startsWith("1.1.5 "),"Pinned Infernal1.1.5 required");version=true;}if(s.equals("# [ok]")){done=true;}continue;}
			f.set(line);require(f.eq(0,r.id) && f.eq(4,"cm"),"Hit must address this exact mask's window using a CM");
			final Model m=models.get(f.text(3));require(m!=null && m.family.equals(r.family),"Scan must cover declared same-fold models only");
			final Hit h=new Hit(r,m,lineNo,f);r.hits.add(h);
			if(h.overlap==0){continue;}
			if(h.strand.equals(r.strand)){
				if(better(h,r.same)){r.same=h;}
				if(h.bits<m.ga && better(h,r.subga)){r.subga=h;}
				if(h.reciprocal && h.bits>=m.ga && better(h,r.ga)){r.ga=h;}
				if(h.reciprocal && h.full && h.trunc.equals("no") && h.bits<m.ga && h.evalue<=0.01 && better(h,r.support)){r.support=h;}
			}else{
				if(better(h,r.opposite)){r.opposite=h;}
				if(h.reciprocal && (h.bits>=m.ga || h.full && h.trunc.equals("no") && h.evalue<=0.01) && better(h,r.oppositeSupport)){r.oppositeSupport=h;}
			}
		}}finally{require(!in.close(),"Sensitive table read failed: "+path);}
		require(program && version && done,"Sensitive table lacks version/program/completion evidence: "+path);
	}
	static boolean better(Hit a,Hit b){
		assert(a!=null):"Ranking requires actual evidence; null means no previous candidate";
		if(b==null){return true;}if(a.overlap!=b.overlap){return a.overlap>b.overlap;}
		final long ca=(long)a.span*b.model.clen,cb=(long)b.span*a.model.clen;if(ca!=cb){return ca>cb;}
		if(a.bits!=b.bits){return a.bits>b.bits;}if(a.evalue!=b.evalue){return a.evalue<b.evalue;}
		final int c=a.model.acc.compareTo(b.model.acc);return c!=0?c<0:a.line<b.line;
	}
	void write(String scans,String prefix,String hash){
		assert(!rows.isEmpty()):"The candidate never drops unsupported original annotations";
		RrnaResourceIO.fresh(prefix,".truth_v2_candidate.tsv.gz",".summary.tsv");
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".truth_v2_candidate.tsv.gz"),sum=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ObjectMap<String,long[]> bins=new ObjectMap<String,long[]>(String.class,long[].class);
		final ObjectMap<String,String> groups=new ObjectMap<String,String>(String.class,String.class);final ByteBuilder b=new ByteBuilder();
		try{
			out.println(EukPanelHistoricalReview.WINDOW_HEADER+"\tannotation_length\tproduct\tannotation_group\trecommended_class\tproposed_action\treason\topposite_strand_flag\texisting_loci\traw_table\tbaseline_truth_sha80\tcandidate_status"+columns("same")+columns("subga")+columns("opposite")+columns("recommended")+columns("opposite_support")+"\tsame_max_overlap_ties\topposite_max_overlap_ties");
			for(Row r:rows){
				final String path=scans+'/'+r.id+"/hits.tbl";hits(r,path);
				final Hit recommended=r.ga!=null?r.ga:r.support;
				final String classification=r.ga!=null?"SUPPORTED_GA":r.support!=null?"SUPPORTED_SUBGA":"UNSUPPORTED";
				final boolean opposite=r.oppositeSupport!=null;
				final String reason=r.ga!=null?"GA_PASS_RECIPROCAL_HALF":r.support!=null?"FULL_MODEL_UNTRUNCATED_RECIPROCAL_HALF_E_LE_0.01_AT_1MB":r.same==null?"NO_SAME_STRAND_OVERLAPPING_HIT_AT_T0":"SAME_STRAND_EVIDENCE_FAILS_RECOMMENDED_SUPPORT_RULE";
				final String overlaps=existing(r,recommended);
				final String action=recommended!=null?(overlaps.equals("NONE")?"PROPOSE_ADDITION_REQUIRES_ADJUDICATION":"LINK_EXISTING_REQUIRES_ADJUDICATION"):(opposite?"REVIEW_STRAND_KEEP_MASK":"KEEP_UNRESOLVED_MASK");
				final String geometry=key(r.assembly,r.seqid,r.strand,r.family)+'\t'+r.start+'\t'+r.end;
				String group=groups.get(geometry);if(group==null){group="A"+(groups.size()+1);groups.put(geometry,group);}
				b.clear().append(r.raw).tab().append(r.end-r.start+1).tab().append(product(r.attrs)).tab().append(group).tab().append(classification).tab().append(action).tab().append(reason).tab().append(opposite?1:0).tab().append(overlaps).tab().append(path).tab().append(hash).tab().append("TRUTH_V2_CANDIDATE_ADDENDUM_NOT_FROZEN");
				append(b,r.same);append(b,r.subga);append(b,r.opposite);append(b,recommended);append(b,r.oppositeSupport);
				int sameTies=0,oppositeTies=0;for(Hit h:r.hits){if(h.overlap>0){
					if(h.strand.equals(r.strand) && r.same!=null && h.overlap==r.same.overlap){sameTies++;}
					else if(!h.strand.equals(r.strand) && r.opposite!=null && h.overlap==r.opposite.overlap){oppositeTies++;}
				}}
				out.print(b.tab().append(sameTies).tab().append(oppositeTies).nl());
				EukPanelCmProposal.increment(bins,r.family+'\t'+r.compartment+'\t'+r.history+'\t'+classification+'\t'+action,1);
			}
			sum.println("fold_family\tcompartment\thistory_class\trecommended_class\tproposed_action\tannotations");long total=0;
			for(String key:EukPanelCmProposal.sorted(bins)){final long n=bins.get(key)[0];sum.print(b.clear().append(key).tab().append(n).nl());total+=n;}
			require(total==rows.size(),"Summary bins must conserve every annotation");
			System.err.println("HISTORICAL_REVIEW_COMPLETE annotations="+rows.size()+" exact_annotation_geometries="+groups.size()+" frozen_truth_changed=false");
		}finally{RrnaResourceIO.finish(out,sum);}
	}
	String existing(Row r,Hit h){
		assert(r!=null):"Locus matching is assembly/fold/strand-qualified, never bare feature ID";
		if(h==null){return "NONE";}final ArrayList<Locus> list=truth.get(key(r.assembly,r.seqid,r.strand,r.family));if(list==null){return "NONE";}
		final ByteBuilder b=new ByteBuilder();for(Locus l:list){if(reciprocal(h.start,h.end,l.start,l.end)){if(b.length()>0){b.append(',');}b.append(l.id);}}
		return b.length()==0?"NONE":b.toString();
	}
	static String product(String attrs){
		assert(attrs!=null):"Original attributes remain alongside decoded product navigation";
		for(String term:attrs.split(";")){if(term.startsWith("product=")){return term.substring(8);}}return "NA";
	}
	static String columns(String p){return "\t"+p+"_model\t"+p+"_table_line\t"+p+"_start\t"+p+"_end\t"+p+"_strand\t"+p+"_model_from\t"+p+"_model_to\t"+p+"_clen\t"+p+"_trunc\t"+p+"_bits\t"+p+"_evalue_1mb\t"+p+"_ga\t"+p+"_ga_margin\t"+p+"_overlap_bp\t"+p+"_annotation_fraction\t"+p+"_hit_fraction";}
	static void append(ByteBuilder b,Hit h){
		assert(b!=null):"Every optional hit occupies the same16 output columns";
		if(h==null){for(int i=0;i<16;i++){b.append("\tNA");}return;}
		b.tab().append(h.model.acc).tab().append(h.line).tab().append(h.start).tab().append(h.end).tab().append(h.strand).tab().append(h.mf).tab().append(h.mt).tab().append(h.model.clen).tab().append(h.trunc).tab().append(h.bits,4).tab().append(h.evalueText).tab().append(h.model.ga,4).tab().append(h.bits-h.model.ga,4).tab().append(h.overlap).tab().append(h.sourceFraction,9).tab().append(h.hitFraction,9);
	}
	public static boolean reciprocal(long a,long b,long c,long d){require(a>0 && b>=a && c>0 && d>=c,"Positive one-based inclusive intervals required");final long n=Math.max(0,Math.min(b,d)-Math.max(a,c)+1);return n>0 && n/(double)(b-a+1)>=0.5 && n/(double)(d-c+1)>=0.5;}
	static String key(String a,String s,String strand,String f){return a+'\t'+s+'\t'+strand+'\t'+f;}
	static void read(String path,String header,java.util.function.Consumer<LineParser1> body){
		assert(path!=null && header!=null && body!=null):"A declared schema governs each native byte-file read";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,header,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);body.accept(p);}}finally{require(!in.close(),"Table read failed: "+path);}
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Model{Model(LineParser1 p){acc=p.parseString(1);family=p.parseString(2);clen=p.parseInt(4);ga=p.parseDouble(5);}final String acc,family;final int clen;final double ga;}
	static final class Locus{Locus(String i,long a,long b){id=i;start=a;end=b;}final String id;final long start,end;}
	static final class Row{
		Row(LineParser1 p){raw=p.line().clone();id=p.parseString(0);assembly=p.parseString(1);seqid=p.parseString(2);strand=p.parseString(3);start=p.parseLong(4);end=p.parseLong(5);family=p.parseString(6);history=p.parseString(7);attrs=p.parseString(13);compartment=p.parseString(16);left=p.parseLong(17);right=p.parseLong(18);}
		final String id,assembly,seqid,strand,family,history,attrs,compartment;final byte[] raw;final long start,end,left,right;final ArrayList<Hit> hits=new ArrayList<Hit>();Hit same,subga,opposite,oppositeSupport,ga,support;
	}
	static final class Hit{
		Hit(Row r,Model m,long l,RrnaCmTruthBuilder.Fields f){model=m;line=l;mf=f.number(5);mt=f.number(6);span=mt-mf+1;full=mf==1 && mt==m.clen;
			final int from=f.number(7),to=f.number(8);strand=f.text(9);trunc=f.text(10);bits=Double.parseDouble(f.text(14));evalueText=f.text(15);evalue=Double.parseDouble(evalueText);
			require(mf>=1 && mt>=mf && mt<=m.clen && from>0 && to>0 && Math.max(from,to)<=r.right-r.left+1,"Hit must fit declared CM/window bounds");
			require(strand.equals("+") && from<=to || strand.equals("-") && from>=to,"CM orientation inconsistent with coordinates");
			require(Double.isFinite(bits) && bits>=0 && Double.isFinite(evalue) && evalue>=0,"Reported CM scores must respect T0 and finite E-values");
			require(trunc.equals("no") || trunc.equals("5'") || trunc.equals("3'") || trunc.equals("5'&3'"),"Unknown truncation flag");
			start=r.left+Math.min(from,to)-1;end=r.left+Math.max(from,to)-1;overlap=Math.max(0,Math.min(end,r.end)-Math.max(start,r.start)+1);
			sourceFraction=overlap/(double)(r.end-r.start+1);hitFraction=overlap/(double)(end-start+1);reciprocal=reciprocal(start,end,r.start,r.end);
		}
		final Model model;final long line,start,end,overlap;final int mf,mt,span;final String strand,trunc,evalueText;final double bits,evalue,sourceFraction,hitFraction;final boolean full,reciprocal;
	}
	final ObjectMap<String,Model> models=new ObjectMap<String,Model>(String.class,Model.class);
	final ObjectMap<String,ArrayList<Locus>> truth=new ObjectMap<String,ArrayList<Locus>>(String.class,(Class)ArrayList.class);
	final ArrayList<Row> rows=new ArrayList<Row>();
}
