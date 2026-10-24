package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.LongList;

/** Joins source annotations to exact known query/reference aligned segments.
 * CM status describes the source annotation, not a new CM scan of the call.
 * Named SSU evidence remains independent, including contradictory CM labels.
 * @author Brian, Ganyu */
public final class NcrnaMitoSsuAnnotationJoin {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoSsuAnnotationJoin.class,false);
		try{final NcrnaMitoSsuAnnotationJoin x=new NcrnaMitoSsuAnnotationJoin();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		for(String arg:args){final String[] p=Parse.splitOnFirst(arg,'=');require(p.length==2,"Use flag=value");final String key=p[0].toLowerCase(Locale.ROOT);
			require(KNOWN.contains("|"+key+"|") && flags.put(key,p[1])==null,"Unknown/repeated flag: "+key);}
		for(String key:new String[]{"calls","alignments","loci","verified","outprefix"}){require(flags.get(key)!=null,"Missing source binding: "+key);}
		expected=flags.get("expected")==null?55:Parse.parseIntKMG(flags.get("expected"));require(expected>0,"Require the complete positive call denominator");
		pins=EukPanelGffValidation.hashes(flags.get("verified"));for(String key:new String[]{"calls","alignments","loci"}){require(pins.get(flags.get(key))!=null,"Producer must verify exact input: "+key);}
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
	}
	private void run(){
		final String prefix=flags.get("outprefix");RrnaResourceIO.fresh(prefix,".calls.tsv",".annotations.tsv.gz",".summary.tsv");
		loadCalls();loadAnnotations();join();write(prefix);
	}
	private void loadCalls(){
		final ByteFile in=RrnaResourceIO.open(flags.get("calls"));final LineParser1 p=new LineParser1('\t');
		try{final byte[] header=in.nextLine();require(header!=null,"Call table header required");p.set(header);
			require(p.terms()==32 && p.termEquals("case_id",0) && p.termEquals("any_supported_bp",10) && p.termEquals("alignment_identity_min",31),"Expected hard-filtered mapping call schema");callHeader=header;
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==32 && p.parseDouble(31)==.76,"Only the explicit measured .76 alignment filter is supported");
				final Call c=new Call(p,row);require(calls.put(c.id,c)==null,"Every call ID must occur exactly once: "+c.id);}
		}finally{require(!in.close(),"Call table read failed");}
		require(calls.size()==expected,"All expected calls, including zero-hit calls, must remain");
	}
	private void loadAnnotations(){
		final ByteFile in=RrnaResourceIO.open(flags.get("loci"));final LineParser1 p=new LineParser1('\t');final ObjectSet<String> ids=new ObjectSet<String>(String.class);long line=1;
		try{p.set(in.nextLine());require(p.terms()==28 && p.termEquals("locus_id",0) && p.termEquals("status",8) && p.termEquals("source_segments_in_transcript_order",23) && p.termEquals("source_gbff",26),"Expected accepted population-evidence locus schema");
			for(byte[] row;(row=in.nextLine())!=null;){line++;p.set(row);require(p.terms()==28 && p.termEquals("EXTRACTED",16),"Every accepted source annotation must retain extracted geometry");
				final Annotation a=new Annotation(p,line);require(ids.add(a.id) && p.termEquals(a.id,9) && p.termEquals(a.reference,10),"Annotation identity must bind the original source locus/accession");
				ArrayList<Annotation> list=annotations.get(a.reference);if(list==null){list=new ArrayList<Annotation>();annotations.put(a.reference,list);}list.add(a);annotationCount++;}
		}finally{require(!in.close(),"Annotation evidence read failed");}
		require(annotationCount>0,"The bound annotation census cannot be empty");
	}
	private void join(){
		final ByteFile in=RrnaResourceIO.open(flags.get("alignments"));final LineParser1 p=new LineParser1('\t');long line=1;
		try{p.set(in.nextLine());require(p.terms()==19 && p.termEquals("case_id",0) && p.termEquals("original_supported_intervals0",14) && p.termEquals("reference_supported_intervals1",16) && p.termEquals("passes_alignment_identity_076",18),"Expected exact paired call-alignment intervals and whole-alignment filter flag");
			for(byte[] row;(row=in.nextLine())!=null;){line++;alignmentRows++;p.set(row);require(p.terms()==19,"Truncated call alignment row");
				final Call c=calls.get(p.parseString(0));require(c!=null && p.termEquals(c.context,1),"Alignment row must bind an original call/context");final int flag=p.parseInt(4);
				require((flag&(1|4|2048))==0,"Only mapped single-ended primary/secondary alignments are supported");
				final boolean primary=(flag&256)==0;require(p.termEquals(primary?"true":"false",5),"Primary flag disagrees with SAM flag");
				final boolean passes=bool(p,18);if(!passes){rejectedRows++;continue;}passingRows++;
				intervals(p.parseString(14),query,false);intervals(p.parseString(16),reference,true);
				require(query.size==reference.size,"Known query and reference segments must be paired in original CIGAR order");
				final boolean reverse=(flag&16)!=0;long paired=0;
				for(int i=0;i<query.size;i++){final int qlo=lo(query.array[i]),qhi=hi(query.array[i]),rlo=lo(reference.array[i]),rhi=hi(reference.array[i]);
					require(qlo>=c.start && qhi<=c.stop && qhi-qlo==rhi-rlo,"Paired segments must have equal lengths and fit their original call");paired+=(long)qhi-qlo+1;
					if(i>0){require(lo(reference.array[i])>hi(reference.array[i-1]) && (reverse?qhi<lo(query.array[i-1]):qlo>hi(query.array[i-1])),"Reference segments ascend while query segments follow SAM orientation");}}
				require(paired==p.parseLong(6)+p.parseLong(7),"Known paired intervals must equal restricted =/X counts; I/D/M cannot become supported bases");c.mapped.addAll(query);
				final ArrayList<Annotation> list=annotations.get(p.parseString(3));if(list==null){continue;}
				final char mappedStrand=(c.strand=='-')^reverse?'-':'+';
				for(Annotation a:list){overlap.clear();for(int i=0;i<query.size;i++){
					final int qlo=lo(query.array[i]),qhi=hi(query.array[i]),rlo=lo(reference.array[i]),rhi=hi(reference.array[i]);
					for(NcrnaMitoCmEvidence.Span segment:a.segments){final int from=Math.max(rlo,segment.lo),to=Math.min(rhi,segment.hi);if(from>to){continue;}
						final int start=reverse?qhi-(to-rlo):qlo+(from-rlo),stop=reverse?qhi-(from-rlo):qlo+(to-rlo);overlap.add(pack(start,stop));}}
					if(overlap.size==0){continue;}Overlap o=c.overlaps.get(a.id);if(o==null){o=new Overlap(a);c.overlaps.put(a.id,o);}o.rows++;if(primary){o.primaryRows++;}o.firstLine=Math.min(o.firstLine,line);o.lastLine=line;
					final boolean strandMatch=a.strand==mappedStrand;
					if(!strandMatch){o.opposite.addAll(overlap);continue;}o.query.addAll(overlap);c.support[a.category].addAll(overlap);
					if(a.named){c.support[4].addAll(overlap);if(a.category==1 || a.category==2){c.support[5].addAll(overlap);}}
				}
			}
		}finally{require(!in.close(),"Call alignment evidence read failed");}
		for(Call c:calls.values()){if(c!=null){require(NcrnaMitoContigMappingSummary.union(c.mapped)==c.known,"Passing aligned query union must reproduce the independently bound mapping call coverage: "+c.id);}}
	}
	private void write(String prefix){
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".calls.tsv"),detail=RrnaResourceIO.writer(prefix+".annotations.tsv.gz"),summary=RrnaResourceIO.writer(prefix+".summary.tsv");final ByteBuilder b=new ByteBuilder();final long[] counts=new long[SUPPORT.length];
		try{b.clear().append(callHeader);for(String name:SUPPORT){b.tab().append(name).append("_any\t").append(name).append("_query_bp\t").append(name).append("_call_fraction");}
			out.print(b.append("\tcalls_sha80\talignments_sha80\tannotation_table_sha80\tinterpretation\n"));
			detail.println("case_id\tcontext_id\treference\tsource_locus_id\tsource_table_line\tsource_segments\tsource_strand\tsource_product\tsource_CM_status\tCM_category\tnamed_SSU\tsame_strand_query_bp\tsame_strand_call_fraction\topposite_or_mixed_strand_query_bp\tpassing_alignment_rows\tprimary_alignment_rows\tfirst_alignment_table_line\tlast_alignment_table_line\toriginal_query_intervals0\tsource_table\tsource_table_sha80\tsource_gbff\tsource_partial\tsource_pseudo\tsource_raw_location");
			for(String id:keys(calls)){final Call c=calls.get(id);final int length=c.stop-c.start+1;b.clear().append(c.raw);
				for(int i=0;i<SUPPORT.length;i++){final long bp=NcrnaMitoContigMappingSummary.union(c.support[i]);require(bp<=c.known && bp<=length,"Annotation-supported query bases cannot exceed passing known call coverage");if(bp>0){counts[i]++;}b.tab().append(bp>0?"yes":"no").tab().append(bp).tab().append(bp/(double)length,9);}
				out.print(b.tab().append(pins.get(flags.get("calls"))).tab().append(pins.get(flags.get("alignments"))).tab().append(pins.get(flags.get("loci"))).append("\tSOURCE_ANNOTATION_CM_STATUS_NOT_NEW_CALL_CM_SCAN_OR_FUNCTION\n"));
				for(String aid:keys(c.overlaps)){final Overlap o=c.overlaps.get(aid);final Annotation a=o.annotation;final long bp=NcrnaMitoContigMappingSummary.union(o.query),opposite=NcrnaMitoContigMappingSummary.union(o.opposite);
					b.clear().append(c.id).tab().append(c.context).tab().append(a.reference).tab().append(a.id).tab().append(a.line).tab().append(a.geometry).tab().append(a.strandText).tab().append(a.product).tab().append(a.status).tab().append(CATEGORIES[a.category]).tab().append(a.named?"yes":"no")
						.tab().append(bp).tab().append(bp/(double)length,9).tab().append(opposite).tab().append(o.rows).tab().append(o.primaryRows).tab().append(o.firstLine).tab().append(o.lastLine).tab();appendIntervals(b,o.query);
					detail.print(b.tab().append(flags.get("loci")).tab().append(pins.get(flags.get("loci"))).tab().append(a.gbff).tab().append(a.partial).tab().append(a.pseudo).tab().append(a.rawLocation).nl());}
			}
			b.clear().append("calls\tsource_annotations\talignment_rows\tpassing_alignment_rows\trejected_alignment_rows");for(String name:SUPPORT){b.tab().append(name).append("_calls");}summary.print(b.nl());
			b.clear().append(calls.size()).tab().append(annotationCount).tab().append(alignmentRows).tab().append(passingRows).tab().append(rejectedRows);for(long n:counts){b.tab().append(n);}summary.print(b.nl());
		}finally{RrnaResourceIO.finish(out,detail,summary);}
		System.err.println("MITO_SSU_ANNOTATION_JOIN_PASS calls="+calls.size()+" annotations="+annotationCount+" no_new_alignment_or_CM_scan=true");
	}
	/** Exact 12S product or explicit small-subunit wording; never a generic 16S name. */
	public static boolean namedSsu(String raw){
		require(raw!=null,"Preserve the actual source product");String value=raw.trim();if(value.length()>=2 && value.charAt(0)=='"' && value.charAt(value.length()-1)=='"'){value=value.substring(1,value.length()-1).trim();}
		value=value.toLowerCase(Locale.ROOT);return value.equals("12s ribosomal rna") || value.contains("small subunit") || value.contains("small-subunit");
	}
	private static int category(String status){for(int i=0;i<STATUSES.length;i++){if(status.equals(STATUSES[i])){return i;}}throw new IllegalArgumentException("Unrecognized accepted annotation CM status: "+status);}
	private static void intervals(String value,LongList list,boolean oneBased){
		list.clear();if(value.isEmpty()){return;}int pos=0;
		while(pos<value.length()){long a=0,z=0;final int first=pos;while(pos<value.length() && Character.isDigit(value.charAt(pos))){a=a*10+value.charAt(pos++)-'0';require(a<=Integer.MAX_VALUE,"Interval coordinate overflow");}
			require(pos>first && pos<value.length() && value.charAt(pos++)=='-',"Expected inclusive start-stop interval");final int after=pos;
			while(pos<value.length() && Character.isDigit(value.charAt(pos))){z=z*10+value.charAt(pos++)-'0';require(z<=Integer.MAX_VALUE,"Interval coordinate overflow");}
			require(pos>after && a>=(oneBased?1:0) && z>=a,"Invalid paired interval bounds");list.add(pack((int)a,(int)z));if(pos<value.length()){require(value.charAt(pos++)==',' && pos<value.length(),"Intervals require comma separators without trailing empties");}}
	}
	private static long pack(int a,int z){require(a>=0 && z>=a,"Projected overlap must be a valid inclusive query interval");return ((long)a<<32)|(z&0xffffffffL);}
	private static int lo(long pair){return (int)(pair>>>32);}private static int hi(long pair){return (int)pair;}
	private static void appendIntervals(ByteBuilder b,LongList list){for(int i=0;i<list.size;i++){if(i>0){b.append(',');}b.append(lo(list.array[i])).append('-').append(hi(list.array[i]));}}
	private static boolean bool(LineParser1 p,int n){require(p.termEquals("true",n) || p.termEquals("false",n),"Explicit recorded boolean required");return p.termEquals("true",n);}
	private static <V> ArrayList<String> keys(ObjectMap<String,V> map){final ArrayList<String> out=new ArrayList<String>();for(String key:map.keys()){if(key!=null){out.add(key);}}Collections.sort(out);return out;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Call{
		Call(LineParser1 p,byte[] row){id=p.parseString(0);context=p.parseString(1);start=p.parseInt(5)-1;stop=p.parseInt(6)-1;final String s=p.parseString(7);require(s.equals("+") || s.equals("-"),"Call strand must be explicit");strand=s.charAt(0);known=p.parseLong(10);raw=row.clone();
			require(start>=0 && stop>=start && p.parseInt(8)==stop-start+1 && known>=0 && known<=stop-start+1,"Original call bounds and passing coverage denominator must agree");for(int i=0;i<support.length;i++){support[i]=new LongList();}}
		final String id,context;final int start,stop;final char strand;final long known;final byte[] raw;final LongList mapped=new LongList();final LongList[] support=new LongList[SUPPORT.length];final ObjectMap<String,Overlap> overlaps=new ObjectMap<String,Overlap>(String.class,Overlap.class);
	}
	private static final class Annotation{
		Annotation(LineParser1 p,long line_){line=line_;id=p.parseString(0);reference=p.parseString(1);status=p.parseString(8);category=category(status);product=p.parseString(15);named=namedSsu(product);geometry=p.parseString(23);strandText=p.parseString(20);require(strandText.equals("+") || strandText.equals("-") || strandText.equals("MIXED"),"Source annotation strand must remain explicit");strand=strandText.equals("MIXED")?'?':strandText.charAt(0);segments=NcrnaMitoCmEvidence.parseSegments(geometry,Integer.MAX_VALUE);
			long length=0;int min=Integer.MAX_VALUE,max=0;for(NcrnaMitoCmEvidence.Span s:segments){require(strand=='?' || s.strand==strand,"Retained source segment strand disagrees with annotation strand");length+=(long)s.hi-s.lo+1;min=Math.min(min,s.lo);max=Math.max(max,s.hi);}require(length==p.parseInt(25) && min==p.parseInt(21) && max==p.parseInt(22),"Retained source geometry must conserve the annotated sequence length and bounds");
			gbff=p.parseString(26);partial=p.parseString(17);pseudo=p.parseString(18);rawLocation=p.parseString(14);}
		final String id,reference,status,product,geometry,strandText,gbff,partial,pseudo,rawLocation;final int category;final long line;final char strand;final boolean named;final ArrayList<NcrnaMitoCmEvidence.Span> segments;
	}
	private static final class Overlap{Overlap(Annotation a){annotation=a;}final Annotation annotation;final LongList query=new LongList(),opposite=new LongList();long rows,primaryRows,firstLine=Long.MAX_VALUE,lastLine;}
	private static final String[] STATUSES={"SSU_SUPPORT_REVIEW","SSU_LSU_CONFLICT_REVIEW","LSU_SUPPORT_REVIEW","NO_SAME_STRAND_GA_HIT_UNRESOLVED"};
	private static final String[] CATEGORIES={"CM_SSU_ONLY","CM_SSU_LSU_CONFLICT","CM_LSU_ONLY","UNRESOLVED_ANNOTATION"};
	private static final String[] SUPPORT={"cm_ssu_only","cm_ssu_lsu_conflict","cm_lsu_only","unresolved_annotation","named_ssu","named_ssu_conflicting_cm"};
	private static final String KNOWN="|calls|alignments|loci|verified|outprefix|expected|";
	private final ObjectMap<String,String> flags=new ObjectMap<String,String>(String.class,String.class);private ObjectMap<String,String> pins;private int expected;private byte[] callHeader;
	private final ObjectMap<String,Call> calls=new ObjectMap<String,Call>(String.class,Call.class);private final ObjectMap<String,ArrayList<Annotation>> annotations=new ObjectMap(String.class,ArrayList.class);
	private final LongList query=new LongList(),reference=new LongList(),overlap=new LongList();private long annotationCount,alignmentRows,passingRows,rejectedRows;
}
