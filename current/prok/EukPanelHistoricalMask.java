package prok;

import java.util.ArrayList;
import java.util.Arrays;
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

/** Conserves historical non-mRNA annotations and marks eligible unscored masks.
 * No source edits or truth admission. A complete raw link to a surviving admitted
 * seed is evidence of actual fixed-locus overlap; a member-only extension is not.
 * @author Ganyu */
public final class EukPanelHistoricalMask {
	public static final String TRUTH_HEADER=EukPanelTruthDisposition.OUTPUT_HEADER;
	public static final String REVIEW_HEADER=EukPanelCmProposal.SOURCE_HEADER+"\tsource_srp_subtype\tproposal_status";
	public static final String LINK_HEADER=EukPanelTruthCandidate.LINK_HEADER;
	public static final String OUTPUT_HEADER="assembly\tnative_seqid\tstrand\tstart\tend\tfold_family\thistory_class\tmask_eligible\tsource_gff\tsource_sha80\tsource_line\toriginal_type\toriginal_attributes\tmask_reason";
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,EukPanelHistoricalMask.class,false);
		String truth=null,review=null,links=null,models=null,sources=null,out=null;int nt=-1,nm=-1;long nr=-1,nl=-1;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("truth")){truth=v;}else if(k.equals("review")){review=v;}else if(k.equals("links")){links=v;}else if(k.equals("models")){models=v;}else if(k.equals("sources")){sources=v;}
			else if(k.equals("outprefix")){out=v;}else if(k.equals("expectedloci")){nt=Integer.parseInt(v);}else if(k.equals("expectedmodels")){nm=Integer.parseInt(v);}
			else if(k.equals("expectedreview")){nr=Long.parseLong(v);}else if(k.equals("expectedlinks")){nl=Long.parseLong(v);}else{throw new IllegalArgumentException("Unknown flag: "+k);}}
		require(truth!=null && review!=null && links!=null && models!=null && sources!=null && out!=null && nt>0 && nm>0 && nr>=0 && nl>=0,"Frozen truth, original evidence, model map, source bindings and denominators required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);final EukPanelHistoricalMask x=new EukPanelHistoricalMask();
		x.hashes=EukPanelGffValidation.hashes(sources);x.loadModels(models,nm);x.loadTruth(truth,nt);x.loadReview(review,nr);x.loadLinks(links,nl);x.write(out);Shared.closeStream(pp.outstream);
	}
	void loadModels(String path,int expected){
		assert(classifier.modelFamilies.isEmpty()):"Only the pinned model map may assign RFAM folds";
		final EukPanelLoci roster=new EukPanelLoci();roster.loadModels(path,expected);
		for(EukPanelLoci.Model m:roster.models.values()){if(m!=null){require(NcrnaClippingComparison.known(m.family),"Unknown mapped fold");
			final String old=classifier.modelFamilies.put(m.accession,m.family);require(old==null || old.equals(m.family),"Conflicting accession-to-fold map");}}
	}
	void loadTruth(String path,int expected){
		assert(seeds.isEmpty()):"Only one frozen surviving-locus universe may govern historical masks";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,TRUTH_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==35 && p.termEquals("TRUTH_V1",32),"Historical masks require frozen35-column truth, not an unfrozen candidate");
			final Seed s=new Seed(p);require(ids.add(p.parseString(0)) && seeds.put(p.parseString(1)+'\t'+p.parseString(10)+'\t'+p.parseLong(12),s)==null,"Duplicate frozen locus or seed identity");truthRows++;
		}}finally{require(!in.close(),"Frozen truth read failed");}require(truthRows==expected,"Frozen locus denominator changed");
	}
	void loadReview(String path,long expected){
		assert(reviewRows==0):"The original review file is consumed once without editing its annotations";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,REVIEW_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==28 && p.termEquals("UNRESOLVED_SOURCE_REVIEW_ONLY",27),"Original review schema/status changed");source(p);reviewRows++;
		}}finally{require(!in.close(),"Historical review read failed");}require(reviewRows==expected,"Historical review denominator changed");
	}
	Source source(LineParser1 p){
		assert(p.terms()>=26):"Historical source identity and support come from the original26 columns";
		final String type=p.parseString(16);if(!type.endsWith("RNA") || type.equals("mRNA")){return null;}
		final String assembly=p.parseString(1),gff=p.parseString(2),hash=p.parseString(3);final long line=p.parseLong(4);
		require(line>0 && hash.equals(hashes.get(gff)),"Historical RNA source lacks its exact file/line/hash binding");
		final String key=assembly+'\t'+gff+'\t'+line;final Source old=bySource.get(key);p.setBounds(25);final int length=p.b();
		if(old!=null){require(old.original.length==length,"Repeated source metadata length differs");for(int i=0;i<length;i++){require(old.original[i]==p.line()[i],"Conflicting original annotation for one assembly/GFF/line");}return old;}
		final Source s=new Source(p,key,Arrays.copyOf(p.line(),length));bySource.put(key,s);rows.add(s);
		try{s.family=classifier.sourceFamily(s.type,s.attrs);}catch(IllegalArgumentException e){s.family="UNKNOWN";s.malformedFamily=true;}
		return s;
	}
	void loadLinks(String path,long expected){
		assert(linkRows==0):"Physical links are consumed once; duplicate links cannot inflate support";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,LINK_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==38 && p.termEquals("PHYSICAL_LINK_REVIEW_ONLY",37),"Original link schema/status changed");linkRows++;
			final Source s=source(p);if(s==null){continue;}
			final long hitLine=p.parseLong(28);final String key=s.assembly+'\t'+p.parseString(27)+'\t'+hitLine;
			require(hitLine>0 && seen.add(s.key+'\t'+key),"Duplicate or invalid historical RNA/hit link");
			final long a=Math.min(p.parseLong(31),p.parseLong(32)),b=Math.max(p.parseLong(31),p.parseLong(32));
			final long ov=overlap(s.start,s.end,a,b);
			require(s.valid && s.stranded && p.termEquals(s.strand,33) && a>0 && b>=a && ov>0 && ov==p.parseLong(34),"Raw RNA link must have valid same-strand positive overlap");
			s.observedLinks++;final Seed seed=seeds.get(key);
			if(seed!=null){
				require(seed.seqid.equals(s.seqid) && seed.strand.equals(s.strand) && seed.start==a && seed.end==b && p.termEquals(seed.model,29) && p.termEquals(seed.family,30),"Surviving seed link identity/geometry changed");
				if(seed.admitted){s.admittedLink=true;}
			}
		}}finally{require(!in.close(),"Historical links read failed");}
		require(linkRows==expected,"Physical-link denominator changed");
		for(Source s:rows){require(s.observedLinks==s.expectedLinks,"Historical RNA links must be complete for seed-based admitted-overlap classification: "+s.key);}
	}
	void write(String prefix){
		assert(rows.size()==bySource.size()):"Each original non-mRNA annotation must survive exactly once";
		RrnaResourceIO.fresh(prefix,".historical_masks.tsv.gz",".historical_counts.tsv",".historical_summary.tsv");
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".historical_masks.tsv.gz"),counts=RrnaResourceIO.writer(prefix+".historical_counts.tsv"),summary=RrnaResourceIO.writer(prefix+".historical_summary.tsv");
		final ByteBuilder b=new ByteBuilder(),reason=new ByteBuilder();final ObjectMap<String,long[]> bins=new ObjectMap<String,long[]>(String.class,long[].class);long admitted=0,linked=0,only=0,eligible=0;
		try{out.println(OUTPUT_HEADER);
			for(Source s:rows){
				final String category=s.admittedLink?"ADMITTED_LINK":s.observedLinks>0?"HIST_LINKED_UNSCORED":"HIST_ONLY";
				final boolean mask=!s.admittedLink && s.valid && s.stranded && NcrnaClippingComparison.known(s.family);
				if(s.admittedLink){admitted++;}else if(s.observedLinks>0){linked++;}else{only++;}if(mask){eligible++;}
				reason.clear();
				if(s.admittedLink){reason.append("ADMITTED_LOCUS_OVERLAP");}
				if(!s.valid){reason(reason,"INVALID_INTERVAL:"+s.validation);}
				if(!s.stranded){reason(reason,"UNRESOLVED_STRAND");}
				if(!NcrnaClippingComparison.known(s.family)){reason(reason,s.malformedFamily?"MALFORMED_FAMILY_ANNOTATION":"UNRESOLVED_FAMILY");}
				if(reason.length()==0){reason.append(s.observedLinks>0?"RAW_LINKS_WITHOUT_ADMITTED_SEED_OVERLAP":"NO_RAW_HIT_LINK");}
				out.print(b.clear().append(s.assembly).tab().append(s.seqid).tab().append(s.strand).tab().append(s.startText).tab().append(s.endText).tab().append(s.family).tab().append(category).tab().append(mask?1:0)
					.tab().append(s.gff).tab().append(s.hash).tab().append(s.line).tab().append(s.type).tab().append(s.attrs).tab().append(reason).nl());
				EukPanelCmProposal.increment(bins,s.family+'\t'+category+'\t'+(mask?1:0),1);
			}
			counts.println("fold_family\thistory_class\tmask_eligible\tannotations");long total=0;
			for(String k:EukPanelCmProposal.sorted(bins)){final long n=bins.get(k)[0];counts.print(b.clear().append(k).tab().append(n).nl());total+=n;}
			require(total==rows.size() && admitted+linked+only==total,"Historical classifications must conserve all unique RNA annotations");
			summary.println("truth_loci\tsource_review_rows\tsource_link_rows\tunique_historical_rna\tadmitted_link\thist_linked_unscored\thist_only\tmask_eligible\tstatus");
			summary.print(b.clear().append(truthRows).tab().append(reviewRows).tab().append(linkRows).tab().append(total).tab().append(admitted).tab().append(linked).tab().append(only).tab().append(eligible).tab().append("HISTORICAL_MASKS_ONLY_NO_TRUTH_ADMISSION").nl());
		}finally{RrnaResourceIO.finish(out,counts,summary);}
		System.err.println("HISTORICAL_MASK_COMPLETE truth_loci="+truthRows+" review_rows="+reviewRows+" link_rows="+linkRows+" historical_rna="+rows.size()+" admitted_link="+admitted+" linked_unscored="+linked+" hist_only="+only+" eligible="+eligible+" no_truth_change=true");
	}
	static void reason(ByteBuilder b,String text){if(b.length()>0){b.append(';');}b.append(text);}
	static long overlap(long a,long b,long c,long d){return Math.max(0,Math.min(b,d)-Math.max(a,c)+1);}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Seed{
		Seed(LineParser1 p){seqid=p.parseString(2);strand=p.parseString(3);start=p.parseLong(4);end=p.parseLong(5);model=p.parseString(9);family=p.parseString(6);admitted=p.termEquals("ADMITTED",26);
			require(start>0 && end>=start && (strand.equals("+") || strand.equals("-")) && (admitted || p.termEquals("BELOW_GA",26)),"Invalid frozen seed geometry/admission");
			require(admitted==EukPanelLocusDiagnostics.bit(p,14),"Frozen admission must retain any-GA policy");}
		final String seqid,strand,model,family;final long start,end;final boolean admitted;
	}
	static final class Source{
		Source(LineParser1 p,String k,byte[] raw){key=k;original=raw;assembly=p.parseString(1);gff=p.parseString(2);hash=p.parseString(3);line=p.parseLong(4);validation=p.parseString(5);seqid=p.parseString(8);
			type=p.parseString(16);startText=p.parseString(17);endText=p.parseString(18);strand=p.parseString(20);attrs=p.parseString(22);expectedLinks=p.parseLong(24);
			valid=validation.equals("VALID_INTERVAL");stranded=strand.equals("+") || strand.equals("-");
			if(valid){start=p.parseLong(17);end=p.parseLong(18);require(start>0 && end>=start && end<=p.parseLong(9) && !seqid.equals("NA"),"Valid source status contradicts original coordinates");}
			else{start=end=-1;require(validation.equals("OUT_OF_RANGE") || validation.equals("UNRESOLVED_SEQID") || validation.equals("INVALID_COORDINATES"),"Unknown historical interval status");}
			require(expectedLinks>=0,"Historical link count cannot be negative");}
		final String key,assembly,gff,hash,validation,seqid,type,startText,endText,strand,attrs;final long line,start,end,expectedLinks;final byte[] original;final boolean valid,stranded;
		String family;long observedLinks;boolean admittedLink,malformedFamily;
	}
	final EukPanelTruthCandidate classifier=new EukPanelTruthCandidate();
	final ObjectMap<String,Seed> seeds=new ObjectMap<String,Seed>(String.class,Seed.class);
	final ObjectMap<String,Source> bySource=new ObjectMap<String,Source>(String.class,Source.class);
	final ArrayList<Source> rows=new ArrayList<Source>();ObjectMap<String,String> hashes;long truthRows,reviewRows,linkRows;
}
