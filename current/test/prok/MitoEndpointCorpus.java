package test.prok;

import prok.RrnaResourceIO;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import tax.TaxTree;

/** Reconstruct every mitochondrial construction member from its original genome,
 * assign whole hosts to an 80/20 split, then emit a deterministic cost-pilot sample.
 * No model fitting, deduplication, locus cap or evaluation-FN selection occurs here.
 * Uses the ByteFile template's parsing/writer idioms and the native sequence Streamer.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointCorpus {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, MitoEndpointCorpus.class, false);
		try{final MitoEndpointCorpus x=new MitoEndpointCorpus();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	private void parse(String[] args){
		final Parser p=new Parser();
		for(String arg:args){
			final String[] kv=Parse.splitOnFirst(arg, '=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("ssu")){ssu=v;}else if(k.equals("insect")){insect=v;}else if(k.equals("metadata")){metadata=v;}
			else if(k.equals("enrollment")){enrollment=v;}else if(k.equals("insectinput")){insectInput=Parse.parseIntKMG(v);}
			else if(k.equals("genomes")){genomes=v;}else if(k.equals("outprefix")){prefix=v;}
			else if(k.equals("expectedssu")){expectedSsu=Parse.parseIntKMG(v);}else if(k.equals("expectedinsect")){expectedInsect=Parse.parseIntKMG(v);}
			else if(k.equals("sample")){sample=Parse.parseIntKMG(v);}else if(k.equals("flank")){flank=Parse.parseIntKMG(v);}
			else if(k.equals("maxspan")){maxSpan=Parse.parseIntKMG(v);}else if(k.equals("seed")){seed=Long.parseLong(v);}
			else if(!p.parse(arg, k, v)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		require(ssu!=null && insect!=null && metadata!=null && genomes!=null && prefix!=null, "Declare both construction arms, metadata, genomic source and output prefix");
		require(expectedSsu>0 && expectedInsect>0 && sample>=2 && sample<=expectedSsu+expectedInsect && flank>=0 && maxSpan>0 && seed>=0 && !p.append && p.maxReads<0,
			"Complete denominators, bounded sample/window, explicit seed and no read truncation are required");
		require(enrollment==null?insectInput<0:insectInput>=expectedInsect,"A complete saved enrollment may select the insect subset from an explicitly sized larger source");
		Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;
		ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;ReadWrite.ZIPLEVEL=9;ReadWrite.setZipThreads(1);
		RrnaResourceIO.fresh(prefix, SUFFIXES);
	}

	private void run(){
		if(enrollment!=null){loadEnrollment();}
		loadArm(ssu, "SSU", expectedSsu, expectedSsu);loadArm(insect, "INSECT", expectedInsect, enrollment==null?expectedInsect:insectInput);
		if(enrollment!=null){for(Enrollment e:enrolled.values()){if(e!=null){require(e.seen,"Missing declared construction member: "+e.id);}}}
		loadMetadata();
		splitHosts();select("SSU", (int)Math.round(sample*(expectedSsu/(double)(expectedSsu+expectedInsect))));
		select("INSECT", sample-selected);
		require(selected==sample, "Pilot selection must conserve its declared member budget");
		final ByteStreamWriter fasta=RrnaResourceIO.writer(prefix+SUFFIXES[0]);final ByteBuilder b=new ByteBuilder();
		try{
			RrnaResourceIO.stream(genomes, r->{
				final Members group=byAccession.get(RrnaResourceIO.first(r.id));if(group==null){return;}
				final float gc=Tools.calcGC(r.bases);
				for(Member m:group){
					require(!m.seen, "A source accession must occur once in the genomic FASTA: "+m.accession);m.seen=true;
					m.geometry.verify(r.bases, m.exact, m.id);m.genomeLength=r.length();m.gc=gc;verified++;
					final MitoEndpointGeometry.Window w=m.selected?m.geometry.window(r.bases, flank, maxSpan):null;
					if(!m.selected){m.disposition="NOT_SAMPLED";continue;}
					if(w==null){m.disposition="UNSUPPORTED_WINDOW";continue;}
					m.start=w.start;m.stop=w.stop;m.windowLength=w.bases.length;
					m.disposition=m.partial||m.pseudo?"UNCERTAIN_ENDS":"EXACT_SOURCE_ENDS";written++;
					// Explicit zero-based inclusive bounds; joined annotations retain
					// introns in this genomic window even though exact parity used splicing.
					fasta.print(b.clear().append('>').append(m.id).append(" assembly=").append(m.accession).append(" host=").append(m.host.id)
						.append(" split=").append(m.host.validation?"VALIDATION":"TRAIN").append(" arm=").append(m.arm)
						.append(" start0=").append(m.start).append(" stop0=").append(m.stop).append(" gccontig=").append(gc, 8)
						.append(" partial=").append(m.partial).append(" pseudo=").append(m.pseudo).append(" segments=").append(m.geometry.lo.length)
						.append(" crossesorigin=").append(m.geometry.crossesOrigin).nl().append(w.bases).nl());
				}
			});
		}finally{RrnaResourceIO.finish(fasta);}
		for(Member m:members){require(m.seen, "Missing original genomic source for "+m.id);}
		require(verified==expectedSsu+expectedInsect, "Every construction copy must pass exact genomic reconstruction before pilot use");
		writeLedgers();System.err.println("MITO_ENDPOINT_CORPUS_PASS members="+verified+" hosts="+hosts.size()+" selected="+selected+" windows="+written+" output_files=4 no_training=true");
	}

	/** Read the accepted full coverage ledger, never derive enrollment from recall. */
	private void loadEnrollment(){
		final ByteFile in=RrnaResourceIO.open(enrollment);final LineParser1 p=new LineParser1('\t');int ssuCount=0,insectCount=0,panelCount=0;
		try{RrnaResourceIO.header(in,ENROLLMENT_HEADER,enrollment);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==19,"Accepted coverage-member schema changed");
				if(p.termEquals("PANEL_INCLUDED",0)){panelCount++;continue;}
				require(p.termEquals("CONSTRUCTION",0),"Only explicitly declared construction rows enter the endpoint corpus");
				final String arm=p.termEquals("ssu",7)?"SSU":p.termEquals("nohit",7)?"INSECT":null;
				require(arm!=null,"Construction provenance must identify SSU or nohit source");
				final Enrollment e=new Enrollment(p.parseString(1),arm,p.parseInt(3),flag(p.parseString(5)),flag(p.parseString(6)));
				require(e.taxid>0 && enrolled.put(e.id,e)==null,"Full construction enrollment must retain unique source IDs and host taxids");
				if(arm.equals("SSU")){ssuCount++;}else{insectCount++;}
			}
		}finally{require(!in.close(),"Accepted construction roster read failed");}
		require(ssuCount==expectedSsu && insectCount==expectedInsect && panelCount==22,"The accepted coverage ledger must conserve11420+4690construction copies and exclude its22panel rows");
	}

	private void loadArm(String path, String arm, int expected, int inputExpected){
		assert(path!=null && expected>0):"Both accepted construction arms need complete populations";
		final int before=members.size();final ObjectSet<String> allIds=new ObjectSet<String>(String.class);
		RrnaResourceIO.stream(path, r->{
			final Member m=new Member();m.id=RrnaResourceIO.first(r.id);m.arm=arm;m.exact=r.bases;
			require(allIds.add(m.id),"The complete original source cannot repeat a locus: "+m.id);
			final Enrollment selected=enrollment==null?null:enrolled.get(m.id);
			if(enrollment!=null && selected==null){return;}
			m.partial=flag(RrnaResourceIO.field(r.id, "partial"));m.pseudo=flag(RrnaResourceIO.field(r.id, "pseudo"));m.taxid=TaxTree.parseTaxID(r.id);
			if(selected!=null){require(!selected.seen && selected.arm.equals(arm) && selected.taxid==m.taxid && selected.partial==m.partial && selected.pseudo==m.pseudo,
				"Saved construction enrollment must match the original source/flags: "+m.id);selected.seen=true;}
			require(m.exact!=null && m.exact.length>0 && m.taxid>0 && byId.put(m.id, m)==null, "Construction IDs must be unique across both arms: "+m.id);
			members.add(m);
		});require(members.size()-before==expected && allIds.size()==inputExpected, "Construction/source arm count changed: "+arm);
	}

	private void loadMetadata(){
		assert(!members.isEmpty()):"Source metadata joins the frozen full construction population";
		final ByteFile in=RrnaResourceIO.open(metadata);final LineParser1 p=new LineParser1('\t');int joined=0;
		try{
			final byte[] header=in.nextLine();require(header!=null, "Missing source-evidence schema");p.set(header);
			final ObjectIntMap<String> columns=new ObjectIntMap<String>(String.class);
			for(int i=0; i<p.terms(); i++){final String name=p.parseString(i);require(columns.get(name)<0, "Duplicate source-evidence column");columns.put(name,i);}
			final int id=column(columns,"locus_id"),sid=column(columns,"source_locus_id"),acc=column(columns,"source_accession"),tid=column(columns,"source_taxid"),
				part=column(columns,"source_partial"),pseudo=column(columns,"source_pseudo"),status=column(columns,"source_status"),
				strand=column(columns,"source_strand"),segments=column(columns,"source_segments_in_transcript_order"),cross=column(columns,"source_crosses_origin"),length=column(columns,"source_extracted_length");
			final int width=p.terms();
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==width, "Source metadata row width changed");final Member m=byId.get(p.parseString(id));if(m==null){continue;}
				require(m.accession==null && p.termEquals(m.id,sid) && p.termEquals("EXTRACTED",status) && p.parseInt(length)==m.exact.length && p.parseInt(tid)==m.taxid,
					"Source metadata must bind each exact member once with matching length/taxid: "+m.id);
				require(flag(p.parseString(part))==m.partial && flag(p.parseString(pseudo))==m.pseudo, "Source/member flags disagree: "+m.id);
				m.accession=p.parseString(acc);m.strand=p.parseString(strand);m.segmentText=p.parseString(segments);
				m.geometry=new MitoEndpointGeometry(m.segmentText, m.strand, flag(p.parseString(cross)));
				final String host=p.parseString(tid);m.host=hostMap.get(host);
				if(m.host==null){m.host=new Host(host);hostMap.put(host,m.host);hosts.add(m.host);}
				Members group=byAccession.get(m.accession);if(group==null){group=new Members();byAccession.put(m.accession,group);}
				else{require(group.get(0).host==m.host, "A source genome cannot cross host split identities");}group.add(m);joined++;
			}
		}finally{require(!in.close(), "Source metadata read failed");}
		require(joined==members.size(), "Every construction member requires original coordinate/host metadata");
	}

	private void splitHosts(){
		require(hosts.size()>=5, "80/20 source-host split needs at least five hosts");
		Collections.sort(hosts, (a,b)->compareKey(a.id,b.id,seed));
		final int n=hosts.size()/5;for(int i=0; i<n; i++){hosts.get(i).validation=true;}
	}
	private void select(String arm, int count){
		require(count>0, "The cost pilot must represent both construction arms");
		final ArrayList<Member> pool=new ArrayList<Member>();for(Member m:members){if(m.arm.equals(arm)){pool.add(m);}}
		require(count<=pool.size(), "Sample exceeds its declared construction arm");
		Collections.sort(pool, (a,b)->compareKey(a.id,b.id,seed^0x6a09e667f3bcc909L));
		for(int i=0; i<count; i++){pool.get(i).selected=true;selected++;}
	}
	private void writeLedgers(){
		assert(verified==members.size()):"No complete cohort receipt before every genomic parity check passes";
		final ByteStreamWriter roster=RrnaResourceIO.writer(prefix+SUFFIXES[1]),pilot=RrnaResourceIO.writer(prefix+SUFFIXES[2]),summary=RrnaResourceIO.writer(prefix+SUFFIXES[3]);
		final ByteBuilder b=new ByteBuilder();
		try{
			roster.println(HEADER);pilot.println(HEADER);
			int partial=0,pseudo=0,joined=0,circular=0,validation=0,eligible=0;
			for(Member m:members){
				partial+=m.partial?1:0;pseudo+=m.pseudo?1:0;joined+=m.geometry.lo.length>1?1:0;circular+=m.geometry.crossesOrigin?1:0;validation+=m.host.validation?1:0;
				eligible+=m.disposition.equals("EXACT_SOURCE_ENDS")?1:0;
				b.clear().append(m.id).tab().append(m.arm).tab().append(m.accession).tab().append(m.host.id).tab().append(m.host.validation?"VALIDATION":"TRAIN")
					.tab().append(m.exact.length).tab().append(m.genomeLength).tab().append(m.strand).tab().append(m.segmentText).tab().append(m.partial).tab().append(m.pseudo)
					.tab().append(m.geometry.crossesOrigin).tab().append(m.selected).tab().append(m.disposition).tab().append(m.start).tab().append(m.stop).tab().append(m.windowLength).tab().append(m.gc,8).nl();
				roster.print(b);if(m.selected){pilot.print(b);}
			}
			summary.println("members\thosts\tvalidation_hosts\tvalidation_members\tpartial_members\tpseudo_members\tjoined_members\torigin_crossing_members\tselected\twindows\tpilot_exact_end_windows\tflank\tmax_span\tseed\tclaim");
			summary.print(b.clear().append(members.size()).tab().append(hosts.size()).tab().append(hosts.size()/5).tab().append(validation).tab().append(partial).tab().append(pseudo)
				.tab().append(joined).tab().append(circular).tab().append(selected).tab().append(written).tab().append(eligible).tab().append(flank).tab().append(maxSpan).tab().append(seed).append("\tINPUT_PARITY_AND_COST_PILOT_ONLY").nl());
		}finally{RrnaResourceIO.finish(roster,pilot,summary);}
	}
	private static int column(ObjectIntMap<String> columns, String name){final int i=columns.get(name);require(i>=0,"Missing source-evidence column: "+name);return i;}
	private static boolean flag(String s){require(s.equals("true") || s.equals("false"),"Expected explicit source boolean: "+s);return s.equals("true");}
	/** Stable sampling mixer, not a provenance checksum. Ties use the source ID. */
	public static int compareKey(String a, String b, long seed){final int c=Long.compareUnsigned(sampleKey(a,seed),sampleKey(b,seed));return c==0?a.compareTo(b):c;}
	private static long sampleKey(String s, long seed){
		assert(s!=null && !s.isEmpty()):"Deterministic sampling must key a stable source identity";
		long x=seed;for(int i=0; i<s.length(); i++){x=(x^s.charAt(i))*0x100000001b3L;}x=(x^(x>>>30))*0xbf58476d1ce4e5b9L;x=(x^(x>>>27))*0x94d049bb133111ebL;return x^(x>>>31);
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Host{Host(String s){id=s;}final String id;boolean validation;}
	private static final class Enrollment{
		Enrollment(String i,String a,int t,boolean p,boolean ps){id=i;arm=a;taxid=t;partial=p;pseudo=ps;}
		final String id,arm;final int taxid;final boolean partial,pseudo;boolean seen;
	}
	private static final class Member{
		String id,arm,accession,strand,segmentText,disposition;byte[] exact;Host host;MitoEndpointGeometry geometry;
		int taxid,genomeLength,start=-1,stop=-1,windowLength;float gc;boolean partial,pseudo,selected,seen;
	}
	private static final class Members extends ArrayList<Member>{private static final long serialVersionUID=1L;}
	private String ssu,insect,metadata,genomes,prefix,enrollment;private int expectedSsu,expectedInsect,insectInput=-1,sample,flank=-1,maxSpan=-1,selected,written,verified;private long seed=-1;
	private final ArrayList<Member> members=new ArrayList<Member>();private final ArrayList<Host> hosts=new ArrayList<Host>();
	private final ObjectMap<String,Member> byId=new ObjectMap<String,Member>(String.class,Member.class);
	private final ObjectMap<String,Host> hostMap=new ObjectMap<String,Host>(String.class,Host.class);
	private final ObjectMap<String,Members> byAccession=new ObjectMap<String,Members>(String.class,Members.class);
	private final ObjectMap<String,Enrollment> enrolled=new ObjectMap<String,Enrollment>(String.class,Enrollment.class);
	private static final String[] SUFFIXES={".fa.gz",".members.tsv.gz",".pilot.tsv.gz",".summary.tsv.gz"};
	public static final String HEADER="locus_id\tarm\taccession\thost_taxid\tsplit\texact_length\tgenome_length\tstrand\tsegments\tpartial\tpseudo\tcrosses_origin\tselected\tdisposition\tstart0\tstop0\twindow_length\tgc";
	private static final String ENROLLMENT_HEADER="role\tlocus_id\tunique_id\ttaxid\tphylum\tpartial\tpseudo\tstratum\tsource_model\tassigned_cover_model\tbest_selected_model\tselected_identity\tselected_trace_identity\tmember_paired_fraction\tmodel_paired_fraction\tbest_pool_model\tbest_pool_identity\tstatus\tprovenance";
}
