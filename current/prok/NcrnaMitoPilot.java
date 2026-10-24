package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Deterministic, annotation-independent-of-score mitochondrial CM pilot targets.
 * Keeps the complete census ledgers and every selected rRNA label. Missing annotations
 * are unknown, not negatives. Circular targets append at most one partial genome copy.
 * @author Ganyu */
public final class NcrnaMitoPilot {

	public static void main(String[] args){
		final NcrnaMitoPilot p=new NcrnaMitoPilot(args);
		try{p.process();}finally{Shared.closeStream(p.log);}
	}
	private NcrnaMitoPilot(String[] args){
		final PreParser pp=new PreParser(args,getClass(),false);log=pp.outstream;
		final Parser parser=new Parser();
		for(String arg:pp.args){
			final String[] kv=Parse.splitOnFirst(arg,'=');final String key=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(key.equals("records")){records=v;}
			else if(key.equals("loci")){loci=v;}
			else if(key.equals("genomes")){genomes=v;}
			else if(key.equals("candidates")){candidates=v;}
			else if(key.equals("outprefix")){prefix=v;}
			else if(key.equals("expectedrecords")){expectedRecords=Parse.parseIntKMG(v);}
			else if(key.equals("expectedloci")){expectedLoci=Parse.parseIntKMG(v);}
			else if(key.equals("expectedkinetoplast")){expectedKinetoplast=Parse.parseIntKMG(v);}
			else if(key.equals("pergroup")){perGroup=Parse.parseIntKMG(v);}
			else if(key.equals("maxbases")){maxBases=Parse.parseKMG(v);}
			else if(key.equals("window")){window=Parse.parseIntKMG(v);}
			else if(key.equals("selection")){require("quartile".equals(v) || "longtail".equals(v),"selection must be quartile or longtail");longTail="longtail".equals(v);}
			else if(!parser.parse(arg,key,v)){throw new IllegalArgumentException("Unknown option: "+arg);}
		}
		require(records!=null && loci!=null && genomes!=null && candidates!=null && prefix!=null,"Require records= loci= genomes= candidates= outprefix=");
		require(expectedRecords>0 && expectedLoci>0 && perGroup>0 && maxBases>0 && window>0,"Require positive expectedrecords/expectedloci/pergroup/maxbases/window; window is largest pinned CM W");
		require(expectedKinetoplast>=0,"Require expectedkinetoplast=N>=0 from the frozen source roster; absence must not silently defeat a declared rare stratum");
		require(!parser.append && parser.maxReads<0,"Pilot selection cannot append or truncate the source census");
		Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;ReadWrite.setZipThreads(Math.min(4,Shared.threads()));
	}
	private void process(){
		for(String suffix:SUFFIXES){require(!new File(prefix+suffix).exists(),"Refusing existing pilot output: "+prefix+suffix);}
		readRecords();readLoci();select();writeTables();
		emitGenomes();emitCandidates();writePopulation();
		final ByteStreamWriter out=table(prefix+".summary.tsv");final ByteBuilder b=new ByteBuilder();
		try{
			out.println("selected_genomes\tannotated\tno_rrna\toriginal_bases\tsearch_bases\tcandidate_loci\tcandidate_bases\tcanary_genomes\tcanary_search_bases\tcanary_candidate_loci\tcanary_candidate_bases\twindow");
			out.print(b.append(selectedGenomes).tab().append(countGroup(0)).tab().append(countGroup(1)).tab().append(selectedBases).tab().append(searchBases)
				.tab().append(selectedLoci).tab().append(candidateBases).tab().append(canaryGenomes).tab().append(canarySearchBases)
				.tab().append(canaryLoci).tab().append(canaryCandidateBases).tab().append(window).nl());
		}finally{require(!out.poisonAndWait(),"Pilot summary write failed");}
		log.println("MITO_PILOT_PASS genomes="+selectedGenomes+" bases="+selectedBases+" loci="+selectedLoci+" canary_genomes="+canaryGenomes);
	}
	private void readRecords(){
		final ByteFile in=ByteFile.makeByteFile(records,true);final LineParser1 p=new LineParser1('\t');
		try{
			recordHeader=in.nextLine();require(recordHeader!=null,"Missing records header");p.set(recordHeader);
			require(p.terms()==14 && p.termEquals("accession",0) && p.termEquals("sequence_source",11),"Expected native 14-column census records schema");
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==14,"Malformed census record");
				final Entry e=new Entry();e.id=p.parseString(0);e.taxid=p.parseInt(1);e.organelle=p.parseString(2);e.length=p.parseInt(5);e.features=p.parseInt(7);
				require(e.taxid>0 && e.length>0 && e.features>=0 && p.parseInt(8)==e.features && p.parseInt(9)==0,"Pilot requires a fully extracted census; unsupported loci need explicit new handling: "+e.id);
				require(p.termEquals("circular",6) || p.termEquals("linear_or_unspecified",6),"Unknown topology: "+e.id);
				require(p.termEquals("FNA_CONTIG",11) || p.termEquals("GBFF_ORIGIN",11),"Unknown sequence provenance: "+e.id);
				e.circular=p.termEquals("circular",6);e.fallback=p.termEquals("FNA_CONTIG",11);e.raw=line.clone();
				populationBases+=e.length;populationSearchBases+=(long)e.length+padding(e.length,window,e.circular);maxGenomeLength=Math.max(maxGenomeLength,e.length);
				require(byId.get(e.id)==null,"Duplicate census accession: "+e.id);byId.put(e.id,e);entries.add(e);
			}
		}finally{require(!in.close(),"Census records read failed");}
		require(entries.size()==expectedRecords,"Census record count mismatch");
	}
	private void readLoci(){
		final ByteFile in=ByteFile.makeByteFile(loci,true);final LineParser1 p=new LineParser1('\t');
		try{
			locusHeader=in.nextLine();require(locusHeader!=null,"Missing loci header");p.set(locusHeader);
			require(p.terms()==18 && p.termEquals("locus_id",0) && p.termEquals("extracted_length",16),"Expected native 18-column census loci schema");
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==18 && p.termEquals("EXTRACTED",7),"Unextracted/malformed locus requires explicit handling");
				final Locus l=new Locus();l.id=p.parseString(0);l.entry=byId.get(p.parseString(1));l.length=p.parseInt(16);l.raw=line.clone();
				require(l.entry!=null && p.parseInt(2)==l.entry.taxid && p.termEquals(l.entry.organelle,3) && l.length>0,"Locus provenance/length mismatch: "+l.id);
				l.entry.partial|=bool(p,8);l.entry.pseudo|=bool(p,9);l.entry.loci++;
				populationCandidateBases+=l.length;maxCandidateLength=Math.max(maxCandidateLength,l.length);
				if(longestLocus==null || l.length>longestLocus.length || l.length==longestLocus.length && l.id.compareTo(longestLocus.id)<0){longestLocus=l;}
				require(byLocus.get(l.id)==null,"Duplicate locus: "+l.id);byLocus.put(l.id,l);allLoci.add(l);
			}
		}finally{require(!in.close(),"Census loci read failed");}
		require(allLoci.size()==expectedLoci,"Census locus count mismatch");
		for(Entry e:entries){require(e.features==e.loci,"Records and loci must conserve per-accession annotation counts: "+e.id);}
	}
	/** Smallest sequence in each length quartile, accession tie break. Mandatory strata
	 * are selected first and count toward the same group quotas and quartile balance.
	 * No identity/CM outcomes are inputs. A mandatory or quota budget conflict fails loud. */
	private void select(){
		entries.sort(Comparator.comparingInt((Entry e)->e.length).thenComparing(e->e.id));
		for(int group=0;group<2;group++){
			final ArrayList<Entry> pool=new ArrayList<Entry>();for(Entry e:entries){if(group(e)==group){pool.add(e);}}
			require(pool.size()>=(longTail?1:perGroup),"Insufficient records for pilot group "+group);
			for(int i=0;i<pool.size();i++){pool.get(i).quartile=(int)(4L*i/pool.size());}
			if(longTail){Entry longest=null;for(Entry e:pool){if(longest==null || e.length>longest.length){longest=e;}}choose(longest);}
		}
		for(int s=0;s<STRATA.length;s++){
			Entry first=null;for(Entry e:entries){if(special(e,s)){available[s]++;if(first==null || (longTail?e.length>first.length:!first.selected && e.selected)){first=e;}}}
			if(s==0){require(available[s]==expectedKinetoplast,"Kinetoplast source count differs from expectedkinetoplast: "+available[s]+" versus "+expectedKinetoplast);}
			if(first!=null){choose(first);first.canary=true;}
		}
		if(longTail){
			// Actual longest group/stratum members, never shortened or skipped to fit a budget.
			// At most two groups, four special strata and the longest candidate's genome.
			// Preserve the candidate long tail even when it lies on a shorter genome.
			require(longestLocus!=null,"A positive expectedloci census must provide a longest candidate");choose(longestLocus.entry);
			for(Entry e:entries){if(e.selected){e.canary=true;canaryGenomes++;}}
			require(selectedGenomes>0 && selectedGenomes<=7,"Longtail selection is bounded by two groups, four strata and longest candidate");
			return;
		}
		for(int group=0;group<2;group++){
			while(countGroup(group)<perGroup){
				Entry choice=null;int bestCount=Integer.MAX_VALUE;
				for(Entry e:entries){if(group(e)!=group || e.selected || selectedBases+e.length>maxBases){continue;}
					final int qCount=quartileCounts[group][e.quartile];
					if(qCount<bestCount){choice=e;bestCount=qCount;}
				}
				require(choice!=null,"Budget prevents filling pilot group "+group+"; selected bases="+selectedBases+" max="+maxBases);
				choose(choice);
			}
			for(Entry e:entries){if(e.selected && group(e)==group){e.canary=true;break;}}
		}
		for(Entry e:entries){if(e.canary){canaryGenomes++;}}
		assert(selectedBases<=maxBases):"Every selection is charged once against original-genome budget";
	}
	private void choose(Entry e){
		if(e.selected){return;}
		require(longTail || countGroup(group(e))<perGroup,"Mandatory strata exceed annotation-group quota: "+e.id);
		require(selectedBases+e.length<=maxBases,"Mandatory pilot member exceeds remaining original-base budget: "+e.id);
		e.selected=true;selectedGenomes++;selectedBases+=e.length;quartileCounts[group(e)][e.quartile]++;
	}
	private int countGroup(int group){int n=0;for(Entry e:entries){if(e.selected && group(e)==group){n++;}}return n;}
	private static int group(Entry e){return e.features>0?0:1;}
	private static boolean special(Entry e,int s){
		// RefSeq NC_000894.1 uses the complete qualifier mitochondrion:kinetoplast.
		// Preserve that raw class in outputs; a bare-label-only match loses the rare stratum.
		return s==0?(e.organelle.equals("kinetoplast") || e.organelle.endsWith(":kinetoplast")):s==1?e.partial:s==2?e.pseudo:e.fallback;
	}
	private void writeTables(){
		final ByteStreamWriter selected=table(prefix+".selection.tsv"),rs=table(prefix+".records.tsv"),ls=table(prefix+".loci.tsv.gz"),strata=table(prefix+".strata.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			selected.println("accession\ttaxid\torganelle_class\tannotation_group\tlength_quartile0\toriginal_length\tprefix_added\tsearch_length\tpartial\tpseudo\tcontig_fallback\tcanary");
			rs.println(recordHeader);ls.println(locusHeader);
			for(Entry e:entries){if(e.selected){final int pad=padding(e.length,window,e.circular);
				selected.print(b.clear().append(e.id).tab().append(e.taxid).tab().append(e.organelle).tab().append(e.features>0?"annotated":"no_rrna")
					.tab().append(e.quartile).tab().append(e.length).tab().append(pad).tab().append((long)e.length+pad).tab().append(e.partial).tab().append(e.pseudo).tab().append(e.fallback).tab().append(e.canary).nl());
				rs.println(e.raw);
			}}
			for(Locus l:allLoci){if(l.entry.selected){ls.println(l.raw);selectedLoci++;if(l.entry.canary){canaryLoci++;}}}
			strata.println("stratum\tavailable_genomes\tselected_genomes\tstatus");
			for(int s=0;s<STRATA.length;s++){int n=0;for(Entry e:entries){if(e.selected && special(e,s)){n++;}}
				require(available[s]==0 || n>0,"An available mandatory stratum cannot silently disappear: "+STRATA[s]);
				strata.print(b.clear().append(STRATA[s]).tab().append(available[s]).tab().append(n).tab().append(n>0?"PRESENT":"ABSENT_IN_SOURCE").nl());
			}
		}finally{require(!(selected.poisonAndWait()|rs.poisonAndWait()|ls.poisonAndWait()|strata.poisonAndWait()),"Pilot tables write failed");}
	}
	private void emitGenomes(){
		final Sink original=new Sink(prefix+".genomes.fa.gz"),search=new Sink(prefix+".search.fa.gz"),canary=new Sink(prefix+".canary.search.fa.gz");
		try{stream(genomes,r->{
			final Entry e=byId.get(first(r.id));require(e!=null && !e.seen && r.length()==e.length,"Genome FASTA membership/length/duplicate mismatch: "+r.id);e.seen=true;
			if(e.selected){
				original.add(r.bases.clone(),r.id);final byte[] expanded=expand(r.bases,window,e.circular);
				search.add(expanded,r.id);searchBases+=expanded.length;
				if(e.canary){canary.add(expanded,r.id);canarySearchBases+=expanded.length;}
			}
		});}finally{closeAll(original,search,canary);}
		for(Entry e:entries){require(e.seen,"Missing source genome: "+e.id);}
		require(original.count==selectedGenomes && search.count==original.count && canary.count==canaryGenomes,"Every selected genome must reach each applicable target stream");
	}
	private void emitCandidates(){
		final Sink all=new Sink(prefix+".candidates.fa.gz"),canary=new Sink(prefix+".canary.candidates.fa.gz");
		try{stream(candidates,r->{
			final Locus l=byLocus.get(first(r.id));require(l!=null && !l.seen && r.length()==l.length,"Candidate FASTA membership/length/duplicate mismatch: "+r.id);l.seen=true;
			if(l.entry.selected){final byte[] bases=r.bases.clone();all.add(bases,r.id);candidateBases+=bases.length;
				if(l.entry.canary){canary.add(bases,r.id);canaryCandidateBases+=bases.length;}
			}
		});}finally{closeAll(all,canary);}
		for(Locus l:allLoci){require(l.seen,"Missing extracted source locus: "+l.id);}
		require(all.count==selectedLoci && canary.count==canaryLoci,"Every selected annotated locus must survive without product/score/partial filtering");
	}
	/** Counts describe the complete validated census, not only selected targets.
	 * Infernal cmsearch1.1.5 dbsize_and_seq_lengths counts both strands; -Z is Mb.
	 * Freeze these same universe denominators across canaries and future shards. */
	private void writePopulation(){
		require(populationBases>0 && populationSearchBases>=populationBases && populationCandidateBases>0,"Full-census denominators must exist before planning search size");
		final ByteStreamWriter out=table(prefix+".population.tsv");final ByteBuilder b=new ByteBuilder();
		try{
			out.println("universe\trecords\toriginal_bases\tsearch_bases\tboth_strands_Mb\tmax_original_length");
			out.print(b.append("genomic").tab().append(entries.size()).tab().append(populationBases).tab().append(populationSearchBases).tab()
				.append(2.0*populationSearchBases/1000000,6).tab().append(maxGenomeLength).nl());
			out.print(b.clear().append("candidate").tab().append(allLoci.size()).tab().append(populationCandidateBases).tab().append(populationCandidateBases).tab()
				.append(2.0*populationCandidateBases/1000000,6).tab().append(maxCandidateLength).nl());
		}finally{require(!out.poisonAndWait(),"Full-census denominator table write failed");}
	}
	public static int padding(int length,int window,boolean circular){require(length>0 && window>0,"Positive genomic length and CM W required");return circular?Math.min(window-1,length-1):0;}
	public static byte[] expand(byte[] bases,int window,boolean circular){
		require(bases!=null && bases.length>0,"Cannot expand an empty genome");final int pad=padding(bases.length,window,circular);
		require((long)bases.length+pad<Integer.MAX_VALUE,"Circular target exceeds Java array capacity");
		final byte[] result=Arrays.copyOf(bases,bases.length+pad);System.arraycopy(bases,0,result,bases.length,pad);return result;
	}
	/** Infernal seq_from/to are 1-based, strand oriented. Retain raw coordinates;
	 * project only spans of at most one genome. Hits wholly in the appended prefix
	 * are appended representations, not automatically discarded hits: a lone valid
	 * hit there must survive physical-coordinate arbitration. Multi-wrap spans fail loud. */
	public static MappedHit mapHit(int from,int to,char strand,int length,int pad){
		require(length>0 && pad>=0 && pad<length && from>0 && to>0 && Math.max(from,to)<=(long)length+pad,"Hit outside declared circular search target");
		require(strand=='+' && from<=to || strand=='-' && from>=to,"Infernal strand must agree with raw oriented coordinates, including single-base hits");
		final int lo=Math.min(from,to),hi=Math.max(from,to);require((long)hi-lo+1<=length,"Hit spans more than one physical genome; retain for explicit adjudication");
		return new MappedHit((from-1)%length+1,(to-1)%length+1,strand,lo<=length && hi>length,lo>length);
	}
	public static final class MappedHit {
		MappedHit(int f,int t,char s,boolean c,boolean d){from=f;to=t;strand=s;crossesOrigin=c;appendedDuplicate=d;}
		public final int from,to;public final char strand;public final boolean crossesOrigin,appendedDuplicate;
	}
	private static boolean bool(LineParser1 p,int term){require(p.termEquals("true",term) || p.termEquals("false",term),"Invalid census boolean");return p.termEquals("true",term);}
	private static String first(String id){require(id!=null && !id.isEmpty(),"FASTA identifier required");int n=0;while(n<id.length() && !Character.isWhitespace(id.charAt(n))){n++;}return id.substring(0,n);}
	private static void stream(String path,java.util.function.Consumer<Read> action){
		final Streamer st=StreamerFactory.makeStreamer(FileFormat.testInput(path,FileFormat.FASTA,null,true,true),null,true,-1,false,true,1);
		st.start();try{for(ListNum<Read> batch;(batch=st.nextList())!=null;){for(Read r:batch.list){action.accept(r);}st.returnList(batch);}require(!st.errorState(),"FASTA stream failed: "+path);}
		finally{require(!ReadWrite.closeStream(st),"FASTA close failed: "+path);}
	}
	private static ByteStreamWriter table(String path){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));w.start();return w;}
	private static void closeAll(Sink... sinks){
		Throwable failure=null;
		for(Sink sink:sinks){try{sink.close();}catch(Throwable t){if(failure==null){failure=t;}else{failure.addSuppressed(t);}}}
		if(failure!=null){throw new IllegalStateException("Pilot FASTA finalization failed after attempting every sink",failure);}
	}
	private static final class Sink {
		Sink(String path){writer=WriterFactory.makeWriter(FileFormat.testOutput(path,FileFormat.FASTA,null,true,false,false,true),null,1);writer.start();}
		void add(byte[] bases,String id){batch.add(new Read(bases,null,id,count++,false));if(batch.size()==64){flush();}}
		void flush(){if(!batch.isEmpty()){writer.add(batch,batchId++);batch=new ArrayList<Read>(64);}}
		void close(){flush();require(!writer.poisonAndWait() && writer.readsWritten()==count,"FASTA output count/write failure");}
		final Writer writer;ArrayList<Read> batch=new ArrayList<Read>(64);long count,batchId;
	}
	private static final class Entry {String id,organelle;int taxid,length,features,loci,quartile;boolean circular,fallback,partial,pseudo,selected,canary,seen;byte[] raw;}
	private static final class Locus {String id;Entry entry;int length;boolean seen;byte[] raw;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private String records,loci,genomes,candidates,prefix;private int expectedRecords,expectedLoci,expectedKinetoplast=-1,perGroup=8,window;
	private long maxBases=5000000,selectedBases,searchBases,candidateBases,canarySearchBases,canaryCandidateBases,populationBases,populationSearchBases,populationCandidateBases;
	private int selectedLoci,canaryLoci,canaryGenomes,selectedGenomes,maxGenomeLength,maxCandidateLength;private boolean longTail;
	private byte[] recordHeader,locusHeader;private final PrintStream log;
	private Locus longestLocus;
	private final ArrayList<Entry> entries=new ArrayList<Entry>();private final ArrayList<Locus> allLoci=new ArrayList<Locus>();
	private final ObjectMap<String,Entry> byId=new ObjectMap<String,Entry>(String.class,Entry.class);
	private final ObjectMap<String,Locus> byLocus=new ObjectMap<String,Locus>(String.class,Locus.class);
	private final int[] available=new int[4];private static final String[] STRATA={"kinetoplast","partial","pseudo","CONTIG_fallback"};
	private final int[][] quartileCounts=new int[2][4];
	private static final String[] SUFFIXES={".selection.tsv",".records.tsv",".loci.tsv.gz",".strata.tsv",".summary.tsv",".population.tsv",".genomes.fa.gz",".search.fa.gz",".candidates.fa.gz",".canary.search.fa.gz",".canary.candidates.fa.gz"};
}
