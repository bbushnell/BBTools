package prok;

import java.io.File;
import java.util.ArrayList;
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

/** Base-balanced whole-accession CM inputs from the complete retained census.
 * Genomic and candidate partitions are independent; all loci of an accession
 * stay together within each universe. Circular padding uses the reviewed pilot rule.
 * @author Ganyu, Raiden */
public final class NcrnaMitoPartitions {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoPartitions.class,false);final Parser parser=new Parser();
		final NcrnaMitoPartitions x=new NcrnaMitoPartitions();
		try{
			for(String arg:pp.args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
				if(k.equals("records")){x.records=v;}else if(k.equals("loci")){x.loci=v;}else if(k.equals("genomes")){x.genomes=v;}
				else if(k.equals("candidates")){x.candidates=v;}else if(k.equals("outprefix")){x.prefix=v;}
				else if(k.equals("partitions")){x.parts=Parse.parseIntKMG(v);}else if(k.equals("window")){x.window=Parse.parseIntKMG(v);}
				else if(k.equals("genomicmaxusable")){x.maxUsable=Parse.parseIntKMG(v);}
				else if(k.equals("genomiconly")){x.genomicOnly=Parse.parseBoolean(v);}else if(k.equals("metadataonly")){x.metadataOnly=Parse.parseBoolean(v);}
				else if(k.equals("expectedrecords")){x.expectedRecords=Parse.parseIntKMG(v);}else if(k.equals("expectedloci")){x.expectedLoci=Parse.parseIntKMG(v);}
				else if(k.equals("expectedsearchbases")){x.expectedSearch=Parse.parseKMG(v);}else if(k.equals("expectedcandidatebases")){x.expectedCandidate=Parse.parseKMG(v);}
				else if(!parser.parse(arg,k,v)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			require(x.records!=null && x.loci!=null && x.prefix!=null,"Require records/loci/outprefix for the complete census");
			require(x.maxUsable>=-1 && (x.maxUsable>=0)==x.genomicOnly && (!x.metadataOnly || x.genomicOnly),"Restricted genomic output needs explicit genomicmaxusable and genomiconly; metadataonly is only for that census");
			require(x.metadataOnly || (x.genomes!=null && (x.genomicOnly || x.candidates!=null)),"Sequence emission requires genomes and, unless genomiconly, candidates");
			require(x.parts>0 && x.parts<=128 && x.window>0 && x.expectedRecords>0 && x.expectedLoci>0 && x.expectedSearch>0 && x.expectedCandidate>0,"Explicit complete population, positive window and 1..128 requested partitions required");
			require(!parser.append && parser.maxReads<0,"Partitioning cannot append or truncate the census");
			Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;ReadWrite.setZipThreads(Math.min(4,Shared.threads()));
			x.process();
			if(x.genomicOnly){pp.outstream.println("MITO_GENOMIC_SUBSET_PASS census_genomes="+x.entries.size()+" census_loci="+x.annotations.size()+" selected_genomes="+x.selectedRecords+" selected_no_rrna="+x.selectedNoRrna+" original_bases="+x.selectedOriginal+" search_bases="+x.selectedSearch+" partitions="+x.parts+" max_usable="+x.maxUsable+" metadata_only="+x.metadataOnly);}
			else{pp.outstream.println("MITO_PARTITIONS_PASS genomes="+x.entries.size()+" loci="+x.annotations.size()+" partitions_per_universe="+x.parts+" search_bases="+x.searchBases+" candidate_bases="+x.candidateBases);}
		}finally{Shared.closeStream(pp.outstream);}
	}
	void process(){
		RrnaResourceIO.fresh(prefix,".selection.tsv",".membership.tsv",".partitions.tsv",".population.tsv");
		if(genomicOnly){RrnaResourceIO.fresh(prefix,".eligibility.tsv");}
		load();
		if(genomicOnly){parts=Math.min(parts,selectedRecords);require(parts>0,"No genomic records satisfy the declared usable-annotation rule");}
		for(String universe:genomicOnly?new String[]{"genomic"}:new String[]{"genomic","candidate"}){for(int i=0;i<parts;i++){require(!new File(path(universe,i)).exists(),"Refusing an existing partition");}}
		assign(false);if(!genomicOnly){assign(true);}tables();
		if(metadataOnly){return;}
		emit(false);if(!genomicOnly){emit(true);}partitionTable();
	}
	void load(){
		assert(entries.isEmpty() && annotations.isEmpty()):"One frozen census defines both complete target universes";
		read(records,14,RECORD_HEADER,p->{
			final Entry e=new Entry();e.id=p.parseString(0);e.taxid=p.parseInt(1);e.organelle=p.parseString(2);e.length=p.parseInt(5);e.features=p.parseInt(7);
			require(e.taxid>0 && e.length>0 && e.features>=0 && e.features==p.parseInt(8) && p.parseInt(9)==0,"Only fully extracted census records may be partitioned");
			require(p.termEquals("circular",6) || p.termEquals("linear_or_unspecified",6),"Unknown topology");
			require(p.termEquals("GBFF_ORIGIN",11) || p.termEquals("FNA_CONTIG",11),"Unknown sequence provenance");
			e.circular=p.termEquals("circular",6);e.fallback=p.termEquals("FNA_CONTIG",11);e.pad=NcrnaMitoPilot.padding(e.length,window,e.circular);
			require(byId.put(e.id,e)==null,"Duplicate genome accession: "+e.id);entries.add(e);originalBases+=e.length;searchBases+=(long)e.length+e.pad;
		});
		read(loci,18,LOCUS_HEADER,p->{
			final Annotation a=new Annotation();a.id=p.parseString(0);a.entry=byId.get(p.parseString(1));a.length=p.parseInt(16);
			require(a.entry!=null && p.termEquals("EXTRACTED",7) && a.length>0 && p.parseInt(2)==a.entry.taxid && p.termEquals(a.entry.organelle,3),"Annotation identity/extraction mismatch");
			final boolean partial=bit(p,8),pseudo=bit(p,9);
			a.entry.partial|=partial;a.entry.pseudo|=pseudo;if(!partial && !pseudo){a.entry.usable++;}
			a.entry.count++;a.entry.candidateBases+=a.length;candidateBases+=a.length;
			require(byAnnotation.put(a.id,a)==null,"Duplicate annotation ID: "+a.id);annotations.add(a);
		});
		require(entries.size()==expectedRecords && annotations.size()==expectedLoci && searchBases==expectedSearch && candidateBases==expectedCandidate,"Complete census counts/base totals must equal the reviewed population");
		for(Entry e:entries){require(e.count==e.features,"Per-accession annotation conservation failed: "+e.id);
			e.eligible=maxUsable<0 || e.usable<=maxUsable;
			if(e.eligible){selectedRecords++;selectedOriginal+=e.length;selectedSearch+=(long)e.length+e.pad;if(e.count==0){selectedNoRrna++;}}
		}
		entries.sort(Comparator.comparing(e->e.id));annotations.sort(Comparator.comparing(a->a.id));
		for(boolean annotated:new boolean[]{false,true}){
			final ArrayList<Entry> group=new ArrayList<Entry>();for(Entry e:entries){if((e.features>0)==annotated){group.add(e);}}
			group.sort(Comparator.comparingInt((Entry e)->e.length).thenComparing(e->e.id));
			for(int i=0;i<group.size();i++){group.get(i).quartile=(int)(4L*i/group.size());}
		}
	}
	/** Longest-first assignment to the lightest bin; ties are accession then bin index.
	 * Metadata order cannot affect membership, and an accession is never split. */
	void assign(boolean candidate){
		assert(!entries.isEmpty()):"Whole-accession lengths are known before balancing";
		final ArrayList<Entry> sorted=new ArrayList<Entry>();for(Entry e:entries){if(candidate?e.count>0:e.eligible){sorted.add(e);}}
		require(sorted.size()>=parts,"Need at least one accession per partition in each universe");
		sorted.sort(Comparator.comparingLong((Entry e)->weight(e,candidate)).reversed().thenComparing(e->e.id));
		final long[] loads=new long[parts];
		for(Entry e:sorted){int best=0;for(int i=1;i<parts;i++){if(loads[i]<loads[best]){best=i;}}
			loads[best]+=weight(e,candidate);if(candidate){e.candidatePart=best;}else{e.genomicPart=best;}
		}
		for(long n:loads){require(n>0,"An empty partition cannot become a search job");}
	}
	static long weight(Entry e,boolean candidate){return candidate?e.candidateBases:(long)e.length+e.pad;}
	void tables(){
		assert(entries.size()==expectedRecords):"Metadata output retains every genome, including unannotated ones";
		final ByteStreamWriter selection=RrnaResourceIO.writer(prefix+".selection.tsv"),membership=RrnaResourceIO.writer(prefix+".membership.tsv"),population=RrnaResourceIO.writer(prefix+".population.tsv");final ByteBuilder b=new ByteBuilder();
		try{
			selection.println("accession\ttaxid\torganelle_class\tannotation_group\tlength_quartile0\toriginal_length\tprefix_added\tsearch_length\tpartial\tpseudo\tcontig_fallback\tcanary");
			membership.println("universe\tpartition\taccession\ttarget\toriginal_length\tsearch_length");
			for(Entry e:entries){if(!e.eligible){continue;}selection.print(b.clear().append(e.id).tab().append(e.taxid).tab().append(e.organelle).tab().append(e.count>0?"annotated":"no_rrna").tab().append(e.quartile).tab().append(e.length).tab().append(e.pad).tab().append((long)e.length+e.pad).tab().append(e.partial).tab().append(e.pseudo).tab().append(e.fallback).append("\tfalse\n"));
				membership.print(b.clear().append("genomic").tab().append(e.genomicPart).tab().append(e.id).tab().append(e.id).tab().append(e.length).tab().append((long)e.length+e.pad).nl());}
			if(!genomicOnly){for(Annotation a:annotations){membership.print(b.clear().append("candidate").tab().append(a.entry.candidatePart).tab().append(a.entry.id).tab().append(a.id).tab().append(a.length).tab().append(a.length).nl());}}
			population.println("universe\trecords\toriginal_bases\tsearch_bases\tboth_strands_Mb");
			population.print(b.clear().append("genomic").tab().append(selectedRecords).tab().append(selectedOriginal).tab().append(selectedSearch).tab().append(2.0*selectedSearch/1000000,6).nl());
			if(!genomicOnly){population.print(b.clear().append("candidate").tab().append(annotations.size()).tab().append(candidateBases).tab().append(candidateBases).tab().append(2.0*candidateBases/1000000,6).nl());}
		}finally{RrnaResourceIO.finish(selection,membership,population);}
		if(genomicOnly){eligibility();}
	}
	/** Reports every included and excluded accession from source flags, never product names. */
	void eligibility(){
		assert(genomicOnly && maxUsable>=0):"An explicit usable-count rule is required before reporting genomic exclusions";
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".eligibility.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("accession\trrna_features\tusable_rrna\tgenomic_selected\toriginal_length\tprefix_added\tsearch_length");
			for(Entry e:entries){out.print(b.clear().append(e.id).tab().append(e.features).tab().append(e.usable).tab().append(e.eligible).tab().append(e.length).tab().append(e.pad).tab().append((long)e.length+e.pad).nl());}
		}finally{RrnaResourceIO.finish(out);}
	}
	void emit(boolean candidate){
		final Sink[] sinks=new Sink[parts];
		final String input=candidate?candidates:genomes;
		Streamer stream=null;
		try{
			for(int i=0;i<parts;i++){sinks[i]=new Sink(path(candidate?"candidate":"genomic",i));}
			stream=StreamerFactory.makeStreamer(FileFormat.testInput(input,FileFormat.FASTA,null,true,true),null,true,-1,false,true,1);stream.start();
			for(ListNum<Read> batch;(batch=stream.nextList())!=null;){for(Read r:batch.list){
				final String id=RrnaResourceIO.first(r.id);
				if(candidate){final Annotation a=byAnnotation.get(id);require(a!=null && !a.seen && a.length==r.length(),"Candidate FASTA membership/length/duplicate mismatch: "+id);a.seen=true;sinks[a.entry.candidatePart].add(r.bases.clone(),r.id);}
				else{final Entry e=byId.get(id);require(e!=null && !e.seen && e.length==r.length(),"Genome FASTA membership/length/duplicate mismatch: "+id);e.seen=true;if(e.eligible){sinks[e.genomicPart].add(NcrnaMitoPilot.expand(r.bases,window,e.circular),r.id);}}
			}stream.returnList(batch);}require(!stream.errorState(),"Census FASTA stream failed");
		}finally{
			final boolean error=stream!=null && ReadWrite.closeStream(stream);Throwable failure=null;
			for(Sink sink:sinks){if(sink==null){continue;}try{sink.close();}catch(Throwable t){if(failure==null){failure=t;}else{failure.addSuppressed(t);}}}
			if(failure!=null){throw new IllegalStateException("Partition finalization failed",failure);}require(!error,"Census FASTA close failed");
		}
		long count=0,bases=0;for(Sink s:sinks){count+=s.count;bases+=s.bases;require(s.count>0,"No empty search partition");}
		if(candidate){for(Annotation a:annotations){require(a.seen,"Missing candidate FASTA: "+a.id);}require(count==expectedLoci && bases==candidateBases,"Candidate output must conserve every complete sequence");candidateSinks=sinks;}
		else{for(Entry e:entries){require(e.seen,"Missing genomic FASTA: "+e.id);}require(count==selectedRecords && bases==selectedSearch,"Genome output must conserve every eligible complete padded sequence");genomicSinks=sinks;}
	}
	void partitionTable(){
		assert(genomicSinks!=null && (genomicOnly || candidateSinks!=null)):"Publish partitions only after all requested universes finish and conserve counts";
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".partitions.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("universe\tpartition\trecords\tbases\tfasta");for(int u=0;u<(genomicOnly?1:2);u++){final String universe=u==0?"genomic":"candidate";final Sink[] sinks=u==0?genomicSinks:candidateSinks;
			for(int i=0;i<parts;i++){out.print(b.clear().append(universe).tab().append(i).tab().append(sinks[i].count).tab().append(sinks[i].bases).tab().append(path(universe,i)).nl());}
		}}finally{RrnaResourceIO.finish(out);}
	}
	String path(String universe,int part){return prefix+'.'+universe+'.'+part+".fa";}
	static boolean bit(LineParser1 p,int column){require(p.termEquals("true",column) || p.termEquals("false",column),"Invalid census boolean");return p.termEquals("true",column);}
	static void read(String path,int columns,String header,java.util.function.Consumer<LineParser1> row){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,header,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==columns,"Census row width changed");row.accept(p);}}
		finally{require(!in.close(),"Census table read failed");}
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Sink{
		Sink(String path){writer=WriterFactory.makeWriter(FileFormat.testOutput(path,FileFormat.FASTA,null,true,false,false,true),null,1);writer.start();}
		void add(byte[] sequence,String id){batch.add(new Read(sequence,null,id,count++,false));bases+=sequence.length;buffered+=sequence.length;if(batch.size()>=64 || buffered>=1048576){flush();}}
		void flush(){if(!batch.isEmpty()){writer.add(batch,batchId++);batch=new ArrayList<Read>(64);buffered=0;}}
		void close(){flush();require(!writer.poisonAndWait() && writer.readsWritten()==count,"Partition FASTA write/count failure");}
		final Writer writer;ArrayList<Read> batch=new ArrayList<Read>(64);long count,bases,batchId,buffered;
	}
	static final class Entry{String id,organelle;int taxid,length,pad,features,count,usable,quartile,genomicPart=-1,candidatePart=-1;long candidateBases;boolean circular,fallback,partial,pseudo,seen,eligible;}
	static final class Annotation{String id;Entry entry;int length;boolean seen;}
	String records,loci,genomes,candidates,prefix;int parts=32,window,expectedRecords,expectedLoci,maxUsable=-1,selectedRecords,selectedNoRrna;
	long expectedSearch,expectedCandidate,originalBases,searchBases,candidateBases,selectedOriginal,selectedSearch;boolean genomicOnly,metadataOnly;
	final ArrayList<Entry> entries=new ArrayList<Entry>();final ArrayList<Annotation> annotations=new ArrayList<Annotation>();
	final ObjectMap<String,Entry> byId=new ObjectMap<String,Entry>(String.class,Entry.class);final ObjectMap<String,Annotation> byAnnotation=new ObjectMap<String,Annotation>(String.class,Annotation.class);
	Sink[] genomicSinks,candidateSinks;
	static final String RECORD_HEADER="accession\ttaxid\torganelle_class\torganism\tstatus\tbases\ttopology\trrna_features\textracted\tunsupported\tgbff\tsequence_source\tsequence_path\tcontig_location";
	static final String LOCUS_HEADER="locus_id\taccession\ttaxid\torganelle_class\trrna_index\traw_location\tproduct\tstatus\tpartial\tpseudo\tpseudogene\tstrand\tmin_coordinate\tmax_coordinate\tsegments_in_transcript_order\tcrosses_origin\textracted_length\tgbff";
}
