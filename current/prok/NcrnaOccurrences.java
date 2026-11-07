package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.LongIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaOccurrenceData.Event;
import prok.NcrnaOccurrenceData.Hit;
import prok.NcrnaOccurrenceData.Original;
import prok.NcrnaOccurrenceData.Overlap;
import prok.NcrnaOccurrenceData.Source;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Resolves primary-model occurrences from an accepted observation census.
 * Forked from template.A_SampleByteFile. One metadata worker applies genomic
 * AND per-event model-span containment; no sequence search, score threshold,
 * synthetic union or biological subtype classification is performed.
 * Other models are linked evidence and never bridge primary occurrences.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaOccurrences{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final NcrnaOccurrences tool=new NcrnaOccurrences(args);
		try{
			tool.process(timer);
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/** Native template parsing; only a fresh output directory is accepted. */
	public NcrnaOccurrences(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		try{
			ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
			final Parser parser=parse(args);
			in=parser.in1;
			out=parser.out1;
			if(parser.append || parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1){
				throw new IllegalArgumentException("Occurrences require complete census metadata, without append or sampling");
			}
			if(in==null || out==null || primary==null || primary.isEmpty()){
				throw new IllegalArgumentException("Required: in=census_directory primary=model_id out=fresh_directory");
			}
			if(!new File(in).isDirectory() || new File(out).exists()){
				throw new IllegalArgumentException("Input must be a census directory and output must not exist");
			}
			ReadWrite.setZipThreads(Shared.threads());
		}catch(RuntimeException | Error failure){
			Shared.closeStream(outstream);
			throw failure;
		}
	}

	private Parser parse(final String[] args){
		assert(args!=null) : "PreParser supplies the native argument array";
		final Parser parser=new Parser();
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("primary")){
				primary=b;
			}else if(parser.parse(arg, a, b)){
				// Standard BBTools options.
			}else{
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		return parser;
	}

	/** All semantic checks, including ambiguous parents, precede output creation. */
	private void process(final Timer timer){
		assert(data==null) : "An exporter instance owns one census generation and may only run once";
		data=new NcrnaOccurrenceData(in, primary);
		resolve();
		linkEvidence();
		if(!new File(out).mkdir()){throw new IllegalArgumentException("Cannot create fresh output directory: "+out);}
		writeOccurrences();
		writePhysicalLinks();
		writeEventLinks();
		writeOriginalLinks();
		writeSummary();
		timer.stop();
		outstream.println("NCRNA_OCCURRENCES_PASS primary_hits="+primaryHits+" occurrences="+roots.size()
			+" collapsed="+collapsed+" raw_events="+data.events.size()+" time="+timer);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Containment          ----------------*/
	/*--------------------------------------------------------------*/

	/** A source-local sweep keeps maximal containers, not overlap components.
	 * Parents sort before genomic subintervals. The every-event witness relation
	 * is transitive, so a collapsed parent may be replaced by its unique maximum.
	 */
	private void resolve(){
		assert(roots.isEmpty()) : "A census generation is resolved once before assigning occurrence IDs";
		final ArrayList<Hit> active=new ArrayList<Hit>();
		for(final Source source:data.sources){
			active.clear();
			char strand=0;
			for(final Hit hit:source.hits){
				if(hit.model!=data.primary){continue;}
				primaryHits++;
				if(hit.strand!=strand){
					active.clear();
					strand=hit.strand;
				}
				Hit parent=null;
				int retained=0;
				for(int i=0; i<active.size(); i++){
					final Hit candidate=active.get(i);
					if(candidate.stop<hit.start){continue;}
					active.set(retained++, candidate);
					if(!contains(candidate, hit)){continue;}
					if(parent!=null){
						throw new IllegalArgumentException("Incomparable maximal containers for physical hit "+hit.id
							+": "+parent.id+" and "+candidate.id+"; refusing to choose a parent or merge occurrences");
					}
					parent=candidate;
				}
				while(active.size()>retained){active.remove(active.size()-1);}
				if(parent==null){
					hit.representative=hit;
					roots.add(hit);
					source.roots.add(hit);
					active.add(hit);
				}else{
					hit.representative=parent;
					collapsed++;
				}
			}
		}
		roots.sort(this::compareHits);
		for(int i=0; i<roots.size(); i++){roots.get(i).occurrence=i;}
		assert(primaryHits==roots.size()+collapsed) : "Every primary physical hit is exactly one representative or collapsed observation";
	}

	/** Every child observation requires one actual containing parent observation.
	 * Several parent observations may not collectively cover one child model span.
	 */
	private boolean contains(final Hit parent, final Hit child){
		assert(parent!=child && parent.source==child.source && parent.model==child.model
			&& parent.strand==child.strand) : "Containment candidates share one source/strand/model and exclude self";
		if(parent.start>child.start || parent.stop<child.stop){return false;}
		for(int i=0; i<child.events.size; i++){
			if(witness(parent, data.events.get(child.events.get(i)))==null){return false;}
		}
		return true;
	}

	/** Lexically first real event makes witnesses stable under input reordering. */
	private Event witness(final Hit parent, final Event child){
		assert(parent.model==child.hit.model) : "Model positions are comparable only within one declared model/version";
		Event best=null;
		for(int i=0; i<parent.events.size; i++){
			final Event event=data.events.get(parent.events.get(i));
			if(event.from<=child.from && event.to>=child.to && (best==null || compareEvents(event, best)<0)){
				best=event;
			}
		}
		return best;
	}

	/** Secondary intervals link by same-strand overlap only; they never make roots. */
	private void linkEvidence(){
		assert(primaryHits==roots.size()+collapsed) : "Evidence cannot change the resolved occurrence population";
		for(final Source source:data.sources){
			source.roots.sort((a, b)->Integer.compare(a.start, b.start));
			final int[] maxEnd=new int[source.roots.size()];
			int end=0;
			for(int i=0; i<maxEnd.length; i++){
				end=Math.max(end, source.roots.get(i).stop);
				maxEnd[i]=end;
			}
			for(final Hit hit:source.hits){
				if(hit.model==data.primary){
					hit.links.add(hit.representative.occurrence);
					continue;
				}
				secondaryHits++;
				int low=0, high=maxEnd.length;
				while(low<high){
					final int mid=(low+high)>>>1;
					if(maxEnd[mid]<hit.start){low=mid+1;}else{high=mid;}
				}
				for(int i=low; i<source.roots.size(); i++){
					final Hit root=source.roots.get(i);
					if(root.start>hit.stop){break;}
					if(root.stop>=hit.start && root.strand==hit.strand){hit.links.add(root.occurrence);}
				}
				hit.links.sort();
				if(hit.links.size==0){unmappedSecondary++;}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeOccurrences(){
		assert(data!=null) : "Occurrence output follows complete census validation and parent resolution";
		try(Output output=output("occurrences.tsv")){
			output.writer.println("occurrence_id\tgenome\tassembly\ttaxid\tseqid\tstrand\tstart1\tstop1\tmodel_id\tmodel_sha80\trepresentative_physical_id\tcontig_length");
			for(final Hit hit:roots){
				row.clear();
				appendId(hit.occurrence);
				row.tab().append(data.symbols.bytes(hit.source.genome)).tab().append(data.symbols.bytes(hit.source.assembly))
					.tab().append(hit.source.taxid).tab().append(data.symbols.bytes(hit.source.seqid)).tab().append(hit.strand)
					.tab().append(hit.start).tab().append(hit.stop).tab().append(data.symbols.bytes(hit.model.id))
					.tab().append(data.symbols.bytes(hit.model.sha)).tab().append(hit.id).tab().append(hit.source.length).nl();
				output.writer.print(row);
			}
		}
	}

	private void writePhysicalLinks(){
		assert(data.hits.size()==primaryHits+secondaryHits) : "Physical links preserve every primary and secondary interval";
		try(Output output=output("physical_links.tsv")){
			output.writer.println("physical_id\toccurrence_id\trelation");
			for(final Hit hit:data.hits){
				for(int i=0; i<Math.max(1, hit.links.size); i++){
					row.clear().append(hit.id).tab();
					appendId(hit.links.size==0 ? -1 : hit.links.get(i));
					row.tab().append(relation(hit)).nl();
					output.writer.print(row);
				}
			}
		}
	}

	private void writeEventLinks(){
		assert(data.events.size()>=data.hits.size()) : "Every physical interval has at least one validated event";
		final Event[] ordered=data.events.toArray(new Event[data.events.size()]);
		Arrays.sort(ordered, this::compareEvents);
		try(Output output=output("event_links.tsv")){
			output.writer.println("table_id\ttable_line\tphysical_id\toccurrence_id\trelation\twitness_table\twitness_line\tcontext_id\torigin_genome\torigin_id");
			for(final Event event:ordered){
				final Hit hit=event.hit;
				final Event parent=hit.model==data.primary ? witness(hit.representative, event) : null;
				assert(hit.model!=data.primary || parent!=null) : "Every collapsed event must retain its actual maximal-parent witness";
				for(int i=0; i<Math.max(1, hit.links.size); i++){
					row.clear().append(data.symbols.bytes(event.search.id)).tab().append(event.line).tab().append(hit.id).tab();
					appendId(hit.links.size==0 ? -1 : hit.links.get(i));
					row.tab().append(relation(hit)).tab();
					if(parent==null){row.append(".\t.");}
					else{row.append(data.symbols.bytes(parent.search.id)).tab().append(parent.line);}
					row.tab().append(data.symbols.bytes(event.context)).tab().append(data.symbols.bytes(event.origin.source.genome))
						.tab().append(data.symbols.bytes(event.origin.id)).nl();
					output.writer.print(row);
					eventLinks++;
				}
			}
		}
	}

	/** Secondary annotation overlap is retained but never transferred across its broad span. */
	private void writeOriginalLinks(){
		final Overlap[] ordered=data.overlaps.toArray(new Overlap[data.overlaps.size()]);
		Arrays.sort(ordered, (a, b)->{
			final int c=compareOriginals(a.original, b.original);
			return c==0 ? Integer.compare(a.hit.id, b.hit.id) : c;
		});
		final LongIntMap seen=new LongIntMap();
		try(Output output=output("original_links.tsv")){
			output.writer.println("genome\toriginal_id\tphysical_id\toccurrence_id\toverlap_bases\tsame_strand\trelation");
			for(final Overlap link:ordered){
				assert(link.original.source==link.hit.source) : "Annotation transfer requires the full source identity validated by loadOverlaps";
				final int occurrence=link.hit.model==data.primary ? link.hit.representative.occurrence : -1;
				row.clear().append(data.symbols.bytes(link.original.source.genome)).tab().append(data.symbols.bytes(link.original.id))
					.tab().append(link.hit.id).tab();
				appendId(occurrence);
				row.tab().append(link.bases).tab().append(link.same ? "true" : "false").tab()
					.append(occurrence<0 ? "SECONDARY_ANNOTATION_EVIDENCE" : "PRIMARY_ANNOTATION_OVERLAP").nl();
				output.writer.print(row);
				if(occurrence>=0){
					final long key=NcrnaOccurrenceData.pair(link.original.index, occurrence);
					if(seen.get(key)<0){
						seen.put(key, 1);
						link.original.occurrences++;
					}
				}
			}
		}
		final Original[] originals=data.originals.toArray(new Original[data.originals.size()]);
		Arrays.sort(originals, this::compareOriginals);
		try(Output output=output("original_status.tsv")){
			output.writer.println("genome\toriginal_id\tobserved_physical_overlaps\toccurrences_with_primary_overlap");
			for(final Original original:originals){
				row.clear().append(data.symbols.bytes(original.source.genome)).tab().append(data.symbols.bytes(original.id))
					.tab().append(original.physicalHits).tab().append(original.occurrences).nl();
				output.writer.print(row);
			}
		}
	}

	private void writeSummary(){
		assert(primaryHits==collapsed+roots.size() && eventLinks>=data.events.size())
			: "Occurrence and event-link accounting must conserve every input observation";
		try(Output output=output("summary.tsv")){
			output.writer.println("metric\tcount");
			metric(output.writer, "original_annotations", data.originals.size());
			metric(output.writer, "searches", data.searches.size());
			metric(output.writer, "physical_hits", data.hits.size());
			metric(output.writer, "raw_events", data.events.size());
			metric(output.writer, "annotation_links", data.overlaps.size());
			metric(output.writer, "primary_physical_hits", primaryHits);
			metric(output.writer, "collapsed_primary_hits", collapsed);
			metric(output.writer, "occurrences", roots.size());
			metric(output.writer, "secondary_physical_hits", secondaryHits);
			metric(output.writer, "unmapped_secondary_hits", unmappedSecondary);
			metric(output.writer, "event_link_rows", eventLinks);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Formatting/order       ----------------*/
	/*--------------------------------------------------------------*/

	private Output output(final String name){return new Output(NcrnaOccurrenceData.path(out, name));}

	private void appendId(final int occurrence){
		assert(occurrence>=-1 && occurrence<roots.size()) : "Only resolved occurrence IDs or explicit unassigned evidence may be emitted";
		if(occurrence<0){row.append('.');}
		else{row.append(data.symbols.bytes(data.primary.id)).append("_occurrence_").append(occurrence+1);}
	}

	private static String relation(final Hit hit){
		if(hit.representative==hit){return "REPRESENTATIVE";}
		if(hit.representative!=null){return "CONTAINED";}
		return hit.links.size==0 ? "UNMAPPED_SECONDARY" : "SECONDARY_OVERLAP";
	}

	private void metric(final ByteStreamWriter writer, final String name, final long value){
		assert(value>=0) : "Summary counts represent observed rows, not missing-value sentinels: "+name+"="+value;
		row.clear().append(name).tab().append(value).nl();
		writer.print(row);
	}

	private int compareHits(final Hit a, final Hit b){
		int c=compareSymbols(a.source.genome, b.source.genome);
		if(c==0){c=compareSymbols(a.source.assembly, b.source.assembly);}
		if(c==0){c=compareSymbols(a.source.seqid, b.source.seqid);}
		if(c==0){c=Character.compare(a.strand, b.strand);}
		if(c==0){c=Integer.compare(a.start, b.start);}
		if(c==0){c=Integer.compare(a.stop, b.stop);}
		return c==0 ? Integer.compare(a.id, b.id) : c;
	}

	private int compareEvents(final Event a, final Event b){
		final int c=compareSymbols(a.search.id, b.search.id);
		return c==0 ? Integer.compare(a.line, b.line) : c;
	}

	private int compareOriginals(final Original a, final Original b){
		final int c=compareSymbols(a.source.genome, b.source.genome);
		return c==0 ? compareSymbols(a.id, b.id) : c;
	}

	/** Native byte order makes output independent of input symbol numbering. */
	private int compareSymbols(final int a, final int b){
		return Tools.compare(data.symbols.bytes(a), data.symbols.bytes(b));
	}

	private String in, out, primary;
	private PrintStream outstream;
	private NcrnaOccurrenceData data;
	private final ArrayList<Hit> roots=new ArrayList<Hit>();
	private final ByteBuilder row=new ByteBuilder();
	private int primaryHits, secondaryHits, collapsed, unmappedSecondary;
	private long eventLinks;
}
