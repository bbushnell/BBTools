package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Comparator;

import map.IntHashMap;
import map.LongIntMap;
import parse.LineParser1;
import prok.NcrnaObservationData.Symbols;
import prok.NcrnaObservationData.Table;
import structures.ByteBuilder;
import structures.IntList;

/** The metadata needed to resolve occurrences from an accepted observation census.
 * Retains original evidence and source-local identities without rereading genomes,
 * query FASTAs or covariance-model search output. One metadata thread owns it.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaOccurrenceData{

	/** Loads complete tables and validates their identity/count joins before output. */
	NcrnaOccurrenceData(final String directory, final String primaryName){
		originalTable=new Table(path(directory, "originals.tsv"), NcrnaCmLocusAudit.HEADER);
		searchTable=new Table(path(directory, "searches.tsv"), NcrnaCmHitProjector.TABLE_HEADER);
		hitTable=new Table(path(directory, "hits.tsv"), NcrnaObservationData.HIT_HEADER);
		eventTable=new Table(path(directory, "events.tsv"), NcrnaObservationData.EVENT_HEADER);
		overlapTable=new Table(path(directory, "overlaps.tsv"), NcrnaObservationData.OVERLAP_HEADER);
		loadOriginals();
		loadSearches(primaryName);
		loadHits();
		loadEvents();
		loadOverlaps();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Metadata joins        ----------------*/
	/*--------------------------------------------------------------*/

	/** Full source keys distinguish equally named contigs in different assemblies. */
	private void loadOriginals(){
		assert(originals.isEmpty()) : "Original identities are established once before event joins";
		for(final byte[] row:originalTable.rows){
			parse(row, 21);
			final int id=symbol(0), genome=symbol(1), assembly=symbol(2), seqid=symbol(6);
			final int key=sourceKey(genome, assembly, seqid);
			final int length=parser.parseInt(15), taxid=parser.parseInt(3);
			require(length>0 && taxid>=0, "Invalid original contig length or taxid");
			int index=bySource.get(key);
			if(index<0){
				index=sources.size();
				bySource.put(key, index);
				sources.add(new Source(genome, assembly, seqid, length, taxid));
			}
			final Source source=sources.get(index);
			require(source.length==length && source.taxid==taxid, "Conflicting metadata for one source contig");
			final long originalKey=pair(genome, id);
			require(byOriginal.get(originalKey)<0, "Repeated original genome/locus ID");
			byOriginal.put(originalKey, originals.size());
			originals.add(new Original(originals.size(), id, source, strand(7)));
		}
	}

	/** A selected model has one declared version; other models remain evidence only. */
	private void loadSearches(final String name){
		assert(searches.isEmpty()) : "Search identities precede event interpretation";
		final byte[] requested=name.getBytes(StandardCharsets.UTF_8);
		final int requestedId=symbols.intern(requested, 0, requested.length);
		for(final byte[] row:searchTable.rows){
			parse(row, 7);
			final int id=symbol(0), modelId=symbol(1), sha=symbol(2), clen=parser.parseInt(3);
			final double ga=parser.parseDouble(4);
			require(clen>0 && Double.isFinite(ga), "Invalid declared model length or gathering threshold");
			final byte[] hash=symbols.bytes(sha);
			require(hash.length==20, "Model versions require the census sha80 identity");
			for(final byte b:hash){
				require(b>='0' && b<='9' || b>='a' && b<='f', "Invalid model sha80 text");
			}
			final long key=pair(modelId, sha);
			int index=byModel.get(key);
			if(index<0){
				index=models.size();
				byModel.put(key, index);
				final Model model=new Model(index, modelId, sha, clen, ga);
				models.add(model);
				if(modelId==requestedId){
					require(primary==null, "Select a census with one version of the primary model: "+name);
					primary=model;
				}
			}
			final Model model=models.get(index);
			require(model.clen==clen && model.ga==ga, "Conflicting metadata for a model/version");
			require(bySearch.get(id)<0, "Repeated search table ID");
			bySearch.put(id, searches.size());
			searches.add(new Search(searches.size(), id, model));
		}
		require(primary!=null, "Primary model is absent from searches.tsv: "+name);
	}

	/** Dense physical IDs are preserved, independent of the input row order. */
	private void loadHits(){
		assert(hits.isEmpty() && primary!=null) : "Physical hits require declared model identities";
		for(final byte[] row:hitTable.rows){
			parse(row, 11);
			final int id=parser.parseInt(0);
			final int source=bySource.get(sourceKey(symbol(1), symbol(2), symbol(3)));
			final int model=byModel.get(pair(symbol(5), symbol(6)));
			require(id>0 && source>=0 && model>=0, "Physical hit has an unknown source/model or invalid ID");
			final int start=parser.parseInt(7), stop=parser.parseInt(8);
			final int count=parser.parseInt(9), overlaps=parser.parseInt(10);
			final Source owner=sources.get(source);
			if(start<1 || stop<start || stop>owner.length || count<1 || overlaps<0){
				throw new IllegalArgumentException("Invalid physical coordinates or evidence counts: "+id);
			}
			final Hit hit=new Hit(id, owner, models.get(model), strand(4), start, stop, count, overlaps);
			hits.add(hit);
			owner.hits.add(hit);
		}
		hits.sort((a, b)->Integer.compare(a.id, b.id));
		for(int i=0; i<hits.size(); i++){
			require(hits.get(i).id==i+1, "Physical IDs must be unique and dense from one");
		}
		for(final Source source:sources){
			source.hits.sort(HIT_ORDER);
			for(int i=1; i<source.hits.size(); i++){
				final Hit a=source.hits.get(i-1), b=source.hits.get(i);
				require(a.model!=b.model || a.strand!=b.strand || a.start!=b.start || a.stop!=b.stop,
					"Repeated physical source/model/strand/coordinate tuple");
			}
		}
	}

	/** Model spans remain individual observations; no synthetic span is constructed. */
	private void loadEvents(){
		assert(events.isEmpty()) : "Raw event keys are registered once";
		final LongIntMap seen=new LongIntMap();
		for(final byte[] row:eventTable.rows){
			parse(row, 15);
			final Hit hit=hit(parser.parseInt(0));
			final int searchIndex=bySearch.get(symbol(1)), line=parser.parseInt(2);
			require(searchIndex>=0 && line>0, "Unknown event search or invalid original table line");
			final Search search=searches.get(searchIndex);
			final long key=pair(searchIndex, line);
			require(seen.get(key)<0, "Repeated original search-table/line event");
			seen.put(key, events.size());
			final Original origin=original(symbol(4), symbol(5));
			require(origin.source==hit.source && search.model==hit.model,
				"Event origin/model differs from the physical hit's source/version");
			final int from=parser.parseInt(6), to=parser.parseInt(7);
			require(from>=1 && to>=from && to<=hit.model.clen, "Observed model span exceeds CLEN");
			final Event event=new Event(events.size(), hit, search, line, symbol(3), origin, from, to);
			hit.events.add(event.index);
			events.add(event);
		}
		for(final Hit hit:hits){
			if(hit.events.size!=hit.expectedEvents){
				throw new IllegalArgumentException("Physical hit event count mismatch: "+hit.id);
			}
		}
	}

	/** Preserve annotation relationships separately from event search origins. */
	private void loadOverlaps(){
		assert(overlaps.isEmpty()) : "Original annotation links load once after events";
		final LongIntMap seen=new LongIntMap();
		for(final byte[] row:overlapTable.rows){
			parse(row, 6);
			final Hit hit=hit(parser.parseInt(0));
			final Original original=original(symbol(1), symbol(2));
			final int bases=parser.parseInt(3);
			final boolean same=bool(4);
			bool(5);
			require(original.source==hit.source && bases>0 && same==(original.strand==hit.strand),
				"Annotation link has inconsistent source, overlap or strand");
			final long key=pair(hit.id, original.index);
			require(seen.get(key)<0, "Repeated physical-hit/original link");
			seen.put(key, overlaps.size());
			overlaps.add(new Overlap(hit, original, bases, same));
			hit.observedOverlaps++;
			original.physicalHits++;
		}
		for(final Hit hit:hits){
			if(hit.observedOverlaps!=hit.expectedOverlaps){
				throw new IllegalArgumentException("Physical hit annotation count mismatch: "+hit.id);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Byte parsing         ----------------*/
	/*--------------------------------------------------------------*/

	private void parse(final byte[] row, final int fields){
		assert(row!=null && fields>0) : "Table supplies a retained row and the reader supplies its schema width";
		parser.set(row);
		current=row;
		if(parser.terms()!=fields){
			throw new IllegalArgumentException("Unexpected census metadata field count; expected "+fields);
		}
	}

	private int symbol(final int column){
		parser.setBounds(column);
		if(parser.a()>=parser.b()){throw new IllegalArgumentException("Empty census identity field: "+column);}
		return symbols.intern(current, parser.a(), parser.b());
	}

	private int sourceKey(final int genome, final int assembly, final int seqid){
		assert(genome>=0 && assembly>=0 && seqid>=0) : "Full source keys require three interned identities, never absent-map sentinels";
		key.clear().append(symbols.bytes(genome)).tab().append(symbols.bytes(assembly))
			.tab().append(symbols.bytes(seqid));
		return symbols.intern(key.array, 0, key.length());
	}

	private char strand(final int column){
		if(parser.termEquals('+', column)){return '+';}
		if(parser.termEquals('-', column)){return '-';}
		throw new IllegalArgumentException("Expected explicit genomic strand");
	}

	private boolean bool(final int column){
		if(parser.termEquals("true", column)){return true;}
		if(parser.termEquals("false", column)){return false;}
		throw new IllegalArgumentException("Expected true/false census flag");
	}

	private Original original(final int genome, final int id){
		final int index=byOriginal.get(pair(genome, id));
		require(index>=0, "Unknown original genome/locus identity");
		return originals.get(index);
	}

	Hit hit(final int id){
		if(id<1 || id>hits.size()){throw new IllegalArgumentException("Physical ID outside census: "+id);}
		return hits.get(id-1);
	}

	static long pair(final int a, final int b){return ((long)a<<32)|(b&0xffffffffL);}
	static String path(final String directory, final String name){return new File(directory, name).getPath();}
	static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Retained data        ----------------*/
	/*--------------------------------------------------------------*/

	static final class Source{
		Source(int g, int a, int s, int length_, int taxid_){
			genome=g;
			assembly=a;
			seqid=s;
			length=length_;
			taxid=taxid_;
		}
		final int genome, assembly, seqid, length, taxid;
		final ArrayList<Hit> hits=new ArrayList<Hit>();
		final ArrayList<Hit> roots=new ArrayList<Hit>();
	}

	static final class Model{
		Model(int index_, int id_, int sha_, int clen_, double ga_){
			index=index_;
			id=id_;
			sha=sha_;
			clen=clen_;
			ga=ga_;
		}
		final int index, id, sha, clen;
		final double ga;
	}

	static final class Search{
		Search(int index_, int id_, Model model_){
			index=index_;
			id=id_;
			model=model_;
		}
		final int index, id;
		final Model model;
	}

	static final class Original{
		Original(int index_, int id_, Source source_, char strand_){
			index=index_;
			id=id_;
			source=source_;
			strand=strand_;
		}
		final int index, id;
		final Source source;
		final char strand;
		int physicalHits, occurrences;
	}

	static final class Hit{
		Hit(int id_, Source source_, Model model_, char strand_, int start_, int stop_, int events_, int overlaps_){
			id=id_;
			source=source_;
			model=model_;
			strand=strand_;
			start=start_;
			stop=stop_;
			expectedEvents=events_;
			expectedOverlaps=overlaps_;
		}
		final int id, start, stop, expectedEvents, expectedOverlaps;
		final Source source;
		final Model model;
		final char strand;
		final IntList events=new IntList(1);
		final IntList links=new IntList(1);
		Hit representative;
		int occurrence=-1, observedOverlaps;
	}

	static final class Event{
		Event(int index_, Hit hit_, Search search_, int line_, int context_, Original origin_, int from_, int to_){
			index=index_;
			hit=hit_;
			search=search_;
			line=line_;
			context=context_;
			origin=origin_;
			from=from_;
			to=to_;
		}
		final int index, line, context, from, to;
		final Hit hit;
		final Search search;
		final Original origin;
	}

	static final class Overlap{
		Overlap(Hit hit_, Original original_, int bases_, boolean same_){
			hit=hit_;
			original=original_;
			bases=bases_;
			same=same_;
		}
		final Hit hit;
		final Original original;
		final int bases;
		final boolean same;
	}

	/** Parents precede genomic subintervals; equal physical tuples are rejected. */
	static final Comparator<Hit> HIT_ORDER=(a, b)->{
		int c=Integer.compare(a.model.index, b.model.index);
		if(c==0){c=Character.compare(a.strand, b.strand);}
		if(c==0){c=Integer.compare(a.start, b.start);}
		if(c==0){c=Integer.compare(b.stop, a.stop);}
		return c==0 ? Integer.compare(a.id, b.id) : c;
	};

	final Symbols symbols=new Symbols();
	final ArrayList<Source> sources=new ArrayList<Source>();
	final ArrayList<Model> models=new ArrayList<Model>();
	final ArrayList<Search> searches=new ArrayList<Search>();
	final ArrayList<Original> originals=new ArrayList<Original>();
	final ArrayList<Hit> hits=new ArrayList<Hit>();
	final ArrayList<Event> events=new ArrayList<Event>();
	final ArrayList<Overlap> overlaps=new ArrayList<Overlap>();
	Model primary;
	private final Table originalTable, searchTable, hitTable, eventTable, overlapTable;
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder key=new ByteBuilder();
	private final IntHashMap bySource=new IntHashMap(), bySearch=new IntHashMap();
	private final LongIntMap byOriginal=new LongIntMap(), byModel=new LongIntMap();
	private byte[] current;
}
