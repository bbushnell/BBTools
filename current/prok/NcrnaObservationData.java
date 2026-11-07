package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import map.IntHashMap;
import map.LongIntMap;
import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.SequenceInput;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/** Retained observations and source-local identities for the ncRNA census.
 * Loads accepted projector tables without changing their rows. Byte symbols are
 * interned only once; subsequent field lookups borrow a reusable slice key.
 * One metadata thread owns this object and every mutable parser/index.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaObservationData{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads and cross-checks all populations before any census output is opened. */
	NcrnaObservationData(final String directory, final String originals,
			final String contexts, final String queryFiles){
		originalTable=new Table(originals, NcrnaCmLocusAudit.HEADER);
		contextTable=new Table(contexts, NcrnaCmHitProjector.CONTEXT_HEADER);
		hitTable=new Table(path(directory, "hits.tsv"), HIT_HEADER);
		eventTable=new Table(path(directory, "events.tsv"), EVENT_HEADER);
		overlapTable=new Table(path(directory, "overlaps.tsv"), OVERLAP_HEADER);
		queryTable=new Table(path(directory, "queries.tsv"), QUERY_HEADER);
		neighborTable=new Table(path(directory, "neighbors.tsv"), NEIGHBOR_HEADER);
		searchTable=new Table(path(directory, "searches.tsv"), NcrnaCmHitProjector.TABLE_HEADER);
		queryFileTable=new Table(queryFiles, "table_id\tquery_fasta");
		loadOriginals();
		loadContexts();
		loadSearches();
		loadMembership();
		loadHits();
		loadEvents();
		loadOverlaps();
		checkQueries();
		checkHistoricalNeighbors();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Population joins       ----------------*/
	/*--------------------------------------------------------------*/

	/** Originals retain all flags and the exact source/assembly/contig identity. */
	private void loadOriginals(){
		assert(loci.isEmpty()) : "Original keys must be fixed before context or hit joins";
		final ByteBuilder key=new ByteBuilder();
		final IntHashMap sourceIndex=new IntHashMap();
		for(final byte[] line:originalTable.rows){
			parse(line, 21);
			final int id=symbol(0), genome=symbol(1), assembly=symbol(2), seqid=symbol(6);
			key.clear().append(symbols.bytes(genome)).tab().append(symbols.bytes(assembly)).tab().append(symbols.bytes(seqid));
			final int sourceKey=symbols.intern(key.array, 0, key.length());
			int source=sourceIndex.get(sourceKey);
			final int length=parser.parseInt(15);
			require(length>0, "Original contig length must be positive");
			if(source<0){
				source=sources.size();
				sourceIndex.put(sourceKey, source);
				sources.add(new Source(genome, assembly, seqid, length));
			}
			final Source contig=sources.get(source);
			require(contig.length==length, "Conflicting source contig lengths");
			parser.setBounds(8);
			int dash=parser.a();
			while(dash<parser.b() && line[dash]!='-'){dash++;}
			require(dash>parser.a() && dash<parser.b()-1, "Expected one contiguous original interval");
			final int start=positive(line, parser.a(), dash), stop=positive(line, dash+1, parser.b());
			require(start<=stop && stop<=length && (long)stop-start+1==parser.parseLong(9), "Original interval/length mismatch");
			final long locusKey=pair(genome, id);
			require(byLocus.get(locusKey)<0, "Duplicate original genome/ID");
			byLocus.put(locusKey, loci.size());
			loci.add(new Locus(loci.size(), id, contig, strand(7), start, stop));
		}
		require(!loci.isEmpty(), "The original population cannot be empty");
	}

	/** Actual flanks determine context bounds on either original strand. */
	private void loadContexts(){
		assert(!loci.isEmpty()) : "Every context must resolve to a retained original";
		for(final byte[] line:contextTable.rows){
			parse(line, 10);
			final int id=symbol(0);
			final Locus locus=locus(symbol(2), symbol(1));
			require(byQuery.get(id)<0, "Duplicate context ID");
			require(symbol(3)==locus.source.seqid && strand(4)==locus.strand
				&& parser.parseLong(5)==locus.start && parser.parseLong(6)==locus.stop
				&& parser.parseLong(7)==locus.source.length, "Context origin differs from original annotation");
			final int left=parser.parseInt(8), right=parser.parseInt(9);
			final long low=(long)locus.start-(locus.strand=='+' ? left : right);
			final long high=(long)locus.stop+(locus.strand=='+' ? right : left);
			require(left>=0 && right>=0 && low>=1 && high<=locus.source.length, "Actual context flanks exceed the contig");
			final Query query=new Query(queries.size(), id, locus, (int)low, (int)high, left, right);
			byQuery.put(id, query.index);
			queries.add(query);
			locus.queries.add(query.index);
			locus.source.queries.add(query);
		}
		require(!queries.isEmpty(), "The searched-context population cannot be empty");
	}

	/** Model/version is a pair; each search keeps its own query population and scores. */
	private void loadSearches(){
		assert(searches.isEmpty()) : "Search IDs are assigned before query membership";
		final LongIntMap byModel=new LongIntMap();
		for(final byte[] line:searchTable.rows){
			parse(line, 7);
			final int id=symbol(0), model=symbol(1), sha=symbol(2);
			final byte[] hash=symbols.bytes(sha);
			require(hash.length==20, "Model provenance requires a 20-character sha80");
			for(final byte b:hash){require((b>='0' && b<='9') || (b>='a' && b<='f'), "Invalid model sha80");}
			final int clen=parser.parseInt(3);
			final double ga=parser.parseDouble(4);
			require(clen>0 && Double.isFinite(ga), "Invalid model length or gathering threshold");
			int index=byModel.get(pair(model, sha));
			if(index<0){
				index=models.size();
				byModel.put(pair(model, sha), index);
				models.add(new Model(index, model, sha, clen, ga));
			}
			final Model definition=models.get(index);
			require(definition.clen==clen && definition.ga==ga, "Conflicting metadata for one model/version");
			require(bySearch.get(id)<0, "Duplicate search table ID");
			bySearch.put(id, searches.size());
			searches.add(new Search(searches.size(), id, definition, parser.parseString(6)));
		}
		require(!searches.isEmpty(), "At least one search table is required, including zero-hit tables");
	}

	/** Reads each distinct query FASTA once and binds it to actual search targets. */
	private void loadMembership(){
		assert(!searches.isEmpty()) : "Membership cannot invent undeclared search IDs";
		final ObjectIntMap<String> files=new ObjectIntMap<String>(String.class);
		final ArrayList<IntList> populations=new ArrayList<IntList>();
		for(final byte[] line:queryFileTable.rows){
			parse(line, 2);
			final Search search=search(symbol(0));
			require(search.members==null, "Duplicate query-file binding for a search");
			final String path=parser.parseString(1);
			checkTarget(search.table, path);
			int index=files.get(path);
			if(index<0){
				index=populations.size();
				files.put(path, index);
				populations.add(readMembership(path));
			}
			search.queryFile=path;
			search.members=populations.get(index);
			for(int i=0; i<search.members.size; i++){
				queries.get(search.members.get(i)).searches.add(search.index);
			}
		}
		for(final Search search:searches){require(search.members!=null, "Search missing its query FASTA binding");}
		for(final Query query:queries){require(query.searches.size>0, "Context has no actual search membership");}
	}

	/** Native FASTA reading; only retained context indices survive each read batch. */
	private IntList readMembership(final String path){
		final IntList members=new IntList();
		final boolean[] seen=new boolean[queries.size()];
		final ByteBuilder id=new ByteBuilder();
		final FileFormat format=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		require(format.fasta(), "Search membership requires the recorded query FASTA");
		// StreamerFactory's zero-thread FASTA branch avoids an asynchronous producer
		// surviving a consumer-side rejection; FastaStreamerST.close does not cancel it.
		final Streamer stream=StreamerFactory.makeStreamer(format, null, false, -1, false, true, 0);
		try(SequenceInput input=new SequenceInput(stream)){
			stream.start();
			for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
				for(final Read read:batch){
					id.clear();
					for(int i=0; i<read.id.length() && read.id.charAt(i)>' '; i++){
						final char c=read.id.charAt(i);
						require(c<128, "Native query identifiers must be ASCII; metadata bytes remain unchanged");
						id.append(c);
					}
					final int symbol=symbols.find(id.array, 0, id.length());
					final int index=symbol<0 ? -1 : byQuery.get(symbol);
					require(index>=0, "Query FASTA contains an ID absent from the accepted context map");
					require(!seen[index], "Repeated context ID within one query FASTA");
					final Query query=queries.get(index);
					require(read.length()==query.stop-query.start+1, "Query sequence length differs from actual context bounds");
					seen[index]=true;
					members.add(index);
				}
			}
		}
		require(members.size>0, "Search query FASTA is empty");
		members.sort();
		return members;
	}

	/** Confirms the manifest against Infernal's own recorded target, including zero-hit tables. */
	private static void checkTarget(final String table, final String fasta){
		int found=0;
		final String prefix="# Target file:";
		try(Input input=new Input(table)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				if(!startsWith(line, prefix)){continue;}
				int from=prefix.length(), to=line.length;
				while(from<to && line[from]<=' '){from++;}
				while(to>from && line[to-1]<=' '){to--;}
				final String recorded=new String(line, from, to-from, StandardCharsets.UTF_8);
				require(recorded.equals(fasta), "Query FASTA differs from recorded search target: "+table);
				found++;
			}
		}
		require(found==1, "Expected one recorded Infernal target file: "+table);
	}

	/** Physical IDs are dense projector IDs, checked rather than inferred from row order. */
	private void loadHits(){
		assert(hits.isEmpty()) : "Physical intervals load once before event joins";
		final LongIntMap byModel=new LongIntMap();
		for(final Model model:models){byModel.put(pair(model.id, model.sha), model.index);}
		final ByteBuilder key=new ByteBuilder();
		final IntHashMap bySource=new IntHashMap();
		for(int i=0; i<sources.size(); i++){
			final Source source=sources.get(i);
			key.clear().append(symbols.bytes(source.genome)).tab().append(symbols.bytes(source.assembly)).tab().append(symbols.bytes(source.seqid));
			bySource.put(symbols.intern(key.array, 0, key.length()), i);
		}
		for(final byte[] line:hitTable.rows){
			parse(line, 11);
			final int id=parser.parseInt(0), genome=symbol(1), assembly=symbol(2), seqid=symbol(3);
			key.clear().append(symbols.bytes(genome)).tab().append(symbols.bytes(assembly)).tab().append(symbols.bytes(seqid));
			final int source=bySource.get(symbols.intern(key.array, 0, key.length()));
			final int model=byModel.get(pair(symbol(5), symbol(6)));
			require(source>=0 && model>=0 && id>0, "Physical interval has unknown source/model or invalid ID");
			final int start=parser.parseInt(7), stop=parser.parseInt(8);
			final int events=parser.parseInt(9), overlaps=parser.parseInt(10);
			require(start>0 && stop>=start && stop<=sources.get(source).length && events>0 && overlaps>=0, "Invalid physical interval or counts");
			final Hit hit=new Hit(id, sources.get(source), strand(4), start, stop, models.get(model), events, overlaps);
			hits.add(hit);
			hit.source.hits.add(hit);
		}
		hits.sort((a, b)->Integer.compare(a.id, b.id));
		for(int i=0; i<hits.size(); i++){require(hits.get(i).id==i+1, "Projector physical IDs must be unique and dense from one");}
		for(final Source source:sources){
			source.hits.sort((a, b)->{
				int c=Character.compare(a.strand, b.strand);
				if(c==0){c=Integer.compare(a.model.index, b.model.index);}
				if(c==0){c=Integer.compare(a.start, b.start);}
				return c==0 ? Integer.compare(a.stop, b.stop) : c;
			});
			for(int i=1; i<source.hits.size(); i++){
				final Hit a=source.hits.get(i-1), b=source.hits.get(i);
				require(a.strand!=b.strand || a.model!=b.model || a.start!=b.start || a.stop!=b.stop,
					"Repeated physical source/strand/model/bounds tuple");
			}
		}
	}

	/** Keeps every raw row; checks origin membership, model and genomic projection. */
	private void loadEvents(){
		assert(!searches.isEmpty()) : "Raw events require declared searches even when no hit is present";
		final LongIntMap keys=new LongIntMap();
		final int[] projected=new int[3];
		for(final byte[] line:eventTable.rows){
			parse(line, 15);
			final Hit hit=hit(parser.parseInt(0));
			final Search search=search(symbol(1));
			final int tableLine=parser.parseInt(2);
			final Query query=query(symbol(3));
			require(tableLine>0 && keys.get(pair(search.index, tableLine))<0, "Repeated or invalid raw event table/line");
			keys.put(pair(search.index, tableLine), 1);
			require(Arrays.binarySearch(search.members.array, 0, search.members.size, query.index)>=0, "Hit context was not in its search query FASTA");
			require(symbol(4)==query.origin.source.genome && symbol(5)==query.origin.id, "Event search origin differs from context map");
			require(hit.source==query.origin.source && hit.model==search.model, "Event model/source differs from its physical interval");
			final int modelFrom=parser.parseInt(6), modelTo=parser.parseInt(7);
			require(modelFrom>=1 && modelTo>=modelFrom && modelTo<=search.model.clen, "Raw model bounds exceed declared CLEN");
			NcrnaCmHitProjector.project(query.start, query.stop, query.origin.strand,
				parser.parseInt(8), parser.parseInt(9), strand(10), projected);
			require(projected[0]==hit.start && projected[1]==hit.stop && projected[2]==hit.strand, "Event bounds differ from physical interval");
			hit.eventRows.add(eventRows++);
			eventQueries.add(query.index);
			query.events++;
			search.events++;
			final int previous=search.eventCounts.get(query.index);
			search.eventCounts.put(query.index, previous<0 ? 1 : previous+1);
		}
		for(final Hit hit:hits){require(hit.eventRows.size==hit.expectedEvents, "Physical event count does not conserve raw observations");}
	}

	/** All annotation links are retained separately from the search-origin links. */
	private void loadOverlaps(){
		assert(eventRows==eventTable.rows.size()) : "Overlap joins follow complete raw-event conservation";
		final LongIntMap pairs=new LongIntMap();
		for(final byte[] line:overlapTable.rows){
			parse(line, 6);
			final Hit hit=hit(parser.parseInt(0));
			final Locus locus=locus(symbol(1), symbol(2));
			final int bases=Math.min(hit.stop, locus.stop)-Math.max(hit.start, locus.start)+1;
			require(hit.source==locus.source && bases>0 && parser.parseInt(3)==bases, "Annotation overlap geometry differs from physical interval");
			require(bool(4)==(hit.strand==locus.strand) && bool(5)==(locus.queries.size>0), "Annotation overlap strand/query state differs");
			final long key=pair(hit.id, locus.index);
			require(pairs.get(key)<0, "Repeated physical/annotation overlap pair");
			pairs.put(key, 1);
			locus.hits.add(hit.id);
			hit.overlaps++;
		}
		for(final Hit hit:hits){require(hit.overlaps==hit.expectedOverlaps, "Physical annotation-link count differs from retained overlap rows");}
	}

	/** The accepted projector's zero-event rows must survive and agree with the joins. */
	private void checkQueries(){
		final boolean[] seen=new boolean[queries.size()];
		for(final byte[] line:queryTable.rows){
			parse(line, 11);
			final Query query=query(symbol(0));
			require(!seen[query.index], "Repeated projector query row");
			seen[query.index]=true;
			require(symbol(1)==query.origin.source.genome && symbol(2)==query.origin.id
				&& parser.parseInt(3)==query.start && parser.parseInt(4)==query.stop
				&& parser.parseInt(5)==query.left && parser.parseInt(6)==query.right
				&& bool(7)==(query.start==1) && bool(8)==(query.stop==query.origin.source.length)
				&& parser.parseInt(9)==query.events && parser.parseInt(10)==query.origin.hits.size,
				"Projector query row differs from complete input joins");
		}
		require(queryTable.rows.size()==queries.size(), "Projector queries do not conserve context IDs");
	}

	/** Historical neighbor rows are validated and copied, never replaced by hit neighbors. */
	private void checkHistoricalNeighbors(){
		assert(!loci.isEmpty()) : "Historical neighbor identity is source-local";
		for(final byte[] line:neighborTable.rows){
			parse(line, 9);
			final Query query=query(symbol(0));
			final Locus a=locus(symbol(1), symbol(2)), b=locus(symbol(1), symbol(4));
			require(query.origin==a && a.source==b.source && a!=b, "Invalid historical neighbor source/origin");
			final int overlap=Math.max(0, Math.min(a.stop, b.stop)-Math.max(a.start, b.start)+1);
			require(parser.parseInt(6)==overlap && bool(7)==(a.strand==b.strand)
				&& bool(8)==(b.queries.size>0), "Historical neighbor geometry/flags differ");
			if(parser.termEquals("OVERLAP", 3)){
				require(overlap>0 && parser.termEquals("NA", 5), "Invalid historical overlap relationship");
			}else{
				final boolean previous=parser.termEquals("PREDECESSOR", 3);
				require(previous || parser.termEquals("SUCCESSOR", 3), "Unknown historical neighbor relation");
				final long gap=previous ? (long)a.start-b.stop-1 : (long)b.start-a.stop-1;
				require(a.strand==b.strand && parser.parseLong(5)==gap, "Historical neighbor gap or strand differs");
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Parsing utilities      ----------------*/
	/*--------------------------------------------------------------*/

	private void parse(final byte[] line, final int fields){
		parser.set(line);
		current=line;
		require(parser.terms()==fields, "Unexpected census input field count");
	}

	private int symbol(final int field){
		parser.setBounds(field);
		require(parser.a()<parser.b(), "Identity fields cannot be empty");
		return symbols.intern(current, parser.a(), parser.b());
	}

	private char strand(final int field){
		if(parser.termEquals('+', field)){return '+';}
		if(parser.termEquals('-', field)){return '-';}
		throw new IllegalArgumentException("Expected explicit + or - strand");
	}

	private boolean bool(final int field){
		if(parser.termEquals("true", field)){return true;}
		if(parser.termEquals("false", field)){return false;}
		throw new IllegalArgumentException("Expected true or false field");
	}

	private Locus locus(final int genome, final int id){
		final int index=byLocus.get(pair(genome, id));
		require(index>=0, "Unknown original genome/ID");
		return loci.get(index);
	}

	private Query query(final int id){
		final int index=byQuery.get(id);
		require(index>=0, "Unknown context ID");
		return queries.get(index);
	}

	private Search search(final int id){
		final int index=bySearch.get(id);
		require(index>=0, "Unknown search table ID");
		return searches.get(index);
	}

	Hit hit(final int id){
		require(id>0 && id<=hits.size(), "Physical ID is outside the accepted population");
		return hits.get(id-1);
	}

	/** Event-origin lookup uses retained primitive indices, without reparsing identity strings. */
	Query eventQuery(final int row){return queries.get(eventQueries.get(row));}

	/** Resolves a source-qualified original from two borrowed TSV fields. */
	Locus originalByFields(final byte[] line, final LineParser1 fields, final int genome, final int id){
		assert(line!=null && fields!=null) : "Historical relationship lookup borrows the caller's parsed row";
		fields.setBounds(genome);
		final int g=symbols.find(line, fields.a(), fields.b());
		fields.setBounds(id);
		final int key=symbols.find(line, fields.a(), fields.b());
		require(g>=0 && key>=0, "Unknown historical relationship identity");
		return locus(g, key);
	}

	static long pair(final int a, final int b){return ((long)a<<32)|(b&0xffffffffL);}
	static String path(final String directory, final String name){return new File(directory, name).getPath();}
	static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}

	private static int positive(final byte[] line, final int from, final int to){
		long value=0;
		for(int i=from; i<to; i++){
			require(line[i]>='0' && line[i]<='9', "Original interval must be contiguous decimal start-stop");
			value=value*10+line[i]-'0';
			require(value<=Integer.MAX_VALUE, "Original interval coordinate exceeds integer range");
		}
		require(value>0, "Original interval coordinates are one-based");
		return (int)value;
	}

	private static boolean startsWith(final byte[] line, final String prefix){
		if(line.length<prefix.length()){return false;}
		for(int i=0; i<prefix.length(); i++){if(line[i]!=prefix.charAt(i)){return false;}}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Data owners         ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns unchanged rows, including flags and raw floating-point text. */
	static final class Table{
		Table(final String path_, final String header_){
			path=path_;
			header=header_;
			try(Input input=new Input(path)){
				final byte[] first=input.reader.nextLine();
				require(first!=null && first.length==header.length() && startsWith(first, header), "Unexpected header: "+path);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					require(line.length>0, "Empty metadata row: "+path);
					rows.add(line);
				}
			}
		}
		final String path, header;
		final ArrayList<byte[]> rows=new ArrayList<byte[]>();
	}

	/** A shared full source identity, not a bare accession or contig name. */
	static final class Source{
		Source(final int g, final int a, final int s, final int length_){
			genome=g;
			assembly=a;
			seqid=s;
			length=length_;
		}
		/** Prefix maximum ends permit bounded overlap lookup even with nested contexts. */
		void indexQueries(){
			queries.sort((a, b)->Integer.compare(a.start, b.start));
			maxQueryEnd=new int[queries.size()];
			int end=0;
			for(int i=0; i<queries.size(); i++){
				end=Math.max(end, queries.get(i).stop);
				maxQueryEnd[i]=end;
			}
			assert(maxQueryEnd.length==queries.size()) : "Context index has one bound per retained query";
		}
		/** First query whose prefix reaches start; caller still tests each actual interval. */
		int firstQuery(final int start){
			assert(maxQueryEnd!=null) : "Coverage lookup follows source context indexing";
			int low=0, high=maxQueryEnd.length;
			while(low<high){
				final int middle=(low+high)>>>1;
				if(maxQueryEnd[middle]<start){low=middle+1;}else{high=middle;}
			}
			return low;
		}
		final int genome, assembly, seqid, length;
		final ArrayList<Hit> hits=new ArrayList<Hit>();
		final ArrayList<Query> queries=new ArrayList<Query>();
		private int[] maxQueryEnd;
	}

	/** Original metadata is retained in originalTable at index, including all flags. */
	static final class Locus{
		Locus(final int index_, final int id_, final Source source_, final char strand_, final int start_, final int stop_){
			index=index_;
			id=id_;
			source=source_;
			strand=strand_;
			start=start_;
			stop=stop_;
		}
		final int index, id, start, stop;
		final Source source;
		final char strand;
		final IntList queries=new IntList(1), hits=new IntList(1);
	}

	/** Actual search window, with its full per-table membership. */
	static final class Query{
		Query(final int index_, final int id_, final Locus origin_, final int start_, final int stop_, final int left_, final int right_){
			index=index_;
			id=id_;
			origin=origin_;
			start=start_;
			stop=stop_;
			left=left_;
			right=right_;
		}
		boolean searched(final Model model, final ArrayList<Search> tables){
			assert(model!=null) : "Coverage queries name a declared model/version";
			for(int i=0; i<searches.size; i++){if(tables.get(searches.get(i)).model==model){return true;}}
			return false;
		}
		final int index, id, start, stop, left, right;
		final Locus origin;
		final IntList searches=new IntList(2);
		int events;
	}

	static final class Model{
		Model(final int index_, final int id_, final int sha_, final int clen_, final double ga_){
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
		Search(final int index_, final int id_, final Model model_, final String table_){
			index=index_;
			id=id_;
			model=model_;
			table=table_;
		}
		final int index, id;
		final Model model;
		final String table;
		String queryFile;
		IntList members;
		int events;
		final IntHashMap eventCounts=new IntHashMap();
	}

	/** Physical bounds never change during census; counts refer to observations only. */
	static final class Hit{
		Hit(final int id_, final Source source_, final char strand_, final int start_, final int stop_, final Model model_, final int events_, final int overlaps_){
			id=id_;
			source=source_;
			strand=strand_;
			start=start_;
			stop=stop_;
			model=model_;
			expectedEvents=events_;
			expectedOverlaps=overlaps_;
		}
		final int id, start, stop, expectedEvents, expectedOverlaps;
		final Source source;
		final char strand;
		final Model model;
		final IntList eventRows=new IntList(1);
		int overlaps;
	}

	/** Reusable byte-slice probe plus immutable copies only for new symbols. */
	static final class Symbols{
		int intern(final byte[] bytes, final int from, final int to){
			probe.set(bytes, from, to);
			int id=index.get(probe);
			if(id<0){
				id=values.size();
				final Key stored=new Key();
				stored.set(Arrays.copyOfRange(bytes, from, to), 0, to-from);
				values.add(stored.bytes);
				index.put(stored, id);
			}
			return id;
		}
		int find(final byte[] bytes, final int from, final int to){
			probe.set(bytes, from, to);
			return index.get(probe);
		}
		byte[] bytes(final int id){return values.get(id);}
		private final ObjectIntMap<Key> index=new ObjectIntMap<Key>(Key.class);
		private final ArrayList<byte[]> values=new ArrayList<byte[]>();
		private final Key probe=new Key();
	}

	/** Content equality for the native map; a stored key is never mutated. */
	private static final class Key{
		void set(final byte[] bytes_, final int from_, final int to_){
			assert(bytes_!=null && from_>=0 && to_>=from_ && to_<=bytes_.length) : "Symbol lookup borrows one validated field slice";
			bytes=bytes_;
			from=from_;
			to=to_;
			hash=1;
			for(int i=from; i<to; i++){hash=31*hash+(bytes[i]&255);}
		}
		@Override
		public int hashCode(){return hash;}
		@Override
		public boolean equals(final Object other){
			if(!(other instanceof Key)){return false;}
			final Key key=(Key)other;
			if(to-from!=key.to-key.from){return false;}
			for(int i=0; i<to-from; i++){if(bytes[from+i]!=key.bytes[key.from+i]){return false;}}
			return true;
		}
		private byte[] bytes;
		private int from, to, hash;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final Symbols symbols=new Symbols();
	final ArrayList<Source> sources=new ArrayList<Source>();
	final ArrayList<Locus> loci=new ArrayList<Locus>();
	final ArrayList<Query> queries=new ArrayList<Query>();
	final ArrayList<Model> models=new ArrayList<Model>();
	final ArrayList<Search> searches=new ArrayList<Search>();
	final ArrayList<Hit> hits=new ArrayList<Hit>();
	final Table originalTable, contextTable, hitTable, eventTable, overlapTable;
	final Table queryTable, neighborTable, searchTable, queryFileTable;
	private final LineParser1 parser=new LineParser1('\t');
	private final LongIntMap byLocus=new LongIntMap();
	private final IntHashMap byQuery=new IntHashMap(), bySearch=new IntHashMap();
	private byte[] current;
	private int eventRows;
	private final IntList eventQueries=new IntList();
	static final String HIT_HEADER="physical_id\tgenome\tassembly\tseqid\tstrand\tmodel_id\tmodel_sha80\tstart1\tstop1\tevents\toverlapping_annotations";
	static final String EVENT_HEADER="physical_id\ttable_id\ttable_line\tcontext_id\torigin_genome\torigin_id\tmodel_from\tmodel_to\tlocal_from\tlocal_to\tlocal_strand\ttruncation\tinclusion\tbit_score\tevalue";
	static final String OVERLAP_HEADER="physical_id\tgenome\toriginal_id\toverlap_bases\tsame_strand\tannotation_queried";
	static final String QUERY_HEADER="context_id\tgenome\toriginal_id\tcontext_start1\tcontext_stop1\tlflank\trflank\ttouches_contig_left\ttouches_contig_right\torigin_events\tannotation_overlapping_hit_intervals";
	static final String NEIGHBOR_HEADER="context_id\tgenome\toriginal_id\trelationship\tneighbor_id\tsigned_gap\toverlap_bases\tsame_strand\tneighbor_queried";
}
