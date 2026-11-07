package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Projects every saved Infernal hit without assigning biological admission.
 * Forked from template.A_SampleByteFile. One metadata worker owns the joins;
 * native byte readers and writers own compressed I/O. Historical annotations,
 * queried contexts, raw events and physical hit intervals remain separate axes.
 * All coordinates are one-based inclusive. Only contiguous historical intervals
 * are supported; multipart inputs fail instead of being flattened silently.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCmHitProjector{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs projection, preserving the diagnostic stream's ownership on failure. */
	public static void main(final String[] args){
		final Timer timer=new Timer();
		final NcrnaCmHitProjector tool=new NcrnaCmHitProjector(args);
		try{
			tool.process(timer);
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/** Parses native options; complete named input populations cannot be sampled. */
	public NcrnaCmHitProjector(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		try{
			final Parser parser=parse(args);
			in=Tools.fixExtension(parser.in1);
			out=parser.out1;
			if(parser.append || parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1){
				throw new IllegalArgumentException("Projection requires complete unpaired metadata inputs, without sampling or append");
			}
			if(in==null || contexts==null || tables==null || out==null){
				throw new IllegalArgumentException("Required: in=historical.tsv contexts=crosswalk.tsv tables=searches.tsv out=fresh_directory");
			}
			if(!Tools.testInputFiles(false, true, in, contexts, tables)){
				throw new IllegalArgumentException("Cannot read historical annotations, context map or search manifest");
			}
			if(new File(out).exists()){throw new IllegalArgumentException("Output directory already exists: "+out);}
			ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
			ReadWrite.setZipThreads(Shared.threads());
		}catch(RuntimeException | Error failure){
			Shared.closeStream(outstream);
			throw failure;
		}
	}

	/** Standard template a/b parsing; expected counts are optional independent checks. */
	private Parser parse(final String[] args){
		assert(args!=null) : "PreParser supplies the complete native argument array";
		final Parser parser=new Parser();
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("contexts")){
				contexts=b;
			}else if(a.equals("tables")){
				tables=b;
			}else if(a.equals("expectedloci")){
				expectedLoci=expected(a, b);
			}else if(a.equals("expectedcontexts")){
				expectedContexts=expected(a, b);
			}else if(a.equals("expectedhits")){
				expectedHits=expected(a, b);
			}else if(parser.parse(arg, a, b)){
				// Standard BBTools options.
			}else{
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		return parser;
	}

	private static int expected(final String name, final String value){
		final long count=Parse.parseKMG(value);
		if(count<0 || count>Integer.MAX_VALUE){throw new IllegalArgumentException(name+" must be a nonnegative integer");}
		return (int)count;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Loading            ----------------*/
	/*--------------------------------------------------------------*/

	/** Completes all input validation before creating the output directory. */
	private void process(final Timer timer){
		loadLoci();
		loadContexts();
		loadHits();
		indexHits();
		write();
		timer.stop();
		outstream.println("Historical annotations: "+loci.size());
		outstream.println("Queried contexts: "+queries.size());
		outstream.println("Raw hit events: "+events.size());
		outstream.println("Physical hit intervals: "+physical.size());
		outstream.println("No subtype, arrangement, admission or biological-negative labels assigned.");
		outstream.println("Time: "+timer);
	}

	/** Reads the complete annotated-corpus schema, preserving exact assembly keys. */
	private void loadLoci(){
		assert(loci.isEmpty()) : "Historical identities must be loaded once before context and hit joins";
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(in)){
			header(input.reader, NcrnaCmLocusAudit.HEADER, in);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.terms()!=21){throw new IllegalArgumentException("Expected 21 historical fields in "+in);}
				final String id=parser.parseString(0), genome=parser.parseString(1);
				final String assembly=parser.parseString(2), seqid=parser.parseString(6);
				require(!id.isEmpty() && !genome.isEmpty() && !assembly.isEmpty() && !seqid.isEmpty(), "Historical identity fields cannot be empty");
				final char strand=strand(parser, 7);
				parser.setBounds(8);
				final int from=parser.a(), to=parser.b();
				int dash=from;
				while(dash<to && line[dash]!='-'){dash++;}
				if(dash<=from || dash>=to-1){throw new IllegalArgumentException("Expected one contiguous start-stop interval for "+id);}
				final int start=positive(line, from, dash), stop=positive(line, dash+1, to);
				final int length=parser.parseInt(15);
				if(start>stop || stop>length || (long)stop-start+1!=parser.parseLong(9)){
					throw new IllegalArgumentException("Historical interval/length mismatch for "+id);
				}
				final String key=genome+'\t'+id;
				if(byLocus.contains(key)){throw new IllegalArgumentException("Duplicate historical source/ID: "+key);}
				final String contigKey=genome+'\t'+seqid;
				Contig contig=byContig.get(contigKey);
				if(contig==null){
					contig=new Contig(genome, assembly, seqid, length);
					byContig.put(contigKey, contig);
					contigs.add(contig);
				}
				if(contig.length!=length || !contig.assembly.equals(assembly)){
					throw new IllegalArgumentException("Conflicting assembly/contig length for "+contigKey);
				}
				final Locus locus=new Locus(id, contig, strand, start, stop);
				byLocus.put(key, locus);
				loci.add(locus);
				contig.loci.add(locus);
			}
		}
		checkCount("historical annotations", loci.size(), expectedLoci);
		require(!loci.isEmpty(), "The complete historical population cannot be empty");
		for(final Contig contig:contigs){contig.index();}
	}

	/** Joins each queried context to an immutable original feature; neighbors may be unqueried. */
	private void loadContexts(){
		assert(!loci.isEmpty()) : "Context origins must resolve against the complete historical population";
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(contexts)){
			header(input.reader, CONTEXT_HEADER, contexts);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.terms()!=10){throw new IllegalArgumentException("Expected ten context fields in "+contexts);}
				final String id=parser.parseString(0);
				final String original=parser.parseString(1), genome=parser.parseString(2);
				final Locus locus=byLocus.get(genome+'\t'+original);
				if(locus==null){throw new IllegalArgumentException("Missing historical origin for context "+id);}
				if(id.isEmpty() || byQuery.contains(id)){throw new IllegalArgumentException("Context IDs must be nonempty and globally unique: "+id);}
				if(!parser.termEquals(locus.contig.seqid, 3) || strand(parser, 4)!=locus.strand
						|| parser.parseLong(5)!=locus.start || parser.parseLong(6)!=locus.stop
						|| parser.parseLong(7)!=locus.contig.length){
					throw new IllegalArgumentException("Context disagrees with historical origin: "+id);
				}
				final long left=parser.parseLong(8), right=parser.parseLong(9);
				if(left<0 || right<0 || left>Integer.MAX_VALUE || right>Integer.MAX_VALUE){throw new IllegalArgumentException("Invalid actual flanks: "+id);}
				final long low=(long)locus.start-(locus.strand=='+' ? left : right);
				final long high=(long)locus.stop+(locus.strand=='+' ? right : left);
				if(low<1 || high>locus.contig.length){throw new IllegalArgumentException("Context exceeds its original contig: "+id);}
				final Query query=new Query(id, locus, (int)low, (int)high, (int)left, (int)right);
				queries.add(query);
				byQuery.put(id, query);
				locus.contexts++;
			}
		}
		checkCount("queried contexts", queries.size(), expectedContexts);
		require(!queries.isEmpty(), "At least one queried context is required");
		Collections.sort(queries, (a, b)->a.id.compareTo(b.id));
	}

	/** Reuses the all-hit Infernal reader without filtering its qualification status. */
	private void loadHits(){
		assert(!queries.isEmpty()) : "Every raw hit target must resolve to a declared queried context";
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		final ObjectMap<String,Model> models=new ObjectMap<String,Model>(String.class, Model.class);
		final LineParser1 parser=new LineParser1('\t');
		final int[] projected=new int[3];
		try(Input input=new Input(tables)){
			header(input.reader, TABLE_HEADER, tables);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.terms()!=7){throw new IllegalArgumentException("Expected seven search-manifest fields in "+tables);}
				final String id=parser.parseString(0), accession=parser.parseString(1), hash=parser.parseString(2);
				final int clen=parser.parseInt(3);
				final double ga=parser.parseDouble(4);
				final String arm=parser.parseString(5), path=parser.parseString(6);
				require(!id.isEmpty() && ids.add(id) && !arm.isEmpty() && !path.isEmpty() && paths.add(path), "Empty or repeated search identity/path: "+id);
				require(hash.matches("[0-9a-f]{20}"), "Model identity requires declared sha80: "+id);
				Model model=models.get(accession);
				if(model==null){
					model=new Model(accession, hash, clen, ga);
					models.put(accession, model);
				}
				require(model.hash.equals(hash) && model.clen==clen && model.ga==ga, "Conflicting model version/metadata: "+accession);
				final Search search=new Search(id, arm, path, model);
				searches.add(search);
				final RrnaCmTruthBuilder.Genome source=new RrnaCmTruthBuilder.Genome(arm, id, path, null);
				for(final RrnaCmTruthBuilder.Hit hit:RrnaCmTruthBuilder.readHits(source, model.criteria)){
					final Query query=byQuery.get(hit.seqid);
					if(query==null){throw new IllegalArgumentException("Unknown context in "+path+":"+hit.line+": "+hit.seqid);}
					project(query.low, query.high, query.locus.strand, hit.seqFrom, hit.seqTo, hit.strand, projected);
					events.add(new Event(search, query, hit, projected[0], projected[1], (char)projected[2]));
					query.events++;
				}
			}
		}
		require(!searches.isEmpty(), "Search manifest cannot be empty");
		checkCount("raw hit events", events.size(), expectedHits);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Geometry and index     ----------------*/
	/*--------------------------------------------------------------*/

	/** Maps an oriented one-based local hit into genomic coordinates and strand.
	 * The reusable result holds start, stop and the genomic strand character.
	 * Public for finite independent geometry fixtures; no allocation per call.
	 */
	public static void project(final int low, final int high, final char orientation,
			final int from, final int to, final char hitStrand, final int[] result){
		require(result!=null && result.length>=3, "Projection needs a three-element result buffer");
		require(low>=1 && high>=low && (orientation=='+' || orientation=='-'), "Invalid genomic context bounds/orientation");
		require((hitStrand=='+' && from<=to) || (hitStrand=='-' && from>=to), "Hit strand disagrees with local coordinate order");
		final int first=Math.min(from, to), last=Math.max(from, to);
		require(first>=1 && (long)last<=(long)high-low+1, "Hit exceeds the declared context length");
		result[0]=orientation=='+' ? low+first-1 : high-last+1;
		result[1]=orientation=='+' ? low+last-1 : high-first+1;
		result[2]=orientation=='+' ? hitStrand : hitStrand=='+' ? '-' : '+';
		assert(result[0]>=low && result[1]<=high && result[0]<=result[1]) : "Projection must preserve a bounded inclusive interval";
	}

	/** Exact physical keys collapse observations only; every event keeps its own scores. */
	private void indexHits(){
		Collections.sort(events, (a, b)->{
			int compare=comparePhysical(a, b);
			if(compare==0){compare=a.search.id.compareTo(b.search.id);}
			return compare==0 ? Integer.compare(a.hit.line, b.hit.line) : compare;
		});
		Physical current=null;
		for(final Event event:events){
			if(current==null || comparePhysical(current.first, event)!=0){
				current=new Physical(physical.size()+1, event);
				physical.add(current);
			}
			event.physical=current;
			current.events++;
		}
		assert(events.isEmpty()==physical.isEmpty()) : "Every raw event must retain exactly one physical interval identity";
	}

	private static int comparePhysical(final Event a, final Event b){
		assert(a!=null && b!=null) : "Physical ordering compares retained validated hit events";
		final Contig x=a.query.locus.contig, y=b.query.locus.contig;
		int compare=x.genome.compareTo(y.genome);
		if(compare==0){compare=x.assembly.compareTo(y.assembly);}
		if(compare==0){compare=x.seqid.compareTo(y.seqid);}
		if(compare==0){compare=Character.compare(a.strand, b.strand);}
		if(compare==0){compare=a.search.model.id.compareTo(b.search.model.id);}
		if(compare==0){compare=a.search.model.hash.compareTo(b.search.model.hash);}
		if(compare==0){compare=Integer.compare(a.start, b.start);}
		return compare==0 ? Integer.compare(a.stop, b.stop) : compare;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes separate event, physical, overlap, query, neighbor and provenance tables. */
	private void write(){
		assert(!queries.isEmpty() && !searches.isEmpty()) : "Output must retain an explicitly queried population and search provenance, even with zero events";
		final File directory=new File(out);
		if(!directory.mkdir()){throw new IllegalArgumentException("Could not create fresh output directory: "+out);}
		final ByteBuilder row=new ByteBuilder();
		try(Output hits=output(directory, "hits.tsv"); Output overlaps=output(directory, "overlaps.tsv")){
			hits.writer.println("physical_id\tgenome\tassembly\tseqid\tstrand\tmodel_id\tmodel_sha80\tstart1\tstop1\tevents\toverlapping_annotations");
			overlaps.writer.println("physical_id\tgenome\toriginal_id\toverlap_bases\tsame_strand\tannotation_queried");
			for(final Physical item:physical){
				final Event event=item.first;
				final Contig contig=event.query.locus.contig;
				int linked=0;
				for(int i=contig.firstOverlap(event.start); i<contig.loci.size(); i++){
					final Locus locus=contig.loci.get(i);
					if(locus.start>event.stop){break;}
					if(locus.stop<event.start){continue;}
					final int bases=Math.min(locus.stop, event.stop)-Math.max(locus.start, event.start)+1;
					locus.overlappingHits++;
					linked++;
					overlaps.writer.print(row.clear().append(item.id).tab().append(contig.genome).tab().append(locus.id)
						.tab().append(bases).tab().append(locus.strand==event.strand).tab().append(locus.contexts>0).nl());
				}
				hits.writer.print(row.clear().append(item.id).tab().append(contig.genome).tab().append(contig.assembly)
					.tab().append(contig.seqid).tab().append(event.strand).tab().append(event.search.model.id)
					.tab().append(event.search.model.hash).tab().append(event.start).tab().append(event.stop)
					.tab().append(item.events).tab().append(linked).nl());
			}
		}
		try(Output output=output(directory, "events.tsv")){
			output.writer.println("physical_id\ttable_id\ttable_line\tcontext_id\torigin_genome\torigin_id\tmodel_from\tmodel_to\tlocal_from\tlocal_to\tlocal_strand\ttruncation\tinclusion\tbit_score\tevalue");
			for(final Event event:events){
				final RrnaCmTruthBuilder.Hit hit=event.hit;
				output.writer.print(row.clear().append(event.physical.id).tab().append(event.search.id).tab().append(hit.line)
					.tab().append(event.query.id).tab().append(event.query.locus.contig.genome).tab().append(event.query.locus.id)
					.tab().append(hit.mdlFrom).tab().append(hit.mdlTo).tab().append(hit.seqFrom).tab().append(hit.seqTo)
					.tab().append(hit.strand).tab().append(hit.trunc).tab().append(hit.included ? '!' : '?')
					.tab().append(hit.score).tab().append(hit.evalue).nl());
			}
		}
		writeQueries(directory, row);
		Collections.sort(searches, (a, b)->a.id.compareTo(b.id));
		try(Output output=output(directory, "searches.tsv")){
			output.writer.println(TABLE_HEADER);
			for(final Search search:searches){
				output.writer.print(row.clear().append(search.id).tab().append(search.model.id).tab().append(search.model.hash)
					.tab().append(search.model.clen).tab().append(search.model.gaText).tab().append(search.arm).tab().append(search.path).nl());
			}
		}
		try(Output output=output(directory, "summary.tsv")){
			output.writer.println("historical_annotations\tqueried_contexts\tsearch_tables\traw_events\tphysical_intervals");
			output.writer.print(row.clear().append(loci.size()).tab().append(queries.size()).tab().append(searches.size())
				.tab().append(events.size()).tab().append(physical.size()).nl());
		}
	}

	/** Retains zero-event queries and distinguishes unqueried historical neighbors. */
	private void writeQueries(final File directory, final ByteBuilder row){
		assert(row!=null) : "One output buffer is reused after each writer copies its current bytes";
		try(Output output=output(directory, "queries.tsv"); Output neighbors=output(directory, "neighbors.tsv")){
			output.writer.println("context_id\tgenome\toriginal_id\tcontext_start1\tcontext_stop1\tlflank\trflank\ttouches_contig_left\ttouches_contig_right\torigin_events\tannotation_overlapping_hit_intervals");
			neighbors.writer.println("context_id\tgenome\toriginal_id\trelationship\tneighbor_id\tsigned_gap\toverlap_bases\tsame_strand\tneighbor_queried");
			for(final Query query:queries){
				final Locus locus=query.locus;
				output.writer.print(row.clear().append(query.id).tab().append(locus.contig.genome).tab().append(locus.id)
					.tab().append(query.low).tab().append(query.high).tab().append(query.left).tab().append(query.right)
					.tab().append(query.low==1).tab().append(query.high==locus.contig.length).tab().append(query.events)
					.tab().append(locus.overlappingHits).nl());
				if(locus.previous!=null){neighbor(neighbors, row, query, locus.previous, "PREDECESSOR", (long)locus.start-locus.previous.stop-1);}
				if(locus.next!=null){neighbor(neighbors, row, query, locus.next, "SUCCESSOR", (long)locus.next.start-locus.stop-1);}
				for(int i=locus.contig.firstOverlap(locus.start); i<locus.contig.loci.size(); i++){
					final Locus other=locus.contig.loci.get(i);
					if(other.start>locus.stop){break;}
					if(other!=locus && other.stop>=locus.start){neighbor(neighbors, row, query, other, "OVERLAP", Long.MIN_VALUE);}
				}
			}
		}
	}

	private static void neighbor(final Output output, final ByteBuilder row, final Query query,
			final Locus other, final String relation, final long gap){
		final Locus origin=query.locus;
		assert(origin.contig==other.contig) : "Neighbor relations cannot cross exact source contigs";
		row.clear().append(query.id).tab().append(origin.contig.genome).tab().append(origin.id).tab().append(relation).tab().append(other.id).tab();
		if(gap==Long.MIN_VALUE){row.append("NA");}else{row.append(gap);}
		row.tab().append(Math.max(0, Math.min(origin.stop, other.stop)-Math.max(origin.start, other.start)+1))
			.tab().append(origin.strand==other.strand).tab().append(other.contexts>0).nl();
		output.writer.print(row);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Utilities           ----------------*/
	/*--------------------------------------------------------------*/

	private static Output output(final File directory, final String name){
		return new Output(FileFormat.testOutput(new File(directory, name).getPath(), FileFormat.TEXT, null, true, false, false, false));
	}

	private static void header(final ByteFile input, final String expected, final String path){
		assert(input!=null && expected!=null) : "Schema validation requires the native reader and a literal header";
		final byte[] line=input.nextLine();
		if(line==null || line.length!=expected.length()){throw new IllegalArgumentException("Unexpected header in "+path);}
		for(int i=0; i<line.length; i++){
			if(line[i]!=expected.charAt(i)){throw new IllegalArgumentException("Unexpected header in "+path);}
		}
	}

	private static char strand(final LineParser1 parser, final int column){
		if(parser.termEquals("+", column)){return '+';}
		if(parser.termEquals("-", column)){return '-';}
		throw new IllegalArgumentException("Expected explicit + or - strand");
	}

	private static int positive(final byte[] line, final int from, final int to){
		assert(line!=null && from>=0 && to<=line.length) : "Coordinate parsing borrows a validated field slice";
		require(from<to, "Empty coordinate token");
		long value=0;
		for(int i=from; i<to; i++){
			require(line[i]>='0' && line[i]<='9', "Expected one contiguous numeric interval");
			value=value*10+line[i]-'0';
			require(value<=Integer.MAX_VALUE, "Coordinate exceeds supported integer range");
		}
		require(value>0, "Coordinates are one-based positive integers");
		return (int)value;
	}

	private static void checkCount(final String name, final int actual, final int expected){
		if(expected>=0 && actual!=expected){throw new IllegalArgumentException(name+": "+actual+" != expected "+expected);}
	}

	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Retained data        ----------------*/
	/*--------------------------------------------------------------*/

	/** Sorted starts and prefix maximum ends support bounded overlap lookup, including nesting. */
	private static final class Contig{
		Contig(final String genome_, final String assembly_, final String seqid_, final int length_){
			genome=genome_;
			assembly=assembly_;
			seqid=seqid_;
			length=length_;
		}
		void index(){
			Collections.sort(loci, (a, b)->{
				int c=Integer.compare(a.start, b.start);
				if(c==0){c=Integer.compare(a.stop, b.stop);}
				return c==0 ? a.id.compareTo(b.id) : c;
			});
			maxEnd=new int[loci.size()];
			final Locus[] previous=new Locus[2];
			for(int i=0; i<loci.size(); i++){
				final Locus locus=loci.get(i);
				final int strand=locus.strand=='+' ? 0 : 1;
				locus.previous=previous[strand];
				if(previous[strand]!=null){previous[strand].next=locus;}
				previous[strand]=locus;
				maxEnd[i]=Math.max(i==0 ? 0 : maxEnd[i-1], locus.stop);
			}
			assert(maxEnd.length==loci.size()) : "Overlap search requires one prefix bound per sorted historical interval";
		}
		int firstOverlap(final int start){
			assert(maxEnd!=null && start>0) : "Overlap lookup follows indexing of valid one-based intervals";
			int low=0, high=maxEnd.length;
			while(low<high){
				final int middle=(low+high)>>>1;
				if(maxEnd[middle]<start){low=middle+1;}else{high=middle;}
			}
			return low;
		}
		final String genome, assembly, seqid;
		final int length;
		final ArrayList<Locus> loci=new ArrayList<Locus>();
		int[] maxEnd;
	}

	/** One retained historical interval; context/hit counts never imply admission. */
	private static final class Locus{
		Locus(final String id_, final Contig contig_, final char strand_, final int start_, final int stop_){
			id=id_;
			contig=contig_;
			strand=strand_;
			start=start_;
			stop=stop_;
		}
		final String id;
		final Contig contig;
		final char strand;
		final int start, stop;
		Locus previous, next;
		int contexts, overlappingHits;
	}

	/** Actual sense-oriented context bounds and the original annotation that supplied it. */
	private static final class Query{
		Query(final String id_, final Locus locus_, final int low_, final int high_, final int left_, final int right_){
			id=id_;
			locus=locus_;
			low=low_;
			high=high_;
			left=left_;
			right=right_;
		}
		final String id;
		final Locus locus;
		final int low, high, left, right;
		int events;
	}

	/** Declared model version and parser metadata; bytes are pinned externally. */
	private static final class Model{
		Model(final String id_, final String hash_, final int clen_, final double ga_){
			id=id_;
			hash=hash_;
			clen=clen_;
			ga=ga_;
			gaText=Double.toString(ga);
			criteria=new RrnaCmTruthBuilder.Criteria(id, id, clen, 1, ga);
		}
		final String id, hash, gaText;
		final int clen;
		final double ga;
		final RrnaCmTruthBuilder.Criteria criteria;
	}

	/** One source table; each event retains its own line within this table. */
	private static final class Search{
		Search(final String id_, final String arm_, final String path_, final Model model_){
			id=id_;
			arm=arm_;
			path=path_;
			model=model_;
		}
		final String id, arm, path;
		final Model model;
	}

	/** Raw observation and projected bounds; no score averaging or winner selection. */
	private static final class Event{
		Event(final Search search_, final Query query_, final RrnaCmTruthBuilder.Hit hit_, final int start_, final int stop_, final char strand_){
			search=search_;
			query=query_;
			hit=hit_;
			start=start_;
			stop=stop_;
			strand=strand_;
		}
		final Search search;
		final Query query;
		final RrnaCmTruthBuilder.Hit hit;
		final int start, stop;
		final char strand;
		Physical physical;
	}

	/** Exact coordinate/model equivalence class with a stable sorted-run ID. */
	private static final class Physical{
		Physical(final int id_, final Event first_){
			id=id_;
			first=first_;
		}
		final int id;
		final Event first;
		int events;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in, contexts, tables, out;
	private int expectedLoci=-1, expectedContexts=-1, expectedHits=-1;
	private PrintStream outstream=System.err;
	private final ArrayList<Locus> loci=new ArrayList<Locus>();
	private final ArrayList<Contig> contigs=new ArrayList<Contig>();
	private final ArrayList<Query> queries=new ArrayList<Query>();
	private final ArrayList<Search> searches=new ArrayList<Search>();
	private final ArrayList<Event> events=new ArrayList<Event>();
	private final ArrayList<Physical> physical=new ArrayList<Physical>();
	private final ObjectMap<String,Locus> byLocus=new ObjectMap<String,Locus>(String.class, Locus.class);
	private final ObjectMap<String,Contig> byContig=new ObjectMap<String,Contig>(String.class, Contig.class);
	private final ObjectMap<String,Query> byQuery=new ObjectMap<String,Query>(String.class, Query.class);
	public static final String CONTEXT_HEADER="context_id\toriginal_id\tgenome\tseqid\tstrand\tstart1\tstop1\tcontig_length\tlflank\trflank";
	public static final String TABLE_HEADER="table_id\tmodel_id\tmodel_sha80\tclen\tga_bits\tarm\ttblout";
}
