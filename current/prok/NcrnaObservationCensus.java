package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaObservationData.Hit;
import prok.NcrnaObservationData.Locus;
import prok.NcrnaObservationData.Model;
import prok.NcrnaObservationData.Query;
import prok.NcrnaObservationData.Search;
import prok.NcrnaObservationData.Source;
import prok.NcrnaObservationData.Table;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Describes saved ncRNA observations without assigning biological labels.
 * Forked from template.A_SampleByteFile. A single metadata worker joins retained
 * rows and sweeps source-local intervals; native byte I/O handles compression.
 * One-based inclusive coordinates and individual observations remain unchanged.
 * No sequence extraction, alignment, score arbitration or admission is performed.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaObservationCensus{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the census and closes the redirected diagnostic stream on every exit. */
	public static void main(final String[] args){
		final Timer timer=new Timer();
		final NcrnaObservationCensus tool=new NcrnaObservationCensus(args);
		try{
			tool.process(timer);
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/** Native template parsing followed by path validation; output must be fresh. */
	public NcrnaObservationCensus(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		try{
			final Parser parser=parse(args);
			in=parser.in1;
			out=parser.out1;
			if(parser.append || parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1){
				throw new IllegalArgumentException("Census requires complete unpaired metadata, without append or sampling");
			}
			if(in==null || originals==null || contexts==null || queryFiles==null || out==null){
				throw new IllegalArgumentException("Required: in=projector_directory originals=loci.tsv contexts=crosswalk.tsv queryfiles=manifest.tsv out=fresh_directory");
			}
			if(!new File(in).isDirectory() || !Tools.testInputFiles(false, true, originals, contexts, queryFiles)){
				throw new IllegalArgumentException("Cannot read census input directory or metadata");
			}
			if(new File(out).exists()){throw new IllegalArgumentException("Output directory already exists: "+out);}
			ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
			ReadWrite.setZipThreads(Shared.threads());
		}catch(RuntimeException | Error failure){
			Shared.closeStream(outstream);
			throw failure;
		}
	}

	/** No biological thresholds: options locate complete inputs and output only. */
	private Parser parse(final String[] args){
		assert(args!=null) : "PreParser supplies the native argument array";
		final Parser parser=new Parser();
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("originals")){
				originals=b;
			}else if(a.equals("contexts")){
				contexts=b;
			}else if(a.equals("queryfiles")){
				queryFiles=b;
			}else if(parser.parse(arg, a, b)){
				// Standard BBTools options.
			}else{
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		return parser;
	}

	/** Validates input populations before opening outputs and summarizes observations. */
	private void process(final Timer timer){
		data=new NcrnaObservationData(in, originals, contexts, queryFiles);
		modelStrandCounts=new long[2][data.models.size()*2];
		if(!new File(out).mkdir()){throw new IllegalArgumentException("Cannot create fresh output directory: "+out);}
		copyEvidence();
		writeEvents();
		writeMembership();
		for(final Source source:data.sources){source.indexQueries();}
		try(Output relationOutput=output("relation_contexts.tsv");
				Output coverageOutput=output("relation_coverage.tsv");
				Output gapOutput=output("coverage_gaps.tsv")){
			relationContexts=relationOutput.writer;
			coverage=coverageOutput.writer;
			gaps=gapOutput.writer;
			relationContexts.println("relation_id\tcontext_id\tsearched_model_a\tsearched_model_b\tcontains_a\tcontains_b\tjoint_search_contains_both\ttouches_contig_left\ttouches_contig_right\tartificial_left_inside_span\tartificial_right_inside_span");
			coverage.println("relation_id\tmodel_id\tmodel_sha80\tspan_start1\tspan_stop1\tcovered_bases\tuncovered_bases\tuncovered_segments\tcomplete_single_contexts\tcontext_evidence");
			gaps.println("relation_id\tmodel_id\tmodel_sha80\tstart1\tstop1");
			writePairs();
			writeAnnotations();
			writeNeighbors();
			writeHistoricalCoverage();
		}
		writeSummary();
		timer.stop();
		outstream.println("NCRNA_OBSERVATION_CENSUS_PASS originals="+data.loci.size()+" contexts="+data.queries.size()
			+" searches="+data.searches.size()+" events="+data.eventTable.rows.size()+" physical_intervals="+data.hits.size());
		outstream.println("No biological admission, subtype, tandem or negative labels assigned.");
		outstream.println("Time: "+timer);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Retained evidence       ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies row bytes through native I/O; schemas and all raw text are preserved. */
	private void copyEvidence(){
		assert(data!=null) : "Only validated complete observations may be copied to census output";
		copy(data.originalTable, "originals.tsv");
		copy(data.contextTable, "context_map.tsv");
		copy(data.hitTable, "hits.tsv");
		copy(data.eventTable, "events.tsv");
		copy(data.overlapTable, "overlaps.tsv");
		copy(data.queryTable, "queries.tsv");
		copy(data.neighborTable, "historical_neighbors.tsv");
		copy(data.searchTable, "searches.tsv");
		copy(data.queryFileTable, "query_files.tsv");
	}

	private void copy(final Table table, final String name){
		assert(table!=null) : "Evidence copies retain a validated table including its header";
		try(Output output=output(name)){
			output.writer.println(table.header);
			for(final byte[] line:table.rows){output.writer.println(line);}
		}
	}

	/** Descriptive per-event geometry; raw score/E-value/truncation remain in events.tsv. */
	private void writeEvents(){
		assert(data.eventTable.rows.size()>=data.hits.size()) : "Every physical interval has at least one retained event";
		final LineParser1 parser=new LineParser1('\t');
		try(Output output=output("event_geometry.tsv")){
			output.writer.println("physical_id\ttable_id\ttable_line\tcontext_id\tmodel_span\tmodel_length\tmodel_coverage\tgenomic_start1\tgenomic_stop1\tgenomic_strand\tcontained_in_origin_context\ttouches_context_left\ttouches_context_right\ttouches_contig_left\ttouches_contig_right");
			for(final Hit hit:data.hits){
				for(int i=0; i<hit.eventRows.size; i++){
					final byte[] line=data.eventTable.rows.get(hit.eventRows.get(i));
					parser.set(line);
					final int span=parser.parseInt(7)-parser.parseInt(6)+1;
					final int localFrom=parser.parseInt(8), localTo=parser.parseInt(9);
					final int low=Math.min(localFrom, localTo), high=Math.max(localFrom, localTo);
					final Query query=data.eventQuery(hit.eventRows.get(i));
					row.clear();
					appendFields(row, line, parser, 0, 4);
					row.tab().append(span).tab().append(hit.model.clen).tab().append((double)span/hit.model.clen, 8)
						.tab().append(hit.start).tab().append(hit.stop).tab().append(hit.strand).tab().append(true)
						.tab().append(hit.start==query.start).tab().append(hit.stop==query.stop)
						.tab().append(hit.start==1).tab().append(hit.stop==hit.source.length).nl();
					assert(low>=1 && high<=query.stop-query.start+1) : "Validated local events must remain contained in their originating query";
					output.writer.print(row);
					count(0, hit.model, hit.strand);
				}
			}
		}
	}

	/** Membership rows include searches with no hits and every searched query. */
	private void writeMembership(){
		assert(!data.searches.isEmpty()) : "Zero-hit searches remain an explicit population";
		try(Output output=output("search_membership.tsv")){
			output.writer.println("table_id\tmodel_id\tmodel_sha80\tcontext_id\tgenome\tassembly\tseqid\tstart1\tstop1\tcontig_length\tlflank\trflank\torigin_strand\tevents");
			for(final Search search:data.searches){
				for(int i=0; i<search.members.size; i++){
					final Query query=data.queries.get(search.members.get(i));
					row.clear().append(bytes(search.id)).tab().append(bytes(search.model.id)).tab().append(bytes(search.model.sha))
						.tab().append(bytes(query.id)).tab();
					appendSource(row, query.origin.source);
					output.writer.print(row.tab().append(query.start).tab().append(query.stop).tab().append(query.origin.source.length)
						.tab().append(query.left).tab().append(query.right).tab().append(query.origin.strand)
						.tab().append(Math.max(0, search.eventCounts.get(query.index))).nl());
					memberships++;
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Interval relations      ----------------*/
	/*--------------------------------------------------------------*/

	/** Sorted-start sweep emits every overlapping unordered pair exactly once. */
	private void writePairs(){
		assert(data!=null) : "Physical-pair sweeps require validated source identities";
		try(Output output=output("physical_pairs.tsv")){
			output.writer.println("relation_id\tphysical_a\tphysical_b\tgenome\tassembly\tseqid\tsame_model_version\tsame_strand\toverlap_bases\tfraction_a\tfraction_b\tstart_b_minus_a\tstop_b_minus_a\tstrand_a\tmodel_a\tsha80_a\tstart_a\tstop_a\tstrand_b\tmodel_b\tsha80_b\tstart_b\tstop_b");
			for(final Source source:data.sources){
				source.hits.sort(START_ORDER);
				for(int i=0; i<source.hits.size(); i++){
					final Hit first=source.hits.get(i);
					for(int j=i+1; j<source.hits.size(); j++){
						final Hit second=source.hits.get(j);
						if(second.start>first.stop){break;}
						final Hit a=compareTuple(first, second)<=0 ? first : second;
						final Hit b=a==first ? second : first;
						final int overlap=overlap(a.start, a.stop, b.start, b.stop);
						assert(overlap>0) : "Sorted-start pair sweep only emits intersecting inclusive intervals";
						final long relation=relate(source, a.start, a.stop, a.model, b.start, b.stop, b.model);
						row.clear().append(relation).tab().append(a.id).tab().append(b.id).tab();
						appendSource(row, source);
						row.tab().append(a.model==b.model).tab().append(a.strand==b.strand)
							.tab().append(overlap).tab().append((double)overlap/(a.stop-a.start+1), 8)
							.tab().append((double)overlap/(b.stop-b.start+1), 8)
							.tab().append((long)b.start-a.start).tab().append((long)b.stop-a.stop).tab();
						appendHit(row, a);
						appendHit(row.tab(), b);
						output.writer.print(row.nl());
						physicalPairs++;
					}
				}
			}
		}
	}

	/** Preserves each original and all child pairs; nonoverlap is geometry only. */
	private void writeAnnotations(){
		assert(!data.loci.isEmpty()) : "Annotation summary includes originals with no observations";
		try(Output summary=output("annotations.tsv"); Output links=output("annotation_relations.tsv");
				Output pairs=output("annotation_child_pairs.tsv")){
			summary.writer.println(data.originalTable.header+"\tcontexts\torigin_events\tphysical_overlaps\tnonoverlapping_child_pairs");
			links.writer.println("relation_id\tgenome\toriginal_id\tphysical_id\toverlap_bases\tsame_strand");
			pairs.writer.println("relation_id\tgenome\toriginal_id\tphysical_a\tphysical_b\tsame_model_version\tsame_strand\toverlap_bases\tgap_bases");
			for(final Locus locus:data.loci){
				long nonoverlap=0, originEvents=0;
				for(int i=0; i<locus.queries.size; i++){originEvents+=data.queries.get(locus.queries.get(i)).events;}
				for(int i=0; i<locus.hits.size; i++){
					final Hit a=data.hit(locus.hits.get(i));
					final long link=relate(locus.source, locus.start, locus.stop, a.model, a.start, a.stop, a.model);
					links.writer.print(row.clear().append(link).tab().append(bytes(locus.source.genome)).tab().append(bytes(locus.id))
						.tab().append(a.id).tab().append(overlap(locus.start, locus.stop, a.start, a.stop)).tab().append(locus.strand==a.strand).nl());
					for(int j=i+1; j<locus.hits.size; j++){
						final Hit b=data.hit(locus.hits.get(j));
						final Hit first=compareTuple(a, b)<=0 ? a : b, second=first==a ? b : a;
						final int overlap=overlap(a.start, a.stop, b.start, b.stop);
						if(overlap==0){nonoverlap++;}
						final long relation=relate(locus.source, first.start, first.stop, first.model, second.start, second.stop, second.model);
						pairs.writer.print(row.clear().append(relation).tab().append(bytes(locus.source.genome)).tab().append(bytes(locus.id))
							.tab().append(first.id).tab().append(second.id).tab().append(a.model==b.model).tab().append(a.strand==b.strand)
							.tab().append(overlap).tab().append(gap(a.start, a.stop, b.start, b.stop)).nl());
						childPairs++;
					}
				}
				if(nonoverlap>0){originalsWithNonoverlap++;}
				summary.writer.print(row.clear().append(data.originalTable.rows.get(locus.index)).tab().append(locus.queries.size)
					.tab().append(originEvents).tab().append(locus.hits.size).tab().append(nonoverlap).nl());
			}
		}
	}

	/** Binary searches by start/end retain all equally nearest neighbors on each strand. */
	private void writeNeighbors(){
		assert(data!=null) : "Observed neighbors never cross a complete source identity";
		try(Output output=output("observed_neighbors.tsv")){
			output.writer.println("relation_id\tphysical_id\tside\tstrand_relation\tneighbor_id\tgap_bases\tobservation_state");
			for(final Source source:data.sources){
				for(int strandIndex=0; strandIndex<2; strandIndex++){
					final char strand=strandIndex==0 ? '+' : '-';
					int count=0;
					for(final Hit hit:source.hits){if(hit.strand==strand){count++;}}
					final Hit[] starts=new Hit[count];
					int next=0;
					for(final Hit hit:source.hits){if(hit.strand==strand){starts[next++]=hit;}}
					Arrays.sort(starts, START_ORDER);
					final Hit[] stops=starts.clone();
					Arrays.sort(stops, STOP_ORDER);
					for(final Hit hit:source.hits){
						writeNeighborSide(output.writer, hit, strand, starts, false);
						writeNeighborSide(output.writer, hit, strand, stops, true);
					}
				}
			}
		}
	}

	/** Emits all distance ties, or an explicit no-observed-neighbor row, never absence. */
	private void writeNeighborSide(final ByteStreamWriter output, final Hit hit, final char strand,
			final Hit[] candidates, final boolean left){
		assert(candidates!=null) : "Neighbor lookup uses a source-local array sorted by the compared bound";
		int low=0, high=candidates.length;
		while(low<high){
			final int middle=(low+high)>>>1;
			final boolean before=left ? candidates[middle].stop<hit.start : candidates[middle].start<=hit.stop;
			if(before){low=middle+1;}else{high=middle;}
		}
		int index=left ? low-1 : low;
		final boolean same=hit.strand==strand;
		if(index<0 || index>=candidates.length){
			output.print(row.clear().append("NA").tab().append(hit.id).tab().append(left ? "LEFT" : "RIGHT")
				.tab().append(same ? "SAME" : "OPPOSITE").append("\tNA\tNA\tUNRESOLVED_NO_OBSERVED_NEIGHBOR\n"));
			missingNeighbors++;
			return;
		}
		final int bound=left ? candidates[index].stop : candidates[index].start;
		while(index>=0 && index<candidates.length && (left ? candidates[index].stop : candidates[index].start)==bound){
			final Hit other=candidates[index];
			final long relation=relate(hit.source, hit.start, hit.stop, hit.model, other.start, other.stop, other.model);
			output.print(row.clear().append(relation).tab().append(hit.id).tab().append(left ? "LEFT" : "RIGHT")
				.tab().append(same ? "SAME" : "OPPOSITE").tab().append(other.id)
				.tab().append(gap(hit.start, hit.stop, other.start, other.stop)).append("\tOBSERVED\n"));
			neighborPairs++;
			index+=left ? -1 : 1;
		}
	}

	/** Links unchanged historical neighbor rows to model-specific context evidence. */
	private void writeHistoricalCoverage(){
		assert(data.neighborTable!=null) : "Historical relations are retained even when no observed-hit neighbor exists";
		final LineParser1 parser=new LineParser1('\t');
		try(Output output=output("historical_neighbor_coverage.tsv")){
			output.writer.println("historical_row1\trelation_id\tmodel_id\tmodel_sha80");
			for(int i=0; i<data.neighborTable.rows.size(); i++){
				final byte[] line=data.neighborTable.rows.get(i);
				parser.set(line);
				final Locus a=data.originalByFields(line, parser, 1, 2), b=data.originalByFields(line, parser, 1, 4);
				for(final Model model:data.models){
					final long relation=relate(a.source, a.start, a.stop, model, b.start, b.stop, model);
					output.writer.print(row.clear().append(i+1).tab().append(relation).tab().append(bytes(model.id)).tab().append(bytes(model.sha)).nl());
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Search coverage        ----------------*/
	/*--------------------------------------------------------------*/

	/** Assigns a reporting key and emits actual shared-context and union-gap evidence. */
	private long relate(final Source source, final int aStart, final int aStop, final Model a,
			final int bStart, final int bStop, final Model b){
		assert(source!=null && a!=null && b!=null) : "Every relationship has one source and explicit relevant model versions";
		final long id=++relations;
		final int start=Math.min(aStart, bStart), stop=Math.max(aStop, bStop);
		for(int q=source.firstQuery(start); q<source.queries.size(); q++){
			final Query query=source.queries.get(q);
			if(query.start>stop){break;}
			if(query.stop<start){continue;}
			final boolean searchedA=query.searched(a, data.searches), searchedB=query.searched(b, data.searches);
			if(!searchedA && !searchedB){continue;}
			final boolean containsA=query.start<=aStart && query.stop>=aStop;
			final boolean containsB=query.start<=bStart && query.stop>=bStop;
			relationContexts.print(contextRow.clear().append(id).tab().append(bytes(query.id)).tab().append(searchedA).tab().append(searchedB)
				.tab().append(containsA).tab().append(containsB).tab().append(searchedA && searchedB && containsA && containsB)
				.tab().append(query.start==1).tab().append(query.stop==source.length)
				.tab().append(query.start>1 && query.start>start && query.start<=stop)
				.tab().append(query.stop<source.length && query.stop>=start && query.stop<stop).nl());
		}
		writeUnion(id, source, start, stop, a);
		if(a!=b){writeUnion(id, source, start, stop, b);}
		return id;
	}

	/** Sorted actual intervals form a union, with every uncovered segment emitted. */
	private void writeUnion(final long relation, final Source source, final int start, final int stop, final Model model){
		assert(start>=1 && stop>=start && stop<=source.length) : "Coverage span is bounded by its complete source contig";
		long cursor=start, covered=0;
		int missing=0, whole=0;
		for(int q=source.firstQuery(start); q<source.queries.size(); q++){
			final Query query=source.queries.get(q);
			if(query.start>stop){break;}
			if(query.stop<start || !query.searched(model, data.searches)){continue;}
			if(query.start<=start && query.stop>=stop){whole++;}
			final int from=Math.max(start, query.start), to=Math.min(stop, query.stop);
			if(from>cursor){
				writeGap(relation, model, cursor, from-1L);
				missing++;
				cursor=from;
			}
			if(to>=cursor){
				covered+=to-cursor+1;
				cursor=(long)to+1;
			}
		}
		if(cursor<=stop){
			writeGap(relation, model, cursor, stop);
			missing++;
		}
		final long uncovered=(long)stop-start+1-covered;
		assert(covered>=0 && uncovered>=0 && (uncovered==0)==(missing==0)) : "Actual interval unions must conserve the compared span and its gaps";
		coverage.print(contextRow.clear().append(relation).tab().append(bytes(model.id)).tab().append(bytes(model.sha))
			.tab().append(start).tab().append(stop).tab().append(covered).tab().append(uncovered).tab().append(missing).tab().append(whole)
			.tab().append(uncovered>0 ? "UNRESOLVED_CONTEXT_LIMITED" : whole>0 ? "SINGLE_CONTEXT_COVERED" : "UNION_ONLY_COVERED").nl());
		if(uncovered>0){contextLimited++;}
	}

	private void writeGap(final long relation, final Model model, final long start, final long stop){
		assert(start<=stop) : "Uncovered-segment output excludes empty spans";
		gaps.print(contextRow.clear().append(relation).tab().append(bytes(model.id)).tab().append(bytes(model.sha))
			.tab().append(start).tab().append(stop).nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------       Summary and math       ----------------*/
	/*--------------------------------------------------------------*/

	/** Includes all physical intervals, including RF03166 claims without cross-model overlap. */
	private void writeSummary(){
		assert(data!=null) : "Summary denominators come from complete retained populations";
		for(final Hit hit:data.hits){count(1, hit.model, hit.strand);}
		try(Output output=output("summary.tsv")){
			output.writer.println("metric\tcount");
			summary(output.writer, "original_annotations", data.loci.size());
			summary(output.writer, "contexts", data.queries.size());
			summary(output.writer, "searches", data.searches.size());
			summary(output.writer, "events", data.eventTable.rows.size());
			summary(output.writer, "physical_intervals", data.hits.size());
			summary(output.writer, "annotation_links", data.overlapTable.rows.size());
			summary(output.writer, "historical_neighbor_rows", data.neighborTable.rows.size());
			summary(output.writer, "search_context_memberships", memberships);
			summary(output.writer, "physical_overlap_pairs", physicalPairs);
			summary(output.writer, "annotation_child_pairs", childPairs);
			summary(output.writer, "originals_with_nonoverlapping_children", originalsWithNonoverlap);
			summary(output.writer, "observed_neighbor_pairs", neighborPairs);
			summary(output.writer, "no_observed_neighbor_rows", missingNeighbors);
			summary(output.writer, "context_limited_relation_model_rows", contextLimited);
			summary(output.writer, "relationships", relations);
		}
		try(Output output=output("model_strand_counts.tsv")){
			output.writer.println("population\tmodel_id\tmodel_sha80\tstrand\tcount");
			for(int population=0; population<2; population++){
				for(final Model model:data.models){
					for(int strand=0; strand<2; strand++){
						output.writer.print(row.clear().append(population==0 ? "events" : "physical_intervals")
							.tab().append(bytes(model.id)).tab().append(bytes(model.sha)).tab().append(strand==0 ? '+' : '-')
							.tab().append(modelStrandCounts[population][model.index*2+strand]).nl());
					}
				}
			}
		}
		writeFlagCounts();
	}

	/** Preserves exact original flag combinations rather than selecting thresholds. */
	private void writeFlagCounts(){
		assert(data.originalTable.rows.size()==data.loci.size()) : "Flag distributions include every original annotation";
		final NcrnaObservationData.Symbols flags=new NcrnaObservationData.Symbols();
		final ArrayList<long[]> counts=new ArrayList<long[]>();
		final LineParser1 parser=new LineParser1('\t');
		try(Output output=output("original_flag_counts.tsv")){
			output.writer.println("partial\tpseudo\tstatus\tedge\toriginals\twithout_overlapping_hits\twith_nonoverlapping_children");
			for(final Locus locus:data.loci){
				final byte[] line=data.originalTable.rows.get(locus.index);
				parser.set(line);
				row.clear();
				appendFields(row, line, parser, 12, 15);
				parser.setBounds(16);
				row.tab().append(line, parser.a(), parser.b()-parser.a());
				final int index=flags.intern(row.array, 0, row.length());
				if(index==counts.size()){counts.add(new long[3]);}
				final long[] bucket=counts.get(index);
				bucket[0]++;
				if(locus.hits.size==0){bucket[1]++;}
				boolean nonoverlap=false;
				for(int i=0; i<locus.hits.size && !nonoverlap; i++){
					final Hit a=data.hit(locus.hits.get(i));
					for(int j=i+1; j<locus.hits.size; j++){
						final Hit b=data.hit(locus.hits.get(j));
						if(overlap(a.start, a.stop, b.start, b.stop)==0){nonoverlap=true; break;}
					}
				}
				if(nonoverlap){bucket[2]++;}
			}
			for(int i=0; i<counts.size(); i++){
				final long[] bucket=counts.get(i);
				output.writer.print(row.clear().append(flags.bytes(i)).tab().append(bucket[0]).tab().append(bucket[1]).tab().append(bucket[2]).nl());
			}
		}
	}

	/** Returns inclusive overlap bases, zero for disjoint intervals. */
	public static int overlap(final int aStart, final int aStop, final int bStart, final int bStop){
		assert(aStart>0 && bStart>0 && aStop>=aStart && bStop>=bStart) : "Census coordinates are valid one-based inclusive intervals";
		return Math.max(0, Math.min(aStop, bStop)-Math.max(aStart, bStart)+1);
	}

	/** Gap between left/right intervals; negative values denote overlapping bounds. */
	public static long gap(final int aStart, final int aStop, final int bStart, final int bStop){
		assert(aStart>0 && bStart>0 && aStop>=aStart && bStop>=bStart) : "Genomic gaps use inclusive physical bounds";
		return aStart<=bStart ? (long)bStart-aStop-1 : (long)aStart-bStop-1;
	}

	/** Full tuple order: source is common, then strand/model/version/start/stop. */
	private int compareTuple(final Hit a, final Hit b){
		assert(a.source==b.source) : "Pair orientation is only defined within one full source key";
		int c=Character.compare(a.strand, b.strand);
		if(c==0){c=compareBytes(bytes(a.model.id), bytes(b.model.id));}
		if(c==0){c=compareBytes(bytes(a.model.sha), bytes(b.model.sha));}
		if(c==0){c=Integer.compare(a.start, b.start);}
		return c==0 ? Integer.compare(a.stop, b.stop) : c;
	}

	private static int compareBytes(final byte[] a, final byte[] b){
		for(int i=0; i<Math.min(a.length, b.length); i++){
			final int c=(a[i]&255)-(b[i]&255);
			if(c!=0){return c;}
		}
		return Integer.compare(a.length, b.length);
	}

	private void count(final int population, final Model model, final char strand){
		assert(population>=0 && population<2) : "Population axis distinguishes raw events from physical intervals";
		modelStrandCounts[population][model.index*2+(strand=='+' ? 0 : 1)]++;
	}

	private void summary(final ByteStreamWriter writer, final String metric, final long value){writer.print(row.clear().append(metric).tab().append(value).nl());}
	private byte[] bytes(final int symbol){return data.symbols.bytes(symbol);}

	private static void appendFields(final ByteBuilder target, final byte[] line, final LineParser1 parser, final int from, final int to){
		assert(to>from) : "A copied field range includes at least one original field";
		parser.setBounds(from);
		final int start=parser.a();
		parser.setBounds(to-1);
		target.append(line, start, parser.b()-start);
	}

	private void appendSource(final ByteBuilder target, final Source source){
		target.append(bytes(source.genome)).tab().append(bytes(source.assembly)).tab().append(bytes(source.seqid));
	}

	private void appendHit(final ByteBuilder target, final Hit hit){
		target.append(hit.strand).tab().append(bytes(hit.model.id)).tab().append(bytes(hit.model.sha))
			.tab().append(hit.start).tab().append(hit.stop);
	}

	private Output output(final String name){
		return new Output(FileFormat.testOutput(NcrnaObservationData.path(out, name), FileFormat.TEXT, null, true, false, false, false));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in, originals, contexts, queryFiles, out;
	private PrintStream outstream=System.err;
	private NcrnaObservationData data;
	private ByteStreamWriter relationContexts, coverage, gaps;
	private final ByteBuilder row=new ByteBuilder(), contextRow=new ByteBuilder();
	private long memberships, physicalPairs, childPairs, originalsWithNonoverlap;
	private long neighborPairs, missingNeighbors, relations, contextLimited;
	private long[][] modelStrandCounts;
	private static final Comparator<Hit> START_ORDER=(a, b)->{
		int c=Integer.compare(a.start, b.start);
		if(c==0){c=Integer.compare(a.stop, b.stop);}
		return c==0 ? Integer.compare(a.id, b.id) : c;
	};
	private static final Comparator<Hit> STOP_ORDER=(a, b)->{
		int c=Integer.compare(a.stop, b.stop);
		return c==0 ? Integer.compare(a.id, b.id) : c;
	};
}
