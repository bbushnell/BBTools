package test.prok;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import prok.NcrnaCmHitProjector;
import prok.NcrnaCmLocusAudit;
import prok.NcrnaObservationCensus;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Literal census fixtures with independent geometry and membership expectations.
 * Includes equal-bound cross-model hits, nesting, an overlap chain, adjacency,
 * reverse-context projection, tied neighbors, disjoint source identities,
 * uncovered gaps, union-only coverage and complete zero-hit populations.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaObservationCensusTest{

	/** Requires a fresh output root; retains all finite inputs and inspected outputs. */
	public static void main(final String[] args){
		String out=null;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0], b=split.length>1 ? split[1] : null;
			if(a.equals("out")){out=b;}else{throw new IllegalArgumentException("Unknown test option: "+arg);}
		}
		check(out!=null && new File(out).mkdir(), "A fresh fixture root is required");
		final NcrnaObservationCensusTest test=new NcrnaObservationCensusTest(out);
		test.run(false);
		test.run(true);
		test.badMembership();
		test.badQueryIds();
		check(NcrnaObservationCensus.overlap(10, 30, 20, 40)==11, "Inclusive overlap oracle");
		check(NcrnaObservationCensus.gap(10, 30, 31, 50)==0, "Adjacent interval gap oracle");
		check(NcrnaObservationCensus.gap(120, 150, 151, 165)==0, "Reverse-context adjacency oracle");
		System.err.println("NCRNA_OBSERVATION_CENSUS_FINITE_PASS");
	}

	private NcrnaObservationCensusTest(final String root_){root=root_;}

	/** End-to-end input/output fixture; expected counts do not call the tested geometry. */
	private void run(final boolean zero){
		final String dir=path(root, zero ? "zero" : "geometry");
		check(new File(dir).mkdir(), "Fresh case directory");
		writeFixture(dir, zero);
		invoke(dir, "queryfiles.tsv", "result");
		final String result=path(dir, "result");
		check(metric(result, "original_annotations")==4, "All four originals must survive");
		check(metric(result, "contexts")==4, "All four contexts must survive");
		check(metric(result, "searches")==3, "All three searches must survive");
		check(metric(result, "search_context_memberships")==8, "Repeated FASTA membership is per search, not per file");
		check(metric(result, "events")== (zero ? 0 : 11), "Raw observations must not collapse into physical intervals");
		check(metric(result, "physical_intervals")== (zero ? 0 : 10), "Distinct source/model/strand/bounds identities");
		check(metric(result, "physical_overlap_pairs")== (zero ? 0 : 7), "Seven literal overlapping physical pairs");
		check(metric(result, "annotation_links")== (zero ? 0 : 10), "Every physical interval overlaps its literal original");
		check(metric(result, "historical_neighbor_rows")==2, "Historical neighbors remain separate from physical neighbors");
		check(metric(result, "annotation_child_pairs")== (zero ? 0 : 13), "Five choose two plus three choose two");
		check(metric(result, "originals_with_nonoverlapping_children")== (zero ? 0 : 2), "Two originals have disjoint child bounds");
		for(final String name:new String[]{"hits.tsv", "events.tsv", "overlaps.tsv", "queries.tsv", "searches.tsv"}){
			identical(path(dir, name), path(result, name));
		}
		identical(path(dir, "neighbors.tsv"), path(result, "historical_neighbors.tsv"));
		identical(path(dir, "originals.tsv"), path(result, "originals.tsv"));
		if(!zero){inspectGeometry(result);}
		final ArrayList<byte[]> membership=rows(path(result, "search_membership.tsv"));
		int zeroMembership=0;
		for(final byte[] line:membership){
			parser.set(line);
			if(parser.termEquals("q4", 3)){
				zeroMembership++;
				if(zero){check(parser.parseInt(13)==0, "Zero-hit search membership keeps explicit zero events");}
			}
		}
		check(zeroMembership==2, "q4 remains explicitly searched in both A tables");
	}

	/** Checks literal pairs, model-specific union gaps, ties and reverse-context boundaries. */
	private void inspectGeometry(final String result){
		int ties=0;
		long relationTo6=-1;
		for(final byte[] line:rows(path(result, "observed_neighbors.tsv"))){
			parser.set(line);
			if(parser.parseInt(1)==5 && parser.termEquals("RIGHT", 2) && parser.termEquals("SAME", 3)){
				final int neighbor=parser.parseInt(4);
				check(neighbor==6 || neighbor==7, "Both equal-start right-neighbor ties survive");
				check(parser.parseLong(5)==49, "Neighbor gap is in genomic bases");
				if(neighbor==6){relationTo6=parser.parseLong(0);}
				ties++;
			}
		}
		check(ties==2 && relationTo6>0, "Exactly two nearest same-strand right neighbors for interval5");
		int foundGap=0;
		for(final byte[] line:rows(path(result, "coverage_gaps.tsv"))){
			parser.set(line);
			if(parser.parseLong(0)==relationTo6){
				check(parser.termEquals("A", 1) && parser.parseInt(3)==101 && parser.parseInt(4)==110, "Actual uncovered gap must not become covered bounding span");
				foundGap++;
			}
		}
		check(foundGap==1, "One ten-base uncovered segment between A's contexts");
		boolean unionOnly=false;
		for(final byte[] line:rows(path(result, "relation_coverage.tsv"))){
			parser.set(line);
			if(parser.termEquals("B", 1) && parser.parseInt(3)==10 && parser.parseInt(4)==110){
				check(parser.parseInt(5)==101 && parser.parseInt(6)==0 && parser.parseInt(8)==0,
					"B union covers historical pair without one complete shared context");
				check(parser.termEquals("UNION_ONLY_COVERED", 9), "Union coverage must not become a joint-search claim");
				unionOnly=true;
			}
		}
		check(unionOnly, "Historical neighbor coverage exposes the union-only fixture");
		int reverse=0;
		for(final byte[] line:rows(path(result, "event_geometry.tsv"))){
			parser.set(line);
			if(parser.parseInt(0)==6){
				check(parser.parseInt(7)==120 && parser.parseInt(8)==150 && parser.termEquals("+", 9), "Reverse-query/negative-local hit becomes positive genomic120-150");
				reverse++;
			}
		}
		check(reverse==1, "Reverse-context event preserved once");
		int equalBounds=0;
		for(final byte[] line:rows(path(result, "physical_pairs.tsv"))){
			parser.set(line);
			check(parser.parseInt(1)!=10 && parser.parseInt(2)!=10, "Equal contig names from separate genome/assembly sources must never join");
			if(parser.parseInt(1)==1 && parser.parseInt(2)==3){
				check(parser.termEquals("false", 6) && parser.parseInt(8)==21
					&& parser.parseDouble(9)==1 && parser.parseDouble(10)==1, "Equal bounds remain a cross-model pair with two own-span denominators");
				equalBounds++;
			}
		}
		check(equalBounds==1, "One unique cross-model equal-bound pair");
	}

	/** Missing and duplicate bindings fail before output creation. */
	private void badMembership(){
		final String dir=path(root, "geometry");
		write(path(dir, "missing.tsv"), "table_id\tquery_fasta\na1\t"+path(dir, "a.fa")+"\n");
		expectFailure(dir, "missing.tsv", "missing_result");
		write(path(dir, "duplicate.tsv"), "table_id\tquery_fasta\na1\t"+path(dir, "a.fa")+"\na1\t"+path(dir, "a.fa")+"\n");
		expectFailure(dir, "duplicate.tsv", "duplicate_result");
	}

	/** Public t=16 must still reject malformed query membership without a stranded producer. */
	private void badQueryIds(){
		for(final String id:new String[]{"q1", "unknown"}){
			final String dir=path(root, "bad_id_"+id);
			check(new File(dir).mkdir(), "Fresh malformed-query fixture");
			writeFixture(dir, true);
			final FileFormat append=FileFormat.testOutput(path(dir, "a.fa"), FileFormat.FASTA, null, true, false, true, false);
			try(Output output=new Output(append)){
				final ByteBuilder row=new ByteBuilder().append('>').append(id).nl();
				for(int i=0; i<100; i++){row.append('A');}
				output.writer.print(row.nl());
			}
			boolean failed=false;
			try{
				NcrnaObservationCensus.main(new String[]{"in="+dir, "originals="+path(dir, "originals.tsv"),
					"contexts="+path(dir, "contexts.tsv"), "queryfiles="+path(dir, "queryfiles.tsv"),
					"out="+path(dir, "result"), "t=16"});
			}catch(IllegalArgumentException expected){failed=true;}
			check(failed && !new File(path(dir, "result")).exists(), "Duplicate/unknown FASTA IDs fail before output and return to the caller");
		}
	}

	private static void expectFailure(final String dir, final String manifest, final String result){
		boolean failed=false;
		try{
			invoke(dir, manifest, result);
		}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !new File(path(dir, result)).exists(), "Invalid membership must fail before output creation");
	}

	private static void invoke(final String dir, final String queryFiles, final String result){
		NcrnaObservationCensus.main(new String[]{"in="+dir, "originals="+path(dir, "originals.tsv"),
			"contexts="+path(dir, "contexts.tsv"), "queryfiles="+path(dir, queryFiles), "out="+path(dir, result), "t=1"});
	}

	/** Writes literal normalized metadata; no inference from the implementation under test. */
	private void writeFixture(final String dir, final boolean zero){
		final StringBuilder originals=new StringBuilder(NcrnaCmLocusAudit.HEADER+"\n");
		original(originals, "L1", "gA", "asmA", '+', 10, 80, 200, "false");
		original(originals, "L2", "gA", "asmA", '-', 120, 160, 200, "true");
		original(originals, "L3", "gA", "asmA", '+', 90, 110, 200, "false");
		original(originals, "Z", "gB", "asmB", '+', 20, 30, 50, "false");
		write(path(dir, "originals.tsv"), originals.toString());
		write(path(dir, "contexts.tsv"), NcrnaCmHitProjector.CONTEXT_HEADER+"\n"
			+"q1\tL1\tgA\tsame\t+\t10\t80\t200\t9\t20\n"
			+"q2\tL2\tgA\tsame\t-\t120\t160\t200\t20\t9\n"
			+"q3\tL3\tgA\tsame\t+\t90\t110\t200\t9\t10\n"
			+"q4\tZ\tgB\tsame\t+\t20\t30\t50\t19\t20\n");
		fasta(path(dir, "a.fa"), new String[]{"q1", "q2", "q4"}, new int[]{100, 70, 50});
		fasta(path(dir, "b.fa"), new String[]{"q1", "q3"}, new int[]{100, 40});
		write(path(dir, "a1.tbl"), "# Target file: "+path(dir, "a.fa")+"\n");
		write(path(dir, "a2.tbl"), "# Target file: "+path(dir, "a.fa")+"\n");
		write(path(dir, "b.tbl"), "# Target file: "+path(dir, "b.fa")+"\n");
		write(path(dir, "searches.tsv"), NcrnaCmHitProjector.TABLE_HEADER+"\n"
			+"a1\tA\t"+SHA_A+"\t100\t1\tfinite\t"+path(dir, "a1.tbl")+"\n"
			+"a2\tA\t"+SHA_A+"\t100\t1\tfinite\t"+path(dir, "a2.tbl")+"\n"
			+"b\tB\t"+SHA_B+"\t100\t2\tfinite\t"+path(dir, "b.tbl")+"\n");
		write(path(dir, "queryfiles.tsv"), "table_id\tquery_fasta\na1\t"+path(dir, "a.fa")
			+"\na2\t"+path(dir, "a.fa")+"\nb\t"+path(dir, "b.fa")+"\n");
		final StringBuilder hits=new StringBuilder(HIT_HEADER+"\n");
		final StringBuilder events=new StringBuilder(EVENT_HEADER+"\n");
		final StringBuilder overlaps=new StringBuilder(OVERLAP_HEADER+"\n");
		if(!zero){
			hit(hits, 1, '+', "A", 10, 30, 2);
			hit(hits, 2, '+', "A", 31, 50, 1);
			hit(hits, 3, '+', "B", 10, 30, 1);
			hit(hits, 4, '-', "A", 20, 40, 1);
			hit(hits, 5, '+', "A", 45, 70, 1);
			hit(hits, 6, '+', "A", 120, 150, 1);
			hit(hits, 7, '+', "A", 120, 155, 1);
			hit(hits, 8, '-', "B", 90, 110, 1);
			hit(hits, 9, '-', "A", 151, 165, 1);
			hits.append("10\tgB\tasmB\tsame\t+\tA\t").append(SHA_A).append("\t20\t30\t1\t1\n");
			event(events, 1, "a1", 1, "q1", "L1", 10, 30, '+');
			event(events, 1, "a2", 1, "q1", "L1", 10, 30, '+');
			event(events, 2, "a1", 2, "q1", "L1", 31, 50, '+');
			event(events, 3, "b", 1, "q1", "L1", 10, 30, '+');
			event(events, 4, "a1", 3, "q1", "L1", 40, 20, '-');
			event(events, 5, "a1", 4, "q1", "L1", 45, 70, '+');
			event(events, 6, "a1", 5, "q2", "L2", 61, 31, '-');
			event(events, 7, "a1", 6, "q2", "L2", 61, 26, '-');
			event(events, 8, "b", 2, "q3", "L3", 30, 10, '-');
			event(events, 9, "a1", 7, "q2", "L2", 16, 30, '+');
			events.append("10\ta1\t8\tq4\tgB\tZ\t1\t11\t20\t30\t+\tno\t!\t40.0\t1e-8\n");
			overlaps.append("1\tgA\tL1\t21\ttrue\ttrue\n2\tgA\tL1\t20\ttrue\ttrue\n")
				.append("3\tgA\tL1\t21\ttrue\ttrue\n4\tgA\tL1\t21\tfalse\ttrue\n5\tgA\tL1\t26\ttrue\ttrue\n")
				.append("6\tgA\tL2\t31\tfalse\ttrue\n7\tgA\tL2\t36\tfalse\ttrue\n")
				.append("8\tgA\tL3\t21\tfalse\ttrue\n9\tgA\tL2\t10\ttrue\ttrue\n10\tgB\tZ\t11\ttrue\ttrue\n");
		}
		write(path(dir, "hits.tsv"), hits.toString());
		write(path(dir, "events.tsv"), events.toString());
		write(path(dir, "overlaps.tsv"), overlaps.toString());
		write(path(dir, "queries.tsv"), QUERY_HEADER+"\n"
			+"q1\tgA\tL1\t1\t100\t9\t20\ttrue\tfalse\t"+(zero ? 0 : 6)+"\t"+(zero ? 0 : 5)+"\n"
			+"q2\tgA\tL2\t111\t180\t20\t9\tfalse\tfalse\t"+(zero ? 0 : 3)+"\t"+(zero ? 0 : 3)+"\n"
			+"q3\tgA\tL3\t81\t120\t9\t10\tfalse\tfalse\t"+(zero ? 0 : 1)+"\t"+(zero ? 0 : 1)+"\n"
			+"q4\tgB\tZ\t1\t50\t19\t20\ttrue\ttrue\t"+(zero ? 0 : 1)+"\t"+(zero ? 0 : 1)+"\n");
		write(path(dir, "neighbors.tsv"), NEIGHBOR_HEADER+"\n"
			+"q1\tgA\tL1\tSUCCESSOR\tL3\t9\t0\ttrue\ttrue\n"
			+"q3\tgA\tL3\tPREDECESSOR\tL1\t9\t0\ttrue\ttrue\n");
	}

	private static void original(final StringBuilder text, final String id, final String genome, final String assembly,
			final char strand, final int start, final int stop, final int length, final String partial){
		text.append(id).append('\t').append(genome).append('\t').append(assembly)
			.append("\t1\tbacteria\tfinite\tsame\t").append(strand).append('\t').append(start).append('-').append(stop)
			.append('\t').append(stop-start+1).append("\t1\tx\t").append(partial).append("\tfalse\tfinite\t")
			.append(length).append("\tfalse\t0\t0\t0\t1\n");
	}

	private static void hit(final StringBuilder text, final int id, final char strand, final String model,
			final int start, final int stop, final int events){
		text.append(id).append("\tgA\tasmA\tsame\t").append(strand).append('\t').append(model).append('\t')
			.append(model.equals("A") ? SHA_A : SHA_B).append('\t').append(start).append('\t').append(stop)
			.append('\t').append(events).append("\t1\n");
	}

	private static void event(final StringBuilder text, final int id, final String table, final int line,
			final String query, final String origin, final int from, final int to, final char strand){
		text.append(id).append('\t').append(table).append('\t').append(line).append('\t').append(query)
			.append("\tgA\t").append(origin).append("\t1\t").append(Math.abs(to-from)+1).append('\t').append(from)
			.append('\t').append(to).append('\t').append(strand).append("\tno\t!\t40.0\t1e-8\n");
	}

	private static void fasta(final String path, final String[] ids, final int[] lengths){
		try(Output output=new Output(path)){
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<ids.length; i++){
				row.clear().append('>').append(ids[i]).nl();
				for(int j=0; j<lengths[i]; j++){row.append('A');}
				output.writer.print(row.nl());
			}
		}
	}

	private static long metric(final String result, final String metric){
		final LineParser1 parser=new LineParser1('\t');
		for(final byte[] line:rows(path(result, "summary.tsv"))){
			parser.set(line);
			if(parser.termEquals(metric, 0)){return parser.parseLong(1);}
		}
		throw new AssertionError("Missing summary metric: "+metric);
	}

	private static ArrayList<byte[]> rows(final String path){
		final ArrayList<byte[]> rows=new ArrayList<byte[]>();
		try(Input input=new Input(path)){
			check(input.reader.nextLine()!=null, "Table header must be present");
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){rows.add(line);}
		}
		return rows;
	}

	private static void identical(final String first, final String second){
		try(Input a=new Input(first); Input b=new Input(second)){
			for(byte[] line=a.reader.nextLine(); line!=null; line=a.reader.nextLine()){
				check(java.util.Arrays.equals(line, b.reader.nextLine()), "Evidence rows must be byte-identical");
			}
			check(b.reader.nextLine()==null, "Evidence copies must have no extra rows");
		}
	}

	private static void write(final String path, final String text){
		try(Output output=new Output(path)){output.writer.print(text);}
	}

	private static String path(final String root, final String name){return new File(root, name).getPath();}
	private static void check(final boolean value, final String why){if(!value){throw new AssertionError(why);}}
	private final String root;
	private final LineParser1 parser=new LineParser1('\t');
	private static final String SHA_A="11111111111111111111", SHA_B="22222222222222222222";
	private static final String HIT_HEADER="physical_id\tgenome\tassembly\tseqid\tstrand\tmodel_id\tmodel_sha80\tstart1\tstop1\tevents\toverlapping_annotations";
	private static final String EVENT_HEADER="physical_id\ttable_id\ttable_line\tcontext_id\torigin_genome\torigin_id\tmodel_from\tmodel_to\tlocal_from\tlocal_to\tlocal_strand\ttruncation\tinclusion\tbit_score\tevalue";
	private static final String OVERLAP_HEADER="physical_id\tgenome\toriginal_id\toverlap_bases\tsame_strand\tannotation_queried";
	private static final String QUERY_HEADER="context_id\tgenome\toriginal_id\tcontext_start1\tcontext_stop1\tlflank\trflank\ttouches_contig_left\ttouches_contig_right\torigin_events\tannotation_overlapping_hit_intervals";
	private static final String NEIGHBOR_HEADER="context_id\tgenome\toriginal_id\trelationship\tneighbor_id\tsigned_gap\toverlap_bases\tsame_strand\tneighbor_queried";
}
