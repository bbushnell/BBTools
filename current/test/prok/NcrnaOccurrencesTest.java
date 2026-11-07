package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import prok.NcrnaCmHitProjector;
import prok.NcrnaCmLocusAudit;
import prok.NcrnaOccurrences;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Finite public-CLI tests with literal containment and occurrence expectations.
 * Inputs and outputs remain beneath a fresh caller-supplied directory. No search,
 * sequence extraction, model construction or production resource is involved.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaOccurrencesTest{

	/** Runs finite cases; the developer launcher supplies assertions and JVM flags. */
	public static void main(final String[] args){
		String out=null;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("out")){out=b;}
			else{throw new IllegalArgumentException("Unknown fixture parameter: "+arg);}
		}
		check(out!=null && new File(out).mkdir(), "A fresh fixture root is required");
		new NcrnaOccurrencesTest(out).run();
		System.err.println("NCRNA_OCCURRENCES_FINITE_PASS");
	}

	private NcrnaOccurrencesTest(final String root_){root=root_;}

	/** Expected representatives are literal oracles, never computed by containment. */
	private void run(){
		pass("genomic_only", hits(hit(100, 220, 10, 80), hit(120, 180, 1, 94)), 1, 2);
		pass("model_only", hits(hit(100, 180, 1, 94), hit(160, 240, 9, 86)), 1, 2);
		pass("one_unwitnessed", hits(hit(100, 220, 1, 60), hit(120, 200, 10, 50, 55, 90)), 1, 2);
		final String two=pass("two_witnesses", hits(hit(100, 220, 1, 60, 50, 94),
			hit(120, 200, 10, 50, 55, 90)), 1, 1);
		expectLine(two, "event_links.tsv", "p\t3\t2\tRF00504_occurrence_1\tCONTAINED\tp\t1\tq\tgA\tL1");
		expectLine(two, "event_links.tsv", "p\t4\t2\tRF00504_occurrence_1\tCONTAINED\tp\t2\tq\tgA\tL1");
		pass("union_trap", hits(hit(100, 220, 1, 55, 40, 94), hit(120, 200, 1, 94)), 1, 2);
		pass("nested", hits(hit(100, 250, 1, 94), hit(120, 230, 9, 86), hit(140, 180, 20, 40)), 1, 1, 1);
		fail("ambiguous", hits(hit(100, 200, 1, 94), hit(150, 250, 1, 94), hit(160, 180, 20, 40)),
			"none", "Incomparable maximal containers");
		pass("adjacent_and_gap14", hits(hit(1, 40, 7, 88), hit(41, 140, 1, 94), hit(155, 254, 1, 94)), 1, 2, 3);
		pass("chain_minus_1", hits(minus(46458, 46557, 1, 94), minus(46466, 46549, 9, 86),
			minus(46555, 46651, 1, 94)), 1, 1, 3);
		pass("chain_minus_2", hits(minus(8610, 8709, 1, 94), minus(8618, 8701, 9, 86, 9, 86),
			minus(8707, 8803, 1, 94, 1, 94)), 1, 1, 3);
		pass("chain_plus", hits(hit(2217, 2313, 1, 94, 1, 94), hit(2311, 2410, 1, 94),
			hit(2319, 2402, 9, 86, 9, 86)), 1, 2, 2);
		final String bridge=pass("secondary_bridge", hits(hit(185146, 185253, 1, 94), hit(185255, 185338, 7, 88),
			new HitSpec("gA", '+', true, 185142, 185338, 1, 94)), 1, 2, 0);
		expectRows(bridge, "physical_links.tsv", 4);
		expectLine(bridge, "physical_links.tsv", "3\tRF00504_occurrence_1\tSECONDARY_OVERLAP");
		expectLine(bridge, "physical_links.tsv", "3\tRF00504_occurrence_2\tSECONDARY_OVERLAP");
		expectLine(bridge, "original_links.tsv", "gA\tL1\t3\t.\t197\ttrue\tSECONDARY_ANNOTATION_EVIDENCE");
		expectLine(bridge, "original_status.tsv", "gA\tL1\t3\t2");
		final String secondary=pass("unmapped_secondary", hits(new HitSpec("gA", '+', true, 100, 200, 1, 94)), 0);
		expectLine(secondary, "physical_links.tsv", "1\t.\tUNMAPPED_SECONDARY");
		expectLine(secondary, "summary.tsv", "unmapped_secondary_hits\t1");
		pass("sources_and_strands", hits(hit(100, 200, 1, 94), minus(100, 200, 1, 94),
			new HitSpec("gB", '+', false, 100, 200, 1, 94)), 1, 2, 3);
		final String zero=pass("zero", hits());
		expectLine(zero, "summary.tsv", "raw_events\t0");
		expectRows(zero, "event_links.tsv", 0);
		fail("duplicate_tuple", hits(hit(100, 200, 1, 94), hit(100, 200, 9, 86)), "none", "Repeated physical");
		fail("bad_id", hits(hit(100, 200, 1, 94)), "id", "dense");
		fail("bad_span", hits(hit(100, 200, 0, 95)), "none", "model span");
		fail("bad_count", hits(hit(100, 200, 1, 94)), "count", "event count mismatch");
		fail("bad_overlap_count", hits(hit(100, 200, 1, 94)), "overlap_count", "annotation count mismatch");
		fail("duplicate_event", hits(hit(100, 200, 1, 94, 1, 94)), "duplicate_event", "Repeated original search");
		fail("unknown_origin", hits(hit(100, 200, 1, 94)), "origin", "Unknown original");
		fail("wrong_model", hits(hit(100, 200, 1, 94)), "model", "Event origin/model differs");
	}

	/** Runs each valid case in both row orders and compares all six output tables. */
	private String pass(final String name, final HitSpec[] hits, final int... representatives){
		check(hits.length==representatives.length, "Each physical interval has a literal expected representative");
		final String first=path(root, name), second=path(root, name+"_reversed");
		fixture(first, hits, false, "none");
		fixture(second, hits, true, "none");
		invoke(first);
		invoke(second);
		final String result=path(first, "result");
		for(final String table:OUTPUTS){identical(path(result, table), path(path(second, "result"), table));}
		final ArrayList<byte[]> occurrences=rows(path(result, "occurrences.tsv"));
		int expectedCount=0;
		for(int i=0; i<representatives.length; i++){
			if(representatives[i]==i+1){expectedCount++;}
		}
		check(occurrences.size()==expectedCount, "Literal occurrence count for "+name);
		final int[] occurrenceByPhysical=new int[hits.length+1];
		int nextOccurrence=0;
		for(final byte[] line:occurrences){
			parser.set(line);
			final int id=parser.parseInt(10);
			check(id>0 && id<=hits.length && representatives[id-1]==id, "Only a literal expected representative may supply coordinates");
			final HitSpec expected=hits[id-1];
			check(parser.parseInt(6)==expected.start && parser.parseInt(7)==expected.stop,
				"Original representative bounds must survive without union or clipping");
			check(parser.termEquals(expected.genome, 1) && parser.termEquals(expected.strand, 5), "Full source and strand survive");
			check(parser.termEquals("RF00504", 8), "Secondary model must never create an occurrence");
			check(occurrenceByPhysical[id]==0, "Each representative occurs exactly once");
			occurrenceByPhysical[id]=occurrenceId(line);
			check(occurrenceByPhysical[id]==++nextOccurrence, "Occurrence IDs must be dense and in output order");
		}
		for(int i=0; i<hits.length; i++){
			if(representatives[i]==0){continue;}
			final int representative=representatives[i];
			check(occurrenceByPhysical[representative]>0, "Expected representative must appear");
			expectLine(result, "physical_links.tsv", (i+1)+"\tRF00504_occurrence_"+occurrenceByPhysical[representative]
				+"\t"+(representative==i+1 ? "REPRESENTATIVE" : "CONTAINED"));
		}
		expectRows(result, "original_status.tsv", 4);
		expectLine(result, "original_status.tsv", "gA\tZ\t0\t0");
		expectLine(result, "summary.tsv", "physical_hits\t"+hits.length);
		expectLine(result, "summary.tsv", "occurrences\t"+expectedCount);
		checkEvents(result, hits, representatives, occurrenceByPhysical);
		return result;
	}

	/** Every original event key survives; primary events have one resolved link. */
	private void checkEvents(final String result, final HitSpec[] hits, final int[] representatives, final int[] ids){
		int total=0;
		for(final HitSpec hit:hits){total+=hit.spans.length/2;}
		final int[] eventHit=new int[total+1], seen=new int[total+1];
		int event=0;
		for(int i=0; i<hits.length; i++){
			for(int j=0; j<hits[i].spans.length; j+=2){eventHit[++event]=i;}
		}
		for(final byte[] line:rows(path(result, "event_links.tsv"))){
			parser.set(line);
			final int key=parser.parseInt(1);
			check(key>0 && key<=total, "Output event keys must come from the literal input");
			final int index=eventHit[key];
			check(parser.parseInt(2)==index+1 && parser.termEquals(hits[index].secondary ? "s" : "p", 0),
				"Event identity and original physical interval must survive");
			seen[key]++;
			if(!hits[index].secondary){
				check(parser.termEquals("RF00504_occurrence_"+ids[representatives[index]], 3), "Primary event must link to its literal expected representative");
			}
		}
		for(int key=1; key<=total; key++){
			check(seen[key]>0 && (hits[eventHit[key]].secondary || seen[key]==1), "Every event survives; primary events occur once");
		}
		expectLine(result, "summary.tsv", "raw_events\t"+total);
	}

	/** Invalid cases must throw the intended error and create no result directory. */
	private void fail(final String name, final HitSpec[] hits, final String fault, final String message){
		for(final boolean reverse:new boolean[]{false, true}){
			final String dir=path(root, name+(reverse ? "_reversed" : ""));
			fixture(dir, hits, reverse, fault);
			boolean failed=false;
			try{invoke(dir);}
			catch(IllegalArgumentException expected){
				check(expected.getMessage().contains(message), "Wrong failure for "+name+": "+expected);
				failed=true;
			}
			check(failed && !new File(path(dir, "result")).exists(), "Invalid case must fail before output: "+name);
		}
	}

	/** Fixture serialization alone; it does not apply the tested containment relation. */
	private void fixture(final String dir, final HitSpec[] hits, final boolean reverse, final String fault){
		check(new File(dir).mkdir(), "Fresh case directory required");
		final ArrayList<byte[]> originals=new ArrayList<byte[]>(), search=new ArrayList<byte[]>();
		final ArrayList<byte[]> physical=new ArrayList<byte[]>(), events=new ArrayList<byte[]>(), overlaps=new ArrayList<byte[]>();
		original(originals, "L1", "gA", '+', 1, 199000);
		original(originals, "Lm", "gA", '-', 1, 199000);
		original(originals, "LB", "gB", '+', 1, 199000);
		original(originals, "Z", "gA", '+', 199990, 200000);
		search.add(new ByteBuilder().append("p\tRF00504\t").append(SHA_PRIMARY).append("\t94\t1\tfixture\tunread-primary.tbl").toBytes());
		search.add(new ByteBuilder().append("s\tRF03166\t").append(SHA_SECONDARY).append("\t94\t1\tfixture\tunread-secondary.tbl").toBytes());
		int line=0;
		final ByteBuilder row=new ByteBuilder();
		for(int i=0; i<hits.length; i++){
			final HitSpec hit=hits[i];
			final int id=i+1, count=hit.spans.length/2;
			final String origin=hit.genome.equals("gB") ? "LB" : (hit.strand=='+' ? "L1" : "Lm");
			row.clear().append(fault.equals("id") ? id+1 : id).tab().append(hit.genome).tab()
				.append(hit.genome.equals("gB") ? "asmB" : "asmA").append("\tsame\t").append(hit.strand).tab()
				.append(hit.secondary ? "RF03166" : "RF00504").tab().append(hit.secondary ? SHA_SECONDARY : SHA_PRIMARY)
				.tab().append(hit.start).tab().append(hit.stop).tab().append(count+(fault.equals("count") ? 1 : 0))
				.tab().append(fault.equals("overlap_count") ? 0 : 1);
			physical.add(row.toBytes());
			for(int j=0; j<hit.spans.length; j+=2){
				line++;
				row.clear().append(id).tab().append(hit.secondary || fault.equals("model") ? "s" : "p").tab()
					.append(fault.equals("duplicate_event") ? 1 : line).append("\tq\t").append(hit.genome).tab()
					.append(fault.equals("origin") ? "unknown" : origin).tab().append(hit.spans[j]).tab().append(hit.spans[j+1])
					.tab().append(hit.strand=='+' ? hit.start : hit.stop).tab().append(hit.strand=='+' ? hit.stop : hit.start)
					.tab().append(hit.strand).append("\tno\t!\t40\t1e-8");
				events.add(row.toBytes());
			}
			row.clear().append(id).tab().append(hit.genome).tab().append(origin).tab()
				.append(hit.stop-hit.start+1).append("\ttrue\ttrue");
			overlaps.add(row.toBytes());
		}
		write(dir, "originals.tsv", NcrnaCmLocusAudit.HEADER, originals, reverse);
		write(dir, "searches.tsv", NcrnaCmHitProjector.TABLE_HEADER, search, reverse);
		write(dir, "hits.tsv", HIT_HEADER, physical, reverse);
		write(dir, "events.tsv", EVENT_HEADER, events, reverse);
		write(dir, "overlaps.tsv", OVERLAP_HEADER, overlaps, reverse);
	}

	private static void original(final ArrayList<byte[]> rows, final String id, final String genome,
			final char strand, final int start, final int stop){
		rows.add(new ByteBuilder().append(id).tab().append(genome).tab().append(genome.equals("gB") ? "asmB" : "asmA")
			.append("\t1\tbacteria\tfixture\tsame\t").append(strand).tab().append(start).append('-').append(stop)
			.tab().append(stop-start+1).append("\t1\tx\tfalse\tfalse\tfixture\t200000\tfalse\t0\t0\t0\t1").toBytes());
	}

	private static void invoke(final String dir){
		NcrnaOccurrences.main(new String[]{"in="+dir, "primary=RF00504", "out="+path(dir, "result"), "t=1"});
	}

	private static void write(final String dir, final String file, final String header, final ArrayList<byte[]> rows, final boolean reverse){
		try(Output output=new Output(path(dir, file))){
			output.writer.println(header);
			for(int i=0; i<rows.size(); i++){output.writer.println(rows.get(reverse ? rows.size()-1-i : i));}
		}
	}

	private static ArrayList<byte[]> rows(final String path){
		final ArrayList<byte[]> rows=new ArrayList<byte[]>();
		try(Input input=new Input(path)){
			check(input.reader.nextLine()!=null, "Table header required: "+path);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){rows.add(line);}
		}
		return rows;
	}

	private static void expectRows(final String dir, final String file, final int count){
		check(rows(path(dir, file)).size()==count, "Unexpected row count for "+file+"; expected "+count);
	}

	private static void expectLine(final String dir, final String file, final String expected){
		final byte[] bytes=expected.getBytes(java.nio.charset.StandardCharsets.UTF_8);
		int found=0;
		for(final byte[] row:rows(path(dir, file))){if(Arrays.equals(row, bytes)){found++;}}
		check(found==1, "Expected exactly one row in "+file+": "+expected+"; found "+found);
	}

	private static void identical(final String first, final String second){
		try(Input a=new Input(first); Input b=new Input(second)){
			for(byte[] line=a.reader.nextLine(); line!=null; line=a.reader.nextLine()){
				check(Arrays.equals(line, b.reader.nextLine()), "Row order must not change output: "+first);
			}
			check(b.reader.nextLine()==null, "Reordered input produced extra rows: "+second);
		}
	}

	/** Parses the documented output ID only; it does not derive biological representatives. */
	private int occurrenceId(final byte[] row){
		parser.set(row);
		final String text=parser.parseString(0);
		check(text.startsWith("RF00504_occurrence_"), "Stable occurrence prefix required");
		return Integer.parseInt(text.substring("RF00504_occurrence_".length()));
	}

	private static HitSpec hit(final int start, final int stop, final int... spans){return new HitSpec("gA", '+', false, start, stop, spans);}
	private static HitSpec minus(final int start, final int stop, final int... spans){return new HitSpec("gA", '-', false, start, stop, spans);}
	private static HitSpec[] hits(final HitSpec... hits){return hits;}
	private static String path(final String root, final String name){return new File(root, name).getPath();}
	private static void check(final boolean value, final String why){if(!value){throw new AssertionError(why);}}

	/** Retained literal fixture description, not production evidence or a model result. */
	private static final class HitSpec{
		HitSpec(final String genome_, final char strand_, final boolean secondary_, final int start_, final int stop_, final int... spans_){
			check(spans_.length>0 && spans_.length%2==0, "Each fixture event needs two literal model positions");
			genome=genome_;
			strand=strand_;
			secondary=secondary_;
			start=start_;
			stop=stop_;
			spans=spans_;
		}
		final String genome;
		final char strand;
		final boolean secondary;
		final int start, stop;
		final int[] spans;
	}

	private final String root;
	private final LineParser1 parser=new LineParser1('\t');
	/** Synthetic model-version identities for finite metadata, not hashes of real files. */
	private static final String SHA_PRIMARY="11111111111111111111", SHA_SECONDARY="22222222222222222222";
	private static final String[] OUTPUTS={"occurrences.tsv", "physical_links.tsv", "event_links.tsv",
		"original_links.tsv", "original_status.tsv", "summary.tsv"};
	private static final String HIT_HEADER="physical_id\tgenome\tassembly\tseqid\tstrand\tmodel_id\tmodel_sha80\tstart1\tstop1\tevents\toverlapping_annotations";
	private static final String EVENT_HEADER="physical_id\ttable_id\ttable_line\tcontext_id\torigin_genome\torigin_id\tmodel_from\tmodel_to\tlocal_from\tlocal_to\tlocal_strand\ttruncation\tinclusion\tbit_score\tevalue";
	private static final String OVERLAP_HEADER="physical_id\tgenome\toriginal_id\toverlap_bases\tsame_strand\tannotation_queried";
}
