package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import prok.NcrnaCmHitProjector;
import prok.NcrnaCmLocusAudit;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Independent finite geometry and relationship checks for the saved-hit projector.
 * Uses synthetic metadata and Infernal tables, never scientific model provenance.
 * All files stay in the explicitly supplied fresh work directory for inspection.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCmHitProjectorTest{

	public static void main(final String[] args){
		boolean enabled=false;
		assert(enabled=true);
		check(enabled, "Run these fixtures with assertions enabled");
		check(args.length==1 && args[0].startsWith("work="), "Required work=fresh_directory");
		final File work=new File(args[0].substring(5));
		check(work.mkdir(), "Fixture directory must be new");
		geometry();
		final String[] loci={locus("A", '+', 100, 109), locus("B", '+', 120, 129),
			locus("N", '+', 102, 105), locus("O", '-', 104, 106), locus("M", '-', 200, 209),
			locus("Z", '+', 1, 10), locus("W", '+', 90, 140)};
		final String[] contexts={context("qA", "A", '+', 100, 109, 10, 30),
			context("qB", "B", '+', 120, 129, 30, 10), context("qM", "M", '-', 200, 209, 7, 13),
			context("qZ", "Z", '+', 1, 10, 0, 5), context("qW", "W", '+', 90, 140, 10, 10)};
		write(new File(work, "loci.tsv"), NcrnaCmLocusAudit.HEADER, loci);
		write(new File(work, "contexts.tsv"), NcrnaCmHitProjector.CONTEXT_HEADER, contexts);
		final String model1="RF00504", model2="RF03166";
		writeCm(new File(work, "one.tbl"), new String[]{hit("qA", model1, 11, 20, '+', "no", "!", "40"),
			hit("qA", model1, 31, 40, '+', "no", "!", "40"), hit("qB", model1, 11, 20, '+', "no", "!", "41"),
			hit("qB", model1, 31, 40, '+', "5'", "?", "1"), hit("qM", model1, 8, 17, '+', "no", "!", "40"),
			hit("qM", model1, 17, 8, '-', "no", "!", "40"), hit("qW", model1, 21, 30, '+', "no", "!", "40"),
			hit("qW", model1, 41, 50, '+', "no", "!", "40")});
		writeCm(new File(work, "two.tbl"), new String[]{hit("qA", model2, 11, 20, '+', "no", "?", "1")});
		// Synthetic identifiers are declared fixture metadata, not hashes of real models.
		final String[] tables={"one\tRF00504\t00000000000000000000\t94\t32\tTEST\t"+new File(work, "one.tbl").getPath(),
			"two\tRF03166\t11111111111111111111\t177\t110\tTEST\t"+new File(work, "two.tbl").getPath()};
		write(new File(work, "tables.tsv"), NcrnaCmHitProjector.TABLE_HEADER, tables);
		run(work, "loci.tsv", "contexts.tsv", "tables.tsv", "first");
		check(rows(new File(work, "first/summary.tsv")).get(1).equals("7\t5\t2\t9\t5"), "Nine events must become five physical intervals, retaining all populations");
		final ArrayList<String> hits=rows(new File(work, "first/hits.tsv"));
		check(has(hits, "\t+\tRF00504\t00000000000000000000\t100\t109\t3\t"), "Three observations of the first aptamer collapse without changing its bounds");
		check(has(hits, "\t+\tRF00504\t00000000000000000000\t120\t129\t3\t"), "The adjacent aptamer remains a separate physical interval");
		check(has(hits, "\t-\tRF00504\t00000000000000000000\t200\t209\t1\t"), "Minus context plus hit must project to genomic minus");
		check(has(hits, "\t+\tRF00504\t00000000000000000000\t200\t209\t1\t"), "Minus context reverse hit must project to genomic plus");
		check(has(hits, "\t+\tRF03166\t11111111111111111111\t100\t109\t1\t"), "Cross-model evidence is a separate physical key");
		check(has(rows(new File(work, "first/events.tsv")), "\t5'\t?\t1\t"), "Truncated below-GA excluded hits remain raw evidence");
		check(has(rows(new File(work, "first/queries.tsv")), "qZ\tg1\tZ\t1\t15\t0\t5\ttrue\tfalse\t0\t0"), "Clipped zero-hit queried contexts must remain explicit");
		check(has(rows(new File(work, "first/queries.tsv")), "qM\tg1\tM\t187\t216\t7\t13\tfalse\tfalse\t2\t2"), "Unequal sense flanks must swap on a minus context");
		check(has(rows(new File(work, "first/queries.tsv")), "qW\tg1\tW\t80\t150\t10\t10\tfalse\tfalse\t2\t3"), "One source spanning two aptamers retains both events and independent cross-model overlaps");
		final ArrayList<String> neighbors=rows(new File(work, "first/neighbors.tsv"));
		check(has(neighbors, "qA\tg1\tA\tOVERLAP\tN\tNA\t4\ttrue\tfalse"), "Nested unqueried same-strand neighbor must survive immediate-neighbor lookup");
		check(has(neighbors, "qA\tg1\tA\tOVERLAP\tO\tNA\t3\tfalse\tfalse"), "Opposite-strand overlap must remain unqueried evidence");
		check(has(neighbors, "qA\tg1\tA\tPREDECESSOR\tW\t-41\t10\ttrue\ttrue"), "Genomic predecessor uses signed inclusive gap, even for nesting");
		write(new File(work, "loci-reversed.tsv"), NcrnaCmLocusAudit.HEADER, reversed(loci));
		write(new File(work, "contexts-reversed.tsv"), NcrnaCmHitProjector.CONTEXT_HEADER, reversed(contexts));
		write(new File(work, "tables-reversed.tsv"), NcrnaCmHitProjector.TABLE_HEADER, reversed(tables));
		run(work, "loci-reversed.tsv", "contexts-reversed.tsv", "tables-reversed.tsv", "reordered");
		for(final String name:new String[]{"summary.tsv", "hits.tsv", "events.tsv", "queries.tsv", "overlaps.tsv", "neighbors.tsv", "searches.tsv"}){
			check(rows(new File(work, "first/"+name)).equals(rows(new File(work, "reordered/"+name))), "Reordering input records changed "+name);
		}
		write(new File(work, "missing.tsv"), NcrnaCmHitProjector.CONTEXT_HEADER,
			new String[]{context("qX", "ABSENT", '+', 100, 109, 10, 10)});
		expectFailure(()->run(work, "loci.tsv", "missing.tsv", "tables.tsv", "missing-output"), "Unknown context origin must fail");
		write(new File(work, "duplicate.tsv"), NcrnaCmHitProjector.CONTEXT_HEADER, new String[]{contexts[0], contexts[0]});
		expectFailure(()->run(work, "loci.tsv", "duplicate.tsv", "tables.tsv", "duplicate-output"), "Duplicate context IDs must fail");
		System.out.println("PASS NcrnaCmHitProjectorTest geometry_all_events_physical_keys_neighbors_zero_hits_reordering_invalid_maps");
	}

	private static void geometry(){
		final int[] out=new int[3];
		NcrnaCmHitProjector.project(91, 130, '+', 10, 19, '+', out);
		check(Arrays.equals(out, new int[]{100, 109, '+'}), "Plus inclusive geometry");
		NcrnaCmHitProjector.project(190, 219, '-', 11, 20, '+', out);
		check(Arrays.equals(out, new int[]{200, 209, '-'}), "Minus inclusive geometry");
		NcrnaCmHitProjector.project(1, 12, '+', 1, 2, '+', out);
		check(Arrays.equals(out, new int[]{1, 2, '+'}), "Clipped physical edge must remain one-based");
		expectFailure(()->NcrnaCmHitProjector.project(1, 12, '+', 0, 2, '+', out), "Zero local coordinate must fail");
		expectFailure(()->NcrnaCmHitProjector.project(1, 12, '+', 1, 13, '+', out), "Out-of-context endpoint must fail");
		expectFailure(()->NcrnaCmHitProjector.project(1, 12, '+', 2, 1, '+', out), "Hit strand/order mismatch must fail");
	}

	private static void run(final File work, final String loci, final String contexts, final String tables, final String out){
		NcrnaCmHitProjector.main(new String[]{"in="+new File(work, loci), "contexts="+new File(work, contexts),
			"tables="+new File(work, tables), "out="+new File(work, out), "expectedloci=7", "expectedcontexts=5", "expectedhits=9", "t=1"});
	}

	private static String locus(final String id, final char strand, final int start, final int stop){
		return id+"\tg1\tGCF_TEST.1\t1\tbacteria\ttest\tchr\t"+strand+'\t'+start+'-'+stop+'\t'+(stop-start+1)
			+"\t1\t"+id+"\tfalse\tfalse\tANNOTATION_FULL_COPY\t1000\tfalse\t0\t0\t0\t1";
	}

	private static String context(final String id, final String original, final char strand, final int start, final int stop, final int left, final int right){
		return id+'\t'+original+"\tg1\tchr\t"+strand+'\t'+start+'\t'+stop+"\t1000\t"+left+'\t'+right;
	}

	private static String hit(final String query, final String model, final int from, final int to,
			final char strand, final String trunc, final String inclusion, final String score){
		return query+" - model "+model+" cm 1 10 "+from+' '+to+' '+strand+' '+trunc+" 1 0.5 0 "+score+" 1e-5 "+inclusion+" fixture";
	}

	private static void writeCm(final File file, final String[] hits){
		try(Output output=new Output(file.getPath())){
			for(final String hit:hits){output.writer.println(hit);}
			output.writer.println("# Program: cmsearch\n# Version: 1.1.5 fixture\n# [ok]");
		}
	}

	private static void write(final File file, final String header, final String[] rows){
		try(Output output=new Output(file.getPath())){
			output.writer.println(header);
			for(final String row:rows){output.writer.println(row);}
		}
	}

	private static ArrayList<String> rows(final File file){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input input=new Input(file.getPath())){
			final ByteFile reader=input.reader;
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				result.add(new String(line, java.nio.charset.StandardCharsets.UTF_8));
			}
		}
		return result;
	}

	private static String[] reversed(final String[] rows){
		final String[] result=rows.clone();
		for(int i=0; i<rows.length; i++){result[i]=rows[rows.length-1-i];}
		return result;
	}

	private static boolean has(final ArrayList<String> rows, final String fragment){
		for(final String row:rows){if(row.contains(fragment)){return true;}}
		return false;
	}

	private static void expectFailure(final Runnable action, final String message){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, message);
	}

	private static void check(final boolean ok, final String message){
		if(!ok){throw new AssertionError(message);}
	}
}
