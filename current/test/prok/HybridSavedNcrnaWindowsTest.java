package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.HybridSavedNcrnaWindows;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.RankSavedNcrnaWindows;
import shared.Shared;
import structures.ByteBuilder;

/** Literal saved native-count binding and hybrid integration, without genomes.
 * @author Brian, Raiden
 */
public final class HybridSavedNcrnaWindowsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HybridSavedNcrnaWindowsTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown hybrid fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_fixture_directory and normal runtime options");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Hybrid fixture invariants require assertions");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Hybrid fixtures must not overwrite prior evidence");
			DynamicDemiLog.setExponent(4);
			run(root.getPath());
			pp.outstream.println("NCRNA_SAVED_HYBRID_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	private static void run(final String root){
		final RankSavedNcrnaWindowsTest.Fixture good=new RankSavedNcrnaWindowsTest.Fixture(root, "good");
		grid(good);
		HybridSavedNcrnaWindows.run(config(good, "one", 1));
		HybridSavedNcrnaWindows.run(config(good, "two", 2));
		for(final String suffix:new String[]{".methods.tsv", ".summary.tsv", ".loci.tsv.gz", ".input_summary.tsv", ".last_losses.tsv"}){
			equal(good.dir+"/one"+suffix, good.dir+"/two"+suffix);
		}
		try(Input in=new Input(good.dir+"/one.summary.tsv")){
			check(in.reader.nextLine()!=null, "A hybrid summary requires its schema header");
			for(final String stratum:new String[]{"NOT_SEEN","SEEN","UNKNOWN"}){
				for(final String method:new String[]{"NATIVE","SKETCH","HYBRID"}){
					final String row=text(in.reader.nextLine());
					final String values=stratum.equals("UNKNOWN") ? "\t0\t0\tNA\t0\t0\t0\t-1\t-1\t-1\t0\t0\t0"
						: "\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0";
					check(row.equals(stratum+'\t'+method+values+"\tSAVED_LOCUS_MINIMA_NOT_CALLER_WORK"),
						"Literal model0 A windows rank first in all orders; unrankable/native-zero loci stay outside common counts");
				}
			}
			check(in.reader.nextLine()==null, "Hybrid summaries cannot append pooled or duplicate strata");
		}
		try(Input in=new Input(good.dir+"/one.loci.tsv.gz")){
			check(in.reader.nextLine()!=null, "Per-parent hybrid evidence requires a header");
			final LineParser1 p=new LineParser1('\t');
			int rows=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.terms()==9, "Hybrid locus rows must retain all three ranks and source disposition");
				final boolean common=p.termEquals("p0", 3), zero=p.termEquals("p4", 3);
				check(p.parseInt(5)==(common ? 1 : -1) && p.parseInt(6)==(common ? 1 : zero ? 2 : -1)
					&& p.parseInt(7)==(common ? 1 : zero ? 2 : -1), "Literal event minima and native-zero reachability must survive the complete pipeline");
				check(p.termEquals("beta", 0) || Arrays.equals(p.parseByteArray(0), PANEL), "Non-ASCII panel names must retain exact UTF-8");
				rows++;
			}
			check(rows==6, "All six selected targets survive, including ambiguous/gap/no-window records");
		}
		try(Input in=new Input(good.dir+"/two.work.tsv")){
			boolean workers=false, nativeRows=false, nativeTuples=false, canonical=false;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				final String value=text(row);
				workers|=value.equals("workers_with_scored_windows\t2");
				nativeRows|=value.equals("cached_native_model_rows\t22");
				nativeTuples|=value.equals("native_event_tuples_checked\t11");
				canonical|=value.equals("canonical_event_tuples_checked\t176");
			}
			check(workers && nativeRows && nativeTuples && canonical, "Both scoring workers must check the complete native/cache/control fixture");
		}
		try(Input in=new Input(good.dir+"/one.last_losses.tsv")){
			check(in.reader.nextLine()!=null && in.reader.nextLine()==null, "The two baseline-reachable fixture inputs retain their only locus");
		}
		final String[] faults={"native_missing","native_duplicate","native_name","native_count","native_tuple","tuple","grid_minimum","grid_binding","grid_method"};
		final String[] reasons={"lose zero-support models","original ordinal/full-header","original ordinal/full-header","MAX_COUNT=65535",
			"Cached native rank/tie/model parity","Canonical hybrid-control","selected new sketch must reproduce","exact original input","supported forward"};
		check(faults.length==reasons.length, "Every malformed cache must have one intended rejection");
		for(int i=0; i<faults.length; i++){
			final RankSavedNcrnaWindowsTest.Fixture bad=new RankSavedNcrnaWindowsTest.Fixture(root, faults[i]);
			grid(bad);
			try{HybridSavedNcrnaWindows.run(config(bad, "rejected", 2));}
			catch(IllegalArgumentException | IllegalStateException failure){
				Throwable cause=failure;
				while(cause.getCause()!=null){cause=cause.getCause();}
				check(cause.getMessage().contains(reasons[i]), "Wrong hybrid rejection for "+faults[i]+": "+cause);
				check(!new File(bad.dir+"/rejected.summary.tsv").exists(), "A failed native/control/minimum binding cannot publish a summary");
				continue;
			}
			throw new AssertionError("Malformed native/grid evidence was accepted: "+faults[i]);
		}
	}

	private static HybridSavedNcrnaWindows.Config config(final RankSavedNcrnaWindowsTest.Fixture f, final String output, final int workers){
		assert(workers>0) : "Hybrid fixture workers must actually score the shared-parent initial windows";
		final RankSavedNcrnaWindows.Config original=f.config(output, workers);
		final HybridSavedNcrnaWindows.Config c=new HybridSavedNcrnaWindows.Config();
		c.panels=original.panels;
		c.selection=original.selection;
		c.grid=f.dir+"/grid";
		c.output=original.output;
		c.members=original.members;
		c.ledger=original.ledger;
		c.assignments=original.assignments;
		c.census=original.census;
		c.summary=original.summary;
		c.references=original.references;
		c.method=64;
		c.models=2;
		c.groups=3;
		c.threads=workers;
		return c;
	}

	/** Independent selected-method metadata/minima; never produced by a scorer. */
	private static void grid(final RankSavedNcrnaWindowsTest.Fixture f){
		assert(f!=null) : "Independent grid witnesses borrow the initialized cache fixture's panel/input paths";
		try(Output out=new Output(f.dir+"/grid.methods.tsv")){
			out.writer.println("method\tk\tbuckets\torientation\tmiddle_mask\tscore\trole");
			out.writer.println("0\t9\t0\tFORWARD_NATIVE_INDEX\tfalse\tPOSITIVE_SHARED_KEYS\tSAVED_BASELINE");
			out.writer.println("64\t"+(f.fault.equals("grid_method") ? 6 : 10)+"\t2048\tFORWARD\ttrue\tRAW_ALL_QUERY_KEYS\tNEW_GRID");
		}
		try(Output out=new Output(f.dir+"/grid.loci.tsv.gz")){
			out.writer.println("panel\tinput_slot\tinput_id\tparent\tseen_stratum\tmethod\tbest_rank\tsource_disposition\tinterpretation");
			final ByteBuilder b=new ByteBuilder();
			for(int panel=0; panel<2; panel++){
				for(int parent=0; parent<(panel==0 ? 5 : 1); parent++){
					final String stratum=panel==1 || parent==1 ? "SEEN" : parent==2 ? "UNKNOWN" : "NOT_SEEN";
					final String reason=parent==0 ? "RANKED" : parent==1 ? "AMBIGUOUS_SPANS_ONLY" : parent==2 ? "NO_CAPTURED_WINDOW"
						: parent==3 ? "NO_COMPATIBLE_MODEL" : "ZERO_NATIVE_INDEX_SUPPORT";
					for(final int method:new int[]{0,64}){
						final int rank=panel==0 && parent==0 && method==64 && f.fault.equals("grid_minimum") ? 2
							: parent==0 ? 1 : parent==4 && method>0 ? 2 : -1;
						b.clear().append(panel==0 ? PANEL : BETA).append("\t0\t")
							.append(panel==0 && method==64 && f.fault.equals("grid_binding") ? "wrong" : panel==0 ? "input_a" : "input_b")
							.append("\tp").append(parent).tab().append(stratum).tab().append(method).tab().append(rank).tab().append(reason)
							.append("\tSAVED_LOCUS_MINIMUM_NOT_CALLER_WORK\n");
						out.writer.print(b);
					}
				}
			}
		}
	}
	private static void equal(final String a, final String b){
		try(Input left=new Input(a); Input right=new Input(b)){
			for(byte[] row=left.reader.nextLine(); row!=null; row=left.reader.nextLine()){
				check(Arrays.equals(row, right.reader.nextLine()), "Worker count changed a deterministic hybrid table: "+a);
			}
			check(right.reader.nextLine()==null, "The second worker may not append an extra hybrid row");
		}
	}
	private static String text(final byte[] row){
		check(row!=null, "Hybrid fixture output ended before its literal expected rows");
		return new String(row, StandardCharsets.UTF_8);
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final byte[] PANEL="α-panel".getBytes(StandardCharsets.UTF_8), BETA="beta".getBytes(StandardCharsets.UTF_8);
}
