package test.prok;

import prok.Orf;
import prok.ProkObject;
import prok.ProteinSupport;
import prok.ProteinSupportFate;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;

import prot.HbmTailGate;
import prot.ProteinTailAssigner;

/** Literal accounting/coordinate fixtures; real caller-hook parity is a separate test.
 * @author Sayu
 */
public final class ProteinSupportFateTest {

	public static void main(String[] args) throws Exception{
		boolean assertions=false; assert(assertions=true);
		check(assertions, "Run diagnostic fixtures with assertions enabled");
		coordinatesAndStages();
		conservation();
		check(ProteinSupportFate.contigToken("fixture:x\tcomment").equals("fixture:x"), "GFF token stops at tab");
		expect(new Action(){public void run(){ProteinSupportFate.contigToken(" comment");}}, "empty GFF token");
		System.err.println("PROTEIN_SUPPORT_FATE_TEST_PASS strands=true stages=true identity=true conservation=true overlaps=true stable_keys=true");
	}

	private static void coordinatesAndStages() throws Exception{
		final ProteinSupportFate.Entry empty=entry(0, 0, 299);
		empty.rejectedLength=2; empty.rejectedStartScore=3; empty.rejectedLimit=1;
		empty.retained(list());
		String[] row=row(empty, list());
		eq(row, "stop_id", "fixture%3Ax:0:300"); eq(row, "fate", "no_starts");
		eq(row, "rejected_length", "2"); eq(row, "rejected_start_score", "3"); eq(row, "rejected_limit", "1");
		eq(row, "admitted", "0"); eq(row, "best_native", ".");
		eq(row, "identity", "80"); eq(row, "family_id", "17"); eq(row, "score64", "1234");
		eq(row, "raw_stop_score", "-9.0"); // Inactive raw stop score must not invent a rejection stage.

		final ProteinSupportFate.Entry pruned=entry(0, 0, 299);
		final Orf p=orf(0, 30, 299, ProkObject.CDS);
		add(pruned, p, 10); pruned.rejected(p, "middle_score"); pruned.retained(list());
		row=row(pruned, list());
		eq(row, "fate", "soft_pruned"); eq(row, "soft_rejected", "1");
		eq(row, "first_pruning_predicates", "middle_score:1");

		final ProteinSupportFate.Entry dp=entry(0, 0, 299);
		final Orf a=orf(0, 30, 299, ProkObject.CDS), b=orf(0, 60, 299, ProkObject.CDS);
		a.orfScore=70; a.proteinBonus=50;
		b.orfScore=60; b.proteinBonus=50;
		add(dp, a, 20); add(dp, b, 10); dp.retained(list(a, b));
		final Orf rna=orf(1, 700, 849, ProkObject.RNA); // Genomic151..300, opposite strand.
		row=row(dp, list(rna));
		eq(row, "admitted", "2"); eq(row, "retained", "2"); eq(row, "selected", "0");
		eq(row, "fate", "dp_not_selected"); eq(row, "best_start", "31"); eq(row, "best_native", "20.0");
		eq(row, "best_bonus", "50.0"); eq(row, "best_total", "70.0");
		eq(row, "observed_selected_overlaps", "fixture%3Ax:1:151:300:7:150:1");
		eq(row, "overlap_basis", "best_retained");
		eq(row(dp, list()), "observed_selected_overlaps", ".");
		eq(row(dp, list()), "fate", "dp_not_selected"); // No local overlap is not evidence of selection.
		eq(row(dp, list(b)), "fate", "selected"); eq(row(dp, list(b)), "selected", "1");
		// A later mutable DP score must not rewrite admission evidence.
		a.orfScore=999; eq(row(dp, list()), "best_total", "70.0");

		final ProteinSupportFate.Entry reverse=entry(1, 0, 299);
		final Orf r=orf(1, 30, 299, ProkObject.CDS);
		add(reverse, r, 10); reverse.retained(list(r));
		row=row(reverse, list(r));
		eq(row, "stop_id", "fixture%3Ax:1:701"); eq(row, "terminal", "701");
		eq(row, "longest_start", "701"); eq(row, "longest_end", "1000");
		eq(row, "best_start", "701"); eq(row, "best_end", "970");
		final Orf plusRna=orf(0, 699, 799, ProkObject.RNA);
		eq(row(reverse, list(plusRna)), "observed_selected_overlaps", "fixture%3Ax:0:700:800:7:100:1");
		final String ordered=ProteinSupportFate.render(entries(empty, reverse), list(r), 60, 30).toString();
		check(ordered.equals(ProteinSupportFate.render(entries(reverse, empty), list(r), 60, 30).toString()),
			"Stable stop ordering must not depend on entry encounter order");
	}

	private static void conservation() throws Exception{
		final ProteinSupportFate.Entry e=entry(0, 0, 299);
		final Orf a=orf(0, 30, 299, ProkObject.CDS);
		expect(new Action(){public void run(){e.admitted(a);}}, "missing annotation evidence");
		add(e, a, 10);
		expect(new Action(){public void run() throws Exception{add(e, a, 10);}}, "duplicate admission");
		expect(new Action(){public void run(){e.retained(list());}}, "unaccounted candidate");
		e.rejected(a, "orf_score");
		expect(new Action(){public void run(){e.rejected(a, "rank_limit");}}, "double pruning");
		expect(new Action(){public void run(){e.retained(list(a));}}, "retaining a rejected candidate");
		e.retained(list());
		expect(new Action(){public void run(){row(e, list(a));}}, "selected pruned candidate");
		expect(new Action(){public void run(){e.retained(list());}}, "second seal");
		expect(new Action(){public void run(){ProteinSupportFate.render(entries(e, e), list(), 60, 30);}}, "duplicate stop key");
		final ProteinSupportFate.Entry unsealed=entry(0, 0, 299);
		expect(new Action(){public void run(){row(unsealed, list());}}, "unsealed entry");
		final ProteinSupportFate.Entry selected=entry(0, 0, 299);
		final Orf s=orf(0, 30, 299, ProkObject.CDS), impostor=orf(0, 30, 299, ProkObject.CDS);
		add(selected, s, 10); selected.retained(list(s));
		expect(new Action(){public void run(){row(selected, list(impostor));}}, "equal coordinates cannot substitute another object");
		expect(new Action(){public void run(){row(selected, list(s, s));}}, "duplicate final selection");
		final ProteinSupportFate.Entry wrong=entry(1, 0, 299);
		expect(new Action(){public void run() throws Exception{add(wrong, s, 10);}}, "wrong strand/stop");
	}

	private static ProteinSupportFate.Entry entry(int strand, int start, int end) throws Exception{
		final Orf longest=orf(strand, start, end, ProkObject.CDS);
		if(strand==1){longest.flip();}
		final HbmTailGate.Metrics metrics=new HbmTailGate.Metrics(); metrics.identity=80f;
		// Retain reflective construction to exercise the original constructor signature.
		final Constructor<ProteinTailAssigner.Hit> ctor=ProteinTailAssigner.Hit.class.getDeclaredConstructor(
			int.class, int.class, int.class, String.class, int.class, int.class, int.class,
			HbmTailGate.Metrics.class, int.class, int[].class);
		ctor.setAccessible(true);
		final ProteinTailAssigner.Hit hit=ctor.newInstance(0, 2, 17, "rep", 1234, -1, -1, metrics, 1, new int[0]);
		return new ProteinSupportFate.Entry(longest, new ProteinSupport.StopEvidence(hit, null, 11, 90), -9f, 80, 0.1f);
	}

	/** Injects an already-scored fixture without loading a real model pack.
	 * Integration parity must separately exercise admitted() after real annotate(). */
	private static void add(ProteinSupportFate.Entry e, Orf o, float nativeScore) throws Exception{
		final Method method=ProteinSupportFate.Entry.class.getDeclaredMethod("add", Orf.class, float.class);
		method.setAccessible(true);
		try{method.invoke(e, o, nativeScore);}
		catch(InvocationTargetException ex){
			if(ex.getCause() instanceof RuntimeException){throw (RuntimeException)ex.getCause();}
			throw ex;
		}
	}
	private static Orf orf(int strand, int start, int end, int type){
		final byte[] bases=new byte[1000]; Arrays.fill(bases, (byte)'A');
		final Orf o=new Orf("fixture:x comment", start, end, strand, 0, bases, strand==1, type, false, false);
		o.startScore=0.5f; o.stopScore=-9; o.kmerScore=20; o.orfScore=10;
		return o;
	}
	private static String[] row(ProteinSupportFate.Entry e, ArrayList<Orf> selected){
		final String[] fields=ProteinSupportFate.render(entries(e), selected, 60, 30).toString().trim().split("\t", -1);
		check(fields.length==ProteinSupportFate.HEADER.split("\t").length, "Every TSV row must match the header width");
		return fields;
	}
	private static void eq(String[] row, String column, String value){
		final String[] columns=ProteinSupportFate.HEADER.split("\t");
		for(int i=0; i<columns.length; i++){if(columns[i].equals(column)){
			check(row[i].equals(value), column+": expected "+value+", observed "+row[i]); return;
		}}
		throw new AssertionError("Unknown fixture column "+column);
	}
	private static ArrayList<Orf> list(Orf... orfs){return new ArrayList<Orf>(Arrays.asList(orfs));}
	private static ArrayList<ProteinSupportFate.Entry> entries(ProteinSupportFate.Entry... es){return new ArrayList<ProteinSupportFate.Entry>(Arrays.asList(es));}
	private interface Action{void run() throws Exception;}
	private static void expect(Action action, String why) throws Exception{
		try{action.run();}catch(IllegalStateException expected){return;}
		throw new AssertionError("Expected rejection: "+why);
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private ProteinSupportFateTest(){}
}
