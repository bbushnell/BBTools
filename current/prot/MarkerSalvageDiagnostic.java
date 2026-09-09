package prot;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/**
 * Read-only threshold/yield diagnostic for the marker-salvage rule (D33/D34/D38: relax the
 * single-copy/presence criteria to rescue phylum groups with too few tracked families -- floor
 * 20, aim for up to 40, without dropping below 70%). Reuses {@link MarkerSelector}'s existing
 * per-org accumulation ({@link MarkerSelector#accumulateGroups}) UNCHANGED -- same two inputs
 * (sparse per-org cache + taxpgm.tsv), same {@link MarkerSelector.Group} present/single counts --
 * and reports ONLY counts. It never writes a family-rank list, never picks a threshold, and never
 * caps or selects a "winner" per group (Yoimiya's review, 2026-09-08: reporting counts above 40
 * is required, not truncated; the aim-40/floor-20/70%-floor language is decision CONTEXT printed
 * in the header, not logic this tool applies).
 *
 * <p><b>Two output tables, not one -- because "the threshold" is itself ambiguous (Yoimiya's
 * review, 2026-09-08).</b> Brian's own words name relaxing presence (P/N, {@code minPrev}) and/or
 * relaxing single-copy-among-present (S/P, {@code minSc}) as the two knobs, but his phrase "more
 * than 40 with &gt;=0.95 having 1 copy" can also be read as a THIRD ratio: single-copy-of-the-
 * WHOLE-population, S/N -- mathematically the product (P/N)*(S/P), a composite this tool does not
 * currently threshold on at all. Rather than commit to one reading, this diagnostic emits:</p>
 * <ol>
 * <li><b>Raw per-(group,rank) table</b> ({@code <out>_raw.tsv}): `n_orgs` (N), `present` (P),
 * `single` (S) for every family rank in every accumulated group. This is the "enough denominator/
 * count information to resolve the interpretation later" artifact -- P/N, S/P, and S/N are all
 * exactly recoverable from these three raw integers, at ANY threshold, without rerunning
 * anything.</li>
 * <li><b>Grid-count table</b> ({@code <out>_grid.tsv}): the approved 6x6 sweep of
 * `n_markers(group, minPrev, minSc)` under the CURRENT (P/N AND S/P as independent constraints)
 * reading only -- a convenience summary, explicitly labeled as one interpretation among the three
 * above, not a resolution of which one is intended.</li>
 * </ol>
 *
 * <p>Every accumulated group appears in both tables, explicitly labeled by `kind`:
 * `phylum` (a true phylum group with `n_orgs&gt;=minOrgs`), `all_domain` (the `ALL_&lt;domain&gt;`
 * cross-phylum summary rows -- never a shipped subnet, {@link MarkerSelector#writeMarkerSets}
 * skips writing a rank file for these), or `excluded_minorgs` (`n_orgs&lt;minOrgs` -- these ARE
 * accumulated by {@link MarkerSelector#accumulateGroups} exactly like any other group, and are
 * only ever filtered OUT at {@link MarkerSelector}'s own output stage; this diagnostic does not
 * drop them, corrected after Yoimiya's review from an earlier plan draft that wrongly said "no
 * group ever formed" for these).</p>
 *
 * <p>Usage: {@code java prot.MarkerSalvageDiagnostic perorg=&lt;perorg_sparse.tsv&gt;
 * taxpgm=&lt;taxpgm.tsv&gt; out=&lt;prefix&gt; [minorgs=3] [prevgrid=0.95,0.90,0.85,0.80,0.75,0.70]
 * [scgrid=0.95,0.90,0.85,0.80,0.75,0.70] [ow=t]}</p>
 *
 * <p><b>Scope (2026-09-08, Yoimiya's assignment):</b> synthetic-fixture increment only -- no real
 * corpus, no cluster execution, no marker set (re)generated. No threshold chosen; the aim-40/
 * floor-20/70%-floor numbers are printed as header context only.</p>
 *
 * @author Eru
 */
public final class MarkerSalvageDiagnostic {

	public static void main(String[] args){
		Timer t=new Timer();
		MarkerSalvageDiagnostic x=new MarkerSalvageDiagnostic(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public MarkerSalvageDiagnostic(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		String prevGridStr=DEFAULT_GRID, scGridStr=DEFAULT_GRID;
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("perorg") || a.equals("in")){perorgFile=b;}
			else if(a.equals("taxpgm")){taxpgmFile=b;}
			else if(a.equals("out")){outPrefix=b;}
			else if(a.equals("minorgs")){minOrgs=Integer.parseInt(b);}
			else if(a.equals("prevgrid")){prevGridStr=b;}
			else if(a.equals("scgrid")){scGridStr=b;}
			else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
			else{outstream.println("Unknown parameter "+arg); assert(false) : "Unknown parameter "+arg;}
		}
		assert(perorgFile!=null) : "perorg= is required.";
		assert(taxpgmFile!=null) : "taxpgm= is required.";
		assert(outPrefix!=null) : "out= is required.";
		prevGrid=parseGrid(prevGridStr);
		scGrid=parseGrid(scGridStr);
	}

	static final String DEFAULT_GRID="0.95,0.90,0.85,0.80,0.75,0.70";

	static double[] parseGrid(String s){
		final String[] parts=s.split(",");
		final double[] out=new double[parts.length];
		for(int i=0; i<parts.length; i++){
			out[i]=Double.parseDouble(parts[i]);
			//Double.parseDouble("NaN") succeeds, and EVERY comparison with NaN (including <0 and
			//>1) is false, so a range check alone silently admits NaN (Yoimiya/Root's review,
			//2026-09-08). Reject all nonfinite values (NaN, +-Infinity) explicitly, first.
			if(!Double.isFinite(out[i])){throw new RuntimeException("Nonfinite grid value: "+parts[i]);}
			if(out[i]<0 || out[i]>1){throw new RuntimeException("Grid value out of [0,1]: "+parts[i]);}
		}
		return out;
	}

	void process(Timer t){
		final HashMap<Integer, String> tid2phylum=MarkerSelector.loadTaxpgm(taxpgmFile);
		outstream.println("taxpgm tids: "+tid2phylum.size());

		final HashMap<String, MarkerSelector.Group> groups=
			MarkerSelector.accumulateGroups(perorgFile, tid2phylum, outstream);

		final ArrayList<String> ordered=new ArrayList<String>(groups.keySet());
		ordered.sort(Comparator.naturalOrder());

		writeRawTable(outPrefix+"_raw.tsv", ordered, groups);
		writeGridTable(outPrefix+"_grid.tsv", ordered, groups);

		t.stop();
		outstream.println("Time: \t"+t);
	}

	static String classify(String name, int nOrgs, int minOrgs){
		if(name.startsWith("ALL_")){return "all_domain";}
		if(nOrgs<minOrgs){return "excluded_minorgs";}
		return "phylum";
	}

	/** Raw per-(group,rank) N/P/S -- the artifact sufficient to recompute P/N, S/P, or S/N (the
	 *  product) at any threshold later, without rerunning anything (Yoimiya's review, 2026-09-08). */
	void writeRawTable(String path, ArrayList<String> ordered, HashMap<String, MarkerSelector.Group> groups){
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.TXT, null, true, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder(1<<16);
		bb.append("#DECISION CONTEXT ONLY, not applied by this tool: D33/D34/D38 -- floor 20, aim").nl();
		bb.append("#for up to 40, without dropping the threshold below 0.70. No selection performed here.").nl();
		bb.append("#group\tkind\tn_orgs\trank\tP_present\tS_single").nl();
		bsw.print(bb);
		bb.clear();
		for(String name : ordered){
			final MarkerSelector.Group g=groups.get(name);
			final String kind=classify(name, g.nOrgs, minOrgs);
			final ArrayList<Integer> ranks=new ArrayList<Integer>(g.present.keySet());
			ranks.sort(null);
			for(int rank : ranks){
				final int p=g.present.get(rank);
				final int s=g.single.getOrDefault(rank, 0);
				bb.append(name).tab().append(kind).tab().append(g.nOrgs).tab()
					.append(rank).tab().append(p).tab().append(s).nl();
			}
			bsw.print(bb);
			bb.clear();
		}
		if(bsw.poisonAndWait()){throw new RuntimeException("ByteStreamWriter reported an I/O error writing "+path);}
	}

	/** Grid-count table: n_markers(group, minPrev, minSc) under the P/N-AND-S/P-independent
	 *  reading only, swept over the full grid, uncapped -- counts above 40 are reported, never
	 *  truncated or selected as a "winner" (Yoimiya's review, 2026-09-08). */
	void writeGridTable(String path, ArrayList<String> ordered, HashMap<String, MarkerSelector.Group> groups){
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.TXT, null, true, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder(1<<16);
		bb.append("#DECISION CONTEXT ONLY, not applied by this tool: D33/D34/D38 -- floor 20, aim").nl();
		bb.append("#for up to 40, without dropping the threshold below 0.70. No selection performed here.").nl();
		bb.append("#Counts under ONE reading (P/N=minPrev AND S/P=minSc, both independently required);").nl();
		bb.append("#see the raw table for the S/N (product) reading. Counts above 40 are NOT truncated.").nl();
		bb.append("#group\tkind\tn_orgs\tminPrev\tminSc\tn_markers").nl();
		bsw.print(bb);
		bb.clear();
		for(String name : ordered){
			final MarkerSelector.Group g=groups.get(name);
			final String kind=classify(name, g.nOrgs, minOrgs);
			for(double minPrev : prevGrid){
				for(double minSc : scGrid){
					final int n=countQualifying(g, minPrev, minSc);
					//Emit the EXACT double's round-trip text (Double.toString is specified to
					//produce the shortest decimal that parses back to the identical double), not
					//a fixed 2-decimal rounding -- a rounded LABEL next to an UNROUNDED computation
					//let e.g. 0.949 print as "0.95" while countQualifying still used 0.949
					//(Yoimiya/Root's review, 2026-09-08: label and computation must agree exactly).
					bb.append(name).tab().append(kind).tab().append(g.nOrgs).tab()
						.append(Double.toString(minPrev)).tab().append(Double.toString(minSc)).tab().append(n).nl();
				}
			}
			bsw.print(bb);
			bb.clear();
		}
		if(bsw.poisonAndWait()){throw new RuntimeException("ByteStreamWriter reported an I/O error writing "+path);}
	}

	/** Same qualification rule as {@link MarkerSelector#writeMarkerSets} (P&gt;=minPrev*n AND
	 *  S&gt;=minSc*P), reapplied at arbitrary thresholds instead of one fixed pair -- counts only,
	 *  no rank retained, no rank file ever written by this method. */
	static int countQualifying(MarkerSelector.Group g, double minPrev, double minSc){
		final int n=g.nOrgs;
		int count=0;
		for(Map.Entry<Integer, Integer> e : g.present.entrySet()){
			final int rank=e.getKey(), p=e.getValue();
			if(p<=0){continue;}
			if(p<minPrev*n){continue;}
			final int s=g.single.getOrDefault(rank, 0);
			if(s<minSc*p){continue;}
			count++;
		}
		return count;
	}

	private String perorgFile=null;
	private String taxpgmFile=null;
	private String outPrefix=null;
	private int minOrgs=3;
	private final double[] prevGrid;
	private final double[] scGrid;
	private boolean overwrite=true;
	private java.io.PrintStream outstream=System.err;
}
