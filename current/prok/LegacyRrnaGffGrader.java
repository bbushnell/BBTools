package prok;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import fileIO.ByteFile;
import fileIO.FileFormat;

import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;

/** Coordinate-based, family-aware, one-to-one legacy rRNA (5S/16S/18S/23S) GFF grader.
 *
 * <p>Sibling of {@link NcrnaGffFamilyGrader}, authorized 2026-09-02 (Citan) after this
 * project's own audit found no existing tool could grade a real {@code callgenes.sh}
 * legacy-rRNA GFF: {@code comparegff.sh}/{@code CompareGff} classifies subtypes
 * internally ({@code gff.GffLine.prokType()}) but only ever prints pooled aggregate
 * {@code rRNA} metrics, never per-subtype, and carries a known non-one-to-one matching
 * defect (accepts the first qualifying reference per query without consuming it);
 * {@link NcrnaGffFamilyGrader} has the right one-to-one matching engine but its parser
 * is built for the ncRNA-scavenger {@code model:family_...} attribute convention, which
 * a real legacy-rRNA GFF row simply does not use at all (verified against a real sample,
 * {@code ~/Noire/eval_bact/tid_100177_Borreliella_lusitaniae.gff}, this project session).
 *
 * <p><b>Deliberately reuses {@link NcrnaGffFamilyGrader}'s proven one-to-one matcher
 * ({@link NcrnaGffFamilyGrader#matchOneToOne}, {@link NcrnaGffFamilyGrader#Locus},
 * {@link NcrnaGffFamilyGrader#MatchResult}, {@link NcrnaGffFamilyGrader#overlapsReciprocal})
 * unchanged, via package-private cross-class access</b> -- no BBTools file outside this
 * grader, its test class, and its launcher script was modified to build this. The matching
 * invariants (same seqid+strand, reciprocal overlap &gt;=50% of BOTH lengths,
 * maximum-cardinality bipartite matching via Kuhn's algorithm rather than greedy,
 * deterministic input-order processing with overlap-descending edge visitation, and the
 * DUP-vs-FP distinction for an unmatched call overlapping an already-claimed truth locus)
 * are exactly {@link NcrnaGffFamilyGrader}'s own, not re-derived or simplified here.
 *
 * <p>What genuinely differs from {@link NcrnaGffFamilyGrader}, per this session's
 * interface specification ({@code legacy_rrna_comparator_interface_spec_20260902.md}):
 * <ul>
 * <li>The GFF parser: type column must be {@code rRNA} (not the generic-ncRNA {@code RNA}),
 * and the subtype token is matched via the SAME substring convention
 * {@code gff.GffLine.prokType()} already uses -- {@code contains("16S"/"23S"/"18S")}, or
 * {@code contains("5S")} gated by {@code length()<300} -- not a {@code model:}-prefix
 * search. An {@code rRNA}-typed row matching none of these fails loud (an attribution gap,
 * not a different family this grader legitimately ignores -- unlike the generic-ncRNA
 * case, there is no other legitimate family an {@code rRNA}-typed row could belong to).</li>
 * <li>Four in-scope families ({@code 5S}/{@code 16S}/{@code 18S}/{@code 23S}) -- so the
 * per-family counters are a real {@code Map<String,FamilyCounts>} keyed by family name,
 * not two hardcoded field pairs; this does not fall out of the two-family design by simple
 * rename, and is called out as such in the interface spec.</li>
 * <li>The cross-label diagnostic generalizes from a single pairwise A-vs-B count to "was
 * this unmatched in-scope truth locus geometrically overlapped by a call of ANY OTHER
 * in-scope family," reported per family (4 counters) -- not a full pairwise matrix, since
 * nothing in the governing plan asks for pairwise attribution beyond the call-equivalence
 * diff's "subtype change" classification.</li>
 * </ul>
 *
 * <p>Housing decision (real BBTools tool vs. benchmark-local code) explicitly left open by
 * the interface spec; this file is the "real BBTools tool, dev-grade" answer Citan chose.
 *
 * @author G11
 */
public class LegacyRrnaGffGrader {

	public static void main(String[] args) throws Exception {
		String gffPath=null, truthPath=null, negPath=null, out=null;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<0){throw new IllegalArgumentException("Bad argument (expected flag=value): "+arg);}
			final String a=arg.substring(0, eq).toLowerCase(), b=arg.substring(eq+1);
			if(a.equals("gff")){gffPath=b;}
			else if(a.equals("truth")){truthPath=b;}
			else if(a.equals("negseqids")){negPath=b;}
			else if(a.equals("out")){out=b;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(gffPath==null || truthPath==null || out==null){
			throw new IllegalArgumentException("Usage: gff=<calls.gff> truth=<truth.tsv> [negseqids=<file>] out=<prefix>");
		}
		final ArrayList<Locus> truth=loadTruthTsv(truthPath);
		final HashSet<String> negSeqids=(negPath==null ? new HashSet<String>() : NcrnaGffFamilyGrader.loadSeqidSet(negPath));
		final ArrayList<Locus> calls=loadCallsFromGff(gffPath);
		final GradeResult r=grade(truth, calls, negSeqids);
		writeReport(out, r);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Core Data Types        ----------------*/
	/*--------------------------------------------------------------*/

	static final String FAM_5S="5S", FAM_16S="16S", FAM_23S="23S", FAM_18S="18S";
	/** In-scope, scored families; iteration order fixes deterministic report order. */
	static final String[] IN_SCOPE_FAMILIES={FAM_5S, FAM_16S, FAM_18S, FAM_23S};
	static final float MIN_OVERLAP_FRAC=0.5f;
	static final String NA="."; //no-match / not-applicable marker in audit output
	/** Matches gff.GffLine.prokType()'s own 5S-vs-longer-feature disambiguation gate,
	 * carried forward unchanged (not re-derived) per the interface spec. */
	static final int MAX_5S_LEN=300;

	/*--------------------------------------------------------------*/
	/*----------------      One-to-One Matching      ----------------*/
	/*--------------------------------------------------------------*/
	// Deliberately NOT reimplemented here -- see NcrnaGffFamilyGrader.matchOneToOne,
	// NcrnaGffFamilyGrader.overlapsReciprocal, and NcrnaGffFamilyGrader.MatchResult,
	// called directly below via package-private cross-class access.

	/*--------------------------------------------------------------*/
	/*----------------            Grading            ----------------*/
	/*--------------------------------------------------------------*/

	static final class FamilyCounts {
		long tp, fn, fp, dup;
	}

	/** One row of the per-locus audit trail -- either a truth locus or a query call. */
	static final class AuditRow {
		final String locusId;
		final Locus locus;
		final String role; //TRUTH or CALL
		String subtypeOutcome=NA, subtypeMatchId=NA;
		String unionOutcome=NA, unionMatchId=NA;
		String crossLabelOutcome="NONE";
		final ArrayList<String> crossLabelMatchIds=new ArrayList<>();
		AuditRow(String locusId_, Locus locus_, String role_){locusId=locusId_; locus=locus_; role=role_;}
	}

	/** Aggregate grading result, per-family (not two hardcoded field pairs -- see class
	 * javadoc). union = pooled {5S,16S,18S,23S} truth; crossLabelCounts are diagnostic-only,
	 * never folded into union or subtype
	 * TP/FN/FP. auditRows (truth rows then call rows, each input-order-stable) backs the
	 * per-locus .rawcalls.tsv output and is sufficient on its own to reconstruct every
	 * field below by scanning/grouping. */
	static final class GradeResult {
		final Map<String, FamilyCounts> byFamily=new LinkedHashMap<>();
		long unionTP, unionFN, unionFP;
		long negFP;
		final Map<String, Long> crossLabelCounts=new LinkedHashMap<>();
		final ArrayList<AuditRow> auditRows=new ArrayList<>();
		GradeResult(){
			for(String f : IN_SCOPE_FAMILIES){byFamily.put(f, new FamilyCounts()); crossLabelCounts.put(f, 0L);}
		}
	}

	/** Pure, file-I/O-free grading core -- directly fixture-testable, mirroring
	 * NcrnaGffFamilyGrader.grade's own "pure so it can be fixture-tested" design. */
	static GradeResult grade(List<Locus> truth, List<Locus> calls, Set<String> negSeqids){
		final GradeResult r=new GradeResult();

		for(Locus t : truth){
			if(negSeqids.contains(t.seqid)){
				throw new IllegalStateException("Error - truth locus present on a negseqids-listed "
					+"seqid: "+t.seqid+" -- data inconsistency, refusing to grade.");
			}
		}

		final LinkedHashMap<Locus, AuditRow> truthRows=new LinkedHashMap<>();
		for(int i=0; i<truth.size(); i++){
			final Locus t=truth.get(i);
			truthRows.put(t, new AuditRow("T"+i, t, "TRUTH"));
		}
		final LinkedHashMap<Locus, AuditRow> callRows=new LinkedHashMap<>();
		for(int i=0; i<calls.size(); i++){
			final Locus c=calls.get(i);
			callRows.put(c, new AuditRow("C"+i, c, "CALL"));
		}

		// --- neg-seqid calls (any in-scope family) -> NEG_FP ---
		final ArrayList<Locus> negCalls=new ArrayList<>();
		final ArrayList<Locus> posCalls=new ArrayList<>();
		for(Locus c : calls){
			if(negSeqids.contains(c.seqid)){negCalls.add(c);}else{posCalls.add(c);}
		}
		for(Locus c : negCalls){
			final FamilyCounts fc=r.byFamily.get(c.family);
			fc.fp++; r.unionFP++; r.negFP++;
			final AuditRow cr=callRows.get(c);
			cr.subtypeOutcome="NEG_FP"; cr.unionOutcome="FP";
		}

		// --- per-family subtype passes ---
		final Map<String, ArrayList<Locus>> truthByFam=new LinkedHashMap<>();
		final Map<String, ArrayList<Locus>> callsByFam=new LinkedHashMap<>();
		final Map<String, MatchResult> matchByFam=new LinkedHashMap<>();
		for(String f : IN_SCOPE_FAMILIES){
			final ArrayList<Locus> truthFam=filter(truth, f), callsFam=filter(posCalls, f);
			truthByFam.put(f, truthFam); callsByFam.put(f, callsFam);
			final MatchResult mr=NcrnaGffFamilyGrader.matchOneToOne(callsFam, truthFam, MIN_OVERLAP_FRAC);
			matchByFam.put(f, mr);
			r.byFamily.get(f).tp=mr.matched;
			applySubtypePass(callsFam, truthFam, mr, callRows, truthRows, r, f);
		}

		// --- union pass: pooled scored-family truth and calls ---
		final ArrayList<Locus> unionTruth=new ArrayList<>();
		unionTruth.addAll(truth);
		final MatchResult mrU=NcrnaGffFamilyGrader.matchOneToOne(posCalls, unionTruth, MIN_OVERLAP_FRAC);
		r.unionTP=mrU.matched;
		applyUnionPass(posCalls, unionTruth, mrU, callRows, truthRows, r);

		// --- cross-label diagnostic: unmatched in-scope truth vs ANY other in-scope family's calls ---
		for(String f : IN_SCOPE_FAMILIES){
			final MatchResult mr=matchByFam.get(f);
			final ArrayList<Locus> truthFam=truthByFam.get(f);
			for(int ti=0; ti<truthFam.size(); ti++){
				if(mr.truthConsumed[ti]){continue;}
				final Locus t=truthFam.get(ti);
				final AuditRow tr=truthRows.get(t);
				boolean found=false;
				for(String other : IN_SCOPE_FAMILIES){
					if(other.equals(f)){continue;}
					for(Locus c : callsByFam.get(other)){
						if(NcrnaGffFamilyGrader.overlapsReciprocal(c, t, MIN_OVERLAP_FRAC)){
							tr.crossLabelOutcome="CROSS_DETECTED"; tr.crossLabelMatchIds.add(callRows.get(c).locusId);
							final AuditRow cr=callRows.get(c);
							cr.crossLabelOutcome="CROSS_MATCH"; cr.crossLabelMatchIds.add(tr.locusId);
							found=true;
						}
					}
				}
				if(found){r.crossLabelCounts.merge(f, 1L, Long::sum);}
			}
		}

		r.auditRows.addAll(truthRows.values());
		r.auditRows.addAll(callRows.values());
		return r;
	}

	private static void applySubtypePass(List<Locus> callsFam, List<Locus> truthFam, MatchResult mr,
			Map<Locus, AuditRow> callRows, Map<Locus, AuditRow> truthRows, GradeResult r, String family){
		final FamilyCounts fc=r.byFamily.get(family);
		for(int ti=0; ti<truthFam.size(); ti++){
			final Locus t=truthFam.get(ti);
			final AuditRow tr=truthRows.get(t);
			if(mr.truthConsumed[ti]){
				final Locus c=callsFam.get(mr.truthMatchedCallIdx[ti]);
				tr.subtypeOutcome="TP"; tr.subtypeMatchId=callRows.get(c).locusId;
			}else{
				tr.subtypeOutcome="FN"; fc.fn++;
			}
		}
		for(int ci=0; ci<callsFam.size(); ci++){
			final Locus c=callsFam.get(ci);
			final AuditRow cr=callRows.get(c);
			if(mr.callConsumed[ci]){
				final Locus t=truthFam.get(mr.callMatchedTruthIdx[ci]);
				cr.subtypeOutcome="TP"; cr.subtypeMatchId=truthRows.get(t).locusId;
			}else{
				final boolean dup=mr.callIsDuplicate[ci];
				cr.subtypeOutcome=(dup ? "DUP" : "FP");
				fc.fp++; if(dup){fc.dup++;}
			}
		}
	}

	private static void applyUnionPass(List<Locus> posCalls, List<Locus> unionTruth, MatchResult mrU,
			Map<Locus, AuditRow> callRows, Map<Locus, AuditRow> truthRows, GradeResult r){
		for(int ti=0; ti<unionTruth.size(); ti++){
			final Locus t=unionTruth.get(ti);
			final AuditRow tr=truthRows.get(t);
			if(mrU.truthConsumed[ti]){
				final Locus c=posCalls.get(mrU.truthMatchedCallIdx[ti]);
				tr.unionOutcome="TP"; tr.unionMatchId=callRows.get(c).locusId;
			}else{
				tr.unionOutcome="FN"; r.unionFN++;
			}
		}
		for(int ci=0; ci<posCalls.size(); ci++){
			final Locus c=posCalls.get(ci);
			final AuditRow cr=callRows.get(c);
			if(mrU.callConsumed[ci]){
				final Locus t=unionTruth.get(mrU.callMatchedTruthIdx[ci]);
				cr.unionOutcome="TP"; cr.unionMatchId=truthRows.get(t).locusId;
			}else{
				cr.unionOutcome="FP"; r.unionFP++;
			}
		}
	}

	private static ArrayList<Locus> filter(List<Locus> list, String family){
		final ArrayList<Locus> out=new ArrayList<>();
		for(Locus l : list){if(family.equals(l.family)){out.add(l);}}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------              I/O              ----------------*/
	/*--------------------------------------------------------------*/

	/** Strips a single trailing '\r' if present. Small, deliberate duplication of
	 * NcrnaGffFamilyGrader's private helper of the same name (that one is `private`, not
	 * package-private, so it is not reusable across classes even within this package) --
	 * see the class javadoc: only the matching ENGINE is reused, trivial private
	 * utilities are not accessible to duplicate around. */
	private static String stripCR(String s){
		final int n=s.length();
		return (n>0 && s.charAt(n-1)=='\r') ? s.substring(0, n-1) : s;
	}

	private static char parseStrand(String s, String context){
		if(s.length()!=1){
			throw new IllegalArgumentException("Error - strand field must be exactly one character "
				+"('+', '-', or '.'), got '"+s+"': "+context);
		}
		final char c=s.charAt(0);
		if(c!='+' && c!='-' && c!='.'){
			throw new IllegalArgumentException("Error - invalid strand character '"+c+"': "+context);
		}
		return c;
	}

	private static int parseCoord(String s, String fieldName, String context){
		final int v;
		try{
			v=Integer.parseInt(s);
		}catch(NumberFormatException e){
			throw new IllegalArgumentException("Error - unparseable "+fieldName+" '"+s+"': "+context);
		}
		if(v<1){
			throw new IllegalArgumentException("Error - "+fieldName+" must be positive (1-based), got "+v+": "+context);
		}
		return v;
	}

	/** Tolerates leading '#'-prefixed comment lines and blank lines before the header
	 * row; the header itself is the first non-comment, non-blank line, and is skipped.
	 * Accepts family in {5S, 16S, 23S, 18S} (18S allowed here so a truth manifest may
	 * legitimately carry any scored legacy rRNA family).
	 * Rejects an exact-duplicate truth row (identical seqid+family+start+stop+strand). */
	static ArrayList<Locus> loadTruthTsv(String path) throws Exception {
		final ArrayList<Locus> out=new ArrayList<>();
		final HashSet<String> seenKeys=new HashSet<>();
		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, true, true));
		byte[] lineBytes=bf.nextLine();
		boolean sawHeader=false;
		while(lineBytes!=null){
			final String raw=stripCR(new String(lineBytes, StandardCharsets.US_ASCII));
			if(raw.isEmpty() || raw.charAt(0)=='#'){lineBytes=bf.nextLine(); continue;}
			if(!sawHeader){sawHeader=true; lineBytes=bf.nextLine(); continue;}//header row
			final String[] f=raw.split("\t", -1);
			if(f.length!=5){
				throw new IllegalArgumentException("Error - malformed truth row (need EXACTLY 5 "
					+"tab-delimited fields: seqid,family,start,stop,strand; got "+f.length+"): "+raw);
			}
			final String seqid=f[0], family=f[1];
			if(seqid.isEmpty()){throw new IllegalArgumentException("Error - empty seqid in truth row: "+raw);}
			if(!FAM_5S.equals(family) && !FAM_16S.equals(family) && !FAM_23S.equals(family) && !FAM_18S.equals(family)){
				throw new IllegalArgumentException("Error - unknown truth family '"+family
					+"' (must be 5S, 16S, 23S, or 18S): "+raw);
			}
			final int start=parseCoord(f[2], "start", raw), stop=parseCoord(f[3], "stop", raw);
			final char strand=parseStrand(f[4], raw);
			final String key=seqid+"\t"+family+"\t"+start+"\t"+stop+"\t"+strand;
			if(!seenKeys.add(key)){
				throw new IllegalArgumentException("Error - exact-duplicate truth row (same "
					+"seqid+family+start+stop+strand appears twice): "+raw);
			}
			out.add(new Locus(seqid, family, start, stop, strand));
			lineBytes=bf.nextLine();
		}
		bf.close();
		return out;
	}

	/** Parses a real callgenes.sh legacy-rRNA GFF. Only type=="rRNA" rows are inspected
	 * (NOT the generic-ncRNA "RNA" type NcrnaGffFamilyGrader reads). The subtype token is
	 * resolved via the exact substring convention gff.GffLine.prokType() already uses --
	 * see extractLegacyFamily. An rRNA-typed row matching none of {16S,23S,18S,5S} is a
	 * genuine attribution gap and fails loud: unlike NcrnaGffFamilyGrader's generic-ncRNA
	 * silent-skip (safe there because an unmatched RNA-typed row legitimately belongs to
	 * a different tracked-elsewhere family), there is no such legitimate other family for
	 * an rRNA-typed row in this grader's scope. */
	static ArrayList<Locus> loadCallsFromGff(String path) throws Exception {
		final ArrayList<Locus> out=new ArrayList<>();
		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, true, true));
		byte[] lineBytes=bf.nextLine();
		while(lineBytes!=null){
			if(lineBytes.length==0 || lineBytes[0]=='#'){lineBytes=bf.nextLine(); continue;}
			final String line=stripCR(new String(lineBytes, StandardCharsets.US_ASCII));
			if(line.isEmpty()){lineBytes=bf.nextLine(); continue;}
			final String[] f=line.split("\t", -1);
			if(f.length!=9){
				throw new IllegalArgumentException("Error - malformed GFF row (need EXACTLY 9 "
					+"tab-delimited fields per GFF3; got "+f.length+"): "+line);
			}
			final String seqid=f[0], type=f[2];
			if(seqid.isEmpty()){throw new IllegalArgumentException("Error - empty seqid in GFF row: "+line);}
			if(type.equals("rRNA")){
				final String attrs=f[8];
				final int start=parseCoord(f[3], "start", line), stop=parseCoord(f[4], "stop", line);
				final int len=stop-start+1;
				final String family=extractLegacyFamily(attrs, len);
				if(family==null){
					throw new IllegalStateException("Error - rRNA-typed GFF row with no recognized "
						+"subtype token (16S/23S/18S/5S) -- attribution gap, refusing to silently "
						+"ignore. Line: "+line);
				}
				final char strand=parseStrand(f[6], line);
				out.add(new Locus(seqid, family, start, stop, strand));
			}
			lineBytes=bf.nextLine();
		}
		bf.close();
		return out;
	}

	/** Mirrors gff.GffLine.prokType()'s exact substring convention (verified against
	 * that source this session), including its check ORDER (16S, 23S, 18S, then 5S gated
	 * by length) and its use of plain .contains(...) rather than isolating the
	 * attribute's first token specifically -- matching that method's real behavior, not a
	 * stricter reading of its intent. The 5S length gate (MAX_5S_LEN=300) is carried
	 * forward unchanged, not re-derived. */
	static String extractLegacyFamily(String attrs, int length){
		if(attrs.contains("16S")){return FAM_16S;}
		else if(attrs.contains("23S")){return FAM_23S;}
		else if(attrs.contains("18S")){return FAM_18S;}
		else if(attrs.contains("5S") && length<MAX_5S_LEN){return FAM_5S;}
		return null;
	}

	/** Writes the .rawcalls.tsv/.summary.tsv split with the same reconstructability
	 * guarantee as NcrnaGffFamilyGrader.writeReport. */
	static void writeReport(String out, GradeResult r) throws Exception {
		final PrintStream rawOut=new PrintStream(out+".rawcalls.tsv");
		rawOut.println("locusId\tseqid\tfamily\tstart\tstop\tstrand\trole\tsubtypeOutcome\tsubtypeMatchId"
			+"\tunionOutcome\tunionMatchId\tcrossLabelOutcome\tcrossLabelMatchId");
		for(AuditRow ar : r.auditRows){
			final String crossIds=(ar.crossLabelMatchIds.isEmpty() ? NA : String.join(",", ar.crossLabelMatchIds));
			rawOut.println(ar.locusId+"\t"+ar.locus.seqid+"\t"+ar.locus.family+"\t"+ar.locus.start+"\t"
				+ar.locus.stop+"\t"+ar.locus.strand+"\t"+ar.role+"\t"+ar.subtypeOutcome+"\t"+ar.subtypeMatchId
				+"\t"+ar.unionOutcome+"\t"+ar.unionMatchId+"\t"+ar.crossLabelOutcome+"\t"+crossIds);
		}
		rawOut.close();

		final PrintStream sumOut=new PrintStream(out+".summary.tsv");
		sumOut.println("scope\tTP\tFN\tFP\tprecision\trecall");
		for(String f : IN_SCOPE_FAMILIES){
			final FamilyCounts fc=r.byFamily.get(f);
			writeGradeRow(sumOut, f, fc.tp, fc.fn, fc.fp);
		}
		writeGradeRow(sumOut, "union", r.unionTP, r.unionFN, r.unionFP);
		sumOut.println();
		for(String f : IN_SCOPE_FAMILIES){
			sumOut.println("duplicateQueryLoci_"+f+"\t"+r.byFamily.get(f).dup);
		}
		sumOut.println("negSeqidFP\t"+r.negFP);
		sumOut.println();
		sumOut.println("# cross-label diagnostic ONLY -- never folded into union or subtype TP/FN/FP above");
		for(String f : IN_SCOPE_FAMILIES){
			sumOut.println("crossLabel_"+f+"truth_calledAsOther_count\t"+r.crossLabelCounts.get(f));
		}
		sumOut.close();

		for(String f : IN_SCOPE_FAMILIES){
			final FamilyCounts fc=r.byFamily.get(f);
			System.err.println(f+": TP="+fc.tp+" FN="+fc.fn+" FP="+fc.fp+" dup="+fc.dup);
		}
		System.err.println("union (pooled 5S/16S/18S/23S truth): TP="+r.unionTP+" FN="+r.unionFN+" FP="+r.unionFP);
		System.err.println("negSeqidFP="+r.negFP);
	}

	private static void writeGradeRow(PrintStream out, String scope, long tp, long fn, long fp){
		final double precision=(tp+fp==0 ? Double.NaN : (double)tp/(tp+fp));
		final double recall=(tp+fn==0 ? Double.NaN : (double)tp/(tp+fn));
		out.println(scope+"\t"+tp+"\t"+fn+"\t"+fp+"\t"
			+String.format(Locale.ROOT, "%.6f", precision)+"\t"+String.format(Locale.ROOT, "%.6f", recall));
	}
}
