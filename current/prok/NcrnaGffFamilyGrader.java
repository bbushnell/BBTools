package prok;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import fileIO.ByteFile;
import fileIO.FileFormat;

/** Coordinate-based, family-aware, one-to-one 6S (RF00013/RF01685) GFF grader.
 *
 * <p>Design (Dori manifest, CALIBRATION_TREE_MANIFEST_20260901.md, amendments
 * 2026-09-02 -- "concise implementation design" and two correction rounds, see below):
 * grades a real {@code callgenes.sh}-produced GFF against a coordinate-level truth
 * manifest, filling the gap that {@code NcrnaCombinedGradingDriver.gradeRecords}
 * cannot -- that method grades one pre-flanked single-locus record at a time (exactly
 * one truth locus per record by construction) and cannot express a real genome
 * carrying MULTIPLE true loci on one contig (e.g. the portable-fixture Bacillus tandem
 * RF00013 pair). This class does genuine per-seqid, per-family, coordinate-based
 * one-to-one matching instead.
 *
 * <p>Deliberately does NOT reuse {@code gff.CompareGff.overlapMetrics}: that method is
 * documented (this project's {@code ncrna} skill) as accepting the first qualifying
 * reference per query WITHOUT consuming it, so duplicate predictions can silently
 * multiply TPs -- exactly the defect {@link #matchOneToOne} is built to avoid.
 *
 * <p><b>Matching is MAXIMUM-CARDINALITY bipartite matching (Kuhn's augmenting-path
 * algorithm), not greedy highest-overlap-first</b> (Citan's review correction,
 * round 1): a first draft used greedy weighted assignment, which can under-count true
 * matches. Counterexample (now {@link NcrnaGffFamilyGraderTest#
 * testMaximumCardinalityCounterexample}, coordinates hand-verified before writing the
 * fixture): call C1 has edges to T1(overlap8) and T2(overlap7); call C2 has an edge
 * ONLY to T1(overlap5). Greedy assigns the single highest-overlap pair first (C1-T1),
 * stranding C2 -- 1 total match, even though C1-T2+C2-T1 gives 2. Kuhn's algorithm
 * finds the true maximum (2) by re-routing C1 onto T2 when C2 needs T1. Edges out of
 * each call are visited overlap-descending (comparator-safe {@code Integer.compare},
 * not overflow-prone subtraction) so that AMONG multiple maximum-cardinality solutions
 * the search is biased toward higher total overlap -- a practical, deterministic
 * heuristic, not a proof of globally overlap-maximal selection (a full weighted
 * assignment/Hungarian algorithm was not implemented; instance sizes here -- loci per
 * seqid+family -- are small enough that this is not expected to matter in practice,
 * per Citan's "if practical" qualifier).
 *
 * <p>Reuses {@code NcrnaCombinedGradingDriver}'s OUTPUT PHILOSOPHY (union = pooled
 * biological 6S truth; a wrong-subtype call is a real union TP, simultaneously a
 * subtype FP+FN; fail-loud on any attribution gap) AND its {@code .rawcalls.tsv}/
 * {@code .summary.tsv} file split.
 *
 * <p><b>Round-2 correction (auditability):</b> {@code .rawcalls.tsv} now carries, per
 * locus, its OWN subtype outcome, its OWN pooled-union outcome, its cross-label
 * diagnostic status, AND stable call&harr;truth match IDs ({@code locusId} per row,
 * {@code *MatchId} columns linking a matched pair) -- every aggregate number in
 * {@code .summary.tsv} is mechanically reconstructable by scanning {@code .rawcalls.tsv}
 * alone (grouping by role/family/outcome), not just human-eyeballable. Verified
 * end-to-end, not just claimed: {@link NcrnaGffFamilyGraderTest#
 * testEndToEndReconstruction} writes real {@code .rawcalls.tsv}/{@code .summary.tsv}
 * files, re-reads them from disk, and reconstructs every subtype/union/neg/dup/cross
 * aggregate purely from the parsed file content. {@code subtypeMatchId}/
 * {@code unionMatchId} are single-valued (a subtype/union match is genuinely
 * one-to-one, by construction of {@link #matchOneToOne}); {@code crossLabelMatchId} is
 * NOT -- an unmatched locus can legitimately cross-relate to more than one
 * opposite-family locus, so it is a comma-joined LIST of locusIds (round-2b
 * correction, Citan's review: an earlier draft used a single-valued field there and
 * silently lost all but the last link on a many-to-one relationship). See
 * {@link #writeReport} for the exact column list and {@link AuditRow}.
 *
 * <p>GFF attribute parsing note: BBTools' generic-ncRNA GFF attribute column is NOT
 * standard GFF3 {@code key=value;key2=value2} syntax -- it is comma-joined free text
 * (e.g. {@code RNA,startScr:0,stopScr:0,innerScr:0,len:186,model:6S_RF00013_unknown_c1113 n=5}).
 * {@link gff.GffLine} leaves this field as a raw String for exactly this reason. This
 * class extracts family by a prefix match on {@code model:6S_RF00013}/
 * {@code model:6S_RF01685} within that raw text, WITH a token-boundary check (round-2
 * correction): the character immediately after the prefix must be {@code _} or a
 * field/token boundary ({@code ,} or end-of-string) -- {@code 6S_RF000135_x} must NOT
 * be misattributed to RF00013 merely because it shares a numeric prefix.
 *
 * <p>Input validation (round-2 hardening, on top of round-1's CRLF/strand/field-count
 * work): truth TSV tolerates leading {@code #}-prefixed comment lines before its
 * header; every coordinate must be positive ({@code >=1}, 1-based); an exact-duplicate
 * truth row (identical seqid+family+start+stop+strand) fails loud, not silently
 * accepted twice.
 *
 * <p>Usage: {@code gff=<calls.gff> truth=<truth.tsv> [negseqids=<file>] out=<prefix>}.
 * {@code truth.tsv}: optional leading {@code #} comment lines, then a header row, then
 * EXACTLY {@code seqid  family  start  stop  strand} (1-based inclusive, family is
 * exactly {@code RF00013} or {@code RF01685}). {@code negseqids}: one seqid per line
 * (no internal whitespace), known to carry zero true 6S loci.
 *
 * @author G11
 */
public class NcrnaGffFamilyGrader {

	public static void main(String[] args) throws Exception {
		String gffPath=null, truthPath=null, negPath=null, out=null, familyMapPath=null;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<0){throw new IllegalArgumentException("Bad argument (expected flag=value): "+arg);}
			final String a=arg.substring(0, eq).toLowerCase(), b=arg.substring(eq+1);
			if(a.equals("gff")){gffPath=b;}
			else if(a.equals("truth")){truthPath=b;}
			else if(a.equals("negseqids")){negPath=b;}
			else if(a.equals("out")){out=b;}
			else if(a.equals("familymap")){familyMapPath=b;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(gffPath==null || truthPath==null || out==null){
			throw new IllegalArgumentException("Usage: gff=<calls.gff> truth=<truth.tsv> [negseqids=<file>] out=<prefix>");
		}
		final HashSet<String> negSeqids=(negPath==null ? new HashSet<String>() : loadSeqidSet(negPath));
		if(familyMapPath==null){
			final ArrayList<Locus> truth=loadTruthTsv(truthPath);
			final ArrayList<Locus> calls=loadCallsFromGff(gffPath);
			final GradeResult r=grade(truth, calls, negSeqids);
			writeReport(out, r);
		}else{
			final FamilyMap families=loadFamilyMap(familyMapPath);
			final ArrayList<Locus> truth=loadTruthTsv(truthPath, families);
			final ArrayList<Locus> calls=loadCallsFromGff(gffPath, families);
			writeGenericReport(out, gradeGeneric(truth, calls, negSeqids, families));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Core Data Types        ----------------*/
	/*--------------------------------------------------------------*/

	static final String FAM_A="RF00013", FAM_B="RF01685";
	static final String MODEL_PREFIX_A="6S_RF00013", MODEL_PREFIX_B="6S_RF01685";
	static final float MIN_OVERLAP_FRAC=0.5f;
	static final String NA="."; //no-match / not-applicable marker in audit output

	/** One genomic feature -- either a truth locus or an extracted query call. 1-based
	 * inclusive [start,stop], matching GFF convention throughout this codebase. Input
	 * validation, not an internal invariant -- IllegalArgumentException, not assert, so
	 * malformed data fails loud regardless of -ea/-da. */
	static final class Locus {
		final String seqid, family;
		final int start, stop;
		final char strand;
		Locus(String seqid_, String family_, int start_, int stop_, char strand_){
			if(start_<1){
				throw new IllegalArgumentException("Invalid locus (start<1, coordinates must be "
					+"positive 1-based): "+seqid_+" "+start_+"-"+stop_);
			}
			if(stop_<start_){
				throw new IllegalArgumentException("Invalid locus (stop<start): "+seqid_+" "+start_+"-"+stop_);
			}
			if(strand_!='+' && strand_!='-' && strand_!='.'){
				throw new IllegalArgumentException("Invalid strand character '"+strand_+"' for "+seqid_+" "+start_+"-"+stop_);
			}
			seqid=seqid_; family=family_; start=start_; stop=stop_; strand=strand_;
		}
		int length(){return stop-start+1;}
	}

	/*--------------------------------------------------------------*/
	/*----------------      One-to-One Matching      ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum-cardinality one-to-one bipartite match via Kuhn's augmenting-path
	 * algorithm -- see the class javadoc for why this replaced greedy weighted
	 * assignment and the exact counterexample. Deterministic: calls are processed in
	 * their input (load) order, and each call's candidate truth edges are visited
	 * overlap-descending via a comparator-safe {@code Integer.compare} -- never
	 * subtraction. Exposes the actual pairing ({@link MatchResult#callMatchedTruthIdx}/
	 * {@link MatchResult#truthMatchedCallIdx}), not just consumed/duplicate flags, so
	 * callers can build stable match-ID links for the audit trail. */
	static MatchResult matchOneToOne(List<Locus> calls, List<Locus> truth, float minOverlapFrac){
		final int nc=calls.size(), nt=truth.size();
		final MatchResult mr=new MatchResult(nc, nt);

		@SuppressWarnings("unchecked")
		final ArrayList<int[]>[] adj=new ArrayList[nc];
		final boolean[] callHasAnyEdge=new boolean[nc];
		for(int ci=0; ci<nc; ci++){
			final Locus c=calls.get(ci);
			final ArrayList<int[]> edges=new ArrayList<>();
			for(int ti=0; ti<nt; ti++){
				final Locus t=truth.get(ti);
				if(!c.seqid.equals(t.seqid) || c.strand!=t.strand){continue;}
				final int overlap=Math.min(c.stop, t.stop)-Math.max(c.start, t.start)+1;
				if(overlap<=0){continue;}
				if(overlap>=minOverlapFrac*c.length() && overlap>=minOverlapFrac*t.length()){
					edges.add(new int[]{ti, overlap});
					callHasAnyEdge[ci]=true;
				}
			}
			edges.sort((x, y)->Integer.compare(y[1], x[1]));//overlap descending, comparator-safe
			adj[ci]=edges;
		}

		for(int ci=0; ci<nc; ci++){
			final boolean[] visited=new boolean[nt];
			if(tryAugment(ci, adj, visited, mr.truthMatchedCallIdx, mr.callMatchedTruthIdx)){mr.matched++;}
		}

		for(int ci=0; ci<nc; ci++){
			if(mr.callMatchedTruthIdx[ci]>=0){mr.callConsumed[ci]=true;}
			else if(callHasAnyEdge[ci]){mr.callIsDuplicate[ci]=true;}
		}
		for(int ti=0; ti<nt; ti++){
			if(mr.truthMatchedCallIdx[ti]>=0){mr.truthConsumed[ti]=true;}
		}
		return mr;
	}

	/** One augmenting-path step of Kuhn's algorithm: try to give call `ci` a match,
	 * re-routing an already-matched truth locus's call onto an alternate edge if that
	 * frees up a valid assignment for `ci`. Returns true iff an augmenting path was
	 * found and truthMatchedTo/callMatchedTo were updated along it. */
	private static boolean tryAugment(int ci, ArrayList<int[]>[] adj, boolean[] visited,
			int[] truthMatchedTo, int[] callMatchedTo){
		for(int[] edge : adj[ci]){
			final int ti=edge[0];
			if(visited[ti]){continue;}
			visited[ti]=true;
			if(truthMatchedTo[ti]<0 || tryAugment(truthMatchedTo[ti], adj, visited, truthMatchedTo, callMatchedTo)){
				truthMatchedTo[ti]=ci;
				callMatchedTo[ci]=ti;
				return true;
			}
		}
		return false;
	}

	static final class MatchResult {
		int matched=0;
		final boolean[] callConsumed;
		final boolean[] truthConsumed;
		final boolean[] callIsDuplicate;
		/** callMatchedTruthIdx[ci] = index into the TRUTH list passed to matchOneToOne
		 * that call ci is matched to, or -1. truthMatchedCallIdx is the mirror. Exposed
		 * (not just consumed booleans) so callers can build stable audit-trail match IDs. */
		final int[] callMatchedTruthIdx;
		final int[] truthMatchedCallIdx;
		MatchResult(int nc, int nt){
			callConsumed=new boolean[nc];
			truthConsumed=new boolean[nt];
			callIsDuplicate=new boolean[nc];
			callMatchedTruthIdx=new int[nc]; Arrays.fill(callMatchedTruthIdx, -1);
			truthMatchedCallIdx=new int[nt]; Arrays.fill(truthMatchedCallIdx, -1);
		}
	}

	/** True if two loci share seqid+strand and pass reciprocal overlap -- the same
	 * predicate matchOneToOne uses internally, exposed for the cross-label diagnostic
	 * below (which needs a plain existence check, not a consuming match). */
	static boolean overlapsReciprocal(Locus a, Locus b, float minOverlapFrac){
		if(!a.seqid.equals(b.seqid) || a.strand!=b.strand){return false;}
		final int overlap=Math.min(a.stop, b.stop)-Math.max(a.start, b.start)+1;
		return overlap>0 && overlap>=minOverlapFrac*a.length() && overlap>=minOverlapFrac*b.length();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Grading            ----------------*/
	/*--------------------------------------------------------------*/

	/** One row of the per-locus audit trail -- either a truth locus or a query call,
	 * carrying its own subtype outcome, pooled-union outcome, cross-label diagnostic
	 * status, and stable match IDs linking it to whatever it matched (round-2
	 * auditability correction). locusId is stable within one grade() call: "T0..Tn-1"
	 * for truth in input-list order, "C0..Cm-1" for calls in input-list order --
	 * assigned before any filtering/partitioning, so it never depends on family or
	 * neg-seqid membership. */
	static final class AuditRow {
		final String locusId;
		final Locus locus;
		final String role; //TRUTH or CALL
		String subtypeOutcome=NA, subtypeMatchId=NA;
		String unionOutcome=NA, unionMatchId=NA;
		/** Round-2b correction (Citan's review): a single locus can legitimately
		 * cross-relate to MORE than one locus of the other family/role (e.g. one
		 * unmatched truth locus geometrically overlapped by two different unmatched
		 * opposite-family calls, or vice versa) -- a single-valued field would silently
		 * overwrite earlier links. Accumulates ALL related locusIds, in the fixed
		 * deterministic iteration order they were found (same input list order every
		 * run, never re-sorted -- no extra sort needed for determinism since the
		 * source lists are themselves already stable/ordered). Formatted as a
		 * comma-joined list at write time ("." if empty). */
		String crossLabelOutcome="NONE";
		final ArrayList<String> crossLabelMatchIds=new ArrayList<>();
		AuditRow(String locusId_, Locus locus_, String role_){locusId=locusId_; locus=locus_; role=role_;}
	}

	/** Pure, file-I/O-free grading core -- directly fixture-testable, mirroring
	 * NcrnaCombinedGradingDriver.gradeRecords/diagnoseCompetition's own "pure so it can
	 * be fixture-tested" design principle.
	 *
	 * <p>Three independent matching passes, not one threaded pass: (1) family-A calls
	 * vs family-A truth, strictly same-family, gives subtype-A TP/FN/FP directly; (2)
	 * the symmetric family-B pass; (3) one family-BLIND pooled pass over ALL non-neg
	 * calls and ALL truth answers "was 6S detected here at all" for the authoritative
	 * union scope. Cross-label counts (crossAasB/crossBasA) are a lightweight POST-HOC
	 * geometric diagnostic over each pass's own leftover unmatched truth -- informational
	 * only, never load-bearing for the TP/FN/FP arithmetic above.
	 *
	 * <p>Every AuditRow is initialized up front (one per truth locus, one per call, in
	 * input order via a LinkedHashMap keyed by object identity -- default Object
	 * equals/hashCode is exactly identity semantics here, which is what's wanted: two
	 * loci with identical field values are still distinct list entries) and then filled
	 * in by each pass -- so a row NOT touched by some pass keeps its NA/NONE defaults,
	 * never a stale or fabricated value. */
	static GradeResult grade(List<Locus> truth, List<Locus> calls, Set<String> negSeqids){
		final GradeResult r=new GradeResult();

		for(Locus t : truth){
			if(negSeqids.contains(t.seqid)){
				throw new IllegalStateException("Error - truth locus present on a negseqids-listed "
					+"seqid: "+t.seqid+" -- data inconsistency, refusing to grade.");
			}
		}

		//Keyed by object identity (default Object equals/hashCode) -- two loci with
		//identical field values are still distinct list entries and get distinct rows.
		//A REPEATED object reference (the exact same Locus instance appearing twice in
		//`truth` or `calls`) would collide and silently lose a row -- Citan's review,
		//2026-09-02: this is loader-impossible (loadTruthTsv/loadCallsFromGff always
		//construct a fresh `new Locus(...)` per parsed row, never reuse a reference) and
		//is not defended against here; documented, not handled, per that review.
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

		final ArrayList<Locus> negCalls=new ArrayList<>();
		final ArrayList<Locus> posCalls=new ArrayList<>();
		for(Locus c : calls){
			if(negSeqids.contains(c.seqid)){negCalls.add(c);}else{posCalls.add(c);}
		}
		for(Locus c : negCalls){
			if(FAM_A.equals(c.family)){r.fpA++;}else if(FAM_B.equals(c.family)){r.fpB++;}
			r.unionFP++; r.negFP++;
			final AuditRow cr=callRows.get(c);
			cr.subtypeOutcome="NEG_FP"; cr.unionOutcome="FP";
		}

		final ArrayList<Locus> truthA=filter(truth, FAM_A), truthB=filter(truth, FAM_B);
		final ArrayList<Locus> callsA=filter(posCalls, FAM_A), callsB=filter(posCalls, FAM_B);

		final MatchResult mrA=matchOneToOne(callsA, truthA, MIN_OVERLAP_FRAC);
		final MatchResult mrB=matchOneToOne(callsB, truthB, MIN_OVERLAP_FRAC);
		r.tpA=mrA.matched; r.tpB=mrB.matched;
		applySubtypePass(callsA, truthA, mrA, callRows, truthRows, r, true);
		applySubtypePass(callsB, truthB, mrB, callRows, truthRows, r, false);

		final MatchResult mrU=matchOneToOne(posCalls, truth, MIN_OVERLAP_FRAC);
		r.unionTP=mrU.matched;
		applyUnionPass(posCalls, truth, mrU, callRows, truthRows, r);

		//Round-2b correction: NO `break` on the first hit -- a single unmatched truth
		//locus can be geometrically overlapped by MULTIPLE unmatched opposite-family
		//calls (and a single call can likewise cross-relate to more than one truth
		//locus across separate iterations of this loop); every such relationship is
		//recorded via crossLabelMatchIds.add, never overwritten. r.crossAasB/crossBasA
		//still count ONCE per truth locus (the diagnostic question is "how many truth
		//loci got a cross-label rescue," not "how many pairwise relationships exist").
		for(int ti=0; ti<truthA.size(); ti++){
			if(mrA.truthConsumed[ti]){continue;}
			final Locus t=truthA.get(ti);
			final AuditRow tr=truthRows.get(t);
			for(Locus c : callsB){
				if(overlapsReciprocal(c, t, MIN_OVERLAP_FRAC)){
					tr.crossLabelOutcome="CROSS_DETECTED"; tr.crossLabelMatchIds.add(callRows.get(c).locusId);
					final AuditRow cr=callRows.get(c);
					cr.crossLabelOutcome="CROSS_MATCH"; cr.crossLabelMatchIds.add(tr.locusId);
				}
			}
			if(!tr.crossLabelMatchIds.isEmpty()){r.crossAasB++;}
		}
		for(int ti=0; ti<truthB.size(); ti++){
			if(mrB.truthConsumed[ti]){continue;}
			final Locus t=truthB.get(ti);
			final AuditRow tr=truthRows.get(t);
			for(Locus c : callsA){
				if(overlapsReciprocal(c, t, MIN_OVERLAP_FRAC)){
					tr.crossLabelOutcome="CROSS_DETECTED"; tr.crossLabelMatchIds.add(callRows.get(c).locusId);
					final AuditRow cr=callRows.get(c);
					cr.crossLabelOutcome="CROSS_MATCH"; cr.crossLabelMatchIds.add(tr.locusId);
				}
			}
			if(!tr.crossLabelMatchIds.isEmpty()){r.crossBasA++;}
		}

		r.auditRows.addAll(truthRows.values());
		r.auditRows.addAll(callRows.values());
		return r;
	}

	/** Fills subtypeOutcome/subtypeMatchId for one family's pass, and increments the
	 * matching FN/FP/dup counters on r. isFamilyA selects which of r's fnA/fpA/dupA vs
	 * fnB/fpB/dupB fields to increment. */
	private static void applySubtypePass(List<Locus> callsFam, List<Locus> truthFam, MatchResult mr,
			Map<Locus, AuditRow> callRows, Map<Locus, AuditRow> truthRows, GradeResult r, boolean isFamilyA){
		for(int ti=0; ti<truthFam.size(); ti++){
			final Locus t=truthFam.get(ti);
			final AuditRow tr=truthRows.get(t);
			if(mr.truthConsumed[ti]){
				final Locus c=callsFam.get(mr.truthMatchedCallIdx[ti]);
				tr.subtypeOutcome="TP"; tr.subtypeMatchId=callRows.get(c).locusId;
			}else{
				tr.subtypeOutcome="FN";
				if(isFamilyA){r.fnA++;}else{r.fnB++;}
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
				if(isFamilyA){r.fpA++; if(dup){r.dupA++;}}else{r.fpB++; if(dup){r.dupB++;}}
			}
		}
	}

	/** Fills unionOutcome/unionMatchId for the family-blind pooled pass, and increments
	 * r.unionFN/unionFP. Neg-seqid calls are handled separately in grade() (they never
	 * enter posCalls/this pass) since they are unmatched-by-construction. */
	private static void applyUnionPass(List<Locus> posCalls, List<Locus> truth, MatchResult mrU,
			Map<Locus, AuditRow> callRows, Map<Locus, AuditRow> truthRows, GradeResult r){
		for(int ti=0; ti<truth.size(); ti++){
			final Locus t=truth.get(ti);
			final AuditRow tr=truthRows.get(t);
			if(mrU.truthConsumed[ti]){
				final Locus c=posCalls.get(mrU.truthMatchedCallIdx[ti]);
				tr.unionOutcome="TP"; tr.unionMatchId=callRows.get(c).locusId;
			}else{
				tr.unionOutcome="FN";
				r.unionFN++;
			}
		}
		for(int ci=0; ci<posCalls.size(); ci++){
			final Locus c=posCalls.get(ci);
			final AuditRow cr=callRows.get(c);
			if(mrU.callConsumed[ci]){
				final Locus t=truth.get(mrU.callMatchedTruthIdx[ci]);
				cr.unionOutcome="TP"; cr.unionMatchId=truthRows.get(t).locusId;
			}else{
				cr.unionOutcome="FP";
				r.unionFP++;
			}
		}
	}

	private static ArrayList<Locus> filter(List<Locus> list, String family){
		final ArrayList<Locus> out=new ArrayList<>();
		for(Locus l : list){if(family.equals(l.family)){out.add(l);}}
		return out;
	}

	/** Aggregate grading result. union = pooled biological 6S truth (authoritative);
	 * crossAasB/crossBasA are an explicitly-diagnostic-only cross-label count, never
	 * folded into union or subtype TP/FN/FP -- mirrors
	 * NcrnaCombinedGradingDriver.GradeResult's union-vs-diagnostic separation exactly.
	 * auditRows (truth rows then call rows, each input-order-stable) backs the
	 * per-locus .rawcalls.tsv output and is sufficient on its own to reconstruct every
	 * field below by scanning/grouping -- see the class javadoc. */
	static final class GradeResult {
		long tpA, fpA, fnA, dupA;
		long tpB, fpB, fnB, dupB;
		long unionTP, unionFP, unionFN;
		long crossAasB, crossBasA;
		long negFP;
		final ArrayList<AuditRow> auditRows=new ArrayList<>();
	}

	/*--------------------------------------------------------------*/
	/*----------------              I/O              ----------------*/
	/*--------------------------------------------------------------*/

	/** Strips a single trailing '\r' if present -- CRLF line endings are a convention
	 * difference, not malformed data, and are tolerated everywhere a line is read. */
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
	 * row (round-2 hardening); the header itself is the first non-comment, non-blank
	 * line, and is skipped. Rejects an exact-duplicate truth row (identical
	 * seqid+family+start+stop+strand) -- virtually certainly a data/extraction bug,
	 * not a genuine feature of biology; fails loud rather than silently double-counting. */
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
			if(!family.equals(FAM_A) && !family.equals(FAM_B)){
				throw new IllegalArgumentException("Error - unknown truth family '"+family
					+"' (must be "+FAM_A+" or "+FAM_B+"): "+raw);
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

	static HashSet<String> loadSeqidSet(String path) throws Exception {
		final HashSet<String> out=new HashSet<>();
		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, true, true));
		byte[] lineBytes=bf.nextLine();
		while(lineBytes!=null){
			final String line=stripCR(new String(lineBytes, StandardCharsets.US_ASCII));
			if(!line.isEmpty()){
				if(line.indexOf('\t')>=0 || line.indexOf(' ')>=0){
					throw new IllegalArgumentException("Error - negseqids line contains internal "
						+"whitespace (expected one bare seqid per line): '"+line+"'");
				}
				out.add(line);
			}
			lineBytes=bf.nextLine();
		}
		bf.close();
		return out;
	}

	/** Parses a real callgenes.sh GFF. Only type=="RNA" rows are inspected (the shared
	 * generic-ncRNA type -- see current-6s.md). A row with no "model:" token at all is a
	 * genuine attribution gap and fails loud (Orf.java always sets trnaModel for a
	 * generic-ncRNA hit, so its absence means something is structurally wrong, not that
	 * this is a different family). A row WITH a model: token whose prefix is neither
	 * sixS family (with a valid token boundary -- see extractFamily) is silently
	 * skipped -- it belongs to a different generic-ncRNA family (rnasep/srp_small/
	 * srp_large/tmrna) this grader does not track, not an error. */
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
			if(type.equals("RNA")){
				final String attrs=f[8];
				final int mi=attrs.indexOf("model:");
				if(mi<0){
					throw new IllegalStateException("Error - RNA-typed GFF row with no model: attribute "
						+"-- attribution gap, refusing to silently ignore. Line: "+line);
				}
				final String family=extractFamily(attrs, mi);
				if(family!=null){
					final int start=parseCoord(f[3], "start", line), stop=parseCoord(f[4], "stop", line);
					final char strand=parseStrand(f[6], line);
					out.add(new Locus(seqid, family, start, stop, strand));
				}
				//else: a different (non-sixS) generic ncRNA family's model, or a
				//boundary-violating near-miss -- not this grader's concern, skip silently.
			}
			lineBytes=bf.nextLine();
		}
		bf.close();
		return out;
	}

	/** Round-2 hardening: a bare prefix match is not enough -- "6S_RF000135_x" must NOT
	 * be misattributed to RF00013 merely because it shares the same leading digits.
	 * After the prefix matches, the very next character (if the string continues at
	 * all) must be a real token boundary: '_' (the model-name's own internal separator,
	 * e.g. "6S_RF00013_unknown_c1113") or ',' (end of this comma-joined attribute
	 * field). Anything else (a digit, a letter) means the prefix match was coincidental,
	 * not a real family identifier -- rejected. */
	static String extractFamily(String attrs, int modelIdx){
		final int valStart=modelIdx+"model:".length();
		final String fam=matchPrefixWithBoundary(attrs, valStart, MODEL_PREFIX_A) ? FAM_A
			: matchPrefixWithBoundary(attrs, valStart, MODEL_PREFIX_B) ? FAM_B : null;
		return fam;
	}

	private static boolean matchPrefixWithBoundary(String attrs, int valStart, String prefix){
		if(!attrs.regionMatches(valStart, prefix, 0, prefix.length())){return false;}
		final int afterPrefix=valStart+prefix.length();
		if(afterPrefix>=attrs.length()){return true;}//prefix runs to end of string -- valid boundary
		final char next=attrs.charAt(afterPrefix);
		return next=='_' || next==',';
	}

	/** Writes the {@code .rawcalls.tsv}/{@code .summary.tsv} split. rawcalls is the
	 * FULL per-locus audit trail (round-2 auditability correction): one row per truth
	 * locus AND one row per query call, each carrying its own subtype outcome, pooled-
	 * union outcome, cross-label diagnostic status, and stable match IDs -- every
	 * .summary.tsv aggregate is mechanically reconstructable from this file alone by
	 * grouping on role/family/outcome (e.g. tpA == count of role=="TRUTH" AND
	 * family=="RF00013" AND subtypeOutcome=="TP"; unionFP == count of role=="CALL" AND
	 * unionOutcome=="FP"; crossAasB == count of role=="TRUTH" AND family=="RF00013" AND
	 * crossLabelOutcome=="CROSS_DETECTED"). Unlike NcrnaCombinedGradingDriver's
	 * rawcalls.tsv (one row per input RECORD -- that unit doesn't exist here, since
	 * this grader works at genome/coordinate scale, not pre-flanked-single-locus-record
	 * scale), the natural unit is one row per LOCUS. */
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
		writeGradeRow(sumOut, FAM_A, r.tpA, r.fnA, r.fpA);
		writeGradeRow(sumOut, FAM_B, r.tpB, r.fnB, r.fpB);
		writeGradeRow(sumOut, "union", r.unionTP, r.unionFN, r.unionFP);
		sumOut.println();
		sumOut.println("duplicateQueryLoci_"+FAM_A+"\t"+r.dupA);
		sumOut.println("duplicateQueryLoci_"+FAM_B+"\t"+r.dupB);
		sumOut.println("negSeqidFP\t"+r.negFP);
		sumOut.println();
		sumOut.println("# cross-label diagnostic ONLY -- never folded into union or subtype TP/FN/FP above");
		sumOut.println("crossLabel_"+FAM_A+"truth_calledAs"+FAM_B+"_count\t"+r.crossAasB);
		sumOut.println("crossLabel_"+FAM_B+"truth_calledAs"+FAM_A+"_count\t"+r.crossBasA);
		sumOut.close();

		System.err.println(FAM_A+": TP="+r.tpA+" FN="+r.fnA+" FP="+r.fpA+" dup="+r.dupA);
		System.err.println(FAM_B+": TP="+r.tpB+" FN="+r.fnB+" FP="+r.fpB+" dup="+r.dupB);
		System.err.println("union (pooled biological 6S truth): TP="+r.unionTP+" FN="+r.unionFN+" FP="+r.unionFP);
		System.err.println("negSeqidFP="+r.negFP);
		System.err.println("crossLabel "+FAM_A+"->"+FAM_B+": "+r.crossAasB+"  "+FAM_B+"->"+FAM_A+": "+r.crossBasA);
	}

	private static void writeGradeRow(PrintStream out, String scope, long tp, long fn, long fp){
		final double precision=(tp+fp==0 ? Double.NaN : (double)tp/(tp+fp));
		final double recall=(tp+fn==0 ? Double.NaN : (double)tp/(tp+fn));
		out.println(scope+"\t"+tp+"\t"+fn+"\t"+fp+"\t"
			+String.format(Locale.ROOT, "%.6f", precision)+"\t"+String.format(Locale.ROOT, "%.6f", recall));
	}

	/*--------------------------------------------------------------*/
	/*----------------   Generic family-map grading  ---------------*/
	/*--------------------------------------------------------------*/

	/** Ordered, validated family-to-model-prefix configuration for development
	 * benchmarks.  The legacy two-6S-family behavior above intentionally remains
	 * byte-for-byte independent when familymap is omitted. */
	static final class FamilyMap {
		final LinkedHashMap<String, String> prefixByFamily=new LinkedHashMap<>();
	}

	static FamilyMap loadFamilyMap(String path) throws Exception {
		final FamilyMap out=new FamilyMap();
		final HashSet<String> prefixes=new HashSet<>();
		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, true, true));
		for(byte[] lineBytes=bf.nextLine(); lineBytes!=null; lineBytes=bf.nextLine()){
			final String line=stripCR(new String(lineBytes, StandardCharsets.US_ASCII));
			if(line.isEmpty() || line.charAt(0)=='#'){continue;}
			final String[] f=line.split("\\t", -1);
			if(f.length!=2 || f[0].isEmpty() || f[1].isEmpty()){
				throw new IllegalArgumentException("Malformed familymap row (need family<TAB>gff_model_prefix): "+line);
			}
			if(out.prefixByFamily.put(f[0], f[1])!=null || !prefixes.add(f[1])){
				throw new IllegalArgumentException("Duplicate family or model prefix in familymap: "+line);
			}
		}
		bf.close();
		if(out.prefixByFamily.isEmpty()){throw new IllegalArgumentException("familymap is empty: "+path);}
		return out;
	}

	static ArrayList<Locus> loadTruthTsv(String path, FamilyMap families) throws Exception {
		final ArrayList<Locus> out=new ArrayList<>();
		final HashSet<String> seenKeys=new HashSet<>();
		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, true, true));
		boolean sawHeader=false;
		for(byte[] lineBytes=bf.nextLine(); lineBytes!=null; lineBytes=bf.nextLine()){
			final String raw=stripCR(new String(lineBytes, StandardCharsets.US_ASCII));
			if(raw.isEmpty() || raw.charAt(0)=='#'){continue;}
			if(!sawHeader){sawHeader=true; continue;}
			final String[] f=raw.split("\\t", -1);
			if(f.length!=5){throw new IllegalArgumentException("Malformed truth row (need 5 tab-delimited fields): "+raw);}
			if(!families.prefixByFamily.containsKey(f[1])){
				throw new IllegalArgumentException("Truth family is not configured in familymap: '"+f[1]+"': "+raw);
			}
			final int start=parseCoord(f[2], "start", raw), stop=parseCoord(f[3], "stop", raw);
			final char strand=parseStrand(f[4], raw);
			final String key=f[0]+"\t"+f[1]+"\t"+start+"\t"+stop+"\t"+strand;
			if(!seenKeys.add(key)){throw new IllegalArgumentException("Exact-duplicate truth row: "+raw);}
			out.add(new Locus(f[0], f[1], start, stop, strand));
		}
		bf.close();
		return out;
	}

	static ArrayList<Locus> loadCallsFromGff(String path, FamilyMap families) throws Exception {
		final ArrayList<Locus> out=new ArrayList<>();
		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, true, true));
		for(byte[] lineBytes=bf.nextLine(); lineBytes!=null; lineBytes=bf.nextLine()){
			if(lineBytes.length==0 || lineBytes[0]=='#'){continue;}
			final String line=stripCR(new String(lineBytes, StandardCharsets.US_ASCII));
			if(line.isEmpty()){continue;}
			final String[] f=line.split("\\t", -1);
			if(f.length!=9){throw new IllegalArgumentException("Malformed GFF row (need 9 tab-delimited fields): "+line);}
			if(!f[2].equals("RNA")){continue;}
			final int mi=f[8].indexOf("model:");
			if(mi<0){throw new IllegalStateException("RNA GFF row has no model: attribute: "+line);}
			final int valStart=mi+"model:".length();
			String family=null;
			for(Map.Entry<String,String> e : families.prefixByFamily.entrySet()){
				if(matchPrefixWithBoundary(f[8], valStart, e.getValue())){
					if(family!=null){throw new IllegalStateException("Ambiguous model-prefix attribution: "+line);}
					family=e.getKey();
				}
			}
			if(family!=null){out.add(new Locus(f[0], family, parseCoord(f[3], "start", line), parseCoord(f[4], "stop", line), parseStrand(f[6], line)));}
		}
		bf.close();
		return out;
	}

	static final class FamilyCounts { long tp, fp, fn, dup; }

	/** One matched TP pair with strand-aware boundary offsets (G11, 2026-09-05; ncrna
	 * skill: detection and boundary accuracy are separate reports). Offsets are in
	 * TRANSCRIPT orientation: d5 = call's 5' end minus truth's 5' end, d3 = call's 3'
	 * end minus truth's 3' end, each positive when the call end lies DOWNSTREAM of the
	 * truth end (toward the 3' direction). Plus strand: d5=call.start-truth.start,
	 * d3=call.stop-truth.stop; minus strand: d5=truth.stop-call.stop,
	 * d3=truth.start-call.start. */
	static final class MatchedPair {
		final String family; final Locus truth, call; final int d5, d3;
		MatchedPair(String family_, Locus t, Locus c){
			family=family_; truth=t; call=c;
			final boolean minus=(t.strand=='-');
			d5=(minus ? t.stop-c.stop : c.start-t.start);
			d3=(minus ? t.start-c.start : c.stop-t.stop);
		}
	}

	static final class GenericGradeResult {
		final FamilyMap families;
		final LinkedHashMap<String, FamilyCounts> counts=new LinkedHashMap<>();
		final ArrayList<MatchedPair> pairs=new ArrayList<>();//per-family TP pairs, load order
		long unionTP, unionFP, unionFN, negFP;
		GenericGradeResult(FamilyMap f){families=f; for(String family : f.prefixByFamily.keySet()){counts.put(family, new FamilyCounts());}}
	}

	static GenericGradeResult gradeGeneric(List<Locus> truth, List<Locus> calls, Set<String> negSeqids, FamilyMap families){
		final GenericGradeResult r=new GenericGradeResult(families);
		for(Locus t : truth){if(negSeqids.contains(t.seqid)){throw new IllegalStateException("Truth locus present on negseqids seqid: "+t.seqid);}}
		final ArrayList<Locus> posCalls=new ArrayList<>();
		for(Locus c : calls){
			if(negSeqids.contains(c.seqid)){r.counts.get(c.family).fp++; r.negFP++; r.unionFP++;}
			else{posCalls.add(c);}
		}
		for(String family : families.prefixByFamily.keySet()){
			final ArrayList<Locus> tf=filter(truth, family), cf=filter(posCalls, family);
			final MatchResult mr=matchOneToOne(cf, tf, MIN_OVERLAP_FRAC);
			final FamilyCounts c=r.counts.get(family);
			c.tp=mr.matched;
			for(int i=0;i<tf.size();i++){if(!mr.truthConsumed[i]){c.fn++;}}
			for(int i=0;i<cf.size();i++){if(!mr.callConsumed[i]){c.fp++; if(mr.callIsDuplicate[i]){c.dup++;}}}
			for(int ti=0; ti<tf.size(); ti++){
				final int ci=mr.truthMatchedCallIdx[ti];
				if(ci>=0){r.pairs.add(new MatchedPair(family, tf.get(ti), cf.get(ci)));}
			}
		}
		final MatchResult union=matchOneToOne(posCalls, truth, MIN_OVERLAP_FRAC);
		r.unionTP=union.matched;
		for(int i=0;i<truth.size();i++){if(!union.truthConsumed[i]){r.unionFN++;}}
		for(int i=0;i<posCalls.size();i++){if(!union.callConsumed[i]){r.unionFP++;}}
		return r;
	}

	static void writeGenericReport(String out, GenericGradeResult r) throws Exception {
		final PrintStream sumOut=new PrintStream(out+".summary.tsv");
		sumOut.println("scope\tTP\tFN\tFP\tprecision\trecall");
		for(Map.Entry<String, FamilyCounts> e : r.counts.entrySet()){
			final FamilyCounts c=e.getValue(); writeGradeRow(sumOut, e.getKey(), c.tp, c.fn, c.fp);
		}
		writeGradeRow(sumOut, "union", r.unionTP, r.unionFN, r.unionFP);
		sumOut.println();
		for(Map.Entry<String, FamilyCounts> e : r.counts.entrySet()){sumOut.println("duplicateQueryLoci_"+e.getKey()+"\t"+e.getValue().dup);}
		sumOut.println("negSeqidFP\t"+r.negFP);
		//Boundary-accuracy block (separate from detection, per the ncRNA methodology):
		//per family over matched TPs -- exact-end counts and signed/absolute offset stats.
		sumOut.println();
		sumOut.println("# boundary accuracy over matched TPs; d5/d3 in transcript orientation, positive=call end downstream of truth end");
		sumOut.println("family\tnTP\texact5\texact3\tbothExact\tmeanAbsD5\tmeanAbsD3\tmaxAbsD5\tmaxAbsD3");
		for(String family : r.counts.keySet()){
			long n=0, e5=0, e3=0, both=0, s5=0, s3=0, m5=0, m3=0;
			for(MatchedPair mp : r.pairs){
				if(!mp.family.equals(family)){continue;}
				n++;
				if(mp.d5==0){e5++;}
				if(mp.d3==0){e3++;}
				if(mp.d5==0 && mp.d3==0){both++;}
				final long a5=Math.abs(mp.d5), a3=Math.abs(mp.d3);
				s5+=a5; s3+=a3; m5=Math.max(m5, a5); m3=Math.max(m3, a3);
			}
			sumOut.println(family+"\t"+n+"\t"+e5+"\t"+e3+"\t"+both
				+"\t"+(n==0 ? "." : String.format("%.2f", s5/(double)n))
				+"\t"+(n==0 ? "." : String.format("%.2f", s3/(double)n))
				+"\t"+(n==0 ? "." : String.valueOf(m5))+"\t"+(n==0 ? "." : String.valueOf(m3)));
		}
		sumOut.close();
		final PrintStream mOut=new PrintStream(out+".matches.tsv");
		mOut.println("#family\tseqid\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3");
		for(MatchedPair mp : r.pairs){
			mOut.println(mp.family+"\t"+mp.truth.seqid+"\t"+mp.truth.strand+"\t"+mp.truth.start+"\t"+mp.truth.stop
				+"\t"+mp.call.start+"\t"+mp.call.stop+"\t"+mp.d5+"\t"+mp.d3);
		}
		mOut.close();
	}
}
