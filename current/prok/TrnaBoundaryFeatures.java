package prok;

import consensus.BaseGraph;
import idaligner.AlignmentStats;
import idaligner.ScrabbleAligner;
import stream.Read;

/**
 * Feature extraction for the tRNA boundary-precision NN scorer (Brian's idea #3,
 * trna repo plans/boundary_precision_ideas.md — a small NN that scores "is this a
 * correct tRNA start/stop" at a candidate boundary and its +-1/+-2bp neighbors,
 * generalizing the existing acceptor-stem snap and a planned 9-mer enrichment rule
 * into learned feature weights instead of hand-picked cutoffs).
 *
 * Assembled while Dori is down (Aug 19, 2026): two of the ~10 planned features
 * (9-mer enrichment table, model-tip fuzziness) need corpus-scale data that does
 * not exist offline. Their outputs are PENDING_DORI sentinels here — not
 * fabricated placeholder numbers — until that data is built on the cluster.
 * @author Neptune
 */
public class TrnaBoundaryFeatures {

	/** Sentinel for a feature that needs corpus-scale data not yet available. */
	public static final float PENDING_DORI=-1f;

	/**
	 * Acceptor-stem palindrome strength at a candidate [s,e) span, normalized to
	 * 0-1 (raw score is 0-14). Reuses TrnaCaller's existing, shipped
	 * acceptorStemScore — ready now, no corpus needed.
	 */
	public static float stemFeature(byte[] seq, int s, int e){
		if(s<0 || e<0 || e>=seq.length || e-s<20){return 0;}
		return TrnaCaller.acceptorStemScore(seq, s, e)/14f;
	}

	/** Mean mature tRNA length (bp), used only to normalize lengthRatioFeature. */
	private static final float MEAN_TRNA_LEN=76f;

	/**
	 * Candidate length ratio (Brian+Noire, Aug 21): this candidate's own length
	 * divided by the mean mature tRNA length. Computed directly from the
	 * candidate span [s,e] already available at every call site -- unlike
	 * contig-GC, this does NOT need a header field or the source corpus at all,
	 * since it varies per shifted candidate (a +-1/+-2bp shift changes e-s by
	 * exactly that much) and is already fully determined by the span itself.
	 */
	public static float lengthRatioFeature(int s, int e){
		return (e-s+1)/MEAN_TRNA_LEN;
	}

	/**
	 * Consensus-alignment identity (0-1) of a candidate span against a specific
	 * library model's consensus sequence. Ready now — same aligner call
	 * TrnaCaller.alignToModel uses (idaligner.ScrabbleAligner.alignAndTraceStatic).
	 *
	 * <p>LEAKAGE FIX (Brian, Aug 21 2026): this must be called with EACH
	 * candidate's own actual [s,e] span, never a span reused from elsewhere (e.g.
	 * the locus's true, unshifted boundary) — a +-1/+-2bp shift genuinely changes
	 * identity, and reusing one span's value across all shifted candidates leaks
	 * true-boundary information into every candidate's training vector. Both real
	 * callers (TrnaBoundaryVectorGen.emitBoundaryVectors,
	 * TrnaBoundaryScorer.score) recompute this per candidate; an earlier version
	 * of this file's own docstring described the OLD, leaky "per-locus constant,
	 * not recomputed per candidate" behavior -- that was the bug, not a design
	 * feature, and is fixed at both call sites now.
	 */
	public static float aniFeature(byte[] candidate, byte[] consensus){
		AlignmentStats stats=new AlignmentStats(true);
		stats.doTrace=true;
		if(candidate.length>=consensus.length){ScrabbleAligner.alignAndTraceStatic(consensus, candidate, stats);}
		else{ScrabbleAligner.alignAndTraceStatic(candidate, consensus, stats);}
		return stats.identity;
	}

	/**
	 * Injectable boundary-kmer enrichment table (idea #2: count of the kmer at
	 * a boundary, tabulated over a known-boundary corpus). Despite the name
	 * (kept for historical continuity with the k=9 design), the table's kmer
	 * length is whatever k it was BUILT with (see prok.TrnaNinemerTableBuilder's
	 * insidecount=/outsidecount= flags) — a table now carries its own k/type/
	 * insideCount/outsideCount (TrnaNinemerTableBuilder.LoadedTable), so a caller
	 * reads those FROM the table rather than separately specifying them (kills the
	 * silent-mismatch risk Noire flagged in review, Aug 21). Null until a real
	 * corpus-derived table is built. A caller MAY inject a locally-derived table
	 * (e.g. for a plumbing smoke test), but that is explicitly NOT the real
	 * enrichment signal and must be labeled as such by the caller.
	 */
	public interface NinemerTable {
		/** RAW observed count of the kmer at [pos,pos+k) in seq (0 if never observed,
		 * out of bounds, or containing a non-ACGT base), for whatever k the table was
		 * built with. Deliberately raw, not log-smoothed -- callers (enrichmentProfile)
		 * apply their own leave-one-out and smoothing, which needs the true integer
		 * count, not a pre-transformed value. */
		int rawCount(byte[] seq, int pos);
	}

	/**
	 * Which end of the tRNA a boundaryPos marks. NOT cosmetic -- inside/outside
	 * flip direction depending on this: for START, bases BEFORE boundaryPos are
	 * outside (upstream flank) and bases FROM boundaryPos onward are inside; for
	 * STOP, boundaryPos IS the last inside base, so bases UP TO AND INCLUDING it
	 * are inside and bases AFTER it are outside (downstream flank). A raw signed
	 * "offset" parameter that ignores this (the ORIGINAL kmerWindowStart(pos,k,out)
	 * design, removed Aug 21) silently flips its own meaning between the two
	 * boundary types -- caught empirically when a k-offset grid swept BOTH types
	 * with the same "out" and the STOP results read backwards from expectation
	 * until traced to this exact bug. kmerWindowStart below takes semantic
	 * insideCount/outsideCount instead, so a caller who says "9 inside, 0 outside"
	 * gets the same real window regardless of which boundary type it's for --
	 * the direction flip is computed HERE, once, and can't be gotten wrong at a
	 * call site.
	 */
	public enum BoundaryType { START, STOP }

	/**
	 * The single shared window-start formula for a boundary kmer lookup, used by
	 * BOTH the table builder (prok.TrnaNinemerTableBuilder) and this class's
	 * ninemerFeatures -- factored out (Noire's review, Aug 21) so the two sides
	 * cannot silently diverge. Hardened a second time the same day (Noire, after
	 * the inside/outside flip was caught) to take semantic insideCount/outsideCount
	 * + boundary type instead of a raw signed offset -- see BoundaryType's javadoc
	 * for why the raw-offset form was retired rather than merely renamed.
	 *
	 * <p>For START: window=[boundaryPos-outsideCount, boundaryPos+insideCount).
	 * For STOP: window=[boundaryPos-(insideCount-1), boundaryPos+outsideCount+1)
	 * -- boundaryPos itself is the last INSIDE base, so it's counted toward
	 * insideCount, not outsideCount.
	 *
	 * <p>Locked production values (discrimination grid on the 2.01M-genome dedup
	 * corpus, Aug 21): k=9, START insideCount=9/outsideCount=0 (fully inside),
	 * STOP insideCount=8/outsideCount=1 (one base peeking into the flank).
	 */
	public static int kmerWindowStart(int boundaryPos, BoundaryType type, int insideCount, int outsideCount){
		if(type==BoundaryType.START){
			return boundaryPos-outsideCount;
		}else{
			return boundaryPos-(insideCount-1);
		}
	}

	/** Sentinel meaning "no true-boundary position to leave-one-out" -- used by the
	 * INFERENCE overload of enrichmentProfile, where there is no ground truth (the
	 * candidate is on a genome that generally isn't even part of the training
	 * corpus, so there is no self-contribution to exclude). Any real sequence
	 * position is a small non-negative int, so a very negative sentinel can never
	 * collide with one. */
	private static final int NO_TRUE_POS=Integer.MIN_VALUE;

	/**
	 * Boundary enrichment-profile feature -- v1 MINIMAL design (Brian+Noire, Aug 21,
	 * superseding a same-day 6-value "full profile" version: the NN scores each
	 * candidate in the search range independently and picks the strongest, so
	 * the full 5-position ratio profile is redundant; magnitude + center-vs-best-
	 * rival is all a per-candidate score needs, and fewer inputs protects a tiny NN
	 * from overfit). For a CANDIDATE boundary position, reads the RAW table counts
	 * at 5 LOCAL positions spanning candidateBoundaryPos-2..+2 (A,B,C,D,E; C is
	 * the candidate's own position) -- this local radius-2 smoothing window is
	 * independent of the per-locus candidate SEARCH range (TrnaBoundaryVectorGen's
	 * START_OFFSETS/STOP_OFFSETS, which as of Aug 22 2026 extend to -3/+4 -- a
	 * per-candidate local neighborhood profile, not a statement about how far
	 * apart candidates are generated), then emits 3 values:
	 * <pre>  F = 1 + max(A,B,C,D,E)
	 *   bestOther = max(A,B,D,E)   (excludes C)
	 *   [ log2(F)/8, (C+1)/F, (bestOther+1)/F ]</pre>
	 * Index 0 (magnitude) says how populated the strongest nearby position is (lets
	 * the NN down-weight sparse/novel regions). Index 1 (center ratio) says how
	 * enriched the candidate itself is relative to the local max. Index 2 (best
	 * rival ratio) says how enriched the single strongest NEIGHBOR is -- together,
	 * index1 near 1 / index2 near 0 says "this candidate is correct"; index1 low /
	 * index2 near 1 says "a neighbor is more likely the true boundary." DIVISOR=8
	 * is Brian's chosen scaling constant (not recalibrated from the observed corpus
	 * max, which runs higher -- log2(~41689)/8~=1.9 at the locked-table's actual
	 * max count -- deliberately left as specified rather than re-tuned; the NN
	 * normalizes what it's given).
	 *
	 * <p>TRAINING overload (this one): trueBoundaryPos is the record's OWN true
	 * boundary (unshifted). If the candidate is within +-2 of the true position,
	 * exactly one of the 5 LOCAL window positions coincides with it and gets
	 * LEAVE-ONE-OUT: that position's raw count is reduced by 1, because this
	 * record's own true kmer was itself counted into the table during the corpus
	 * build (same window formula), so its own count always includes a trivial
	 * self-match that would otherwise bias its own training label. For a
	 * candidate further than +-2 from true (possible since Aug 22 2026's
	 * asymmetric search ranges extend to -3/+4), none of the 5 local positions
	 * coincide with trueBoundaryPos, so leave-one-out simply never fires -- correct,
	 * not a gap: there is no self-contribution in this local window to exclude.
	 *
	 * <p>Use the INFERENCE overload (no trueBoundaryPos) when scoring a real,
	 * unknown candidate (TrnaBoundaryScorer) -- there is no ground truth to leave
	 * out there, and the query genome is generally not part of the training corpus
	 * at all, so there is no self-contribution to exclude.
	 */
	public static float[] enrichmentProfile(byte[] seq, int candidateBoundaryPos, int trueBoundaryPos,
			BoundaryType type, int insideCount, int outsideCount, NinemerTable table){
		return enrichmentProfileCore(seq, candidateBoundaryPos, trueBoundaryPos, type, insideCount, outsideCount, table);
	}

	/** Inference-mode overload: no ground truth, no leave-one-out. See the training
	 * overload's javadoc for the full feature design. */
	public static float[] enrichmentProfile(byte[] seq, int candidateBoundaryPos,
			BoundaryType type, int insideCount, int outsideCount, NinemerTable table){
		return enrichmentProfileCore(seq, candidateBoundaryPos, NO_TRUE_POS, type, insideCount, outsideCount, table);
	}

	private static final double LOG2=Math.log(2);
	private static final double MAGNITUDE_DIVISOR=8.0;

	private static float[] enrichmentProfileCore(byte[] seq, int candidateBoundaryPos, int trueBoundaryPos,
			BoundaryType type, int insideCount, int outsideCount, NinemerTable table){
		if(table==null){
			float[] out=new float[3];
			java.util.Arrays.fill(out, PENDING_DORI);
			return out;
		}
		final int[] raw=new int[5];
		for(int i=-2; i<=2; i++){
			final int pos=candidateBoundaryPos+i;
			final int windowStart=kmerWindowStart(pos, type, insideCount, outsideCount);
			int c=table.rawCount(seq, windowStart);
			assert(c>=0) : "NinemerTable.rawCount contract violated (must be >=0, never a sentinel): "+c;
			if(pos==trueBoundaryPos && c>=1){
				//Leave-one-out: exclude this record's own contribution to its own true kmer. Only
				//fires when c>=1 -- do NOT assume c is always >=1 here. rawCount's own documented
				//contract gives THREE legitimate reasons for c==0 even at a record's own true
				//position: never-observed, OUT OF BOUNDS (a k-mer window can extend past a
				//record's own sequence -- e.g. STOP's window peeks 1bp into the flank, which can
				//run out at a contig-edge-clipped record whose rflank was trimmed below the
				//window's needs), or a non-ACGT base in the window. All three mean the SAME
				//window-computation formula would have skipped this exact record during the
				//table's own build pass too, so there is nothing to subtract -- c==0 is the
				//correct, already-leave-one-out'd value, not a violated invariant. Discovered via
				//a real AssertionError at corpus scale (50K-record prefix of train_gc.fa, Aug 21)
				//-- the smaller smoke/regression fixtures never happened to hit this edge case.
				c=c-1;
			}
			raw[i+2]=c;
		}
		final int C=raw[2];//candidate's own position (A,B,C,D,E = raw[0..4])
		final int bestOther=Math.max(Math.max(raw[0], raw[1]), Math.max(raw[3], raw[4]));//max(A,B,D,E), C excluded
		final int max=Math.max(C, bestOther);//max(A,B,C,D,E)
		final float F=1+max;
		return new float[]{
			(float)(Math.log(F)/LOG2/MAGNITUDE_DIVISOR),
			(C+1)/F,
			(bestOther+1)/F
		};
	}

	/**
	 * Sentinel for "the fuzziness machinery works and models are loaded, but THIS
	 * candidate didn't align confidently to any of them" -- distinct from
	 * PENDING_DORI, which means the whole feature category is unavailable (no
	 * models loaded at all). Deliberately a value no real computation could ever
	 * produce (unlike e.g. 0.5, which sits inside the plausible 0-1ish range a
	 * genuine coverage ratio can take and so could be mistaken for real data) --
	 * same "sentinel must be out-of-band" principle PENDING_DORI already follows
	 * in this file, applied to a semantically different situation (Noire's
	 * review, Aug 21: the NN needs to tell "uncertain because novel" apart from
	 * "feature not computed at all").
	 */
	public static final float NO_MODEL_MATCH=-2f;

	/**
	 * Aligns a candidate tRNA span against one HBM model, returning the alignment
	 * in MODEL coordinates (r.start/r.stop index into model.ref[]/model.del[]),
	 * or null if no confident alignment exists. Deliberately mirrors
	 * TrnaConsensusBuilder.scoreAgainstModel's exact logic -- same 0.2 identity
	 * floor, same length-branch choosing which side is the aligner's query -- so
	 * this picks the identical alignment production model-scoring would use; a
	 * different floor or branch here would silently diverge from what TrnaCaller
	 * actually does at inference time. toAlignedRead/invertToModelFrame are
	 * package-private in TrnaConsensusBuilder but reachable here (same package),
	 * so no changes to TrnaConsensusBuilder.java are needed to reuse them.
	 */
	public static Read alignToModelFrame(byte[] candidate, BaseGraph model){
		return alignToModelFrame(candidate, model, new AlignmentStats(true));
	}

	/** Overload exposing the underlying AlignmentStats (matches/subs/ins/dels/identity) to the
	 * caller -- needed by the tip-adjustment approximation below, which derives its running
	 * (matches,alnLen) counters from the ONE real alignment instead of realigning per shift. */
	public static Read alignToModelFrame(byte[] candidate, BaseGraph model, AlignmentStats stats){
		if(model==null || candidate==null){return null;}
		final byte[] cons=model.original;
		stats.doTrace=true;
		final Read r;
		if(candidate.length<=cons.length){
			float id=ScrabbleAligner.alignAndTraceStatic(candidate, cons, stats);
			if(stats.matchString==null || id<0.2f){return null;}
			r=TrnaConsensusBuilder.toAlignedRead(candidate, stats, "candidate", 0);
		}else{
			float id=ScrabbleAligner.alignAndTraceStatic(cons, candidate, stats);
			if(stats.matchString==null || id<0.2f){return null;}
			r=TrnaConsensusBuilder.invertToModelFrame(candidate, stats);
		}
		return r;
	}

	/** Sum of ref+del training-read counts at one model position -- the same
	 * formula as BaseGraph.countSum(int), reimplemented here because that method
	 * is package-private in consensus (not reachable from prok) and widening a
	 * shared/high-blast-radius class for a 2-line accessor isn't worth it when
	 * BaseNode.countSum and BaseGraph.ref[]/del[] are already public. Out-of-range
	 * positions return 0 (not an error): a position beyond the model genuinely has
	 * no training coverage to report, which is honest information, not a guess. */
	private static int modelCoverage(BaseGraph model, int rpos){
		if(rpos<0 || rpos>=model.ref.length){return 0;}
		return model.ref[rpos].countSum+model.del[rpos].countSum;
	}

	/** True max training coverage anywhere in the model -- the external reimplementation
	 * of BaseGraph.maxDepth()'s exact formula (package-private in consensus, not reachable
	 * from prok). NOT the same as "coverage at the mid-model position": some shipped
	 * consensus pivots contain a literal internal N-run (an assembly gap in the source
	 * genome that survived pivot selection -- confirmed by inspection, e.g.
	 * tRNA_consensus_GUA_c0's pivot is 76 real bases + ~220 N's + 30 real bases), and
	 * buildBaseGraph never adds the pivot's own read to the graph (it's the reference,
	 * skipped in the cluster loop), so a position landing inside that N-gap can have
	 * ZERO coverage even at the geometric middle of a 297bp model. A true max-scan finds
	 * the model's actual depth plateau regardless of where in the model it falls. */
	private static int modelMaxCoverage(BaseGraph model){
		int max=0;
		for(int i=0; i<model.ref.length; i++){
			final int c=modelCoverage(model, i);
			if(c>max){max=c;}
		}
		return max;
	}

	/**
	 * Model-tip fuzziness (Brian+Noire, Aug 21): how well the HBM's training
	 * coverage supports THIS candidate's boundary position -- both in absolute
	 * terms (near the model's plateau depth vs past its coverage cliff) and
	 * locally (how steep the cliff is right at this candidate's own endpoint).
	 * "If there's a sharp cliff but the tRNA in question doesn't end right at the
	 * cliff, that too is useful to know" (Brian) -- the 3-value local profile
	 * captures both the cliff's sharpness AND where this candidate sits on it.
	 *
	 * <p>3 values: countSum at alignPos-1, alignPos, alignPos+1 (model
	 * coordinates), each divided by the model's TRUE max coverage anywhere
	 * (modelMaxCoverage -- NOT a fixed mid-model position; some shipped
	 * consensus pivots have an internal N-gap that a fixed position could land
	 * inside, see modelMaxCoverage's javadoc). alignPos is r.start for a
	 * START-boundary row, r.stop for a STOP-boundary row.
	 *
	 * <p>PENDING_DORI when model==null: the whole capability is off (no HBM
	 * library loaded). NO_MODEL_MATCH when alignment fails for this specific
	 * candidate (model is loaded and other candidates may align fine) -- see
	 * NO_MODEL_MATCH's own javadoc for why this needs a distinct sentinel from
	 * PENDING_DORI rather than reusing it.
	 *
	 * @param candidate The candidate span's bases (already trimmed to [s,e], same
	 *   convention as stemFeature/aniFeature).
	 * @param model The single best-matching HBM model for this locus (same model
	 *   chosen for the ANI feature via the production TrnaKmerIndex shortlist --
	 *   not re-shortlisted here).
	 * @param useStart True for a START-boundary row (use the alignment's start
	 *   endpoint), false for a STOP-boundary row (use its stop endpoint).
	 */
	/**
	 * Tip-adjustment approximation for ANI under a +-1/+-2bp boundary shift, given ONE real
	 * alignment of the unshifted candidate (Brian, Aug 21 2026 -- the CallGenes-facing
	 * shipping design: one alignment per locus, not one per shifted candidate, so the NN
	 * integration stays at CallGenes' real per-genome speed budget). O(1) per shift: at most
	 * 2 byte comparisons plus one division, vs a full realignment.
	 *
	 * <p>NOT used for TrnaBoundaryVectorGen's training vectors -- those brute-force realign
	 * every candidate (correct, already generated). This is a deliberate train/inference
	 * approximation gap Brian accepted for v1 ("~1-2% error is well within what the NN can
	 * absorb -- it's learning relative differences between candidates, not absolute values").
	 *
	 * <p>Model: under the approximation that the alignment near the boundary is a simple 1:1
	 * (ungapped) correspondence, shifting the candidate boundary by a signed offset shifts the
	 * model-coordinate alignment position by the SAME offset (newAlignPos=alignPos+offset --
	 * verified algebraically identical whether the shift GROWS or SHRINKS the candidate, see
	 * the derivation this method's implementation follows). Growing the candidate (extending
	 * past the base alignment's own edge) pairs each newly-included window base against the
	 * model position it would extrapolate to; shrinking it (trimming from the edge) pairs each
	 * newly-EXCLUDED window base against the model position it ACTUALLY occupied in the base
	 * alignment (exact, not extrapolated) and reverses its contribution.
	 *
	 * @param window Genomic window containing the candidate.
	 * @param boundaryPos The base alignment's own boundary position (the candidate's actual s
	 *   for START, or e for STOP) -- i.e. the position BEFORE evaluating any shift.
	 * @param alignPos The base alignment's model-coordinate endpoint at boundaryPos (r.start
	 *   for START, r.stop for STOP, from the ONE real alignToModelFrame call).
	 * @param growDirection +1 if a POSITIVE offset grows the candidate at this boundary (STOP:
	 *   increasing e grows it), -1 if a NEGATIVE offset grows it (START: decreasing s grows it).
	 * @param offset The signed shift to evaluate -- one of 0,-1,1,-2,2.
	 * @param modelBases The SAME reference sequence the base alignment was computed against
	 *   (model.original -- the BaseGraph's pivot, NOT library[bestModel]'s cleaned consensus;
	 *   see this method's class-level design note on why ANI here measures identity-to-pivot,
	 *   a deliberate choice to let one alignment serve both ani and fuzziness).
	 * @param baseMatches @param baseAlnLen The base alignment's own match count and total
	 *   aligned length (matches+subs+ins+dels), from the AlignmentStats of that ONE alignment.
	 * @return {ani, newAlignPos} -- ani is Float.NaN if the resulting alnLen is non-positive
	 *   (degenerate; should not occur for any real candidate at these small offsets).
	 */
	public static float[] tipAdjustAni(byte[] window, int boundaryPos, int alignPos, int growDirection,
			int offset, byte[] modelBases, int baseMatches, int baseAlnLen){
		int matches=baseMatches, alnLen=baseAlnLen;
		if(offset!=0){
			final boolean growing=(Integer.signum(offset)==growDirection);
			if(growing){
				for(int i=1; i<=Math.abs(offset); i++){
					final int step=growDirection*i;
					final int wPos=boundaryPos+step, mPos=alignPos+step;
					if(wPos<0 || wPos>=window.length || mPos<0 || mPos>=modelBases.length){continue;}
					alnLen++;
					if(window[wPos]==modelBases[mPos] && window[wPos]!='N'){matches++;}
				}
			}else{
				for(int i=0; i<Math.abs(offset); i++){
					final int step=(-growDirection)*i;
					final int wPos=boundaryPos+step, mPos=alignPos+step;
					if(wPos<0 || wPos>=window.length || mPos<0 || mPos>=modelBases.length){continue;}
					alnLen--;
					if(window[wPos]==modelBases[mPos] && window[wPos]!='N'){matches--;}
				}
			}
		}
		final int newAlignPos=alignPos+offset;
		final float ani=(alnLen>0) ? matches/(float)alnLen : Float.NaN;
		return new float[]{ani, newAlignPos};
	}

	/** Tip fuzziness directly from a (possibly tip-adjusted) model alignment position -- no
	 * realignment. Companion to tipAdjustAni: given newAlignPos (from tipAdjustAni's second
	 * return value), reads countSum at newAlignPos-1/newAlignPos/newAlignPos+1, plateau-
	 * normalized, identically to tipFuzzinessFeature's real-alignment version. */
	public static float[] fuzzinessFromAlignPos(BaseGraph model, int alignPos){
		if(model==null){
			float[] out=new float[3];
			java.util.Arrays.fill(out, PENDING_DORI);
			return out;
		}
		final float plateau=modelMaxCoverage(model);
		assert(plateau>0) : "Model '"+model.name+"' (ref.length="+model.ref.length+") has ZERO training "
			+"coverage anywhere -- see tipFuzzinessFeature's identical assert for the full rationale.";
		return new float[]{
			modelCoverage(model, alignPos-1)/plateau,
			modelCoverage(model, alignPos)/plateau,
			modelCoverage(model, alignPos+1)/plateau
		};
	}

	public static float[] tipFuzzinessFeature(byte[] candidate, BaseGraph model, boolean useStart){
		if(model==null){
			float[] out=new float[3];
			java.util.Arrays.fill(out, PENDING_DORI);
			return out;
		}
		assert(model.pad==0) : "tipFuzzinessFeature assumes unpadded models -- "
			+"TrnaConsensusBuilder.loadTextModels hardcodes pad=0 for every model loaded from the shipped "
			+".hbm text format, so a nonzero pad here means a model was built/loaded a different way and the "
			+"countSum position mapping below (which assumes ref[i] IS the i'th model base, no N-padding "
			+"offset) would silently read the wrong position. model="+model.name+", pad="+model.pad;
		final Read r=alignToModelFrame(candidate, model);
		if(r==null){
			float[] out=new float[3];
			java.util.Arrays.fill(out, NO_MODEL_MATCH);
			return out;
		}
		final int alignPos=(useStart ? r.start : r.stop);
		final float plateau=modelMaxCoverage(model);
		assert(plateau>0) : "Model '"+model.name+"' (ref.length="+model.ref.length+") has ZERO training "
			+"coverage anywhere -- buildBaseGraph should have added at least one non-pivot cluster member's "
			+"alignment (TrnaConsensusBuilder minClusterSize gates this), so an all-zero model means either "
			+"a genuinely degenerate cluster or a loader bug, not a normal model shape.";
		return new float[]{
			modelCoverage(model, alignPos-1)/plateau,
			modelCoverage(model, alignPos)/plateau,
			modelCoverage(model, alignPos+1)/plateau
		};
	}
}
