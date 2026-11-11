package assemble;

import ukmer.Kmer;
import structures.IntList;

/**
 * Non-mutating +/-1 homopolymer candidate verification against an immutable count
 * table. This is evidence for a candidate, NOT a sequencing-error/allele classifier.
 * In particular a sufficiently rare true allele can have the same evidence as an
 * error. Do not apply proposals without a separately validated correction policy.
 * One instance belongs to one worker; the lookup must not retain or change its Kmer.
 * @author Fischl
 */
public final class HomopolymerIndelProposal{

	/** Borrowed immutable table access using this worker's reusable query word. */
	public interface CountLookup{
		/**
		 * Return accumulated depth, or -1 for an absent key, without changing the
		 * table or the query's sequence arrays. Use the TABLE's canonicalization
		 * and layout, not an independently chosen hash/orientation. Kmer packing
		 * and core-mask settings must remain fixed for this worker/table lifetime.
		 * Reverse complements must have equal depth for a canonical table. A
		 * strand-specific table needs separate validation before use here.
		 */
		int count(Kmer key);
	}

	/** Outcome of the latest proposal attempt; rejection reasons do not classify truth. */
	enum Status{NOT_EVALUATED, NOT_A_RUN, EDGE, UNOBSERVABLE, UNDEFINED,
		ORIGINAL_SUPPORTED, NO_SUPPORTED_ALTERNATIVE, AMBIGUOUS, PROPOSED, UNSUPPORTED_ACTION, COMPETING_INSERTION, COMPETING_DELETION}

	/** Uses runs of at least two bases, without neighboring-edit competition vetoes. */
	HomopolymerIndelProposal(final int k_, final CountLookup lookup_, final int minSupport_,
			final int maxOriginal_, final float ratio_){
		this(k_, lookup_, minSupport_, maxOriginal_, ratio_, false);
	}

	/** Optionally allows singleton extension; competing-edit vetoes remain disabled. */
	HomopolymerIndelProposal(final int k_, final CountLookup lookup_, final int minSupport_,
			final int maxOriginal_, final float ratio_, final boolean singletons_){
		this(k_, lookup_, minSupport_, maxOriginal_, ratio_, singletons_, false);
	}

	/** Adds optional interruption-insertion competition; adjacent-deletion vetoes remain off. */
	HomopolymerIndelProposal(final int k_, final CountLookup lookup_, final int minSupport_,
			final int maxOriginal_, final float ratio_, final boolean singletons_, final boolean competing_){
		this(k_, lookup_, minSupport_, maxOriginal_, ratio_, singletons_, competing_, false);
	}

	/**
	 * Creates one worker's verifier and optional competing-profile scratch.
	 * @param k_ Exact, unrounded table word length, at least five
	 * @param lookup_ Immutable table access; must not retain the borrowed query word
	 * @param minSupport_ Inclusive minimum count for every candidate-context window
	 * @param maxOriginal_ Inclusive maximum allowed original run-spanning count
	 * @param ratio_ Finite contrast ratio, at least one; candidate support must strictly exceed it times the original maximum
	 * @param singletons_ Whether one-base runs may propose extension only
	 * @param competing_ Whether supported internal non-run insertions veto extension
	 * @param deletionCompeting_ Whether supported adjacent deletions veto shortening
	 */
	HomopolymerIndelProposal(final int k_, final CountLookup lookup_, final int minSupport_,
			final int maxOriginal_, final float ratio_, final boolean singletons_, final boolean competing_, final boolean deletionCompeting_){
		if(k_<5 || lookup_==null || minSupport_<1 || maxOriginal_<0 ||
				!(ratio_>=1) || Float.isInfinite(ratio_)){
			throw new IllegalArgumentException("Invalid homopolymer proposal count policy.");
		}
		k=k_;
		lookup=lookup_;
		minSupport=minSupport_;
		maxOriginal=maxOriginal_;
		ratio=ratio_;
		singletons=singletons_;
		key=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Kmer layout rounded requested k="+k+" to "+key.kbig);}
		competitor=competing_ ? new HomopolymerInsertionEvidence(k, lookup) : null;
		competingDepths=competing_ ? new IntList() : null;
		deletionCompetitor=deletionCompeting_ ? new HomopolymerDeletionEvidence(k, lookup) : null;
		deletionDepths=deletionCompeting_ ? new IntList() : null;
	}

	/**
	 * Return -1 (shorten), +1 (lengthen), or 0 (abstain) for [start,end), in
	 * zero-based coordinates of the supplied oriented read. This deliberately
	 * does not choose a physical edit position or an inserted-base quality;
	 * sequence-equivalent placements can carry different quality/provenance.
	 * The interval must be a maximal uppercase A/C/G/T run of length>=2, or an
	 * explicitly enabled singleton eligible only for +1 (still testing -1 evidence).
	 * Both +/-1 alternatives must be observable across both distinct run flanks.
	 * Every kmer in the candidate context, including two untouched kmers on each
	 * side, must pass the support floor and contrast against ALL original kmers
	 * spanning the run and its two distinct adjacent bases. No input bytes change.
	 * Resets diagnostics on entry; callers must not consume a result after an exception.
	 * @return One-base length delta, or zero when the current evidence abstains
	 */
	int propose(final byte[] bases, final int start, final int end){
		status=Status.NOT_EVALUATED;
		originalMax=shorterMin=longerMin=-1;
		leftDeletionMin=rightDeletionMin=-1;
		originalComplete=false;
		queries=0;
		if(bases==null || start<0 || end<start || end>bases.length){
			throw new IllegalArgumentException("Invalid run interval: "+start+".."+end);
		}
		final int run=end-start;
		if(run<1 || (run==1 && !singletons)){status=Status.NOT_A_RUN; return 0;}
		final byte base=bases[start];
		if(!defined(base)){status=Status.UNDEFINED; return 0;}
		for(int i=start+1; i<end; i++){
			if(bases[i]!=base){status=Status.NOT_A_RUN; return 0;}
		}
		if((start>0 && bases[start-1]==base) || (end<bases.length && bases[end]==base)){
			status=Status.NOT_A_RUN; return 0;
		}
		if(run>k-3){status=Status.UNOBSERVABLE; return 0;}
		final int flank=k+1;
		if(start<flank || bases.length-end<flank){status=Status.EDGE; return 0;}
		final int from=start-flank, to=end+flank;
		for(int i=from; i<to; i++){
			if(!defined(bases[i])){status=Status.UNDEFINED; return 0;}
		}
		assert(from>=0 && to<=bases.length && run+3<=k) :
			"Proposal windows require two untouched flank kmers and a kmer spanning the longer run plus both distinct neighbors.";
		originalMax=score(bases, from, run, 0, true);
		if(originalMax>maxOriginal){status=Status.ORIGINAL_SUPPORTED; return 0;}
		assert(originalComplete) : "Candidate decisions require the complete low-support original profile.";
		shorterMin=score(bases, from, run, -1, false);
		longerMin=score(bases, from, run, 1, false);
		final boolean shorter=supported(shorterMin), longer=supported(longerMin);
		if(shorter && longer){status=Status.AMBIGUOUS; return 0;}
		if(!shorter && !longer){status=Status.NO_SUPPORTED_ALTERNATIVE; return 0;}
		// Singleton -1 remains an evidence hypothesis, never an allowed action.
		if(run==1 && shorter){status=Status.UNSUPPORTED_ACTION; return 0;}
		if(shorter && deletionCompetitor!=null && competingDeletion(bases, start, end)){
			status=Status.COMPETING_DELETION; return 0;
		}
		if(longer && competitor!=null && competingInsertion(bases, start, end)){
			status=Status.COMPETING_INSERTION; return 0;
		}
		status=Status.PROPOSED;
		return shorter ? -1 : 1;
	}

	/** A nearby non-HP insertion can mimic an extra run base in pooled repeat counts.
	 * Compare both adjacent deletions in the same original span; abstain only.
	 * Vetoing a proposal may release other batch edits, so validate whole reads. */
	private boolean competingDeletion(final byte[] bases, final int start, final int end){
		assert(end-start>=2 && supported(shorterMin) && bases[start-1]!=bases[start] && bases[end]!=bases[start]) :
			"Only an observable maximal run with an accepted HP-1 hypothesis reaches adjacent-deletion evidence.";
		deletionCompetitor.profile(bases, start-k-1, end+k+1, start-1, deletionDepths);
		leftDeletionMin=deletionCompetitor.minimum; queries+=deletionDepths.size;
		deletionCompetitor.profile(bases, start-k-1, end+k+1, end, deletionDepths);
		rightDeletionMin=deletionCompetitor.minimum; queries+=deletionDepths.size;
		return supported(leftDeletionMin) || supported(rightDeletionMin);
	}

	/** Extend ambiguity detection, not the set of permitted corrections.
	 * AA may be missing an intervening G (truth AGA), not an A from AAA.
	 * Evaluate non-run-base insertions strictly INSIDE the run in the SAME context;
	 * a supported internal interruption vetoes +1. Boundary alternatives describe
	 * neighboring insertions, not missing run interruptions (T4 reads170/173).
	 * This narrower opt-in policy deliberately leaves boundary errors unresolved.
	 * Never choose or apply that other base. Singletons have no interior gaps.
	 * Proposal-level veto can change later batch/overlap selection, so this opt-in
	 * policy requires whole-read validation; it is not automatically monotone. */
	private boolean competingInsertion(final byte[] bases, final int start, final int end){
		assert(competitor!=null && competingDepths!=null && supported(longerMin)) :
			"Competing insertion checks only follow a supported HP extension and require worker-local scratch.";
		for(int gap=start+1; gap<end; gap++){
			for(byte base : ALPHABET){
				if(base==bases[start]){continue;}
				competitor.profile(bases, start-k-1, end+k+1, gap, base, competingDepths);
				queries+=competingDepths.size;
				if(supported(competitor.minimum)){return true;}
			}
		}
		return false;
	}

	/** Streams a virtual edited context without candidate arrays or per-window allocations.
	 * Original scoring takes the maximum over windows spanning the run and both
	 * neighbors; candidate scoring takes the minimum over the entire context.
	 * Singleton original scoring may stop early once support already rules out an edit.
	 * Updates query counts and marks originalComplete only after a full original scan. */
	private int score(final byte[] bases, final int from, final int run, final int delta,
			final boolean original){
		assert(delta>=-1 && delta<=1 && (!original || delta==0)) :
			"The original support profile must be unedited; candidates change only one run base.";
		final int flank=k+1, changedRun=run+delta;
		final int size=2*flank+changedRun;
		key.clearFast();
		int score=original ? 0 : Integer.MAX_VALUE, measured=0;
		for(int i=0; i<size; i++){
			final int source=i<flank ? from+i :
				i<flank+changedRun ? from+flank : from+i-delta;
			key.addRight(bases[source]);
			if(i<k-1){continue;}
			if(original && !(i-k+1<flank && i>=flank+run)){continue;}
			final int raw=lookup.count(key);
			assert(raw>=-1) : "CountLookup permits only nonnegative counts or -1 for an absent key; got "+raw;
			final int depth=Math.max(0, raw);
			queries++;
			measured++;
			score=original ? Math.max(score, depth) : Math.min(score, depth);
			// Any supported original spanning kmer already proves abstention.
			// Only the opt-in singleton path short-circuits; ordinary profiles
			// and every emitted candidate retain their exact previous semantics.
			if(original && run==1 && score>maxOriginal){return score;}
		}
		assert(measured>0) : "Both distinct run flanks must fit in a measured kmer.";
		if(original){originalComplete=true;}
		return score;
	}

	/** Requires the absolute floor and strict contrast against the current original maximum. */
	private boolean supported(final int depth){
		assert(originalMax>=0) : "Candidate contrast requires the original run-spanning profile first.";
		return depth>=minSupport && depth>originalMax*ratio;
	}

	/** Accepts only uppercase DNA bases in the full verification context. */
	private static boolean defined(final byte b){return b=='A' || b=='C' || b=='G' || b=='T';}

	/** Last proposal outcome, reset before validating each new request. */
	Status status=Status.NOT_EVALUATED;
	/** False means originalMax may be only a lower bound (or not evaluated).
	 * Emitted proposals always have a complete original profile. */
	boolean originalComplete=false;
	/** Original maximum and candidate minima (-1 until measured), plus lookups in this call. */
	int originalMax=-1, shorterMin=-1, longerMin=-1, queries=0;
	/** Minima for deleting either adjacent flank base; -1 until the competition check runs. */
	int leftDeletionMin=-1, rightDeletionMin=-1;
	/** Exact word length, candidate support floor and original support ceiling. */
	private final int k, minSupport, maxOriginal;
	/** Strict multiplicative contrast threshold applied to the original maximum. */
	private final float ratio;
	/** Whether a length-one run may propose +1; deletion remains evidence only. */
	private final boolean singletons;
	/** Immutable borrowed table access shared with optional competing-profile helpers. */
	private final CountLookup lookup;
	/** Reused rolling query word; never retained by the lookup. */
	private final Kmer key;
	/** Optional internal interruption-insertion profiler; null when disabled. */
	private final HomopolymerInsertionEvidence competitor;
	/** Worker-owned insertion-profile buffer, allocated only with competitor. */
	private final IntList competingDepths;
	/** Optional adjacent-base deletion profiler; null when disabled. */
	private final HomopolymerDeletionEvidence deletionCompetitor;
	/** Worker-owned deletion-profile buffer, allocated only with deletionCompetitor. */
	private final IntList deletionDepths;
	/** Defined insertion alternatives in deterministic A/C/G/T order. */
	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
}
