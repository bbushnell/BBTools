package prok;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;
import structures.IntList;

/**
 * Inverted k-mer index over a tRNA consensus library plus a reusable per-query counting shortlister.
 *
 * Replaces TrnaCaller's per-query allocations (int[nModels] hits + int[nModels][2] scored, ~nModels+1
 * allocations per candidate window) and its boxed-comparator Arrays.sort with: a reused int[] count
 * buffer cleared by touched-index, char[][] postings (model IDs 0..nModels-1 fit char), and an
 * O(nModels) counting sort (count desc, ties ascending id) instead of an O(n log n) comparator sort.
 *
 * The output ORDERING rule (count desc, ties by ascending model id) is unchanged from the original
 * shortlist. The COUNT VALUES are no longer bit-identical to that original as of the unique-query-k-mer
 * change below (2026-08-27): shortlist() now dedupes query k-mer occurrences before touching postings,
 * so a query region that repeats a k-mer contributes that k-mer's postings ONCE, not once per raw
 * occurrence. This is a deliberate semantics change, not a performance-only refactor -- see that
 * section's javadoc for why. The qFrac denominator (qKmers) follows the SAME dedup: it is the unique
 * query-k-mer count (qk.size), confirmed decision (Citan/Brian via Noire, 2026-08-27) -- not the raw
 * positional count the pre-dedup code used. See the qKmers assignment in shortlist() for the exact
 * asymmetric-formula note (numerator stays model-multiplicity-weighted; only the denominator is unique).
 *
 * One instance per calling thread: the count buffer is mutable and NOT thread-safe, matching the old
 * per-instance kmerIndex (TrnaCaller is constructed per calling thread).
 *
 * <p>The per-model count buffer was originally a byte[] (max 255), sized for tRNA-scale windows
 * (~120 k-mers against ~76bp mostly-unique refs, so max shared << 255). NcrnaScavenger reusing this
 * same class for a ~400bp+ family (RNase P) hit the overflow for real: a long window sharing more than
 * 255 k-mer occurrences with a similarly-long, closely-related model. Widened to int[] here (Neptune,
 * 2026-08-23) rather than saturating at 255, per the original overflow assert's own guidance ("widen
 * the counter type") -- saturating would silently flatten resolution among the top matches instead of
 * crashing loud, which is exactly the failure mode this codebase's assertion philosophy rejects. The
 * counting-sort scratch (hist/pos) is sized generously (MAX_COUNT+1) rather than unboundedly, with its
 * own loud assert if a query ever exceeds that -- see MAX_COUNT's javadoc for why a fixed bound is safe.
 *
 * @author Noire
 */
final class TrnaKmerIndex {

	/** Shortlist index k-mer length (INDEX_K).  numKmers=1<<(2*indexK). */
	private final int indexK;
	private final int numKmers;
	private final int nModels;
	/** postings[kmer] = model IDs whose sequence contains that k-mer, once PER OCCURRENCE (multiplicity
	 * preserved so the per-query count is occurrence-weighted, exactly as the old int[][] index did). */
	private final char[][] postings;
	private static final char[] EMPTY=new char[0];
	private static final int[] EMPTY_INT=new int[0];

	/*-- Per-query reusable scratch (mutable; per-thread) --*/
	/** counts[model] = shared k-mer count with the current query, summed over the query's UNIQUE k-mer
	 * types (2026-08-27 change), each weighted by that model's OWN occurrence multiplicity for that
	 * k-mer type (postings[kmer] is untouched/undeduped). So this is no longer occurrence-weighted on
	 * the query side, only on the model side -- bounded by sum over the query's unique k-mer types of
	 * that model's occurrence count for each, which is <= that model's total index-k-mer occurrences.
	 * int (not byte/short): see the class javadoc for why a fixed-width counter is unsafe once callers
	 * aren't all tRNA-scale. */
	private final int[] counts;
	/** Model IDs touched (count>0) this query, so the buffer clears in O(touched), not O(nModels). */
	private final IntList touched=new IntList();
	/** Raw query k-mer values for the current query (one entry per position), sorted+condensed to the
	 * UNIQUE set in-place before postings are touched (see shortlist()). Reused across queries via
	 * clear() -- deliberately never shrunk, to keep this path allocation-free after warmup. */
	private final IntList queryKmers=new IntList();
	/** Upper bound on a single query's shared-count with one model, for sizing hist/pos (the counting-sort
	 * scratch) as a fixed array rather than reallocating per query. A count this high would mean a single
	 * candidate window shares 65536+ k-mer occurrences with one model -- at k=INDEX_K (5-7bp) that needs a
	 * query with tens of thousands of valid bases, far beyond any conserved ncRNA family (even 16S/23S-scale
	 * callers don't route through this shortlist). Enforced by a loud assert in shortlist(), not assumed. */
	private static final int MAX_COUNT=65535;
	/** Counting-sort scratch indexed by count value 0..MAX_COUNT (reused across queries). */
	private final int[] hist=new int[MAX_COUNT+1];
	private final int[] pos=new int[MAX_COUNT+1];
	/** Models sorted count-desc / id-asc for the current query (reused; valid entries = touched.size). */
	private final int[] sortedOut;

	/*-- Cutoff policy (constant per run; resolved by TrnaCaller from its static flags) --*/
	private final boolean adaptive;
	private final float floor, topFrac, qFrac;
	private final int fixedMinHits;

	/*-- Reporting counters + last-query stats (for REFHIST / SHORTLIST_STATS) --*/
	private long queries=0, totalShortlisted=0;
	private int lastMaxShared=0;

	TrnaKmerIndex(byte[][] library, int indexK_, boolean adaptive_, float floor_, float topFrac_,
			float qFrac_, int fixedMinHits_){
		indexK=indexK_;
		numKmers=1<<(2*indexK);
		nModels=library.length;
		adaptive=adaptive_; floor=floor_; topFrac=topFrac_; qFrac=qFrac_; fixedMinHits=fixedMinHits_;
		counts=new int[nModels];
		sortedOut=new int[nModels];
		postings=build(library);
	}

	/** Builds the inverted index in two linear passes (size, then fill) -- no ArrayList<Integer> boxing.
	 * refId is added once per OCCURRENCE of each k-mer, matching the old buildKmerIndex. */
	private char[][] build(byte[][] library){
		final int kmask=numKmers-1;
		final byte[] bton=AminoAcid.baseToNumber;
		final int[] sizes=new int[numKmers];
		for(int refId=0; refId<library.length; refId++){
			final byte[] seq=library[refId];
			int kmer=0, len=0;
			for(int j=0; j<seq.length; j++){
				final int x=bton[seq[j]];
				if(x>=0){
					kmer=((kmer<<2)|x)&kmask; len++;
					if(len>=indexK){sizes[kmer]++;}
				}else{len=0; kmer=0;}
			}
		}
		final char[][] idx=new char[numKmers][];
		for(int i=0; i<numKmers; i++){idx[i]=(sizes[i]==0 ? EMPTY : new char[sizes[i]]);}
		final int[] fill=new int[numKmers];
		for(int refId=0; refId<library.length; refId++){
			final byte[] seq=library[refId];
			int kmer=0, len=0;
			for(int j=0; j<seq.length; j++){
				final int x=bton[seq[j]];
				if(x>=0){
					kmer=((kmer<<2)|x)&kmask; len++;
					if(len>=indexK){final char[] p=idx[kmer]; p[fill[kmer]++]=(char)refId;}
				}else{len=0; kmer=0;}
			}
		}
		return idx;
	}

	/** Returns the shortlisted model IDs (count desc / id asc), cut to the resolved minHits + topN with
	 * the keep-top-1 fallback. Ordering rule is unchanged from the original shortlistByKmer; COUNT VALUES
	 * are not bit-identical to it for queries with repeated k-mers -- see the class javadoc and the
	 * unique-query-k-mer block below. */
	int[] shortlist(byte[] seq, int topN){
		queries++;
		//Clear the previous query's counts in O(touched).
		final int[] cnt=counts;
		for(int t=0, ts=touched.size; t<ts; t++){cnt[touched.get(t)]=0;}
		touched.clear();

		//Pass 1: collect every valid query k-mer (one entry per position, same enumeration as before),
		//then sort+dedup to the UNIQUE set. Brian's directive (2026-08-27): postings must be traversed
		//once per unique query k-mer, not once per raw occurrence -- the old per-position traversal let a
		//repetitive QUERY region revisit the same postings list arbitrarily many times, so a model's
		//shared count was bounded only by seq.length, not by the model's own k-mer occurrences/reference
		//length. Deduping the query side (this block) while leaving postings[kmer] itself untouched
		//(model-side multiplicity preserved -- a model whose sequence repeats a k-mer still contributes
		//that many entries from ONE unique-query-k-mer visit) restores the intended bound.
		final IntList qk=queryKmers;
		qk.clear();
		final int kmask=numKmers-1;
		final byte[] bton=AminoAcid.baseToNumber;
		int kmer=0, len=0;
		for(int j=0; j<seq.length; j++){
			final int x=bton[seq[j]];
			if(x>=0){
				kmer=((kmer<<2)|x)&kmask; len++;
				if(len>=indexK){qk.add(kmer);}
			}else{len=0; kmer=0;}
		}
		qk.sort();
		//condense(), NOT shrinkToUnique(): shrinkToUnique() calls shrink(), which reallocates the backing
		//array whenever size!=array.length -- true on nearly every call once deduped, which would put a
		//per-query allocation on this hot reused-scratch path. condense() alone leaves the sorted-unique
		//values in the existing array's prefix [0,qk.size) with no reallocation. DEVIATION FLAGGED for
		//Brian/Citan: the directive named shrinkToUnique() explicitly; this substitutes the non-allocating
		//half of it for the class's stated zero-alloc-after-warmup design. Functionally identical sorted-
		//unique output, no established API contract, and one already-verified sort() call.
		qk.condense();

		//Pass 2: traverse postings ONCE per unique query k-mer.
		for(int qi=0, qs=qk.size; qi<qs; qi++){
			final char[] refIds=postings[qk.get(qi)];
			for(int i=0; i<refIds.length; i++){
				final int m=refIds[i];
				if(cnt[m]==0){touched.add(m);}
				cnt[m]++;
				assert(cnt[m]<=MAX_COUNT) : KillSwitch.assertDie("Shared k-mer count for model "+m
					+" exceeded MAX_COUNT="+MAX_COUNT+" (TrnaKmerIndex.MAX_COUNT) -- the counting-sort "
					+"scratch (hist/pos) is sized for this bound; a query this repetitive against one "
					+"model is unexpected for any conserved ncRNA family this class currently serves. "
					+"Runs on a CallGenes worker thread -- a plain AssertionError here would only kill "
					+"this thread and could hang the pipeline's producer/consumer, so this must halt "
					+"the whole process instead (bbtools skill assertDie idiom).");
			}
		}

		final int distinct=touched.size;
		if(distinct==0){lastMaxShared=0; return EMPTY_INT;}

		//maxShared over touched.
		int maxShared=0;
		for(int t=0; t<distinct; t++){final int c=cnt[touched.get(t)]; if(c>maxShared){maxShared=c;}}
		lastMaxShared=maxShared;

		//Resolve the shortlist cutoff. qKmers is the UNIQUE valid query-k-mer count (qk.size after
		//condense), NOT the raw positional count (seq.length-indexK+1) the pre-dedup code used --
		//confirmed decision (Citan/Brian via Noire, 2026-08-27): the qFrac denominator follows the same
		//query-side dedup as the numerator's counting pass, so qFrac*qKmers means "fraction of the
		//query's DISTINCT k-mer types", not "fraction of raw query positions". This IS an asymmetric
		//formula -- the numerator (maxShared/cnt[m]) is still weighted by MODEL-side postings
		//multiplicity, only the query side is deduped -- documented here since it's not a simple set
		//fraction. (Earlier version of this file computed qKmers positionally; that was wrong per the
		//confirmed spec, not an open question -- see current/prok/old/TrnaKmerIndex_archive_20260827_
		//intermediate-dedup-positional-denom.java for the superseded intermediate state.)
		final int minHits;
		if(adaptive){
			final int qKmers=Tools.max(1, qk.size);
			final double v=Math.max(floor, Math.max(topFrac*maxShared, qFrac*qKmers));
			minHits=(int)Math.ceil(v);
		}else{
			minHits=fixedMinHits;
		}

		//Counting sort of the touched models: count desc, ties ascending id.
		final int[] h=hist, p=pos;
		for(int c=1; c<=maxShared; c++){h[c]=0;}
		for(int t=0; t<distinct; t++){h[cnt[touched.get(t)]]++;}
		int acc=0;
		for(int c=maxShared; c>=1; c--){p[c]=acc; acc+=h[c];}//higher counts occupy lower output indices
		final int[] out=sortedOut;
		for(int m=0; m<nModels; m++){final int c=cnt[m]; if(c>0){out[p[c]++]=m;}}
		//out[0..acc-1] == touched models in count-desc / id-asc order (acc==distinct).

		//Apply topN + minHits cutoff with the keep-top-1 fallback.
		final int limit=Tools.min(topN, nModels);
		int count=0;
		while(count<limit && count<acc && cnt[out[count]]>=minHits){count++;}
		if(count<1 && acc>0){count=1;}//keep the single best if anything was shared (old: scored[0][1]>0)

		totalShortlisted+=count;
		final int[] result=new int[count];
		System.arraycopy(out, 0, result, 0, count);
		return result;
	}

	/*-- Accessors for TrnaCaller instrumentation + reporting --*/
	/** Shared k-mer count of a model with the LAST query (valid until the next shortlist call). */
	int lastSharedCount(int model){return (model>=0 && model<nModels) ? counts[model] : 0;}
	int lastMaxShared(){return lastMaxShared;}
	int lastDistinctRefs(){return touched.size;}
	long queriesProcessed(){return queries;}
	/** Total models placed on shortlists = upper bound on alignments issued (the true align count is
	 * TrnaCaller.alignmentCount, which the borderline early-exit can make smaller). */
	long totalShortlisted(){return totalShortlisted;}
}
