package prok;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;
import structures.IntList;

/**
 * Inverted k-mer index over a tRNA consensus library plus a reusable per-query counting shortlister.
 *
 * Replaces TrnaCaller's per-query allocations (int[nModels] hits + int[nModels][2] scored, ~nModels+1
 * allocations per candidate window) and its boxed-comparator Arrays.sort with: a reused byte[] count
 * buffer cleared by touched-index, char[][] postings (model IDs 0..nModels-1 fit char), and an
 * O(nModels) counting sort (count desc, ties ascending id) instead of an O(n log n) comparator sort.
 *
 * The output order is BIT-IDENTICAL to the old shortlist (count desc, ties by ascending model id), so
 * a caller's order-sensitive early-exit (the borderline id>=ID_PASS break in alignWindow) selects the
 * same model and issues the same number of alignments -- this is a pure performance change with zero
 * behavioral difference (validated by identical GFF + identical alignment count).
 *
 * One instance per calling thread: the count buffer is mutable and NOT thread-safe, matching the old
 * per-instance kmerIndex (TrnaCaller is constructed per calling thread).
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
	/** counts[model] = occurrence-weighted shared k-mer count with the current query (read via &0xFF).
	 * A byte suffices for the realistic tRNA range (window ~120 k-mers, refs ~76bp mostly-unique, so
	 * max shared << 255); an overflow past 255 wraps to 0 and is caught loud by an assert under -ea. */
	private final byte[] counts;
	/** Model IDs touched (count>0) this query, so the buffer clears in O(touched), not O(nModels). */
	private final IntList touched=new IntList();
	/** Counting-sort scratch indexed by count value 0..255 (reused across queries). */
	private final int[] hist=new int[256];
	private final int[] pos=new int[256];
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
		counts=new byte[nModels];
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
	 * the keep-top-1 fallback.  Byte-identical to the old shortlistByKmer. */
	int[] shortlist(byte[] seq, int topN){
		queries++;
		//Clear the previous query's counts in O(touched).
		final byte[] cnt=counts;
		for(int t=0, ts=touched.size; t<ts; t++){cnt[touched.get(t)]=0;}
		touched.clear();

		//Count occurrence-weighted shared k-mers into the byte buffer, tracking touched models.
		final int kmask=numKmers-1;
		final byte[] bton=AminoAcid.baseToNumber;
		int kmer=0, len=0;
		for(int j=0; j<seq.length; j++){
			final int x=bton[seq[j]];
			if(x>=0){
				kmer=((kmer<<2)|x)&kmask; len++;
				if(len>=indexK){
					final char[] refIds=postings[kmer];
					for(int i=0; i<refIds.length; i++){
						final int m=refIds[i];
						if(cnt[m]==0){touched.add(m);}
						cnt[m]++;//branch-free; wraps at 256, and the assert catches that under -ea
						assert(cnt[m]!=0) : KillSwitch.assertDie("tRNA shortlist byte count overflowed a "
							+"byte (model "+m+"); widen the counter type.");
					}
				}
			}else{len=0; kmer=0;}
		}

		final int distinct=touched.size;
		if(distinct==0){lastMaxShared=0; return EMPTY_INT;}

		//maxShared over touched (true value; saturation asserted impossible above).
		int maxShared=0;
		for(int t=0; t<distinct; t++){final int c=cnt[touched.get(t)]&0xFF; if(c>maxShared){maxShared=c;}}
		lastMaxShared=maxShared;

		//Resolve the shortlist cutoff (identical to the old fixed/adaptive logic).
		final int minHits;
		if(adaptive){
			final int qKmers=Tools.max(1, seq.length-indexK+1);
			final double v=Math.max(floor, Math.max(topFrac*maxShared, qFrac*qKmers));
			minHits=(int)Math.ceil(v);
		}else{
			minHits=fixedMinHits;
		}

		//Counting sort of the touched models: count desc, ties ascending id.
		final int[] h=hist, p=pos;
		for(int c=1; c<=maxShared; c++){h[c]=0;}
		for(int t=0; t<distinct; t++){h[cnt[touched.get(t)]&0xFF]++;}
		int acc=0;
		for(int c=maxShared; c>=1; c--){p[c]=acc; acc+=h[c];}//higher counts occupy lower output indices
		final int[] out=sortedOut;
		for(int m=0; m<nModels; m++){final int c=cnt[m]&0xFF; if(c>0){out[p[c]++]=m;}}
		//out[0..acc-1] == touched models in count-desc / id-asc order (acc==distinct).

		//Apply topN + minHits cutoff with the keep-top-1 fallback.
		final int limit=Tools.min(topN, nModels);
		int count=0;
		while(count<limit && count<acc && (cnt[out[count]]&0xFF)>=minHits){count++;}
		if(count<1 && acc>0){count=1;}//keep the single best if anything was shared (old: scored[0][1]>0)

		totalShortlisted+=count;
		final int[] result=new int[count];
		System.arraycopy(out, 0, result, 0, count);
		return result;
	}

	/*-- Accessors for TrnaCaller instrumentation + reporting --*/
	/** Shared k-mer count of a model with the LAST query (valid until the next shortlist call). */
	int lastSharedCount(int model){return (model>=0 && model<nModels) ? (counts[model]&0xFF) : 0;}
	int lastMaxShared(){return lastMaxShared;}
	int lastDistinctRefs(){return touched.size;}
	long queriesProcessed(){return queries;}
	/** Total models placed on shortlists = upper bound on alignments issued (the true align count is
	 * TrnaCaller.alignmentCount, which the borderline early-exit can make smaller). */
	long totalShortlisted(){return totalShortlisted;}
}
