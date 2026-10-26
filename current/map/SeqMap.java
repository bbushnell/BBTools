package map;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicLong;

import dna.AminoAcid;
import fileIO.FileFormat;
import ml.CellNet;
import ml.ScoreSequence;
import parse.Parse;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;
import structures.Range;

/**
 * K-mer based sequence mapping structure with middle-masking support.
 * Extends LongArrayListHashMap to store sequences indexed by masked k-mers,
 * enabling approximate candidate retrieval for sequence overlap detection.
 * Build and mutate with exclusive access, then publish for concurrent queries.
 * Stored sequence arrays are borrowed and must remain unchanged while indexed.
 * Each thread reuses its result list; copy the list before another query on the
 * same map and thread if the results need to be retained. Queries return candidate
 * placements, not verified alignments: geometry and a GC heuristic filter the
 * hits, but no base-by-base mismatch count is computed here. Forward kmers are
 * indexed as supplied; reverse complements are added only when explicitly requested.
 * @author Brian Bushnell
 */
public class SeqMap extends LongArrayListHashMap<SeqPos>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an empty index with k=11, maskMiddle=1 and minCount=0. */
	public SeqMap(){this(11, 1, 0);}
	/**
	 * Creates a SeqMap with specified k-mer parameters and middle-masking.
	 * Computes a mask that gives the middle bases a common all-ones pattern.
	 *
	 * @param k_ Packed DNA kmer length, from 1 through 32
	 * @param maskMiddle_ Nonnegative number of middle bases to mask; positive values require k_&gt;maskMiddle_+1
	 * @param minCount_ Nonnegative threshold recorded for loading; add() does not enforce it
	 */
	public SeqMap(int k_, int maskMiddle_, int minCount_){
		super();

		k=k_;
		maskMiddle=maskMiddle_;
		minCount=minCount_;
		long mask=0;
		if(maskMiddle>0){
			assert(k>maskMiddle+1);
			int bits=maskMiddle*2;
			int shift=(k-maskMiddle)&(~1);//Equivalent to (x/2)*2
			mask=~((-1L)<<bits);
			mask<<=(shift);
		}
		midMask=mask;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Indexing          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds one placement for every fully defined forward kmer, applying the middle mask.
	 * Borrows the same sequence array for every placement. Undefined symbols reset the
	 * valid-base run; no kmer spans one, even when that position would be masked.
	 * Null or shorter-than-k sequences are ignored. Repeated calls add repeated entries;
	 * minCount is not applied here. Requires exclusive access to the inherited map.
	 *
	 * @param s Borrowed DNA bytes addressable by AminoAcid.baseToNumber, or null
	 * @param count Nonnegative occurrence count for every indexed placement
	 * @param score Ranking score; -1 conventionally denotes an unscored sequence
	 */
	//CrisprFinder builds once with load(), then queries the frozen index with per-thread scratch.
	//Parent constructor/resize synchronization does not make concurrent mutation safe.
	public void add(byte[] s, int count, final float score){
		if(s==null || s.length<k){return;}

		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		int len=0;

		final int code=Tools.hash(s, 22);
		final float gc=Tools.calcGC(s);

		/* Loop through s, maintaining a forward kmer via bitshifts */
		for(int i=0; i<s.length; i++){
			byte b=s[i];
			long x=symbolToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){
				len=0;
			}else{len++;}
			if(len>=k){
				long mmKmer=kmer|midMask;
				SeqPos sp=new SeqPos(s, i, count, code, gc, score);
				put(mmKmer, sp);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Candidate Queries     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Retrieves distinct sequence/start placements sharing a defined masked kmer with a region.
	 * Requires a frozen index and valid inclusive bounds for a non-null, long-enough
	 * query. Candidates must pass overlap, endpoint and GC filters; mismatches are not
	 * counted. A GC difference times overlap must be strictly below max(2.8,maxMM+0.5).
	 * Returned positions are projected sequence starts in the query array's coordinates,
	 * possibly negative. Increments queries once after the null/short-input guards.
	 *
	 * @param query DNA bytes addressable by the base lookup tables; null returns null
	 * @param a1 Inclusive region start, at least zero
	 * @param b1 Inclusive region end, at least a1 and less than query.length
	 * @param minOverlap0 Nonnegative absolute overlap floor in bases
	 * @param maxMM Mismatch-budget hint for the GC heuristic, not an alignment mismatch limit
	 * @param maxTrim Maximum uncovered region-end length on either side, in bases
	 * @param maxLopsidedness Maximum Range.lopsidedness of the region and clipped candidate
	 * @param minOverlapFractionQ Fraction of region length; multiplied by region length and truncated to int
	 * @param sort Whether to rank results by descending count and length
	 * @return Borrowed thread-local result list, or null for no matches/null/short input;
	 * the next eligible query on this map/thread clears the list. Candidate objects
	 * are distinct per result but share indexed sequence arrays
	 */
	public ArrayList<SeqPosM> fetch(byte[] query, final int a1, final int b1, final int minOverlap0, int maxMM,
			final int maxTrim, final int maxLopsidedness, float minOverlapFractionQ, boolean sort){
		final int qlen=b1-a1+1;
		//Fixed guard ordering [map/SeqMap#001]: reject null/short input before calculating region GC.
		if(query==null || query.length<k || qlen<k){return null;}
		final float qgc=Tools.calcGC(query, a1, b1);
		queries.incrementAndGet();

		final int minOverlapQ=Tools.max(minOverlap0, (int)(minOverlapFractionQ*qlen));
		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		int len=0;
		long lookups=0;

		HashSet<SeqPosM> set=localSet.get();
		ArrayList<SeqPosM> retList=localList.get();
		SeqPosM temp=localSP.get();

		assert(set.isEmpty());
		retList.clear();
		final float maxGCO=Tools.max(2.8f, maxMM+0.5f); //gcdif*overlap, proxy for mismatches

		/* Loop through s, maintaining a forward kmer via bitshifts */
		for(int i=a1; i<=b1; i++){
			byte b=query[i];
			long x=symbolToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){
				len=0;
			}else{len++;}
			if(len>=k){
				long mmKmer=kmer|midMask;
				ArrayList<SeqPos> candidates=get(mmKmer);
				if(candidates!=null){
					lookups+=addCandidates(candidates, query, a1, b1, qgc,
							minOverlapQ, maxTrim, maxLopsidedness, maxGCO, i, temp, set, retList);
				}
			}
		}

		if(!set.isEmpty()){
			assert(set.size()==retList.size());
			if(sort){Collections.sort(retList);}
		}else{retList=null;}
		set.clear();
		setQueries.addAndGet(lookups);
		return retList;
	}

	/**
	 * Scans two equal-length regions in lockstep on the calling thread, retrieving
	 * candidates from either region. Returned positions remain relative to the
	 * first region's coordinate system, so its caller can perform one subsequent alignment.
	 * Each region supplies its own GC fraction; equal kmers skip duplicate work only
	 * when both windows are valid and the GC filters agree. Bounds must be valid,
	 * inclusive and equally long. Uses the same geometry/GC heuristics as fetch(),
	 * without counting mismatches, and increments queries once after input guards.
	 *
	 * @param query DNA bytes addressable by the base lookup tables; null returns null
	 * @param a1 Inclusive first-region start, at least zero
	 * @param b1 Inclusive first-region end, at least a1 and less than query.length
	 * @param a2 Inclusive second-region start, at least zero
	 * @param b2 Inclusive second-region end below query.length, with b2-a2 equal to b1-a1
	 * @param minOverlap0 Nonnegative absolute overlap floor in bases
	 * @param maxMM Mismatch-budget hint for the GC heuristic, not an alignment mismatch limit
	 * @param maxTrim Maximum uncovered first-region-end length on either side, in bases
	 * @param maxLopsidedness Maximum Range.lopsidedness relative to the first region
	 * @param minOverlapFractionQ Fraction of region length, converted to a truncated integer overlap floor
	 * @param sort Whether to rank results by descending count and length
	 * @return Borrowed thread-local result list, or null for no matches/null/short input;
	 * another eligible query on this map/thread reuses it. Candidates are deduplicated
	 * by sequence and projected placement, irrespective of which region found them
	 */
	public ArrayList<SeqPosM> doubleFetch(byte[] query, final int a1, final int b1, final int a2, final int b2,
			final int minOverlap0, int maxMM,
			final int maxTrim, final int maxLopsidedness, float minOverlapFractionQ, boolean sort){
		final int qlen=b1-a1+1;
		assert(qlen==b2-a2+1);
		//[map/SeqMap#001] FIXED: guard before calcGC (which derefs query); was below -> null query NPE'd.
		if(query==null || query.length<k || qlen<k){return null;}
		final float qgc=Tools.calcGC(query, a1, b1);
		final float qgc2=Tools.calcGC(query, a2, b2);
		final boolean sameGC=(qgc==qgc2);
		queries.incrementAndGet();

		final int minOverlapQ=Tools.max(minOverlap0, (int)(minOverlapFractionQ*qlen));
		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, kmer2=0;
		int len=0, len2=0;
		long lookups=0;

		HashSet<SeqPosM> set=localSet.get();
		ArrayList<SeqPosM> retList=localList.get();
		SeqPosM temp=localSP.get();

		assert(set.isEmpty());
		retList.clear();
		final float maxGCO=Tools.max(2.8f, maxMM+0.5f); //gcdif*overlap, proxy for mismatches

		/* Loop through s, maintaining a forward kmer via bitshifts */
		for(int i=a1, j=a2; i<=b1; i++, j++){
			byte b=query[i], c=query[j];
			long x=symbolToNumber[b], y=symbolToNumber[c];
			kmer=((kmer<<2)|x)&mask;
			kmer2=((kmer2<<2)|y)&mask;
			if(x<0){len=0;}else{len++;}
			if(y<0){len2=0;}else{len2++;}
			final long mmKmer=kmer|midMask, mmKmer2=kmer2|midMask;
			if(len>=k){
				ArrayList<SeqPos> candidates=get(mmKmer);
				if(candidates!=null){
					lookups+=addCandidates(candidates, query, a1, b1, qgc,
							minOverlapQ, maxTrim, maxLopsidedness, maxGCO, i, temp, set, retList);
				}
			}
			//Equal keys only duplicate a lookup when the first window was valid and used the same GC filter.
			if(len2>=k && (len<k || mmKmer!=mmKmer2 || !sameGC)){
				ArrayList<SeqPos> candidates=get(mmKmer2);
				if(candidates!=null){//This a1 and b1 is intentional, to keep pos left-relative
					lookups+=addCandidates(candidates, query, a1, b1, qgc2,
							minOverlapQ, maxTrim, maxLopsidedness, maxGCO, i, temp, set, retList);
				}
			}
		}

		if(!set.isEmpty()){
			assert(set.size()==retList.size());
			if(sort){Collections.sort(retList);}
		}else{retList=null;}
		set.clear();
		setQueries.addAndGet(lookups);
		return retList;
	}

	/**
	 * Filters one kmer's candidates, then appends previously unseen sequence/start placements.
	 * Projects each sequence start as qpos-sp.pos and clips the compared interval to
	 * the query array. temp is mutated for membership tests but never inserted; only
	 * its shallow clones are retained in set/list. Both containers must represent
	 * the same accumulated results on entry and remain consistent on return.
	 *
	 * @param candidates Non-null, stable list of indexed placements with non-null sequences
	 * @param query The query sequence bytes
	 * @param a1 Start position of query region (inclusive)
	 * @param b1 End position of query region (inclusive)
	 * @param qgc GC content of the query region
	 * @param minOverlapQ Minimum overlap required with query
	 * @param maxTrim Maximum bases that can be trimmed
	 * @param maxLopsidedness Maximum alignment asymmetry allowed
	 * @param maxGCO Maximum GC difference times overlap (mismatch proxy)
	 * @param qpos Matched kmer end, projected into first-region coordinates for doubleFetch
	 * @param temp Thread-confined scratch, not a live member of set
	 * @param set Thread-confined deduplication set, mutated by adding new placements
	 * @param retList Thread-confined result list, appended in discovery order
	 * @return Candidate occurrences passing all filters before deduplication, not map lookups
	 */
	private int addCandidates(ArrayList<SeqPos> candidates, byte[] query, final int a1, final int b1, float qgc,
			int minOverlapQ, int maxTrim, int maxLopsidedness, float maxGCO, int qpos,
			SeqPosM temp, HashSet<SeqPosM> set, ArrayList<SeqPosM> retList){
		int lookups=0;
		//geometric+GC filter per candidate; temp (reused SeqPosM) is keyed into `set` to DEDUP by
		//(offset pos, seq-hash, seq) so the same alignment found via multiple k-mers is added once.
		for(SeqPos sp : candidates){
			final float gcdif=Tools.absdif(qgc, sp.gc);
			final int rlen0=sp.seq().length, pos=sp.pos();
			int a3=qpos-pos;
			int b3=Tools.min(query.length-1, a3+rlen0-1);
			a3=Tools.max(0, a3);
			final int overlap=Range.overlap(a1, b1, a3, b3);
			final int lopsidedness=Range.lopsidedness(a1, b1, a3, b3);
			final int trim=Tools.max(a3-a1, b1-b3);

			if(overlap>=minOverlapQ && trim<=maxTrim && lopsidedness<=maxLopsidedness
					&& (gcdif*overlap<maxGCO)){
				temp.setFrom(sp);
				temp.setPos(qpos-pos);
				lookups++;
				if(!set.contains(temp)){
					SeqPosM clone=temp.clone();
					set.add(clone);
					retList.add(clone);
				}
			}
		}
		return lookups;
	}

	/**
	 * Sorts indexed lists by count, length and score, then finds maximum count.
	 * Mutates the indexed lists and therefore requires exclusive access.
	 * Iterates through all stored sequence lists, sorts each by count/score,
	 * and returns the highest count value found.
	 * @return Maximum count across indexed placements, or zero for an empty index
	 */
	public int sort(){
		int max=0;
		for(ArrayList<SeqPos> list : values()){
			if(list!=null && !list.isEmpty()){
				Collections.sort(list);
				max=Tools.max(max, list.get(0).count);
			}
		}
		return max;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Loading           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads sequences from a FASTA file into a new SeqMap.
	 * Convenience method that reads the file and delegates to the main load method.
	 *
	 * @param ref Path to reference FASTA file
	 * @param k Packed kmer length in [1,32]
	 * @param mm Nonnegative middle-mask width; positive widths require k&gt;mm+1
	 * @param minCount Nonnegative explicit-header count threshold; see load(ArrayList,...)
	 * @param rcomp Whether to include reverse complements
	 * @param net Neural network for scoring sequences (may be null)
	 * @return New SeqMap containing the loaded sequences
	 */
	public static SeqMap load(String ref, int k, int mm, int minCount, boolean rcomp, CellNet net){
		ArrayList<Read> reads=FastaReadInputStream.toReads(ref, FileFormat.FASTA, -1);
		return load(reads, k, mm, minCount, rcomp, net);
	}

	/**
	 * Loads sequences from Read objects into a new SeqMap with k-mer indexing.
	 * Parses count and score information from read headers, applies neural network
	 * scoring if provided, and optionally includes reverse complements.
	 * A missing count defaults to one and bypasses filtering. The highest observed
	 * count, including rejected records, must reach minCount when assertions are enabled.
	 * A missing score or score=-1 is replaced by a network score when a network is supplied;
	 * this happens before count filtering. Forward arrays are borrowed; reverse complements
	 * are new arrays. Read objects, their IDs and the list are not retained or modified.
	 * The result is not sorted. Network scoring uses the supplied network on this thread.
	 *
	 * @param reads Non-null list of non-null reads with non-null IDs and score-compatible bases if net is used
	 * @param k Packed kmer length in [1,32]
	 * @param mm Nonnegative middle-mask width; positive widths require k&gt;mm+1
	 * @param minCount Minimum explicit header count; records without count= are retained with count 1
	 * @param rcomp Whether to include reverse complement sequences
	 * @param net Neural network for sequence scoring (may be null)
	 * @return New index of accepted sequences' defined kmers, without an implicit sort
	 */
	public static SeqMap load(ArrayList<Read> reads, int k, int mm, int minCount, boolean rcomp, CellNet net){
		SeqMap map=new SeqMap(k, mm, minCount);
		int maxCount=0;
		float[] vec=(net==null ? null : new float[net.numInputs()]);
		synchronized(map){
			for(Read r : reads){
				String id=r.id;
				int x=id.indexOf("count=");
				int count=(x>=0 ? Parse.parseInt(id, x+6) : 1);
				int y=id.indexOf("score=");
				float score=(y>=0 ? Parse.parseFloat(id, y+6) : -1);
				if((y<0 || score==-1) && net!=null){
					score=ScoreSequence.score(r.bases, vec, 0, net);
				}
				maxCount=Tools.max(count, maxCount);
				if(x<0 || count>=minCount){
					map.add(r.bases, count, score);
					if(rcomp){map.add(AminoAcid.reverseComplementBases(r.bases), count, score);}
				}
			}
		}
		assert(maxCount>=minCount) : "Error: minRefCount="+minCount+
		" was set too high for the reference, whose highest count was "+maxCount+
		"\nOnly set minRefCount if your ref headers have a 'count=X' term.";
		return map;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed kmer width, expected in [1,32]. */
	public final int k;
	/** Number of middle bases to mask in k-mers for fuzzy matching */
	public final int maskMiddle;
	/** Loader count threshold, retained as metadata; direct add() does not enforce it. */
	public final int minCount;
	/** Mask that replaces middle bases with an all-ones pattern via OR. */
	private final long midMask;
	/** Borrowed DNA symbol table: defined bases map to two-bit codes, undefined symbols to negatives. */
	private static final byte[] symbolToNumber=AminoAcid.baseToNumber;

	/**
	 * Reusable per-thread deduplication set; accepted candidate objects are still allocated.
	 */
	private final ThreadLocal<HashSet<SeqPosM>> localSet=new ThreadLocal<HashSet<SeqPosM>>(){
		/** Creates this thread's empty deduplication set.
		 * @return New set owned by this ThreadLocal/thread pair */
		@Override
		protected HashSet<SeqPosM> initialValue(){return new HashSet<SeqPosM>();}
	};
	/**
	 * Reusable per-thread result list; returned lists are borrowed until the next query.
	 */
	private final ThreadLocal<ArrayList<SeqPosM>> localList=new ThreadLocal<ArrayList<SeqPosM>>(){
		/** Creates this thread's reusable result list.
		 * @return New list owned by this ThreadLocal/thread pair */
		@Override
		protected ArrayList<SeqPosM> initialValue(){return new ArrayList<SeqPosM>();}
	};
	/**
	 * Thread-local storage for reusable SeqPosM object to avoid allocations during comparisons
	 */
	private final ThreadLocal<SeqPosM> localSP=new ThreadLocal<SeqPosM>(){
		/** Creates an uninitialized scratch key whose fields are replaced before lookup.
		 * @return New thread-confined scratch placement */
		@Override
		protected SeqPosM initialValue(){return new SeqPosM(null, 0);}
	};
	/** Atomic count of fetch/doubleFetch calls that passed null/short-input guards. */
	public final AtomicLong queries=new AtomicLong(0);
	/** Atomic count of candidate occurrences passing filters, counted before deduplication. */
	public final AtomicLong setQueries=new AtomicLong(0);

}
