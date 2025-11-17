package structures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicLong;

import dna.AminoAcid;
import fileIO.FileFormat;
import ml.CellNet;
import ml.ScoreSequence;
import shared.Parse;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;

/**
 * K-mer based sequence mapping structure with middle-masking support.
 * Extends LongArrayListHashMap to store sequences indexed by masked k-mers,
 * enabling efficient sequence similarity search and overlap detection.
 * @author Brian Bushnell
 */
public class SeqMap extends LongArrayListHashMap<SeqPos> {

	/** Creates a SeqMap with default parameters: k=11, maskMiddle=1, minCount=0 */
	public SeqMap(){this(11,1,0);}
	/**
	 * Creates a SeqMap with specified k-mer parameters and middle-masking.
	 * Computes bit mask to zero out middle bases for fuzzy k-mer matching.
	 *
	 * @param k_ K-mer length for indexing
	 * @param maskMiddle_ Number of middle bases to mask (0 disables masking)
	 * @param minCount_ Minimum count threshold for sequences
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
	
	public void add(byte[] s, int count, final float score) {
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
	
	public ArrayList<SeqPosM> fetch(byte[] query, final int a1, final int b1, final int minOverlap0, int maxMM,
			final int maxTrim, final int maxLopsidedness, float minOverlapFractionQ, boolean sort){
		final float qgc=Tools.calcGC(query, a1, b1);
		final int qlen=b1-a1+1;
		if(query==null || query.length<k || qlen<k){return null;}
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
//			System.err.println("i="+i);
			byte b=query[i];
			long x=symbolToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){
				len=0;
			}else{len++;}
			if(len>=k){
				long mmKmer=kmer|midMask;
				ArrayList<SeqPos> candidates=get(mmKmer);
//				if(candidates!=null) {
//					for(SeqPos sp : candidates) {
//						final float gcdif=Tools.absdif(qgc, sp.gc);
////						if(Tools.absdif(qgc, sp.gc)<0.24f) {//Intended to increase speed, but not very useful
//						final int rlen0=sp.seq().length, pos=sp.pos();
//						int a2=i-pos;
//						int b2=Tools.min(query.length-1, a2+rlen0-1);
//						a2=Tools.max(0, a2);
////						final int rlen=b2-a2+1;
//						final int overlap=Range.overlap(a1, b1, a2, b2);
//						final int lopsidedness=Range.lopsidedness(a1, b1, a2, b2);
////						final int minOverlapR=Tools.max(minOverlap0, (int)(minOverlapFractionR*rlen));
//						final int trim=Tools.max(a2-a1, b1-b2);
////						System.err.println("rlen="+rlen+", minOverlapR="+minOverlapR+", overlap="+overlap+", "+
////								new String(query, a2, rlen));
////						System.err.println("("+a1+"-"+b1+", "+a2+"-"+b2);
//
//						if(overlap>=minOverlapQ && trim<=maxTrim && lopsidedness<=maxLopsidedness
//								&& (gcdif*overlap<maxGCO)) {
////							assert(a2<=a1 && b2>=b1) : a1+"-"+b1+", "+a2+"-"+b2+"; trim="+trim+"; max="+maxTrim;
//							temp.setFrom(sp);
//							temp.setPos(i-pos);
//							lookups++;
//							if(!set.contains(temp)) {
//								SeqPosM clone=temp.clone();
//								set.add(clone);
//								retList.add(clone);
//							}
////							assert(a1==a2 && b1==b2) : "("+a1+"-"+b1+", "+a2+"-"+b2+"), "+
////								overlap+", "+minOverlapQ+", "+minOverlapR;
//						}
////						}
//					}
//				}
				if(candidates!=null) {//More concise.  Probably same speed.
					lookups+=addCandidates(candidates, query, a1, b1, qgc, 
							minOverlapQ, maxTrim, maxLopsidedness, maxGCO, i, temp, set, retList);
				}
			}
		}
		
//		ArrayList<SeqPosM> retList=null;
		if(!set.isEmpty()) {
			assert(set.size()==retList.size());
//			list=new ArrayList<SeqPosM>(set.size());
//			list.addAll(set);
			if(sort) {Collections.sort(retList);}
		}else {retList=null;}
		set.clear();
//		temp.seq=null;
		setQueries.addAndGet(lookups);
		return retList;
	}
	
	public ArrayList<SeqPosM> doubleFetch(byte[] query, final int a1, final int b1, final int a2, final int b2,
			final int minOverlap0, int maxMM,
			final int maxTrim, final int maxLopsidedness, float minOverlapFractionQ, boolean sort){
		final float qgc=Tools.calcGC(query, a1, b1);
		final float qgc2=Tools.calcGC(query, a2, b2);
		final int qlen=b1-a1+1;
		assert(qlen==b2-a2+1);
		if(query==null || query.length<k || qlen<k){return null;}
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
//			System.err.println("i="+i);
			byte b=query[i], c=query[j];
			long x=symbolToNumber[b], y=symbolToNumber[c];
			kmer=((kmer<<2)|x)&mask;
			kmer2=((kmer2<<2)|y)&mask;
			if(x<0){len=0;}else{len++;}
			if(y<0){len2=0;}else{len2++;}
			final long mmKmer=kmer|midMask, mmKmer2=kmer2|midMask;
			if(len>=k){
				ArrayList<SeqPos> candidates=get(mmKmer);
				if(candidates!=null) {
					lookups+=addCandidates(candidates, query, a1, b1, qgc, 
							minOverlapQ, maxTrim, maxLopsidedness, maxGCO, i, temp, set, retList);
				}
			}
			if(len2>=k && mmKmer!=mmKmer2){
				ArrayList<SeqPos> candidates=get(mmKmer2);
				if(candidates!=null) {//This a1 and b1 is intentional, to keep pos left-relative
					lookups+=addCandidates(candidates, query, a1, b1, qgc2, 
							minOverlapQ, maxTrim, maxLopsidedness, maxGCO, i, temp, set, retList);
				}
			}
		}
		
//		ArrayList<SeqPosM> retList=null;
		if(!set.isEmpty()) {
			assert(set.size()==retList.size());
//			list=new ArrayList<SeqPosM>(set.size());
//			list.addAll(set);
			if(sort) {Collections.sort(retList);}
		}else {retList=null;}
		set.clear();
//		temp.seq=null;
		setQueries.addAndGet(lookups);
		return retList;
	}
	
	private int addCandidates(ArrayList<SeqPos> candidates, byte[] query, final int a1, final int b1, float qgc,
			int minOverlapQ, int maxTrim, int maxLopsidedness, float maxGCO, int qpos, 
			SeqPosM temp, HashSet<SeqPosM> set, ArrayList<SeqPosM> retList) {
		int lookups=0;
		for(SeqPos sp : candidates) {
			final float gcdif=Tools.absdif(qgc, sp.gc);
			final int rlen0=sp.seq().length, pos=sp.pos();
			int a3=qpos-pos;
			int b3=Tools.min(query.length-1, a3+rlen0-1);
			a3=Tools.max(0, a3);
			final int overlap=Range.overlap(a1, b1, a3, b3);
			final int lopsidedness=Range.lopsidedness(a1, b1, a3, b3);
			final int trim=Tools.max(a3-a1, b1-b3);

			if(overlap>=minOverlapQ && trim<=maxTrim && lopsidedness<=maxLopsidedness
					&& (gcdif*overlap<maxGCO)) {
				temp.setFrom(sp);
				temp.setPos(qpos-pos);
				lookups++;
				if(!set.contains(temp)) {
					SeqPosM clone=temp.clone();
					set.add(clone);
					retList.add(clone);
				}
			}
		}
		return lookups;
	}
	
	/**
	 * Sorts all sequence lists by their natural ordering and finds maximum count.
	 * Iterates through all stored sequence lists, sorts each by count/score,
	 * and returns the highest count value found.
	 * @return Maximum count value across all stored sequences
	 */
	public int sort() {
		int max=0;
		for(ArrayList<SeqPos> list : values()) {
			if(list!=null && !list.isEmpty()) {
				Collections.sort(list);
				max=Tools.max(max, list.get(0).count);
			}
//			assert(false) : max;
		}
		return max;
	}
	
	/**
	 * Loads sequences from a FASTA file into a new SeqMap.
	 * Convenience method that reads the file and delegates to the main load method.
	 *
	 * @param ref Path to reference FASTA file
	 * @param k K-mer length for indexing
	 * @param mm Middle masking length
	 * @param minCount Minimum count threshold for sequences
	 * @param rcomp Whether to include reverse complements
	 * @param net Neural network for scoring sequences (may be null)
	 * @return New SeqMap containing the loaded sequences
	 */
	public static SeqMap load(String ref, int k, int mm, int minCount, boolean rcomp, CellNet net) {
		ArrayList<Read> reads=FastaReadInputStream.toReads(ref, FileFormat.FASTA, -1);
		return load(reads, k, mm, minCount, rcomp, net);
	}
	
	/**
	 * Loads sequences from Read objects into a new SeqMap with k-mer indexing.
	 * Parses count and score information from read headers, applies neural network
	 * scoring if provided, and optionally includes reverse complements.
	 *
	 * @param reads List of Read objects containing sequences and metadata
	 * @param k K-mer length for indexing
	 * @param mm Middle masking length (number of middle bases to ignore)
	 * @param minCount Minimum count threshold - only sequences with count >= minCount are added
	 * @param rcomp Whether to include reverse complement sequences
	 * @param net Neural network for sequence scoring (may be null)
	 * @return New SeqMap containing indexed sequences meeting the count threshold
	 */
	public static SeqMap load(ArrayList<Read> reads, int k, int mm, int minCount, boolean rcomp, CellNet net) {
		SeqMap map=new SeqMap(k, mm, minCount);
		int maxCount=0;
		float[] vec=(net==null ? null : new float[net.numInputs()]);
		synchronized(map) {
			for(Read r : reads) {
				String id=r.id;
				int x=id.indexOf("count=");
				int count=(x>=0 ? Parse.parseInt(id, x+6) : 1);
				int y=id.indexOf("score=");
				float score=(y>=0 ? Parse.parseFloat(id, y+6) : -1);
				if((y<0 || score==-1) && net!=null) {
					score=ScoreSequence.score(r.bases, vec, 0, net);
				}
				maxCount=Tools.max(count, maxCount);
				if(x<0 || count>=minCount) {
					map.add(r.bases, count, score);
					if(rcomp) {map.add(AminoAcid.reverseComplementBases(r.bases), count, score);}
				}
			}
		}
		assert(maxCount>=minCount) : "Error: minRefCount="+minCount+
		" was set too high for the reference, whose highest count was "+maxCount+
		"\nOnly set minRefCount if your ref headers have a 'count=X' term.";
		return map;
	}
	
	/** K-mer length used for sequence indexing */
	public final int k;
	/** Number of middle bases to mask in k-mers for fuzzy matching */
	public final int maskMiddle;
	/** Minimum count threshold for sequences to be included */
	public final int minCount;
	/** Bit mask for zeroing out middle bases in k-mers */
	private final long midMask;
	/** Lookup table for converting DNA bases to numeric values */
	private static final byte[] symbolToNumber=AminoAcid.baseToNumber;
	
	/**
	 * Thread-local storage for reusable HashSet to avoid allocations during searches
	 */
	private final ThreadLocal<HashSet<SeqPosM>> localSet=new ThreadLocal<HashSet<SeqPosM>>(){
        @Override protected HashSet<SeqPosM> initialValue() {return new HashSet<SeqPosM>();}
    };
	/**
	 * Thread-local storage for reusable ArrayList to avoid allocations during searches
	 */
	private final ThreadLocal<ArrayList<SeqPosM>> localList=new ThreadLocal<ArrayList<SeqPosM>>(){
        @Override protected ArrayList<SeqPosM> initialValue() {return new ArrayList<SeqPosM>();}
    };
	/**
	 * Thread-local storage for reusable SeqPosM object to avoid allocations during comparisons
	 */
	private final ThreadLocal<SeqPosM> localSP=new ThreadLocal<SeqPosM>(){
        @Override protected SeqPosM initialValue() {return new SeqPosM(null, 0);}
    };
	/** Counter for total number of search queries performed */
	public final AtomicLong queries=new AtomicLong(0);
	/** Counter for total number of k-mer lookups performed during searches */
	public final AtomicLong setQueries=new AtomicLong(0);
	
}
