package prok;

import dna.AminoAcid;
import map.LongHashSet;
import map.LongIntMap;
import structures.IntList;

/**
 * Shared 17-mer front end for conserved-RNA scavengers.
 *
 * <p>Each registered slot owns one bit in an int mask.  A seed present in
 * several families retains every bit, so one oriented-sequence scan fans the
 * hit out to every downstream family without changing family-local window,
 * shortlist, alignment, or acceptance logic.</p>
 */
final class ConservedRnaSeedIndex {

	static final int K=17;
	static final int MAX_SLOTS=32;

	ConservedRnaSeedIndex(LongHashSet[] seedSets){
		if(seedSets==null){throw new IllegalArgumentException("Conserved-RNA seed sets are null");}
		if(seedSets.length>MAX_SLOTS){
			throw new IllegalArgumentException("Conserved-RNA 17-mer index has "+seedSets.length
				+" family/type slots; the int bit-vector capacity is "+MAX_SLOTS);
		}
		slotCount=seedSets.length;
		long seedUpperBound=0;
		for(int slot=0; slot<slotCount; slot++){
			LongHashSet set=seedSets[slot];
			if(set==null || set.size()<1){
				throw new IllegalArgumentException("Conserved-RNA 17-mer slot "+slot+" has no seeds");
			}
			seedUpperBound+=set.size();
		}
		final int initialSize=(int)Math.min(1_000_000_000L, Math.max(16L, seedUpperBound*2L));
		seedToTypeBits=new LongIntMap(initialSize);
		for(int slot=0; slot<slotCount; slot++){
			final int bit=1<<slot;
			for(long seed : seedSets[slot].toArray()){
				final int old=seedToTypeBits.getOrDefault(seed, 0);
				seedToTypeBits.put(seed, old|bit);
			}
		}
	}

	ScanResult scan(byte[] bases){
		if(bases==null || bases.length<K){return new ScanResult(new int[slotCount][]);}
		final IntList[] lists=new IntList[slotCount];
		final long mask=~((-1L)<<(2*K));
		final byte[] bton=AminoAcid.baseToNumber;
		long kmer=0;
		int len=0;
		for(int i=0; i<bases.length; i++){
			final int x=bton[bases[i]];
			if(x<0){
				len=0;
				kmer=0;
				continue;
			}
			kmer=((kmer<<2)|x)&mask;
			len++;
			if(len<K){continue;}
			int typeBits=seedToTypeBits.getOrDefault(kmer, 0);
			while(typeBits!=0){
				final int slot=Integer.numberOfTrailingZeros(typeBits);
				IntList list=lists[slot];
				if(list==null){lists[slot]=list=new IntList();}
				list.add(i-K/2);
				typeBits&=typeBits-1;
			}
		}
		final int[][] hits=new int[slotCount][];
		for(int slot=0; slot<slotCount; slot++){
			hits[slot]=(lists[slot]==null ? EMPTY : lists[slot].toArray());
		}
		return new ScanResult(hits);
	}

	int slotCount(){return slotCount;}

	static final class ScanResult {
		ScanResult(int[][] hits_){hits=hits_;}
		int[] hits(int slot){
			if(slot<0 || slot>=hits.length){
				throw new IllegalArgumentException("Invalid conserved-RNA seed slot "+slot+" for "+hits.length+" slots");
			}
			return hits[slot]==null ? EMPTY : hits[slot];
		}
		private final int[][] hits;
	}

	private final LongIntMap seedToTypeBits;
	private final int slotCount;
	private static final int[] EMPTY=new int[0];
}
