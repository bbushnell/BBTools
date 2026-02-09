package ifa;

import map.IntLongHashMap2;
import dna.AminoAcid;

/**
 * Hash-Backed Compressed Sparse Row Index.
 * Stores k-mer positions in a single contiguous array (positions).
 * * Singleton Optimization:
 * K-mers appearing only once store their position directly in the HashMap value,
 * bypassing the positions array entirely.
 * * Value Layout:
 * - Missing: -1L
 * - Singleton: (RefPos << 32) | 1
 * - Multi-Hit: (Offset << 32) | Count (where Count > 1)
 * * @author Brian Bushnell
 * @contributor Amber
 * @date February 4, 2026
 */
public class PackedIndex {

	public PackedIndex(byte[] ref, int k, int midMaskLen, int rStep){
		build(ref, k, midMaskLen, rStep);
	}

	private void build(byte[] ref, int k, int midMaskLen, int rStep){
		final int len=ref.length;
		if(len<k){return;}

		// 1. Estimation & Allocation
		final int defined=Math.max(k-midMaskLen, 2);
		final int kSpace=(1<<(2*defined));
		final int initialSize=(int)Math.min(kSpace, len);
		map=new IntLongHashMap2(initialSize, 0.7);

		final int shift=2*k, mask=~((-1)<<shift);
		final int stepMask=rStep-1;
		final int stepTarget=(k-1)&stepMask;
		final int midMask=Query.makeMidMask(k, midMaskLen);

		// 2. Count Pass
		int kmer=0, clen=0;
		for(int i=0; i<len; i++){
			final byte b=ref[i];
			final int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){clen=0; kmer=0;}else{clen++;}
			if(clen>=k && ((i&stepMask)==stepTarget)){
				int key=(kmer&midMask);
				map.increment(key);
			}
		}

		// 3. Prefix Sum (Packing)
		int[] keys=map.keys();
		long[] values=map.values();
		long totalHits=0;
		int invalid=map.invalid();

		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){
				long count=values[i];
				if(count==1){
					// Flag as singleton. 
					// We use -1L because no valid packed value (offset<<32 | count) 
					// can ever be -1 (since count must be > 0).
					values[i]=-1L; 
				}else{
					// Multi-hit: store offset in high bits, reset count to 0.
					values[i]=(totalHits<<32); 
					totalHits+=count;
				}
			}
		}

		positions=new int[(int)totalHits];

		// 4. Fill Pass
		kmer=0; clen=0;
		for(int i=0; i<len; i++){
			final byte b=ref[i];
			final int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){clen=0; kmer=0;}else{clen++;}
			if(clen>=k && ((i&stepMask)==stepTarget)){
				int key=(kmer&midMask);

				// Lookup packed state
				long packed=map.get(key);

				if(packed==-1L){
					// Singleton Case:
					// Store (RefPos << 32) | 1
					// Note: RefPos is guaranteed positive.
					long val=((long)(i-k+1)<<32) | 1L;
					map.set(key, val);
				}else{
					// Multi-Hit Case:
					// Standard CSR fill
					int offset=(int)(packed>>>32);
					int count=(int)packed;

					positions[offset+count]=i-k+1;
					map.set(key, packed+1);
				}
			}
		}
	}

	/** Returns packed (offset << 32) | count, or -1 if missing */
	public long get(int key){
		return map.get(key);
	}

	public IntLongHashMap2 map;
	public int[] positions;
}