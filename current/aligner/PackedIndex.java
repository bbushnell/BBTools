package aligner;

import map.IntLongHashMap2;
import dna.AminoAcid;

/**
 * Hash-Backed Compressed Sparse Row Index.
 * Stores k-mer positions in a single contiguous array (positions)
 * and uses a hash map to map k-mers to (offset, count) pairs packed into a long.
 * Drastically reduces memory overhead and pointer chasing compared to IntListHashMap.
 * @author Brian Bushnell
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
		map=new IntLongHashMap2(initialSize, 0.7f);

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
				values[i]=(totalHits<<32); // Store offset in high bits, reset count to 0 in low bits
				totalHits+=count;
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
				long packed=map.increment(key);
				int offset=(int)(packed>>>32);
				int count=(int)packed-1;
				assert(count>=0);
				positions[offset+count]=i-k+1;
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