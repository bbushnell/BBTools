package ifa;

import map.IntLongHashMap2;
import dna.AminoAcid;

/**
 * Forward masked k-mer index with packed counts and contiguous position lists.
 * Missing keys return {@code -1L}. A singleton stores its zero-based reference
 * start in the high 32 bits and 1 in the low bits. For repeated keys, the high
 * bits store an offset into {@link #positions} and the low bits store the count.
 * Each repeated key's positions are ascending. Masked keys are not canonicalized
 * against their reverse complements. Ambiguous bases invalidate the whole window,
 * even when the ambiguous position would have been masked.
 * <p>
 * Construction reads but does not retain or change the reference. The public map
 * and array are mutable; this class provides no synchronization.
 * @author Brian Bushnell
 * @contributor Amber
 * @date February 4, 2026
 */
public class PackedIndex{

	/**
	 * Builds an index of sampled reference windows.
	 * References shorter than k produce an empty map and positions array.
	 * @param ref Non-null ASCII reference bases
	 * @param k K-mer length, expected to be 1 through 15
	 * @param midMaskLen Number of central bases masked by Query.makeMidMask; 0 disables masking
	 *        and positive values require {@code k > midMaskLen + 1}
	 * @param rStep Positive power-of-two stride; sampled window starts are multiples of it
	 */
	public PackedIndex(byte[] ref, int k, int midMaskLen, int rStep){
		build(ref, k, midMaskLen, rStep);
	}

	/** Counts sampled keys, allocates repeated-key slices, then fills them in reference order. */
	private void build(byte[] ref, int k, int midMaskLen, int rStep){
		final int len=ref.length;
		if(len<k){
			//FIXED IFA-002: IFA3 prescan reads map directly, so even a seedless reference needs an empty map.
			map=new IntLongHashMap2(1, 0.7);
			positions=new int[0];
			return;
		}

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
					// RefPos is nonnegative; a singleton at position 0 packs as 1L.
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

	/**
	 * Looks up an already masked forward key without modifying it.
	 * @param key Two-bit encoded key with the same middle mask used at construction
	 * @return -1L if absent; otherwise high bits hold a singleton position or list offset,
	 *         and low bits hold the count (1 identifies the singleton form)
	 */
	public long get(int key){
		return map.get(key);
	}

	/** Owned map of packed values; empty for references shorter than k. */
	public IntLongHashMap2 map;
	/** Owned storage for repeated keys only; each slice is ascending; empty for a short reference. */
	public int[] positions;
}
