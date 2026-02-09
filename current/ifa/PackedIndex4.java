package ifa;

import map.IntHashMap2;
import dna.AminoAcid;

/**
 * Hash-Backed Compressed Sparse Row Index for k-mers.
 * Stores k-mer positions in a single contiguous int array.
 * <p>
 * Key Features:
 * <ul>
 * <li><b>Implicit State Construction:</b> Uses the zero-initialized state of the array to detect list ends during the fill pass, avoiding sentinel passes.</li>
 * <li><b>Backwards Fill:</b> Fills lists from End to Start. The final map values naturally point to the list heads.</li>
 * <li><b>Stop Bit Encoding:</b> The sign bit (negative value) in the positions array indicates the end of a list.</li>
 * <li><b>Singleton Optimization:</b> K-mers appearing once are stored directly in the HashMap value (negative).</li>
 * </ul>
 * @author Brian Bushnell
 * @author Amber
 * @date Feb 5, 2026
 */
public class PackedIndex4{

	public PackedIndex4(byte[] ref, int k, int midMaskLen, int rStep){
		build(ref, k, midMaskLen, rStep);
	}

	/**
	 * Builds the index from a reference sequence.
	 * @param ref Reference sequence bytes
	 * @param k K-mer length
	 * @param midMaskLen Length of middle mask for hash reduction
	 * @param rStep Step size for reference indexing
	 */
	private void build(byte[] ref, int k, int midMaskLen, int rStep){
		final int len=ref.length;
		if(len<k){return;}
		
		// 1. Setup Map
		final int initialSize=(int)Math.min(1<<(2*Math.max(k-midMaskLen, 2)), len);
		map=new IntHashMap2(initialSize);
		
		final int shift=2*k, mask=~((-1)<<shift);
		final int stepMask=rStep-1, stepTarget=(k-1)&stepMask;
		final int midMask=Query.makeMidMask(k, midMaskLen);
		
		// 2. Count Pass (Forward)
		int kmer=0, clen=0;
		for(int i=0; i<len; i++){
			final byte b=ref[i];
			final int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){clen=0; kmer=0;}else{clen++;}
			if(clen>=k && ((i&stepMask)==stepTarget)){
				map.increment(kmer&midMask);
			}
		}

		// 3. Pack Pass (Linear Scan of Map)
		int[] keys=map.keys();
		int[] values=map.values();
		int totalHits=0;
		final int invalid=map.invalid();
		
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){
				int count=values[i];
				if(count==1){
					// Mark as Singleton (MAX_VALUE temp flag)
					values[i]=Integer.MAX_VALUE;
				}else{
					// Multi-hit: Store END INDEX.
					// We will fill backwards from here.
					totalHits+=count;
					values[i]=totalHits-1; 
				}
			}
		}
		
		// Allocate positions array (Zero-Initialized).
		positions=new int[totalHits];
		
		// 4. Fill Pass (Backwards)
		final int backShift=2*(k-1);
		kmer=0; clen=0;
		
		for(int j=len-k; j>=0; j--){
			final byte b=ref[j];
			final int x=AminoAcid.baseToNumber[b];
			
			if(x<0){
				clen=0; kmer=0; 
			}else{
				clen++;
				kmer=(kmer>>>2)|(x<<backShift);
			}
			
			int endPos=j+k-1;
			if(clen>=k && ((endPos&stepMask)==stepTarget)){
				int key=kmer&midMask;
				
				// Direct array access for speed
				int cell=map.findCell(key);
				int val=values[cell];
				int refPos=j; 
				
				if(val==Integer.MAX_VALUE){
					// Singleton: Encode (RefPos | MIN_VALUE)
					values[cell]=refPos|Integer.MIN_VALUE;
				}else if(val<0){
					// Already filled singleton
				}else{
					// Multi-Hit: 'val' is the current pointer.
					// Check if this slot is empty (0) to detect First Write vs Subsequent Write.
					
					if(positions[val]==0){
						// Case A: Slot is 0. This is the End Index (First Write).
						// Write Stop Bit. Do NOT decrement pointer.
						positions[val]=refPos|Integer.MIN_VALUE;
					}else{
						// Case B: Slot is non-zero. We have visited this list before.
						// Decrement pointer to new empty slot.
						val--;
						positions[val]=refPos; // Write positive RefPos
						values[cell]=val;      // Update map pointer
					}
				}
			}
		}
		
		// No Cleanup Pass needed. 
		// Map values for Multi-Hits now point to the Start Index (Head).
	}

	/** Returns value from map: -1 (Missing), < -1 (Singleton), or >= 0 (List Pointer) */
	public int get(int key){
		return map.get(key);
	}

	public IntHashMap2 map;
	public int[] positions;
}