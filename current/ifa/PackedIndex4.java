package ifa;

import map.IntHashMap2;
import dna.AminoAcid;

/**
 * Forward masked k-mer index with singleton values and sign-terminated position lists.
 * A missing key returns -1; a singleton's value encodes its zero-based reference
 * start with the sign bit set. A nonnegative value is a list head in
 * {@link #positions}; the last position in that list carries the sign bit.
 * Other list positions are unflagged. Clear the sign bit to recover a position.
 * Forward keys are not canonicalized against reverse complements. Ambiguous
 * bases invalidate windows even at masked positions.
 * <p>
 * Key Features:
 * <ul>
 * <li><b>Implicit State Construction:</b> Uses the zero-initialized state of the array to detect list ends during the fill pass, avoiding sentinel passes.</li>
 * <li><b>Backwards Fill:</b> Fills lists from End to Start. The final map values naturally point to the list heads.</li>
 * <li><b>Stop Bit Encoding:</b> The sign bit (negative value) in the positions array indicates the end of a list.</li>
 * <li><b>Singleton Optimization:</b> K-mers appearing once are stored directly in the HashMap value (negative).</li>
 * </ul>
 * Construction reads but does not retain or change the reference. The public map
 * and array are mutable and no synchronization is provided. Position lists
 * are filled backward so their final order is ascending.
 * @author Brian Bushnell
 * @author Amber
 * @date Feb 5, 2026
 */
public class PackedIndex4{

	/**
	 * Builds the count map and packed position storage.
	 * References shorter than k produce an empty map and positions array.
	 * @param ref Non-null ASCII reference bases
	 * @param k K-mer length, expected to be 1 through 15
	 * @param midMaskLen Central bases masked by Query.makeMidMask; 0 disables masking
	 *        and positive values require {@code k > midMaskLen + 1}
	 * @param rStep Positive power-of-two stride; sampled starts are multiples of it
	 */
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
		if(len<k){
			//FIXED IFA-002: IFA4 prescan reads map directly, so even a seedless reference needs an empty map.
			map=new IntHashMap2(1);
			positions=new int[0];
			return;
		}
		
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
		
		//FIXED IFA-001: scan from the last base so the final k-1 bases prime the reverse rolling kmer.
		//Starting at len-k leaves counted tail hits unfilled; IFA4's decoders require every count to be filled.
		for(int j=len-1; j>=0; j--){
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
					
					//The positions[val]==0 "empty slot" test is safe even though refPos can be 0: refPos==0 only when j==0,
					//the globally-lowest j, which is processed last (backwards scan) and thus is always a key's terminal HEAD
					//write — no later write re-reads that slot. Non-head slots always get a j>0 (nonzero) refPos.
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

	/**
	 * Looks up an already masked forward key.
	 * @param key Two-bit encoded key with the construction-time middle mask
	 * @return -1 if absent, less than -1 for a sign-flagged singleton position,
	 *         or a nonnegative head offset for a sign-terminated positions list
	 */
	public int get(int key){
		return map.get(key);
	}

	/** Owned packed-value map; empty for references shorter than k. */
	public IntHashMap2 map;
	/** Owned multi-hit storage; sign bit marks a slice's last position; empty for a short reference. */
	public int[] positions;
}
