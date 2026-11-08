package bloom;

import java.util.Arrays;
import java.util.BitSet;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/**
 * Implements long-kmer correction primitives using multi-word {@link Kmer} state.
 * Normally used for k&gt;31. Each stored window holds its canonical xor() hash;
 * dual-hash mode also maintains xor2() in a parallel thread-local list. A hash
 * cannot reconstruct the original sequence, so candidate generation rebuilds
 * a rolling Kmer from bases. This implementation was ported from assemble.Tadpole2;
 * the separate short-kmer implementation retains packed sequence values instead.
 *
 * Build and query must use the same effective length after {@link Kmer#getKbig}
 * normalization and must keep Kmer.MASK_CORE false. The normal wrapper establishes
 * these settings; direct callers must do so themselves. A second, smaller window
 * size is not supported. Counts retain the table's collision and saturation limits.
 *
 * The filter is borrowed. Superclass entry points initialize thread-local scratch;
 * direct use of primitives that need scratch requires that initialization first.
 * Mutable sequence/count buffers belong to the caller. Dual-hash stored-position
 * queries must use the primary list paired with the current thread's secondary
 * list; a later fill replaces that pairing. Concurrent callers need separate
 * buffers and stable filter/global settings. Independent operations on the same
 * thread share this scratch and secondary-list state.
 *
 * @author Brian Bushnell, Nepgear
 * @date August 19, 2026
 */
public class BloomFilterCorrector2 extends BloomFilterCorrector{

	/**
	 * Stores the borrowed table and checks the long-kmer geometry assertions.
	 * Does not initialize thread-local scratch or normalize the requested length.
	 * @param filter_ Borrowed count table, required before a lookup
	 * @param k_ Effective k-mer length, already normalized for Kmer storage
	 * @param ksmall_ Must equal k_; smaller constituent-window counting is unsupported
	 */
	public BloomFilterCorrector2(BloomFilter filter_, int k_, int ksmall_){
		super(filter_, k_, ksmall_);
		assert(k_==ksmall_) : "BFC2 does not support the ksmall sliding-window trick; k must equal ksmall for true k>31 mode. k="+k_+", ksmall="+ksmall_;
		assert(!Kmer.MASK_CORE) : "MASK_CORE must be false for build/query hash agreement (plan's crux pin); was true.";
	}

	/**
	 * Diagnostic-only hook: {perWordK, mult}, or null in production. When set,
	 * forces every Kmer this instance constructs to use this explicit split
	 * instead of the natural Kmer.getK/getMult(k) derivation -- lets a test
	 * harness exercise the multi-word carry/combine path (addLeft/addRight
	 * across words, xor's rotate+combine) at a k value BFC1 can also run
	 * (k&lt;=31), for direct side-by-side comparison isolated from any
	 * k-length confound. newKmer reads this field when scratch is initialized;
	 * normal operation leaves it null. Tests must not change it while a caller
	 * uses this class, and must restore it afterward.
	 */
	public static int[] FORCE_MULT_SPLIT=null;

	/** Constructs scratch geometry from the diagnostic split, or the normal effective k. */
	@Override
	protected Kmer newKmer(){
		return FORCE_MULT_SPLIT!=null ? new Kmer(FORCE_MULT_SPLIT[0], FORCE_MULT_SPLIT[1]) : super.newKmer();
	}

	/*--------------------------------------------------------------*/
	/*----------------      Abstract Seam: Fill      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Replaces the primary window hashes and, in dual mode, this thread's secondary list.
	 * For reads of at least k bases, writes length-k+1 entries indexed by window
	 * start; undefined windows receive -1 in both lists. Shorter reads return zero
	 * before clearing either list. Uses initialized localKmer scratch, not bases as
	 * mutable workspace.
	 *
	 * Historical port failure: recording before the first complete window produced
	 * one entry per base, desynchronized counts and crashed smooth's edge lookup.
	 * Retain window-start indexing; it is shared by all inherited correction paths.
	 * @param bases Non-null sequence, read without modification
	 * @param kmers Caller-owned primary list to replace when the read is long enough
	 * @return Number of valid windows, excluding undefined-window entries
	 */
	@Override
	public int fillKmers(byte[] bases, LongList kmers){
		final int blen=bases.length;
		if(blen<k){return 0;}
		final int min=k-1;
		kmers.clear();
		final boolean dualHash=dualHash();
		final LongList kmers2=(dualHash ? xor2List() : null);
		if(dualHash){kmers2.clear();}
		final Kmer kmer=localKmer.get();
		kmer.clear();
		int valid=0;
		for(int i=0; i<blen; i++){
			kmer.addRight(bases[i]);
			if(i>=min){
				if(kmer.len>=k){
					kmers.add(kmer.xor());
					if(dualHash){kmers2.add(kmer.xor2());}
					valid++;
				}else{
					kmers.add(-1);
					if(dualHash){kmers2.add(-1);}
				}
			}
		}
		return valid;
	}

	/**
	 * Rebuilds the unchanged anchor starting at a, then refreshes windows after a substitution at a+k.
	 * Updates up to k following starts, clipped to counts.size, in the primary list,
	 * count list and current secondary list. Undefined windows receive -1 hashes
	 * and zero counts. Requires initialized scratch and aligned lists.
	 * @param bases Non-null modified sequence, read without modification
	 * @param kmers Primary window list to update in place
	 * @param counts Full count list aligned with kmers, updated in place
	 * @param a In-range start of the unchanged anchor immediately preceding the changed base
	 */
	@Override
	public void regenerateKmers(byte[] bases, LongList kmers, IntList counts, final int a){
		final int loc=a+k;
		final int lim=Tools.min(counts.size, a+k+1);
		final Kmer kmer=localKmer.get();
		final boolean dualHash=dualHash();
		final LongList kmers2=(dualHash ? xor2List() : null);
		kmer.clear();
		for(int i=a; i<a+k; i++){kmer.addRight(bases[i]);}

		for(int i=loc, j=a+1; j<lim; i++, j++){
			kmer.addRight(bases[i]);
			if(kmer.len>=k){
				final long xor=kmer.xor();
				final long xor2=(dualHash ? kmer.xor2() : 0);
				kmers.set(j, xor);
				if(dualHash){kmers2.set(j, xor2);}
				final int count=(dualHash ? filter.filter.read(xor, xor2) : filter.filter.read(xor));
				counts.set(j, count);
			}else{
				kmers.set(j, -1);
				if(dualHash){kmers2.set(j, -1);}
				counts.set(j, 0);
			}
		}
	}

	/**
	 * Recounts windows spanning the first through last changed base, including intervening gaps.
	 * Uses initialized localKmer scratch and writes zero for undefined windows.
	 * When smoothing is enabled, suppresses strict interior peaks across the whole
	 * count list afterward. Does not refresh stored hash lists or modify bases.
	 * @param bases Non-null modified sequence
	 * @param counts Existing window-count list, updated in place
	 * @param dummy Unused compatibility argument; may be null
	 * @param changed Nonempty set of in-range modified base indices, read without modification
	 * @return Number of valid windows recounted before smoothing
	 */
	@Override
	public int regenerateCounts(byte[] bases, IntList counts, final Kmer dummy, BitSet changed){
		assert(!changed.isEmpty());
		final int firstBase=changed.nextSetBit(0), lastBase=changed.length()-1;
		final int start=Tools.max(0, firstBase-k+1);
		final int lim=Tools.min(lastBase+k-1, bases.length-1);

		final Kmer kmer=localKmer.get();
		kmer.clear();
		int valid=0;

		for(int i=start; i<=lim; i++){
			kmer.addRight(bases[i]);
			final int c=i-k+1;
			if(i>=firstBase){
				if(kmer.len>=k){
					valid++;
					final int count=getCount(kmer);
					counts.set(c, count);
				}else if(c>=0){
					counts.set(c, 0);
				}
			}
		}

		if(smooth){
			int[] array=counts.array;
			for(int i=1; i<counts.size-1; i++){
				int a2=array[i-1], b=array[i], c2=array[i+1];
				if(b>a2 && b>c2){array[i]=Tools.max(a2, c2);}
			}
		}

		return valid;
	}

	/**
	 * Regenerates k-mer counts for a specific region after base changes.
	 * Mirrors BFC1's window (starts at ca+k-1-k+1=ca, extends to lim), but rolls
	 * a live Kmer instead of a bit-packed long. Must use localKmer2, NOT
	 * localKmer -- this is called mid-loop from reassemble_inner, which keeps
	 * its own rolling window live in localKmer across the whole scan; clearing
	 * that shared object here (as a prior version did) corrupted the caller's
	 * position tracking for every base after the first correction in a read,
	 * producing systematic under-correction (root cause of the BFC2-vs-BFC1
	 * divergence found via TraceRead2.java: kmer.len desynced from the loop's
	 * position `a` immediately after the first regenerateCounts call).
	 * Recounts the anchor at ca and up to k following windows, clipped at the read
	 * end; undefined windows receive zero. Does not smooth or refresh stored hashes.
	 * @param bases Non-null modified sequence, read without modification
	 * @param counts Existing full count list, updated without changing its size
	 * @param ca In-range anchor-window start preceding the changed base at ca+k
	 * @return Number of valid windows recounted
	 */
	@Override
	public int regenerateCounts(byte[] bases, IntList counts, final int ca){
		final int b=ca+k-1;
		final int lim=Tools.min(bases.length, b+k+1);
		final int start=b-k+1;

		final Kmer kmer=localKmer2.get();
		kmer.clear();
		int valid=0;
		final int clen=counts.size;

		for(int i=start; i<lim; i++){
			kmer.addRight(bases[i]);
			if(i>=b){
				if(kmer.len>=k){
					valid++;
					final int count=getCount(kmer);
					counts.set(i-k+1, count);
				}else{
					counts.set(i-k+1, 0);
				}
			}
		}

		assert((counts.size==0 && bases.length<k) || counts.size==bases.length-k+1) : bases.length+", "+k+", "+counts.size;
		assert(clen==counts.size);

		return valid;
	}

	/*--------------------------------------------------------------*/
	/*----------------   Abstract Seam: Reassemble   ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Scans left to right and attempts substitutions at count drops.
	 * Changes bb's bases and refreshes counts in place; does not extend the builder
	 * or modify qualities. The rolling state uses localKmer, while candidate and
	 * recount operations use localKmer2 so they cannot reset the active scan.
	 * Candidate identity is compared as a base code, independently of its support count (BLM203).
	 * @param bb Mutable read builder with backing capacity equal to logical length; recounting uses array length
	 * @param quals Read-only qualities aligned to bases, or null to use quality 20
	 * @param rightCounts Scratch array with exactly four entries, overwritten
	 * @param counts Full window-count list aligned to bb, updated in place
	 * @param errorExtension Number of preceding windows used for the similarity check
	 * @return Number of corrections recorded by this directional scan
	 */
	@Override
	public int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension){
		final int length=bb.length();
		if(length<k+1+deadZone){return 0;}
		final byte[] bases=bb.array;

		final Kmer kmer=localKmer.get();
		kmer.clear();

		int detected=0;
		int corrected=0;

		for(int a=0, lim=length-deadZone-1; a<lim; a++){
			final byte aBase=bases[a];
			kmer.addRight(aBase);

			if(kmer.len>=k){
				final int b=a+1;
				final int ca=a-k+1;
				final int cb=ca+1;

				final int aCount=counts.get(ca);
				final int bCount=counts.get(cb);
				final byte qb=(quals==null ? 20 : quals[b]);

				if(verbose){
					System.err.println("BFC2 ca="+ca+", cb="+cb+"; aCount="+aCount+", bCount="+bCount);
					System.err.println(isError(aCount, bCount, qb)+", "+isSimilar(aCount, ca-errorExtension, ca-1, counts)+
							", "+isError(aCount, ca+2, ca+k, counts));
				}

				if(isSubstitution(ca, errorExtension, qb, counts)){
					if(verbose){System.err.println("BFC2 ***Found error: "+aCount+", "+bCount);}
					final int rightMaxPos=fillRightCounts(kmer, rightCounts);
					final int rightMax=rightCounts[rightMaxPos];
					final int rightSecondPos=Tools.secondHighestPosition(rightCounts);
					final int rightSecond=rightCounts[rightSecondPos];

					final byte base=bases[b];
					final byte num=AminoAcid.baseToNumber[base];

					if(verbose){
						System.err.println("BFC2 kmer.xor()="+kmer.xor()+", selfCount="+getCount(kmer)+" vs aCount="+aCount+", kmer.len="+kmer.len);
						System.err.println("BFC2 Counts: "+aCount+", "+Arrays.toString(rightCounts));
						System.err.println("BFC2 rightMaxPos="+rightMaxPos+", rightMax="+rightMax+", rightSecondPos="+rightSecondPos+", rightSecond="+rightSecond);
						System.err.println("BFC2 base="+(char)base+", num="+num+"; num==rightMaxPos? "+(num==rightMaxPos));
					}

					if(rightMax>=minCountExtend){
						detected++;
						//BLM203: compare base codes, not a base code with its support count.
						//An unchanged winning base is not a correction.
						if(num==rightMaxPos){
							detected--;
							if(verbose){System.err.println("BFC2 rejected: num==rightMaxPos");}
						}else if((isError(rightMax, rightSecond, qb) || !isJunction(rightMax, rightSecond)) && isSimilar(aCount, rightMax)){
							bases[b]=AminoAcid.numberToBase[rightMaxPos];
							corrected++;
							regenerateCounts(bases, counts, ca);
							if(verbose){System.err.println("BFC2 Corrected error: "+num+"->"+rightMaxPos+". New counts:\n"+counts);}
							//`kmer` currently ends at position a; bases[b]=bases[a+1] hasn't been
							//rolled in yet, so no correction to the live rolling state is needed here.
						}else if(verbose){
							System.err.println("BFC2 rejected: isError(rightMax,rightSecond,qb)="+isError(rightMax, rightSecond, qb)+
									", isJunction="+isJunction(rightMax, rightSecond)+", isSimilar(aCount,rightMax)="+isSimilar(aCount, rightMax));
						}
					}else if(verbose){
						System.err.println("BFC2 rejected: rightMax<minCountExtend: "+rightMax+"<"+minCountExtend);
					}
				}
			}
		}

		return corrected;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Abstract Seam: Extend      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Derives the trailing window from bb and attempts rightward extension.
	 * Returns zero for a short builder or undefined trailing window. Requires
	 * initialized scratch; appends to bb without changing the count table.
	 * @param bb Caller-owned sequence builder to extend
	 * @param leftCounts Four-entry left-branch scratch array, or null to omit those checks
	 * @param rightCounts Scratch array with exactly four entries
	 * @param distance Maximum number of bases to append
	 * @param includeJunctionBase Allow the terminal junction base when the forward strand is canonical
	 * @return Number of bases appended
	 */
	@Override
	public int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase){
		final int initialLength=bb.length();
		if(initialLength<k){return 0;}
		final byte[] bases=bb.array;

		final Kmer kmer=localKmer.get();
		kmer.clear();
		for(int i=initialLength-k; i<initialLength; i++){
			kmer.addRight(bases[i]);
		}
		if(kmer.len<k){return 0;}

		return extendCore(bb, leftCounts, rightCounts, distance, includeJunctionBase, kmer);
	}

	/**
	 * Extends rightward starting from the kmer at list position originPos,
	 * derived fresh from bases (the stored xor hash cannot be un-hashed back
	 * into a rollable kmer). If fromRcomp, flips to the reverse-complement
	 * orientation first via Kmer.rcomp(), mirroring BFC1's rkmer-as-starting-
	 * point trick for the pincer's right-side anchor.
	 * @param bb Caller-owned builder to append to
	 * @param leftCounts Four-entry left-branch scratch array, or null
	 * @param rightCounts Scratch array with exactly four entries
	 * @param distance Maximum number of bases to append
	 * @param includeJunctionBase Allow the terminal junction base on the canonical forward strand
	 * @param bases Sequence containing a valid anchor at originPos
	 * @param kmers Unused hash list; the anchor must be reconstructed from bases
	 * @param originPos Start base index of the anchor window
	 * @param fromRcomp Whether to reverse-complement the anchor before extending
	 * @return Number of bases appended
	 */
	@Override
	protected int extendFromStored(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance,
			boolean includeJunctionBase, final byte[] bases, final LongList kmers, final int originPos, final boolean fromRcomp){
		final Kmer kmer=localKmer.get();
		kmer.clear();
		for(int i=originPos; i<originPos+k; i++){
			kmer.addRight(bases[i]);
		}
		assert(kmer.len==k) : "expected a fully valid stored kmer at originPos="+originPos+"; len="+kmer.len+", k="+k;
		if(fromRcomp){kmer.rcomp();}

		return extendCore(bb, leftCounts, rightCounts, distance, includeJunctionBase, kmer);
	}

	/**
	 * Shared core of both extension entry points above -- mirrors BFC1's
	 * extendToRight2raw structurally, one base at a time via Kmer.addRight
	 * instead of shift/mask, with fillLeftCounts/fillRightCounts/getCount
	 * translated to the Kmer/xor equivalents.
	 * Stops for insufficient support, junctions, or conflicting left predecessors
	 * when leftCounts is supplied. May append a junction base only on the forward
	 * canonical strand. Mutates the live anchor and caller buffers, not the table.
	 * @param bb Caller-owned builder to append to
	 * @param leftCounts Four-entry left-branch scratch array, or null
	 * @param rightCounts Scratch array with exactly four entries
	 * @param distance Maximum append length in bases
	 * @param includeJunctionBase Allow the terminal junction base on the canonical forward strand
	 * @param kmer Valid live anchor, advanced in place; must not be localKmer2
	 * @return Number of bases appended
	 */
	private int extendCore(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance,
			boolean includeJunctionBase, final Kmer kmer){
		final int initialLength=bb.length();

		int count=getCount(kmer);
		if(count<minCountCorrect){return 0;}

		int leftMaxPos=0;
		int leftMax=minCountExtend;
		int leftSecondPos=1;
		int leftSecond=0;

		if(leftCounts!=null){
			leftMaxPos=fillLeftCounts(kmer, leftCounts);
			leftMax=leftCounts[leftMaxPos];
			leftSecondPos=Tools.secondHighestPosition(leftCounts);
			leftSecond=leftCounts[leftSecondPos];
		}

		int rightMaxPos=fillRightCounts(kmer, rightCounts);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];

		if(rightMax<minCountExtend){return 0;}
		if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){return 0;}

		final int maxLen=bb.length()+distance;

		while(bb.length()<maxLen){
			final byte b=AminoAcid.numberToBase[rightMaxPos];

			//The base that falls off the left end as we roll right -- needed for
			//the "hidden branch" junction check below, mirroring BFC1's `evicted`.
			final byte evictedBase=kmer.addRight(b);
			final byte evictedNum=AminoAcid.baseToNumber[evictedBase];

			count=rightMax;

			if(leftCounts!=null){
				leftMaxPos=fillLeftCounts(kmer, leftCounts);
				leftMax=leftCounts[leftMaxPos];
				leftSecondPos=Tools.secondHighestPosition(leftCounts);
				leftSecond=leftCounts[leftSecondPos];
			}

			rightMaxPos=fillRightCounts(kmer, rightCounts);
			rightMax=rightCounts[rightMaxPos];
			rightSecondPos=Tools.secondHighestPosition(rightCounts);
			rightSecond=rightCounts[rightSecondPos];

			if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){
				//Mirrors BFC1's `kmer>rkmer` orientation check (append only when the
				//forward strand is canonical) -- see BFC1 for the reasoning; Tadpole2's
				//analogous check is array2-based (opposite sense), but BFC2 mirrors BFC1's
				//OWN semantic since that's the algorithm being ported. Flagged for the
				//STEP 4 Tadpole2 cross-check to catch if this guess is wrong.
				if(includeJunctionBase && kmer.key()==kmer.array1()){
					bb.append(b);
				}
				break;
			}

			if(leftCounts!=null && leftMaxPos!=evictedNum){
				if(includeJunctionBase && kmer.key()==kmer.array1()){
					bb.append(b);
				}
				break;
			}

			bb.append(b);

			if(rightMax<minCountExtend){break;}
		}
		return bb.length()-initialLength;
	}

	/**
	 * Fills array with counts for all possible left extensions of a live kmer.
	 * Tests all four candidate preceding bases via Kmer.addLeft on a scratch
	 * clone (addLeft prepends and evicts the current rightmost base -- exactly
	 * the candidate BFC1's shift-based version constructs).
	 * @param kmer Valid live anchor, preserved; must not alias localKmer2
	 * @param counts Four-entry output array in A/C/G/T order, overwritten
	 * @return Base code of the first maximum-count candidate, with ties in A/C/G/T order
	 */
	private int fillLeftCounts(final Kmer kmer, final int[] counts){
		final Kmer scratch=localKmer2.get();
		int max=-1, maxPos=0;
		for(int i=0; i<=3; i++){
			scratch.setFrom(kmer);
			final byte b=AminoAcid.numberToBase[i];
			scratch.addLeft(b);
			int count=Tools.max(getCount(scratch), 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}

	/**
	 * Fills array with counts for all possible right extensions of a live kmer.
	 * Tests all four candidate following bases via Kmer.addRight on a scratch
	 * clone (does not mutate the caller's live rolling kmer).
	 * @param kmer Valid live anchor, preserved; must not alias localKmer2
	 * @param counts Four-entry output array in A/C/G/T order, overwritten
	 * @return Base code of the first maximum-count candidate, with ties in A/C/G/T order
	 */
	private int fillRightCounts(final Kmer kmer, final int[] counts){
		final Kmer scratch=localKmer2.get();
		int max=-1, maxPos=0;
		for(int i=0; i<=3; i++){
			scratch.setFrom(kmer);
			final byte b=AminoAcid.numberToBase[i];
			scratch.addRight(b);
			int count=Tools.max(getCount(scratch), 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}

	/*--------------------------------------------------------------*/
	/*----------------   Abstract Seam: Extend edges ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Max count among all possible left extensions of the kmer at list position
	 * pos, rebuilt from bases (the stored hash cannot be un-hashed). Called only
	 * from smooth() for the sequence's leading edge.
	 * @param bases Sequence with an in-range k-base window at pos, read without modification
	 * @param kmers Unused stored hash list
	 * @param pos Window start index in bases
	 * @return Maximum nonnegative left-extension count, or zero for an undefined anchor
	 */
	@Override
	protected int maxLeftCount(byte[] bases, LongList kmers, int pos){
		final Kmer kmer=localKmer.get();
		kmer.clear();
		for(int i=pos; i<pos+k; i++){kmer.addRight(bases[i]);}
		//2026-08-21 (Nepgear+Amber): an N anywhere in bases[pos..pos+k-1] is a legitimate
		//condition, not a bug -- addRightNumeric resets len on an invalid base, so the
		//window can finish with kmer.len<k even though exactly k in-bounds bases were fed
		//in. Real Illumina data has N-calls, especially near read ends, which is exactly
		//where smooth() calls this (the read's own edges) -- confirmed live on real data
		//(BFC2, k=62), previously unexercised because all prior testing (Item 1 and Item
		//3a) used only randomreads.sh-generated data, which never inserts N's. Treat a
		//short window the same as fillKmers already does for one (no valid extension).
		if(kmer.len<k){return 0;}

		final Kmer scratch=localKmer2.get();
		int max=-1;
		for(int i=0; i<=3; i++){
			scratch.setFrom(kmer);
			scratch.addLeft(AminoAcid.numberToBase[i]);
			int count=Tools.max(getCount(scratch), 0);
			if(count>max){max=count;}
		}
		return max;
	}

	/**
	 * Max count among all possible right extensions of the kmer at list position
	 * pos, rebuilt from bases. Called only from smooth() for the trailing edge.
	 * @param bases Sequence with an in-range k-base window at pos, read without modification
	 * @param kmers Unused stored hash list
	 * @param pos Window start index in bases
	 * @return Maximum nonnegative right-extension count, or zero for an undefined anchor
	 */
	@Override
	protected int maxRightCount(byte[] bases, LongList kmers, int pos){
		final Kmer kmer=localKmer.get();
		kmer.clear();
		for(int i=pos; i<pos+k; i++){kmer.addRight(bases[i]);}
		//See maxLeftCount's comment above -- same N-at-the-edge condition, same fix.
		if(kmer.len<k){return 0;}

		final Kmer scratch=localKmer2.get();
		int max=-1;
		for(int i=0; i<=3; i++){
			scratch.setFrom(kmer);
			scratch.addRight(AminoAcid.numberToBase[i]);
			int count=Tools.max(getCount(scratch), 0);
			if(count>max){max=count;}
		}
		return max;
	}

	/*--------------------------------------------------------------*/
	/*----------------  Abstract Seam: Merge/Similar ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Tests count discontinuities around the two inferred overlap boundaries.
	 * Rebuilds stored hashes and reads them through getCountAt, preserving both
	 * hash lanes in dual mode. Returns true without testing when the read fails
	 * the initial length requirements. This is a merge heuristic, not proof that
	 * a read is nonchimeric. Each tested range must fit the dense stored-window list.
	 *
	 * Historical port notes (Amber) retained this override rather than hoisting
	 * the similar loop into the base during the port.
	 * @param merged Proposed merged read; bases are not modified
	 * @param len1 First original mate's length in bases
	 * @param len2 Second original mate's length in bases
	 * @param kmers Caller-owned primary hash list, replaced when the test runs
	 * @param width Number of neighboring windows on each side to examine, positive
	 * @param thresh Minimum supporting count for rejection against counts at most one
	 * @param highMult Count-ratio multiplier for rejection against positive counts
	 * @return false on a detected discontinuity; true when accepted or not tested
	 */
	@Override
	public boolean mergeOK(Read merged, int len1, int len2, LongList kmers, final int width, final int thresh, final long highMult){
		final int len=merged.length();
		final byte[] bases=merged.bases;
		if(len<len1+width+1 || len<len2+width+1 || len<k+width){return true;}
		fillKmers(bases, kmers);
		final int kmerCount=kmers.size();
		final int a=len-len2-1;
		final int b=a+1;
		final int c=len1-1;
		final int d=c+1;

		final int ak=a-k+1;
		final int bk=b;
		final int ck=c-k+1;
		final int dk=d;

		//BLM204: use dense window positions, not base length or the count of valid windows.
		if(ak-width>=0 && ak+width<kmerCount){
			int min=getCountAt(kmers, ak);
			for(int i=ak-width+1; i<ak; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, ak+1);
			for(int i=ak+2; i<=ak+width; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(ck-width>=0 && ck+width<kmerCount){
			int min=getCountAt(kmers, ck);
			for(int i=ck-width+1; i<ck; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, ck+1);
			for(int i=ck+2; i<=ck+width; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(bk-width>=0 && bk+width<kmerCount){
			int min=getCountAt(kmers, bk);
			for(int i=bk+1; i<bk+width+1; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, bk-1);
			for(int i=bk-width; i<bk-1; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(dk-width>=0 && dk+width<kmerCount){
			int min=getCountAt(kmers, dk);
			for(int i=dk+1; i<dk+width+1; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, dk-1);
			for(int i=dk-width; i<dk-1; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		return true;
	}

	/**
	 * Rebuilds the window starting at a, drops its first base and appends newBase.
	 * Compares that candidate's count with counts[a] using inherited similarity
	 * thresholds. Uses initialized localKmer scratch without modifying the inputs;
	 * the hash list cannot reconstruct the candidate and is not read.
	 * @param bases Sequence containing a valid k-base window at a
	 * @param a Start base index of the anchor window
	 * @param newBase Defined DNA base to append to the candidate
	 * @param kmers Unused stored primary hash list
	 * @param counts Read-only count list indexed by window start
	 * @return Whether the candidate count meets the configured similarity thresholds
	 */
	@Override
	protected boolean isSimilar(byte[] bases, int a, byte newBase, LongList kmers, IntList counts){
		final Kmer kmer=localKmer.get();
		kmer.clear();
		for(int i=a; i<a+k; i++){kmer.addRight(bases[i]);}
		assert(kmer.len==k) : "expected a fully valid kmer at a="+a;
		kmer.addRight(newBase);

		final int count=getCount(kmer);
		final int aCount=counts.get(a);
		return isSimilar(aCount, count);
	}

	/*--------------------------------------------------------------*/
	/*----------------  Abstract Seam: Primitives    ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns this thread's reusable secondary-hash list, allocating it on first use. */
	private LongList xor2List(){
		LongList list=localXor2List.get();
		if(list==null){
			list=new LongList();
			localXor2List.set(list);
		}
		return list;
	}

	/** Queries a valid live k-mer using the borrowed filter's one- or two-hash mode. */
	private int getCount(final Kmer kmer){
		return dualHash() ? filter.filter.read(kmer.xor(), kmer.xor2()) : filter.filter.read(kmer.xor());
	}

	/**
	 * Queries a stored window using its aligned secondary hash when dual mode is enabled.
	 * @param kmers Primary list paired with this thread's current secondary list
	 * @param pos In-range window index in both lists
	 * @return Table count, or zero for a negative undefined-window primary hash
	 */
	@Override
	protected int getCountAt(final LongList kmers, final int pos){
		final long key=kmers.get(pos);
		return key<0 ? 0 : (dualHash() ? filter.filter.read(key, xor2List().get(pos)) : filter.filter.read(key));
	}

	/** Whether the currently assigned filter requests the secondary hash lane. */
	private boolean dualHash(){return filter!=null && filter.dualHash;}

	/** Not a real kmer to decode (it's a hash proxy) -- formats the hash value itself, debug-only. */
	@Override
	protected StringBuilder toText(long kmer){
		return new StringBuilder("xor:").append(kmer);
	}

	/**
	 * Identity: Kmer.xor() is already strand-canonical (fwd and rcomp give the
	 * same value), so there is no separate
	 * "reverse complement" of a stored hash to compute.
	 * @param kmer Stored primary hash, not packed sequence
	 * @return The supplied hash unchanged
	 */
	@Override
	protected long rcomp(long kmer){
		return kmer;
	}

	/**
	 * Queries a hash pair; these arguments are not forward/reverse sequence encodings.
	 * @param kmer Primary canonical hash
	 * @param rkmer Secondary canonical hash, ignored in single-hash mode
	 * @return Table count using the filter's configured hash mode
	 */
	@Override
	public int getCount(long kmer, long rkmer){
		return dualHash() ? filter.filter.read(kmer, rkmer) : filter.filter.read(kmer);
	}

	/**
	 * Performs the legacy single-hash table lookup, even when dual mode is configured.
	 * Callers requiring dual-hash lookup must supply both hashes or a stored position.
	 * @param key Primary hash for a single-hash table; no sentinel handling
	 * @return Count from the single-hash lookup path
	 */
	@Override
	public int getCount(long key){
		return filter.filter.read(key);
	}

	/**
	 * Performs a single-hash lookup with negative-sentinel handling.
	 * This compatibility method does not recover the secondary hash; inherited
	 * stored-position processing uses the getCountAt override instead.
	 * @param kmer Primary hash for a single-hash table, or a negative undefined sentinel
	 * @return Zero for a negative value, otherwise the single-hash table count
	 */
	@Override
	public int getCount2(long kmer){
		return kmer<0 ? 0 : filter.filter.read(kmer);
	}

	/**
	 * Retains the primary canonical hash without combining it with the secondary one.
	 * @param kmer Primary canonical hash
	 * @param rkmer Unused secondary hash
	 * @return The supplied primary hash unchanged
	 */
	@Override
	public long toValue(long kmer, long rkmer){
		return kmer;
	}

	/** Secondary window hashes aligned with the most recently filled primary list on each thread. */
	private final ThreadLocal<LongList> localXor2List=new ThreadLocal<LongList>();

}
