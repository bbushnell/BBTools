package bloom;

import java.util.Arrays;
import java.util.BitSet;

import assemble.ErrorTracker;
import dna.AminoAcid;
import kmer.AbstractKmerTable;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/**
 * Long-primitive (k&lt;=31) implementation of Bloom filter-based sequencing
 * error correction. Encodes kmers as bit-packed longs via shift/mask rolling.
 * All representation-agnostic orchestration lives in the abstract superclass
 * {@link BloomFilterCorrector}; this class implements only the kmer-representation
 * primitives (the abstract seam).
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class BloomFilterCorrector1 extends BloomFilterCorrector {

	/**
	 * Constructs a BloomFilterCorrector1 with specified parameters.
	 * @param filter_ The Bloom filter containing reference k-mer counts
	 * @param k_ K-mer length for error correction operations
	 * @param ksmall_ Small k-mer length for efficient lookups (must be <= k_)
	 */
	public BloomFilterCorrector1(BloomFilter filter_, int k_, int ksmall_) {
		super(filter_, k_, ksmall_);
	}

	/**
	 * Regenerates k-mer counts for positions that were changed during correction.
	 * Updates count array after bases have been modified to reflect new k-mer frequencies.
	 *
	 * @param bases Modified sequence bases
	 * @param counts Count array to update
	 * @param dummy Unused Kmer parameter
	 * @param changed BitSet indicating which positions were modified
	 * @return Number of valid k-mers processed
	 */
	@Override
	public int regenerateCounts(byte[] bases, IntList counts, final Kmer dummy, BitSet changed){
		assert(!changed.isEmpty());
		final int firstBase=changed.nextSetBit(0), lastBase=changed.length()-1;
		final int ca=firstBase-k;
		final int firstCount=Tools.max(firstBase-k+1, 0), lastCount=Tools.min(counts.size-1, lastBase); //count limit

		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		int valid=0;

		for(int i=Tools.max(0, firstBase-k+1), lim=Tools.min(lastBase+k-1, bases.length-1); i<=lim; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			final long x2=AminoAcid.baseToComplementNumber[base];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
			}

			final int c=i-k+1;
			if(i>=firstBase){
				if(len>=k){
					valid++;
					int count=getCount(kmer, rkmer);
					counts.set(c, count);
				}else if(c>=0){
					counts.set(c, 0);
				}
			}
		}

		if(smooth){
			int[] array=counts.array;
			for(int i=1; i<counts.size-1; i++){
				int a=array[i-1], b=array[i], c=array[i+1];
				if(b>a && b>c){array[i]=Tools.max(a, c);}
			}
		}

		return valid;
	}

	/**
	 * Regenerates k-mer counts for a specific region after base changes.
	 * Updates count array for positions affected by sequence modifications.
	 *
	 * @param bases Modified sequence bases
	 * @param counts Count array to update
	 * @param ca Starting count array index
	 * @return Number of valid k-mers processed
	 */
	@Override
	public int regenerateCounts(byte[] bases, IntList counts, final int ca){
		final int b=ca+k-1;
		final int lim=Tools.min(bases.length, b+k+1);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		int valid=0;

		final int clen=counts.size;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts.
		 * i is an index in the base array, j is an index in the count array. */
		for(int i=b-k+1; i<lim; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			final long x2=AminoAcid.baseToComplementNumber[base];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
			}

			if(i>=b){
				if(len>=k){
					valid++;
					int count=getCount(kmer, rkmer);
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

	/**
	 * Inner reassembly algorithm for single-direction extension.
	 * Extends sequence rightward looking for error patterns and corrections.
	 * Identifies substitution errors by analyzing k-mer count transitions.
	 *
	 * @param bb ByteBuilder containing sequence to extend
	 * @param quals Quality scores
	 * @param rightCounts Buffer for extension counts
	 * @param counts K-mer count values
	 * @param errorExtension Maximum extension distance
	 * @return Number of errors corrected
	 */
	@Override
	public int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension){
		final int length=bb.length();
		if(length<k+1+deadZone){return 0;}
		final byte[] bases=bb.array;

		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;

		int detected=0;
		int corrected=0;
		int len=0;

		//a is the index of the rightmost base of the kmer
		//b=a+1 is the index of the next base
		//ca=a-k+1 is the index of the count of the kmer
		//cb=a-k+2 is the index of the count of the next kmer
		for(int a=0, lim=length-deadZone-1; a<lim; a++){

			//Generate the new kmer
			final byte aBase=bases[a];
			final long x=AminoAcid.baseToNumber[aBase];
			final long x2=AminoAcid.baseToComplementNumber[aBase];

			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{


				//Now consider the next kmer
				kmer=((kmer<<2)|(long)x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				len++;

				if(verbose){
					System.err.println("len: "+len+" vs "+k+"; a="+a);
				}

				if(len>=k){

					final int b=a+1;
					final int ca=a-k+1;
					final int cb=ca+1;

					final int aCount=counts.get(ca);
					final int bCount=counts.get(cb);
					final byte qb=(quals==null ? 20 : quals[b]);

					if(verbose){
						System.err.println("ca="+ca+", cb="+cb+"; aCount="+aCount+", bCount="+bCount);
						System.err.println(isError(aCount, bCount, qb)+", "+isSimilar(aCount, ca-errorExtension, ca-1, counts)+
								", "+isError(aCount, ca+2, ca+k, counts));
					}

					if(isSubstitution(ca, errorExtension, qb, counts)){
						if(verbose){
							System.err.println("***Found error: "+aCount+", "+bCount);
						}
						//Assume a 1bp substitution; attempt to correct.

						int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
						int rightMax=rightCounts[rightMaxPos];
						int rightSecondPos=Tools.secondHighestPosition(rightCounts);
						int rightSecond=rightCounts[rightSecondPos];

						byte base=bases[b];
						byte num=AminoAcid.baseToNumber[base];

						if(rightMax>=minCountExtend){
							detected++;
							if(num==rightMax){
								detected--;
							}else if((isError(rightMax, rightSecond, qb) || !isJunction(rightMax, rightSecond)) && isSimilar(aCount, rightMax)){
								bases[b]=AminoAcid.numberToBase[rightMaxPos];
								corrected++;
								regenerateCounts(bases, counts, ca);
								if(verbose){System.err.println("Corrected error: "+num+"->"+rightMaxPos+". New counts:\n"+counts);}
							}
						}

					}else{
						if(verbose){
							System.err.println("Not an error: "+aCount+", "+bCount+
									";  "+isError(aCount, bCount, qb)+", "+isSimilar(aCount, a-errorExtension, a-1, counts)+", "+isError(aCount, a+2, a+k, counts));
						}
					}
				}
			}
		}

		return corrected;
	}

	/**
	 * Extends sequence rightward by generating k-mer from existing sequence.
	 * Wrapper method that computes the rightmost k-mer and calls main extension.
	 *
	 * @param bb ByteBuilder containing initial sequence
	 * @param leftCounts Buffer for left neighbor analysis (may be null)
	 * @param rightCounts Buffer for right neighbor counts
	 * @param distance Maximum extension distance
	 * @param includeJunctionBase Whether to include base at junction positions
	 * @return Number of bases added to sequence
	 */
	@Override
	public int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase){
		if(verbose || verbose2){outstream.println("Entering extendToRight2 (no kmers).");}
		final int initialLength=bb.length();
		if(initialLength<k){return 0;}
		final int k2=k-1;
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		long rkmer=0;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			int len=0;
			final byte[] bases=bb.array;
			for(int i=initialLength-k; i<initialLength; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("B: Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
			if(len<k){
				if(verbose || verbose2){outstream.println("Returning because len<k: "+len+"<"+k);}
				return 0;
			}
			else{assert(len==k);}
		}
		return extendToRight2raw(bb, leftCounts, rightCounts, distance, includeJunctionBase, kmer, rkmer);
	}

	/**
	 * Extends rightward starting from the kmer stored at list position originPos
	 * (or its reverse complement, if fromRcomp). BFC1 already has the actual
	 * kmer content in the LongList, so this is a cheap O(1) lookup + optional
	 * rcomp, then delegates to the same raw extension core as the 5-arg overload.
	 */
	@Override
	protected int extendFromStored(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance,
			boolean includeJunctionBase, final byte[] bases, final LongList kmers, final int originPos, final boolean fromRcomp){
		final long stored=kmers.get(originPos);
		final long kmer, rkmer;
		if(!fromRcomp){
			kmer=stored;
			rkmer=rcomp(stored);
		}else{
			rkmer=stored;
			kmer=rcomp(stored);
		}
		return extendToRight2raw(bb, leftCounts, rightCounts, distance, includeJunctionBase, kmer, rkmer);
	}

	/**
	 * Extend these bases to the right by at most 'distance'.
	 * Stops at right junctions only.
	 * Does not claim ownership.
	 * Not part of the abstract contract -- a k&gt;31 kmer cannot be represented as a
	 * single long, so this exact (kmer,rkmer) signature is BFC1-only; it is used
	 * internally by the two abstract entry points below (extendToRight2's 5-arg
	 * form, and extendFromStored).
	 */
	private int extendToRight2raw(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase,
			long kmer, long rkmer){
		if(verbose || verbose2){outstream.println("Entering extendToRight2 (with kmers).");}
		final int initialLength=bb.length();
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));

		/* Now the trailing kmer has been initialized. */

		long key=toValue(kmer, rkmer);
		int count=getCount(key);
		if(count<minCountCorrect){
			if(verbose || verbose2){outstream.println("Returning because count was too low: "+count+"<"+minCountCorrect);}
			return 0;
		}

		int leftMaxPos=0;
		int leftMax=minCountExtend;
		int leftSecondPos=1;
		int leftSecond=0;

		if(leftCounts!=null){
			leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
			leftMax=leftCounts[leftMaxPos];
			leftSecondPos=Tools.secondHighestPosition(leftCounts);
			leftSecond=leftCounts[leftSecondPos];
		}

		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];

		if(verbose){
			outstream.println("kmer: "+toText(kmer)+", "+toText(rkmer));
			outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
			outstream.println("rightMaxPos="+rightMaxPos);
			outstream.println("rightMax="+rightMax);
			outstream.println("rightSecondPos="+rightSecondPos);
			outstream.println("rightSecond="+rightSecond);
		}

		if(rightMax<minCountExtend){
			if(verbose || verbose2){outstream.println("Returning because rightMax was too low: "+rightMax+"<"+minCountExtend+"\n"+count+", "+Arrays.toString(rightCounts));}
			return 0;
		}
		if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){
			if(verbose || verbose2){outstream.println("Returning because isJunction: "+rightMax+", "+rightSecond+"; "+leftMax+", "+leftSecond);}
			return 0;
		}

		final int maxLen=bb.length()+distance;

		while(bb.length()<maxLen){

			//Generate the new kmer
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[(int)x];

			final long evicted=(kmer>>>shift2); //Binary value that falls off the end.

			//Now consider the next kmer
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			key=toValue(kmer, rkmer);

			assert(getCount(key)==rightMax || rightMax==0);
			count=rightMax;

			assert(count>=minCountExtend) : count;

			if(leftCounts!=null){
				leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
				leftMax=leftCounts[leftMaxPos];
				leftSecondPos=Tools.secondHighestPosition(leftCounts);
				leftSecond=leftCounts[leftSecondPos];
			}

			rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
			rightMax=rightCounts[rightMaxPos];
			rightSecondPos=Tools.secondHighestPosition(rightCounts);
			rightSecond=rightCounts[rightSecondPos];

			if(verbose){
				outstream.println("kmer: "+toText(kmer)+", "+toText(rkmer));
				outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}

			if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){
				if(verbose){outstream.println("B: Breaking because isJunction("+rightMax+", "+rightSecond+", "+leftMax+", "+leftSecond+")");}
				if(includeJunctionBase && kmer>rkmer){
					bb.append(b);
					if(verbose){outstream.println("Added base "+(char)b);}
				}
				break;
			}

			if(leftCounts!=null && leftMaxPos!=evicted){
				if(verbose){outstream.println("B: Breaking because of hidden branch: leftMaxPos!=evicted ("+leftMaxPos+"!="+evicted+")" +
						"\nleftMaxPos="+leftMaxPos+", leftMax="+leftMax+", leftSecondPos="+leftSecondPos+", leftSecond="+leftSecond);}
				if(includeJunctionBase && kmer>rkmer){
					bb.append(b);
					if(verbose){outstream.println("Added base "+(char)b);}
				}
				break;
			}

			bb.append(b);
			if(verbose){outstream.println("Added base "+(char)b);}

			if(rightMax<minCountExtend){
				if(verbose || verbose2){outstream.println("C: Breaking because highest right was too low: "+rightMax+"<"+minCountExtend);}
				break;
			}
		}
		if(verbose || verbose2){System.err.println("Extended by "+(bb.length()-initialLength));}
		return bb.length()-initialLength;
	}

	/**
	 * Regenerates k-mers and counts after sequence modification.
	 * Updates k-mer list and count array for positions affected by base changes.
	 *
	 * @param bases Modified sequence bases
	 * @param kmers K-mer list to update
	 * @param counts Count array to update
	 * @param a Starting index for regeneration
	 */
	@Override
	public void regenerateKmers(byte[] bases, LongList kmers, IntList counts, final int a){
		final int loc=a+k;
		final int lim=Tools.min(counts.size, a+k+1);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=kmers.get(a);
		long rkmer=rcomp(kmer);
		int len=k;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=loc, j=a+1; j<lim; i++, j++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}

			if(len>=k){
				assert(kmers.get(j)!=kmer);
				kmers.set(j, kmer);
				int count=getCount(kmer, rkmer);
				counts.set(j, count);
			}else{
				kmers.set(j, -1);
				counts.set(j, 0);
			}
		}
	}

	/**
	 * Cheap O(1) form using the stored kmer directly. Called by the abstract
	 * bases-based override below; not part of the abstract contract itself
	 * (BFC2 cannot accept a raw long).
	 * @param kmer K-mer to extend leftward
	 * @return Maximum count among left extensions
	 */
	private int maxLeftCount(long kmer){
		long rkmer=rcomp(kmer);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		rkmer=(rkmer<<2)&mask;
		kmer=(kmer>>>2);
		int max=-1, maxPos=0;
		for(int i=0; i<=3; i++){
			long rkmer2=rkmer|((long)AminoAcid.numberToComplement[i]);
			long kmer2=kmer|(((long)i)<<shift2);
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			long key=toValue(rkmer2, kmer2);
			int count=getCount(key);
			count=Tools.max(count, 0);
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return max;
	}

	/**
	 * Cheap O(1) form using the stored kmer directly. Called by the abstract
	 * bases-based override below; not part of the abstract contract itself.
	 * @param kmer K-mer to extend rightward
	 * @return Maximum count among right extensions
	 */
	private int maxRightCount(long kmer){
		long rkmer=rcomp(kmer);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		if(verbose){outstream.println("fillRightCounts:   "+toText(kmer)+",   "+toText(rkmer));}
		kmer=(kmer<<2)&mask;
		rkmer=(rkmer>>>2);
		int max=-1, maxPos=0;

		for(int i=0; i<=3; i++){
			long kmer2=kmer|((long)i);
			long rkmer2=rkmer|(((long)AminoAcid.numberToComplement[i])<<shift2);
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			assert(kmer2==rcomp(rkmer2));
			long key=toValue(kmer2, rkmer2);
			int count=getCount(key);
			count=Tools.max(count, 0);
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return max;
	}

	/** Abstract-contract entry: O(1) lookup from the stored kmer, delegates to the raw helper above. */
	@Override
	protected int maxLeftCount(byte[] bases, LongList kmers, int pos){
		return maxLeftCount(kmers.get(pos));
	}

	/** Abstract-contract entry: O(1) lookup from the stored kmer, delegates to the raw helper above. */
	@Override
	protected int maxRightCount(byte[] bases, LongList kmers, int pos){
		return maxRightCount(kmers.get(pos));
	}

	/**
	 * Fills array with counts for all possible left extensions of a k-mer.
	 * Computes counts for all four bases that could precede the given k-mer.
	 * Not part of the abstract contract -- BFC1-internal, used only by
	 * extendToRight2raw and reassemble_inner's own rolling state.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param counts Output array for extension counts
	 * @param mask Bit mask for k-mer operations
	 * @param shift2 Bit shift value for operations
	 * @return Index of extension with maximum count
	 */
	private int fillLeftCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){
		assert(kmer==rcomp(rkmer));
		rkmer=(rkmer<<2)&mask;
		kmer=(kmer>>>2);
		int max=-1, maxPos=0;
		for(int i=0; i<=3; i++){
			long rkmer2=rkmer|((long)AminoAcid.numberToComplement[i]);
			long kmer2=kmer|(((long)i)<<shift2);
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			long key=toValue(rkmer2, kmer2);
			int count=getCount(key);
			count=Tools.max(count, 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}

	/**
	 * Fills array with counts for all possible right extensions of a k-mer.
	 * Computes counts for all four bases that could follow the given k-mer.
	 * Not part of the abstract contract -- BFC1-internal, same reasoning as
	 * fillLeftCounts above.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param counts Output array for extension counts
	 * @param mask Bit mask for k-mer operations
	 * @param shift2 Bit shift value for operations
	 * @return Index of extension with maximum count
	 */
	private int fillRightCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){
		assert(kmer==rcomp(rkmer));
		if(verbose){outstream.println("fillRightCounts:   "+toText(kmer)+",   "+toText(rkmer));}
		kmer=(kmer<<2)&mask;
		rkmer=(rkmer>>>2);
		int max=-1, maxPos=0;

		for(int i=0; i<=3; i++){
			long kmer2=kmer|((long)i);
			long rkmer2=rkmer|(((long)AminoAcid.numberToComplement[i])<<shift2);
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			assert(kmer2==rcomp(rkmer2));
			long key=toValue(kmer2, rkmer2);
			int count=getCount(key);
			count=Tools.max(count, 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}

	/** Returns number of valid kmers */
	@Override
	public int fillKmers(byte[] bases, LongList kmers){
		final int blen=bases.length;
		if(blen<k){return 0;}
		final int min=k-1;
		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		int len=0;
		int valid=0;

		kmers.clear();

		/* Loop through the bases, maintaining a forward kmer via bitshifts */
		for(int i=0; i<blen; i++){
			final byte b=bases[i];
			assert(b>=0) : Arrays.toString(bases);
			final long x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){
				len=0;
				kmer=0;
			}else{len++;}
			if(i>=min){
				if(len>=k){
					kmers.add(kmer);
					valid++;
				}else{
					kmers.add(-1);
				}
			}
		}
		return valid;
	}

	/**
	 * Tests if a proposed base change results in similar k-mer count.
	 * Verifies that the new k-mer has a count similar to the original
	 * to avoid introducing artifacts during correction.
	 *
	 * @param a K-mer index
	 * @param newBase Proposed replacement base
	 * @param kmers K-mer list
	 * @param counts K-mer count values
	 * @return true if the new k-mer has similar count
	 */
	@Override
	protected boolean isSimilar(byte[] bases, int a, byte newBase, LongList kmers, IntList counts){
		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=kmers.get(a);

		final long x=AminoAcid.baseToNumber[newBase];
		kmer=((kmer<<2)|x)&mask;
		long rkmer=rcomp(kmer);
		int count=getCount(kmer, rkmer);
		int aCount=counts.get(a);
		boolean similar=isSimilar(aCount, count);
		return similar;
	}

	/** Converts encoded k-mer to text representation */
	@Override
	protected StringBuilder toText(long kmer){return AbstractKmerTable.toText(kmer, k);}
	/** Computes reverse complement of encoded k-mer */
	@Override
	protected long rcomp(long kmer){return AminoAcid.reverseComplementBinaryFast(kmer, k);}
	/**
	 * Gets count for k-mer using appropriate method based on k size.
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return K-mer count from Bloom filter
	 */
	@Override
	public int getCount(long kmer, long rkmer){
		return (k==ksmall ? filter.getCount(kmer, rkmer) : filter.getCountBig(kmer));
	}
	/**
	 * Gets count for k-mer key using appropriate method based on k size.
	 * @param key K-mer key value
	 * @return K-mer count from Bloom filter
	 */
	@Override
	public int getCount(long key){
		return (k==ksmall ? filter.getCount(key) : filter.getCountBig(key));
	}
	/**
	 * Gets count for k-mer with invalid k-mer handling.
	 * Returns 0 for invalid k-mers (negative values).
	 * @param kmer K-mer to look up
	 * @return K-mer count, or 0 if k-mer is invalid
	 */
	@Override
	public int getCount2(long kmer){
		return kmer<0 ? 0 : (k==ksmall ? filter.getCount(toValue(kmer, rcomp(kmer))) : filter.getCountBig(kmer));
	}
	/**
	 * Converts k-mer pair to lookup key value.
	 * Uses canonical representation if reverse complement mode is enabled.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return Key value for Bloom filter lookup
	 */
	@Override
	public long toValue(long kmer, long rkmer){
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}

}
