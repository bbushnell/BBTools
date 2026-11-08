package bloom;

import java.util.Arrays;
import java.util.BitSet;

import dna.AminoAcid;
import kmer.AbstractKmerTable;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/**
 * Implements short-kmer correction primitives using two-bit packed longs.
 * Supports k-mer lengths up to 31; {@link BloomFilterCorrector} supplies the
 * correction orchestration. Stored values retain forward sequence orientation,
 * with -1 for windows containing an undefined base. Queries honor the attached
 * filter's strand policy: canonical selection or the forward packed value.
 * When ksmall is smaller than k, the filter combines smaller-window counts
 * through its getCountBig path.
 *
 * The count table is borrowed. Sequence builders, count lists and extension
 * buffers are caller-owned and may be modified as documented below. Concurrent
 * callers must supply separate mutable buffers and keep filter/settings stable.
 * Counts remain subject to the filter's collisions and saturation; they are
 * evidence for correction heuristics, not proof of the true sequence.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class BloomFilterCorrector1 extends BloomFilterCorrector{

	/**
	 * Stores the filter and requested k-mer geometry without allocating a table.
	 * @param filter_ Borrowed filter, required before any count lookup
	 * @param k_ Packed k-mer length in bases, from 1 through 31
	 * @param ksmall_ Filter's small-window length, positive and at most k_
	 */
	public BloomFilterCorrector1(BloomFilter filter_, int k_, int ksmall_){
		super(filter_, k_, ksmall_);
	}

	/**
	 * Recounts windows spanning the interval from the first through last changed base.
	 * Rewrites the whole affected interval, including gaps between changed bits;
	 * undefined windows receive zero. When smoothing is enabled, also suppresses
	 * strict interior count peaks across the entire list. Does not modify bases.
	 *
	 * @param bases Non-null modified sequence
	 * @param counts Existing count list indexed by window start, updated in place
	 * @param dummy Unused compatibility parameter; may be null
	 * @param changed Nonempty set of in-range base indices, read without modification
	 * @return Number of valid windows recounted before smoothing
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
	 * Recounts the anchor window and up to k following windows after a base change.
	 * Starts at window ca, clips the scan at the read end and writes zero for
	 * undefined windows. Unlike the changed-bit overload, does not smooth counts.
	 *
	 * @param bases Non-null sequence, read without modification
	 * @param counts Existing full-length count list, updated without changing its size
	 * @param ca In-range anchor-window start preceding the changed base at ca+k
	 * @return Number of valid windows recounted
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
		 * i is a base index; i-k+1 is the corresponding count-list index. */
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
	 * Scans left to right and attempts substitutions at count drops.
	 * Changes bb's bases in place and refreshes affected counts after each recorded
	 * correction; it neither extends the builder nor changes quality scores.
	 * Reads too short for a candidate beyond the anchor and dead zone return zero.
	 * Candidate identity is compared as a base code, independently of its support count (BLM203).
	 *
	 * @param bb Mutable read builder with backing capacity equal to logical length; recounting uses array length
	 * @param quals Read-only qualities aligned to bases, or null to use quality 20
	 * @param rightCounts Caller-owned scratch array with exactly four entries
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
							//BLM203: compare base codes, not a base code with its support count.
							//An unchanged winning base is not a correction.
							if(num==rightMaxPos){
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
	 * Builds an anchor from the trailing k bases and attempts rightward extension.
	 * Returns zero without extension if the builder is too short or the anchor
	 * contains an undefined base; otherwise delegates to the packed extension core.
	 *
	 * @param bb Caller-owned sequence builder to append to
	 * @param leftCounts Four-entry scratch buffer for left-branch checks, or null to omit them
	 * @param rightCounts Scratch buffer with exactly four entries, overwritten
	 * @param distance Maximum number of bases to append
	 * @param includeJunctionBase Allow the terminal junction base when its forward strand is canonical
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
			}else{assert(len==k);}
		}
		return extendToRight2raw(bb, leftCounts, rightCounts, distance, includeJunctionBase, kmer, rkmer);
	}

	/**
	 * Uses a stored forward window, optionally reverse-complemented, as the anchor.
	 * No sequence reconstruction is needed for this packed representation.
	 * @param bb Caller-owned builder to append to
	 * @param leftCounts Four-entry left-branch scratch buffer, or null
	 * @param rightCounts Scratch buffer with exactly four entries
	 * @param distance Maximum number of bases to append
	 * @param includeJunctionBase Allow the terminal junction base on the canonical forward strand
	 * @param bases Unused in this representation; supplied by the superclass contract
	 * @param kmers Read-only list containing a valid packed anchor at originPos
	 * @param originPos In-range anchor index in kmers
	 * @param fromRcomp Whether to reverse-complement the anchor before extending
	 * @return Number of bases appended
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
	 * Appends at most distance bases from an already-encoded anchor.
	 * Stops at insufficient support, junctions, or a conflicting left predecessor
	 * when leftCounts is supplied. Junction-base inclusion is conditional on the
	 * forward packed strand being larger than its reverse complement. Does not
	 * claim graph ownership or modify the count table.
	 * @param bb Caller-owned builder to append to
	 * @param leftCounts Four-entry left-branch scratch buffer, or null
	 * @param rightCounts Scratch buffer with exactly four entries
	 * @param distance Maximum append length in bases
	 * @param includeJunctionBase Allow the terminal junction base on the canonical forward strand
	 * @param kmer Valid packed forward anchor
	 * @param rkmer Reverse complement of kmer
	 * @return Number of bases appended
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
	 * Refreshes windows affected by a substitution at base a+k using stored anchor a.
	 * Leaves anchor a unchanged and rolls through the next k window starts, clipped
	 * to counts.size. Writes -1/zero for undefined windows. The assertion on each
	 * valid replacement assumes that the corresponding old packed value changed.
	 *
	 * @param bases Non-null modified sequence, read without modification
	 * @param kmers Full aligned window list with a valid unchanged anchor at a; updated in place
	 * @param counts Full aligned count list, updated in place
	 * @param a Start index of the unchanged window immediately preceding the changed base
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
	 * Counts the four windows formed by prepending a base and dropping the last base.
	 * Uses packed values directly; does not allocate a sequence or modify the table.
	 * @param kmer Valid forward packed k-mer, not the undefined-window sentinel
	 * @return Maximum nonnegative count among the four candidate windows
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

	/**
	 * Counts the four windows formed by appending a base and dropping the first base.
	 * Uses packed values directly; does not allocate a sequence or modify the table.
	 * @param kmer Valid forward packed k-mer, not the undefined-window sentinel
	 * @return Maximum nonnegative count among the four candidate windows
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

	/**
	 * Reads the packed anchor at pos and finds its strongest left extension.
	 * @param bases Unused for this representation
	 * @param kmers Read-only packed-window list
	 * @param pos In-range index of a valid packed window
	 * @return Maximum nonnegative count among the four left candidate windows
	 */
	@Override
	protected int maxLeftCount(byte[] bases, LongList kmers, int pos){
		return maxLeftCount(kmers.get(pos));
	}

	/**
	 * Reads the packed anchor at pos and finds its strongest right extension.
	 * @param bases Unused for this representation
	 * @param kmers Read-only packed-window list
	 * @param pos In-range index of a valid packed window
	 * @return Maximum nonnegative count among the four right candidate windows
	 */
	@Override
	protected int maxRightCount(byte[] bases, LongList kmers, int pos){
		return maxRightCount(kmers.get(pos));
	}

	/**
	 * Replaces counts[0..3] with nonnegative counts for the A/C/G/T left extensions.
	 * Each candidate prepends one base and evicts the current last base. Ties for
	 * the largest count select the first base in A/C/G/T order.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param counts Caller-owned output array with at least four entries
	 * @param mask Mask retaining the low 2*k bits
	 * @param shift2 Position of the first base's low bit, 2*(k-1)
	 * @return Base code of the first maximum-count extension, from 0 through 3
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

	/**
	 * Replaces counts[0..3] with nonnegative counts for the A/C/G/T right extensions.
	 * Each candidate appends one base and evicts the current first base. Ties for
	 * the largest count select the first base in A/C/G/T order.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param counts Caller-owned output array with at least four entries
	 * @param mask Mask retaining the low 2*k bits
	 * @param shift2 Position of the first base's low bit, 2*(k-1)
	 * @return Base code of the first maximum-count extension, from 0 through 3
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

	/**
	 * Stores each forward packed window, using -1 for windows with undefined bases.
	 * For a read of at least k bases, clears kmers and writes length-k+1 entries.
	 * A shorter read returns zero before clearing the list; callers must not treat
	 * its old contents as windows from that read. Does not modify bases.
	 * @param bases Non-null sequence of nonnegative encoded character bytes
	 * @param kmers Caller-owned list to replace when the read is long enough
	 * @return Number of valid windows, excluding the -1 entries
	 */
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
	 * Advances the stored window at a by one base and compares its candidate count.
	 * Drops the first base, appends newBase, then applies the inherited count-
	 * similarity heuristic against counts[a]. Does not mutate any supplied buffer.
	 *
	 * @param bases Unused for this representation
	 * @param a In-range index of a valid stored window
	 * @param newBase Defined DNA base to append to the candidate window
	 * @param kmers Read-only forward packed-window list
	 * @param counts Read-only count list aligned with kmers
	 * @return Whether the candidate count satisfies the configured similarity thresholds
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

	/** Formats a packed k-mer as k DNA bases in a new builder, for diagnostics. */
	@Override
	protected StringBuilder toText(long kmer){return AbstractKmerTable.toText(kmer, k);}
	/** Returns the reverse complement of the low 2*k bits, without canonicalizing. */
	@Override
	protected long rcomp(long kmer){return AminoAcid.reverseComplementBinaryFast(kmer, k);}
	/**
	 * Queries a strand pair using the filter's strand policy when k equals ksmall.
	 * Otherwise queries its constituent smaller windows.
	 * The smaller-window path derives all queries from kmer and ignores rkmer.
	 * @param kmer Valid forward packed k-mer
	 * @param rkmer Its packed reverse complement
	 * @return Direct count, or minimum count across the constituent small windows
	 */
	@Override
	public int getCount(long kmer, long rkmer){
		return (k==ksmall ? filter.getCount(kmer, rkmer) : filter.getCountBig(kmer));
	}
	/**
	 * Queries an already-selected key, or the small windows encoded within it.
	 * Does not recognize undefined-window sentinels; use getCount2 for that case.
	 * @param key Selected table key when k equals ksmall, otherwise a valid packed k-mer
	 * @return Direct count, or minimum count across the constituent small windows
	 */
	@Override
	public int getCount(long key){
		return (k==ksmall ? filter.getCount(key) : filter.getCountBig(key));
	}
	/**
	 * Returns zero for negative values; otherwise applies the filter's strand policy.
	 * With ksmall smaller than k, the filter applies its strand policy to each small window.
	 * @param kmer Stored forward packed k-mer, or a negative undefined-window sentinel
	 * @return Direct or minimum constituent-window count; zero for a negative input
	 */
	@Override
	public int getCount2(long kmer){
		return kmer<0 ? 0 : (k==ksmall ? filter.getCount(toValue(kmer, rcomp(kmer))) : filter.getCountBig(kmer));
	}
	/**
	 * Selects a packed strand using the attached filter's policy.
	 * Without a filter, retains the original canonical key-conversion behavior.
	 *
	 * @param kmer Valid forward packed k-mer
	 * @param rkmer Its packed reverse complement
	 * @return Packed strand value to query; no hashing or table access occurs here
	 */
	@Override
	public long toValue(long kmer, long rkmer){
		//BLM206: stored and extension queries must use the same strand policy as the table.
		long value=(filter==null || filter.rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}

}
