package bloom;

import java.io.PrintStream;
import java.util.BitSet;

import assemble.ErrorTracker;
import assemble.Rollback;
import shared.Tools;
import simd.Vector;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/**
 * Bloom filter-based sequencing error correction system using k-mer frequency analysis.
 * Identifies and corrects sequencing errors by analyzing k-mer abundance patterns
 * with memory-efficient Bloom filters. Supports multiple correction strategies including
 * pincer-based correction, tail correction, and reassembly-based approaches.
 *
 * Holds all kmer-representation-agnostic logic. The kmer-representation primitives
 * (encoding, rolling, and count lookup) are abstract; {@link BloomFilterCorrector1}
 * implements them with bit-packed longs (k&lt;=31) and {@link BloomFilterCorrector2}
 * with {@link ukmer.Kmer} objects (k&gt;31).
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public abstract class BloomFilterCorrector {

	/**
	 * Constructs a BloomFilterCorrector with specified parameters.
	 * @param filter_ The Bloom filter containing reference k-mer counts
	 * @param k_ K-mer length for error correction operations
	 * @param ksmall_ Small k-mer length for efficient lookups (must be <= k_)
	 */
	public BloomFilterCorrector(BloomFilter filter_, int k_, int ksmall_) {
		filter=filter_;
		k=k_;
		ksmall=ksmall_;
		assert(ksmall<=k);
	}

	/**
	 * Performs error correction on a read using thread-local storage.
	 * Initializes thread-local data structures and delegates to the main correction method.
	 * @param r The read to error correct
	 * @return Number of errors corrected
	 */
	public int errorCorrect(Read r){
		initializeThreadLocals();
		int corrected=errorCorrect(r, localLeftCounts.get(), localRightCounts.get(), localLongList.get(),
				localIntList.get(), localIntList2.get(), localByteBuilder.get(), localByteBuilder2.get(), localTracker.get(), localBitSet.get());
		return corrected;
	}

	/**
	 * Wrapper method that delegates to the main error correction implementation.
	 *
	 * @param r The read to correct
	 * @param leftCounts Buffer for left extension counts
	 * @param rightCounts Buffer for right extension counts
	 * @param kmers List of k-mers from the read
	 * @param counts K-mer count values
	 * @param counts2 Secondary count buffer
	 * @param bb Primary byte buffer
	 * @param bb2 Secondary byte buffer
	 * @param tracker Error tracking statistics
	 * @param bs BitSet for marking changed positions
	 * @param kmer Kmer object for operations
	 * @param kmer2 Secondary Kmer object
	 * @return Number of errors corrected
	 */
	public int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, LongList kmers, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs, Kmer kmer, Kmer kmer2){
		return errorCorrect(r, leftCounts, rightCounts, kmers, counts, counts2, bb, bb2, tracker, bs);
	}

	/**
	 * Fast error detection by sampling k-mers at regular intervals.
	 * Checks k-mer counts and transitions to identify potential errors without
	 * performing full correction analysis.
	 *
	 * @param kmers List of k-mers from the sequence
	 * @return true if errors are likely present, false otherwise
	 */
	public boolean hasErrorsFast(LongList kmers){
		if(kmers.size<1){return false;}
		int prev=-1;

		final int incr=Tools.mid(1, k/2, 9), mcc=minCountCorrect();
		for(int i=0; i<kmers.size; i+=incr){
			long kmer=kmers.get(i);
			if(kmer<0){
				return true;
			}
			long rkmer=rcomp(kmer);
			int count=getCount(kmer, rkmer);
			final int min=Tools.min(count, prev), max=Tools.max(count, prev);
			if(count<mcc || (i>0 && (isError(max+1, min-1)))){return true;}
			prev=count;
		}

		long kmer=kmers.get(kmers.size()-1);
		if(kmer<0){return true;}
		long rkmer=rcomp(kmer);
		int count=getCount(kmer, rkmer);
		final int min=Tools.min(count, prev), max=Tools.max(count, prev);
		return count<mcc || isError(max+1, min-1);
	}

	/**
	 * Main error correction method implementing multiple correction strategies.
	 * Applies pincer correction, tail correction, and reassembly approaches
	 * with rollback capability for quality control.
	 *
	 * @param r The read to correct
	 * @param leftCounts Buffer for left extension counts
	 * @param rightCounts Buffer for right extension counts
	 * @param kmers List of k-mers from the read
	 * @param counts K-mer count values
	 * @param counts2 Secondary count buffer
	 * @param bb Primary byte buffer
	 * @param bb2 Secondary byte buffer
	 * @param tracker Error tracking statistics
	 * @param bs BitSet for marking changed positions
	 * @return Number of errors corrected
	 */
	public int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, LongList kmers, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs){

		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		tracker.clear();
		int valid=fillKmers(bases, kmers);
		if(valid<2){return 0;}

		fillCounts(bases, kmers, counts);
		final int possibleErrors=tracker.suspected=countErrors(counts, quals);
		if(possibleErrors<0){return 0;}
		final float expectedErrors=r.expectedErrors(true, r.length());
		final Rollback roll=ECC_ROLLBACK ? new Rollback(r, counts) : null;

		assert(counts.size>0) : counts.size;

		int correctedPincer=0;
		int correctedTail=0;
		int correctedBrute=0;
		int correctedReassemble=0;

		if(ECC_PINCER){
			correctedPincer+=errorCorrectPincer(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, errorExtensionPincer);
		}

		if(ECC_TAIL || ECC_ALL){
			int start=(ECC_ALL ? 0 : counts.size-k-1);
			correctedTail+=errorCorrectTail(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, start, errorExtensionTail);
			r.reverseComplementFast();
			valid=fillKmers(bases, kmers);
			counts.reverse();
			correctedTail+=errorCorrectTail(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, start, errorExtensionTail);
			r.reverseComplementFast();
			counts.reverse();
		}

		if(ECC_REASSEMBLE){
			if(verbose){System.err.println("Correcting "+possibleErrors+" errors.  Counts:\n"+counts);}
			if((correctedPincer<1 && correctedTail<1) || countErrors(counts, quals)>0){
				correctedReassemble=reassemble(bases, quals, rightCounts, counts, counts2, tracker, errorExtensionReassemble, bb, bb2, null, null, bs);
			}
			if(verbose){System.err.println("Corrected  "+correctedReassemble+" errors.  Counts:\n"+counts);}
		}
		assert(counts.size>0);

		assert(correctedPincer+correctedTail+correctedReassemble+correctedBrute==tracker.corrected())
			: correctedPincer+", "+correctedTail+", "+correctedReassemble+", "+correctedBrute+", "+tracker;

		if(ECC_ROLLBACK && (tracker.corrected()>0 || tracker.rollback)){

			if(!tracker.rollback && quals!=null && tracker.corrected()>3){
				float mult=Tools.max(1, 0.5f*(0.5f+0.01f*r.length()));//1 for a 150bp read.
				if(countErrors(counts, quals)>0 && tracker.corrected()>mult+expectedErrors){tracker.rollback=true;}
				else if(tracker.corrected()>2.5f*mult+expectedErrors){tracker.rollback=true;}
			}

			IntList counts0=roll.counts0;
			for(int i=0; !tracker.rollback && i<counts.size; i++){
				int a=Tools.max(0, counts0.get(i));
				int b=Tools.max(0, counts.get(i));
				if(b<a-1 && !isSimilar(a, b)){
					if(verbose){System.err.println("Y: RID="+r.numericID+"; "+a+"->"+b+"\n"+counts0+"\n"+counts);}
					tracker.rollback=true;
				}
			}

			if(tracker.rollback){
				roll.rollback(r, counts);
				tracker.clearCorrected();
				return 0;
			}
		}

		if(MARK_BAD_BASES>0 && (!MARK_ERROR_READS_ONLY || countErrors(counts, quals)>0 ||
				r.expectedErrors(false, r.length())>3)){
			int marked=markBadBases(bases, quals, counts, bs, MARK_BAD_BASES, MARK_DELTA_ONLY, MARK_QUALITY);
			tracker.marked=marked;
		}

		return tracker.corrected();
	}

	/** Changes to N any base covered strictly by kmers with count below minCount */
	public final int markBadBases(final byte[] bases, final byte[] quals, final IntList counts, final BitSet bs,
			final int minCount, boolean deltaOnly, final byte markQuality){
		if(counts.size<1){return 0;}

		bs.clear();
		assert(counts.size==bases.length-k+1) : counts.size+", "+bases.length;

		for(int i=0; i<counts.size;){
			final int count=counts.get(i);
			if(count>=minCount){
				bs.set(i, i+k);
				i+=k;
			}else{
				i++;
			}
		}
		{//Last cycle
			final int i=counts.size-1;
			final int count=counts.get(i);
			if(count>=minCount){
				bs.set(i, i+k);
			}
		}

		final int card=bs.cardinality();
		final int toMark=bases.length-card;
		int marked=0;
		assert(card<=bases.length);

		int consecutiveBad=0;
		for(int i=0; i<bases.length; i++){
			if(bs.get(i)){
				consecutiveBad=0;
			}else{
				consecutiveBad++;
				boolean mark=((quals!=null && quals[i]>markQuality) || bases[i]!='N');
				if(mark && deltaOnly){
					mark=(consecutiveBad>=k) || bs.get(i+1) || (i>0 && bs.get(i-1));
				}
				if(mark){
					marked++;

					if(markQuality<1){
						bases[i]='N';
					}
					if(quals!=null){
						quals[i]=Tools.min(quals[i], (bases[i]=='N' ? 0 : markQuality));
					}
				}
				if(bases[i]=='N' || (quals!=null && quals[i]<=markQuality)){consecutiveBad=0;}
			}
		}

		return marked;
	}

	/**
	 * Fills k-mer count array for the given sequence.
	 * Uses either small k-mer lookup or big k-mer counting depending on k size.
	 * Applies smoothing if enabled to reduce noise in count profiles.
	 *
	 * @param bases Sequence bases
	 * @param kmers K-mer list (used for small k)
	 * @param counts Output count array
	 */
	public void fillCounts(byte[] bases, LongList kmers, IntList counts){
		counts.clear();

		if(k==ksmall){
			fillCountsFromKmers(kmers, counts);
		}else{
			filter.fillCountsBig(bases, counts);
		}

		if(smooth){
			smooth(bases, kmers, counts, smoothWidth);
		}
	}

	/**
	 * Fills counts by looking up each k-mer individually.
	 * Used when k equals ksmall for direct k-mer lookup.
	 * @param kmers List of encoded k-mers
	 * @param counts Output count array
	 */
	protected void fillCountsFromKmers(LongList kmers, IntList counts){
		for(int i=0; i<kmers.size; i++){
			long kmer=kmers.get(i);
			if(kmer>=0){
				long rkmer=rcomp(kmer);
				int count=getCount(kmer, rkmer);
				counts.add(count);
			}else{
				counts.add(0);
			}
		}
	}

	/**
	 * Smooths k-mer count profiles to reduce noise and artifacts.
	 * Applies multi-pass smoothing to eliminate isolated peaks and valleys
	 * that likely represent sequencing errors rather than real abundance changes.
	 *
	 * @param bases Sequence bases (needed only for the two edge positions, which
	 *   must derive a candidate k-mer rather than merely look one up -- see
	 *   {@link #maxLeftCount})
	 * @param kmerList List of k-mers
	 * @param countList Count values to smooth
	 * @param width Maximum width of peaks to smooth
	 */
	public void smooth(byte[] bases, LongList kmerList, IntList countList, int width){
		final int size=countList.size;
		if(size<3){return;}
		final int[] counts=countList.array;
		final long[] kmers=kmerList.array;

		counts[0]=Tools.min(counts[0], Tools.max(counts[1], maxLeftCount(bases, kmerList, 0)));
		counts[size-1]=Tools.min(counts[size-1], Tools.max(counts[size-2], maxRightCount(bases, kmerList, size-1)));

		for(int i=1, max=size-1; i<max; i++){
			long ka=kmers[i-1], kc=kmers[i+1];
			if(ka>=0 && kc>=0){
				int a=counts[i-1], b=counts[i], c=counts[i+1];
				int maxCount=Tools.max(a, c);
				counts[i]=Tools.min(b, maxCount);
			}
		}

		//Smooth peaks 2-wide.
		if(width>=2){
			for(int i=1, max=size-2; i<max; i++){
				long ka=kmers[i-1], kd=kmers[i+2];
				if(ka>=0 && kd>=0){
					int a=counts[i-1], b=counts[i], c=counts[i+1], d=counts[i+2];
					int maxCount=Tools.max(a, d);
					counts[i]=Tools.min(b, maxCount);
					counts[i+1]=Tools.min(c, maxCount);
				}
			}
		}

		//Smooth peaks 3-wide.
		if(width>=3){
			boolean changed=false;
			for(int i=2, max=size-2; i<max; i++){
				long ka=kmers[i-2], kc=kmers[i+2];
				if(ka>=0 && kc>=0){
					int a=counts[i-2], b=counts[i], c=counts[i+2];
					int maxCount=Tools.max(a, c);
					if(maxCount<b){
						counts[i]=Tools.min(b, maxCount);
						changed=true;
					}
				}
			}
			if(changed){smooth(bases, kmerList, countList, 2);}
		}
	}

	/**
	 * Reassembly-based error correction using bidirectional extension.
	 * Extends sequences from both ends to identify and correct errors
	 * through consensus agreement between forward and reverse extensions.
	 *
	 * @param bases Sequence bases to correct
	 * @param quals Quality scores
	 * @param rightCounts Buffer for right extension counts
	 * @param counts K-mer count values
	 * @param counts2 Secondary count buffer
	 * @param tracker Error tracking statistics
	 * @param errorExtension Extension distance for error detection
	 * @param bb Primary byte buffer
	 * @param bb2 Secondary byte buffer
	 * @param kmer Kmer object for operations
	 * @param regenKmer Secondary Kmer object
	 * @param bs BitSet for marking changed positions
	 * @return Number of errors corrected
	 */
	public final int reassemble(final byte[] bases, final byte[] quals, final int[] rightCounts, final IntList counts, final IntList counts2,
			final ErrorTracker tracker, final int errorExtension, final ByteBuilder bb, final ByteBuilder bb2, final Kmer kmer, final Kmer regenKmer, BitSet bs){
		if(bases.length<k+1+deadZone){return 0;}
		final ByteBuilder fromLeft=new ByteBuilder(bases.length);
		final ByteBuilder fromRight=new ByteBuilder(bases.length);

		int detected0=tracker.detectedReassemble;
		int corrected=reassemble_pass(bases, quals, fromLeft, fromRight, rightCounts, counts, counts2, tracker, errorExtension, kmer, regenKmer, bs);

		int correctedIncr=corrected;
		int detectedIncr=tracker.detectedReassemble-detected0;
		int uncorrected=detectedIncr-correctedIncr;

		for(int passes=1; passes<6 && correctedIncr>0 && uncorrected>0; passes++){//Without a pass limit this could, in rare cases, make an infinite loop
			tracker.detectedReassemble-=uncorrected;
			detected0=tracker.detectedReassemble;
			correctedIncr=reassemble_pass(bases, quals, fromLeft, fromRight, rightCounts, counts, counts2, tracker, errorExtension, kmer, regenKmer, bs);

			corrected+=correctedIncr;
			detectedIncr=tracker.detectedReassemble-detected0;
			uncorrected=detectedIncr-correctedIncr;
		}

		return corrected;
	}

	/**
	 * Single pass of the reassembly algorithm.
	 * Performs bidirectional extension and compares results to identify corrections.
	 * Uses windowing to prevent over-correction in high-error regions.
	 *
	 * @param bases Original sequence bases
	 * @param quals Quality scores
	 * @param fromLeft Extension from left end
	 * @param fromRight Extension from right end
	 * @param rightCounts Buffer for extension counts
	 * @param counts K-mer count values
	 * @param counts2 Secondary count buffer
	 * @param tracker Error tracking statistics
	 * @param errorExtension Extension distance
	 * @param kmer Primary Kmer object
	 * @param kmer2 Secondary Kmer object
	 * @param bs BitSet for marking changes
	 * @return Number of errors corrected
	 */
	public final int reassemble_pass(final byte[] bases, final byte[] quals, final ByteBuilder fromLeft, final ByteBuilder fromRight,
			final int[] rightCounts, final IntList counts, final IntList counts2, final ErrorTracker tracker, final int errorExtension,
			final Kmer kmer, final Kmer kmer2, final BitSet bs){
		if(bases.length<k+1+deadZone){return 0;}

		fromLeft.clear();
		fromRight.clear();
		for(byte b : bases){
			fromLeft.append(b);
			fromRight.append(b);
		}

		assert(counts.size>0) : counts+", "+bases.length;

		counts2.clear();
		counts2.addAll(counts);
		reassemble_inner(fromLeft, quals, rightCounts, counts2, errorExtension, kmer, kmer2);

		fromRight.reverseComplementInPlace();
		counts2.clear();
		counts2.addAll(counts);
		counts2.reverse();

		reassemble_inner(fromRight, quals, rightCounts, counts2, errorExtension, kmer, kmer2);
		fromRight.reverseComplementInPlace();

		int correctedInner=0;
		int correctedOuter=0;
		int detectedInner=0;
		int detectedOuter=0;
		boolean rollback=false;

		for(int i=0; i<bases.length; i++){
			byte a=bases[i];
			byte b=fromLeft.get(i);
			byte c=fromRight.get(i);
			if(a!=b || a!=c){
				if(b==c){detectedInner++;}
				else{
					detectedOuter++;
					if(a!=b && a!=c){
						assert(b!=c);
						rollback=true;
					}
				}
			}
			if(b==a){fromLeft.set(i, (byte)0);}
			if(c==a){fromRight.set(i, (byte)0);}
		}

		final int detected=detectedInner+detectedOuter;
		tracker.detectedReassemble+=detected;
		if(rollback || detected==0){return 0;}
		bs.clear();

		int clearedLeft=clearWindow2(fromLeft, quals, windowLen, windowCount, windowQualSum);
		fromRight.reverseInPlace();
		Vector.reverseInPlace(quals);
		int clearedRight=clearWindow2(fromRight, quals, windowLen, windowCount, windowQualSum);
		fromRight.reverseInPlace();
		Vector.reverseInPlace(quals);

		for(int i=0; i<bases.length; i++){
			byte a=bases[i];
			byte b=fromLeft.get(i);
			byte c=fromRight.get(i);
			byte d=a;
			if(b==0 && c==0){
				//do nothing
			}else if(b==c){
				d=b;
			}else if(b==0){
				d=c;
			}else if(c==0){
				d=b;
			}else if(b!=c){
				//no action
			}

			if(ECC_REQUIRE_BIDIRECTIONAL && b!=c && i>=k && i<bases.length-k){d=a;}//Clause to force pincer mode in the middle

			if(d!=a){
				byte q=(quals==null ? 30 : quals[i]);
				if(b==c){
					correctedInner++;
					q=(byte)Tools.mid(q+qIncreasePincer, qMinPincer, qMaxPincer);
				}else{
					correctedOuter++;
					q=(byte)Tools.mid(q+qIncreaseTail, qMinTail, qMaxTail);
				}
				if(!rollback){
					bs.set(i);
					bases[i]=d;
					if(quals!=null){quals[i]=q;}
				}
			}
		}

		if(rollback && correctedInner+correctedOuter>0){
			tracker.rollback=true;
			return 0;
		}

		{
			tracker.correctedReassembleInner+=correctedInner;
			tracker.correctedReassembleOuter+=correctedOuter;
		}
		int corrected=correctedOuter+correctedInner;

		if(corrected>0){
			regenerateCounts(bases, counts, kmer, bs);
			assert(counts.size>0);
		}

		return corrected;
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
	public abstract int regenerateCounts(byte[] bases, IntList counts, final Kmer dummy, BitSet changed);

	/**
	 * Clears correction candidates in windows with too many changes.
	 * Prevents over-correction by removing proposed changes when too many
	 * corrections are concentrated in a small region.
	 *
	 * @param bb ByteBuilder containing proposed changes
	 * @param quals Quality scores
	 * @param window Window size for analysis
	 * @param limit Maximum changes allowed in window
	 * @param qsumLimit Maximum quality sum allowed in window
	 * @return Number of positions cleared
	 */
	private static int clearWindow2(final ByteBuilder bb, final byte[] quals, final int window,
			final int limit, final int qsumLimit){
		final int len=bb.length;
		final byte[] array=bb.array;

		int cleared=0;
		int count=0, countHQ=0, qsum=0;
		for(int i=0, prev=-window; i<len; i++, prev++){
			byte b=array[i];

			if(b!=0 && (quals==null || quals[i]>0)){
				count++;
				if(quals!=null){
					qsum+=quals[i];
					countHQ++;
				}
				if(count>limit || qsum>qsumLimit){
					for(int j=Tools.max(0, i-window), lim=bb.length(); j<lim; j++){
						if(array[j]!=0){
							array[j]=0;
							cleared++;
						}
					}
					return cleared;
				}
			}
			if(prev>=0 && array[prev]>0 && (quals==null || quals[prev]>0)){
				count--;
				if(quals!=null){
					qsum-=quals[prev];
					countHQ--;
				}
			}
		}
		return cleared;
	}

	/**
	 * Wrapper for reassemble_inner with Kmer parameters.
	 *
	 * @param bb ByteBuilder containing sequence
	 * @param quals Quality scores
	 * @param rightCounts Buffer for right extension counts
	 * @param counts K-mer count values
	 * @param errorExtension Extension distance for error detection
	 * @param kmer Primary Kmer object
	 * @param regenKmer Secondary Kmer object
	 * @return Number of errors corrected
	 */
	public int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension, final Kmer kmer, final Kmer regenKmer){
		return reassemble_inner(bb, quals, rightCounts, counts, errorExtension);
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
	public abstract int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension);

	/**
	 * Determines if a position represents a substitution error.
	 * Analyzes k-mer count patterns around a position to identify
	 * the characteristic signature of single-base substitutions.
	 *
	 * @param ca Count array index for the position
	 * @param errorExtension Extension distance for context analysis
	 * @param qb Quality score at the position
	 * @param counts K-mer count values
	 * @return true if position appears to be a substitution error
	 */
	protected final boolean isSubstitution(int ca, int errorExtension, byte qb, IntList counts){
		final int cb=ca+1;
		final int aCount=counts.get(ca);
		final int bCount=counts.get(cb);
		if(isError(aCount, bCount, qb) && isSimilar(aCount, ca-errorExtension, ca-1, counts) &&
				isError(aCount, ca+2, ca+k, counts)){
			final int cc=ca+k;
			final int cd=cc+1;
			if(cd<counts.size){
				final int cCount=counts.get(cc);
				final int dCount=counts.get(cd);
				if(isError(aCount, dCount) || isError(dCount, cCount, qb)){
					return true;
				}
			}else{return true;}
		}
		return false;
	}

	/**
	 * Counts potential errors in a sequence based on k-mer count transitions.
	 * Identifies positions where k-mer counts suggest sequencing errors
	 * using bidirectional analysis and quality score thresholds.
	 *
	 * @param counts K-mer count values
	 * @param quals Quality scores (may be null)
	 * @return Number of potential errors detected
	 */
	public final int countErrors(IntList counts, byte[] quals){
		int possibleErrors=0;
		for(int i=1; i<counts.size; i++){
			final int a=counts.get(i-1), b=counts.get(i);
			boolean error;
			if(quals!=null){
				error=isErrorBidirectional(a, b, quals[i-1], quals[i+k-1]);
			}else{
				error=isErrorBidirectional(a, b, (byte)20, (byte)20);
			}
			if(error){
				possibleErrors++;
				i+=k;
			}
		}
		return possibleErrors;
	}

	/**
	 * Pincer-based error correction using flanking k-mer analysis.
	 * Identifies errors by comparing count patterns of adjacent k-mers
	 * and corrects single-base substitutions when flanking k-mers agree.
	 *
	 * @param bases Sequence bases to correct
	 * @param quals Quality scores
	 * @param leftBuffer Buffer for left extension analysis
	 * @param rightBuffer Buffer for right extension analysis
	 * @param kmers K-mer list from sequence
	 * @param counts K-mer count values
	 * @param bb ByteBuilder for operations
	 * @param tracker Error tracking statistics
	 * @param errorExtension Extension distance for verification
	 * @return Number of errors corrected
	 */
	public int errorCorrectPincer(final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final ErrorTracker tracker, final int errorExtension){

		int detected=0;
		int corrected=0;

		//a is the index of the left kmer
		//b is a+1 (right-extension of left kmer)
		//c is d-1 (left-extension of right kmer)
		//d is the index of the right kmer
		//the base between the kmers is at a+k
		for(int a=0, d=k+1; d<counts.size; a++, d++){
			final int aCount=counts.get(a);
			final int bCount=counts.get(a+1);
			final int cCount=counts.get(d-1);
			final int dCount=counts.get(d);
			final byte qb=(quals==null ? 20 : quals[a+k]);
			if(isError(aCount, bCount, qb) && isError(dCount, cCount, qb) && isSimilar(aCount, dCount)){
				if(verbose){
					System.err.println("Found error: "+aCount+", "+bCount+", "+cCount+", "+dCount);
				}
				//Looks like a 1bp substitution; attempt to correct.
				detected++;
				int ret=correctSingleBasePincer(a, d, bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, errorExtension);
				corrected+=ret;
				if(verbose){
					System.err.println("Corrected error.");
				}
			}else{
				if(verbose){
					System.err.println("Not an error: "+aCount+", "+bCount+", "+cCount+", "+dCount+
							";  "+isError(aCount, bCount, qb)+", "+isError(dCount, cCount, qb)+", "+isSimilar(aCount, dCount));
				}
			}
		}

		{
			tracker.detectedPincer+=detected;
			tracker.correctedPincer+=corrected;
		}

		return corrected;
	}

	/**
	 * Tail-based error correction for sequence ends.
	 * Corrects errors near sequence ends where pincer correction is not possible.
	 * Uses unidirectional extension to verify and correct substitutions.
	 *
	 * @param bases Sequence bases to correct
	 * @param quals Quality scores
	 * @param leftBuffer Buffer for left extension analysis
	 * @param rightBuffer Buffer for right extension analysis
	 * @param kmers K-mer list from sequence
	 * @param counts K-mer count values
	 * @param bb ByteBuilder for operations
	 * @param tracker Error tracking statistics
	 * @param startPos Starting position for correction
	 * @param errorExtension Extension distance for verification
	 * @return Number of errors corrected
	 */
	public int errorCorrectTail(final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final ErrorTracker tracker, final int startPos, final int errorExtension){
		if(bases.length<k+2+errorExtension+deadZone){return 0;}
		int detected=0;
		int corrected=0;

		//a is the index of the left kmer
		//b is a+1
		//the base between the kmers is at a+k
		for(int a=Tools.max(startPos, errorExtension), lim=counts.size-deadZone-1; a<lim; a++){//errorExtension-1
			final int aCount=counts.get(a);
			final int bCount=counts.get(a+1);
			final byte qb=(quals==null ? 20 : quals[a+k]);
			if(isError(aCount, bCount, qb) && isSimilar(aCount, a-errorExtension, a-1, counts) && isError(aCount, a+2, a+k, counts)){
				if(verbose){
					System.err.println("Found error: "+aCount+", "+bCount);
				}
				//Assume like a 1bp substitution; attempt to correct.
				detected++;
				int ret=correctSingleBaseRight(a, bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, errorExtension);
				corrected+=ret;
				if(verbose){
					System.err.println("Corrected error.");
				}
			}else{
				if(verbose){
					System.err.println("Not an error: "+aCount+", "+bCount+
							";  "+isError(aCount, bCount, qb)+", "+isSimilar(aCount, a-errorExtension, a-1, counts)+", "+isError(aCount, a+2, a+k, counts));
				}
			}
		}

		{
			tracker.detectedTail+=detected;
			tracker.correctedTail+=corrected;
		}

		return corrected;
	}

	/**
	 * Corrects a single base using pincer analysis with bidirectional verification.
	 * Extends from both flanking k-mers to confirm the correct base at the error position.
	 *
	 * @param a Left k-mer index
	 * @param d Right k-mer index
	 * @param bases Sequence bases
	 * @param quals Quality scores
	 * @param leftBuffer Buffer for left extension
	 * @param rightBuffer Buffer for right extension
	 * @param kmers K-mer list
	 * @param counts K-mer count values
	 * @param bb ByteBuilder for extension
	 * @param errorExtension Extension distance required
	 * @return 1 if corrected, 0 if not corrected
	 */
	private int correctSingleBasePincer(final int a, final int d, final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final int errorExtension){
		final byte leftReplacement, rightReplacement;
		final int loc=a+k;
		{
			bb.clear();
			int extension=extendFromStored(bb, null, rightBuffer, errorExtension, true, bases, kmers, a, false);
			if(extension<errorExtension){return 0;}
			for(int i=1; i<extension; i++){
				if(bb.get(i)!=bases[loc+i]){
					return 0;
				}
			}
			leftReplacement=bb.get(0);
		}
		{
			bb.clear();
			int extension=extendFromStored(bb, null, rightBuffer, errorExtension, true, bases, kmers, d, true);
			if(extension<errorExtension){return 0;}
			bb.reverseComplementInPlace();
			for(int i=0; i<extension-1; i++){
				if(bb.get(i)!=bases[loc+i+1-extension]){
					return 0;
				}
			}
			rightReplacement=bb.get(extension-1);
		}
		if(leftReplacement!=rightReplacement){return 0;}
		if(bases[loc]==leftReplacement){return 0;}
		if(!isSimilar(bases, a, leftReplacement, kmers, counts)){return 0;}

		bases[loc]=leftReplacement;
		assert(d==a+k+1);
		regenerateKmers(bases, kmers, counts, a);
		return 1;
	}

	/**
	 * Corrects a single base using rightward extension only.
	 * Used for tail correction where only one direction is available for verification.
	 *
	 * @param a K-mer index for extension start
	 * @param bases Sequence bases
	 * @param quals Quality scores
	 * @param leftBuffer Buffer for left extension (unused)
	 * @param rightBuffer Buffer for right extension
	 * @param kmers K-mer list
	 * @param counts K-mer count values
	 * @param bb ByteBuilder for extension
	 * @param errorExtension0 Extension distance required
	 * @return 1 if corrected, 0 if not corrected
	 */
	private int correctSingleBaseRight(final int a, final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final int errorExtension0){
		final byte leftReplacement;
		final int loc=a+k;
		final int errorExtension=Tools.min(errorExtension0, bases.length-loc);
		{
			bb.clear();
			int extension=extendFromStored(bb, null, rightBuffer, errorExtension, true, bases, kmers, a, false);
			if(extension<errorExtension){return 0;}
			for(int i=1; i<extension; i++){
				if(bb.get(i)!=bases[loc+i]){
					return 0;
				}
			}
			leftReplacement=bb.get(0);
		}

		if(bases[loc]==leftReplacement){return 0;}
		if(!isSimilar(bases, a, leftReplacement, kmers, counts)){return 0;}

		bases[loc]=leftReplacement;
		regenerateKmers(bases, kmers, counts, a);
		return 1;
	}

	/**
	 * Tests if a proposed base change results in similar k-mer count.
	 * Verifies that the new k-mer has a count similar to the original
	 * to avoid introducing artifacts during correction. Takes the underlying
	 * bases (not just the stored kmer list) because the candidate kmer must
	 * be DERIVED (old kmer with one base substituted), not merely looked up --
	 * a k&gt;31 implementation's stored per-position value is a one-way hash
	 * proxy and cannot be un-hashed, so it must rebuild the candidate from
	 * the live bases instead (mirrors assemble/Tadpole2.isSimilar).
	 *
	 * @param bases Sequence bases (unmodified at the position being tested)
	 * @param a K-mer index
	 * @param newBase Proposed replacement base
	 * @param kmers K-mer list
	 * @param counts K-mer count values
	 * @return true if the new k-mer has similar count
	 */
	protected abstract boolean isSimilar(byte[] bases, int a, byte newBase, LongList kmers, IntList counts);

	/**
	 * Tests if all counts in a range are similar to a reference count.
	 *
	 * @param a Reference count value
	 * @param loc1 Start of range (inclusive)
	 * @param loc2 End of range (inclusive)
	 * @param counts Count array
	 * @return true if all counts in range are similar to reference
	 */
	protected final boolean isSimilar(final int a, int loc1, int loc2, final IntList counts){
		loc1=Tools.max(loc1, 0);
		loc2=Tools.min(loc2, counts.size-1);
		for(int i=loc1; i<=loc2; i++){
			if(!isSimilar(a, counts.get(i))){return false;}
		}
		return true;
	}

	/**
	 * Tests if two count values are similar within tolerance.
	 * Uses both absolute and relative thresholds for similarity testing.
	 *
	 * @param a First count value
	 * @param b Second count value
	 * @return true if counts are considered similar
	 */
	protected final boolean isSimilar(final int a, final int b){
		int min=Tools.min(a, b);
		int max=Tools.max(a, b);
		int dif=max-min;
		assert(dif>=0);
		return (dif<pathSimilarityConstant || dif<max*pathSimilarityFraction);
	}

	/**
	 * Tests if all positions in a range represent errors relative to reference count.
	 *
	 * @param a Reference count value
	 * @param loc1 Start of range (inclusive)
	 * @param loc2 End of range (inclusive)
	 * @param counts Count array
	 * @return true if all positions in range are errors
	 */
	protected final boolean isError(final int a, int loc1, int loc2, final IntList counts){
		loc1=Tools.max(loc1, 0);
		loc2=Tools.min(loc2, counts.size-1);
		for(int i=loc1; i<=loc2; i++){
			if(!isError(a, counts.get(i))){return false;}
		}
		return true;
	}

	/**
	 * Bidirectional error test considering both count values and quality scores.
	 *
	 * @param a First count value
	 * @param b Second count value
	 * @param qa Quality score for first position
	 * @param qb Quality score for second position
	 * @return true if counts represent an error transition
	 */
	protected final boolean isErrorBidirectional(final int a, final int b, final byte qa, final byte qb){
		return (a>=b ? isError(a, b, qb) : isError(b, a, qa));
	}

	/**
	 * Tests if count transition represents a sequencing error.
	 * Uses configurable thresholds and multipliers to identify error patterns.
	 *
	 * @param high Higher count value
	 * @param low Lower count value
	 * @return true if transition indicates an error
	 */
	protected final boolean isError(final int high, final int low){
		float em1;
		if(errorPath==1){
			em1=errorMult1;
		}else if(errorPath==2){
			em1=Tools.mid(errorMult1, 4, low*1.6f-3);
		}else{throw new RuntimeException(""+errorPath);}

		return (low*em1<high || (low<=errorLowerConst && high>=Tools.max(minCountCorrect, low*errorMult2)));
	}

	/**
	 * Tests if count transition represents an error, incorporating quality score.
	 * Quality scores modify error detection thresholds for more accurate identification.
	 *
	 * @param high Higher count value
	 * @param low Lower count value
	 * @param q Quality score at the position
	 * @return true if transition indicates an error
	 */
	protected final boolean isError(final int high, final int low, final byte q){
		float em1;
		if(errorPath==1){
			em1=errorMult1*(1+q*errorMultQFactor);
		}else if(errorPath==2){
			if(low<minCountCorrect && high>=minCountCorrect && q<20 && high>2*low){return true;}
			em1=Tools.mid(errorMult1, 4, low*(q<=10 ? 1.6f : 2f)-3);
			em1=em1*(1+q*errorMultQFactor);
		}else{throw new RuntimeException(""+errorPath);}

		return (low*em1<high || (low<=errorLowerConst && high>=Tools.max(minCountCorrect, low*errorMult2)));
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
	public abstract int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase);

	/**
	 * Extends rightward starting from the kmer stored at list position {@code originPos}.
	 * If {@code fromRcomp} is true, starts from the reverse-complement orientation of
	 * that kmer (used when extending backward from a right-side pincer anchor -- see
	 * correctSingleBasePincer's right-side call). Takes {@code bases} and derives the
	 * starting representation itself, rather than accepting a raw pre-rolled kmer,
	 * because a k&gt;31 implementation's stored per-position value is a one-way hash
	 * proxy that cannot be un-hashed back into a rollable kmer (same reasoning as
	 * {@link #isSimilar}; a k&gt;31 kmer also cannot be represented as a single long
	 * at all, so no long-based signature could work here regardless).
	 *
	 * @param bb ByteBuilder to extend into (cleared and populated by this call)
	 * @param leftCounts Buffer for left neighbor analysis (may be null)
	 * @param rightCounts Buffer for right neighbor counts
	 * @param distance Maximum extension distance
	 * @param includeJunctionBase Whether to include base at junction positions
	 * @param bases Sequence bases (unmodified at the position being tested)
	 * @param kmers K-mer list (used by k&lt;=31 implementations for O(1) lookup)
	 * @param originPos Index into kmers/bases of the anchor kmer
	 * @param fromRcomp Whether to start from the reverse-complement orientation
	 * @return Number of bases added to sequence
	 */
	protected abstract int extendFromStored(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance,
			boolean includeJunctionBase, final byte[] bases, final LongList kmers, final int originPos, final boolean fromRcomp);

	/**
	 * Regenerates k-mers and counts after sequence modification.
	 * Updates k-mer list and count array for positions affected by base changes.
	 *
	 * @param bases Modified sequence bases
	 * @param kmers K-mer list to update
	 * @param counts Count array to update
	 * @param a Starting index for regeneration
	 */
	public abstract void regenerateKmers(byte[] bases, LongList kmers, IntList counts, final int a);

	/**
	 * Finds maximum count among all possible left extensions of the k-mer at
	 * list position {@code pos} (spanning bases[pos..pos+k-1]). Tests all four
	 * possible bases that could precede it. Bases-based (not a raw long) for the
	 * same reason as {@link #isSimilar}: this DERIVES a new candidate k-mer
	 * (prepend a base, drop the last), which a stored hash proxy cannot support.
	 * Called only from {@link #smooth} for the sequence's leading edge.
	 *
	 * @param bases Sequence bases
	 * @param kmers K-mer list (used by k&lt;=31 implementations for O(1) lookup)
	 * @param pos Index into kmers/bases of the k-mer being extended
	 * @return Maximum count among left extensions
	 */
	protected abstract int maxLeftCount(byte[] bases, LongList kmers, int pos);

	/**
	 * Finds maximum count among all possible right extensions of the k-mer at
	 * list position {@code pos} (spanning bases[pos..pos+k-1]). Tests all four
	 * possible bases that could follow it. Same bases-based reasoning as
	 * {@link #maxLeftCount}. Called only from {@link #smooth} for the
	 * sequence's trailing edge.
	 *
	 * @param bases Sequence bases
	 * @param kmers K-mer list (used by k&lt;=31 implementations for O(1) lookup)
	 * @param pos Index into kmers/bases of the k-mer being extended
	 * @return Maximum count among right extensions
	 */
	protected abstract int maxRightCount(byte[] bases, LongList kmers, int pos);

	/**
	 * Tests if position represents a junction using both left and right analysis.
	 *
	 * @param rightMax Highest count in right extensions
	 * @param rightSecond Second highest count in right extensions
	 * @param leftMax Highest count in left extensions
	 * @param leftSecond Second highest count in left extensions
	 * @return true if position is a junction
	 */
	protected final boolean isJunction(int rightMax, int rightSecond, int leftMax, int leftSecond){
		if(isJunction(rightMax, rightSecond)){return true;}
		return isJunction(leftMax, leftSecond);
	}

	/**
	 * Tests if count pattern indicates a junction point.
	 * Uses configurable thresholds to identify positions where multiple
	 * high-count paths diverge, indicating assembly complexity.
	 *
	 * @param max Highest count among extensions
	 * @param second Second highest count among extensions
	 * @return true if counts indicate a junction
	 */
	protected final boolean isJunction(int max, int second){
		if(second<1 || second*branchMult1<max || (second<=branchLowerConst && max>=Tools.max(minCountExtend, second*branchMult2))){
			return false;
		}
		if(verbose){outstream.println("Breaking because second-highest was too high:\n" +
				"max="+max+", second="+second+", branchMult1="+branchMult1+"\n" +
				"branchLowerConst="+branchLowerConst+", minCountExtend="+minCountExtend+", branchMult2="+branchMult2+"\n" +
				(second*branchMult1<max)+", "+(second<=branchLowerConst)+", "+(max>=Tools.max(minCountExtend, second*branchMult2)));}
		return true;
	}

	/** Returns number of valid kmers */
	public abstract int fillKmers(byte[] bases, LongList kmers);

	/** Examines kmer counts around the merge borders to ensure the merge was not chimeric */
	public boolean mergeOK(Read merged, int len1, int len2, LongList kmers, final int width, final int thresh, final long highMult){
		final int len=merged.length();
		final int overlap=Tools.min(len, len1+len2-len);
		final byte[] bases=merged.bases;
		if(len<len1+width+1 || len<len2+width+1 || len<k+width){return true;}
		int valid=fillKmers(bases, kmers);
		final int a=len-len2-1; //Base to left of first boundary
		final int b=a+1; //Base to right of first boundary
		final int c=len1-1;
		final int d=c+1;

		final int ak=a-k+1; //kmer to left of first boundary
		final int bk=b; //kmer to right of first boundary
		final int ck=c-k+1; //kmer to left of second boundary
		final int dk=d; //kmer to right of second boundary

		//This is faster since fewer counts are looked up
		if(ak-width>=0 && ak+width<len){
			int min=getCount2(kmers.get(ak));
			for(int i=ak-width+1; i<ak; i++){min=Tools.min(min, getCount2(kmers.get(i)));}
			int min2=getCount2(kmers.get(ak+1));
			for(int i=ak+2; i<=ak+width; i++){min2=Tools.min(min2, getCount2(kmers.get(i)));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(ck-width>=0 && ck+width<len){
			int min=getCount2(kmers.get(ck));
			for(int i=ck-width+1; i<ck; i++){min=Tools.min(min, getCount2(kmers.get(i)));}
			int min2=getCount2(kmers.get(ck+1));
			for(int i=ck+2; i<=ck+width; i++){min2=Tools.min(min2, getCount2(kmers.get(i)));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}

		if(bk-width>=0 && bk+width<len){
			int min=getCount2(kmers.get(bk));
			for(int i=bk+1; i<bk+width+1; i++){min=Tools.min(min, getCount2(kmers.get(i)));}
			int min2=getCount2(kmers.get(bk-1));
			for(int i=bk-width; i<bk-1; i++){min2=Tools.min(min2, getCount2(kmers.get(i)));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(dk-width>=0 && dk+width<len){
			int min=getCount2(kmers.get(dk));
			for(int i=dk+1; i<dk+width+1; i++){min=Tools.min(min, getCount2(kmers.get(i)));}
			int min2=getCount2(kmers.get(dk-1));
			for(int i=dk-width; i<dk-1; i++){min2=Tools.min(min2, getCount2(kmers.get(i)));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		return true;
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
	public abstract int regenerateCounts(byte[] bases, IntList counts, final int ca);

	/** Converts encoded k-mer to text representation */
	protected abstract StringBuilder toText(long kmer);
	/** Computes reverse complement of encoded k-mer */
	protected abstract long rcomp(long kmer);
	/**
	 * Gets count for k-mer using appropriate method based on k size.
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return K-mer count from Bloom filter
	 */
	public abstract int getCount(long kmer, long rkmer);
	/**
	 * Gets count for k-mer key using appropriate method based on k size.
	 * @param key K-mer key value
	 * @return K-mer count from Bloom filter
	 */
	public abstract int getCount(long key);
	/**
	 * Gets count for k-mer with invalid k-mer handling.
	 * Returns 0 for invalid k-mers (negative values).
	 * @param kmer K-mer to look up
	 * @return K-mer count, or 0 if k-mer is invalid
	 */
	public abstract int getCount2(long kmer);
	/**
	 * Converts k-mer pair to lookup key value.
	 * Uses canonical representation if reverse complement mode is enabled.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return Key value for Bloom filter lookup
	 */
	public abstract long toValue(long kmer, long rkmer);

	/*--------------------------------------------------------------*/
	/*----------------       ThreadLocal Temps      ----------------*/
	/*--------------------------------------------------------------*/

	/** Initializes thread-local storage for correction operations */
	protected final void initializeThreadLocals(){
		if(localLeftCounts.get()!=null){return;}
		localLeftCounts.set(new int[4]);
		localRightCounts.set(new int[4]);
		localLongList.set(new LongList());
		localIntList.set(new IntList());
		localIntList2.set(new IntList());
		localByteBuilder.set(new ByteBuilder());
		localByteBuilder2.set(new ByteBuilder());
		localBitSet.set(new BitSet(300));
		localKmer.set(newKmer());
		localKmer2.set(newKmer());
		localTracker.set(new ErrorTracker());
	}

	/**
	 * Constructs a new Kmer instance for thread-local use. Overridable so a
	 * diagnostic harness can force a specific (perWordK,mult) split
	 * independent of the natural getK/getMult derivation -- e.g. exercising
	 * the multi-word carry/combine path at a k value BFC1 can also run
	 * (k&lt;=31), for direct side-by-side comparison. Null override in
	 * production; see BloomFilterCorrector2.FORCE_MULT_SPLIT.
	 */
	protected Kmer newKmer(){
		return new Kmer(k);
	}

	/** Thread-local storage for left extension count arrays */
	protected ThreadLocal<int[]> localLeftCounts=new ThreadLocal<int[]>();
	/** Thread-local storage for right extension count arrays */
	protected ThreadLocal<int[]> localRightCounts=new ThreadLocal<int[]>();
	/** Thread-local storage for k-mer lists */
	protected ThreadLocal<LongList> localLongList=new ThreadLocal<LongList>();
	/** Thread-local storage for count arrays */
	protected ThreadLocal<IntList> localIntList=new ThreadLocal<IntList>();
	/** Thread-local storage for secondary count arrays */
	protected ThreadLocal<IntList> localIntList2=new ThreadLocal<IntList>();
	/** Thread-local storage for primary byte builders */
	protected ThreadLocal<ByteBuilder> localByteBuilder=new ThreadLocal<ByteBuilder>();
	/** Thread-local storage for secondary byte builders */
	protected ThreadLocal<ByteBuilder> localByteBuilder2=new ThreadLocal<ByteBuilder>();
	/** Thread-local storage for position tracking bit sets */
	protected ThreadLocal<BitSet> localBitSet=new ThreadLocal<BitSet>();
	/** Thread-local storage for primary Kmer objects */
	protected ThreadLocal<Kmer> localKmer=new ThreadLocal<Kmer>();
	/** Thread-local storage for secondary Kmer objects */
	protected ThreadLocal<Kmer> localKmer2=new ThreadLocal<Kmer>();
	/** Thread-local storage for error tracking statistics */
	protected ThreadLocal<ErrorTracker> localTracker=new ThreadLocal<ErrorTracker>();

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Enable pincer-based error correction mode */
	protected boolean ECC_PINCER=false;
	/** Enable tail-based error correction mode */
	protected boolean ECC_TAIL=false;
	/** Enable correction across entire sequence length */
	protected boolean ECC_ALL=false;
	/** Enable reassembly-based error correction mode */
	protected boolean ECC_REASSEMBLE=true;
	/** Enable aggressive error correction parameters */
	protected boolean ECC_AGGRESSIVE=false;
	/** Enable conservative error correction parameters */
	protected boolean ECC_CONSERVATIVE=false;
	/** Enable rollback of questionable corrections */
	protected boolean ECC_ROLLBACK=true;
	/** Require bidirectional confirmation for corrections in middle regions */
	protected boolean ECC_REQUIRE_BIDIRECTIONAL=true;

	/** Mark bases as bad if they are completely covered by kmers with a count below this */
	protected int MARK_BAD_BASES=0;
	/** Only mark bad bases that are adjacent to good bases */
	protected boolean MARK_DELTA_ONLY=true;
	/** Only mark bad bases in reads that appear to have errors */
	protected boolean MARK_ERROR_READS_ONLY=true;
	/** Assign this quality score to marked bases */
	protected byte MARK_QUALITY=0;

	/*--------------------------------------------------------------*/

	/** The Bloom filter containing k-mer counts for correction */
	BloomFilter filter;

	/** K-mer length for error correction operations */
	int k=31;
	/** Small k-mer length for efficient lookups */
	int ksmall=31;
	/** Use reverse complement canonical k-mer representation */
	final boolean rcomp=true;

	/** Minimum count required for k-mer extension operations */
	int minCountExtend=2;
	/** Primary multiplier for branch detection thresholds */
	float branchMult1=20;
	/** Secondary multiplier for branch detection at low counts */
	float branchMult2=3;
	/** Constant threshold for low-count branch detection */
	int branchLowerConst=3;

	/** Error detection algorithm variant selection */
	int errorPath=1;
	/** Primary multiplier for error detection thresholds */
	float errorMult1=16;
	/** Secondary multiplier for error detection at low counts */
	float errorMult2=2.6f;
	/** Quality score factor for adjusting error detection sensitivity */
	float errorMultQFactor=0.002f;
	/** Constant threshold for low-count error detection */
	int errorLowerConst=4;//3 seems fine
	/** Minimum count required for error correction operations */
	int minCountCorrect=3;//5 is more conservative...
	/** Gets the minimum count threshold for error correction */
	int minCountCorrect(){return minCountCorrect;}
	/** Constant threshold for path similarity testing */
	int pathSimilarityConstant=3;
	/** Fractional threshold for path similarity testing */
	float pathSimilarityFraction=0.45f;//0.3
	/** Extension distance for reassembly-based error detection */
	int errorExtensionReassemble=3;//default 2; higher is more conservative
	/** Extension distance for pincer-based error detection */
	int errorExtensionPincer=3;//default 5; higher is more conservative
	/** Extension distance for tail-based error detection */
	int errorExtensionTail=8;//default 9; higher is more conservative
	/** Dead zone size at sequence ends where correction is avoided */
	int deadZone=0;
	/** Window length for correction density analysis */
	int windowLen=12;
	/** Maximum corrections allowed within a window */
	int windowCount=6;
	/** Maximum quality sum allowed for corrections within a window */
	int windowQualSum=80;

	/** Quality score increase for pincer-corrected bases */
	byte qIncreasePincer=8;
	/** Minimum quality score for pincer-corrected bases */
	byte qMinPincer=24;
	/** Maximum quality score for pincer-corrected bases */
	byte qMaxPincer=32;

	/** Quality score increase for tail-corrected bases */
	byte qIncreaseTail=4;
	/** Minimum quality score for tail-corrected bases */
	byte qMinTail=20;
	/** Maximum quality score for tail-corrected bases */
	byte qMaxTail=28;

	/** Enable verbose debugging output */
	boolean verbose=false;
	/** Enable additional verbose debugging output */
	boolean verbose2=false;
	/** Enable smoothing of k-mer count profiles */
	boolean smooth=true;
	/** Width parameter for count profile smoothing */
	int smoothWidth=3;

	/** Output stream for debugging and verbose messages */
	PrintStream outstream=System.err;

}
