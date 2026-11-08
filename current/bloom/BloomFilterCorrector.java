package bloom;

import java.io.PrintStream;
import java.util.BitSet;

import assemble.ErrorTracker;
import assemble.Rollback;
import shared.Tools;
import shared.ErrorRunMarker;
import simd.Vector;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/**
 * Shared substitution-correction and base-marking heuristics over a borrowed count filter.
 * Pincer and tail passes test extensions; reassembly compares directional proposals
 * without changing read length. Counts are probabilistic support, not proof of a base.
 *
 * Holds all kmer-representation-agnostic logic. The kmer-representation primitives
 * (encoding, rolling, and count lookup) are abstract; {@link BloomFilterCorrector1}
 * implements them with bit-packed longs (k&lt;=31) and {@link BloomFilterCorrector2}
 * with {@link ukmer.Kmer} objects (k&gt;31).
 * Stored lists have one position per read window, including undefined-window sentinels;
 * valid-window counts must not be used as dense list bounds. Long implementations can
 * maintain an aligned secondary-hash list on the current thread.
 *
 * The convenience entry points reuse per-instance, per-thread scratch. Explicit-buffer
 * callers own all supplied arrays/lists/trackers and must initialize subclass scratch
 * where required. Calls may mutate read bases and qualities; ordinary return restores
 * temporary orientation changes, but exceptions are not a transactional rollback.
 * Concurrent calls require separate reads/scratch and stable settings/filter contents;
 * the same thread's scratch is not reentrant. This class owns no input/output streams.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public abstract class BloomFilterCorrector{

	/** Independent counters returned by markErrors; mirrors ErrorRunMarker's diagnostic units. */
	public static final class MarkErrorStats{
		/** Non-N bases changed to N, including accepted edge marks. */
		public int marked;
		/** Interior marks attributed to exactly k low-count windows. */
		public int exactKRuns;
		/** Interior marks attributed to exactly k-1 low-count windows. */
		public int exactKMinus1Runs;
		/** Marks attributed to accepted one-sided edge runs. */
		public int edgeRunsMarked;
		/** Qualified exact-length candidates without an accepted interior/edge mark. */
		public int skippedEndRuns;
		/** Qualified isolated runs longer than k that were skipped. */
		public int skippedMergedRuns;
	}

	/**
	 * Retains geometry and a borrowed filter without initializing thread-local scratch.
	 * @param filter_ Count table wrapper, required before table queries; not closed here
	 * @param k_ Positive read-window length supported by the concrete implementation
	 * @param ksmall_ Positive filter-window length, at most k_
	 */
	public BloomFilterCorrector(BloomFilter filter_, int k_, int ksmall_){
		filter=filter_;
		k=k_;
		ksmall=ksmall_;
		assert(ksmall<=k);
	}

	/**
	 * Initializes this thread's scratch and corrects one read in place; does not visit its mate.
	 * @param r Read with non-null bases and optional aligned qualities
	 * @return Retained substitution count; marking is reported separately in the tracker
	 */
	public int errorCorrect(Read r){
		initializeThreadLocals();
		int corrected=errorCorrect(r, localLeftCounts.get(), localRightCounts.get(), localLongList.get(),
				localIntList.get(), localIntList2.get(), localByteBuilder.get(), localByteBuilder2.get(), localTracker.get(), localBitSet.get());
		return corrected;
	}

	/** Marks bounded low-count runs through ErrorRunMarker, mutating bases/available qualities.
	 * Counts may first be smoothed. Fewer than two valid windows returns new zero counters.
	 * Resets this thread's tracker and records marked bases; does not correct or process mates.
	 * @param r Read with non-null bases and optional aligned qualities
	 * @return New counters independent of the reusable thread-local tracker */
	public MarkErrorStats markErrors(Read r){
		initializeThreadLocals();
		final ErrorTracker tracker=localTracker.get();
		tracker.clear();
		final LongList kmers=localLongList.get();
		final IntList counts=localIntList.get();
		final MarkErrorStats out=new MarkErrorStats();
		final int valid=fillKmers(r.bases, kmers);
		if(valid<2){return out;}
		fillCounts(r.bases, kmers, counts);
		final ErrorRunMarker.Result result=ErrorRunMarker.mark(r.bases, r.quality, counts.array, counts.size, k, errorMult1,
				new ErrorRunMarker.ErrorPredicate(){
					@Override
					public boolean isError(final int high, final int low){return BloomFilterCorrector.this.isError(high, low);}
				});
		out.marked=result.marked;
		out.exactKRuns=result.exactKRuns;
		out.exactKMinus1Runs=result.exactKMinus1Runs;
		out.edgeRunsMarked=result.edgeRunsMarked;
		out.skippedEndRuns=result.skippedEndRuns;
		out.skippedMergedRuns=result.skippedMergedRuns;
		tracker.marked=out.marked;
		return out;
	}

	/**
	 * Compatibility overload; the two Kmer arguments are unused by this adapter.
	 * Delegates with the same caller-owned scratch and mutation contract as the shorter overload.
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
	 * @param kmer Unused compatibility argument
	 * @param kmer2 Unused compatibility argument
	 * @return Number of errors corrected
	 */
	public int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, LongList kmers, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs, Kmer kmer, Kmer kmer2){
		return errorCorrect(r, leftCounts, rightCounts, kmers, counts, counts2, bb, bb2, tracker, bs);
	}

	/**
	 * Samples stored-window counts every clamp(k/2,1,9) positions and checks the last.
	 * Negative stored values, low depth or an abrupt sampled transition report suspicion.
	 * Empty lists return false; a negative result does not establish correctness of skipped windows.
	 *
	 * @param kmers Dense stored windows, with any subclass secondary hashes aligned on this thread
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
			int count=getCountAt(kmers, i);
			final int min=Tools.min(count, prev), max=Tools.max(count, prev);
			//BLM213: preserve a saturated high count instead of wrapping it negative.
			if(count<mcc || (i>0 && (isError(max<Integer.MAX_VALUE ? max+1 : max, min-1)))){return true;}
			prev=count;
		}

		long kmer=kmers.get(kmers.size()-1);
		if(kmer<0){return true;}
		int count=getCountAt(kmers, kmers.size()-1);
		final int min=Tools.min(count, prev), max=Tools.max(count, prev);
		//BLM213: the final-window comparison needs the same saturation guard.
		return count<mcc || isError(max<Integer.MAX_VALUE ? max+1 : max, min-1);
	}

	/**
	 * Corrects one read in place, resetting tracker and filling its dense window/count scratch.
	 * Fewer than two valid windows returns zero. Enabled pincer, tail and reassembly passes
	 * update bases/counts and tracker; the tail pass visits both orientations. Reassembly
	 * can change qualities. Optional rollback restores saved read/count state and clears
	 * correction counters; it does not provide exception cleanup or refresh every scratch list.
	 * Optional final marking is excluded from the returned substitution count. Does not process mates.
	 *
	 * @param r Mutable read with non-null bases and optional aligned qualities
	 * @param leftCounts Caller-owned scratch with exactly four entries
	 * @param rightCounts Caller-owned scratch with exactly four entries
	 * @param kmers Reusable dense stored-window list; its contents are scratch, not final output
	 * @param counts Reusable full-window count profile, filled and changed in place
	 * @param counts2 Reusable independent count profile for directional proposals
	 * @param bb Reusable extension builder
	 * @param bb2 Compatibility scratch; reassemble allocates its own exact-size builders
	 * @param tracker Caller-owned counters reset at entry, retained for inspection afterward
	 * @param bs Reusable changed-position/support bit set
	 * @return Retained substitutions, excluding final marking; zero after rollback
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

	/** Marks bases outside the union of windows whose counts reach minCount.
	 * Clears/rebuilds bs as that full union; an empty count list returns before clearing it.
	 * With markQuality below one, selected bases become N. Otherwise bases are retained;
	 * existing qualities are capped at markQuality (or zero for N), never raised.
	 * deltaOnly selects coverage boundaries or a sufficiently long unmarked bad run;
	 * N bases and qualities at or below the cap reset that distance after each position.
	 * @param bases Mutable read bases, with counts.size equal to bases.length-k+1
	 * @param quals Optional aligned mutable qualities
	 * @param counts Read-only dense window counts
	 * @param bs Caller-owned support scratch, replaced when counts is nonempty
	 * @param minCount Minimum supporting count
	 * @param deltaOnly Restrict marking by boundary/distance rules rather than mark every unsupported base
	 * @param markQuality Quality cap; values below one also replace selected bases with N
	 * @return Number of selected marking operations, not necessarily bases whose byte changed */
	public final int markBadBases(final byte[] bases, final byte[] quals, final IntList counts, final BitSet bs,
			final int minCount, boolean deltaOnly, final byte markQuality){
		if(counts.size<1){return 0;}

		bs.clear();
		assert(counts.size==bases.length-k+1) : counts.size+", "+bases.length;

		//BLM202: grow the union through the furthest supported start touching it.
		//Earlier starts cannot extend farther. Starts after that choice were checked
		//low, so the next backward search examines only the newly reached starts.
		for(int i=0; i<counts.size; i++){
			if(counts.get(i)>=minCount){
				final int start=i;
				int end=i+k;
				while(i+1<counts.size){
					final int limit=Tools.min(end, counts.size-1);
					int next=limit;
					while(next>i && counts.get(next)<minCount){next--;}
					if(next==i){i=limit; break;}
					i=limit;
					end=next+k;
				}
				bs.set(start, end);
			}
		}

		final int card=bs.cardinality();
		final int toMark=bases.length-card;
		int marked=0;
		assert(card<=bases.length);
		if(toMark==0){return 0;}

		int consecutiveBad=0;
		for(int i=0; i<bases.length; i++){
			if(bs.get(i)){
				consecutiveBad=0;
				i=bs.nextClearBit(i)-1;//Supported runs neither mutate bases nor accumulate bad-base distance.
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
	 * Clears counts and fills one value per stored window, or uses constituent-small-window counts.
	 * When k equals ksmall, positional dispatch preserves subclass secondary-hash handling.
	 * Applies downward smoothing when enabled. Does not mutate bases or stored keys.
	 *
	 * @param bases Current sequence matching the supplied stored windows
	 * @param kmers Dense stored-window list and any aligned subclass secondary hashes
	 * @param counts Caller-owned count list to clear/fill and optionally smooth
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
	 * Appends positional lookups; does not clear counts. Called by fillCounts when k equals ksmall.
	 * @param kmers List of encoded k-mers
	 * @param counts Caller-owned destination; clear beforehand when replacing a profile
	 */
	protected void fillCountsFromKmers(LongList kmers, IntList counts){
		for(int i=0; i<kmers.size; i++){
			counts.add(getCountAt(kmers, i));
		}
	}

	/** Looks up a stored kmer by position; long-k implementations may carry a parallel fingerprint list. */
	protected int getCountAt(LongList kmers, int pos){return getCount2(kmers.get(pos));}

	/**
	 * Lowers edge and interior peaks in place; never raises valleys or changes bases/keys.
	 * Fewer than three counts is a no-op. Interior checks require valid neighboring stored
	 * keys. Width at least two enables paired peak clipping; width at least three also
	 * clips centers between two-step flanks, then repeats width-two smoothing if changed.
	 *
	 * @param bases Sequence bases (needed only for the two edge positions, which
	 *   must derive a candidate k-mer rather than merely look one up -- see
	 *   {@link #maxLeftCount})
	 * @param kmerList Dense stored windows aligned with counts and subclass secondary hashes
	 * @param countList Mutable populated count profile
	 * @param width Selects the one-, two- and three-position passes; values above three add no pass
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
	 * Reconciles substitutions proposed by leftward and rightward scans on fixed-length copies.
	 * Allocates two exact-size builders regardless of supplied bb/bb2. Runs at most six passes,
	 * continuing only while the last pass corrected something and left detected candidates.
	 * Mutates bases/optional qualities/counts and accumulates tracker counters; does not reset tracker.
	 *
	 * @param bases Mutable original-orientation bases
	 * @param quals Optional aligned mutable qualities
	 * @param rightCounts Scratch with exactly four entries
	 * @param counts Mutable complete profile aligned to bases
	 * @param counts2 Independent scratch profile copied/reversed for proposals
	 * @param tracker Existing counters to accumulate into
	 * @param errorExtension Number of preceding windows checked for consistent anchor depth
	 * @param bb Unused compatibility builder
	 * @param bb2 Unused compatibility builder
	 * @param kmer Compatibility argument forwarded to recount/scan hooks
	 * @param regenKmer Compatibility argument forwarded to scan hooks
	 * @param bs Scratch for accepted changed-base positions
	 * @return Sum of accepted substitutions across passes, not distinct base positions
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
	 * Builds directional substitution proposals, rejects conflicts/dense changes and applies survivors.
	 * Proposal builders become sparse zero/nonzero change masks. Disagreement where both
	 * directions change the same base to different values aborts this pass before mutation.
	 * Counts and bases retain original orientation on normal return; qualities are temporarily
	 * reversed for the reverse proposal scan and right-side density filtering. The proposal
	 * scan restores quality orientation in finally, including when that scan throws.
	 * Agreement gets pincer quality bounds; a surviving one-sided proposal gets tail bounds.
	 * Recounts accepted positions and adds tracker counters. No whole-read transactional rollback here.
	 *
	 * @param bases Mutable sequence, replaced only at accepted proposal positions
	 * @param quals Optional aligned qualities; may be raised/capped at corrected positions
	 * @param fromLeft Reused left buffer with backing capacity exactly bases.length, required by recounting
	 * @param fromRight Reused right buffer with backing capacity exactly bases.length, required by recounting
	 * @param rightCounts Scratch with exactly four entries
	 * @param counts Complete original-orientation count profile, updated for accepted changes
	 * @param counts2 Independent scratch profile for each directional scan
	 * @param tracker Counters to accumulate; not cleared here
	 * @param errorExtension Number of preceding windows used to check anchor depth
	 * @param kmer Compatibility scratch forwarded to scan/recount hooks
	 * @param kmer2 Compatibility scratch forwarded to scan hooks
	 * @param bs Changed-position scratch, cleared only after nonconflicting candidates are found
	 * @return Accepted substitutions in this pass; zero for insufficient length, no candidates or conflicts
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

		//BLM212: qualities must follow reversed bases, then return to the caller's orientation.
		Vector.reverseInPlace(quals);
		try{
			reassemble_inner(fromRight, quals, rightCounts, counts2, errorExtension, kmer, kmer2);
		}finally{
			Vector.reverseInPlace(quals);
		}
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
	 * Recounts windows affected by a nonempty changed-base set, preserving count-list size.
	 * Implementations may also smooth counts; callers must not assume only those slots changed.
	 *
	 * @param bases Modified sequence bases
	 * @param counts Count array to update
	 * @param dummy Compatibility scratch; current implementations use their own rolling state
	 * @param changed Read-only nonempty set of changed base indices
	 * @return Valid windows visited by the recount, not number of changed bases
	 */
	public abstract int regenerateCounts(byte[] bases, IntList counts, final Kmer dummy, BitSet changed);

	/**
	 * Clears the proposal suffix when a sliding density/quality budget is exceeded.
	 * Zero entries mean no proposal. Nonzero entries with positive quality (or absent
	 * qualities) contribute; the current position is tested before evicting i-window,
	 * so the test can span window+1 coordinates. On failure, clears all nonzero entries
	 * from max(0,i-window) to the end and returns immediately, not just a local window.
	 *
	 * @param bb Mutable zero/nonzero proposal mask, logical length bounds the scan
	 * @param quals Optional read-only qualities in the same direction as bb
	 * @param window Nonnegative coordinate distance retained before each eviction
	 * @param limit Maximum contributing proposals before rejection
	 * @param qsumLimit Maximum contributing quality sum; no quality contribution when quals is null
	 * @return Nonzero proposal slots cleared, including zero-quality proposals in the rejected suffix
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
	 * Compatibility adapter; ignores both Kmer arguments and delegates without allocating.
	 *
	 * @param bb Sequence builder whose backing capacity equals its logical length; recounting uses array length
	 * @param quals Optional aligned read-only qualities for this scan direction
	 * @param rightCounts Caller-owned scratch with exactly four entries
	 * @param counts Mutable dense window-count profile aligned to bb
	 * @param errorExtension Number of preceding windows used for the similarity check
	 * @param kmer Unused compatibility argument
	 * @param regenKmer Unused compatibility argument
	 * @return Number of errors corrected
	 */
	public int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension, final Kmer kmer, final Kmer regenKmer){
		return reassemble_inner(bb, quals, rightCounts, counts, errorExtension);
	}

	/**
	 * Scans count drops in one direction, proposing substitutions without extending bb.
	 * Changes bases and affected counts in place, leaving qualities and logical length unchanged.
	 *
	 * @param bb Sequence builder whose backing capacity equals its logical length; recounting uses array length
	 * @param quals Optional aligned read-only qualities in the same direction as bb
	 * @param rightCounts Caller-owned scratch with exactly four entries
	 * @param counts Mutable complete window-count profile aligned to bb
	 * @param errorExtension Number of preceding windows used for the similarity check
	 * @return Number of substitutions made by this directional scan
	 */
	public abstract int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension);

	/**
	 * Tests a depth drop after anchor ca against preceding similar and following low windows.
	 * If a complete post-error anchor exists, also tests its transition; otherwise the
	 * clipped end can qualify. This is a heuristic predicate, not base validation.
	 *
	 * @param ca Anchor-window start, with ca+1 in the populated profile
	 * @param errorExtension Number of preceding windows to compare, clipped at profile start
	 * @param qb Quality of the proposed error base at ca+k
	 * @param counts Read-only complete profile
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
	 * Counts sampled abrupt depth transitions as suspected errors, without changing the profile.
	 * After a hit, skips k subsequent transition positions; the result is not an exact
	 * error count. Uses qualities of the departing/entering bases, or quality twenty.
	 *
	 * @param counts Read-only dense count profile; fewer than two entries returns zero
	 * @param quals Optional full read qualities, aligned to the profile's sequence
	 * @return Nonnegative number of suspected transitions retained by the skip rule
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
	 * Scans flanking anchors separated by k+1 windows and tries their intervening base.
	 * Accepted substitutions update bases, stored windows and counts, but not qualities.
	 * Adds detected/corrected counts to tracker without clearing existing counters.
	 *
	 * @param bases Sequence bases to correct
	 * @param quals Optional aligned read-only qualities for depth-transition tests
	 * @param leftBuffer Compatibility scratch; current pincer helper does not use it
	 * @param rightBuffer Four-entry mutable extension scratch
	 * @param kmers Dense stored windows and aligned subclass state, refreshed after changes
	 * @param counts Mutable full window-count profile
	 * @param bb Reused extension builder, cleared for each side
	 * @param tracker Existing counters to accumulate into
	 * @param errorExtension Positive number of extension bases required on each side
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
	 * Scans from a selected window start and tests rightward extension substitutions.
	 * startPos determines whether only a tail or a broader region is examined. The caller
	 * handles the reverse direction separately. Mutates bases/windows/counts, not qualities,
	 * and accumulates tracker counts. Insufficient read length returns zero unchanged.
	 *
	 * @param bases Sequence bases to correct
	 * @param quals Optional aligned read-only qualities for transition tests
	 * @param leftBuffer Unused compatibility scratch in the current single-base helper
	 * @param rightBuffer Four-entry mutable extension scratch
	 * @param kmers Dense stored windows and aligned subclass state, refreshed after changes
	 * @param counts Mutable full window-count profile
	 * @param bb Reused extension builder
	 * @param tracker Existing counters to accumulate into
	 * @param startPos First anchor-window index to consider, clamped to at least errorExtension
	 * @param errorExtension Positive extension distance and preceding-depth comparison span
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
	 * Replaces base a+k only if both extensions agree, their other bases match, and depth is similar.
	 * Refreshes affected stored windows/counts after a substitution; does not change qualities.
	 *
	 * @param a Left k-mer index
	 * @param d Right k-mer index
	 * @param bases Sequence bases
	 * @param quals Unused compatibility argument
	 * @param leftBuffer Unused compatibility scratch
	 * @param rightBuffer Four-entry mutable extension scratch, reused for both directions
	 * @param kmers K-mer list
	 * @param counts K-mer count values
	 * @param bb ByteBuilder for extension
	 * @param errorExtension Positive number of extension bases required in each direction
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
	 * Replaces base a+k when right extension, subsequent matching bases and depth support it.
	 * Clips the requested extension at the sequence end. Recounts/rebuilds affected
	 * windows on success; does not change qualities or tracker counters itself.
	 *
	 * @param a K-mer index for extension start
	 * @param bases Sequence bases
	 * @param quals Unused compatibility argument
	 * @param leftBuffer Buffer for left extension (unused)
	 * @param rightBuffer Buffer for right extension
	 * @param kmers K-mer list
	 * @param counts K-mer count values
	 * @param bb ByteBuilder for extension
	 * @param errorExtension0 Positive requested extension distance, clipped to available bases
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
	 * @param a Anchor-window start preceding the candidate base at a+k
	 * @param newBase Proposed replacement at a+k
	 * @param kmers Dense stored-window scratch; consulted according to representation
	 * @param counts Read-only profile providing the anchor's original count
	 * @return true if the new k-mer has similar count
	 */
	protected abstract boolean isSimilar(byte[] bases, int a, byte newBase, LongList kmers, IntList counts);

	/**
	 * Tests every count in an inclusive range clipped to populated indices; empty ranges pass.
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
	 * Accepts a strict absolute difference or a strict fraction of the larger nonnegative count.
	 *
	 * @param a First nonnegative count
	 * @param b Second nonnegative count
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
	 * Requires the reference to dominate every count in a clipped inclusive range; empty ranges pass.
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
	 * Orders the two depths and passes the lower-depth side's quality to the error predicate.
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
	 * Applies the configured high/low depth heuristic, not a verified error classification.
	 * errorPath selects fixed or depth-dependent primary multiplication; unsupported values throw.
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
	 * Applies the high/low depth heuristic using quality of the lower-depth candidate.
	 * errorPath two includes a low-quality shortcut; unsupported errorPath values throw.
	 *
	 * @param high Higher count value
	 * @param low Lower count value
	 * @param q Quality at the lower-depth candidate position
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
	 * Derives an anchor from bb's trailing k bases and appends supported rightward extensions.
	 * May append fewer than distance bases; does not modify any separate read/quality object.
	 *
	 * @param bb ByteBuilder containing initial sequence
	 * @param leftCounts Four-entry left-neighbor scratch, or null to omit left checks
	 * @param rightCounts Four-entry right-neighbor scratch
	 * @param distance Nonnegative maximum number of bases to append
	 * @param includeJunctionBase Permit the implementation's terminal junction-base rule
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
	 * Refreshes windows affected by the changed base at a+k, with read-end clipping.
	 * Mutates stored keys/counts and any aligned subclass secondary-hash scratch.
	 *
	 * @param bases Modified sequence bases
	 * @param kmers K-mer list to update
	 * @param counts Count array to update
	 * @param a Anchor-window start immediately preceding the changed base
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
		if(verbose){outstream.println("Breaking because second-highest was too high:\n"+
				"max="+max+", second="+second+", branchMult1="+branchMult1+"\n"+
				"branchLowerConst="+branchLowerConst+", minCountExtend="+minCountExtend+", branchMult2="+branchMult2+"\n"+
				(second*branchMult1<max)+", "+(second<=branchLowerConst)+", "+(max>=Tools.max(minCountExtend, second*branchMult2)));}
		return true;
	}

	/** Fills dense stored positions, including negative sentinels for undefined windows.
	 * A sequence shorter than k may return zero without clearing old scratch; honor the
	 * return value before consuming it. The valid count is not the list's dense length.
	 * @param bases Read-only sequence
	 * @param kmers Caller-owned destination; long variants also fill aligned thread-local hashes
	 * @return Number of valid windows, excluding undefined-window sentinels */
	public abstract int fillKmers(byte[] bases, LongList kmers);

	/** Heuristically checks depth discontinuities around the two original mate boundaries.
	 * A true result means no tested boundary rejected the merge, not proof against chimerism.
	 * Insufficient base context returns true without filling scratch; otherwise each range
	 * must fit dense stored-window positions, not merely base length or valid-window count.
	 * @param merged Read-only merged read with non-null bases
	 * @param len1 Original first-mate length in bases
	 * @param len2 Original second-mate length in bases
	 * @param kmers Reused dense window scratch, replaced when enough context exists
	 * @param width Positive number of neighboring windows examined on each side
	 * @param thresh High-side minimum depth for rejecting an adjacent depth at most one
	 * @param highMult Maximum allowed high/positive-low depth ratio
	 * @return false for a qualifying discontinuity; true when all available checks pass or are skipped */
	public boolean mergeOK(Read merged, int len1, int len2, LongList kmers, final int width, final int thresh, final long highMult){
		final int len=merged.length();
		final int overlap=Tools.min(len, len1+len2-len);
		final byte[] bases=merged.bases;
		if(len<len1+width+1 || len<len2+width+1 || len<k+width){return true;}
		int valid=fillKmers(bases, kmers);
		final int kmerCount=kmers.size();
		final int a=len-len2-1; //Base to left of first boundary
		final int b=a+1; //Base to right of first boundary
		final int c=len1-1;
		final int d=c+1;

		final int ak=a-k+1; //kmer to left of first boundary
		final int bk=b; //kmer to right of first boundary
		final int ck=c-k+1; //kmer to left of second boundary
		final int dk=d; //kmer to right of second boundary

		//BLM204: bound by dense window positions, including undefined windows.
		//The valid-window count and base length are not bounds for stored-position lookups.
		//This is faster since fewer counts are looked up
		if(ak-width>=0 && ak+width<kmerCount){
			int min=getCountAt(kmers, ak);
			for(int i=ak-width+1; i<ak; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, ak+1);
			for(int i=ak+2; i<=ak+width; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(ck-width>=0 && ck+width<kmerCount){
			int min=getCountAt(kmers, ck);
			for(int i=ck-width+1; i<ck; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, ck+1);
			for(int i=ck+2; i<=ck+width; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}

		if(bk-width>=0 && bk+width<kmerCount){
			int min=getCountAt(kmers, bk);
			for(int i=bk+1; i<bk+width+1; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, bk-1);
			for(int i=bk-width; i<bk-1; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(dk-width>=0 && dk+width<kmerCount){
			int min=getCountAt(kmers, dk);
			for(int i=dk+1; i<dk+width+1; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, dk-1);
			for(int i=dk-width; i<dk-1; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		return true;
	}

	/**
	 * Recounts the anchor ca and up to k following windows after changing base ca+k.
	 * Clips at read end, preserves list size, and does not refresh stored keys or smooth counts.
	 *
	 * @param bases Modified sequence bases
	 * @param counts Count array to update
	 * @param ca In-range anchor-window start preceding the changed base
	 * @return Number of valid k-mers processed
	 */
	public abstract int regenerateCounts(byte[] bases, IntList counts, final int ca);

	/** Allocates diagnostic text: DNA for a packed k-mer, or hash text for a long representation. */
	protected abstract StringBuilder toText(long kmer);
	/** Representation hook: packed reverse complement for short k, identity for a canonical hash proxy. */
	protected abstract long rcomp(long kmer);
	/**
	 * Queries the concrete representation; argument meanings differ between implementations.
	 * @param kmer Forward packed k-mer for short k, primary canonical hash for long k
	 * @param rkmer Packed reverse complement for short k, secondary hash for long k
	 * @return Direct or constituent-window count according to the concrete implementation
	 */
	public abstract int getCount(long kmer, long rkmer);
	/**
	 * Compatibility lookup without sentinel handling; long dual-hash callers need a pair or position.
	 * @param key Selected key, or packed big window on the short implementation's small-window path
	 * @return Implementation-specific lookup count
	 */
	public abstract int getCount(long key);
	/**
	 * Compatibility stored-value lookup, returning zero for a negative sentinel.
	 * Short implementations apply strand selection; long single-value lookup cannot recover
	 * the secondary hash. Internal positional consumers should use getCountAt.
	 * @param kmer Stored forward packed value or primary hash, or a negative sentinel
	 * @return Lookup count, or zero for a negative input
	 */
	public abstract int getCount2(long kmer);
	/**
	 * Selects a primary lookup value without querying the table.
	 * Short implementations apply filter strand policy; long implementations retain the primary hash.
	 *
	 * @param kmer Forward packed k-mer or primary canonical hash
	 * @param rkmer Packed reverse complement or secondary canonical hash
	 * @return Key value for Bloom filter lookup
	 */
	public abstract long toValue(long kmer, long rkmer);

	/*--------------------------------------------------------------*/
	/*----------------       ThreadLocal Temps      ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates this instance's scratch once per calling thread, using left counts as the initialized flag.
	 * Buffers are reused across reads and retained for the thread's lifetime; this does not clear them.
	 * Subclasses may extend initialization with their own representation-specific scratch. */
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
	 * @return New rolling k-mer scratch for this corrector's window length
	 */
	protected Kmer newKmer(){
		return new Kmer(k);
	}

	/** Four left-extension counts per thread; also the initialization sentinel. */
	protected ThreadLocal<int[]> localLeftCounts=new ThreadLocal<int[]>();
	/** Four right-extension counts per thread, reused by candidate scans. */
	protected ThreadLocal<int[]> localRightCounts=new ThreadLocal<int[]>();
	/** Dense stored-window values, including negative undefined-window sentinels. */
	protected ThreadLocal<LongList> localLongList=new ThreadLocal<LongList>();
	/** Primary dense window-count profile for the current read. */
	protected ThreadLocal<IntList> localIntList=new ThreadLocal<IntList>();
	/** Independent copy/reversal of the primary profile for directional proposals. */
	protected ThreadLocal<IntList> localIntList2=new ThreadLocal<IntList>();
	/** Reusable extension scratch; reassembly allocates separate exact-size builders. */
	protected ThreadLocal<ByteBuilder> localByteBuilder=new ThreadLocal<ByteBuilder>();
	/** Secondary compatibility builder retained for explicit-buffer adapters. */
	protected ThreadLocal<ByteBuilder> localByteBuilder2=new ThreadLocal<ByteBuilder>();
	/** Reused changed-base or supported-base set, depending on the current operation. */
	protected ThreadLocal<BitSet> localBitSet=new ThreadLocal<BitSet>();
	/** Primary rolling-k-mer scratch available to subclasses. */
	protected ThreadLocal<Kmer> localKmer=new ThreadLocal<Kmer>();
	/** Independent rolling-k-mer scratch available to subclasses. */
	protected ThreadLocal<Kmer> localKmer2=new ThreadLocal<Kmer>();
	/** Mutable per-read counters; convenience entry points clear these before use. */
	protected ThreadLocal<ErrorTracker> localTracker=new ThreadLocal<ErrorTracker>();

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Enable pincer-based error correction mode */
	protected boolean ECC_PINCER=false;
	/** Enable tail-based error correction mode */
	protected boolean ECC_TAIL=false;
	/** Run the tail algorithm from profile start instead of restricting its starting window. */
	protected boolean ECC_ALL=false;
	/** Enable reassembly-based error correction mode */
	protected boolean ECC_REASSEMBLE=true;
	/** Retained compatibility flag; current base and concrete correctors do not consult it. */
	protected boolean ECC_AGGRESSIVE=false;
	/** Retained compatibility flag; current base and concrete correctors do not consult it. */
	protected boolean ECC_CONSERVATIVE=false;
	/** Enable rollback of questionable corrections */
	protected boolean ECC_ROLLBACK=true;
	/** Require bidirectional confirmation for corrections in middle regions */
	protected boolean ECC_REQUIRE_BIDIRECTIONAL=true;

	/** Positive minimum supporting-window count for final marking; zero disables it. */
	protected int MARK_BAD_BASES=0;
	/** Restrict marking to support boundaries or sufficiently long unmarked unsupported runs. */
	protected boolean MARK_DELTA_ONLY=true;
	/** Require suspected count transitions or more than three expected errors before final marking. */
	protected boolean MARK_ERROR_READS_ONLY=true;
	/** Cap marked qualities at this value; values below one also replace selected bases with N. */
	protected byte MARK_QUALITY=0;

	/*--------------------------------------------------------------*/

	/** Borrowed count filter; concurrent correction requires stable table contents and settings. */
	BloomFilter filter;

	/** K-mer length for error correction operations */
	int k=31;
	/** Filter-window length; smaller constituent windows can support a larger correction window. */
	int ksmall=31;
	/** Legacy constant; concrete lookup methods determine strand policy independently. */
	final boolean rcomp=true;

	/** Minimum count required for k-mer extension operations */
	int minCountExtend=2;
	/** Primary multiplier for branch detection thresholds */
	float branchMult1=20;
	/** Secondary multiplier for branch detection at low counts */
	float branchMult2=3;
	/** Constant threshold for low-count branch detection */
	int branchLowerConst=3;

	/** Error-threshold formula selector; only paths 1 and 2 are implemented. */
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
	/** Preceding-window distance checked for consistent anchor depth during reassembly. */
	int errorExtensionReassemble=3;//default 2; higher is more conservative
	/** Extension distance for pincer-based error detection */
	int errorExtensionPincer=3;//default 5; higher is more conservative
	/** Extension distance for tail-based error detection */
	int errorExtensionTail=8;//default 9; higher is more conservative
	/** Terminal bases excluded from directional reassembly candidate scanning. */
	int deadZone=0;
	/** Preceding-coordinate span for proposal density; the current position adds one coordinate. */
	int windowLen=12;
	/** Maximum contributing proposals before clearing the remaining directional proposal suffix. */
	int windowCount=6;
	/** Maximum positive-quality sum for contributing proposals; absent qualities contribute zero. */
	int windowQualSum=80;

	/** Quality increase for agreeing bidirectional reassembly proposals, not the pincer helper. */
	byte qIncreasePincer=8;
	/** Lower quality bound for agreeing bidirectional reassembly proposals. */
	byte qMinPincer=24;
	/** Upper quality bound for agreeing bidirectional reassembly proposals. */
	byte qMaxPincer=32;

	/** Quality increase for accepted one-sided reassembly proposals, not the tail helper. */
	byte qIncreaseTail=4;
	/** Lower quality bound for accepted one-sided reassembly proposals. */
	byte qMinTail=20;
	/** Upper quality bound for accepted one-sided reassembly proposals. */
	byte qMaxTail=28;

	/** Enable verbose debugging output */
	boolean verbose=false;
	/** Enable additional verbose debugging output */
	boolean verbose2=false;
	/** Lower isolated count peaks after profile fills/recounts. */
	boolean smooth=true;
	/** Enable peak-width passes through 1, 2 or 3; values above 3 add no further pass. */
	int smoothWidth=3;

	/** Borrowed diagnostic stream; never closed here. Some diagnostics write directly to System.err. */
	PrintStream outstream=System.err;

}
