package bloom;

import java.util.Arrays;
import java.util.BitSet;

import assemble.ErrorTracker;
import dna.AminoAcid;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/**
 * ukmer.Kmer-based (k&gt;31) implementation of Bloom filter-based sequencing
 * error correction. Encodes kmers via {@link ukmer.Kmer}'s multi-word rolling
 * (addRight/addLeft), and stores each position's canonical xor() hash in the
 * per-position LongList plus an independent xor2() in a parallel thread-local
 * list. Both hashes are strand-canonical; the abstract superclass performs
 * stored-position lookups through getCountAt(), overridden here to use both.
 * Anything that must DERIVE a new candidate kmer (not
 * merely look one up) rebuilds a live Kmer from the underlying bases, since
 * xor() is a one-way hash and cannot be un-hashed -- mirrors assemble/Tadpole2,
 * which this was ported from (see the plan's Port recipe).
 *
 * Correctness pins (see BloomFilterCorrectorWrapper's k&gt;31 dispatch and the
 * plan's §The crux): ukmer.Kmer.MASK_CORE must be false on both filter BUILD
 * and QUERY, and both sides must agree on k after {@link ukmer.Kmer#getKbig}
 * normalization -- the Wrapper handles both before this class is ever touched.
 *
 * @author Brian Bushnell, Nepgear
 * @date August 19, 2026
 */
public class BloomFilterCorrector2 extends BloomFilterCorrector {

	public BloomFilterCorrector2(BloomFilter filter_, int k_, int ksmall_) {
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
	 * k-length confound. Not read anywhere in the normal dispatch path.
	 */
	public static int[] FORCE_MULT_SPLIT=null;

	@Override
	protected Kmer newKmer(){
		return FORCE_MULT_SPLIT!=null ? new Kmer(FORCE_MULT_SPLIT[0], FORCE_MULT_SPLIT[1]) : super.newKmer();
	}

	/*--------------------------------------------------------------*/
	/*----------------      Abstract Seam: Fill      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns number of valid kmers. Rolls a live Kmer over bases, storing each
	 * position's xor() hash. Only starts RECORDING once i>=min (mirroring BFC1
	 * exactly): the list has blen-k+1 entries, kmers[j] spans bases[j..j+k-1] --
	 * NOT one entry per base. Everything downstream (counts, errorCorrectPincer/
	 * Tail, etc., all inherited unchanged from the base) assumes this exact
	 * length; recording from i=0 instead of i=min produced a list of length
	 * blen with a bogus leading run of -1s, which desynced counts.size from
	 * bases.length and crashed smooth()'s edge lookup (caught on the first
	 * real k>31 run, ArrayIndexOutOfBounds at maxRightCount).
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
	 * Regenerates k-mers and counts after sequence modification. Since the stored
	 * value is a hash proxy, rebuilds the kmer ending at position `a` fresh from
	 * bases[a..a+k-1] (unmodified by the edit at a+k), then rolls forward.
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
	 * Regenerates k-mer counts for positions changed during reassembly. The Kmer
	 * parameter is always null at BFC1's call site too (vestigial scaffolding,
	 * per the plan's Gotchas) -- ignored here as well, using our own scratch Kmer.
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
	 * Inner reassembly algorithm for single-direction extension. Rolls a live
	 * Kmer rightward over bb's bases (always fresh from position 0 -- no stored
	 * value to invert here, same as BFC1's own fresh roll).
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
						System.err.println("BFC2 base="+(char)base+", num="+num+"; num==rightMax? "+(num==rightMax));
					}

					if(rightMax>=minCountExtend){
						detected++;
						if(num==rightMax){
							detected--;
							if(verbose){System.err.println("BFC2 rejected: num==rightMax");}
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

	/** Derives the initial kmer from bb's own trailing k bases, then extends. */
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
	 * Examines kmer counts around merge borders to ensure the merge was not
	 * chimeric. Pure count lookup through the abstract getCount2 primitive --
	 * no kmer derivation, so this could stay in the shared base, but the plan
	 * lists it under the seam and duplicating a lookup-only loop costs nothing
	 * for correctness. (Amber: this one COULD be hoisted; flagging rather than
	 * silently deviating a third time on the same class of judgment call.)
	 */
	@Override
	public boolean mergeOK(Read merged, int len1, int len2, LongList kmers, final int width, final int thresh, final long highMult){
		final int len=merged.length();
		final byte[] bases=merged.bases;
		if(len<len1+width+1 || len<len2+width+1 || len<k+width){return true;}
		fillKmers(bases, kmers);
		final int a=len-len2-1;
		final int b=a+1;
		final int c=len1-1;
		final int d=c+1;

		final int ak=a-k+1;
		final int bk=b;
		final int ck=c-k+1;
		final int dk=d;

		if(ak-width>=0 && ak+width<len){
			int min=getCountAt(kmers, ak);
			for(int i=ak-width+1; i<ak; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, ak+1);
			for(int i=ak+2; i<=ak+width; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(ck-width>=0 && ck+width<len){
			int min=getCountAt(kmers, ck);
			for(int i=ck-width+1; i<ck; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, ck+1);
			for(int i=ck+2; i<=ck+width; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(bk-width>=0 && bk+width<len){
			int min=getCountAt(kmers, bk);
			for(int i=bk+1; i<bk+width+1; i++){min=Tools.min(min, getCountAt(kmers, i));}
			int min2=getCountAt(kmers, bk-1);
			for(int i=bk-width; i<bk-1; i++){min2=Tools.min(min2, getCountAt(kmers, i));}
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(dk-width>=0 && dk+width<len){
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
	 * Tests if a proposed base change results in similar k-mer count. Rebuilds
	 * the candidate kmer from bases[a..a+k-1] with newBase substituted at the
	 * end (via addRight, which evicts bases[a]) -- the stored hash at position a
	 * cannot be un-hashed to construct this candidate, per the plan's Port
	 * recipe / Tadpole2.isSimilar precedent.
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

	private LongList xor2List(){
		LongList list=localXor2List.get();
		if(list==null){
			list=new LongList();
			localXor2List.set(list);
		}
		return list;
	}

	private int getCount(final Kmer kmer){
		return dualHash() ? filter.filter.read(kmer.xor(), kmer.xor2()) : filter.filter.read(kmer.xor());
	}

	@Override
	protected int getCountAt(final LongList kmers, final int pos){
		final long key=kmers.get(pos);
		return key<0 ? 0 : (dualHash() ? filter.filter.read(key, xor2List().get(pos)) : filter.filter.read(key));
	}

	private boolean dualHash(){return filter!=null && filter.dualHash;}

	/** Not a real kmer to decode (it's a hash proxy) -- formats the hash value itself, debug-only. */
	@Override
	protected StringBuilder toText(long kmer){
		return new StringBuilder("xor:").append(kmer);
	}

	/**
	 * Identity: Kmer.xor() is already strand-canonical (fwd and rcomp give the
	 * same value -- see the plan's §Bridge resolution), so there is no separate
	 * "reverse complement" of a stored hash to compute.
	 */
	@Override
	protected long rcomp(long kmer){
		return kmer;
	}

	/** The two arguments are the primary and secondary hashes for long kmers. */
	@Override
	public int getCount(long kmer, long rkmer){
		return dualHash() ? filter.filter.read(kmer, rkmer) : filter.filter.read(kmer);
	}

	@Override
	public int getCount(long key){
		return filter.filter.read(key);
	}

	@Override
	public int getCount2(long kmer){
		return kmer<0 ? 0 : filter.filter.read(kmer);
	}

	/** Already canonical -- no max(kmer,rkmer) selection needed. */
	@Override
	public long toValue(long kmer, long rkmer){
		return kmer;
	}

	private final ThreadLocal<LongList> localXor2List=new ThreadLocal<LongList>();

}
