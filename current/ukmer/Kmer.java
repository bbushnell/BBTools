package ukmer;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Represents a k-mer (DNA sequence of length k) with support for large k-values.
 * Uses multi-word arrays to handle k-mers longer than 31 bases efficiently.
 * Maintains forward and reverse-complement words and lazily selects a lookup strand.
 * Supports rolling operations for sequence sliding windows and hash-based comparisons.
 * Instances are mutable and thread-confined. Returned arrays are borrowed; direct
 * array1 edits require fillArray2 before another lookup. Keep global packing,
 * canonicalization and mixing settings stable while related keys/tables are in use.
 * Lookup/hash operations require a complete valid window; addRight builds that
 * window, whereas left rolling and in-place edits require it already to exist.
 * @author Brian Bushnell
 * @date Jul 9, 2015
 */
public class Kmer implements Cloneable{

	/*--------------------------------------------------------------*/
	/*----------------     Geometry and Construction    ------------*/
	/*--------------------------------------------------------------*/

	/** Copies geometry and rolling sequence state into independent arrays. */
	public Kmer(Kmer o){
		this(o.k, o.mult, o.kbig, o.hollow);
		setFrom(o);
	}

	/**
	 * Convenience constructor. When PACKED is false, retains the legacy behavior
	 * (symmetric even-split; kbig may round DOWN from kbig_ if it doesn't
	 * divide evenly -- see getKbig). When PACKED is true, uses the packed
	 * layout (leading words full at FULL_WORD_K=32, last word partial) which
	 * represents supported positive lengths exactly, without rounding.
	 */
	public Kmer(int kbig_){
		this(kbig_, false);
	}

	/** Retains the full rolling span, but optionally hashes only its outer 32-base arms. */
	public Kmer(int kbig_, boolean hollow_){
		this(PACKED ? FULL_WORD_K : getK(kbig_),
			PACKED ? getMultPacked(kbig_) : getMult(kbig_),
			PACKED ? kbig_ : getK(kbig_)*getMult(kbig_), hollow_);
	}

	/** Copies geometry, rolling words and valid-run length into independent arrays.
	 * @return New key whose lookup caches will be computed lazily */
	@Override
	public Kmer clone(){
		return new Kmer(this);
	}

	/**
	 * Legacy 2-arg constructor -- unchanged public behavior (symmetric, kbig
	 * derived as k_*mult_). Still used directly by test/diagnostic code that
	 * wants explicit control over the split (e.g.
	 * BloomFilterCorrector2.FORCE_MULT_SPLIT), independent of PACKED.
	 */
	public Kmer(int k_, int mult_){
		this(k_, mult_, k_*mult_);
	}

	/**
	 * Master constructor. kbig_ is taken as given, NOT recomputed as k_*mult_
	 * -- required so the packed convenience constructor above can pass the
	 * true (possibly non-uniform-word-implied) kbig. lastWordK/lastShift2/
	 * lastMask describe word[maxindex] (in array1's indexing -- the word
	 * addRight/addLeft roll into/out of first), which is the only word that
	 * can differ from k_ under packing; the formula collapses to exactly k_
	 * whenever kbig_==k_*mult_ (the old symmetric invariant), so it is a
	 * no-op change for every existing (non-packed) caller.
	 */
	private Kmer(int k_, int mult_, int kbig_){
		this(k_, mult_, kbig_, false);
	}

	/** Instance-local hash geometry prevents a hollow invocation changing later solid callers. */
	private Kmer(int k_, int mult_, int kbig_, boolean hollow_){
		if(hollow_ && (k_!=FULL_WORD_K || kbig_<=64)){
			throw new IllegalArgumentException("Hollow kmers require packed storage and a span greater than 64 bases.");
		}
		hollow=hollow_;
		k=k_;
		mult=mult_;
		maxindex=mult-1;
		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));

		kbig=kbig_;
		array1=new long[mult];
		array2=new long[mult];
		key=null;

		lastWordK=kbig-k*(mult-1);
		assert(lastWordK>0 && lastWordK<=k) : "lastWordK="+lastWordK+", k="+k+", mult="+mult+", kbig="+kbig;
		final int ls2=2*lastWordK-2;
		lastShift2=ls2;
		lastMask=(ls2+2>63 ? -1L : ~((-1L)<<(ls2+2)));

		//Item 1b (2026-08-19): coreMask deliberately IGNORES the two GLOBAL
		//ends of the kmer (the 5'-most base of word[0] and the 3'-most base
		//of word[maxindex]) so that, in explicit-key/assembly mode, all four
		//{A,C,T,G} extensions of a kmer hash/collide to the same table slot
		//(resolved by linear probing -- 1 random access instead of 4). Under
		//the old uniform per-word coreMask, EVERY word's ends got masked,
		//which was correct only because every word WAS an end (k<=31,
		//single word). With multi-word packing that's wrong -- only the
		//truly-first and truly-last words have a "global end"; middle words
		//are unmasked (-1L). word[0]'s width is perWordK(0), not k, so
		//mult==1 (word[0]==word[maxindex]) is handled correctly by
		//perWordK's own i==maxindex check, not a special case here.
		firstWordCoreMask=(MASK_CORE ? ~(3L<<(2*perWordK(0)-2)) : -1L);
		lastWordCoreMask=(MASK_CORE ? ~3L : -1L);
	}

	/** Effective coreMask for word index i: the global-first-base mask if
	 * i==0, the global-last-base mask if i==maxindex (both, ANDed together,
	 * when mult==1 and they're the same word), -1L (no masking) otherwise. */
	private long coreMaskFor(int i){
		long m=-1L;
		if(i==0){m&=firstWordCoreMask;}
		if(i==maxindex){m&=lastWordCoreMask;}
		return m;
	}

	/** Returns the ceiling number of32-base words; requires a positive supported span. */
	private static int getMultPacked(int kbig){
		return (kbig+FULL_WORD_K-1)/FULL_WORD_K;
	}

	/** Legacy single-word mask excluding both end bases; intended for widths1..31.
	 * For32, Java's shift wrapping makes this return zero when MASK_CORE is enabled. */
	public static final long toCoreMask(int k){
//		System.err.println(k+", "+MASK_CORE);
		return MASK_CORE ? ((~((-1L)<<(2*k)))>>2)&~(3L) : -1L;
	}

	/** Returns the number of words for a positive requested span under current PACKED. */
	public static int getMult(int kbig){
		//PACKED-aware (Amber, 2026-08-19): same bug class as getKbig -- any
		//external code deriving mult/k from these (KmerBufferU, KmerTableSetU)
		//must see the packed values under PACKED, or it silently builds
		//exact-table structures sized for the OLD symmetric mult while real
		//Kmers use the packed one. Diverges from old getMult0 at kbig>=63
		//(e.g. getMult0(64)=3, packed=2) -- k=51 in this session's testing
		//happened not to diverge, which is why it wasn't caught until now.
		if(PACKED){return getMultPacked(kbig);}
		final int mult=getMult0(kbig);
//		assert(mult==getMult0(kbig*(mult/mult))) : mult+", "+getMult0(mult*(kbig/mult));
		return mult;
	}

	/** Returns the represented length; legacy equal-width packing may round down. */
	public static int getKbig(int kbig){
		//Under packing every kbig>0 is exactly representable -- no rounding.
		//Callers like BloomFilterCorrectorWrapper use this to decide the
		//ACTUAL k a request normalizes to, so it must be PACKED-aware or a
		//packed=t run would still silently round down to the old symmetric
		//value (found 2026-08-19 testing k=51 under packed=t).
		if(PACKED){return kbig;}
		int x=getMult(kbig)*getK(kbig);
		assert(x<=kbig) : x+", "+kbig;
		assert(kbig>31 || x==kbig);
		return x;
	}

	/** Chooses the equal-width legacy layout retaining the most requested bases. */
	private static int getMult0(int kbig){
//		if(true){return 2;}//TODO: 123 //Enable to allow multi-word arrays for k<32
		final int word=31;

		final int mult1=(kbig+word-1)/word;
		final int mult2=Tools.max(1, kbig/word);
		if(mult1==mult2){return mult1;}

		final int k1=Tools.min(word, kbig/mult1);
		final int k2=Tools.min(word, kbig/mult2);

		final int kbig1=k1*mult1;
		final int kbig2=k2*mult2;

		assert(kbig1<=kbig);
		assert(kbig2<=kbig);
		assert(mult2<=mult1);

//		assert(false) : mult1+", "+mult2+", "+k1+", "+k2;

		final int mult=kbig2>=kbig1 ? mult2 : mult1;

		return mult;
	}

	/** Returns leading-word width, which is32 in packed mode even for a short final word. */
	public static int getK(int kbig){
		//PACKED-aware, matching what the constructor actually sets `.k` to
		//(FULL_WORD_K -- the leading-word width, NOT kbig/mult, which would
		//be wrong whenever lastWordK!=FULL_WORD_K, i.e. whenever kbig isn't
		//an exact multiple of 32). Same fix reasoning as getMult above.
		if(PACKED){return FULL_WORD_K;}
		int mult=getMult(kbig);
		int x=kbig/mult;
		assert(x*mult<=kbig) : x+", "+kbig;
		assert(x<=31) : kbig+", "+mult+", "+x;
		return x;
	}

	/**
	 * Returns word i's base-count under whichever layout is active. Extended
	 * (Nepgear, 2026-08-19) from Amber's skeleton: word[maxindex] uses
	 * lastWordK (== k whenever kbig==k*mult, i.e. always under the old
	 * symmetric layout, so this is a no-op there); every other word is k.
	 */
	public int perWordK(int i){
		return i==maxindex ? lastWordK : k;
	}

	/**
	 * Static form for dump/table code that has raw kbig/mult without a live
	 * Kmer instance (see ukmer/KmerTableSetU/HashArrayU dump paths). Mirrors
	 * the instance form's logic without a Kmer to read lastWordK from.
	 */
	public static int perWordK(int kbig, int mult, int i){
		if(!PACKED){return kbig/mult;}
		return i==mult-1 ? kbig-FULL_WORD_K*(mult-1) : FULL_WORD_K;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Rolling State and Edits     ------------*/
	/*--------------------------------------------------------------*/

	/** Copies rolling arrays and valid-run length from a key with identical geometry.
	 * Invalidates cached lookup state and retains this object's own arrays/masks. */
	public Kmer setFrom(Kmer o){
		assert(hollow==o.hollow && kbig==o.kbig && k==o.k) :
			"Copying rolling state must preserve span, packing and hollow hash semantics.";
		for(int i=0; i<mult; i++){
			array1[i]=o.array1[i];
			array2[i]=o.array2[i];
			len=o.len;
		}
		incarnation++;
		return this;
	}

	/** Copies mult forward words, rebuilds the reverse strand, and marks a complete window.
	 * The source must supply correctly packed words for this instance's geometry. */
	public Kmer setFrom(long[] array){
		for(int i=0; i<mult; i++){
			array1[i]=array[i];
		}
		fillArray2();
		incarnation++;
		return this;
	}

	/** Resets valid-run length and caches, and zeros both arrays for a fresh rolling window. */
	public void clear(){
		len=0;
		for(int i=0; i<mult; i++){
			array1[i]=0;
			array2[i]=0;
		}
		lastIncarnation=-1;
		lastXor2Incarnation=-1;
		incarnation=0;
		//incarnation++;
	}

	/** Resets valid-run length/caches without clearing words; roll in a full window before lookup. */
	public void clearFast(){
		len=0;
		lastIncarnation=-1;
		lastXor2Incarnation=-1;
		incarnation=0;
		//incarnation++;
	}

	/** Checks complementary sequence content and, under assertions, cache freshness.
	 * @param update Refresh caches first; with assertions enabled this requires a complete window.
	 * If false, a complete window must already have current caches.
	 * @return After optional refresh, true for an incomplete window or matching strands */
	public boolean verify(boolean update){
//		boolean b=verify();
//		if(b){
//			if(update){update();}
//			b=verify();
//			assert(len<kbig || incarnation==lastIncarnation);
//		}
		if(update){
			update();
			assert(len<kbig || incarnation==lastIncarnation) : "incarnation="+incarnation+", last="+lastIncarnation+", len="+len+", kbig="+kbig;
		}
		boolean b=verify();
		return b;
	}

	/** Diagnostic sequence check; assertions also require current cached state for full windows. */
	private boolean verify(){
		if(len<kbig){return true;}
		//Brian's design (2026-08-19): array2 is F's own reverse-complement,
		//independently packed leading-full/partial-last -- NOT a per-word
		//mirror of array1, so there is no fixed word-index relationship
		//between the two arrays to check word-by-word once mult>1 and the
		//last word is partial. Decode both to flat base sequences and compare
		//at the base level; verify() is a diagnostic method (assertion-gated
		//call sites), not hot-path, so the O(kbig) cost here is fine.
		if(!forwardSequence().equals(AminoAcid.reverseComplementBases(reverseSequence()))){
//			assert(false) : toString()+" vs rc of "+reverseSequence();
			return false;
		}
		assert(incarnation==lastIncarnation) : "incarnation="+incarnation+", last="+lastIncarnation+", len="+len+", kbig="+kbig;
		return true;
	}

	/** Decodes array2 (the reverse-complement strand) to a base-sequence
	 * String, the same way toString() decodes array1. Diagnostic-path only. */
	private String reverseSequence(){
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<mult; i++){
			bb.appendKmer(array2[i], perWordK(i));
		}
		return bb.toString();
	}

	/** Rolls in a table-supported base from the right and returns the evicted base letter.
	 * A base mapped to a negative code resets the valid run and inserts an A placeholder. */
	public byte addRight(final byte b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addRightNumeric(x)];
	}

	/** Character form of right rolling; the character must index the alphabet lookup table. */
	public byte addRight(final char b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addRightNumeric(x)];
	}

	/** Rolls a defined base into the left of an already complete window; returns the evicted base. */
	public byte addLeft(final byte b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addLeftNumeric(x)];
	}

	/** Item 3c (Fix #3, Noelle's perf review, Amber directing): substitutes the base at
	 * position `pos` (0-indexed from array1's own 5' end -- word[0]'s first base is
	 * position 0, same convention forwardSequence()/toString() decode) with `newBase`
	 * (numeric, 0-3), updating BOTH array1 (F) and array2 (R=RC(F), at the mirrored
	 * position kbig-1-pos with complement(newBase), same leading-full/partial-last
	 * packing convention as array1 per this class's own doc) in place -- no O(kbig)
	 * rebuild via repeated addRight calls. Only the single affected word in each array
	 * is touched; word index and within-word bit offset are computed directly since
	 * every word except word[maxindex] holds exactly k=FULL_WORD_K bases under PACKED.
	 * <p>
	 * Bumps incarnation so the next key()/xor() call's update() recomputes lazily --
	 * same invalidation contract as addRightNumeric/addLeftNumeric/fillArray2/setFrom.
	 * len is untouched (this pokes an already-built kmer; it is not a roll/extend).
	 * <p>
	 * Substitution is its own inverse: to restore the original base, call this again
	 * with the value this call returns. The caller MUST restore before reusing this
	 * Kmer for anything else that assumes the original content -- same "temporarily
	 * mutated" contract as any other in-place Kmer op.
	 * @param pos Absolute base position, 0-indexed from the 5' end (array1's word[0] start)
	 * @param newBase Numeric base (0-3) to place at that position
	 * @return The numeric base (0-3) previously at that position, for restoration */
	public long substituteBase(final int pos, final long newBase){
		assert(pos>=0 && pos<kbig) : "pos="+pos+", kbig="+kbig;
		assert(newBase>=0 && newBase<4) : newBase;
		assert(len>=kbig) : "Kmer must be fully built before an in-place substitution; len="+len+", kbig="+kbig;

		final int wi=pos/k;
		final int w=perWordK(wi);
		final int withinWordIndex=pos-wi*k;
		final int shiftAmt=2*(w-1-withinWordIndex);
		final long oldBase=(array1[wi]>>>shiftAmt)&3L;
		array1[wi]=(array1[wi]&~(3L<<shiftAmt))|(newBase<<shiftAmt);

		final int qPos=kbig-1-pos;
		final int wi2=qPos/k;
		final int w2=perWordK(wi2);
		final int withinWordIndex2=qPos-wi2*k;
		final int shiftAmt2=2*(w2-1-withinWordIndex2);
		final long newComp=AminoAcid.numberToComplement[(int)newBase];
		array2[wi2]=(array2[wi2]&~(3L<<shiftAmt2))|(newComp<<shiftAmt2);

		incarnation++;
		return oldBase;
	}

	/** Delete the base at pos, shift the suffix left, and append incomingBase.
	 * Fixed kbig and len; returns the deleted numeric base for insertBase restore.
	 * Updates F and RC(F) in O(affected words), without allocation or base replay.
	 * Numeric A/C/G/T=0/1/2/3 (AminoAcid.baseToNumber), so complement(x)=3-x. */
	public long deleteBase(final int pos, final long incomingBase){
		assert(pos>=0 && pos<kbig && incomingBase>=0 && incomingBase<4 && len>=kbig) :
			"Fixed-window deletion needs a complete kmer, valid position and incoming base: pos="+pos+", base="+incomingBase+", len="+len+", K="+kbig;
		final long removed=deleteAndAppend(array1, pos, incomingBase);
		final long reverseRemoved=deleteAndPrepend(array2, kbig-1-pos, 3-incomingBase);
		assert(reverseRemoved==3-removed) : "Deleting mirrored complementary bases must preserve array2=RC(array1).";
		incarnation++;
		return removed;
	}

	/** Insert before pos, shift the suffix right, and drop the final base.
	 * Fixed kbig and len; returns the dropped numeric base for deleteBase restore.
	 * Like substituteBase, invalidates canonical/hash caches lazily. */
	public long insertBase(final int pos, final long newBase){
		assert(pos>=0 && pos<kbig && newBase>=0 && newBase<4 && len>=kbig) :
			"Fixed-window insertion needs a complete kmer, valid position and base: pos="+pos+", base="+newBase+", len="+len+", K="+kbig;
		final long dropped=insertAndDropLast(array1, pos, newBase);
		final long reverseDropped=dropFirstAndInsert(array2, kbig-1-pos, 3-newBase);
		assert(reverseDropped==3-dropped) : "Dropping opposite complementary ends must preserve array2=RC(array1).";
		incarnation++;
		return dropped;
	}

	/** Deletes one packed base, shifts its suffix and appends a base across word boundaries. */
	private long deleteAndAppend(final long[] a, final int pos, final long incoming){
		assert(a.length==mult) : "Word shifts must use this kmer's leading-full/partial-last geometry.";
		final int wi=pos/k, bit=2*(perWordK(wi)-1-pos%k);
		final long old=a[wi], removed=(old>>>bit)&3;
		long carry=incoming;
		for(int i=maxindex; i>wi; i--){
			final long word=a[i];
			a[i]=((word<<2)|carry)&(i==maxindex ? lastMask : mask);
			carry=(word>>>(i==maxindex ? lastShift2 : shift2))&3;
		}
		// Inclusive suffix mask has 2..64 bits; never shift a Java long by 64.
		final long suffix=-1L>>>(62-bit);
		a[wi]=(old&~suffix)|(((old<<2)|carry)&suffix);
		return removed;
	}

	/** Reverse-strand counterpart: deletes a base, shifts its prefix and prepends a base. */
	private long deleteAndPrepend(final long[] a, final int pos, final long incoming){
		assert(a.length==mult) : "RC prefix shifts use the same word widths as F, not reversed widths.";
		final int wi=pos/k, top=2*perWordK(wi)-2, bit=top-2*(pos%k);
		final long old=a[wi], removed=(old>>>bit)&3;
		long carry=incoming;
		for(int i=0; i<wi; i++){
			final long word=a[i];
			a[i]=(word>>>2)|(carry<<shift2);
			carry=word&3;
		}
		final long suffix=(1L<<bit)-1;
		a[wi]=(old&suffix)|(((old>>>2)|(carry<<top))&~suffix);
		return removed;
	}

	/** Inserts at pos while shifting its suffix right; returns the evicted final base. */
	private long insertAndDropLast(final long[] a, final int pos, final long inserted){
		assert(a.length==mult) : "Insertion carries propagate from the changed word into the partial last word.";
		final int wi=pos/k, bit=2*(perWordK(wi)-1-pos%k);
		final long old=a[wi], inclusive=-1L>>>(62-bit), suffix=(1L<<bit)-1;
		a[wi]=(old&~inclusive)|(inserted<<bit)|((old>>>2)&suffix);
		long carry=old&3;
		for(int i=wi+1; i<=maxindex; i++){
			final long word=a[i];
			a[i]=(word>>>2)|(carry<<(i==maxindex ? lastShift2 : shift2));
			carry=word&3;
		}
		return carry;
	}

	/** Drops the global first base and inserts at pos while preserving its unchanged suffix. */
	private long dropFirstAndInsert(final long[] a, final int pos, final long inserted){
		assert(a.length==mult) : "Reverse insertion drops the global first base and retains the unchanged suffix.";
		final int wi=pos/k, top=2*perWordK(wi)-2, bit=top-2*(pos%k);
		final long old=a[wi], inclusive=-1L>>>(62-bit), suffix=(1L<<bit)-1;
		a[wi]=(old&suffix)|((old<<2)&~inclusive&(wi==maxindex ? lastMask : mask))|(inserted<<bit);
		long carry=(old>>>top)&3;
		for(int i=wi-1; i>=0; i--){
			final long word=a[i];
			a[i]=((word<<2)|carry)&mask;
			carry=(word>>>shift2)&3;
		}
		return carry;
	}

	/**
	 * Brian's design (2026-08-19): array2 is NOT a per-word mirror of
	 * array1 -- it is F's reverse-complement, maintained as its OWN
	 * independent sequence in the SAME leading-full/partial-last convention
	 * as array1. Extending F at its 3' end (this method) is equivalent to
	 * prepending comp(x) at R's 5' end (RC(F+x) = comp(x)+RC(F)) -- so F
	 * rolls right-entry (unchanged from before packing) while R rolls
	 * left-entry, evicting from word[maxindex] (R's own partial word) on the
	 * far end. Because both arrays share the same per-word width pattern,
	 * this is still O(mult), not O(kbig).
	 * @param x Numeric base0..3; any negative value resets len and rolls an A placeholder
	 * @return Evicted numeric base from the forward array's left end
	 */
	public long addRightNumeric(long x){
		long x2;

		if(x<0){
			x=0;
			x2=3;
			len=0;
		}else{
			x2=AminoAcid.numberToComplement[(int)x];
			len++;
		}

		long carryF=x;
		for(int i=maxindex; i>=0; i--){
			final int s2=(i==maxindex ? lastShift2 : shift2);
			final long m=(i==maxindex ? lastMask : mask);
			final long evicted=(array1[i]>>>s2)&3L;
			array1[i]=((array1[i]<<2)|carryF)&m;
			carryF=evicted;
		}

		rollRLeft(x2);

		incarnation++;
		return carryF;
	}

	/** See addRightNumeric's doc. F extends at its 5' end here, which is
	 * equivalent to appending comp(x) at R's 3' end (RC(x+F) = RC(F)+comp(x))
	 * -- R rolls right-entry, the mirror image of addRightNumeric's case.
	 * Requires a complete window and x in0..3. Does not change len.
	 * @return Evicted numeric base from the forward array's right end */
	public long addLeftNumeric(long x){
		assert(x>=0 && x<4) : x;
		long x2=AminoAcid.numberToComplement[(int)x];

		assert(x>=0);
		assert(len>=kbig);

		long carryF=x;
		for(int i=0; i<=maxindex; i++){
			final int s2=(i==maxindex ? lastShift2 : shift2);
			final long evicted=array1[i]&3L;
			array1[i]=(array1[i]>>>2)|(carryF<<s2);
			carryF=evicted;
		}

		rollRRight(x2);

		incarnation++;
		return carryF;
	}

	/** R gains base x2 at its 5' end (word[0]/MSB), carry ripples toward the
	 * partial word[maxindex] where the oldest base is evicted and discarded.
	 * Shared by addRightNumeric and fillArray2 (see fillArray2's derivation). */
	private void rollRLeft(long x2){
		long carryR=x2;
		for(int i=0; i<=maxindex; i++){
			final int s2=(i==maxindex ? lastShift2 : shift2);
			final long evicted=array2[i]&3L;
			array2[i]=(array2[i]>>>2)|(carryR<<s2);
			carryR=evicted;
		}
	}

	/** R gains base x2 at its 3' end (word[maxindex]/LSB, the partial word),
	 * carry ripples toward word[0] where the oldest base is evicted. */
	private void rollRRight(long x2){
		long carryR=x2;
		for(int i=maxindex; i>=0; i--){
			final int s2=(i==maxindex ? lastShift2 : shift2);
			final long m=(i==maxindex ? lastMask : mask);
			final long evicted=(array2[i]>>>s2)&3L;
			array2[i]=((array2[i]<<2)|carryR)&m;
			carryR=evicted;
		}
	}

	/**
	 * Rebuilds R(array2) from scratch as F(array1)'s reverse-complement, in
	 * R's own leading-full/partial-last packing -- NOT the old per-word
	 * mirror derivation, which produced a REVERSED width pattern (proven
	 * wrong via Amber's canonical-invariance test: array2 came out
	 * [partial,...,full] against array1's [full,...,partial], so a kmer and
	 * its own reverse-complement disagreed on key()/xor()).
	 * <p>
	 * Derivation: feeding F's bases through {@link #rollRLeft(long)} (the
	 * same left-entry R-update addRightNumeric uses) in F's own 5'-&gt;3'
	 * order correctly reconstructs R=RC(F). Worked by hand on F="AB":
	 * feeding comp(A) then comp(B) through a left-entry roll yields R
	 * (read 5'-&gt;3') = comp(B)+comp(A), which is exactly RC("AB").
	 * Marks the window complete and invalidates lookup caches. Callers must supply
	 * correctly packed forward words; this method does not validate their unused bits.
	 */
	public void fillArray2(){
		Arrays.fill(array2, 0L);
		for(int i=0; i<mult; i++){
			final int pk=perWordK(i);
			for(int b=pk-1; b>=0; b--){
				final long base=(array1[i]>>>(2*b))&3L;
				rollRLeft(AminoAcid.numberToComplement[(int)base]);
			}
		}
		len=kbig;
		incarnation++;
	}

	/**
	 * Returns string representation of the k-mer sequence.
	 * Shows the forward orientation as concatenated bases from all words.
	 * Under assertions, first refreshes and verifies a complete window. Otherwise
	 * decodes the stored words directly, including any zero-filled positions.
	 * @return DNA sequence string representation
	 */
	@Override
	public String toString(){
//		update();
		assert(verify(true));
		return forwardSequence();
	}

	/** Decodes array1 (the forward strand) to a base-sequence String,
	 * WITHOUT the verify() assertion toString() carries -- needed so
	 * verify()'s own base-level comparison doesn't recurse into itself via
	 * toString(). */
	private String forwardSequence(){
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<mult; i++){
			bb.appendKmer(array1[i], perWordK(i));
//			bb.append(" ");
		}
////		bb.append("~");
//		for(int i=0; i<mult; i++){
//			bb.appendKmer(array2[i], k);
////			bb.append(" ");
//		}
		return bb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------      Key Selection and Hashing   ------------*/
	/*--------------------------------------------------------------*/

	/** Compares lookup hashes and complete selected words; requires a nonnull compatible key.
	 * Even hollow keys compare full retained words here, unlike hash-only table matching.
	 * This typed overload does not override Object.equals or define a hashCode contract. */
	public boolean equals(Kmer x){
		if(xor()!=x.xor()){return false;}
		return AbstractKmerTableU.equals(key(), x.key());
	}

	/** Compares hashes and forward words for nonnull compatible keys. */
	public boolean sameOrientation(Kmer x){
		if(xor()!=x.xor()){return false;}
		return Tools.equals(array1, x.array1);
	}

	/** Compares current key words using AbstractKmerTableU's ordering; requires compatible keys.
	 * This is distinct from strand selection's unsigned ordering; see that comparator's UKM106 note. */
	public int compareTo(Kmer x){
		return compare(key(), x.key());
	}

	/** Legacy raw-key overload guarded as unsupported by an unconditional assertion. */
	public int compareTo(long[] key2){
		assert(false);
		return compare(key(), key2);
	}

	/** Delegates equal-length raw-key comparison to AbstractKmerTableU. */
	public static int compare(long[] key1, long[] key2){
//		assert(false); //Why was this here?
		return AbstractKmerTableU.compare(key1, key2);
	}

	/** Legacy raw-key equality overload guarded as unsupported by an unconditional assertion. */
	public static boolean equals(long[] key1, long[] key2){
		assert(false);
		return AbstractKmerTableU.equals(key1, key2);
	}

	/** Returns borrowed forward words; use fillArray2 after direct mutation. */
	public long[] array1(){return array1;}

	/** Returns borrowed reverse-complement words; callers must not mutate them independently. */
	public long[] array2(){return array2;}

	/**
	 * Returns the canonical key array for this k-mer.
	 * Updates internal state if needed before returning the key.
	 * In solid canonical mode, chooses the unsigned lexicographically greater masked
	 * core, breaking core ties with complete words. Hollow mode compares the two arms;
	 * equal arms may retain either orientation. Forward-only mode selects array1.
	 * @return Borrowed internal key array; copy it before retaining across mutation
	 */
	public long[] key(){
		update();
//		assert(verify(false));
		return key;
	}

	/** Returns whether strand-selection cores/arms tie, after refreshing cached state. */
	public boolean corePalindrome(){//TODO: This can be set as a flag from setKey0
		update();
		return corePalindrome;
	}

	/** Selects a strand using instance masks or hollow arms; records a selection-core tie. */
	private void setKey0(){
		corePalindrome=false;
		key=array1;
		if(!rcomp){return;}
		if(hollow){
			// Interior bases cannot choose the strand: that would leak masked errors
			// back into both fingerprints when the selected arms are unchanged.
			int c=Long.compareUnsigned(array1[0], array2[0]);
			if(c==0){c=Long.compareUnsigned(rightArm(array1), rightArm(array2));}
			if(c<0){key=array2;}
			corePalindrome=(c==0);
			return;
		}
		for(int i=0; i<mult; i++){
			final long cm=coreMaskFor(i);
			final long a=array1[i]&cm, b=array2[i]&cm;
			//Unsigned compare (Amber, 2026-08-19): before packing, every word held
			//<=31 bases, so bit 63 was always 0 and signed compare behaved
			//identically to unsigned. A full 32-base word under PACKED can
			//set bit 63 (negative long) -- signed compare would then pick
			//the wrong strand as canonical. Required alongside packing, not
			//an optional follow-up; it preserves ordering for the legacy 31-base words.
			final int c=Long.compareUnsigned(a, b);
			if(c>0){return;}
			else if(c<0){
				key=array2;
				return;
			}
		}
		corePalindrome=true;
		setKey0safe();
	}

	/** Breaks a solid core tie using the complete unsigned words. */
	private void setKey0safe(){
		key=array1;
		for(int i=0; i<mult; i++){
			final long a=array1[i], b=array2[i];
			final int c=Long.compareUnsigned(a, b);
			if(c>0){break;}
			else if(c<0){
				key=array2;
				break;
			}
		}
	}

	/** Folds a nonempty word array with the same caller-supplied mask on every word.
	 * FULL_MIX controls first-word mixing. Returns a nonnegative placement hash;
	 * instance solid keys use the two-end-mask overload instead. */
	public static long xor(long[] key, long coreMask){
		long xor=(FULL_MIX ? Tools.hash64shift(key[0]&coreMask) : key[0]&coreMask);
		for(int i=1; i<key.length; i++){
			xor=(Long.rotateLeft(xor, 25))^Tools.hash64shift(key[i]&coreMask);
		}
		return xor&mask63;
	}

	/**
	 * Item 1b dual-mask overload: firstMask covers word[0]'s global-first-
	 * base, lastMask covers word[key.length-1]'s global-last-base; every
	 * word strictly between them is unmasked (raw). When key.length==1,
	 * word[0] IS the last word, so both masks apply to it together (ANDed) --
	 * matching the single-word toCoreMask formula exactly for that case.
	 */
	public static long xor(long[] key, long firstMask, long lastMask){
		final int maxi=key.length-1;
		final long w0mask=(maxi==0 ? (firstMask&lastMask) : firstMask);
		long xor=(FULL_MIX ? Tools.hash64shift(key[0]&w0mask) : key[0]&w0mask);
		for(int i=1; i<key.length; i++){
			final long m=(i==maxi ? lastMask : -1L);
			xor=(Long.rotateLeft(xor, 25))^Tools.hash64shift(key[i]&m);
		}
		return xor&mask63;
	}

	/**
	 * Second 63-bit fingerprint for supplied canonical words; an empty array yields zero. Unlike
	 * xor(), this deliberately retains both terminal bases: xor() is also the
	 * placement hash and may mask the two ends so extensions share a probe
	 * neighborhood, while xor2() is used to distinguish exact keys in compact
	 * hash-only tables.  Folding right with a different rotation and mixer
	 * supplies a separate fingerprint from xor(); neither fingerprint is collision-free.
	 */
	public static long xor2(long[] key){
		long xor=0;
		for(long word : key){
			xor=Tools.splitMix64(Long.rotateRight(xor, 23)^word);
		}
		return xor&mask63;
	}

	/**
	 * When true, the first word is put through the same hash64shift mixer as
	 * every other word before folding, instead of being used raw. Initially disabled
	 * to preserve prior behavior (including the mult=1 pass-through
	 * identity xor()==key[0] some test harnesses depend on). Item 1c fix
	 * (2026-08-19, Amber's roadmap): with FULL_MIX off, the final {@code &mask63}
	 * truncation drops bit 63 -- and because rotateLeft(xor,25) moves LOWER,
	 * UNMIXED bits of the raw first word into position 63, that truncation can
	 * silently discard real distinguishing information in bloom-filter mode
	 * (no stored key to fall back on). Mixing word 0 first means the bit lost
	 * at position 63 is one bit of a well-avalanched value, not a targeted raw
	 * bit -- the same reasoning already applied to words 1+ when the k>31
	 * BloomFilterCorrector2 weak-hash-fold theory was confirmed (small but
	 * real effect, see bloomfiltercorrector_ukmer_k31.plan).
	 * Brian changed the default to true on 2026-08-21. Configure this before using
	 * related solid keys/tables; changing it does not invalidate already cached hashes. */
	public static boolean FULL_MIX=true;

	/**
	 * Returns the cached XOR hash value for this k-mer.
	 * Updates internal state if needed before returning the hash.
	 * @return 63-bit hash value for efficient k-mer comparison and table indexing
	 */
	public long xor(){
		update();
		return lastXor;
	}

	/** Returns the lazily computed second fingerprint: full words for solid keys,
	 * outer arms and span for hollow keys. The retained interior is not hashed in hollow mode. */
	public long xor2(){
		update();
		if(lastXor2Incarnation!=incarnation){
			lastXor2=hollow ? Tools.splitMix64(Long.rotateRight(Tools.splitMix64(key[0]), 23)
					^rightArm(key)^Tools.splitMix64(kbig))&mask63 : xor2(key);
			lastXor2Incarnation=incarnation;
		}
		return lastXor2;
	}

	/**
	 * Returns this k-mer's hash value modulo the specified divisor.
	 * Useful for hash table indexing and partitioning operations.
	 * @param divisor Positive shard count or modulus
	 * @return Hash value modulo divisor
	 */
	public int mod(int divisor){
		int x=(int)(xor()%divisor);
//		System.err.println(xor()+"%"+value+"="+x);
		return x;
	}

	/** Exchanges the forward and reverse arrays without copying; valid-run length is unchanged. */
	public void rcomp(){
		long[] temp=array1;
		array1=array2;
		array2=temp;
		//Canonical fingerprints are orientation-invariant; forward-only keys must be refreshed.
		incarnation+=(rcomp ? 0 : 1);
	}

	/** Lazily refreshes strand selection and placement hash for the current incarnation. */
	private void update(){
		if(verbose){System.err.println("update() - len="+len);}
		assert(TESTMODE || len>=kbig) : len+", "+kbig+", "+array1[0];
		if(incarnation==lastIncarnation){return;}
		setKey0();
		lastXor=xor0();
		lastIncarnation=incarnation;
		if(verbose){System.err.println("After update - kmer "+this+"; key="+Arrays.toString(key)+"; a1="+Arrays.toString(array1())+"; a2="+Arrays.toString(array2()));}
	}

	/** Computes the placement hash from the selected key and this instance's hash geometry. */
	private long xor0(){
		if(hollow){
			return (Long.rotateLeft(Tools.hash64shift(key[0]), 25)
					^Tools.hash64shift(rightArm(key))^Tools.hash64shift(kbig))&mask63;
		}
		return xor(key, firstWordCoreMask, lastWordCoreMask);
	}

	/** Extracts exactly 32 suffix bases even when the packed last word is partial. */
	private long rightArm(final long[] strand){
		assert(hollow && k==FULL_WORD_K && mult>=3) :
			"Hollow suffix extraction needs nonoverlapping full-width arms.";
		if(lastWordK==FULL_WORD_K){return strand[maxindex];}
		return (strand[maxindex-1]<<(2*lastWordK))|strand[maxindex];
	}

	/** Displays raw arrays and the last selected key without refreshing cached state. */
	public String arraysToString(){
		return "key="+Arrays.toString(key)+", a1="+Arrays.toString(array1)+", a2="+Arrays.toString(array2);
	}

	/** Counts C/G encodings in forward storage; use only with a complete valid window. */
	public final int gc(){
		int gc=0;
		for(long kmer : array1){
			//Amber, 2026-08-19: was `kmer>0` -- a full 32-base word can set bit
			//63 (negative long), which would end the loop early and undercount.
			//The fix was introduced before packed storage became the default.
			//>>>= already shifts unsigned, so `!=0` terminates correctly
			//in both regimes.
			while(kmer!=0){
				long x=kmer&3;
				kmer>>>=2;
				if(x==1 || x==2){gc++;}
			}
		}
		return gc;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	/** Global strand-selection mode; configure before creating/using related tables and keys. */
	static boolean rcomp=true;

	/** Cached nonnegative placement hash. */
	private long lastXor=-1;
	/** Cached second fingerprint. */
	private long lastXor2=-1;
	/** Mutation generation, incremented by rolling/edit/rebuild operations. */
	private long incarnation=0;
	/** Generation represented by key, corePalindrome and lastXor. */
	private long lastIncarnation=-1;
	/** Generation represented by lastXor2. */
	private long lastXor2Incarnation=-1;
	/** Whether compared selection cores/arms tie; not necessarily full-sequence palindromy. */
	private boolean corePalindrome=false;
	/** Borrowed alias of one of this instance's strand arrays. */
	private long[] key=null;

	/** Owned forward words, exchanged with array2 by rcomp. */
	private long[] array1;
	/** Owned reverse-complement words with the same leading-full/partial-tail layout. */
	private long[] array2;
	/** Effective represented length in bases. */
	public final int kbig;
	/** Leading-word width in bases; the final word can be shorter. */
	public final int k;
	/** True only for an end-anchor key; explicit sequence arrays still retain the entire span. */
	public final boolean hollow;
	/** Word count and final valid array index (mult-1), respectively. */
	final int mult, maxindex;

	/** Bit width of a leading word. */
	private final int shift;
	/** Position of the leading base in a full word. */
	private final int shift2;
	/** Mask retaining the valid bits of a full word. */
	private final long mask;
	/** Masks for the two GLOBAL kmer ends (see the master constructor's
	 * comment) -- replaces the old uniform per-word coreMask, which was only
	 * correct when every word was an end (single-word/{@code k<=31}). */
	private final long firstWordCoreMask;
	private final long lastWordCoreMask;

	/** Width/shift/mask for word[maxindex] specifically -- see the master
	 * constructor. Equal to k/shift2/mask whenever kbig==k*mult (every
	 * pre-packing caller); genuinely different only under PACKED. */
	private final int lastWordK;
	private final int lastShift2;
	private final long lastMask;

	/** Length of the valid right-rolling run; may exceed the represented window length. */
	public int len=0; //TODO: Make private; use getter.
	/** Returns valid-run length, not storage capacity. */
	public final int len(){return len;}

	/** Whether new instances exclude global terminal bases from their selection/placement core. */
	public static boolean MASK_CORE=false;
	/** Item 1a (2026-08-19): when true, Kmer(int kbig_) uses the packed
	 * layout (leading words full at 32 bases, last word partial) instead
	 * of the old symmetric even-split. Default flipped to true 2026-08-21
	 * (Brian). Convenience constructors and static geometry helpers consult this;
	 * existing instances retain their captured widths. The explicit k,mult constructor
	 * bypasses it. Keep the setting consistent with code that derives table/buffer
	 * geometry through the static helpers. */
	public static boolean PACKED=true;
	/** Base count of a full leading word under the packed layout (64 bits
	 * / 2 bits-per-base). */
	private static final int FULL_WORD_K=32;
	/** Drops the sign bit from fingerprints used for nonnegative table routing. */
	private static final long mask63=Long.MAX_VALUE;
	/** Allows incomplete-window lookup only in explicitly rebuilt diagnostic code. */
	private static final boolean TESTMODE=false; //123
	/** Enables internal lookup diagnostics; disabled in ordinary builds. */
	private static final boolean verbose=false;
}
