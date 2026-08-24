package ukmer;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Represents a k-mer (DNA sequence of length k) with support for large k-values.
 * Uses multi-word arrays to handle k-mers longer than 31 bases efficiently.
 * Maintains both forward and reverse complement representations with canonical ordering.
 * Supports rolling operations for sequence sliding windows and hash-based comparisons.
 * @author Brian Bushnell
 * @date Jul 9, 2015
 */
public class Kmer implements Cloneable {
	
	public Kmer(Kmer o){
		this(o.k, o.mult, o.kbig);
		setFrom(o);
	}

	/**
	 * Convenience constructor. When PACKED is false, exactly today's behavior
	 * (symmetric even-split; kbig may round DOWN from kbig_ if it doesn't
	 * divide evenly -- see getKbig). When PACKED is true, uses the packed
	 * layout (leading words full at FULL_WORD_K=32, last word partial) which
	 * represents every kbig_&gt;0 exactly, no rounding.
	 */
	public Kmer(int kbig_){
		this(PACKED ? FULL_WORD_K : getK(kbig_),
			PACKED ? getMultPacked(kbig_) : getMult(kbig_),
			PACKED ? kbig_ : getK(kbig_)*getMult(kbig_));
	}

	/** Creates and returns a copy of this Kmer.
	 * @return A new Kmer instance identical to this one */
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

	private static int getMultPacked(int kbig){
		return (kbig+FULL_WORD_K-1)/FULL_WORD_K;
	}
	
	public static final long toCoreMask(int k){
//		System.err.println(k+", "+MASK_CORE);
		return MASK_CORE ? ((~((-1L)<<(2*k)))>>2)&~(3L) : -1L;
	}
	
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

	public Kmer setFrom(Kmer o){
		for(int i=0; i<mult; i++){
			array1[i]=o.array1[i];
			array2[i]=o.array2[i];
			len=o.len;
		}
		incarnation++;
		return this;
	}
	
	public Kmer setFrom(long[] array){
		for(int i=0; i<mult; i++){
			array1[i]=array[i];
		}
		fillArray2();
		incarnation++;
		return this;
	}
	
	public void clear() {
		len=0;
		for(int i=0; i<mult; i++){
			array1[i]=0;
			array2[i]=0;
		}
		lastIncarnation=-1;
		incarnation=0;
		//incarnation++;
	}
	
	public void clearFast() {
		len=0;
		lastIncarnation=-1;
		incarnation=0;
		//incarnation++;
	}
	
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
	
	public byte addRight(final byte b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addRightNumeric(x)];
	}
	
	public byte addRight(final char b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addRightNumeric(x)];
	}
	
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
	 * -- R rolls right-entry, the mirror image of addRightNumeric's case. */
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
	 */
	public void fillArray2() {
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
	
	public boolean equals(Kmer x){
		if(xor()!=x.xor()){return false;}
		return AbstractKmerTableU.equals(key(), x.key());
	}
	
	public boolean sameOrientation(Kmer x){
		if(xor()!=x.xor()){return false;}
		return Tools.equals(array1, x.array1);
	}
	
	public int compareTo(Kmer x){
		return compare(key(), x.key());
	}
	
	public int compareTo(long[] key2){
		assert(false);
		return compare(key(), key2);
	}
	
	public static int compare(long[] key1, long[] key2){
//		assert(false); //Why was this here?
		return AbstractKmerTableU.compare(key1, key2);
	}
	
	public static boolean equals(long[] key1, long[] key2){
		assert(false);
		return AbstractKmerTableU.equals(key1, key2);
	}
	
	public long[] array1(){return array1;}
		
	public long[] array2(){return array2;}
	
	/**
	 * Returns the canonical key array for this k-mer.
	 * Updates internal state if needed before returning the key.
	 * The key is either the forward or reverse array, whichever is lexicographically smaller.
	 * @return Canonical key array (either array1 or array2)
	 */
	public long[] key(){
		update();
//		assert(verify(false));
		return key;
	}
	
	public boolean corePalindrome(){//TODO: This can be set as a flag from setKey0
		update();
		return corePalindrome;
	}
	
	private void setKey0(){
		corePalindrome=false;
		key=array1;
		if(!rcomp){return;}
		for(int i=0; i<mult; i++){
			final long cm=coreMaskFor(i);
			final long a=array1[i]&cm, b=array2[i]&cm;
			//Unsigned compare (Amber, 2026-08-19): today every word holds
			//<=31 bases so bit 63 is always 0 and signed compare behaves
			//identically to unsigned. A full 32-base word under PACKED can
			//set bit 63 (negative long) -- signed compare would then pick
			//the wrong strand as canonical. Required alongside packing, not
			//an optional follow-up, since it's a no-op change today.
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
	
	/** Uniform-mask overload -- unchanged, still used directly by
	 * HashForestU (single-mask context, not the multi-word packed layout).
	 * Kmer's own instance path uses the dual-mask overload below. */
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
	 * When true, the first word is put through the same hash64shift mixer as
	 * every other word before folding, instead of being used raw. Default OFF
	 * to preserve exact prior behavior (including the mult=1 pass-through
	 * identity xor()==key[0] some test harnesses depend on). Item 1c fix
	 * (2026-08-19, Amber's roadmap): with FULL_MIX off, the final `&mask63`
	 * truncation drops bit 63 -- and because rotateLeft(xor,25) moves LOWER,
	 * UNMIXED bits of the raw first word into position 63, that truncation can
	 * silently discard real distinguishing information in bloom-filter mode
	 * (no stored key to fall back on). Mixing word 0 first means the bit lost
	 * at position 63 is one bit of a well-avalanched value, not a targeted raw
	 * bit -- the same reasoning already applied to words 1+ when the k>31
	 * BloomFilterCorrector2 weak-hash-fold theory was confirmed (small but
	 * real effect, see bloomfiltercorrector_ukmer_k31.plan).
	 * Default flipped to true 2026-08-21 (Brian) -- the more accurate mode,
	 * pending a k=64 regression check confirming results barely change. */
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

	/**
	 * Returns this k-mer's hash value modulo the specified divisor.
	 * Useful for hash table indexing and partitioning operations.
	 * @param divisor The divisor for modulo operation
	 * @return Hash value modulo divisor
	 */
	public int mod(int divisor) {
		int x=(int)(xor()%divisor);
//		System.err.println(xor()+"%"+value+"="+x);
		return x;
	}
	
	public void rcomp() {
		long[] temp=array1;
		array1=array2;
		array2=temp;
	}
	
	private void update(){
		if(verbose){System.err.println("update() - len="+len);}
		assert(TESTMODE || len>=kbig) : len+", "+kbig+", "+array1[0];
		if(incarnation==lastIncarnation){return;}
		setKey0();
		lastXor=xor0();
		lastIncarnation=incarnation;
		if(verbose){System.err.println("After update - kmer "+this+"; key="+Arrays.toString(key)+"; a1="+Arrays.toString(array1())+"; a2="+Arrays.toString(array2()));}
	}
	
	private long xor0(){
		return xor(key, firstWordCoreMask, lastWordCoreMask);
	}
	
	public String arraysToString() {
		return "key="+Arrays.toString(key)+", a1="+Arrays.toString(array1)+", a2="+Arrays.toString(array2);
	}
	
	public final int gc(){
		int gc=0;
		for(long kmer : array1){
			//Amber, 2026-08-19: was `kmer>0` -- a full 32-base word can set bit
			//63 (negative long), which would end the loop early and undercount.
			//Inert today (no word holds 32 bases yet) but real once 1a packing
			//lands. >>>= already shifts unsigned, so `!=0` terminates correctly
			//in both regimes.
			while(kmer!=0){
				long x=kmer&3;
				kmer>>>=2;
				if(x==1 || x==2){gc++;}
			}
		}
		return gc;
	}
	
	static boolean rcomp=true;
	
	private long lastXor=-1;
	private long incarnation=0;
	private long lastIncarnation=-1;
	private boolean corePalindrome=false;
	private long[] key=null;
	
	private long[] array1;
	private long[] array2;
	public final int kbig;
	public final int k;
	/** Maximum valid index in the storage arrays (mult - 1) */
	final int mult, maxindex;
	
	private final int shift;
	private final int shift2;
	private final long mask;
	/** Masks for the two GLOBAL kmer ends (see the master constructor's
	 * comment) -- replaces the old uniform per-word coreMask, which was only
	 * correct when every word was an end (single-word/k<=31). */
	private final long firstWordCoreMask;
	private final long lastWordCoreMask;

	/** Width/shift/mask for word[maxindex] specifically -- see the master
	 * constructor. Equal to k/shift2/mask whenever kbig==k*mult (every
	 * pre-packing caller); genuinely different only under PACKED. */
	private final int lastWordK;
	private final int lastShift2;
	private final long lastMask;

	public int len=0; //TODO: Make private; use getter.
	public final int len(){return len;}

	public static boolean MASK_CORE=false;
	/** Item 1a (2026-08-19): when true, Kmer(int kbig_) uses the packed
	 * layout (leading words full at 32 bases, last word partial) instead
	 * of the old symmetric even-split. Default flipped to true 2026-08-21
	 * (Brian), pending a k=64 regression check confirming results barely
	 * change. Only Kmer(int kbig_) consults this directly -- every other
	 * method derives its per-word geometry from lastWordK/lastShift2/
	 * lastMask, which are set correctly for either regime at construction
	 * time. */
	public static boolean PACKED=true;
	/** Base count of a full leading word under the packed layout (64 bits
	 * / 2 bits-per-base). */
	private static final int FULL_WORD_K=32;
	private static final long mask63=Long.MAX_VALUE;
	private static final boolean TESTMODE=false; //123
	private static final boolean verbose=false;
}
