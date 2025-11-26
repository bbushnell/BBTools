package ukmer;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date Jul 9, 2015
 *
 */
public class Kmer implements Cloneable {
	
	/** Copy constructor that creates a new Kmer identical to the provided one.
	 * @param o The Kmer to copy from */
	public Kmer(Kmer o){
		this(o.k, o.mult);
		setFrom(o);
	}
	
	/**
	 * Constructs a Kmer with the specified total k-mer size.
	 * Automatically calculates the optimal word size and multiplicity.
	 * @param kbig_ Total k-mer length in bases
	 */
	public Kmer(int kbig_){
		this(getK(kbig_), getMult(kbig_));
	}
	
	@Override
	public Kmer clone(){
		return new Kmer(this);
	}
	
	/**
	 * Constructs a Kmer with specified word size and multiplicity.
	 * Initializes bit manipulation constants and storage arrays.
	 * @param k_ K-mer size per word (maximum 31)
	 * @param mult_ Number of words needed to store the complete k-mer
	 */
	public Kmer(int k_, int mult_){
		k=k_;
		mult=mult_;
		maxindex=mult-1;
		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		coreMask=toCoreMask(k);
		
		kbig=k*mult;
		array1=new long[mult];
		array2=new long[mult];
		key=null;
	}
	
	/**
	 * Creates a core mask for k-mer comparison based on MASK_CORE setting.
	 * When core masking is enabled, masks the first and last bases for fuzzy matching.
	 * @param k K-mer length to create mask for
	 * @return Bit mask for core comparison, or -1L if core masking is disabled
	 */
	public static final long toCoreMask(int k){
//		System.err.println(k+", "+MASK_CORE);
		return MASK_CORE ? ((~((-1L)<<(2*k)))>>2)&~(3L) : -1L;
	}
	
	/**
	 * Calculates the number of words needed to store a k-mer of given size.
	 * Optimizes for memory efficiency while ensuring adequate storage.
	 * @param kbig Total k-mer length in bases
	 * @return Number of long words needed to store the k-mer
	 */
	public static int getMult(int kbig){
		final int mult=getMult0(kbig);
//		assert(mult==getMult0(kbig*(mult/mult))) : mult+", "+getMult0(mult*(kbig/mult));
		return mult;
	}
	
	/**
	 * Returns the actual k-mer size that can be stored with the given target size.
	 * May be smaller than requested due to word alignment requirements.
	 * @param kbig Requested k-mer size
	 * @return Actual k-mer size that will be used
	 */
	public static int getKbig(int kbig){
		int x=getMult(kbig)*getK(kbig);
		assert(x<=kbig) : x+", "+kbig;
		assert(kbig>31 || x==kbig);
		return x;
	}
	
	/**
	 * Internal calculation of word multiplicity for k-mer storage.
	 * Balances between minimizing words used and maximizing k-mer size stored.
	 * @param kbig Target k-mer size
	 * @return Optimal number of words for storage
	 */
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
	
	/**
	 * Calculates the k-mer size per word for a given total k-mer size.
	 * @param kbig Total k-mer length
	 * @return K-mer size per individual word (maximum 31)
	 */
	public static int getK(int kbig){
		int mult=getMult(kbig);
		int x=kbig/mult;
		assert(x*mult<=kbig) : x+", "+kbig;
		assert(x<=31) : kbig+", "+mult+", "+x;
		return x;
	}
	
	/**
	 * Copies the contents of another Kmer into this one.
	 * Updates incarnation number to track modifications.
	 * @param o The Kmer to copy from
	 * @return This Kmer instance for method chaining
	 */
	public Kmer setFrom(Kmer o){
		for(int i=0; i<mult; i++){
			array1[i]=o.array1[i];
			array2[i]=o.array2[i];
			len=o.len;
		}
		incarnation++;
		return this;
	}
	
	/**
	 * Sets this Kmer's forward array from the provided long array.
	 * Automatically generates the reverse complement array.
	 * @param array Forward k-mer representation as long array
	 * @return This Kmer instance for method chaining
	 */
	public Kmer setFrom(long[] array){
		for(int i=0; i<mult; i++){
			array1[i]=array[i];
		}
		fillArray2();
		incarnation++;
		return this;
	}
	
	/** Resets the Kmer to empty state, clearing all arrays and counters.
	 * Zeroes both forward and reverse complement arrays completely. */
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
	
	/** Quickly resets the Kmer by only clearing length and incarnation counters.
	 * Does not zero arrays for performance, assuming they will be overwritten. */
	public void clearFast() {
		len=0;
		lastIncarnation=-1;
		incarnation=0;
		//incarnation++;
	}
	
	/**
	 * Verifies that forward and reverse complement arrays are properly synchronized.
	 * Optionally updates internal state before verification.
	 * @param update Whether to update internal state before verification
	 * @return true if arrays are properly synchronized, false otherwise
	 */
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
	
	/**
	 * Internal verification that forward and reverse arrays are reverse complements.
	 * Checks each word pair for proper reverse complement relationship.
	 * @return true if verification passes, false otherwise
	 */
	private boolean verify(){
		if(len<kbig){return true;}
		for(int i=maxindex, j=0; j<mult; j++, i--){
			long kmer=array1[i];
			long rkmer=array2[j];
			if(kmer!=AminoAcid.reverseComplementBinaryFast(rkmer, k)){
//				assert(false) : Arrays.toString(array1);
				return false;
			}
		}
		assert(incarnation==lastIncarnation) : "incarnation="+incarnation+", last="+lastIncarnation+", len="+len+", kbig="+kbig;
		return true;
	}
	
	/**
	 * Adds a base to the right end of the k-mer in a rolling fashion.
	 * Shifts existing bases left and drops the leftmost base.
	 * @param b The base to add (A, C, G, T)
	 * @return The base that was dropped from the left end
	 */
	public byte addRight(final byte b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addRightNumeric(x)];
	}
	
	/**
	 * Adds a base to the right end of the k-mer in a rolling fashion.
	 * Character version of addRight for convenience.
	 * @param b The base character to add (A, C, G, T)
	 * @return The base that was dropped from the left end
	 */
	public byte addRight(final char b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addRightNumeric(x)];
	}
	
	/**
	 * Adds a base to the left end of the k-mer in a rolling fashion.
	 * Shifts existing bases right and drops the rightmost base.
	 * @param b The base to add (A, C, G, T)
	 * @return The base that was dropped from the right end
	 */
	public byte addLeft(final byte b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addLeftNumeric(x)];
	}
	
	/**
	 * Adds a numeric base (0-3) to the right end using rolling window.
	 * Updates both forward and reverse complement arrays simultaneously.
	 * Handles ambiguous bases by resetting length counter.
	 *
	 * @param x Numeric base value (0=A, 1=C, 2=G, 3=T, <0=ambiguous)
	 * @return Numeric value of the base dropped from the left end
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
		
		for(int i=maxindex, j=0; j<mult; j++, i--){
			
			long y=(array1[i]>>>shift2)&3L;
			long y2=array2[j]&3L;
			
			//Update kmers
			array1[i]=((array1[i]<<2)|x)&mask;
			array2[j]=((array2[j]>>>2)|(x2<<shift2))&mask;
			
			x=y;
			x2=y2;
		}
		incarnation++;
		return x;
	}
	
	/**
	 * Adds a numeric base to the left end using rolling window.
	 * Requires the k-mer to be full-length before calling.
	 * @param x Numeric base value (0=A, 1=C, 2=G, 3=T)
	 * @return Numeric value of the base dropped from the right end
	 */
	public long addLeftNumeric(long x){
		assert(x>=0 && x<4) : x;
		long x2=AminoAcid.numberToComplement[(int)x];
		
		assert(x>=0);
		assert(len>=kbig);
		
		for(int i=0, j=maxindex; i<mult; i++, j--){

			long y=array1[i]&3L;
			long y2=(array2[j]>>>shift2)&3L;
			
			//Update kmers
			array1[i]=(array1[i]>>>2)|(x<<shift2);
			array2[j]=((array2[j]<<2)|x2)&mask;
			
			x=y;
			x2=y2;
		}
		incarnation++;
		return x;
	}
	
	/** Generates the reverse complement array from the forward array.
	 * Sets length to full k-mer size and updates incarnation counter. */
	public void fillArray2() {
		for(int i=maxindex, j=0; j<mult; j++, i--){
			array2[j]=AminoAcid.reverseComplementBinaryFast(array1[i], k);
		}
		len=kbig;
		incarnation++;
	}
	
	@Override
	public String toString(){
//		update();
		assert(verify(true));
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<mult; i++){
			bb.appendKmer(array1[i], k);
//			bb.append(" ");
		}
////		bb.append("~");
//		for(int i=0; i<mult; i++){
//			bb.appendKmer(array2[i], k);
////			bb.append(" ");
//		}
		return bb.toString();
	}
	
	/**
	 * Tests equality with another Kmer using canonical representation.
	 * First compares hash values, then performs detailed comparison if needed.
	 * @param x The Kmer to compare with
	 * @return true if k-mers represent the same sequence (any orientation)
	 */
	public boolean equals(Kmer x){
		if(xor()!=x.xor()){return false;}
		return AbstractKmerTableU.equals(key(), x.key());
	}
	
	/**
	 * Tests if two k-mers have the same orientation (not just canonical equality).
	 * Compares forward arrays directly rather than canonical keys.
	 * @param x The Kmer to compare with
	 * @return true if k-mers have identical forward orientation
	 */
	public boolean sameOrientation(Kmer x){
		if(xor()!=x.xor()){return false;}
		return Tools.equals(array1, array2); //Possible bug: should be Tools.equals(array1, x.array1)
	}
	
	/**
	 * Compares this k-mer to another using lexicographic ordering of canonical keys.
	 * @param x The Kmer to compare with
	 * @return Negative, zero, or positive for less than, equal, or greater than
	 */
	public int compareTo(Kmer x){
		return compare(key(), x.key());
	}
	
	/**
	 * Compares this k-mer's canonical key to a provided key array.
	 * @param key2 The key array to compare with
	 * @return Negative, zero, or positive for less than, equal, or greater than
	 */
	public int compareTo(long[] key2){
		assert(false);
		return compare(key(), key2);
	}
	
	/**
	 * Static method to compare two k-mer key arrays lexicographically.
	 * @param key1 First key array
	 * @param key2 Second key array
	 * @return Negative, zero, or positive for less than, equal, or greater than
	 */
	public static int compare(long[] key1, long[] key2){
//		assert(false); //Why was this here?
		return AbstractKmerTableU.compare(key1, key2);
	}
	
	/**
	 * Static method to test equality of two k-mer key arrays.
	 * @param key1 First key array
	 * @param key2 Second key array
	 * @return true if the key arrays are identical
	 */
	public static boolean equals(long[] key1, long[] key2){
		assert(false);
		return AbstractKmerTableU.equals(key1, key2);
	}
	
	/** Returns the forward k-mer array */
	public long[] array1(){return array1;}
		
	/** Returns the reverse complement k-mer array */
	public long[] array2(){return array2;}
	
	/** WARNING!
	 * Do not confuse this with xor()! */
	public long[] key(){
		update();
//		assert(verify(false));
		return key;
	}
	
	/**
	 * Checks if the core (masked) regions of forward and reverse are identical.
	 * Indicates whether this k-mer is a palindrome when considering only core bases.
	 * @return true if core regions are palindromic
	 */
	public boolean corePalindrome(){//TODO: This can be set as a flag from setKey0
		update();
		return corePalindrome;
	}
	
	/** Determines the canonical key by comparing core-masked forward and reverse arrays.
	 * Sets corePalindrome flag and falls back to full comparison if cores are equal. */
	private void setKey0(){
		corePalindrome=false;
		key=array1;
		if(!rcomp){return;}
		for(int i=0; i<mult; i++){
			final long a=array1[i]&coreMask, b=array2[i]&coreMask;
			if(a>b){return;}
			else if(a<b){
				key=array2;
				return;
			}
		}
		corePalindrome=true;
		setKey0safe();
	}
	
	/** Fallback method for key selection when core regions are identical.
	 * Compares full arrays without masking to determine canonical orientation. */
	private void setKey0safe(){
		key=array1;
		for(int i=0; i<mult; i++){
			final long a=array1[i], b=array2[i];
			if(a>b){break;}
			else if(a<b){
				key=array2;
				break;
			}
		}
	}
	
	/**
	 * Computes XOR hash of a k-mer key array using rotating XOR operations.
	 * Applies core masking to focus hash on significant regions.
	 *
	 * @param key The k-mer key array to hash
	 * @param coreMask Mask to apply before hashing
	 * @return 63-bit hash value for the k-mer
	 */
	public static long xor(long[] key, long coreMask){
		long xor=key[0]&coreMask;
		for(int i=1; i<key.length; i++){
			xor=(Long.rotateLeft(xor, 25))^(key[i]&coreMask);
		}
		return xor&mask63;
	}
	
	/** WARNING!
	 * Do not confuse this with key()! */
	public long xor(){
		update();
		return lastXor;
	}

	/**
	 * @param divisor
	 * @return This kmer's xor modulo the divisor
	 */
	public int mod(int divisor) {
		int x=(int)(xor()%divisor);
//		System.err.println(xor()+"%"+value+"="+x);
		return x;
	}
	
	/** Swaps the forward and reverse complement arrays.
	 * Effectively converts the k-mer to its reverse complement orientation. */
	public void rcomp() {
		long[] temp=array1;
		array1=array2;
		array2=temp;
	}
	
	/**
	 * Updates internal cached values if the k-mer has been modified.
	 * Recalculates canonical key and hash value as needed.
	 * Uses incarnation counter to track when updates are necessary.
	 */
	private void update(){
		if(verbose){System.err.println("update() - len="+len);}
		assert(TESTMODE || len>=kbig) : len+", "+kbig+", "+array1[0];
		if(incarnation==lastIncarnation){return;}
		setKey0();
		lastXor=xor0();
		lastIncarnation=incarnation;
		if(verbose){System.err.println("After update - kmer "+this+"; key="+Arrays.toString(key)+"; a1="+Arrays.toString(array1())+"; a2="+Arrays.toString(array2()));}
	}
	
	/** Internal method to compute XOR hash using current key and core mask.
	 * @return XOR hash value for the current canonical key */
	private long xor0(){
		return xor(key, coreMask);
	}
	
	/**
	 * Returns debug string showing all internal arrays.
	 * Displays canonical key, forward array, and reverse complement array.
	 * @return Formatted string with array contents for debugging
	 */
	public String arraysToString() {
		return "key="+Arrays.toString(key)+", a1="+Arrays.toString(array1)+", a2="+Arrays.toString(array2);
	}
	
	/**
	 * Counts GC bases (C and G) in the forward k-mer representation.
	 * Iterates through all words and counts bases with numeric values 1 or 2.
	 * @return Number of G and C bases in the k-mer
	 */
	public final int gc(){
		int gc=0;
		for(long kmer : array1){
			while(kmer>0){
				long x=kmer&3;
				kmer>>>=2;
				if(x==1 || x==2){gc++;}
			}
		}
		return gc;
	}
	
	/** Whether to consider reverse complement for canonical representation */
	static boolean rcomp=true;
	
	/** Cached XOR hash value to avoid recomputation */
	private long lastXor=-1;
	/** Counter tracking modifications to detect when updates are needed */
	private long incarnation=0;
	/** Last incarnation when internal state was updated */
	private long lastIncarnation=-1;
	/** Whether the core-masked regions form a palindrome */
	private boolean corePalindrome=false;
	/** Canonical key array (reference to either array1 or array2) */
	private long[] key=null;
	
	/** Forward k-mer representation as array of longs */
	private long[] array1;
	/** Reverse complement k-mer representation as array of longs */
	private long[] array2;
	/** Total k-mer length in bases (k * mult) */
	public final int kbig;
	/** K-mer length per word (maximum 31 bases) */
	public final int k;
	/** Maximum valid index in the storage arrays (mult - 1) */
	/** Number of long words needed to store the complete k-mer */
	final int mult, maxindex;
	
	/** Bit shift amount for full k-mer (2 * k) */
	private final int shift;
	/** Bit shift amount for k-mer minus one base (shift - 2) */
	private final int shift2;
	/** Bit mask to retain only k-mer bits after shifting operations */
	private final long mask;
	/** Mask for core region comparison, excluding terminal bases if enabled */
	private final long coreMask;
	
	/** Current number of bases accumulated in the k-mer */
	public int len=0; //TODO: Make private; use getter.
	/** Returns the current length of accumulated bases in the k-mer */
	public final int len(){return len;}
	
	/** Whether to enable core masking for fuzzy k-mer comparisons */
	public static boolean MASK_CORE=false;
	/** Mask to ensure hash values fit in 63 bits (Long.MAX_VALUE) */
	private static final long mask63=Long.MAX_VALUE;
	/** Debug flag for testing mode with relaxed assertions */
	private static final boolean TESTMODE=false; //123
	/** Debug flag for verbose logging of operations */
	private static final boolean verbose=false;
}
