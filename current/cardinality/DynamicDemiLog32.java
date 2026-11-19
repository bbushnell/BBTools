package cardinality;

import java.util.Arrays;

import map.LongHashSet;
import parse.Parser;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Mutable DemiLog with unsigned 32-bit registers and an instance-bound layout.
 * With exponent6 each register has 26 inverted mantissa bits. Zero alone means
 * empty; a populated register may have a zero high16 prefix. Hash bucket bits
 * and seed convention match DynamicDemiLog. Overflow ranks saturate rather than
 * wrapping, and unavailable low hash bits are zero-padded before inversion.
 * No 16-bit calibration table or process-global exponent is reused.
 * Use one instance per worker and merge after workers finish. Borrowed register,
 * count, GC and witness arrays must be treated as read-only; direct writes bypass
 * floor maintenance and cardinality-cache invalidation. A shared blacklist must
 * remain unchanged while any sketch uses it.
 * @author Sayu
 */
public final class DynamicDemiLog32 extends AbstractDemiLog{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates the default 2048-bucket, k31, exponent6 sketch with counts and no witnesses. */
	DynamicDemiLog32(){
		this(2048, 31, -1, 0, true, false, 6);
	}

	/** Uses the parser's hashing settings with counts enabled and exponent6. */
	DynamicDemiLog32(Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob, true, false, 6);
	}

	/** Creates an exponent6 sketch with counts and no literal-kmer witnesses. */
	DynamicDemiLog32(int buckets, int k, long seed, float minProb){
		this(buckets, k, seed, minProb, true, false, 6);
	}

	/** Allocates owned arrays; counts also enables the matching GC payload. */
	private DynamicDemiLog32(int buckets_, int k_, long seed_, float minProb_, boolean counts, boolean kmers, int exponent_){
		super(buckets_, k_, seed_, minProb_);
		if(exponent_<1 || exponent_>8){throw new IllegalArgumentException("DDL32 exponent must be 1..8");}
		exponent=exponent_;
		mantissaBits=32-exponent;
		mantissaMask=(1L<<mantissaBits)-1;
		maxRank=Math.min(64, (1<<exponent)-1);
		maxArray=new int[buckets];
		countArray=counts ? new char[buckets] : null;
		gcArray=counts ? new byte[buckets] : null;
		kmerArray=kmers ? new long[buckets] : null;
		floorCount=buckets;
		assert(mantissaBits>0 && mantissaBits<32) : "32-bit score packing reserves exponent bits and an unsigned mantissa";
	}

	/** Creates an empty exponent6 sketch without witnesses, optionally retaining counts and GC. */
	public static DynamicDemiLog32 create(int buckets, int k, long seed, float minProb, boolean counts){
		return create(buckets, k, seed, minProb, counts, false, 6);
	}
	/**
	 * Creates an empty sketch with independently selected count and witness payloads.
	 * @param buckets Requested bucket count; the superclass rounds/clamps it to a power of two
	 * @param k Requested kmer length, normalized by CardinalityTracker
	 * @param seed XOR value used before hashing; -1 selects the process default seed
	 * @param minProb Minimum kmer correctness probability used by sequence hashing
	 * @param counts Whether to allocate both count and GC arrays
	 * @param kmers Whether to retain the winning literal kmer in each bucket
	 * @param exponent Exponent width in bits, from 1 through 8
	 * @return New empty sketch with owned storage
	 */
	public static DynamicDemiLog32 create(int buckets, int k, long seed, float minProb, boolean counts, boolean kmers, int exponent){
		return new DynamicDemiLog32(buckets, k, seed, minProb, counts, kmers, exponent);
	}

	/**
	 * Copies absolute wire values and optional witnesses into owned storage.
	 * The caller must supply the matching seed/exponent; no global layout is consulted.
	 * Counts/GC, added count, micro-index and blacklist are not restored; minProb is zero.
	 * @param values Nonempty power-of-two array of absolute unsigned register bit patterns
	 * @param k Kmer length for the reconstructed sketch
	 * @param seed Original effective hash seed; -1 still resolves to the process default
	 * @param exponent Original exponent width, from 1 through 8
	 * @param kmers Optional literal witnesses with the same length as values
	 * @return New sketch with recounted occupancy and early-exit floor
	 */
	public static DynamicDemiLog32 fromArray(int[] values, int k, long seed, int exponent, long[] kmers){
		if(values==null || values.length<1 || Integer.bitCount(values.length)!=1 || (kmers!=null && kmers.length!=values.length)){
			throw new IllegalArgumentException("DDL32 requires power-of-two registers and a matching optional kmer payload");
		}
		final DynamicDemiLog32 out=create(values.length, k, seed, 0, false, kmers!=null, exponent);
		for(int i=0; i<values.length; i++){
			if((values[i]>>>out.mantissaBits)>out.maxRank){throw new IllegalArgumentException("DDL32 rank exceeds the 64-bit hash at bucket "+i);}
			out.maxArray[i]=values[i];
		}
		if(kmers!=null){System.arraycopy(kmers, 0, out.kmerArray, 0, kmers.length);}
		out.recount();
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Worker Lifecycle       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns an empty worker copy with the same recipe and payload flags, not a data clone. */
	@Override
	public DynamicDemiLog32 copy(){
		return emptyCopy(countArray!=null, kmerArray!=null);
	}
	/**
	 * Creates an empty worker with owned arrays and the same hashing/filter recipe.
	 * The blacklist object is shared read-only; existing observations are not copied.
	 */
	public DynamicDemiLog32 emptyCopy(boolean counts, boolean kmers){
		final DynamicDemiLog32 out=create(buckets, k, hashXor, minProb, counts, kmers, exponent);
		out.blacklist=blacklist;
		return out;
	}

	/**
	 * Packs an already hashed 64-bit value into a nonzero unsigned register bit pattern.
	 * Excess rank saturates; unavailable low hash bits are zero-padded before inversion.
	 * This pure operation does not select a bucket or change a sketch.
	 */
	public static int scoreHash(long hash, int exponent){
		if(exponent<1 || exponent>8){throw new IllegalArgumentException("DDL32 exponent must be 1..8");}
		final int bits=32-exponent, limit=(1<<exponent)-1;
		final int zeros=Long.numberOfLeadingZeros(hash), rank=Math.min(zeros, limit);
		final long mask=(1L<<bits)-1;
		final long inv= zeros>limit || zeros>=63 ? mask : (~(hash<<(zeros+1)))>>>(64-bits);
		final long score=((long)rank<<bits)|inv;
		assert((score&~0xffffffffL)==0) : "The exponent plus mantissa must fit exactly one uint32 register";
		return (int)Math.max(1L, score);
	}

	/**
	 * Updates the winning register and optional payloads for a literal kmer.
	 * Equal scores increment the saturating count without replacing the witness.
	 * Use inherited add(long) externally to update the separate added counter as well.
	 */
	@Override
	public void hashAndStore(long number){
		final long hash=Tools.hash64shift(number^hashXor);
		if(Long.compareUnsigned(hash, earlyExitMask)>0){return;}
		final int bucket=(int)(hash&bucketMask), score=scoreHash(hash, exponent), old=maxArray[bucket];
		final int comparison=Integer.compareUnsigned(score, old);
		if(comparison<0 || (blacklist!=null && blacklist.contains(number))){return;}
		lastCardinality=-1;
		if(comparison>0){
			maxArray[bucket]=score;
			if(old==0){filled++;}
			if(kmerArray!=null){kmerArray[bucket]=number;}
			if(gcArray!=null){gcArray[bucket]=(byte)Long.bitCount(number&0x5555555555555555L);}
			if(countArray!=null){countArray[bucket]=1;}
			if((old>>>mantissaBits)==floor && (score>>>mantissaBits)>floor && --floorCount==0){advanceFloor();}
		}else if(countArray!=null){countArray[bucket]=(char)Math.min(65535, countArray[bucket]+1);}
		assert(filled<=buckets && floorCount>0) : "Early exit requires a nonempty minimum-rank tier after every successful update";
	}

	/** Borrows words before hashing; the caller must keep it unchanged while sketches use it. */
	public void setBlacklist(LongHashSet words){
		if(filled!=0 || added!=0){throw new IllegalStateException("Set the DDL32 blacklist before hashing");}
		blacklist=words;
	}

	/** Finds the next occupied rank after the last register at the old floor advances. */
	private void advanceFloor(){
		assert(floorCount==0 && floor<maxRank) : "Only exhausting the current minimum rank may advance the early-exit floor";
		do{
			floor++;
			for(int value : maxArray){if((value>>>mantissaBits)==floor){floorCount++;}}
		}while(floorCount==0 && floor<maxRank);
		earlyExitMask=floor==64 ? 0 : -1L>>>floor;
		assert(floorCount>0) : "A finite register array must have a minimum occupied rank";
	}

	/** Rebuilds occupancy/floor state after loading or merging and invalidates the cached estimate. */
	private void recount(){
		filled=0;
		floor=maxRank;
		floorCount=0;
		for(int value : maxArray){
			if(value!=0){filled++;}
			final int rank=value>>>mantissaBits;
			if(rank<floor){floor=rank; floorCount=1;}else if(rank==floor){floorCount++;}
		}
		earlyExitMask=floor==64 ? 0 : -1L>>>floor;
		lastCardinality=-1;
		assert(floorCount>0 && filled<=buckets) : "Deserialized and merged sketches must rebuild the same early-exit invariant as live sketches";
	}

	/**
	 * Merges corresponding unsigned maxima after workers have stopped mutating either sketch.
	 * Layout, payload flags and blacklist object identity must agree; self-merge is a no-op.
	 * Equal scores sum saturating counts and retain this sketch's witness/GC value.
	 */
	@Override
	public void add(CardinalityTracker source){
		if(!(source instanceof DynamicDemiLog32)){throw new IllegalArgumentException("Cannot merge a different tracker into DDL32");}
		final DynamicDemiLog32 other=(DynamicDemiLog32)source;
		requireCompatible(other);
		if(blacklist!=other.blacklist || (countArray!=null)!=(other.countArray!=null) || (kmerArray!=null)!=(other.kmerArray!=null)){
			throw new IllegalArgumentException("DDL32 merge requires matching counts, witness and blacklist configuration");
		}
		if(other==this){return;}
		added=Math.addExact(added, other.added);
		microIndex|=other.microIndex;
		for(int i=0; i<buckets; i++){
			final int cmp=Integer.compareUnsigned(maxArray[i], other.maxArray[i]);
			if(cmp<0){
				maxArray[i]=other.maxArray[i];
				if(countArray!=null){countArray[i]=other.countArray[i];gcArray[i]=other.gcArray[i];}
				if(kmerArray!=null){kmerArray[i]=other.kmerArray[i];}
			}else if(cmp==0 && countArray!=null){countArray[i]=(char)Math.min(65535, countArray[i]+other.countArray[i]);}
		}
		recount();
	}

	/** Folds with this sketch's blacklist preference when requested; see the filter overload. */
	@Override
	public DynamicDemiLog32 condenseSketch(int newBuckets, boolean useBlacklist){
		return condenseSketch(newBuckets, useBlacklist ? blacklist : null);
	}

	/**
	 * Folds source buckets modulo a smaller power-of-two width without changing this sketch.
	 * With witnesses and a filter, a nonblacklisted witness takes precedence; if all
	 * candidates are blacklisted, their unsigned maximum remains. Without witnesses,
	 * filter has no effect. This is a preference, not removal of every blacklisted word.
	 * @return This instance when width is unchanged, with no filtering; otherwise a new
	 * sketch retaining added/micro-index/blacklist and witnesses, but no counts or GC,
	 * and minProb reset to zero by fromArray
	 */
	public DynamicDemiLog32 condenseSketch(int newBuckets, LongHashSet filter){
		if(newBuckets<1 || Integer.bitCount(newBuckets)!=1 || newBuckets>buckets){throw new IllegalArgumentException("DDL32 can fold only to a smaller power-of-two width");}
		if(newBuckets==buckets){return this;}
		final int[] values=new int[newBuckets];
		final long[] words=kmerArray==null ? null : new long[newBuckets];
		final boolean[] bad=filter!=null && kmerArray!=null ? new boolean[newBuckets] : null;
		for(int i=0; i<buckets; i++){
			final int value=maxArray[i], target=i&(newBuckets-1);
			if(value==0){continue;}
			final boolean rejected=bad!=null && filter.contains(kmerArray[i]);
			if(values[target]==0 || (bad!=null && bad[target] && !rejected) ||
				((bad==null || bad[target]==rejected) && Integer.compareUnsigned(value, values[target])>0)){
				values[target]=value;
				if(words!=null){words[target]=kmerArray[i];}
				if(bad!=null){bad[target]=rejected;}
			}
		}
		final DynamicDemiLog32 out=fromArray(values, k, hashXor, exponent, words);
		out.added=added;
		out.microIndex=microIndex;
		out.blacklist=blacklist;
		assert(out.buckets==newBuckets) : "Folding must preserve each source bucket modulo the requested width";
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/

	/** Decodes registers into statistics, allocating packed values and reusing rank-count scratch. */
	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}else{Arrays.fill(nlzCounts, 0);}
		final long[] packed=new long[buckets];
		for(int i=0; i<buckets; i++){
			final int value=maxArray[i];
			if(value==0){nlzCounts[0]++;continue;}
			final int rank=value>>>mantissaBits;
			nlzCounts[rank+1]++;
			packed[i]=((long)(rank+1)<<mantissaBits)|(Integer.toUnsignedLong(value)&mantissaMask);
		}
		assert(nlzCounts[0]+filled==buckets) : "Statistics must preserve populated registers even when their upper16 bits are zero";
		return CardStats.forDemiLog32(packed, nlzCounts, exponent, mantissaBits, buckets, microIndex, added);
	}

	/** Returns the cached or newly computed estimate, optionally bounded by a positive added count. */
	@Override
	public long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats stats=summarize();
		long count=Math.max((long)stats.hybridDDL(), stats.microCardinality());
		if(clampToAdded && added>0){count=Math.min(count, added);}
		lastCardinalityStatic=lastCardinality=count;
		assert(count>=0) : "Cardinality estimates cannot be negative after unsigned-register decoding";
		return count;
	}
	/** Returns newly computed estimates in the legacy result-array layout. */
	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		return AbstractCardStats.buildLegacyArray(s, s.hybridDDL());
	}

	/** Returns null because this implementation does not use the legacy correction table. */
	@Override
	public float[] compensationFactorLogBucketsArray(){return null;}

	/*--------------------------------------------------------------*/
	/*----------------       Layout and Arrays      ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the full stored width, including exponent and mantissa. */
	@Override
	public int registerBits(){return 32;}
	/** Returns the register storage width used by cardinality clients. */
	@Override
	public int bitsPerWord(){return 32;}
	/** Returns this instance's exponent width. */
	@Override
	public int exponentWidth(){return exponent;}
	/** Returns the minimum rank across all registers, including empty buckets. */
	@Override
	public int getGlobalNLZ(){return floor;}
	/** Returns the number of nonzero registers. */
	@Override
	public int filledBuckets(){return filled;}
	/** Returns the occupied fraction of the bucket array. */
	@Override
	public double occupancy(){return (double)filled/buckets;}
	/** Returns borrowed register storage; direct edits invalidate cached state. */
	@Override
	public int[] maxArray32(){return maxArray;}
	/** Returns a detached copy of all unsigned register bit patterns. */
	@Override
	public int[] bucketValues(){return maxArray.clone();}
	/** Returns one unsigned register bit pattern; zero denotes an empty bucket. */
	@Override
	public int registerValue(int bucket){
		assert(bucket>=0 && bucket<buckets) : "Register lookup must stay within this sketch";
		return maxArray[bucket];
	}
	/** Returns borrowed saturating count storage, or null when counts are disabled. */
	@Override
	public char[] counts16(){return countArray;}
	/** Returns borrowed winner-GC storage, or null when counts are disabled. */
	@Override
	public byte[] gcArray(){return gcArray;}
	/** Returns borrowed literal-kmer witnesses, or null when not retained. */
	@Override
	public long[] kmerArray(){return kmerArray;}
	/** Returns whether a literal witness array was allocated. */
	@Override
	public boolean hasKmers(){return kmerArray!=null;}

	/*--------------------------------------------------------------*/
	/*----------------         Serialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Appends unsigned registers as a tab-separated A48 row without a newline; returns bb. */
	@Override
	public ByteBuilder toBytes(ByteBuilder bb){
		assert(bb!=null) : "A48 output requires an owned reusable byte buffer";
		for(int i=0; i<buckets; i++){
			if(i>0){bb.tab();}
			bb.appendA48(Integer.toUnsignedLong(maxArray[i]));
		}
		return bb;
	}
	/** Appends witnesses in bucket order, or nothing when absent; returns bb without a newline. */
	@Override
	public ByteBuilder toBytesKmers(ByteBuilder bb){
		assert(bb!=null) : "Optional witness output follows the same register order as the A48 row";
		if(kmerArray!=null){
			for(int i=0; i<buckets; i++){
				if(i>0){bb.tab();}
				bb.appendA48(kmerArray[i]);
			}
		}
		return bb;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Instance layout and largest rank representable from a 64-bit hash. */
	private final int exponent, mantissaBits, maxRank;
	/** Low mantissa bits used while decoding unsigned registers. */
	private final long mantissaMask;
	/** Owned unsigned register bit patterns; zero alone denotes an empty bucket. */
	private final int[] maxArray;
	/** Optional equal-score counts, saturated at 65535. */
	private final char[] countArray;
	/** Optional GC counts associated with the winning literal kmer. */
	private final byte[] gcArray;
	/** Optional winning literal kmer for each bucket, meaningful only for nonzero registers. */
	private final long[] kmerArray;
	/** Minimum rank, its population, and the total number of nonzero registers. */
	private int floor, floorCount, filled;
	/** Largest unsigned hash admitted by the current rank floor, including possible ties. */
	private long earlyExitMask=-1L;
	/** Optional borrowed literal-word filter, shared unchanged by worker copies. */
	private LongHashSet blacklist;
}
