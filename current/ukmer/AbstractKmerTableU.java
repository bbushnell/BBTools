package ukmer;

import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import kmer.KmerTableInterface;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Shared long-kmer table contract and serialization helpers.
 * Implementations store packed {@link Kmer} words or long-key hash pairs,
 * optionally retain multiple positive values per key, expose optional owner
 * arrays for parallel claiming, and use -1 as the missing value / list
 * terminator sentinel.
 *
 * @author Brian Bushnell
 * @date Oct 23, 2013
 */
public abstract class AbstractKmerTableU implements KmerTableInterface{

	/*--------------------------------------------------------------*/
	/*----------------         Kmer methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Increments this k-mer's scalar count, inserting it if absent.
	 * @param kmer canonical key to increment
	 * @return updated count reported by the concrete table
	 */
	public abstract int increment(Kmer kmer);

	/**
	 * Increments this k-mer and reports whether a new key was created.
	 * @param kmer canonical key to increment
	 * @return 1 when the key was inserted, 0 when it was already present
	 */
	public abstract int incrementAndReturnNumCreated(final Kmer kmer);

	/**
	 * Sets the value associated with a k-mer.
	 * @param kmer The k-mer key
	 * @param value The value to associate with the k-mer
	 * @return 1 when a new key was created, 0 when an existing key was updated
	 */
	public abstract int set(Kmer kmer, int value);

	/**
	 * Sets multiple values associated with a k-mer.
	 * @param kmer The k-mer key
	 * @param vals Array of values to associate with the k-mer
	 * @return 1 when a new key was created, 0 when an existing key was updated
	 */
	public abstract int set(Kmer kmer, int[] vals);

	/**
	 * Stores a scalar value only when the k-mer is currently absent.
	 * @param kmer The k-mer key
	 * @param value The value to associate with a newly inserted k-mer
	 * @return number of k-mers added
	 */
	public abstract int setIfNotPresent(Kmer kmer, int value);

	/**
	 * Fetches the primary scalar value associated with a k-mer.
	 * @param kmer The k-mer key
	 * @return value, or -1 when the k-mer is not present
	 */
	public abstract int getValue(Kmer kmer);

	/**
	 * Fetches all values associated with a k-mer.
	 * @param kmer The k-mer key
	 * @param singleton single-element scratch space for scalar tables
	 * @return null, singleton with slot 0 set to -1, or a value array with -1
	 * marking unused tail slots; returned arrays may be borrowed live storage
	 */
	public abstract int[] getValues(Kmer kmer, int[] singleton);

	/**
	 * Tests whether a k-mer is present in the table.
	 * @param kmer The k-mer to test
	 * @return true if the k-mer is present, false otherwise
	 */
	public abstract boolean contains(Kmer kmer);

//	public abstract boolean contains(Kmer kmer, int v);
//
//	public abstract boolean contains(Kmer kmer, int[] vals);
//
//	public abstract Object get(Kmer kmer);

	/**
	 * Compares two k-mer keys represented as long arrays.
	 * Intends lexicographic comparison element by element; UKM106 documents the
	 * signed-word ordering defect.
	 *
	 * @param key1 First k-mer key
	 * @param key2 Second k-mer key
	 * @return Negative, zero, or positive value indicating relative ordering
	 */
	public static final int compare(long[] key1, long[] key2){
		//TODO: Probable bug (latent) - signed subtraction gives the wrong SIGN when words span
		//the signed-long boundary, which packed word[0] (32 bases, bit 63 live) can do; the
		//induced order is circular/non-transitive. Currently harmless: equality (dif==0) is
		//exact, and every BST consumer (HashForestU/KmerNodeU insert+find+resize-reinsert)
		//navigates with THIS same deterministic comparator, so paths stay self-consistent.
		//It becomes REAL if anything assumes a true total order: KmerNodeU.rebalance() (rebuilds
		//assuming sorted inorder; currently uncalled on these paths) or any sorted dump/merge.
		//Fix would be Long.compareUnsigned per word (matches Kmer.setKey0's convention) -- not
		//changed now to avoid altering tree shapes mid-validation. (Amber, 2026-08-19)
		for(int i=0; i<key1.length; i++){
			long dif=key1[i]-key2[i];
			if(dif!=0){return (int)Tools.mid(-1, dif, 1);}
		}
		return 0;
	}

	/**
	 * Tests equality of two k-mer keys represented as long arrays.
	 * @param key1 First k-mer key
	 * @param key2 Second k-mer key
	 * @return true if keys are equal, false otherwise
	 */
	public static final boolean equals(long[] key1, long[] key2){
		return compare(key1, key2)==0;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Fetches the value associated with a raw key and bucket-placement hash.
	 * @param key The k-mer key array
	 * @param xor placement hash used to choose the bucket or probe chain
	 * @return The value associated with the k-mer, or -1 if not present
	 */
	public abstract int getValue(long[] key, long xor);

//	/** Returns count */
//	public final int increment(long[] key){throw new RuntimeException();}
//
//	/** Returns number of entries created */
//	public final int incrementAndReturnNumCreated(final long[] key){throw new RuntimeException();}
//
//	public final int set(long[] key, int value){throw new RuntimeException();}
//
//	public final int set(long[] key, int[] vals){throw new RuntimeException();}
//
//	/** Returns number of kmers added */
//	public final int setIfNotPresent(long[] key, int value){throw new RuntimeException();}
//
//	/**
//	 * Fetch the value associated with a kmer.
//	 * @param kmer
//	 * @return A value.  -1 means the kmer was not present.
//	 */
//	final int getValue(long[] key){throw new RuntimeException();}
//
//	/**
//	 * Fetch the values associated with a kmer.
//	 * @param kmer
//	 * @param singleton A blank array of length 1.
//	 * @return An array filled with values.  Values of -1 are invalid.
//	 */
//	public final int[] getValues(long[] key, int[] singleton){throw new RuntimeException();}
//
//	public final boolean contains(long[] key){throw new RuntimeException();}

	/**
	 * Tests whether a k-mer contains a specific value.
	 * Only active in test mode.
	 *
	 * @param kmer The k-mer to test
	 * @param v The value to search for
	 * @return true if the k-mer contains the specified value
	 */
	public final boolean contains(Kmer kmer, int v){
		assert(TESTMODE);
		int[] set=getValues(kmer, new int[] {-1});
		if(set==null){return false;}
		for(int s : set){
			if(s==-1){break;}
			if(s==v){return true;}
		}
		return false;
	}

	/**
	 * Tests whether a k-mer contains all specified values.
	 * Only active in test mode.
	 *
	 * @param kmer The k-mer to test
	 * @param vals Array of values to search for
	 * @return true if the k-mer contains all specified values
	 */
	public final boolean contains(Kmer kmer, int[] vals){
		assert(TESTMODE);
		int[] set=getValues(kmer, new int[] {-1});
		if(set==null){return false;}
		boolean success=true;
		for(int v : vals){
			if(v==-1){break;}
			success=false;
			for(int s : set){
				if(s==v){
					success=true;
					break;
				}
			}
			if(!success){break;}
		}
		return success;
	}

	/** Rebalances the table when the implementation supports tree reshaping. */
	public abstract void rebalance();

	/** Returns this implementation's occupied key count, which may exclude victims. */
	public abstract long size();
	/** Returns the length of the implementation's primary storage. */
	public abstract int arrayLength();
	/** Returns true if the table supports rebalancing operations. */
	public abstract boolean canRebalance();

	/**
	 * Emits k-mers to a text stream writer with representation-specific filtering.
	 *
	 * @param tsw Text stream writer for output
	 * @param k K-mer length
	 * @param mincount Minimum count threshold (inclusive)
	 * @param maxcount Maximum count threshold (inclusive)
	 * @return implementation-specific emission status, not writer completion proof
	 */
	public abstract boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount);
	/**
	 * Emits k-mers to a byte stream writer with representation-specific filtering.
	 *
	 * @param bsw Byte stream writer for output
	 * @param k K-mer length
	 * @param mincount Minimum count threshold (inclusive)
	 * @param maxcount Maximum count threshold (inclusive)
	 * @param remaining Atomic counter tracking remaining items to dump
	 * @return implementation-specific emission status, not writer completion proof
	 */
	public abstract boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining);
	/**
	 * Emits k-mers as byte records, appending through a caller-owned
	 * thread-exclusive buffer.
	 *
	 * @param bsw Byte stream writer for output
	 * @param bb caller-thread-exclusive builder; the caller submits any final tail
	 * @param k K-mer length
	 * @param mincount Minimum count threshold (inclusive)
	 * @param maxcount Maximum count threshold (inclusive)
	 * @param remaining Atomic counter tracking remaining items to dump
	 * @return implementation-specific emission status, not writer completion proof
	 */
	public abstract boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, final int maxcount, AtomicLong remaining);

	/**
	 * Fills a histogram array with k-mer count frequencies.
	 * @param ca Count array to fill with histogram data
	 * @param max Maximum count value to include
	 */
	public abstract void fillHistogram(long[] ca, int max);
	/** Fills a SuperLongList with k-mer count histogram data.
	 * @param sll SuperLongList to populate with histogram data */
	public abstract void fillHistogram(SuperLongList sll);
	/**
	 * Counts GC content distribution across k-mers.
	 * @param gcCounts Array to fill with GC count data
	 * @param max Maximum GC count to track
	 */
	public abstract void countGC(long[] gcCounts, int max);

	/**
	 * Calculates GC count in a k-mer represented as a long.
	 * Counts bases with values 1 (C) or 2 (G).
	 * @param kmer K-mer encoded as bit-packed long
	 * @return Number of G and C bases in the k-mer
	 */
	public static final int gc(long kmer){
		return Long.bitCount((kmer^(kmer<<1))&0xAAAAAAAAAAAAAAAAL); //May be faster
//		int gc=0;
//		while(kmer>0){
//			long x=kmer&3;
//			kmer>>>=2;
//			if(x==1 || x==2){gc++;}
//		}
//		return gc;
	}

	/**
	 * Gets the object associated with a k-mer.
	 * @param kmer The k-mer key
	 * @return Object associated with the k-mer
	 */
	Object get(Kmer kmer){return get(kmer.key());}
	/**
	 * Gets the object associated with a k-mer key array.
	 * @param key The k-mer key array
	 * @return Object associated with the key
	 */
	abstract Object get(long[] key);
	/** Resizes the table to accommodate more entries. */
	abstract void resize();
	/** Returns true if the table supports resizing operations. */
	abstract boolean canResize();



	/**
	 * Removes entries with a value of {@code limit} or less.
	 * Rehashes the remainder.
	 * @return Number removed.
	 */
	abstract long regenerate(final int limit);

	/** Acquires this table's concrete lock. */
	final void lock(){getLock().lock();}
	/** Releases this table's concrete lock. */
	final void unlock(){getLock().unlock();}
	/**
	 * Attempts to acquire this table's concrete lock without blocking.
	 * @return true if the lock was acquired, false otherwise
	 */
	final boolean tryLock(){return getLock().tryLock();}
	/**
	 * Gets the lock object for this table.
	 * Default implementation throws RuntimeException.
	 * @return The lock object
	 * @throws RuntimeException if not implemented by subclass
	 */
	Lock getLock(){
		throw new RuntimeException("Unimplemented.");
	}

	/*--------------------------------------------------------------*/
	/*---------------       Allocation Methods      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Allocates an AtomicIntegerArray with memory management.
	 * @param len Length of the array to allocate
	 * @return Newly allocated AtomicIntegerArray
	 */
	final static AtomicIntegerArray allocAtomicInt(int len){
		return KillSwitch.allocAtomicInt(len);
	}

	/**
	 * Allocates a one-dimensional long array with memory management.
	 * @param len Length of the array to allocate
	 * @return Newly allocated long array
	 */
	final static long[] allocLong1D(int len){
		return KillSwitch.allocLong1D(len);
	}

	/**
	 * Allocates a two-dimensional long array with memory management.
	 * @param mult Number of sub-arrays
	 * @param len Length of each sub-array
	 * @return Newly allocated two-dimensional long array
	 */
	final static long[][] allocLong2D(int mult, int len){
		return KillSwitch.allocLong2D(mult, len);
	}

	/**
	 * Allocates a one-dimensional int array with memory management.
	 * @param len Length of the array to allocate
	 * @return Newly allocated int array
	 */
	final static int[] allocInt1D(int len){
		return KillSwitch.allocInt1D(len);
	}

	/**
	 * Allocates a two-dimensional int array with memory management.
	 * @param len Number of sub-arrays to allocate
	 * @return Newly allocated two-dimensional int array
	 */
	final static int[][] allocInt2D(int len){
		return KillSwitch.allocInt2D(len);
	}

	/**
	 * Allocates an array of KmerNodeU objects with out-of-memory handling.
	 * @param len Length of the array to allocate
	 * @return Newly allocated KmerNodeU array, or null if out of memory
	 */
	final static KmerNodeU[] allocKmerNodeArray(int len){
		KmerNodeU[] ret=null;
		try{
			ret=new KmerNodeU[len];
		}catch(OutOfMemoryError e){
			synchronized(killMessage){
				e.printStackTrace();
				System.err.println(killMessage);
				KillSwitch.killSilent();
			}
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	/*---------------       Ownership Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Claims this k-mer for a thread when {@code newOwner} exceeds its current owner.
	 * @param kmer The k-mer key
	 * @param newOwner candidate owner thread id
	 * @return resulting owner id
	 */
	public abstract int setOwner(Kmer kmer, int newOwner);

	/**
	 * Clears the owner when it matches the supplied id.
	 * @param kmer The k-mer key
	 * @param owner owner id to clear
	 * @return true when the owner was cleared
	 */
	public abstract boolean clearOwner(Kmer kmer, int owner);

	/**
	 * Returns the owner thread id for this key.
	 * @param kmer The k-mer key
	 * @return owner id, or -1 when unowned
	 */
	public abstract int getOwner(Kmer kmer);

	/** Creates any auxiliary storage needed for ownership tracking. */
	public abstract void initializeOwnership();

	/** Clears or drops ownership storage; later lookups are implementation-defined. */
	public abstract void clearOwnership();

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Converts a k-mer to its text representation.
	 * @param kmer The k-mer to convert
	 * @return StringBuilder containing the k-mer sequence
	 */
	public static final StringBuilder toText(Kmer kmer){
		return toText(kmer.key(), kmer.kbig);
	}

	/**
	 * Converts a k-mer array to text representation.
	 * @param array Long array containing k-mer data
	 * @param kbig K-mer length
	 * @return StringBuilder containing the k-mer sequence
	 */
	//Gate-1 packed routing (Amber, 2026-08-19): the length parameter of every multi-word
	//toText/toBytes below is the TOTAL kmer length (kbig), and each word is decoded at its own
	//width via Kmer.perWordK(kbig, array.length, i) -- under the symmetric layout
	//perWordK==kbig/mult==the per-word k these methods used to receive, so output is
	//byte-identical for every pre-packing caller; under PACKED the last word is partial and a
	//uniform width would emit phantom leading bases. Callers updated to pass kbig
	//(KmerTableSetU dump entry points, toText(Kmer) above).
	public static final StringBuilder toText(long[] array, int kbig){
		StringBuilder sb=new StringBuilder(kbig+2);
		for(int pos=0; pos<array.length; pos++){
			long kmer=array[pos];
			for(int i=Kmer.perWordK(kbig, array.length, pos)-1; i>=0; i--){
				int x=(int)((kmer>>(2*i))&3);
				sb.append((char)AminoAcid.numberToBase[x]);
			}
		}
		return sb;
	}

	/**
	 * Converts a k-mer array with count to text representation.
	 *
	 * @param array Long array containing k-mer data
	 * @param count Count value associated with the k-mer
	 * @param k K-mer length
	 * @return StringBuilder containing formatted k-mer and count
	 */
	static final StringBuilder toText(long[] array, int count, int k){
		StringBuilder sb=new StringBuilder(k+10);
		return toText(array, count, k, sb);
	}

	/**
	 * Converts a k-mer array with count to byte representation.
	 *
	 * @param array Long array containing k-mer data
	 * @param count Count value associated with the k-mer
	 * @param k K-mer length
	 * @return ByteBuilder containing formatted k-mer and count
	 */
	static final ByteBuilder toBytes(long[] array, int count, int k){
		ByteBuilder bb=new ByteBuilder(k+10);
		return toBytes(array, count, k, bb);
	}

	/**
	 * Converts a k-mer array with multiple values to text representation.
	 *
	 * @param array Long array containing k-mer data
	 * @param values Array of values associated with the k-mer
	 * @param k K-mer length
	 * @return StringBuilder containing formatted k-mer and values
	 */
	static final StringBuilder toText(long[] array, int[] values, int k){
		StringBuilder sb=new StringBuilder(k+10);
		return toText(array, values, k, sb);
	}

	/**
	 * Converts a k-mer array with multiple values to byte representation.
	 *
	 * @param array Long array containing k-mer data
	 * @param values Array of values associated with the k-mer
	 * @param k K-mer length
	 * @return ByteBuilder containing formatted k-mer and values
	 */
	static final ByteBuilder toBytes(long[] array, int[] values, int k){
		ByteBuilder bb=new ByteBuilder(k+10);
		return toBytes(array, values, k, bb);
	}

	/**
	 * Converts k-mer array with count to text using provided StringBuilder.
	 * Supports both FASTA and tabular output formats.
	 *
	 * @param array Long array containing k-mer data
	 * @param count Count value associated with the k-mer
	 * @param k K-mer length
	 * @param sb StringBuilder to append output to
	 * @return The StringBuilder with appended k-mer data
	 */
	static final StringBuilder toText(long[] array, int count, int k, StringBuilder sb){
		if(FASTA_DUMP){
			sb.append('>');
			sb.append(count);
			sb.append('\n');
			for(int i=0; i<array.length; i++){
				append(array[i], Kmer.perWordK(k, array.length, i), sb);//k here is kbig; see toText header comment
			}
		}else{
			for(int i=0; i<array.length; i++){
				append(array[i], Kmer.perWordK(k, array.length, i), sb);//k here is kbig; see toText header comment
			}
			sb.append('\t');
			sb.append(count);
		}
		return sb;
	}

	/**
	 * Converts k-mer array with values to text using provided StringBuilder.
	 * Supports both FASTA and tabular output formats.
	 *
	 * @param array Long array containing k-mer data
	 * @param values Array of values associated with the k-mer
	 * @param k K-mer length
	 * @param sb StringBuilder to append output to
	 * @return The StringBuilder with appended k-mer data
	 */
	static final StringBuilder toText(long[] array, int[] values, int k, StringBuilder sb){
		if(FASTA_DUMP){
			sb.append('>');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
			sb.append('\n');
			for(int i=0; i<array.length; i++){
				append(array[i], Kmer.perWordK(k, array.length, i), sb);//k here is kbig; see toText header comment
			}
		}else{
			for(int i=0; i<array.length; i++){
				append(array[i], Kmer.perWordK(k, array.length, i), sb);//k here is kbig; see toText header comment
			}
			sb.append('\t');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
		}
		return sb;
	}

	/**
	 * Appends a k-mer sequence to a StringBuilder.
	 * Decodes bit-packed k-mer into nucleotide characters.
	 *
	 * @param kmer Bit-packed k-mer as long
	 * @param k K-mer length
	 * @param sb StringBuilder to append to
	 */
	private static final void append(long kmer, int k, StringBuilder sb){
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(2*i))&3);
			sb.append((char)AminoAcid.numberToBase[x]);
		}
	}

	/**
	 * Converts k-mer array with count to bytes using provided ByteBuilder.
	 * Supports both FASTA and tabular output formats.
	 *
	 * @param array Long array containing k-mer data
	 * @param count Count value associated with the k-mer
	 * @param k K-mer length
	 * @param sb ByteBuilder to append output to
	 * @return The ByteBuilder with appended k-mer data
	 */
	public static final ByteBuilder toBytes(long[] array, long count, int k, ByteBuilder sb){
		if(FASTA_DUMP){
			sb.append('>');
			sb.append(count);
			sb.append('\n');
			for(int i=0; i<array.length; i++){
				append(array[i], Kmer.perWordK(k, array.length, i), sb);//k here is kbig; see toText header comment
			}
		}else{
			for(int i=0; i<array.length; i++){
				append(array[i], Kmer.perWordK(k, array.length, i), sb);//k here is kbig; see toText header comment
			}
			sb.append('\t');
			sb.append(count);
		}
		return sb;
	}

	/**
	 * Converts k-mer array with values to bytes using provided ByteBuilder.
	 * Supports both FASTA and tabular output formats.
	 *
	 * @param array Long array containing k-mer data
	 * @param values Array of values associated with the k-mer
	 * @param k K-mer length
	 * @param sb ByteBuilder to append output to
	 * @return The ByteBuilder with appended k-mer data
	 */
	public static final ByteBuilder toBytes(long[] array, int[] values, int k, ByteBuilder sb){
		if(FASTA_DUMP){
			sb.append('>');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
			sb.append('\n');
			for(int i=0; i<array.length; i++){
				append(array[i], Kmer.perWordK(k, array.length, i), sb);//k here is kbig; see toText header comment
			}
		}else{
			for(int i=0; i<array.length; i++){
				append(array[i], Kmer.perWordK(k, array.length, i), sb);//k here is kbig; see toText header comment
			}
			sb.append('\t');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
		}
		return sb;
	}

	/**
	 * Appends a k-mer sequence to a ByteBuilder.
	 * Decodes bit-packed k-mer into nucleotide characters.
	 *
	 * @param kmer Bit-packed k-mer as long
	 * @param k K-mer length
	 * @param sb ByteBuilder to append to
	 */
	private static final void append(long kmer, int k, ByteBuilder sb){
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(2*i))&3);
			sb.append((char)AminoAcid.numberToBase[x]);
		}
	}


//	static void appendKmerText(long kmer, int count, int k, StringBuilder sb){
//		sb.setLength(0);
//		toText(kmer, count, k, sb);
//		sb.append('\n');
//	}

	/**
	 * Appends k-mer text with count to ByteBuilder and adds newline.
	 * Resets ByteBuilder length before appending.
	 *
	 * @param array Long array containing k-mer data
	 * @param count Count value associated with the k-mer
	 * @param k K-mer length
	 * @param bb ByteBuilder to use for output
	 */
	static void appendKmerText(long[] array, int count, int k, ByteBuilder bb){
		bb.setLength(0);
		toBytes(array, count, k, bb);
		bb.nl();
	}


	/** Unsupported base hook for flushing buffered inserts into this table. */
	long flush(){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * This allocates the data structures in multiple threads.  Unfortunately,
	 * it does not lead to any speedup, at least for ARRAY type.
	 * @param ways number of tables to allocate
	 * @param tableType implementation constant to instantiate
	 * @param schedule size schedule for every table
	 * @param k standard k-mer length
	 * @param kbig extended k-mer length for long-key tables
	 * @return preallocated tables
	 */
	public static final AbstractKmerTableU[] preallocate(int ways, int tableType, int[] schedule, int k, int kbig){

		final AbstractKmerTableU[] tables=new AbstractKmerTableU[ways];

		{
			final int t=Tools.max(1, Tools.min(Shared.threads(), 2, ways));
			final AllocThread[] allocators=new AllocThread[t];
			for(int i=0; i<t; i++){
				allocators[i]=new AllocThread(tableType, schedule, i, t, k, kbig, tables);
			}
			for(AllocThread at : allocators){at.start();}
			for(AllocThread at : allocators){
				while(at.getState()!=Thread.State.TERMINATED){
					try{
						at.join();
					}catch(InterruptedException e){
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}

		synchronized(tables){
			for(int i=0; i<tables.length; i++){
				final AbstractKmerTableU akt=tables[i];
				if(akt==null){
					throw new RuntimeException("KmerTable allocation failed, probably due to lack of RAM: "+i+", "+tables.length);
				}
			}
		}

		return tables;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Worker thread for parallel allocation of k-mer tables.
	 * Each thread allocates a subset of tables based on modulo assignment.
	 * Supports different table types including arrays, forests, and hybrid structures.
	 */
	private static class AllocThread extends Thread{

		/**
		 * Creates an allocation thread for k-mer table creation.
		 *
		 * @param type_ Table type constant (ARRAY1D, FOREST1D, etc.)
		 * @param schedule_ Size schedule for table growth
		 * @param mod_ Modulo offset for this thread's work assignment
		 * @param div_ Thread count divisor for work distribution
		 * @param k_ Standard k-mer length
		 * @param kbig_ Extended k-mer length
		 * @param tables_ Array to store allocated tables
		 */
		AllocThread(int type_, int[] schedule_, int mod_, int div_,
				int k_, int kbig_, AbstractKmerTableU[] tables_){
			type=type_;
			schedule=schedule_;
			size=schedule[0];
			mod=mod_;
			div=div_;
			growable=schedule.length>1;
			tables=tables_;
			k=k_;
			kbig=kbig_;
		}

		@Override
		public void run(){
			for(int i=mod; i<tables.length; i+=div){
//				System.err.println("T"+i+" allocating "+i);
				final AbstractKmerTableU akt;
				if(type==FOREST1D){
					akt=new HashForestU(size, k, kbig, growable, false);
				}else if(type==ARRAY1D){
					akt=new HashArrayU1D(schedule, k, kbig);
//					akt=new HashArrayU1D(size, k, kbig, growable);
				}else if(type==NODE1D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
//					akt=new KmerNode2(-1, 0);
				}else if(type==FOREST2D){
					akt=new HashForestU(size, k, kbig, growable, true);
				}else if(type==ARRAY2D){
					akt=new HashArrayU2D(schedule, k, kbig);
//					akt=new HashArrayU2D(size, k, kbig, growable);
				}else if(type==NODE2D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
//					akt=new KmerNode(-1, 0);
				}else if(type==ARRAYH){
					akt=new HashArrayUHybrid(schedule, k, kbig);
//					akt=new HashArrayUHybrid(size, k, kbig, growable);
				}else if(type==ARRAYHASH1D){
					akt=new HashArrayH1D(schedule);
				}else if(type==ARRAYHASH1D_FIXED){
					akt=new HashArrayH1D(schedule, false);
				}else{
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
				}
				synchronized(tables){
					tables[i]=akt;
				}
//				System.err.println("T"+i+" allocated "+i);
			}
		}

		/** Table type constant for this allocation thread */
		private final int type;
		/** Size schedule array defining table growth parameters */
		private final int[] schedule;
		/** Initial size for table allocation from schedule[0] */
		private final int size;
		/** Modulo offset for work assignment among threads */
		private final int mod;
		/** Thread count divisor for distributing table allocation work */
		private final int div;
		/** Standard k-mer length for table creation */
		private final int k;
		/** Extended k-mer length for certain table implementations */
		private final int kbig;
		/** Whether allocated tables support dynamic growth */
		private final boolean growable;
		/** Array to store the allocated k-mer tables */
		final AbstractKmerTableU[] tables;

	}

	/**
	 * Creates a walker for iterating through table contents.
	 * Default implementation throws RuntimeException.
	 * @return WalkerU for table iteration
	 * @throws RuntimeException if not implemented by subclass
	 */
	public WalkerU walk(){
		throw new RuntimeException("Unimplemented");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Controls whether k-mer output uses FASTA format instead of tabular output. */
	public static boolean FASTA_DUMP=true;
	/** Legacy dump flag; serializers in this class ignore it. */
	public static boolean NUMERIC_DUMP=false;

	/** Enables verbose output mode. */
	public static final boolean verbose=false; //slow
	/** Enables test-mode assertion checks. */
	public static final boolean TESTMODE=false; //slow

	/**
	 * Table type constants: unknown, one-dimensional array/forest/node,
	 * two-dimensional array/forest/node, hybrid array, and long-key hash-pair arrays.
	 */
	public static final int UNKNOWN=0, ARRAY1D=1, FOREST1D=2, NODE1D=4, ARRAY2D=5, FOREST2D=6, NODE2D=8,
			ARRAYH=9, ARRAYHASH1D=10, ARRAYHASH1D_FIXED=11;

	/** Sentinel return values for absent keys and primary hash collisions. */
	public static final int NOT_PRESENT=-1, HASH_COLLISION=-2;
	/** Sentinel owner id used when no thread owns the k-mer. */
	public static final int NO_OWNER=-1;

	/** Error message displayed when program runs out of memory */
	private static final String killMessage=new String("\nThis program ran out of memory.  Try increasing the -Xmx flag and setting prealloc.");

}
