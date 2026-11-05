package kmer;

import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Common operations for short-key tables, tree nodes and buffered table adapters.
 * Keys are nonnegative packed longs. Count tables store one integer per key;
 * multi-value tables merge integer values instead. Not every implementation
 * supports counting, resizing, ownership or every output operation.
 * Mutating a shared table requires caller coordination; the available table lock
 * is not acquired automatically by each operation. Borrowed arrays and walkers
 * must not be used across concurrent mutations or resizing.
 * @author Brian Bushnell
 * @date Oct 23, 2013
 *
 */
public abstract class AbstractKmerTable implements KmerTableInterface{

	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/

//	/** Returns count */
//	public final int increment(long kmer){return increment(kmer, 1);}

	/** Adds a positive increment and returns the resulting count, where supported. */
	public abstract int increment(final long kmer, final int incr);

//	/** Returns number of entries created */
//	public final int incrementAndReturnNumCreated(final long kmer){return incrementAndReturnNumCreated(kmer, 1);}

	/** Adds a positive increment and reports newly created keys; buffers report flushed batches. */
	public abstract int incrementAndReturnNumCreated(final long kmer, final int incr);

	/**
	 * Stores a value for a k-mer. Count tables replace the count; multi-value
	 * tables add the value to the key's list instead.
	 * @param kmer The k-mer key (bit-packed)
	 * @param value The value to set
	 * @return Number of new keys reported by the operation; buffers report flushed batches
	 */
	public abstract int set(long kmer, int value);

//	public abstract int set(long kmer, int[] vals);

	/**
	 * Stores a value list where supported, including IntList3-backed hybrid tables.
	 * Implementations may retain the supplied array; callers must preserve its contents.
	 * @param kmer Packed key
	 * @param vals Backing array containing the supplied values
	 * @param vlen Logical number of supplied values; some implementations instead scan sentinel padding
	 * @return Number of new keys reported by the operation
	 */
	public abstract int set(long kmer, int[] vals, int vlen);

	/** Inserts only an absent key where supported; returns the implementation's new-key count. */
	public abstract int setIfNotPresent(long kmer, int value);

	/**
	 * Fetches the count or first stored value associated with a k-mer.
	 * @param kmer Packed key
	 * @return Stored value, or -1 when absent according to this interface's contract
	 */
	public abstract int getValue(long kmer);

	/**
	 * Fetches the values associated with a k-mer. The result may be the caller's
	 * singleton buffer or a borrowed internal array; it is not a defensive copy.
	 * @param kmer Packed key
	 * @param singleton Reusable array of length one for single-value results
	 * @return Stored or singleton-buffer array; negative entries are invalid, and absence may return null
	 */
	public abstract int[] getValues(long kmer, int[] singleton);

	/**
	 * Tests whether the table contains a k-mer.
	 * @param kmer The k-mer key (bit-packed)
	 * @return true if the k-mer is present
	 */
	public abstract boolean contains(long kmer);

	/**
	 * Tests whether the table contains a k-mer with a specific value.
	 * Diagnostic helper requiring TESTMODE when assertions are enabled.
	 *
	 * @param kmer The k-mer key (bit-packed)
	 * @param v The value to search for
	 * @return true if the k-mer is present with the specified value
	 */
	public final boolean contains(long kmer, int v){
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
	 * Tests whether the table contains a k-mer with all specified values.
	 * Diagnostic helper requiring TESTMODE when assertions are enabled.
	 *
	 * @param kmer The k-mer key (bit-packed)
	 * @param vals Array of values to search for
	 * @return true if the k-mer is present with all specified values
	 */
	public final boolean contains(long kmer, int[] vals){
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

	/** Requests implementation-specific rebalancing; unsupported implementations may throw. */
	public abstract void rebalance();

	/** Returns this structure's key count; HashArray excludes its separately counted victims. */
	public abstract long size();
	/** Returns implementation-specific array capacity, including any allocated probe tail. */
	public abstract int arrayLength();
	/** Returns true if this table supports rebalancing operations */
	public abstract boolean canRebalance();

	/**
	 * Dumps k-mers as text format to a stream writer.
	 *
	 * @param tsw Text stream writer for output
	 * @param k K-mer length for sequence reconstruction
	 * @param mincount Minimum count threshold
	 * @param maxcount Maximum count threshold
	 * @return Implementation status; does not establish asynchronous writer completion or success
	 */
	public abstract boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount);
	/**
	 * Dumps k-mers as bytes to a stream writer.
	 *
	 * @param bsw Byte stream writer for output
	 * @param k K-mer length for sequence reconstruction
	 * @param mincount Minimum count threshold
	 * @param maxcount Maximum count threshold
	 * @param remaining Optional shared output-record allowance, consumed by supported implementations
	 * @return Implementation status; does not establish asynchronous writer completion or success
	 */
	public abstract boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining);
	/**
	 * Appends k-mer records using a caller-owned buffer, flushing blocks to the writer.
	 * This operation supports externally partitioned work; it does not itself launch workers.
	 *
	 * @param bsw Byte stream writer for output
	 * @param bb Exclusive caller buffer; caller submits the final partial block and finishes the writer
	 * @param k K-mer length for sequence reconstruction
	 * @param mincount Minimum count threshold
	 * @param maxcount Maximum count threshold
	 * @param remaining Optional shared output-record allowance, consumed by supported implementations
	 * @return Implementation status; caller retains any unflushed buffer contents
	 */
	public abstract boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining);

	/**
	 * Fills a histogram array with k-mer count frequencies.
	 * @param ca Histogram array to fill
	 * @param max Highest histogram bin; larger counts accumulate in this bin
	 */
	public abstract void fillHistogram(long[] ca, int max);
	/** Adds stored counts to the supplied histogram accumulator.
	 * @param sll Accumulator receiving counts */
	public abstract void fillHistogram(SuperLongList sll);
	/**
	 * Adds each visited key's number of G/C bases to its stored-count bin.
	 * @param gcCounts Accumulated G/C-base totals, indexed by count rather than GC content
	 * @param max Highest count bin; larger counts accumulate in this bin
	 */
	public abstract void countGC(long[] gcCounts, int max);

	/**
	 * Counts G/C symbols in a two-bit nucleotide encoding; not an amino-acid calculation.
	 * There is no length argument, so unused bits must be zero (which encode A).
	 * @param kmer The k-mer to analyze (bit-packed)
	 * @return Number of G and C bases
	 */
	public static final int gc(long kmer){
		//Comprehension (clever, verified): for 2-bit symbols A=00,C=01,G=10,T=11, a base is G or C iff its two bits DIFFER
		//(b1^b0: C=01->1, G=10->1, A/T->0). (kmer^(kmer<<1)) puts b1^b0 at each symbol's high bit; the 0xAA... mask keeps
		//one such bit per symbol; bitCount = #GC. NUCLEOTIDE-only (assumes 2-bit); meaningless under AMINO_IN 5-bit, but GC
		//is a nucleotide concept; callers must supply nucleotide encodings.
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
	 * Gets an implementation-specific live node/object where supported; other types throw.
	 * @param kmer The k-mer key (bit-packed)
	 * @return The raw object, or null if not present
	 */
	abstract Object get(long kmer);
	/** Requests implementation-specific growth; may leave capacity unchanged or reject exhaustion. */
	abstract void resize();
	/** Returns true if this table supports resizing operations */
	abstract boolean canResize();

	/**
	 * Removes entries with a value of the limit or less.
	 * Rehashes the remainder.
	 * @param limit Inclusive maximum value to remove
	 * @return Number removed where supported; unsupported implementations may reject the operation
	 */
	abstract long regenerate(int limit);

	/** Acquires the lock for this table */
	final void lock(){getLock().lock();}
	/** Releases the lock for this table */
	final void unlock(){getLock().unlock();}
	/** Attempts to acquire the lock without blocking.
	 * @return true if the lock was acquired */
	final boolean tryLock(){return getLock().tryLock();}
	/**
	 * Gets the lock object for this table.
	 * @return The lock object
	 * @throws RuntimeException If not implemented by subclass
	 */
	Lock getLock(){
		throw new RuntimeException("Unimplemented.");
	}

	/*--------------------------------------------------------------*/
	/*---------------       Allocation Methods      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Delegates atomic-array allocation and out-of-memory handling to KillSwitch.
	 * @param len Array length
	 * @return Allocated AtomicIntegerArray
	 */
	final static AtomicIntegerArray allocAtomicInt(int len){
		return KillSwitch.allocAtomicInt(len);
	}

	/**
	 * Delegates long-array allocation and out-of-memory handling to KillSwitch.
	 * @param len Array length
	 * @return Allocated long array
	 */
	final static long[] allocLong1D(int len){
		return KillSwitch.allocLong1D(len);
	}

	/**
	 * Delegates rectangular long-array allocation and out-of-memory handling to KillSwitch.
	 * @param mult Number of sub-arrays
	 * @param len Length of each sub-array
	 * @return Allocated 2D long array
	 */
	final static long[][] allocLong2D(int mult, int len){
		return KillSwitch.allocLong2D(mult, len);
	}

	/**
	 * Delegates int-array allocation and out-of-memory handling to KillSwitch.
	 * @param len Array length
	 * @return Allocated int array
	 */
	final static int[] allocInt1D(int len){
		return KillSwitch.allocInt1D(len);
	}

	/**
	 * Allocates only the outer int-array dimension through KillSwitch; rows remain null.
	 * @param len Array length for first dimension
	 * @return Allocated 2D int array
	 */
	final static int[][] allocInt2D(int len){
		return KillSwitch.allocInt2D(len);
	}

	/**
	 * Allocates a KmerNode array with out-of-memory protection.
	 * @param len Array length
	 * @return Allocated node array; out-of-memory handling prints diagnostics and halts the VM
	 */
	final static KmerNode[] allocKmerNodeArray(int len){
		KmerNode[] ret=null;
		try{
			ret=new KmerNode[len];
		}catch(OutOfMemoryError e){
			synchronized(killMessage){
				e.printStackTrace();
				System.err.println(killMessage);
//				Shared.printMemory();
				KillSwitch.killSilent();
			}
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	/*---------------       Ownership Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Raises an existing key's owner where supported, returning an observed/installed owner.
	 * Requires ownership initialization and stable structure; synchronization varies by implementation. */
	public abstract int setOwner(long kmer, int newOwner);

	/** Clears a matching owner of an existing key, returning whether the comparison succeeded. */
	public abstract boolean clearOwner(long kmer, int owner);

	/** Returns an existing key's owner or -1 for unowned; this is not an absent-key lookup. */
	public abstract int getOwner(long kmer);

	/** Initializes or resets ownership under exclusive access; implementation support varies. */
	public abstract void initializeOwnership();

	/** Releases ownership storage or resets owner values to -1 under exclusive access. */
	public abstract void clearOwnership();

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a sequence builder using Shared.AMINO_IN to choose two- or five-bit symbols.
	 * Appends no record terminator; k must fit the selected long encoding.
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param k K-mer length
	 * @return StringBuilder containing the k-mer sequence
	 */
	public static final StringBuilder toText(long kmer, int k){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		StringBuilder sb=new StringBuilder(k);
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(bits*i))&mask);
			sb.append((char)lookup[x]);
		}
		return sb;
	}

	/**
	 * Converts a k-mer and count to text format.
	 *
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param count Count value associated with k-mer
	 * @param k K-mer length
	 * @return StringBuilder containing formatted output
	 */
	static final StringBuilder toText(long kmer, int count, int k){
		StringBuilder sb=new StringBuilder(k+10);
		return toText(kmer, count, k, sb);
	}

	/**
	 * Converts a k-mer and count to byte format.
	 *
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param count Count value associated with k-mer
	 * @param k K-mer length
	 * @return ByteBuilder containing formatted output
	 */
	static final ByteBuilder toBytes(long kmer, int count, int k){
		ByteBuilder bb=new ByteBuilder(k+10);
		return toBytes(kmer, count, k, bb);
	}

	/**
	 * Converts a k-mer and multiple values to text format.
	 *
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param values Nonnull value array; formatting stops at the first -1 or array end
	 * @param k K-mer length
	 * @return StringBuilder containing formatted output
	 */
	static final StringBuilder toText(long kmer, int[] values, int k){
		StringBuilder sb=new StringBuilder(k+10);
		return toText(kmer, values, k, sb);
	}

	/**
	 * Converts a k-mer and multiple values to byte format.
	 *
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param values Nonnull value array; formatting stops at the first -1 or array end
	 * @param k K-mer length
	 * @return ByteBuilder containing formatted output
	 */
	static final ByteBuilder toBytes(long kmer, int[] values, int k){
		ByteBuilder bb=new ByteBuilder(k+10);
		return toBytes(kmer, values, k, bb);
	}

	/**
	 * Converts a k-mer and count to text format using provided StringBuilder.
	 * Appends FASTA when FASTA_DUMP is set, otherwise sequence and tab-separated values.
	 * NUMERIC_DUMP does not affect this StringBuilder path. Appends no final newline.
	 *
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param count Count value associated with k-mer
	 * @param k K-mer length
	 * @param sb StringBuilder to append to
	 * @return The provided StringBuilder with formatted output
	 */
	static final StringBuilder toText(long kmer, int count, int k, StringBuilder sb){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		if(FASTA_DUMP){
			sb.append('>');
			sb.append(count);
			sb.append('\n');
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				sb.append((char)lookup[x]);
			}
		}else{
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				sb.append((char)lookup[x]);
			}
			sb.append('\t');
			sb.append(count);
		}
		return sb;
	}

	/**
	 * Converts a k-mer and multiple values to text format using provided StringBuilder.
	 * Appends FASTA when FASTA_DUMP is set, otherwise sequence and tab-separated values.
	 * NUMERIC_DUMP does not affect this StringBuilder path. Appends no final newline.
	 *
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param values Nonnull value array; formatting stops at the first -1 or array end
	 * @param k K-mer length
	 * @param sb StringBuilder to append to
	 * @return The provided StringBuilder with formatted output
	 */
	static final StringBuilder toText(long kmer, int[] values, int k, StringBuilder sb){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		if(FASTA_DUMP){
			sb.append('>');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
			sb.append('\n');
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				sb.append((char)lookup[x]);
			}
		}else{
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				sb.append((char)lookup[x]);
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
	 * Converts a k-mer and count to byte format using provided ByteBuilder.
	 * Appends FASTA first, otherwise hexadecimal keys when NUMERIC_DUMP is set,
	 * otherwise sequence and tab-separated values. Appends no final newline.
	 *
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param count Count value associated with k-mer
	 * @param k K-mer length
	 * @param bb ByteBuilder to append to
	 * @return The provided ByteBuilder with formatted output
	 */
	public static final ByteBuilder toBytes(long kmer, long count, int k, ByteBuilder bb){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		if(FASTA_DUMP){
			bb.append('>');
			bb.append(count);
			bb.nl();
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				bb.append(lookup[x]);
			}
//			assert(false) : kmer+"->\n"+bb+"\n"+AminoAcid.kmerToStringAA(kmer, k);
		}else if(NUMERIC_DUMP){
			bb.append(Long.toHexString(kmer));
			bb.tab();
			bb.append(count);
		}else{
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				bb.append(lookup[x]);
			}
			bb.tab();
			bb.append(count);
		}
		return bb;
	}

	/**
	 * Converts a k-mer and multiple values to byte format using provided ByteBuilder.
	 * Appends FASTA first, otherwise hexadecimal keys when NUMERIC_DUMP is set,
	 * otherwise sequence and tab-separated values. Appends no final newline.
	 *
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param values Nonnull value array; formatting stops at the first -1 or array end
	 * @param k K-mer length
	 * @param bb ByteBuilder to append to
	 * @return The provided ByteBuilder with formatted output
	 */
	public static final ByteBuilder toBytes(long kmer, int[] values, int k, ByteBuilder bb){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		if(FASTA_DUMP){
			bb.append('>');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){bb.append(',');}
				bb.append(x);
			}
			bb.nl();
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				bb.append(lookup[x]);
			}
		}else if(NUMERIC_DUMP){
			bb.append(Long.toHexString(kmer));
			bb.tab();
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){bb.append(',');}
				bb.append(x);
			}
		}else{
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				bb.append(lookup[x]);
			}
			bb.tab();
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){bb.append(',');}
				bb.append(x);
			}
		}
		return bb;
	}

//	static void appendKmerText(long kmer, int count, int k, StringBuilder sb){
//		sb.setLength(0);
//		toText(kmer, count, k, sb);
//		sb.append('\n');
//	}

	/**
	 * Appends k-mer text with newline to a ByteBuilder.
	 * Clears the ByteBuilder before appending.
	 *
	 * @param kmer The k-mer to convert (bit-packed)
	 * @param count Count value associated with k-mer
	 * @param k K-mer length
	 * @param bb ByteBuilder to use for output
	 */
	static void appendKmerText(long kmer, int count, int k, ByteBuilder bb){
		bb.setLength(0);
		toBytes(kmer, count, k, bb);
		bb.nl();
	}

	/** Flushes a buffered implementation where supported; the base implementation throws. */
	long flush(){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Allocates independent table shards with at most two workers and joins them all.
	 * Workers take disjoint strided indices. A missing shard after termination causes
	 * an exception; worker exceptions themselves are not forwarded with their cause.
	 * @param ways Number of table shards to allocate
	 * @param tableType Supported table-family constant; node-only/unknown types are rejected
	 * @param schedule Nonempty borrowed capacity schedule; callers must not mutate it
	 * @param mask Hash mask passed to array-backed implementations
	 * @return Fully populated table array, possibly empty when ways is zero
	 */
	public static final AbstractKmerTable[] preallocate(int ways, int tableType, int[] schedule, long mask){

		final AbstractKmerTable[] tables=new AbstractKmerTable[ways];

		{
			shared.Timer tm=new shared.Timer();
			final int t=Tools.max(1, Tools.min(Shared.threads(), 2, ways)); //More than 2 still improves allocation time, but only slightly; ~25% faster at t=4.
			final AllocThread[] allocators=new AllocThread[t];
			for(int i=0; i<t; i++){
				allocators[i]=new AllocThread(tableType, schedule, i, t, mask, tables);
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
			tm.stop();
			if(AbstractKmerTableSet.DISPLAY_PROGRESS){System.err.println(tm);}
		}

		synchronized(tables){
			for(int i=0; i<tables.length; i++){
				final AbstractKmerTable akt=tables[i];
				if(akt==null){
					throw new RuntimeException("KmerTable allocation failed, probably due to lack of RAM: "+i+", "+tables.length);
				}
			}
		}

		return tables;
	}

	/**
	 * Creates a bitmask that zeros out the middle symbol(s) of a k-mer.
	 * Used for fuzzy k-mer matching by ignoring the center base(s).
	 *
	 * @param k Positive k-mer length fitting the selected packed encoding
	 * @param amino true for amino acid sequences, false for nucleotides
	 * @return Bitmask with middle bits zeroed
	 */
	public static final long makeMiddleMask(int k, boolean amino){
		final boolean odd=((k&1)==1);
		final int bitsPerSymbol=(amino ? 5 : 2);
		final int shift=bitsPerSymbol*((k-1)/2);
		final long middle=~((-1L)<<(odd ? bitsPerSymbol : 2*bitsPerSymbol));
		return ~(middle<<shift);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Worker thread for parallel k-mer table allocation.
	 * Allocates a disjoint strided subset and publishes each result under the array monitor.
	 * @author Brian Bushnell
	 * @date Oct 23, 2013
	 */
	private static class AllocThread extends Thread{

		/**
		 * Creates allocation thread for specific table subset.
		 *
		 * @param type_ Table type to create
		 * @param schedule_ Size schedule for table growth
		 * @param mod_ Starting index modulo for this thread
		 * @param div_ Thread count divisor for work distribution
		 * @param mask_ Bit mask for table operations
		 * @param tables_ Shared array to store created tables
		 */
		AllocThread(int type_, int[] schedule_, int mod_, int div_,
				long mask_, AbstractKmerTable[] tables_){
			type=type_;
			schedule=schedule_;
			size=schedule[0];
			mod=mod_;
			div=div_;
			mask=mask_;
			growable=schedule.length>1;
			tables=tables_;
		}

		/** Constructs supported table families; failed slots remain null for the coordinator check. */
		@Override
		public void run(){
			//Initialize tables
			long sum=0;

//			Shared.printMemory();}
			for(int i=mod; i<tables.length; i+=div){
//				System.err.println("T"+i+" allocating "+i);
				final AbstractKmerTable akt;
				if(type==FOREST1D){
					akt=new HashForest(size, growable, false);
				}else if(type==TABLE){
					akt=new KmerTable(size, growable);
				}else if(type==ARRAY1D){
					akt=new HashArray1D(schedule, mask);
//					long len=((HashArray1D)akt).arrayLength();
//					long mem=((HashArray1D)akt).calcMem();
//					sum+=mem;
//					akt=new HashArray1D(size, -1, mask, growable);//TODO: Set maxSize
				}else if(type==NODE1D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
//					akt=new KmerNode2(-1, 0);
				}else if(type==FOREST2D){
					akt=new HashForest(size, growable, true);
				}else if(type==TABLE2D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
				}else if(type==ARRAY2D){
					akt=new HashArray2D(schedule, mask);
//					akt=new HashArray2D(size, -1, mask, growable);//TODO: Set maxSize
				}else if(type==NODE2D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
//					akt=new KmerNode(-1, 0);
				}else if(type==ARRAYH){
					akt=new HashArrayHybrid(schedule, mask);
//					akt=new HashArrayHybrid(size, -1, mask, growable);//TODO: Set maxSize
				}else if(type==ARRAYHF){
					akt=new HashArrayHybridFast(schedule, mask);
//					akt=new HashArrayHybrid(size, -1, mask, growable);//TODO: Set maxSize
				}else{
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
				}
				synchronized(tables){
					tables[i]=akt;
				}
//				System.err.println("T"+i+" allocated "+i);

//				if(i%100==0){
//					System.err.println("Allocated "+(sum/1000000)+"MB");
//					Shared.printMemory();
//				}
			}
//			Shared.printMemory();
//			assert(false) : ("Allocated "+(sum/1000000)+"MB");
		}

		/** Table type to allocate */
		private final int type;
		/** Size schedule for table growth */
		private final int[] schedule;
		/** Initial table size from schedule */
		private final int size;
		/** Starting index modulo for this thread */
		private final int mod;
		/** Thread count divisor for work distribution */
		private final int div;
		/** Bit mask for table operations */
		private final long mask;
		/** Whether tables can grow beyond initial size */
		private final boolean growable;
		/** Shared array to store created tables */
		final AbstractKmerTable[] tables;

	}

	/**
	 * Creates a walker for iterating through k-mer entries.
	 * @return Walker instance for this table
	 * @throws RuntimeException If not implemented by subclass
	 */
	public Walker walk(){
		throw new RuntimeException("Unimplemented");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Controls whether k-mer output uses FASTA format */
	public static boolean FASTA_DUMP=true;
	/** Selects hexadecimal keys in byte formatters when FASTA_DUMP is false. */
	public static boolean NUMERIC_DUMP=false;
	/** Enables high-count-first reinsertion in implementations that consult this flag. */
	public static boolean TWO_PASS_RESIZE=false;

	/** Enables verbose output for debugging */
	public static final boolean verbose=false;
	/** Enables test mode with additional validation (slow) */
	public static final boolean TESTMODE=false; //123 SLOW!

	/** Table-family identifiers. ARRAYH/ARRAYHF use IntList2/IntList3 respectively;
	 * UNKNOWN, NODE1D, TABLE2D and NODE2D are not supported by preallocate. */
	public static final int UNKNOWN=0, ARRAY1D=1, FOREST1D=2, TABLE=3, NODE1D=4, ARRAY2D=5, FOREST2D=6, TABLE2D=7, NODE2D=8, ARRAYH=9, ARRAYHF=10;

	/** Lookup sentinels: absent key, or exhausted primary probe window requiring victim lookup. */
	public static final int NOT_PRESENT=-1, HASH_COLLISION=-2;
	/** Constant indicating no thread owns the k-mer */
	public static final int NO_OWNER=-1;

	/** Error message displayed when program runs out of memory */
	private static final String killMessage=new String("\nThis program ran out of memory.  Try increasing the -Xmx flag and setting prealloc.");

}
