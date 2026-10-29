package tax;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import shared.Timer;
import structures.IntList;
import structures.LongList;

/**
 * Disk-backed accession-to-TaxID table using a memory-mapped
 * open-addressing hash table, partitioned for build efficiency.
 *
 * Each slot is 12 bytes: [8-byte long key][4-byte int value].
 * Empty slots have key=0 (EMPTY sentinel).
 *
 * Partitioned file format:
 *   Prefix (24 bytes):
 *     [long totalCapacity] [long entryCount] [int numPartitions] [int reserved]
 *   Partition table (numPartitions x 16 bytes):
 *     [long partCapacity] [long partOffset]
 *   Padding (8 bytes):
 *     [long unused]
 *   Body:
 *     numPartitions x (partCapacity x 12-byte slots)
 *
 * Lookup: hash the key, select partition, probe within that partition.
 *
 * @author Brian Bushnell, Chloe
 */
public class DiskAccessionTable{

	/*--------------------------------------------------------------*/
	/*----------------        Construction         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens an existing accession table and maps each partition for lookup.
	 *
	 * @param path table file path
	 * @throws IOException if the table cannot be opened or mapped
	 */
	public DiskAccessionTable(String path) throws IOException{
		file=new File(path);
		if(!file.exists()){throw new IOException("DiskAccessionTable file not found: "+path);}

		raf=new RandomAccessFile(file, "r");
		channel=raf.getChannel();
		boolean success=false;
		try{
			long totalCapacity=raf.readLong();
			entryCount=raf.readLong();
			numPartitions=raf.readInt();
			raf.readInt(); //reserved

			if(numPartitions<1 || numPartitions>(1<<20)){
				//Legacy non-partitioned format
				numPartitions=1;
				partCapacities=new long[]{totalCapacity};
				partOffsets=new long[]{HEADER_SIZE_LEGACY};
				partSegments=new MappedByteBuffer[1][];
				partSegments[0]=mapSegments(channel, HEADER_SIZE_LEGACY, totalCapacity*SLOT_SIZE);
			}else{
				partCapacities=new long[numPartitions];
				partOffsets=new long[numPartitions];
				for(int i=0; i<numPartitions; i++){
					partCapacities[i]=raf.readLong();
					partOffsets[i]=raf.readLong();
				}
				partSegments=new MappedByteBuffer[numPartitions][];
				for(int i=0; i<numPartitions; i++){
					partSegments[i]=mapSegments(channel, partOffsets[i], partCapacities[i]*SLOT_SIZE);
				}
			}

			partMasks=new long[numPartitions];
			for(int i=0; i<numPartitions; i++){
				partMasks[i]=partCapacities[i]-1;
			}
			partitionBits=Integer.numberOfTrailingZeros(Integer.highestOneBit(numPartitions));
			success=true;
		}finally{
			if(!success){closeQuietly(raf);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Lookup              ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Looks up a digitized accession key.
	 *
	 * @param key digitized accession key
	 * @return matching TaxID, or {@code -1} when absent
	 */
	public int get(long key){
		if(key<=0){return -1;}
		long h=hash(key);
		int part=(int)((h>>>(63-partitionBits))&(numPartitions-1));
		long mask=partMasks[part];
		MappedByteBuffer[] segs=partSegments[part];

		long slot=h&mask;

		for(long probe=0; probe<=mask; probe++){
			long idx=(slot+probe)&mask;
			long bytePos=idx*SLOT_SIZE;
			int segIdx=(int)(bytePos/SEGMENT_SIZE);
			int segOff=(int)(bytePos%SEGMENT_SIZE);

			long storedKey=segs[segIdx].getLong(segOff);
			if(storedKey==EMPTY){return -1;}
			if(storedKey==key){return segs[segIdx].getInt(segOff+8);}
		}
		return -1;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Hashing             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Mixes a digitized accession key into a nonnegative probe hash.
	 *
	 * @param key digitized accession key
	 * @return nonnegative hash value
	 */
	static long hash(long key){
		key^=(key>>>33);
		key*=0xff51afd7ed558ccdL;
		key^=(key>>>33);
		key*=0xc4ceb9fe1a85ec53L;
		key^=(key>>>33);
		return key&Long.MAX_VALUE;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Building             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds a partitioned disk accession table from digitized accession keys.
	 *
	 * @param outPath output table path
	 * @param keys digitized accession keys
	 * @param values TaxIDs parallel to {@code keys}
	 * @param count number of keys and values to read
	 * @throws IOException if the output table cannot be written
	 */
	public static void build(String outPath, long[] keys, int[] values, int count) throws IOException{
		final Timer t=new Timer();
		System.err.println("Building DiskAccessionTable: entries="+count+", file="+outPath);
		long added=0;
		try(Builder builder=new Builder(outPath, Math.max(1, count))){
			for(int i=0; i<count; i++){
				final long key=keys[i];
				if(key<=0){continue;}
				if(builder.add(key, values[i])){added++;}
			}
		}

		t.stop();
		final long fileSize=new File(outPath).length();
		System.err.println("DiskAccessionTable built: "+fileSize/(1024*1024)+"MB, entries="+added+" in "+t);
	}

	/**
	 * Starts an incremental table build sized for an expected entry count.
	 *
	 * @param outPath output table path
	 * @param expectedEntries expected number of positive keys
	 * @return builder for streaming entries into {@code outPath}
	 * @throws IOException if the output table cannot be created or mapped
	 */
	public static Builder build(String outPath, long expectedEntries) throws IOException{
		return new Builder(outPath, expectedEntries);
	}

	/**
	 * Inserts or updates one slot in a mapped partition.
	 *
	 * @param segs mapped segments for one partition
	 * @param mask partition capacity mask
	 * @param slot initial probe slot
	 * @param key digitized accession key
	 * @param value TaxID to store
	 * @return {@code true} if a new key was inserted
	 * @throws IOException if the partition has no empty slot
	 */
	private static boolean insert(MappedByteBuffer[] segs, long mask, long slot, long key, int value) throws IOException{
		for(long probe=0; probe<=mask; probe++){
			final long idx=(slot+probe)&mask;
			final long bytePos=idx*SLOT_SIZE;
			final int segIdx=(int)(bytePos/SEGMENT_SIZE);
			final int segOff=(int)(bytePos%SEGMENT_SIZE);
			final long storedKey=segs[segIdx].getLong(segOff);
			if(storedKey==EMPTY){
				segs[segIdx].putLong(segOff, key);
				segs[segIdx].putInt(segOff+8, value);
				return true;
			}
			if(storedKey==key){
				segs[segIdx].putInt(segOff+8, value);
				return false;
			}
		}
		throw new IOException("DiskAccessionTable partition is full; mask="+mask+" key="+key);
	}

	/**
	 * Chooses a power-of-two table capacity for the expected entry count.
	 *
	 * @param entries expected number of entries
	 * @return total slot count across all build partitions
	 */
	private static long tableCapacity(long entries){
		if(entries>(Long.MAX_VALUE/4)*3){throw new RuntimeException("Too many entries: "+entries);}
		long capacity=BUILD_PARTITIONS*MIN_PART_CAPACITY;
		final long target=entries+(entries/3)+(entries%3==0 ? 0 : 1);
		while(capacity<target){
			if(capacity>(Long.MAX_VALUE/2)){throw new RuntimeException("Too many entries: "+entries);}
			capacity*=2;
		}
		return capacity;
	}

	/**
	 * Calculates the partitioned-file header length.
	 *
	 * @param numPartitions number of partition descriptors
	 * @return bytes occupied by the prefix, partition table, and padding
	 */
	private static long headerSize(int numPartitions){
		return HEADER_SIZE_LEGACY+(long)numPartitions*16;
	}

	/** Closes a partially constructed table file during constructor failure cleanup.
	 * @param raf Open table file to close */
	private static void closeQuietly(RandomAccessFile raf){
		try{
			raf.close();
		}catch(IOException e){/* Failure is intentionally suppressed. */}
	}

	/**
	 * Incremental writer for a partitioned disk accession table. {@link BufferSet}
	 * instances are externally owned: callers must quiesce all adder threads and
	 * flush every buffer set before closing this builder.
	 */
	public static class Builder implements AutoCloseable{

		/**
		 * Creates and memory-maps a table file sized for the expected entry count.
		 *
		 * @param outPath output table path
		 * @param expectedEntries expected number of positive keys
		 * @throws IOException if the output table cannot be created or mapped
		 */
		public Builder(String outPath, long expectedEntries) throws IOException{
			final long expected=Math.max(1, expectedEntries);
			final long totalCapacity=tableCapacity(expected);
			numPartitions=BUILD_PARTITIONS;
			partitionBits=Integer.numberOfTrailingZeros(numPartitions);
			final long partCapacity=totalCapacity/numPartitions;
			partMask=partCapacity-1;
			final long partLength=partCapacity*SLOT_SIZE;
			final long headerSize=headerSize(numPartitions);
			final long fileLength=headerSize+totalCapacity*SLOT_SIZE;

			//TODO: Probable bug - length, header, or mapping failures can leave
			// the output handle open and a partially initialized table on disk.
			raf=new RandomAccessFile(outPath, "rw");
			raf.setLength(0);
			raf.setLength(fileLength);
			raf.seek(0);
			raf.writeLong(totalCapacity);
			raf.writeLong(0);
			raf.writeInt(numPartitions);
			raf.writeInt(0);
			for(int i=0; i<numPartitions; i++){
				raf.writeLong(partCapacity);
				raf.writeLong(headerSize+i*partLength);
			}
			final FileChannel channel=raf.getChannel();
			parts=new MappedByteBuffer[numPartitions][];
			for(int i=0; i<numPartitions; i++){
				parts[i]=mapSegments(channel, headerSize+i*partLength, partLength, FileChannel.MapMode.READ_WRITE);
			}
		}

		/**
		 * Adds or replaces one digitized accession mapping. Direct calls are
		 * serialized only by the caller; use one shared {@link BufferSet} for
		 * multiple producer threads.
		 *
		 * @param key digitized accession key
		 * @param value TaxID to store
		 * @return {@code true} if a new key was inserted
		 * @throws IOException if the target partition is full
		 */
		public boolean add(long key, int value) throws IOException{
			if(key<=0){return false;}
			final long h=hash(key);
			final int part=partition(h);
			final boolean added=insert(parts[part], partMask, h&partMask, key, value);
			if(added){addEntries(1);}
			return added;
		}

		/**
		 * Creates a per-partition buffered inserter. One shared {@link BufferSet}
		 * synchronizes each target partition's buffer; independent buffer sets do
		 * not protect their same-partition writes from each other.
		 * {@code entriesPerPartition} must be positive.
		 *
		 * @param entriesPerPartition maximum buffered entries per partition
		 * @return buffered inserter backed by this builder
		 */
		public BufferSet bufferSet(int entriesPerPartition){
			return new BufferSet(this, entriesPerPartition);
		}

		/**
		 * Selects a partition from the high bits of a mixed hash.
		 *
		 * @param h nonnegative hash value
		 * @return partition index
		 */
		private int partition(final long h){
			return (int)((h>>>(63-partitionBits))&(numPartitions-1));
		}

		/**
		 * Sorts one partition buffer by slot and inserts it into its partition.
		 * The slot sort is unstable, so conflicting duplicate keys in one buffered
		 * batch need not preserve arrival order. Parallel builds are deterministic
		 * only when duplicate accession keys agree on TaxID.
		 *
		 * @param part target partition index
		 * @param slots initial probe slots
		 * @param keys digitized accession keys
		 * @param values TaxIDs parallel to {@code keys}
		 * @param size number of buffered entries
		 * @return number of newly inserted keys
		 * @throws IOException if the target partition is full
		 */
		private long addSorted(int part, long[] slots, long[] keys, int[] values, int size) throws IOException{
			sort(slots, keys, values, size);
			long added=0;
			final MappedByteBuffer[] segs=parts[part];
			for(int i=0; i<size; i++){
				if(insert(segs, partMask, slots[i], keys[i], values[i])){added++;}
			}
			addEntries(added);
			return added;
		}

		/**
		 * Adds to the total entry count after a batch insert.
		 *
		 * @param added newly inserted entry count
		 */
		private synchronized void addEntries(long added){
			entryCount+=added;
		}

		/**
		 * Forces mapped partitions, writes the final entry count, and closes the
		 * output file. This does not flush any {@link BufferSet}; callers must
		 * quiesce and flush all buffer sets before closing.
		 * File closure is attempted even if forcing or header writing fails.
		 * A later close failure is suppressed behind the earlier failure.
		 */
		@Override
		public void close() throws IOException{
			if(closed){return;}
			closed=true;
			try(final RandomAccessFile r=raf){
				for(int i=0; i<numPartitions; i++){
					for(MappedByteBuffer mbb : parts[i]){mbb.force();}
				}
				r.seek(8);
				r.writeLong(entryCount);
			}
		}

		/** Output table file handle. */
		private final RandomAccessFile raf;
		/** Mapped write segments by partition. */
		private final MappedByteBuffer[][] parts;
		/** Number of build partitions. */
		private final int numPartitions;
		/** High hash bits used to choose a partition. */
		private final int partitionBits;
		/** Mask for probing within one partition. */
		private final long partMask;
		/** Number of newly inserted keys. */
		private long entryCount=0;
		/** True after {@link #close()} has been entered; suppresses retries. */
		private boolean closed=false;

	}

	/**
	 * Buffered inserter that groups entries by target partition. Share one
	 * instance across producer threads, then flush it after those threads finish.
	 */
	public static class BufferSet{

		/**
		 * Allocates one insertion buffer per builder partition.
		 * {@code entriesPerPartition} must be positive.
		 *
		 * @param builder_ destination builder
		 * @param entriesPerPartition maximum buffered entries per partition
		 */
		private BufferSet(Builder builder_, int entriesPerPartition){
			builder=builder_;
			buffers=new Buffer[builder.numPartitions];
			for(int i=0; i<buffers.length; i++){
				buffers[i]=new Buffer(builder, i, entriesPerPartition);
			}
		}

		/**
		 * Adds one mapping to its target partition buffer. Calls on one
		 * {@link BufferSet} are safe because each partition's {@link Buffer} is
		 * synchronized.
		 *
		 * @param key digitized accession key
		 * @param value TaxID to store
		 * @throws IOException if flushing the target buffer fails
		 */
		public void add(long key, int value) throws IOException{
			if(key<=0){return;}
			buffers[builder.partition(hash(key))].add(key, value);
		}

		/**
		 * Flushes all buffered mappings into the destination table. Call this
		 * after all producer threads using this buffer set have quiesced.
		 *
		 * @throws IOException if any target partition is full
		 */
		public void flush() throws IOException{
			for(int i=0; i<buffers.length; i++){buffers[i].flush();}
		}

		/**
		 * Counts the new keys inserted by all buffers. Call this after all producer
		 * threads have quiesced and {@link #flush()} has returned.
		 *
		 * @return total newly inserted keys
		 */
		public long added(){
			long sum=0;
			for(int i=0; i<buffers.length; i++){sum+=buffers[i].added;}
			return sum;
		}

		/** Destination builder. */
		private final Builder builder;
		/** Per-partition insertion buffers. */
		private final Buffer[] buffers;

	}

	/**
	 * Accumulates keys and values for one partition before sorted insertion.
	 */
	private static class Buffer{

		/**
		 * Allocates a buffer for one target partition.
		 *
		 * @param builder_ destination builder
		 * @param part_ target partition index
		 * @param capacity maximum buffered entries
		 */
		Buffer(Builder builder_, int part_, int capacity){
			builder=builder_;
			part=part_;
			keys=new LongList(capacity);
			values=new IntList(capacity);
			sortSlots=new long[capacity];
			maxSize=capacity;
		}

		/**
		 * Adds a mapping and flushes when the buffer reaches capacity.
		 *
		 * @param key digitized accession key
		 * @param value TaxID to store
		 * @throws IOException if the flush fails
		 */
		synchronized void add(long key, int value) throws IOException{
			keys.add(key);
			values.add(value);
			if(keys.size>=maxSize){flushInner();}
		}

		/**
		 * Flushes pending mappings into the destination table.
		 *
		 * @throws IOException if the target partition is full
		 */
		synchronized void flush() throws IOException{
			flushInner();
		}

		/**
		 * Inserts the currently buffered mappings and clears the lists.
		 *
		 * @throws IOException if the target partition is full
		 */
		private void flushInner() throws IOException{
			final int size=keys.size;
			if(size>0){
				final long[] ka=keys.array;
				final int[] va=values.array;
				final long mask=builder.partMask;
				for(int i=0; i<size; i++){sortSlots[i]=hash(ka[i])&mask;}
				added+=builder.addSorted(part, sortSlots, ka, va, size);
				keys.size=0;
				values.size=0;
			}
		}

		/** Destination builder. */
		private final Builder builder;
		/** Target partition index. */
		private final int part;
		/** Buffered digitized accession keys. */
		private final LongList keys;
		/** Buffered TaxIDs. */
		private final IntList values;
		/** Initial probe slots parallel to {@link #keys}. */
		private final long[] sortSlots;
		/** Flush threshold for this buffer. */
		private final int maxSize;
		/** Number of new keys inserted by this buffer. */
		private long added=0;

	}

	/**
	 * Sorts buffered entries by initial probe slot.
	 *
	 * @param slots initial probe slots
	 * @param keys digitized accession keys parallel to {@code slots}
	 * @param values TaxIDs parallel to {@code slots}
	 * @param size number of elements to sort
	 */
	private static void sort(long[] slots, long[] keys, int[] values, int size){
		if(size>1){quickSort(slots, keys, values, 0, size-1);}
	}

	/**
	 * Sorts a slot range with quicksort and an insertion-sort cutoff.
	 *
	 * @param slots initial probe slots
	 * @param keys digitized accession keys parallel to {@code slots}
	 * @param values TaxIDs parallel to {@code slots}
	 * @param left inclusive left index
	 * @param right inclusive right index
	 */
	private static void quickSort(long[] slots, long[] keys, int[] values, int left, int right){
		while(right-left>SORT_SMALL){
			final long pivot=slots[(left+right)>>>1];
			int a=left, b=right;
			while(a<=b){
				while(slots[a]<pivot){a++;}
				while(slots[b]>pivot){b--;}
				if(a<=b){
					swap(slots, keys, values, a, b);
					a++;
					b--;
				}
			}
			if(b-left<right-a){
				if(left<b){quickSort(slots, keys, values, left, b);}
				left=a;
			}else{
				if(a<right){quickSort(slots, keys, values, a, right);}
				right=b;
			}
		}
		insertionSort(slots, keys, values, left, right);
	}

	/**
	 * Insertion-sorts a small inclusive slot range.
	 *
	 * @param slots initial probe slots
	 * @param keys digitized accession keys parallel to {@code slots}
	 * @param values TaxIDs parallel to {@code slots}
	 * @param left inclusive left index
	 * @param right inclusive right index
	 */
	private static void insertionSort(long[] slots, long[] keys, int[] values, int left, int right){
		for(int i=left+1; i<=right; i++){
			final long slot=slots[i], key=keys[i];
			final int value=values[i];
			int j=i;
			while(j>left && slots[j-1]>slot){
				slots[j]=slots[j-1];
				keys[j]=keys[j-1];
				values[j]=values[j-1];
				j--;
			}
			slots[j]=slot;
			keys[j]=key;
			values[j]=value;
		}
	}

	/**
	 * Swaps two parallel slot/key/value entries.
	 *
	 * @param slots initial probe slots
	 * @param keys digitized accession keys parallel to {@code slots}
	 * @param values TaxIDs parallel to {@code slots}
	 * @param a first index
	 * @param b second index
	 */
	private static void swap(long[] slots, long[] keys, int[] values, int a, int b){
		final long slot=slots[a], key=keys[a];
		final int value=values[a];
		slots[a]=slots[b];
		keys[a]=keys[b];
		values[a]=values[b];
		slots[b]=slot;
		keys[b]=key;
		values[b]=value;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Mapping              ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Maps a table region read-only.
	 *
	 * @param channel open table channel
	 * @param offset first byte to map
	 * @param length number of bytes to map
	 * @return mapped segments covering the requested region
	 * @throws IOException if a segment cannot be mapped
	 */
	private static MappedByteBuffer[] mapSegments(FileChannel channel, long offset, long length) throws IOException{
		return mapSegments(channel, offset, length, FileChannel.MapMode.READ_ONLY);
	}

	/**
	 * Maps a table region using chunks small enough for buffer indexes.
	 *
	 * @param channel open table channel
	 * @param offset first byte to map
	 * @param length number of bytes to map
	 * @param mode requested mapping mode
	 * @return mapped segments covering the requested region
	 * @throws IOException if a segment cannot be mapped
	 */
	private static MappedByteBuffer[] mapSegments(FileChannel channel, long offset, long length, FileChannel.MapMode mode) throws IOException{
		int numSegments=(int)((length+SEGMENT_SIZE-1)/SEGMENT_SIZE);
		MappedByteBuffer[] segs=new MappedByteBuffer[numSegments];
		for(int i=0; i<numSegments; i++){
			long segStart=offset+((long)i*SEGMENT_SIZE);
			long segLen=Math.min(SEGMENT_SIZE, offset+length-segStart);
			segs[i]=channel.map(mode, segStart, segLen);
		}
		return segs;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Closing              ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Closes the table channel and file.
	 *
	 * @throws IOException if either close fails
	 */
	public void close() throws IOException{
		channel.close();
		raf.close();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Fields              ----------------*/
	/*--------------------------------------------------------------*/

	/** Backing table file. */
	private final File file;
	/** Read-only random-access handle for the table. */
	private final RandomAccessFile raf;
	/** Read-only channel used for mapping table segments. */
	private final FileChannel channel;
	/** Mapped read segments by partition. */
	private final MappedByteBuffer[][] partSegments;
	/** Slot capacity of each partition. */
	private final long[] partCapacities;
	/** Slot masks for wrapping probes within each partition. */
	private final long[] partMasks;
	/** File offsets for partition bodies. */
	private final long[] partOffsets;
	/** Number of entries recorded in the table header. */
	private final long entryCount;
	/** Number of mapped table partitions. */
	private int numPartitions;
	/** High hash bits used to choose a lookup partition. */
	private final int partitionBits;

	/** Bytes in the fixed legacy header. */
	private static final int HEADER_SIZE_LEGACY=32;
	/** Bytes in one key/value slot. */
	private static final long SLOT_SIZE=12;
	/** Empty-slot sentinel key. */
	private static final long EMPTY=0;
	/** Maximum mapped segment length, aligned to whole slots. */
	private static final long SEGMENT_SIZE=((long)Integer.MAX_VALUE/SLOT_SIZE)*SLOT_SIZE;
	/** Number of partitions written by new builds. */
	static final int BUILD_PARTITIONS=1<<14;
	/** Minimum slots per build partition. */
	private static final int MIN_PART_CAPACITY=16;
	/** Quicksort range size handled by insertion sort. */
	private static final int SORT_SMALL=32;

}
