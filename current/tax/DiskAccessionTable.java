package tax;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import shared.Timer;

/**
 * Disk-backed accession-to-TaxID table using a memory-mapped
 * open-addressing hash table, partitioned for build efficiency.
 *
 * Each slot is 12 bytes: [8-byte long key][4-byte int value].
 * Empty slots have key=0 (EMPTY sentinel).
 *
 * Partitioned file format:
 *   Header (32 bytes):
 *     [long totalCapacity] [long entryCount] [int numPartitions] [int reserved]
 *   Partition table (numPartitions x 16 bytes):
 *     [long partCapacity] [long partOffset]
 *   Body:
 *     numPartitions x (partCapacity x 12-byte slots)
 *
 * Lookup: hash the key, select partition, probe within that partition.
 *
 * @author Brian Bushnell, Chloe
 */
public class DiskAccessionTable {

	/*--------------------------------------------------------------*/
	/*----------------        Construction         ----------------*/
	/*--------------------------------------------------------------*/

	public DiskAccessionTable(String path) throws IOException {
		file=new File(path);
		if(!file.exists()){throw new IOException("DiskAccessionTable file not found: "+path);}

		raf=new RandomAccessFile(file, "r");
		channel=raf.getChannel();

		long totalCapacity=raf.readLong();
		entryCount=raf.readLong();
		numPartitions=raf.readInt();
		raf.readInt(); //reserved

		if(numPartitions<1 || numPartitions>64){
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
		for(int i=0; i<numPartitions; i++){partMasks[i]=partCapacities[i]-1;}
		partitionBits=Integer.numberOfTrailingZeros(Integer.highestOneBit(numPartitions));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Lookup              ----------------*/
	/*--------------------------------------------------------------*/

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
	 */
	public static void build(String outPath, long[] keys, int[] values, int count) throws IOException {
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

	public static Builder build(String outPath, long expectedEntries) throws IOException {
		return new Builder(outPath, expectedEntries);
	}

	private static boolean insert(MappedByteBuffer[] segs, long mask, long slot, long key, int value) throws IOException {
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

	private static long headerSize(int numPartitions){
		return HEADER_SIZE_LEGACY+(long)numPartitions*16;
	}

	public static class Builder implements AutoCloseable {

		public Builder(String outPath, long expectedEntries) throws IOException {
			final long expected=Math.max(1, expectedEntries);
			final long totalCapacity=tableCapacity(expected);
			numPartitions=BUILD_PARTITIONS;
			partitionBits=Integer.numberOfTrailingZeros(numPartitions);
			final long partCapacity=totalCapacity/numPartitions;
			partMask=partCapacity-1;
			final long partLength=partCapacity*SLOT_SIZE;
			final long headerSize=headerSize(numPartitions);
			final long fileLength=headerSize+totalCapacity*SLOT_SIZE;

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
			raf.writeLong(0);//Pad the body to byte 160 for 8 partitions, matching existing tables.

			final FileChannel channel=raf.getChannel();
			parts=new MappedByteBuffer[numPartitions][];
			for(int i=0; i<numPartitions; i++){
				parts[i]=mapSegments(channel, headerSize+i*partLength, partLength, FileChannel.MapMode.READ_WRITE);
			}
		}

		public boolean add(long key, int value) throws IOException {
			if(key<=0){return false;}
			final long h=hash(key);
			final int part=partition(h);
			final boolean added=insert(parts[part], partMask, h&partMask, key, value);
			if(added){addEntries(1);}
			return added;
		}

		public BufferSet bufferSet(int entriesPerPartition){
			return new BufferSet(this, entriesPerPartition);
		}

		private int partition(final long h){
			return (int)((h>>>(63-partitionBits))&(numPartitions-1));
		}

		private long addSorted(int part, long[] slots, long[] keys, int[] values, int size) throws IOException {
			sort(slots, keys, values, size);
			long added=0;
			final MappedByteBuffer[] segs=parts[part];
			for(int i=0; i<size; i++){
				if(insert(segs, partMask, slots[i], keys[i], values[i])){added++;}
			}
			addEntries(added);
			return added;
		}

		private synchronized void addEntries(long added){
			entryCount+=added;
		}

		@Override
		public void close() throws IOException {
			if(closed){return;}
			closed=true;
			for(int i=0; i<numPartitions; i++){
				for(MappedByteBuffer mbb : parts[i]){mbb.force();}
			}
			raf.seek(8);
			raf.writeLong(entryCount);
			raf.close();
		}

		private final RandomAccessFile raf;
		private final MappedByteBuffer[][] parts;
		private final int numPartitions;
		private final int partitionBits;
		private final long partMask;
		private long entryCount=0;
		private boolean closed=false;

	}

	public static class BufferSet {

		private BufferSet(Builder builder_, int entriesPerPartition){
			builder=builder_;
			buffers=new Buffer[builder.numPartitions];
			for(int i=0; i<buffers.length; i++){
				buffers[i]=new Buffer(builder, i, entriesPerPartition);
			}
		}

		public void add(long key, int value) throws IOException {
			if(key<=0){return;}
			final long h=hash(key);
			buffers[builder.partition(h)].add(h&builder.partMask, key, value);
		}

		public void flush() throws IOException {
			for(int i=0; i<buffers.length; i++){buffers[i].flush();}
		}

		public long added(){
			long sum=0;
			for(int i=0; i<buffers.length; i++){sum+=buffers[i].added;}
			return sum;
		}

		private final Builder builder;
		private final Buffer[] buffers;

	}

	private static class Buffer {

		Buffer(Builder builder_, int part_, int capacity){
			builder=builder_;
			part=part_;
			slots=new long[capacity];
			keys=new long[capacity];
			values=new int[capacity];
		}

		synchronized void add(long slot, long key, int value) throws IOException {
			slots[size]=slot;
			keys[size]=key;
			values[size]=value;
			size++;
			if(size>=slots.length){flushInner();}
		}

		synchronized void flush() throws IOException {
			flushInner();
		}

		private void flushInner() throws IOException {
			if(size>0){
				added+=builder.addSorted(part, slots, keys, values, size);
				size=0;
			}
		}

		private final Builder builder;
		private final int part;
		private final long[] slots;
		private final long[] keys;
		private final int[] values;
		private int size=0;
		private long added=0;

	}

	private static void sort(long[] slots, long[] keys, int[] values, int size){
		if(size>1){quickSort(slots, keys, values, 0, size-1);}
	}

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

	private static MappedByteBuffer[] mapSegments(FileChannel channel, long offset, long length) throws IOException {
		return mapSegments(channel, offset, length, FileChannel.MapMode.READ_ONLY);
	}

	private static MappedByteBuffer[] mapSegments(FileChannel channel, long offset, long length, FileChannel.MapMode mode) throws IOException {
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

	public void close() throws IOException {
		channel.close();
		raf.close();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Fields              ----------------*/
	/*--------------------------------------------------------------*/

	private final File file;
	private final RandomAccessFile raf;
	private final FileChannel channel;
	private final MappedByteBuffer[][] partSegments;
	private final long[] partCapacities;
	private final long[] partMasks;
	private final long[] partOffsets;
	private final long entryCount;
	private int numPartitions;
	private final int partitionBits;

	private static final int HEADER_SIZE_LEGACY=32;
	private static final long SLOT_SIZE=12;
	private static final long EMPTY=0;
	private static final long SEGMENT_SIZE=((long)Integer.MAX_VALUE/SLOT_SIZE)*SLOT_SIZE;
	static final int BUILD_PARTITIONS=8;
	private static final int MIN_PART_CAPACITY=16;
	private static final int SORT_SMALL=32;

}
