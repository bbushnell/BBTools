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
	/*----------------        Mapping              ----------------*/
	/*--------------------------------------------------------------*/

	private static MappedByteBuffer[] mapSegments(FileChannel channel, long offset, long length) throws IOException {
		int numSegments=(int)((length+SEGMENT_SIZE-1)/SEGMENT_SIZE);
		MappedByteBuffer[] segs=new MappedByteBuffer[numSegments];
		for(int i=0; i<numSegments; i++){
			long segStart=offset+((long)i*SEGMENT_SIZE);
			long segLen=Math.min(SEGMENT_SIZE, offset+length-segStart);
			segs[i]=channel.map(FileChannel.MapMode.READ_ONLY, segStart, segLen);
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

}
