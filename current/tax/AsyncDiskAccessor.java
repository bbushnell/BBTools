package tax;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.nio.file.*;
import java.util.concurrent.*;

/**
 * Async disk accessor — uses AsynchronousFileChannel to submit all
 * first-probe reads at once. Subsequent probes use the same channel
 * but are nearly always page-cache hits (12B slots, 4KB pages).
 * @author Brian Bushnell, Chloe
 */
public class AsyncDiskAccessor extends DiskAccessor {

	public AsyncDiskAccessor(String path) throws IOException {
		afc=AsynchronousFileChannel.open(Paths.get(path), StandardOpenOption.READ);

		// Read header synchronously via a temporary channel
		RandomAccessFile raf=new RandomAccessFile(path, "r");
		long totalCapacity=raf.readLong();
		entryCount=raf.readLong();
		numPartitions=raf.readInt();
		raf.readInt(); //reserved

		partCapacities=new long[numPartitions];
		partOffsets=new long[numPartitions];
		partMasks=new long[numPartitions];
		for(int i=0; i<numPartitions; i++){
			partCapacities[i]=raf.readLong();
			partOffsets[i]=raf.readLong();
			partMasks[i]=partCapacities[i]-1;
		}
		partitionBits=Integer.numberOfTrailingZeros(Integer.highestOneBit(numPartitions));
		raf.close();
	}

	@Override
	public int[] getBatch(long[] keys){
		final int n=keys.length;
		int[] results=new int[n];

		// Phase 1: submit all first-probe reads
		ByteBuffer[] bufs=new ByteBuffer[n];
		@SuppressWarnings("unchecked")
		Future<Integer>[] futures=(Future<Integer>[])new Future[n];
		long[] probePositions=new long[n];
		int[] parts=new int[n];
		long[] masks=new long[n];

		for(int i=0; i<n; i++){
			if(keys[i]<=0){results[i]=-1; continue;}
			long h=hash(keys[i]);
			int part=(int)((h>>>(63-partitionBits))&(numPartitions-1));
			parts[i]=part;
			masks[i]=partMasks[part];
			long slot=h&masks[i];
			long bytePos=partOffsets[part]+slot*SLOT_SIZE;
			probePositions[i]=bytePos;
			bufs[i]=ByteBuffer.allocate(12);
			futures[i]=afc.read(bufs[i], bytePos);
		}

		// Phase 2: collect results, probe further if needed
		for(int i=0; i<n; i++){
			if(futures[i]==null){continue;}
			try{
				futures[i].get();
			}catch(Exception e){results[i]=-1; continue;}
			bufs[i].flip();
			long storedKey=bufs[i].getLong();
			if(storedKey==EMPTY){results[i]=-1; continue;}
			if(storedKey==keys[i]){results[i]=bufs[i].getInt(); continue;}

			// Collision — sequential probes (same page, fast)
			long slot=(probePositions[i]-partOffsets[parts[i]])/SLOT_SIZE;
			long mask=masks[i];
			results[i]=-1;
			for(long probe=1; probe<=mask; probe++){
				long idx=(slot+probe)&mask;
				long bytePos=partOffsets[parts[i]]+idx*SLOT_SIZE;
				ByteBuffer buf=ByteBuffer.allocate(12);
				try{
					afc.read(buf, bytePos).get();
				}catch(Exception e){break;}
				buf.flip();
				storedKey=buf.getLong();
				if(storedKey==EMPTY){break;}
				if(storedKey==keys[i]){results[i]=buf.getInt(); break;}
			}
		}
		return results;
	}

	@Override public String name(){return "Async (AsynchronousFileChannel)";}

	@Override
	public void close() throws IOException {
		afc.close();
	}

	private final AsynchronousFileChannel afc;
	private final long[] partCapacities;
	private final long[] partOffsets;
	private final long[] partMasks;
	private final long entryCount;
	private final int numPartitions;
	private final int partitionBits;
}
