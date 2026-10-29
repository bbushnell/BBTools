package tax;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.nio.file.*;
import java.util.concurrent.*;

/**
 * Async disk accessor — uses AsynchronousFileChannel to submit all
 * first-probe reads at once. Sequential collision probes reuse the same channel.
 *
 * @author Brian Bushnell, Chloe
 */
public class AsyncDiskAccessor extends DiskAccessor{

	/** Opens a disk accession table and reads its partition header.
	 * The asynchronous channel remains open for lookup batches; the temporary
	 * {@link RandomAccessFile} used for header loading is closed before return.
	 *
	 * @param path Path to the binary accession table
	 * @throws IOException if the table or its header cannot be opened or read
	 */
	public AsyncDiskAccessor(String path) throws IOException{
		afc=AsynchronousFileChannel.open(Paths.get(path), StandardOpenOption.READ);
		boolean success=false;
		RandomAccessFile raf=null;
		try{

			// Read header synchronously via a temporary channel
			raf=new RandomAccessFile(path, "r");
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
			success=true;
		}finally{
			if(!success){
				if(raf!=null){closeQuietly(raf);}
				closeQuietly(afc);
			}
		}
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
				//TODO: Probable bug - asynchronous read failures are caught and
				// converted to absent results in both first and collision probes.
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
	public void close() throws IOException{
		afc.close();
	}

	/** Closes a partially constructed handle during constructor failure cleanup.
	 *
	 * @param c Handle to close
	 */
	private static void closeQuietly(final Closeable c){
		try{
			c.close();
		}catch(IOException e){/* Failure is intentionally suppressed. */}
	}

	/** Asynchronous channel shared by all submitted slot reads */
	private final AsynchronousFileChannel afc;
	/** Slot count for each table partition, parsed from the file header */
	private final long[] partCapacities;
	/** Byte offset for the first slot of each partition */
	private final long[] partOffsets;
	/** Probe-mask for wrapping within each power-of-two partition */
	private final long[] partMasks;
	/** Number of accession entries recorded in the table header */
	private final long entryCount;
	/** Number of table partitions recorded in the table header */
	private final int numPartitions;
	/** High hash bits consumed to choose a partition */
	private final int partitionBits;
}
