package tax;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * NIO disk accessor — uses FileChannel positioned reads.
 * Thread count is configurable; 1 thread = pure serial FileChannel.
 * @author Brian Bushnell, Chloe
 */
public class NIODiskAccessor extends DiskAccessor {

	public NIODiskAccessor(String path, int maxThreads) throws IOException {
		this.maxThreads=maxThreads;
		raf=new RandomAccessFile(path, "r");
		channel=raf.getChannel();

		ByteBuffer hdr=ByteBuffer.allocate(32);
		channel.read(hdr, 0);
		hdr.flip();
		long totalCapacity=hdr.getLong();
		entryCount=hdr.getLong();
		numPartitions=hdr.getInt();
		hdr.getInt();

		partCapacities=new long[numPartitions];
		partOffsets=new long[numPartitions];
		partMasks=new long[numPartitions];
		ByteBuffer ptable=ByteBuffer.allocate(numPartitions*16);
		channel.read(ptable, 24);
		ptable.flip();
		for(int i=0; i<numPartitions; i++){
			partCapacities[i]=ptable.getLong();
			partOffsets[i]=ptable.getLong();
			partMasks[i]=partCapacities[i]-1;
		}
		partitionBits=Integer.numberOfTrailingZeros(Integer.highestOneBit(numPartitions));
	}

	@Override
	public int[] getBatch(long[] keys){
		final int n=keys.length;
		int[] results=new int[n];
		int threads=Math.min(maxThreads, Math.max(1, n/MIN_PER_THREAD));
		if(threads<=1){
			for(int i=0; i<n; i++){results[i]=getOne(keys[i]);}
			return results;
		}
		int perThread=(n+threads-1)/threads;
		Thread[] workers=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int start=t*perThread;
			final int end=Math.min(start+perThread, n);
			workers[t]=new Thread(()->{
				for(int i=start; i<end; i++){
					results[i]=getOne(keys[i]);
				}
			});
			workers[t].start();
		}
		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){}
		}
		return results;
	}

	private int getOne(long key){
		if(key<=0){return -1;}
		long h=hash(key);
		int part=(int)((h>>>(63-partitionBits))&(numPartitions-1));
		long mask=partMasks[part];
		long baseOffset=partOffsets[part];
		long slot=h&mask;

		ByteBuffer buf=ByteBuffer.allocate(12);
		for(long probe=0; probe<=mask; probe++){
			long idx=(slot+probe)&mask;
			long bytePos=baseOffset+idx*SLOT_SIZE;
			buf.clear();
			try{
				int bytesRead=0;
				while(bytesRead<12){
					int r=channel.read(buf, bytePos+bytesRead);
					if(r<0){return -1;}
					bytesRead+=r;
				}
			}catch(IOException e){return -1;}
			buf.flip();
			long storedKey=buf.getLong();
			if(storedKey==EMPTY){return -1;}
			if(storedKey==key){return buf.getInt();}
		}
		return -1;
	}

	@Override public String name(){return "NIO (FileChannel, "+maxThreads+" thread"+(maxThreads>1?"s":"")+")";}

	@Override
	public void close() throws IOException {
		channel.close();
		raf.close();
	}

	private final RandomAccessFile raf;
	private final FileChannel channel;
	private final long[] partCapacities;
	private final long[] partOffsets;
	private final long[] partMasks;
	private final long entryCount;
	private final int numPartitions;
	private final int partitionBits;
	private final int maxThreads;

	private static final int MIN_PER_THREAD=50;
}
