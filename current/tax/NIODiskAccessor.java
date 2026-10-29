package tax;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * Disk accessor using positioned {@link FileChannel} reads.
 * A single accessor may split a large batch across worker threads; positioned
 * reads let those workers share the channel while reading independent slots.
 * @author Brian Bushnell, Chloe
 */
public class NIODiskAccessor extends DiskAccessor{

	/**
	 * Opens an immutable disk table and loads its partition metadata.
	 * @param path Table file built by DiskAccessionTable
	 * @param maxThreads Maximum worker threads per batch
	 * @throws IOException If the table cannot be opened or its header cannot be read
	 */
	public NIODiskAccessor(final String path, final int maxThreads) throws IOException{
		this.maxThreads=maxThreads;
		raf=new RandomAccessFile(path, "r");
		channel=raf.getChannel();
		boolean success=false;
		try{
			final ByteBuffer hdr=ByteBuffer.allocate(32);
			channel.read(hdr, 0);
			hdr.flip();
			long totalCapacity=hdr.getLong();
			entryCount=hdr.getLong();
			numPartitions=hdr.getInt();
			hdr.getInt();

			partCapacities=new long[numPartitions];
			partOffsets=new long[numPartitions];
			partMasks=new long[numPartitions];
			final ByteBuffer ptable=ByteBuffer.allocate(numPartitions*16);
			channel.read(ptable, 24);
			ptable.flip();
			for(int i=0; i<numPartitions; i++){
				partCapacities[i]=ptable.getLong();
				partOffsets[i]=ptable.getLong();
				partMasks[i]=partCapacities[i]-1;
			}
			partitionBits=Integer.numberOfTrailingZeros(Integer.highestOneBit(numPartitions));
			success=true;
		}finally{
			if(!success){closeQuietly(raf);}
		}
	}

	/** Splits large batches over worker threads and writes results by input index.
	 * Joins every worker before returning, restores caller interruption after
	 * joining, and propagates any worker failure on the caller thread. */
	@Override
	public int[] getBatch(final long[] keys){
		final int n=keys.length;
		final int[] results=new int[n];
		final int threads=Math.min(maxThreads, Math.max(1, n/MIN_PER_THREAD));
		if(threads<=1){
			try{
				for(int i=0; i<n; i++){results[i]=getOne(keys[i]);}
			}catch(IOException e){
				throw new RuntimeException("NIO disk slot read failed.", e);
			}
			return results;
		}
		final int perThread=(n+threads-1)/threads;
		final LookupTask[] tasks=new LookupTask[threads];
		final Thread[] workers=new Thread[threads];
		Throwable startupFailure=null;
		try{
			for(int t=0; t<threads; t++){
				final int start=t*perThread;
				final int end=Math.min(start+perThread, n);
				final LookupTask task=new LookupTask(keys, results, start, end);
				tasks[t]=task;
				final Thread worker=makeWorker(task);
				workers[t]=worker;
				worker.start();
			}
		}catch(Throwable e){
			startupFailure=e;
		}
		boolean interrupted=false;
		for(Thread w : workers){
			if(w==null || w.getState()==Thread.State.NEW){continue;}
			boolean joined=false;
			while(!joined){
				try{
					w.join();
					joined=true;
				}catch(InterruptedException e){interrupted=true;}
			}
		}
		final Throwable workerFailure=firstWorkerFailure(tasks);
		if(interrupted){Thread.currentThread().interrupt();}
		if(startupFailure!=null){
			suppressDistinct(startupFailure, workerFailure);
			rethrowStartupFailure(startupFailure);
		}
		if(workerFailure!=null){
			rethrowWorkerFailure(workerFailure);
		}
		if(interrupted){
			throw new RuntimeException("Interrupted while waiting for NIO disk lookup workers.");
		}
		return results;
	}

	/** Creates one worker thread for a slice of a positioned-read batch.
	 * @param task Slice reader to run in the returned thread
	 * @return New unstarted worker thread */
	protected Thread makeWorker(final Runnable task){
		return new Thread(task);
	}

	/** Selects the first worker failure in slice order after workers join.
	 * @param tasks Lookup tasks by slice order
	 * @return First worker failure, or null if all started workers completed */
	private static Throwable firstWorkerFailure(final LookupTask[] tasks){
		for(LookupTask task : tasks){
			if(task!=null && task.failure!=null){return task.failure;}
		}
		return null;
	}

	/** Retains a distinct secondary failure on a primary failure.
	 * @param primary Primary failure to retain
	 * @param secondary Secondary failure to suppress, if present and distinct */
	private static void suppressDistinct(final Throwable primary, final Throwable secondary){
		if(primary!=null && secondary!=null && primary!=secondary){
			primary.addSuppressed(secondary);
		}
	}

	/** Lookup slice with join-visible per-worker failure storage. */
	private final class LookupTask implements Runnable{

		/** Stores the slice bounds and shared batch arrays.
		 * @param keys_ Query keys shared across workers
		 * @param results_ Result array partitioned by worker
		 * @param start_ Inclusive first result index
		 * @param end_ Exclusive last result index */
		LookupTask(final long[] keys_, final int[] results_, final int start_,
				final int end_){
			keys=keys_;
			results=results_;
			start=start_;
			end=end_;
		}

		/** Reads one slice and retains any failure for the caller after join. */
		@Override
		public void run(){
			try{
				for(int i=start; i<end; i++){results[i]=getOne(keys[i]);}
			}catch(Throwable e){failure=e;}
		}

		/** Query keys shared across all workers. */
		private final long[] keys;
		/** Result array partitioned by worker-owned index ranges. */
		private final int[] results;
		/** Inclusive first index in the worker-owned slice. */
		private final int start;
		/** Exclusive last index in the worker-owned slice. */
		private final int end;
		/** First failure from this slice; deliberately nonvolatile because the
		 * caller reads it only after a completed Thread.join(). */
		private Throwable failure;
	}

	/** Finds one key by linear probing inside its selected partition.
	 * @param key Digitized accession key
	 * @return TaxID for the key, or -1 on absence
	 * @throws IOException If a matching slot cannot be read fully */
	private int getOne(final long key) throws IOException{
		if(key<=0){return -1;}
		final long h=hash(key);
		final int part=(int)((h>>>(63-partitionBits))&(numPartitions-1));
		final long mask=partMasks[part];
		final long baseOffset=partOffsets[part];
		final long slot=h&mask;

		final ByteBuffer buf=ByteBuffer.allocate(12);
		for(long probe=0; probe<=mask; probe++){
			final long idx=(slot+probe)&mask;
			final long bytePos=baseOffset+idx*SLOT_SIZE;
			buf.clear();
			int bytesRead=0;
			while(bytesRead<12){
				final int r=channel.read(buf, bytePos+bytesRead);
				if(r<0){
					throw new EOFException("Unexpected EOF reading NIO disk slot at "+
							bytePos+" for key "+key+" after "+bytesRead+" bytes.");
				}
				bytesRead+=r;
			}
			buf.flip();
			final long storedKey=buf.getLong();
			if(storedKey==EMPTY){return -1;}
			if(storedKey==key){return buf.getInt();}
		}
		return -1;
	}

	/** Throws a worker failure on the calling thread after all workers join.
	 * @param failure First worker failure, or null if all workers completed */
	private static void rethrowWorkerFailure(final Throwable failure){
		if(failure==null){return;}
		if(failure instanceof Error){throw (Error)failure;}
		throw new RuntimeException("NIO disk slot read failed.", failure);
	}

	/** Throws an exact startup failure on the calling thread.
	 * @param failure Failure thrown while creating or starting a worker */
	private static void rethrowStartupFailure(final Throwable failure){
		if(failure instanceof Error){throw (Error)failure;}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		throw new RuntimeException("NIO disk worker startup failed.", failure);
	}

	/** Closes a partially constructed table file during constructor failure cleanup.
	 * @param raf Open table file to close */
	private static void closeQuietly(final RandomAccessFile raf){
		try{
			raf.close();
		}catch(IOException e){/* Failure is intentionally suppressed. */}
	}

	@Override public String name(){return "NIO (FileChannel, "+maxThreads+" thread"+(maxThreads>1?"s":"")+")";}

	/** Closes the FileChannel and its RandomAccessFile. */
	@Override
	public void close() throws IOException{
		channel.close();
		raf.close();
	}

	/** Owning file handle for the disk table. */
	private final RandomAccessFile raf;
	/** Channel used for positioned slot reads by one or more threads. */
	private final FileChannel channel;
	/** Number of slots per partition. */
	private final long[] partCapacities;
	/** Byte offsets of partition starts. */
	private final long[] partOffsets;
	/** Capacity-1 masks for partition-local linear probing. */
	private final long[] partMasks;
	/** Number of occupied slots recorded in the table header. */
	private final long entryCount;
	/** Number of hash partitions in the table. */
	private final int numPartitions;
	/** Log2 partition count used to extract high hash bits. */
	private final int partitionBits;
	/** Maximum worker count for one batch. */
	private final int maxThreads;

	/** Minimum query count per worker before enabling another worker. */
	private static final int MIN_PER_THREAD=50;
}
