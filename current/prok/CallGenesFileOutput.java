package prok;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;
import stream.Read;
import structures.ByteBuilder;

/** Bounded output handoff for one file, consumed in input-file order.
 * The file worker transfers ownership of every submitted buffer/read list.
 * @author Brian, Raiden */
public final class CallGenesFileOutput {

	/** All active files share one failure signal so a failed producer cannot strand publication. */
	public static final class Failure{
		public void fail(Throwable cause){
			assert(cause!=null) : "An output failure must retain the originating worker exception";
			error=cause;
		}
		void check(){if(error!=null){throw new IllegalStateException("Multi-file caller failed", error);}}
		private volatile Throwable error;
	}
	public CallGenesFileOutput(Failure failure_){
		this(failure_, MIN_CHUNKS);
	}
	public CallGenesFileOutput(Failure failure_, int chunks){
		if(failure_==null){throw new IllegalArgumentException("Shared failure state is required");}
		if(chunks<1){throw new IllegalArgumentException("Output queue capacity must be positive: "+chunks);}
		failure=failure_;
		queue=new ArrayBlockingQueue<Chunk>(chunks);
	}
	/** Budget ordinary queued payload across all active files. Preserve the established
	 * four-MiB minimum; larger heaps permit later nuclear files to finish before publication. */
	public static int queueCapacity(long heapBytes, int files){
		if(heapBytes<=0 || files<=0){throw new IllegalArgumentException("Output budgeting requires a positive heap and active-file count");}
		final long budget=heapBytes/8/files/CHUNK_BYTES;
		return (int)Math.max(MIN_CHUNKS, Math.min(MAX_CHUNKS, budget));
	}
	/** GFF is a byte stream, so bounded pieces may split a line without changing its bytes. */
	public void gff(ByteBuilder bytes){
		assert(bytes!=null) : "Every GFF batch, including empty output, must be published";
		if(bytes.length()==0){put(new Chunk(0, bytes, null));return;}
		for(int start=0; start<bytes.length(); start+=CHUNK_BYTES){
			final int length=Math.min(CHUNK_BYTES, bytes.length()-start);
			final ByteBuilder part=new ByteBuilder(length);
			part.append(bytes.array, start, length);
			put(new Chunk(0, part, null));
		}
	}
	/** Keep complete reads together; one oversized read is the unavoidable lower memory bound. */
	public void reads(int channel, ArrayList<Read> reads){
		assert(channel>=1 && channel<=4) : "Sequence channels are amino,16S,18S,ITS";
		if(reads==null || reads.isEmpty()){put(new Chunk(channel, null, new ArrayList<Read>(0)));return;}
		ArrayList<Read> part=new ArrayList<Read>();
		long bases=0;
		for(Read read:reads){
			part.add(read);
			bases+=read.length()+read.mateLength();
			if(bases>=CHUNK_BYTES){put(new Chunk(channel, null, part));part=new ArrayList<Read>();bases=0;}
		}
		if(!part.isEmpty()){put(new Chunk(channel, null, part));}
	}
	/** A terminal follows all output even for an empty input file. */
	public void finish(){put(END);}
	public void fail(Throwable cause){failure.fail(cause);}
	/** Returns the next immutable handoff, or null at this file's terminal. */
	public Chunk take(){
		try{
			while(true){
				failure.check();
				final Chunk chunk=queue.poll(100, TimeUnit.MILLISECONDS);
				if(chunk!=null){return chunk==END ? null : chunk;}
			}
		}catch(InterruptedException e){Thread.currentThread().interrupt();failure.fail(e);throw new IllegalStateException(e);}
	}
	private void put(Chunk chunk){
		assert(chunk!=null) : "Queue entries must be a payload or the explicit terminal";
		try{
			while(true){failure.check();if(queue.offer(chunk, 100, TimeUnit.MILLISECONDS)){return;}}
		}catch(InterruptedException e){Thread.currentThread().interrupt();failure.fail(e);throw new IllegalStateException(e);}
	}
	public static final class Chunk{
		private Chunk(int channel_, ByteBuilder gff_, ArrayList<Read> reads_){channel=channel_;gff=gff_;reads=reads_;}
		public final int channel;
		public final ByteBuilder gff;
		public final ArrayList<Read> reads;
	}
	private final Failure failure;
	// Ordinary queued payload is roughly bounded by one-eighth of the JVM heap,
	// subject to the established4MiB/file floor and64MiB/file ceiling. Read/caller
	// buffers, object overhead, headers and oversized individual reads are separate.
	private final ArrayBlockingQueue<Chunk> queue;
	private static final int MIN_CHUNKS=64, MAX_CHUNKS=1024;
	private static final int CHUNK_BYTES=65536;
	private static final Chunk END=new Chunk(-1, null, null);
}
