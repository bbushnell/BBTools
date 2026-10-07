package stream.bam;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

import stream.OrderedQueueSystem2;

/**
 * BGZF output with producer-side buffering and OrderedQueueSystem2 job routing.
 * Construction starts one writer and the requested number of compression workers;
 * both roles use ProcessThread, while OQS2 supplies the queues and completion flag.
 * One caller owns write/flush/close and the accumulation buffer. Submitted buffers
 * become jobs with ascending IDs; workers produce raw deflate payloads and footer
 * metadata, which the writer emits in order.
 *
 * flush() submits pending bytes and flushes the sink without waiting for queued work.
 * close() signals OQS2 and waits for its completion flag, not thread joins. The writer
 * normally emits EOF and closes the retained sink; any helper's caught failure can
 * also signal completion. No caller method rejects writes merely because closed is set.
 * @author Brian Bushnell
 * @contributor Collei
 * @date December 6, 2025
 */
public class BgzfOutputStreamMT2 extends OutputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts one compressor and one writer, using level 6 and 65280-byte input blocks. */
	public BgzfOutputStreamMT2(OutputStream out){this(out, 1, 6, DEFAULT_BLOCK_SIZE);}

	/** Starts the requested compressors and one writer with the default input-block capacity. */
	public BgzfOutputStreamMT2(OutputStream out, int threads, int compressionLevel){this(out, threads, compressionLevel, DEFAULT_BLOCK_SIZE);}

	/** Retains configuration, allocates the input buffer/queues and immediately starts helpers.
	 * Only nonnull out and positive threads are explicitly asserted here; other arguments
	 * are passed to allocation/compression code without additional constructor validation.
	 * @param out Sink used for data, flushes and normal writer-side closure
	 * @param threads Compressor count, excluding the separate writer
	 * @param compressionLevel Deflater level passed to each worker
	 * @param blockSize Accumulation-buffer capacity used as supplied */
	public BgzfOutputStreamMT2(OutputStream out, int threads, int compressionLevel, int blockSize){
		assert(out!=null);
		assert(threads>0);

		this.out=out;
		this.workerThreads=threads;
		this.compressionLevel=compressionLevel;
		this.maxBlockSize=blockSize;
		this.buffer=new byte[maxBlockSize];

		// Create OQS
		// Retain prototype instances for BgzfJob's poison/last marker factories.
		BgzfJob inputProto=new BgzfJob(0, null, null, false);
		BgzfJob outputProto=new BgzfJob(0, null, null, false);
		
		this.oqs=new OrderedQueueSystem2<>(threads, true, inputProto, outputProto);
		
		startThreads();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts workerThreads+1 helpers: ID zero writes, positive IDs compress. */
	private void startThreads(){
		alpt=new ArrayList<>(workerThreads+1);
		for(int i=0; i<=workerThreads; i++){
			alpt.add(new ProcessThread(i));
		}
		for(ProcessThread pt : alpt){pt.start();}
	}

	/** Checks the recorded error, buffers the low byte and submits a full block.
	 * Does not check the closed flag or wait for completion of submitted work.
	 * @throws IOException If a helper has recorded an error */
	@Override
	public void write(int b) throws IOException{
		checkError();
		buffer[bufferPos++]=(byte)b;
		if(bufferPos>=maxBlockSize){flushBlock();}
	}

	/** Copies a range into owned blocks, checking the recorded error once on entry.
	 * Bounds are enforced by arraycopy when reached; a nonpositive len skips copying.
	 * The source array is not retained, and the closed flag is not checked.
	 * @param b Source array
	 * @param off First source position
	 * @param len Number of bytes to copy
	 * @throws IOException If a helper has recorded an error */
	@Override
	public void write(byte[] b, int off, int len) throws IOException{
		checkError();
		while(len>0){
			int available=maxBlockSize-bufferPos;
			int toWrite=Math.min(available, len);
			System.arraycopy(b, off, buffer, bufferPos, toWrite);
			bufferPos+=toWrite;
			off+=toWrite;
			len-=toWrite;
			if(bufferPos>=maxBlockSize){flushBlock();}
		}
	}

	/** Transfers a nonempty buffer into an OQS input job, then allocates its replacement.
	 * The job retains the old array without copying. Empty buffers return immediately;
	 * this helper does not check errorState or closed itself. */
	private void flushBlock() throws IOException{
		if(bufferPos==0){return;}
		
		// Create job
		BgzfJob job=new BgzfJob(nextJobId++, buffer, null, false);
		synchronized(job){job.decompressedSize=bufferPos;}
		
		// Add to OQS (handles blocking if queue is full)
		oqs.addInput(job);

		// Reset buffer
		buffer=new byte[maxBlockSize];
		bufferPos=0;
	}

	/** Submits pending bytes and flushes the sink without awaiting compression or writes.
	 * Does not check the stored error or closed flag directly.
	 * @throws IOException If flushing the underlying sink fails */
	@Override
	public void flush() throws IOException{
		flushBlock();
		out.flush();
	}

	/** Submits pending bytes, signals OQS termination and waits for its completion flag.
	 * Sets closed after that wait, then checks errorState. Already-closed calls return
	 * without checking it again. The wait is not a helper-thread join; a writer or worker
	 * failure can also set the completion flag. The writer owns normal sink closure.
	 * @throws IOException If a recorded helper error is observed */
	@Override
	public void close() throws IOException{
		if(closed){return;}
		
		// Flush final data
		flushBlock();
		
		// Signal OQS to shutdown
		// This injects the Last job for the writer and Poison for the workers.
		oqs.poison();
		
		// Wait for OQS's externally set completion flag; this does not join helpers.
		oqs.waitForFinish();
		closed=true;
		
		// Check for errors one last time
		checkError();
	}
	
	/** Throws the currently recorded helper error, without checking the close latch. */
	private void checkError() throws IOException{
		if(errorState!=null){throw errorState;}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Helper role selected by tid: zero writes and positive IDs compress. */
	private class ProcessThread extends Thread{
		
		/** Retains the role ID and names the thread; daemon status remains inherited. */
		ProcessThread(int tid_){
			tid=tid_;
			setName("BGZF-"+(tid==0 ? "Writer" : "Worker"+tid));
		}
		
		/** Runs the assigned role; records caught failures and signals forced OQS completion. */
		@Override
		public void run(){
			try{
				if(tid==0){writerLoop();}else{workerLoop();}
			}catch(Throwable e){
				errorState=new IOException(getName()+" failed", e);
				oqs.setFinished(true); // Force shutdown
			}
		}
		
		/** Compresses ordered input jobs and submits payload/footer metadata to the output queue.
		 * Reuses a local Deflater and CRC32; ordinary loop exit ends the deflater. */
		private void workerLoop(){
			Deflater deflater=new Deflater(compressionLevel, true);
			if(FILTERED_BGZF){deflater.setStrategy(Deflater.FILTERED);}
			CRC32 crc=new CRC32();
			
			BgzfJob job=oqs.getInput();
			while(job!=null && !job.poison()){
				
				// Compression Logic
				synchronized(job){
					// CRC32
					crc.reset();
					crc.update(job.decompressed, 0, job.decompressedSize);
					long crcValue=crc.getValue();
					
					// Deflate
					deflater.reset();
					deflater.setInput(job.decompressed, 0, job.decompressedSize);
					deflater.finish();
					
					// Dynamic buffer sizing (simplified for brevity)
					int cap=maxBlockSize+1024;
					job.compressed=new byte[cap];
					int compressedSize=0;
					while(!deflater.finished()){
						int n=deflater.deflate(job.compressed, compressedSize, cap-compressedSize);
						compressedSize+=n;
						if(compressedSize==cap && !deflater.finished()){
							// Grow buffer
							cap*=2;
							job.compressed=Arrays.copyOf(job.compressed, cap);
						}
					}
					job.compressedSize=compressedSize;
					
					// Footer metadata
					byte[] meta=new byte[8];
					int mpos=0;
					mpos=writeInt32(meta, mpos, (int)crcValue);
					mpos=writeInt32(meta, mpos, job.decompressedSize);
					job.decompressed=meta; // Reuse field for footer
					job.decompressedSize=8;
				}
				
				oqs.addOutput(job);
				job=oqs.getInput();
			}
			
			deflater.end();
			// Historical cascade rationale follows. Current JobQueue retains poison and maps it
			//to null, so getInput normally takes the null branch here; the explicit re-add below
			//only applies to a nonnull marker. This comment is not a new liveness certification.
			// Re-inject poison for other workers
			//Each worker exits its loop on a poison job (getInput()==null OR job.poison()) and re-injects it
			//so the NEXT worker also wakes -> the poison cascades through all N workers (OQS2's poison
			//protocol, reviewed clean). This is why MT2 needs no manual per-worker poison accounting like the
			//original MT engine - OQS2 owns ordering+shutdown+backpressure, which is the whole point of the
			//rewrite (clean vs MT's hand-rolled poison bookkeeping). job is non-null here only on the poison
			//path (the while-guard exits on null OR poison; null can't be re-added, poison is).
			if(job!=null){oqs.addInput(job);}
		}
		
		/** Writes ordered payload jobs until null or LAST, then emits EOF and closes the sink.
		 * Successful finalization signals OQS completion; exceptions are handled by run(). */
		private void writerLoop() throws IOException{
			BgzfJob job=oqs.getOutput();
			while(job!=null && !job.last()){
				
				// Writer Logic
				synchronized(job){
					//Historical #001 rationale below describes the pre-guard version. The BSIZE assertion
					//is now present and the default input-block size is 65280; custom blockSize is still
					//used as supplied. No new runtime validation is claimed by this documentation update.
					//Historical bug [stream/bam/BgzfOutputStreamMT2#001] - BSIZE 16-bit overflow on an
					//incompressible full block. The worker's grow-loop (L181-189) correctly captures the full
					//compressed bytes (no truncation, unlike the ST writer), but for a 65536-byte incompressible
					//block compressedSize~=65556 -> bsize=18+65556+8-1=65581 > 65535, and writerHeader[16/17]
					//take only the low 16 bits -> BSIZE wraps -> corrupt BGZF framing, SILENTLY. The original
					//BgzfOutputStreamMT guards this with assert(bsize>=27 && bsize<=65535) (crash-loud under -ea);
					//MT2 has NO such guard. LOW (generic-.gz MT path, USE_BGZFOS_MT2=t; needs incompressible
					//full block, which BBTools' own compressible data never produces). Root cause shared with
					//#BgzfOutputStream#001: maxBlockSize=65536 leaves no BSIZE headroom (spec/samtools cap at
					//0xff00). Fix: cap uncompressed block <=0xff00, or at least assert(bsize<=65535) loud.
					int bsize=18+job.compressedSize+8-1;
					assert(bsize>=27 && bsize<=65535) : "BSIZE overflow: "+bsize+" (block too large for 16-bit "+
						"BSIZE; uncompressed block must be <=0xff00). [#001 - restored the guard the original MT "+
						"engine has; with DEFAULT_BLOCK_SIZE=0xff00 this can't fire on valid input)";

					// Header
					if(writerHeader==null){writerHeader=new byte[18];}
					System.arraycopy(GZIP_HEADER_TEMPLATE, 0, writerHeader, 0, 18);
					writerHeader[16]=(byte)(bsize&0xFF);
					writerHeader[17]=(byte)((bsize>>8)&0xFF);
					out.write(writerHeader, 0, 18);
					
					// Payload
					out.write(job.compressed, 0, job.compressedSize);
					
					// Footer
					out.write(job.decompressed, 0, 8);
				}
				
				job=oqs.getOutput();
			}
			
			// Finish
			writeEOFMarker();
			out.close();
			oqs.setFinished(true); // Signal completion to OQS
		}
		
		/** Role ID; zero selects writerLoop. */
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Writes the fixed 28-byte empty BGZF block and flushes the retained sink. */
	private void writeEOFMarker() throws IOException{
		byte[] eof=new byte[]{
			0x1f, (byte)0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
			0x00, (byte)0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
			0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
			0x00, 0x00, 0x00, 0x00
		};
		out.write(eof);
		out.flush();
	}

	/** Writes 32 bits in little-endian order and returns the position after four bytes. */
	private int writeInt32(byte[] buf, int pos, int val){
		buf[pos++]=(byte)(val&0xFF);
		buf[pos++]=(byte)((val>>8)&0xFF);
		buf[pos++]=(byte)((val>>16)&0xFF);
		buf[pos++]=(byte)((val>>24)&0xFF);
		return pos;
	}

	/** Fixed 18-byte header with BC subfield; the writer patches the two BSIZE bytes. */
	private static final byte[] GZIP_HEADER_TEMPLATE=new byte[]{
		31, (byte)139, 8, 4, 0, 0, 0, 0, 0, (byte)255, 6, 0, 66, 67, 2, 0, 0, 0
	};
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained sink; the writer writes/closes it and caller-side flush also flushes it. */
	private final OutputStream out;
	/** Number of compressors, excluding the writer. */
	private final int workerThreads;
	/** Deflater level supplied to each compressor. */
	private final int compressionLevel;
	/** Configured accumulation-buffer capacity, used without clamping. */
	private final int maxBlockSize;
	
	/** Input/output ordering, backpressure and external completion signaling. */
	private final OrderedQueueSystem2<BgzfJob, BgzfJob> oqs;
	/** Retained helper threads; close does not join this list. */
	private ArrayList<ProcessThread> alpt;
	
	/** Caller-owned accumulation buffer, transferred to a job on submission. */
	private byte[] buffer;
	/** Number of accumulated bytes. */
	private int bufferPos=0;
	/** Next input job ID, assigned sequentially from zero. */
	private long nextJobId=0;
	
	/** Writer-owned reusable 18-byte header buffer. */
	private byte[] writerHeader;
	/** Most recently assigned helper failure, initially null. */
	private volatile IOException errorState;
	/** Close-call latch set after waiting for OQS completion; write/flush do not check it. */
	private volatile boolean closed=false;
	
	/** Default uncompressed block capacity, 65280 bytes. */
	//[block-size lead FIXED 2026-06-20 (greenlit)] 0xff00 (65280), NOT 65536 - see BgzfOutputStreamMT: caps
	//the uncompressed block so compressed+overhead fits the 16-bit BSIZE (avoids the incompressible-block wrap).
	public static final int DEFAULT_BLOCK_SIZE=65280;
	/** Strategy preference sampled by each compressor at startup; set before construction. */
	public static boolean FILTERED_BGZF=false;
}
