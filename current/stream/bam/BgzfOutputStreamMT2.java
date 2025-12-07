package stream.bam;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

import stream.OrderedQueueSystem2;

/**
 * Multithreaded BGZF output stream using OrderedQueueSystem2.
 * Significant simplification of the original architecture:
 * - OQS2 manages all synchronization, ordering, and backpressure.
 * - Workers and Writer are unified into ProcessThreads.
 * - Explicit shutdown logic is replaced by OQS2's poison protocol.
 * @author Brian Bushnell
 * @contributor Collei
 * @date December 6, 2025
 */
public class BgzfOutputStreamMT2 extends OutputStream {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public BgzfOutputStreamMT2(OutputStream out){
		this(out, 1, 6, DEFAULT_BLOCK_SIZE);
	}

	public BgzfOutputStreamMT2(OutputStream out, int threads, int compressionLevel){
		this(out, threads, compressionLevel, DEFAULT_BLOCK_SIZE);
	}

	public BgzfOutputStreamMT2(OutputStream out, int threads, int compressionLevel, int blockSize){
		assert(out!=null);
		assert(threads>0);

		this.out=out;
		this.workerThreads=threads;
		this.compressionLevel=compressionLevel;
		this.maxBlockSize=blockSize;
		this.buffer=new byte[maxBlockSize];

		// Create OQS
		// Note: We need prototype instances for Poison/Last generation.
		// Assuming BgzfJob has a constructor or factory for this.
		BgzfJob inputProto=new BgzfJob(0, null, null, false);
		BgzfJob outputProto=new BgzfJob(0, null, null, false);
		
		this.oqs=new OrderedQueueSystem2<>(threads, true, inputProto, outputProto);
		
		startThreads();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void startThreads(){
		alpt=new ArrayList<>(workerThreads+1);
		for(int i=0; i<=workerThreads; i++){
			alpt.add(new ProcessThread(i));
		}
		for(ProcessThread pt : alpt){pt.start();}
	}

	@Override
	public void write(int b) throws IOException{
		checkError();
		buffer[bufferPos++]=(byte)b;
		if(bufferPos>=maxBlockSize){flushBlock();}
	}

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

	/** Submit current buffer to OQS. */
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

	@Override
	public void flush() throws IOException{
		flushBlock();
		out.flush();
	}

	@Override
	public void close() throws IOException{
		if(closed){return;}
		
		// Flush final data
		flushBlock();
		
		// Signal OQS to shutdown
		// This injects the Last job for the writer and Poison for the workers.
		oqs.poison();
		
		// Wait for completion
		oqs.waitForFinish();
		closed=true;
		
		// Check for errors one last time
		checkError();
	}
	
	private void checkError() throws IOException {
		if(errorState!=null) {throw errorState;}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		ProcessThread(int tid_){
			tid=tid_;
			setName("BGZF-"+(tid==0 ? "Writer" : "Worker"+tid));
		}
		
		@Override
		public void run(){
			try {
				if(tid==0) {writerLoop();}
				else {workerLoop();}
			} catch (Throwable e) {
				errorState=new IOException(getName()+" failed", e);
				oqs.setFinished(true); // Force shutdown
			}
		}
		
		private void workerLoop(){
			Deflater deflater=new Deflater(compressionLevel, true);
			if(FILTERED_BGZF) {deflater.setStrategy(Deflater.FILTERED);}
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
			// Re-inject poison for other workers
			if(job!=null) {oqs.addInput(job);}
		}
		
		private void writerLoop() throws IOException{
			BgzfJob job=oqs.getOutput();
			while(job!=null && !job.last()){
				
				// Writer Logic
				synchronized(job){
					int bsize=18+job.compressedSize+8-1;
					
					// Header
					if(writerHeader==null){writerHeader=new byte[18];}
					System.arraycopy(GZIP_HEADER_TEMPLATE, 0, writerHeader, 0, 18);
					writerHeader[16]=(byte)(bsize & 0xFF);
					writerHeader[17]=(byte)((bsize>>8) & 0xFF);
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
		
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
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

	private int writeInt32(byte[] buf, int pos, int val){
		buf[pos++]=(byte)(val & 0xFF);
		buf[pos++]=(byte)((val>>8) & 0xFF);
		buf[pos++]=(byte)((val>>16) & 0xFF);
		buf[pos++]=(byte)((val>>24) & 0xFF);
		return pos;
	}

	private static final byte[] GZIP_HEADER_TEMPLATE=new byte[]{
		31, (byte)139, 8, 4, 0, 0, 0, 0, 0, (byte)255, 6, 0, 66, 67, 2, 0, 0, 0
	};
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final OutputStream out;
	private final int workerThreads;
	private final int compressionLevel;
	private final int maxBlockSize;
	
	private final OrderedQueueSystem2<BgzfJob, BgzfJob> oqs;
	private ArrayList<ProcessThread> alpt;
	
	private byte[] buffer;
	private int bufferPos=0;
	private long nextJobId=0;
	
	private byte[] writerHeader;
	private volatile IOException errorState;
	private volatile boolean closed=false;
	
	public static final int DEFAULT_BLOCK_SIZE=65536;
	public static boolean FILTERED_BGZF=false;
}