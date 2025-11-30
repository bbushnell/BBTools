package stream.bam;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.GZIPInputStream;
import java.util.zip.Inflater;

import stream.OrderedQueueSystem;
import structures.BinaryByteWrapperLE;
import template.ThreadWaiter;

/**
 * Multithreaded BGZF (Blocked GZIP Format) input stream using OrderedQueueSystem.
 *
 * Architecture:
 * - Producer thread: Reads BGZF blocks from file, creates jobs with ascending IDs
 * - Worker thread(s): Decompresses blocks in parallel
 * - Consumer (main thread): Uses OQS to retrieve jobs in sequential order
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 15, 2025
 */
public class BgzfInputStreamMT3 extends InputStream {
	
	public static void main(String[] args) throws IOException{
		if(args.length<1){
			System.err.println("Usage: BgzfInputStreamMT3 <file.gz>");
			System.exit(1);
		}
		int threads=(args.length>1 ? Integer.parseInt(args[1]) : BgzfSettings.READ_THREADS);
		boolean write=args.length>2;
		long maxBlocks=(args.length>3 ? Long.parseLong(args[3]) : Long.MAX_VALUE);
		
		
		String filename=args[0];
		byte[] buffer=new byte[131072];
		
		long totalReads=0;
		long totalBytes=0;
		long startTime=System.nanoTime();
		
		try(InputStream fis=new java.io.FileInputStream(filename);
			InputStream bis=new BufferedInputStream(fis, 131072);
			InputStream bgzf=new BgzfInputStreamMT3(bis, threads)){
			
			int bytesRead;
			while((bytesRead=bgzf.read(buffer))>=0 && totalReads<maxBlocks){
				totalReads++;
				totalBytes+=bytesRead;
				if(write){System.out.write(buffer, 0, bytesRead);}
			}
		}
		long endTime=System.nanoTime();
		float seconds=(endTime-startTime)/1e9f;
		
		System.err.println("Read operations: "+totalReads);
		System.err.println("Total bytes:     "+totalBytes);
		System.err.println("Time:            "+String.format("%.3f", seconds)+" seconds");
		System.err.println("Throughput:      "+String.format("%.2f", totalBytes/seconds/1e6)+" MB/s");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public BgzfInputStreamMT3(InputStream in){
		this(in, BgzfSettings.READ_THREADS);
	}

	public BgzfInputStreamMT3(InputStream in, int threads){
		assert(in!=null) : "Null input stream";
		assert(threads>0 && threads<=32) : "Invalid thread count: "+threads;

		this.in=in;
		this.workerThreads=threads;

		// Create OQS with prototypes
		BgzfInputJob inputPrototype=new BgzfInputJob(0, null, 0, 0, false);
		BgzfInputJob outputPrototype=new BgzfInputJob(0, null, 0, 0, false);
		this.oqs=new OrderedQueueSystem<BgzfInputJob, BgzfInputJob>(
			threads, true, inputPrototype, outputPrototype);

		startThreads();

		assert(repOK()) : "Constructor postcondition failed";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Start producer and worker threads. */
	private void startThreads(){
		assert(producer==null) : "Threads already started";
		assert(workers==null) : "Workers already started";

		//Start producer thread
		producer=new Thread(new Runnable(){
			public void run(){producerLoop();}
		}, "BGZF-InputProducer");
		producer.setDaemon(true);
		producer.start();

		//Start worker threads
		workers=new Thread[workerThreads];
		for(int i=0; i<workerThreads; i++){
			final int threadNum=i;
			workers[i]=new Thread(new Runnable(){
				public void run(){workerLoop();}
			}, "BGZF-InputWorker-"+threadNum);
			workers[i].setDaemon(true);
			workers[i].start();
		}
	}

	/** Producer thread: Read BGZF blocks and create jobs. */
	private void producerLoop(){
		// Declare temp arrays outside loop for reuse
		final byte[] fields=new byte[10];
		final byte[] xlenBytes=new byte[2];
		final byte[] extra=new byte[1024];
		
		try{
			while(!closed){
				BgzfInputJob job=readNextBlock(fields, xlenBytes, extra);
				if(job==null){
					break;
				}

				assert(job.compressed!=null || job.decompressed!=null) : 
					"Producer created job without data";
				if(job.compressed!=null){
					assert(job.compressed.length>0) : 
						"Producer created job with zero compressed size";
				}else{
					assert(job.decompressedSize>0) : 
						"Producer created job with zero decompressed size";
				}
				assert(!DEBUG || job.repOK()) : "Producer created invalid job";

				oqs.addInput(job);
			}
		}catch(IOException e){
			workerError=e;
		}finally{
			producerFinished=true;
			oqs.poison();
		}
	}

	/** Read next BGZF block from input stream and create job. Returns null on EOF. */
	private BgzfInputJob readNextBlock(byte[] fields, byte[] xlenBytes, byte[] extra) 
			throws IOException{
		while(true){
			if(readingPlainGzip){
				BgzfInputJob job=readPlainGzipChunk();
				if(job!=null){return job;}
				closePlainStream();
				continue;
			}

			//Read gzip header (minimum 10 bytes)
			int bytesRead=readFully(fields, 0, fields.length);
			if(bytesRead==0){return null;} //EOF
			if(bytesRead<fields.length){
				throw new EOFException("Truncated BGZF block header");
			}

			//Verify gzip signature
			assert((fields[0] & 0xFF)==31 && (fields[1] & 0xFF)==139) : 
				"Not a gzip file: "+(fields[0] & 0xFF)+", "+(fields[1] & 0xFF);
			if((fields[0] & 0xFF)!=31 || (fields[1] & 0xFF)!=139){
				throw new IOException("Not a gzip file");
			}

			//Check compression method (should be 8=DEFLATE)
			if(fields[2]!=8){
				throw new IOException("Unsupported compression method: "+fields[2]);
			}

			//Check flags - FEXTRA must be set for BGZF
			int flags=fields[3] & 0xFF;
			boolean fextra=(flags & 0x04)!=0;
			if(!fextra){
				startPlainGzip(fields, null, null, 0);
				continue;
			}

			//Read XLEN (2 bytes, little-endian)
			if(readFully(xlenBytes, 0, 2)<2){
				throw new EOFException("Truncated XLEN");
			}
			int xlen=((xlenBytes[1] & 0xFF)<<8) | (xlenBytes[0] & 0xFF);

			//Read extra field and find BC subfield
			if(xlen>extra.length) {extra=new byte[xlen];}
			if(readFully(extra, 0, xlen)<xlen){
				throw new EOFException("Truncated extra field");
			}

			int bsize=findBsizeInExtra(extra, xlen);
			if(bsize<0){
				startPlainGzip(fields, xlenBytes, extra, xlen);
				continue;
			}

			//Calculate compressed data length
			int alreadyRead=10+2+xlen;
			int remaining=(bsize+1)-alreadyRead;
			if(remaining<8){throw new IOException("Invalid BSIZE: "+bsize);}

			final int compressedSize=remaining-8;

			//Read compressed data
			final byte[] compressed=new byte[compressedSize];
			if(readFully(compressed, 0, compressedSize)<compressedSize){
				throw new EOFException("Truncated compressed data");
			}

			//Read trailer (CRC32+ISIZE)
			if(readFully(fields, 0, 8)<8){
				throw new EOFException("Truncated block trailer");
			}

			//Extract expected CRC and size using wrapper
			BinaryByteWrapperLE wrapper=new BinaryByteWrapperLE(fields);
			long expectedCrc=wrapper.getInt() & 0xFFFFFFFFL;
			int expectedSize=wrapper.getInt();

			boolean isLast=(compressedSize==0 && expectedSize==0);

			BgzfInputJob job=new BgzfInputJob(nextJobId++, compressed, 
				expectedCrc, expectedSize, isLast);
			
			if(verbose && job.lastJob){
				System.err.println("Producer: job "+job.id+" marked LAST (compressed="+
					compressedSize+", expected="+expectedSize+")");
			}

			assert(!DEBUG || job.repOK()) : "readNextBlock created invalid job";
			return job;
		}
	}
	
	/** Find BC subfield in gzip extra field and extract BSIZE. */
	private final int findBsizeInExtra(final byte[] extra, final int xlen){
		for(int pos=0, lim=xlen-6; pos<=lim;){
			final int slen=((extra[pos+3] & 0xFF)<<8) | (extra[pos+2] & 0xFF);
			if(extra[pos]=='B' && extra[pos+1]=='C' && slen==2){
				return ((extra[pos+5] & 0xFF)<<8) | (extra[pos+4] & 0xFF);
			}
			pos+=4+slen;
		}
		return -1;
	}

	private void startPlainGzip(byte[] header, byte[] xlenBytes, byte[] extra, int xlen) 
		throws IOException{
		if(readingPlainGzip){return;}
		byte[] prefix=buildPrefix(header, xlenBytes, extra, xlen);
		plainGzipStream=new GZIPInputStream(new PrefixedInputStream(prefix, in), 65536);
		readingPlainGzip=true;
	}

	private byte[] buildPrefix(byte[] header, byte[] xlenBytes, byte[] extra, final int xlen){
		int prefixLen=header.length;
		if(xlenBytes!=null){prefixLen+=xlenBytes.length;}
		prefixLen+=xlen;

		byte[] prefix=new byte[prefixLen];
		int pos=0;
		System.arraycopy(header, 0, prefix, pos, header.length);
		pos+=header.length;
		if(xlenBytes!=null){
			System.arraycopy(xlenBytes, 0, prefix, pos, xlenBytes.length);
			pos+=xlenBytes.length;
		}
		if(xlen>0){
			System.arraycopy(extra, 0, prefix, pos,xlen);
		}
		return prefix;
	}

	private BgzfInputJob readPlainGzipChunk() throws IOException{
		if(plainGzipStream==null){return null;}

		final byte[] buffer=new byte[65536];
		int total=0;
		while(total<buffer.length){
			int n=0;
			try {//Protects from a closing race condition
				n=plainGzipStream.read(buffer, total, buffer.length-total);
			} catch (Exception e) {
				if(closed) {return null;} // Expected shutdown error
				throw e; // Real error
			}
			if(n<0){break;}
			if(n==0){continue;}
			total+=n;
			break;
		}

		if(total<=0){return null;}

		BgzfInputJob job=new BgzfInputJob(nextJobId++, null, 0, 0, false);
		synchronized(job){
			job.decompressed=buffer;
			job.decompressedSize=total;
		}
		return job;
	}

	private void closePlainStream() throws IOException{
		if(plainGzipStream!=null){
			plainGzipStream.close();
			plainGzipStream=null;
		}
		readingPlainGzip=false;
	}

	private static final class PrefixedInputStream extends InputStream{
		private final byte[] prefix;
		private int position=0;
		private final InputStream tail;

		PrefixedInputStream(byte[] prefix, InputStream tail){
			this.prefix=prefix;
			this.tail=tail;
		}

		@Override
		public int read() throws IOException{
			if(position<prefix.length){return prefix[position++] & 0xFF;}
			return tail.read();
		}

		@Override
		public int read(byte[] b, int off, int len) throws IOException{
			if(position<prefix.length){
				int toCopy=Math.min(len, prefix.length-position);
				System.arraycopy(prefix, position, b, off, toCopy);
				position+=toCopy;
				return toCopy;
			}
			return tail.read(b, off, len);
		}

		@Override
		public void close(){
			//Do not close tail stream; caller manages lifecycle
		}
	}

	/** Worker thread: Decompress BGZF blocks. */
	private void workerLoop(){
		final Inflater inflater=new Inflater(true);

		try{
			while(!closed){
				BgzfInputJob job=oqs.getInput();

				if(verbose){
					System.err.println("Worker: dequeued job "+job.id+
						(job.isPoisonPill() ? " (POISON)" : ""));
				}

				if(job.isPoisonPill()){
					if(verbose){
						System.err.println("Worker: got POISON, terminating");
					}
					//Re-inject for other workers
					oqs.addInput(job);
					break;
				}

				if(job.compressed==null){
					assert(job.decompressed!=null || job.last()) : 
						"Job "+job.id+" missing decompressed data";
					if(verbose){
						System.err.println("Worker: job "+job.id+
							" already decompressed ("+job.decompressedSize+" bytes)");
					}
				}else{
					assert(!DEBUG || job.repOK()) : "Worker received invalid job";

					//Decompress
					inflater.reset();
					inflater.setInput(job.compressed, 0, job.compressed.length);

					if(verbose){System.err.println("Worker: inflating job "+job.id);}
					final byte[] decompressed=new byte[65536];
					
					try{
						synchronized(job){
							job.decompressedSize=inflater.inflate(decompressed, 0, 65536);
							job.decompressed=decompressed;
						}
						if(verbose){
							System.err.println("Worker: inflated job "+job.id+" to "+
								job.decompressedSize+" bytes");
						}
					}catch(DataFormatException e){
						job.error=new IOException("Decompression failed for job "+
							job.id, e);
						oqs.addOutput(job);
						continue;
					}

					//Validate decompressed size
					assert(job.decompressedSize==job.expectedSize) : 
						"Size mismatch for job "+job.id+": expected "+job.expectedSize+
						", got "+job.decompressedSize;
					if(job.decompressedSize!=job.expectedSize){
						job.error=new IOException("Uncompressed size mismatch for job "+
							job.id+": expected "+job.expectedSize+", got "+
							job.decompressedSize);
						oqs.addOutput(job);
						continue;
					}

					//Verify CRC32
					CRC32 crc=new CRC32();
					crc.update(job.decompressed, 0, job.decompressedSize);
					long actualCrc=crc.getValue();

					assert(actualCrc==job.expectedCrc) : 
						"CRC32 mismatch for job "+job.id+": expected "+job.expectedCrc+
						", got "+actualCrc;
					if(actualCrc!=job.expectedCrc){
						job.error=new IOException("CRC32 mismatch for job "+job.id);
						oqs.addOutput(job);
						continue;
					}

					assert(!DEBUG || job.repOK()) : "Worker produced invalid job";
				}

				//Add to output queue
				oqs.addOutput(job);
			}
		}catch(Exception e){
			e.printStackTrace();
		}finally{
			inflater.end();
		}
	}

	@Override
	public int read() throws IOException{
		byte[] b=new byte[1];
		int n=read(b, 0, 1);
		return n<0 ? -1 : (b[0] & 0xFF);
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException{
		assert(b!=null) : "Null buffer";
		assert(off>=0 && len>=0 && len<=b.length-off) : 
			"Invalid offset/length: off="+off+", len="+len+", buf.length="+b.length;

		if(closed){throw new IOException("Stream closed");}
		if(len==0){return 0;}
		if(eofReached){return -1;}

		if(workerError!=null){throw workerError;}

		int totalRead=0;
		while(totalRead<len){
			if(currentBlockPos>=currentBlockSize){
				BgzfInputJob nextJob=oqs.getOutput();
				if(nextJob==null){
					eofReached=true;
					oqs.setFinished(true);
					return totalRead==0 ? -1 : totalRead;
				}

				if(verbose){
					System.err.println("Consumer: received job "+nextJob.id+
						", decompressedSize="+nextJob.decompressedSize+
						", decompressed="+(nextJob.decompressed!=null ? 
							"not null" : "NULL"));
				}

				if(nextJob.error!=null){
					throw new IOException("Decompression failed", nextJob.error);
				}

				assert(nextJob.decompressed!=null || nextJob.lastJob) : 
					"Job "+nextJob.id+" has null decompressed data";

				if(nextJob.lastJob){
					if(verbose){
						System.err.println("Consumer: job "+nextJob.id+
							" marked LAST, returning EOF");
					}
					eofReached=true;
					oqs.setFinished(true);
					return totalRead==0 ? -1 : totalRead;
				}

				if(nextJob.decompressedSize==0){
					if(verbose){
						System.err.println("Consumer: job "+nextJob.id+
							" has 0 decompressed bytes, skipping");
					}
					continue;
				}

				assert(nextJob.decompressedSize>0) : 
					"Job "+nextJob.id+" has zero decompressed size";

				currentBlock=nextJob.decompressed;
				currentBlockSize=nextJob.decompressedSize;
				currentBlockPos=0;
			}

			int available=currentBlockSize-currentBlockPos;
			int toCopy=Math.min(available, len-totalRead);

			assert(toCopy>0) : "toCopy should be positive: "+toCopy;
			assert(currentBlockPos+toCopy<=currentBlockSize) : 
				"Copy would exceed block: pos="+currentBlockPos+", toCopy="+toCopy+
				", size="+currentBlockSize;

			System.arraycopy(currentBlock, currentBlockPos, b, off+totalRead, toCopy);
			currentBlockPos+=toCopy;
			totalRead+=toCopy;
		}

		assert(totalRead==len) : "Read wrong amount: expected "+len+", got "+totalRead;
		return totalRead;
	}

	private int readFully(byte[] b, int off, int len) throws IOException{
		int total=0;
		while(total<len){
			int n=in.read(b, off+total, len-total);
			if(n<0){return total;}
			total+=n;
		}
		return total;
	}

	@Override
	public synchronized void close() throws IOException{
		if(closed){return;}
		if(verbose)
		closed=true;
		if(verbose) {System.err.println("Called close.");}
		try{closePlainStream();}catch(IOException ignore){}
		try{in.close();}catch(IOException ignore){}
		if(verbose) {System.err.println("Calling setFinished.");}
		oqs.setFinished(true);
		if(verbose) {System.err.println("Wiating for workers.");}
		ThreadWaiter.waitForThreadsToFinish(workers);//Not strictly needed for daemons
		if(verbose) {System.err.println("Close finished.");}
	}

	private boolean repOK(){
		if(in==null){return false;}
		if(workerThreads<=0 || workerThreads>32){return false;}
		if(oqs==null){return false;}
		if(currentBlockPos<0 || currentBlockSize<0){return false;}
		if(currentBlockPos>currentBlockSize){return false;}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int workerThreads;
	private final OrderedQueueSystem<BgzfInputJob, BgzfInputJob> oqs;
	private Thread producer;
	private Thread[] workers;
	private long nextJobId=0;
	private byte[] currentBlock;
	private int currentBlockPos=0;
	private int currentBlockSize=0;
	private final InputStream in;
	private volatile IOException workerError=null;
	private volatile boolean producerFinished=false;
	private volatile boolean closed=false;
	private boolean eofReached=false;
	private static final boolean verbose=false;
	private static final boolean DEBUG=false;
	private GZIPInputStream plainGzipStream=null;
	private boolean readingPlainGzip=false;
}