package stream.bam;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.GZIPInputStream;
import java.util.zip.Inflater;

import stream.JobQueue;
import structures.BinaryByteWrapperLE;

/**
 * Multithreaded BGZF (Blocked GZIP Format) input stream.
 *
 * Architecture:
 * - Producer thread: Reads BGZF blocks from file, creates jobs with ascending IDs
 * - Worker thread(s): Decompresses blocks in parallel
 * - Consumer (main thread): Uses JobQueue to retrieve jobs in sequential order
 *
 * Jobs are automatically ordered by JobQueue to maintain sequential output
 * even when workers complete out of order.
 *
 * @author Chloe
 * @contributor Isla
 * @date November 14, 2025
 */
public class BgzfInputStreamMT2 extends InputStream {
	
	public static void main(String[] args) throws IOException{
		if(args.length<1){
			System.err.println("Usage: BgzfInputStreamMT2 <file.gz>");
			System.exit(1);
		}
		int threads=(args.length>1 ? Integer.parseInt(args[1]) : BgzfSettings.READ_THREADS);
		
		String filename=args[0];
		byte[] buffer=new byte[65536];
		
		long totalReads=0;
		long totalBytes=0;
		long startTime=System.nanoTime();
		
		try(InputStream fis=new java.io.FileInputStream(filename);
			InputStream bgzf=new BgzfInputStreamMT2(fis, threads)){
			
			int bytesRead;
			while((bytesRead=bgzf.read(buffer))>=0){
				totalReads++;
				totalBytes+=bytesRead;
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

	public BgzfInputStreamMT2(InputStream in){
		this(in, BgzfSettings.READ_THREADS);
	}

	public BgzfInputStreamMT2(InputStream in, int threads){
		assert(in!=null) : "Null input stream";
		assert(threads>0 && threads<=32) : "Invalid thread count: "+threads;

		this.in=in;
		this.workerThreads=threads;

		//Queue size: 3+workers*2 allows some buffering without excessive memory
		final int queueSize=3+(3*workerThreads)/2;
		this.inputQueue=new ArrayBlockingQueue<>(queueSize);
		this.jobQueue=new JobQueue<BgzfInputJob>(queueSize, true, true, 0);

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
		byte[] header=new byte[10];
		byte[] xlenBytes=new byte[2];
		byte[] trailer=new byte[8];
		
		try{
			while(!closed){
				BgzfInputJob job=readNextBlock(header, xlenBytes, trailer);
				if(job==null){
					if(!lastJobQueued){
						BgzfInputJob eofJob=new BgzfInputJob(nextJobId++, null, 0, 0, 0, true);
						while(eofJob!=null){
							try{
								inputQueue.put(eofJob);
								eofJob=null;
							}catch(InterruptedException e){
								synchronized(BgzfInputStreamMT2.this){
									if(closed){
										inputQueue.offer(eofJob);
										producerFinished=true;
										return;
									}
								}
							}
						}
						lastJobQueued=true;
					}
					break;
				}

				assert(job.compressed!=null || job.decompressed!=null) : 
					"Producer created job without data";
				if(job.compressed!=null){
					assert(job.compressedSize>0) : 
						"Producer created job with zero compressed size";
				}else{
					assert(job.decompressedSize>0) : 
						"Producer created job with zero decompressed size";
				}
				assert(!DEBUG || job.repOK()) : "Producer created invalid job";

				// Submit to input queue
				while(job!=null){
					try{
						inputQueue.put(job);
						job=null;
					}catch(InterruptedException e){
						synchronized(BgzfInputStreamMT2.this){
							if(closed){
								job=null;
								producerFinished=true;
								return;
							}
						}
					}
				}
			}
		}catch(IOException e){
			workerError=e;
		}finally{
			producerFinished=true;
			if(!lastJobQueued){
				inputQueue.offer(BgzfInputJob.POISON_PILL);
			}
		}
	}

	/** Read next BGZF block from input stream and create job. Returns null on EOF. */
	private BgzfInputJob readNextBlock(byte[] header, byte[] xlenBytes, byte[] trailer) throws IOException{
		while(true){
			if(readingPlainGzip){
				BgzfInputJob job=readPlainGzipChunk();
				if(job!=null){return job;}
				closePlainStream();
				continue;
			}

			//Read gzip header (minimum 10 bytes)
			int bytesRead=readFully(header, 0, header.length);
			if(bytesRead==0){return null;} //EOF
			if(bytesRead<header.length){
				throw new EOFException("Truncated BGZF block header");
			}

			//Verify gzip signature
			assert((header[0] & 0xFF)==31 && (header[1] & 0xFF)==139) : 
				"Not a gzip file: "+(header[0] & 0xFF)+", "+(header[1] & 0xFF);
			if((header[0] & 0xFF)!=31 || (header[1] & 0xFF)!=139){
				throw new IOException("Not a gzip file");
			}

			//Check compression method (should be 8=DEFLATE)
			if(header[2]!=8){
				throw new IOException("Unsupported compression method: "+header[2]);
			}

			//Check flags - FEXTRA must be set for BGZF
			int flags=header[3] & 0xFF;
			boolean fextra=(flags & 0x04)!=0;
			if(!fextra){
				startPlainGzip(header, null, null);
				continue;
			}

			//Read XLEN (2 bytes, little-endian)
			if(readFully(xlenBytes, 0, 2)<2){
				throw new EOFException("Truncated XLEN");
			}
			int xlen=((xlenBytes[1] & 0xFF)<<8) | (xlenBytes[0] & 0xFF);

			//Read extra field and find BC subfield
			byte[] extra=new byte[xlen];
			if(readFully(extra, 0, xlen)<xlen){
				throw new EOFException("Truncated extra field");
			}

			int bsize=findBsizeInExtra(extra, xlen);
			if(bsize<0){
				startPlainGzip(header, xlenBytes, extra);
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
			if(readFully(trailer, 0, 8)<8){
				throw new EOFException("Truncated block trailer");
			}

			//Extract expected CRC and size using wrapper
			BinaryByteWrapperLE wrapper=new BinaryByteWrapperLE(trailer);
			long expectedCrc=wrapper.getInt() & 0xFFFFFFFFL;
			int expectedSize=wrapper.getInt();

			if(compressedSize==0 && expectedSize==0){lastJobQueued=true;}

			BgzfInputJob job=new BgzfInputJob(nextJobId++, compressed, compressedSize, 
				expectedCrc, expectedSize, lastJobQueued);
			
			if(DEBUG && job.lastJob){
				System.err.println("Producer: job "+job.id+" marked LAST (compressed="+
					compressedSize+", expected="+expectedSize+")");
			}

			assert(!DEBUG || job.repOK()) : "readNextBlock created invalid job";
			return job;
		}
	}

	/** Find BC subfield in gzip extra field and extract BSIZE. */
	private int findBsizeInExtra(byte[] extra, int xlen){
		int pos=0;
		while(pos+4<=xlen){
			int si1=extra[pos] & 0xFF;
			int si2=extra[pos+1] & 0xFF;
			int slen=((extra[pos+3] & 0xFF)<<8) | (extra[pos+2] & 0xFF);

			if(si1==66 && si2==67){ //'B' 'C'
				if(slen==2 && pos+6<=xlen){
					return ((extra[pos+5] & 0xFF)<<8) | (extra[pos+4] & 0xFF);
				}
			}

			pos+=4+slen;
		}
		return -1;
	}

	private void startPlainGzip(byte[] header, byte[] xlenBytes, byte[] extra) 
		throws IOException{
		if(readingPlainGzip){return;}
		byte[] prefix=buildPrefix(header, xlenBytes, extra);
		plainGzipStream=new GZIPInputStream(new PrefixedInputStream(prefix, in), 65536);
		readingPlainGzip=true;
	}

	private BgzfInputJob readPlainGzipChunk() throws IOException{
		if(plainGzipStream==null){return null;}

		byte[] buffer=new byte[65536];
		int total=0;
		while(total<buffer.length){
			int n=plainGzipStream.read(buffer, total, buffer.length-total);
			if(n<0){break;}
			if(n==0){continue;}
			total+=n;
			break;
		}

		if(total<=0){return null;}

		BgzfInputJob job=new BgzfInputJob(nextJobId++, null, 0, 0, 0, false);
		job.decompressed=buffer;
		job.decompressedSize=total;
		return job;
	}

	private void closePlainStream() throws IOException{
		if(plainGzipStream!=null){
			plainGzipStream.close();
			plainGzipStream=null;
		}
		readingPlainGzip=false;
	}

	private byte[] buildPrefix(byte[] header, byte[] xlenBytes, byte[] extra){
		int prefixLen=header.length;
		if(xlenBytes!=null){prefixLen+=xlenBytes.length;}
		if(extra!=null){prefixLen+=extra.length;}

		byte[] prefix=new byte[prefixLen];
		int pos=0;
		System.arraycopy(header, 0, prefix, pos, header.length);
		pos+=header.length;
		if(xlenBytes!=null){
			System.arraycopy(xlenBytes, 0, prefix, pos, xlenBytes.length);
			pos+=xlenBytes.length;
		}
		if(extra!=null && extra.length>0){
			System.arraycopy(extra, 0, prefix, pos, extra.length);
		}
		return prefix;
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
		Inflater inflater=new Inflater(true);

		try{
			while(!closed){
				BgzfInputJob job=inputQueue.take();

				if(DEBUG){
					System.err.println("Worker: dequeued job "+job.id+
						(job.isPoisonPill() ? " (POISON)" : ""));
				}

				if(job.isPoisonPill()){
					if(DEBUG){
						System.err.println("Worker: got POISON, terminating");
					}
					break;
				}

				if(job.compressed==null){
					assert(job.decompressed!=null) : 
						"Job "+job.id+" missing decompressed data";
					if(DEBUG){
						System.err.println("Worker: job "+job.id+
							" already decompressed ("+job.decompressedSize+" bytes)");
					}
				}else{
					assert(!DEBUG || job.repOK()) : "Worker received invalid job";

					//Decompress
					inflater.reset();
					inflater.setInput(job.compressed, 0, job.compressedSize);

					if(DEBUG){System.err.println("Worker: inflating job "+job.id);}

					// TODO: Brian to move this allocation here from producer
					final byte[] decompressed=new byte[65536];
					
					try{
						job.decompressedSize=inflater.inflate(decompressed, 0, 65536);
						job.decompressed=decompressed;
						if(DEBUG){
							System.err.println("Worker: inflated job "+job.id+" to "+
								job.decompressedSize+" bytes");
						}
					}catch(DataFormatException e){
						job.error=new IOException("Decompression failed for job "+
							job.id, e);
						if(!jobQueue.add(job)){break;}
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
						if(!jobQueue.add(job)){break;}
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
						if(!jobQueue.add(job)){break;}
						continue;
					}

					assert(!DEBUG || job.repOK()) : "Worker produced invalid job";
				}

				//Add to job queue
				if(!jobQueue.add(job)){break;}

				if(job.lastJob){
					if(DEBUG){
						System.err.println("Worker: job "+job.id+
							" marked LAST, injecting POISON to wake others");
					}
					final int toSignal=Math.max(1, workerThreads-1);
					int signaled=0;
					while(signaled<toSignal){
						if(inputQueue.offer(BgzfInputJob.POISON_PILL)){signaled++;}
						else{Thread.yield();}
					}
				}
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
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
				BgzfInputJob nextJob=jobQueue.take();
				if(nextJob==null){
					eofReached=true;
					return totalRead==0 ? -1 : totalRead;
				}

				if(DEBUG){
					System.err.println("Consumer: received job "+nextJob.id+
						", decompressedSize="+nextJob.decompressedSize+
						", decompressed="+(nextJob.decompressed!=null ? 
							"not null" : "NULL"));
				}

				if(nextJob.error!=null){
					throw new IOException("Decompression failed", nextJob.error);
				}

				assert(nextJob.decompressed!=null) : 
					"Job "+nextJob.id+" has null decompressed data";

				if(nextJob.lastJob){
					if(DEBUG){
						System.err.println("Consumer: job "+nextJob.id+
							" marked LAST, returning EOF");
					}
					eofReached=true;
					return totalRead==0 ? -1 : totalRead;
				}

				if(nextJob.decompressedSize==0){
					if(DEBUG){
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
	public void close() throws IOException{
		if(closed){return;}

		closed=true;

		try{closePlainStream();}catch(IOException ignore){}
		try{in.close();}catch(IOException ignore){}

		inputQueue.offer(BgzfInputJob.POISON_PILL);

		if(producer!=null){producer.interrupt();}
		if(workers!=null){
			for(Thread worker : workers){worker.interrupt();}
		}

		try{
			if(producer!=null){producer.join(10);}
			if(workers!=null){
				for(Thread worker : workers){worker.join(10);}
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
		}
	}

	private boolean repOK(){
		if(in==null){return false;}
		if(workerThreads<=0 || workerThreads>32){return false;}
		if(inputQueue==null || jobQueue==null){return false;}
		if(currentBlockPos<0 || currentBlockSize<0){return false;}
		if(currentBlockPos>currentBlockSize){return false;}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int workerThreads;
	private final ArrayBlockingQueue<BgzfInputJob> inputQueue;
	private final JobQueue<BgzfInputJob> jobQueue;
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
	private boolean lastJobQueued=false;
	private boolean eofReached=false;
	private static final boolean DEBUG=false;
	private GZIPInputStream plainGzipStream=null;
	private boolean readingPlainGzip=false;
}