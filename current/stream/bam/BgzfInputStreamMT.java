package stream.bam;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.PriorityQueue;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * Multithreaded BGZF (Blocked GZIP Format) input stream.
 *
 * Architecture:
 * - Producer thread: Reads BGZF blocks from file, creates jobs with ascending IDs
 * - Worker thread(s): Decompresses blocks in parallel
 * - Consumer (main thread): Waits for next sequential job, returns decompressed data
 *
 * Jobs are ordered by ID in a PriorityQueue heap to maintain sequential output
 * even when workers complete out of order.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BgzfInputStreamMT extends InputStream {

	/** Number of worker threads (start with 1 for correctness, then scale up) */
	private final int workerThreads;

	/** Input queue for jobs to be processed */
	private final ArrayBlockingQueue<BgzfJob> inputQueue;

	/** Output heap maintaining sequential order */
	private final PriorityQueue<BgzfJob> outputHeap;

	/** Maximum number of jobs allowed in output heap before workers must wait */
	private final int maxHeapSize;

	/** Producer thread reading BGZF blocks */
	private Thread producer;

	/** Worker threads decompressing blocks */
	private Thread[] workers;

	/** Next job ID expected by consumer */
	private long nextExpectedId = 0;

	/** Next job ID to assign by producer */
	private long nextJobId = 0;

	/** Current decompressed block being read from */
	private byte[] currentBlock;

	/** Position in current block */
	private int currentBlockPos = 0;

	/** Size of current block */
	private int currentBlockSize = 0;

	/** Underlying input stream */
	private final InputStream in;

	/** Error state from worker threads */
	private volatile IOException workerError = null;

	/** Producer finished flag */
	private volatile boolean producerFinished = false;

	/** Stream closed flag */
	private volatile boolean closed = false;

	/** Whether a last-job marker has already been queued */
	private boolean lastJobQueued = false;

	/** Lock protecting last-job state */
	private final Object lastJobLock = new Object();

	/** Whether EOF has already been delivered to the caller */
	private boolean eofReached = false;

	/** Debug flag (enable with -Dbgzf.debug=true) */
	private static final boolean DEBUG = false;//Boolean.getBoolean("bgzf.debug");

	public BgzfInputStreamMT(InputStream in) {
		this(in, 1); // Start with 1 worker for correctness
	}

	public BgzfInputStreamMT(InputStream in, int threads) {
		assert in != null : "Null input stream";
		assert threads > 0 && threads <= 32 : "Invalid thread count: " + threads;

		this.in = in;
		this.workerThreads = threads;

		// Queue size: workers * 2 allows some buffering without excessive memory
		int queueSize = Math.max(workerThreads * 2, 2);
		this.inputQueue = new ArrayBlockingQueue<>(queueSize);
		this.outputHeap = new PriorityQueue<>(queueSize);
		this.maxHeapSize = queueSize;

		startThreads();

		assert repOK() : "Constructor postcondition failed";
	}

	/**
	 * Start producer and worker threads.
	 */
	private void startThreads() {
		assert producer == null : "Threads already started";
		assert workers == null : "Workers already started";

		// Start producer thread
		producer = new Thread(this::producerLoop, "BGZF-Producer");
		producer.setDaemon(true);
		producer.start();

		// Start worker threads
		workers = new Thread[workerThreads];
		for (int i = 0; i < workerThreads; i++) {
			workers[i] = new Thread(this::workerLoop, "BGZF-Worker-" + i);
			workers[i].setDaemon(true);
			workers[i].start();
		}

		// Don't assert isAlive() immediately - thread might not be scheduled yet
	}

	/**
	 * Producer thread: Read BGZF blocks and create jobs.
	 */
	private void producerLoop() {
		try {
			while (!closed) {
				// Read next BGZF block
				BgzfJob job = readNextBlock();
				if (job == null) {
					// EOF reached
					break;
				}

				assert job.compressed != null : "Producer created job with null compressed data";
				assert job.compressedSize > 0 : "Producer created job with zero compressed size";
				assert(!DEBUG || job.repOK()) : "Producer created invalid job";

				// Submit to input queue (blocks if queue full)
				inputQueue.put(job);
			}
		} catch (IOException e) {
			workerError = e;
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
		} finally {
			producerFinished = true;
			if (!lastJobQueued) {
				try {
					inputQueue.put(BgzfJob.POISON_PILL);
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
				}
			}
		}
	}

	/**
	 * Read next BGZF block from input stream and create job.
	 * @return Job with compressed data, or null on EOF
	 */
	private BgzfJob readNextBlock() throws IOException {
		// Read gzip header (minimum 10 bytes)
		byte[] header = new byte[12];
		int bytesRead = readFully(header, 0, 10);
		if (bytesRead == 0) {
			return null; // EOF
		}
		if (bytesRead < 10) {
			throw new EOFException("Truncated BGZF block header");
		}

		// Verify gzip signature
		assert((header[0] & 0xFF) == 31 && (header[1] & 0xFF) == 139) :
			"Not a gzip file: " + (header[0] & 0xFF) + ", " + (header[1] & 0xFF);
		if ((header[0] & 0xFF) != 31 || (header[1] & 0xFF) != 139) {
			throw new IOException("Not a gzip file");
		}

		// Check compression method (should be 8 = DEFLATE)
		if (header[2] != 8) {
			throw new IOException("Unsupported compression method: " + header[2]);
		}

		// Check flags - FEXTRA must be set for BGZF
		int flags = header[3] & 0xFF;
		boolean fextra = (flags & 0x04) != 0;
		assert fextra : "BGZF block missing FEXTRA flag";
		if (!fextra) {
			throw new IOException("BGZF block missing FEXTRA flag");
		}

		// Read XLEN (2 bytes, little-endian)
		if (readFully(header, 0, 2) < 2) {
			throw new EOFException("Truncated XLEN");
		}
		int xlen = ((header[1] & 0xFF) << 8) | (header[0] & 0xFF);

		// Read extra field and find BC subfield
		byte[] extra = new byte[xlen];
		if (readFully(extra, 0, xlen) < xlen) {
			throw new EOFException("Truncated extra field");
		}

		int bsize = findBsizeInExtra(extra, xlen);
		assert bsize >= 0 : "BGZF block missing BC subfield";
		if (bsize < 0) {
			throw new IOException("BGZF block missing BC subfield");
		}

		// Calculate compressed data length
		int alreadyRead = 10 + 2 + xlen;
		int remaining = (bsize + 1) - alreadyRead;
		assert remaining >= 8 : "Invalid BSIZE: " + bsize;
		if (remaining < 8) {
			throw new IOException("Invalid BSIZE: " + bsize);
		}

		int compressedSize = remaining - 8; // Subtract CRC32 and ISIZE

		// Read complete block into job
		BgzfJob job = new BgzfJob(nextJobId++);
		job.compressed = new byte[compressedSize];
		job.compressedSize = compressedSize;

		if (readFully(job.compressed, 0, compressedSize) < compressedSize) {
			throw new EOFException("Truncated compressed data");
		}

		// Read trailer (CRC32 + ISIZE)
		byte[] trailer = new byte[8];
		if (readFully(trailer, 0, 8) < 8) {
			throw new EOFException("Truncated block trailer");
		}

		// Store trailer in job for worker validation
		ByteBuffer bb = ByteBuffer.wrap(trailer).order(ByteOrder.LITTLE_ENDIAN);
		long expectedCrc = bb.getInt() & 0xFFFFFFFFL;
		int expectedSize = bb.getInt();

		// Store expected values for worker validation
		// We'll use the first 8 bytes of decompressed array to store these temporarily
		job.decompressed = new byte[8 + 65536]; // 8 bytes for metadata + max block size
		ByteBuffer meta = ByteBuffer.wrap(job.decompressed).order(ByteOrder.LITTLE_ENDIAN);
		meta.putInt((int)expectedCrc);
		meta.putInt(expectedSize);

		if (compressedSize == 0 && expectedSize == 0) {
			synchronized(lastJobLock) {
				if (!lastJobQueued) {
					job.lastJob = true;
					lastJobQueued = true;
				}
			}
		}

		assert(!DEBUG || job.repOK()) : "readNextBlock created invalid job";

		return job;
	}

	/**
	 * Find BC subfield in gzip extra field and extract BSIZE.
	 */
	private int findBsizeInExtra(byte[] extra, int xlen) {
		int pos = 0;
		while (pos + 4 <= xlen) {
			int si1 = extra[pos] & 0xFF;
			int si2 = extra[pos + 1] & 0xFF;
			int slen = ((extra[pos + 3] & 0xFF) << 8) | (extra[pos + 2] & 0xFF);

			if (si1 == 66 && si2 == 67) { // 'B' 'C'
				if (slen == 2 && pos + 6 <= xlen) {
					return ((extra[pos + 5] & 0xFF) << 8) | (extra[pos + 4] & 0xFF);
				}
			}

			pos += 4 + slen;
		}
		return -1;
	}

	/**
	 * Worker thread: Decompress BGZF blocks.
	 */
	private void workerLoop() {
		Inflater inflater = new Inflater(true); // true = nowrap mode for raw deflate

		try {
			while (!closed) {
				BgzfJob job = inputQueue.take();

				if (DEBUG) {
					System.err.println("Worker: dequeued job " + job.id +
						(job.isPoisonPill() ? " (POISON)" : ""));
				}

				if (job.isPoisonPill()) {
					if (DEBUG) {
						System.err.println("Worker: got POISON, re-injecting and terminating");
					}
					inputQueue.put(BgzfJob.POISON_PILL);
					break;
				}

				assert job.compressed != null : "Worker received job with null compressed data";
				assert(!DEBUG || job.repOK()) : "Worker received invalid job";

				// Extract expected CRC and size from metadata
				ByteBuffer meta = ByteBuffer.wrap(job.decompressed).order(ByteOrder.LITTLE_ENDIAN);
				long expectedCrc = meta.getInt() & 0xFFFFFFFFL;
				int expectedSize = meta.getInt();

				// Decompress into same array (after metadata)
				inflater.reset();
				inflater.setInput(job.compressed, 0, job.compressedSize);

				if (DEBUG) {
					System.err.println("Worker: inflating job " + job.id);
				}

				try {
					job.decompressedSize = inflater.inflate(job.decompressed, 8, 65536);
					if (DEBUG) {
						System.err.println("Worker: inflated job " + job.id + " to " + job.decompressedSize + " bytes");
					}
				} catch (DataFormatException e) {
					job.error = new IOException("Decompression failed for job " + job.id, e);
					if (!offerJobToHeap(job)) {
						break;
					}
					continue;
				}

				// Validate decompressed size
				assert job.decompressedSize == expectedSize :
					"Size mismatch for job " + job.id + ": expected " + expectedSize + ", got " + job.decompressedSize;
				if (job.decompressedSize != expectedSize) {
					job.error = new IOException("Uncompressed size mismatch for job " + job.id +
						": expected " + expectedSize + ", got " + job.decompressedSize);
					if (!offerJobToHeap(job)) {
						break;
					}
					continue;
				}

				// Verify CRC32
				CRC32 crc = new CRC32();
				crc.update(job.decompressed, 8, job.decompressedSize);
				long actualCrc = crc.getValue();

				assert actualCrc == expectedCrc :
					"CRC32 mismatch for job " + job.id + ": expected " + expectedCrc + ", got " + actualCrc;
				if (actualCrc != expectedCrc) {
					job.error = new IOException("CRC32 mismatch for job " + job.id);
					if (!offerJobToHeap(job)) {
						break;
					}
					continue;
				}

				// Move decompressed data to start of array (remove metadata)
				System.arraycopy(job.decompressed, 8, job.decompressed, 0, job.decompressedSize);

				// Clear compressed data to free memory
				job.compressed = null;
				job.compressedSize = 0;

				assert(!DEBUG || job.repOK()) : "Worker produced invalid job";

				// Add to output heap
				if (!offerJobToHeap(job)) {
					break;
				}

				if (job.lastJob) {
					if (DEBUG) {
						System.err.println("Worker: job " + job.id + " marked LAST, injecting POISON");
					}
					inputQueue.put(BgzfJob.POISON_PILL);
				}
			}
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
		} finally {
			inflater.end();
		}
	}

	/**
	 * Add job to output heap with proper synchronization.
	 */
	private void addToOutputHeap(BgzfJob job) throws InterruptedException {
		synchronized(outputHeap) {
			while (true) {
				if (outputHeap.size() < maxHeapSize ||
					outputHeap.isEmpty() ||
					job.id <= outputHeap.peek().id) {
					outputHeap.add(job);
					if (DEBUG) {
						System.err.println("Worker: added job " + job.id +
							" to heap (heap size now " + outputHeap.size() + ")");
					}
					outputHeap.notifyAll();
					return;
				}
				if (DEBUG) {
					long smallest = outputHeap.isEmpty() ? -1 : outputHeap.peek().id;
					System.err.println("Worker: waiting to add job " + job.id +
						" (heap size " + outputHeap.size() + ", smallest " + smallest + ")");
				}
				outputHeap.wait();
			}
		}
	}

	/**
	 * Helper that adds a job to the heap, returning false if interrupted.
	 */
	private boolean offerJobToHeap(BgzfJob job) {
		try {
			addToOutputHeap(job);
			return true;
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			return false;
		}
	}

	@Override
	public int read() throws IOException {
		byte[] b = new byte[1];
		int n = read(b, 0, 1);
		return n < 0 ? -1 : (b[0] & 0xFF);
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException {
		assert b != null : "Null buffer";
		assert off >= 0 && len >= 0 && len <= b.length - off :
			"Invalid offset/length: off=" + off + ", len=" + len + ", buf.length=" + b.length;

		if (closed) {
			throw new IOException("Stream closed");
		}
		if (len == 0) {
			return 0;
		}

		if (eofReached) {
			return -1;
		}

		// Check for worker errors
		if (workerError != null) {
			throw workerError;
		}

		int totalRead = 0;
		while (totalRead < len) {
			// Need new block?
			if (currentBlockPos >= currentBlockSize) {
				BgzfJob nextJob = waitForNextJob();
				if (nextJob == null) {
					// EOF
					eofReached = true;
					return totalRead == 0 ? -1 : totalRead;
				}

				if (DEBUG) {
					System.err.println("Consumer: received job " + nextJob.id +
						", decompressedSize=" + nextJob.decompressedSize +
						", decompressed=" + (nextJob.decompressed != null ? "not null" : "NULL"));
				}

				// Check for job error
				if (nextJob.error != null) {
					throw new IOException("Decompression failed", nextJob.error);
				}

				assert nextJob.decompressed != null : "Job " + nextJob.id + " has null decompressed data";

				// EOF marker blocks have 0 decompressed bytes - this is normal
				if (nextJob.decompressedSize == 0 || nextJob.lastJob) {
					if (DEBUG) {
						System.err.println("Consumer: job " + nextJob.id + " has 0 decompressed bytes, treating as EOF");
					}
					// EOF block encountered
					nextExpectedId++;
					eofReached = true;
					return totalRead == 0 ? -1 : totalRead;
				}

				assert nextJob.decompressedSize > 0 : "Job " + nextJob.id + " has zero decompressed size";

				currentBlock = nextJob.decompressed;
				currentBlockSize = nextJob.decompressedSize;
				currentBlockPos = 0;
				nextExpectedId++;
			}

			// Copy from current block
			int available = currentBlockSize - currentBlockPos;
			int toCopy = Math.min(available, len - totalRead);

			assert toCopy > 0 : "toCopy should be positive: " + toCopy;
			assert currentBlockPos + toCopy <= currentBlockSize :
				"Copy would exceed block: pos=" + currentBlockPos + ", toCopy=" + toCopy + ", size=" + currentBlockSize;

			System.arraycopy(currentBlock, currentBlockPos, b, off + totalRead, toCopy);
			currentBlockPos += toCopy;
			totalRead += toCopy;
		}

		assert totalRead == len : "Read wrong amount: expected " + len + ", got " + totalRead;
		return totalRead;
	}

	/**
	 * Wait for next sequential job from output heap.
	 * @return Next job, or null on EOF
	 */
	private BgzfJob waitForNextJob() throws IOException {
		synchronized(outputHeap) {
			while (true) {
				// Check for errors
				if (workerError != null) {
					throw workerError;
				}

				// Check if next job is available
				if (!outputHeap.isEmpty() && outputHeap.peek().id == nextExpectedId) {
					BgzfJob job = outputHeap.poll();
					assert job != null : "poll() returned null when peek() was non-null";
					assert job.id == nextExpectedId :
						"Wrong job ID: expected " + nextExpectedId + ", got " + job.id;
					if (DEBUG) {
						System.err.println("Consumer: taking job " + job.id +
							" (heap size now " + outputHeap.size() + ")");
					}
					outputHeap.notifyAll();
					return job;
				}

				// Check if we're done
				if (producerFinished) {
					// All jobs submitted, wait for workers to finish processing
					boolean allWorkersFinished = true;
					for (Thread worker : workers) {
						if (worker.isAlive()) {
							allWorkersFinished = false;
							break;
						}
					}

					if (allWorkersFinished && outputHeap.isEmpty()) {
						if (DEBUG) {
							System.err.println("Consumer: no more jobs, returning EOF");
						}
						return null; // EOF
					}
				}

				if (DEBUG) {
					long smallest = outputHeap.isEmpty() ? -1 : outputHeap.peek().id;
					System.err.println("Consumer: waiting for job " + nextExpectedId +
						" (heap size " + outputHeap.size() + ", smallest " + smallest + ")");
				}
				// Wait for next job
				try {
					outputHeap.wait();
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
					throw new IOException("Interrupted while waiting for data");
				}
			}
		}
	}

	/**
	 * Read exactly n bytes from input stream.
	 * @return number of bytes read (0 on immediate EOF, n on success)
	 */
	private int readFully(byte[] b, int off, int len) throws IOException {
		int total = 0;
		while (total < len) {
			int n = in.read(b, off + total, len - total);
			if (n < 0) {
				return total;
			}
			total += n;
		}
		return total;
	}

	@Override
	public void close() throws IOException {
		if (closed) {
			return;
		}

		closed = true;

		// Interrupt threads
		if (producer != null) {
			producer.interrupt();
		}
		if (workers != null) {
			for (Thread worker : workers) {
				worker.interrupt();
			}
		}

		// Wait for threads to finish
		try {
			if (producer != null) {
				producer.join(1000);
			}
			if (workers != null) {
				for (Thread worker : workers) {
					worker.join(1000);
				}
			}
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
		}

		// Close underlying stream
		in.close();
	}

	/**
	 * Validate internal state for debugging.
	 */
	private boolean repOK() {
		if (in == null) return false;
		if (workerThreads <= 0 || workerThreads > 32) return false;
		if (inputQueue == null || outputHeap == null) return false;
		if (maxHeapSize < 1) return false;
		if (outputHeap.size() > maxHeapSize) return false;
		if (nextExpectedId < 0 || nextJobId < 0) return false;
		if (nextExpectedId > nextJobId) return false; // Can't expect a job that hasn't been created yet
		if (currentBlockPos < 0 || currentBlockSize < 0) return false;
		if (currentBlockPos > currentBlockSize) return false;
		return true;
	}
}
