package stream;

import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.FileFormat;
import shared.KillSwitch;
import stream.bam.BamReader;
import stream.bam.BamToSamConverter;
import stream.bam.BgzfInputStream;
import structures.ListNum;

/**
 * Loads BAM files rapidly with multiple threads.
 * Thread 0 reads BAM binary format and converts to SAM text.
 * Worker threads convert SAM text byte[] to SamLine objects.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamLineStreamer extends SamStreamer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructor.
	 */
	public BamLineStreamer(String fname_, int threads_, boolean saveHeader_, long maxReads_) {
		this(FileFormat.testInput(fname_, FileFormat.BAM, null, true, false), threads_, saveHeader_, maxReads_);
	}

	/**
	 * Constructor.
	 */
	public BamLineStreamer(FileFormat ffin_, int threads_, boolean saveHeader_, long maxReads_) {
		super(ffin_, threads_, saveHeader_, maxReads_);
		outq = new ArrayBlockingQueue<ListNum<SamLine>>(threads_ + 2);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public ListNum<SamLine> nextLines() {
		ListNum<SamLine> list = null;
		while (list == null) {
			try {
				list = outq.take();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			if (verbose) {
				outstream.println("a. Got list size " + list.size());
			}
		}
		while (list == POISON_LINES) {
			if (verbose) {
				outstream.println("b. Got poison.");
			}
			try {
				outq.put(list);
				list = null;
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		if (verbose) {
			outstream.println("c. done.");
		}
		return list;
	}

	@Override
	public ListNum<Read> nextReads() {
		KillSwitch.kill("Unsupported.");
		return null;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Spawn process threads */
	@Override
	void spawnThreads() {

		// Determine how many threads may be used
		final int threads = this.threads + 1;

		// Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt = new ArrayList<ProcessThread>(threads);
		for (int i = 0; i < threads; i++) {
			alpt.add(new ProcessThread(i, alpt));
		}
		if (verbose) {
			outstream.println("Spawned threads.");
		}

		// Start the threads
		for (ProcessThread pt : alpt) {
			pt.start();
		}
		if (verbose) {
			outstream.println("Started threads.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {

		// Constructor
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_) {
			tid = tid_;
			alpt = (tid == 0 ? alpt_ : null);
		}

		// Called by start()
		@Override
		public void run() {
			// Process the reads
			if (tid == 0) {
				processBamBytes();
			} else {
				makeReads();
			}

			// Indicate successful exit status
			success = true;
		}

		void processBamBytes() {
			try {
				FileInputStream fis = new FileInputStream(fname);
				BgzfInputStream bgzf = new BgzfInputStream(fis);
				BamReader reader = new BamReader(bgzf);

				// Read BAM magic
				byte[] magic = reader.readBytes(4);
				if (!Arrays.equals(magic, new byte[]{'B', 'A', 'M', 1})) {
					throw new RuntimeException("Not a BAM file: " + fname);
				}

				// Read header text
				long l_text = reader.readUint32();
				byte[] text = reader.readBytes((int)l_text);

				// Parse header if requested
				if (saveHeader && header != null) {
					synchronized (header) {
						// Split by newline and add to header
						int start = 0;
						for (int i = 0; i < text.length; i++) {
							if (text[i] == '\n') {
								if (i > start) {
									byte[] line = Arrays.copyOfRange(text, start, i);
									header.add(line);
								}
								start = i + 1;
							}
						}
						// Add last line if not ending with newline
						if (start < text.length) {
							byte[] line = Arrays.copyOfRange(text, start, text.length);
							header.add(line);
						}
						// Set shared header for SamLine processing
						SamReadInputStream.setSharedHeader(header);
					}
				}

				// Read reference sequence dictionary
				int n_ref = reader.readInt32();
				String[] refNames = new String[n_ref];
				for (int i = 0; i < n_ref; i++) {
					long l_name = reader.readUint32();
					refNames[i] = reader.readString((int)l_name - 1); // Exclude NUL
					reader.readUint8(); // Skip NUL terminator
					long l_ref = reader.readUint32(); // Reference length (unused here)
				}

				BamToSamConverter converter = new BamToSamConverter(refNames);

				// Read alignment records
				ArrayList<byte[]> list = new ArrayList<byte[]>(LIST_SIZE);
				long listNumber = 0;

				try {
					while (true) {
						long block_size = reader.readUint32();
						byte[] bamRecord = reader.readBytes((int)block_size);
						byte[] samLine = converter.convertAlignment(bamRecord);

						list.add(samLine);

						if (list.size() >= LIST_SIZE) {
							putBytes(new ListNum<byte[]>(list, listNumber));
							listNumber++;
							list = new ArrayList<byte[]>(LIST_SIZE);
						}
					}
				} catch (EOFException e) {
					// Normal end of file
				}

				if (list.size() > 0) {
					putBytes(new ListNum<byte[]>(list, listNumber));
				}
				putBytes(POISON_BYTES);

				bgzf.close();
				fis.close();

			} catch (IOException e) {
				throw new RuntimeException("Error reading BAM file: " + fname, e);
			}

			success = true;

			// Wait for completion of all threads
			boolean allSuccess = true;
			for (ProcessThread pt : alpt) {

				// Wait until this thread has terminated
				if (pt != this) {
					if (verbose) {
						outstream.println("Waiting for thread " + pt.tid);
					}
					while (pt.getState() != Thread.State.TERMINATED) {
						try {
							pt.join();
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
					}

					// Accumulate per-thread statistics
					readsProcessed += pt.readsProcessedT;
					basesProcessed += pt.basesProcessedT;
					allSuccess &= pt.success;
				}
			}

			putReads(POISON_LINES);
			if (verbose || verbose2) {
				outstream.println("tid " + tid + " done poisoning reads.");
			}

			// Track whether any threads failed
			if (!allSuccess) {
				errorState = true;
			}
			if (verbose || verbose2) {
				outstream.println("tid " + tid + " finished!");
			}
		}

		void putReads(ListNum<SamLine> list) {
			if (verbose) {
				outstream.println("tid " + tid + " putting rlist size " + list.size());
			}
			while (list != null) {
				try {
					outq.put(list);
					list = null;
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			if (verbose) {
				outstream.println("tid " + tid + " done putting rlist");
			}
		}

		/** Iterate through the reads */
		void makeReads() {
			if (verbose) {
				outstream.println("tid " + tid + " started makeReads.");
			}

			ListNum<byte[]> list = takeBytes();
			while (list != POISON_BYTES) {
				ListNum<SamLine> reads = new ListNum<SamLine>(new ArrayList<SamLine>(list.size()), list.id);
				for (byte[] line : list) {
					if (line[0] == '@') {
						// Ignore header lines
					} else {
						SamLine sl = new SamLine(line);
						reads.add(sl);

						readsProcessedT++;
						basesProcessedT += (sl.seq == null ? 0 : sl.length());
					}
				}
				if (reads.size() > 0) {
					putReads(reads);
				}
				list = takeBytes();
			}
			if (verbose || verbose2) {
				outstream.println("tid " + tid + " done making reads.");
			}

			putBytes(POISON_BYTES);
			if (verbose || verbose2) {
				outstream.println("tid " + tid + " done poisoning bytes.");
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT = 0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT = 0;

		/** True only if this thread has completed successfully */
		boolean success = false;

		/** Thread ID */
		final int tid;

		ArrayList<ProcessThread> alpt;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	final ArrayBlockingQueue<ListNum<SamLine>> outq;

}
