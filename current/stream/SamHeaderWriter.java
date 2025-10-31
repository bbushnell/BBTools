package stream;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import structures.ByteBuilder;
import structures.ListNum;
import structures.StringNum;

/**
 * Thread-safe writer for SAM headers using ordered job queue pattern.
 * Accepts sequence names in batches via ListNum, maintains insertion order,
 * and writes properly formatted SAM header with @HD, @SQ, @RG, and @PG lines.
 * 
 * @author Isla
 * @date October 30, 2025
 */
public class SamHeaderWriter {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates header writer with output file configuration.
	 * Initializes ByteStreamWriter, JobQueue, and starts consumer thread.
	 * @param ff Output file format for header file
	 */
	public SamHeaderWriter(FileFormat ff){
		this(ff, 128);
	}

	/**
	 * Creates header writer with specified queue capacity.
	 * @param ff Output file format for header file
	 * @param queueSize Maximum pending jobs before blocking on add()
	 */
	public SamHeaderWriter(FileFormat ff, int queueSize){
		bsw=ByteStreamWriter.makeBSW(ff);
		queue=new JobQueue<ListNum<StringNum>>(queueSize);
		writerThread=new WriterThread();
		writerThread.start();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds batch of sequence entries to header in ordered fashion.
	 * Submits ListNum directly to JobQueue which handles blocking and ordering.
	 * @param ln Batch of sequences where StringNum.s is name and StringNum.n is length
	 */
	public void add(ListNum<StringNum> ln){
		queue.add(ln);
	}

	/**
	 * Poisons the queue to signal no more jobs coming.
	 * Thread-safe, only poisons once even if called multiple times.
	 */
	public synchronized void poison(){
		if(!closed){
			ListNum<StringNum> last=new ListNum<StringNum>(null, queue.maxSeen()+1, ListNum.LAST);
			queue.add(last);
			closed=true;
		}
	}

	/**
	 * Waits for writer thread to finish processing all jobs.
	 * Does not poison the queue - call poison() first.
	 */
	public synchronized void waitForFinish(){
		try {
			writerThread.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Poisons queue and waits for completion.
	 * Convenience method combining poison() and waitForFinish().
	 */
	public synchronized void poisonAndWait(){
		poison();
		waitForFinish();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Consumer thread that reads ordered jobs and writes complete SAM header.
	 * Writes @HD line, @SQ lines for each sequence, then @RG/@PG lines.
	 */
	private class WriterThread extends Thread {

		@Override
		public void run(){
			ByteBuilder bb=new ByteBuilder(4096);

			// Write @HD header line
			SamHeader.header0B(bb);
			bb.nl();

			// Process @SQ sequence dictionary entries
			ListNum<StringNum> ln;
			while((ln=queue.take())!=null){
				if(ln.list!=null){
					for(StringNum sn : ln.list){
						bb.append("@SQ\tSN:");
						bb.append(sn.s);
						bb.append("\tLN:");
						bb.append(sn.n);
						bb.nl();

						if(bb.length()>=16384){
							bsw.addJob(bb);
							bb=new ByteBuilder(4096);
						}
					}
				}
			}

			// Write @RG and @PG lines
			SamHeader.header2B(bb);
			bb.nl();

			// Flush remaining data
			if(!bb.isEmpty()){
				bsw.addJob(bb);
			}

			// Close output stream
			bsw.poisonAndWait();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Output writer for batched text writing */
	private final ByteStreamWriter bsw;
	/** Ordered job queue maintaining insertion sequence */
	private final JobQueue<ListNum<StringNum>> queue;
	/** Consumer thread for writing header lines */
	private final WriterThread writerThread;
	/** Flag indicating queue has been poisoned */
	private boolean closed=false;
}