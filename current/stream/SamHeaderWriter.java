package stream;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.KillSwitch;
import structures.ByteBuilder;
import structures.ListNum;
import structures.StringNum;

/** Asynchronous SAM header writer consuming sequence batches in numeric ID order.
 * Submit each contiguous batch ID once, starting at zero, then finish all producers
 * before poisoning. Entries retain their order within a batch. Submitted lists and
 * entries are shared, not copied; do not mutate them after submission.
 * Construction starts both the byte writer and header consumer. Header settings
 * are read during execution, not captured as one constructor snapshot. Sequence
 * names/lengths are appended without validation or escaping. A finalization error
 * reported by the byte writer terminates the JVM immediately. An interrupted wait
 * can still return before completion; this class does not validate all output failures.
 * @author Isla, Shinobu
 * @date October 30, 2025
 */
public class SamHeaderWriter{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens output and starts both writer threads with queue capacity 128.
	 * @param ff Nonnull usable output descriptor for the header file
	 */
	public SamHeaderWriter(final FileFormat ff){
		this(ff, 128);
	}

	/** Opens output, creates the ordered queue and starts the header consumer.
	 * There is no separate start call. Output initialization precedes queue creation.
	 * @param ff Nonnull usable output descriptor for the header file
	 * @param queueSize JobQueue backpressure capacity, greater than one; not a strict pending-job maximum
	 */
	public SamHeaderWriter(final FileFormat ff, final int queueSize){
		bsw=ByteStreamWriter.makeBSW(ff);
		queue=new JobQueue<ListNum<StringNum>>(queueSize);
		writerThread=new WriterThread();
		writerThread.start();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Submits an unchanged batch reference to JobQueue, which may block for admission.
	 * IDs determine output order, not call order; supply contiguous IDs from zero,
	 * including empty batches when needed. Does not reject calls after poison locally;
	 * callers must finish all submissions before poisoning and must not mutate data.
	 * @param ln Batch with StringNum.s sequence names and StringNum.n sequence lengths
	 */
	public void add(final ListNum<StringNum> ln){
		queue.add(ln);
	}

	/** Posts one LAST marker at maxSeen()+1; repeated calls do not post another.
	 * Call only after all producers have finished submitting. This method may block
	 * on queue admission; it does not cancel producers or wait for writer completion.
	 */
	public synchronized void poison(){
		if(!closed){
			ListNum<StringNum> last=new ListNum<StringNum>(null, queue.maxSeen()+1, ListNum.LAST);
			queue.add(last);
			closed=true;
		}
	}

	/** Joins the header thread once; poison after producers finish before calling this.
	 * Does not itself poison. Interruption prints a trace and returns without retrying
	 * or restoring interrupt status, so return alone does not guarantee completion.
	 * The worker aborts the JVM if the byte writer reports a finalization error.
	 */
	public synchronized void waitForFinish(){
		try{
			writerThread.join();
		}catch(InterruptedException e){
			e.printStackTrace();
		}
	}

	/** Calls poison then waitForFinish, with their submission/interruption requirements.
	 * A reported finalization error causes the worker to abort the JVM.
	 */
	public synchronized void poisonAndWait(){
		poison();
		waitForFinish();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes configured header text and ID-ordered sequence entries through ByteStreamWriter.
	 * Read-group output is conditional; program output comes from SamHeader.header2B.
	 */
	private class WriterThread extends Thread{

		/** Appends header helpers and unvalidated sequence entries, then finishes the byte writer.
		 * Enqueues text after an SQ entry reaches the 16384-byte threshold, then enqueues
		 * remaining text at completion. This does not itself flush the underlying stream.
		 * Does not catch worker exceptions. A reported finalization error terminates the
		 * JVM through KillSwitch without an additional cleanup phase.
		 */
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
			//STR-027: callers use this void completion API before reporting success. Reject the
			//byte writer's final error flag so failed flush/close cannot look like successful completion.
			if(bsw.poisonAndWait()){KillSwitch.kill("Error writing SAM header to "+bsw.fname);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Output writer created and started before the header consumer. */
	private final ByteStreamWriter bsw;
	/** Queue consuming contiguous numeric batch IDs starting at zero. */
	private final JobQueue<ListNum<StringNum>> queue;
	/** Header-building consumer, started during construction. */
	private final WriterThread writerThread;
	/** Whether this wrapper has submitted its LAST marker; not proof of output completion. */
	private boolean closed=false;
}
