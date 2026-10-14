package fileIO;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;

import dna.Data;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/**
 * Queues text for one background PrintWriter, opening the output during construction.
 * Start the worker before submitting output. Ordinary print/println calls use an
 * unsynchronized buffer and require one producer. Ordered submissions are synchronized,
 * but must not overlap ordinary print calls or shutdown.
 *
 * <p>Submitted CharSequences, lists and wrapped byte arrays are not copied. Their
 * contents must remain unchanged until the writer has finished. The worker converts
 * each CharSequence to text, then PrintWriter encodes it. Plain file output uses
 * the platform default charset.
 * Byte arrays are wrapped in ByteBuilder, whose toString decodes US-ASCII; this is
 * not a byte-preserving output API. Read formatting occurs in the producer.
 *
 * <p>The queue holds at most five batches. The ordinary buffer flushes after reaching
 * either its entry or length threshold; neither threshold bounds one large entry.
 * Ordered output retains an unbounded map of entries waiting for contiguous keys
 * starting at zero. Supply empty text for a key with no output.
 *
 * <p>Finish all producers, then call poisonAndWait and inspect its error result.
 * A missing ordered key currently throws before queuing termination (FIO-035).
 * Normal shutdown collects retained PrintWriter, stream-close and subprocess errors;
 * it does not provide recovery from arbitrary worker exceptions.
 *
 * @author Brian Bushnell
 * @date Aug 23, 2010
 */
public class TextStreamWriter extends Thread{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves the output with TEXT as fallback and no explicit format code.
	 * @param fname_ Output path
	 * @param overwrite_ Whether existing output may be replaced
	 * @param append_ Whether output is appended; mutually exclusive with overwrite
	 * @param allowSubprocess_ Policy passed to output selection and final closure */
	public TextStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_){
		this(fname_, overwrite_, append_, allowSubprocess_, 0);
	}

	/** Resolves FileFormat and delegates construction without starting the worker.
	 * @param fname_ Output path
	 * @param overwrite_ Whether existing output may be replaced
	 * @param append_ Whether output is appended
	 * @param allowSubprocess_ Policy passed to output selection and final closure
	 * @param format Format code passed to FileFormat.testOutput */
	public TextStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_, int format){
		this(FileFormat.testOutput(fname_, FileFormat.TEXT, format, 0, allowSubprocess_, overwrite_, append_, false));
	}

	/**
	 * Opens the output and creates the queue/buffer without starting this thread.
	 * Generic TEXT selects FASTQ when printing Read objects; plain text is unaffected.
	 * Supported BAM output uses ReadWrite's BAM stream, otherwise the generic stream
	 * path is used. Compression and subprocess behavior are delegated to ReadWrite.
	 * @param ff Non-null output configuration; append is allowed only for raw/gzip
	 */
	public TextStreamWriter(FileFormat ff){
		FASTQ=ff.fastq() || ff.text();
		FASTA=ff.fasta();
		BREAD=ff.bread();
		SAM=ff.samOrBam();
		BAM=ff.bam();
		SITES=ff.sites();
		INFO=ff.attachment();
		OTHER=(!FASTQ && !FASTA && !BREAD && !SAM && !BAM && !SITES && !INFO);

		fname=ff.name();
		overwrite=ff.overwrite();
		append=ff.append();
		allowSubprocess=ff.allowSubprocess();
		assert(!(overwrite&append));
		assert(ff.canWrite()) : "File "+fname+" exists and overwrite=="+overwrite;
		if(append && !(ff.raw() || ff.gzip())){throw new RuntimeException("Can't append to compressed files.");}

		if(!BAM || !Data.BAM_SUPPORT_OUT()){
			myOutstream=ReadWrite.getOutputStream(fname, append, true, allowSubprocess);
		}else{
			myOutstream=ReadWrite.getBamOutputStream(fname, append);
		}
		myWriter=new PrintWriter(myOutstream);
		if(verbose){System.err.println("Created PrintWriter for "+myOutstream);}

		queue=new ArrayBlockingQueue<ArrayList<CharSequence>>(5);
		buffer=new ArrayList<CharSequence>(buffersize);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Primary Method        ----------------*/
	/*--------------------------------------------------------------*/

	/** Drains batches until the identity terminal, then closes output and records errors. */
	@Override
	public void run(){
		if(verbose){System.err.println("running");}
		assert(open) : fname;

		synchronized(this){
			started=true;
			this.notify();
		}

		ArrayList<CharSequence> job=null;

		if(verbose){System.err.println("waiting for jobs");}
		while(job==null){
			try{
				job=queue.take();
				if(verbose){System.err.println("grabbed first job of size "+job.size());}
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		if(verbose){System.err.println("processing jobs");}
		while(job!=null && job!=POISON2){
			if(!job.isEmpty()){
				for(final CharSequence cs : job){
					assert(cs!=POISON);
					myWriter.print(cs);
				}
			}

			job=null;
			while(job==null){
				try{
					job=queue.take();
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		if(verbose){System.err.println("null/poison job");}
		open=false;
		if(verbose){System.err.println("call finish writing");}
		errorState=ReadWrite.finishWriting(myWriter, myOutstream, fname, allowSubprocess)|errorState;
		if(verbose){System.err.println("finished writing");}
		synchronized(this){notifyAll();}
		if(verbose){System.err.println("done");}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Control and Helpers     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Starts the worker and waits for its startup handshake. No buffered text may
	 * precede this call. Interrupted waits are logged and retried without restoring
	 * the interrupt status.
	 */
	@Override
	public synchronized void start(){
		//ADDED [TextStreamWriter start-guard]: catch the write-before-start bug class. Buffered print/println before start()
		//force an addJob that asserts(started) (-> -ea crash once buffered), or are lost if the writer is never started.
		//TextStreamWriter has NO direct-write path (no forcePrint), so it must ALWAYS be started before any use.
		//Motivated by tax/ImgRecord#001 (never started/closed) + tax/AnalyzeAccession#001 (started after the writes).
		assert(buffer==null || buffer.isEmpty()) : "TextStreamWriter.start() called after "+buffer.size()+" entr(ies) were already buffered via print/println; start() must precede all writes.";
		super.start();
		if(verbose){System.err.println(this.getState());}
		synchronized(this){
			while(!started){
				try{
					this.wait(20);
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Queues remaining buffered text and the terminal after all producers finish.
	 * Requires a separately started worker; waits for startup but does not start it.
	 * Normal return requests shutdown without waiting for output closure.
	 * @throws RuntimeException If ordered text remains behind a missing key; currently
	 * this happens before shutdown is queued and can leave the worker waiting (FIO-035)
	 */
	public synchronized void poison(){
		//Don't allow thread to shut down before it has started
		while(!started || this.getState()==Thread.State.NEW){
			if(verbose){System.err.println("waiting for start.");}
			try{
				this.wait(20);
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		if(verbose){System.err.println("testing if open.");}
		if(!open){return;}
		//Fail loud rather than silently drop: a non-empty ordered map means a writeOrdered key was skipped.
		//TODO: Probable bug - FIO-035: this exception precedes queue termination, leaving
		//the worker waiting for a job and never reaching normal output closure.
		if(!map.isEmpty()){
			throw new RuntimeException("TextStreamWriter: ordered output incomplete - key "+nextKey+
				" was never supplied; "+map.size()+" buffered entr(ies) would be lost. "+
				"Every key passed to writeOrdered() must be contiguous from 0 and supplied exactly once.");
		}
		addJob(buffer);
		buffer=null;

		open=false;
		addJob(POISON2);
	}

	/** Joins a started worker after shutdown is requested; does not start or poison it.
	 * Interrupted joins are logged and retried without restoring interrupt status. */
	public void waitForFinish(){
		if(verbose){System.err.println("waiting for finish.");}
		while(this.getState()!=Thread.State.TERMINATED){
			if(verbose){System.err.println("attempting join.");}
			try{
				this.join(1000);
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
	}

	/**
	 * Requests shutdown and joins the worker. All producers must have finished.
	 * A poison exception propagates without joining; normal return observes final status.
	 * @return Retained output/closure/subprocess errors or a previously caller-set error
	 */
	public boolean poisonAndWait(){
		poison();
		waitForFinish();
		assert(buffer==null || buffer.isEmpty());
		return errorState;
	}

	//Synchronized to share this monitor with writeOrdered()/poison() and to gate on the started handshake; queue.put() is itself thread-safe, so the lock is belt-and-suspenders.
	/**
	 * Transfers a batch to the bounded queue, waiting for capacity as needed.
	 * Requires completed startup. Does not flush pending print-buffer entries, so
	 * callers must coordinate ordering when mixing entry points. No payload is copied.
	 * Interrupted queue waits are logged and retried.
	 * @param j Non-null list whose contents and mutable entries must no longer change
	 */
	public synchronized void addJob(ArrayList<CharSequence> j){
		if(verbose){System.err.println("Got job "+(j==null ? "null" : j.size()));}

		assert(started) : "Wait for start() to return before using the writer.";

		boolean success=false;
		while(!success){
			try{
				queue.put(j);
				success=true;
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
				assert(!queue.contains(j)); //Hopefully it was not added.
			}
		}
		if(verbose){System.err.println("Put job in queue: "+success);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Print             ----------------*/
	/*--------------------------------------------------------------*/

	/** Buffers text without a newline; null becomes the literal text "null".
	 * Single-producer operation. Retains the supplied object until worker conversion.
	 * @param cs Text whose contents must remain unchanged until writing finishes */
	public void print(CharSequence cs){
		if(cs==null){cs="null";}
		assert(open) : cs;
		buffer.add(cs);
		bufferLen+=cs.length();
		if(buffer.size()>=buffersize || bufferLen>=maxBufferLen){
			addJob(buffer);
			buffer=new ArrayList<CharSequence>(buffersize);
			bufferLen=0;
		}
	}

	/** Buffers a ByteBuilder wrapping the array without copying or adding a newline.
	 * Null uses "null".getBytes(); ByteBuilder later decodes the bytes as US-ASCII.
	 * @param cs Bytes to retain unchanged until writing finishes */
	public void print(byte[] cs){
		if(cs==null){cs="null".getBytes();}
		assert(open) : cs;
		buffer.add(new ByteBuilder(cs));
		bufferLen+=cs.length;
		if(buffer.size()>=buffersize || bufferLen>=maxBufferLen){
			addJob(buffer);
			buffer=new ArrayList<CharSequence>(buffersize);
			bufferLen=0;
		}
	}

	/** Buffers a decimal long without a newline. */
	public void print(long number){
		print(Long.toString(number));
	}

	/**
	 * Formats one Read using the selected format and buffers the resulting builder.
	 * Does not append a newline or explicitly serialize its mate. INFO can retain
	 * the Read's mutable ByteBuilder attachment directly; do not mutate it afterward.
	 * @param r Non-null Read with fields required by its selected formatter
	 */
	public void print(Read r){
		assert(!OTHER);
		ByteBuilder sb=(FASTQ ? r.toFastq() : FASTA ? r.toFasta(FASTA_WRAP) : SAM ? r.toSam() :
			SITES ? r.toSites() : INFO ? r.toInfo() : r.toText(true));
		print(sb);
	}

	/**
	 * Retains out-of-order entries and buffers each available contiguous prefix.
	 * Keys start at zero and must be supplied exactly once, including empty text
	 * for a key with no output. The pending map has no capacity limit. Finish all
	 * producers before shutdown; do not overlap ordinary print calls with this method.
	 * A gap makes poison throw before queue termination, leaving a waiting worker (FIO-035).
	 * @param cs Non-null text retained without copying until worker conversion
	 * @param key Dense sequence number no smaller than the next expected key
	 */
	public synchronized void writeOrdered(CharSequence cs, long key){
		assert(cs!=null) : key;
		assert(key>=nextKey) : key+", "+nextKey;
		assert(!map.containsKey(key));

		map.put(key, cs);

		while(map.containsKey(nextKey)){
			CharSequence value=map.remove(nextKey);
			print(value);
			nextKey++;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Println            ----------------*/
	/*--------------------------------------------------------------*/

	/** Buffers one LF character. */
	public void println(){
		print("\n");
	}

	/** Buffers the text, then a separate LF entry; this is not an atomic multi-producer call. */
	public void println(CharSequence cs){
		print(cs);
		print("\n");
	}

	/** Buffers the wrapped bytes, then LF, with the same ownership/encoding as print(byte[]). */
	public void println(byte[] cs){
		print(cs);
		print("\n");
	}

	/** Formats one Read and appends LF to its builder before submission.
	 * INFO may return the attached ByteBuilder itself, so this call can append LF
	 * to that attachment. Requires the same ownership as print(Read). */
	public void println(Read r){
		assert(!OTHER);
		ByteBuilder sb=(FASTQ ? r.toFastq() : FASTA ? r.toFasta(FASTA_WRAP) : SAM ? r.toSam() :
			SITES ? r.toSites() : INFO ? r.toInfo() : r.toText(true)).append('\n');
		print(sb);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Single-producer batch, transferred to the queue on flush and nulled at shutdown. */
	private ArrayList<CharSequence> buffer;

	/** Entry threshold checked after appending; also the next batch's initial capacity. */
	public int buffersize=100;
	/** Sum-of-lengths threshold checked after appending, not an encoded-byte bound. */
	public int maxBufferLen=60000;
	/** Sum of lengths at submission time for the current ordinary buffer. */
	private int bufferLen=0;
	/** Resolved overwrite setting used during construction. */
	public final boolean overwrite;
	/** Resolved append setting used during construction. */
	public final boolean append;
	/** Policy passed to ReadWrite output selection and final closure. */
	public final boolean allowSubprocess;
	/** Resolved output name and subprocess-registry key. */
	public final String fname;
	/** Output selected during construction; closure is delegated to ReadWrite. */
	private final OutputStream myOutstream;
	/** Worker-owned text encoder, retaining checked write errors for final inspection. */
	private final PrintWriter myWriter;
	/** Five-batch queue; the static POISON2 list is recognized by identity. */
	private final ArrayBlockingQueue<ArrayList<CharSequence>> queue;
	/** Submission/shutdown state, not a thread-safe public liveness API. */
	private boolean open=true;
	/** Worker-start handshake published before start returns. */
	private volatile boolean started=false;

	/** Final output/closure/subprocess errors, ORed with any caller-set error.
	 * Read after joining through poisonAndWait; this is not live per-write status. */
	public boolean errorState=false;

	/** Unbounded pending ordered entries; mutable values are retained without copying. */
	private HashMap<Long, CharSequence> map=new HashMap<Long, CharSequence>();
	/** Next ordered key, initially zero. */
	private long nextKey=0;

	/*--------------------------------------------------------------*/

	/** Requests the BAM output stream when BAM support is available. */
	private final boolean BAM;
	/** Selects SAM text formatting for both SAM and BAM output. */
	private final boolean SAM;
	/** Selects FASTQ Read formatting, including a generic TEXT FileFormat. */
	private final boolean FASTQ;
	/** Selects FASTA Read formatting. */
	private final boolean FASTA;
	/** Native Read-text format; handled by the formatting fallback. */
	private final boolean BREAD;
	/** Selects candidate-site formatting. */
	private final boolean SITES;
	/** Selects Read attachment formatting, which can share a mutable builder. */
	private final boolean INFO;
	/** Unsupported Read formatting; raw text calls remain usable. */
	private final boolean OTHER;

	/** FASTA wrapping sampled from Shared during construction. */
	private final int FASTA_WRAP=Shared.FASTA_WRAP;

	/*--------------------------------------------------------------*/

	/** Legacy assertion-only reserved text marker; not the queue shutdown marker. */
	private static final String POISON=new String("POISON_TextStreamWriter");
	/** Identity-only queue terminal; other empty lists are ordinary batches. */
	private static final ArrayList<CharSequence> POISON2=new ArrayList<CharSequence>(1);

	/** Enables diagnostic logging. */
	public static boolean verbose=false;

}
