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
 * @author Brian Bushnell
 * @date Aug 23, 2010
 *
 */
public class TextStreamWriter extends Thread {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a TextStreamWriter with default format (TEXT).
	 *
	 * @param fname_ Output file name
	 * @param overwrite_ Whether to overwrite existing files
	 * @param append_ Whether to append to existing files
	 * @param allowSubprocess_ Whether to allow subprocess execution
	 */
	public TextStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_){
		this(fname_, overwrite_, append_, allowSubprocess_, 0);
	}
	
	/**
	 * Creates a TextStreamWriter with specified format.
	 *
	 * @param fname_ Output file name
	 * @param overwrite_ Whether to overwrite existing files
	 * @param append_ Whether to append to existing files
	 * @param allowSubprocess_ Whether to allow subprocess execution
	 * @param format Output format type
	 */
	public TextStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_, int format){
		this(FileFormat.testOutput(fname_, FileFormat.TEXT, format, 0, allowSubprocess_, overwrite_, append_, false));
	}
	
	/**
	 * Creates a TextStreamWriter from a FileFormat configuration.
	 * Initializes all format flags, output streams, and buffer settings based on
	 * the provided FileFormat. Handles BAM files through samtools subprocess if
	 * available, otherwise uses standard output streams.
	 *
	 * @param ff FileFormat containing output configuration and format settings
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
	
	
	@Override
	public void run() {
		if(verbose){System.err.println("running");}
		assert(open) : fname;
		
		synchronized(this){
			started=true;
			this.notify();
		}
		
		ArrayList<CharSequence> job=null;

		if(verbose){System.err.println("waiting for jobs");}
		while(job==null){
			try {
				job=queue.take();
				if(verbose){System.err.println("grabbed first job of size "+job.size());}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(verbose){System.err.println("processing jobs");}
		while(job!=null && job!=POISON2){
			if(!job.isEmpty()){
//				if(verbose){System.err.println("writing job of size "+job.size());}
				for(final CharSequence cs : job){
//					if(verbose){System.err.println("writing cs of size "+cs.length());}
					assert(cs!=POISON);
					myWriter.print(cs);
//					if(verbose){System.err.println("printing "+cs);}
				}
			}
			
			job=null;
			while(job==null){
				try {
					job=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		if(verbose){System.err.println("null/poison job");}
//		assert(false);
		open=false;
		if(verbose){System.err.println("call finish writing");}
		ReadWrite.finishWriting(myWriter, myOutstream, fname, allowSubprocess);
		if(verbose){System.err.println("finished writing");}
		synchronized(this){notifyAll();}
		if(verbose){System.err.println("done");}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Control and Helpers     ----------------*/
	/*--------------------------------------------------------------*/
	
	
	@Override
	public synchronized void start(){
		super.start();
		if(verbose){System.err.println(this.getState());}
		synchronized(this){
			while(!started){
				try {
					this.wait(20);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}

	
	/**
	 * Signals the writer thread to terminate gracefully.
	 * Waits for thread to start if not already running, adds any remaining
	 * buffer contents to the queue, and sends a poison job to terminate
	 * the writing loop. Sets the writer to closed state.
	 */
	public synchronized void poison(){
		//Don't allow thread to shut down before it has started
		while(!started || this.getState()==Thread.State.NEW){
			if(verbose){System.err.println("waiting for start.");}
			try {
				this.wait(20);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(verbose){System.err.println("testing if open.");}
		if(!open){return;}
//		if(verbose){System.err.println("adding buffer: "+buffer.size());}
		addJob(buffer);
		buffer=null;
//		System.err.println("Poisoned!");
//		assert(false);
		
//		assert(false) : open+", "+this.getState()+", "+started;
		open=false;
		addJob(POISON2);
	}
	
	/**
	 * Blocks until the writer thread has completely terminated.
	 * Repeatedly attempts to join the thread with 1-second timeouts
	 * until the thread state becomes TERMINATED.
	 */
	public void waitForFinish(){
		if(verbose){System.err.println("waiting for finish.");}
		while(this.getState()!=Thread.State.TERMINATED){
			if(verbose){System.err.println("attempting join.");}
			try {
				this.join(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * @return true if there was an error, false otherwise
	 */
	public boolean poisonAndWait(){
		poison();
		waitForFinish();
		assert(buffer==null || buffer.isEmpty());
		return errorState;
	}
	
	//TODO Why is this synchronized?
	/**
	 * Adds a writing job to the queue for processing by the writer thread.
	 * Blocks until the job is successfully added to the ArrayBlockingQueue.
	 * Requires the thread to be started before accepting jobs.
	 * @param j ArrayList of CharSequences to write
	 */
	public synchronized void addJob(ArrayList<CharSequence> j){
		if(verbose){System.err.println("Got job "+(j==null ? "null" : j.size()));}
		
		assert(started) : "Wait for start() to return before using the writer.";
//		while(!started || this.getState()==Thread.State.NEW){
//			try {
//				this.wait(20);
//			} catch (InterruptedException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
		
		boolean success=false;
		while(!success){
			try {
				queue.put(j);
				success=true;
			} catch (InterruptedException e) {
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
	

	/**
	 * Buffers a CharSequence for writing to the output stream.
	 * Adds the sequence to the current buffer and flushes the buffer
	 * if it reaches the maximum size or length threshold.
	 * @param cs CharSequence to write (null converted to "null")
	 */
	public void print(CharSequence cs){
		if(cs==null){cs="null";}
//		System.err.println("Added line '"+cs+"'");
//		System.err.println("Adding "+cs.length()+" chars.");
		assert(open) : cs;
		buffer.add(cs);
		bufferLen+=cs.length();
		if(buffer.size()>=buffersize || bufferLen>=maxBufferLen){
			addJob(buffer);
			buffer=new ArrayList<CharSequence>(buffersize);
			bufferLen=0;
		}
	}
	

	public void print(byte[] cs){
		if(cs==null){cs="null".getBytes();}
//		System.err.println("Added line '"+cs+"'");
//		System.err.println("Adding "+cs.length()+" chars.");
		assert(open) : cs;
		buffer.add(new ByteBuilder(cs));
		bufferLen+=cs.length;
		if(buffer.size()>=buffersize || bufferLen>=maxBufferLen){
			addJob(buffer);
			buffer=new ArrayList<CharSequence>(buffersize);
			bufferLen=0;
		}
	}
	
	/** Prints a long number by converting to string.
	 * @param number The number to print */
	public void print(long number){
		print(Long.toString(number));
	}
	
	/**
	 * Prints a Read object in the appropriate format.
	 * Converts the Read to the configured output format (FASTQ, FASTA, SAM, etc.)
	 * and buffers it for writing. Cannot be used with OTHER format type.
	 * @param r The Read object to print
	 */
	public void print(Read r){
		assert(!OTHER);
		ByteBuilder sb=(FASTQ ? r.toFastq() : FASTA ? r.toFasta(FASTA_WRAP) : SAM ? r.toSam() :
			SITES ? r.toSites() : INFO ? r.toInfo() : r.toText(true));
		print(sb);
	}
	
	/**
	 * Writes CharSequences in order based on key values.
	 * Stores sequences in a HashMap until all preceding keys are available,
	 * then writes them in sequential order. Ensures ordered output even
	 * when sequences arrive out of order.
	 *
	 * @param cs CharSequence to write
	 * @param key Ordering key (must be >= nextKey)
	 */
	public synchronized void writeOrdered(CharSequence cs, long key){
		assert(cs!=null) : key;
		assert(key>=nextKey) : key+", "+nextKey;
		assert(!map.containsKey(key));

		map.put(key, cs);

		while(map.containsKey(nextKey)){
			//		System.err.println("Writing list "+first.get(0).numericID);
			CharSequence value=map.remove(nextKey);
			print(value);
			nextKey++;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Println            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Prints a newline character */
	public void println(){
		print("\n");
	}
	
	/** Prints a CharSequence followed by a newline.
	 * @param cs CharSequence to print */
	public void println(CharSequence cs){
		print(cs);
		print("\n");
	}
	
	public void println(byte[] cs){
		print(cs);
		print("\n");
	}
	
	/**
	 * Prints a Read object followed by a newline.
	 * Converts the Read to the appropriate format and appends a newline.
	 * @param r The Read object to print
	 */
	public void println(Read r){
		assert(!OTHER);
		ByteBuilder sb=(FASTQ ? r.toFastq() : FASTA ? r.toFasta(FASTA_WRAP) : SAM ? r.toSam() :
			SITES ? r.toSites() : INFO ? r.toInfo() : r.toText(true)).append('\n');
		print(sb);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<CharSequence> buffer;
	
	/** Number of CharSequences to buffer before flushing (default 100) */
	public int buffersize=100;
	/** Maximum total character length before flushing buffer (default 60000) */
	public int maxBufferLen=60000;
	private int bufferLen=0;
	/** Whether this writer will overwrite existing files */
	public final boolean overwrite;
	/** Whether this writer will append to existing files */
	public final boolean append;
	/** Whether this writer allows subprocess execution for compression */
	public final boolean allowSubprocess;
	/** Output file name */
	public final String fname;
	private final OutputStream myOutstream;
	private final PrintWriter myWriter;
	private final ArrayBlockingQueue<ArrayList<CharSequence>> queue;
	private boolean open=true;
	private volatile boolean started=false;
	
	/** TODO */
	public boolean errorState=false;
	
	private HashMap<Long, CharSequence> map=new HashMap<Long, CharSequence>();
	private long nextKey=0;
	
	/*--------------------------------------------------------------*/
	
	private final boolean BAM;
	private final boolean SAM;
	private final boolean FASTQ;
	private final boolean FASTA;
	private final boolean BREAD;
	private final boolean SITES;
	private final boolean INFO;
	private final boolean OTHER;
	
	private final int FASTA_WRAP=Shared.FASTA_WRAP;
	
	/*--------------------------------------------------------------*/

	private static final String POISON=new String("POISON_TextStreamWriter");
	private static final ArrayList<CharSequence> POISON2=new ArrayList<CharSequence>(1);
	
	/** Global verbose output flag for debugging */
	public static boolean verbose=false;
	
}
