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
	
	public TextStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_){
		this(fname_, overwrite_, append_, allowSubprocess_, 0);
	}
	
	public TextStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_, int format){
		this(FileFormat.testOutput(fname_, FileFormat.TEXT, format, 0, allowSubprocess_, overwrite_, append_, false));
	}
	
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
		
		if(!BAM || !Data.SAMTOOLS() || !Data.SH()){
			myOutstream=ReadWrite.getOutputStream(fname, append, true, allowSubprocess);
			if(verbose){System.err.println("Created output stream for "+fname+", "+append+", "+true+", "+allowSubprocess);}
		}else{
			myOutstream=ReadWrite.getOutputStreamFromProcess(fname, "samtools view -S -b -h - ", true, append, true, true);
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
	
	public void print(long number){
		print(Long.toString(number));
	}
	
	public void print(Read r){
		assert(!OTHER);
		ByteBuilder sb=(FASTQ ? r.toFastq() : FASTA ? r.toFasta(FASTA_WRAP) : SAM ? r.toSam() :
			SITES ? r.toSites() : INFO ? r.toInfo() : r.toText(true));
		print(sb);
	}
	
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
	
	public void println(){
		print("\n");
	}
	
	public void println(CharSequence cs){
		print(cs);
		print("\n");
	}
	
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
	
	public int buffersize=100;
	public int maxBufferLen=60000;
	private int bufferLen=0;
	public final boolean overwrite;
	public final boolean append;
	public final boolean allowSubprocess;
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
	
	public static boolean verbose=false;
	
}
