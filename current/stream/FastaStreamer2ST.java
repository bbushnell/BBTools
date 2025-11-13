package stream;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.ByteFile1Fc;
import fileIO.FileFormat;
import shared.KillSwitch;
import shared.Shared;
import structures.IntList;
import structures.ListNum;

/**
 * Single-threaded FASTA file loader using ByteFile1Fc (pre-stripped records).
 * Not suitable for interleaved files.
 * Runs in its own thread with output queue.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 12, 2025
 */
public class FastaStreamer2ST implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public FastaStreamer2ST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}

	/** Constructor. */
	public FastaStreamer2ST(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(!interleaved) : "FastaStreamer2ST does not support interleaved files";
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		// Simple output queue
		outputQueue=new ArrayBlockingQueue<ListNum<Read>>(QUEUE_SIZE);

		if(verbose){outstream.println("Made FastaStreamer2ST");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void start(){
		if(verbose){outstream.println("FastaStreamer2ST.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Start processing thread
		thread=new ProcessThread();
		thread.start();

		if(verbose){outstream.println("FastaStreamer2ST started.");}
	}

	@Override
	public void close(){
		if(bf!=null) {bf.close(); bf=null;}
	}

	@Override
	public String fname(){return fname;}

	@Override
	public boolean hasMore(){
		return !finished;
	}

	@Override
	public boolean errorState(){return errorState;}

	@Override
	public boolean paired(){return false;}

	@Override
	public int pairnum(){return pairnum;}

	@Override
	public long readsProcessed(){return readsProcessed;}

	@Override
	public long basesProcessed(){return basesProcessed;}

	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	@Override
	public ListNum<Read> nextList(){
		try{
			ListNum<Read> list=outputQueue.take();
			if(verbose){
				if(list==null || list.last()) {outstream.println("Consumer got terminal list.");}
				else {outstream.println("Consumer got list "+list.id());}
			}
			if(list!=null && list.last()){
				finished=true;
				readsProcessed=thread.readsProcessedT;
				basesProcessed=thread.basesProcessedT;
				errorState=!thread.success;
				return null;
			}
			return list;
		}catch(InterruptedException e){
			errorState=true;
			return null;
		}
	}

	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTA does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private class ProcessThread extends Thread {

		/** Constructor */
		ProcessThread(){
			setName("FastaStreamer2ST-Worker");
		}

		/** Called by start() */
		@Override
		public void run(){
			try{
				processSingle();
				success=true;
			}catch(Exception e){
				e.printStackTrace();
				errorState=true;
			}finally{
				// Send terminal list
				try{
					ListNum<Read> terminal=new ListNum<Read>(null, -1, ListNum.LAST);
					outputQueue.put(terminal);
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("ProcessThread terminated.");}
		}

		void processSingle() throws InterruptedException{
			if(verbose){outstream.println("Started processSingle.");}

			bf=new ByteFile1Fc(ffin);
			IntList newlines=new IntList(256);

			long listNumber=0;

			while(readsProcessedT<maxReads){
				// Get block of records with newline positions
				byte[] block=bf.nextLine(newlines);
				if(block==null || block.length==0){break;}

				ArrayList<Read> readList=new ArrayList<Read>();
				ListNum<Read> reads=new ListNum<Read>(readList, listNumber++);
				reads.firstRecordNum=readsProcessedT;

				for(int i=0, nl0=-1; i<newlines.size() && readsProcessedT<maxReads; i++){
					int nl1=newlines.get(i);
					int nl2=(newlines.size()>i+1 ? newlines.get(i+1) : nl1);
					assert(block[nl0+1]=='>') : nl0+", "+(char)block[nl0+1];
					final byte[] header=KillSwitch.copyOfRange(block, nl0+2, nl1);
					final byte[] bases=(nl2>nl1 ? KillSwitch.copyOfRange(block, nl1+1, nl2) : null);
					
					if(samplerate>=1f || randy.nextFloat()<samplerate){
						Read r=new Read(bases, null, new String(header, StandardCharsets.US_ASCII), readsProcessedT);
						r.setPairnum(pairnum);
						readList.add(r);
						readsProcessedT++;
						basesProcessedT+=r.length();
					}
					
					if(bases!=null) {
						i++;
						nl0=nl2;
					}else {
						nl0=nl1;
					}
				}

				if(readList.size()>0){
					outputQueue.put(reads);
				}
			}

			bf.close();
			if(verbose){outstream.println("Finished processSingle.");}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** True only if this thread has completed successfully */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	public final String fname;

	/** Primary input file */
	final FileFormat ffin;

	/** Output queue */
	final ArrayBlockingQueue<ListNum<Read>> outputQueue;

	/** Processing thread */
	private ProcessThread thread;
	
	/** Input source */
	private ByteFile1Fc bf;

	final int pairnum;
	final boolean interleaved;

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads */
	final long maxReads;

	/** Set when terminal list is received */
	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static final int QUEUE_SIZE=2;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	private float samplerate=1f;
	private java.util.Random randy=null;

}
