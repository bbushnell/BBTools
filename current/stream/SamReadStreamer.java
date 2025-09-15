package stream;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.FileFormat;
import shared.KillSwitch;
import structures.ListNum;

/**
 * Loads sam files rapidly with multiple threads.
 * 
 * @author Brian Bushnell
 * @date November 4, 2016
 *
 */
public class SamReadStreamer extends SamStreamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 */
	public SamReadStreamer(String fname_, int threads_, boolean saveHeader_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false), threads_, saveHeader_, maxReads_);
	}
	
	/**
	 * Constructor.
	 */
	public SamReadStreamer(FileFormat ffin_, int threads_, boolean saveHeader_, long maxReads_){
		super(ffin_, threads_, saveHeader_, maxReads_);
		outq=new ArrayBlockingQueue<ListNum<Read>>(threads+1);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public ListNum<Read> nextReads(){
		ListNum<Read> list=null;
		while(list==null){
			try {
				list=outq.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(verbose){outstream.println("a. Got list size "+list.size());}
		}
		while(list==POISON_READS){
			if(verbose){outstream.println("b. Got poison.");}
			try {
				outq.put(list);
				list=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(verbose){outstream.println("c. done.");}
		return list;
	}
	
	@Override
	public ListNum<SamLine> nextLines(){
		KillSwitch.kill("Unsupported.");
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** Spawn process threads */
	@Override
	void spawnThreads(){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=this.threads+1;
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, alpt));
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_){
			tid=tid_;
			alpt=(tid==0 ? alpt_ : null);
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			if(tid==0){
				processBytes();
			}else{
				makeReads();
			}
			
			//Indicate successful exit status
			success=true;
		}
		
		void processBytes(){
			processBytes0(tid);
			
			success=true;
			
			//Wait for completion of all threads
			boolean allSuccess=true;
			for(ProcessThread pt : alpt){
				
				//Wait until this thread has terminated
				if(pt!=this){
					if(verbose){outstream.println("Waiting for thread "+pt.tid);}
					while(pt.getState()!=Thread.State.TERMINATED){
						try {
							//Attempt a join operation
							pt.join();
						} catch (InterruptedException e) {
							//Potentially handle this, if it is expected to occur
							e.printStackTrace();
						}
					}

					//Accumulate per-thread statistics
					readsProcessed+=pt.readsProcessedT;
					basesProcessed+=pt.basesProcessedT;
					allSuccess&=pt.success;
				}
			}
			
			putReads(POISON_READS);
			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning reads.");}
			
			//Track whether any threads failed
			if(!allSuccess){errorState=true;}
			if(verbose || verbose2){outstream.println("tid "+tid+" finished!");}
			
		}
		
		void putReads(ListNum<Read> list){
			if(verbose){outstream.println("tid "+tid+" putting rlist size "+list.size());}
			while(list!=null){
				try {
					outq.put(list);
					list=null;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("tid "+tid+" done putting rlist");}
		}
		
		/** Iterate through the reads */
		void makeReads(){
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}
			
			ListNum<byte[]> list=takeBytes();
			while(list!=POISON_BYTES){
				ListNum<Read> reads=new ListNum<Read>(new ArrayList<Read>(list.size()), list.id);
				for(byte[] line : list){
					if(line.length==0 || line[0]=='@'){
						//ignore;
					}else{
						SamLine sl=new SamLine(line);
						Read r=sl.toRead(false);
//						assert(!r.mapped()) : sl+"\n"+r;
						if(!r.validated()){r.validate(true);}
						r.samline=sl;
						
						reads.add(r);
						
						readsProcessedT++;
						basesProcessedT+=r.length();
					}
				}
				if(reads.size()>0){putReads(reads);}
				list=takeBytes();
			}
			if(verbose || verbose2){outstream.println("tid "+tid+" done making reads.");}

			putBytes(POISON_BYTES);
			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning bytes.");}
			
//			putReads(POISON_READS);
//			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning reads.");}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Thread ID */
		final int tid;
		
		ArrayList<ProcessThread> alpt;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final ArrayBlockingQueue<ListNum<Read>> outq;
	
}
