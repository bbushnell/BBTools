package stream;

import java.util.ArrayList;
import fileIO.FileFormat;
import shared.LineParser1;
import shared.Tools;
import structures.ListNum;

/**
 * Loads sam files rapidly with multiple threads.
 * 
 * @author Brian Bushnell
 * @date November 4, 2016
 *
 */
public class SamLineStreamer extends SamStreamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	SamLineStreamer(String fname_, int threads_, boolean saveHeader_, boolean ordered_, 
			long maxReads_, boolean makeReads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false), threads_, 
			saveHeader_, ordered_, maxReads_, makeReads_);
	}
	
	/** Constructor. */
	SamLineStreamer(FileFormat ffin_, int threads_, boolean saveHeader_, boolean ordered_, 
			long maxReads_, boolean makeReads_){
		super(ffin_, threads_, saveHeader_, ordered_, maxReads_, makeReads_);
		final int queueSize=3+(3*threads)/2;
		outq=new JobQueue<ListNum<SamLine>>(queueSize, ordered, true, 0);
		if(verbose){outstream.println("Made sls-"+threads);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean hasMore() {return outq.hasMore();}
	
	@Override
	public boolean paired(){return false;}

	@Override
	public int pairnum(){return 0;}
	
	@Override
	public long readsProcessed() {return readsProcessed;}
	
	@Override
	public long basesProcessed() {return basesProcessed;}
	
	@Override
	public ListNum<SamLine> nextLines(){
		ListNum<SamLine> list=outq.take();
		if(list==null || list.last()) {
			assert(list==null || list.isEmpty());
			assert(!outq.hasMore());
			return null;
		}
		if(verbose && list!=null){outstream.println("Got list size "+list.size());}
		return list;
	}

	@Override
	public ListNum<Read> nextReads(){
		assert(makeReads);
		ListNum<SamLine> lines=nextLines();
		if(lines==null){return null;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		if(!lines.isEmpty()) {
			for(SamLine line : lines){
				assert(line.obj!=null);
				reads.add((Read)line.obj);
			}
		}
		ListNum<Read> ln=new ListNum<Read>(reads, lines.id);
		return ln;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	@Override
	void spawnThreads(){
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
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		/** Constructor */
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_){
			tid=tid_;
			setName("SamLineStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}
		
		/** Called by start() */
		@Override
		public void run(){
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
			if(verbose || verbose2){outstream.println("tid "+tid+" done with processBytes0.");}
			
			success=true;
			
			//Wait for completion of all threads
			boolean allSuccess=true;
			for(ProcessThread pt : alpt){
				
				//Wait until this thread has terminated
				if(pt!=this){
					if(verbose){outstream.println("Waiting for thread "+pt.tid);}
					while(pt.getState()!=Thread.State.TERMINATED){
						try{
							pt.join(); //Attempt a join operation
						}catch(InterruptedException e){
							e.printStackTrace();
						}
					}
					if(verbose || verbose2){outstream.println("tid "+tid+" joined tid "+pt.tid);}

					//Accumulate per-thread statistics
					readsProcessed+=pt.readsProcessedT;
					basesProcessed+=pt.basesProcessedT;
					allSuccess&=pt.success;
//					assert(pt.success) : pt.tid;
				}
			}
			if(verbose || verbose2){outstream.println("tid "+tid+" noted all process threads finished.");}
			
			ListNum<byte[]> list=takeBytes();//Poison
			putReads(new ListNum<SamLine>(null, list.id(), false, true));//Last
			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning reads.");}
			
			//Track whether any threads failed
			if(!allSuccess){errorState=true;}
			if(verbose || verbose2){outstream.println("tid "+tid+" finished! Error="+errorState);}
		}
		
		void putReads(ListNum<SamLine> list){
			if(verbose){outstream.println("tid "+tid+" putting rlist "+list.id+" size "+list.size());}
			outq.add(list);
			if(verbose){outstream.println("tid "+tid+" done putting rlist");}
		}
		
		/** Iterate through the reads */
		void makeReads(){
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}
			
			final LineParser1 lp=new LineParser1('\t');
			ListNum<byte[]> list=takeBytes();
			while(!list.poison()){
				if(verbose || verbose2){outstream.println("tid "+tid+" grabbed blist "+list.id);}
				
				// Apply subsampling if needed
				if(samplerate<1f && randy!=null){
					int nulled=0;
					for(int i=0; i<list.size(); i++){
						if(randy.nextFloat()>=samplerate){
							list.list.set(i, null);
							nulled++;
						}
					}
					if(nulled>0) {Tools.condenseStrict(list.list);}
				}
				
				ListNum<SamLine> reads=new ListNum<SamLine>(
					new ArrayList<SamLine>(list.size()), list.id);
				long readID=list.firstRecordNum;
				for(byte[] line : list){
					if(line[0]=='@'){
						//Ignore header lines
					}else{
						SamLine sl=new SamLine(lp.set(line));
						reads.add(sl);
						if(makeReads){
							Read r=sl.toRead(FASTQ.PARSE_CUSTOM);
							sl.obj=r;
							r.samline=sl;
							r.numericID=readID++;
							if(!r.validated()){r.validate(true);}
						}
						readsProcessedT++;
						basesProcessedT+=(sl.seq==null ? 0 : sl.length());
					}
				}
//				if(reads.size()>0){putReads(reads);}
				putReads(reads);//Must put them no matter what
				list=takeBytes();
			}
			if(verbose || verbose2){outstream.println("tid "+tid+" done making reads.");}

			putBytes(list);
			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning bytes.");}
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
	
	final JobQueue<ListNum<SamLine>> outq;
	
}