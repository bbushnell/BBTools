package stream;

import java.io.File;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.FileFormat;
import fileIO.ReadWrite;

/**
 * Thread-safe output stream for writing sequence reads to files concurrently.
 * Extends ConcurrentReadOutputStream with generic file format support and ordering.
 * Manages dual ReadStreamByteWriter instances for paired-end output with buffering.
 *
 * @author Brian Bushnell
 * @date Jan 26, 2015
 */
public final class ConcurrentGenericReadOutputStream extends ConcurrentReadOutputStream {
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	ConcurrentGenericReadOutputStream(FileFormat ff1_, FileFormat ff2_, String qf1, String qf2, int rswBuffers, CharSequence header, boolean useSharedHeader){
		super(ff1_, ff2_);
		
		if(verbose){
			System.err.println("ConcurrentGenericReadOutputStream("+ff1+", "+ff2+", "+qf1+", "+qf2+", "+rswBuffers+", "+useSharedHeader+")");
		}
		
		assert(ff1!=null);
		assert(!ff1.text() && !ff1.unknownFormat()) : "Unknown format for "+ff1;
		
		if(ff1.hasName() && ff1.devnull()){
			File f=new File(ff1.name());
			assert(ff1.overwrite() || !f.exists() || ff1.name().equals("/dev/null")) : f.getAbsolutePath()+" already exists; please delete it.";
			if(ff2!=null){assert(!ff1.name().equals(ff2.name())) : ff1.name()+"=="+ff2.name();}
		}
		
		if(ff1.samOrBam() && ReadWrite.USE_READ_STREAM_SAM_WRITER) {
			readstream1=new ReadStreamSamWriter(ff1, rswBuffers, header, useSharedHeader);
			readstream2=null;
		}else {
			readstream1=new ReadStreamByteWriter(ff1, qf1, true, rswBuffers, header, useSharedHeader);
			readstream2=ff1.stdio() || ff2==null ? null : new ReadStreamByteWriter(ff2, qf2, false, rswBuffers, header, useSharedHeader);
		}
		
		if(readstream2==null && readstream1!=null){
//			System.out.println("ConcurrentReadOutputStream detected interleaved output.");
			readstream1.OUTPUT_INTERLEAVED=true;
		}
		
		table=(ordered ? new HashMap<Long, ArrayList<Read>>(MAX_CAPACITY) : null);
		
		assert(readstream1==null || readstream1.read1==true);
		assert(readstream2==null || (readstream2.read1==false));
	}
	
	@Override
	public synchronized void start(){
		if(started){
			System.err.println("Resetting output stream.");
			nextListID=0;
			throw new RuntimeException();
		}else{
			started=true;
			if(readstream1!=null){readstream1.start();}
			if(readstream2!=null){readstream2.start();}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds a list of reads to the output queue with optional ordering.
	 * For ordered output, blocks if buffer becomes full waiting for sequential IDs.
	 * For unordered output, writes immediately without buffering.
	 *
	 * @param list Read list to add to output queue
	 * @param listnum Sequential identifier for ordering (ignored if unordered)
	 */
	@Override
	public synchronized void add(ArrayList<Read> list, long listnum){
		
		if(ordered){
			int size=table.size();
//			System.err.print(size+", ");
			final boolean flag=(size>=HALF_LIMIT);
			//BACKPRESSURE (deadlock-free by construction): only a list STRICTLY AHEAD of order (listnum>nextListID) blocks when the buffer is full. The list AT nextListID - the one whose write advances the order and drains the buffer - is NEVER caught by this guard, so it always gets through, writes, and bumps nextListID. wait() releases 'this', so a sibling add() carrying nextListID runs while this one waits. NOTE (#002, LOW latency, RESOLVED): the notifyAll coverage is INCOMPLETE - a drain that drops the buffer below HALF_LIMIT without emptying it, when flag was false at entry, signals neither waiter path. Rather than add fragile notify-tracking, the wait timeout below was tightened 20000ms->500ms, so a missed signal costs at most ~0.5s instead of ~20s (and the poll only runs during active backpressure). Not a deadlock/correctness bug; a bounded latency tail.
			if(listnum>nextListID && size>=ADD_LIMIT){
				if(printBufferNotification){
					System.err.println("Output buffer became full; key "+listnum+" waiting on "+nextListID+".");
					printBufferNotification=false;
				}
				while(listnum>nextListID && size>=HALF_LIMIT){
					try {
						this.wait(500);//#002 fix: 500ms (was 20000) bounds the worst-case spurious stall from the incomplete notifyAll coverage to ~0.5s instead of ~20s. Cheap: this only polls while a producer is actually blocked (active backpressure - rare).
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
					size=table.size();
				}
				if(printBufferNotification){
					System.err.println("Output buffer became clear for key "+listnum+"; next="+nextListID+", size="+size);
				}
			}
			addOrdered(list, listnum);
			assert(listnum!=nextListID);
			if(flag && listnum<nextListID){this.notifyAll();}
		}else{
			addDisordered(list, listnum);
		}
	}
	
	/**
	 * Closes the output stream and terminates writer threads.
	 * Sets error state if unfinished lists remain in buffer.
	 * Poisons ReadStreamByteWriter instances to trigger shutdown.
	 */
	@Override
	public synchronized void close(){
		
		if(table!=null && !table.isEmpty()){
			errorState=true;
			System.err.println("Error: An unfinished ReadOutputStream was closed.");
		}
		//assert(table==null || table.isEmpty()); //TODO Seems like a race condition.  Probably, I should wait at this point until the condition is true before proceeding.
		//REASONED (not a race under valid usage): add() and close() are BOTH synchronized(this), so they cannot run concurrently; under the contract (all add()s by the producer, THEN close()), a non-empty table at close means a real GAP in list numbering (some listnum never arrived) -> those buffered lists are genuinely unwritten -> the errorState above is correct gap-detection, not a transient race. The disabled assert would false-fire on that legitimate-error path. (Full confirmation of the writer-thread side is gated on ReadStreamByteWriter's own V2 review.)
		
//		readstream1.addList(null);
//		if(readstream2!=null){readstream2.addList(null);}
		readstream1.poison();
		if(readstream2!=null){readstream2.poison();}
	}
	
	/**
	 * Waits for all writer threads to complete before returning.
	 * Ensures proper cleanup by joining both ReadStreamByteWriter instances.
	 * Sets finishedSuccessfully flag upon completion.
	 */
	@Override
	public void join(){
		while(readstream1!=null && readstream1.getState()!=Thread.State.TERMINATED){
			try {
				readstream1.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		while(readstream2!=null && readstream2.getState()!=Thread.State.TERMINATED){
			try {
				if(readstream2!=null){readstream2.join();}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		assert(table==null || table.isEmpty());
		finishedSuccessfully=true;
	}
	
	/**
	 * Resets the next list ID counter to zero after clearing buffers.
	 * Waits up to ~66 minutes (2000 iterations x wait(2000ms)) for the table to clear, warns once past that, then waits indefinitely.
	 * (#001 doc fix: the javadoc previously said "4 minutes" but 2000 x 2000ms = 4000s; the loop bound may be higher than intended - flagged, code left unchanged.)
	 * Issues a warning if the table doesn't clear within the timeout period.
	 */
	@Override
	public synchronized void resetNextListID(){
		for(int i=0; i<2000 && !table.isEmpty(); i++){
			try {this.wait(2000);}
			catch (InterruptedException e) {e.printStackTrace();}
		}
		if(!table.isEmpty()){
			System.err.println("WARNING! resetNextListID() waited a long time and the table never cleared.  Process may have stalled.");
		}
		while(!table.isEmpty()){
			try {this.wait(2000);}
			catch (InterruptedException e) {e.printStackTrace();}
		}
		nextListID=0;
	}
	
	/** Gets the filename of the primary output stream.
	 * @return Filename from primary ReadStreamByteWriter */
	@Override
	public final String fname(){
//		if(STANDARD_OUT){return "stdout";}
		return readstream1.fname();
	}
	
	/** Checks if any component is in an error state.
	 * @return true if this stream or either ReadStreamByteWriter has errors */
	@Override
	public boolean errorState(){
		return errorState || (readstream1!=null && readstream1.errorState()) || (readstream2!=null && readstream2.errorState());
	}
	
	/** Checks if all components finished without errors.
	 * @return true if this stream and both ReadStreamByteWriter instances completed successfully */
	@Override
	public boolean finishedSuccessfully(){
		return finishedSuccessfully && (readstream1==null || readstream1.finishedSuccessfully()) && (readstream2==null || readstream2.finishedSuccessfully());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	private synchronized void addOrdered(ArrayList<Read> list, long listnum){
//		System.err.println("RTOS got "+listnum+" of size "+(list==null ? "null" : list.size())+
//				" with first read id "+(list==null || list.isEmpty() || list.get(0)==null ? "null" : ""+list.get(0).numericID));
		assert(list!=null) : listnum;
		assert(listnum>=nextListID) : listnum+", "+nextListID;
//		assert(list.isEmpty() || list.get(0)==null || list.get(0).numericID>=nextReadID) : list.get(0).numericID+", "+nextReadID;
		assert(!table.containsKey(listnum));
		
		table.put(listnum, new ArrayList<Read>(list));//defensive COPY: the caller may reuse/recycle its list after add() returns

		//Flush every now-contiguous list starting at nextListID, in strict order, advancing nextListID until a gap (a missing listnum) stops it.
		while(table.containsKey(nextListID)){
//			System.err.println("Writing list "+first.get(0).numericID);
			ArrayList<Read> value=table.remove(nextListID);
			write(value);
			nextListID++;
		}
		//#002 (LOW latency, RESOLVED 2026-06-18): this notifies ONLY on a full drain (table empty), and add()'s other notify is gated on flag, so a partial drain below HALF_LIMIT can miss waking a blocked producer. Rather than make this notify edge-perfect (a waiter-count or size-crossing test - more shared state, more to test), the backpressure wait above was tightened 20000ms->500ms, bounding any missed-signal stall to ~0.5s. This notify left as-is.
		if(table.isEmpty()){notifyAll();}
	}
	
	private synchronized void addDisordered(ArrayList<Read> list, long listnum){
		assert(list!=null);
		assert(table==null);
		write(new ArrayList<Read>(list));
	}
	
	private synchronized void write(ArrayList<Read> list){
		//Crash-loud guard: writing to an already-terminated writer thread would silently drop the lists, so fail loudly instead.
		if(readstream1!=null){
			if(readstream1.getState()==State.TERMINATED){throw new RuntimeException("Writing to a terminated thread.");}
			readstream1.addList(list);
		}
		if(readstream2!=null){
			if(readstream2.getState()==State.TERMINATED){throw new RuntimeException("Writing to a terminated thread.");}
			readstream2.addList(list);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets the primary ReadStreamWriter instance.
	 * @return Primary ReadStreamWriter for first-in-pair or single-end reads */
	@Override
	public final ReadStreamWriter getRS1(){return readstream1;}
	@Override
	public final ReadStreamWriter getRS2(){return readstream2;}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	private final ReadStreamWriter readstream1;
	private final ReadStreamWriter readstream2;
	private long nextListID=0;
	
	private final int MAX_CAPACITY=256;
	private final int ADD_LIMIT=MAX_CAPACITY-2;
	private final int HALF_LIMIT=ADD_LIMIT/2;
	
	private final HashMap<Long, ArrayList<Read>> table;
	
	{if(HALF_LIMIT<1){throw new RuntimeException("Capacity too low.");}}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean printBufferNotification=true;
	
}
