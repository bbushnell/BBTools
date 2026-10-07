package stream;

import java.io.File;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.FileFormat;
import fileIO.ReadWrite;

/**
 * Submits read lists to one or two background writers, optionally restoring list order.
 * Synchronized submissions copy the outer list; Read objects and mates remain borrowed.
 * Ordered callers supply every ID from zero, including empty lists for filtered batches.
 * Unordered submission bypasses only the ordering table, not backend buffering.
 * For normal shutdown, stop producers, call close, then join. Coordinate ordering
 * resets with producers. Status accessors do not close or join the writers.
 *
 * @author Brian Bushnell
 * @date Jan 26, 2015
 */
public final class ConcurrentGenericReadOutputStream extends ConcurrentReadOutputStream{
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructs the selected backend writers and optional ordered-list table.
	 * Ordering is captured from the primary descriptor by the base constructor.
	 * SAM/BAM uses ReadStreamSamWriter when its ReadWrite selector is enabled; that
	 * adapter starts its delegate during construction, before this stream's start().
	 * Otherwise, byte writers are created, with a separate mate writer only when the
	 * primary is not standard I/O and a secondary descriptor is present.
	 * @param ff1_ Nonnull primary output descriptor
	 * @param ff2_ Optional mate descriptor; unused by the selected SAM/BAM adapter
	 * @param qf1 Optional primary quality filename for the byte-writer route
	 * @param qf2 Optional mate quality filename for the byte-writer route
	 * @param rswBuffers Buffering capacity forwarded to each selected backend
	 * @param header Optional header text forwarded to the backend
	 * @param useSharedHeader Whether to request shared-header handling */
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
		
		if(ff1.samOrBam() && ReadWrite.USE_READ_STREAM_SAM_WRITER){
			readstream1=new ReadStreamSamWriter(ff1, rswBuffers, header, useSharedHeader);
			readstream2=null;
		}else{
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
	
	/** Starts the outer backend threads once, before any submissions.
	 * @throws RuntimeException On a repeated call, after resetting the next ordered ID to zero */
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
	 * Submits a shallow list copy; callers may reuse the list after this call returns.
	 * Read payloads must remain stable until downstream writing completes.
	 * Ordered IDs start at zero and must be unique and contiguous, including empty lists.
	 * Future IDs can wait for ordering capacity; the currently required ID bypasses
	 * that wait. Either mode may wait for capacity in the backend writers.
	 * @param list Nonnull list of borrowed Read references
	 * @param listnum Ordered batch ID, ignored by the unordered route
	 * @throws RuntimeException If this stream has been aborted or a backend has terminated
	 */
	@Override
	public synchronized void add(ArrayList<Read> list, long listnum){
		if(aborted){throw abortedException(listnum);}
		if(ordered){
			int size=table.size();
//			System.err.print(size+", ");
			final boolean flag=(size>=HALF_LIMIT);
			//Only future IDs wait for ordering capacity; the required ID can enter and drain
			//the table while wait() releases this monitor. This assumes every required ID arrives.
			//Historical #002: a partial drain can miss notification when flag was false at entry.
			//The timed wait was reduced from 20000ms to 500ms to recheck that condition sooner;
			//it does not bound scheduler delays, backend I/O or the total time spent in add().
			if(listnum>nextListID && size>=ADD_LIMIT){
				if(printBufferNotification){
					System.err.println("Output buffer became full; key "+listnum+" waiting on "+nextListID+".");
					printBufferNotification=false;
				}
				while(!aborted && listnum>nextListID && size>=HALF_LIMIT){
					try{
						this.wait(500);//#002: timed recheck only while ordering capacity blocks this future ID.
					}catch(InterruptedException e){
						e.printStackTrace();
					}
					size=table.size();
				}
				if(aborted){throw abortedException(listnum);}
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
	 * Requests writer termination without joining; all producers must have stopped.
	 * A nonempty ordering table marks an error because some required IDs never arrived.
	 * Sends poison to the selected backends; an already aborted stream is left alone.
	 */
	@Override
	public synchronized void close(){
		if(aborted){return;}
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
	 * Discards invalid buffered output and wakes this stream's ordered-list waiters
	 * after an upstream failure. Downstream writer I/O performed while another
	 * synchronized caller holds this monitor can still delay entry into abort().
	 */
	@Override
	public synchronized void abort(){
		if(aborted){return;}
		aborted=true;
		errorState=true;
		finishedSuccessfully=false;
		if(table!=null){table.clear();}
		if(readstream1!=null){readstream1.abortNow();}
		if(readstream2!=null){readstream2.abortNow();}
		notifyAll();
	}
	
	/**
	 * Waits for previously started outer writer threads, after close or abort was requested.
	 * Does not request termination itself. Interrupted waits print a trace and retry.
	 * After joining, asserts that ordered pending data is empty and records whether
	 * this stream was aborted; finishedSuccessfully() also checks child/error flags.
	 */
	@Override
	public void join(){
		while(readstream1!=null && readstream1.getState()!=Thread.State.TERMINATED){
			try{
				readstream1.join();
			}catch(InterruptedException e){
				//Report the interruption and retry the join.
				e.printStackTrace();
			}
		}
		while(readstream2!=null && readstream2.getState()!=Thread.State.TERMINATED){
			try{
				if(readstream2!=null){readstream2.join();}
			}catch(InterruptedException e){
				//Report the interruption and retry the join.
				e.printStackTrace();
			}
		}
		synchronized(this){
			assert(table==null || table.isEmpty());
			finishedSuccessfully=!aborted;
		}
	}
	
	/**
	 * Waits for ordered pending data to drain, then resets its next ID to zero.
	 * Requires an ordered instance and coordinated producers; does not clear data or
	 * wait for already submitted backend writes. Unordered instances have no table.
	 * The initial loop permits 2000 waits of 2000ms, then warns and continues waiting
	 * without an iteration limit. Notifications and scheduling affect elapsed time.
	 * Historical #001 corrected the former "4 minutes" claim; the wait counts remain.
	 */
	@Override
	public synchronized void resetNextListID(){
		for(int i=0; i<2000 && !table.isEmpty(); i++){
			try{this.wait(2000);}catch(InterruptedException e){e.printStackTrace();}
		}
		if(!table.isEmpty()){
			System.err.println("WARNING! resetNextListID() waited a long time and the table never cleared.  Process may have stalled.");
		}
		while(!table.isEmpty()){
			try{this.wait(2000);}catch(InterruptedException e){e.printStackTrace();}
		}
		nextListID=0;
	}
	
	/** @return Primary backend filename */
	@Override
	public final String fname(){
//		if(STANDARD_OUT){return "stdout";}
		return readstream1.fname();
	}
	
	/** Combines currently observed error flags without waiting for completion.
	 * @return true if this stream or either selected backend reports an error */
	@Override
	public boolean errorState(){
		return errorState || (readstream1!=null && readstream1.errorState()) || (readstream2!=null && readstream2.errorState());
	}
	
	/** Checks stored completion/error flags without closing or joining anything.
	 * @return true after a successful join when this stream and both present backends
	 * report success, with no local error or abort */
	@Override
	public synchronized boolean finishedSuccessfully(){
		return !errorState && !aborted && finishedSuccessfully &&
				(readstream1==null || readstream1.finishedSuccessfully()) &&
				(readstream2==null || readstream2.finishedSuccessfully());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** Copies one unique ordered list, then drains consecutive IDs from the expected ID.
	 * @param list Nonnull list whose Read payloads remain borrowed
	 * @param listnum Pending ID at or after the currently expected ID */
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
		//Historical #002 (2026-06-18): this notifies only on a full drain; add()'s other
		//notification depends on its entry flag. The unchanged 500ms polling wait covers
		//a missed partial-drain notification without promising a wall-clock completion bound.
		if(table.isEmpty()){notifyAll();}
	}
	
	/** Copies the outer list and submits it directly to the backend queues.
	 * @param list Nonnull list of borrowed Read references
	 * @param listnum Unused compatibility argument */
	private synchronized void addDisordered(ArrayList<Read> list, long listnum){
		assert(list!=null);
		assert(table==null);
		write(new ArrayList<Read>(list));
	}

	/** Creates a diagnostic rejection for a submission after abort.
	 * @param listnum Rejected input batch ID
	 * @return Exception naming the attempted batch and primary output */
	private RuntimeException abortedException(long listnum){
		return new RuntimeException("Cannot add list "+listnum+" to aborted output stream "+fname()+".");
	}
	
	/** Enqueues the supplied list in each present backend, which may wait for capacity.
	 * @param list Already copied list; shared by both mate writers when present
	 * @throws RuntimeException If a backend thread has already terminated */
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
	
	/** @return Live primary backend; includes interleaved mates when there is no secondary */
	@Override
	public final ReadStreamWriter getRS1(){return readstream1;}
	/** @return Live secondary backend, or null for single/interleaved or selected SAM/BAM output */
	@Override
	public final ReadStreamWriter getRS2(){return readstream2;}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary writer selected at construction. */
	private final ReadStreamWriter readstream1;
	/** Optional separate-mate writer; null for the single-backend routes. */
	private final ReadStreamWriter readstream2;
	/** Next required ordered ID, guarded by this object's monitor. */
	private long nextListID=0;
	
	/** Initial ordered-table capacity and basis for the admission thresholds. */
	private final int MAX_CAPACITY=256;
	/** Future-ID table size at which an add begins waiting. */
	private final int ADD_LIMIT=MAX_CAPACITY-2;
	/** Table size below which an existing future-ID waiter may resume. */
	private final int HALF_LIMIT=ADD_LIMIT/2;
	
	/** Pending shallow-copied lists, guarded by this monitor; null when unordered. */
	private final HashMap<Long, ArrayList<Read>> table;
	/** Abort state guarded by this monitor. */
	private boolean aborted=false;
	
	{if(HALF_LIMIT<1){throw new RuntimeException("Capacity too low.");}}
	
	/*--------------------------------------------------------------*/
	/*----------------         Diagnostics          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Controls the initial ordering-capacity notification; accessed under this monitor. */
	private boolean printBufferNotification=true;
	
}
