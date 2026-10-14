package fileIO;

import java.io.InputStream;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;

import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.bam.BgzfSettings;
import structures.ListNum;

/**
 * Buffers ByteFile1 lines on a background thread for a serialized consumer.
 * Construction opens the wrapped reader and starts its producer. Four reusable
 * arrays carry line references through full/empty queues; individual line arrays
 * are not recycled. Delimiter handling follows ByteFile1. A null slot ends a
 * partial buffer, while a shared zero-length outer array marks queue termination.
 *
 * Direct nextLine and pushBack calls require caller serialization. Inherited
 * nextList serializes its own callers and numbers record origins; nextListAsync
 * uses the base scalar batching implementation without that monitor. Mixing direct
 * calls with concurrent batch consumers is not coordinated by this wrapper.
 * lineNum observes the wrapped reader, which may read ahead or reset its counter
 * on closure; it is not the count of lines delivered to this consumer.
 *
 * close requests shutdown, joins the producer, closes the wrapped reader and
 * reports a sticky error flag. Reset closes and constructs a new reader/thread;
 * it requires a reopenable input and does not clear earlier errors. Interrupted
 * queue/join waits print traces and retry, rather than implementing cancellation.
 *
 * @author Brian Bushnell
 * @contributor Shinobu
 * @date Sep 23, 2013
 */
public final class ByteFile2 extends ByteFile{
	
	/** Prints a line range or reports whole-remaining-input speed statistics.
	 * Retains the diagnostic close/reset/close sequence and reports its error flags.
	 * Print mode also checks borrowed stdout, with read/close errors reported first.
	 * The simd option enables SIMD only when the optional capability probe succeeds.
	 * @param args Positional filename, first line or speedtest, last line or simd,
	 * and optional native decoder request; defaults to stdin and the first 100 lines */
	public static void main(String[] args){
		long first=0, last=100;
		boolean speedtest=false;
		if(args.length>1){
			if(args[1].equalsIgnoreCase("speedtest")){
				speedtest=true;
				first=0;
				last=Long.MAX_VALUE;
			}else{
				first=Integer.parseInt(args[1]);
				last=first+100;
			}
		}
		if(args.length>2){
			if(args[2].equalsIgnoreCase("simd")){
				//FIO-018: preserve scalar fallback when the optional Vector API is unavailable.
				if(simd.Vector.simd256){Shared.SIMD=true;}
			}else{
				last=Integer.parseInt(args[2]);
			}
		}
		if(args.length>3){
			if(args[3].equalsIgnoreCase("native")){
				ReadWrite.ALLOW_NATIVE_BGZF=BgzfSettings.USE_MULTITHREADED_BGZF=true;
			}
		}
		ByteFile2 tf=new ByteFile2(args.length>0 ? args[0] : "stdin", true);
		speedtest(tf, first, last, !speedtest);
		//FIO-019: capture PrintStream failure, but preserve read/close error priority below.
		final boolean outputError=!speedtest && System.out.checkError();
		
		//FIO-017: retain errors from either pass while preserving the diagnostic reset sequence.
		boolean error=tf.close();
		tf.reset();
		error|=tf.close();
		if(error){throw new RuntimeException("ByteFile2: error reading or closing "+tf.name());}
		if(outputError){throw new RuntimeException("ByteFile2: error writing output");}
	}
	
	/** Consumes a scalar prefix, then prints the selected range or drains batch input.
	 * @param bf Reader consumed by the diagnostic
	 * @param first Number of initial scalar lines to skip
	 * @param last Exclusive print limit; batch speed mode ignores this limit
	 * @param reprint Whether to print decoded Strings to borrowed stdout */
	private static void speedtest(ByteFile2 bf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){bf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=bf.nextLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length+1;
				System.out.println(new String(s));
			}
			
			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
//			for(long i=first; i<last; i++){
//				byte[] s=tf.nextLine();
//				if(s==null){break;}
//				lines++;
//				bytes+=s.length+1;
//			}
			for(ListNum<byte[]> ln=bf.nextListAsync(); ln!=null; ln=bf.nextListAsync()){
				for(byte[] line : ln.list){
					lines++;
					bytes+=line.length+1;
				}
			}
		}
		t.stop();
		
		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
		}
	}
	
//	public ByteFile2(String name()){this(name(), false);}
	
	/*--------------------------------------------------------------*/
	/*----------------  Initialization  ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves a text descriptor and starts background reading.
	 * @param fname Input filename or supported stream alias
	 * @param allowSubprocess_ Whether subprocess decoding is allowed */
	public ByteFile2(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}
	
	/** Opens a wrapped ByteFile1 and starts its producer thread.
	 * @param ff Input descriptor retained for the initial open and later resets */
	public ByteFile2(FileFormat ff){
		super(ff);
		if(verbose){System.err.println("ByteFile2("+ff+")");}
		open();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------  Reader API  ----------------*/
	/*--------------------------------------------------------------*/

	/** Closes the current pipeline, creates a new reader/thread and resets base batch counters.
	 * Reopens the same descriptor; it cannot rewind a consumed stdin or other one-shot
	 * source. The wrapper's read/close error flag remains sticky across resets. */
	@Override
	public final void reset(){
		close();
		open();
		superReset();
	}
	
	/** Requests producer shutdown, joins it and incorporates the wrapped close status.
	 * Clears thread/cache/pushback references after the join. Callers stop consuming
	 * before closure; synchronization here does not coordinate direct scalar calls.
	 * @return Sticky read/close error flag; repeated closed calls preserve it */
	@Override
	public synchronized final boolean close(){
		if(verbose){System.err.println("ByteFile2("+name()+").close()");}
		if(isOpen()){
//			errorState|=ReadWrite.killProcess(name());
			thread.shutdown();
			while(thread.getState()!=Thread.State.TERMINATED){
				try{
					thread.join();
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			errorState|=thread.bf1.close();//Fold the wrapped ByteFile1's read-error state into ours (else a truncated input is lost)
		}
		thread=null;
		currentList=null;
		currentLoc=0;
//		assert(numIn==numOut) : numIn+", "+numOut;
		pushBack=null;
		if(verbose){System.err.println("ByteFile2("+name()+").close() returned "+errorState);}
		return errorState;
	}
	
	/** Returns a pending pushback line, then buffered producer lines in order.
	 * The borrowed outer buffer is recycled only after consumption; returned line
	 * arrays remain available to the caller. Requires serialized consumption.
	 * @return Next line array, possibly empty, or null after the poison buffer */
	@Override
	public final byte[] nextLine(){
		
		if(pushBack!=null){//Commenting out does not seem to improve speed here.
			byte[] temp=pushBack;
			pushBack=null;
			return temp;
		}
		
//		if(verbose){System.err.println("Reading line.");}
//		byte[] r=null;
		
		byte[][] temp=currentList;
		int tempLoc=currentLoc;
		
		if(temp==null || tempLoc>=temp.length || temp[tempLoc]==null){
			boolean b=getBuffer();
			if(!b){
				if(verbose2){System.err.println("nextLine()->getBuffer() returned false.");}
				return null;
			}
			temp=currentList;
			tempLoc=currentLoc;
			if(temp==null || temp==poison || temp[tempLoc]==null){
				return null;
			}
		}
		
		//TODO: This is a race condition; currentList can be changed to null.  A defensive copy could be created.
		//Note that I read the above warning and added "temp" and "temploc" but I'm not sure if that fixed anything.
		assert(temp!=null && temp!=poison);
		assert(tempLoc<temp.length);
		assert(temp[tempLoc]!=null);
		byte[] r=temp[tempLoc];
		assert(r!=null);
		currentLoc++;
//		numOut++;
		return r;
	}
	
	/** Clears/recycles the exhausted outer array and waits for the next full buffer.
	 * The poison array is retained as terminal state and is never recycled.
	 * @return Whether a nonterminal buffer was obtained */
	private boolean getBuffer(){
		if(verbose2){System.err.println("Getting new buffer.");}
		currentLoc=0;
		final BF1Thread bft=thread;
		if(bft==null){
			currentList=null;
			if(verbose2){System.err.println("No buffers available.  thread="+thread+", shutdown="+(thread==null ? "X" : ""+thread.shutdown));}
			return false;
		}
		if(currentList==poison){
			if(verbose2){System.err.println("A: Current list is poison.");}
			return false;
		}
		if(currentList!=null){
			Arrays.fill(currentList, null); //MUST be done or lines get recycled at end of file.
			while(currentList!=null){
				try{
					if(verbose2){System.err.println("adding to qEmpty list size "+currentList.length+"\n"+toString(currentList));}
					bft.qEmpty.put(currentList);
					currentList=null;
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		assert(currentList==null);
		while(currentList==null){
			try{
				assert(bft!=null);
				if(verbose2){System.err.println("C: qFull.size()="+bft.qFull.size());}
				currentList=bft.qFull.take();
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(verbose2){
			if(currentList==poison){
				System.err.println("B: Current list is poison.");
			}else{
				System.err.println("getBuffer fetched a new buffer of size "+currentList.length);
			}
		}
		return currentList!=poison;
	}
	
	/** Creates a new wrapped reader and starts a producer; requires no current thread.
	 * @return Newly started producer, whose queues own the reusable outer arrays */
	private final synchronized BF1Thread open(){
		if(verbose2){System.err.println("ByteFile2("+name()+").open()");}
		assert(thread==null);
		currentList=null;
		currentLoc=0;
//		numIn=0;
//		numOut=0;
		thread=new BF1Thread(ff);
		thread.start();
		return thread;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------  Producer thread  ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads through ByteFile1 and transfers line references in reusable outer arrays. */
	private class BF1Thread extends Thread{
		
//		public BF1Thread(String fname){
//			bf1=new ByteFile1(fname, false, allowSubprocess);
//			qFull=new ArrayBlockingQueue<byte[][]>(buffs+2);
//			qEmpty=new ArrayBlockingQueue<byte[][]>(buffs+2);
//			for(int i=0; i<buffs; i++){
//				try{
//					qEmpty.put(new byte[bufflen][]);
//				}catch(InterruptedException e){
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
//			}
//		}
		
		/** Opens the wrapped reader and seeds the empty queue with four arrays.
		 * Both queues reserve two additional slots beyond the ordinary buffer count.
		 * @param ff Descriptor used by the wrapped reader */
		public BF1Thread(FileFormat ff){
			bf1=new ByteFile1(ff);
			qFull=new ArrayBlockingQueue<byte[][]>(buffs+2);
			qEmpty=new ArrayBlockingQueue<byte[][]>(buffs+2);
			for(int i=0; i<buffs; i++){
				try{
					qEmpty.put(new byte[bufflen][]);
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		/** Transfers batches until input ends or a shutdown marker is observed.
		 * Publishes a buffer at 1000 lines or 256000 payload bytes, allowing the last
		 * line to exceed the byte threshold. Publishes any partial tail, then poison. */
		@Override
		public void run(){
			if(verbose){System.err.println("ByteFile2("+name()+").run()");}
			byte[] s=null;
			byte[][] list=null;
			while(list==null){
				try{
					list=qEmpty.take();
				}catch(InterruptedException e1){
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}
			synchronized(this){
				if(list==poison || shutdown){
					shutdown();
					return;
				}
			}
			
			int loc=0;
			long bases=0;
			
			//At this point, list is not null
			for(s=bf1.nextLine(); s!=null; s=bf1.nextLine()){
				bases+=s.length;
				assert(list!=null) : "Somehow the list became null for "+bf1.name()+" at line "+cntr;
				list[loc]=s;
				loc++;
//				numIn++;
//				if(verbose){System.err.println("Added line "+numIn);}
				if(loc>=bufflen || bases>=buffcapacity){
					if(verbose2){System.err.println("Capacity exceeded.");}
					while(list!=null){
						try{
//							synchronized(this){
//								if(!shutdown){
									if(verbose2){
										System.err.println("A: Adding to qFull list of size "+loc);
										System.err.println(ByteFile2.toString(list));
									}
									cntr+=list.length;
									qFull.put(list);
									if(verbose2){System.err.println("A: qFull.size()="+qFull.size());}
//								}
//							}
							list=null;
							loc=0;
						}catch(InterruptedException e){
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					//At this point, list is null
					if(shutdown){
						if(verbose2){System.err.println("Break 1");}
						break;
					}
					while(list==null){
						if(verbose2){System.err.println("Taking empty list.");}
						try{
							list=qEmpty.take();
						}catch(InterruptedException e1){
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}
					}
					//At this point, list is not null
					bases=0;
					if(list==poison){
						if(verbose2){System.err.println("Break 2");}
						break;
					}
					//At this point, list is not null
				}
			}
			if(verbose2){System.err.println("Run loop exit.");}
			
			while(list!=null && loc>0){
				try{
//					synchronized(this){
//						if(!shutdown){
							if(verbose2){System.err.println("B: Adding list of size "+loc);}
							qFull.put(list);
							if(verbose2){System.err.println("B: qFull.size()="+qFull.size());}
//						}
//					}
					list=null;
					loc=0;
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			//At this point, list is null
			shutdown();
			
			if(verbose){System.err.println("ByteFile2("+name()+").run() finished");}
		}
		
		/** Sets shutdown once and adds poison to both queues; does not close the stream.
		 * The wrapper's close joins this thread before closing its ByteFile1. */
		synchronized void shutdown(){
			if(verbose || verbose2){System.err.println("ByteFile2("+name()+").shutdown()");}
			if(shutdown){return;}
			shutdown=true;
			if(verbose2){System.err.println("Adding poison.");}
			qFull.add(poison);
			qEmpty.add(poison);
			if(verbose2){System.err.println("D: qFull.size()="+qFull.size());}
			if(verbose || verbose2){System.err.println("ByteFile2("+name()+").shutdown() finished");}
		}
		
		/** Shutdown request/marker-publication latch; not a thread-termination test. */
		private boolean shutdown=false;
		/** Wrapped reader opened before this thread starts. */
		final ByteFile1 bf1;
		/** Filled outer buffers or the shared poison marker, consumed by the wrapper. */
		final ArrayBlockingQueue<byte[][]> qFull;
		/** Cleared outer buffers or poison, consumed by the producer. */
		final ArrayBlockingQueue<byte[][]> qEmpty;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------  State and pushback  ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks retained consumer-buffer/thread references without waiting.
	 * Producer completion or exhausted input does not clear the thread reference.
	 * @return True while a cached line or producer reference remains; not an EOF test */
	@Override
	public boolean isOpen(){
		final byte[][] list=currentList;
		final int loc=currentLoc;
		if(list!=null && loc<list.length && list[loc]!=null){return true;}
		final BF1Thread bft=thread;
		if(bft==null){
			return false;
		}
		return true;
//		synchronized(bft){
//		//NOTE!!!  This cannot be used because qFull.size() will not return a correctly synchronized value.  Poll() may work.
//			assert(bft.bf1.isOpen() || !bft.qFull.isEmpty()) : bft.bf1.isOpen()+", "+bft.qFull.isEmpty()+", "+bft.qFull.size();
//			return (bft.bf1.isOpen() || !bft.qFull.isEmpty());
//		}
	}
	
	/** Retains one borrowed line for the next scalar read, before buffered input.
	 * Requires no pending pushback; null leaves the slot empty. Does not change
	 * the wrapped producer's line counter. Requires caller serialization.
	 * @param line Line reference to return without copying */
	@Override
	public final void pushBack(byte[] line){
		assert(pushBack==null);
		pushBack=line;
	}

//	@Override
//	public void pushBack(byte[] line){
//		if(bstart>line.length){
//			bstart--;
//			buffer[bstart]='\n';
//			for(int i=0, j=bstart-line.length; i<line.length; i++, j++){
//				buffer[j]=line[i];
//			}
//			bstart=bstart-line.length;
//			return;
//		}
//		
//		int bLen=bstop-bstart;
//		int newLen=bLen+line.length+1;
//		int rShift=line.length+1-bstart;
//		assert(rShift>0) : bstop+", "+bstart+", "+line.length;
//		while(newLen>buffer.length){
//			//This could get big if pushback is used often,
//			//unless special steps are taken to prevent it, like leaving extra space for pushbacks.
//			buffer=Arrays.copyOf(buffer, buffer.length*2);
//		}
//		
//		Tools.shiftRight(buffer, rShift);
//		
//		for(int i=0; i<line.length; i++){
//			buffer[i]=line[i];
//		}
//		buffer[line.length]='\n';
//		bstart=0;
//		bstop=newLen;
//	}
	
	/** Renders every outer-array slot for diagnostic logging, including null entries.
	 * @param x Nonnull array of line references
	 * @return Decoded lines/null markers, each followed by LF */
	private static String toString(byte[][] x){
		StringBuilder sb=new StringBuilder();
		for(byte[] z : x){
			sb.append(z==null ? "null" : new String(z)).append('\n');
		}
		return sb.toString();
	}
	
	/** @return Borrowed wrapped stream, or null if the producer reference is cleared */
	@Override
	public final InputStream is(){return thread==null ? null : thread.bf1.is();}
	
	/** @return Wrapped producer counter (possibly ahead of consumption), or -1 after close */
	@Override
	public final long lineNum(){return thread==null ? -1 : thread.bf1.lineNum();}

	/*--------------------------------------------------------------*/
	/*----------------  Fields and constants  ----------------*/
	/*--------------------------------------------------------------*/

	/** Sum of capacities at threshold publication; excludes the final tail and survives reset. */
	long cntr;
	/** Current producer; retained after EOF until close clears it. */
	private BF1Thread thread=null;
	/** Borrowed consumer buffer or the shared poison array. */
	private byte[][] currentList=null;
	/** Next slot in the current consumer buffer. */
	private int currentLoc=0;
//	private int currentSize=0;
	
//	private long numIn=0, numOut=0;
	
	/** One borrowed line returned before the next buffered line. */
	private byte[] pushBack=null;
	
	/** Shared terminal marker, compared by identity. */
	static final byte[][] poison=new byte[0][];
	/** Lifecycle diagnostic switch. */
	public static boolean verbose=false;
	/** Compile-time detailed queue diagnostic switch. */
	private static final boolean verbose2=false;
	/** Maximum line references per ordinary outer buffer. */
	private static final int bufflen=1000;
	/** Number of reusable ordinary outer buffers. */
	private static final int buffs=4;
	/** Soft payload-byte threshold before publishing a buffer. */
	private static final int buffcapacity=256000;
	
	/** Sticky wrapped read/close error flag, preserved across resets. */
	private boolean errorState=false;
	
}
