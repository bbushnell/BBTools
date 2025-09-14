package stream;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;

import cardinality.CardinalityTracker;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Uses a retire queue with independent threads.
 * 
 * @author Brian Bushnell
 * @date April 5, 2024
 *
 */
public class MultiCros4 extends BufferedMultiCros {
	
	/** 
	 * For testing.<br>
	 * args should be:
	 * {input file, output pattern, names...}
	 */
	public static void main(String[] args){
		
		//Do some basic parsing
		String in=args[0];
		String pattern=args[1];
		ArrayList<String> names=new ArrayList<String>();
		for(int i=2; i<args.length; i++){names.add(args[i]);}
		
		//Create an mcros for this pattern
		MultiCros4 mcros=new MultiCros4(pattern, null, false, false, false, false, FileFormat.FASTQ, false, 4);
		
		//Create the input stream
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, false, in);
		cris.start();
		
		//Fetch the first list
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		//Process all remaining lists
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			
			//Add the reads by barcode
			for(Read r1 : reads){
				mcros.add(r1, r1.barcode(true));
			}
			
			//Return the old list and get a new one
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		
		//Close the streams
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
		mcros.close();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** @See Details in superclass constructor */
	public MultiCros4(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, 
			boolean useSharedHeader_, int defaultFormat_, boolean threaded_, int maxStreams_){
		super(pattern1_, pattern2_, overwrite_, append_, allowSubprocess_, 
				useSharedHeader_, defaultFormat_, threaded_, maxStreams_);
		
		bufferMap=new LinkedHashMap<String, Buffer>();
		streamQueue=new ArrayDeque<String>(maxStreams);
		freeTokens=new ArrayBlockingQueue<Token>(maxStreams);
		for(int i=0; i<maxStreams; i++) {freeTokens.add(new Token(i));}
		maxRetireThreads=Tools.max(1, (int)(maxStreams*0.25f));
		retireQueue=new ArrayBlockingQueue<Buffer>(maxRetireThreads);
		retireThreads=new ArrayList<RetireThread>(maxRetireThreads);
		maxOpenStreams=Tools.max(1, maxStreams-maxRetireThreads);
		for(int i=0; i<8 && i<maxRetireThreads; i++) {
			RetireThread rt=new RetireThread();
			rt.start();
			retireThreads.add(rt);
		}
		System.err.println("maxStreams="+maxStreams+", maxRetireThreads="+
				maxRetireThreads+", maxOpenStreams="+maxOpenStreams);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public boolean finishedSuccessfully(){
		return !errorState;
	}
	
	@Override
	public void add(Read r, String name){
		Buffer b=bufferMap.get(name);
		if(b==null){
			b=new Buffer(name);
			bufferMap.put(name, b);
			//Note: I could adjust bytesPerBuffer threshold here in response to the number of buffers.
		}
		b.add(r);
	}
	
	@Override
	public long dumpResidual(ConcurrentReadOutputStream rosu){
		//For each Buffer, check if it contains residual reads
		//If so, dump it into the stream
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			Buffer b=e.getValue();
			assert((b.readsIn<minReadsToDump) == (b.list!=null && !b.list.isEmpty()));
			if(b.readsIn>0 && b.readsIn<minReadsToDump){
				assert(b.list!=null && !b.list.isEmpty());
				residualReads+=b.readsIn;
				residualBases+=b.basesIn;
				if(rosu!=null){rosu.add(b.list, 0);}
			}
			b.list=null;
		}
		return residualReads;
	}
	
	@Override
	public ByteBuilder report(){
		ByteBuilder bb=new ByteBuilder(1024);
		
		//Add a line for residual reads dumped 
		if(minReadsToDump>0){
			bb.append("Residual").tab().append(residualReads).tab().append(residualBases).nl();
		}
		
		//Add a line for each Buffer
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			Buffer buffer=e.getValue();
			if(buffer.numDumps>0){//Only add a line if this Buffer actually created a file
				buffer.appendTo(bb);
			}
		}
		return bb;
	}
	
	@Override
	public Set<String> getKeys(){return bufferMap.keySet();}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	long closeInner() {
		//First dump everything
		final long x=dumpAll();
		//Then, retire any active streams
		while(!streamQueue.isEmpty()){retire(1);}
		try {
			for(boolean success=false; !success;) {
				retireQueue.put(POISON_BUFFER);
				success=true;
			}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		waitForFinishInner();
		return x;
	}
	
	/** Wait for this object's thread to terminate */
	public final void waitForFinishInner(){
		if(verbose){System.err.println("Waiting for finish.");}
		for(RetireThread rt : retireThreads) {
			synchronized(rt) {
				while(rt.getState()!=Thread.State.TERMINATED){
					if(verbose){System.err.println("Attempting join: state="+rt.getState());}
					try {
						rt.join(1000);
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
			}
		}
	}
	
	@Override
	long dumpAll(){
		if(verbose) {
			System.err.println("before dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight);
		}else {
//			System.err.println("dumpAll triggered due to memory pressure.");//Also happens at the end
		}
		long dumped=0;
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			dumped+=e.getValue().dump();
		}
		if(verbose) {
			System.err.println("after dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight+", dumped="+dumped);
		}
		return dumped;
	}
	
	/** Close the least-recently-used stream */
	private void retire(int retCount){
		if(verbose){System.err.println("Enter retire(); streamQueue="+streamQueue);}
		assert(retCount==1); //For now
		final long time0=System.nanoTime(), time1, time2, time3, time4;
		
		//Select the first name in the queue, which is the least-recently-used.
		String name=streamQueue.removeFirst();
		Buffer b=bufferMap.get(name);
		
		if(verbose){
			System.err.println("retire("+name+"); state="+b.getState()+", streamQueue="+streamQueue);
		}
		assert(b!=null);
		assert(b.getState()==OPEN) : name+", "+b.getState()+"\n"+streamQueue+"\n";
		b.dump();
		b.setState(RETIRING);
		time1=System.nanoTime();
		
//		if(retireThreads.size()<maxRetireThreads && retireQueue.size()>0) {
//			RetireThread rt=new RetireThread();
//			rt.start();
//			retireThreads.add(rt);
//		}
		
		try {
			for(boolean success=false; !success; ) {
				retireQueue.put(b);
				success=true;
			}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		time2=System.nanoTime();
		if(verbose) {System.err.println("Added "+name+" to retire queue.");}
		retireTime1+=(time1-time0);
		retireTime2+=(time2-time1);
//		retireTime3+=(time3-time2);
//		retireTime4+=(time4-time3);
		retireCount+=retCount;
		retireCalls++;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Profiling           ----------------*/
	/*--------------------------------------------------------------*/
	
	public String printRetireTime() {
		ByteBuilder bb=new ByteBuilder();
		float mult=0.001f/retireCount;
		bb.append("Max Streams:\t").append(maxStreams).nl();
		bb.append("Retire Count:\t").append(retireCount).nl();
		bb.append("Retires Per Call:\t").append(retireCount/(float)retireCalls, 2).nl();
		bb.append("streamsToRetire:\t").append(streamsToRetire).nl();
		bb.append("Retire Time 1:\t").append(retireTime1*mult, 2).append(" us").nl();
		bb.append("Retire Time 2:\t").append(retireTime2*mult, 2).append(" us").nl();
//		bb.append("Retire Time 3:\t").append(retireTime3*mult, 2).append(" us").nl();
//		bb.append("Retire Time 4:\t").append(retireTime4*mult, 2).append(" us").nl();
		bb.append("Total:        \t").append((retireTime1+retireTime2+retireTime3+retireTime4)/1000000000.0, 3).append(" s").nl();
		return bb.toString();
	}

	private long retireTime1=0;
	private long retireTime2=0;
	private long retireTime3=0;
	private long retireTime4=0;
	private long retireCount=0;
	private long retireCalls=0;
	
	public String printCreateTime() {
		ByteBuilder bb=new ByteBuilder();
		float mult=0.001f/retireCount;
		bb.append("Create Time 1:\t").append(createTime1*mult, 2).append(" us").nl();
		bb.append("Create Time 2:\t").append(createTime2*mult, 2).append(" us").nl();
		bb.append("Create Time 3:\t").append(createTime3*mult, 2).append(" us").nl();
		bb.append("Create Time 4:\t").append(createTime4*mult, 2).append(" us").nl();
		bb.append("Create Time 5:\t").append(createTime5*mult, 2).append(" us").nl();
		bb.append("Total:        \t").append((createTime2+createTime3+createTime4+createTime5)/1000000000.0, 3).append(" s").nl();
		return bb.toString();
	}
	
	private long createTime1=0;
	private long createTime2=0;
	private long createTime3=0;
	private long createTime4=0;
	private long createTime5=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * A Buffer holds reads destined for to a specific file.
	 * When sufficient reads are present, it opens a stream and writes them.
	 * If too many streams are open, it closes another stream first.
	 */
	private class Buffer {
		
		Buffer(String name_){
			name=name_;
			String s1=pattern1.replaceFirst("%", name);
			String s2=pattern2==null ? null : pattern2.replaceFirst("%", name);
			
			//These are created with overwrite=false append=true because 
			//the files will be appended to if the stream gets prematurely retired.
			//Therefore, files must be explicitly deleted first.
			//Alternative would be to create a new FileFormat each time.
			ff1=FileFormat.testOutput(s1, defaultFormat, null, allowSubprocess, false, true, false);
			ff2=FileFormat.testOutput(s2, defaultFormat, null, allowSubprocess, false, true, false);
			
			list=new ArrayList<Read>(readsPerBuffer);
			if(trackCardinality){loglog=CardinalityTracker.makeTracker();}
//			if(verbose){System.err.println("Made buffer for "+name);}
		}
		
		/** 
		 * Add a read to this buffer, and update all the tracking variables.
		 * This may trigger a dump.
		 */
		void add(Read r){
			//Add the read
			list.add(r);
			
			//Gather statistics
			long size=r.countPairBytes();
			int count=r.pairCount();
			currentBytes+=size;
			bytesInFlight+=size;
			basesIn+=r.pairLength();
			readsInFlight+=count;
			readsIn+=count;
			if(trackCardinality){loglog.hash(r);}
			
			//Decide whether to dump
			handleLoad();
		}
		
		/**
		 * Determine whether to dump this buffer based on its current size.
		 * Then, determine whether to dump all buffers based on their combined size.
		 */
		private void handleLoad(){
			//3rd term allows preemptive dumping
			//More generally, this triggers a dump if the reads in this buffer exceed
			//the maximum allowed reads or bytes 
			final int size=list.size();
			if(size>=readsPerBuffer || currentBytes>=bytesPerBuffer) {
				if(verbose){
					System.err.println("list.size="+list.size()+"/"+readsPerBuffer+
							", bytes="+currentBytes+"/"+bytesPerBuffer+", bytesInFlight="+bytesInFlight+"/"+memLimitUpper);
				}
				dump();
			}else if((size>=200 || currentBytes>=400000) && state==OPEN) {
				assert(getState()==OPEN);//Synchronization should not be needed here 
				if(verbose){
					//					System.err.println("list.size="+list.size()+"/"+readsPerBuffer+
					//							", bytes="+currentBytes+"/"+bytesPerBuffer+", bytesInFlight="+bytesInFlight+"/"+memLimit);
				}
				dump();
			}
			
			//Too much buffered data in ALL buffers; dump everything.
			if(bytesInFlight>=memLimitUpper){
				long dumped=dumpAll();
				if(dumped<1 && Shared.EA()){//Dump failed; exit
					KillSwitch.kill("\nThis program ran out of memory."
							+ "\nTry increasing the -Xmx flag or get rid of the minreads flag,"
							+ "\nor disable assertions to skip this message and try anyway.");
				}
			}
		}
		
		/** Dump buffered reads, creating a stream if needed */
		long dump(){
			if(list.isEmpty() || readsIn<minReadsToDump){return 0;}
			final int state=getState();
			if(state==RETIRING) {return 0;}
			ConcurrentReadOutputStream ros=getStream();
			return dump(ros);
		}
		
		/** 
		 * Dump buffered reads to the stream.
		 * If the buffer is empty, nothing happens. */
		long dump(final ConcurrentReadOutputStream ros){
			if(verbose && list.size()>400){System.err.println("Dumping "+name);}
			if(list.isEmpty()){return 0;}
			final long size0=list.size();
			
			//Send the list to the output stream
			//This part is async
			ros.add(list, numDumps);
			
			//Create a new list, since the old one is busy
			list=new ArrayList<Read>(400);
			
			//Manage statistics
			bytesInFlight-=currentBytes;
			readsInFlight-=size0;
			readsWritten+=size0;
			currentBytes=0;
			numDumps++;
			return size0;
		}
		
		/** Fetch the stream for this buffer, creating a new one if needed */
		private synchronized ConcurrentReadOutputStream getStream(){
			if(state==OPEN) {return currentRos;}
			else if(state==CLOSED) {assert(currentRos==null);}
			else {assert(false);}
			if(verbose){System.err.println("Enter getStream("+name+"); ros="+(currentRos!=null)+", +streamQueue="+streamQueue);}
			
			if(currentRos!=null){//The stream already exists
//				assert(streamQueue.contains(name));//slow
				if(streamQueue.peekLast()!=name){
					//Move to end to prevent early retire
					boolean b=streamQueue.remove(name);
					assert(b) : "streamQueue did not contain "+name+", but the ros was open.";
					streamQueue.addLast(name);
				}
			}else{//The stream does not exist, so create it
				createStream();
			}
			
			assert(currentRos!=null) : "The stream for "+name+" was not created.";
			assert(streamQueue.peekLast()==name) : "The stream for "+name+" was not placed in the queue.";
			if(verbose){System.err.println("Exit  getStream("+name+"); ros="+(currentRos!=null)+", +streamQueue="+streamQueue);}
			return currentRos;
		}
		
		/** Create a stream for this buffer, and stick it in the queue */
		private synchronized ConcurrentReadOutputStream createStream(){
			assert(state==CLOSED);
			assert(currentRos==null) : "This should never be called if there is an existing stream.";
			if(numDumps==0 && overwrite){
				//First time, an existing file must be deleted first, because the ff is set to append mode
				if(verbose){System.err.println("Deleting "+name+" ; exists? "+ff1.exists());}
				delete(ff1);
				delete(ff2);
			}
			final long time0=System.nanoTime(), time1, time2, time3, time4, time5;
			
//			assert(!streamQueue.contains(name));//slow
			
			//Active streams should never exceed maxStreams
			assert(streamQueue.size()<=maxOpenStreams) : "Too many streams: "+streamQueue+", "+maxStreams;
			if(streamQueue.size()>=maxOpenStreams){
				//Too many open streams; retire one.
				retire(1);
			}
			time1=System.nanoTime();
			//After retire, open streams must be less than maxStreams
			assert(streamQueue.size()<maxOpenStreams) : "Too many streams: "+streamQueue+", "+maxStreams;
			assert(token==null);
			Token t=null;
			if(verbose) {System.err.println("Fetching token for "+name);}
			while(t==null) {
				try {
					t=freeTokens.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			time2=System.nanoTime();
			if(verbose) {System.err.println("Got token for "+name);}
			giveToken(t);
			
			//Create a stream
			currentRos=ConcurrentReadOutputStream.getStream(ff1, ff2, rswBuffers, null, useSharedHeader && numDumps==0);
			time3=System.nanoTime();
			currentRos.start();
			time4=System.nanoTime();
			
			if(verbose){System.err.println("Created ros "+name+"; ow="+ff1.overwrite()+", append="+ff1.append());}
			setState(OPEN);
			
			//Add it to the queue
			streamQueue.addLast(name);
			time5=System.nanoTime();

			createTime1+=(time1-time0);
			createTime2+=(time2-time1);
			createTime3+=(time3-time2);
			createTime4+=(time4-time3);
			createTime5+=(time5-time4);
			return currentRos;
		}
		
		/** Delete this file if it exists */
		private void delete(FileFormat ff){
			if(ff==null){return;}
			assert(overwrite || !ff.exists()) : "Trying to delete file "+ff.name()+", but overwrite=f.  Please add the flag overwrite=t.";
			ff.deleteIfPresent();
		}
		
		/** 
		 * Format this buffer's summary as a line of text.
		 * @param bb ByteBuilder to append the text
		 * @return The modified ByteBuilder
		 */
		ByteBuilder appendTo(ByteBuilder bb) {
			bb.append(name).tab().append(readsIn).tab().append(basesIn);
			if(trackCardinality){bb.tab().append(loglog.cardinality());}
			return bb.nl();
		}
		
		@Override
		public String toString(){
			if(verbose) {return toString2();}
			return appendTo(new ByteBuilder()).toString();
		}
		
		public String toString2(){
			return name+" "+state+" "+(currentRos!=null);
		}
		
		synchronized void setState(int newState) {
			if(verbose){
				System.err.println("setState "+name+" "+state+" -> "+newState+"; ros="+(currentRos!=null));
			}
			assert(currentRos!=null);
			assert(token!=null);
			
			int x=(state+1)%3;
			assert(state!=newState);
			assert(x==newState);
			state=newState;
			if(state==CLOSED) {currentRos=null;}
			else if(state==RETIRING) {assert(list.isEmpty());}
		}
		
		synchronized int getState() {
			return state;
		}
		
		synchronized void giveToken(Token t) {
			assert(token==null);
			assert(state==CLOSED);
			assert(t!=null);
			token=t;
		}
		
		synchronized Token takeToken() {
			assert(token!=null);
			assert(state==CLOSED) : state;
			Token t=token;
			token=null;
			return t;
		}
		
		/** Stream name, which is the variable part of the file pattern */
		private final String name;
		/** Output file 1 */
		private final FileFormat ff1;
		/** Output file 2 */
		private final FileFormat ff2;
		
		/** Once created, the stream sticks around to be re-used unless it is retired. */
		private ConcurrentReadOutputStream currentRos;
		
		
		private int state=CLOSED;
		private Token token;
		
		/** Current list of buffered reads */
		private ArrayList<Read> list;
		/** List of buffered reads to dump */
		private ArrayList<Read> dumpList;
		
		/** Number of reads entering the buffer */
		private long readsIn=0;
		/** Number of bases entering the buffer */
		private long basesIn=0;
		/** Number of reads written to disk */
		@SuppressWarnings("unused")
		private long readsWritten=0;//This does not count read2!
		/** Number of bytes currently in this buffer (estimated) */
		private long currentBytes=0;
		/** Number of dumps executed */
		private long numDumps=0;
		/** Optional, for tracking cardinality */
		private CardinalityTracker loglog;
		
	}
	
	private class RetireThread extends Thread {
		
		@Override
		public void run() {
			while(true) {
				Buffer b=null;
				if(verbose) {System.err.println("Retire thread fetching buffer; retQueue="+retireQueue);}
				while(b==null) {
					try {
						b=retireQueue.take();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				if(verbose) {System.err.println("Retire thread fetched "+b.name+"; retQueue="+retireQueue);}
				if(b==POISON_BUFFER) {
					for(boolean success=false; success==false;) {
						try {
							retireQueue.put(b);
							success=true;
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					if(verbose) {System.err.println("Retire thread terminated.");}
					return;
				}
				retire(b);
			}
		}
		
		void retire(Buffer b) {
			if(verbose) {System.err.println("Retire thread retiring "+b.name+".");}
			
			assert(b.currentRos!=null) : b.state+", "+b.name;
			errorState=ReadWrite.closeStream(b.currentRos) | errorState;//Traditional synchronous close-and-wait
			
//			ros.close();
//			ros.join();
//			errorState|=(ros.errorState() || !ros.finishedSuccessfully());
			if(verbose){System.err.println("Exit retire("+b.name+"); ros="+(b.currentRos!=null)+
					", streamQueue="+streamQueue);}
			b.setState(CLOSED);//Delete the pointer to output stream
			Token t=b.takeToken();
			freeTokens.add(t);

			if(verbose) {System.err.println("Retire thread retired "+b.name+".");}
		}
		
	}
	
	private static class Token {
		Token(int id_){id=id_;}
		final int id;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	private final int maxRetireThreads;

	/** Allow this many open streams */
	public final int maxOpenStreams;
	
	/** Open stream names */
	private final ArrayDeque<String> streamQueue;
	
	/** Tokens available for use */
	private final ArrayBlockingQueue<Token> freeTokens;
	
	/** Buffers waiting to retire */
	private final ArrayBlockingQueue<Buffer> retireQueue;
	
	/** Buffers waiting to retire */
	private final ArrayList<RetireThread> retireThreads;
	
	/** Map of names to buffers */
	public final LinkedHashMap<String, Buffer> bufferMap;
	
	private final Buffer POISON_BUFFER=new Buffer("POISON_BUFFER_NOT_A_FILE");
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final int CLOSED=0, OPEN=1, RETIRING=2;

}
