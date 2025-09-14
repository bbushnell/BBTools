package stream;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.Set;

import cardinality.CardinalityTracker;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.HeapLoc;
import structures.ListNum;
import structures.SetLoc;

/**
 * Timestamps and sorts open buffers prior to retirement.
 * Also uses a heap to open the largest closed buffer.
 * This is the most efficient at minimizing file opens/closes,
 * and seems to be the fastest too.
 * 
 * @author Brian Bushnell
 * @date April 9, 2024
 *
 */
public class MultiCros6 extends BufferedMultiCros {
	
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
		MultiCros6 mcros=new MultiCros6(pattern, null, false, false, false, false, FileFormat.FASTQ, false, 4);
		
		//Create the input stream
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, false, in);
		cris.start();
		
		//Fetch the first list
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		//Process all remaining lists
		while(reads!=null && reads.size()>0){
			
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
	public MultiCros6(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, boolean threaded_, int maxStreams_){
		super(pattern1_, pattern2_, overwrite_, append_, allowSubprocess_, useSharedHeader_, defaultFormat_, threaded_, maxStreams_);
		
		bufferMap=new LinkedHashMap<String, Buffer>();
		streamQueue=new ArrayDeque<String>(maxStreams);
		heap=new HeapLoc<Buffer>(2047, false);
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
		
		//In threaded mode this is called automatically by addToBuffers.
//		if(!threaded) {handleLoad0();}
		//Actually, this is always called by addToBuffers.
	}
	
	@Override
	void handleLoad0() {
		
		//Dump opportunistically because there are free streams
		final int streams=streamQueue.size();
		final int freeStreams=maxStreams-streams;
		if(freeStreams>0 && !heap.isEmpty()) {
			final float mult=((streams+3f)/(maxStreams+2f));
			final long mll=(long)(memLimitLower*mult);
			final int bpb=(int)(bytesPerBuffer*mult);
			
			Buffer b=heap.peek();
			boolean dump=(b.currentBytes>=bpb && bytesInFlight>=mll);
			dump|=(b.currentBytes>=bpb*2);
			if(dump) {
//				System.err.println("Triggered mll dump: "+
//						(b.currentBytes>=bpb && bytesInFlight>=mll)
//						+", "+(b.currentBytes*2>=bpb)+", "+
//						"cb="+b.currentBytes+", bpb="+bpb+", bif="+bytesInFlight+
//						", mll="+mll+", fs="+freeStreams);
				heap.poll();
				b.dump(true);
			}
		}
		
		//TODO: Getting rid of the oldest buffers is also a good idea since that will free up memory for a longer time.
		
		//Dump biggest buffers to free memory
		while(bytesInFlight>=memLimitMid && !heap.isEmpty()) {
			Buffer b=heap.poll();

//			System.err.println("Triggered mlm dump: "+
//					"cb="+b.currentBytes+", bpb="+bytesPerBuffer+", bif="+bytesInFlight+
//					", mlm="+memLimitMid+", fs="+freeStreams);
			b.dump(true);
		}
		
		//Dump everything.
		if(bytesInFlight>=memLimitUpper){
			assert(heap.isEmpty());
			long dumped=dumpAll();
			if(dumped<1 && Shared.EA()){//Dump failed; exit
				KillSwitch.kill("\nThis program ran out of memory."
						+ "\nTry increasing the -Xmx flag or get rid of the minreads flag,"
						+ "\nor disable assertions to skip this message and try anyway.");
			}
		}
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
		assert(heap.isEmpty());
		//Then, retire any active streams
//		while(!streamQueue.isEmpty()){retire(4);}
		if(!streamQueue.isEmpty()){retire(streamQueue.size());}
		assert(streamQueue.isEmpty());
		return x;
	}
	
	@Override
	long dumpAll(){
		if(verbose) {
			System.err.println("before dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight);
		}
		long dumped=0;
		while(!heap.isEmpty()) {
			dumped+=heap.poll().dump(true);
		}
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			dumped+=e.getValue().dump(true);
		}
		if(verbose) {
			System.err.println("after dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight+", dumped="+dumped);
		}
		return dumped;
	}
	
	/** Close the least-recently-used streams */
	private void retire(int retCount){
		if(verbose){System.err.println("Enter retire("+retCount+"); streamQueue="+streamQueue);}
		final long time0=System.nanoTime(), time1, time2, time3, time4;
		retCount=Tools.min(streamQueue.size(), retCount);
		final int sortCount=Tools.min(streamQueue.size(), retCount*2+1, retCount+4);
		
		ArrayList<Buffer> rlist=new ArrayList<Buffer>(sortCount);
		//Select the first name in the queue, which is the least-recently-used.
		for(int i=0; i<sortCount; i++) {
			String name=streamQueue.removeFirst();
			Buffer b=bufferMap.get(name);
			rlist.add(b);
			if(verbose){
				System.err.println("rlist.add("+name+","+b.timestamp+"); ros="+(b.currentRos!=null));
			}
		}
		Collections.sort(rlist, TimestampComparator.instance);
//		Collections.reverse(rlist);//This is to test if any of this logic is even useful
		time1=System.nanoTime();
		for(int i=0; i<sortCount; i++) {
			Buffer b=rlist.get(i);
			if(i<retCount) {
				b.dump(b.currentRos);
				if(verbose){System.err.println("retire("+b.name+","+b.timestamp+")");}
			}else{
				streamQueue.add(b.name);
				rlist.set(i, null);
			}
		}
		
		time2=System.nanoTime();
		for(int i=0; i<retCount; i++) {rlist.get(i).currentRos.close();}
		time3=System.nanoTime();
		for(int i=0; i<retCount; i++) {
			Buffer b=rlist.get(i);
			ConcurrentReadOutputStream ros=b.currentRos;
			ros.join();
			errorState|=(ros.errorState() || !ros.finishedSuccessfully());
			//Delete the pointer to output stream
			b.currentRos=null;
			if(verbose){System.err.println("Exit retire("+b.name+"); ros="+(b.currentRos!=null)+
					", streamQueue="+streamQueue);}
		}
			
		time4=System.nanoTime();
		retireTime1+=(time1-time0);
		retireTime2+=(time2-time1);
		retireTime3+=(time3-time2);
		retireTime4+=(time4-time3);
		retireCount+=retCount;
		retireCalls++;
	}
	
	private void addToHeap(Buffer b) {
		assert(b.loc()<0);
		assert(b.list.size()>0);
		if(!heap.hasRoom()) {
			heap=heap.resizeNew(heap.CAPACITY*2+1);
		}
		heap.add(b);
		assert(b.loc()>=0);
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
		bb.append("Retire Time 3:\t").append(retireTime3*mult, 2).append(" us").nl();
		bb.append("Retire Time 4:\t").append(retireTime4*mult, 2).append(" us").nl();
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
	private class Buffer implements SetLoc<Buffer> {
		
		Buffer(String name_){
			name=name_;
			timestamp=(bufferTimer++);
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
			if(verbose){System.err.println("Made buffer for "+name);}
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
			basesIn+=size;//should be pairlen...
			readsInFlight+=count;
			readsIn+=count;
//			bufferTimer+=8*count;
			if(trackCardinality){loglog.hash(r);}
			
			//Decide whether to dump
			handleLoadB();
		}
		
		/**
		 * Determine whether to dump this buffer based on its current size.
		 * Then, determine whether to dump all buffers based on their combined size.
		 */
		private void handleLoadB(){
			final int size=list.size();
			
			if(currentRos!=null) {
				assert(heapLoc<0);
				if(size>=200 || currentBytes>400000) {
					dump(false);
				}
				return;
			}
			if(heapLoc>=0) {
				heap.jiggleDown(this);
			}else if(currentBytes>=10000 && readsIn>=minReadsToDump) {
				//Add eventually, but don't pollute the heap with tiny buffers
				addToHeap(this);
			}
			
		}
		
		/** Dump buffered reads, creating a stream if needed */
		long dump(boolean force){
			if(list.isEmpty() || readsIn<minReadsToDump){return 0;}
			ConcurrentReadOutputStream ros=getStream(force);
			return ros==null ? 0 : dump(ros);
		}
		
		/** 
		 * Dump buffered reads to the stream.
		 * If the buffer is empty, nothing happens. */
		long dump(final ConcurrentReadOutputStream ros){
			if(verbose){System.err.println("Dumping "+name);}
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
			timestamp=(bufferTimer++);
			return size0;
		}
		
		/** Fetch the stream for this buffer, creating a new one if needed */
		private ConcurrentReadOutputStream getStream(boolean force){
			if(verbose){System.err.println("Enter getStream("+name+"); ros="+(currentRos!=null)+", +streamQueue="+streamQueue);}
			
			if(currentRos!=null){//The stream already exists
//				assert(streamQueue.contains(name));//slow
				if(streamQueue.peekLast()!=name){
					//Move to end to prevent early retirement
					boolean b=streamQueue.remove(name);
					assert(b) : "streamQueue did not contain "+name+", but the ros was open.";
					streamQueue.addLast(name);
				}//TODO: This is no longer necessary
			}else{//The stream does not exist, so create it
				return createStream(force);
			}
			
			assert(currentRos!=null) : "The stream for "+name+" was not created.";
			assert(streamQueue.peekLast()==name) : "The stream for "+name+" was not placed in the queue.";
			if(verbose){System.err.println("Exit  getStream("+name+"); ros="+(currentRos!=null)+", +streamQueue="+streamQueue);}
			return currentRos;
		}
		
		/** Create a stream for this buffer, and stick it in the queue */
		private ConcurrentReadOutputStream createStream(boolean force){
			if(!force && streamQueue.size()>=maxStreams && bytesInFlight<memLimitLower) {return null;}
			assert(currentRos==null) : "This should never be called if there is an existing stream.";
			if(!deleted) {
				if(numDumps==0 && overwrite){
					//First time, an existing file must be deleted first, because the ff is set to append mode
					if(verbose){System.err.println("Deleting "+name+" ; exists? "+ff1.exists());}
					delete(ff1);
					delete(ff2);
				}
				deleted=true;
			}
			final long time0=System.nanoTime(), time1, time2, time3, time4, time5;
			
//			assert(!streamQueue.contains(name));//slow
			
			//Active streams should never exceed maxStreams
			assert(streamQueue.size()<=maxStreams) : "Too many streams: "+streamQueue+", "+maxStreams;
			if(streamQueue.size()>=maxStreams){
				//Too many open streams; retire one.
				retire(streamsToRetire);
			}
			//After retirement, open streams must be less than maxStreams
			assert(streamQueue.size()<maxStreams) : "Too many streams: "+streamQueue+", "+maxStreams;
			time1=System.nanoTime();
			time2=time1;
			
			//Create a stream
			currentRos=ConcurrentReadOutputStream.getStream(ff1, ff2, rswBuffers, null, useSharedHeader && numDumps==0);
			time3=System.nanoTime();
			currentRos.start();
			if(verbose){System.err.println("Created ros "+name+"; ow="+ff1.overwrite()+", append="+ff1.append());}
			time4=System.nanoTime();
			
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
			return appendTo(new ByteBuilder()).toString();
		}
		
		@Override
		public int compareTo(Buffer b) {
			long dif=b.currentBytes-currentBytes;
			//Using the timestamp does reduce retire count, but only by ~1% on a shuffled file.
			//This might be more noticable on a non-shuffled file; it's worth testing.
//			long dif=b.currentBytes-currentBytes+1024*(timestamp-b.timestamp);
//			long dif=((b.currentBytes*(bufferTimer-b.timestamp))-(currentBytes*(bufferTimer-timestamp)));
//			long dif=(long)(((b.currentBytes*Math.sqrt(bufferTimer-b.timestamp))-
//					(currentBytes*Math.sqrt(bufferTimer-timestamp))));
			return dif<0 ? -1 : dif>0 ? 1 : 0;
		}
		
		@Override
		public void setLoc(int newLoc) {
			heapLoc=newLoc;
		}

		@Override
		public int loc() {
			return heapLoc;
		}
		
		/** Stream name, which is the variable part of the file pattern */
		private final String name;
		/** Output file 1 */
		private final FileFormat ff1;
		/** Output file 2 */
		private final FileFormat ff2;
		
		/** Once created, the stream sticks around to be re-used unless it is retired. */
		private ConcurrentReadOutputStream currentRos;
		
		/** Current list of buffered reads */
		private ArrayList<Read> list;
		
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
		/** Whether the existing files have been checked or deleted yet */
		private boolean deleted=false;
		/** Time of last dump */
		private long timestamp=-1;
		/** Location in heap; -1 means not in heap */
		private int heapLoc=-1;
		/** Optional, for tracking cardinality */
		private CardinalityTracker loglog;
		
	}
	
	private static final class TimestampComparator implements Comparator<Buffer>{

		private TimestampComparator() {}
		
		@Override
		public final int compare(Buffer a, Buffer b) {
			assert(a.timestamp!=b.timestamp);
			return a.timestamp<b.timestamp ? -1 : 1;
		}
		
		static final TimestampComparator instance=new TimestampComparator();
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Essentially the number of dumps.  Does not distinguish by dump size. */
	private long bufferTimer=0;
	
	private HeapLoc<Buffer> heap;
	
	/** Open stream names */
	private final ArrayDeque<String> streamQueue;
	
	/** Map of names to buffers */
	public final LinkedHashMap<String, Buffer> bufferMap;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

}
