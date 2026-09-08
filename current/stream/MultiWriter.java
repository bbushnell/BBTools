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
 * A modern, Writer-based successor to MultiCros6.
 * Structurally identical to MultiCros6 (timestamp-sorted LRU retirement + a size heap for
 * opportunistic dumps), but each per-name output is a lightweight single-threaded {@link Writer}
 * obtained from {@link WriterFactory} with threads=1 (ST2 writers in threaded mode: one background
 * format+write thread per open stream, so dumps stay asynchronous like MultiCros6's CROS path but
 * with far less machinery), NOT a ConcurrentReadOutputStream. The 4-phase close+join of MultiCros6's
 * retire() collapses to a single poisonAndWait() per writer (drain queue, join thread, close stream,
 * return errorState). Retirement is still essential: there can be tens of thousands of
 * output files, so only maxStreams writers are kept open at once; a retired name's file is reopened in
 * append mode on reactivation (concatenated gzip members, exactly as MultiCros6 relies on).
 * Measured on an 800k-read 48-group demux (streams=8): wall-time parity with MultiCros6,
 * byte-identical outputs across unpaired, twin-file paired, gzipped, nonthreaded, and minreads modes.
 *
 * Registered as mcrostype=7 in BufferedMultiCros.make(). It is a BufferedMultiCros subclass for
 * backwards-compatible drop-in into NovaDemux/DemuxByName2 even though it never uses CROS.
 *
 * @author Brian Bushnell
 * @contributor Furina
 * @date June 25, 2026
 */
public class MultiWriter extends BufferedMultiCros {

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
		MultiWriter mcros=new MultiWriter(pattern, null, false, false, false, false, FileFormat.FASTQ, false, 4);

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
	public MultiWriter(String pattern1_, String pattern2_,
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
		}
		b.add(r);
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
				heap.poll();
				b.dump(true);
			}
		}

		//Dump biggest buffers to free memory
		while(bytesInFlight>=memLimitMid && !heap.isEmpty()) {
			Buffer b=heap.poll();
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
		//If so, dump it into the stream. rosu is the CALLER's unmatched-reads catch-all (a CROS),
		//not a per-name MultiWriter stream - forwarded as-is, exactly like MultiCros6.
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
		//NOTE (not a bug - do not "simplify"): buffers that were in the heap are visited AGAIN here (they
		//remain in bufferMap), so dump(true) hits them twice - but the second call short-circuits on the
		//list.isEmpty() guard in dump() (no double-write, no double-count). This bufferMap pass exists to
		//also catch buffers that were never heaped (the small ones, below the addToHeap threshold).
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
		//Select the first names in the queue, which are the least-recently-used.
		for(int i=0; i<sortCount; i++) {
			String name=streamQueue.removeFirst();
			Buffer b=bufferMap.get(name);
			rlist.add(b);
			if(verbose){
				System.err.println("rlist.add("+name+","+b.timestamp+"); writer="+(b.currentWriter!=null));
			}
		}
		Collections.sort(rlist, TimestampComparator.instance);
		time1=System.nanoTime();
		//Flush remaining buffered reads to the (still-open) writers for the retCount oldest; re-queue the rest.
		for(int i=0; i<sortCount; i++) {
			Buffer b=rlist.get(i);
			if(i<retCount) {
				b.dump(b.currentWriter);
				if(verbose){System.err.println("retire("+b.name+","+b.timestamp+")");}
			}else{
				streamQueue.add(b.name);
				rlist.set(i, null);
			}
		}
		time2=System.nanoTime();
		//Close each retiring writer. poisonAndWait() drains the writer's queue, joins its worker
		//thread, closes the stream, and returns errorState (true=error) - so the 4-phase CROS
		//close()+join()+fold collapses to this one call.
		for(int i=0; i<retCount; i++) {
			Buffer b=rlist.get(i);
			Writer w=b.currentWriter;
			errorState|=w.poisonAndWait();
			b.currentWriter=null;
			if(verbose){System.err.println("Exit retire("+b.name+"); writer="+(b.currentWriter!=null)+
					", streamQueue="+streamQueue);}
		}
		time3=System.nanoTime();
		time4=time3;//Phase 4 (separate join) is folded into phase 3 for ZT writers.
		retireTime1+=(time1-time0);
		retireTime2+=(time2-time1);
		retireTime3+=(time3-time2);
		retireTime4+=(time4-time3);
		retireCount+=retCount;
		retireCalls++;
	}

	/**
	 * Adds buffer to the priority heap for size-based processing.
	 * Resizes heap if necessary to accommodate new buffer.
	 * @param b Buffer to add to heap (must not already be in heap)
	 */
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

	@Override
	public String printRetireTime() {
		ByteBuilder bb=new ByteBuilder();
		float mult=0.001f/Tools.max(1, retireCount);//guard: /retireCount = Infinity when retireCount==0 (no streams retired). Profiling-only.
		bb.append("Max Streams:\t").append(maxStreams).nl();
		bb.append("Retire Count:\t").append(retireCount).nl();
		bb.append("Retires Per Call:\t").append(retireCount/(float)Tools.max(1, retireCalls), 2).nl();
		bb.append("streamsToRetire:\t").append(streamsToRetire).nl();
		bb.append("Retire Time 1:\t").append(retireTime1*mult, 2).append(" us").nl();
		bb.append("Retire Time 2:\t").append(retireTime2*mult, 2).append(" us").nl();
		bb.append("Retire Time 3:\t").append(retireTime3*mult, 2).append(" us").nl();
		bb.append("Total:        \t").append((retireTime1+retireTime2+retireTime3+retireTime4)/1000000000.0, 3).append(" s").nl();
		return bb.toString();
	}

	private long retireTime1=0;
	private long retireTime2=0;
	private long retireTime3=0;
	private long retireTime4=0;
	private long retireCount=0;
	private long retireCalls=0;

	@Override
	public String printCreateTime() {
		ByteBuilder bb=new ByteBuilder();
		float mult=0.001f/Tools.max(1, retireCount);//guard: /retireCount = Infinity when retireCount==0. Profiling-only (approximate normalizer; no createCount field).
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
	 * A Buffer holds reads destined for a specific file.
	 * When sufficient reads are present, it opens a Writer and dumps them to it.
	 * If too many streams are open, it closes another stream first.
	 */
	private class Buffer implements SetLoc<Buffer> {

		/**
		 * Constructs buffer for specified output name.
		 * Files configured for append mode to handle retirement and recreation (the first open is
		 * preceded by an explicit delete via the deleted flag when overwrite is set).
		 * @param name_ Buffer identifier used in file pattern substitution
		 */
		Buffer(String name_){
			name=name_;
			timestamp=(bufferTimer++);
			String s1=pattern1.replaceFirst("%", name);
			String s2=pattern2==null ? null : pattern2.replaceFirst("%", name);

			//Created overwrite=false append=true because the files are appended to if a stream is
			//prematurely retired then reopened; the first-time truncate is the explicit delete below.
			//ordered=true is REQUIRED: WriterFactory.makeWriter asserts ordered() for twin files, and
			//an ordered ZT writer's JobQueue gives deterministic R1/R2 pairing from this single caller thread.
			ff1=FileFormat.testOutput(s1, defaultFormat, null, allowSubprocess, false, true, true);
			ff2=FileFormat.testOutput(s2, defaultFormat, null, allowSubprocess, false, true, true);

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
			basesIn+=r.pairLength();//pairLength (bases), NOT countPairBytes: residualBases feeds DemuxByName2's Bases Out subtraction (the bytes version went negative; replicated 2026-09-05)
			readsInFlight+=count;
			readsIn+=count;
			if(trackCardinality){loglog.hash(r);}

			//Decide whether to dump
			handleLoadB();
		}

		/**
		 * Determine whether to dump this buffer based on its current size,
		 * and manage its position in the size heap.
		 */
		private void handleLoadB(){
			final int size=list.size();

			if(currentWriter!=null) {
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
			Writer ros=getStream(force);
			return ros==null ? 0 : dump(ros);
		}

		/**
		 * Dump buffered reads to the writer (asynchronous handoff to its worker thread).
		 * If the buffer is empty, nothing happens. */
		long dump(final Writer ros){
			if(verbose){System.err.println("Dumping "+name);}
			if(list.isEmpty()){return 0;}
			final long size0=list.size();

			//Hand the list to the writer; its worker thread formats and writes it.
			//Ids must be a dense ascending sequence STARTING AT 0 for each writer instance (JobQueue
			//firstID=0), so nextJobId resets in createStream; numDumps spans retirements and cannot be used.
			ros.add(list, nextJobId);
			nextJobId++;

			//Create a new list, since the old one is now owned by the writer/output.
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

		/** Fetch the writer for this buffer, creating a new one if needed */
		private Writer getStream(boolean force){
			if(verbose){System.err.println("Enter getStream("+name+"); writer="+(currentWriter!=null)+", +streamQueue="+streamQueue);}

			if(currentWriter!=null){//The stream already exists
				if(streamQueue.peekLast()!=name){
					//Move to end to prevent early retirement
					boolean b=streamQueue.remove(name);
					assert(b) : "streamQueue did not contain "+name+", but the writer was open.";
					streamQueue.addLast(name);
				}
			}else{//The stream does not exist, so create it
				return createStream(force);
			}

			assert(currentWriter!=null) : "The stream for "+name+" was not created.";
			assert(streamQueue.peekLast()==name) : "The stream for "+name+" was not placed in the queue.";
			if(verbose){System.err.println("Exit  getStream("+name+"); writer="+(currentWriter!=null)+", +streamQueue="+streamQueue);}
			return currentWriter;
		}

		/** Create a ZT writer for this buffer, and stick it in the queue */
		private Writer createStream(boolean force){
			if(!force && streamQueue.size()>=maxStreams && bytesInFlight<memLimitLower) {return null;}
			assert(currentWriter==null) : "This should never be called if there is an existing stream.";
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

			//Create a lightweight writer via the drop-in factory entry point. threads=1 selects the
			//ST2 writers in THREADED mode (one background format+write thread per open stream), which
			//keeps dumps asynchronous like MultiCros6's CROS path; threads=0 (synchronous ZT) was
			//measured ~60% slower here because every dump then formats+writes on this thread.
			currentWriter=WriterFactory.getStream(ff1, ff2, rswBuffers, null, useSharedHeader && numDumps==0, 1);
			nextJobId=0;//Fresh writer, fresh JobQueue: ids restart at 0
			time3=System.nanoTime();
			currentWriter.start();
			if(verbose){System.err.println("Created writer "+name+"; ow="+ff1.overwrite()+", append="+ff1.append());}
			time4=System.nanoTime();

			//Add it to the queue
			streamQueue.addLast(name);
			time5=System.nanoTime();

			createTime1+=(time1-time0);
			createTime2+=(time2-time1);
			createTime3+=(time3-time2);
			createTime4+=(time4-time3);
			createTime5+=(time5-time4);
			return currentWriter;
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

		/** Once created, the writer sticks around to be re-used unless it is retired. */
		private Writer currentWriter;

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
		/** Next job id for the current writer; dense from 0 per writer instance */
		private long nextJobId=0;
		/** Whether the existing files have been checked or deleted yet */
		private boolean deleted=false;
		/** Time of last dump */
		private long timestamp=-1;
		/** Location in heap; -1 means not in heap */
		private int heapLoc=-1;
		/** Optional, for tracking cardinality */
		private CardinalityTracker loglog;

	}

	/**
	 * Comparator for sorting buffers by timestamp for retirement ordering.
	 * Ensures oldest streams are retired first to maintain LRU behavior.
	 */
	private static final class TimestampComparator implements Comparator<Buffer>{

		private TimestampComparator() {}

		@Override
		public final int compare(Buffer a, Buffer b) {
			assert(a.timestamp!=b.timestamp);
			return a.timestamp<b.timestamp ? -1 : 1;
		}

		/** Singleton instance of timestamp comparator */
		static final TimestampComparator instance=new TimestampComparator();

	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	/** Essentially the number of dumps.  Does not distinguish by dump size. */
	private long bufferTimer=0;

	/** Priority heap containing buffers ordered by size for dump prioritization */
	private HeapLoc<Buffer> heap;

	/** Open stream names */
	private final ArrayDeque<String> streamQueue;

	/** Map of names to buffers */
	public final LinkedHashMap<String, Buffer> bufferMap;

}
