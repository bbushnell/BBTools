package stream;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
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
import structures.ListNum;

/**
 * Buffers named reads and retires a subset of open outputs by logical dump recency.
 * Retirement sorts a prefix of the active queue, not every open destination.
 * The stream limit counts logical destinations; each can own separate mate files.
 * Direct routing requires one owner; threaded mode uses the inherited transfer thread.
 * Read payloads remain borrowed until downstream output completes.
 *
 * @author Brian Bushnell
 * @date April 8, 2024
 */
public class MultiCros5 extends BufferedMultiCros{
	
	/**
	 * For testing.
	 * Creates a MultiCros5 instance and processes reads from input file,
	 * distributing them by barcode to separate output files.
	 * @param args Positional input/output pattern; trailing names are parsed but unused
	 */
	public static void main(String[] args){
		
		//Do some basic parsing
		String in=args[0];
		String pattern=args[1];
		ArrayList<String> names=new ArrayList<String>();
		for(int i=2; i<args.length; i++){names.add(args[i]);}
		
		//Create an mcros for this pattern
		MultiCros5 mcros=new MultiCros5(pattern, null, false, false, false, false, FileFormat.FASTQ, false, 4);
		
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
	
	/**
	 * Creates routing state without starting the inherited outer transfer thread.
	 * Destination streams are created lazily, with append-mode descriptors for reopening.
	 *
	 * @param pattern1_ Primary output file pattern with % placeholder for variable substitution
	 * @param pattern2_ Secondary output file pattern (may be null for single-end)
	 * @param overwrite_ Attempt deletion before a destination's first append-mode open
	 * @param append_ Stored base option; destination descriptors always use append mode
	 * @param allowSubprocess_ Whether subprocess execution is allowed
	 * @param useSharedHeader_ Request a shared header only on each destination's first open
	 * @param defaultFormat_ Default file format for outputs
	 * @param threaded_ Enable outer list transfer; callers must start that thread
	 * @param maxStreams_ Positive limit on simultaneously open logical destinations
	 */
	public MultiCros5(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, boolean threaded_, int maxStreams_){
		super(pattern1_, pattern2_, overwrite_, append_, allowSubprocess_, useSharedHeader_, defaultFormat_, threaded_, maxStreams_);
		
		bufferMap=new LinkedHashMap<String, Buffer>();
		streamQueue=new ArrayDeque<String>(maxStreams);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns cached success without waiting for routing or output completion.
	 * @return true if the cached error indication is clear */
	@Override
	public boolean finishedSuccessfully(){
		return !errorState;
	}
	
	/**
	 * Routes a borrowed root on the single routing owner, bypassing the outer transfer queue.
	 * Creates buffers as needed, then applies per-buffer and aggregate load checks.
	 * @param r Nonnull root; read/base counts include its mate
	 * @param name Replacement string for the first percent in each output pattern
	 */
	@Override
	public void add(Read r, String name){
		Buffer b=bufferMap.get(name);
		if(b==null){
			b=new Buffer(name);
			bufferMap.put(name, b);
			//Note: I could adjust bytesPerBuffer threshold here in response to the number of buffers.
		}
		b.add(r);
		//TODO: Move load-handling logic here?
		if(bytesInFlight>=memLimitUpper){handleLoad0();}
	}
	
	/**
	 * Handles memory overload condition by dumping all buffers.
	 * Called when total bytes in flight exceeds the upper memory limit.
	 * Under assertions, exits if that dump hands off no root entries.
	 */
	@Override
	void handleLoad0(){
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
	
	/**
	 * Performs single-use residual disposal after eligible data has been drained.
	 * Counts below-minimum reads even without a sink and clears every list reference.
	 * @param rosu Optional sink; every list uses ID zero, requiring support for repeated IDs
	 * @return Cumulative residual individual-read count, including mates
	 */
	@Override
	public long dumpResidual(ConcurrentReadOutputStream rosu){
		//For each Buffer, check if it contains residual reads
		//If so, dump it into the stream
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			Buffer b=e.getValue();
			assert((b.readsIn<minReadsToDump)==(b.list!=null && !b.list.isEmpty()));
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
	
	/**
	 * Reports dumped destinations in first-observation order without draining them.
	 * Rows contain name, individual reads, bases and optional cardinality. A positive
	 * minimum adds a three-column residual row; observe after routing for stable totals.
	 * @return ByteBuilder containing formatted statistics report
	 */
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
	
	/** @return Live key-set view of all observed names, including never-opened destinations */
	@Override
	public Set<String> getKeys(){return bufferMap.keySet();}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Dumps eligible roots after submissions stop, then closes and joins active outputs.
	 * For an unchanged minimum, below-minimum roots remain for terminal residual handling.
	 * @return Root-list entries handed off by dumpAll, not mate-inclusive reads
	 */
	@Override
	long closeInner(){
		//First dump everything
		final long x=dumpAll();
		//Then, retire any active streams
		while(!streamQueue.isEmpty()){retire(4);}
		return x;
	}
	
	/**
	 * Attempts every buffer, bypassing stream-opening deferral but retaining minreads.
	 * Used during memory pressure and final cleanup; may open outputs and retire others.
	 * @return Root-list entries handed off, excluding separate mate counts
	 */
	@Override
	long dumpAll(){
		if(verbose){
			System.err.println("before dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight);
		}
		long dumped=0;
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			dumped+=e.getValue().dump(true);
		}
		if(verbose){
			System.err.println("after dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight+", dumped="+dumped);
		}
		return dumped;
	}
	
	/**
	 * Sorts an active-queue prefix by logical creation/last-dump ordinal and retires its oldest.
	 * Candidates number min(queue size, 2*retCount, retCount+4); unselected names rejoin
	 * the queue's end. Dumps pending data, requests all selected closes, then joins them.
	 * @param retCount Requested retirements, clamped to the active queue size
	 */
	private void retire(int retCount){
		if(verbose){System.err.println("Enter retire("+retCount+"); streamQueue="+streamQueue);}
		final long time0=System.nanoTime(), time1, time2, time3, time4;
		retCount=Tools.min(streamQueue.size(), retCount);
		final int sortCount=Tools.min(streamQueue.size(), retCount*2, retCount+4);
		
		ArrayList<Buffer> rlist=new ArrayList<Buffer>(sortCount);
		//Select the first name in the queue, which is the least-recently-used.
		for(int i=0; i<sortCount; i++){
			String name=streamQueue.removeFirst();
			Buffer b=bufferMap.get(name);
			rlist.add(b);
			if(verbose){
				System.err.println("rlist.add("+name+","+b.timestamp+"); ros="+(b.currentRos!=null));
			}
		}
		Collections.sort(rlist);
		time1=System.nanoTime();
		for(int i=0; i<sortCount; i++){
			Buffer b=rlist.get(i);
			if(i<retCount){
				b.dump(b.currentRos);
				if(verbose){System.err.println("retire("+b.name+","+b.timestamp+")");}
			}else{
				streamQueue.add(b.name);
				rlist.set(i, null);
			}
		}
		
		time2=System.nanoTime();
		for(int i=0; i<retCount; i++){rlist.get(i).currentRos.close();}
		time3=System.nanoTime();
		for(int i=0; i<retCount; i++){
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
	
	/*--------------------------------------------------------------*/
	/*----------------          Profiling           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Formats retirement timings normalized by at least one recorded retirement.
	 * The retires-per-call ratio remains unguarded when no requests occurred.
	 * @return Diagnostic text, not a completion or benchmark result */
	public String printRetireTime(){
		ByteBuilder bb=new ByteBuilder();
		float mult=0.001f/Tools.max(1, retireCount);//#001 guard: was /retireCount = Infinity when retireCount==0 (no streams retired, e.g. fewer barcodes than maxStreams). Profiling-only. Family twin of MultiCros6#001.
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

	/** Accumulated candidate selection/sort nanoseconds. */
	private long retireTime1=0;
	/** Accumulated pending-dump/queue-return nanoseconds. */
	private long retireTime2=0;
	/** Accumulated close-request nanoseconds. */
	private long retireTime3=0;
	/** Accumulated join/status-check nanoseconds. */
	private long retireTime4=0;
	/** Total selected destinations across retirement calls. */
	private long retireCount=0;
	/** Number of retirement calls, including calls that select no destination. */
	private long retireCalls=0;
	
	/** Formats creation timings using retirement count as the existing denominator.
	 * The displayed total excludes phase one; phase two currently has zero duration.
	 * @return Diagnostic text, not an independently measured throughput result */
	public String printCreateTime(){
		ByteBuilder bb=new ByteBuilder();
		float mult=0.001f/Tools.max(1, retireCount);//#001 guard: was /retireCount = Infinity when retireCount==0 (no streams retired, e.g. fewer barcodes than maxStreams). Profiling-only. Family twin of MultiCros6#001.
		bb.append("Create Time 1:\t").append(createTime1*mult, 2).append(" us").nl();
		bb.append("Create Time 2:\t").append(createTime2*mult, 2).append(" us").nl();
		bb.append("Create Time 3:\t").append(createTime3*mult, 2).append(" us").nl();
		bb.append("Create Time 4:\t").append(createTime4*mult, 2).append(" us").nl();
		bb.append("Create Time 5:\t").append(createTime5*mult, 2).append(" us").nl();
		bb.append("Total:        \t").append((createTime2+createTime3+createTime4+createTime5)/1000000000.0, 3).append(" s").nl();
		return bb.toString();
	}
	
	/** Creation setup/retirement nanoseconds, excluding first-open deletion. */
	private long createTime1=0;
	/** Retained phase-two counter; time2 currently equals time1. */
	private long createTime2=0;
	/** Output construction nanoseconds. */
	private long createTime3=0;
	/** Output startup and optional diagnostic nanoseconds. */
	private long createTime4=0;
	/** Active-queue insertion nanoseconds. */
	private long createTime5=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Retains one destination's roots, counters, append descriptors and logical recency. */
	private class Buffer implements Comparable<Buffer>{
		
		/** Creates descriptors and storage without opening files.
		 * @param name_ Replacement string used in both output patterns */
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
		
		/** Retains a borrowed root and updates pair-aware input statistics before checking load.
		 * @param r Nonnull root whose payload remains stable through downstream completion */
		void add(Read r){
			//Add the read
			list.add(r);
			
			//Gather statistics
			long size=r.countPairBytes();
			int count=r.pairCount();
			currentBytes+=size;
			bytesInFlight+=size;
			basesIn+=r.pairLength();//was +=size (countPairBytes): bytes-as-bases made DemuxByName2 subtract oversized residualBases -> negative Bases Out (replicated 2026-09-05); MultiCros4 already used pairLength
			readsInFlight+=count;
			readsIn+=count;
			if(trackCardinality){loglog.hash(r);}
			
			//Decide whether to dump
			handleLoadB();
		}
		
		/** Attempts an ordinary dump at size thresholds, including smaller thresholds for open outputs. */
		private void handleLoadB(){
			//3rd term allows preemptive dumping
			//More generally, this triggers a dump if the reads in this buffer exceed
			//the maximum allowed reads or bytes 
			final int size=list.size();
			if(size>=readsPerBuffer || currentBytes>=bytesPerBuffer || 
					(currentRos!=null && (size>=200 || currentBytes>400000))){
				if(verbose){
					System.err.println("list.size="+list.size()+"/"+readsPerBuffer+
							", bytes="+currentBytes+"/"+bytesPerBuffer+", bytesInFlight="+bytesInFlight+"/"+memLimitUpper);
				}
				dump(false);
			}
		}
		
		/** Dumps eligible data, allowing ordinary calls to defer opening a new destination.
		 * @param force Bypass opening deferral, not the minimum-read requirement
		 * @return Root entries handed off, or zero for empty/ineligible/deferred buffers */
		long dump(boolean force){
			if(list.isEmpty() || readsIn<minReadsToDump){return 0;}
			ConcurrentReadOutputStream ros=getStream(force);
			return ros==null ? 0 : dump(ros);
		}
		
		/** Transfers the old root list downstream and replaces it, advancing the logical clock.
		 * Bypasses minreads and opening policy; Read payloads remain borrowed.
		 * @param ros Existing output stream
		 * @return Root-list entries handed off, excluding separate mates */
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
			//TODO: Probable bug - STR233: add uses pairCount(), while this subtracts list.size().
			//Keep the existing counter units pending their separate review.
			readsInFlight-=size0;
			readsWritten+=size0;
			currentBytes=0;
			numDumps++;
			timestamp=(bufferTimer++);
			return size0;
		}
		
		/** Returns existing output, promoting its queue position, or attempts a lazy open.
		 * @param force Bypass the new-stream deferral rule
		 * @return Active stream, or null when opening is deferred */
		private ConcurrentReadOutputStream getStream(boolean force){
			if(verbose){System.err.println("Enter getStream("+name+"); ros="+(currentRos!=null)+", +streamQueue="+streamQueue);}
			
			if(currentRos!=null){//The stream already exists
//				assert(streamQueue.contains(name));//slow
				if(streamQueue.peekLast()!=name){
					//Move to end to prevent early retirement
					boolean b=streamQueue.remove(name);
					assert(b) : "streamQueue did not contain "+name+", but the ros was open.";
					streamQueue.addLast(name);
				}//Historical TODO: promotion may be redundant, but it changes the candidate prefix.
			}else{//The stream does not exist, so create it
				return createStream(force);
			}
			
			assert(currentRos!=null) : "The stream for "+name+" was not created.";
			assert(streamQueue.peekLast()==name) : "The stream for "+name+" was not placed in the queue.";
			if(verbose){System.err.println("Exit  getStream("+name+"); ros="+(currentRos!=null)+", +streamQueue="+streamQueue);}
			return currentRos;
		}
		
		/** Opens append-mode output, retiring candidates first when the destination limit is reached.
		 * Ordinary calls defer if the limit is reached and buffered bytes remain below memLimitLower.
		 * @param force Bypass that deferral, retaining the stream-count limit
		 * @return Newly started stream, or null if deferred */
		private ConcurrentReadOutputStream createStream(boolean force){
			if(!force && streamQueue.size()>=maxStreams && bytesInFlight<memLimitLower){return null;}
			assert(currentRos==null) : "This should never be called if there is an existing stream.";
			if(!deleted){
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
		
		/** Attempts first-open removal before append-mode output.
		 * @param ff Descriptor to remove; null is ignored */
		private void delete(FileFormat ff){
			if(ff==null){return;}
			assert(overwrite || !ff.exists()) : "Trying to delete file "+ff.name()+", but overwrite=f.  Please add the flag overwrite=t.";
			//TODO: Probable bug - STR238: deleteIfPresent ignores File.delete's result,
			//so failed removal can leave old content for subsequent append output.
			ff.deleteIfPresent();
		}
		
		/** Appends cumulative input counts and optional cardinality.
		 * @param bb Destination builder
		 * @return The same builder */
		ByteBuilder appendTo(ByteBuilder bb){
			bb.append(name).tab().append(readsIn).tab().append(basesIn);
			if(trackCardinality){bb.tab().append(loglog.cardinality());}
			return bb.nl();
		}
		
		/** @return This destination's ordinary report row */
		@Override
		public String toString(){
			return appendTo(new ByteBuilder()).toString();
		}
		
		/** Orders distinct retirement candidates by their unique logical ordinal.
		 * @param b Another candidate with a different timestamp
		 * @return Negative for older, positive for newer */
		@Override
		public int compareTo(Buffer b){
			assert(timestamp!=b.timestamp);
			return timestamp<b.timestamp ? -1 : 1;
		}
		
		/** Stable destination name. */
		private final String name;
		/** Primary append-mode output descriptor. */
		private final FileFormat ff1;
		/** Optional mate append-mode output descriptor. */
		private final FileFormat ff2;
		
		/** Open output, or null after joined retirement/before first open. */
		private ConcurrentReadOutputStream currentRos;
		
		/** Current list of buffered reads awaiting output */
		private ArrayList<Read> list;
		
		/** Cumulative individual input reads, including mates. */
		private long readsIn=0;
		/** Number of bases that have entered this buffer */
		private long basesIn=0;
		/** Root entries handed downstream, not a durable-output or mate-inclusive count. */
		@SuppressWarnings("unused")
		private long readsWritten=0;//This does not count read2!
		/** Estimated bytes retained in the pending root list. */
		private long currentBytes=0;
		/** Number of dump operations executed for this buffer */
		private long numDumps=0;
		/** Whether first-open deletion policy has been processed, including overwrite=false. */
		private boolean deleted=false;
		/** Logical ordinal assigned on construction and each nonempty dump, not wall-clock time. */
		private long timestamp=-1;
		/** Optional tracker estimating unique k-mers across reads and mates. */
		private CardinalityTracker loglog;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Next logical ordinal for buffer creation or nonempty dumps, owned by the routing thread. */
	private long bufferTimer=0;
	
	/** Active names; reuse promotes to the end and retirement sorts only a prefix. */
	private final ArrayDeque<String> streamQueue;
	
	/** Map of buffer names to Buffer objects for managing multiple outputs */
	public final LinkedHashMap<String, Buffer> bufferMap;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Unused retained option; changing it does not alter this class's retirement behavior.
	 */
	private static final boolean closeFast=false;

}
