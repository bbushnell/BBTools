package stream;

import java.util.ArrayDeque;
import java.util.ArrayList;
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
 * Buffered named-read routing with a bounded number of open CROS destinations.
 * Each destination retains a buffer and optional cardinality tracker after retirement.
 * The stream limit counts destinations, not physical mate files or internal threads.
 * Eligible buffers flush by size or memory pressure; cumulative individual-read
 * counts below minReadsToDump remain buffered even during normal close.
 * Destination descriptors append on reopen; overwrite deletes files on the first open.
 * The inherited threaded flag controls outer list transfer, independently of CROS threads.
 * Direct add(Read, String) calls require one routing owner and must not bypass a
 * running outer transfer thread. Lists handed to CROS are replaced, not copied;
 * retain Read payloads until downstream output has completed.
 *
 * @author Brian Bushnell
 * @date May 14, 2019
 */
public class MultiCros3 extends BufferedMultiCros{
	
	/**
	 * Testing main method that demonstrates basic usage.
	 * Creates a MultiCros3 instance and processes reads from input file,
	 * distributing them to output files based on barcode.
	 * @param args Positional input and output pattern; trailing names are parsed but unused
	 */
	public static void main(String[] args){
		
		//Do some basic parsing
		String in=args[0];
		String pattern=args[1];
		ArrayList<String> names=new ArrayList<String>();
		for(int i=2; i<args.length; i++){names.add(args[i]);}
		
		//Create an mcros for this pattern
		MultiCros3 mcros=new MultiCros3(pattern, null, false, false, false, false, FileFormat.FASTQ, false, 4);
		
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
	 * Constructs a MultiCros3 with specified configuration parameters.
	 * Initializes buffer map and stream queue for managing multiple concurrent streams.
	 *
	 * @param pattern1_ Primary pattern; the first percent uses replaceFirst name substitution
	 * @param pattern2_ Secondary output file pattern (may be null)
	 * @param overwrite_ Whether to overwrite existing output files
	 * @param append_ Stored base option; this implementation opens destination descriptors in append mode
	 * @param allowSubprocess_ Whether to allow subprocess spawning for compression
	 * @param useSharedHeader_ Whether to use shared header across files
	 * @param defaultFormat_ Default file format for outputs
	 * @param threaded_ Whether the inherited list path uses an outer routing thread
	 * @param maxStreams_ Maximum number of open destination streams, not physical files
	 */
	public MultiCros3(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, boolean threaded_, int maxStreams_){
		super(pattern1_, pattern2_, overwrite_, append_, allowSubprocess_, useSharedHeader_, defaultFormat_, threaded_, maxStreams_);
		
		bufferMap=new LinkedHashMap<String, Buffer>();
		streamQueue=new ArrayDeque<String>(maxStreams);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the cached error flag without checking completion or waiting.
	 * @return true when errorState is clear */
	@Override
	public boolean finishedSuccessfully(){
		return !errorState;
	}
	
	/**
	 * Adds a read to the buffer associated with the given name.
	 * Creates a new buffer if one doesn't exist for this name.
	 * May trigger buffer dumps if thresholds are exceeded.
	 *
	 * This direct method runs on the calling routing owner; it does not enqueue work.
	 * @param r Nonnull read root; its mate contributes to individual-read/base counts
	 * @param name Destination name used as a replaceFirst replacement string
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
	}
	
	/**
	 * Disposes of below-minimum buffers after ordinary close has drained eligible ones.
	 * All buffer list references are cleared, so this is a terminal, single-use operation.
	 * Residual counts include mates and accumulate even when the output is null.
	 * @param rosu Optional residual output; each submitted list uses ID zero, so the
	 * sink must tolerate repeated IDs rather than require a dense increasing sequence
	 * @return Cumulative residual individual-read count
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
	 * Reports destinations that have dumped, in first-observation order.
	 * Rows contain name, cumulative individual reads and bases, plus optional cardinality.
	 * A three-column residual row is included when the configured minimum is positive.
	 * Read after routing/finalization to obtain stable totals.
	 * @return Newly allocated report; this does not itself drain output
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
	
	/** Exposes all observed names, including destinations that never opened a file.
	 * @return Live key-set view backed by bufferMap; callers must preserve routing state */
	@Override
	public Set<String> getKeys(){return bufferMap.keySet();}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Dumps eligible buffers, then closes and joins all open destination streams.
	 * Below-minimum buffers remain available for a later terminal dumpResidual call.
	 * @return Number of root-list entries handed off by dumpAll, excluding separate mate counts
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
	 * Tries every buffer without bypassing the cumulative minimum-read requirement.
	 * Used under memory pressure and during close; zero may leave ineligible data buffered.
	 * @return Number of root-list entries handed to destination streams
	 */
	@Override
	long dumpAll(){
		if(verbose){
			System.err.println("before dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight);
		}
		long dumped=0;
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			dumped+=e.getValue().dump();
		}
		if(verbose){
			System.err.println("after dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight+", dumped="+dumped);
		}
		return dumped;
	}
	
	/**
	 * Retained older retirement path; no current in-class caller uses it.
	 * Removes the first stream-queue entry to free resources.
	 * Dumps any remaining buffered reads before closing the stream.
	 * Uses either fast or traditional close based on closeFast setting.
	 */
	private void retire_old(){
		if(verbose){System.err.println("Enter retire(); streamQueue="+streamQueue);}
		
		//Select the first name in the queue, which is the least-recently-used.
		String name=streamQueue.removeFirst();
		Buffer b=bufferMap.get(name);
		
		if(verbose){System.err.println("retire("+name+"); ros="+(b.currentRos!=null)+", streamQueue="+streamQueue);}
		assert(b!=null);
		assert(b.currentRos!=null) : name+"\n"+streamQueue+"\n";
		
		//Purge any remaining reads first
		b.dump(b.currentRos);
		
		//Then close the stream
		if(closeFast){
			b.currentRos.close();//Faster, but the error state is not caught
			//Also the stream could be re-opened before writing is done, which is very bad, so this is unsafe as implemented.
		}else{
			//Unfortunately, this is very slow.
			errorState=ReadWrite.closeStream(b.currentRos) | errorState;//Traditional synchronous close-and-wait
		}
		
//		ros.close();
//		ros.join();
//		errorState|=(ros.errorState() || !ros.finishedSuccessfully());
		
		//Delete the pointer to output stream
		b.currentRos=null;
		if(verbose){System.err.println("Exit retire("+name+"); ros="+(b.currentRos!=null)+
				", streamQueue="+streamQueue);}
//		assert(!streamQueue.contains(name)); //Slow
	}
	
	/**
	 * Removes up to count oldest stream-queue entries, flushes them, then closes and joins.
	 * Recency is refreshed when a buffer requests its stream, not on every incoming read.
	 * Close requests are issued together before joins; errors are accumulated afterward.
	 * @param count Nonnegative maximum number of destination streams to retire
	 */
	private void retire(int count){
		if(verbose){System.err.println("Enter retire("+count+"); streamQueue="+streamQueue);}
		final long time0=System.nanoTime(), time1, time2, time3, time4;
		count=Tools.min(streamQueue.size(), count);
		
		ArrayList<Buffer> rlist=new ArrayList<Buffer>(count);
		//Select the first name in the queue, which is the least-recently-used.
		for(int i=0; i<count; i++){
			String name=streamQueue.removeFirst();
			Buffer b=bufferMap.get(name);
			rlist.add(b);
		}
		time1=System.nanoTime();
		for(Buffer b : rlist){
			if(verbose){
				System.err.println("retire("+b.name+"); ros="+
						(b.currentRos!=null)+", streamQueue="+streamQueue);
			}
			assert(b!=null);
			assert(b.currentRos!=null) : b.name+"\n"+streamQueue+"\n";
			//Purge any remaining reads first
			b.dump(b.currentRos);
		}
		time2=System.nanoTime();
		for(Buffer b : rlist){b.currentRos.close();}
		time3=System.nanoTime();
		for(Buffer b : rlist){
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
		retireCount+=count;
		retireCalls++;
	}
	
	/** Formats four retirement phase timings, normalized by at least one retired stream.
	 * The separate retires-per-call ratio remains unguarded when there were no calls.
	 * @return Diagnostic timings; reading this does not wait for output completion */
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

	/** Nanoseconds spent selecting, flushing, requesting close, and joining retired streams. */
	private long retireTime1=0;
	private long retireTime2=0;
	private long retireTime3=0;
	private long retireTime4=0;
	/** Total retired destinations and number of retirement invocations. */
	private long retireCount=0;
	private long retireCalls=0;
	
	/** Formats creation phase timings using the existing retired-stream normalization.
	 * The displayed total excludes phase one, which includes any retirement cost.
	 * @return Diagnostic timing text; not a benchmark or completion barrier */
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
	
	/** Accumulated creation phases; phase one includes capacity-driven retirement. */
	private long createTime1=0;
	private long createTime2=0;
	private long createTime3=0;
	private long createTime4=0;
	private long createTime5=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Per-name retained statistics, pending roots and optional active output stream. */
	private class Buffer implements Comparable<Buffer>{
		
		/** Creates append-mode descriptors and optional cardinality state, without opening output.
		 * @param name_ Replacement string used for the first percent in each pattern */
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
		
		/** Retains a root, updates pair-aware counters/cardinality and checks load thresholds.
		 * @param r Nonnull root whose payload remains owned by the caller */
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
			handleLoad();
		}
		
		/** Attempts local/global dumps after size thresholds; minimum-read gating still applies. */
		private void handleLoad(){
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
		
		/** Flushes a nonempty eligible buffer, opening or refreshing its destination stream.
		 * @return Root-list entries handed off, or zero for an empty/ineligible buffer */
		long dump(){
			if(list.isEmpty() || readsIn<minReadsToDump){return 0;}
			ConcurrentReadOutputStream ros=getStream();
			return dump(ros);
		}
		
		/** Hands off this list and replaces it; the supplied stream owns asynchronous output.
		 * This overload bypasses minimum-read gating and is used during retirement.
		 * @param ros Existing destination stream
		 * @return Root-list entries handed off, excluding separate mate counts */
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
			//TODO: Probable bug - STR233: add increases readsInFlight by pairCount(), but
			//this subtracts root-list size. Paired dumps leave an in-flight counter remainder.
			//Keep existing units until a separate counter review; memory gating uses bytes.
			readsInFlight-=size0;
			readsWritten+=size0;
			currentBytes=0;
			numDumps++;
			timestamp=(bufferTimer++);
			return size0;
		}
		
		/** Opens the destination if needed and moves its name to the stream queue's end.
		 * @return Active destination output, newest in the retirement queue */
		private ConcurrentReadOutputStream getStream(){
			if(verbose){System.err.println("Enter getStream("+name+"); ros="+(currentRos!=null)+", +streamQueue="+streamQueue);}
			
			if(currentRos!=null){//The stream already exists
//				assert(streamQueue.contains(name));//slow
				if(streamQueue.peekLast()!=name){
					//Move to end to prevent early retirement
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
		
		/** Deletes first-open outputs when overwrite is set, retires capacity, then opens append output.
		 * Shared headers are requested only before the first dump for this name.
		 * @return Newly started destination stream */
		private ConcurrentReadOutputStream createStream(){
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
		
		/** Attempts removal before the first append-mode open when overwrite is requested.
		 * @param ff Output descriptor; null is ignored */
		private void delete(FileFormat ff){
			if(ff==null){return;}
			assert(overwrite || !ff.exists()) : "Trying to delete file "+ff.name()+", but overwrite=f.  Please add the flag overwrite=t.";
			//TODO: Probable bug - STR238: FileFormat.deleteIfPresent ignores File.delete's
			//boolean result. Failed removal can leave old content for the append-mode open.
			ff.deleteIfPresent();
		}
		
		/** Appends this name's cumulative individual reads/bases and optional cardinality.
		 * @param bb Destination builder
		 * @return The supplied builder with one newline-terminated row added */
		ByteBuilder appendTo(ByteBuilder bb){
			bb.append(name).tab().append(readsIn).tab().append(basesIn);
			if(trackCardinality){bb.tab().append(loglog.cardinality());}
			return bb.nl();
		}
		
		/** @return A new String containing this buffer's report row */
		@Override
		public String toString(){
			return appendTo(new ByteBuilder()).toString();
		}
		
		/** Compares creation/last-dump timestamps, which are expected to be distinct.
		 * @param b Other buffer
		 * @return Negative when older, positive otherwise */
		@Override
		public int compareTo(Buffer b){
			assert(timestamp!=b.timestamp);
			return timestamp<b.timestamp ? -1 : 1;
		}
		
		/** Destination replacement name and permanently append-mode format descriptors. */
		private final String name;
		private final FileFormat ff1;
		private final FileFormat ff2;
		
		/** Currently open stream, or null after retirement. */
		private ConcurrentReadOutputStream currentRos;
		
		/** Pending roots; handed off and replaced on dump, cleared by terminal residual handling. */
		private ArrayList<Read> list;
		
		/** Cumulative individual reads and sequence bases, including mates. */
		private long readsIn=0;
		private long basesIn=0;
		@SuppressWarnings("unused")
		private long readsWritten=0;//This does not count read2!
		/** Estimated bytes in the currently pending roots. */
		private long currentBytes=0;
		/** Submitted batch count, also used as the next CROS batch ID across reopenings. */
		private long numDumps=0;
		/** Monotonic stamp assigned at creation and after each dump. */
		private long timestamp=-1;
		/** Optional tracker created with the buffer; configuration must remain stable. */
		private CardinalityTracker loglog;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Source of creation and dump timestamps for retained buffers. */
	private long bufferTimer=0;
	
	/** Active names in stream-request recency order; paired files share one entry. */
	private final ArrayDeque<String> streamQueue;
	
	/** All observed names and retained state, in first-observation order; preserve ownership. */
	public final LinkedHashMap<String, Buffer> bufferMap;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Disabled historical fast-close option, used only by the old retirement method. */
	private static final boolean closeFast=false;

}
