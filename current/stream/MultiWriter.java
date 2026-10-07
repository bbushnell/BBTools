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
 * Writer-based multi-destination output selected as mcrostype=7 by BufferedMultiCros.
 * Buffers reads by name, uses a size heap for opportunistic dumps, and bounds the number
 * of open per-name Writer instances. Retirement sorts a prefix of queue candidates by
 * logical timestamp; it is not a global least-recently-used search.
 * Reopened outputs use append mode, with an optional first-open deletion for overwrite.
 * Each writer is requested from WriterFactory with threads=1. Its actual implementation
 * and background-thread policy depend on format and available threads; paired outputs
 * may wrap two writers. The inherited threaded flag controls the outer transfer queue,
 * independently of those per-output writer choices.
 *
 * Historical implementation note: an 800k-read, 48-group demux with streams=8 was
 * reported to match MultiCros6 wall time and produce byte-identical unpaired, twin-file,
 * gzipped, nonthreaded and minreads outputs. This is retained prior evidence, not a
 * benchmark of the current implementation.
 *
 * @author Brian Bushnell
 * @contributor Furina
 * @date June 25, 2026
 */
public class MultiWriter extends BufferedMultiCros{

	/**
	 * Diagnostic driver that routes input reads by their barcode to the output pattern.
	 * Takes input and output-pattern arguments; additional names are parsed but unused.
	 * Uses direct single-read additions and closes both input and output after iteration.
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
		//TODO: Probable load-policy bug (STR383) - this direct-add driver omits the per-list
		//handleLoad0 check used by batch submission; lazy outputs can retain input until final close.
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

	/**
	 * Initializes inherited output settings, the name map, retirement queue and size heap.
	 * Does not start this outer thread or open any per-name writer.
	 * @param pattern1_ Primary pattern containing a percent placeholder; the parent expands # for twin files
	 * @param pattern2_ Optional second pattern containing a percent placeholder
	 * @param overwrite_ Whether to delete existing outputs before the first writer is created
	 * @param append_ Stored by the parent; per-name descriptors always use append mode
	 * @param allowSubprocess_ Whether output descriptors permit subprocesses
	 * @param useSharedHeader_ Whether to request shared headers on the first writer for each name
	 * @param defaultFormat_ Fallback output format
	 * @param threaded_ Whether the inherited list-submission API uses an outer transfer thread
	 * @param maxStreams_ Maximum number of open per-name Writer instances
	 */
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

	/** Returns the inverse of the accumulated error flag; this query does not wait for completion. */
	@Override
	public boolean finishedSuccessfully(){return !errorState;}

	/**
	 * Adds a nonnull read and its mate to the named logical buffer without copying them.
	 * Called by inherited list routing, or directly by the single owner in nonthreaded mode.
	 * This overload updates the buffer but does not perform the list API's global load check.
	 * @param r Submitted first read, with an optional mate
	 * @param name Destination name used in pattern replacement
	 */
	@Override
	public void add(Read r, String name){
		Buffer b=bufferMap.get(name);
		if(b==null){
			b=new Buffer(name);
			bufferMap.put(name, b);
		}
		b.add(r);
	}

	/** Dumps eligible buffers according to available writer slots and estimated-byte thresholds. */
	@Override
	void handleLoad0(){

		//Dump opportunistically because there are free streams
		final int streams=streamQueue.size();
		final int freeStreams=maxStreams-streams;
		if(freeStreams>0 && !heap.isEmpty()){
			final float mult=((streams+3f)/(maxStreams+2f));
			final long mll=(long)(memLimitLower*mult);
			final int bpb=(int)(bytesPerBuffer*mult);

			Buffer b=heap.peek();
			boolean dump=(b.currentBytes>=bpb && bytesInFlight>=mll);
			dump|=(b.currentBytes>=bpb*2);
			if(dump){
				heap.poll();
				b.dump(true);
			}
		}

		//Dump biggest buffers to free memory
		while(bytesInFlight>=memLimitMid && !heap.isEmpty()){
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

	/**
	 * Performs a terminal residual pass after normal output has been drained.
	 * Buffers below minReadsToDump contribute individual-read and base totals, and their
	 * lists are optionally handed to the caller's residual CROS with ID 0. Every buffer's
	 * list is then set to null; this is not a repeatable flush or a preparation for more adds.
	 * @param rosu Optional residual destination supplied by the caller
	 * @return Accumulated residual individual-read count
	 */
	@Override
	public long dumpResidual(ConcurrentReadOutputStream rosu){
		//For each Buffer, check if it contains residual reads
		//If so, dump it into the stream. rosu is the CALLER's unmatched-reads catch-all (a CROS),
		//not a per-name MultiWriter stream - forwarded as-is, exactly like MultiCros6.
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
	 * Returns tab-separated input read/base totals for destinations that have dumped at least once.
	 * Preserves name insertion order, optionally includes cardinality, and prepends residual
	 * totals when minReadsToDump is positive. Totals are not a durable-write confirmation.
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

	/** Returns the live key-set view of the insertion-ordered destination map. */
	@Override
	public Set<String> getKeys(){return bufferMap.keySet();}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Dumps eligible buffers, retires open writers and returns the number of list entries dumped. */
	@Override
	long closeInner(){
		//First dump everything
		final long x=dumpAll();
		assert(heap.isEmpty());
		//Then, retire any active streams
		if(!streamQueue.isEmpty()){retire(streamQueue.size());}
		assert(streamQueue.isEmpty());
		return x;
	}

	/** Empties the size heap and visits all buffers; destinations below minReadsToDump remain buffered. */
	@Override
	long dumpAll(){
		if(verbose){
			System.err.println("before dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight);
		}
		long dumped=0;
		while(!heap.isEmpty()){
			dumped+=heap.poll().dump(true);
		}
		//NOTE (not a bug - do not "simplify"): buffers that were in the heap are visited AGAIN here (they
		//remain in bufferMap), so dump(true) hits them twice - but the second call short-circuits on the
		//list.isEmpty() guard in dump() (no double-write, no double-count). This bufferMap pass exists to
		//also catch buffers that were never heaped (the small ones, below the addToHeap threshold).
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
	 * Retires up to retCount writers from a timestamp-sorted prefix of the open-name queue.
	 * Flushes selected buffers, returns unselected candidates to the queue, and folds each
	 * selected writer's poisonAndWait result into errorState before dropping its reference.
	 * @param retCount Requested number of open writers to retire
	 */
	private void retire(int retCount){
		if(verbose){System.err.println("Enter retire("+retCount+"); streamQueue="+streamQueue);}
		final long time0=System.nanoTime(), time1, time2, time3, time4;
		retCount=Tools.min(streamQueue.size(), retCount);
		final int sortCount=Tools.min(streamQueue.size(), retCount*2+1, retCount+4);

		ArrayList<Buffer> rlist=new ArrayList<Buffer>(sortCount);
		//Select a prefix of retirement candidates; timestamps order this subset, not the whole queue.
		for(int i=0; i<sortCount; i++){
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
		for(int i=0; i<sortCount; i++){
			Buffer b=rlist.get(i);
			if(i<retCount){
				b.dump(b.currentWriter);
				if(verbose){System.err.println("retire("+b.name+","+b.timestamp+")");}
			}else{
				streamQueue.add(b.name);
				rlist.set(i, null);
			}
		}
		time2=System.nanoTime();
		//Delegate normal finalization to each concrete writer and fold its reported error state.
		//This replaces the separate CROS close/join/fold phases; writer implementations differ.
		for(int i=0; i<retCount; i++){
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
	private void addToHeap(Buffer b){
		assert(b.loc()<0);
		assert(b.list.size()>0);
		if(!heap.hasRoom()){
			heap=heap.resizeNew(heap.CAPACITY*2+1);
		}
		heap.add(b);
		assert(b.loc()>=0);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Profiling           ----------------*/
	/*--------------------------------------------------------------*/

	/** Formats retirement profiling totals, with per-retired-writer microseconds and a guarded divisor. */
	@Override
	public String printRetireTime(){
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

	/** Nanoseconds spent selecting and sorting retirement candidates. */
	private long retireTime1=0;
	/** Nanoseconds spent submitting retiring buffers and requeuing other candidates. */
	private long retireTime2=0;
	/** Nanoseconds spent finalizing selected writers. */
	private long retireTime3=0;
	/** Retained fourth-phase timing field; the folded phase currently adds zero. */
	private long retireTime4=0;
	/** Total number of per-name writers retired. */
	private long retireCount=0;
	/** Number of retirement calls. */
	private long retireCalls=0;

	/** Formats creation profiling using retirement count as an approximate guarded normalizer. */
	@Override
	public String printCreateTime(){
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

	/** Nanoseconds spent freeing writer slots before creation. */
	private long createTime1=0;
	/** Retained second-phase creation timer; the folded phase currently adds zero. */
	private long createTime2=0;
	/** Nanoseconds spent constructing writers and resetting their batch IDs. */
	private long createTime3=0;
	/** Nanoseconds spent starting writers. */
	private long createTime4=0;
	/** Nanoseconds spent adding new writer names to the retirement queue. */
	private long createTime5=0;

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * A Buffer holds reads destined for a specific file.
	 * When sufficient reads are present, it opens a Writer and dumps them to it.
	 * If too many streams are open, it closes another stream first.
	 */
	private class Buffer implements SetLoc<Buffer>{

		/**
		 * Constructs buffer for specified output name.
		 * Replaces the first percent placeholder using String.replaceFirst and creates ordered,
		 * append-mode descriptors. First-open overwrite deletion is deferred until writer creation.
		 * Initializes a logical timestamp, a read list and optional cardinality tracking.
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

			if(currentWriter!=null){
				assert(heapLoc<0);
				if(size>=200 || currentBytes>400000){
					dump(false);
				}
				return;
			}
			if(heapLoc>=0){
				heap.jiggleDown(this);
			}else if(currentBytes>=10000 && readsIn>=minReadsToDump){
				//Add eventually, but don't pollute the heap with tiny buffers
				addToHeap(this);
			}
		}

		/**
		 * Submits a nonempty eligible buffer, creating a writer if permitted.
		 * @param force Whether writer creation may retire another output below the lower memory threshold
		 * @return Number of submitted list entries, or zero if ineligible or no writer is available
		 */
		long dump(boolean force){
			if(list.isEmpty() || readsIn<minReadsToDump){return 0;}
			Writer ros=getStream(force);
			return ros==null ? 0 : dump(ros);
		}

		/**
		 * Hands a nonempty list to the supplied writer with its next dense batch ID.
		 * Replaces the submitted list, clears estimated buffered bytes and advances dump/timestamp
		 * statistics. Whether submission formats inline or enqueues work belongs to the writer.
		 * @param ros Destination writer
		 * @return Number of list entries handed off, excluding separate mate counts
		 */
		long dump(final Writer ros){
			if(verbose){System.err.println("Dumping "+name);}
			if(list.isEmpty()){return 0;}
			final long size0=list.size();

			//Hand the list to the writer under its concrete synchronous/asynchronous submission policy.
			//Ids must be a dense ascending sequence STARTING AT 0 for each writer instance (JobQueue
			//firstID=0), so nextJobId resets in createStream; numDumps spans retirements and cannot be used.
			ros.add(list, nextJobId);
			nextJobId++;

			//Create a new list, since the old one is now owned by the writer/output.
			list=new ArrayList<Read>(400);

			//Manage statistics
			bytesInFlight-=currentBytes;
			//TODO: Probable bug (STR233) - add counts individual reads via pairCount, but this subtracts only list entries.
			readsInFlight-=size0;
			readsWritten+=size0;
			currentBytes=0;
			numDumps++;
			timestamp=(bufferTimer++);
			return size0;
		}

		/** Returns an existing writer and moves its name to the queue tail, or attempts writer creation. */
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

		/**
		 * Creates and starts a factory-selected writer, resetting its dense batch IDs to zero.
		 * May decline creation when not forced, writer slots are full and estimated bytes are low.
		 * Otherwise performs first-open overwrite deletion and retires writers as needed.
		 * Shared headers are requested only before the destination's first successful handoff.
		 * @param force Whether to permit retirement below the lower memory threshold
		 * @return The new writer, or null when creation is deferred
		 */
		private Writer createStream(boolean force){
			if(!force && streamQueue.size()>=maxStreams && bytesInFlight<memLimitLower){return null;}
			assert(currentWriter==null) : "This should never be called if there is an existing stream.";
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

			//Request threads=1; the factory chooses format-specific writers and may run FASTQ/FASTA
			//inline when Shared.threads()<4. Preserve this request and per-instance dense IDs.
			//Historical implementation note: threads=0 synchronous output was reported about 60%
			//slower on the original workload; this is not a current performance measurement.
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
			//TODO: Probable bug (STR238) - deleteIfPresent ignores File.delete's result;
			//failed removal can leave old content for append-mode output.
			ff.deleteIfPresent();
		}

		/**
		 * Format this buffer's summary as a line of text.
		 * @param bb ByteBuilder to append the text
		 * @return The modified ByteBuilder
		 */
		ByteBuilder appendTo(ByteBuilder bb){
			bb.append(name).tab().append(readsIn).tab().append(basesIn);
			if(trackCardinality){bb.tab().append(loglog.cardinality());}
			return bb.nl();
		}

		@Override
		public String toString(){return appendTo(new ByteBuilder()).toString();}

		/** Orders larger estimated byte counts first in the min-heap. */
		@Override
		public int compareTo(Buffer b){
			long dif=b.currentBytes-currentBytes;
			return dif<0 ? -1 : dif>0 ? 1 : 0;
		}

		/** Stores the heap-maintained array position, or -1 when absent. */
		@Override
		public void setLoc(int newLoc){heapLoc=newLoc;}

		/** Returns the heap-maintained array position. */
		@Override
		public int loc(){return heapLoc;}

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

		/** Individual reads received, including mates; used for minReadsToDump and reporting. */
		private long readsIn=0;
		/** Bases received, including mate bases. */
		private long basesIn=0;
		/** List entries handed to writers, excluding mate counts; not a durable-write count. */
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
		/** Logical timestamp assigned on buffer creation and each nonempty dump. */
		private long timestamp=-1;
		/** Location in heap; -1 means not in heap */
		private int heapLoc=-1;
		/** Optional, for tracking cardinality */
		private CardinalityTracker loglog;

	}

	/**
	 * Comparator for sorting buffers by timestamp for retirement ordering.
	 * Orders distinct logical timestamps within the selected retirement candidate subset.
	 */
	private static final class TimestampComparator implements Comparator<Buffer>{

		/** Creates the shared timestamp comparator. */
		private TimestampComparator(){}

		/** Compares distinct logical timestamps; equal timestamps violate the caller invariant. */
		@Override
		public final int compare(Buffer a, Buffer b){
			assert(a.timestamp!=b.timestamp);
			return a.timestamp<b.timestamp ? -1 : 1;
		}

		/** Singleton instance of timestamp comparator */
		static final TimestampComparator instance=new TimestampComparator();

	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	/** Monotonic logical clock advanced for each buffer creation and nonempty dump. */
	private long bufferTimer=0;

	/** Priority heap containing buffers ordered by size for dump prioritization */
	private HeapLoc<Buffer> heap;

	/** Open destination names in retirement-candidate order. */
	private final ArrayDeque<String> streamQueue;

	/** Insertion-ordered map retaining buffers after their writers are retired. */
	public final LinkedHashMap<String, Buffer> bufferMap;

}
