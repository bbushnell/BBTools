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
 * Buffered named-read routing with token-limited outputs and background retirement.
 * Construction starts up to eight retirement helpers; the inherited threaded flag
 * separately controls the outer list-routing thread and does not start that thread.
 * Tokens count logical destinations, including retiring streams, not physical mate files.
 * Direct routing requires one owner; Read payloads remain borrowed by downstream output.
 * Below-minimum buffers remain pending, and a retiring buffer skips ordinary dump attempts.
 * After submissions stop, close waits for outstanding retirements and drains eligible
 * pending roots before retiring the final outputs and joining helpers.
 * 
 * @author Brian Bushnell
 * @date April 5, 2024
 *
 */
public class MultiCros4 extends BufferedMultiCros{
	
	/** 
	 * For testing.<br>
	 * args should be:
	 * {input file, output pattern, names...}
	 * @param args Positional input and output pattern; trailing names are parsed but unused
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
	
	/** Creates routing state and immediately starts the retirement helpers.
	 * Actual helpers are capped at eight; the configured allowance is
	 * {@code max(1, (int)(maxStreams*0.25f))}, which also sizes the retirement queue
	 * and active-stream reserve.
	 * @param pattern1_ Primary pattern containing a percent placeholder
	 * @param pattern2_ Optional mate pattern; the inherited constructor may expand a hash
	 * @param overwrite_ Attempt deletion before each destination's first append-mode open
	 * @param append_ Stored base option; destination descriptors always use append mode
	 * @param allowSubprocess_ Output descriptor subprocess permission
	 * @param useSharedHeader_ Request shared headers on a destination's first open
	 * @param defaultFormat_ Fallback output format
	 * @param threaded_ Enable inherited outer list transfer; callers must start that thread
	 * @param maxStreams_ Positive logical-destination token budget, including retiring outputs */
	public MultiCros4(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, 
			boolean useSharedHeader_, int defaultFormat_, boolean threaded_, int maxStreams_){
		super(pattern1_, pattern2_, overwrite_, append_, allowSubprocess_, 
				useSharedHeader_, defaultFormat_, threaded_, maxStreams_);
		
		bufferMap=new LinkedHashMap<String, Buffer>();
		streamQueue=new ArrayDeque<String>(maxStreams);
		freeTokens=new ArrayBlockingQueue<Token>(maxStreams);
		for(int i=0; i<maxStreams; i++){freeTokens.add(new Token(i));}
		maxRetireThreads=Tools.max(1, (int)(maxStreams*0.25f));
		retireQueue=new ArrayBlockingQueue<Buffer>(maxRetireThreads);
		retireThreads=new ArrayList<RetireThread>(maxRetireThreads);
		maxOpenStreams=Tools.max(1, maxStreams-maxRetireThreads);
		for(int i=0; i<8 && i<maxRetireThreads; i++){
			RetireThread rt=new RetireThread();
			rt.start();
			retireThreads.add(rt);
		}
		if(verbose){System.err.println("maxStreams="+maxStreams+", maxRetireThreads="+
				maxRetireThreads+", maxOpenStreams="+maxOpenStreams);}//#002 gate: was unconditional stderr on every mcrostype=4 construction; all other diagnostics here are verbose-gated.
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the cached error indication without waiting or checking completion.
	 * @return true when errorState is clear */
	@Override
	public boolean finishedSuccessfully(){
		return !errorState;
	}
	
	/** Routes a root directly on the routing owner; does not use the outer transfer queue.
	 * @param r Nonnull root whose mate contributes to read/base counts
	 * @param name Name used as a replaceFirst replacement string in output patterns */
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
	
	/** Terminal, single-use disposal of below-minimum buffers after eligible data is drained.
	 * Clears all list references and accumulates pair-aware residual counts even without a sink.
	 * @param rosu Optional sink; each list uses ID zero, so repeated IDs must be supported
	 * @return Cumulative residual individual-read count */
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
	
	/** Reports dumped destinations in first-observation order without draining them.
	 * Rows contain name, individual reads, bases and optional cardinality; a positive
	 * minimum adds a three-column residual row. Observe after routing for stable totals.
	 * @return Newly allocated report builder */
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

	/** Drains eligible roots after submissions stop, then retires outputs and joins helpers.
	 * Waits for earlier retirements before the final dump; below-minimum roots remain pending.
	 * @return Root-list entries handed off by both dump passes, not mate-inclusive reads */
	@Override
	long closeInner(){
		long x=dumpAll();
		while(!streamQueue.isEmpty()){retire(1);}
		//STR244: a retiring destination may retain roots added after retirement began.
		//Start the final pass with every destination CLOSED and no earlier retirement in flight.
		awaitRetirements();
		x+=dumpAll();
		//With submissions stopped, this pass can retire only destinations already drained.
		while(!streamQueue.isEmpty()){retire(1);}
		try{
			for(boolean success=false; !success;){
				retireQueue.put(POISON_BUFFER);
				success=true;
			}
		}catch(InterruptedException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		waitForFinishInner();
		return x;
	}

	/** Waits for every outstanding output retirement without terminating the helpers.
	 * Called by the routing owner after all active names have been queued for retirement.
	 * Holding every token establishes completion; tokens are returned before reopening output. */
	private void awaitRetirements(){
		assert(streamQueue.isEmpty()) : "Active outputs must be retired before waiting for their tokens";
		final Token[] tokens=new Token[maxStreams];
		for(int i=0; i<tokens.length; i++){
			while(tokens[i]==null){
				try{
					tokens[i]=freeTokens.take();
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}
		for(Token token : tokens){freeTokens.add(token);}
	}
	
	/** Joins retirement helpers, not the inherited outer routing thread.
	 * Does not signal termination; helpers must first receive their terminal marker.
	 * Interrupted joins are logged and retried. */
	public final void waitForFinishInner(){
		if(verbose){System.err.println("Waiting for finish.");}
		for(RetireThread rt : retireThreads){
			synchronized(rt){
				while(rt.getState()!=Thread.State.TERMINATED){
					if(verbose){System.err.println("Attempting join: state="+rt.getState());}
					try{
						rt.join(1000);
					}catch(InterruptedException e){
						e.printStackTrace();
					}
				}
			}
		}
	}
	
	/** Attempts every buffer's ordinary dump, retaining minimum-read and RETIRING checks.
	 * @return Root-list entries handed off, excluding separate mate counts */
	@Override
	long dumpAll(){
		if(verbose){
			System.err.println("before dumpAll: bytesInFlight="+bytesInFlight+
					", limit="+memLimitUpper+", readsInFlight="+readsInFlight);
		}else{
//			System.err.println("dumpAll triggered due to memory pressure.");//Also happens at the end
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
	
	/** Removes the oldest active queue entry and requests asynchronous retirement.
	 * Existing OPEN streams return early from getStream, so access does not refresh this order.
	 * @param retCount Required to be one by the current implementation */
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
		
//		if(retireThreads.size()<maxRetireThreads && retireQueue.size()>0){
//			RetireThread rt=new RetireThread();
//			rt.start();
//			retireThreads.add(rt);
//		}
		
		try{
			for(boolean success=false; !success;){
				retireQueue.put(b);
				success=true;
			}
		}catch(InterruptedException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		time2=System.nanoTime();
		if(verbose){System.err.println("Added "+name+" to retire queue.");}
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
	
	/**
	 * Formats retirement-request timings, normalized by at least one recorded retirement.
	 * The retires-per-call ratio remains unguarded when no requests occurred.
	 * @return Diagnostic text; not a completion barrier or measured benchmark conclusion
	 */
	public String printRetireTime(){
		ByteBuilder bb=new ByteBuilder();
		float mult=0.001f/Tools.max(1, retireCount);//#001 guard: was /retireCount = Infinity when retireCount==0 (no streams retired, e.g. fewer barcodes than maxStreams). Profiling-only. Family twin of MultiCros6#001.
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

	/** Timing for retire phase 1: preparation and queue selection */
	private long retireTime1=0;
	/** Timing for retire phase 2: queue operations */
	private long retireTime2=0;
	/** Timing for retire phase 3: unused */
	private long retireTime3=0;
	/** Timing for retire phase 4: unused */
	private long retireTime4=0;
	/** Destination count recorded by retirement requests, not an independent completion count. */
	private long retireCount=0;
	/** Total number of retire() method calls */
	private long retireCalls=0;
	
	/**
	 * Formats creation timings using the existing retired-stream normalization.
	 * The displayed total excludes phase one, which includes retirement requests.
	 * @return Diagnostic timing text; reading does not wait for output completion
	 */
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
	
	/** Timing for creation phase 1: pre-stream setup */
	private long createTime1=0;
	/** Timing for creation phase 2: token acquisition */
	private long createTime2=0;
	/** Timing for creation phase 3: stream object creation */
	private long createTime3=0;
	/** Timing for creation phase 4: stream startup */
	private long createTime4=0;
	/** Timing for creation phase 5: queue management */
	private long createTime5=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Retains one destination's pending roots, cumulative statistics, stream and token.
	 * Minimum-read eligibility and retirement state determine whether ordinary dumps proceed.
	 */
	private class Buffer{
		
		/**
		 * Creates buffer for specific output file pattern.
		 * Initializes file formats and read list for buffering.
		 * No output is opened here; descriptors use append mode for later reopening.
		 * @param name_ Replacement string used for the first percent in each pattern
		 */
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
		 * @param r Nonnull borrowed root; counts and cardinality include its mate
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
			if(size>=readsPerBuffer || currentBytes>=bytesPerBuffer){
				if(verbose){
					System.err.println("list.size="+list.size()+"/"+readsPerBuffer+
							", bytes="+currentBytes+"/"+bytesPerBuffer+", bytesInFlight="+bytesInFlight+"/"+memLimitUpper);
				}
				dump();
			}else if((size>=200 || currentBytes>=400000) && state==OPEN){
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
		
		/** Dumps eligible data unless the buffer is RETIRING, opening output if needed.
		 * @return Root-list entries handed off, or zero for empty/ineligible/retiring state */
		long dump(){
			if(list.isEmpty() || readsIn<minReadsToDump){return 0;}
			final int state=getState();
			if(state==RETIRING){return 0;}
			ConcurrentReadOutputStream ros=getStream();
			return dump(ros);
		}
		
		/** 
		 * Dump buffered reads to the stream.
		 * If the buffer is empty, nothing happens. Otherwise ownership of the old list
		 * passes downstream and a fresh list replaces it; Read payloads remain borrowed.
		 * This overload bypasses minimum-read and retirement checks.
		 * @param ros Existing destination stream
		 * @return Root-list entries handed off, not individual-read count */
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
			//TODO: Probable bug - STR233: add increases readsInFlight by pairCount(),
			//but this subtracts root-list size. Preserve existing units pending separate review.
			readsInFlight-=size0;
			readsWritten+=size0;
			currentBytes=0;
			numDumps++;
			return size0;
		}
		
		/** Returns an OPEN stream immediately or creates output from CLOSED state.
		 * The early OPEN return bypasses the later queue-promotion code.
		 * @return Active destination output */
		private synchronized ConcurrentReadOutputStream getStream(){
			if(state==OPEN){return currentRos;}else if(state==CLOSED){assert(currentRos==null);}else{assert(false);}
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
		
		/** Acquires a token, starts append-mode output and appends this active name.
		 * May request another retirement or wait for a token. Overwrite deletion and
		 * shared-header forwarding apply only before this destination's first dump.
		 * @return Newly started stream */
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
			if(verbose){System.err.println("Fetching token for "+name);}
			while(t==null){
				try{
					t=freeTokens.take();
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			time2=System.nanoTime();
			if(verbose){System.err.println("Got token for "+name);}
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
		
		/** Attempts first-open removal before append-mode output.
		 * @param ff Descriptor to remove; null is ignored */
		private void delete(FileFormat ff){
			if(ff==null){return;}
			assert(overwrite || !ff.exists()) : "Trying to delete file "+ff.name()+", but overwrite=f.  Please add the flag overwrite=t.";
			//TODO: Probable bug - STR238: FileFormat.deleteIfPresent ignores File.delete's
			//boolean result, so an unsuccessful removal can leave old data for append output.
			ff.deleteIfPresent();
		}
		
		/** 
		 * Appends name, cumulative individual reads/bases and optional cardinality.
		 * @param bb ByteBuilder to append the text
		 * @return The modified ByteBuilder
		 */
		ByteBuilder appendTo(ByteBuilder bb){
			bb.append(name).tab().append(readsIn).tab().append(basesIn);
			if(trackCardinality){bb.tab().append(loglog.cardinality());}
			return bb.nl();
		}
		
		/** @return The ordinary report row, or state diagnostics when verbose is enabled */
		@Override
		public String toString(){
			if(verbose){return toString2();}
			return appendTo(new ByteBuilder()).toString();
		}
		
		/** @return Name, state value and whether a stream reference is present */
		public String toString2(){
			return name+" "+state+" "+(currentRos!=null);
		}
		
		/**
		 * Changes buffer state with validation.
		 * Requires the CLOSED-to-OPEN-to-RETIRING-to-CLOSED cycle under assertions.
		 * Entering CLOSED clears currentRos but does not detach or return the token.
		 * @param newState The new state to set
		 */
		synchronized void setState(int newState){
			if(verbose){
				System.err.println("setState "+name+" "+state+" -> "+newState+"; ros="+(currentRos!=null));
			}
			assert(currentRos!=null);
			assert(token!=null);
			
			int x=(state+1)%3;
			assert(state!=newState);
			assert(x==newState);
			state=newState;
			if(state==CLOSED){currentRos=null;}else if(state==RETIRING){assert(list.isEmpty());}
		}
		
		/** @return Synchronized state snapshot; this does not transfer token ownership */
		synchronized int getState(){
			return state;
		}
		
		/** Assigns token to this buffer for stream creation.
		 * @param t The token to assign */
		synchronized void giveToken(Token t){
			assert(token==null);
			assert(state==CLOSED);
			assert(t!=null);
			token=t;
		}
		
		/** Removes and returns this buffer's token.
		 * @return The token that was assigned to this buffer */
		synchronized Token takeToken(){
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
		
		
		/** Current buffer state: CLOSED, OPEN, or RETIRING */
		private int state=CLOSED;
		/** Token required for creating streams, limits concurrent stream creation */
		private Token token;
		
		/** Current list of buffered reads */
		private ArrayList<Read> list;
		/** Unused retained field; current dumping operates directly on list. */
		private ArrayList<Read> dumpList;
		
		/** Cumulative individual reads, including mates. */
		private long readsIn=0;
		/** Cumulative sequence bases, including mates. */
		private long basesIn=0;
		/** Root-list entries handed downstream, not a durable-output or mate-inclusive count. */
		@SuppressWarnings("unused")
		private long readsWritten=0;//This does not count read2!
		/** Number of bytes currently in this buffer (estimated) */
		private long currentBytes=0;
		/** Number of dumps executed */
		private long numDumps=0;
		/** Optional, for tracking cardinality */
		private CardinalityTracker loglog;
		
	}
	
	/** Constructor-started helper that closes queued outputs and returns their tokens. */
	private class RetireThread extends Thread{
		
		/** Consumes retirement requests; forwards the shared terminal marker before exiting. */
		@Override
		public void run(){
			while(true){
				Buffer b=null;
				if(verbose){System.err.println("Retire thread fetching buffer; retQueue="+retireQueue);}
				while(b==null){
					try{
						b=retireQueue.take();
					}catch(InterruptedException e){
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				if(verbose){System.err.println("Retire thread fetched "+b.name+"; retQueue="+retireQueue);}
				if(b==POISON_BUFFER){
					for(boolean success=false; success==false;){
						try{
							retireQueue.put(b);
							success=true;
						}catch(InterruptedException e){
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					if(verbose){System.err.println("Retire thread terminated.");}
					return;
				}
				retire(b);
			}
		}
		
		/** Closes a stream, records its result, publishes CLOSED and returns its token.
		 * @param b The buffer to retire */
		void retire(Buffer b){
			if(verbose){System.err.println("Retire thread retiring "+b.name+".");}
			
			assert(b.currentRos!=null) : b.state+", "+b.name;
			//TODO: Possible bug - STR246: multiple helpers update shared errorState with
			//an unsynchronized read-modify-write; a false result can overwrite a concurrent true.
			errorState=ReadWrite.closeStream(b.currentRos) | errorState;//Traditional synchronous close-and-wait
			
//			ros.close();
//			ros.join();
//			errorState|=(ros.errorState() || !ros.finishedSuccessfully());
			if(verbose){System.err.println("Exit retire("+b.name+"); ros="+(b.currentRos!=null)+
					", streamQueue="+streamQueue);}
			//STR245: reopening must observe CLOSED and an absent old token together.
			final Token t;
			synchronized(b){
				b.setState(CLOSED);
				t=b.takeToken();
			}
			freeTokens.add(t);

			if(verbose){System.err.println("Retire thread retired "+b.name+".");}
		}
		
	}
	
	/** One logical-destination reservation held through output retirement. */
	private static class Token{
		/** Creates token with specified ID.
		 * @param id_ Unique identifier for this token */
		Token(int id_){id=id_;}
		/** Unique identifier for this token */
		final int id;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Configured retirement allowance/queue capacity; actual helper count is capped at eight. */
	private final int maxRetireThreads;

	/** Active-name queue limit, reserving part of the token budget for retiring streams. */
	public final int maxOpenStreams;
	
	/** Active names in creation/reopen order; OPEN access does not refresh recency. */
	private final ArrayDeque<String> streamQueue;
	
	/** Tokens available for use */
	private final ArrayBlockingQueue<Token> freeTokens;
	
	/** Buffers waiting to retire */
	private final ArrayBlockingQueue<Buffer> retireQueue;
	
	/** Retirement helper instances started by the constructor. */
	private final ArrayList<RetireThread> retireThreads;
	
	/** Map of names to buffers */
	public final LinkedHashMap<String, Buffer> bufferMap;
	
	/** Special buffer used to signal retire threads to terminate */
	private final Buffer POISON_BUFFER=new Buffer("POISON_BUFFER_NOT_A_FILE");
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Cyclic stream states: closed, active, then queued/ongoing retirement. */
	private static final int CLOSED=0, OPEN=1, RETIRING=2;

}
