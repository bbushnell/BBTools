package stream;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.Parse;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;
import template.ThreadWaiter;

/** Parallel FASTQ reader with one byte-input thread and ordered conversion workers.
 * Call start once per instance, then consume batches through nextList. Interleaved
 * input links adjacent records as mates; maxReads and sampling positions count pairs
 * in that mode. Output order follows input batch order, including empty sampled batches.
 * Configure shared parser options and sampling before starting. Record decoding uses
 * nucleotide flags; this reader does not capture Shared.AMINO_IN.
 * EOF and close do not join the input coordinator or wait for its counter aggregation.
 * Error reporting retains the documented deferred shared-state concern below; this
 * class does not promise general thread safety for lifecycle/configuration changes.
 * @author Isla, Shinobu
 * @date October 30, 2025
 */
public class FastqStreamer implements Streamer{

	/** Legacy throughput driver that delegates reader selection to StreamerFactory.
	 * Positional arguments select input, default workers, the legacy SIMD toggle and
	 * vector validation. A present SIMD argument requests capability-gated SIMD.
	 * @param args Input filename followed by optional driver settings
	 */
	public static void main(final String[] args){
		Timer t=new Timer();
		String fname=args[0];
		if(args.length>1){DEFAULT_THREADS=Integer.parseInt(args[1]);}
		//STR-006: Preserve the argument-presence trigger while honoring JVM/hardware capability.
		if(args.length>2){Shared.SIMD=simd.Vector.simd256;}
		if(args.length>3){Read.VALIDATE_VECTOR=Parse.parseBoolean(args[3]);}

		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, true);
		Streamer st=StreamerFactory.makeStreamer(ff, 0, true, -1, true, true);
		st.start();
		long reads=0, bases=0;
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
			for(Read r : ln){
				reads+=r.pairCount();
				bases+=r.pairLength();
			}
		}
		t.stop();
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves a FASTQ descriptor with subprocess input allowed, then delegates.
	 * @param fname_ Input path or standard-input name
	 * @param threads_ Conversion workers; values below one select DEFAULT_THREADS
	 * @param pairnum_ Pair marker for noninterleaved reads: 0 or 1
	 * @param maxReads_ Input-record limit, or pair limit when interleaved; negative means unlimited
	 */
	public FastqStreamer(final String fname_, final int threads_, final int pairnum_, final long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTQ, null, true, false), threads_, pairnum_, maxReads_);
	}

	/** Captures the descriptor/configuration and creates queues without opening input.
	 * Adds one input coordinator beyond the conversion-worker count. Assertions require
	 * a valid pair marker and marker zero for interleaved input.
	 * @param ffin_ Nonnull input descriptor, including pairing and subprocess settings
	 * @param threads_ Worker request, defaulted below one and clamped to 1..Shared.threads()
	 * @param pairnum_ Pair marker for noninterleaved reads: 0 or 1
	 * @param maxReads_ Input-record or interleaved-pair limit before sampling; negative means unlimited
	 */
	public FastqStreamer(final FileFormat ffin_, final int threads_, final int pairnum_, final long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(pairnum==0 || !interleaved);
//		if(interleaved && maxReads_<Long.MAX_VALUE/2) {maxReads_*=2;}
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		// Create OQS with prototypes for LAST/POISON generation
		ListNum<byte[][]> inputPrototype=new ListNum<byte[][]>(null, 0, ListNum.PROTO);
		ListNum<Read> outputPrototype=new ListNum<Read>(null, 0, ListNum.PROTO);
		oqs=new OrderedQueueSystem<ListNum<byte[][]>, ListNum<Read>>(
			threads, true, inputPrototype, outputPrototype);

		if(verbose){outstream.println("Made FastqStreamer-"+threads);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts one input coordinator and the configured conversion workers.
	 * Resets aggregate counters, but does not reset queues or error state. Call once;
	 * repeated starts are not a supported restart operation. Returns without joining.
	 */
	@Override
	public void start(){
		if(verbose){outstream.println("FastqStreamer.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Process the reads in separate threads
		spawnThreads();

		if(verbose){outstream.println("FastqStreamer started.");}
	}

	/** Closes an available byte backend, then forces output-queue termination.
	 * Does not join this streamer's coordinator or conversion workers, or wait for
	 * their counter aggregation. Backend close may itself wait for backend threads.
	 * Does not fold the backend close result into errorState here. Use as early/emergency
	 * shutdown, not restart.
	 */
	@Override
	public void close(){
		if(bf!=null){bf.close(); bf=null;}
		//Emergency-abort completeness: also force-finish the OQS so a close() from a dying/early-exiting
		//consumer frees blocked workers (outq capacity-wait) and lets the input thread drain to its own
		//poison. Without this, a consumer death left every non-daemon pipeline thread blocked and the JVM
		//alive forever (replicated + jstack-proven via the paired-sampling crash, 2026-09-05). Idempotent
		//and harmless on the normal path (nextList already called setFinished on LAST).
		oqs.setFinished(true);
	}

	/** Returns the captured descriptor name.
	 * @return Input filename or standard-input name
	 */
	@Override
	public String fname(){return fname;}

	/** Returns the delegated queue hint without consuming a batch.
	 * Use nextList to detect EOF; this hint is not a termination guarantee after forced close.
	 * @return Current OrderedQueueSystem availability hint
	 */
	@Override
	public boolean hasMore(){
		return oqs.hasMore();
	}

	/** Returns the shared error flag without waiting for producer/worker completion.
	 * Retains the separately recorded shared-state limitation; not a final-status barrier.
	 * @return Currently observed error flag
	 */
	@Override
	public boolean errorState(){return errorState;}

	/** Returns the captured interleaved-input mode.
	 * @return true when adjacent input records are processed as mates
	 */
	@Override
	public boolean paired(){return interleaved;}

	/** Returns the captured marker used for noninterleaved reads.
	 * Interleaved reads receive their individual 0/1 markers during conversion.
	 * @return Configured pair marker
	 */
	@Override
	public int pairnum(){return pairnum;}

	/** Returns the aggregate count of converted individual reads, including both mates.
	 * Sampled-out records are excluded. The input coordinator aggregates after worker
	 * joins; EOF and close do not wait for that step or guarantee its visibility.
	 * @return Currently observed aggregate converted-read count
	 */
	@Override
	public long readsProcessed(){return readsProcessed;}

	/** Returns the aggregate base count for converted reads, including both mates.
	 * Aggregation and visibility have the same limits as readsProcessed.
	 * @return Currently observed aggregate converted-base count
	 */
	@Override
	public long basesProcessed(){return basesProcessed;}

	/** Configures positional sampling; call before start, not concurrently with workers.
	 * Workers keep every record for rates at least one; lower rates use sampleKeep.
	 * Input positions supply numeric IDs even when other records are skipped. Interleaved
	 * mates are kept or dropped together. This setter does not validate the rate.
	 * @param rate Requested fraction, normally in [0,1]
	 * @param seed Sampling seed; negative requests a newly chosen random seed
	 */
	@Override
	public void setSampleRate(final float rate, final long seed){
		samplerate=rate;
		sampleSeed=Streamer.resolveSampleSeed(seed);
	}

	/** Waits for the next ordered batch, which may be empty after sampling.
	 * Consuming LAST forces queue termination. A null/LAST result with an observed
	 * error flag invokes the existing fatal helper rather than returning clean EOF.
	 * This method does not wait for the coordinator to aggregate worker counters.
	 * @return Next batch, or null at terminal input when no error flag is observed
	 */
	@Override
	public ListNum<Read> nextList(){
		ListNum<Read> list=oqs.getOutput();
		if(verbose){
			if(list==null){outstream.println("Consumer got null.");}
			else {outstream.println("Consumer got list "+list.id()+" type "+list.type);}
		}
		if(list==null || list.last()){
			if(list!=null && list.last()){
				oqs.setFinished(true);
			}
			//[stream/FastqStreamer#001 FIXED] crash LOUD on a producer/worker error instead of silently truncating: a thread death
			//set errorState + force-poisoned the OQS (so we got null/last here). A bare `return null` would look like clean EOF ->
			//wrong/partial results. KillSwitch.kill exits loudly (BBTools contract: crash, never silently wrong).
			if(errorState){KillSwitch.kill("Error reading FASTQ file (corrupt or truncated): "+fname);}
			return null;
		}
		return list;
	}

	/** Rejects the SAM-line API for this FASTQ reader.
	 * @return Never returns normally
	 * @throws UnsupportedOperationException Always, because FASTQ has no SamLine output
	 */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTQ does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates and starts a shared thread list: coordinator at index zero, then workers. */
	void spawnThreads(){
		//Determine how many threads may be used
		final int threads=this.threads+1;

		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, alpt));
		}
		if(verbose){outstream.println("Spawned threads.");}

		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Input coordinator for tid zero, otherwise a conversion worker sharing the queues. */
	private class ProcessThread extends Thread{

		/** Assigns a diagnostic name and retains the join list only for the coordinator.
		 * @param tid_ Zero for byte input, positive for a conversion worker
		 * @param alpt_ Shared thread list, filled before any member is started
		 */
		ProcessThread(final int tid_, final ArrayList<ProcessThread> alpt_){
			tid=tid_;
			setName("FastqStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}

		/** Dispatches the assigned role and marks success only when that role returns.
		 * Worker failures rethrow; input-reading failures are caught inside processBytes.
		 */
		@Override
		public void run(){
			//Process the reads
			if(tid==0){
				processBytes();
			}else{
				if(interleaved){
					makeReadsInterleaved();
				}else {
					makeReadsSingle();
				}
			}

			//Indicate successful exit status
			success=true;
			if(verbose){outstream.println("tid "+tid+" terminated.");}
		}

		/** Reads/enqueues byte batches, signals completion in finally, then joins workers.
		 * Joins skip this coordinator itself. Aggregates converted-read/base counters and
		 * unsuccessful-worker status after joining; terminal output can be consumed earlier.
		 */
		void processBytes(){
			//[stream/FastqStreamer#001 FIXED 2026-06-20]: try/finally GUARANTEES oqs.poison() even when processBytes0() throws (an I/O
			//error from bf.nextLine, OOM, etc.), so the workers in getInput() + the consumer in getOutput() wake instead of hanging;
			//the catch records errorState so nextList crashes LOUD via KillSwitch instead of silently truncating. Input-thread death
			//leaves NO ordered gap (chunks delivered in order, LAST reachable) so plain poison() suffices here (the worker side needs
			//setFinished(true), see makeReadsSingle/makeReadsInterleaved). Same OQS-thread-death fix as the greenlit SamStreamer/BamWriter#001.
			try{
				processBytes0();
			}catch(Throwable t){
				errorState=true;
				outstream.println("FastqStreamer: error reading "+fname+": "+t);
			}finally{
				oqs.poison();// Signal completion via OQS -- ALWAYS, so workers (getInput) + consumer (getOutput) wake
			}
			if(verbose){outstream.println("tid "+tid+" done with processBytes0 + poisoning.");}

			//Wait for completion of all threads
			boolean allSuccess=true;
			ThreadWaiter.waitForThreadsToFinish(alpt);
			for(ProcessThread pt : alpt){
				//Wait until this thread has terminated
				if(pt!=this){
					//Accumulate per-thread statistics
					readsProcessed+=pt.readsProcessedT;
					basesProcessed+=pt.basesProcessedT;
					allSuccess&=pt.success;
				}
			}
			if(verbose){outstream.println("tid "+tid+" noted all process threads finished.");}

			//Track whether any threads failed
			if(!allSuccess){errorState=true;}
			if(verbose){outstream.println("tid "+tid+" finished! Error="+errorState);}
		}

		/** Opens input and enqueues four-line records with original record/pair positions.
		 * Samples no records here. Checks incomplete records/pairs, normalizes only real
		 * '+' separator lines, and leaves other validation to conversion of retained records.
		 * Batch thresholds count quad records and an approximate two-times-bases byte size;
		 * interleaved pairs stay together and can exceed a threshold. Header bytes are omitted.
		 * Reads batch thresholds once on entry and folds the final backend-close result.
		 */
		private void processBytes0(){
			if(verbose){outstream.println("tid "+tid+" started processBytes.");}

			bf=ByteFile.makeByteFile(ffin);

			long listNumber=0;
			long reads=0;
			int bytes=0;
			final int sections=(interleaved ? 2 : 1);

			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<byte[][]> ln=new ListNum<byte[][]>(new ArrayList<byte[][]>(slimit), listNumber++);
			ln.firstRecordNum=reads;

			while(reads<maxReads){
				{
					// Read 4 lines per FASTQ record
					byte[] header=bf.nextLine();
					if(header==null){break;}
					byte[] bases=bf.nextLine();
					byte[] plus=bf.nextLine();
					byte[] quals=bf.nextLine();

					if(bases==null || plus==null || quals==null){
						System.err.println("Incomplete record at end of file: "+new String(header));
						errorState=true;
						break;
					}
					//STR-003: retain malformed separators for quadToReadVec instead of disguising them as '+'.
					if(plus.length>1 && plus[0]=='+'){plus=PLUS;}
					bytes+=2*bases.length;//Ignore header, usually short

					byte[][] record=new byte[][]{header, bases, plus, quals};
					ln.add(record);
				}
				if(interleaved){
					// Read 4 lines per FASTQ record
					byte[] header=bf.nextLine();
					if(header==null){
						System.err.println("Incomplete pair at end of file, pairnum "+reads);
						errorState=true;
						break;
					}
					byte[] bases=bf.nextLine();
					byte[] plus=bf.nextLine();
					byte[] quals=bf.nextLine();

					if(bases==null || plus==null || quals==null){
						System.err.println("Incomplete record at end of file: "+new String(header));
						errorState=true;
						break;
					}
					//STR-003: apply the same preservation rule to the second mate's separator.
					if(plus.length>1 && plus[0]=='+'){plus=PLUS;}
					bytes+=2*bases.length;//Ignore header, usually short

					byte[][] record=new byte[][]{header, bases, plus, quals};
					ln.add(record);
				}
				reads++;

				if(ln.size()>=slimit || bytes>=blimit){
					oqs.addInput(ln);
					ln=new ListNum<byte[][]>(new ArrayList<byte[][]>(slimit), listNumber++);
					ln.firstRecordNum=reads;
					bytes=0;
				}
			}

			if(verbose){outstream.println("tid "+tid+" ran out of input.");}
			if(ln.size()>0){
				oqs.addInput(ln);
			}
			ln=null;
			if(verbose){outstream.println("tid "+tid+" done reading bytes.");}
			//TODO: Probable bug - STR-004: this read-modify-write can overwrite a worker's concurrent true; trace queue publication before choosing a fix.
			errorState|=bf.close();//Fold the reader's error state (truncated/corrupt input) so it isn't silently dropped at the streamer boundary
			if(verbose){outstream.println("tid "+tid+" closed stream.");}
		}

		/** Converts selected single records with positional IDs and the configured pair marker.
		 * Publishes an output batch even if sampling retains nothing. Re-enqueues terminal
		 * poison for other workers; failures mark error, force queue completion and rethrow.
		 */
		void makeReadsSingle(){
			//[stream/FastqStreamer#001 FIXED 2026-06-20]: worker death (a throw in quadToRead/quadToReadVec) used to leave its ordered
			//job UNDELIVERED -> the ordered consumer blocked forever on the gap. The catch force-poisons outq via oqs.setFinished(true)
			//(JobQueue.poisoned -> take() returns null past the gap) + records errorState -> consumer wakes + crashes LOUD in nextList.
			//NOT oqs.poison() (it sets lastSeen, tripping addInput's assert in the still-reading input thread). Rethrow so run() skips success=true.
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}

			try{
				ListNum<byte[][]> list=oqs.getInput();
				while(list!=null && !list.poison()){
					if(verbose){outstream.println("tid "+tid+" grabbed blist "+list.id());}

					ListNum<Read> reads=new ListNum<Read>(new ArrayList<Read>(list.size()), list.id());
					long readID=list.firstRecordNum;

					if(samplerate>=1f){
						for(byte[][] quad : list){
							Read r=quadToRead(quad, pairnum, readID++);
							reads.add(r);
						}
					}else{
						//Positional sampling (Streamer.sampleKeep): decision + numericID both come from the
						//record's file position, so R1/R2 twins keep matching subsets with matching ids
						//regardless of thread scheduling. The old shared-PRNG decision desynced pairs under MT.
						for(byte[][] quad : list){
							if(Streamer.sampleKeep(readID, sampleSeed, samplerate)){
								Read r=quadToRead(quad, pairnum, readID);
								reads.add(r);
							}
							readID++;
						}
					}

					oqs.addOutput(reads);
					list=oqs.getInput();
				}
				if(verbose){outstream.println("tid "+tid+" done making reads.");}
				//Re-inject poison for other workers
				if(list!=null){oqs.addInput(list);}
			}catch(Throwable t){
				errorState=true;
				oqs.setFinished(true);//force-poison outq -> release the consumer past THIS worker's undelivered-job gap
				throw new RuntimeException("FastqStreamer worker "+tid+" failed: "+fname, t);
			}
		}

		/** Converts selected adjacent pairs, links mates and publishes their first reads.
		 * Both mates share the input-pair numeric ID. An odd trailing quad records an error;
		 * complete pairs before it are still processed. Sampling is by pair position.
		 * Terminal and failure handling match the single-record worker.
		 */
		void makeReadsInterleaved(){
			//[stream/FastqStreamer#001 FIXED 2026-06-20]: same worker-death gap fix as makeReadsSingle — a throw in quadToRead leaves an
			//undelivered ordered job -> consumer hang; the catch force-poisons outq via setFinished(true) + records errorState so the
			//consumer wakes + crashes LOUD in nextList. (The odd-quads branch below already sets errorState without throwing.)
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}

			try{
				ListNum<byte[][]> list=oqs.getInput();
				while(list!=null && !list.poison()){
					if(verbose){outstream.println("tid "+tid+" grabbed blist "+list.id());}

					ListNum<Read> reads=new ListNum<Read>(new ArrayList<Read>((list.size()+1)/2), list.id());
					//firstRecordNum already counts PAIRS in interleaved mode (processBytes0 increments reads
					//once per pair); the old /2 halved it again, giving overlapping numericID ranges between
					//batches (batch at pair 100 restarted ids at 50).
					long readID=list.firstRecordNum;
					ArrayList<byte[][]> quads=list.list;
//					assert((quads.size()&1)==0) : "Odd number of quads for interleaved list: "+quads.size();
					if((quads.size()&1)!=0){
						System.err.println("Incomplete pair near pairnum "+readID);
						errorState=true;
					}

					final int lim=(quads.size()|1)-1;
					if(samplerate>=1f){
						for(int i=0; i<lim; i+=2){
							byte[][] quad1=quads.get(i);
							byte[][] quad2=quads.get(i+1);
							Read r1=quadToRead(quad1, 0, readID);
							Read r2=quadToRead(quad2, 1, readID++);
							r1.mate=r2;
							r2.mate=r1;
							reads.add(r1);
						}
					}else{
						//Positional sampling by PAIR index; the pair is kept or dropped as a unit.
						for(int i=0; i<lim; i+=2){
							if(Streamer.sampleKeep(readID, sampleSeed, samplerate)){
								byte[][] quad1=quads.get(i);
								byte[][] quad2=quads.get(i+1);
								Read r1=quadToRead(quad1, 0, readID);
								Read r2=quadToRead(quad2, 1, readID);
								r1.mate=r2;
								r2.mate=r1;
								reads.add(r1);
							}
							readID++;
						}
					}

					oqs.addOutput(reads);
					list=oqs.getInput();
				}
				if(verbose){outstream.println("tid "+tid+" done making reads.");}
				//Re-inject poison for other workers
				if(list!=null){oqs.addInput(list);}
			}catch(Throwable t){
				errorState=true;
				oqs.setFinished(true);//force-poison outq -> release the consumer past THIS worker's undelivered-job gap
				throw new RuntimeException("FastqStreamer worker "+tid+" failed: "+fname, t);
			}
		}

		/** Decodes one retained nucleotide record, sets its pair marker and counts it.
		 * Explicitly validates when the Read constructor has not already done so.
		 * @param quad Header, bases, separator and encoded qualities
		 * @param pairnum Pair marker to assign
		 * @param id Original record position, or shared pair position when interleaved
		 * @return Converted and validated read
		 */
		private Read quadToRead(final byte[][] quad, final int pairnum, final long id){
//			Read r=FASTQ.quadToRead_slow(quad, false, null, readID, 0);

			Read r=FASTQ.quadToReadVec(quad, id, 0, fname);
			r.setPairnum(pairnum);

			if(!r.validated()){r.validate(true);}

			readsProcessedT++;
			basesProcessedT+=r.length();
			return r;
		}

		/** Converted individual reads retained by sampling; both mates are counted. */
		protected long readsProcessedT=0;
		/** Bases in successfully converted reads retained by this worker. */
		protected long basesProcessedT=0;
		/** True when the assigned role returns; coordinator catches input failures internally. */
		boolean success=false;
		/** Zero for the input coordinator; positive for conversion workers. */
		final int tid;

		/** Coordinator-only join list; null in workers and fully populated before startup. */
		ArrayList<ProcessThread> alpt;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	public final String fname;

	/** Primary input file */
	final FileFormat ffin;

	/** Byte backend opened by the coordinator; close may clear this reference. */
	private ByteFile bf;

	/** Shared input jobs and output batches, with output ordered by batch ID. */
	final OrderedQueueSystem<ListNum<byte[][]>, ListNum<Read>> oqs;

	/** Conversion-worker count, excluding the additional byte-input coordinator. */
	final int threads;
	/** Marker for noninterleaved input; zero is required when interleaved. */
	final int pairnum;
	/** Descriptor pairing mode captured at construction. */
	final boolean interleaved;

	/** Converted individual reads, aggregated by the coordinator after worker joins. */
	protected long readsProcessed=0;
	/** Converted bases, aggregated by the coordinator after worker joins. */
	protected long basesProcessed=0;

	/** Input-record or interleaved-pair limit before sampling; negative requests become Long.MAX_VALUE. */
	final long maxReads;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Soft quad-record count per input batch, read once when byte input begins. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Soft approximate byte threshold per input batch; counts twice each base length. */
	public static int TARGET_LIST_BYTES=262144;
	/** Default conversion-worker request before the constructor's global-thread clamp. */
	public static int DEFAULT_THREADS=2;
	/** Shared canonical separator replacing valid plus lines with trailing text. */
	private static final byte[] PLUS=new byte[]{(byte)'+'};

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	/** Shared producer/worker error flag; retains the deferred STR-004 update concern. */
	public boolean errorState=false;
	/** Live positional sampling rate; configure before starting workers. */
	private float samplerate=1f;
	/** Seed for positional sampling (Streamer.sampleKeep); resolved from setSampleRate's seed */
	private long sampleSeed=17;

}
