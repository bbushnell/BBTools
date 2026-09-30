package stream;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/** Single-worker FASTQ reader using buffered four-line record groups.
 * Builds Reads, optionally pairs adjacent records, then samples by input-position ID.
 * Configure sampling and shared parsing settings before start; consume batches until null
 * and inspect errorState. Sampling can leave empty nonterminal batches. Public statistics
 * describe scanned/constructed input, not necessarily the records emitted after sampling.
 * The fastqscan.sh launcher uses FastqScan, not this class's separate diagnostic main.
 * @author Brian Bushnell
 * @contributor Collei
 * @contributor Shinobu (documentation and formatting)
 * @date November 29, 2025
 */
public final class FastqScanStreamer implements Streamer{

	/** Drains input and prints statistics plus structural corruption diagnostics.
	 * The second argument controls slow Read validation, not thread count.
	 * Missing/extra arguments print usage but do not cause an immediate return.
	 * Checks both structural diagnostics and reader error status before returning.
	 * @param args Input path (optionally prefixed by in=) and optional validation boolean
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		if(args.length<1 || args.length>2){
			System.err.println("Usage: FastqScanStreamer filename");
		}
		String fname=args[0];
		if(args.length>1){Read.SKIP_SLOW_VALIDATION=!Parse.parseBoolean(args[1]);}
		while(fname.startsWith("-")){fname=fname.substring(1);}
		if(fname.startsWith("in=")){fname=fname.substring(3);}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		if(ff.stdin()){
			//Do nothing
		}else{
			File f=new File(fname);
			if(!f.isFile() || !f.canRead()){
				throw new RuntimeException("Can't read "+fname);
			}
		}
		FastqScanStreamer fqs=new FastqScanStreamer(ff, 0, -1);
		fqs.start();
		for(ListNum<Read> ln=fqs.nextList(); ln!=null; ln=fqs.nextList()){
			//Do nothing
		}
		fqs.close();
		t.stop("Time:   \t");
		System.err.println("Records:\t"+fqs.totalRecords);
		System.err.println("Bases:  \t"+fqs.totalBases);
		if(ff.samOrBam()){System.err.println("Headers:\t"+fqs.totalHeaders);}
		ByteBuilder bb=fqs.corruption();
		//[FastqScanStreamer#003; STR-046] Worker failures need not set a structural corruption flag.
		//Include the reader error status in the exit decision after printing the existing diagnostics.
		if(fqs.slashrLines>0){
			System.err.println("Contained Windows-style \r\n");
		}
		if(bb!=null){System.err.print(bb);}
		if(bb!=null || fqs.errorState()){System.exit(1);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an unstarted reader, inferring input format and interleaving.
	 * @param fname Input name
	 * @param pairnum_ Pair number assigned before optional adjacent-record pairing
	 * @param maxReads_ Construction limit in reads or interleaved pairs; negative for unlimited
	 */
	FastqScanStreamer(String fname, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false), pairnum_, maxReads_);
	}

	/** Captures the input descriptor and construction limit without opening input.
	 * @param ff_ Input descriptor; interleaving is captured now
	 * @param pairnum_ Pair number assigned before optional adjacent-record pairing
	 * @param maxReads_ Construction limit in reads or interleaved pairs; negative for unlimited
	 */
	FastqScanStreamer(FileFormat ff_, int pairnum_, long maxReads_){
		ff=ff_;
		interleaved=ff.interleaved();
		pairnum=pairnum_;
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens input, groups buffered lines, constructs/samples reads and publishes batches.
	 * Counts complete buffered groups before limiting construction; residues survive fills.
	 * A final missing LF is synthesized and recorded. Caller handles errors and termination.
	 * @throws IOException If an input read fails
	 */
	void readFastq() throws IOException{
		//FIXED [stream/FastqScanStreamer bbnorm-race]: open the stream, then publish 'is' under the same monitor close()
		//uses. If close() already ran early (it saw is==null), we own cleanup of the just-opened stream and bail;
		//otherwise publish it so the synchronized close() closes it exactly once. No NPE, no leak, no double-close.
		final InputStream localIs=ReadWrite.getInputStream(ff.name(), false, false);
		synchronized(this){
			//[FastqScanStreamer#002; STR-045] Preserve reported cleanup errors even after early closure.
			if(closed){errorState|=ReadWrite.finishReading(localIs, ff.name(), ff.allowSubprocess()); return;}
			is=localIs;
		}
		IntList newlines=new IntList(8192);
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)){
			assert(bstart==0);
			r=Math.max(r, 0);
			bstop+=r;
			//bstop>0 here, so buffer[bstop-1] never underflows: the for-condition (r>0 || bstop>0) plus bstop+=max(r,0) guarantee it
			//(an empty file exits the loop before the body). On a 0-length/EOF read with no trailing '\n', append a synthetic terminal
			//newline so the last record still parses (expand() first if the buffer is exactly full).
			if(r==0 && buffer[bstop-1]!='\n'){
				if(bstop>=buffer.length){expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			final int records=(!interleaved ? newlines.size/4 : (newlines.size/4)&(~1));
			totalRecords+=records;
			final ArrayList<Read> reads;
			if(interleaved){
				reads=makeReadsInterleaved(records, newlines);
			}else{
				reads=makeReadsSingle(records, newlines);
			}
			if(reads!=null && !reads.isEmpty()){
				if(samplerate<1){sample(reads);}
				ListNum<Read> ln=new ListNum<Read>(reads, nextLID++);
				boolean b=add(ln);
				if(!b){break;}
			}
			
			final int residue=bstop-bstart;
			if(residue>0){
				if(bstart>0){
					System.arraycopy(buffer, bstart, buffer, 0, residue);
//				}else if(r>0){
				}else if(r>0 && bstop>=buffer.length){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1){
				if(residue>0){partialRecords++;}
				break;
			}
			if(nextRID>=maxReads){break;}
		}
	}
	
	/** Builds up to the remaining single-read limit and advances the consumed buffer offset.
	 * Converts qualities and applies configured Read validation before collecting structural flags.
	 * @param records Complete buffered record groups available
	 * @param newlines Offsets of LF bytes, four per record
	 * @return Constructed reads in input order, before sampling
	 */
	private ArrayList<Read> makeReadsSingle(int records, IntList newlines){
		int recordStart=0;
		long maxRecords=Math.min(records, maxReads-nextRID);
		records=(int)maxRecords;
		ArrayList<Read> reads=new ArrayList<Read>(records);
		final int offset=-FASTQ.ASCII_OFFSET;
		for(int i=0, j=0; i<records; i++, j+=4){
			final int headerEnd=newlines.get(j);
			final int basesEnd=newlines.get(j+1);
			final int plusEnd=newlines.get(j+2);
			final int recordEnd=newlines.get(j+3);
			final int slashr0=(buffer[headerEnd-1]=='\r') ? 1 : 0;
			final int slashr1=(buffer[basesEnd-1]=='\r') ? 1 : 0;
			final int slashr2=(buffer[recordEnd-1]=='\r') ? 1 : 0;
			final int headerLen=headerEnd-recordStart-1-slashr0;
			final int basesLen=basesEnd-headerEnd-1-slashr1;
			final int qualsLen=recordEnd-plusEnd-1-slashr2;
			
			final String header=new String(buffer, recordStart+1, headerLen);
			final byte[] bases=Arrays.copyOfRange(buffer, headerEnd+1, basesEnd-slashr1);
			final byte[] quals=Arrays.copyOfRange(buffer, plusEnd+1, recordEnd-slashr2);
			FASTQ.applyQualityOffset(quals, bases, offset);
			final Read r=new Read(bases, quals, header, nextRID);
			r.setPairnum(pairnum);
			reads.add(r);
			nextRID++;
			
			slashrLines+=slashr1+slashr2;
			totalBases+=basesLen;
			bstart=recordEnd+1;
			qualMismatch|=(qualsLen!=basesLen);
			missingAt|=(buffer[recordStart]!='@');
			missingPlus|=(buffer[basesEnd+1]!='+');
			recordStart=recordEnd+1;
		}
		return reads;
	}
	
	/** Constructs adjacent mates within the remaining pair limit, then links them.
	 * Mates share an input-position ID; names are not checked for pair agreement here.
	 * @param records Even number of complete buffered record groups available
	 * @param newlines Offsets of LF bytes, four per record
	 * @return First-mate list with linked second mates, before sampling
	 */
	private ArrayList<Read> makeReadsInterleaved(int records, IntList newlines){
		int recordStart=0;
		final long remainingRecords=(maxReads-nextRID);
		long maxRecords=(records<remainingRecords ? records : Math.min(records, remainingRecords*2));
		records=(int)maxRecords;
		ArrayList<Read> reads=new ArrayList<Read>(records);
		final int offset=-FASTQ.ASCII_OFFSET;
		for(int i=0, j=0; i<records; i++, j+=4){
			final int headerEnd=newlines.get(j);
			final int basesEnd=newlines.get(j+1);
			final int plusEnd=newlines.get(j+2);
			final int recordEnd=newlines.get(j+3);
			final int slashr0=(buffer[headerEnd-1]=='\r') ? 1 : 0;
			final int slashr1=(buffer[basesEnd-1]=='\r') ? 1 : 0;
			final int slashr2=(buffer[recordEnd-1]=='\r') ? 1 : 0;
			final int headerLen=headerEnd-recordStart-1-slashr0;
			final int basesLen=basesEnd-headerEnd-1-slashr1;
			final int qualsLen=recordEnd-plusEnd-1-slashr2;
			
			final String header=new String(buffer, recordStart+1, headerLen);
			final byte[] bases=Arrays.copyOfRange(buffer, headerEnd+1, basesEnd-slashr1);
			final byte[] quals=Arrays.copyOfRange(buffer, plusEnd+1, recordEnd-slashr2);
			FASTQ.applyQualityOffset(quals, bases, offset);
			final Read r=new Read(bases, quals, header, nextRID);
			r.setPairnum(pairnum);
			reads.add(r);
			nextRID+=(i&1);//pairs share a numericID: increments only after the odd (second) mate, so reads 0,1→id 0; 2,3→id 1; pair() then copies r1.numericID onto r2
			
			slashrLines+=slashr1+slashr2;
			totalBases+=basesLen;
			bstart=recordEnd+1;
			qualMismatch|=(qualsLen!=basesLen);
			missingAt|=(buffer[recordStart]!='@');
			missingPlus|=(buffer[basesEnd+1]!='+');
			recordStart=recordEnd+1;
		}
		pair(reads);
		return reads;
	}
	
	/** Links adjacent Reads and compacts second mates out of the top-level list.
	 * An odd leftover remains unpaired and produces a warning.
	 * @param reads Mutable list of individual reads
	 */
	private void pair(ArrayList<Read> reads){
		final int lim=reads.size()&~1;
		if(reads.size()!=lim){
			System.err.println("Warning: File "+ff.name()+" was processed"
				+ " as interleaved but had an odd number of reads.");
		}
		if(lim<1){return;}
		for(int i=0; i<lim; i+=2){
			Read r1=reads.get(i), r2=reads.set(i+1, null);
			r1.mate=r2;
			r2.mate=r1;
			r2.numericID=r1.numericID;
			r2.setPairnum(1);
		}
		Tools.condenseStrict(reads);
	}
	
	/** Retains input-position samples and compacts rejected top-level reads in place.
	 * @param reads Single reads or first mates with linked second mates
	 */
	private void sample(ArrayList<Read> reads){
		assert(samplerate<1);
		//Positional sampling by numericID (Streamer.sampleKeep): reproducible across runs and
		//thread counts, and pair-safe automatically since mates share a numericID.
		for(int i=0; i<reads.size(); i++){
			if(!Streamer.sampleKeep(reads.get(i).numericID, sampleSeed, samplerate)){
				reads.set(i, null);
			}
		}
		Tools.condenseStrict(reads);
	}
	
	/** Enqueues a batch while open; may block until a consumer makes room.
	 * @param ln Data or terminal batch; null is a no-op
	 * @return False only when an interrupted put observes closure; true also if already closed
	 * @throws RuntimeException If a put is interrupted while still open
	 */
	private boolean add(ListNum<Read> ln){
		while(ln!=null && !closed){
			try{
				queue.put(ln);
				ln=null;
			}catch(InterruptedException e){
				if(closed){return false;}
				throw new RuntimeException(e);
			}
		}
		return true;
	}

	/** Doubles the byte buffer up to Shared.MAX_ARRAY_LEN, asserting that growth is possible. */
	private void expand(){
		long newlen=Math.min(buffer.length*2L, Shared.MAX_ARRAY_LEN);
		assert(newlen>buffer.length) : "Record "+totalRecords+" is too long.";
		buffer=Arrays.copyOf(buffer, (int)newlen);
	}
	
	/** Describes structural anomalies, including a missing terminal LF.
	 * Excludes worker/cleanup failures and CRLF alone; consult errorState too.
	 * Counts are lower bounds, not a complete malformed-record census.
	 * @return New diagnostic buffer, or null when none of these flags are set
	 */
	public ByteBuilder corruption(){
		if(partialRecords<1 && !qualMismatch && !missingTerminalNewline && !missingPlus && !missingAt){
			return null;
		}
		ByteBuilder bb=new ByteBuilder();
		if(partialRecords>0 || missingAt || missingPlus || qualMismatch){
			bb.appendln("At least "+Math.max(partialRecords, 1)+" corrupt records.");
		}
		if(partialRecords>0){bb.appendln("At least "+partialRecords+" incomplete records.");}
		if(qualMismatch){bb.appendln("At least "+1+" base/quality mismatches.");}
		if(missingAt){bb.appendln("At least "+1+" missing @ symbols.");}
		if(missingPlus){bb.appendln("At least "+1+" missing + symbols.");}
		if(missingTerminalNewline){bb.appendln("Missing terminal newline.");}
		assert(bb.length()>0);
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Reader Thread         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Runs input processing then attempts terminal publication/closure, recording thrown failures. */
	private class ReaderRunnable implements Runnable{
		/** Prints caught failures and sets the outer error flag before attempting cleanup. */
		@Override
		public void run(){
			try{
				readFastq();
			}catch(Throwable e){
				e.printStackTrace();
				synchronized(FastqScanStreamer.this){errorState=true;}
			}
			try{
				poisonAndClose();
			}catch(Throwable e){
				e.printStackTrace();
				synchronized(FastqScanStreamer.this){errorState=true;}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Thread Control        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts one worker; subsequent calls are ignored. Does not reset a closed reader. */
	@Override
	public synchronized void start(){
		if(started){return;}
		started=true;
		readerThread=new Thread(new ReaderRunnable());
		readerThread.start();
		if(verbose){System.err.println("Started "+getClass().getName());}
	}

	/** Marks this reader closed and interrupts its worker; closes published input except bare System.in.
	 * Preserves earlier failures and folds any failure reported by input cleanup into errorState.
	 * Does not join the worker or independently publish a terminal batch.
	 */
	@Override
	public synchronized void close(){
		if(closed){return;}
		//FIXED [stream/FastqScanStreamer bbnorm-race]: set 'closed' FIRST so an interrupted reader's add() sees it and
		//returns cleanly (instead of rethrowing InterruptedException as a RuntimeException), and so a reader still
		//starting up cleans up its own stream. finishReading is now null-safe (ReadWrite#006); if 'is' isn't published
		//yet (early close races ahead of readFastq), readFastq's synchronized publish closes it instead — exactly once.
		closed=true;
		if(readerThread!=null){readerThread.interrupt();}
		//[FastqScanStreamer#002; STR-045] Always perform cleanup and retain its reported failure.
		errorState|=ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
	}
	
	/** Attempts to enqueue a terminal marker, then closes input; enqueueing may block. */
	public synchronized void poisonAndClose(){
		poison();
		close();
	}
	
	/** Attempts terminal publication once; closure can suppress the enqueue in add(). */
	private synchronized void poison(){
		if(poisoned){return;}
		add(new ListNum<Read>(null, nextLID++, ListNum.POISON));
		poisoned=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Overrides           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the configured input name. */
	@Override
	public String fname(){return ff.name();}

	/** Returns the interleaving setting captured at construction. */
	@Override
	public boolean paired(){return interleaved;}

	/** Returns the configured pair number before interleaved pairing. */
	@Override
	public int pairnum(){return pairnum;}

	/** Returns scanned record groups, divided by two for interleaved input.
	 * Groups are counted before limits and sampling; this need not equal output size.
	 */
	@Override
	public long readsProcessed(){return totalRecords/(interleaved ? 2 : 1);}

	/** Returns bases in constructed individual records before sampling, including both mates. */
	@Override
	public long basesProcessed(){return totalBases;}

	/** Configures positional sampling; call before starting the worker.
	 * @param rate Retention fraction; values below one activate sampling
	 * @param seed Sampling seed, resolved through Streamer.resolveSampleSeed
	 */
	@Override
	public synchronized void setSampleRate(float rate, long seed){
		samplerate=rate;
		sampleSeed=Streamer.resolveSampleSeed(seed);
	}

	/** Takes the next batch, possibly empty after sampling, or null at terminal.
	 * Normally requeues terminal markers for repeated EOF calls. Interrupted takes retry
	 * unless already drained; interrupted terminal reinsertion is not retried once ln is set.
	 * Call start before consuming.
	 * @return Next data batch, or null after a terminal marker is consumed
	 */
	@Override
	public ListNum<Read> nextList(){
		ListNum<Read> ln=null;
		while(ln==null){
			try{
				ln=queue.take();
				//re-inject the poison so it stays in the queue: makes nextList() idempotent across repeated end-of-stream calls
				//(and would release a second consumer); the poison's id sorts nowhere — it's a sentinel, not an ordered job.
				if(ln.poison()){queue.put(ln);}
			}catch(InterruptedException e){
				if(drained){return null;}
			}
		}
		if(ln.poison()){
			synchronized(this){
				drained=true;
				close();
			}
		}
		return ln.poison() ? null : ln;
	}

	/** Unsupported SAM-line view.
	 * @return Never returns normally
	 * @throws RuntimeException Unless an enabled interleaved-input assertion fails first
	 * @throws AssertionError If assertions are enabled and input is interleaved
	 */
	@Override
	public ListNum<SamLine> nextLines(){
		assert(!interleaved);
		throw new RuntimeException("Not supported");
	}

	/** Returns whether terminal consumption is still pending; does not inspect the queue. */
	@Override
	public boolean hasMore(){return !drained;}

	/** Returns caught failures, reported cleanup failures or structural record flags.
	 * Missing terminal LF alone and CRLF alone are excluded; this is not exhaustive validation.
	 */
	@Override
	public synchronized boolean errorState(){
		return errorState || partialRecords>0 || qualMismatch || missingPlus || missingAt;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Final Fields          ----------------*/
	/*--------------------------------------------------------------*/

	/** Input descriptor captured at construction. */
	private final FileFormat ff;
	/** Whether adjacent records form pairs. */
	private final boolean interleaved;
	/** Construction limit in single reads or pairs; independent of buffered group statistics. */
	private final long maxReads;
	/** Pair number initially assigned to each constructed Read. */
	private final int pairnum;
	/** Bounded handoff of data batches and a terminal marker. */
	private final ArrayBlockingQueue<ListNum<Read>> queue=new ArrayBlockingQueue<ListNum<Read>>(4);

	/*--------------------------------------------------------------*/
	/*----------------        Mutable Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Closure flag observed by input publication and enqueueing. */
	private volatile boolean closed=false;
	/** Whether start has launched a worker; never reset. */
	private volatile boolean started=false;
	/** Whether terminal publication has been attempted; add may skip it after closure. */
	private volatile boolean poisoned=false;
	/** Whether a consumer has observed the terminal marker. */
	private volatile boolean drained=false;
	/** Caught or reported worker/cleanup failure flag, accessed under this instance's monitor. */
	private boolean errorState=false;
	/** Initial byte-buffer capacity. */
	private int bufferLen=262144;
	/** Input bytes, with any unconsumed residue compacted to the front between fills. */
	private byte[] buffer=new byte[bufferLen];
	/** Published input stream; publication and closure coordinate through this instance's monitor. */
	private volatile InputStream is;//FIXED [stream/FastqScanStreamer bbnorm-race]: was non-volatile; close() could read a stale null mid-shutdown
	/** Single input/parsing worker, assigned by start. */
	private Thread readerThread;
	
	/** Consumed offset and exclusive valid-byte limit within buffer. */
	private int bstart=0, bstop=0;
	/** Next input-position ID, incremented per single read or interleaved pair. */
	private long nextRID=0;
	/** Next monotonically increasing batch or terminal ID. */
	private long nextLID=0;
	
	/** Sampling rate configured before start; values at least one bypass sampling. */
	private float samplerate=1f;
	/** Seed for positional sampling (Streamer.sampleKeep); resolved from setSampleRate's seed */
	private long sampleSeed=17;
	
	/*--------------------------------------------------------------*/
	/*----------------            Stats             ----------------*/
	/*--------------------------------------------------------------*/
	
	//TODO: Possible bug [stream/FastqScanStreamer#001] - vestigial: never incremented (this is a FASTQ-only scanner), so always 0;
	//main()'s `if(ff.samOrBam()) print(totalHeaders)` would print "Headers: 0" if ever reached. LOW/DOC. Remove the field + the
	//samOrBam print, or wire it up if SAM-header counting is intended. (Public field — grep for external readers before deleting.)
	/** Legacy header counter; this FASTQ implementation never increments it. */
	public long totalHeaders;
	/** Complete buffered four-line groups before construction limits or sampling; even groups when interleaved. */
	public long totalRecords;
	/** Bases in constructed records before sampling; includes both mates. */
	public long totalBases;
	
	/** Number of terminal nonempty-residue observations, not the exact number of incomplete records. */
	public long partialRecords;
	/** Sequence/quality lines ending in CRLF; excludes header and separator lines. */
	public long slashrLines;
	/** Whether an inspected record had differing sequence and quality lengths. */
	public boolean qualMismatch;
	/** Whether a terminal LF was synthesized for residual input. */
	public boolean missingTerminalNewline;
	/** Whether an inspected separator line lacked its leading plus sign. */
	public boolean missingPlus;
	/** Whether an inspected header lacked its leading at sign. */
	public boolean missingAt;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Compile-time optional startup diagnostics. */
	private static final boolean verbose=false;
	
}
