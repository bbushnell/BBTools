package stream;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import structures.ByteBuilder;
import structures.ListNum;

/** FASTA reader that parses synchronously on the thread calling nextList.
 * Read names honor Shared.TRIM_READ_DESCRIPTION using legacy byte-whitespace semantics.
 * This class has no parsing worker or output queue; its ByteFile backend may use
 * input threads. StreamerFactory uses this fallback when its worker or SIMD
 * alternatives are not selected. Configure before start, start once on a fresh
 * instance, consume to terminal, then close. Restart is not implemented.
 *
 * Headers are decoded as ASCII and nonheader lines are concatenated into Reads
 * with null qualities and configured Read validation. Limits count fragments:
 * unpaired reads or complete interleaved pairs. Counters count constructed individual
 * reads/bases, including data not returned when a parsing call throws. Rejected
 * headers do not count. Parsing exceptions and assertions propagate synchronously;
 * they are not automatically latched into errorState. Close after such a failure.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @contributor Shinobu (limit guards, documentation and formatting)
 * @date November 10, 2025
 */
public class FastaStreamerZT implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an unstarted reader using FASTA input-format defaults.
	 * @param fname_ Input filename
	 * @param pairnum_ 0 for unpaired/R1, or 1 for a separate R2 input
	 * @param maxReads_ Fragment limit; negative for unlimited
	 */
	public FastaStreamerZT(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}

	/** Captures input format, pairing, amino flags and the fragment limit.
	 * Input opens later in start. Interleaved input requires pair number zero.
	 * @param ffin_ Input descriptor, including interleaving/amino settings
	 * @param pairnum_ 0 for unpaired/interleaved/R1, or 1 for a separate R2 input
	 * @param maxReads_ Unpaired-read or interleaved-pair limit; negative for unlimited
	 */
	public FastaStreamerZT(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		flag=(ffin.amino() || Shared.AMINO_IN ? Read.AAMASK : 0);
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(pairnum==0 || !interleaved);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		if(verbose){outstream.println("Made FastaStreamerZT");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resets counters and opens the ByteFile; call once before reading a fresh instance.
	 * Does not reset finished, header, list number or error state. Repeated start also
	 * replaces the reader reference without closing any previously open reader.
	 */
	@Override
	public synchronized void start(){
		if(verbose){outstream.println("FastaStreamerZT.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Open the file
		bf=ByteFile.makeByteFile(ffin);

		if(verbose){outstream.println("FastaStreamerZT started.");}
	}

	/** Closes the current ByteFile, folds its reported error status and clears the reference.
	 * Does not mark finished or clear pending header state. Close after terminal
	 * consumption or a thrown parsing failure; repeated close is a no-op.
	 */
	@Override
	public synchronized void close(){
		if(bf!=null){
			errorState|=bf.close();//Fold the reader's error state (truncated/corrupt input) so it isn't silently dropped at the streamer boundary
			bf=null;
		}
	}

	/** Returns the captured input name. */
	@Override
	public String fname(){return fname;}

	/** Reports whether terminal consumption is still pending; not a ByteFile availability check. */
	@Override
	public boolean hasMore(){return !finished;}

	/** Returns accumulated ByteFile-close status; thrown parsing failures are separate. */
	@Override
	public boolean errorState(){return errorState;}

	/** Returns the interleaving setting captured at construction. */
	@Override
	public boolean paired(){return interleaved;}

	/** Returns the configured input-side pair number. */
	@Override
	public int pairnum(){return pairnum;}

	/** Returns constructed individual reads, including both mates and unreturned data on failure. */
	@Override
	public synchronized long readsProcessed(){return readsProcessed;}

	/** Returns bases in constructed reads, including unreturned data on failure. */
	@Override
	public synchronized long basesProcessed(){return basesProcessed;}

	/** Configures PRNG header selection; call before start for the ordinary lifecycle.
	 * All interleaved rates below one, including zero, are rejected with assertions
	 * enabled. With assertions disabled, this mode does not guarantee pair alignment.
	 * This is header-level random selection, not positional Streamer.sampleKeep sampling.
	 * @param rate Requested retention threshold, normally in [0,1]
	 * @param seed Seed passed to Shared.random when sampling is active
	 * @throws AssertionError If assertions are enabled and an interleaved rate is below one
	 */
	@Override
	public synchronized void setSampleRate(float rate, long seed){
		//[stream/FastaStreamerZT#003] Crash-loud guard (Brian 2026-06-22): fractional sampling of INTERLEAVED FASTA desyncs read pairs (the start-read1 roll has no file-parity guard). Unsupported weird corner — best-effort only. -ea (default) crash loud with the workaround; -da silent best-effort. Shared with FastaStreamerST#001; crash-don't-corrupt.
		assert(!(interleaved && rate<1f)) : "Fractional sampling of interleaved FASTA is unsupported (read pairs would desync). Workaround: convert to FASTQ, subsample, then convert back to FASTA. ["+fname+"]";
		samplerate=rate;
		randy=(rate>=1f ? null : shared.Shared.random(seed));
	}

	/** Parses one batch inline, or closes and marks finished on a null/empty result.
	 * Repeated terminal calls return null. Calling before start also reaches terminal.
	 * Data batches advance listNum; fragment origins and constructed counters may differ
	 * from returned data after a thrown parsing failure. Exceptions/assertions escape
	 * without automatically setting errorState or marking finished; close after failure.
	 * @return Next nonempty batch, or null after terminal consumption
	 */
	@Override
	public synchronized ListNum<Read> nextList(){
		if(finished){return null;}

		ListNum<Read> list=interleaved ? nextListInterleaved() : nextListSingle();

		if(list==null || list.size()==0){
			finished=true;
			close();
			return null;
		}

		listNum++;
		return list;
	}

	/** Rejects the unsupported SAM-line view.
	 * @return Never returns normally
	 * @throws UnsupportedOperationException Always
	 */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTA does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs selected unpaired reads using a persistent header and a per-call builder.
	 * Batch boundaries retain the next header; sequence accumulation resumes next call.
	 * Stops at the read target or after exceeding the base target. Limits/IDs count
	 * constructed reads, not rejected headers; loop lookahead may consume beyond the cap.
	 * @return Batch, possibly empty, or null when no ByteFile is open
	 */
	private ListNum<Read> nextListSingle(){
		if(bf==null){return null;}

		ArrayList<Read> readList=new ArrayList<Read>(TARGET_LIST_SIZE);
		ListNum<Read> reads=new ListNum<Read>(readList, listNum);
		reads.firstRecordNum=readsProcessed;

		final ByteBuilder bb=new ByteBuilder(4096);

		int readsInList=0;
		int bytesInList=0;
		byte[] line=null;
		for(line=bf.nextLine(); line!=null && readsProcessed<maxReads; line=bf.nextLine()){

			if(line.length>0 && line[0]=='>'){
				if(header!=null){
					Read r=new Read(bb.toBytes(), null, new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1, 
						StandardCharsets.US_ASCII), readsProcessed, flag);
					r.setPairnum(pairnum);
					readList.add(r);
					readsProcessed++;
					basesProcessed+=r.length();
					readsInList++;
					bytesInList+=r.length();
				}
				header=null;
				bb.clear();
				
				if(samplerate>=1f || randy.nextFloat()<samplerate){header=line;}
				if(readsInList>=TARGET_LIST_SIZE || bytesInList>TARGET_LIST_BYTES){break;}
			}else if(header!=null){
				bb.append(line);
			}
		}
		//[STR-052] A carried lookahead header must not flush after the read limit, even on a later nextList call.
		if(line==null && header!=null && readsProcessed<maxReads) {//EOF
			Read r=new Read(bb.toBytes(), null, new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1, 
				StandardCharsets.US_ASCII), readsProcessed, flag);
			r.setPairnum(pairnum);
			readList.add(r);
			readsProcessed++;
			basesProcessed+=r.length();
			readsInList++;
			bytesInList+=r.length();
			header=null;
			bb.clear();
		}
		return reads;
	}

	/** Constructs adjacent records and returns first-mate entries with reciprocal links.
	 * With assertions disabled, EOF may also return an unmatched final record.
	 * Names are not checked for pair agreement. Limits and IDs count pairs, counters
	 * count individual reads, and batching thresholds are checked after complete pairs.
	 * With assertions enabled, an unmatched record before the limit throws synchronously.
	 * The sampling setter guards unsupported fractional interleaving before parsing.
	 * @return Batch, possibly empty, or null when no ByteFile is open
	 */
	private ListNum<Read> nextListInterleaved(){
		if(bf==null){return null;}

		ArrayList<Read> readList=new ArrayList<Read>(TARGET_LIST_SIZE);
		ListNum<Read> reads=new ListNum<Read>(readList, listNum);
		reads.firstRecordNum=readsProcessed/2;

		final ByteBuilder bb=new ByteBuilder(4096);
		long readID=readsProcessed/2;

		int readsInList=0;
		int bytesInList=0;

		Read pending=null;
		byte[] line=null;
		//[STR-051] Limit complete pairs while keeping readsProcessed in individual-read units.
		for(line=bf.nextLine(); line!=null && readID<maxReads; line=bf.nextLine()){

			if(line.length>0 && line[0]=='>'){
				if(header!=null){
					// Finish current read
					Read r=new Read(bb.toBytes(), null, new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1, 
						StandardCharsets.US_ASCII), 0, flag);
					readsProcessed++;
					basesProcessed+=r.length();
					bb.clear();

					if(pending==null){
						// This is read1
						pending=r;
						pending.setPairnum(0);
						header=line; // Start building read2
					}else{
						// This is read2
						r.setPairnum(1);
						pending.mate=r;
						r.mate=pending;
						pending.numericID=readID;
						r.numericID=readID++;

						readList.add(pending);
						readsInList+=2;
						bytesInList+=pending.length()+r.length();
						pending=null;

						// Decide whether to start building next pair
						header=(samplerate>=1f || randy.nextFloat()<samplerate) ? line : null;

						// Check if we should ship current list
						if(readsInList>=TARGET_LIST_SIZE || bytesInList>=TARGET_LIST_BYTES){
							break;
						}
					}
				}else{
					// Not currently building - decide whether to start
					header=(samplerate>=1f || randy.nextFloat()<samplerate) ? line : null;
					bb.clear();
				}
			}else if(header!=null){
				// Accumulate bases
				bb.append(line);
			}
		}
		
		if(line==null && header!=null && readID<maxReads) {//EOF; do not flush a lookahead header beyond the pair limit
			// Finish current read
			Read r=new Read(bb.toBytes(), null, new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1,
				StandardCharsets.US_ASCII), 0, flag);
			readsProcessed++;
			basesProcessed+=r.length();
			bb.clear();

			if(pending==null){
				// File ends on read1 - incomplete pair
				pending=r;
				pending.setPairnum(0);
				readList.add(pending);
			}else{
				// This is read2 - complete the pair
				r.setPairnum(1);
				pending.mate=r;
				r.mate=pending;
				pending.numericID=readID;
				r.numericID=readID++;

				readList.add(pending);
				readsInList+=2;
				bytesInList+=pending.length()+r.length();
				pending=null;
			}
			header=null;
		}

		// Handle incomplete pair at end
		//Historical intentional dual-mode behavior: the EOF pending==null branch adds a lone read1, then this
		//assert rejects the odd input under -ea. With assertions disabled that lone read can be returned best-effort.
		//Unlike ST's worker path, the assertion propagates directly to the caller; close in its failure handler.
		assert(pending==null) : "Odd number of reads in interleaved FASTA file: "+fname;

		return reads;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Captured input name. */
	public final String fname;

	/** Captured input descriptor used when start opens the reader. */
	final FileFormat ffin;

	/** Reader opened by start and cleared by close. */
	private ByteFile bf;

	/** Input-side pair marker; interleaving requires zero under assertions. */
	final int pairnum;
	/** Whether adjacent constructed records are paired. */
	final boolean interleaved;

	/** Constructed individual-read total, including unreturned data on failure. */
	protected long readsProcessed=0;
	/** Bases in constructed reads, including unreturned data on failure. */
	protected long basesProcessed=0;

	/** Fragment limit: unpaired reads or complete interleaved pairs. */
	final long maxReads;
	/** Flags passed to Read construction, initially including the requested amino flag. */
	public int flag;

	/** ID of the next nonempty returned batch; not reset by start. */
	private long listNum=0;

	/** Set only when nextList obtains a null/empty result; not set by close or a parse failure. */
	private boolean finished=false;
	
	/** Header line carried ACROSS nextList() calls: when a list fills mid-file, the next record's '>' header is already consumed into this field, so the following call resumes accumulating that record's sequence. bb is fresh per call; this header field is what persists the cross-call parse state (the whole reason a host-driven streamer can chunk output without losing a record boundary). */
	private byte[] header=null;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructed individual-read batching target, initialized from Shared. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Base target: unpaired batching checks greater than; paired batching checks greater than or equal. */
	public static int TARGET_LIST_BYTES=262144;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Destination for optional verbose messages. */
	protected PrintStream outstream=System.err;
	/** Compile-time optional diagnostics switch. */
	public static final boolean verbose=false;
	/** Accumulated ByteFile-close error status; synchronous parsing failures propagate separately. */
	public boolean errorState=false;
	/** Header-retention threshold configured before start. */
	private float samplerate=1f;
	/** PRNG used below rate one, assigned by setSampleRate. */
	private shared.Random randy=null;

}
