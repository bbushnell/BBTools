package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Serializes SAM records on the submitting thread with per-batch buffering.
 * This class has no formatting queue or batch-ID reordering; WriterFactory selects
 * other implementations. Serialize submissions and lifecycle operations at the caller.
 * addReads/addLines, header initialization and counters are not synchronized. The lock
 * around individual byte writes does not make concurrent submissions safe.
 * Output opens during construction, with headers emitted by start or the first batch.
 * The selected output backend may buffer, compress or convert SAM text to BAM.
 * Do not concurrently mutate shared serialization settings or supplied records.
 *
 * @author Isla
 * @date November 10, 2025
 */
public class SamWriterST implements Writer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens output with both list entries and linked mates enabled.
	 * @param ffout_ Nonnull output descriptor
	 * @param header_ Optional explicit header lines, retained by reference
	 * @param useSharedHeader_ Prefer the shared input header over the explicit header
	 */
	public SamWriterST(FileFormat ffout_, ArrayList<byte[]> header_, boolean useSharedHeader_){
		this(ffout_, header_, useSharedHeader_, true, true);
	}

	/**
	 * Captures header-suppression settings and opens the selected output backend.
	 * BAM uses ReadWrite's BAM backend; other output uses the descriptor's append and
	 * subprocess settings. SAM header lines are submitted by start or a nonnull batch submission.
	 * Selection flags apply to Read conversion, not direct SamLine submission.
	 * @param ffout_ Nonnull output descriptor, opened immediately
	 * @param header_ Optional explicit header lines, retained by reference
	 * @param useSharedHeader_ Select shared input header before considering other sources
	 * @param writeR1_ Include each supplied list entry and its selected secondary records
	 * @param writeR2_ Include each entry's linked mate and its selected secondary records
	 */
	public SamWriterST(FileFormat ffout_, ArrayList<byte[]> header_, boolean useSharedHeader_,
			boolean writeR1_, boolean writeR2_){
		ffout=ffout_;
		fname=ffout.name();
		header=header_;
		useSharedHeader=useSharedHeader_;
		writeR1=writeR1_;
		writeR2=writeR2_;
		supressHeader=(ReadStreamWriter.NO_HEADER || (ffout.append() && ffout.exists()));
		supressHeaderSequences=(ReadStreamWriter.NO_HEADER_SEQUENCES || supressHeader);
		
		if(ffout.bam()){
			outstream=ReadWrite.getBamOutputStream(fname, ffout.append());
		}else{
			outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
		}
		if(verbose){outstream2.println("Made SamWriterST");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Emits the header if needed, then marks submissions started.
	 * Shared-header selection may wait for publication by the input side.
	 */
	@Override
	public void start(){
		writeHeader();
		started=true;
	}
	
	/** Wraps a nonnull Read list and submits it immediately.
	 * @param list Read entries; null entries are skipped
	 * @param id Accepted for the Writer API but ignored for ordering
	 */
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}

	/** Converts selected entries/mates to SamLines and writes the resulting batch.
	 * A null wrapper is ignored; an empty list can still trigger initial header output.
	 * @param reads Optional wrapper containing a nonnull list; its ID is ignored
	 */
	@Override
	public final void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		ArrayList<SamLine> lines=toSamLines(reads.list, writeR1, writeR2);
		writeLines(lines);
	}

	/** Writes the supplied SamLines directly, bypassing Read mate-selection flags.
	 * A null wrapper is ignored; null entries are skipped and the batch ID is ignored.
	 * @param lines Optional wrapper containing a nonnull list of SAM records
	 */
	@Override
	public final void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		writeLines(lines.list);
	}

	/** Starts lazily and serializes each nonnull record, flushing at the soft size threshold
	 * and at batch end. Counts records/bases before output succeeds, including secondaries.
	 * @param lines Nonnull records in caller-supplied order
	 */
	private void writeLines(ArrayList<SamLine> lines){
		if(!started){start();}
		
		ByteBuilder bb=new ByteBuilder();
		for(SamLine sl : lines){
			if(sl==null){continue;}
			sl.toBytes(bb);
			bb.nl();
			readsWritten++;
			basesWritten+=sl.length();

			if(bb.length()>=BUFFER_SIZE){
				write(bb);
			}
		}
		if(bb.length()>0){
			write(bb);
		}
	}
	
	/** Writes a copy of buffered bytes under this object's output lock and clears on success.
	 * IOException is propagated as RuntimeException without updating the cached error flag.
	 * @param bb Buffer to copy and clear
	 */
	private void write(ByteBuilder bb){
		if(bb.length()<0){return;}
		byte[] array=bb.toBytes();
		try{
			synchronized(this){outstream.write(array);}
			bb.clear();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}

	/** Marks end-of-input for the success flag; does not flush, close or reject later calls. */
	public final synchronized void poison(){
		poisoned=true;
	}

	/** Finalizes output once and caches the returned status; does not set the poison flag.
	 * Delegates stream/backend finalization to ReadWrite and does not emit an unwritten header.
	 * @return Cached error state, including the output finalizer's result
	 */
	public final synchronized boolean waitForFinish(){
		if(closed){return errorState;}
		boolean b=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
		closed=true;
		return (errorState|=b);
	}

	/** Marks end-of-input and finalizes output.
	 * @return Cached error state after finalization
	 */
	public final synchronized boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	/** Marks an external error, then uses the same poison/finalization path.
	 * This class has no own formatting backlog, but backend finalization can still block.
	 * Repeated calls retain the error flag and skip already-completed finalization.
	 */
	@Override
	public final synchronized void finishError(){
		errorState=true;
		poisonAndWait();
	}

	/** Returns the number of SamLines serialized, including secondary records.
	 * The count is updated before buffered output succeeds.
	 */
	@Override
	public long readsWritten(){return readsWritten;}

	/** Returns the sum of serialized SamLine.length values, including secondary records.
	 * The count is updated before buffered output succeeds.
	 */
	@Override
	public long basesWritten(){return basesWritten;}

	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Converts both list entries and linked mates using current shared SAM settings.
	 * @param reads Nonnull Read list; null entries are skipped
	 * @return New list whose SamLines may be reused attached records
	 */
	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads){
		return toSamLines(reads, true, true);
	}

	/**
	 * Creates or reuses primary SamLines, then appends the selected entries and secondaries.
	 * Selection is by list-entry/linked-mate position, not by filtering Read.pairnum.
	 * Both primary lines are prepared before selection. When KEEP_NAMES is false,
	 * mate qname normalization can mutate an attached SamLine even if that mate is excluded.
	 * @param reads Nonnull Read list; null entries are skipped
	 * @param writeR1 Include each list entry and its optional secondary alignments
	 * @param writeR2 Include each linked mate and its optional secondary alignments
	 * @return New ordered list; records are not guaranteed to be detached copies
	 */
	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads, boolean writeR1, boolean writeR2){
		ArrayList<SamLine> samLines=new ArrayList<SamLine>();

		for(final Read r1 : reads){
			if(r1==null){continue;}
			Read r2=(r1==null ? null : r1.mate);

			SamLine sl1=(r1==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE 
				&& r1.samline!=null ? r1.samline : new SamLine(r1, 0)));
			SamLine sl2=(r2==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE 
				&& r2.samline!=null ? r2.samline : new SamLine(r2, 1)));

			if(!SamLine.KEEP_NAMES && sl1!=null && sl2!=null && ((sl2.qname==null) || 
				!sl2.qname.equals(sl1.qname))){
				sl2.qname=sl1.qname;
			}
			assert(sl1!=null) : r1;
			if(writeR1){addSamLine(r1, sl1, samLines);}
			if(writeR2){addSamLine(r2, sl2, samLines);}
		}
		return samLines;
	}

	/** Appends a primary record and optional secondary records for sites at indexes one onward.
	 * Null Read/primary inputs are ignored; existing CIGAR/secondary assertions are preserved.
	 * @param r Read supplying secondary sites
	 * @param primary Primary SAM record
	 * @param samLines Destination list, appended in place
	 */
	private static void addSamLine(Read r, SamLine primary, ArrayList<SamLine> samLines){
		if(r==null || primary==null){return;}
		
		assert(!ReadStreamWriter.ASSERT_CIGAR || !r.mapped() || primary.cigar!=null) : r;
		samLines.add(primary);

		// Handle secondary alignments
		ArrayList<SiteScore> list=r.sites;
		if(ReadStreamWriter.OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
			final Read clone=r.clone();
			for(int i=1; i<list.size(); i++){
				SiteScore ss=list.get(i);
				clone.match=null;
				clone.setFromSite(ss);
				clone.setSecondary(true);
				SamLine secondary=new SamLine(clone, r.pairnum());
				assert(!secondary.nonSecondary());
				assert(!ReadStreamWriter.USE_ATTACHED_SAMLINE || secondary.cigar!=null) : r;
				samLines.add(secondary);
			}
		}
	}
	
	/** Selects shared, explicit or generated header lines in that precedence order.
	 * Shared lookup requests waiting. A null selected header warns and becomes empty;
	 * it does not fall through to another source. Sequence suppression affects generation only.
	 * @return Selected header reference or a newly generated/empty list
	 */
	ArrayList<byte[]> getHeader(){
		if(verbose){outstream2.println("Fetching header: "+useSharedHeader+","+(header!=null));}
		ArrayList<byte[]> headerLines;
		if(useSharedHeader){
			headerLines=SamReadInputStream.getSharedHeader(true);
		}else if(header!=null){
			headerLines=header;
		}else{
			headerLines=SamHeader.makeHeaderList(supressHeaderSequences, 
				ReadStreamWriter.MINCHROM, ReadStreamWriter.MAXCHROM);
		}
		if(headerLines==null){
			outstream2.println("Warning: Header was null, creating empty header");
			headerLines=new ArrayList<byte[]>();
		}
		if(verbose){outstream2.println("Fetched header: "+(headerLines==null ? "null" : headerLines.size()));}
		return headerLines;
	}

	/** Emits the selected header once unless all headers are suppressed.
	 * Uses its own 16KiB soft flush threshold. Empty selected headers still mark handling done.
	 * IOException is wrapped without setting the cached error flag; callers serialize this method.
	 */
	protected void writeHeader(){
		if(headerWritten || supressHeader){return;}
		ArrayList<byte[]> headerLines=getHeader();
		
		ByteBuilder bb=new ByteBuilder();
		try{
			for(byte[] line : headerLines){
				bb.append(line).nl();
				if(bb.length()>=16384){
					outstream.write(bb.toBytes());
					bb.clear();
				}
			}
			if(bb.length()>=1){
				outstream.write(bb.toBytes());
				bb.clear();
			}
		}catch(IOException e){
			throw new RuntimeException(e);
		}
		headerWritten=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Getters and Setters      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the cached finalization status or explicit error marker.
	 * Propagated write/header exceptions are not automatically folded into this field.
	 */
	@Override
	public boolean errorState(){return errorState;}
	
	/** Returns {@code poisoned && !errorState}; this does not check whether output is closed. */
	@Override
	public boolean finishedSuccessfully(){return !errorState && poisoned;}
	
	/** Returns the retained output name. */
	@Override
	public final String fname(){return fname;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Output name retained from the descriptor. */
	final String fname;
	/** Retained output descriptor used for backend selection and finalization. */
	final FileFormat ffout;
	/** Captured preference for shared input header lookup. */
	final boolean useSharedHeader;
	/** Captured full-header suppression, including append to an existing destination. */
	final boolean supressHeader;
	/** Captured suppression of generated reference-sequence header lines. */
	final boolean supressHeaderSequences;
	/** True after nonsuppressed header handling completes, even if it emitted no lines. */
	boolean headerWritten=false;
	/** Optional explicit header retained by reference; ignored when shared mode is selected. */
	final ArrayList<byte[]> header;
	/** Backend opened by the constructor and finalized by waitForFinish. */
	final OutputStream outstream;
	
	/** Serialized SAM records, including secondaries, counted before buffer output. */
	private long readsWritten=0;
	/** Sum of serialized record lengths, counted before buffer output. */
	private long basesWritten=0;
	/** Cached finalization status or marker set explicitly by finishError. */
	private boolean errorState=false;
	/** Include supplied list entries during Read conversion. */
	private final boolean writeR1;
	/** Include linked mates during Read conversion. */
	private final boolean writeR2;
	/** True after start completes normally. */
	private boolean started=false;
	/** End-of-input marker; not a submission guard. */
	private boolean poisoned=false;
	/** True after output finalization completes normally. */
	private boolean closed=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Soft data flush threshold checked after each serialized record. */
	private static final int BUFFER_SIZE=65536;
	/** Enables optional construction/header diagnostics. */
	public static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	protected PrintStream outstream2=System.err;

}
