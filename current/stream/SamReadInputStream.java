package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import shared.KillSwitch;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/** ReadInputStream adapter that constructs and immediately starts a SAM/BAM Streamer.
 * Returns the delegate's Read lists without copying them; null and empty batches both
 * map to null. The wrapper exposes no sampling setter and does not link mates itself;
 * paired() always returns false. Reader selection and parsing options remain delegated.
 * Header publication uses JVM-wide shared state, not the unused per-instance header
 * field. Configure readers before use; this adapter does not synchronize its ordinary
 * iteration/close methods or provide a general concurrent-lifecycle guarantee.
 * @author Brian Bushnell, Shinobu
 * @contributor Isla
 * @date Original, refactored October 23, 2025
 */
public class SamReadInputStream extends ReadInputStream{

	/** Legacy throughput driver: counts list entries and sums each entry's pairLength.
	 * Uses the filename constructor with header publication disabled, then closes
	 * the adapter and aborts on reported errors before printing completion statistics.
	 * @param args Input filename in args[0]
	 */
	public static void main(final String[] args){
		SamReadInputStream sris=new SamReadInputStream(args[0], false, true, -1, -1);

		Timer t=new Timer();
		long reads=0, bases=0;
		for(ArrayList<Read> ln=sris.nextList(); ln!=null; ln=sris.nextList()){
			for(Read r : ln){bases+=r.pairLength();}
			reads+=ln.size();
		}
		//STR-021: a delegate may report failure only through its final error flag and terminal batch.
		//Reject that status before printing a normal completion summary.
		if(sris.close()){
			KillSwitch.kill("Error reading SAM/BAM file: "+sris.fname());
		}
		t.stop();
		System.err.println();
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves a SAM-fallback input descriptor and delegates with the default thread hint.
	 * Construction starts the selected reader; header publication may finish later.
	 * @param fname Input filename or standard-input name
	 * @param loadHeader_ Request publication of this input's header in shared state
	 * @param allowSubprocess_ Permit subprocess-assisted input through the descriptor
	 * @param maxReads_ Limit forwarded to the selected reader; negative means unlimited
	 */
	public SamReadInputStream(final String fname, final boolean loadHeader_,
			final boolean allowSubprocess_, final long maxReads_){
		this(fname, loadHeader_, allowSubprocess_, -1, maxReads_);
	}

	/** Resolves an input descriptor without content probing, then starts its reader.
	 * Header retention controls publication, not all parsing required by the format.
	 * @param fname Input filename or standard-input name; SAM is the fallback format
	 * @param loadHeader_ Request shared-header publication; false does not clear an old header
	 * @param allowSubprocess_ Permit subprocess-assisted input
	 * @param threads_ Reader-selection hint; negative selects the factory default
	 * @param maxReads_ Limit forwarded unchanged; units follow the selected reader
	 */
	public SamReadInputStream(final String fname, final boolean loadHeader_,
			final boolean allowSubprocess_, final int threads_, final long maxReads_){
		this(FileFormat.testInput(fname, FileFormat.SAM, null, allowSubprocess_, false),
			loadHeader_, threads_, maxReads_);
	}

	/** Captures configuration, constructs the delegate and starts it immediately.
	 * Requests ordered Read output. Unexpected formats warn but are still passed to
	 * the factory. Sets the shared input-present marker before constructing the reader,
	 * even when header publication is disabled; that marker is not proof of a publisher.
	 * @param ff Nonnull input descriptor, including format and subprocess permissions
	 * @param loadHeader_ Request shared-header publication by the selected reader
	 * @param threads_ Reader-selection hint; negative selects the factory default
	 * @param maxReads_ Limit forwarded unchanged; negative means unlimited
	 */
	public SamReadInputStream(final FileFormat ff, final boolean loadHeader_,
			final int threads_, final long maxReads_){
		SAM_INPUT_PRESENT=true; //Preserve the input-present marker even when loadHeader_ is false; publication is separate.
		loadHeader=loadHeader_;
		stdin=ff.stdio();

		if(!ff.samOrBam()){
			System.err.println("Warning: Did not find expected sam file extension for filename "+
				ff.name());
		}

		//Create streamer with appropriate thread count
		streamer=StreamerFactory.makeSamOrBamStreamer(ff, threads_, loadHeader_, true, maxReads_, true);

//		//Extract header if requested
//		if(loadHeader){
//			header=streamer.header;
//			if(header!=null){setSharedHeader(header);}
//		}
		streamer.start();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the delegate's availability hint without consuming a batch.
	 * Use nextList for iteration; this hint is not a universal EOF or shutdown guarantee.
	 * @return Currently reported availability hint
	 */
	@Override
	public boolean hasMore(){
		return streamer.hasMore();
	}

	/** Returns the next nonempty delegate batch's backing list without copying it.
	 * A null or empty batch maps directly to null; this method does not skip empty
	 * batches to search for later data. Normal construction leaves delegate sampling
	 * disabled, and this adapter exposes no sampling setter.
	 * @return Delegate list, or null for a null/empty delegate batch
	 */
	@Override
	public ArrayList<Read> nextList(){
		ListNum<Read> ln=streamer.nextList();
		return ln==null || ln.isEmpty() ? null : ln.list;
	}

	/** Closes the delegate, then folds its observed error flag into the local flag.
	 * Completion/join guarantees depend on the delegate; this wrapper adds no join.
	 * Does not clear shared headers or the shared input-present marker.
	 * @return Local error flag after including the delegate's reported state
	 */
	@Override
	public boolean close(){
		streamer.close();
		errorState|=streamer.errorState();//#001 fix: was returning the never-set LOCAL errorState, dropping the streamer's. SamStreamer/BamStreamer set errorState=true when a reader thread fails (e.g. a truncated/corrupt BAM), so without this a failed SAM/BAM read was reported as SUCCESS.
		return errorState;
	}

	/** Reports the local flag or the delegate's current error state, without waiting.
	 * #001 fix: the inherited base getter masked delegate errors in the CRIS error chain.
	 * @return true if either currently observed flag reports an error
	 */
	@Override
	public boolean errorState(){return errorState || streamer.errorState();}

	/** Rejects restart; construct a new adapter to read the source again.
	 * @throws RuntimeException Always, because restart is not implemented
	 */
	@Override
	public synchronized void restart(){
		throw new RuntimeException("SamReadInputStream does not support restart.");
	}

	/** Returns the input filename reported by the underlying streamer.
	 * @return Input filename */
	@Override
	public String fname(){return streamer.fname();}

	/** Returns false; this adapter does not advertise linked paired batches.
	 * SAM pairing metadata in individual Read/SamLine objects is a separate concern.
	 * @return false
	 */
	@Override
	public boolean paired(){return false;}

	/*--------------------------------------------------------------*/
	/*----------------      Shared Header Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the shared header reference, optionally waiting for a nonnull publication.
	 * No copy is made. Returns immediately when wait is false, a header already exists,
	 * or no input-present marker has been set. An empty nonnull list counts as published.
	 * Otherwise waits in 100ms intervals, printing caught interruptions and continuing.
	 * The caller must arrange a publisher: constructing with loadHeader=false, opening
	 * an input unsuccessfully or merely setting the marker does not guarantee a header.
	 * There is no timeout when the marker is set but no nonnull header arrives.
	 * @param wait Request waiting when input has been marked present and no header exists
	 * @return Shared header, possibly empty; null on an immediate unavailable return
	 */
	public static synchronized ArrayList<byte[]> getSharedHeader(final boolean wait){
		if(!wait || SHARED_HEADER!=null){return SHARED_HEADER;}
		//Crash-loud, never hang [SamReadInputStream#002]: if no SAM/BAM input stream was ever opened, no
		//reader will ever call setSharedHeader, so waiting here deadlocks (e.g. fastq/fasta -> sam, the
		//common case). Return null (unavailable) and let the caller generate a header instead. This also
		//resolves the original //TODO for the no-shared-header case.
		if(!SAM_INPUT_PRESENT){return SHARED_HEADER;}
		if(printHeaderWait){System.err.println("Waiting on header to be read from a sam file.");}
		while(SHARED_HEADER==null){
			try{
				SamReadInputStream.class.wait(100);
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
		return SHARED_HEADER;
	}

	/** Replaces the JVM-wide shared header reference and notifies all header waiters.
	 * Makes no defensive copy and accepts null, which clears the reference without
	 * resetting the input-present marker. Waiters continue while the header is null.
	 * @param list Header lines to publish, an empty published header, or null to clear
	 */
	public static synchronized void setSharedHeader(final ArrayList<byte[]> list){
		SHARED_HEADER=list;
		SamReadInputStream.class.notifyAll();
	}

	/** Marks that a SAM/BAM input source now exists in this JVM, so getSharedHeader(true) may legitimately
	 * block until that input parses and sets the shared header. Historically this was set in the
	 * SamReadInputStream constructor, but the native BamStreamer/SamStreamer streamers bypass this
	 * class entirely, leaving the flag false; getSharedHeader(true) then hit the #002 no-hang gate and
	 * returned null instead of waiting, racing callers like var2/ScafMap.loadSamHeader (assert header!=null).
	 * StreamerFactory calls this synchronously, on the constructing thread, whenever it builds a native
	 * SAM/BAM streamer with saveHeader=true (requesting publication by that streamer), so the
	 * gate invariant is visible before the streamer's worker thread starts. This sets only
	 * the marker; it neither opens input nor guarantees a later successful publication. */
	public static void markSamInputPresent(){SAM_INPUT_PRESENT=true;}

	/** Truncates the first tab-delimited SN: field at its first non-tab whitespace.
	 * Requires only an @SQ prefix. Preserves subsequent tab-delimited fields and does
	 * not mutate the input array. An unchanged line, including null/non-@SQ input,
	 * retains its original reference. Missing SN: asserts, then returns the original
	 * line when assertions are disabled; this helper is not a full SAM header validator.
	 * @param line Header bytes, or null
	 * @return Original reference if unchanged, otherwise an exactly sized trimmed array
	 * @throws AssertionError With assertions enabled, if an @SQ-prefixed line lacks SN:
	 */
	public static byte[] trimHeaderSQ(final byte[] line){
		if(line==null || !Tools.startsWith(line, "@SQ")){return line;}

		final int idx=Tools.indexOfDelimited(line, "SN:", 2, (byte)'\t');
		if(idx<0){
			assert(false) : "Bad header: "+new String(line);
			return line;
		}

		int trimStart=-1;
		for(int i=idx; i<line.length; i++){
			final byte b=line[i];
			if(b=='\t'){return line;}
			if(Character.isWhitespace(b)){
				trimStart=i;
				break;
			}
		}
		if(trimStart<0){return line;}

		final int trimStop=Tools.indexOf(line, (byte)'\t', trimStart+1);
		final int bbLen=trimStart+(trimStop<0 ? 0 : line.length-trimStop);
		final ByteBuilder bb=new ByteBuilder(bbLen);
		for(int i=0; i<trimStart; i++){bb.append(line[i]);}
		if(trimStop>=0){
			for(int i=trimStop; i<line.length; i++){bb.append(line[i]);}
		}
		assert(bb.length==bbLen) : bbLen+", "+bb.length+", idx="+idx+", trimStart="+
			trimStart+", trimStop="+trimStop+"\n\n"+new String(line)+"\n\n"+bb+"\n\n";

		return bb.array;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Unused legacy per-instance placeholder; active header publication uses SHARED_HEADER. */
	private ArrayList<byte[]> header=null;

	/** Factory-selected delegate, started during construction; threading follows its implementation. */
	private final Streamer streamer;
	/** Captured publication request; forwarded at construction and otherwise unused locally. */
	private final boolean loadHeader;

	/** Descriptor standard-stream flag captured at construction; does not imply restart support. */
	public final boolean stdin;

	/*--------------------------------------------------------------*/
	/*----------------        Static Settings       ----------------*/
	/*--------------------------------------------------------------*/

	/** JVM-wide shared header reference, not owned by any one input; contents are not copied. */
	private static volatile ArrayList<byte[]> SHARED_HEADER;
	/** Set by construction or markSamInputPresent to enable getSharedHeader's wait path.
	 * Does not prove successful opening or header publication. The false value preserves
	 * the #002 immediate-return gate for applications that have not marked a SAM input.
	 */
	static volatile boolean SAM_INPUT_PRESENT=false;
	/** Prints one diagnostic before a blocking shared-header wait begins. */
	public static boolean printHeaderWait=false;

}
