package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;

/**
 * Buffered legacy SCARF reader backed by a ByteFile and the shared FASTQ parser.
 * Opens the source during construction and reads batches on demand. Pairing is
 * captured from FASTQ.FORCE_INTERLEAVED; paired batches contain first mates with
 * their second mates attached. Counters and numeric IDs count those batch entries.
 * Returned lists and reads are handed to the caller without defensive copies.
 * <p>Use one owning thread: synchronization on selected methods does not protect
 * every lifecycle operation. Parser configuration/error state is shared with FASTQ.
 * Consumed incomplete pairs use the parser's process-fatal error path.
 */
public class ScarfReadInputStream extends ReadInputStream{

	/** Prints the first read from a nonempty SCARF file supplied as args[0]. */
	public static void main(final String[] args){
		final ScarfReadInputStream fris=new ScarfReadInputStream(args[0], true);
		final Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves the path with SCARF as the default format and opens the byte source.
	 * @param fname Input path
	 * @param allowSubprocess_ Whether format resolution may permit subprocess input
	 */
	public ScarfReadInputStream(final String fname, final boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.SCARF, null, allowSubprocess_, false));
	}

	/**
	 * Opens the input, captures pairing and buffer settings, and defers parsing until requested.
	 * A non-SCARF descriptor produces a warning but is still read through the SCARF parser.
	 * @param ff Nonnull input descriptor
	 */
	public ScarfReadInputStream(final FileFormat ff){
		if(verbose){System.err.println("ScarfReadInputStream("+ff.name()+")");}

		stdin=ff.stdio();
		if(!ff.scarf()){
			System.err.println("Warning: Did not find expected scarf file extension for filename "+ff.name());
		}

		tf=ByteFile.makeByteFile(ff);

		//SCARF pairing is driven solely by FORCE_INTERLEAVED, without per-file pairing auto-detection.
		interleaved=FASTQ.FORCE_INTERLEAVED;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks buffered entries, filling a batch when needed while the source is open.
	 * @return True when an unread entry is buffered
	 * @throws AssertionError if an exhausted buffer is checked after source closure
	 * with zero generated entries and assertions are enabled
	 */
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}

	/**
	 * Returns the current or next batch and relinquishes this reader's reference to it.
	 * Empty batches become null; consumed counts returned entries, not individual mates.
	 * The retained single-read cursor guard rejects a nonzero cursor.
	 * @return Next list, or null when the parser supplies an empty batch
	 */
	@Override
	public synchronized ArrayList<Read> nextList(){
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> list=buffer;
		buffer=null;
		if(list!=null && list.size()==0){list=null;}
		consumed+=(list==null ? 0 : list.size());
		return list;
	}

	/**
	 * Replaces an exhausted buffer with up to BUF_LEN entries and advances IDs/generated.
	 * A short batch closes the byte source. Parsing and quality handling are delegated
	 * to FASTQ.toScarfReadList; the legacy MAX_DATA byte threshold is not consulted.
	 */
	private synchronized void fillBuffer(){
		assert(buffer==null || next>=buffer.size());

		buffer=null;
		next=0;

		buffer=FASTQ.toScarfReadList(tf, BUF_LEN, nextReadID, interleaved);
		final int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){tf.close();}

		generated+=bsize;
		//defensive/unreachable: FASTQ.toScarfReadList always returns a (possibly empty) list, never null. Clean EOF -> empty buffer (bsize=0<BUF_LEN closes tf), and nextList() maps empty->null as the EOF signal.
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in ScarfReadInputStream.");
			}
		}
	}

	/**
	 * Closes the byte source and folds its close result into the local error flag.
	 * The shared FASTQ parser flag is additionally consulted by errorState(), not here.
	 * @return Accumulated local error flag after closing
	 */
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
	}

	/**
	 * Clears batch/counter/ID state and delegates resetting the source to ByteFile.reset().
	 * Retains pairing/buffer configuration and both local and shared parser error flags.
	 * Whether the source can be reopened is determined by the underlying byte reader.
	 */
	@Override
	public synchronized void restart(){
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		tf.reset();
	}

	/** Returns the FORCE_INTERLEAVED setting captured when this reader was constructed. */
	@Override
	public boolean paired(){return interleaved;}

	/**
	 * Combines the local flag with the process-wide FASTQ parser error flag.
	 * That shared flag can reflect errors from another reader.
	 * @return True when either flag is set
	 */
	@Override
	public boolean errorState(){return errorState || FASTQ.errorState();}//anti-swallow: ORs the static (process-global) FASTQ parser error state so a parse failure isn't dropped -> contrast SamReadInputStream#001. Over-reports rather than under (safe direction).

	/** Returns the underlying byte source's name. */
	@Override
	public String fname(){return tf.name();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Parsed entries waiting to be handed to the caller; null when no batch is held. */
	private ArrayList<Read> buffer=null;
	/** Retained single-read cursor; blockwise access leaves it at zero. */
	private int next=0;//vestigial: never incremented (blockwise-only reader via nextList) -> always 0, so the next-based single-read guards above are constant.

	/** Owned byte source; closed at a short batch or explicit close(). */
	private final ByteFile tf;
	/** Pairing choice captured from the global parser setting. */
	private final boolean interleaved;

	/** Maximum entries per batch, captured during construction. */
	private final int BUF_LEN=Shared.bufferLen();
	/** Legacy byte threshold, currently unused by the batching implementation. */
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long scarf reads.  Must be disabled for paired-ends.

	/** Parsed entry count since construction/restart; pairs count as one entry. */
	public long generated=0;
	/** Entry count handed to callers since construction/restart. */
	public long consumed=0;
	/** Numeric ID assigned to the first entry of the next parsed batch. */
	private long nextReadID=0;

	/** Whether the input descriptor identifies standard input. */
	public final boolean stdin;
	/** Enables shared diagnostic output for instances of this class. */
	public static boolean verbose=false;

}
