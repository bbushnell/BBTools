package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Batch reader that shreds FASTA sequence bytes into overlapping fragments as Reads.
 * Lines beginning with > delimit sequences; original header text is discarded.
 * Reads receive consecutive numeric IDs and matching decimal names, with no supplied
 * qualities or mate links. Read construction applies its configured validation.
 * Minimum/target lengths and batch limits are captured at construction; overlap is
 * read from the shared setting for each emitted fragment. Configure settings before use.
 * Instance state and mutable static settings do not provide a concurrent-use contract.
 *
 * @author Brian Bushnell
 * @date Feb 13, 2013
 */
public class FastaShredInputStream extends ReadInputStream{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Legacy diagnostic that prints sampled reads and timing to stdout.
	 * Takes a filename and optional print/scan counts, minimum length and target length.
	 * Reads consecutive entries through next(), including the complete initial batch (#002).
	 * This driver does not explicitly close the reader when it stops early.
	 * @param args Filename followed by up to four integer options
	 */
	public static void main(String[] args){

		int a=20, b=Integer.MAX_VALUE;
		if(args.length>1){a=Integer.parseInt(args[1]);}
		if(args.length>2){b=Integer.parseInt(args[2]);}
		if(args.length>3){MIN_READ_LEN=Integer.parseInt(args[3]);}
		if(args.length>4){TARGET_READ_LEN=Integer.parseInt(args[4]);}

		Timer t=new Timer();

		FastaShredInputStream fris=new FastaShredInputStream(args[0], false, false, Shared.bufferData());
		//#002-fix: Take the first Read through the cursor; transferring its whole list discarded the remainder.
		Read r=fris.next();
		int i=0;

		while(r!=null){
			if(i<a){System.out.println(r.toText(false));}
			r=fris.next();
			if(++i>=a){break;}
		}
		while(r!=null && i++<b){r=fris.next();}
		t.stop();
		System.out.println("Time: \t"+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves the filename with FASTA format forced, then opens a reader.
	 * @param fname Input name, including supported stdin names
	 * @param amino_ Set the amino-acid flag on generated Reads
	 * @param allowSubprocess_ Allow input decompression subprocesses
	 * @param maxdata Positive batch base target, or nonpositive for Shared.bufferData
	 */
	public FastaShredInputStream(String fname, boolean amino_, boolean allowSubprocess_, long maxdata){
		this(FileFormat.testInput(fname, FileFormat.FASTA, FileFormat.FASTA, 0, allowSubprocess_, false, false), amino_, maxdata);
	}

	/**
	 * Captures length/batch settings and opens input, then checks lengths via an assertion.
	 * Only the descriptor's name and subprocess permission are forwarded to ByteFile.
	 * Adding the last fragment can exceed the batch base target, checked between fragments.
	 * @param ff Nonnull input descriptor
	 * @param amino_ Set Read.AAMASK on generated Reads
	 * @param maxData_ Positive batch base target, or nonpositive for Shared.bufferData
	 */
	public FastaShredInputStream(FileFormat ff, boolean amino_, long maxData_){
		name=ff.name();
		amino=amino_;
		flag=(amino ? Read.AAMASK : 0);

		if(!fileIO.FileFormat.hasFastaExtension(name) && !name.startsWith("stdin")){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+name);
		}

		allowSubprocess=ff.allowSubprocess();
		minLen=MIN_READ_LEN;
		maxLen=TARGET_READ_LEN;
		maxData=maxData_>0 ? maxData_ : Shared.bufferData();

		bf=open();

		assert(settingsOK());
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns one buffered entry, reading ahead if needed and retaining the remaining entries.
	 * Counts only the returned Read as consumed; prefetched entries remain available after close.
	 * @return Next unpaired Read, or null when no entries remain
	 */
	public Read next(){
		//#002-fix: Keep a cursor instead of consuming a whole list and discarding all but its first entry.
		if(!hasMore()){return null;}
		assert(nextIndex>=0 && nextIndex<currentList.size()) : "Single-read selection requires an unread buffered entry; index="+nextIndex;
		final Read r=currentList.set(nextIndex, null);
		nextIndex++;
		consumed++;
		if(nextIndex==currentList.size()){
			currentList=null;
			nextIndex=0;
		}
		return r;
	}

	/**
	 * Fills if needed and transfers an untouched prefetched batch without copying.
	 * After single-read access, copies only remaining Read references into a new list.
	 * Counts every returned list entry as consumed; a batch may still exist after close.
	 * @return Next unpaired batch, or null when no entries remain
	 */
	@Override
	public ArrayList<Read> nextList(){
		if(currentList==null){
			boolean b=fillList();
		}
		ArrayList<Read> list=currentList;
		if(nextIndex>0){
			assert(list!=null && nextIndex<list.size()) : "A partial batch must retain unread entries before list transfer";
			list=new ArrayList<Read>(list.subList(nextIndex, list.size()));
		}
		currentList=null;
		nextIndex=0;
		if(list==null || list.isEmpty()){
			list=null;
		}else{
			consumed+=list.size();
		}
		return list;
	}

	/** Returns whether a batch is available, reading ahead when open and no batch is buffered. */
	@Override
	public boolean hasMore(){
		if(currentList==null || currentList.size()==0){
			if(open){
				fillList();
			}else{
			}
		}
		return (currentList!=null && currentList.size()>0);
	}

	/**
	 * Reopens the input, resets IDs/consumed count and discards prefetched reads/sequence bytes.
	 * Retains captured configuration, sequence-buffer capacity and the inherited error flag.
	 */
	@Override
	public void restart(){
		//#003-fix: Rewinding the file must also discard sequence bytes retained after a partial batch.
		//Otherwise restart emits the old suffix ahead of the reopened input despite resetting IDs.
		if(bf!=null){close();}
		assert(bf==null);
		consumed=0;
		nextReadID=0;
		currentList=null;
		nextIndex=0;
		buffer.clear();

		if(bf==null){
			bf=open();
		}else{
			assert(false) : "bf should be null";
		}
	}

	/**
	 * Closes and detaches the ByteFile, latching its returned status into errorState.
	 * Leaves prefetched reads and sequence bytes in place; repeated calls return immediately.
	 * @return false regardless of the stored error flag (#005)
	 */
	@Override
	public final boolean close(){
		//TODO: Probable bug #005 - the ByteFile close result is latched below but never returned.
		//Callers using only the ReadInputStream.close return value cannot observe that status.
		synchronized(this){
			if(!open){return false;}
			open=false;
			assert(bf!=null);
			errorState|=bf.close();
			bf=null;
		}
		return false;
	}

	/** Indicates whether this stream produces paired reads.
	 * @return false as FASTA shredding produces single-end reads */
	@Override
	public boolean paired(){return false;}

	/**
	 * Builds a batch up to the captured entry limit and soft base target.
	 * Closes the input when generation returns null, retaining any partial batch.
	 * An already-closed reader clears currentList without reading.
	 * @return true if the new batch has entries, false otherwise
	 */
	private final boolean fillList(){
		if(!open){
			currentList=null;
			return false;
		}
		assert(currentList==null);
		currentList=new ArrayList<Read>(BUF_LEN);

		long len=0;
		for(int i=0; i<BUF_LEN && len<maxData; i++){
			Read r=generateRead();
			if(r==null){
				close();
				break;
			}
			currentList.add(r);
			len+=r.length();
			//[stream/FastaShredInputStream#001] nextReadID is already advanced inside generateRead() when the
			//read is created; the extra ++ here double-incremented it, yielding non-contiguous IDs (0,2,4,...).
			if(verbose){System.err.println("Generated a read; i="+i+", BUF_LEN="+BUF_LEN);}
		}

		return currentList.size()>0;
	}

	/**
	 * Accumulates nonheader lines until the target, a header boundary or EOF.
	 * Emits a copied prefix when at least minLen bytes remain. Full target fragments
	 * retain the current overlap plus unused suffix; shorter fragments clear the buffer.
	 * Discards shorter-than-minimum remnants and continues to the next sequence.
	 * @return Generated Read with the next decimal ID, or null after exhaustion
	 */
	private final Read generateRead(){
		Read r=null;
		boolean eof=false;
		while(r==null && !eof){
			while(buffer.length()<maxLen){
				byte[] line=bf.nextLine();
				if(line==null){
					eof=true;
					break;
				}
				if(line.length>0 && line[0]==carrot){break;}
				buffer.append(line);
			}
			if(buffer.length>=minLen){
				byte[] bases=buffer.expelAndShift(Tools.min(maxLen, buffer.length()), TARGET_READ_OVERLAP);
				r=new Read(bases, null, Long.toString(nextReadID), nextReadID, flag);
				nextReadID++;
				if(verbose){System.err.println("Made read:\t"+(r.length()>1000 ? r.id : r.toString()));}
				if(bases.length<maxLen){buffer.clear();}
			}else{
				buffer.clear();
			}
		}
		return r;
	}

	/** Opens via the filename-based ByteFile factory and marks this reader open.
	 * @return Newly created byte reader
	 * @throws RuntimeException If already marked open */
	private final ByteFile open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened fasta file "+name);
		}
		open=true;
		ByteFile bf=ByteFile.makeByteFile(name, allowSubprocess);
		return bf;
	}

	/** Returns the local open flag without reading ahead. */
	public boolean isOpen(){return open;}

	/**
	 * Checks current static minimum/target lengths; constructor invocation is assertion-gated.
	 * Overlap is not validated here (#006); keep it nonnegative, below the target and
	 * no greater than any emitted fragment length. Settings must remain stable during use.
	 * @return Whether the minimum/target length checks pass
	 * @throws RuntimeException For an invalid minimum or target length
	 */
	public static final boolean settingsOK(){
		//TODO: Probable bug #006 - overlap bounds are not checked with the configured read lengths.
		//Overlap equal to a full emitted length retains the buffer unchanged; larger overlap violates
		//ByteBuilder.expelAndShift's assertion. Default settings avoid this; public settings can differ.
		if(MIN_READ_LEN>=Integer.MAX_VALUE-1){
			throw new RuntimeException("Minimum FASTA read length is too long: "+MIN_READ_LEN);
		}
		if(MIN_READ_LEN<1){
			throw new RuntimeException("Minimum FASTA read length is too short: "+MIN_READ_LEN);
		}
		if(TARGET_READ_LEN<1){
			throw new RuntimeException("Target FASTA read length is too short: "+TARGET_READ_LEN);
		}
		if(MIN_READ_LEN>TARGET_READ_LEN){
			throw new RuntimeException("Minimum FASTA read length is longer than maximum read length: "+MIN_READ_LEN+">"+TARGET_READ_LEN);
		}
		if(MIN_READ_LEN>=Integer.MAX_VALUE-1 || MIN_READ_LEN<1){return false;}
		if(TARGET_READ_LEN<1 || MIN_READ_LEN>TARGET_READ_LEN){return false;}
		return true;
	}

	/** Returns the attached ByteFile name; normal close/EOF detaches it (#004).
	 * @return Name supplied by the attached byte reader
	 * @throws NullPointerException After the byte reader has been detached */
	@Override
	public String fname(){
		//TODO: Probable bug #004 - close, including normal EOF, nulls bf although name is retained.
		return bf.name();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained descriptor name, still available after close. */
	public final String name;

	/** Prefetched list; consumed single-read slots are cleared, and close retains unread entries. */
	private ArrayList<Read> currentList=null;
	/** First unread list entry; zero when no entries remain or the list has been transferred. */
	private int nextIndex;

	/** Local input-open flag, cleared by close. */
	private boolean open=false;
	/** Sequence suffix retained between fragments; restart clears its logical contents. */
	private ByteBuilder buffer=new ByteBuilder();
	/** Attached byte reader; publicly exposed and nulled by close. */
	public ByteFile bf;

	/** Numeric ID assigned to the next generated fragment. */
	private long nextReadID=0;
	/** Entries returned through either single-read or list access. */
	private long consumed=0;

	/** Captured input subprocess permission. */
	public final boolean allowSubprocess;
	/** Captured amino-acid mode. */
	public final boolean amino;
	/** Read.AAMASK in amino mode, otherwise zero. */
	public final int flag;
	/** Maximum batch entries captured during instance initialization. */
	private final int BUF_LEN=Shared.bufferLen();
	/** Soft batch base target checked before each generated fragment. */
	private final long maxData;
	/** Captured target and minimum fragment lengths, respectively. */
	private final int maxLen, minLen;


	/** Enables diagnostic messages for fragment generation. */
	public static boolean verbose=false;
	/** First-byte FASTA header delimiter. */
	private static final byte carrot='>';

	/** Target fragment length copied into maxLen at construction. */
	public static int TARGET_READ_LEN=800;
	/** Shared overlap read on each fragment emission; configure before use. */
	public static int TARGET_READ_OVERLAP=31;
	/** Minimum emitted length copied into minLen at construction. */
	public static int MIN_READ_LEN=31;
	/** Legacy compatibility setting, not read by this class. */
	public static boolean FAKE_QUALITY=false;
	/** Legacy compatibility setting, not read by this class. */
	public static boolean WARN_IF_NO_SEQUENCE=true;
	/** Legacy compatibility setting, not read by this class. */
	public static boolean WARN_FIRST_TIME_ONLY=true;
	/** Legacy compatibility setting, not read by this class. */
	public static boolean ABORT_IF_NO_SEQUENCE=false;

}
