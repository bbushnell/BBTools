package stream.bam;

import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import shared.Shared;
import stream.FASTQ;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;
import stream.SamReadInputStream;
import structures.ByteBuilder;

/**
 * Caller-driven BAM reader that opens a named file and returns batches of Read objects.
 * Record parsing runs on the consuming caller; optional BGZF decompression can start
 * background threads. This reader is intended for one consuming caller and supports
 * neither restart nor individual-read access. The interleaved flag is reported as
 * metadata; this class does not attach mate links. Current EOF and close-status
 * limitations are documented at their implementation points below.
 *
 * @author Brian Bushnell
 * @date October 2025
 */
public class BamReadInputStreamST extends ReadInputStream{

	/**
	 * Prints the first read and its SAM representation, then closes on normal completion.
	 * @param args First element is a BAM filename containing at least one alignment
	 */
	public static void main(String[] args){
		BamReadInputStreamST bris=new BamReadInputStreamST(args[0], false, false, true);

		//TODO: Bug [stream/bam/BamReadInputStreamST#004] - nextList returns null for an
		//empty batch. This diagnostic dereferences it before close; handle empty input
		//and cleanup separately from the documentation pass (STR-191).
		Read r=bris.nextList().get(0);
		System.out.println(r.toText(false));
		System.out.println();
		if(r.samline!=null){
			System.out.println(r.samline.toText());
			System.out.println();
		}

		bris.close();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Detects the input format and immediately opens the named file.
	 * @param fname BAM filename passed to FileFormat and then FileInputStream
	 * @param loadHeader_ Whether to retain and publish the textual header
	 * @param interleaved_ Value reported by paired(); does not link mates here
	 * @param allowSubprocess_ Passed to format detection; this reader opens the file directly
	 */
	public BamReadInputStreamST(String fname, boolean loadHeader_, boolean interleaved_, boolean allowSubprocess_){this(FileFormat.testInput(fname, FileFormat.BAM, null, allowSubprocess_, false), loadHeader_, interleaved_);}

	/**
	 * Opens the file, selects BGZF decompression, and consumes the BAM header and dictionary.
	 * The constructor owns the opened streams. Header text is always consumed; when
	 * requested, nonempty newline-terminated lines are collected and optionally passed
	 * through trimHeaderSQ. A final unterminated fragment is retained without that trim.
	 * The collected header is published through SamReadInputStream.setSharedHeader.
	 * @param ff Input descriptor; its name is opened as a file, even if stdio() is true
	 * @param loadHeader_ Whether to retain and publish header text
	 * @param interleaved_ Paired-input metadata reported without mate-link construction
	 * @throws RuntimeException If the magic is not BAM or an opening/header IOException occurs
	 */
	public BamReadInputStreamST(FileFormat ff, boolean loadHeader_, boolean interleaved_){
		loadHeader=loadHeader_;
		interleaved=interleaved_;

		stdin=ff.stdio();
		if(!ff.bam()){System.err.println("Warning: Did not find expected bam file extension for filename "+ff.name());}

		fname=ff.name();
		header=new ArrayList<byte[]>();

		try{
			fis=new FileInputStream(fname);
			//Historical backend/call-site note below; its failure and severity claims have not
			//been revalidated here. This class's main is the only current in-tree constructor caller;
			//public external uses are not ruled out. The selected MT constructor starts threads.
			//Despite the "ST" name (no background threads in THIS class - it reads records inline in
			//fillBuffer), the underlying BGZF decompression is MULTIThreaded when USE_MULTITHREADED_BGZF: this
			//is the SOLE non-test constructor of BgzfInputStreamMT (the 702-line engine). So this reader is the
			//concrete path that would exhibit BgzfInputStreamMT#001 (producer-death → consumer HANG on corrupt/
			//truncated BAM). Both are TEST-ONLY (this reader's only caller is its own main; the live BAM read
			//path is BamStreamer→BgzfInputStreamMT2), so #001 stays LOW here too.
			if(BgzfSettings.USE_MULTITHREADED_BGZF){
				int threads=Math.max(1, BgzfSettings.READ_THREADS);
				bgzf=new BgzfInputStreamMT(fis, threads);
			}else{bgzf=new BgzfInputStream(fis);}
			reader=new BamReader(bgzf);

			// Read BAM magic
			byte[] magic=reader.readBytes(4);
			if(!Arrays.equals(magic, new byte[]{'B', 'A', 'M', 1})){throw new RuntimeException("Not a BAM file: "+fname);}

			// Read header text
			long l_text=reader.readUint32();
			byte[] text=reader.readBytes((int)l_text);

			// Parse header if requested
			if(loadHeader){
				int start=0;
				for(int i=0; i<text.length; i++){
					if(text[i]=='\n'){
						if(i>start){
							byte[] line=Arrays.copyOfRange(text, start, i);
							if(Shared.TRIM_READ_DESCRIPTION){line=SamReadInputStream.trimHeaderSQ(line);}
							header.add(line);
						}
						start=i+1;
					}
				}
				// Add last line if not ending with newline
				if(start<text.length){
					byte[] line=Arrays.copyOfRange(text, start, text.length);
					header.add(line);
				}
				SamReadInputStream.setSharedHeader(header);
			}

			// Read reference sequence dictionary
			int n_ref=reader.readInt32();
			refNames=new String[n_ref];
			for(int i=0; i<n_ref; i++){
				long l_name=reader.readUint32();
				refNames[i]=reader.readString((int)l_name-1); // Exclude NUL
				reader.readUint8(); // Skip NUL terminator
				long l_ref=reader.readUint32(); // Reference length (unused here)
			}

			converter=new BamToSamConverter(refNames);

		}catch(IOException e){throw new RuntimeException("Error opening BAM file: "+fname, e);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reports buffered availability, filling a batch if needed before finished is set.
	 * @return Whether an unread buffered entry is present
	 */
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			if(!finished){fillBuffer();}else{
				//TODO: Bug [stream/bam/BamReadInputStreamST#003] - finished with zero generated
				//records includes empty input. A repeated query then asserts instead of returning
				//false when assertions are enabled (STR-190); no runtime case was run here.
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}

	/**
	 * Transfers the current batch and counts its entries as consumed.
	 * Fills when the buffer is absent/exhausted without checking finished first;
	 * callers should not request another batch after close. Returned lists are no
	 * longer retained in buffer, and this method does not attach mates.
	 * @return Next batch, or null when the filled batch contains no records
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
	 * Reads up to Shared.bufferLen() record bodies and converts them to Read objects.
	 * Each read receives its SamLine and a monotonically increasing numeric ID.
	 * A caught EOFException sets finished regardless of which record read raised it.
	 * Generated counts include buffered entries that have not yet been returned.
	 */
	private synchronized void fillBuffer(){
		assert(buffer==null || next>=buffer.size());

		buffer=null;
		next=0;

		int BUF_LEN=Shared.bufferLen();
		buffer=new ArrayList<Read>(BUF_LEN);
		final ByteBuilder cigar=new ByteBuilder(1024);

		// Read alignment records until buffer full or EOF
		try{
			while(buffer.size()<BUF_LEN){
				long block_size=reader.readUint32();
				byte[] bamRecord=reader.readBytes((int)block_size);
				SamLine sl=converter.toSamLine(bamRecord, cigar);
				Read r=sl.toRead(FASTQ.PARSE_CUSTOM);
				r.samline=sl;
				r.numericID=nextReadID;
				buffer.add(r);

				nextReadID++;
			}
		}catch(EOFException e){
			//TODO: Probable bug [stream/bam/BamReadInputStreamST#005] - BamReader.readFully
			//uses EOFException for both a boundary EOF and a partial read. Catching all such
			//exceptions here does not distinguish an incomplete size/body from normal EOF
			//(STR-192). Source-only finding; no truncated-input probe or behavior repair here.
			//The older "Clean" and "acceptable" assessments below are historical, not new conclusions.
			// Normal end of file
			//EOF detection relies on BamReader.readFully throwing EOFException on ANY end-of-stream (even an
			//immediate 0-byte read). So readUint32(block_size) at the BGZF EOF marker throws → caught here →
			//finished. Clean. (A non-EOF IOException re-throws loud below. The hasMore() assert(generated>0)
			//"Was the file empty?" is a dev sanity guard that fires under -ea on a legitimately-empty BAM;
			//test-only, acceptable.)
			finished=true;
		}catch(IOException e){throw new RuntimeException("Error reading BAM file: "+fname, e);}

		generated+=buffer.size();
	}

	/**
	 * Sets finished and attempts to close the BGZF wrapper followed by the raw file.
	 * Does not discard a buffered batch or update the inherited errorState field.
	 * @return True when both close calls complete, false on a caught IOException;
	 * this is the opposite polarity of the ReadInputStream error-status contract
	 */
	@Override
	public boolean close(){
		//TODO: Bug [stream/bam/BamReadInputStreamST#002] - this returns success rather than
		//the parent API's error status, and caught errors are not latched in errorState.
		//Retain the current behavior until a separately verified correction (STR-189).
		finished=true;
		try{
			if(bgzf!=null){bgzf.close();}
			if(fis!=null){fis.close();}
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
		return true;
	}

	/**
	 * Rejects restart; this implementation does not reopen its file or reset counters.
	 * @throws RuntimeException Always, because this reader does not implement restart
	 */
	@Override
	public synchronized void restart(){throw new RuntimeException("BamReadInputStreamST.restart() not supported - BAM streams cannot be reset");}

	/** Returns the filename retained from the input descriptor. */
	@Override
	public String fname(){return fname;}

	/** Returns the constructor's interleaved flag; does not imply attached mate links. */
	@Override
	public boolean paired(){return interleaved;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Current unreturned batch; cleared when transferred by nextList. */
	private ArrayList<Read> buffer=null;
	/** Collected textual header lines when loadHeader is enabled. */
	private ArrayList<byte[]> header=null;
	/** Buffered-entry index, reset on filling; no individual-read accessor advances it here. */
	private int next=0;

	/** Raw file stream opened by the constructor. */
	private final FileInputStream fis;
	/** Selected BGZF wrapper, which may use background decompression. */
	private final InputStream bgzf;
	/** Little-endian scalar reader over the BGZF wrapper. */
	private final BamReader reader;
	/** Alignment decoder using refNames in dictionary order. */
	private final BamToSamConverter converter;
	/** Reference names read from the BAM binary dictionary; lengths are discarded here. */
	private final String[] refNames;

	/** Paired-input metadata supplied by the caller. */
	private final boolean interleaved;
	/** Whether textual header lines are retained and published. */
	private final boolean loadHeader;
	/** Filename opened directly through FileInputStream. */
	private final String fname;
	/** Set on caught EOF or close; nextList does not consult it before filling. */
	private boolean finished=false;

	/** Decoded record entries, including any still buffered. */
	public long generated=0;
	/** Record entries transferred through nextList. */
	public long consumed=0;
	/** Numeric ID assigned to the next decoded read. */
	private long nextReadID=0;

	/** Descriptor's stdio flag, retained as metadata without selecting a special input stream. */
	public final boolean stdin;

}
