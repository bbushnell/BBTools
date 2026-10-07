package stream.bam;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Objects;

import fileIO.FileFormat;
import shared.Shared;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Presents BAM content as SAM text through the {@link InputStream} interface.
 * The implementation wraps {@link Streamer} so callers can treat a BAM
 * file as a plain-text stream of SAM header lines followed by SAM alignment
 * lines, one per newline.
 *
 * <p>This is intentionally inefficient (it round-trips through {@link SamLine})
 * but convenient for components that expect a text stream.</p>
 * Header text is obtained from JVM-wide shared state, not a per-stream header.
 * Instance buffers require a single reader or external synchronization; synchronized
 * close alone does not make concurrent reads and close safe. Delegate error flags
 * are not inspected by this facade.
 */
public class BamInputStream extends InputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens a BAM-default input with a worker hint capped at six.
	 * @param fname Input path
	 * @param ordered Whether to request ordered delegate batches
	 */
	public BamInputStream(String fname, boolean ordered){
		this(FileFormat.testInput(fname, FileFormat.BAM, null, true, false), ordered,
			Math.min(6, Shared.threads()/2+1));
	}

	/**
	 * Resolves a BAM-default input path and starts its selected reader.
	 * @param fname Input path
	 * @param ordered Whether to request ordered delegate batches
	 * @param threads Worker hint passed to the factory
	 */
	public BamInputStream(String fname, boolean ordered, int threads){
		this(FileFormat.testInput(fname, FileFormat.BAM, null, true, false), ordered, threads);
	}

	/**
	 * Starts a reader requesting shared header publication and SamLine batches,
	 * with no read limit and without conversion to Read objects.
	 * The caller supplies a descriptor supported by the factory's SAM-line interface.
	 * @param ff Nonnull input descriptor
	 * @param ordered Whether to request ordered delegate batches
	 * @param threads Worker hint passed to the factory
	 * @throws NullPointerException If ff is null
	 */
	public BamInputStream(FileFormat ff, boolean ordered, int threads){
		Objects.requireNonNull(ff, "FileFormat must not be null");
		//Header retention requests publication to SamReadInputStream's shared state;
		//makeReads=false requests SAM lines rather than Read objects.
		streamer=StreamerFactory.makeSamOrBamStreamer(ff, threads, true, ordered, -1L, false);
		streamer.start();
		headerQueue=new ArrayDeque<>();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the next unsigned SAM-text byte, or -1 after EOF or explicit close.
	 * @return Byte value from 0 through 255, or -1
	 * @throws IOException Declared by the InputStream API; delegate runtime exceptions propagate
	 */
	@Override
	public int read() throws IOException{
		if(!ensureBuffer()){return -1;}
		return currentBuffer[currentOffset++]&0xFF;
	}

	/**
	 * Copies SAM text across line boundaries until len bytes or delegate EOF.
	 * A valid zero-length request returns zero even after close. A positive request
	 * returns -1 when no byte is available; a final partial copy returns its count.
	 * @param b Destination array
	 * @param off First destination index
	 * @param len Maximum bytes to copy
	 * @return Copied count, zero for len=0, or -1 when exhausted or closed
	 * @throws NullPointerException If b is null
	 * @throws IndexOutOfBoundsException If the destination range is invalid
	 * @throws IOException Declared by the InputStream API; delegate runtime exceptions propagate
	 */
	@Override
	public int read(byte[] b, int off, int len) throws IOException{
		if(b==null){throw new NullPointerException();}
		if(off<0 || len<0 || len>b.length-off){throw new IndexOutOfBoundsException();}
		if(len==0){return 0;}
		if(!ensureBuffer()){return -1;}

		int copied=0;
		while(copied<len){
			if(currentOffset>=currentBuffer.length){
				if(!ensureBuffer()){break;}
			}
			int toCopy=Math.min(len-copied, currentBuffer.length-currentOffset);
			System.arraycopy(currentBuffer, currentOffset, b, off+copied, toCopy);
			currentOffset+=toCopy;
			copied+=toCopy;
		}

		return copied==0 ? -1 : copied;
	}

	/**
	 * Closes the delegate unless EOF or an earlier close has already set closed.
	 * On normal return, marks closed and drops the current batch and byte buffer.
	 * Does not clear the shared header or inspect the delegate's error state.
	 */
	@Override
	public synchronized void close(){
		if(closed){return;}
		streamer.close();
		closed=true;
		currentList=null;
		currentBuffer=null;
	}

	/**
	 * Supplies the next header or alignment line when buffered bytes are exhausted.
	 * Copies shared header references once, polling up to 100 times with 10ms sleeps
	 * before requesting the shared accessor's waiting behavior. That fallback has no
	 * facade timeout. Skips empty batches and null alignments; null batch marks EOF.
	 * @return Whether currentBuffer contains an unread byte
	 */
	private boolean ensureBuffer(){
		if(closed){return false;}
		if(currentBuffer!=null && currentOffset<currentBuffer.length){return true;}
		currentBuffer=null;
		currentOffset=0;

		//Header lines first.
		//TODO: Probable bug - stream/bam/BamInputStream#001: this JVM-wide header
		//has no association with this reader. A previously published nonnull header
		//can be consumed before this reader publishes its own; source concern only,
		//not reproduced here. Per-reader header ownership remains separate work.
		if(!headerLoaded){
			ArrayList<byte[]> sharedHeader=stream.SamReadInputStream.getSharedHeader(false);
			if(sharedHeader==null){
				for(int attempts=0; attempts<100 && sharedHeader==null; attempts++){
					try{
						Thread.sleep(10L);
					}catch(InterruptedException ignored){
						Thread.currentThread().interrupt();
					}
					sharedHeader=stream.SamReadInputStream.getSharedHeader(false);
				}
				if(sharedHeader==null){
					//Historical rationale assumed native readers left the input marker unset.
					//StreamerFactory now marks SAM/BAM inputs when saveHeader is true, so
					//this accessor may wait for publication; the polling loop is not a timeout.
					sharedHeader=stream.SamReadInputStream.getSharedHeader(true);
				}
			}
			if(sharedHeader!=null){headerQueue.addAll(sharedHeader);}
			headerLoaded=true;
		}
		if(!headerQueue.isEmpty()){
			byte[] headerLine=headerQueue.poll();
			if(headerLine!=null){
				currentBuffer=appendNewline(headerLine);
				return true;
			}
		}

		while(true){
			if(currentList==null || samIndex>=currentList.size()){
				currentList=streamer.nextLines();
				samIndex=0;
				if(currentList==null){
					closed=true;
					return false;
				}
				if(currentList.size()==0){continue;}
			}

			ArrayList<SamLine> samLines=currentList.list;
			if(samIndex<samLines.size()){
				SamLine sam=samLines.get(samIndex++);
				if(sam==null){continue;}
				ByteBuilder bb=sam.toText();
				bb.append('\n');
				currentBuffer=bb.toBytes();
				return true;
			}
		}
	}

	/**
	 * Copies a header line after removing trailing CR/LF, then appends exactly one LF.
	 * @param line Nonnull header bytes; the input array is not changed
	 * @return New normalized line bytes
	 */
	private static byte[] appendNewline(byte[] line){
		int len=line.length;
		while(len>0 && (line[len-1]=='\r' || line[len-1]=='\n')){len--;}
		ByteBuilder bb=new ByteBuilder(len+1);
		bb.append(line, 0, len);
		bb.append('\n');
		return bb.toBytes();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Factory-selected reader started during construction. */
	private final Streamer streamer;
	/** Shared header array references waiting for newline normalization. */
	private final ArrayDeque<byte[]> headerQueue;
	/** Whether shared header acquisition has been attempted. */
	private boolean headerLoaded=false;

	/** Current delegate batch, retained until consumed or explicitly closed. */
	private ListNum<SamLine> currentList=null;
	/** Next alignment index in the current batch. */
	private int samIndex=0;

	/** Current normalized header or rendered SAM line. */
	private byte[] currentBuffer=null;
	/** Next unread byte offset. */
	private int currentOffset=0;

	/** EOF observed or explicit close completed; both prevent further positive reads. */
	private boolean closed=false;
}
