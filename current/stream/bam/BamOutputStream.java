package stream.bam;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import parse.LineParser1;
import stream.SamLine;
import structures.ByteBuilder;

/**
 * {@link OutputStream} facade that accepts SAM text and emits BAM output.
 * Buffers leading header lines and parses alignment lines into {@link SamLine}
 * objects, then writes length-prefixed BAM records through a BGZF backend.
 * A newline, {@link #flush()}, or {@link #close()} completes the pending line;
 * callers must not flush in the middle of a SAM record. This mutable facade
 * requires a single caller or external synchronization.
 * <p>Backend close bypasses the requested {@code closeUnderlying} policy.
 * In particular, the ST backend closes the supplied stream even when the flag
 * is false; see the ownership TODO in {@link #close()}.
 */
public class BamOutputStream extends OutputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens a file for replacement using the configured BGZF defaults.
	 * @param filename Destination path; an existing file is truncated
	 * @throws IOException If the destination cannot be opened
	 */
	public BamOutputStream(String filename) throws IOException{
		this(new FileOutputStream(filename), true);
	}

	/**
	 * Opens a file for replacement with explicit compression and worker settings.
	 * @param filename Destination path; an existing file is truncated
	 * @param compression Compression level passed to the selected backend
	 * @param threads Requested workers; MT also requires the global MT flag
	 * @throws IOException If the destination cannot be opened
	 */
	public BamOutputStream(String filename, int compression, int threads) throws IOException{
		this(new FileOutputStream(filename), true, compression, threads);
	}

	/**
	 * Append-capable constructor. When appending to an existing nonempty BAM, the stream must NOT
	 * re-emit the BAM magic/header/ref-dict (a second header block mid-file is corrupt). Incoming
	 * header lines are still consumed to build the converter's reference dictionary.
	 * The caller must provide the same reference names and order as the existing file;
	 * this facade does not read or compare the existing dictionary.
	 * Mirrors stream.BamWriter's suppressHeader handling. (Previously ReadWrite.getBamOutputStream's
	 * native path ignored append entirely and truncated; 2026-09-05.)
	 * @param filename Destination path
	 * @param compression Compression level passed to the selected backend
	 * @param threads Requested workers; MT also requires the global MT flag
	 * @param append Whether to append rather than truncate the file
	 * @throws IOException If the destination cannot be opened
	 */
	public BamOutputStream(String filename, int compression, int threads, boolean append) throws IOException{
		this(new FileOutputStream(filename, append), true, compression, threads);
		suppressHeaderEmit=append && appendTargetNonEmpty(filename);
	}

	/**
	 * Wraps a stream using the configured BGZF compression and worker defaults.
	 * @param out Destination receiving compressed BAM bytes
	 * @param closeUnderlying Requested ownership policy; currently bypassed by backend close
	 */
	public BamOutputStream(OutputStream out, boolean closeUnderlying){
		this(out, closeUnderlying, BgzfSettings.WRITE_COMPRESSION_LEVEL, BgzfSettings.WRITE_THREADS);
	}

	/**
	 * Selects an MT backend only when enabled globally and threads exceeds one.
	 * Otherwise selects the ST backend. The MT block size is clamped to at least one.
	 * No BAM header is written until an alignment is processed or this facade closes.
	 * @param out Destination receiving compressed BAM bytes
	 * @param closeUnderlying Requested ownership policy; currently bypassed by backend close
	 * @param compression Compression level passed to the selected backend
	 * @param threads Requested MT workers
	 */
	public BamOutputStream(OutputStream out, boolean closeUnderlying, int compression, int threads){
		this.closeUnderlying=closeUnderlying;
		if(BgzfSettings.USE_MULTITHREADED_BGZF && threads>1){
			int blockSize=Math.max(1, BgzfSettings.WRITE_BLOCK_SIZE);
			mtOut=new BgzfOutputStreamMT(out, threads, compression, blockSize);
			stOut=null;
			bgzf=mtOut;
//			System.err.println("MT: compression="+compression+", threads="+threads+", blockSize="+blockSize);
		}else{
			stOut=new BgzfOutputStream(out, compression);
			mtOut=null;
			bgzf=stOut;
//			System.err.println("ST: compression="+compression+", threads="+threads);
		}
		helper=new BamWriterHelper(bgzf);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Accepts one SAM byte; LF completes the pending line.
	 * Interior CR bytes are retained by this overload; trailing CR is trimmed later.
	 * @param b Byte value, narrowed when appended to the pending line
	 * @throws IOException If closed or a downstream write fails
	 */
	@Override
	public void write(int b) throws IOException{
		ensureOpen();
		if(b=='\n'){flushLine();}else{lineBuffer.append((byte)b);}
	}

	/**
	 * Accepts SAM bytes, completing lines at LF and discarding all CR bytes.
	 * The caller supplies a valid array range; this loop does not prevalidate it.
	 * @param b Source bytes
	 * @param off First byte to consume
	 * @param len Number of bytes to consume
	 * @throws IOException If closed or a downstream write fails
	 */
	@Override
	public void write(byte[] b, int off, int len) throws IOException{
		ensureOpen();
		int end=off+len;
		for(int i=off; i<end; i++){
			byte ch=b[i];
			if(ch=='\n'){flushLine();}else if(ch!='\r'){lineBuffer.append(ch);}
		}
	}

	/**
	 * Completes any pending line, drains alignment records, and flushes BGZF.
	 * Leading header lines alone remain buffered until an alignment or close.
	 * @throws IOException If closed or a downstream write fails
	 */
	@Override
	public void flush() throws IOException{
		ensureOpen();
		flushLine();
		flushSamBuffer();
		bgzf.flush();
	}

	/**
	 * Completes the pending line and records, emits a header if still needed,
	 * and closes the selected backend with its EOF marker after all pending data.
	 * The ST path flushes newly emitted headers and alignment bytes before its marker.
	 * Repeated calls return immediately. The closed flag is set even if an operation
	 * throws; closing again does not retry the unfinished work.
	 * @throws IOException If a downstream write, flush, or close fails
	 */
	@Override
	public void close() throws IOException{
		if(closed){return;}
		try{
			flushLine();
			flushSamBuffer();
			if(!headerWritten){writeHeader();}
			//Existing per-engine BGZF EOF rationale: the MT engine writes the 28-byte EOF
			//marker INSIDE its own close() (BgzfOutputStreamMT.writerLoop on lastJob), so we just close it;
			//the ST engine's close() deliberately does NOT write EOF, so flush its pending data,
			//writeEOF explicitly, THEN close — exactly once each way. Getting this wrong would either omit the EOF
			//(truncation-undetectable BAM) or double-write it. See BgzfOutputStream#close contract note.
			//TODO: Probable bug - stream/bam/BamOutputStream#001: these backend
			//branches bypass closeUnderlying. In the ST path, BgzfOutputStream.close()
			//calls out.close() even when this facade was given closeUnderlying=false.
			//Confirmed from source only; ownership repair remains separate.
			if(mtOut!=null){mtOut.close();}else if(stOut!=null){
				// Fixed STR208: writeEOF does not drain data, including a header just emitted above.
				stOut.flush();
				stOut.writeEOF();
				stOut.close();
			}else if(closeUnderlying){bgzf.close();}else{bgzf.flush();}
		}finally{
			closed=true;
		}
	}

	/**
	 * Consumes the pending line, trimming trailing CR/LF and ignoring empty lines.
	 * Accumulates leading headers; otherwise ensures the header, parses an alignment,
	 * and drains the SAM batch when it reaches the threshold.
	 * @throws IOException If emitting the header or a record fails
	 */
	private void flushLine() throws IOException{
		if(lineBuffer.isEmpty()){return;}
		byte[] raw=lineBuffer.toBytes();
		lineBuffer.clear();
		int len=raw.length;
		while(len>0 && (raw[len-1]=='\r' || raw[len-1]=='\n')){len--;}
		if(len==0){return;}
		byte[] line=(len==raw.length) ? raw : trimBytes(raw, len);
		//Header accumulation invariant: @-lines are buffered into headerLines UNTIL the first non-@ line,
		//which triggers writeHeader() (emits BAM magic + header text + ref dict + builds the converter) and
		//flips headerWritten. After that, headerWritten is permanent, so a stray @-line arriving post-header
		//(malformed interleaved SAM) would fall through and be parsed as an alignment — out of scope (valid
		//SAM puts all @-lines first). Header-only input (no alignment) is handled by close()'s
		//if(!headerWritten) writeHeader(). Empty headerLines (headerless SAM) → writeHeader builds a 0-ref
		//dict + converter. SamToBamConverter returns refID=-1 only for null/"*" names;
		//a named reference absent from that dictionary throws IllegalArgumentException.
		if(!headerWritten && line[0]=='@'){
			headerLines.add(line);
			return;
		}
		if(!headerWritten){writeHeader();}
		SamLine sam=new SamLine(new LineParser1('\t').set(line));
		samBuffer.add(sam);
		if(samBuffer.size()>=FLUSH_THRESHOLD){flushSamBuffer();}
	}

	/**
	 * Writes each buffered alignment body with its four-byte length prefix.
	 * Returns without work until the header has been handled or if no records wait.
	 * Clears the batch only after all its records have been written.
	 * @throws IOException If a record write fails
	 */
	private void flushSamBuffer() throws IOException{
		if(!headerWritten || samBuffer.isEmpty()){return;}
		for(SamLine sam : samBuffer){
			byte[] record=converter.convertAlignment(sam);
			helper.writeUint32(record.length);
			helper.writeBytes(record);
		}
		samBuffer.clear();
	}

	/**
	 * Collects complete SN/LN entries from leading SQ lines in encounter order.
	 * Emits BAM magic, text and dictionary unless append mode suppresses emission,
	 * and initializes the converter from those names. Marks the header handled last.
	 * Does not compare an append target's existing dictionary with these entries.
	 * @throws IOException If header output fails
	 */
	private void writeHeader() throws IOException{
		if(headerWritten){return;}
		if(!suppressHeaderEmit){
			helper.writeBytes(new byte[]{'B', 'A', 'M', 1});
		}

		StringBuilder textBuilder=new StringBuilder();
		List<String> refs=new ArrayList<>();
		List<Integer> refLengths=new ArrayList<>();
		for(byte[] line : headerLines){
			String lineStr=new String(line, StandardCharsets.US_ASCII);
			textBuilder.append(lineStr).append('\n');
			if(lineStr.startsWith("@SQ")){
				String[] fields=lineStr.split("\\t");
				String sn=null;
				Integer ln=null;
				for(int i=1; i<fields.length; i++){
					String field=fields[i];
					if(field.startsWith("SN:")){sn=field.substring(3);}else if(field.startsWith("LN:")){ln=Integer.parseInt(field.substring(3));}
				}
				if(sn!=null && ln!=null){
					refs.add(sn);
					refLengths.add(ln);
				}
			}
		}

		if(!suppressHeaderEmit){
			byte[] headerText=textBuilder.toString().getBytes(StandardCharsets.US_ASCII);
			helper.writeUint32(headerText.length);
			helper.writeBytes(headerText);

			helper.writeInt32(refs.size());
			for(int i=0; i<refs.size(); i++){
				byte[] nameBytes=refs.get(i).getBytes(StandardCharsets.US_ASCII);
				helper.writeUint32(nameBytes.length+1);
				helper.writeBytes(nameBytes);
				helper.writeUint8(0);
				helper.writeUint32(refLengths.get(i));
			}
		}

		if(converter==null){converter=new SamToBamConverter(refs.toArray(new String[0]));}

		headerWritten=true;
	}

	/** Returns whether the destination length is positive; returns false on any throwable. */
	private static boolean appendTargetNonEmpty(String filename){
		try{
			return new java.io.File(filename).length()>0;
		}catch(Throwable t){
			return false;
		}
	}

	/** @throws IOException If this facade has already been closed */
	private void ensureOpen() throws IOException{
		if(closed){throw new IOException("Stream closed");}
	}

	/** Copies the first len bytes into a new array; the caller supplies a valid prefix length. */
	private static byte[] trimBytes(byte[] data, int len){
		byte[] trimmed=new byte[len];
		System.arraycopy(data, 0, trimmed, 0, len);
		return trimmed;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of parsed alignments that triggers a batch write. */
	private static final int FLUSH_THRESHOLD=256;

	/** Requested ownership policy, currently bypassed by both selected-backend branches. */
	private final boolean closeUnderlying;
	/** Selected MT backend, or null when using ST. */
	private final BgzfOutputStreamMT mtOut;
	/** Selected ST backend, or null when using MT. */
	private final BgzfOutputStream stOut;
	/** Selected BGZF output backend. */
	private final OutputStream bgzf;
	/** Little-endian writer for the selected backend. */
	private final BamWriterHelper helper;

	/** Leading header lines retained for header text and reference dictionary construction. */
	private final ArrayList<byte[]> headerLines=new ArrayList<>();
	/** Parsed alignments awaiting conversion and output. */
	private final ArrayList<SamLine> samBuffer=new ArrayList<>();
	/** Pending SAM line, excluding consumed LF bytes. */
	private final ByteBuilder lineBuffer=new ByteBuilder(4096);

	/** Alignment converter initialized from the incoming reference names. */
	private SamToBamConverter converter;
	/** Whether the header was emitted or deliberately suppressed for append. */
	private boolean headerWritten=false;
	/** Whether close was attempted, including an unsuccessful attempt. */
	private boolean closed=false;
	/** Appending to an existing nonempty BAM: consume header lines for the ref dict but emit no
	 * magic/header/dict (a second header block mid-file is corrupt). */
	private boolean suppressHeaderEmit=false;
}
